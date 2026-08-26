/*
   LKina — AutoGrid4 scoring function engine
   Copyright (C) 2025 LK-Studio1128 and LKina contributors
   SPDX-License-Identifier: GPL-3.0-or-later

   This file implements the AutoGrid 4.2 pairwise interaction scoring
   function (vdW, H-bond, electrostatics, desolvation). It was written
   with reference to the AutoGrid4 source code (GPL-2.0 or later) and
   is therefore distributed under the GNU General Public License v3.

   Algorithm reference:
     Morris et al., J. Comput. Chem. 19:1639-1662 (1998)
     Morris et al., J. Comput. Chem. 30:2785-2791 (2009)
   AutoGrid4 upstream: https://github.com/ccsb-scripps/AutoDock4
   AutoGrid4 license:  GPL-2.0 or later

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program. If not, see <https://www.gnu.org/licenses/>.
*/

#include "ag4_engine.h"
#include "ad4_parameter_data.h"
#include "atom_constants.h"

std::string get_adtype_str(sz& t);

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <array>
#include <cassert>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <map>
#include <algorithm>
#ifdef _OPENMP
#  include <omp.h>
#endif

// ------------------------------------------------------------------
// Constants — AutoGrid 4.2 defaults (no USE_8A_NBCUTOFF)
// ------------------------------------------------------------------
static const double AG4_A_DIV     = 128.0;
static const double AG4_INV_ADIV  = 1.0 / 128.0;
static const int    AG4_NEINT     = 2048;      // NBC(16) * A_DIV(128)
static const int    AG4_NDIEL     = 8192;      // DIELCUT(64) * A_DIV(128)
static const double AG4_SOFTNBC   = 8.0;       // affinity cutoff Ang
static const double AG4_DIELCUT   = 64.0;      // elec/desolv cutoff Ang
static const double AG4_FACTOR    = 332.0;     // Coulomb kcal*Ang/e^2
static const double AG4_DIEL      = 4.0;       // constant dielectric
static const double AG4_SOLPAR_Q  = 0.01097;
static const double AG4_SIGMA     = 3.6;       // desolvation sigma Ang
// AutoGrid 4.2 clamps affinity-map values at ±MAXVALUE (±1000 kcal/mol).
// LKina originally used ±100000, which let unphysical repulsive spikes
// (e.g. a TZ pseudoatom sampled near r=0 on a 0.375 Å grid, measured
// +14630 kcal/mol) propagate into the affinity maps.  ±1000 matches
// AutoGrid4's grid files exactly and bounds all pair energies.
static const double AG4_EINTCLAMP = 1000.0;
static const double AG4_R_SMOOTH  = 0.5;       // smoothing half-width Ang
static const double AG4_APPROX0   = 1.0e-6;

inline int    ag4_lookup  (double r) { return (int)(r * AG4_A_DIV); }
inline double ag4_angstrom(int    i) { return i * AG4_INV_ADIV;     }

// ------------------------------------------------------------------
// Parsed AD4 atom parameters
// ------------------------------------------------------------------
struct ag4_atom_par {
    double Rij;        // vdW radius
    double epsij;      // vdW well depth  (unweighted — multiply by coeff_vdW)
    double vol;        // solvation volume Ang^3
    double solpar;     // solvation parameter
    double Rij_hb;     // H-bond radius
    double epsij_hb;   // H-bond well depth (unweighted — multiply by coeff_hbond)
    int    hbond;      // 0=NON,1=DS,2=D1,3=AS,4=A1,5=A2
    bool   valid;
    ag4_atom_par() : Rij(0),epsij(0),vol(0),solpar(0),
                     Rij_hb(0),epsij_hb(0),hbond(0),valid(false) {}
};

typedef std::map<std::string, ag4_atom_par> ag4_param_table;

// ------------------------------------------------------------------
// Receptor atom
// ------------------------------------------------------------------
struct ag4_rec_atom {
    double x, y, z, charge;
    std::string type_name;
    double vol, solpar;
    int    hbond;       // 0=NON,1=DS,2=D1,3=AS,4=A1,5=A2
    double rvector[3];  // H-bond direction (D1/A1/A2)
    double rvector2[3]; // A2 plane normal
    int    rexp;        // D1: exponent for racc (2=N-H, 4=O-H)
    bool   has_hb_vec;
    ag4_rec_atom()
        : x(0),y(0),z(0),charge(0),vol(0),solpar(0),hbond(0),rexp(2),has_hb_vec(false)
    { rvector[0]=rvector[1]=rvector[2]=rvector2[0]=rvector2[1]=rvector2[2]=0; }
};

// ------------------------------------------------------------------
// Parse builtin AD4 parameter text
// ------------------------------------------------------------------
static void ag4_parse_params(ag4_param_table& table,
                              double& coeff_vdW, double& coeff_hbond,
                              double& coeff_estat, double& coeff_desolv)
{
    coeff_vdW    = ad4_param::FE_coeff_vdW;
    coeff_hbond  = ad4_param::FE_coeff_hbond;
    coeff_estat  = ad4_param::FE_coeff_estat;
    coeff_desolv = ad4_param::FE_coeff_desolv;

    std::istringstream iss(ad4_param::builtin_ad4_parameter_text);
    std::string line;
    while (std::getline(iss, line)) {
        if (line.size() < 9 || line.compare(0, 9, "atom_par ") != 0) continue;
        char name[32] = {};
        double r=0,e=0,v=0,s=0,rh=0,eh=0;
        int hb=0,ri=0,mi=0,bi=0;
        int n = sscanf(line.c_str(),
            "atom_par %31s %lf %lf %lf %lf %lf %lf %d %d %d %d",
            name, &r, &e, &v, &s, &rh, &eh, &hb, &ri, &mi, &bi);
        if (n < 8) continue;
        ag4_atom_par p;
        p.Rij=r; p.epsij=e; p.vol=v; p.solpar=s;
        p.Rij_hb=rh; p.epsij_hb=eh; p.hbond=hb; p.valid=true;
        table[std::string(name)] = p;
    }
}

// ------------------------------------------------------------------
// Minimal receptor PDBQT parser
// ------------------------------------------------------------------
static void ag4_parse_receptor(const std::string& pdbqt_path,
                                const ag4_param_table& params,
                                std::vector<ag4_rec_atom>& atoms)
{
    std::ifstream ifs(pdbqt_path.c_str());
    if (!ifs) throw std::runtime_error("ag4_engine: cannot open receptor: " + pdbqt_path);
    std::string line;
    while (std::getline(ifs, line)) {
        if (line.size() < 6) continue;
        const std::string rec6 = line.substr(0, 6);
        if (rec6 != "ATOM  " && rec6 != "HETATM") continue;
        if (line.size() < 54) continue;
        // Pad to 80 characters for safe column access
        while (line.size() < 80) line.push_back(' ');

        ag4_rec_atom a;
        try {
            a.x = std::atof(line.substr(30, 8).c_str());
            a.y = std::atof(line.substr(38, 8).c_str());
            a.z = std::atof(line.substr(46, 8).c_str());
        } catch (...) { continue; }

        // Partial charge at cols 69-76 (0-indexed)
        {
            std::string cs = line.substr(68, 8);
            a.charge = std::atof(cs.c_str());
        }

        // AD4 type at cols 78-79 (0-indexed 77-78, up to 2 chars)
        {
            std::string ts = line.substr(77, 2);
            while (!ts.empty() && ts[0]  == ' ') ts.erase(0, 1);
            while (!ts.empty() && ts.back()== ' ') ts.erase(ts.size()-1, 1);
            if (ts.empty()) continue;
            a.type_name = ts;
        }

        // Apply same equivalences as atom_equivalence_data in atom_constants.h
        if (a.type_name == "Cr") a.type_name = "Cr1";
        if (a.type_name == "W")  a.type_name = "Tg"; // Tungsten element; "W" alone = water probe (skip)

        ag4_param_table::const_iterator it = params.find(a.type_name);
        if (it == params.end()) continue;   // unknown type — skip

        const ag4_atom_par& par = it->second;
        a.vol    = par.vol;
        a.solpar = par.solpar;
        a.hbond  = par.hbond;
        atoms.push_back(a);
    }
}

static bool ag4_is_donor_type_name(const std::string& tn) {
    return tn == "NA" || tn == "N" || tn == "OA" || tn == "O" || tn == "SA" || tn == "S";
}

static const double AG4_PSEUDO_COS_OCC = 0.766;

static const double AG4_PSEUDO_CANDS[26][3] = {
    { 1, 0, 0},{-1, 0, 0},{ 0, 1, 0},{ 0,-1, 0},{ 0, 0, 1},{ 0, 0,-1},
    { 0.70710678, 0.70710678, 0},{ 0.70710678,-0.70710678, 0},{-0.70710678, 0.70710678, 0},{-0.70710678,-0.70710678, 0},
    { 0.70710678, 0, 0.70710678},{ 0.70710678, 0,-0.70710678},{-0.70710678, 0, 0.70710678},{-0.70710678, 0,-0.70710678},
    { 0, 0.70710678, 0.70710678},{ 0, 0.70710678,-0.70710678},{ 0,-0.70710678, 0.70710678},{ 0,-0.70710678,-0.70710678},
    { 0.57735027, 0.57735027, 0.57735027},{ 0.57735027, 0.57735027,-0.57735027},{ 0.57735027,-0.57735027, 0.57735027},{ 0.57735027,-0.57735027,-0.57735027},
    {-0.57735027, 0.57735027, 0.57735027},{-0.57735027, 0.57735027,-0.57735027},{-0.57735027,-0.57735027, 0.57735027},{-0.57735027,-0.57735027,-0.57735027}
};

static std::vector<ag4_metal_mode> ag4_collect_active_modes(ag4_metal_mode primary_mode,
                                                             const std::vector<ag4_metal_mode>& extra_metal_modes) {
    std::vector<ag4_metal_mode> out;
    if (primary_mode != ag4_metal_mode::none)
        out.push_back(primary_mode);
    for (ag4_metal_mode mm : extra_metal_modes) {
        if (mm == ag4_metal_mode::none) continue;
        if (std::find(out.begin(), out.end(), mm) == out.end())
            out.push_back(mm);
    }
    return out;
}

static std::vector<std::array<double,3>> ag4_collect_coord_dirs(const std::vector<ag4_rec_atom>& rec,
                                                                int metal_index,
                                                                double coord_dist,
                                                                const std::vector<std::string>& pseudo_types) {
    std::vector<std::array<double,3>> occ;
    const ag4_rec_atom& M = rec[metal_index];
    double dmax2 = coord_dist * coord_dist;
    for (int i = 0; i < (int)rec.size(); i++) {
        if (i == metal_index) continue;
        bool accept = ag4_is_donor_type_name(rec[i].type_name);
        if (!accept) {
            for (const std::string& pt : pseudo_types)
                if (rec[i].type_name == pt) { accept = true; break; }
        }
        if (!accept) continue;
        double dx = rec[i].x - M.x;
        double dy = rec[i].y - M.y;
        double dz = rec[i].z - M.z;
        double d2 = dx*dx + dy*dy + dz*dz;
        if (d2 < dmax2 && d2 > 0.01) {
            double len = std::sqrt(d2);
            occ.push_back({dx/len, dy/len, dz/len});
        }
    }
    return occ;
}

static std::vector<std::array<double,3>> ag4_select_vacant_dirs(const std::vector<std::array<double,3>>& occ,
                                                                 int vacancies) {
    std::vector<std::array<double,3>> placed;
    for (int v = 0; v < vacancies; v++) {
        double best_score = 2.0;
        int best_k = -1;
        for (int k = 0; k < 26; k++) {
            double max_sim = -2.0;
            bool blocked = false;
            for (const auto& od : occ) {
                double c = AG4_PSEUDO_CANDS[k][0]*od[0] + AG4_PSEUDO_CANDS[k][1]*od[1] + AG4_PSEUDO_CANDS[k][2]*od[2];
                if (c > AG4_PSEUDO_COS_OCC) { blocked = true; break; }
                if (c > max_sim) max_sim = c;
            }
            if (blocked) continue;
            for (const auto& pd : placed) {
                double c = AG4_PSEUDO_CANDS[k][0]*pd[0] + AG4_PSEUDO_CANDS[k][1]*pd[1] + AG4_PSEUDO_CANDS[k][2]*pd[2];
                if (c > AG4_PSEUDO_COS_OCC) { blocked = true; break; }
                if (c > max_sim) max_sim = c;
            }
            if (blocked) continue;
            if (max_sim < best_score) { best_score = max_sim; best_k = k; }
        }
        if (best_k < 0) break;
        placed.push_back({AG4_PSEUDO_CANDS[best_k][0], AG4_PSEUDO_CANDS[best_k][1], AG4_PSEUDO_CANDS[best_k][2]});
    }
    return placed;
}

static std::array<double,3> ag4_tetragonal_axis(const std::vector<std::array<double,3>>& occ) {
    if (occ.size() >= 2) {
        double best_dot = 1.0;
        int bi = 0, bj = 1;
        for (int i = 0; i < (int)occ.size(); i++)
            for (int j = i + 1; j < (int)occ.size(); j++) {
                double dot = occ[i][0]*occ[j][0] + occ[i][1]*occ[j][1] + occ[i][2]*occ[j][2];
                if (dot < best_dot) { best_dot = dot; bi = i; bj = j; }
            }
        double ax = occ[bi][0] - occ[bj][0];
        double ay = occ[bi][1] - occ[bj][1];
        double az = occ[bi][2] - occ[bj][2];
        double len = std::sqrt(ax*ax + ay*ay + az*az);
        if (len > 1e-6) return {ax/len, ay/len, az/len};
    }
    if (!occ.empty()) {
        double ax = occ[0][0], ay = occ[0][1], az = occ[0][2];
        double len = std::sqrt(ax*ax + ay*ay + az*az);
        if (len > 1e-6) return {ax/len, ay/len, az/len};
    }
    return {0.0, 0.0, 1.0};
}

static std::vector<std::string> ag4_mode_metal_names(ag4_metal_mode mode) {
    if (ag4_needs_jt_injection(mode)) return ag4_jt_injection_params(mode).metal_names;
    if (ag4_needs_sq_injection(mode)) return ag4_sq_injection_params(mode).metal_names;
    if (ag4_needs_mh_injection(mode)) return ag4_mh_injection_params(mode).metal_names;
    if (ag4_needs_tz_injection(mode)) return ag4_tz_injection_params(mode).metal_names;
    switch (mode) {
        case ag4_metal_mode::zn:
            return {"Zn"};
        case ag4_metal_mode::fe:
        case ag4_metal_mode::fe2:
        case ag4_metal_mode::fe3:
        case ag4_metal_mode::fe3_aq:
            return {"Fe"};
        case ag4_metal_mode::mn:
        case ag4_metal_mode::mn2:
        case ag4_metal_mode::mn3:
        case ag4_metal_mode::mn2_aq:
            return {"Mn"};
        case ag4_metal_mode::co:
        case ag4_metal_mode::co2:
        case ag4_metal_mode::co3:
        case ag4_metal_mode::co2_aq:
            return {"Co"};
        case ag4_metal_mode::mg:
        case ag4_metal_mode::mg_aq:
            return {"Mg"};
        case ag4_metal_mode::ca:
        case ag4_metal_mode::ca_aq:
            return {"Ca"};
        case ag4_metal_mode::cu:
        case ag4_metal_mode::cu1:
        case ag4_metal_mode::cu2:
            return {"Cu"};
        case ag4_metal_mode::ru: return {"Ru"};
        case ag4_metal_mode::ir: return {"Ir"};
        case ag4_metal_mode::rh: return {"Rh"};
        case ag4_metal_mode::pt: return {"Pt"};
        case ag4_metal_mode::pd: return {"Pd"};
        case ag4_metal_mode::ni: return {"Ni"};
        case ag4_metal_mode::au: return {"Au"};
        case ag4_metal_mode::na_ion: return {"Na"};
        case ag4_metal_mode::k_ion: return {"K"};
        default: return {};
    }
}

// ------------------------------------------------------------------
// AutoDock4Zn-style: inject TZ pseudoatoms at vacant tetrahedral coord sites
// (Santos-Martins et al. 2014 — zinc_pseudo.py geometry)
// Now generalised: metal_names may include any d10 tetrahedral metal (Zn, Cd, …)
// ------------------------------------------------------------------
static void ag4_inject_tz_pseudoatoms(std::vector<ag4_rec_atom>& rec,
                                       const ag4_param_table& params,
                                       const std::vector<std::string>& metal_names = std::vector<std::string>(1,"Zn"),
                                       double tz_dist   = 2.0,
                                       int    max_coord = 4)
{
    const double ZN_COORD_DIST  = 2.5;  // max coordination radius (Å)
    const double TZ_DIST        = tz_dist;
    const int    MAX_COORD      = max_coord;

    // Locate TZ params (for vol/solpar/hbond fields)
    ag4_param_table::const_iterator tz_par = params.find("TZ");

    // Collect target metal atom indices (snapshot before appending)
    std::vector<int> zn_idx;
    for (int i = 0; i < (int)rec.size(); i++) {
        for (size_t mi = 0; mi < metal_names.size(); mi++)
            if (rec[i].type_name == metal_names[mi]) { zn_idx.push_back(i); break; }
    }

    std::vector<ag4_rec_atom> new_tz;
    for (int zi : zn_idx) {
        const ag4_rec_atom& zn = rec[zi];

        // Find receptor coordination partners within ZN_COORD_DIST
        std::vector<std::array<double,3>> occ_dirs; // occupied unit vectors Zn→partner
        for (int i = 0; i < (int)rec.size(); i++) {
            if (i == zi) continue;
            double dx = rec[i].x - zn.x;
            double dy = rec[i].y - zn.y;
            double dz = rec[i].z - zn.z;
            double d2 = dx*dx + dy*dy + dz*dz;
            if (d2 < ZN_COORD_DIST * ZN_COORD_DIST && d2 > 0.01) {
                double len = std::sqrt(d2);
                occ_dirs.push_back({dx/len, dy/len, dz/len});
            }
        }

        int n_occ = (int)occ_dirs.size();
        if (n_occ >= MAX_COORD) continue;  // already saturated — skip

        int n_tz = MAX_COORD - n_occ;

        // Iteratively place TZ atoms: each new direction is opposite to the
        // centroid of ALL current directions (occupied + already-placed TZ).
        // This approximates the ideal tetrahedral complement geometry.
        std::vector<std::array<double,3>> placed;
        for (int t = 0; t < n_tz; t++) {
            // Collect all current directions
            std::vector<std::array<double,3>> all = occ_dirs;
            all.insert(all.end(), placed.begin(), placed.end());

            double cx = 0, cy = 0, cz = 0;
            for (auto& d : all) { cx -= d[0]; cy -= d[1]; cz -= d[2]; }

            double clen = std::sqrt(cx*cx + cy*cy + cz*cz);
            if (clen < 1e-6) {
                // Degenerate (e.g. bare Zn with no receptor coords): use arbitrary dir
                if (!placed.empty()) {
                    double ax = placed.back()[0], ay = placed.back()[1], az = placed.back()[2];
                    cx = -ay; cy = ax; cz = 0.0;
                    clen = std::sqrt(cx*cx + cy*cy);
                    if (clen < 1e-6) { cx = az; cy = 0; cz = -ax; clen = std::sqrt(cx*cx + cz*cz); }
                } else {
                    cx = 1.0; cy = 0.0; cz = 0.0;
                }
                clen = std::sqrt(cx*cx + cy*cy + cz*cz);
                if (clen < 1e-6) clen = 1.0;
            }

            std::array<double,3> dir = {cx/clen, cy/clen, cz/clen};
            placed.push_back(dir);

            ag4_rec_atom tz;
            tz.x         = zn.x + TZ_DIST * dir[0];
            tz.y         = zn.y + TZ_DIST * dir[1];
            tz.z         = zn.z + TZ_DIST * dir[2];
            tz.charge    = 0.0;
            tz.type_name = "TZ";
            tz.hbond     = 0;
            tz.has_hb_vec = false;
            if (tz_par != params.end()) {
                tz.vol    = tz_par->second.vol;
                tz.solpar = tz_par->second.solpar;
                tz.hbond  = tz_par->second.hbond;
            }
            new_tz.push_back(tz);
        }
    }

    rec.insert(rec.end(), new_tz.begin(), new_tz.end());
}

// ------------------------------------------------------------------
// O1: Generic 3D greedy pseudoatom injection (octahedral, etc.)
// Places `pseudo_type` atoms at the `vacancies` most-separated unoccupied
// directions around each target metal atom.
// Algorithm: greedy max-min-separation on a 26-direction sphere sample.
// References: GPDOCK (Wang et al. 2023, Briefings Bioinformatics bbac620);
//             MBD (Clemente et al. 2024, JCIM acs.jcim.3c01853).
// ------------------------------------------------------------------
static void ag4_inject_3d_pseudoatoms(
        std::vector<ag4_rec_atom>&       rec,
        const ag4_param_table&           params,
        const std::vector<std::string>&  metal_names,
        const std::string&               pseudo_type,
        double                           coord_dist,
        int                              max_coord)
{
    const double COORD_DIST = 3.4;   // donor search radius (Å)

    ag4_param_table::const_iterator par = params.find(pseudo_type);

    std::vector<int> metal_idx;
    for (int i = 0; i < (int)rec.size(); i++)
        for (const auto& mn : metal_names)
            if (rec[i].type_name == mn) { metal_idx.push_back(i); break; }

    std::vector<ag4_rec_atom> new_atoms;
    for (int mi : metal_idx) {
        const ag4_rec_atom& M = rec[mi];

        std::vector<std::array<double,3>> occ = ag4_collect_coord_dirs(rec, mi, COORD_DIST, {pseudo_type});
        if ((int)occ.size() >= max_coord) continue;

        int vacancies = max_coord - (int)occ.size();
        std::vector<std::array<double,3>> placed = ag4_select_vacant_dirs(occ, vacancies);

        for (const auto& dir : placed) {
            ag4_rec_atom pa;
            pa.x=M.x+coord_dist*dir[0]; pa.y=M.y+coord_dist*dir[1]; pa.z=M.z+coord_dist*dir[2];
            pa.charge=0.0; pa.type_name=pseudo_type; pa.hbond=0; pa.has_hb_vec=false;
            if (par != params.end()) {
                pa.vol=par->second.vol; pa.solpar=par->second.solpar;
                pa.hbond=par->second.hbond;
            }
            new_atoms.push_back(pa);
        }
    }
    rec.insert(rec.end(), new_atoms.begin(), new_atoms.end());
}

// ------------------------------------------------------------------
// M4: Inject SQ pseudoatoms at vacant sq-planar / linear coord sites
// (Pt, Pd, Ni, Cu — square planar d8/d9; Au — linear d10)
// Algorithm: determine coordination plane from existing N/O/S donors;
// place SQ at equally-spaced in-plane directions that are unoccupied.
// References: MetalDock 2024 (PMC10751784); analogous to TZ approach.
// ------------------------------------------------------------------
static void ag4_inject_sq_pseudoatoms(
        std::vector<ag4_rec_atom>&       rec,
        const ag4_param_table&           params,
        const std::vector<std::string>&  metal_names,
        double sq_dist,   // ideal coordination bond length (Å)
        int    max_coord) // 4 = square planar, 2 = linear
{
    const double COORD_DIST = 3.0;   // coordination partner search radius (Å)
    const double COS_THR    = 0.766; // within 40° = already occupied

    ag4_param_table::const_iterator sq_par = params.find("SQ");

    // Collect target metal atom indices (snapshot before appending)
    std::vector<int> metal_idx;
    for (int i = 0; i < (int)rec.size(); i++)
        for (const auto& mn : metal_names)
            if (rec[i].type_name == mn) { metal_idx.push_back(i); break; }

    std::vector<ag4_rec_atom> new_sq;
    for (int mi : metal_idx) {
        const ag4_rec_atom& M = rec[mi];

        std::vector<std::array<double,3>> occ = ag4_collect_coord_dirs(rec, mi, COORD_DIST, {"SQ"});
        int n_occ = (int)occ.size();
        if (n_occ >= max_coord) continue; // already saturated

        // ── Determine coordination-plane normal n ──────────────────────────
        double nx=0, ny=0, nz=1; // default: arbitrary z-axis
        if (n_occ >= 2) {
            nx = occ[0][1]*occ[1][2] - occ[0][2]*occ[1][1];
            ny = occ[0][2]*occ[1][0] - occ[0][0]*occ[1][2];
            nz = occ[0][0]*occ[1][1] - occ[0][1]*occ[1][0];
            double nlen = std::sqrt(nx*nx+ny*ny+nz*nz);
            if (nlen < 1e-6) { // dirs are parallel (trans pair) — perp to first
                nx=-occ[0][1]; ny=occ[0][0]; nz=0.0;
                double nl2=std::sqrt(nx*nx+ny*ny);
                if (nl2<1e-6) { nx=0.0; ny=-occ[0][2]; nz=occ[0][1]; nl2=std::sqrt(ny*ny+nz*nz); }
                nlen=nl2; if (nlen<1e-6) nlen=1.0;
            }
            nx/=nlen; ny/=nlen; nz/=nlen;
        } else if (n_occ == 1) {
            nx=-occ[0][1]; ny=occ[0][0]; nz=0.0;
            double nlen=std::sqrt(nx*nx+ny*ny);
            if (nlen<1e-6) { nx=0.0; ny=-occ[0][2]; nz=occ[0][1]; nlen=std::sqrt(ny*ny+nz*nz); }
            if (nlen<1e-6) nlen=1.0;
            nx/=nlen; ny/=nlen; nz/=nlen;
        }

        // ── Build in-plane orthonormal basis u1, u2 ────────────────────────
        double u1x, u1y, u1z;
        if (n_occ >= 1) {
            // Project first occupied dir onto plane → u1
            double dot = occ[0][0]*nx + occ[0][1]*ny + occ[0][2]*nz;
            u1x=occ[0][0]-dot*nx; u1y=occ[0][1]-dot*ny; u1z=occ[0][2]-dot*nz;
        } else {
            // Arbitrary in-plane vector
            u1x=(std::abs(nx)<0.9)?1.0:0.0;
            u1y=(u1x<0.5&&std::abs(ny)<0.9)?1.0:0.0;
            u1z=(u1x<0.5&&u1y<0.5)?1.0:0.0;
            if (u1x<0.5&&u1y<0.5&&u1z<0.5) u1x=1.0;
            double dot=u1x*nx+u1y*ny+u1z*nz;
            u1x-=dot*nx; u1y-=dot*ny; u1z-=dot*nz;
        }
        double u1len=std::sqrt(u1x*u1x+u1y*u1y+u1z*u1z);
        if (u1len<1e-6) { u1x=1.0; u1y=0.0; u1z=0.0; }
        else             { u1x/=u1len; u1y/=u1len; u1z/=u1len; }
        // u2 = n × u1
        double u2x=ny*u1z-nz*u1y, u2y=nz*u1x-nx*u1z, u2z=nx*u1y-ny*u1x;

        // ── Inject SQ at vacant positions ──────────────────────────────────
        int n_sq = max_coord - n_occ;
        int injected = 0;
        for (int k = 0; k < max_coord && injected < n_sq; k++) {
            double theta = k * (2.0 * M_PI / max_coord);
            double cx = std::cos(theta)*u1x + std::sin(theta)*u2x;
            double cy = std::cos(theta)*u1y + std::sin(theta)*u2y;
            double cz = std::cos(theta)*u1z + std::sin(theta)*u2z;
            // Skip if within 40° of an occupied direction
            bool busy = false;
            for (auto& od : occ)
                if (cx*od[0]+cy*od[1]+cz*od[2] > COS_THR) { busy=true; break; }
            if (busy) continue;

            ag4_rec_atom sq;
            sq.x=M.x+sq_dist*cx; sq.y=M.y+sq_dist*cy; sq.z=M.z+sq_dist*cz;
            sq.charge=0.0; sq.type_name="SQ"; sq.hbond=0; sq.has_hb_vec=false;
            if (sq_par != params.end()) {
                sq.vol=sq_par->second.vol; sq.solpar=sq_par->second.solpar;
                sq.hbond=sq_par->second.hbond;
            }
            new_sq.push_back(sq);
            injected++;
        }
    }
    rec.insert(rec.end(), new_sq.begin(), new_sq.end());
}

static void ag4_inject_jt_pseudoatoms(
        std::vector<ag4_rec_atom>&       rec,
        const ag4_param_table&           params,
        const ag4_jt_injection_params_t& jp)
{
    if (jp.metal_names.empty()) return;

    ag4_param_table::const_iterator eq_par = params.find(jp.equatorial_type);
    ag4_param_table::const_iterator ax_par = params.find(jp.axial_type);
    std::vector<int> metal_idx;
    for (int i = 0; i < (int)rec.size(); i++)
        for (const auto& mn : jp.metal_names)
            if (rec[i].type_name == mn) { metal_idx.push_back(i); break; }

    std::vector<ag4_rec_atom> new_atoms;
    for (int mi : metal_idx) {
        const ag4_rec_atom& M = rec[mi];
        std::vector<std::array<double,3>> occ = ag4_collect_coord_dirs(rec, mi, 3.4, {jp.equatorial_type, jp.axial_type});
        std::array<double,3> axis = ag4_tetragonal_axis(occ);

        std::vector<std::array<double,3>> eq_occ;
        int axial_occ = 0;
        for (const auto& od : occ) {
            double dot = std::fabs(od[0]*axis[0] + od[1]*axis[1] + od[2]*axis[2]);
            if (dot > 0.75) axial_occ++;
            else eq_occ.push_back(od);
        }

        int eq_vac = std::max(0, jp.equatorial_coord - (int)eq_occ.size());
        std::vector<std::array<double,3>> eq_dirs = ag4_select_vacant_dirs(eq_occ, eq_vac);
        for (const auto& dir : eq_dirs) {
            if (std::fabs(dir[0]*axis[0] + dir[1]*axis[1] + dir[2]*axis[2]) > 0.65) continue;
            ag4_rec_atom pa;
            pa.x=M.x+jp.equatorial_dist*dir[0]; pa.y=M.y+jp.equatorial_dist*dir[1]; pa.z=M.z+jp.equatorial_dist*dir[2];
            pa.charge=0.0; pa.type_name=jp.equatorial_type; pa.hbond=0; pa.has_hb_vec=false;
            if (eq_par != params.end()) { pa.vol=eq_par->second.vol; pa.solpar=eq_par->second.solpar; pa.hbond=eq_par->second.hbond; }
            new_atoms.push_back(pa);
        }

        if (axial_occ < jp.axial_coord) {
            const double dirs[2][3] = {{ axis[0], axis[1], axis[2]}, {-axis[0], -axis[1], -axis[2]}};
            for (int s = 0; s < 2; s++) {
                bool busy = false;
                for (const auto& od : occ) {
                    double dot = dirs[s][0]*od[0] + dirs[s][1]*od[1] + dirs[s][2]*od[2];
                    if (dot > 0.75) { busy = true; break; }
                }
                if (busy) continue;
                ag4_rec_atom pa;
                pa.x=M.x+jp.axial_dist*dirs[s][0]; pa.y=M.y+jp.axial_dist*dirs[s][1]; pa.z=M.z+jp.axial_dist*dirs[s][2];
                pa.charge=0.0; pa.type_name=jp.axial_type; pa.hbond=0; pa.has_hb_vec=false;
                if (ax_par != params.end()) { pa.vol=ax_par->second.vol; pa.solpar=ax_par->second.solpar; pa.hbond=ax_par->second.hbond; }
                new_atoms.push_back(pa);
            }
        }
    }
    rec.insert(rec.end(), new_atoms.begin(), new_atoms.end());
}

// ------------------------------------------------------------------
// H-bond geometry setup for D1 (hbond=2), A1 (hbond=4), A2 (hbond=5)
// ------------------------------------------------------------------
static void ag4_setup_hbond_geometry(std::vector<ag4_rec_atom>& atoms)
{
    int N = (int)atoms.size();
    for (int ia = 0; ia < N; ia++) {
        ag4_rec_atom& a = atoms[ia];
        if (a.hbond == 0) continue;

        int from = std::max(0,   ia - 20);
        int to   = std::min(N-1, ia + 20);

        // ---- D1: hydrogen bond donor (ia is H) ----
        if (a.hbond == 2) {
            double best_d2 = 1.378 * 1.378 + 0.01;
            int best_ib = -1;
            for (int ib = from; ib <= to; ib++) {
                if (ib == ia) continue;
                double dx = atoms[ib].x - a.x;
                double dy = atoms[ib].y - a.y;
                double dz = atoms[ib].z - a.z;
                double d2 = dx*dx + dy*dy + dz*dz;
                if (d2 <= best_d2) { best_d2 = d2; best_ib = ib; }
            }
            if (best_ib >= 0) {
                // rvector = normalized vector from bonded heavy atom -> H
                double dx = a.x - atoms[best_ib].x;
                double dy = a.y - atoms[best_ib].y;
                double dz = a.z - atoms[best_ib].z;
                double len = std::sqrt(dx*dx+dy*dy+dz*dz);
                if (len > AG4_APPROX0) {
                    a.rvector[0]=dx/len; a.rvector[1]=dy/len; a.rvector[2]=dz/len;
                    a.has_hb_vec = true;
                }
                // rexp: 4 for O-H/S-H, 2 for N-H
                const std::string& bt = atoms[best_ib].type_name;
                if (bt=="OA"||bt=="OS"||bt=="OW"||bt=="SA"||bt=="SS") a.rexp=4;
                else a.rexp=2;
            }
            continue;
        }

        // ---- A2: oxygen acceptor (2 lone pairs) ----
        if (a.hbond == 5) {
            double bvx=0, bvy=0, bvz=0;
            struct BAtom { double dx,dy,dz; };
            std::vector<BAtom> bonds;
            for (int ib = from; ib <= to && (int)bonds.size() < 2; ib++) {
                if (ib == ia) continue;
                double dx = atoms[ib].x - a.x;
                double dy = atoms[ib].y - a.y;
                double dz = atoms[ib].z - a.z;
                double d2 = dx*dx+dy*dy+dz*dz;
                bool isH = (atoms[ib].type_name=="HD"||atoms[ib].type_name=="H");
                double cut = isH ? 1.3*1.3 : 1.9*1.9;
                if (d2 <= cut) { BAtom b={dx,dy,dz}; bonds.push_back(b); }
            }
            if (!bonds.empty()) {
                // rvector = average direction of bonded atoms (points toward bonded)
                for (size_t k=0; k<bonds.size(); k++) { bvx+=bonds[k].dx; bvy+=bonds[k].dy; bvz+=bonds[k].dz; }
                double len = std::sqrt(bvx*bvx+bvy*bvy+bvz*bvz);
                if (len > AG4_APPROX0) {
                    a.rvector[0]=bvx/len; a.rvector[1]=bvy/len; a.rvector[2]=bvz/len;
                }
                if (bonds.size() >= 2) {
                    // rvector2 = normal to the plane spanned by the two bonds
                    double ax=bonds[0].dx, ay=bonds[0].dy, az=bonds[0].dz;
                    double bx=bonds[1].dx, by=bonds[1].dy, bz=bonds[1].dz;
                    double nx=ay*bz-az*by, ny=az*bx-ax*bz, nz=ax*by-ay*bx;
                    double nlen=std::sqrt(nx*nx+ny*ny+nz*nz);
                    if (nlen > AG4_APPROX0) {
                        a.rvector2[0]=nx/nlen; a.rvector2[1]=ny/nlen; a.rvector2[2]=nz/nlen;
                    }
                }
                a.has_hb_vec = true;
            }
            continue;
        }

        // ---- A1: planar N acceptor ----
        if (a.hbond == 4) {
            double vx=0, vy=0, vz=0; int cnt=0;
            for (int ib = from; ib <= to; ib++) {
                if (ib == ia) continue;
                double dx = atoms[ib].x - a.x;
                double dy = atoms[ib].y - a.y;
                double dz = atoms[ib].z - a.z;
                double d2 = dx*dx+dy*dy+dz*dz;
                if (d2 <= 1.9*1.9) { vx+=dx; vy+=dy; vz+=dz; cnt++; }
            }
            if (cnt > 0) {
                vx/=cnt; vy/=cnt; vz/=cnt;
                double len=std::sqrt(vx*vx+vy*vy+vz*vz);
                if (len > AG4_APPROX0) {
                    a.rvector[0]=vx/len; a.rvector[1]=vy/len; a.rvector[2]=vz/len;
                    a.has_hb_vec=true;
                }
            }
        }
    }
}

// ------------------------------------------------------------------
// Build desolvation lookup table
// ------------------------------------------------------------------
static void ag4_build_sol_fn(double coeff_desolv, std::vector<double>& sol_fn)
{
    sol_fn.resize(AG4_NDIEL, 0.0);
    double s2 = AG4_SIGMA * AG4_SIGMA;
    for (int i = 1; i < AG4_NDIEL; i++) {
        double r = ag4_angstrom(i);
        sol_fn[i] = coeff_desolv * std::exp(-(r*r) / (2.0*s2));
    }
}

// ------------------------------------------------------------------
// Build LJ energy lookup table for one probe/receptor type pair
// and apply window-smoothing (min over ±i_smooth bins)
// ------------------------------------------------------------------
static void ag4_build_lj_table(double nbp_r, double nbp_eps,
                                int xA, int xB,
                                std::vector<double>& table)
{
    table.assign(AG4_NEINT, 0.0);
    if (nbp_eps <= 0.0 || nbp_r <= 0.0 || xA == xB) {
        // degenerate: leave all zeros
        return;
    }
    double tmpconst = nbp_eps / (double)(xA - xB);
    double cA = tmpconst * std::pow(nbp_r, (double)xA) * xB;
    double cB = tmpconst * std::pow(nbp_r, (double)xB) * xA;
    for (int i = 1; i < AG4_NEINT; i++) {
        double r  = ag4_angstrom(i);
        double rA = std::pow(r, (double)xA);
        double rB = std::pow(r, (double)xB);
        table[i] = std::max(-AG4_EINTCLAMP,
                            std::min(AG4_EINTCLAMP, cA/rA - cB/rB));
    }
    table[0] = AG4_EINTCLAMP;
    table[AG4_NEINT-1] = 0.0;

    // Smoothing: replace each bin with minimum over ±i_smooth neighbourhood
    int i_smooth = (int)(AG4_R_SMOOTH * AG4_A_DIV / 2.0);
    if (i_smooth > 0) {
        std::vector<double> s(AG4_NEINT, AG4_EINTCLAMP);
        for (int i = 0; i < AG4_NEINT; i++) {
            int lo = std::max(0,          i - i_smooth);
            int hi = std::min(AG4_NEINT,  i + i_smooth + 1);
            for (int j = lo; j < hi; j++) s[i] = std::min(s[i], table[j]);
        }
        table.swap(s);
    }
}

// ------------------------------------------------------------------
// Utility: squared distance
// ------------------------------------------------------------------
inline double ag4_dist2(const ag4_rec_atom& a,
                         double gx, double gy, double gz)
{
    double dx=a.x-gx, dy=a.y-gy, dz=a.z-gz;
    return dx*dx + dy*dy + dz*dz;
}

ag4_metal_state ag4_build_metal_state(const std::string& receptor_pdbqt_path,
                                      ag4_metal_mode primary_mode,
                                      const std::vector<ag4_metal_mode>& extra_metal_modes)
{
    ag4_metal_state state;
    ag4_param_table params;
    double coeff_vdW, coeff_hbond, coeff_estat, coeff_desolv;
    ag4_parse_params(params, coeff_vdW, coeff_hbond, coeff_estat, coeff_desolv);

    std::vector<ag4_rec_atom> rec;
    ag4_parse_receptor(receptor_pdbqt_path, params, rec);
    std::vector<ag4_metal_mode> active_modes = ag4_collect_active_modes(primary_mode, extra_metal_modes);
    for (ag4_metal_mode mode : active_modes) {
        std::vector<std::string> metal_names = ag4_mode_metal_names(mode);
        if (metal_names.empty()) continue;

        double coord_dist = 2.1;
        int max_coord = 6;
        bool jt_enabled = false;
        double jt_axial_dist = 0.0;
        if (ag4_needs_jt_injection(mode)) {
            ag4_jt_injection_params_t jp = ag4_jt_injection_params(mode);
            coord_dist = jp.equatorial_dist;
            max_coord = jp.equatorial_coord + jp.axial_coord;
            jt_enabled = true;
            jt_axial_dist = jp.axial_dist;
        } else if (ag4_needs_mh_injection(mode)) {
            ag4_mh_injection_params_t mp = ag4_mh_injection_params(mode);
            coord_dist = mp.coord_dist;
            max_coord = mp.max_coord;
        } else if (ag4_needs_sq_injection(mode)) {
            ag4_sq_injection_params_t sp = ag4_sq_injection_params(mode);
            coord_dist = sp.coord_dist;
            max_coord = sp.max_coord;
        }

        std::vector<ag4_nbp_override> ov;
        ag4_apply_metal_mode(mode, ov);

        for (int i = 0; i < (int)rec.size(); i++) {
            bool is_target = false;
            for (const std::string& mn : metal_names)
                if (rec[i].type_name == mn) { is_target = true; break; }
            if (!is_target) continue;

            std::vector<std::array<double,3>> occ = ag4_collect_coord_dirs(rec, i, 3.4, {"SQ","MH","JT","TZ"});
            ag4_metal_site_state site;
            site.mode = mode;
            site.metal_type = rec[i].type_name;
            site.metal_xyz = vec(rec[i].x, rec[i].y, rec[i].z);
            site.coord_dist = coord_dist;
            site.max_coord = max_coord;
            site.receptor_cn = (int)occ.size();
            site.direct_overrides = ov;
            site.jt_enabled = jt_enabled;
            site.jt_axis = vec(0,0,1);
            site.jt_axial_dist = jt_axial_dist;
            if (jt_enabled) {
                std::array<double,3> axis = ag4_tetragonal_axis(occ);
                site.jt_axis = vec(axis[0], axis[1], axis[2]);
            }

            int water_sites = std::min(2, std::max(0, max_coord - (int)occ.size()));
            std::vector<std::array<double,3>> waters = ag4_select_vacant_dirs(occ, water_sites);
            for (int w = 0; w < (int)waters.size(); w++) {
                ag4_bridge_water_site ws;
                ws.xyz = vec(rec[i].x + coord_dist*waters[w][0], rec[i].y + coord_dist*waters[w][1], rec[i].z + coord_dist*waters[w][2]);
                ws.weight = (w == 0) ? 0.90 : 0.70;
                ws.target_dist = 2.80;
                site.bridge_waters.push_back(ws);
            }
            state.sites.push_back(site);
        }
    }
    return state;
}

// ------------------------------------------------------------------
// P3-OPT: Near-atom cache for inner grid loop
// Precomputed per near atom (within softnbc=8Å) to share sqrt/lookup
// across all probe types, eliminating O(n_valid_probes) redundant work.
// ------------------------------------------------------------------
struct NearAtomCache {
    int    ia;               // atom index → rec[ia] for H-bond properties
    int    ridx;             // receptor atom type index
    int    indx;             // ag4_lookup(r); valid as both LJ and sol index within 8Å
    float  dnorm_x, dnorm_y, dnorm_z;   // normalised grid→atom direction
    double rec_vol_solval;              // rec[ia].vol * sol_fn[indx]
    double rec_sp_chg_solval;           // (rec[ia].solpar + SOLPAR_Q*|chg|) * sol_fn[indx]
};

// ------------------------------------------------------------------
// Public entry point
// ------------------------------------------------------------------
ag4_inline_result ag4_compute_maps(const std::string& receptor_pdbqt_path,
                                    const vec&         center,
                                    const vec&         box_size,
                                    double             spacing,
                                    const std::vector<sz>& ligand_ad_types,
                                    const std::vector<ag4_nbp_override>& nbp_overrides,
                                    bool zn_mode,
                                    ag4_metal_mode metal_mode,
                                    const std::vector<ag4_metal_mode>& extra_metal_modes)
{
    ag4_inline_result result;

    // ---- 1. Parse parameters ----
    ag4_param_table params;
    double coeff_vdW, coeff_hbond, coeff_estat, coeff_desolv;
    ag4_parse_params(params, coeff_vdW, coeff_hbond, coeff_estat, coeff_desolv);

    // ---- 2. Parse receptor ----
    std::vector<ag4_rec_atom> rec;
    ag4_parse_receptor(receptor_pdbqt_path, params, rec);
    if (rec.empty()) return result;     // nothing to compute

    // ---- 2b. zn_mode: accumulate AD4Zn pairwise overrides + zero Zn charges ----
    std::vector<ag4_nbp_override> effective_overrides = nbp_overrides;
    bool effective_zn_mode = zn_mode || metal_mode == ag4_metal_mode::zn;
    for (ag4_metal_mode mm : extra_metal_modes)
        if (mm == ag4_metal_mode::zn) effective_zn_mode = true;
    if (effective_zn_mode) {
        std::vector<ag4_nbp_override> zn_ovr = ag4_zn_mode_overrides();
        for (const auto& x : zn_ovr)
            ag4_append_nbp_pair_if_missing(effective_overrides, x);
        for (auto& ra : rec) {
            if (ra.type_name == "Zn")
                ra.charge = 0.0;  // suppress electrostatic dominance of bare Zn ion
        }
        // Auto-inject TZ pseudoatoms at vacant Zn tetrahedral coordination sites.
        // Only injects for Zn atoms that don't already have TZ neighbours in the PDBQT.
        bool has_tz = false;
        for (const auto& ra : rec) if (ra.type_name == "TZ") { has_tz = true; break; }
        if (!has_tz)
            ag4_inject_tz_pseudoatoms(rec, params);
    }

    std::vector<ag4_metal_mode> active_modes = ag4_collect_active_modes(metal_mode, extra_metal_modes);
    for (ag4_metal_mode mode : active_modes) {
        if (ag4_needs_sq_injection(mode)) {
            std::vector<ag4_nbp_override> sq_ovr = ag4_sq_nbp_overrides(mode);
            for (const auto& x : sq_ovr)
                ag4_append_nbp_pair_if_missing(effective_overrides, x);
            ag4_sq_injection_params_t sp = ag4_sq_injection_params(mode);
            if (!sp.metal_names.empty())
                ag4_inject_sq_pseudoatoms(rec, params, sp.metal_names, sp.coord_dist, sp.max_coord);
        }
        if (ag4_needs_mh_injection(mode)) {
            std::vector<ag4_nbp_override> mh_ovr = ag4_mh_nbp_overrides(mode);
            for (const auto& x : mh_ovr)
                ag4_append_nbp_pair_if_missing(effective_overrides, x);
            ag4_mh_injection_params_t mp = ag4_mh_injection_params(mode);
            if (!mp.metal_names.empty())
                ag4_inject_3d_pseudoatoms(rec, params, mp.metal_names, "MH", mp.coord_dist, mp.max_coord);
        }
        if (ag4_needs_jt_injection(mode)) {
            std::vector<ag4_nbp_override> jt_ovr = ag4_jt_nbp_overrides(mode);
            for (const auto& x : jt_ovr)
                ag4_append_nbp_pair_if_missing(effective_overrides, x);
            ag4_jt_injection_params_t jp = ag4_jt_injection_params(mode);
            ag4_inject_jt_pseudoatoms(rec, params, jp);
        }
        if (ag4_needs_tz_injection(mode)) {
            std::vector<ag4_nbp_override> tz_ovr = ag4_tz_nbp_overrides(mode);
            for (const auto& x : tz_ovr)
                ag4_append_nbp_pair_if_missing(effective_overrides, x);
            ag4_tz_injection_params_t tp = ag4_tz_injection_params(mode);
            if (!tp.metal_names.empty())
                ag4_inject_tz_pseudoatoms(rec, params, tp.metal_names, tp.coord_dist, tp.max_coord);
        }
    }
    ag4_append_ligand_metal_overrides(effective_overrides);

    // ---- 3. H-bond geometry ----
    ag4_setup_hbond_geometry(rec);

    // ---- 4. Grid dimensions ----
    int nx = (int)std::ceil(box_size[0] / spacing);
    int ny = (int)std::ceil(box_size[1] / spacing);
    int nz = (int)std::ceil(box_size[2] / spacing);
    if (nx % 2 != 0) nx++;
    if (ny % 2 != 0) ny++;
    if (nz % 2 != 0) nz++;

    grid_dims& gd = result.gd;
    for (int d = 0; d < 3; d++) {
        int np = (d==0?nx : d==1?ny : nz);
        double hs = np * spacing / 2.0;
        gd[d].n_voxels = np;
        gd[d].begin    = center[d] - hs;
        gd[d].end      = center[d] + hs;
    }

    // Number of grid SAMPLE POINTS per axis = n_voxels + 1
    int px = nx+1, py = ny+1, pz = nz+1;
    int npts = px * py * pz;

    // ---- 5. Collect active probe types ----
    std::vector<sz> probe_types;
    for (size_t k=0; k<ligand_ad_types.size(); k++) {
        sz t = ligand_ad_types[k];
        if (t==AD_TYPE_G0||t==AD_TYPE_G1||t==AD_TYPE_G2||t==AD_TYPE_G3) continue;
        if (t==AD_TYPE_CG0||t==AD_TYPE_CG1||t==AD_TYPE_CG2||t==AD_TYPE_CG3) continue;
        if (t==AD_TYPE_W) continue;
        probe_types.push_back(t);
    }
    int Np = (int)probe_types.size();
    if (Np == 0) return result;

    // ---- 6. Collect unique receptor types ----
    std::vector<std::string> rec_type_names;
    for (size_t ia=0; ia<rec.size(); ia++) {
        bool found=false;
        for (size_t k=0; k<rec_type_names.size(); k++)
            if (rec_type_names[k]==rec[ia].type_name) { found=true; break; }
        if (!found) rec_type_names.push_back(rec[ia].type_name);
    }
    int Nr = (int)rec_type_names.size();
    std::map<std::string,int> rec_idx_map;
    for (int k=0; k<Nr; k++) rec_idx_map[rec_type_names[k]] = k;
    // Per-atom rec type index
    std::vector<int> atom_ridx(rec.size());
    for (size_t ia=0; ia<rec.size(); ia++)
        atom_ridx[ia] = rec_idx_map[rec[ia].type_name];

    // ---- 7. Build LJ tables [probe_idx][rec_type_idx][NEINT] ----
    // Also store hb_pair flags and probe parameters
    std::vector<std::vector<std::vector<double>>> lj(Np,
        std::vector<std::vector<double>>(Nr));
    std::vector<bool>   valid_probe  (Np, false);  // false if type not in param table
    std::vector<bool>   probe_is_hbonder(Np, false);
    std::vector<bool>   hb_pair_flag(Np * Nr, false); // [pi*Nr + ri]
    std::vector<double> probe_vol   (Np, 0.0);
    std::vector<double> probe_solpar(Np, 0.0);
    std::vector<int>    probe_hbond (Np, 0);

    for (int pi=0; pi<Np; pi++) {
        sz t = probe_types[pi];
        sz tt = t;
        std::string pname = get_adtype_str(tt);
        ag4_param_table::const_iterator pit = params.find(pname);
        if (pit == params.end()) continue;  // type not in param file — skip entirely
        valid_probe[pi] = true;
        const ag4_atom_par& pp = pit->second;
        probe_vol[pi]    = pp.vol;
        probe_solpar[pi] = pp.solpar;
        probe_hbond[pi]  = pp.hbond;
        probe_is_hbonder[pi] = (pp.hbond > 0);

        for (int ri=0; ri<Nr; ri++) {
            ag4_param_table::const_iterator rit = params.find(rec_type_names[ri]);
            if (rit == params.end()) continue;
            const ag4_atom_par& rp = rit->second;

            int xA=12, xB=6;
            double nbp_r, nbp_eps;
            bool is_hb = false;

            // Check for explicit nbp_r_eps overrides first (e.g. from --zn_mode)
            bool overridden = false;
            for (const auto& ovr : effective_overrides) {
                if (ovr.probe == pname && ovr.receptor == rec_type_names[ri]) {
                    nbp_r   = ovr.r_eq;
                    nbp_eps = ovr.eps;
                    xA      = ovr.xA;
                    xB      = ovr.xB;
                    overridden = true;
                    break;
                }
            }

            if (!overridden) {
                // Standard AD4 pairwise: H-bond or vdW
                // Probe is acceptor (AS/A1/A2, hbond>=3) and rec is donor (DS/D1, hbond=1 or 2)
                if (pp.hbond >= 3 && (rp.hbond==1||rp.hbond==2)) {
                    xB=10;
                    nbp_r   = pp.Rij_hb;
                    nbp_eps = pp.epsij_hb * coeff_hbond;
                    is_hb   = true;
                }
                // Probe is donor (DS/D1, hbond=1 or 2) and rec is acceptor (AS/A1/A2, hbond>=3)
                else if ((pp.hbond==1||pp.hbond==2) && rp.hbond >= 3) {
                    xB=10;
                    nbp_r   = rp.Rij_hb;
                    nbp_eps = rp.epsij_hb * coeff_hbond;
                    is_hb   = true;
                } else {
                    nbp_r   = (pp.Rij + rp.Rij) / 2.0;
                    nbp_eps = std::sqrt(pp.epsij * coeff_vdW * rp.epsij * coeff_vdW);
                }
            }

            ag4_build_lj_table(nbp_r, nbp_eps, xA, xB, lj[pi][ri]);
            hb_pair_flag[pi*Nr + ri] = is_hb;
        }
    }

    // ---- 8. Desolvation table ----
    std::vector<double> sol_fn;
    ag4_build_sol_fn(coeff_desolv, sol_fn);

    // ---- 9. Precompute H-bond donor list (D1 atoms only, for closestH) ----
    std::vector<int> hbond_donors;
    for (int ia=0; ia<(int)rec.size(); ia++)
        if (rec[ia].hbond==1||rec[ia].hbond==2) hbond_donors.push_back(ia);

    // ---- 10. Allocate output arrays ----
    result.ad_types = probe_types;
    result.aff_maps.assign(Np, std::vector<double>(npts, 0.0));
    result.elec_map  .assign(npts, 0.0);
    result.desolv_map.assign(npts, 0.0);

    // Electrostatic constant: factor / (diel * coeff_estat) ... actually:
    // E_elec = charge * (1/max(r,0.5)) * (factor/diel) * coeff_estat
    double invdielcal = (AG4_FACTOR / AG4_DIEL) * coeff_estat;

    // P2: SoA pre-packing of receptor coords/charge/vol for SIMD-friendly inner loops
    // Contiguous scalar arrays enable auto-vectorisation + omp simd on the tight loops.
    const int Nrec = (int)rec.size();

    // Print grid info so users can diagnose slow runs (output goes to docking log via stdout redirect)
    {
        long long total_ops = (long long)npts * Nrec;
        std::cout << "  [AD4 grid] " << nx << " x " << ny << " x " << nz
                  << " = " << npts << " pts, receptor=" << Nrec
                  << " atoms, probes=" << Np << ", spacing=" << spacing << " A" << std::endl;
        if (total_ops > 5000000000LL) {
            std::cout << "  [WARNING] Large grid (" << npts << " pts, spacing=" << spacing
                      << " A). Map generation may take several minutes."
                      << " Use --spacing 0.375 for faster results (~"
                      << (int)std::round((double)npts / ((int)(std::ceil(box_size[0]/0.375)+2)
                                                       * (int)(std::ceil(box_size[1]/0.375)+2)
                                                       * (int)(std::ceil(box_size[2]/0.375)+2)))
                      << "x speedup)." << std::endl;
        }
        std::cout.flush();
    }
    std::vector<double> rec_px(Nrec), rec_py(Nrec), rec_pz(Nrec);
    std::vector<double> rec_chg(Nrec), rec_vol2(Nrec);
    for (int ia = 0; ia < Nrec; ia++) {
        rec_px[ia]   = rec[ia].x;
        rec_py[ia]   = rec[ia].y;
        rec_pz[ia]   = rec[ia].z;
        rec_chg[ia]  = rec[ia].charge;
        rec_vol2[ia] = rec[ia].vol;
    }
    const double* __restrict__ soa_x   = rec_px.data();
    const double* __restrict__ soa_y   = rec_py.data();
    const double* __restrict__ soa_z   = rec_pz.data();
    const double* __restrict__ soa_chg = rec_chg.data();
    const double* __restrict__ soa_vol = rec_vol2.data();

    // ---- 9b. P3-OPT: precompute valid probe index list and O(1) H-bond donor flag ----
    // valid_pi: only iterate over probe types actually present in the ligand.
    // is_hbd:   O(1) check for H-bond donor atoms (avoids scanning hbond_donors per grid pt).
    std::vector<int> valid_pi;
    valid_pi.reserve(static_cast<size_t>(Np));
    for (int pi = 0; pi < Np; pi++)
        if (valid_probe[pi]) valid_pi.push_back(pi);

    std::vector<uint8_t> is_hbd(static_cast<size_t>(Nrec), 0);
    for (int ia : hbond_donors) is_hbd[static_cast<size_t>(ia)] = 1;

    // ---- 11. Grid point loop ----
    double softnbc2 = AG4_SOFTNBC * AG4_SOFTNBC;
    double dielcut2 = AG4_DIELCUT * AG4_DIELCUT;

#ifdef _OPENMP
#  pragma omp parallel for schedule(dynamic) default(shared)
#endif
    for (int iz = 0; iz < pz; iz++) {
        // P3-OPT: thread-private near-atom cache; allocated once per z-slice,
        // reused (via clear()) for every grid point within the slice.
        std::vector<NearAtomCache> near_buf;
        near_buf.reserve(static_cast<size_t>(Nrec));

        double gz = gd[2].begin + iz * spacing;
        for (int iy = 0; iy < py; iy++) {
            double gy = gd[1].begin + iy * spacing;
            for (int ix = 0; ix < px; ix++) {
                double gx = gd[0].begin + ix * spacing;
                int flat = iz * py * px + iy * px + ix;

                // ---- P3-OPT: Single merged pass over all receptor atoms ----
                // Computes elec + desolv + H-bond donor tracking +
                // near-atom cache in one traversal, sharing sqrt/lookup
                // across all subsequent probe types.
                near_buf.clear();
                double e_elec = 0.0, desolv_acc = 0.0;
                double rminH = 1e30;
                int    closestH = -1;
                for (int ia = 0; ia < Nrec; ia++) {
                    double dx = soa_x[ia]-gx, dy = soa_y[ia]-gy, dz = soa_z[ia]-gz;
                    double d2 = dx*dx+dy*dy+dz*dz;
                    if (d2 >= dielcut2) continue;
                    double r = std::sqrt(d2);
                    if (r < AG4_APPROX0) r = AG4_APPROX0;
                    // Electrostatics (all atoms within dielcut)
                    e_elec += soa_chg[ia] / std::max(r, 0.5);
                    // Desolv + near-atom cache (atoms within softnbc only)
                    if (d2 < softnbc2) {
                        int    indx    = std::min(ag4_lookup(r), AG4_NEINT-1);
                        double sol_val = sol_fn[indx];
                        desolv_acc    += soa_vol[ia] * sol_val;
                        double inv_r   = 1.0 / r;
                        near_buf.push_back({
                            ia,
                            atom_ridx[ia],
                            indx,
                            (float)(dx * inv_r),
                            (float)(dy * inv_r),
                            (float)(dz * inv_r),
                            rec[ia].vol * sol_val,
                            (rec[ia].solpar + AG4_SOLPAR_Q * std::fabs(rec[ia].charge)) * sol_val
                        });
                        // Closest H-bond donor for Hramp (O(1) flag lookup)
                        if (is_hbd[ia] && r < rminH) { rminH = r; closestH = ia; }
                    }
                }
                result.elec_map[flat]   = e_elec * invdielcal;
                result.desolv_map[flat] = AG4_SOLPAR_Q * desolv_acc;

                // ---- Affinity maps (probe loop over valid_pi only) ----
                for (int pi : valid_pi) {
                    double E = 0.0;
                    double hbondmin =  1e30;
                    double hbondmax = -1e30;
                    bool   hbondflag = false;

                    for (const NearAtomCache& na : near_buf) {
                        const double dnorm_x = na.dnorm_x;
                        const double dnorm_y = na.dnorm_y;
                        const double dnorm_z = na.dnorm_z;
                        int    ri    = na.ridx;
                        int    ia    = na.ia;
                        double lj_e  = lj[pi][ri][na.indx];
                        bool   is_hb = hb_pair_flag[pi*Nr + ri];

                        // ---------- H-bond directional factors ----------
                        double racc = 1.0, rdon = 1.0, Hramp = 1.0;
                        if (is_hb && probe_is_hbonder[pi]) {
                            int phb = probe_hbond[pi];
                            int rhb = rec[ia].hbond;

                            // Probe acceptor, rec donor (D1)
                            if (phb >= 3 && rhb == 2 && rec[ia].has_hb_vec) {
                                // racc = [cos theta]^rexp where theta = angle between
                                // H->grid_point vector and N->H (rvector)
                                double cos_t = -(dnorm_x*rec[ia].rvector[0]
                                              +  dnorm_y*rec[ia].rvector[1]
                                              +  dnorm_z*rec[ia].rvector[2]);
                                if (cos_t <= 0.0) {
                                    racc = 0.0;
                                } else {
                                    switch (rec[ia].rexp) {
                                        case 4: { double t=cos_t*cos_t; racc=t*t; break; }
                                        case 2: racc=cos_t*cos_t; break;
                                        default: racc=cos_t; break;
                                    }
                                }
                                // Hramp: ramp for multiple donors
                                if (closestH == ia) {
                                    Hramp = 1.0;
                                } else if (closestH >= 0 && rec[closestH].has_hb_vec) {
                                    double ct2 = 0.0;
                                    for (int d=0;d<3;d++)
                                        ct2 += rec[closestH].rvector[d]*rec[ia].rvector[d];
                                    ct2 = std::max(-1.0, std::min(1.0, ct2));
                                    double theta = std::acos(ct2);
                                    Hramp = 0.5 - 0.5*std::cos(theta * 120.0/90.0);
                                }
                                double rsph = lj_e / 100.0;
                                rsph = std::max(0.0, std::min(1.0, rsph));
                                E += lj_e * Hramp * (racc + (1.0-racc)*rsph);

                            // Probe donor, rec acceptor (A1): directional N
                            } else if ((phb==1||phb==2) && rhb == 4 && rec[ia].has_hb_vec) {
                                double cos_t = -(dnorm_x*rec[ia].rvector[0]
                                              +  dnorm_y*rec[ia].rvector[1]
                                              +  dnorm_z*rec[ia].rvector[2]);
                                rdon = (cos_t > 0.0) ? cos_t*cos_t : 0.0;
                                double rsph = lj_e / 100.0;
                                rsph = std::max(0.0, std::min(1.0, rsph));
                                double hb_e = lj_e * (rdon + (1.0-rdon)*rsph);
                                hbondmin = std::min(hbondmin, hb_e);
                                hbondmax = std::max(hbondmax, hb_e);
                                hbondflag = true;

                            // Probe donor, rec acceptor (A2): lone-pair oxygen
                            } else if ((phb==1||phb==2) && rhb == 5 && rec[ia].has_hb_vec) {
                                double cos_t = -(dnorm_x*rec[ia].rvector[0]
                                              +  dnorm_y*rec[ia].rvector[1]
                                              +  dnorm_z*rec[ia].rvector[2]);
                                // t0: angle out of lone-pair plane
                                double t0 = dnorm_x*rec[ia].rvector2[0]
                                          + dnorm_y*rec[ia].rvector2[1]
                                          + dnorm_z*rec[ia].rvector2[2];
                                t0 = std::max(-1.0, std::min(1.0, t0));
                                t0 = M_PI/2.0 - std::acos(t0);
                                // ti: angle in lone-pair plane
                                double cross_x = dnorm_y*rec[ia].rvector2[2] - dnorm_z*rec[ia].rvector2[1];
                                double cross_y = dnorm_z*rec[ia].rvector2[0] - dnorm_x*rec[ia].rvector2[2];
                                double cross_z = dnorm_x*rec[ia].rvector2[1] - dnorm_y*rec[ia].rvector2[0];
                                double rd2 = cross_x*cross_x+cross_y*cross_y+cross_z*cross_z;
                                if (rd2 < AG4_APPROX0) rd2 = AG4_APPROX0;
                                double inv_rd = 1.0 / std::sqrt(rd2);
                                double ti = (cross_x*inv_rd*rec[ia].rvector[0]
                                           + cross_y*inv_rd*rec[ia].rvector[1]
                                           + cross_z*inv_rd*rec[ia].rvector[2]);
                                if (cos_t >= 0.0) {
                                    ti = std::max(-1.0, std::min(1.0, ti));
                                    ti = std::acos(ti) - M_PI/2.0;
                                    if (ti < 0.0) ti = -ti;
                                    rdon = (0.9 + 0.1*std::sin(ti+ti)) * std::cos(t0);
                                } else if (cos_t >= -0.34202) {
                                    rdon = 562.25*std::pow(0.116978-cos_t*cos_t,3.0)*std::cos(t0);
                                } else {
                                    rdon = 0.0;
                                }
                                rdon = std::max(0.0, rdon);
                                double rsph = lj_e / 100.0;
                                rsph = std::max(0.0, std::min(1.0, rsph));
                                double hb_e = lj_e * (rdon + (1.0-rdon)*rsph);
                                hbondmin = std::min(hbondmin, hb_e);
                                hbondmax = std::max(hbondmax, hb_e);
                                hbondflag = true;

                            // Probe acceptor, rec non-directional donor (DS)
                            } else if (phb >= 3 && rhb == 1) {
                                double rsph = lj_e / 100.0;
                                rsph = std::max(0.0, std::min(1.0, rsph));
                                E += lj_e * (1.0 + (1.0-1.0)*rsph);   // racc=1 for DS

                            } else {
                                // hbonder probe but non-H-bond interaction with this atom
                                E += lj_e;
                            }
                        } else {
                            // No H-bond: straight vdW
                            E += lj_e;
                        }

                        // ---- Desolvation contribution to affinity map ----
                        // P3-OPT: na.rec_vol_solval and na.rec_sp_chg_solval already
                        // incorporate sol_fn[indx], precomputed in the merged pass.
                        E += probe_solpar[pi] * na.rec_vol_solval
                           + probe_vol[pi]    * na.rec_sp_chg_solval;
                    } // near-atom loop

                    // Apply directional H-bond correction (min+max rule)
                    if (hbondflag) {
                        E += hbondmin + std::max(hbondmax, 0.0);
                    }
                    // AutoGrid4-compatible final clamp: affinity-map values are
                    // bounded to ±AG4_EINTCLAMP (±1000 kcal/mol) at each grid
                    // point, exactly as AutoGrid 4.2's MAXVALUE behaviour.
                    result.aff_maps[pi][flat] = std::max(-AG4_EINTCLAMP,
                                                         std::min(AG4_EINTCLAMP, E));
                } // probe type loop
            } // ix
        } // iy
    } // iz

    result.valid = true;
    return result;
}
