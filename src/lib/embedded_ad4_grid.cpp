/*
   LKina — embedded AD4 grid generation
   Copyright (C) 2025 LK-Studio1128 and LKina contributors
   SPDX-License-Identifier: GPL-3.0-or-later

   Generates AutoGrid4-compatible affinity maps inline via ag4_engine,
   which is derived from AutoGrid4 (GPL-2.0 or later). This file is
   therefore also distributed under the GNU General Public License v3.

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

#include "embedded_ad4_grid.h"
#include "ag4_engine.h"
#include "ad4_parameter_data.h"
#include "atom_constants.h"
#include "file.h"

std::string get_adtype_str(sz& t);

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <fstream>
#include <stdexcept>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

namespace fs = boost::filesystem;

namespace {

std::string make_tmp_dir() {
    char tmpl[] = "/tmp/vina_ad4_XXXXXX";
    char* d = mkdtemp(tmpl);
    if (!d) throw std::runtime_error("embedded_ad4_grid: mkdtemp failed");
    return std::string(d);
}

void write_builtin_parameter_file(const std::string& path) {
    std::ofstream ofs(path.c_str());
    if (!ofs) throw std::runtime_error("embedded_ad4_grid: cannot write parameter file to " + path);
    ofs << ad4_param::builtin_ad4_parameter_text;
}

void write_gpf(const std::string& gpf_path,
               const std::string& receptor_path,
               const std::string& param_path,
               const std::string& map_prefix,
               const vec& center,
               const vec& box_size,
               fl spacing,
               const std::vector<sz>& ligand_ad_types) {

    std::ofstream gpf(gpf_path.c_str());
    if (!gpf) throw std::runtime_error("embedded_ad4_grid: cannot write GPF to " + gpf_path);

    // npts: must be even
    int nx = (int)std::ceil(box_size[0] / spacing);
    int ny = (int)std::ceil(box_size[1] / spacing);
    int nz = (int)std::ceil(box_size[2] / spacing);
    if (nx % 2 != 0) nx++;
    if (ny % 2 != 0) ny++;
    if (nz % 2 != 0) nz++;

    gpf << "parameter_file " << param_path << "\n";
    gpf << "npts " << nx << " " << ny << " " << nz << "\n";
    gpf << "gridfld " << map_prefix << ".maps.fld\n";
    gpf << "spacing " << spacing << "\n";
    gpf << "receptor_types ";
    gpf << "C A N OA HD SA NA Mg Ca Fe Mn Zn\n";
    gpf << "ligand_types ";
    bool first = true;
    for (sz t : ligand_ad_types) {
        if (t >= AD_TYPE_G0 && t <= AD_TYPE_G3) continue;
        if (t >= AD_TYPE_CG0 && t <= AD_TYPE_CG3) continue;
        if (t == AD_TYPE_W) continue;
        if (!first) gpf << " ";
        first = false;
        sz tt = t;
        gpf << get_adtype_str(tt);
    }
    gpf << "\n";
    gpf << "receptor " << receptor_path << "\n";
    gpf << "gridcenter " << center[0] << " " << center[1] << " " << center[2] << "\n";
    gpf << "smooth 0.5\n";
    gpf << "map " << map_prefix << ".\n"; // autogrid4 will append type+.map
    gpf << "elecmap " << map_prefix << ".e.map\n";
    gpf << "dsolvmap " << map_prefix << ".d.map\n";
    gpf << "dielectric -0.1465\n";
}

std::vector<fl> read_flat_map(const std::string& path) {
    std::ifstream ifs(path.c_str());
    if (!ifs) throw std::runtime_error("embedded_ad4_grid: cannot read map: " + path);

    std::vector<fl> data;
    std::string line;
    int line_no = 0;
    while (std::getline(ifs, line)) {
        ++line_no;
        if (line_no <= 6) continue;
        if (line.empty()) continue;
        data.push_back((fl)std::atof(line.c_str()));
    }
    return data;
}

vec read_map_origin_and_spacing(const std::string& path, fl& out_spacing,
                                 sz& out_nx, sz& out_ny, sz& out_nz) {
    std::ifstream ifs(path.c_str());
    if (!ifs) throw std::runtime_error("embedded_ad4_grid: cannot open map: " + path);

    std::string line;
    int line_no = 0;
    fl spacing = 0.375f;
    sz nx = 0, ny = 0, nz = 0;
    fl cx = 0, cy = 0, cz = 0;
    while (std::getline(ifs, line) && line_no < 7) {
        ++line_no;
        if (line_no == 4) { // SPACING
            std::istringstream ss(line);
            std::string tok; ss >> tok >> spacing;
        }
        if (line_no == 5) { // NELEMENTS
            std::istringstream ss(line);
            std::string tok; ss >> tok >> nx >> ny >> nz;
        }
        if (line_no == 6) { // CENTER
            std::istringstream ss(line);
            std::string tok; ss >> tok >> cx >> cy >> cz;
        }
    }
    out_spacing = spacing;
    out_nx = nx; out_ny = ny; out_nz = nz;
    fl hx = nx * spacing / 2.0f;
    fl hy = ny * spacing / 2.0f;
    fl hz = nz * spacing / 2.0f;
    vec origin;
    origin[0] = cx - hx;
    origin[1] = cy - hy;
    origin[2] = cz - hz;
    return origin;
}

} // anonymous namespace

ad4_grid_data generate_embedded_ad4_maps(const embedded_ad4_request& req) {
    ad4_grid_data result;
    result.valid = false;

    // --- Try fully inline computation first (no external executable) ---
    {
        // Merge metal-mode overrides via the central dispatch helper (supports multi-metal)
        std::vector<ag4_nbp_override> combined_overrides = req.nbp_overrides;
        ag4_apply_metal_mode(req.metal_mode, combined_overrides);
        for (ag4_metal_mode mm : req.extra_metal_modes)
            ag4_apply_metal_mode(mm, combined_overrides);

        ag4_inline_result ir = ag4_compute_maps(
            req.receptor_pdbqt_path,
            req.center, req.box_size, req.spacing,
            req.ligand_ad_types,
            combined_overrides,
            req.zn_mode,
            req.metal_mode,
            req.extra_metal_modes);
        if (ir.valid) {
            result.valid            = true;
            result.map_prefix       = "";  // empty: populated in-memory, no files
            result.spacing          = req.spacing;
            result.nx               = (sz)ir.gd[0].n_voxels;
            result.ny               = (sz)ir.gd[1].n_voxels;
            result.nz               = (sz)ir.gd[2].n_voxels;
            result.origin[0]        = ir.gd[0].begin;
            result.origin[1]        = ir.gd[1].begin;
            result.origin[2]        = ir.gd[2].begin;
            result.gd               = ir.gd;
            result.map_ad_types     = ir.ad_types;
            result.affinity_maps.clear();
            for (size_t k = 0; k < ir.aff_maps.size(); k++) {
                const std::vector<double>& src = ir.aff_maps[k];
                result.affinity_maps.push_back(std::vector<fl>(src.begin(), src.end()));
            }
            result.electrostatic_map  = std::vector<fl>(ir.elec_map.begin(),   ir.elec_map.end());
            result.desolvation_map    = std::vector<fl>(ir.desolv_map.begin(), ir.desolv_map.end());
            return result;
        }
    }

    if (req.zn_mode || req.metal_mode != ag4_metal_mode::none ||
        !req.extra_metal_modes.empty() || !req.nbp_overrides.empty()) {
        throw std::runtime_error("embedded_ad4_grid: inline metal/nbp map generation failed; external autogrid4 fallback is disabled because it cannot reproduce LKina metal pseudoatoms and nbp_r_eps overrides");
    }

    // --- Fallback: external autogrid4 subprocess ---
    if (req.autogrid4_exe.empty()) {
        return result;  // no fallback available
    }

    // --- Build working directory ---
    std::string work_dir = req.work_dir.empty() ? make_tmp_dir() : req.work_dir;
    fs::create_directories(work_dir);

    // --- Write builtin parameter file if needed ---
    std::string param_path = req.parameter_file;
    if (param_path.empty()) {
        param_path = work_dir + "/ad4_params.dat";
        write_builtin_parameter_file(param_path);
    }

    // --- Generate map file prefix ---
    std::string map_prefix = work_dir + "/receptor";

    // --- Write GPF ---
    std::string gpf_path = work_dir + "/grid.gpf";
    write_gpf(gpf_path, req.receptor_pdbqt_path, param_path,
              map_prefix, req.center, req.box_size, req.spacing, req.ligand_ad_types);

    // --- Invoke autogrid4 ---
    std::string glg_path = work_dir + "/grid.glg";
    std::string cmd = "\"" + req.autogrid4_exe + "\" -p \"" + gpf_path + "\" -l \"" + glg_path + "\" 2>&1";
    int ret = std::system(cmd.c_str());
    if (ret != 0) {
        throw std::runtime_error("embedded_ad4_grid: autogrid4 returned non-zero exit code");
    }

    // --- Read maps into memory ---
    fl spacing;
    sz nx, ny, nz;

    // find first valid map for geometry
    std::string first_map;
    for (sz t : req.ligand_ad_types) {
        if (t >= AD_TYPE_G0 && t <= AD_TYPE_G3) continue;
        if (t >= AD_TYPE_CG0 && t <= AD_TYPE_CG3) continue;
        if (t == AD_TYPE_W) continue;
        sz tt = t;
        std::string name = get_adtype_str(tt);
        std::string p = map_prefix + "." + name + ".map";
        if (fs::exists(p)) { first_map = p; break; }
    }
    if (first_map.empty()) {
        throw std::runtime_error("embedded_ad4_grid: no affinity maps found after autogrid4 run");
    }

    result.map_prefix = map_prefix;
    result.origin  = read_map_origin_and_spacing(first_map, spacing, nx, ny, nz);
    result.spacing = spacing;
    result.nx = nx; result.ny = ny; result.nz = nz;

    for (sz t : req.ligand_ad_types) {
        if (t >= AD_TYPE_G0 && t <= AD_TYPE_G3) continue;
        if (t >= AD_TYPE_CG0 && t <= AD_TYPE_CG3) continue;
        if (t == AD_TYPE_W) continue;
        sz tt = t;
        std::string name = get_adtype_str(tt);
        std::string p = map_prefix + "." + name + ".map";
        if (fs::exists(p)) {
            result.map_ad_types.push_back(t);
            result.affinity_maps.push_back(read_flat_map(p));
        }
    }

    std::string emap = map_prefix + ".e.map";
    std::string dmap = map_prefix + ".d.map";
    if (fs::exists(emap)) result.electrostatic_map  = read_flat_map(emap);
    if (fs::exists(dmap)) result.desolvation_map    = read_flat_map(dmap);

    result.valid = true;
    return result;
}
