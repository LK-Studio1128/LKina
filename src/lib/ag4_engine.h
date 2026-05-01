/*
   LKina — AutoGrid4 scoring function engine (header)
   Copyright (C) 2025 LK-Studio1128 and LKina contributors
   SPDX-License-Identifier: GPL-3.0-or-later
   See COPYING for the full GNU General Public License v3 text.
*/

#ifndef VINA_AG4_ENGINE_H
#define VINA_AG4_ENGINE_H

#include <string>
#include <vector>
#include "common.h"
#include "grid_dim.h"

// Pairwise nonbonded potential override (from AutoDock4 nbp_r_eps / nbp_coeffs)
// Used by --zn_mode (AutoDock4Zn) and future metal-specific modes.
struct ag4_nbp_override {
    std::string probe;    // ligand probe atom type name  (e.g. "NA")
    std::string receptor; // receptor atom type name      (e.g. "TZ")
    double r_eq;          // equilibrium distance (Å)     (e.g. 0.25)
    double eps;           // well depth (kcal/mol)         (e.g. 23.2135)
    int    xA = 12;       // repulsive exponent
    int    xB = 6;        // attractive exponent
};

// AD4Zn hardcoded pairwise overrides (Santos-Martins et al. 2014, J Chem Inf Model 54:2371)
inline std::vector<ag4_nbp_override> ag4_zn_mode_overrides() {
    return {
        {"NA",  "TZ",  0.25,   23.2135, 12, 6},  // TZ–NA tetrahedral attraction
        {"OA",  "Zn",  2.1,     3.8453, 12, 6},  // Zn–OA coordination
        {"SA",  "Zn",  2.25,    7.5914, 12, 6},  // Zn–SA coordination
        {"HD",  "Zn",  1.0,     0.0,    12, 6},  // Zn–HD repulsion only
        {"NA",  "Zn",  2.0,     0.0060, 12, 6},  // Zn–NA weak direct
        {"N",   "Zn",  2.0,     0.2966, 12, 6},  // Zn–N  weak direct
    };
}

inline bool ag4_is_coord_donor_probe_type(const std::string& t) {
    return t == "OA" || t == "O" || t == "NA" || t == "N" || t == "SA" || t == "S";
}

inline bool ag4_is_coord_pseudo_type(const std::string& t) {
    return t == "TZ" || t == "SQ" || t == "MH" || t == "JT";
}

inline bool ag4_has_nbp_pair(const std::vector<ag4_nbp_override>& ov, const ag4_nbp_override& x) {
    for (const auto& o : ov)
        if (o.probe == x.probe && o.receptor == x.receptor)
            return true;
    return false;
}

inline void ag4_append_nbp_pair_if_missing(std::vector<ag4_nbp_override>& ov, const ag4_nbp_override& x) {
    if (!ag4_has_nbp_pair(ov, x))
        ov.push_back(x);
}

inline bool ag4_metaldock_reverse_override(const ag4_nbp_override& x, ag4_nbp_override& r) {
    double e_na = 0.0, e_oa = 0.0, e_sa = 0.0, e_hd = 0.0;
    bool supported = true;
    if      (x.receptor == "V")   { e_na = 4.696; e_oa = 6.825; e_sa = 5.658; e_hd = 3.984; }
    else if (x.receptor == "Cr1" || x.receptor == "Cr") { e_na = 6.371; e_oa = 1.998; e_sa = 0.144; e_hd = 3.625; }
    else if (x.receptor == "Co")  { e_na = 5.280; e_oa = 0.050; e_sa = 6.673; e_hd = 5.929; }
    else if (x.receptor == "Ni")  { e_na = 0.630; e_oa = 2.732; e_sa = 4.462; e_hd = 2.820; }
    else if (x.receptor == "Cu")  { e_na = 4.696; e_oa = 1.277; e_sa = 6.791; e_hd = 1.114; }
    else if (x.receptor == "Mo")  { e_na = 1.330; e_oa = 0.014; e_sa = 0.168; e_hd = 5.620; }
    else if (x.receptor == "Ru")  { e_na = 6.936; e_oa = 2.796; e_sa = 4.295; e_hd = 6.357; }
    else if (x.receptor == "Rh")  { e_na = 5.559; e_oa = 2.056; e_sa = 0.573; e_hd = 5.471; }
    else if (x.receptor == "Pd")  { e_na = 4.688; e_oa = 0.845; e_sa = 5.574; e_hd = 3.159; }
    else if (x.receptor == "Re")  { e_na = 6.738; e_oa = 0.645; e_sa = 3.309; e_hd = 4.502; }
    else if (x.receptor == "Os")  { e_na = 5.958; e_oa = 0.135; e_sa = 4.102; e_hd = 6.589; }
    else if (x.receptor == "Pt")  { e_na = 6.532; e_oa = 2.020; e_sa = 6.332; e_hd = 1.844; }
    else supported = false;
    if (!supported) return false;

    if (x.probe == "NA" || x.probe == "N") {
        r = {x.receptor, x.probe, 2.20, e_na, 12, 10};
    } else if (x.probe == "OA" || x.probe == "O") {
        r = {x.receptor, x.probe, 2.25, e_oa, 12, 10};
    } else if (x.probe == "SA" || x.probe == "S") {
        r = {x.receptor, x.probe, 2.30, e_sa, 12, 10};
    } else if (x.probe == "HD") {
        r = {x.receptor, x.probe, 1.00, e_hd, 12, 6};
    } else {
        return false;
    }
    return true;
}

inline void ag4_append_ligand_metal_overrides(std::vector<ag4_nbp_override>& ov) {
    std::vector<ag4_nbp_override> rev;
    for (const auto& o : ov) {
        if (!ag4_is_coord_donor_probe_type(o.probe) && o.probe != "HD") continue;
        if (ag4_is_coord_pseudo_type(o.receptor)) continue;
        ag4_nbp_override r;
        if (!ag4_metaldock_reverse_override(o, r)) {
            if (o.eps <= 0.0) continue;
            r = {o.receptor, o.probe, o.r_eq, o.eps, o.xA, o.xB};
        }
        if (!ag4_has_nbp_pair(ov, r) && !ag4_has_nbp_pair(rev, r))
            rev.push_back(r);
    }
    ov.insert(ov.end(), rev.begin(), rev.end());
}

// Metal coordination mode selector
// Covers first-row transition metals, medicinal organometallics, and alkali/toxicological metals.
// Parameters based on:
//   Bazayeva et al. 2024 (PDB <1.5 Å survey, PMC11066882) for r_eq
//   Hakkennes et al. 2024 MetalDock (MC-optimised LJ, PMC10751784) for ε
//   Harding, Acta Cryst D 2006 (CSD survey) for first-row TM distances
//   AutoDock4Zn (Santos-Martins et al. 2014) as reference framework
enum class ag4_metal_mode {
    none,
    // M1: Zn (special, uses TZ pseudoatom — handled by zn_mode flag, not metal_mode)
    zn,
    // Biological first-row transition metals
    mg,   // Mg2+  octahedral, hard Lewis acid, O/N preference
    ca,   // Ca2+  7-8 coord, very hard, O-only
    mn,   // Mn2+  octahedral, hard, O>N>>S
    fe,   // Fe2+/3+  octahedral, heme + non-heme
    co,   // Co2+  tetrahedral/octahedral, borderline
    ni,   // Ni2+  sq-planar/octahedral, borderline-soft
    cu,   // Cu2+  sq-planar + Jahn-Teller, borderline-soft
    // Medicinal organometallics
    pt,   // Pt2+  sq-planar d8, cisplatin family
    pd,   // Pd2+  sq-planar d8, softer than Pt
    ru,   // Ru3+  octahedral d5, NAMI-A/RAPTA anticancer
    ir,   // Ir3+  octahedral d6, organometallic anticancer
    au,   // Au+   linear d10, auranofin / gold drugs — extreme S-philic
    // Toxicological / environmental
    cd,   // Cd2+  tetrahedral d10, Zn-mimic but softer
    hg,   // Hg2+  linear d10, extreme S-philic, metallothionein inhibition
    // Alkali / alkali-earth (crystal contacts, enzyme cofactors)
    na_ion, // Na+  octahedral, hard, O-only
    k_ion,  // K+   7-8 coord, very hard, large radius
    // ── Early transition metals (hard Lewis acids, O>N) ──────────────────────
    v,      // V3+/4+ octahedral, insulin-mimetic (VOSO4)
    cr,     // Cr3+  octahedral d3, glucose tolerance factor (Cr1 AD_TYPE)
    ti,     // Ti4+  octahedral d0, titanocene anticancer
    sc,     // Sc3+  octahedral d0, 44Sc PET imaging
    y_ion,  // Y3+   octahedral d0, 90Y radiotherapy
    zr,     // Zr4+  8-coord, 89Zr PET, hafnium-like
    nb,     // Nb5+  octahedral d0
    hf,     // Hf4+  8-coord, hafnium nanoparticle radiotherapy
    ta,     // Ta5+  octahedral d0
    tg_mode,// W6+   (Tg AD_TYPE) tungstoenzymes
    mo,     // Mo4+/6+ xanthine oxidase / sulfite oxidase
    // ── Second/third row TMs — borderline to soft ────────────────────────────
    rh,     // Rh3+  octahedral d6, anticancer organometallics
    ag,     // Ag+   linear/tetrahedral d10, antimicrobial
    tc,     // Tc5+  99mTc SPECT (most used radionuclide)
    re,     // Re5+  186Re/188Re radiotherapy
    os,     // Os2+/3+ octahedral d6, osmium arene anticancer
    // ── Post-transition metals / metalloids ──────────────────────────────────
    ga,     // Ga3+  tetrahedral d10, 68Ga PET / antimicrobial
    in_ion, // In3+  octahedral d10, 111In SPECT
    sn,     // Sn4+  octahedral d10, organotin anticancer
    sb,     // Sb3+  octahedral, antiparasitic (Pentostam)
    bi,     // Bi3+  irregular d10, anti-H.pylori (bismuth subcitrate)
    tl,     // Tl+   large d10, 201Tl cardiac SPECT / toxic
    pb,     // Pb2+  hemispherical d10, toxicology target
    // ── s-block metals ───────────────────────────────────────────────────────
    li,     // Li+   tetrahedral, lithium carbonate (GSK-3β)
    al,     // Al3+  tetrahedral, AlF4- phosphatase TS analog
    sr,     // Sr2+  8-coord, strontium ranelate bone drug
    ba,     // Ba2+  9-coord, very hard, large
    // ── Lanthanides (all hard O-preferring, lanthanide contraction) ──────────
    la,     // La3+  9-coord, lanthanum carbonate (Fosrenol)
    ce,     // Ce3+/4+ 8-9 coord
    pr,     // Pr3+  8-coord
    nd,     // Nd3+  8-coord
    sm,     // Sm3+  8-coord, 153Sm radiotherapy
    eu,     // Eu3+  8-coord
    gd,     // Gd3+  8-9 coord, MRI contrast agent
    tb,     // Tb3+  8-coord, 149Tb targeted alpha therapy
    dy,     // Dy3+  8-coord
    ho,     // Ho3+  8-coord, 166Ho radiotherapy
    er,     // Er3+  8-coord
    tm,     // Tm3+  8-coord
    yb,     // Yb3+  8-coord
    lu,     // Lu3+  8-coord, 177Lu Lutathera radiotherapy
    // ── Metalloids with pharmacological interest ─────────────────────────────
    se,     // Se0/2- selenoenzymes, selenocysteine
    as_met, // As3+  arsenic trioxide (Trisenox) AML treatment
    ge,     // Ge4+  organic germanium (spirogermanium) anticancer
    // ── Other main-group elements ────────────────────────────────────────────
    b,      // B3+   boron; bortezomib (Velcade) proteasome inhibitor
    be,     // Be2+  beryllium; occupational toxicology / metalloprotein research
    te,     // Te2-  tellurium; organic telluride anticancer / antibacterials
    po,     // Po2+/4+ polonium; research / environmental toxicology
    at,     // At-   astatine; At-211 targeted alpha therapy
    // ── Alkali metals — heavy ────────────────────────────────────────────────
    rb,     // Rb+   rubidium; K+ mimic, Rb-82 cardiac PET
    cs,     // Cs+   cesium; K+ mimic, Cs-131 brachytherapy seeds
    // ── Lanthanide gap: Pm ───────────────────────────────────────────────────
    pm,     // Pm3+  promethium; Pm-149 brachytherapy
    // ── Heavy alkaline-earth ─────────────────────────────────────────────────
    ra,     // Ra2+  radium; Ra-223 dichloride (Xofigo) bone metastases
    // ── Actinides ────────────────────────────────────────────────────────────
    ac,     // Ac3+  actinium; Ac-225 targeted alpha therapy
    th,     // Th4+  thorium; Th-227 targeted alpha therapy (Xofigo successor)
    pa,     // Pa4+/5+ protactinium; research
    u,      // U4+/UO2^2+ uranium; environmental toxicology target
    np,     // Np5+  neptunium; research
    pu,     // Pu3+/4+ plutonium; research
    am,     // Am3+  americium; research / environmental
    cm,     // Cm3+  curium; research
    bk,     // Bk3+  berkelium; research
    cf,     // Cf3+  californium; research
    es,     // Es3+  einsteinium; research (AD_TYPE = E)
    fm,     // Fm3+  fermium; research
    // ── Oxidation-state variants ──────────────────────────────────────────────
    fe2,     // Fe²⁺ high-spin oct. (heme, ferritin, Rieske clusters)
    fe3,     // Fe³⁺ high-spin (transferrin, cytochrome, lactoferrin)
    cu1,     // Cu⁺  soft d10 linear/tetrahedral (blue-copper, metallothionein)
    cu2,     // Cu²⁺ Jahn-Teller sq-planar (plastocyanin, laccase, CuZn-SOD)
    cu2_jt,  // Cu²⁺ explicit 4+2 Jahn-Teller mode (equatorial SQ + weak axial JT)
    mn2,     // Mn²⁺ high-spin octahedral (Mn-SOD, ConA, arginase)
    mn3,     // Mn³⁺ Jahn-Teller distorted (Mn-peroxidase, Mn-catalase)
    mn3_jt,  // Mn³⁺ explicit 4+2 Jahn-Teller mode (equatorial MH + elongated axial JT)
    co2,     // Co²⁺ octahedral (carbonic anhydrase, B12 enzymes)
    co3,     // Co³⁺ inert d6 octahedral (cobalamin, NAMI-A complex)
    ni2,     // Ni²⁺ sq-planar/octahedral d8 (urease, NiSOD, acetyl-CoA synthase)
    ni3,     // Ni³⁺ activated Ni-Fe hydrogenase; d7 (stronger N/S pref. vs Ni2+)
    as3,     // As³⁺ pyramidal Cys-binding (Trisenox/ATO AML treatment)
    as5,     // As⁵⁺ arsenate tetrahedral (phosphate mimic, phosphatase inh.)
    sb3,     // Sb³⁺ antimonite (sodium stibogluconate leishmania treatment)
    sb5,     // Sb⁵⁺ meglumine antimoniate (pentavalent antimonial drug)
    v4,      // VO²⁺ vanadyl equatorial (vanadium SOD mimics, antidiabetic)
    v5,      // VO₄³⁻ vanadate (phosphatase inhibitor, phosphate analogue)
    mo4,     // Mo⁴⁺ sulfido (xanthine oxidase, DMSO reductase, nitrogenase)
    mo6,     // Mo⁶⁺ dioxo (sulfite oxidase; hard, oxo-O dominant)
    uo2,     // UO₂²⁺ uranyl equatorial (5-6 eq. ligands, r_eq≈2.44 Å)
    // ── Lanthanide chelate contrast-agent modes ───────────────────────────────
    gd_dtpa, // Gd-DTPA (Magnevist): 9-coord 3N+5O+H₂O; ε(N) enhanced vs free Gd³⁺
    gd_dota, // Gd-DOTA (Dotarem): 9-coord 4N+4O+H₂O; N≈O (macrocyclic pre-org.)
    // ── Aqua-mode variants (_aq): metal with ≥1 inner-sphere H₂O ────────────
    mg_aq,   // Mg²⁺ + H₂O: 5 open sites vs 6 (reduced ε ~20%)
    ca_aq,   // Ca²⁺ + H₂O: 6 open sites vs 7
    fe3_aq,  // Fe³⁺ + H₂O: 5 open sites vs 6 (common in transferrin-like sites)
    mn2_aq,  // Mn²⁺ + H₂O: 5 open sites vs 6 (ConA / arginase)
    co2_aq   // Co²⁺ + H₂O: 5 open sites vs 6
};

// ── Mg2+ ── octahedral, hard Lewis acid, O>N>>S
// r_eq from Bazayeva 2024: Mg-O ≈2.10 Å, Mg-N(His) ≈2.20 Å
inline std::vector<ag4_nbp_override> ag4_mg_mode_overrides() {
    return {
        {"OA",  "Mg",  2.07,    3.5,    12, 6},  // Mg–OA primary (Asp/Glu/Tyr backbone O)
        {"NA",  "Mg",  2.20,    1.5,    12, 6},  // Mg–NA (His)
        {"N",   "Mg",  2.20,    1.5,    12, 6},  // Mg–N
        {"SA",  "Mg",  2.50,    0.5,    12, 6},  // Mg–SA (rare, Cys)
        {"HD",  "Mg",  1.0,     0.0,    12, 6},  // Mg–HD repulsion only
    };
}

// ── Ca2+ ── 7-8 coord, very hard Lewis acid, O-only
// r_eq from Bazayeva 2024: Ca-O ≈2.38 Å (carboxylate), Ca-backbone ≈2.4 Å
inline std::vector<ag4_nbp_override> ag4_ca_mode_overrides() {
    return {
        {"OA",  "Ca",  2.38,    3.0,    12, 6},  // Ca–OA dominant
        {"NA",  "Ca",  2.50,    1.0,    12, 6},  // Ca–NA (weak)
        {"N",   "Ca",  2.50,    1.0,    12, 6},  // Ca–N
        {"HD",  "Ca",  1.2,     0.0,    12, 6},  // Ca–HD repulsion
    };
}

// ── Mn2+ ── octahedral d5, hard Lewis acid, O>N>>S
// r_eq Bazayeva 2024: Mn-His(Nε2)≈2.20 Å, Mn-OA(Asp/Glu)≈2.15 Å
inline std::vector<ag4_nbp_override> ag4_mn_mode_overrides() {
    return {
        {"OA",  "Mn",  2.15,    2.5,    12, 6},  // Mn–OA carboxylate syn
        {"NA",  "Mn",  2.20,    2.0,    12, 6},  // Mn–NA (His Nε2)
        {"N",   "Mn",  2.20,    2.0,    12, 6},  // Mn–N
        {"SA",  "Mn",  2.47,    0.5,    12, 6},  // Mn–SA (rare)
        {"HD",  "Mn",  1.0,     0.0,    12, 6},  // Mn–HD repulsion
    };
}

// ── Fe2+/3+ ── octahedral, heme and non-heme sites
// r_eq Bazayeva 2024: Fe-His(N)≈2.05 Å, Fe-Cys(S)≈2.30 Å, Fe-OA≈2.10 Å, Fe-Tyr(O)≈2.00 Å
inline std::vector<ag4_nbp_override> ag4_fe_mode_overrides() {
    return {
        {"NA",  "Fe",  2.05,    5.0,    12, 6},  // Fe–NA (porphyrin His)
        {"N",   "Fe",  2.05,    5.0,    12, 6},  // Fe–N
        {"OA",  "Fe",  2.05,    3.0,    12, 6},  // Fe–OA (Tyr/Asp/axial O)
        {"SA",  "Fe",  2.30,    2.5,    12, 6},  // Fe–SA (Cys/Met iron-sulfur)
        {"HD",  "Fe",  1.0,     0.0,    12, 6},  // Fe–HD repulsion
    };
}

// ── Co2+ ── tetrahedral/octahedral d7, borderline Lewis acid
// r_eq Harding 2006 CSD: Co-N≈2.07 Å, Co-O≈2.08 Å, Co-S≈2.25 Å
inline std::vector<ag4_nbp_override> ag4_co_mode_overrides() {
    return {
        {"NA",  "Co",  2.07,    3.0,    12, 6},  // Co–NA
        {"N",   "Co",  2.07,    3.0,    12, 6},  // Co–N
        {"OA",  "Co",  2.08,    2.0,    12, 6},  // Co–OA
        {"SA",  "Co",  2.25,    2.5,    12, 6},  // Co–SA
        {"HD",  "Co",  1.0,     0.0,    12, 6},  // Co–HD repulsion
    };
}

// ── Ni2+ ── sq-planar/octahedral d8, borderline-soft
// r_eq Bazayeva 2024: Ni-His(Nε2)≈2.10 Å, Ni-OA≈2.06 Å, Ni-SA≈2.15 Å
inline std::vector<ag4_nbp_override> ag4_ni_mode_overrides() {
    return {
        {"NA",  "Ni",  2.10,    3.5,    12, 6},  // Ni–NA (His; dominant in Ni-SOD)
        {"N",   "Ni",  2.10,    3.5,    12, 6},  // Ni–N
        {"OA",  "Ni",  2.06,    1.5,    12, 6},  // Ni–OA (Asp/backbone O)
        {"SA",  "Ni",  2.15,    3.0,    12, 6},  // Ni–SA (Cys; urease)
        {"HD",  "Ni",  1.0,     0.0,    12, 6},  // Ni–HD repulsion
    };
}

// ── Cu2+ ── sq-planar + axial Jahn-Teller d9, borderline-soft
// r_eq Bazayeva 2024: Cu-His(N)≈2.00 Å, Cu-Cys(S)≈2.20 Å, Cu-Met(S)≈2.50 Å mapped→SA
inline std::vector<ag4_nbp_override> ag4_cu_mode_overrides() {
    return {
        {"NA",  "Cu",  2.00,    3.5,    12, 6},  // Cu–NA (His; strong)
        {"N",   "Cu",  2.00,    3.5,    12, 6},  // Cu–N
        {"SA",  "Cu",  2.20,    4.5,    12, 6},  // Cu–SA (Cys thiolate; covalent-like)
        {"OA",  "Cu",  2.25,    1.5,    12, 6},  // Cu–OA (Asp/Glu backbone)
        {"HD",  "Cu",  1.0,     0.0,    12, 6},  // Cu–HD repulsion
    };
}

// ── Pt2+ ── sq-planar d8, strong trans-influence, cisplatin family
// r_eq Lippard/Jamieson: Pt-N≈2.03 Å, Pt-S≈2.30 Å, Pt-O≈2.10 Å
inline std::vector<ag4_nbp_override> ag4_pt_mode_overrides() {
    return {
        {"NA",  "Pt",  2.03,    6.0,    12, 6},  // Pt–NA (imidazole, DNA-N7)
        {"N",   "Pt",  2.03,    6.0,    12, 6},  // Pt–N
        {"SA",  "Pt",  2.30,    5.0,    12, 6},  // Pt–SA (Cys/Met; strong)
        {"OA",  "Pt",  2.10,    2.0,    12, 6},  // Pt–OA
        {"HD",  "Pt",  1.0,     0.0,    12, 6},  // Pt–HD repulsion
    };
}

// ── Pd2+ ── sq-planar d8, softer than Pt, MetalDock-parameterised
// r_eq CSD: Pd-N≈2.02 Å, Pd-S≈2.28 Å; ε from MetalDock MC optimisation
inline std::vector<ag4_nbp_override> ag4_pd_mode_overrides() {
    return {
        {"NA",  "Pd",  2.02,    4.0,    12, 6},  // Pd–NA
        {"N",   "Pd",  2.02,    4.0,    12, 6},  // Pd–N
        {"SA",  "Pd",  2.28,    4.5,    12, 6},  // Pd–SA (strong)
        {"OA",  "Pd",  2.02,    1.0,    12, 6},  // Pd–OA (weak)
        {"HD",  "Pd",  1.0,     0.0,    12, 6},  // Pd–HD repulsion
    };
}

// ── Ru3+ ── octahedral d5, NAMI-A / RAPTA anticancer; MetalDock-parameterised
// r_eq from NAMI-A crystal + Bazayeva: Ru-N≈2.12 Å, Ru-O≈2.05 Å, Ru-S≈2.34 Å
inline std::vector<ag4_nbp_override> ag4_ru_mode_overrides() {
    return {
        {"NA",  "Ru",  2.09,    6.0,    12, 6},  // Ru–NA (His; MetalDock MC-opt: r↓0.03, ε↑×1.7)
        {"N",   "Ru",  2.09,    6.0,    12, 6},  // Ru–N
        {"OA",  "Ru",  2.06,    4.0,    12, 6},  // Ru–OA (MetalDock: ε↑×2.0)
        {"SA",  "Ru",  2.34,    4.5,    12, 6},  // Ru–SA (Met; MetalDock: ε↑×1.5)
        {"HD",  "Ru",  1.0,     0.0,    12, 6},  // Ru–HD repulsion
    };
}

// ── Ir3+ ── octahedral d6, strong ligand-field organometallics
// r_eq CSD Ir complexes: Ir-N≈2.05 Å, Ir-O≈2.02 Å, Ir-S≈2.27 Å
inline std::vector<ag4_nbp_override> ag4_ir_mode_overrides() {
    return {
        {"NA",  "Ir",  2.05,    7.0,    12, 6},  // Ir–NA (MetalDock MC-opt: ε⇑×1.75)
        {"N",   "Ir",  2.05,    7.0,    12, 6},  // Ir–N
        {"OA",  "Ir",  2.02,    4.5,    12, 6},  // Ir–OA (MetalDock: ε⇑×2.25)
        {"SA",  "Ir",  2.27,    5.5,    12, 6},  // Ir–SA (MetalDock: ε⇑×1.57)
        {"HD",  "Ir",  1.0,     0.0,    12, 6},  // Ir–HD repulsion
    };
}

// ── Au+ ── linear d10, auranofin / gold drugs; extreme S-philic (HSAB soft acid)
// r_eq CSD Au-thiolate: Au-S≈2.29 Å, Au-N≈2.09 Å, Au-O≈2.5 Å (rare)
inline std::vector<ag4_nbp_override> ag4_au_mode_overrides() {
    return {
        {"SA",  "Au",  2.29,   11.0,    12, 6},  // Au–SA (MetalDock MC-opt: extreme S-philic, ε⇑×1.83)
        {"NA",  "Au",  2.09,    3.5,    12, 6},  // Au–NA (His; MetalDock: ε⇑×1.4)
        {"N",   "Au",  2.09,    3.5,    12, 6},  // Au–N
        {"OA",  "Au",  2.50,    0.5,    12, 6},  // Au–OA (very weak)
        {"HD",  "Au",  1.0,     0.0,    12, 6},  // Au–HD repulsion
    };
}

// ── Cd2+ ── tetrahedral d10, Zn-mimic; prefers S>N>O (softer than Zn)
// r_eq Harding 2001: Cd-S≈2.52 Å, Cd-N(His)≈2.27 Å, Cd-O≈2.35 Å
inline std::vector<ag4_nbp_override> ag4_cd_mode_overrides() {
    return {
        {"SA",  "Cd",  2.52,    3.5,    12, 6},  // Cd–SA (metallothionein Cys)
        {"NA",  "Cd",  2.27,    2.0,    12, 6},  // Cd–NA (His)
        {"N",   "Cd",  2.27,    2.0,    12, 6},  // Cd–N
        {"OA",  "Cd",  2.35,    1.0,    12, 6},  // Cd–OA (weak)
        {"HD",  "Cd",  1.0,     0.0,    12, 6},  // Cd–HD repulsion
    };
}

// ── Hg2+ ── linear d10, most S-philic metal; metallothionein inhibitor
// r_eq CSD: Hg-Cys(S)≈2.34 Å, Hg-N(His)≈2.04 Å, Hg-O≈2.4 Å
inline std::vector<ag4_nbp_override> ag4_hg_mode_overrides() {
    return {
        {"SA",  "Hg",  2.34,    7.0,    12, 6},  // Hg–SA (extreme; Hg-Cys bond ~covalent)
        {"NA",  "Hg",  2.04,    2.0,    12, 6},  // Hg–NA
        {"N",   "Hg",  2.04,    2.0,    12, 6},  // Hg–N
        {"OA",  "Hg",  2.40,    0.3,    12, 6},  // Hg–OA (very weak)
        {"HD",  "Hg",  1.0,     0.0,    12, 6},  // Hg–HD repulsion
    };
}

// ── Na+ ── octahedral, hard Lewis acid; exclusively O-coordinated in proteins
// r_eq Bazayeva 2024: Na-backbone-O≈2.35 Å, Na-Asp/Glu-O≈2.40 Å
inline std::vector<ag4_nbp_override> ag4_na_ion_mode_overrides() {
    return {
        {"OA",  "Na",  2.35,    2.0,    12, 6},  // Na–OA (dominant)
        {"NA",  "Na",  2.50,    0.5,    12, 6},  // Na–NA (very rare)
        {"N",   "Na",  2.50,    0.5,    12, 6},  // Na–N
        {"SA",  "Na",  2.80,    0.1,    12, 6},  // Na–SA (essentially zero)
        {"HD",  "Na",  1.2,     0.0,    12, 6},  // Na–HD repulsion
    };
}

// ── K+ ── 7-8 coord, very large radius, hard but weakly bound
// r_eq Bazayeva 2024: K-O≈2.70 Å (backbone/carboxylate); very broad distributions
inline std::vector<ag4_nbp_override> ag4_k_ion_mode_overrides() {
    return {
        {"OA",  "K",   2.70,    1.5,    12, 6},  // K–OA
        {"NA",  "K",   2.80,    0.3,    12, 6},  // K–NA (very weak)
        {"N",   "K",   2.80,    0.3,    12, 6},  // K–N
        {"SA",  "K",   3.00,    0.1,    12, 6},  // K–SA (essentially zero)
        {"HD",  "K",   1.5,     0.0,    12, 6},  // K–HD repulsion
    };
}

// ════════════════════════════════════════════════════════════════════════════
// ── Early transition metals (hard Lewis acids) ───────────────────────────
// ════════════════════════════════════════════════════════════════════════════

// ── V3+/4+ ── octahedral, borderline-hard; insulin-mimetic VOSO4
// r_eq CSD: V-O ≈2.00 Å, V-N(His) ≈2.17 Å, V-S(Cys) ≈2.50 Å (rare)
inline std::vector<ag4_nbp_override> ag4_v_mode_overrides() {
    return {
        {"OA",  "V",   2.00,    3.0,    12, 6},  // V–OA (carboxylate/Tyr dominant)
        {"NA",  "V",   2.17,    2.0,    12, 6},  // V–NA (His)
        {"N",   "V",   2.17,    2.0,    12, 6},  // V–N
        {"SA",  "V",   2.50,    0.8,    12, 6},  // V–SA (Cys rare)
        {"HD",  "V",   1.0,     0.0,    12, 6},  // V–HD repulsion
    };
}

// ── Cr3+ ── octahedral d3, very inert; glucose tolerance factor
// r_eq CSD Cr(III): Cr-O ≈1.97 Å, Cr-N ≈2.07 Å; hard Lewis acid
inline std::vector<ag4_nbp_override> ag4_cr_mode_overrides() {
    return {
        {"OA",  "Cr1", 1.97,    3.5,    12, 6},  // Cr–OA (dominant; Asp/Glu/acetate)
        {"NA",  "Cr1", 2.07,    2.5,    12, 6},  // Cr–NA (His)
        {"N",   "Cr1", 2.07,    2.5,    12, 6},  // Cr–N
        {"SA",  "Cr1", 2.40,    0.5,    12, 6},  // Cr–SA (rare)
        {"HD",  "Cr1", 1.0,     0.0,    12, 6},  // Cr–HD repulsion
    };
}

// ── Ti4+ ── octahedral d0, hard Lewis acid; titanocene anticancer
// r_eq CSD Ti(IV): Ti-O ≈1.95 Å, Ti-N ≈2.15 Å; very hard
inline std::vector<ag4_nbp_override> ag4_ti_mode_overrides() {
    return {
        {"OA",  "Ti",  1.95,    3.5,    12, 6},  // Ti–OA (dominant; Asp/Glu)
        {"NA",  "Ti",  2.15,    1.5,    12, 6},  // Ti–NA (His; weaker than O)
        {"N",   "Ti",  2.15,    1.5,    12, 6},  // Ti–N
        {"SA",  "Ti",  2.50,    0.3,    12, 6},  // Ti–SA (very rare)
        {"HD",  "Ti",  1.0,     0.0,    12, 6},  // Ti–HD repulsion
    };
}

// ── Sc3+ ── octahedral d0, very hard; Sc-44 PET / enzyme cofactor
// r_eq CSD: Sc-O ≈2.11 Å, Sc-N ≈2.23 Å
inline std::vector<ag4_nbp_override> ag4_sc_mode_overrides() {
    return {
        {"OA",  "Sc",  2.11,    3.5,    12, 6},  // Sc–OA
        {"NA",  "Sc",  2.23,    2.0,    12, 6},  // Sc–NA
        {"N",   "Sc",  2.23,    2.0,    12, 6},  // Sc–N
        {"SA",  "Sc",  2.55,    0.3,    12, 6},  // Sc–SA (essentially none)
        {"HD",  "Sc",  1.0,     0.0,    12, 6},  // Sc–HD repulsion
    };
}

// ── Y3+ ── octahedral d0, hard; 90Y radiotherapy chelates
// r_eq CSD: Y-O ≈2.35 Å, Y-N ≈2.47 Å
inline std::vector<ag4_nbp_override> ag4_y_ion_mode_overrides() {
    return {
        {"OA",  "Y",   2.35,    3.0,    12, 6},  // Y–OA (dominant)
        {"NA",  "Y",   2.47,    1.5,    12, 6},  // Y–NA
        {"N",   "Y",   2.47,    1.5,    12, 6},  // Y–N
        {"SA",  "Y",   2.70,    0.2,    12, 6},  // Y–SA (negligible)
        {"HD",  "Y",   1.2,     0.0,    12, 6},  // Y–HD repulsion
    };
}

// ── Zr4+ ── 8-coord, hard; Zr-89 PET chelates (DFO)
// r_eq CSD: Zr-O ≈2.21 Å, Zr-N ≈2.36 Å
inline std::vector<ag4_nbp_override> ag4_zr_mode_overrides() {
    return {
        {"OA",  "Zr",  2.21,    3.5,    12, 6},  // Zr–OA (dominant)
        {"NA",  "Zr",  2.36,    2.0,    12, 6},  // Zr–NA
        {"N",   "Zr",  2.36,    2.0,    12, 6},  // Zr–N
        {"SA",  "Zr",  2.65,    0.2,    12, 6},  // Zr–SA (negligible)
        {"HD",  "Zr",  1.2,     0.0,    12, 6},  // Zr–HD repulsion
    };
}

// ── Nb5+ ── octahedral d0, hard Lewis acid
// r_eq CSD: Nb-O ≈1.98 Å, Nb-N ≈2.10 Å
inline std::vector<ag4_nbp_override> ag4_nb_mode_overrides() {
    return {
        {"OA",  "Nb",  1.98,    3.5,    12, 6},  // Nb–OA
        {"NA",  "Nb",  2.10,    2.0,    12, 6},  // Nb–NA
        {"N",   "Nb",  2.10,    2.0,    12, 6},  // Nb–N
        {"SA",  "Nb",  2.48,    0.3,    12, 6},  // Nb–SA (rare)
        {"HD",  "Nb",  1.0,     0.0,    12, 6},  // Nb–HD repulsion
    };
}

// ── Hf4+ ── 8-coord, hard; hafnium nanoparticle radiotherapy (NBTXR3)
// r_eq CSD: Hf-O ≈2.22 Å, Hf-N ≈2.38 Å; similar to Zr
inline std::vector<ag4_nbp_override> ag4_hf_mode_overrides() {
    return {
        {"OA",  "Hf",  2.22,    3.5,    12, 6},  // Hf–OA (dominant)
        {"NA",  "Hf",  2.38,    2.0,    12, 6},  // Hf–NA
        {"N",   "Hf",  2.38,    2.0,    12, 6},  // Hf–N
        {"SA",  "Hf",  2.65,    0.2,    12, 6},  // Hf–SA (negligible)
        {"HD",  "Hf",  1.2,     0.0,    12, 6},  // Hf–HD repulsion
    };
}

// ── Ta5+ ── octahedral d0, hard; tantalum complexes
// r_eq CSD: Ta-O ≈1.97 Å, Ta-N ≈2.12 Å
inline std::vector<ag4_nbp_override> ag4_ta_mode_overrides() {
    return {
        {"OA",  "Ta",  1.97,    3.5,    12, 6},  // Ta–OA
        {"NA",  "Ta",  2.12,    2.0,    12, 6},  // Ta–NA
        {"N",   "Ta",  2.12,    2.0,    12, 6},  // Ta–N
        {"SA",  "Ta",  2.50,    0.3,    12, 6},  // Ta–SA (rare)
        {"HD",  "Ta",  1.0,     0.0,    12, 6},  // Ta–HD repulsion
    };
}

// ── W6+ (Tg AD_TYPE) ── octahedral d0; tungstoenzymes (aldehyde oxidoreductase)
// r_eq CSD: W-O ≈1.95 Å, W-N ≈2.20 Å; very hard
inline std::vector<ag4_nbp_override> ag4_tg_mode_overrides() {
    return {
        {"OA",  "Tg",  1.95,    4.0,    12, 6},  // W–OA (dominant; W=O terminal)
        {"NA",  "Tg",  2.20,    1.5,    12, 6},  // W–NA (His ligand)
        {"N",   "Tg",  2.20,    1.5,    12, 6},  // W–N
        {"SA",  "Tg",  2.60,    0.5,    12, 6},  // W–SA (Cys in some W-enzymes)
        {"HD",  "Tg",  1.0,     0.0,    12, 6},  // W–HD repulsion
    };
}

// ── Mo4+/6+ ── octahedral; xanthine oxidase, sulfite oxidase, nitrogenase
// r_eq CSD: Mo-O(oxo) ≈1.70 Å, Mo-O(carboxylate) ≈2.08 Å; here use ligand-binding O
// r_eq Bazayeva 2024: Mo-N ≈2.19 Å, Mo-S(Cys) ≈2.49 Å
inline std::vector<ag4_nbp_override> ag4_mo_mode_overrides() {
    return {
        {"OA",  "Mo",  2.08,    2.5,    12, 6},  // Mo–OA (Asp/Glu; non-oxo)
        {"NA",  "Mo",  2.19,    2.0,    12, 6},  // Mo–NA (His)
        {"N",   "Mo",  2.19,    2.0,    12, 6},  // Mo–N
        {"SA",  "Mo",  2.49,    2.5,    12, 6},  // Mo–SA (Cys/molybdopterin S)
        {"HD",  "Mo",  1.0,     0.0,    12, 6},  // Mo–HD repulsion
    };
}

// ════════════════════════════════════════════════════════════════════════════
// ── 2nd/3rd row TMs — borderline to soft ────────────────────────────────
// ════════════════════════════════════════════════════════════════════════════

// ── Rh3+ ── octahedral d6, borderline-soft; anticancer organometallics
// r_eq CSD: Rh-N ≈2.04 Å, Rh-O ≈2.04 Å, Rh-S ≈2.35 Å; N ≥ O > S priority
inline std::vector<ag4_nbp_override> ag4_rh_mode_overrides() {
    return {
        {"NA",  "Rh",  2.04,    4.0,    12, 6},  // Rh–NA (His; dominant in Rh drugs)
        {"N",   "Rh",  2.04,    4.0,    12, 6},  // Rh–N
        {"OA",  "Rh",  2.04,    2.0,    12, 6},  // Rh–OA (carboxylate)
        {"SA",  "Rh",  2.35,    3.0,    12, 6},  // Rh–SA (Cys/Met)
        {"HD",  "Rh",  1.0,     0.0,    12, 6},  // Rh–HD repulsion
    };
}

// ── Ag+ ── linear/tetrahedral d10, moderate S-philic; antimicrobial
// r_eq CSD: Ag-S ≈2.46 Å, Ag-N(His) ≈2.20 Å, Ag-O ≈2.35 Å
inline std::vector<ag4_nbp_override> ag4_ag_mode_overrides() {
    return {
        {"SA",  "Ag",  2.46,    4.0,    12, 6},  // Ag–SA (Cys; primary in antimicrobial)
        {"NA",  "Ag",  2.20,    2.5,    12, 6},  // Ag–NA (His)
        {"N",   "Ag",  2.20,    2.5,    12, 6},  // Ag–N
        {"OA",  "Ag",  2.35,    0.8,    12, 6},  // Ag–OA (weak)
        {"HD",  "Ag",  1.0,     0.0,    12, 6},  // Ag–HD repulsion
    };
}

// ── Tc5+ ── octahedral d0; 99mTc SPECT (most widely used diagnostic radionuclide)
// r_eq CSD Tc(V): Tc-N ≈2.02 Å, Tc-O ≈2.05 Å, Tc-S ≈2.34 Å
inline std::vector<ag4_nbp_override> ag4_tc_mode_overrides() {
    return {
        {"NA",  "Tc",  2.02,    3.5,    12, 6},  // Tc–NA (His; HYNIC coordination)
        {"N",   "Tc",  2.02,    3.5,    12, 6},  // Tc–N
        {"OA",  "Tc",  2.05,    2.0,    12, 6},  // Tc–OA
        {"SA",  "Tc",  2.34,    3.0,    12, 6},  // Tc–SA (Cys/thiolate)
        {"HD",  "Tc",  1.0,     0.0,    12, 6},  // Tc–HD repulsion
    };
}

// ── Re5+ ── octahedral; 186Re/188Re radiotherapy; congener of Mn/Tc
// r_eq CSD Re(V): Re-N ≈2.10 Å, Re-O ≈2.05 Å, Re-S ≈2.40 Å
inline std::vector<ag4_nbp_override> ag4_re_mode_overrides() {
    return {
        {"NA",  "Re",  2.10,    3.5,    12, 6},  // Re–NA
        {"N",   "Re",  2.10,    3.5,    12, 6},  // Re–N
        {"OA",  "Re",  2.05,    2.0,    12, 6},  // Re–OA
        {"SA",  "Re",  2.40,    3.0,    12, 6},  // Re–SA
        {"HD",  "Re",  1.0,     0.0,    12, 6},  // Re–HD repulsion
    };
}

// ── Os2+/3+ ── octahedral d6; osmium arene anticancer (RAPTA-C family)
// r_eq CSD Os complexes: Os-N ≈2.07 Å, Os-O ≈2.04 Å, Os-S ≈2.30 Å
inline std::vector<ag4_nbp_override> ag4_os_mode_overrides() {
    return {
        {"NA",  "Os",  2.07,    4.0,    12, 6},  // Os–NA (His)
        {"N",   "Os",  2.07,    4.0,    12, 6},  // Os–N
        {"SA",  "Os",  2.30,    3.5,    12, 6},  // Os–SA (Cys/Met)
        {"OA",  "Os",  2.04,    1.5,    12, 6},  // Os–OA (weaker than N)
        {"HD",  "Os",  1.0,     0.0,    12, 6},  // Os–HD repulsion
    };
}

// ════════════════════════════════════════════════════════════════════════════
// ── Post-transition metals / metalloids ─────────────────────────────────
// ════════════════════════════════════════════════════════════════════════════

// ── Ga3+ ── tetrahedral/octahedral d10, borderline-hard; 68Ga PET, antimicrobial
// r_eq CSD: Ga-O ≈1.97 Å, Ga-N ≈2.09 Å, Ga-S ≈2.35 Å; O ≥ N > S
inline std::vector<ag4_nbp_override> ag4_ga_mode_overrides() {
    return {
        {"OA",  "Ga",  1.97,    3.5,    12, 6},  // Ga–OA (Asp/Glu; dominant)
        {"NA",  "Ga",  2.09,    3.0,    12, 6},  // Ga–NA (His; DFO-type chelates)
        {"N",   "Ga",  2.09,    3.0,    12, 6},  // Ga–N
        {"SA",  "Ga",  2.35,    0.8,    12, 6},  // Ga–SA (Cys; weaker)
        {"HD",  "Ga",  1.0,     0.0,    12, 6},  // Ga–HD repulsion
    };
}

// ── In3+ ── octahedral d10, borderline; 111In SPECT / 113In
// r_eq CSD: In-O ≈2.16 Å, In-N ≈2.28 Å, In-S ≈2.53 Å
inline std::vector<ag4_nbp_override> ag4_in_ion_mode_overrides() {
    return {
        {"OA",  "In",  2.16,    3.0,    12, 6},  // In–OA (dominant; DTPA/DOTA)
        {"NA",  "In",  2.28,    2.5,    12, 6},  // In–NA (His; chelate N)
        {"N",   "In",  2.28,    2.5,    12, 6},  // In–N
        {"SA",  "In",  2.53,    1.5,    12, 6},  // In–SA (Cys; modest)
        {"HD",  "In",  1.0,     0.0,    12, 6},  // In–HD repulsion
    };
}

// ── Sn4+ ── octahedral d10; organotin anticancer/antifungal (tributyltin)
// r_eq CSD Sn(IV): Sn-O ≈2.15 Å, Sn-N ≈2.20 Å, Sn-S ≈2.50 Å
inline std::vector<ag4_nbp_override> ag4_sn_mode_overrides() {
    return {
        {"OA",  "Sn",  2.15,    2.5,    12, 6},  // Sn–OA (carboxylate)
        {"NA",  "Sn",  2.20,    2.0,    12, 6},  // Sn–NA
        {"N",   "Sn",  2.20,    2.0,    12, 6},  // Sn–N
        {"SA",  "Sn",  2.50,    2.0,    12, 6},  // Sn–SA (Cys; tri-organotin)
        {"HD",  "Sn",  1.0,     0.0,    12, 6},  // Sn–HD repulsion
    };
}

// ── Sb3+ ── octahedral d10; antiparasitic drugs (sodium stibogluconate/Pentostam)
// r_eq CSD Sb(III): Sb-O ≈2.10 Å, Sb-N ≈2.25 Å, Sb-S ≈2.45 Å
inline std::vector<ag4_nbp_override> ag4_sb_mode_overrides() {
    return {
        {"SA",  "Sb",  2.45,    3.5,    12, 6},  // Sb–SA (Cys/Sb-S; thiol interaction)
        {"OA",  "Sb",  2.10,    2.5,    12, 6},  // Sb–OA (glycolate in Pentostam)
        {"NA",  "Sb",  2.25,    1.5,    12, 6},  // Sb–NA
        {"N",   "Sb",  2.25,    1.5,    12, 6},  // Sb–N
        {"HD",  "Sb",  1.0,     0.0,    12, 6},  // Sb–HD repulsion
    };
}

// ── Bi3+ ── irregular 6-8 coord, d10; anti-H.pylori bismuth subcitrate
// r_eq CSD: Bi-O ≈2.30 Å, Bi-S(Cys) ≈2.52 Å, Bi-N ≈2.40 Å
inline std::vector<ag4_nbp_override> ag4_bi_mode_overrides() {
    return {
        {"SA",  "Bi",  2.52,    4.0,    12, 6},  // Bi–SA (Cys; metallothionein binding)
        {"OA",  "Bi",  2.30,    2.5,    12, 6},  // Bi–OA (citrate oxygen)
        {"NA",  "Bi",  2.40,    1.0,    12, 6},  // Bi–NA (His; secondary)
        {"N",   "Bi",  2.40,    1.0,    12, 6},  // Bi–N
        {"HD",  "Bi",  1.0,     0.0,    12, 6},  // Bi–HD repulsion
    };
}

// ── Tl+ ── soft d10, large; 201Tl cardiac SPECT / toxic
// r_eq CSD: Tl-S ≈2.67 Å, Tl-O ≈2.71 Å, Tl-N ≈2.72 Å
inline std::vector<ag4_nbp_override> ag4_tl_mode_overrides() {
    return {
        {"SA",  "Tl",  2.67,    3.0,    12, 6},  // Tl–SA (soft S-philic)
        {"OA",  "Tl",  2.71,    1.0,    12, 6},  // Tl–OA
        {"NA",  "Tl",  2.72,    0.8,    12, 6},  // Tl–NA
        {"N",   "Tl",  2.72,    0.8,    12, 6},  // Tl–N
        {"HD",  "Tl",  1.5,     0.0,    12, 6},  // Tl–HD repulsion
    };
}

// ── Pb2+ ── hemispherical 4-8 coord d10; toxicology / environmental target
// r_eq CSD Pb(II): Pb-S ≈2.73 Å, Pb-O ≈2.58 Å, Pb-N ≈2.64 Å
inline std::vector<ag4_nbp_override> ag4_pb_mode_overrides() {
    return {
        {"SA",  "Pb",  2.73,    3.0,    12, 6},  // Pb–SA (Cys/Met; soft character)
        {"OA",  "Pb",  2.58,    1.5,    12, 6},  // Pb–OA (Asp/Glu)
        {"NA",  "Pb",  2.64,    1.2,    12, 6},  // Pb–NA (His)
        {"N",   "Pb",  2.64,    1.2,    12, 6},  // Pb–N
        {"HD",  "Pb",  1.5,     0.0,    12, 6},  // Pb–HD repulsion
    };
}

// ════════════════════════════════════════════════════════════════════════════
// ── s-block metals ───────────────────────────────────────────────────────
// ════════════════════════════════════════════════════════════════════════════

// ── Li+ ── tetrahedral, very hard; lithium carbonate (GSK-3β inhibition)
// r_eq CSD: Li-O ≈2.00 Å, Li-N ≈2.13 Å; charge-charge dominant
inline std::vector<ag4_nbp_override> ag4_li_mode_overrides() {
    return {
        {"OA",  "Li",  2.00,    2.5,    12, 6},  // Li–OA (dominant)
        {"NA",  "Li",  2.13,    1.0,    12, 6},  // Li–NA (rare)
        {"N",   "Li",  2.13,    1.0,    12, 6},  // Li–N
        {"SA",  "Li",  2.70,    0.1,    12, 6},  // Li–SA (negligible)
        {"HD",  "Li",  1.0,     0.0,    12, 6},  // Li–HD repulsion
    };
}

// ── Al3+ ── tetrahedral, extremely hard; AlF4-/AlF3 phosphatase TS analogs
// r_eq CSD Al(III): Al-O ≈1.88 Å, Al-N ≈2.05 Å; very short bonds
inline std::vector<ag4_nbp_override> ag4_al_mode_overrides() {
    return {
        {"OA",  "Al",  1.88,    4.5,    12, 6},  // Al–OA (very strong; AlF4- binding)
        {"NA",  "Al",  2.05,    2.0,    12, 6},  // Al–NA
        {"N",   "Al",  2.05,    2.0,    12, 6},  // Al–N
        {"SA",  "Al",  2.37,    0.2,    12, 6},  // Al–SA (negligible)
        {"HD",  "Al",  1.0,     0.0,    12, 6},  // Al–HD repulsion
    };
}

// ── Sr2+ ── 8-coord, hard; strontium ranelate (bone disease drugs)
// r_eq Bazayeva 2024: Sr-O ≈2.59 Å, Sr-N ≈2.72 Å
inline std::vector<ag4_nbp_override> ag4_sr_mode_overrides() {
    return {
        {"OA",  "Sr",  2.59,    2.5,    12, 6},  // Sr–OA (dominant; Asp/Glu/phosphate)
        {"NA",  "Sr",  2.72,    0.8,    12, 6},  // Sr–NA (rare)
        {"N",   "Sr",  2.72,    0.8,    12, 6},  // Sr–N
        {"SA",  "Sr",  3.00,    0.1,    12, 6},  // Sr–SA (negligible)
        {"HD",  "Sr",  1.5,     0.0,    12, 6},  // Sr–HD repulsion
    };
}

// ── Ba2+ ── 9-coord, very hard, largest alkaline-earth
// r_eq Bazayeva 2024: Ba-O ≈2.74 Å, Ba-N ≈2.85 Å
inline std::vector<ag4_nbp_override> ag4_ba_mode_overrides() {
    return {
        {"OA",  "Ba",  2.74,    2.0,    12, 6},  // Ba–OA (dominant)
        {"NA",  "Ba",  2.85,    0.5,    12, 6},  // Ba–NA (very rare)
        {"N",   "Ba",  2.85,    0.5,    12, 6},  // Ba–N
        {"SA",  "Ba",  3.20,    0.1,    12, 6},  // Ba–SA (negligible)
        {"HD",  "Ba",  1.5,     0.0,    12, 6},  // Ba–HD repulsion
    };
}

// ════════════════════════════════════════════════════════════════════════════
// ── Lanthanides (hard O-preferring, lanthanide contraction r_eq trend) ──
// r_eq sources: Bazayeva 2024 (PMC11066882), CSD lanthanide surveys
// All: O > N >> S priority; ε(OA)=3.0, ε(NA)=1.5, ε(SA)=0.2 (hard Lewis acids)
// ════════════════════════════════════════════════════════════════════════════
inline std::vector<ag4_nbp_override> ag4_la_mode_overrides() {
    return {{"OA","La",2.55,3.0,12,6},{"NA","La",2.72,1.5,12,6},{"N","La",2.72,1.5,12,6},
            {"SA","La",2.95,0.2,12,6},{"HD","La",1.5,0.0,12,6}};
}
inline std::vector<ag4_nbp_override> ag4_ce_mode_overrides() {
    return {{"OA","Ce",2.46,3.0,12,6},{"NA","Ce",2.63,1.5,12,6},{"N","Ce",2.63,1.5,12,6},
            {"SA","Ce",2.88,0.2,12,6},{"HD","Ce",1.5,0.0,12,6}};
}
inline std::vector<ag4_nbp_override> ag4_pr_mode_overrides() {
    return {{"OA","Pr",2.45,3.0,12,6},{"NA","Pr",2.61,1.5,12,6},{"N","Pr",2.61,1.5,12,6},
            {"SA","Pr",2.85,0.2,12,6},{"HD","Pr",1.5,0.0,12,6}};
}
inline std::vector<ag4_nbp_override> ag4_nd_mode_overrides() {
    return {{"OA","Nd",2.44,3.0,12,6},{"NA","Nd",2.60,1.5,12,6},{"N","Nd",2.60,1.5,12,6},
            {"SA","Nd",2.83,0.2,12,6},{"HD","Nd",1.5,0.0,12,6}};
}
inline std::vector<ag4_nbp_override> ag4_sm_mode_overrides() {
    return {{"OA","Sm",2.41,3.0,12,6},{"NA","Sm",2.57,1.5,12,6},{"N","Sm",2.57,1.5,12,6},
            {"SA","Sm",2.80,0.2,12,6},{"HD","Sm",1.5,0.0,12,6}};
}
inline std::vector<ag4_nbp_override> ag4_eu_mode_overrides() {
    return {{"OA","Eu",2.40,3.0,12,6},{"NA","Eu",2.56,1.5,12,6},{"N","Eu",2.56,1.5,12,6},
            {"SA","Eu",2.78,0.2,12,6},{"HD","Eu",1.5,0.0,12,6}};
}
inline std::vector<ag4_nbp_override> ag4_gd_mode_overrides() {
    return {{"OA","Gd",2.40,3.0,12,6},{"NA","Gd",2.55,1.5,12,6},{"N","Gd",2.55,1.5,12,6},
            {"SA","Gd",2.78,0.2,12,6},{"HD","Gd",1.5,0.0,12,6}};
}
inline std::vector<ag4_nbp_override> ag4_tb_mode_overrides() {
    return {{"OA","Tb",2.38,3.0,12,6},{"NA","Tb",2.53,1.5,12,6},{"N","Tb",2.53,1.5,12,6},
            {"SA","Tb",2.75,0.2,12,6},{"HD","Tb",1.5,0.0,12,6}};
}
inline std::vector<ag4_nbp_override> ag4_dy_mode_overrides() {
    return {{"OA","Dy",2.36,3.0,12,6},{"NA","Dy",2.51,1.5,12,6},{"N","Dy",2.51,1.5,12,6},
            {"SA","Dy",2.73,0.2,12,6},{"HD","Dy",1.5,0.0,12,6}};
}
inline std::vector<ag4_nbp_override> ag4_ho_mode_overrides() {
    return {{"OA","Ho",2.34,3.0,12,6},{"NA","Ho",2.49,1.5,12,6},{"N","Ho",2.49,1.5,12,6},
            {"SA","Ho",2.70,0.2,12,6},{"HD","Ho",1.5,0.0,12,6}};
}
inline std::vector<ag4_nbp_override> ag4_er_mode_overrides() {
    return {{"OA","Er",2.33,3.0,12,6},{"NA","Er",2.48,1.5,12,6},{"N","Er",2.48,1.5,12,6},
            {"SA","Er",2.68,0.2,12,6},{"HD","Er",1.5,0.0,12,6}};
}
inline std::vector<ag4_nbp_override> ag4_tm_mode_overrides() {
    return {{"OA","Tm",2.32,3.0,12,6},{"NA","Tm",2.47,1.5,12,6},{"N","Tm",2.47,1.5,12,6},
            {"SA","Tm",2.65,0.2,12,6},{"HD","Tm",1.5,0.0,12,6}};
}
inline std::vector<ag4_nbp_override> ag4_yb_mode_overrides() {
    return {{"OA","Yb",2.28,3.0,12,6},{"NA","Yb",2.44,1.5,12,6},{"N","Yb",2.44,1.5,12,6},
            {"SA","Yb",2.63,0.2,12,6},{"HD","Yb",1.5,0.0,12,6}};
}
inline std::vector<ag4_nbp_override> ag4_lu_mode_overrides() {
    return {{"OA","Lu",2.24,3.0,12,6},{"NA","Lu",2.40,1.5,12,6},{"N","Lu",2.40,1.5,12,6},
            {"SA","Lu",2.60,0.2,12,6},{"HD","Lu",1.5,0.0,12,6}};
}

// ════════════════════════════════════════════════════════════════════════════
// ── Metalloids with pharmacological interest ─────────────────────────────
// ════════════════════════════════════════════════════════════════════════════

// ── Se0/Se2- ── selenocysteine; selenoenzymes (GPx, TrxR, deiodinases)
// r_eq CSD Se complexes: Se-S ≈2.35 Å, Se-N ≈2.20 Å, Se-O ≈2.05 Å
// Chemistry: chalcogenide, soft-like; S ≥ N > O
inline std::vector<ag4_nbp_override> ag4_se_mode_overrides() {
    return {
        {"SA",  "Se",  2.35,    4.0,    12, 6},  // Se–SA (Se-S selenosulfide; GPx active site)
        {"NA",  "Se",  2.20,    2.0,    12, 6},  // Se–NA (His; TrxR)
        {"N",   "Se",  2.20,    2.0,    12, 6},  // Se–N
        {"OA",  "Se",  2.05,    0.8,    12, 6},  // Se–OA (rare)
        {"HD",  "Se",  1.0,     0.0,    12, 6},  // Se–HD repulsion
    };
}

// ── As3+ ── trigonal pyramidal; arsenic trioxide (Trisenox) AML / trypanosomiasis
// r_eq CSD: As-S(Cys) ≈2.35 Å, As-N ≈2.22 Å, As-O ≈1.95 Å
// Chemistry: class-b soft acid; S >> N > O
inline std::vector<ag4_nbp_override> ag4_as_met_mode_overrides() {
    return {
        {"SA",  "As",  2.35,    5.0,    12, 6},  // As–SA (Cys; dithiol binding site in ATO)
        {"NA",  "As",  2.22,    1.5,    12, 6},  // As–NA (His; secondary)
        {"N",   "As",  2.22,    1.5,    12, 6},  // As–N
        {"OA",  "As",  1.95,    0.8,    12, 6},  // As–OA (arsenite oxide O; very short)
        {"HD",  "As",  1.0,     0.0,    12, 6},  // As–HD repulsion
    };
}

// ── Ge4+ ── tetrahedral; organic germanium (spirogermanium) anticancer
// r_eq CSD Ge(IV): Ge-O ≈1.93 Å, Ge-N ≈2.10 Å, Ge-S ≈2.40 Å; borderline
inline std::vector<ag4_nbp_override> ag4_ge_mode_overrides() {
    return {
        {"OA",  "Ge",  1.93,    2.5,    12, 6},  // Ge–OA (dominant; propane-1,2-diol ligand)
        {"NA",  "Ge",  2.10,    2.0,    12, 6},  // Ge–NA
        {"N",   "Ge",  2.10,    2.0,    12, 6},  // Ge–N
        {"SA",  "Ge",  2.40,    1.0,    12, 6},  // Ge–SA
        {"HD",  "Ge",  1.0,     0.0,    12, 6},  // Ge–HD repulsion
    };
}

// ════════════════════════════════════════════════════════════════════════════
// ── Other main-group elements ────────────────────────────────────────────
// ════════════════════════════════════════════════════════════════════════════

// ── B3+ ── trigonal/tetrahedral; bortezomib (Velcade) proteasome inhibitor
// r_eq CSD boronic acid complexes: B-O ≈1.55 Å, B-N ≈1.62 Å (covalent bond lengths)
// Very strong OA preference — boron acts as electrophilic warhead, B←:O covalent
inline std::vector<ag4_nbp_override> ag4_b_mode_overrides() {
    return {
        {"OA",  "B",   1.55,    6.0,    12, 6},  // B–OA (dominant; Ser/Thr hydroxyl; Velcade mode)
        {"NA",  "B",   1.62,    4.0,    12, 6},  // B–NA (His; imidazole N; secondary)
        {"N",   "B",   1.62,    4.0,    12, 6},  // B–N
        {"SA",  "B",   2.00,    0.5,    12, 6},  // B–SA (Cys; rare)
        {"HD",  "B",   0.8,     0.0,    12, 6},  // B–HD repulsion
    };
}

// ── Be2+ ── tetrahedral, extremely hard, very small cation
// r_eq CSD Be(II): Be-O ≈1.65 Å, Be-N ≈1.83 Å (very short bonds)
inline std::vector<ag4_nbp_override> ag4_be_mode_overrides() {
    return {
        {"OA",  "Be",  1.65,    5.0,    12, 6},  // Be–OA (dominant; very short bond)
        {"NA",  "Be",  1.83,    2.0,    12, 6},  // Be–NA (His)
        {"N",   "Be",  1.83,    2.0,    12, 6},  // Be–N
        {"SA",  "Be",  2.20,    0.2,    12, 6},  // Be–SA (negligible)
        {"HD",  "Be",  0.8,     0.0,    12, 6},  // Be–HD repulsion
    };
}

// ── Te2- ── heavy chalcogenide, soft; organic telluride anticancer
// r_eq CSD: Te-S ≈2.65 Å, Te-N ≈2.40 Å, Te-O ≈2.20 Å; S ≥ N > O (soft)
inline std::vector<ag4_nbp_override> ag4_te_mode_overrides() {
    return {
        {"SA",  "Te",  2.65,    4.5,    12, 6},  // Te–SA (Se-Te / S-Te bond; dominant)
        {"NA",  "Te",  2.40,    2.0,    12, 6},  // Te–NA (His)
        {"N",   "Te",  2.40,    2.0,    12, 6},  // Te–N
        {"OA",  "Te",  2.20,    0.8,    12, 6},  // Te–OA (weak)
        {"HD",  "Te",  1.0,     0.0,    12, 6},  // Te–HD repulsion
    };
}

// ── Po2+/4+ ── heavy chalcogenide/metalloid, very soft; research / env. toxicology
// r_eq estimated: Po-S ≈2.72 Å, Po-N ≈2.55 Å, Po-O ≈2.35 Å; similar to Pb/Bi
inline std::vector<ag4_nbp_override> ag4_po_mode_overrides() {
    return {
        {"SA",  "Po",  2.72,    3.5,    12, 6},  // Po–SA (soft; Te/S-philic)
        {"NA",  "Po",  2.55,    1.5,    12, 6},  // Po–NA
        {"N",   "Po",  2.55,    1.5,    12, 6},  // Po–N
        {"OA",  "Po",  2.35,    1.0,    12, 6},  // Po–OA
        {"HD",  "Po",  1.5,     0.0,    12, 6},  // Po–HD repulsion
    };
}

// ── At- ── radiohalogen; At-211 targeted alpha therapy (like I but softer)
// r_eq estimated CSD: At-C ≈2.35 Å (covalent), At-S ≈2.82 Å (electrostatic binding)
// Soft halide: S > N >> O in non-covalent contacts
inline std::vector<ag4_nbp_override> ag4_at_mode_overrides() {
    return {
        {"SA",  "At",  2.82,    3.0,    12, 6},  // At–SA (electrostatic/halogen bond to S)
        {"NA",  "At",  2.70,    2.0,    12, 6},  // At–NA (halogen bond to N)
        {"N",   "At",  2.70,    2.0,    12, 6},  // At–N
        {"OA",  "At",  2.55,    0.8,    12, 6},  // At–OA (weaker halogen bond)
        {"HD",  "At",  1.5,     0.0,    12, 6},  // At–HD repulsion
    };
}

// ════════════════════════════════════════════════════════════════════════════
// ── Heavy alkali metals ─────────────────────────────────────────────────
// ════════════════════════════════════════════════════════════════════════════

// ── Rb+ ── K+ mimic, hard; Rb-82 cardiac PET imaging
// r_eq CSD: Rb-O ≈2.87 Å, Rb-N ≈2.99 Å (larger than K)
inline std::vector<ag4_nbp_override> ag4_rb_mode_overrides() {
    return {
        {"OA",  "Rb",  2.87,    1.5,    12, 6},  // Rb–OA (dominant; like K but larger)
        {"NA",  "Rb",  2.99,    0.3,    12, 6},  // Rb–NA (very weak)
        {"N",   "Rb",  2.99,    0.3,    12, 6},  // Rb–N
        {"SA",  "Rb",  3.20,    0.1,    12, 6},  // Rb–SA (negligible)
        {"HD",  "Rb",  1.5,     0.0,    12, 6},  // Rb–HD repulsion
    };
}

// ── Cs+ ── largest alkali metal cation, very hard; Cs-131 brachytherapy
// r_eq CSD: Cs-O ≈3.15 Å, Cs-N ≈3.25 Å
inline std::vector<ag4_nbp_override> ag4_cs_mode_overrides() {
    return {
        {"OA",  "Cs",  3.15,    1.0,    12, 6},  // Cs–OA (dominant; very weak due to size)
        {"NA",  "Cs",  3.25,    0.2,    12, 6},  // Cs–NA (essentially none)
        {"N",   "Cs",  3.25,    0.2,    12, 6},  // Cs–N
        {"SA",  "Cs",  3.50,    0.1,    12, 6},  // Cs–SA (negligible)
        {"HD",  "Cs",  1.5,     0.0,    12, 6},  // Cs–HD repulsion
    };
}

// ── Pm3+ ── lanthanide between Nd and Sm; Pm-149 brachytherapy
// r_eq interpolated: Pm-O ≈2.43 Å, Pm-N ≈2.58 Å (between Nd and Sm)
inline std::vector<ag4_nbp_override> ag4_pm_mode_overrides() {
    return {{"OA","Pm",2.43,3.0,12,6},{"NA","Pm",2.58,1.5,12,6},{"N","Pm",2.58,1.5,12,6},
            {"SA","Pm",2.81,0.2,12,6},{"HD","Pm",1.5,0.0,12,6}};
}

// ── Ra2+ ── largest alkaline-earth cation; Ra-223 dichloride (Xofigo) FDA-approved
// r_eq CSD/estimated: Ra-O ≈2.88 Å, Ra-N ≈2.98 Å; hard, Ca/Ba-mimic but larger
inline std::vector<ag4_nbp_override> ag4_ra_mode_overrides() {
    return {
        {"OA",  "Ra",  2.88,    2.0,    12, 6},  // Ra–OA (dominant; phosphate/carboxylate)
        {"NA",  "Ra",  2.98,    0.3,    12, 6},  // Ra–NA (very rare)
        {"N",   "Ra",  2.98,    0.3,    12, 6},  // Ra–N
        {"SA",  "Ra",  3.30,    0.1,    12, 6},  // Ra–SA (negligible)
        {"HD",  "Ra",  1.5,     0.0,    12, 6},  // Ra–HD repulsion
    };
}

// ════════════════════════════════════════════════════════════════════════════
// ── Actinides (all hard O-preferring; similar to heavy lanthanides) ──────
// r_eq from Shannon ionic radii + CSD estimates for early actinides
// All: O >> N >> S; ε(OA)=2.5–3.0, ε(NA)=1.0–1.5, ε(SA)=0.2
// ════════════════════════════════════════════════════════════════════════════
inline std::vector<ag4_nbp_override> ag4_ac_mode_overrides() {
    return {{"OA","Ac",2.78,3.0,12,6},{"NA","Ac",2.92,1.2,12,6},{"N","Ac",2.92,1.2,12,6},
            {"SA","Ac",3.10,0.2,12,6},{"HD","Ac",1.5,0.0,12,6}};
}
inline std::vector<ag4_nbp_override> ag4_th_mode_overrides() {
    return {{"OA","Th",2.52,3.5,12,6},{"NA","Th",2.68,1.5,12,6},{"N","Th",2.68,1.5,12,6},
            {"SA","Th",2.88,0.2,12,6},{"HD","Th",1.5,0.0,12,6}};
}
inline std::vector<ag4_nbp_override> ag4_pa_mode_overrides() {
    return {{"OA","Pa",2.48,3.5,12,6},{"NA","Pa",2.63,1.5,12,6},{"N","Pa",2.63,1.5,12,6},
            {"SA","Pa",2.83,0.2,12,6},{"HD","Pa",1.5,0.0,12,6}};
}
inline std::vector<ag4_nbp_override> ag4_u_mode_overrides() {
    return {{"OA","U", 2.50,3.5,12,6},{"NA","U", 2.65,1.5,12,6},{"N","U", 2.65,1.5,12,6},
            {"SA","U", 2.85,0.2,12,6},{"HD","U", 1.5,0.0,12,6}};
}
inline std::vector<ag4_nbp_override> ag4_np_mode_overrides() {
    return {{"OA","Np",2.48,3.0,12,6},{"NA","Np",2.63,1.2,12,6},{"N","Np",2.63,1.2,12,6},
            {"SA","Np",2.83,0.2,12,6},{"HD","Np",1.5,0.0,12,6}};
}
inline std::vector<ag4_nbp_override> ag4_pu_mode_overrides() {
    return {{"OA","Pu",2.44,3.0,12,6},{"NA","Pu",2.59,1.2,12,6},{"N","Pu",2.59,1.2,12,6},
            {"SA","Pu",2.79,0.2,12,6},{"HD","Pu",1.5,0.0,12,6}};
}
inline std::vector<ag4_nbp_override> ag4_am_mode_overrides() {
    return {{"OA","Am",2.54,3.0,12,6},{"NA","Am",2.69,1.2,12,6},{"N","Am",2.69,1.2,12,6},
            {"SA","Am",2.89,0.2,12,6},{"HD","Am",1.5,0.0,12,6}};
}
inline std::vector<ag4_nbp_override> ag4_cm_mode_overrides() {
    return {{"OA","Cm",2.51,3.0,12,6},{"NA","Cm",2.66,1.2,12,6},{"N","Cm",2.66,1.2,12,6},
            {"SA","Cm",2.86,0.2,12,6},{"HD","Cm",1.5,0.0,12,6}};
}
inline std::vector<ag4_nbp_override> ag4_bk_mode_overrides() {
    return {{"OA","Bk",2.50,3.0,12,6},{"NA","Bk",2.64,1.2,12,6},{"N","Bk",2.64,1.2,12,6},
            {"SA","Bk",2.84,0.2,12,6},{"HD","Bk",1.5,0.0,12,6}};
}
inline std::vector<ag4_nbp_override> ag4_cf_mode_overrides() {
    return {{"OA","Cf",2.48,3.0,12,6},{"NA","Cf",2.62,1.2,12,6},{"N","Cf",2.62,1.2,12,6},
            {"SA","Cf",2.82,0.2,12,6},{"HD","Cf",1.5,0.0,12,6}};
}
inline std::vector<ag4_nbp_override> ag4_es_mode_overrides() {
    // AD_TYPE_E = Einsteinium; using "E" as atom-type string
    return {{"OA","E", 2.46,3.0,12,6},{"NA","E", 2.60,1.2,12,6},{"N","E", 2.60,1.2,12,6},
            {"SA","E", 2.80,0.2,12,6},{"HD","E", 1.5,0.0,12,6}};
}
inline std::vector<ag4_nbp_override> ag4_fm_mode_overrides() {
    return {{"OA","Fm",2.43,3.0,12,6},{"NA","Fm",2.58,1.2,12,6},{"N","Fm",2.58,1.2,12,6},
            {"SA","Fm",2.78,0.2,12,6},{"HD","Fm",1.5,0.0,12,6}};
}

// ── Oxidation-state variants ──────────────────────────────────────────────────
// Fe²⁺ high-spin oct.; Harding 2006: Fe²⁺-O=2.12, Fe²⁺-N=2.20 Å
// HSAB borderline: ε(N)≈ε(O) > ε(S)
inline std::vector<ag4_nbp_override> ag4_fe2_mode_overrides() {
    return {{"OA","Fe",2.12,5.0,12,6},{"NA","Fe",2.20,5.5,12,6},{"N","Fe",2.20,5.5,12,6},
            {"SA","Fe",2.54,1.8,12,6},{"HD","Fe",1.5,0.0,12,6}};
}
// Fe³⁺ high-spin oct.; Harding 2006: Fe³⁺-O=2.00, Fe³⁺-N=2.10 Å
// HSAB harder than Fe²⁺: ε(O) > ε(N) >> ε(S)
inline std::vector<ag4_nbp_override> ag4_fe3_mode_overrides() {
    return {{"OA","Fe",2.00,6.5,12,6},{"NA","Fe",2.10,4.5,12,6},{"N","Fe",2.10,4.5,12,6},
            {"SA","Fe",2.26,0.8,12,6},{"HD","Fe",1.5,0.0,12,6}};
}
// Cu⁺ soft d10; type-1 blue-copper CSD: Cu-S(Cys)=2.30, Cu-N(His)=2.05 Å
// HSAB soft: ε(S) >> ε(N) > ε(O)
inline std::vector<ag4_nbp_override> ag4_cu1_mode_overrides() {
    return {{"SA","Cu",2.30,9.0,12,6},{"NA","Cu",2.05,4.5,12,6},{"N","Cu",2.05,4.5,12,6},
            {"OA","Cu",2.40,1.5,12,6},{"HD","Cu",1.5,0.0,12,6}};
}
// Cu²⁺ Jahn-Teller; CSD: Cu-N(His)_eq=2.02, Cu-O_eq=1.97, Cu-S(Cys)_ax=2.30 Å
// HSAB borderline-soft: ε(N)≈ε(O) > ε(S)
inline std::vector<ag4_nbp_override> ag4_cu2_mode_overrides() {
    return {{"NA","Cu",2.02,5.5,12,6},{"N","Cu",2.02,5.5,12,6},
            {"OA","Cu",1.97,5.0,12,6},{"SA","Cu",2.30,3.0,12,6},
            {"HD","Cu",1.5,0.0,12,6}};
}
// Mn²⁺ high-spin oct.; Harding 2006: Mn²⁺-O=2.18, Mn²⁺-N=2.25 Å
// HSAB hard: ε(O) > ε(N) >> ε(S)
inline std::vector<ag4_nbp_override> ag4_mn2_mode_overrides() {
    return {{"OA","Mn",2.18,5.5,12,6},{"NA","Mn",2.25,3.5,12,6},{"N","Mn",2.25,3.5,12,6},
            {"SA","Mn",2.68,0.5,12,6},{"HD","Mn",1.5,0.0,12,6}};
}
// Mn³⁺ Jahn-Teller; Harding 2006: Mn³⁺-O=2.02, Mn³⁺-N=2.10 Å
// HSAB harder than Mn²⁺: ε(O) >> ε(N) >> ε(S)
inline std::vector<ag4_nbp_override> ag4_mn3_mode_overrides() {
    return {{"OA","Mn",2.02,7.0,12,6},{"NA","Mn",2.10,4.0,12,6},{"N","Mn",2.10,4.0,12,6},
            {"SA","Mn",2.42,0.3,12,6},{"HD","Mn",1.5,0.0,12,6}};
}
// Co²⁺ oct.; Harding 2006: Co²⁺-O=2.10, Co²⁺-N=2.12 Å
// HSAB borderline: ε(N)≈ε(O) > ε(S)
inline std::vector<ag4_nbp_override> ag4_co2_mode_overrides() {
    return {{"OA","Co",2.10,5.0,12,6},{"NA","Co",2.12,5.5,12,6},{"N","Co",2.12,5.5,12,6},
            {"SA","Co",2.52,2.0,12,6},{"HD","Co",1.5,0.0,12,6}};
}
// Co³⁺ inert d6; CSD: Co³⁺-N=1.98 (cobalamin, NAMI-A), Co³⁺-O=1.90 Å
// HSAB N-selective (π back-donation): ε(N) > ε(O) >> ε(S)
inline std::vector<ag4_nbp_override> ag4_co3_mode_overrides() {
    return {{"NA","Co",1.98,7.0,12,6},{"N","Co",1.98,7.0,12,6},
            {"OA","Co",1.90,4.5,12,6},{"SA","Co",2.25,0.8,12,6},
            {"HD","Co",1.5,0.0,12,6}};
}
// As³⁺ pyramidal; CSD: As³⁺-S(Cys)=2.25 (PML triad, Trisenox AML)
// HSAB borderline-soft; thiolate dominant: ε(S) >> ε(N) > ε(O)
inline std::vector<ag4_nbp_override> ag4_as3_mode_overrides() {
    return {{"SA","As",2.25,8.0,12,6},{"NA","As",2.35,3.0,12,6},{"N","As",2.35,3.0,12,6},
            {"OA","As",1.96,1.5,12,6},{"HD","As",1.5,0.0,12,6}};
}
// As⁵⁺ arsenate tetrahedral; r(As⁵⁺-O)≈1.68 Å (phosphate mimic, PO₄ analogue)
// HSAB ultra-hard (pentavalent d0): ε(O) >> ε(N) >> ε(S)
inline std::vector<ag4_nbp_override> ag4_as5_mode_overrides() {
    return {{"OA","As",1.68,7.0,12,6},{"NA","As",2.05,2.0,12,6},{"N","As",2.05,2.0,12,6},
            {"SA","As",2.45,0.3,12,6},{"HD","As",1.5,0.0,12,6}};
}
// Sb³⁺ antimonite; r(Sb³⁺-S)≈2.45 (Pentostam Cys-targeting, leishmania)
// HSAB soft-ish: ε(S) > ε(O) > ε(N)
inline std::vector<ag4_nbp_override> ag4_sb3_mode_overrides() {
    return {{"SA","Sb",2.45,5.5,12,6},{"OA","Sb",2.15,3.5,12,6},
            {"NA","Sb",2.35,2.0,12,6},{"N","Sb",2.35,2.0,12,6},
            {"HD","Sb",1.5,0.0,12,6}};
}
// Sb⁵⁺ meglumine antimoniate; r(Sb⁵⁺-O)≈1.97 Å octahedral
// HSAB hard (pentavalent): ε(O) >> ε(N) > ε(S)
inline std::vector<ag4_nbp_override> ag4_sb5_mode_overrides() {
    return {{"OA","Sb",1.97,6.0,12,6},{"NA","Sb",2.15,3.0,12,6},{"N","Sb",2.15,3.0,12,6},
            {"SA","Sb",2.50,0.5,12,6},{"HD","Sb",1.5,0.0,12,6}};
}
// VO²⁺ vanadyl equatorial; CSD: V-O_eq=2.05, V-N=2.15 Å
// (axial V=O=1.60 Å is spectator; equatorial belt is ligand-binding site)
// HSAB borderline-hard: ε(O) > ε(N) >> ε(S)
inline std::vector<ag4_nbp_override> ag4_v4_mode_overrides() {
    return {{"OA","V",2.05,6.0,12,6},{"NA","V",2.15,4.0,12,6},{"N","V",2.15,4.0,12,6},
            {"SA","V",2.65,0.8,12,6},{"HD","V",1.5,0.0,12,6}};
}
// VO₄³⁻ vanadate phosphatase inhibitor; r(V⁵⁺-O)≈1.73 Å tetrahedral
// HSAB ultra-hard (pentavalent d0): ε(O) >> ε(N) >> ε(S)
inline std::vector<ag4_nbp_override> ag4_v5_mode_overrides() {
    return {{"OA","V",1.73,8.0,12,6},{"NA","V",2.05,2.5,12,6},{"N","V",2.05,2.5,12,6},
            {"SA","V",2.50,0.3,12,6},{"HD","V",1.5,0.0,12,6}};
}
// Ni²⁺ oxidation-state specific; Brese & O'Keeffe 1991: Ni-O r0=1.654, Ni-N r0=1.679, Ni-S r0=1.978
// sq-planar d8: N > S > O; HSAB borderline-soft
inline std::vector<ag4_nbp_override> ag4_ni2_mode_overrides() {
    return {
        {"NA","Ni", 2.06, 5.0, 12, 6},  // Ni²⁺–NA (His; Ni-SOD equatorial; r_eq Harding: 2.06)
        {"N", "Ni", 2.06, 5.0, 12, 6},  // Ni²⁺–N
        {"SA","Ni", 2.15, 4.5, 12, 6},  // Ni²⁺–SA (Cys; NiSOD axial; Harding: 2.18)
        {"OA","Ni", 2.06, 3.0, 12, 6},  // Ni²⁺–OA (Asp/Glu; urease)
        {"HD","Ni", 1.0,  0.0, 12, 6},  // repulsion
    };
}
// Ni³⁺; Brese & O'Keeffe: Ni³⁺-O r0≈1.62, Ni³⁺-N r0≈1.65, Ni³⁺-S r0≈1.95
// Ni-Fe hydrogenase activated state: N>>S>O; shorter bonds than Ni²⁺
inline std::vector<ag4_nbp_override> ag4_ni3_mode_overrides() {
    return {
        {"NA","Ni", 2.00, 6.0, 12, 6},  // Ni³⁺–NA (His; shorter by ~0.06 vs Ni²⁺)
        {"N", "Ni", 2.00, 6.0, 12, 6},  // Ni³⁺–N
        {"SA","Ni", 2.10, 5.5, 12, 6},  // Ni³⁺–SA (Cys; stronger than Ni²⁺)
        {"OA","Ni", 2.00, 3.5, 12, 6},  // Ni³⁺–OA
        {"HD","Ni", 1.0,  0.0, 12, 6},  // repulsion
    };
}
// TODO Mo-Fe heteronuclear (nitrogenase FeMo-co):
//   Correct BVS for Fe in FeMo-co requires knowing Mo oxidation state simultaneously.
//   Placeholder for future multi-metal BVS coupling implementation.
//   Polynuclear Mn (OEC Mn4CaO5): requires 4-centre BVS sum with Ca bridging.
//   These are tracked as future work; single-metal modes mo4/mo6 are used in the interim.

// Mo⁴⁺ sulfido; CSD: Mo-S(xanthine oxidase)=2.35, Mo-O_eq=2.10 Å
// HSAB borderline (Mo⁴⁺): ε(S) > ε(O) > ε(N)
inline std::vector<ag4_nbp_override> ag4_mo4_mode_overrides() {
    return {{"SA","Mo",2.35,6.0,12,6},{"OA","Mo",2.10,4.5,12,6},
            {"NA","Mo",2.25,2.5,12,6},{"N","Mo",2.25,2.5,12,6},
            {"HD","Mo",1.5,0.0,12,6}};
}
// Mo⁶⁺ dioxo; r(Mo⁶⁺=O)≈1.73 Å (sulfite oxidase, molybdate MoO₄²⁻)
// HSAB hard (hexavalent d0): ε(O) >> ε(N) > ε(S)
inline std::vector<ag4_nbp_override> ag4_mo6_mode_overrides() {
    return {{"OA","Mo",1.73,8.0,12,6},{"NA","Mo",2.00,2.5,12,6},{"N","Mo",2.00,2.5,12,6},
            {"SA","Mo",2.45,0.5,12,6},{"HD","Mo",1.5,0.0,12,6}};
}
// UO₂²⁺ uranyl equatorial; CSD: U-O_eq=2.44, U-N_eq=2.61 Å
// (axial U=O=1.76 Å are spectators; equatorial belt is drug-binding site)
// HSAB hard actinide: ε(O) >> ε(N) > ε(S)
inline std::vector<ag4_nbp_override> ag4_uo2_mode_overrides() {
    return {{"OA","U",2.44,4.5,12,6},{"NA","U",2.61,1.8,12,6},{"N","U",2.61,1.8,12,6},
            {"SA","U",2.80,0.3,12,6},{"HD","U",1.5,0.0,12,6}};
}
// Gd-DTPA (Magnevist): 9-coord 3N+5O from DTPA + 1 inner-sphere H₂O
// Gd-O_DTPA=2.43, Gd-N_DTPA=2.54 Å (Burai 2009 CSD survey)
// Chelation reduces Lewis acidity; N selectivity enhanced vs free Gd³⁺
inline std::vector<ag4_nbp_override> ag4_gd_dtpa_mode_overrides() {
    return {{"OA","Gd",2.43,2.8,12,6},{"NA","Gd",2.54,2.0,12,6},{"N","Gd",2.54,2.0,12,6},
            {"SA","Gd",3.10,0.1,12,6},{"HD","Gd",1.5,0.0,12,6}};
}
// Gd-DOTA (Dotarem): 9-coord 4N+4O macrocycle + 1 H₂O
// Gd-O_DOTA=2.43, Gd-N_DOTA=2.56 Å; macrocyclic rigidity pre-organizes N
// Most N-balanced Gd mode: ε(N)≈ε(O)
inline std::vector<ag4_nbp_override> ag4_gd_dota_mode_overrides() {
    return {{"OA","Gd",2.43,2.5,12,6},{"NA","Gd",2.56,2.5,12,6},{"N","Gd",2.56,2.5,12,6},
            {"SA","Gd",3.10,0.1,12,6},{"HD","Gd",1.5,0.0,12,6}};
}

// ── O4: Aqua-mode variants (_aq) ─────────────────────────────────────────────
// Models metal with ≥1 coordinated water in the receptor.
// r_eq unchanged; ε reduced ~20% (water competes for open sites).
inline std::vector<ag4_nbp_override> ag4_mg_aq_mode_overrides() {
    return {{"OA","Mg",2.07,2.8,12,6},{"NA","Mg",2.20,1.2,12,6},{"N","Mg",2.20,1.2,12,6},
            {"SA","Mg",2.50,0.4,12,6},{"HD","Mg",1.0,0.0,12,6}};
}
inline std::vector<ag4_nbp_override> ag4_ca_aq_mode_overrides() {
    return {{"OA","Ca",2.38,2.5,12,6},{"NA","Ca",2.50,0.8,12,6},{"N","Ca",2.50,0.8,12,6},
            {"SA","Ca",2.80,0.2,12,6},{"HD","Ca",1.2,0.0,12,6}};
}
inline std::vector<ag4_nbp_override> ag4_fe3_aq_mode_overrides() {
    return {{"OA","Fe",2.00,5.2,12,6},{"NA","Fe",2.10,3.6,12,6},{"N","Fe",2.10,3.6,12,6},
            {"SA","Fe",2.26,0.6,12,6},{"HD","Fe",1.5,0.0,12,6}};
}
inline std::vector<ag4_nbp_override> ag4_mn2_aq_mode_overrides() {
    return {{"OA","Mn",2.18,4.4,12,6},{"NA","Mn",2.25,2.8,12,6},{"N","Mn",2.25,2.8,12,6},
            {"SA","Mn",2.68,0.4,12,6},{"HD","Mn",1.5,0.0,12,6}};
}
inline std::vector<ag4_nbp_override> ag4_co2_aq_mode_overrides() {
    return {{"OA","Co",2.10,4.0,12,6},{"NA","Co",2.12,4.4,12,6},{"N","Co",2.12,4.4,12,6},
            {"SA","Co",2.52,1.6,12,6},{"HD","Co",1.5,0.0,12,6}};
}

// ── Central dispatch: apply one metal mode's overrides to accumulated list ──
// Replaces the monolithic switch in embedded_ad4_grid.cpp; enables multi-metal.
inline void ag4_apply_metal_mode(ag4_metal_mode mode, std::vector<ag4_nbp_override>& ov) {
    auto append = [&](std::vector<ag4_nbp_override> v) {
        for (const auto& x : v)
            ag4_append_nbp_pair_if_missing(ov, x);
    };
    switch (mode) {
        case ag4_metal_mode::zn:      append(ag4_zn_mode_overrides());      break;
        case ag4_metal_mode::mg:      append(ag4_mg_mode_overrides());      break;
        case ag4_metal_mode::ca:      append(ag4_ca_mode_overrides());      break;
        case ag4_metal_mode::mn:      append(ag4_mn_mode_overrides());      break;
        case ag4_metal_mode::fe:      append(ag4_fe_mode_overrides());      break;
        case ag4_metal_mode::co:      append(ag4_co_mode_overrides());      break;
        case ag4_metal_mode::ni:      append(ag4_ni_mode_overrides());      break;
        case ag4_metal_mode::cu:      append(ag4_cu_mode_overrides());      break;
        case ag4_metal_mode::pt:      append(ag4_pt_mode_overrides());      break;
        case ag4_metal_mode::pd:      append(ag4_pd_mode_overrides());      break;
        case ag4_metal_mode::ru:      append(ag4_ru_mode_overrides());      break;
        case ag4_metal_mode::ir:      append(ag4_ir_mode_overrides());      break;
        case ag4_metal_mode::au:      append(ag4_au_mode_overrides());      break;
        case ag4_metal_mode::cd:      append(ag4_cd_mode_overrides());      break;
        case ag4_metal_mode::hg:      append(ag4_hg_mode_overrides());      break;
        case ag4_metal_mode::na_ion:  append(ag4_na_ion_mode_overrides());  break;
        case ag4_metal_mode::k_ion:   append(ag4_k_ion_mode_overrides());   break;
        case ag4_metal_mode::v:       append(ag4_v_mode_overrides());       break;
        case ag4_metal_mode::cr:      append(ag4_cr_mode_overrides());      break;
        case ag4_metal_mode::ti:      append(ag4_ti_mode_overrides());      break;
        case ag4_metal_mode::sc:      append(ag4_sc_mode_overrides());      break;
        case ag4_metal_mode::y_ion:   append(ag4_y_ion_mode_overrides());   break;
        case ag4_metal_mode::zr:      append(ag4_zr_mode_overrides());      break;
        case ag4_metal_mode::nb:      append(ag4_nb_mode_overrides());      break;
        case ag4_metal_mode::hf:      append(ag4_hf_mode_overrides());      break;
        case ag4_metal_mode::ta:      append(ag4_ta_mode_overrides());      break;
        case ag4_metal_mode::tg_mode: append(ag4_tg_mode_overrides());      break;
        case ag4_metal_mode::mo:      append(ag4_mo_mode_overrides());      break;
        case ag4_metal_mode::rh:      append(ag4_rh_mode_overrides());      break;
        case ag4_metal_mode::ag:      append(ag4_ag_mode_overrides());      break;
        case ag4_metal_mode::tc:      append(ag4_tc_mode_overrides());      break;
        case ag4_metal_mode::re:      append(ag4_re_mode_overrides());      break;
        case ag4_metal_mode::os:      append(ag4_os_mode_overrides());      break;
        case ag4_metal_mode::ga:      append(ag4_ga_mode_overrides());      break;
        case ag4_metal_mode::in_ion:  append(ag4_in_ion_mode_overrides());  break;
        case ag4_metal_mode::sn:      append(ag4_sn_mode_overrides());      break;
        case ag4_metal_mode::sb:      append(ag4_sb_mode_overrides());      break;
        case ag4_metal_mode::bi:      append(ag4_bi_mode_overrides());      break;
        case ag4_metal_mode::tl:      append(ag4_tl_mode_overrides());      break;
        case ag4_metal_mode::pb:      append(ag4_pb_mode_overrides());      break;
        case ag4_metal_mode::li:      append(ag4_li_mode_overrides());      break;
        case ag4_metal_mode::al:      append(ag4_al_mode_overrides());      break;
        case ag4_metal_mode::sr:      append(ag4_sr_mode_overrides());      break;
        case ag4_metal_mode::ba:      append(ag4_ba_mode_overrides());      break;
        case ag4_metal_mode::la:      append(ag4_la_mode_overrides());      break;
        case ag4_metal_mode::ce:      append(ag4_ce_mode_overrides());      break;
        case ag4_metal_mode::pr:      append(ag4_pr_mode_overrides());      break;
        case ag4_metal_mode::nd:      append(ag4_nd_mode_overrides());      break;
        case ag4_metal_mode::sm:      append(ag4_sm_mode_overrides());      break;
        case ag4_metal_mode::eu:      append(ag4_eu_mode_overrides());      break;
        case ag4_metal_mode::gd:      append(ag4_gd_mode_overrides());      break;
        case ag4_metal_mode::tb:      append(ag4_tb_mode_overrides());      break;
        case ag4_metal_mode::dy:      append(ag4_dy_mode_overrides());      break;
        case ag4_metal_mode::ho:      append(ag4_ho_mode_overrides());      break;
        case ag4_metal_mode::er:      append(ag4_er_mode_overrides());      break;
        case ag4_metal_mode::tm:      append(ag4_tm_mode_overrides());      break;
        case ag4_metal_mode::yb:      append(ag4_yb_mode_overrides());      break;
        case ag4_metal_mode::lu:      append(ag4_lu_mode_overrides());      break;
        case ag4_metal_mode::se:      append(ag4_se_mode_overrides());      break;
        case ag4_metal_mode::as_met:  append(ag4_as_met_mode_overrides());  break;
        case ag4_metal_mode::ge:      append(ag4_ge_mode_overrides());      break;
        case ag4_metal_mode::b:       append(ag4_b_mode_overrides());       break;
        case ag4_metal_mode::be:      append(ag4_be_mode_overrides());      break;
        case ag4_metal_mode::te:      append(ag4_te_mode_overrides());      break;
        case ag4_metal_mode::po:      append(ag4_po_mode_overrides());      break;
        case ag4_metal_mode::at:      append(ag4_at_mode_overrides());      break;
        case ag4_metal_mode::rb:      append(ag4_rb_mode_overrides());      break;
        case ag4_metal_mode::cs:      append(ag4_cs_mode_overrides());      break;
        case ag4_metal_mode::pm:      append(ag4_pm_mode_overrides());      break;
        case ag4_metal_mode::ra:      append(ag4_ra_mode_overrides());      break;
        case ag4_metal_mode::ac:      append(ag4_ac_mode_overrides());      break;
        case ag4_metal_mode::th:      append(ag4_th_mode_overrides());      break;
        case ag4_metal_mode::pa:      append(ag4_pa_mode_overrides());      break;
        case ag4_metal_mode::u:       append(ag4_u_mode_overrides());       break;
        case ag4_metal_mode::np:      append(ag4_np_mode_overrides());      break;
        case ag4_metal_mode::pu:      append(ag4_pu_mode_overrides());      break;
        case ag4_metal_mode::am:      append(ag4_am_mode_overrides());      break;
        case ag4_metal_mode::cm:      append(ag4_cm_mode_overrides());      break;
        case ag4_metal_mode::bk:      append(ag4_bk_mode_overrides());      break;
        case ag4_metal_mode::cf:      append(ag4_cf_mode_overrides());      break;
        case ag4_metal_mode::es:      append(ag4_es_mode_overrides());      break;
        case ag4_metal_mode::fm:      append(ag4_fm_mode_overrides());      break;
        // Oxidation-state variants
        case ag4_metal_mode::fe2:     append(ag4_fe2_mode_overrides());     break;
        case ag4_metal_mode::fe3:     append(ag4_fe3_mode_overrides());     break;
        case ag4_metal_mode::cu1:     append(ag4_cu1_mode_overrides());     break;
        case ag4_metal_mode::cu2:     append(ag4_cu2_mode_overrides());     break;
        case ag4_metal_mode::cu2_jt:  append(ag4_cu2_mode_overrides());     break;
        case ag4_metal_mode::mn2:     append(ag4_mn2_mode_overrides());     break;
        case ag4_metal_mode::mn3:     append(ag4_mn3_mode_overrides());     break;
        case ag4_metal_mode::mn3_jt:  append(ag4_mn3_mode_overrides());     break;
        case ag4_metal_mode::co2:     append(ag4_co2_mode_overrides());     break;
        case ag4_metal_mode::co3:     append(ag4_co3_mode_overrides());     break;
        case ag4_metal_mode::as3:     append(ag4_as3_mode_overrides());     break;
        case ag4_metal_mode::as5:     append(ag4_as5_mode_overrides());     break;
        case ag4_metal_mode::sb3:     append(ag4_sb3_mode_overrides());     break;
        case ag4_metal_mode::sb5:     append(ag4_sb5_mode_overrides());     break;
        case ag4_metal_mode::v4:      append(ag4_v4_mode_overrides());      break;
        case ag4_metal_mode::v5:      append(ag4_v5_mode_overrides());      break;
        case ag4_metal_mode::mo4:     append(ag4_mo4_mode_overrides());     break;
        case ag4_metal_mode::mo6:     append(ag4_mo6_mode_overrides());     break;
        case ag4_metal_mode::uo2:     append(ag4_uo2_mode_overrides());     break;
        case ag4_metal_mode::gd_dtpa: append(ag4_gd_dtpa_mode_overrides()); break;
        case ag4_metal_mode::gd_dota: append(ag4_gd_dota_mode_overrides()); break;
        // Aqua-mode variants
        case ag4_metal_mode::mg_aq:   append(ag4_mg_aq_mode_overrides());   break;
        case ag4_metal_mode::ca_aq:   append(ag4_ca_aq_mode_overrides());   break;
        case ag4_metal_mode::fe3_aq:  append(ag4_fe3_aq_mode_overrides());  break;
        case ag4_metal_mode::mn2_aq:  append(ag4_mn2_aq_mode_overrides());  break;
        case ag4_metal_mode::co2_aq:  append(ag4_co2_aq_mode_overrides());  break;
        case ag4_metal_mode::ni2:     append(ag4_ni2_mode_overrides());     break;
        case ag4_metal_mode::ni3:     append(ag4_ni3_mode_overrides());     break;
        default: break;
    }
}

// ── Square-planar / linear SQ pseudoatom nbp_r_eps overrides (M4) ──────────
// Probe → SQ directional attractor: r_eq=0.25 Å (probe sits ON the SQ site).
// ε values encode metal-specific ligand-field preference (HSAB + MetalDock 2024).
// Used for: Pt/Pd/Ni/Cu (square planar d8/d9) and Au (linear d10).
inline std::vector<ag4_nbp_override> ag4_sq_nbp_overrides(ag4_metal_mode mode) {
    switch (mode) {
        case ag4_metal_mode::pt:  // sq-planar d8; N >> S > O (DNA-N7, His; cisplatin)
            return {{"NA","SQ",0.25,20.0,12,6}, {"N","SQ",0.25,20.0,12,6},
                    {"SA","SQ",0.25,18.0,12,6}, {"OA","SQ",0.25, 8.0,12,6},
                    {"HD","SQ", 1.0,  0.0,12,6}};
        case ag4_metal_mode::pd:  // sq-planar d8, softer than Pt; S ≥ N > O
            return {{"NA","SQ",0.25,18.0,12,6}, {"N","SQ",0.25,18.0,12,6},
                    {"SA","SQ",0.25,20.0,12,6}, {"OA","SQ",0.25, 6.0,12,6},
                    {"HD","SQ", 1.0,  0.0,12,6}};
        case ag4_metal_mode::ni:  // sq-planar/oct d8; N > S > O (urease His, NiSOD Cys)
        case ag4_metal_mode::ni2:  // same as ni; explicit oxidation-state variant
            return {{"NA","SQ",0.25,16.0,12,6}, {"N","SQ",0.25,16.0,12,6},
                    {"SA","SQ",0.25,14.0,12,6}, {"OA","SQ",0.25, 6.0,12,6},
                    {"HD","SQ", 1.0,  0.0,12,6}};
        case ag4_metal_mode::cu:  // Jahn-Teller d9; S ≥ N > O (cupredoxin: Cys, 2×His, Met)
        case ag4_metal_mode::cu2_jt:
            return {{"NA","SQ",0.25,15.0,12,6}, {"N","SQ",0.25,15.0,12,6},
                    {"SA","SQ",0.25,18.0,12,6}, {"OA","SQ",0.25, 5.0,12,6},
                    {"HD","SQ", 1.0,  0.0,12,6}};
        case ag4_metal_mode::cu2:  // sq-planar d9 non-JT; N ≈ O > S
            return {{"NA","SQ",0.25,15.0,12,6}, {"N","SQ",0.25,15.0,12,6},
                    {"OA","SQ",0.25,13.0,12,6}, {"SA","SQ",0.25, 6.0,12,6},
                    {"HD","SQ", 1.0,  0.0,12,6}};
        case ag4_metal_mode::au:  // linear d10; extreme S-philic (auranofin → Cys/thioredoxin)
            return {{"SA","SQ",0.25,22.0,12,6},
                    {"NA","SQ",0.25,10.0,12,6}, {"N","SQ",0.25,10.0,12,6},
                    {"OA","SQ",0.25, 3.0,12,6},
                    {"HD","SQ", 1.0,  0.0,12,6}};
        case ag4_metal_mode::hg:  // linear d10; extreme S-philic (metallothionein Cys, softer than Au)
            return {{"SA","SQ",0.25,24.0,12,6},
                    {"NA","SQ",0.25, 8.0,12,6}, {"N","SQ",0.25, 8.0,12,6},
                    {"OA","SQ",0.25, 2.0,12,6},
                    {"HD","SQ", 1.0,  0.0,12,6}};
        case ag4_metal_mode::ag:  // linear d10 Ag⁺; S >> N > O (antimicrobial, Ag-thiolate)
            return {{"SA","SQ",0.25,20.0,12,6},
                    {"NA","SQ",0.25, 9.0,12,6}, {"N","SQ",0.25, 9.0,12,6},
                    {"OA","SQ",0.25, 2.5,12,6},
                    {"HD","SQ", 1.0,  0.0,12,6}};
        default: return {};
    }
}

// Returns true if the given metal_mode needs SQ pseudoatom injection
inline bool ag4_needs_sq_injection(ag4_metal_mode mode) {
    return mode == ag4_metal_mode::pt || mode == ag4_metal_mode::pd
        || mode == ag4_metal_mode::ni  || mode == ag4_metal_mode::ni2
        || mode == ag4_metal_mode::cu  || mode == ag4_metal_mode::cu2
        || mode == ag4_metal_mode::cu2_jt
        || mode == ag4_metal_mode::au
        || mode == ag4_metal_mode::hg   // linear d10; 2-coord like Au
        || mode == ag4_metal_mode::ag;  // linear d10 Ag⁺; 2-coord
}

// ── O1: MH (octahedral) pseudoatom nbp_r_eps overrides ──────────────────────
// Probe → MH directional attractor: r_eq=0.25 Å (probe sits ON the MH site).
// ε values encode metal HSAB preference: hard acid → OA dominant, soft → SA.
// Used for: octahedral metals (Fe2+/3+, Mg, Mn, Co, Ru, Ir, Rh, etc.)
inline std::vector<ag4_nbp_override> ag4_mh_nbp_overrides(ag4_metal_mode mode) {
    switch (mode) {
        case ag4_metal_mode::fe3:  case ag4_metal_mode::fe:
            // Fe3+ hard: OA >> NA > SA
            return {{"OA","MH",0.25,12.0,12,6},{"NA","MH",0.25,8.0,12,6},
                    {"N","MH",0.25,8.0,12,6},  {"SA","MH",0.25,2.0,12,6},
                    {"HD","MH",1.0,0.0,12,6}};
        case ag4_metal_mode::fe2:
            // Fe2+ borderline: OA ≈ NA > SA
            return {{"OA","MH",0.25,10.0,12,6},{"NA","MH",0.25,10.0,12,6},
                    {"N","MH",0.25,10.0,12,6}, {"SA","MH",0.25,4.0,12,6},
                    {"HD","MH",1.0,0.0,12,6}};
        case ag4_metal_mode::mg:  case ag4_metal_mode::na_ion:
            // Mg2+/Na+ very hard: OA dominant
            return {{"OA","MH",0.25,14.0,12,6},{"NA","MH",0.25,5.0,12,6},
                    {"N","MH",0.25,5.0,12,6},  {"SA","MH",0.25,0.5,12,6},
                    {"HD","MH",1.0,0.0,12,6}};
        case ag4_metal_mode::ca:  case ag4_metal_mode::k_ion:
            // Ca2+/K+ ultra-hard O-only
            return {{"OA","MH",0.25,10.0,12,6},{"NA","MH",0.25,2.0,12,6},
                    {"N","MH",0.25,2.0,12,6},  {"SA","MH",0.25,0.2,12,6},
                    {"HD","MH",1.0,0.0,12,6}};
        case ag4_metal_mode::mn:  case ag4_metal_mode::mn2:
            // Mn2+ hard: OA > NA > SA
            return {{"OA","MH",0.25,11.0,12,6},{"NA","MH",0.25,7.0,12,6},
                    {"N","MH",0.25,7.0,12,6},  {"SA","MH",0.25,1.0,12,6},
                    {"HD","MH",1.0,0.0,12,6}};
        case ag4_metal_mode::mn3:
            // Mn3+ harder still: OA >> NA >> SA
            return {{"OA","MH",0.25,14.0,12,6},{"NA","MH",0.25,8.0,12,6},
                    {"N","MH",0.25,8.0,12,6},  {"SA","MH",0.25,0.5,12,6},
                    {"HD","MH",1.0,0.0,12,6}};
        case ag4_metal_mode::co:  case ag4_metal_mode::co2:
            // Co2+ borderline: NA ≈ OA > SA
            return {{"NA","MH",0.25,10.0,12,6},{"N","MH",0.25,10.0,12,6},
                    {"OA","MH",0.25,9.0,12,6}, {"SA","MH",0.25,3.0,12,6},
                    {"HD","MH",1.0,0.0,12,6}};
        case ag4_metal_mode::co3:
            // Co3+ inert d6: NA >> OA >> SA
            return {{"NA","MH",0.25,14.0,12,6},{"N","MH",0.25,14.0,12,6},
                    {"OA","MH",0.25,9.0,12,6}, {"SA","MH",0.25,1.5,12,6},
                    {"HD","MH",1.0,0.0,12,6}};
        case ag4_metal_mode::ru:
            // Ru3+ octahedral, N≈S > O (MetalDock-calibrated)
            return {{"NA","MH",0.25,11.0,12,6},{"N","MH",0.25,11.0,12,6},
                    {"OA","MH",0.25,8.0,12,6}, {"SA","MH",0.25,8.0,12,6},
                    {"HD","MH",1.0,0.0,12,6}};
        case ag4_metal_mode::ir:
            // Ir3+ strong N-preference
            return {{"NA","MH",0.25,13.0,12,6},{"N","MH",0.25,13.0,12,6},
                    {"OA","MH",0.25,8.5,12,6}, {"SA","MH",0.25,10.0,12,6},
                    {"HD","MH",1.0,0.0,12,6}};
        case ag4_metal_mode::rh:
            return {{"NA","MH",0.25,12.0,12,6},{"N","MH",0.25,12.0,12,6},
                    {"OA","MH",0.25,8.0,12,6}, {"SA","MH",0.25,7.0,12,6},
                    {"HD","MH",1.0,0.0,12,6}};
        case ag4_metal_mode::ni3:
            // Ni³⁺ activated (Ni-Fe hydrogenase): N >> S > O; shorter bonds than Ni²⁺
            return {{"NA","MH",0.25,13.0,12,6},{"N","MH",0.25,13.0,12,6},
                    {"SA","MH",0.25,10.0,12,6},{"OA","MH",0.25, 3.5,12,6},
                    {"HD","MH", 1.0, 0.0,12,6}};
        default: return {};
    }
}

// Returns true when this metal mode should have MH octahedral pseudoatom injection
inline bool ag4_needs_mh_injection(ag4_metal_mode mode) {
    switch (mode) {
        case ag4_metal_mode::fe:    case ag4_metal_mode::fe2:  case ag4_metal_mode::fe3:
        case ag4_metal_mode::mg:    case ag4_metal_mode::ca:   case ag4_metal_mode::mn:
        case ag4_metal_mode::mn2:   case ag4_metal_mode::mn3:  case ag4_metal_mode::co:
        case ag4_metal_mode::co2:   case ag4_metal_mode::co3:  case ag4_metal_mode::ru:
        case ag4_metal_mode::ir:    case ag4_metal_mode::rh:   case ag4_metal_mode::na_ion:
        case ag4_metal_mode::k_ion: case ag4_metal_mode::ni3:
            return true;
        default: return false;
    }
}

// Metal names, coordination distance (Å), and max coordination number for MH injection
struct ag4_mh_injection_params_t {
    std::vector<std::string> metal_names;
    double coord_dist;
    int    max_coord;
};
inline ag4_mh_injection_params_t ag4_mh_injection_params(ag4_metal_mode mode) {
    switch (mode) {
        case ag4_metal_mode::fe:  case ag4_metal_mode::fe2:
            return {{"Fe"}, 2.12, 6};
        case ag4_metal_mode::fe3:
            return {{"Fe"}, 2.00, 6};
        case ag4_metal_mode::mg:
            return {{"Mg"}, 2.10, 6};
        case ag4_metal_mode::ca:
            return {{"Ca"}, 2.40, 7};  // Ca often 7-coord
        case ag4_metal_mode::mn:  case ag4_metal_mode::mn2:
            return {{"Mn"}, 2.18, 6};
        case ag4_metal_mode::mn3:
            return {{"Mn"}, 2.02, 6};
        case ag4_metal_mode::ni3:
            return {{"Ni"}, 2.00, 6};  // octahedral activated Ni-Fe hydrogenase
        case ag4_metal_mode::co:  case ag4_metal_mode::co2:
            return {{"Co"}, 2.10, 6};
        case ag4_metal_mode::co3:
            return {{"Co"}, 1.95, 6};
        case ag4_metal_mode::ru:
            return {{"Ru"}, 2.09, 6};
        case ag4_metal_mode::ir:
            return {{"Ir"}, 2.05, 6};
        case ag4_metal_mode::rh:
            return {{"Rh"}, 2.05, 6};
        case ag4_metal_mode::na_ion:
            return {{"Na"}, 2.35, 6};
        case ag4_metal_mode::k_ion:
            return {{"K"},  2.70, 7};
        default: return {{}, 2.10, 6};
    }
}

struct ag4_sq_injection_params_t {
    std::vector<std::string> metal_names;
    double coord_dist;
    int    max_coord;
};
inline ag4_sq_injection_params_t ag4_sq_injection_params(ag4_metal_mode mode) {
    switch (mode) {
        case ag4_metal_mode::pt: return {{"Pt"}, 2.05, 4};
        case ag4_metal_mode::pd: return {{"Pd"}, 2.00, 4};
        case ag4_metal_mode::ni:  case ag4_metal_mode::ni2:
            return {{"Ni"}, 2.06, 4};  // sq-planar, r(Ni-N/O)≈2.06 Å
        case ag4_metal_mode::cu: return {{"Cu"}, 2.05, 4};
        case ag4_metal_mode::cu2: return {{"Cu"}, 2.02, 4};  // sq-planar Cu²⁺, r(Cu-N/O)≈2.02 Å
        case ag4_metal_mode::cu2_jt: return {{"Cu"}, 2.03, 4};
        case ag4_metal_mode::au: return {{"Au"}, 2.30, 2};
        case ag4_metal_mode::hg: return {{"Hg"}, 2.35, 2};  // linear Hg²⁺, r(Hg-S)≈2.35 Å
        case ag4_metal_mode::ag: return {{"Ag"}, 2.30, 2};  // linear Ag⁺, r(Ag-N/S)≈2.30 Å
        default: return {{}, 2.05, 4};
    }
}

// ── TZ (tetrahedral) pseudoatom overrides for Cd — Zn-mimic metals (M4-TZ) ──
// Cd²⁺ is a d10 tetrahedral metal that replaces Zn in metalloproteins;
// it needs TZ injection (same geometry as Zn) but with Cd-tuned ε values.
// Parameters: HSAB soft acid; S >> N > O (metallothionein Cys, carbonic anhydrase).
inline std::vector<ag4_nbp_override> ag4_tz_nbp_overrides(ag4_metal_mode mode) {
    switch (mode) {
        case ag4_metal_mode::cd:  // Cd²⁺ d10 tetrahedral; S >> N > O (softer than Zn)
            return {{"SA","TZ",0.25,20.0,12,6},
                    {"NA","TZ",0.25,10.0,12,6}, {"N","TZ",0.25,10.0,12,6},
                    {"OA","TZ",0.25, 4.0,12,6},
                    {"HD","TZ", 1.0,  0.0,12,6}};
        default: return {};
    }
}

// Returns true when this metal mode needs TZ tetrahedral pseudoatom injection
inline bool ag4_needs_tz_injection(ag4_metal_mode mode) {
    return mode == ag4_metal_mode::cd;
}

// Metal name, coordination distance (Å), and max coordination number for TZ injection
struct ag4_tz_injection_params_t {
    std::vector<std::string> metal_names;
    double coord_dist;
    int    max_coord;
};
inline ag4_tz_injection_params_t ag4_tz_injection_params(ag4_metal_mode mode) {
    switch (mode) {
        case ag4_metal_mode::cd: return {{"Cd"}, 2.52, 4};  // tetrahedral Cd²⁺, r(Cd-S)≈2.52 Å
        default: return {{}, 2.5, 4};
    }
}

inline std::vector<ag4_nbp_override> ag4_jt_nbp_overrides(ag4_metal_mode mode) {
    switch (mode) {
        case ag4_metal_mode::cu2_jt:
            return {{"NA","JT",0.25,3.0,12,6}, {"N","JT",0.25,3.0,12,6},
                    {"SA","JT",0.25,2.5,12,6}, {"OA","JT",0.25,1.5,12,6},
                    {"HD","JT", 1.0,0.0,12,6}};
        case ag4_metal_mode::mn3_jt:
            // JT axial (2 sites at 2.28 Å): Mn³⁺ hard → OA > NA >> SA
            // MH equatorial (4 sites at 1.92 Å): same Mn³⁺ selectivity (OA >> NA >> SA)
            // Both included here because ag4_needs_mh_injection(mn3_jt)=false (JT branch
            // takes precedence in the injection guard), so MH overrides must be carried here.
            return {{"OA","JT",0.25, 5.0,12,6}, {"NA","JT",0.25, 2.5,12,6},
                    {"N", "JT",0.25, 2.5,12,6}, {"SA","JT",0.25, 0.2,12,6},
                    {"HD","JT", 1.0, 0.0,12,6},
                    {"OA","MH",0.25,14.0,12,6}, {"NA","MH",0.25, 8.0,12,6},
                    {"N", "MH",0.25, 8.0,12,6}, {"SA","MH",0.25, 0.5,12,6},
                    {"HD","MH", 1.0, 0.0,12,6}};
        default:
            return {};
    }
}

inline bool ag4_needs_jt_injection(ag4_metal_mode mode) {
    return mode == ag4_metal_mode::cu2_jt || mode == ag4_metal_mode::mn3_jt;
}

struct ag4_jt_injection_params_t {
    std::vector<std::string> metal_names;
    std::string equatorial_type;
    double equatorial_dist;
    int    equatorial_coord;
    std::string axial_type;
    double axial_dist;
    int    axial_coord;
};
inline ag4_jt_injection_params_t ag4_jt_injection_params(ag4_metal_mode mode) {
    switch (mode) {
        case ag4_metal_mode::cu2_jt:
            return {{"Cu"}, "SQ", 2.03, 4, "JT", 2.45, 2};
        case ag4_metal_mode::mn3_jt:
            return {{"Mn"}, "MH", 1.92, 4, "JT", 2.28, 2};
        default:
            return {{}, "SQ", 2.05, 4, "JT", 2.35, 2};
    }
}

struct ag4_bridge_water_site {
    vec xyz;
    double weight;
    double target_dist;
    ag4_bridge_water_site() : xyz(0,0,0), weight(0.0), target_dist(2.80) {}
};

struct ag4_metal_site_state {
    ag4_metal_mode mode;
    std::string metal_type;
    vec metal_xyz;
    double coord_dist;
    int max_coord;
    int receptor_cn;
    bool jt_enabled;
    vec jt_axis;
    double jt_axial_dist;
    std::vector<ag4_nbp_override> direct_overrides;
    std::vector<ag4_bridge_water_site> bridge_waters;
    ag4_metal_site_state()
        : mode(ag4_metal_mode::none), metal_xyz(0,0,0), coord_dist(0.0), max_coord(0), receptor_cn(0),
          jt_enabled(false), jt_axis(0,0,1), jt_axial_dist(0.0) {}
};

struct ag4_metal_state {
    std::vector<ag4_metal_site_state> sites;
    bool enabled() const { return !sites.empty(); }
    void clear() { sites.clear(); }
};

struct ag4_inline_result {
    grid_dims                        gd;
    std::vector<sz>                  ad_types;
    std::vector<std::vector<double>> aff_maps;   // [i][flat_idx] x-fastest
    std::vector<double>              elec_map;
    std::vector<double>              desolv_map;
    bool                             valid;
    ag4_inline_result() : valid(false) {}
};

ag4_inline_result ag4_compute_maps(const std::string& receptor_pdbqt_path,
                                    const vec&         center,
                                    const vec&         box_size,
                                    double             spacing,
                                    const std::vector<sz>& ligand_ad_types,
                                    const std::vector<ag4_nbp_override>& nbp_overrides = {},
                                    bool zn_mode = false,
                                    ag4_metal_mode metal_mode = ag4_metal_mode::none,
                                    const std::vector<ag4_metal_mode>& extra_metal_modes = {});

ag4_metal_state ag4_build_metal_state(const std::string& receptor_pdbqt_path,
                                      ag4_metal_mode primary_mode,
                                      const std::vector<ag4_metal_mode>& extra_metal_modes = {});

#endif
