/*
   LKina — Reactive Covalent Docking: type definitions
   All structs default to "disabled" — zero impact on existing code paths.
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#ifndef VINA_REACTIVE_TYPES_H
#define VINA_REACTIVE_TYPES_H

#include <string>
#include <limits>
#include <cmath>
#include "common.h"

enum class reactive_mode {
    off = 0,
    distance,
    hybrid
};

enum class reactive_geometry_mode {
    none = 0,
    line,
    angle
};

// ── User-facing configuration (populated from CLI / API) ─────────────────────
struct ReactiveOptions {
    reactive_mode          mode          = reactive_mode::off;
    reactive_geometry_mode geometry_mode = reactive_geometry_mode::none;

    std::string receptor_atom_spec;        // "chain:resnum:name"  OR  "x,y,z"
    std::string receptor_frame_atom_spec;  // optional second receptor atom for directional term
    std::string ligand_atom_spec;          // 1-based PDBQT serial (e.g. "5") or atom name (e.g. "SG")

    double bond_length         = 0.0;  // 0 = attractor well centred at anchor
    double attractor_width     = 1.5;  // Gaussian well σ (Å)
    double attractor_strength  = 8.0;  // well depth (kcal/mol)
    double angle_strength      = 4.0;  // angular constraint K (kcal/mol), used with frame atom

    double hybrid_vdw_scale = 0.0;   // fraction of VdW+HB grid to keep in hybrid mode (0=fully suppressed, 1=full)
    double target_angle_deg = 180.0; // ideal angle at reactive atom (degrees)
    double angle_width_deg  = 0.0;   // flat-bottom half-width in degrees (0=pure harmonic, >0=flat zone)
    std::string ligand_frame_atom_spec; // ligand-side frame atom — enables approach-angle constraint

    bool debug = false;
    bool debug_energy = false;
    bool gradient_check = false;
    double gradient_check_eps = 1e-4;

    // ── C3 Two-step strategy ─────────────────────────────────────────────────
    // Phase 1: MC presample without reactive energy.
    // Phase 2: filter poses by ligand-anchor distance ≤ presample_max_dist, then
    //          run quasi-newton with full reactive constraints.
    bool   two_step           = false;  // enable C3 two-step strategy
    double presample_max_dist = 10.0;   // phase-1 distance filter threshold (Å)
    bool   weak_attractor     = false;  // C3b: use broad/weak Gaussian in Phase 1 (Goullieux 2023 "attracting cavities")

    // ── Reaction-type preset name (populated from --reactive_preset CLI arg) ──
    // Valid values: "cys_michael", "cys_sn2", "ser_covalent", "lys_targeting",
    //               "boronic_acid", "tyr_covalent"
    // Preset fills in defaults for mode/bond_length/attractor_*/angle_* fields.
    // Any individually specified CLI options override the preset values.
    std::string preset_name;

    // Apply built-in preset defaults. Returns false if preset_name is unknown.
    // Only fills fields that have not been explicitly overridden (caller must
    // track which flags were set and call this before the override pass).
    bool apply_preset() {
        if (preset_name.empty()) return true;
        // ── Cys Michael addition: C=C warhead attacks Cys SG (sp3 product)
        if (preset_name == "cys_michael") {
            mode              = reactive_mode::hybrid;
            geometry_mode     = reactive_geometry_mode::angle;
            bond_length        = 1.82;
            attractor_width    = 0.8;
            attractor_strength = 10.0;
            target_angle_deg   = 109.5;
            angle_width_deg    = 25.0;
            angle_strength     = 4.0;
            hybrid_vdw_scale   = 0.2;
            return true;
        }
        // ── Cys SN2: electrophile attacks Cys SG from the back (180° approach)
        if (preset_name == "cys_sn2") {
            mode              = reactive_mode::hybrid;
            geometry_mode     = reactive_geometry_mode::angle;
            bond_length        = 1.82;
            attractor_width    = 0.8;
            attractor_strength = 10.0;
            target_angle_deg   = 180.0;
            angle_width_deg    = 15.0;
            angle_strength     = 5.0;
            hybrid_vdw_scale   = 0.2;
            return true;
        }
        // ── Ser covalent: acyl/phospho group reacts with Ser OG (beta-lactam etc.)
        if (preset_name == "ser_covalent") {
            mode              = reactive_mode::hybrid;
            geometry_mode     = reactive_geometry_mode::angle;
            bond_length        = 1.34;
            attractor_width    = 0.8;
            attractor_strength = 10.0;
            target_angle_deg   = 109.5;
            angle_width_deg    = 25.0;
            angle_strength     = 4.0;
            hybrid_vdw_scale   = 0.2;
            return true;
        }
        // ── Lys targeting: aldehyde/activated ester forms Schiff base with Lys NZ
        if (preset_name == "lys_targeting") {
            mode              = reactive_mode::hybrid;
            geometry_mode     = reactive_geometry_mode::angle;
            bond_length        = 1.47;
            attractor_width    = 0.8;
            attractor_strength = 8.0;
            target_angle_deg   = 109.5;
            angle_width_deg    = 30.0;
            angle_strength     = 3.0;
            hybrid_vdw_scale   = 0.3;
            return true;
        }
        // ── Boronic acid: reversible covalent with Ser/Thr/Tyr OH (no angle required)
        if (preset_name == "boronic_acid") {
            mode              = reactive_mode::distance;
            geometry_mode     = reactive_geometry_mode::none;
            bond_length        = 1.47;
            attractor_width    = 1.0;
            attractor_strength = 6.0;
            hybrid_vdw_scale   = 0.5;
            return true;
        }
        // ── Tyr covalent: phenol OH attacks electrophile
        if (preset_name == "tyr_covalent") {
            mode              = reactive_mode::hybrid;
            geometry_mode     = reactive_geometry_mode::angle;
            bond_length        = 1.38;
            attractor_width    = 0.8;
            attractor_strength = 10.0;
            target_angle_deg   = 109.5;
            angle_width_deg    = 25.0;
            angle_strength     = 4.0;
            hybrid_vdw_scale   = 0.2;
            return true;
        }
        return false;  // unknown preset
    }

    bool enabled() const { return mode != reactive_mode::off; }
    void clear()         { *this = ReactiveOptions(); }
};

// ── Map-generation result: resolved receptor geometry ────────────────────────
struct reactive_payload {
    bool                   enabled       = false;
    reactive_mode          mode          = reactive_mode::off;
    reactive_geometry_mode geometry_mode = reactive_geometry_mode::none;

    vec  receptor_atom_xyz;
    bool has_receptor_atom  = false;

    vec  receptor_frame_xyz;
    bool has_receptor_frame = false;

    double bond_length        = 0.0;
    double attractor_width    = 0.0;
    double attractor_strength = 0.0;
    double angle_strength     = 0.0;

    double hybrid_vdw_scale = 0.0;
    double target_angle_deg = 180.0;
    double angle_width_deg  = 0.0;   // flat-bottom half-width in degrees (0=pure harmonic, >0=flat zone)

    bool debug = false;
    bool debug_energy = false;
    bool gradient_check = false;
    double gradient_check_eps = 0.0;

    reactive_payload() : receptor_atom_xyz(0,0,0), receptor_frame_xyz(0,0,0) {}
    void clear() { *this = reactive_payload(); }
};

// ── Scoring-time runtime state: fully resolved, handed to ad4cache ───────────
struct reactive_state {
    bool                   enabled       = false;
    reactive_mode          mode          = reactive_mode::off;
    reactive_geometry_mode geometry_mode = reactive_geometry_mode::none;

    sz   ligand_atom_index;    // index into m.atoms / m.coords
    bool has_ligand_atom  = false;

    vec  receptor_atom_xyz;
    bool has_receptor_atom  = false;

    vec  receptor_frame_xyz;
    bool has_receptor_frame = false;

    double bond_length        = 0.0;
    double attractor_width    = 0.0;
    double attractor_strength = 0.0;
    double angle_strength     = 0.0;

    double hybrid_vdw_scale      = 0.0;
    double cos_target_angle = -1.0;  // cos(target_angle_deg), default 180° → -1
    // Flat-bottom boundaries in cosine space derived from exact angle-space bounds:
    //   cos_angle_lo = cos(target + width)  (lower cosine = upper angle boundary)
    //   cos_angle_hi = cos(target - width)  (upper cosine = lower angle boundary)
    // When width=0: cos_angle_lo == cos_angle_hi == cos_target_angle → pure harmonic.
    double cos_angle_lo = -1.0;
    double cos_angle_hi = -1.0;
    sz     ligand_frame_atom_index;       // ligand-side frame atom for approach-angle constraint
    bool   has_ligand_frame  = false;

    bool debug = false;
    bool debug_energy = false;
    bool gradient_check = false;
    double gradient_check_eps = 0.0;

    reactive_state()
        : ligand_atom_index(std::numeric_limits<sz>::max())
        , ligand_frame_atom_index(std::numeric_limits<sz>::max())
        , receptor_atom_xyz(0,0,0)
        , receptor_frame_xyz(0,0,0)
    {}

    // Ready to compute reactive energy — all required fields resolved
    bool ready() const {
        return enabled && has_ligand_atom && has_receptor_atom
            && ligand_atom_index != std::numeric_limits<sz>::max()
            && attractor_width > 0.0;
    }

    void clear() { *this = reactive_state(); }
};

#endif
