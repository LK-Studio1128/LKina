/*
   LKina — embedded AD4 grid generation (header)
   Copyright (C) 2025 LK-Studio1128 and LKina contributors
   SPDX-License-Identifier: GPL-3.0-or-later
   See COPYING for the full GNU General Public License v3 text.
*/

#ifndef VINA_EMBEDDED_AD4_GRID_H
#define VINA_EMBEDDED_AD4_GRID_H

#include <string>
#include <vector>
#include "common.h"
#include "grid_dim.h"
#include "reactive_types.h"
#include "ag4_engine.h"

struct embedded_ad4_request {
    std::string receptor_pdbqt_path;
    vec         center;
    vec         box_size;          // angstroms
    fl          spacing;           // default 0.375 A
    std::vector<sz> ligand_ad_types;
    std::string parameter_file;    // empty -> use builtin
    std::string autogrid4_exe;     // path for fallback; empty -> skip
    std::string work_dir;          // temp dir for files; empty -> use system tmp
    reactive_payload reactive;     // reactive covalent docking geometry (disabled by default)
    bool zn_mode = false;          // AutoDock4Zn: TZ pseudoatom + AD4Zn pairwise overrides
    ag4_metal_mode metal_mode = ag4_metal_mode::none; // primary metal mode
    std::vector<ag4_metal_mode> extra_metal_modes;    // additional modes (multi-metal sites)
    std::vector<ag4_nbp_override> nbp_overrides; // additional pairwise potential overrides
};

struct ad4_grid_data {
    vec       origin;
    fl        spacing;
    sz        nx, ny, nz;
    grid_dims gd;                           // set for inline path
    std::vector<sz>                  map_ad_types;
    std::vector<std::vector<fl>>     affinity_maps;
    std::vector<fl>                  electrostatic_map;
    std::vector<fl>                  desolvation_map;
    std::string                      map_prefix;  // empty => inline path
    bool      valid;
    reactive_payload reactive;     // forwarded from request; unused unless enabled
};

ad4_grid_data generate_embedded_ad4_maps(const embedded_ad4_request& req);

#endif
