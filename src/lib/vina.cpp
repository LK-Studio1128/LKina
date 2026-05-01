/*

   Copyright (c) 2006-2010, The Scripps Research Institute

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

	   http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   Author: Dr. Oleg Trott <ot14@columbia.edu>, 
		   The Olson Lab,
		   The Scripps Research Institute

*/

#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <iomanip>
#include "vina.h"
#include "scoring_function.h"
#include "precalculate.h"
#include "embedded_ad4_grid.h"
#include "parse_pdbqt.h"

// ── file-local: parse receptor anchor for reactive covalent docking ─────────
// spec formats:
//   "chain:resnum:atom_name"  e.g. "A:145:SG"
//   "x,y,z"                  e.g. "12.3,8.7,1.1"
static bool parse_receptor_reactive_atom(const std::string& pdbqt_path,
                                          const std::string& spec,
                                          vec& out_xyz,
                                          std::string* out_resname = 0,
                                          std::string* out_atom_name = 0) {
    if (spec.empty()) return false;
    // coordinate literal
    if (spec.find(',') != std::string::npos) {
        std::istringstream ss(spec);
        double x, y, z; char c;
        if (ss >> x >> c >> y >> c >> z) {
            out_xyz = vec(x, y, z);
            if (out_resname) out_resname->clear();
            if (out_atom_name) out_atom_name->clear();
            return true;
        }
        return false;
    }
    // chain:resnum:atom_name
    std::vector<std::string> parts;
    std::istringstream ss(spec);
    std::string tok;
    while (std::getline(ss, tok, ':')) parts.push_back(tok);
    if (parts.size() != 3 && parts.size() != 4) return false;
    const std::string& chain = parts[0];
    int resnum = 0;
    try { resnum = std::stoi(parts[1]); } catch (...) { return false; }
    const std::string target_icode = (parts.size() == 4) ? parts[2] : "";
    const std::string& target_name = (parts.size() == 4) ? parts[3] : parts[2];
    std::ifstream f(pdbqt_path.c_str());
    if (!f) return false;
    std::string line;
    while (std::getline(f, line)) {
        if (line.size() < 54) continue;
        std::string rec = line.substr(0, 6);
        if (rec.substr(0,4) != "ATOM" && rec != "HETATM") continue;
        std::string aname = line.substr(12, 4);
        size_t p = aname.find_first_not_of(" ");
        if (p == std::string::npos) continue;
        aname = aname.substr(p);
        size_t q = aname.find_last_not_of(" ");
        if (q != std::string::npos) aname = aname.substr(0, q + 1);
        std::string achain = line.substr(21, 1);
        int aresnum = 0;
        try { aresnum = std::stoi(line.substr(22, 4)); } catch (...) { continue; }
        std::string icode = (line.size() > 26) ? line.substr(26, 1) : "";
        if (icode == " ") icode.clear();
        char altloc = (line.size() > 16) ? line[16] : ' ';
        if (achain == chain && aresnum == resnum && (target_icode.empty() || icode == target_icode) && aname == target_name && (altloc == ' ' || altloc == 'A')) {
            try {
                double x = std::stod(line.substr(30, 8));
                double y = std::stod(line.substr(38, 8));
                double z = std::stod(line.substr(46, 8));
                out_xyz = vec(x, y, z);
                if (out_resname) {
                    std::string rname = line.substr(17, 3);
                    size_t rp = rname.find_first_not_of(" ");
                    if (rp != std::string::npos) rname = rname.substr(rp);
                    size_t rq = rname.find_last_not_of(" ");
                    if (rq != std::string::npos) rname = rname.substr(0, rq + 1);
                    *out_resname = rname;
                }
                if (out_atom_name) *out_atom_name = aname;
                return true;
            } catch (...) { continue; }
        }
    }
    return false;
}

static bool reactive_anchor_matches_preset(const std::string& preset,
                                           const std::string& resname,
                                           const std::string& atom_name) {
    if (preset.empty() || resname.empty() || atom_name.empty()) return true;
    if (preset == "cys_michael" || preset == "cys_sn2")
        return resname == "CYS" && atom_name == "SG";
    if (preset == "ser_covalent")
        return resname == "SER" && atom_name == "OG";
    if (preset == "lys_targeting")
        return resname == "LYS" && atom_name == "NZ";
    if (preset == "boronic_acid")
        return (resname == "SER" && atom_name == "OG") ||
               (resname == "THR" && atom_name == "OG1") ||
               (resname == "TYR" && atom_name == "OH");
    if (preset == "tyr_covalent")
        return resname == "TYR" && atom_name == "OH";
    return true;
}
// ───────────────────────────────────────────────────────────────────────────────


void Vina::cite() {
	const std::string cite_message = "\
#################################################################\n\
#                    LKina Docking Engine                       #\n\
#      Metal-Enhanced AutoDock4 Force Field  (AD4 + inline      #\n\
#      AutoGrid, no external autogrid4 required)                #\n\
#                                                               #\n\
# Based on AutoDock Vina v1.2.7 (Eberhardt et al. 2021)        #\n\
# Extended with inline AG4 grid engine for full-metal docking   #\n\
# Supports 80+ metal atom types with AD4 force field            #\n\
#                                                               #\n\
# Usage: LKina --scoring ad4 --generate_maps --receptor r.pdbqt #\n\
#              --center_x X --center_y Y --center_z Z           #\n\
#              --size_x N  --size_y N  --size_z N               #\n\
#              --ligand metal.pdbqt --out out.pdbqt             #\n\
#################################################################\n";

	std::cout << cite_message << '\n';
}

int Vina::generate_seed(const int seed) {
	// Seed generator, if the global seed (m_seed) was defined to 0
	// it seems that we want to generate random seed, otherwise it means
	// that we want to use a particular seed.
	if (seed == 0) {
		return auto_seed();
	} else {
		return seed;
	}
}

void Vina::set_receptor(const std::string& rigid_name, const std::string& flex_name) {
	// Read the receptor PDBQT file
	/* CONDITIONS:
		- 1. AD4/Vina rigid  NO, flex  NO: FAIL
		- 2. AD4      rigid YES, flex YES: FAIL
		- 3. AD4      rigid YES, flex  NO: FAIL
		- 4. AD4      rigid  NO, flex YES: SUCCESS (need to read maps later)
		- 5. Vina     rigid YES, flex YES: SUCCESS
		- 6. Vina     rigid YES, flex  NO: SUCCESS
		- 7. Vina     rigid  NO, flex YES: SUCCESS (need to read maps later)
	*/
	if (rigid_name.empty() && flex_name.empty() && m_sf_choice == SF_VINA) {
		// CONDITION 1
		std::cerr << "ERROR: No (rigid) receptor or flexible residues were specified. (vina.cpp)\n";
		exit(EXIT_FAILURE);
	} else if (m_sf_choice == SF_AD42 && !rigid_name.empty()) {
		// CONDITIONS 2, 3
		std::cerr << "ERROR: Only flexible residues allowed with the AD4 scoring function. No (rigid) receptor.\n";
		exit(EXIT_FAILURE);
	}

	// CONDITIONS 4, 5, 6, 7 (rigid_name and flex_name are empty strings per default)
	m_receptor = parse_receptor_pdbqt(rigid_name, flex_name, m_scoring_function->get_atom_typing());

	m_model = m_receptor;
	m_receptor_initialized = true;
	// If we are reading another receptor we should not consider the ligand and the map as initialized anymore
	m_ligand_initialized = false;
	m_map_initialized = false;
}

void Vina::set_ligand_from_string(const std::string& ligand_string) {
	// Read ligand PDBQT string and add it to the model
	if (ligand_string.empty()) {
		std::cerr << "ERROR: Cannot read ligand file. Ligand string is empty.\n";
		exit(EXIT_FAILURE);
	}

	atom_type::t atom_typing = m_scoring_function->get_atom_typing();

	if (!m_receptor_initialized) {
		// This situation will happen if we don't need a receptor and we are using affinity maps
		model m(atom_typing);
		m_model = m;
		m_receptor = m;
	} else {
		// Replace current model with receptor and reinitialize poses
		m_model = m_receptor;
	}

	// ... and add ligand to the model
	m_model.append(parse_ligand_pdbqt_from_string(ligand_string, atom_typing));

	// Because we precalculate ligand atoms interactions
	precalculate_byatom precalculated_byatom(*m_scoring_function, m_model);

	// Check that all atom types are in the grid (if initialized)
	if (m_map_initialized) {
		szv atom_types = m_model.get_movable_atom_types(atom_typing);

		if (m_sf_choice == SF_VINA || m_sf_choice == SF_VINARDO) {
			if(!m_grid.are_atom_types_grid_initialized(atom_types))
				exit(EXIT_FAILURE);
		} else {
			if(!m_ad4grid.are_atom_types_grid_initialized(atom_types))
				exit(EXIT_FAILURE);
		}
	}

	// Store in Vina object
	output_container poses;
	m_poses = poses;
	m_precalculated_byatom = precalculated_byatom;
	m_ligand_initialized = true;
	finalize_reactive_state_if_possible();
}

void Vina::set_ligand_from_string(const std::vector<std::string>& ligand_string) {
	// Read ligand PDBQT strings and add them to the model
	if (ligand_string.empty()) {
		std::cerr << "ERROR: Cannot read ligand list. Ligands list is empty.\n";
		exit(EXIT_FAILURE);
	}

	atom_type::t atom_typing = m_scoring_function->get_atom_typing();

	if (!m_receptor_initialized) {
		// This situation will happen if we don't need a receptor and we are using affinity maps
		model m(atom_typing);
		m_model = m;
		m_receptor = m;
	} else {
		// Replace current model with receptor and reinitialize poses
		m_model = m_receptor;
	}

	VINA_RANGE(i, 0, ligand_string.size())
		m_model.append(parse_ligand_pdbqt_from_string(ligand_string[i], atom_typing));

	// Because we precalculate ligand atoms interactions
	precalculate_byatom precalculated_byatom(*m_scoring_function, m_model);

	// Check that all atom types are in the grid (if initialized)
	if (m_map_initialized) {
		szv atom_types = m_model.get_movable_atom_types(atom_typing);

		if (m_sf_choice == SF_VINA || m_sf_choice == SF_VINARDO) {
			if(!m_grid.are_atom_types_grid_initialized(atom_types))
				exit(EXIT_FAILURE);
		} else {
			if(!m_ad4grid.are_atom_types_grid_initialized(atom_types))
				exit(EXIT_FAILURE);
		}
	}

	// Store in Vina object
	output_container poses;
	m_poses = poses;
	m_precalculated_byatom = precalculated_byatom;
	m_ligand_initialized = true;
	finalize_reactive_state_if_possible();
}

void Vina::set_ligand_from_file(const std::string& ligand_name) {
	set_ligand_from_string(get_file_contents(ligand_name));
}

void Vina::set_ligand_from_file(const std::vector<std::string>& ligand_name) {
	std::vector<std::string> ligand_string;

	VINA_RANGE(i, 0, ligand_name.size())
		ligand_string.push_back(get_file_contents(ligand_name[i]));

	set_ligand_from_string(ligand_string);
}

/*
void Vina::set_ligand(OpenBabel::OBMol* mol) {
	// Add OBMol to the model
	OpenBabel::OBConversion conv;
	conv.SetOutFormat("PDBQT");
	set_ligand_from_string(conv.WriteString(mol));
}

void Vina::set_ligand(std::vector<OpenBabel::OBMol*> mol) {
	// Add OBMols to the model
	std::vector<std::string> ligand_string;

	OpenBabel::OBConversion conv;
	conv.SetOutFormat("PDBQT");

	VINA_RANGE(i, 0, ligand_name.size())
		ligand_string.push_back(conv.WriteString(mol[i]));

	set_ligand_from_string(ligand_string);
}
*/

void Vina::set_vina_weights(double weight_gauss1, double weight_gauss2, double weight_repulsion,
							double weight_hydrophobic, double weight_hydrogen, double weight_glue,
							double weight_rot) {
	flv weights;

	if (m_sf_choice == SF_VINA) {
		weights.push_back(weight_gauss1);
		weights.push_back(weight_gauss2);
		weights.push_back(weight_repulsion);
		weights.push_back(weight_hydrophobic);
		weights.push_back(weight_hydrogen);
		weights.push_back(weight_glue);
		weights.push_back(5 * weight_rot / 0.1 - 1);

		// Store in Vina object
		m_weights = weights;

		// Since we set (different) weights, we automatically initialize the forcefield
		set_forcefield();
	}
}

void Vina::set_vinardo_weights(double weight_gauss1, double weight_repulsion,
							   double weight_hydrophobic, double weight_hydrogen, double weight_glue,
							   double weight_rot) {
	flv weights;

	if (m_sf_choice == SF_VINARDO) {
		weights.push_back(weight_gauss1);
		weights.push_back(weight_repulsion);
		weights.push_back(weight_hydrophobic);
		weights.push_back(weight_hydrogen);
		weights.push_back(weight_glue);
		weights.push_back(5 * weight_rot / 0.1 - 1);

		// Store in Vina object
		m_weights = weights;

		// Since we set (different) weights, we automatically initialize the forcefield
		set_forcefield();
	}
}

void Vina::set_ad4_weights(double weight_ad4_vdw , double weight_ad4_hb,
						   double weight_ad4_elec, double weight_ad4_dsolv,
						   double weight_glue, double weight_ad4_rot) {
	flv weights;

	if (m_sf_choice == SF_AD42) {
		weights.push_back(weight_ad4_vdw);
		weights.push_back(weight_ad4_hb);
		weights.push_back(weight_ad4_elec);
		weights.push_back(weight_ad4_dsolv);
		weights.push_back(weight_glue);
		weights.push_back(weight_ad4_rot);

		// Store in Vina object
		m_weights = weights;

		// Since we set (different) weights, we automatically initialize the forcefield
		set_forcefield();
	}
}

void Vina::set_forcefield() {
    // Store in Vina object
    m_scoring_function = std::make_shared<ScoringFunction>(m_sf_choice, m_weights);
}

std::vector<double> Vina::grid_dimensions_from_ligand(double buffer_size) {
	std::vector<double> box_dimensions(6, 0);
	std::vector<double> box_center(3, 0);
	std::vector<double> max_distance(3, 0);

	// The center of the ligand will be the center of the box
	box_center = m_model.center();

	// Get the furthest atom coordinates from the center in each dimensions
	VINA_FOR(i, m_model.num_movable_atoms()) {
		const vec& atom_coords = m_model.get_coords(i);

		VINA_FOR_IN(j, atom_coords) {
			double distance = std::fabs(box_center[j] - atom_coords[j]);

			if (max_distance[j] < distance)
				max_distance[j] = distance;
		}
	}

	// Get the final dimensions of the box
	box_dimensions[0] = box_center[0];
	box_dimensions[1] = box_center[1];
	box_dimensions[2] = box_center[2];
	box_dimensions[3] = std::ceil((max_distance[0] + buffer_size) * 2);
	box_dimensions[4] = std::ceil((max_distance[1] + buffer_size) * 2);
	box_dimensions[5] = std::ceil((max_distance[2] + buffer_size) * 2);

	return box_dimensions;
}

void Vina::compute_vina_maps(double center_x, double center_y, double center_z, double size_x, double size_y, double size_z, double granularity, bool force_even_voxels) {
	// Setup the search box
	// Check first that the receptor was added
	if (m_sf_choice == SF_AD42) {
		throw vina_runtime_error("Cannot compute Vina affinity maps using the AD4 scoring function.");
	} else if (!m_receptor_initialized) {
		// m_model
		throw vina_runtime_error("Cannot compute Vina or Vinardo affinity maps. The (rigid) receptor was not initialized.");
	} else if (size_x <= 0 || size_y <= 0 || size_z <= 0) {
		throw vina_runtime_error("Grid box dimensions must be greater than 0 Angstrom");
	} else if (size_x * size_y * size_z > 27e3) {
		std::cerr << "WARNING: Search space volume is greater than 27000 Angstrom^3 (See FAQ)\n";
	}

	grid_dims gd;
	vec span(size_x, size_y, size_z);
	vec center(center_x, center_y, center_z);
	const fl slope = 1e6; // FIXME: too large? used to be 100
	szv atom_types;
	atom_type::t atom_typing = m_scoring_function->get_atom_typing();

	/* Atom types initialization
	If a ligand was defined before, we only use those present in the ligand
	otherwise we use all the atom types present in the forcefield
	*/
	if (m_ligand_initialized)
		atom_types = m_model.get_movable_atom_types(atom_typing);
	else
		atom_types = m_scoring_function->get_atom_types();

	// Grid dimensions
	VINA_FOR_IN(i, gd) {
		gd[i].n_voxels = sz(std::ceil(span[i] / granularity));

		// If odd n_voxels increment by 1
		if (force_even_voxels && (gd[i].n_voxels % 2 == 1))
			// because sample points (npts) == n_voxels + 1
			gd[i].n_voxels += 1;

		fl real_span = granularity * gd[i].n_voxels;
		gd[i].begin = center[i] - real_span / 2;
		gd[i].end = gd[i].begin + real_span;
	}

	// Initialize the scoring function
	precalculate precalculated_sf(*m_scoring_function);
	// Store it now in Vina object because of non_cache
	m_precalculated_sf = precalculated_sf;

	if (m_sf_choice == SF_VINA)
		doing("Computing Vina grid", m_verbosity, 0);
	else
		doing("Computing Vinardo grid", m_verbosity, 0);

	// Compute the Vina grids
	cache grid(gd, slope);
	grid.populate(m_model, precalculated_sf, atom_types);

	done(m_verbosity, 0);

	// create non_cache for scoring with explicit receptor atoms (instead of grids)
	if (!m_no_refine) {
		non_cache nc(m_model, gd, &m_precalculated_sf, slope);
		m_non_cache = nc;
	}

	// Store in Vina object
	m_grid = grid;
	m_map_initialized = true;
}

void Vina::compute_ad4_maps(const std::string& receptor_pdbqt_path,
                            double center_x, double center_y, double center_z,
                            double size_x, double size_y, double size_z,
                            const std::string& autogrid4_exe,
                            double granularity,
                            const std::string& map_output_dir) {
	if (m_sf_choice != SF_AD42)
		throw vina_runtime_error("compute_ad4_maps requires the AD4 scoring function.");
	if (receptor_pdbqt_path.empty())
		throw vina_runtime_error("compute_ad4_maps: receptor PDBQT path is empty.");
	if (size_x <= 0 || size_y <= 0 || size_z <= 0)
		throw vina_runtime_error("Grid box dimensions must be greater than 0 Angstrom");

	szv ligand_types;
	if (m_ligand_initialized) {
		ligand_types = m_model.get_movable_atom_types(atom_type::AD);
	} else {
		for (sz t = 0; t < AD_TYPE_SIZE; ++t) {
			switch (t) {
				case AD_TYPE_G0: case AD_TYPE_G1: case AD_TYPE_G2: case AD_TYPE_G3:
				case AD_TYPE_CG0: case AD_TYPE_CG1: case AD_TYPE_CG2: case AD_TYPE_CG3:
				case AD_TYPE_W: continue;
				default: break;
			}
			ligand_types.push_back(t);
		}
	}

	embedded_ad4_request req;
	req.receptor_pdbqt_path = receptor_pdbqt_path;
	req.center  = vec(center_x, center_y, center_z);
	req.box_size = vec(size_x, size_y, size_z);
	req.spacing = (fl)granularity;
	req.ligand_ad_types = ligand_types;
	req.autogrid4_exe   = autogrid4_exe;
	req.work_dir        = map_output_dir;
	req.zn_mode           = m_zn_mode;
	req.metal_mode        = m_metal_mode;
	req.extra_metal_modes = m_extra_metal_modes;
	m_ad4_receptor_path   = receptor_pdbqt_path;  // remember for reactive reuse via --maps

	doing("Running embedded AD4 map generation", m_verbosity, 0);
	ad4_grid_data result = generate_embedded_ad4_maps(req);
	done(m_verbosity, 0);

	if (!result.valid)
		throw vina_runtime_error("compute_ad4_maps: map generation failed (autogrid4 not available or failed).");

	m_last_ad4_map_prefix = result.map_prefix;

	doing("Loading AD4 affinity maps", m_verbosity, 0);
	if (result.map_prefix.empty()) {
		// Inline computation path: populate ad4cache directly from memory
		m_ad4grid.populate_from_data(result.gd,
		                             result.map_ad_types,
		                             std::vector<std::vector<double>>(result.affinity_maps.begin(), result.affinity_maps.end()),
		                             std::vector<double>(result.electrostatic_map.begin(), result.electrostatic_map.end()),
		                             std::vector<double>(result.desolvation_map.begin(), result.desolvation_map.end()));
	} else {
		// Autogrid4 fallback path: read from .map files
		m_ad4grid.read(result.map_prefix);
	}
	done(m_verbosity, 0);
	m_ad4grid.set_metal_state(ag4_build_metal_state(receptor_pdbqt_path, m_metal_mode, m_extra_metal_modes));
	if (m_metal_soft_weight > 0.0f) m_ad4grid.set_metal_soft_weight(m_metal_soft_weight);

	if (m_ligand_initialized) {
		szv atom_types = m_model.get_movable_atom_types(atom_type::AD);
		if (!m_ad4grid.are_atom_types_grid_initialized(atom_types))
			exit(EXIT_FAILURE);
	}

	m_map_initialized = true;

	// ── Reactive: resolve receptor anchor coordinates (no-op when not enabled) ──
	if (m_reactive_enabled) {
		reactive_payload pl;
		pl.enabled            = true;
		pl.mode               = m_reactive_options.mode;
		pl.geometry_mode      = m_reactive_options.geometry_mode;
		pl.bond_length        = m_reactive_options.bond_length;
		pl.attractor_width    = m_reactive_options.attractor_width;
		pl.attractor_strength = m_reactive_options.attractor_strength;
		pl.angle_strength     = m_reactive_options.angle_strength;
		pl.hybrid_vdw_scale   = m_reactive_options.hybrid_vdw_scale;
		pl.target_angle_deg   = m_reactive_options.target_angle_deg;
		pl.angle_width_deg    = m_reactive_options.angle_width_deg;
		pl.debug              = m_reactive_options.debug;
		pl.debug_energy       = m_reactive_options.debug_energy;
		pl.gradient_check     = m_reactive_options.gradient_check;
		pl.gradient_check_eps = m_reactive_options.gradient_check_eps;
		vec anchor(0,0,0);
		std::string anchor_resname, anchor_atom_name;
		if (parse_receptor_reactive_atom(receptor_pdbqt_path,
		                                  m_reactive_options.receptor_atom_spec,
		                                  anchor,
		                                  &anchor_resname,
		                                  &anchor_atom_name)) {
			pl.receptor_atom_xyz = anchor;
			pl.has_receptor_atom = true;
			if (!reactive_anchor_matches_preset(m_reactive_options.preset_name, anchor_resname, anchor_atom_name))
				std::cerr << "WARNING: reactive preset '" << m_reactive_options.preset_name
				          << "' is unusual for receptor atom " << anchor_resname << ":" << anchor_atom_name << ".\n";
		} else {
			throw vina_runtime_error("reactive: could not resolve receptor atom spec '" + m_reactive_options.receptor_atom_spec + "'.");
		}
		if (!m_reactive_options.receptor_frame_atom_spec.empty()) {
			vec frame_xyz(0,0,0);
			if (parse_receptor_reactive_atom(receptor_pdbqt_path,
			                                  m_reactive_options.receptor_frame_atom_spec,
			                                  frame_xyz)) {
				pl.receptor_frame_xyz = frame_xyz;
				pl.has_receptor_frame = true;
				pl.geometry_mode      = reactive_geometry_mode::angle;
			} else {
				throw vina_runtime_error("reactive: could not resolve frame atom spec '" + m_reactive_options.receptor_frame_atom_spec + "'.");
			}
		}
		m_reactive_payload       = pl;
		m_reactive_payload_ready = pl.has_receptor_atom;
		finalize_reactive_state_if_possible();
	}
}

void Vina::load_maps(std::string maps) {
	const fl slope = 1e6; // FIXME: too large? used to be 100
	grid_dims gd;

	if (m_sf_choice == SF_VINA || m_sf_choice == SF_VINARDO) {
		doing("Reading Vina maps", m_verbosity, 0);
		cache grid(slope);
		grid.read(maps);
		done(m_verbosity, 0);
		m_grid = grid;
	} else {
		doing("Reading AD4.2 maps", m_verbosity, 0);
		ad4cache grid(slope);
		grid.read(maps);
		done(m_verbosity, 0);
		m_ad4grid = grid;
	}

	// Check that all the affinity map are present for ligands/flex residues (if initialized already)
	if (m_ligand_initialized) {
		atom_type::t atom_typing = m_scoring_function->get_atom_typing();
		szv atom_types = m_model.get_movable_atom_types(atom_typing);

		if (m_sf_choice == SF_VINA || m_sf_choice == SF_VINARDO) {
			if(!m_grid.are_atom_types_grid_initialized(atom_types))
				exit(EXIT_FAILURE);
		} else {
			if(!m_ad4grid.are_atom_types_grid_initialized(atom_types))
				exit(EXIT_FAILURE);
		}
	}

	// Store in Vina object
	m_map_initialized = true;

	// ── C2: Reactive anchor resolution when loading pre-built --maps ──
	// If reactive is enabled and we have a receptor path (from --receptor or a previous
	// compute_ad4_maps call), resolve the anchor now — same logic as compute_ad4_maps.
	if (m_reactive_enabled && m_ad4_receptor_path.empty() && !m_reactive_payload_ready)
		throw vina_runtime_error("reactive docking with pre-built --maps requires --receptor for anchor resolution.");
	if (m_reactive_enabled && !m_ad4_receptor_path.empty() && !m_reactive_payload_ready) {
		reactive_payload pl;
		pl.enabled            = true;
		pl.mode               = m_reactive_options.mode;
		pl.geometry_mode      = m_reactive_options.geometry_mode;
		pl.bond_length        = m_reactive_options.bond_length;
		pl.attractor_width    = m_reactive_options.attractor_width;
		pl.attractor_strength = m_reactive_options.attractor_strength;
		pl.angle_strength     = m_reactive_options.angle_strength;
		pl.hybrid_vdw_scale   = m_reactive_options.hybrid_vdw_scale;
		pl.target_angle_deg   = m_reactive_options.target_angle_deg;
		pl.angle_width_deg    = m_reactive_options.angle_width_deg;
		pl.debug              = m_reactive_options.debug;
		pl.debug_energy       = m_reactive_options.debug_energy;
		pl.gradient_check     = m_reactive_options.gradient_check;
		pl.gradient_check_eps = m_reactive_options.gradient_check_eps;
		vec anchor(0,0,0);
		std::string anchor_resname, anchor_atom_name;
		if (parse_receptor_reactive_atom(m_ad4_receptor_path,
		                                  m_reactive_options.receptor_atom_spec,
		                                  anchor,
		                                  &anchor_resname,
		                                  &anchor_atom_name)) {
			pl.receptor_atom_xyz = anchor;
			pl.has_receptor_atom = true;
			if (!reactive_anchor_matches_preset(m_reactive_options.preset_name, anchor_resname, anchor_atom_name))
				std::cerr << "WARNING: reactive preset '" << m_reactive_options.preset_name
				          << "' is unusual for receptor atom " << anchor_resname << ":" << anchor_atom_name << ".\n";
		} else {
			throw vina_runtime_error("reactive (--maps): could not resolve receptor atom spec '" + m_reactive_options.receptor_atom_spec + "'.");
		}
		if (!m_reactive_options.receptor_frame_atom_spec.empty()) {
			vec frame_xyz(0,0,0);
			if (parse_receptor_reactive_atom(m_ad4_receptor_path,
			                                  m_reactive_options.receptor_frame_atom_spec,
			                                  frame_xyz)) {
				pl.receptor_frame_xyz = frame_xyz;
				pl.has_receptor_frame = true;
				pl.geometry_mode      = reactive_geometry_mode::angle;
			} else {
				throw vina_runtime_error("reactive (--maps): could not resolve frame atom spec '" + m_reactive_options.receptor_frame_atom_spec + "'.");
			}
		}
		m_reactive_payload       = pl;
		m_reactive_payload_ready = pl.has_receptor_atom;
		finalize_reactive_state_if_possible();
	}
}

void Vina::write_maps(const std::string& map_prefix, const std::string& gpf_filename,
					  const std::string& fld_filename, const std::string& receptor_filename) {
	if (!m_map_initialized) {
		throw vina_runtime_error("Cannot write affinity maps. Affinity maps were not initialized.");
	}

	szv atom_types;
	atom_type::t atom_typing = m_scoring_function->get_atom_typing();

	if (m_ligand_initialized)
		atom_types = m_model.get_movable_atom_types(atom_typing);
	else
		atom_types = m_scoring_function->get_atom_types();

	if (m_sf_choice == SF_VINA || m_sf_choice == SF_VINARDO) {
		doing("Writing Vina maps", m_verbosity, 0);
		m_grid.write(map_prefix, atom_types, gpf_filename, fld_filename, receptor_filename);
		done(m_verbosity, 0);
	} else {
		// Add electrostatics and desolvation maps
		atom_types.push_back(AD_TYPE_SIZE);
		atom_types.push_back(AD_TYPE_SIZE + 1);
		doing("Writing AD4.2 maps", m_verbosity, 0);
		m_ad4grid.write(map_prefix, atom_types, gpf_filename, fld_filename, receptor_filename);
		done(m_verbosity, 0);
	}
}

std::vector< std::vector<double> > Vina::get_poses_coordinates(int how_many, double energy_range) {
	int n = 0;
	double best_energy = 0;
	std::vector< std::vector<double> > coordinates;

	if (how_many < 0) {
		throw vina_runtime_error("number of poses asked must be greater than zero.");
	}

	if (energy_range < 0) {
		throw vina_runtime_error("energy range must be greater than zero.");
	}

	if (!m_poses.empty()) {
		// Get energy from the best conf
		best_energy = m_poses[0].e;

		VINA_FOR_IN(i, m_poses) {
			/* Stop if:
				- We wrote the number of conf asked
				- If there is no conf to write
				- The energy of the current conf is superior than best_energy + energy_range
			*/
			if (n >= how_many || !not_max(m_poses[i].e) || m_poses[i].e > best_energy + energy_range)
				break; // check energy_range sanity FIXME

			// Push the current pose to model
			m_model.set(m_poses[i].c);
			coordinates.push_back(m_model.get_ligand_coords());

			n++;
		}

		// Push back the best conf in model
		m_model.set(m_poses[0].c);
	} else {
		std::cerr << "WARNING: Could not find any pose coordinaates.\n";
	}

	return coordinates;
}

std::vector< std::vector<double> > Vina::get_poses_energies(int how_many, double energy_range) {
	int n = 0;
	double best_energy = 0;
	std::vector< std::vector<double> > energies;

	if (how_many < 0) {
		throw vina_runtime_error("number of poses asked must be greater than zero.");
	}

	if (energy_range < 0) {
		throw vina_runtime_error("energy range must be greater than zero.");
	}

	if (!m_poses.empty()) {
		// Get energy from the best conf
		best_energy = m_poses[0].e;

		VINA_FOR_IN(i, m_poses) {
			/* Stop if:
				- We wrote the number of conf asked
				- If there is no conf to write
				- The energy of the current conf is superior than best_energy + energy_range
			*/
			if (n >= how_many || !not_max(m_poses[i].e) || m_poses[i].e > best_energy + energy_range)
				break; // check energy_range sanity FIXME

			// Push the current pose to model
			energies.push_back({m_poses[i].e,
								m_poses[i].inter, m_poses[i].intra,
								m_poses[i].conf_independent, m_poses[i].unbound});

			n++;
		}
	} else {
		std::cerr << "WARNING: Could not find any pose energies.\n";
	}

	return energies;
}

// E1: parse receptor PDBQT once and cache metal atom positions
void Vina::load_metal_cache_if_needed() {
    if (m_metal_cache_loaded || m_ad4_receptor_path.empty()) return;
    static const char* metal_types[] = {
        "Mg","Ca","Mn","Fe","Co","Ni","Cu","Zn",
        "Pt","Pd","Ru","Ir","Au","Cd","Hg","Na","K", nullptr
    };
    std::ifstream f(m_ad4_receptor_path);
    std::string line;
    while (std::getline(f, line)) {
        if (line.size() < 54) continue;
        std::string rec = line.substr(0,6);
        if (rec != "ATOM  " && rec != "HETATM") continue;
        double x=0, y=0, z=0;
        try {
            x = std::stod(line.substr(30,8));
            y = std::stod(line.substr(38,8));
            z = std::stod(line.substr(46,8));
        } catch (...) { continue; }
        std::istringstream iss(line);
        std::string tok, last;
        while (iss >> tok) last = tok;
        for (int k = 0; metal_types[k]; k++)
            if (last == metal_types[k]) { m_metal_cache.push_back({last,x,y,z}); break; }
    }
    m_metal_cache_loaded = true;
}

struct ligand_metal_candidate {
	std::string type;
	sz index;
	vec xyz;
	int expected_cn;
	bool square_planar;
};

struct ligand_coord_candidate {
	std::string type;
	vec unit;
	fl dist;
	fl target;
};

static bool ligand_metal_geometry_info(sz ad, std::string& type, int& expected_cn, bool& square_planar) {
	if (ad == AD_TYPE_Pt) { type = "Pt"; expected_cn = 4; square_planar = true;  return true; }
	if (ad == AD_TYPE_Pd) { type = "Pd"; expected_cn = 4; square_planar = true;  return true; }
	if (ad == AD_TYPE_Ru) { type = "Ru"; expected_cn = 6; square_planar = false; return true; }
	if (ad == AD_TYPE_Os) { type = "Os"; expected_cn = 6; square_planar = false; return true; }
	if (ad == AD_TYPE_Re) { type = "Re"; expected_cn = 6; square_planar = false; return true; }
	return false;
}

static bool ligand_coord_target(sz ad, bool square_planar, std::string& type, fl& target) {
	if (ad == AD_TYPE_NA || ad == AD_TYPE_N) { type = (ad == AD_TYPE_NA ? "NA" : "N"); target = square_planar ? 2.05 : 2.10; return true; }
	if (ad == AD_TYPE_OA || ad == AD_TYPE_O) { type = (ad == AD_TYPE_OA ? "OA" : "O"); target = square_planar ? 2.00 : 2.08; return true; }
	if (ad == AD_TYPE_SA || ad == AD_TYPE_S) { type = (ad == AD_TYPE_SA ? "SA" : "S"); target = square_planar ? 2.28 : 2.35; return true; }
	if (ad == AD_TYPE_P)  { type = "P";  target = square_planar ? 2.25 : 2.30; return true; }
	if (ad == AD_TYPE_F)  { type = "F";  target = square_planar ? 1.95 : 2.05; return true; }
	if (ad == AD_TYPE_Cl) { type = "Cl"; target = square_planar ? 2.32 : 2.40; return true; }
	if (ad == AD_TYPE_Br) { type = "Br"; target = square_planar ? 2.45 : 2.52; return true; }
	if (ad == AD_TYPE_I)  { type = "I";  target = square_planar ? 2.65 : 2.70; return true; }
	return false;
}

static fl ligand_vec_dot(const vec& a, const vec& b) {
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static vec ligand_vec_cross(const vec& a, const vec& b) {
	return vec(a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]);
}

static fl ligand_vec_norm(const vec& a) {
	return std::sqrt(ligand_vec_dot(a, a));
}

static vec ligand_vec_unit(const vec& a) {
	fl n = ligand_vec_norm(a);
	if (n <= (fl)1e-8) return vec(0,0,0);
	return vec(a[0]/n, a[1]/n, a[2]/n);
}

fl Vina::ligand_metal_geometry_penalty(const model& m, std::string* detail) const {
	if (detail) detail->clear();
	std::vector<ligand_metal_candidate> metals;
	VINA_FOR(i, m.num_movable_atoms()) {
		std::string type;
		int expected_cn = 0;
		bool square_planar = false;
		sz ad = m.get_atom(i).get(atom_type::AD);
		if (!ligand_metal_geometry_info(ad, type, expected_cn, square_planar)) continue;
		metals.push_back({type, i, m.movable_coords(i), expected_cn, square_planar});
	}
	if (metals.empty()) return 0;

	fl total_penalty = 0;
	std::ostringstream ds;
	ds.setf(std::ios::fixed, std::ios::floatfield);
	ds.setf(std::ios::showpoint);
	for (const auto& metal : metals) {
		std::vector<ligand_coord_candidate> coord;
		VINA_FOR(j, m.num_movable_atoms()) {
			if (j == metal.index) continue;
			std::string dtype;
			fl target = 0;
			sz ad = m.get_atom(j).get(atom_type::AD);
			if (!ligand_coord_target(ad, metal.square_planar, dtype, target)) continue;
			const vec& xyz = m.movable_coords(j);
			vec d(xyz[0] - metal.xyz[0], xyz[1] - metal.xyz[1], xyz[2] - metal.xyz[2]);
			fl dist = ligand_vec_norm(d);
			if (dist < (fl)0.5 || dist > target + (fl)0.75) continue;
			coord.push_back({dtype, ligand_vec_unit(d), dist, target});
		}
		std::sort(coord.begin(), coord.end(), [](const ligand_coord_candidate& a, const ligand_coord_candidate& b) {
			return a.dist < b.dist;
		});
		std::size_t used_n = std::min<std::size_t>(coord.size(), metal.expected_cn);
		fl site_penalty = 0;
		if ((int)coord.size() < metal.expected_cn)
			site_penalty += (fl)(metal.expected_cn - (int)coord.size()) * (fl)1.50;
		if ((int)coord.size() > metal.expected_cn)
			site_penalty += (fl)((int)coord.size() - metal.expected_cn) * (fl)0.35;
		for (std::size_t i = 0; i < used_n; i++) {
			fl dr = coord[i].dist - coord[i].target;
			fl z = dr / (fl)0.22;
			site_penalty += std::min((fl)4.0, z*z) * (fl)0.20;
		}
		if (metal.square_planar && used_n >= 3) {
			int trans_pairs = 0;
			std::vector<bool> used(used_n, false);
			for (std::size_t i = 0; i < used_n; i++) {
				if (used[i]) continue;
				fl best_dot = 1;
				std::size_t best_j = used_n;
				for (std::size_t j = i + 1; j < used_n; j++) {
					if (used[j]) continue;
					fl dot = ligand_vec_dot(coord[i].unit, coord[j].unit);
					if (dot < best_dot) { best_dot = dot; best_j = j; }
				}
				if (best_j < used_n && best_dot < (fl)-0.72) {
					used[i] = true;
					used[best_j] = true;
					trans_pairs++;
				}
			}
			site_penalty += (fl)std::max(0, 2 - trans_pairs) * (fl)0.90;
			if (used_n >= 4) {
				fl best_plane = 10;
				for (std::size_t i = 0; i < used_n; i++) {
					for (std::size_t j = i + 1; j < used_n; j++) {
						vec n = ligand_vec_cross(coord[i].unit, coord[j].unit);
						fl nn = ligand_vec_norm(n);
						if (nn < (fl)0.15) continue;
						n = ligand_vec_unit(n);
						fl mean_abs = 0;
						for (std::size_t k = 0; k < used_n; k++)
							mean_abs += std::fabs(ligand_vec_dot(n, coord[k].unit));
						mean_abs /= (fl)used_n;
						best_plane = std::min(best_plane, mean_abs);
					}
				}
				if (best_plane < 10) site_penalty += best_plane * (fl)2.0;
			}
		} else if (!metal.square_planar && used_n >= 4) {
			int trans_pairs = 0;
			std::vector<bool> used(used_n, false);
			for (std::size_t i = 0; i < used_n; i++) {
				if (used[i]) continue;
				fl best_dot = 1;
				std::size_t best_j = used_n;
				for (std::size_t j = i + 1; j < used_n; j++) {
					if (used[j]) continue;
					fl dot = ligand_vec_dot(coord[i].unit, coord[j].unit);
					if (dot < best_dot) { best_dot = dot; best_j = j; }
				}
				if (best_j < used_n && best_dot < (fl)-0.72) {
					used[i] = true;
					used[best_j] = true;
					trans_pairs++;
				}
			}
			site_penalty += (fl)std::max(0, 3 - trans_pairs) * (fl)0.75;
			for (std::size_t i = 0; i < used_n; i++) {
				for (std::size_t j = i + 1; j < used_n; j++) {
					fl dot = ligand_vec_dot(coord[i].unit, coord[j].unit);
					if (dot < (fl)-0.72) continue;
					site_penalty += std::max((fl)0.0, std::fabs(dot) - (fl)0.35) * (fl)0.12;
				}
			}
		}
		const char* quality = site_penalty < (fl)1.0 ? "good" : (site_penalty < (fl)3.0 ? "fair" : "poor");
		if (detail) {
			ds << "REMARK LIGAND_METAL_SITE: " << std::left << std::setw(3) << metal.type
			   << "  cn=" << std::right << coord.size() << "/" << metal.expected_cn
			   << "  geom=" << (metal.square_planar ? "square_planar" : "octahedral")
			   << "  penalty=" << std::setw(8) << std::setprecision(3) << site_penalty
			   << "  quality=" << quality << "\n";
		}
		total_penalty += site_penalty;
	}
	if (detail) *detail = ds.str();
	return total_penalty;
}

std::string Vina::vina_remarks(output_type &pose, fl lb, fl ub) {
	std::ostringstream remark;

	remark.setf(std::ios::fixed, std::ios::floatfield);
	remark.setf(std::ios::showpoint);

	remark << "REMARK VINA RESULT: "
		   << std::setw(9) << std::setprecision(3) << pose.e
		   << "  " << std::setw(9) << std::setprecision(3) << lb
		   << "  " << std::setw(9) << std::setprecision(3) << ub
		   << '\n';

	remark << "REMARK INTER + INTRA:    " << std::setw(12) << std::setprecision(3) << pose.total << "\n";
	remark << "REMARK INTER:            " << std::setw(12) << std::setprecision(3) << pose.inter << "\n";
	remark << "REMARK INTRA:            " << std::setw(12) << std::setprecision(3) << pose.intra << "\n";
	if (m_sf_choice == SF_AD42)
		remark << "REMARK CONF_INDEPENDENT: " << std::setw(12) << std::setprecision(3) << pose.conf_independent << "\n";
	remark << "REMARK UNBOUND:          " << std::setw(12) << std::setprecision(3) << pose.unbound << "\n";
	if (m_sf_choice == SF_AD42 && pose.metal_rerank != 0) {
		remark << "REMARK METAL_RERANK:     " << std::setw(12) << std::setprecision(3) << pose.metal_rerank << "\n";
		remark << "REMARK METAL_GEO_E:      " << std::setw(12) << std::setprecision(3) << pose.metal_geo_rerank << "\n";
		remark << "REMARK METAL_WATER_E:    " << std::setw(12) << std::setprecision(3) << pose.metal_water_rerank << "\n";
		remark << "REMARK METAL_JT_E:       " << std::setw(12) << std::setprecision(3) << pose.metal_jt_rerank << "\n";
	}
	if (m_sf_choice == SF_AD42) {
		std::string ligand_metal_detail;
		fl ligand_metal_penalty = ligand_metal_geometry_penalty(m_model, &ligand_metal_detail);
		if (!ligand_metal_detail.empty()) {
			remark << "REMARK LIGAND_METAL_GEOM:" << std::setw(12) << std::setprecision(3) << ligand_metal_penalty << "\n";
			if (pose.ligand_metal_geom_rerank != 0)
				remark << "REMARK LIGAND_METAL_RERANK:" << std::setw(9) << std::setprecision(3) << pose.ligand_metal_geom_rerank << "\n";
			remark << ligand_metal_detail;
		}
	}

	if (m_sf_choice == SF_AD42 && m_ad4grid.has_reactive_state()) {
		fl dist_ang = 0, angle_deg_val = 0;
		fl dist_energy = 0, ang_energy = 0;
		m_ad4grid.get_reactive_geometry(m_model, dist_ang, angle_deg_val);
		m_ad4grid.get_reactive_terms(m_model, dist_energy, ang_energy);
		remark << "REMARK REACTIVE_DIST:    " << std::setw(12) << std::setprecision(3) << dist_ang     << " A\n";
		remark << "REMARK REACTIVE_ANGLE:   " << std::setw(12) << std::setprecision(3) << angle_deg_val << " deg\n";
		remark << "REMARK REACTIVE_DIST_E:  " << std::setw(12) << std::setprecision(3) << dist_energy  << " kcal/mol\n";
		remark << "REMARK REACTIVE_ANGLE_E: " << std::setw(12) << std::setprecision(3) << ang_energy   << " kcal/mol\n";
		remark << "REMARK REACTIVE_TOTAL:   " << std::setw(12) << std::setprecision(3) << (dist_energy+ang_energy) << " kcal/mol\n";
		// NAC: near-attack conformer (dist < 3.5 A AND angle within 30 deg of target)
		bool nac_dist = (dist_ang > 0 && dist_ang < 3.5);
		bool has_angle_constraint = !m_reactive_options.receptor_frame_atom_spec.empty();
		bool nac_angle = !has_angle_constraint ||
		    std::abs(angle_deg_val - m_reactive_options.target_angle_deg) < 30.0;
		remark << "REMARK REACTIVE_NAC:     " << (nac_dist && nac_angle ? "       YES" : "        NO") << "\n";
	}

	// E1: metal coordination geometry — nearest donor atom per metal center
	if (m_sf_choice == SF_AD42 && m_metal_mode != ag4_metal_mode::none) {
		load_metal_cache_if_needed();
		sz nat = m_model.num_movable_atoms();
		for (const auto& ma : m_metal_cache) {
			double best_d = 5.0; // report only if within 5 Å
			const char* best_donor = nullptr;
			VINA_FOR(ai, nat) {
				sz at = m_model.get_atom(ai).get(atom_type::AD);
				const char* dn = nullptr;
				if      (at == AD_TYPE_NA) dn = "NA";
				else if (at == AD_TYPE_OA) dn = "OA";
				else if (at == AD_TYPE_SA) dn = "SA";
				else if (at == AD_TYPE_N)  dn = "N";
				else if (at == AD_TYPE_O)  dn = "O";
				else if (at == AD_TYPE_S)  dn = "S";
				if (!dn) continue;
				const vec& lp = m_model.movable_coords(ai);
				double dx=lp[0]-ma.x, dy=lp[1]-ma.y, dz=lp[2]-ma.z;
				double d = std::sqrt(dx*dx+dy*dy+dz*dz);
				if (d < best_d) { best_d=d; best_donor=dn; }
			}
			if (best_donor) {
				remark << "REMARK METAL_COORD: " << std::left << std::setw(3) << ma.type
				       << "  donor=" << best_donor
				       << "  d=" << std::right << std::setprecision(2) << best_d << " A\n";
			}
		}
	}

	return remark.str();
}

std::string Vina::get_poses(int how_many, double energy_range) {
	int n = 0;
	double best_energy = 0;
	std::ostringstream out;
	std::string remarks;

	if (how_many < 0) {
		throw vina_runtime_error("number of poses written must be greater than zero.");
	}

	if (energy_range < 0) {
		throw vina_runtime_error("energy range must be greater than zero.");
	}

	if (!m_poses.empty()) {
		// Get energy from the best conf
		best_energy = m_poses[0].e;

		VINA_FOR_IN(i, m_poses) {
			/* Stop if:
				- We wrote the number of conf asked
				- If there is no conf to write
				- The energy of the current conf is superior than best_energy + energy_range
			*/
			if (n >= how_many || !not_max(m_poses[i].e) || m_poses[i].e > best_energy + energy_range)
				break; // check energy_range sanity FIXME

			// Push the current pose to model
			m_model.set(m_poses[i].c);

			// Write conf
			remarks = vina_remarks(m_poses[i], m_poses[i].lb, m_poses[i].ub);
			out << m_model.write_model(n + 1, remarks);

			n++;
		}

		// Push back the best conf in model
		m_model.set(m_poses[0].c);

	} else {
		std::cerr << "WARNING: Could not find any poses. No poses were written.\n";
	}

	return out.str();
}

void Vina::write_poses(const std::string& output_name, int how_many, double energy_range) {
	std::string out;

	if (!m_poses.empty()) {
		// Open output file
		ofile f(make_path(output_name));
		out = get_poses(how_many, energy_range);
		f << out;
	} else {
		std::cerr << "WARNING: Could not find any poses. No poses were written.\n";
	}
}

void Vina::write_pose(const std::string& output_name, const std::string& remark) {
	std::ostringstream format_remark;
	format_remark.setf(std::ios::fixed, std::ios::floatfield);
	format_remark.setf(std::ios::showpoint);

	// Add REMARK keyword to be PDB valid
	if(!remark.empty()){
		format_remark << "REMARK " << remark << " \n";
	}
	if (m_sf_choice == SF_AD42 && m_ad4grid.has_reactive_state()) {
		fl dist_ang = 0, angle_deg_val = 0;
		fl dist_energy = 0, ang_energy = 0;
		m_ad4grid.get_reactive_geometry(m_model, dist_ang, angle_deg_val);
		m_ad4grid.get_reactive_terms(m_model, dist_energy, ang_energy);
		format_remark << "REMARK REACTIVE_DIST:    " << std::setw(12) << std::setprecision(3) << dist_ang     << " A\n";
		format_remark << "REMARK REACTIVE_ANGLE:   " << std::setw(12) << std::setprecision(3) << angle_deg_val << " deg\n";
		format_remark << "REMARK REACTIVE_DIST_E:  " << std::setw(12) << std::setprecision(3) << dist_energy  << " kcal/mol\n";
		format_remark << "REMARK REACTIVE_ANGLE_E: " << std::setw(12) << std::setprecision(3) << ang_energy   << " kcal/mol\n";
		format_remark << "REMARK REACTIVE_TOTAL:   " << std::setw(12) << std::setprecision(3) << (dist_energy+ang_energy) << " kcal/mol\n";
		bool nac_dist = (dist_ang > 0 && dist_ang < 3.5);
		bool has_angle_constraint = !m_reactive_options.receptor_frame_atom_spec.empty();
		bool nac_angle = !has_angle_constraint ||
		    std::abs(angle_deg_val - m_reactive_options.target_angle_deg) < 30.0;
		format_remark << "REMARK REACTIVE_NAC:     " << (nac_dist && nac_angle ? "       YES" : "        NO") << "\n";
	}
	if (m_sf_choice == SF_AD42) {
		std::string ligand_metal_detail;
		fl ligand_metal_penalty = ligand_metal_geometry_penalty(m_model, &ligand_metal_detail);
		if (!ligand_metal_detail.empty()) {
			format_remark << "REMARK LIGAND_METAL_GEOM:" << std::setw(12) << std::setprecision(3) << ligand_metal_penalty << "\n";
			format_remark << ligand_metal_detail;
		}
	}

	ofile f(make_path(output_name));
	m_model.write_structure(f, format_remark.str());
}

void Vina::randomize(const int max_steps) {
	// Randomize ligand/flex residues conformation
	// Check the box was defined
	if (!m_ligand_initialized) {
		throw vina_runtime_error("Cannot do ligand randomization. Ligand(s) was(ere) not initialized.");
	} else if (!m_map_initialized) {
		throw vina_runtime_error("Cannot do ligand randomization. Affinity maps were not initialized.");
	}

	conf c;
	int seed = generate_seed();
	double penalty = 0;
	double best_clash_penalty = 0;
	std::stringstream sstm;
	rng generator(static_cast<rng::result_type>(seed));

	// It's okay to take the initial conf since we will randomize it
	conf init_conf = m_model.get_initial_conf();
	conf best_conf = init_conf;

	sstm << "Randomize conformation (random seed: " << std::to_string(seed) << ")";
	doing(sstm.str(), m_verbosity, 0);
	VINA_FOR(i, max_steps) {
		c = init_conf;
		if (m_sf_choice == SF_VINA || m_sf_choice == SF_VINARDO)
			c.randomize(m_grid.corner1(), m_grid.corner2(), generator);
		else
			c.randomize(m_ad4grid.corner1(), m_ad4grid.corner2(), generator);
		m_model.set(c);
		penalty = m_model.clash_penalty();

		if (i == 0 || penalty < best_clash_penalty) {
			best_conf = c;
			best_clash_penalty = penalty;
		}
	}
	done(m_verbosity, 0);

	m_model.set(best_conf);

	if (m_verbosity > 1) {
		std::cout << "Clash penalty: " << best_clash_penalty << "\n";
	}
}

void Vina::show_score(const std::vector<double> energies) {
	std::cout << "Estimated Free Energy of Binding   : " << std::fixed << std::setprecision(3) << energies[0] << " (kcal/mol) [=(1)+(2)+(3)-(4)]\n";
	std::cout << "(1) Final Intermolecular Energy    : " << std::fixed << std::setprecision(3) << energies[1] + energies[2] << " (kcal/mol)\n";
	std::cout << "    Ligand - Receptor              : " << std::fixed << std::setprecision(3) << energies[1] << " (kcal/mol)\n";
	std::cout << "    Ligand - Flex side chains      : " << std::fixed << std::setprecision(3) << energies[2] << " (kcal/mol)\n";
	std::cout << "(2) Final Total Internal Energy    : " << std::fixed << std::setprecision(3) << energies[3] + energies[4] + energies[5] << " (kcal/mol)\n";
	std::cout << "    Ligand                         : " << std::fixed << std::setprecision(3) << energies[5] << " (kcal/mol)\n";
	std::cout << "    Flex   - Receptor              : " << std::fixed << std::setprecision(3) << energies[3] << " (kcal/mol)\n";
	std::cout << "    Flex   - Flex side chains      : " << std::fixed << std::setprecision(3) << energies[4] << " (kcal/mol)\n";
	std::cout << "(3) Torsional Free Energy          : " << std::fixed << std::setprecision(3) << energies[6] << " (kcal/mol)\n";
	if (m_sf_choice == SF_VINA || m_sf_choice == SF_VINARDO) {
		std::cout << "(4) Unbound System's Energy        : " << std::fixed << std::setprecision(3) << energies[7] << " (kcal/mol)\n";
	} else {
		std::cout << "(4) Unbound System's Energy [=(2)] : " << std::fixed << std::setprecision(3) << energies[7] << " (kcal/mol)\n";
	}
	if (m_sf_choice == SF_AD42 && m_ad4grid.has_reactive_state() && (m_reactive_options.debug_energy || m_reactive_options.gradient_check)) {
		fl reactive_distance = 0;
		fl reactive_angle = 0;
		if (m_ad4grid.get_reactive_terms(m_model, reactive_distance, reactive_angle, 0)) {
			if (m_reactive_options.debug_energy) {
				std::cout << "Reactive Distance Energy           : " << std::fixed << std::setprecision(3) << reactive_distance << " (kcal/mol)\n";
				std::cout << "Reactive Angle Energy              : " << std::fixed << std::setprecision(3) << reactive_angle << " (kcal/mol)\n";
				std::cout << "Reactive Total                     : " << std::fixed << std::setprecision(3) << (reactive_distance + reactive_angle) << " (kcal/mol)\n";
			}
			if (m_reactive_options.gradient_check)
				m_ad4grid.debug_check_reactive_gradient(m_model, (fl)m_reactive_options.gradient_check_eps, std::cout);
		}
	}
}

std::vector<double> Vina::score(double intramolecular_energy) {
	// Score the current conf in the model
	double total = 0;
	double inter = 0;
	double intra = 0;
	double all_grids = 0; // ligand & flex
	double lig_grids = 0;
	double flex_grids = 0;
	double lig_intra = 0;
	double conf_independent = 0;
	double inter_pairs = 0;
	double intra_pairs = 0;
	const vec authentic_v(1000, 1000, 1000);
	std::vector<double> energies;

	if (m_sf_choice == SF_VINA || m_sf_choice == SF_VINARDO) {
		// Inter
		if (m_no_refine || !m_receptor_initialized)
			all_grids = m_grid.eval(m_model, authentic_v[1]); // [1] ligand & flex -- grid
		else
			all_grids = m_non_cache.eval(m_model, authentic_v[1]); // [1] ligand & flex -- grid
		inter_pairs = m_model.eval_inter(m_precalculated_byatom, authentic_v); // [1] ligand -- flex
		// Intra
		if (m_no_refine || !m_receptor_initialized)
			flex_grids = m_grid.eval_intra(m_model, authentic_v[1]); // [1] flex -- grid
		else
			flex_grids = m_non_cache.eval_intra(m_model, authentic_v[1]); // [1] flex -- grid
		intra_pairs = m_model.evalo(m_precalculated_byatom, authentic_v); // [1] flex_i -- flex_i and flex_i -- flex_j
		lig_grids = all_grids - flex_grids;
		inter = lig_grids + inter_pairs;
		lig_intra = m_model.evali(m_precalculated_byatom, authentic_v); // [2] ligand_i -- ligand_i
		intra = flex_grids + intra_pairs + lig_intra;
		// Total
		total = m_scoring_function->conf_independent(m_model, inter + intra - intramolecular_energy); // we pass intermolecular energy from the best pose
		// Torsion, we want to know how much torsion penalty was added to the total energy
		conf_independent = total - (inter + intra - intramolecular_energy);
	} else {
		// Inter
		lig_grids = m_ad4grid.eval(m_model, authentic_v[1]); // [1] ligand -- grid
		inter_pairs = m_model.eval_inter(m_precalculated_byatom, authentic_v); // [1] ligand -- flex
		inter = lig_grids + inter_pairs;
		// Intra
		flex_grids = m_ad4grid.eval_intra(m_model, authentic_v[1]); // [1] flex -- grid
		intra_pairs = m_model.evalo(m_precalculated_byatom, authentic_v); // [1] flex_i -- flex_i and flex_i -- flex_j
		lig_intra = m_model.evali(m_precalculated_byatom, authentic_v); // [2] ligand_i -- ligand_i
		intra = flex_grids + intra_pairs + lig_intra;
		// Torsion
		conf_independent = m_scoring_function->conf_independent(m_model, 0); // [3] we can pass e=0 because we do not modify the energy like in vina
		// Total
		total = inter + conf_independent; // (+ intra - intra)
	}

	energies.push_back(total);
	energies.push_back(lig_grids);
	energies.push_back(inter_pairs);
	energies.push_back(flex_grids);
	energies.push_back(intra_pairs);
	energies.push_back(lig_intra);
	energies.push_back(conf_independent);

	if (m_sf_choice == SF_VINA  || m_sf_choice == SF_VINARDO) {
		energies.push_back(intramolecular_energy);
	} else {
		energies.push_back(intra);
	}

	return energies;
}

std::vector<double> Vina::score() {
	// Score the current conf in the model
	// Check if ff and ligand were initialized
	// Check if the ligand is not outside the box
	if (!m_ligand_initialized) {
		throw vina_runtime_error("Cannot score the pose. Ligand(s) was(ere) not initialized.");
	} else if (!m_map_initialized) {
		throw vina_runtime_error("Cannot score the pose. Affinity maps were not initialized.");
	} else if ((m_sf_choice == SF_VINA || m_sf_choice == SF_VINARDO) ? !m_grid.is_in_grid(m_model) : !m_ad4grid.is_in_grid(m_model)) {
		throw vina_runtime_error("The ligand is outside the grid box. Increase the size of the grid box or center it accordingly around the ligand.");
	}

	double intramolecular_energy = 0;
	const vec authentic_v(1000, 1000, 1000);

	if(m_sf_choice == SF_VINA || m_sf_choice == SF_VINARDO) {
		intramolecular_energy = m_model.eval_intramolecular(m_precalculated_byatom, m_grid, authentic_v);
	}

	std::vector<double> energies = score(intramolecular_energy);
	return energies;
}

std::vector<double> Vina::optimize(output_type& out, int max_steps) {
	// Local optimization of the ligand conf
	change g(m_model.get_size());
	quasi_newton quasi_newton_par;
	const fl slope = 1e6;
	const vec authentic_v(1000, 1000, 1000);
	std::vector<double> energies_before_opt;
	std::vector<double> energies_after_opt;
	int evalcount = 0;

	// Define the number minimization steps based on the number moving atoms
	if (max_steps == 0) {
		max_steps = unsigned((25 + m_model.num_movable_atoms()) / 3);
		if (m_verbosity > 1)
			std::cout << "Number of local optimization steps: " << max_steps << "\n";
	}
	quasi_newton_par.max_steps = max_steps;

	if (m_verbosity > 1) {
		std::cout << "Before local optimization:\n";
		energies_before_opt = score();
		show_score(energies_before_opt);
	}

	doing("Performing local search", m_verbosity, 0);
	// Try 5 five times to optimize locally the conformation
	VINA_FOR(p, 5) {
		if (m_sf_choice == SF_VINA || m_sf_choice == SF_VINARDO) {
			quasi_newton_par(m_model, m_precalculated_byatom, m_grid,    out, g, authentic_v, evalcount);
			// Break if we succeed to bring (back) the ligand within the grid
			if (m_grid.is_in_grid(m_model))
				break;
		} else {
			quasi_newton_par(m_model, m_precalculated_byatom, m_ad4grid, out, g, authentic_v, evalcount);
			if (m_ad4grid.is_in_grid(m_model))
				break;
		}
	}
	done(m_verbosity, 0);

	energies_after_opt = score();

	return energies_after_opt;
}

std::vector<double> Vina::optimize(int max_steps) {
	// Local optimization of the ligand conf
	// Check if ff, box and ligand were initialized
	// Check if the ligand is not outside the box
	if (!m_ligand_initialized) {
		throw vina_runtime_error("Cannot do the optimization. Ligand(s) was(ere) not initialized.");
	} else if (!m_map_initialized) {
		throw vina_runtime_error("Cannot do the optimization. Affinity maps were not initialized.");
	} else if ((m_sf_choice == SF_VINA || m_sf_choice == SF_VINARDO) ? !m_grid.is_in_grid(m_model) : !m_ad4grid.is_in_grid(m_model)) {
		throw vina_runtime_error("The ligand is outside the grid box. Increase the size of the grid box or center it accordingly around the ligand.");
	}

	double e = 0;
	conf c;

	if (!m_poses.empty()) {
		// if m_poses is not empty, it means that we did a docking before
		// But it is really that useful to minimize after docking?
		e = m_poses[0].e;
		c = m_poses[0].c;
	} else {
		c = m_model.get_initial_conf();
	}

	output_type out(c, e);

	std::vector<double> energies = optimize(out, max_steps);

	return energies;
}

output_container Vina::remove_redundant(const output_container &in, fl min_rmsd) {
	output_container tmp;
	VINA_FOR_IN(i, in)
	add_to_output_container(tmp, in[i], min_rmsd, in.size());
	return tmp;
}

void Vina::global_search(const int exhaustiveness, const int n_poses, const double min_rmsd, const int max_evals) {
	// Vina search (Monte-carlo and local optimization)
	// Check if ff, box and ligand were initialized
	if (!m_ligand_initialized) {
		throw vina_runtime_error("Cannot do the global search. Ligand(s) was(ere) not initialized.");
	} else if (!m_map_initialized) {
		throw vina_runtime_error("Cannot do the global search. Affinity maps were not initialized.");
	} else if (exhaustiveness < 1) {
		throw vina_runtime_error("Exhaustiveness must be 1 or greater");
	}

	if (exhaustiveness < m_cpu) {
		std::cerr << "WARNING: At low exhaustiveness, it may be impossible to utilize all CPUs.\n";
	}

	double e = 0;
	double intramolecular_energy = 0;
	const vec authentic_v(1000, 1000, 1000);
	model best_model;
	boost::optional<model> ref;
	output_container poses;
	std::stringstream sstm;
	rng generator(static_cast<rng::result_type>(m_seed));

	// ── C3 Two-step: save & suppress reactive state during MC (phase 1) ─────
	const bool two_step_active = (m_sf_choice == SF_AD42)
	                           && m_ad4grid.has_reactive_state()
	                           && m_reactive_enabled
	                           && m_reactive_options.two_step;
	reactive_state saved_rs;
	if (two_step_active) {
		saved_rs = m_ad4grid.get_reactive_state();
		if (m_reactive_options.weak_attractor) {
			// C3b: broad/weak Gaussian bias toward anchor (Goullieux 2023 "attracting cavities")
			reactive_state weak_rs         = saved_rs;
			weak_rs.attractor_width       *= 3.0;   // σ×3: very broad spatial basin
			weak_rs.attractor_strength    *= 0.15;  // ε×0.15: ~1-2 kcal/mol gentle pull
			weak_rs.angle_strength         = 0.0;   // no angular constraint in phase 1
			m_ad4grid.set_reactive_state(weak_rs);
			if (m_verbosity > 0)
				std::cout << "Two-step: phase 1 — MC presample with weak attractor "
				          << "(σ=" << weak_rs.attractor_width << " Å, "
				          << "ε=" << weak_rs.attractor_strength << " kcal/mol)\n";
		} else {
			m_ad4grid.clear_reactive_state();
			if (m_verbosity > 0)
				std::cout << "Two-step: phase 1 — MC presample without reactive constraints\n";
		}
	}

	// Setup Monte-Carlo search
	parallel_mc parallelmc;
	sz heuristic = m_model.num_movable_atoms() + 10 * m_model.get_size().num_degrees_of_freedom();
	parallelmc.mc.global_steps = unsigned(70 * 3 * (50 + heuristic) / 2); // 2 * 70 -> 8 * 20 // FIXME
	parallelmc.mc.local_steps = unsigned((25 + m_model.num_movable_atoms()) / 3);
	parallelmc.mc.max_evals = max_evals;
	parallelmc.mc.min_rmsd = min_rmsd;
	parallelmc.mc.num_saved_mins = n_poses;
	parallelmc.mc.hunt_cap = vec(10, 10, 10);
	parallelmc.num_tasks = exhaustiveness;
	parallelmc.num_threads = m_cpu;
	parallelmc.display_progress = (m_verbosity > 0);

	// Docking search
	sstm << "Performing docking (random seed: " << m_seed << ")";
	doing(sstm.str(), m_verbosity, 0);
	if (m_sf_choice == SF_VINA || m_sf_choice == SF_VINARDO) {
		parallelmc(m_model, poses, m_precalculated_byatom,    m_grid, m_grid.corner1(), m_grid.corner2(), generator, m_progress_callback);
	} else {
		parallelmc(m_model, poses, m_precalculated_byatom, m_ad4grid, m_ad4grid.corner1(), m_ad4grid.corner2(), generator, m_progress_callback);
	}
	done(m_verbosity, 1);

	// ── C3 Two-step: phase 2 — restore reactive, filter, refine ─────────────
	if (two_step_active) {
		m_ad4grid.set_reactive_state(saved_rs);
		if (m_verbosity > 0)
			std::cout << "Two-step: phase 2 — filtering poses within "
			          << m_reactive_options.presample_max_dist << " Å of reactive anchor\n";

		const double max_dist   = m_reactive_options.presample_max_dist;
		const sz     lig_idx    = saved_rs.ligand_atom_index;
		const vec&   anchor_xyz = saved_rs.receptor_atom_xyz;

		// Filter: keep poses where reactive ligand atom is within max_dist of anchor
		output_container cands;
		VINA_FOR_IN(i, poses) {
			m_model.set(poses[i].c);
			const vec lp = m_model.get_coords(lig_idx);
			fl dx = lp[0] - anchor_xyz[0];
			fl dy = lp[1] - anchor_xyz[1];
			fl dz = lp[2] - anchor_xyz[2];
			if (std::sqrt(dx*dx + dy*dy + dz*dz) <= max_dist)
				add_to_output_container(cands, poses[i], 0.0, poses.size());
		}

		if (cands.empty()) {
			m_ad4grid.set_reactive_state(saved_rs);
			throw vina_runtime_error("Two-step reactive docking found no phase-1 poses within the reactive presample distance. Increase --reactive_presample_dist, enable --reactive_weak_attractor, or increase exhaustiveness.");
		} else if (m_verbosity > 0) {
			std::cout << "Two-step: " << cands.size()
			          << " pose(s) passed distance filter — running quasi-newton refinement\n";
		}

		change g(m_model.get_size());
		quasi_newton qn;
		qn.max_steps = unsigned((25 + m_model.num_movable_atoms()) / 3);
		int evalcnt = 0;
		VINA_FOR_IN(i, cands) {
			m_model.set(cands[i].c);
			qn(m_model, m_precalculated_byatom, m_ad4grid, cands[i], g, authentic_v, evalcnt);
			cands[i].coords = m_model.get_heavy_atom_movable_coords();
		}

		poses = remove_redundant(cands, min_rmsd);
		if (m_verbosity > 0)
			std::cout << "Two-step: " << poses.size() << " pose(s) after deduplication\n";
	} else {
		// Docking post-processing and rescoring
		poses = remove_redundant(poses, min_rmsd);
	}
	if (!poses.empty()) {
		// For the Vina scoring function, we take the intramolecular energy from the best pose
		// the order must not change because of non-decreasing g (see paper), but we'll re-sort in case g is non strictly increasing
		if (m_sf_choice == SF_VINA || m_sf_choice == SF_VINARDO) {
			// Refine poses if no_refine is false and got receptor
			if (!m_no_refine & m_receptor_initialized) {
				change g(m_model.get_size());
				quasi_newton quasi_newton_par;
				//std::vector<double> energies_before_opt;
				//std::vector<double> energies_after_opt;
				int evalcount = 0;
				const fl slope = 1e6;
				m_non_cache.slope = slope;
				quasi_newton_par.max_steps = unsigned((25 + m_model.num_movable_atoms()) / 3);

				VINA_FOR_IN(i, poses){
					VINA_FOR(p, 5){
						m_non_cache.slope = 100 * std::pow(10.0, 2.0*p);
						quasi_newton_par(m_model, m_precalculated_byatom, m_non_cache, poses[i], g, authentic_v, evalcount);
						if(m_non_cache.within(m_model))
							break;
					}
					poses[i].coords = m_model.get_heavy_atom_movable_coords();
					m_non_cache.slope = slope;
					// rescoring in case a ligand or flex sidechain atom is outside box
					// ensuring poses will be sorted with the same slope (a.k.a. out of
					// box penalty) that will be used to calculate final energies.
					m_model.set(poses[i].c);
					double all_grids = m_non_cache.eval(m_model, authentic_v[1]);
					double inter_pairs = m_model.eval_inter(m_precalculated_byatom, authentic_v); // ligand -- flex
					double intra_pairs = m_model.evalo(m_precalculated_byatom, authentic_v); // flex_i -- flex_i and flex_i -- flex_j
					double lig_intra = m_model.evali(m_precalculated_byatom, authentic_v); // ligand_i -- ligand_i
					poses[i].e = all_grids + inter_pairs + intra_pairs + lig_intra;
				}
			}

			poses.sort(); // order often changes after non_cache refinement
			m_model.set(poses[0].c);
			if (m_no_refine || !m_receptor_initialized)
				intramolecular_energy = m_model.eval_intramolecular(m_precalculated_byatom, m_grid, authentic_v);
			else
				intramolecular_energy = m_model.eval_intramolecular(m_precalculated_byatom, m_non_cache, authentic_v);
		}

		VINA_FOR_IN(i, poses) {
			if (m_verbosity > 1)
				std::cout << "ENERGY FROM SEARCH: " << poses[i].e << "\n";

			m_model.set(poses[i].c);

			// For AD42 intramolecular_energy is equal to 0
			std::vector<double> energies = score(intramolecular_energy);
			// Store energy components in current pose
			poses[i].e = energies[0]; // specific to each scoring function
			poses[i].inter = energies[1] + energies[2];
			poses[i].intra = energies[3] + energies[4] + energies[5];
			poses[i].total = poses[i].inter + poses[i].intra; // cost function for optimization
			poses[i].conf_independent = energies[6]; // "torsion"
			poses[i].unbound = energies[7]; // specific to each scoring function
			poses[i].metal_rerank = 0;
			poses[i].metal_geo_rerank = 0;
			poses[i].metal_water_rerank = 0;
			poses[i].metal_jt_rerank = 0;
			poses[i].ligand_metal_geom_penalty = 0;
			poses[i].ligand_metal_geom_rerank = 0;
			if (m_sf_choice == SF_AD42 && m_ad4grid.has_metal_state()) {
				fl geo_term = 0, water_term = 0, jt_term = 0;
				fl rerank_delta = m_ad4grid.get_metal_rerank_terms(m_model, geo_term, water_term, jt_term);
				poses[i].metal_rerank = rerank_delta;
				poses[i].metal_geo_rerank = geo_term;
				poses[i].metal_water_rerank = water_term;
				poses[i].metal_jt_rerank = jt_term;
				poses[i].e += rerank_delta;
				poses[i].inter += rerank_delta;
				poses[i].total += rerank_delta;
			}
			if (m_sf_choice == SF_AD42) {
				fl ligand_metal_penalty = ligand_metal_geometry_penalty(m_model);
				poses[i].ligand_metal_geom_penalty = ligand_metal_penalty;
				fl ligand_metal_delta = ligand_metal_penalty * m_ligand_metal_geometry_weight;
				poses[i].ligand_metal_geom_rerank = ligand_metal_delta;
				if (ligand_metal_delta != 0) {
					poses[i].e += ligand_metal_delta;
					poses[i].intra += ligand_metal_delta;
					poses[i].total += ligand_metal_delta;
				}
			}
			if (!m_no_refine && m_receptor_initialized)
				poses[i].coords = m_model.get_heavy_atom_movable_coords();

			if (m_verbosity > 1) {
				std::cout << "FINAL ENERGY: \n";
				show_score(energies);
			}
		}

		// In AD4, the unbound energy is intra for each pose, so order may have changed
		// The order does not change in Vina because unbound is intra of the 1st pose
		if (m_sf_choice == SF_AD42)
			poses.sort();

		// Now compute RMSD from the best model
		// Necessary to do it in two pass for AD4 scoring function
		m_model.set(poses[0].c);
		best_model = m_model;

		if (m_verbosity > 0) {
			std::cout << '\n';
			std::cout << "mode |   affinity | dist from best mode\n";
			std::cout << "     | (kcal/mol) | rmsd l.b.| rmsd u.b.\n";
			std::cout << "-----+------------+----------+----------\n";
		}

		VINA_FOR_IN(i, poses) {
			m_model.set(poses[i].c);

			// Get RMSD between current pose and best_model
			const model &r = ref ? ref.get() : best_model;
			poses[i].lb = m_model.rmsd_lower_bound(r);
			poses[i].ub = m_model.rmsd_upper_bound(r);

			if (m_verbosity > 0) {
				std::cout << std::setw(4) << i + 1 << "    " << std::setw(9) << std::setprecision(4) << poses[i].e;
				std::cout << "  " << std::setw(9) << std::setprecision(4) << poses[i].lb;
				std::cout << "  " << std::setw(9) << std::setprecision(4) << poses[i].ub << "\n";
			}
		}

		// Clean up by putting back the best pose in model
		m_model.set(poses[0].c);
	} else {
		std::cerr << "WARNING: Zero poses in output container after global search. This should not be happening and is likely a bug.\n";
		std::cerr << "WARNING: If possible, please file a bug report with your input files and random seed on GitHub.\n";
	}

	// Store results in Vina object
	m_poses = poses;
}

Vina::~Vina() {
	//OpenBabel::OBMol m_mol;
	// model and poses
	model m_receptor;
	model m_model;
	output_container m_poses;
	bool m_receptor_initialized;
	bool m_ligand_initialized;
	// scoring function
	scoring_function_choice m_sf_choice;
	flv m_weights;
	precalculate_byatom m_precalculated_byatom;
	precalculate m_precalculated_sf;
	// maps
	cache m_grid;
	ad4cache m_ad4grid;
	non_cache m_non_cache;
	bool m_map_initialized;
	// global search
	int m_cpu;
	int m_seed;
	// others
	int m_verbosity;
}

// ── Reactive covalent docking: Vina method implementations ───────────────────

void Vina::set_reactive_options(const ReactiveOptions& opts) {
	if (!opts.enabled()) {
		m_reactive_enabled = false;
		m_reactive_options.clear();
		m_reactive_payload.clear();
		m_reactive_payload_ready = false;
		m_ad4grid.clear_reactive_state();
		return;
	}
	if (m_sf_choice != SF_AD42)
		throw vina_runtime_error("set_reactive_options: reactive covalent docking is only supported with AD4/LKDock scoring.");
	if (opts.receptor_atom_spec.empty())
		throw vina_runtime_error("set_reactive_options: receptor_atom_spec is empty.");
	if (opts.ligand_atom_spec.empty())
		throw vina_runtime_error("set_reactive_options: ligand_atom_spec is empty.");
	if (opts.bond_length < 0.0)
		throw vina_runtime_error("set_reactive_options: bond_length must be >= 0.");
	if (opts.mode == reactive_mode::hybrid && opts.bond_length <= 0.0)
		throw vina_runtime_error("set_reactive_options: hybrid covalent mode requires bond_length > 0.");
	if (opts.attractor_width <= 0.0)
		throw vina_runtime_error("set_reactive_options: attractor_width must be > 0.");
	if (opts.attractor_strength < 0.0)
		throw vina_runtime_error("set_reactive_options: attractor_strength must be >= 0.");
	if (opts.hybrid_vdw_scale < 0.0 || opts.hybrid_vdw_scale > 1.0)
		throw vina_runtime_error("set_reactive_options: hybrid_vdw_scale must be between 0 and 1.");
	if (opts.target_angle_deg < 0.0 || opts.target_angle_deg > 180.0)
		throw vina_runtime_error("set_reactive_options: target_angle must be between 0 and 180 degrees.");
	if (opts.angle_width_deg < 0.0 || opts.angle_width_deg > 180.0)
		throw vina_runtime_error("set_reactive_options: angle_width must be between 0 and 180 degrees.");
	if (opts.gradient_check && opts.gradient_check_eps <= 0.0)
		throw vina_runtime_error("set_reactive_options: gradient_check_eps must be > 0.");
	if (opts.two_step && opts.presample_max_dist <= 0.0)
		throw vina_runtime_error("set_reactive_options: reactive_presample_dist must be > 0 when two_step is enabled.");
	if (opts.geometry_mode == reactive_geometry_mode::angle
	        && opts.receptor_frame_atom_spec.empty())
		std::cerr << "WARNING: reactive geometry_mode=angle is set but receptor_frame_atom "
		             "(--reactive_frame_atom) is not specified — the angle term will not activate.\n";
	if (!opts.ligand_frame_atom_spec.empty() && opts.receptor_frame_atom_spec.empty())
		std::cerr << "WARNING: reactive ligand_frame_atom (--reactive_lig_frame_atom) is set "
		             "but receptor_frame_atom is not — ligand-frame has no effect without receptor frame.\n";
	if (opts.mode == reactive_mode::distance && opts.hybrid_vdw_scale != 0.0)
		std::cerr << "WARNING: reactive hybrid_vdw_scale=" << opts.hybrid_vdw_scale
		          << " is ignored in distance mode (only used in hybrid mode).\n";
	if (opts.weak_attractor && !opts.two_step)
		std::cerr << "WARNING: reactive_weak_attractor has no effect when reactive_two_step is not enabled.\n";
	m_reactive_options       = opts;
	m_reactive_enabled       = true;
	m_reactive_payload.clear();
	m_reactive_payload_ready = false;
}

void Vina::clear_reactive_options() {
	m_reactive_options.clear();
	m_reactive_enabled       = false;
	m_reactive_payload.clear();
	m_reactive_payload_ready = false;
	m_ad4grid.clear_reactive_state();
}

sz Vina::resolve_ligand_reactive_atom(const std::string& spec) const {
	if (!m_ligand_initialized || m_model.num_ligands() == 0) return max_sz;
	if (spec.empty()) return max_sz;

	ligand lig = m_model.get_ligand(0);
	std::string key = spec;
	std::string value = spec;
	size_t colon = spec.find(':');
	if (colon != std::string::npos) {
		key = spec.substr(0, colon);
		value = spec.substr(colon + 1);
	}

	bool is_numeric = !value.empty();
	for (char c : value) {
		if (!std::isdigit((unsigned char)c)) { is_numeric = false; break; }
	}
	if (is_numeric && (key == spec || key == "index")) {
		int local_index = std::stoi(value);
		if (local_index < 1) return max_sz;
		sz idx = lig.begin + (sz)(local_index - 1);
		if (idx < lig.end && idx < m_model.num_atoms()) return idx;
		return max_sz;
	}

	const context& cont = lig.cont;
	if (is_numeric && key == "serial") {
		int serial = std::stoi(value);
		for (const auto& pl : cont) {
			if (!pl.second) continue;
			const std::string& line = pl.first;
			if (line.size() < 11) continue;
			try {
				if (std::stoi(line.substr(6, 5)) == serial) return *pl.second;
			} catch (...) {}
		}
		return max_sz;
	}

	std::string target_name = (key == "name" || key == "atom") ? value : spec;
	for (const auto& pl : cont) {
		if (!pl.second) continue;          // not an atom line
		sz idx = *pl.second;
		const std::string& line = pl.first;
		if (line.size() < 16) continue;
		std::string aname = line.substr(12, 4);
		size_t p = aname.find_first_not_of(" ");
		if (p == std::string::npos) continue;
		aname = aname.substr(p);
		size_t q = aname.find_last_not_of(" ");
		if (q != std::string::npos) aname = aname.substr(0, q + 1);
		if (aname == target_name) return idx;
	}
	return max_sz;
}

void Vina::finalize_reactive_state_if_possible() {
	if (!m_reactive_enabled || !m_reactive_payload_ready || !m_ligand_initialized) return;
	if (m_sf_choice != SF_AD42) {
		m_ad4grid.clear_reactive_state();
		throw vina_runtime_error("reactive covalent docking is only supported with AD4/LKDock scoring.");
	}

	sz lig_idx = resolve_ligand_reactive_atom(m_reactive_options.ligand_atom_spec);
	if (lig_idx == max_sz) {
		m_ad4grid.clear_reactive_state();
		throw vina_runtime_error("reactive: ligand atom spec '" + m_reactive_options.ligand_atom_spec + "' not found.");
	}

	// Compute cos(target_angle_deg) portably
	const double pi_val = std::acos(-1.0);

	reactive_state rs;
	rs.enabled               = true;
	rs.mode                  = m_reactive_payload.mode;
	rs.geometry_mode         = m_reactive_payload.geometry_mode;
	rs.ligand_atom_index     = lig_idx;
	rs.has_ligand_atom       = true;
	rs.receptor_atom_xyz     = m_reactive_payload.receptor_atom_xyz;
	rs.has_receptor_atom     = m_reactive_payload.has_receptor_atom;
	rs.receptor_frame_xyz    = m_reactive_payload.receptor_frame_xyz;
	rs.has_receptor_frame    = m_reactive_payload.has_receptor_frame;
	rs.bond_length           = m_reactive_payload.bond_length;
	rs.attractor_width       = m_reactive_payload.attractor_width;
	rs.attractor_strength    = m_reactive_payload.attractor_strength;
	rs.angle_strength        = m_reactive_payload.angle_strength;
	rs.hybrid_vdw_scale      = m_reactive_payload.hybrid_vdw_scale;
	rs.cos_target_angle = std::cos(m_reactive_payload.target_angle_deg * pi_val / 180.0);
	{
	    double tgt = m_reactive_payload.target_angle_deg;
	    double wid = m_reactive_payload.angle_width_deg;
	    double lo  = std::min(180.0, tgt + wid);   // upper angle boundary (clamped)
	    double hi  = std::max(0.0,   tgt - wid);   // lower angle boundary (clamped)
	    // cos is decreasing: larger angle → smaller cosine
	    rs.cos_angle_lo = std::cos(lo * pi_val / 180.0);
	    rs.cos_angle_hi = std::cos(hi * pi_val / 180.0);
	}
	rs.debug                 = m_reactive_payload.debug;
	rs.debug_energy          = m_reactive_payload.debug_energy;
	rs.gradient_check        = m_reactive_payload.gradient_check;
	rs.gradient_check_eps    = m_reactive_payload.gradient_check_eps;
	// Ligand-side frame atom — enables approach-angle energy term (angle at lig_atom)
	if (!m_reactive_options.ligand_frame_atom_spec.empty()) {
		sz lf_idx = resolve_ligand_reactive_atom(m_reactive_options.ligand_frame_atom_spec);
		if (lf_idx != max_sz) {
			rs.ligand_frame_atom_index = lf_idx;
			rs.has_ligand_frame        = true;
		} else {
			m_ad4grid.clear_reactive_state();
			throw vina_runtime_error("reactive: ligand frame atom spec '" + m_reactive_options.ligand_frame_atom_spec + "' not found.");
		}
	}
	m_ad4grid.set_reactive_state(rs);

	if (m_reactive_options.debug) {
		std::cerr << "[reactive] state finalized: ligand_atom_index=" << lig_idx
		          << "  receptor_anchor=("
		          << m_reactive_payload.receptor_atom_xyz[0] << ","
		          << m_reactive_payload.receptor_atom_xyz[1] << ","
		          << m_reactive_payload.receptor_atom_xyz[2] << ")\n";
	}
}
// ─────────────────────────────────────────────────────────────────────────────
