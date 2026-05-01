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

#ifndef VINA_AD4CACHE_H
#define VINA_AD4CACHE_H

#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <boost/serialization/split_member.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/static_assert.hpp>
#include "igrid.h"
#include "grid.h"
#include "model.h"
#include "file.h"
#include "szv_grid.h"
#include "common.h"
#include "reactive_types.h"
#include "ag4_engine.h"

struct ad4cache : public igrid {
public:
    ad4cache(fl slope=1e6): m_slope(slope), m_grids(AD_TYPE_SIZE + 2) {}
	fl eval      (const model& m, fl v) const; // needs m.coords // clean up
	fl eval_intra(      model& m, fl v) const; // needs m.coords, sets m.minus_forces // clean up
	fl eval_deriv(      model& m, fl v) const; // needs m.coords, sets m.minus_forces // clean up
    grid_dims gd() const { return m_gd; }
    vec corner1() const { vec corner(m_gd[0].begin, m_gd[1].begin, m_gd[2].begin); return corner; }
    vec corner2() const { vec corner(m_gd[0].end, m_gd[1].end, m_gd[2].end); return corner; }
    bool is_in_grid(const model &m, fl margin=0.0001) const;
    bool is_atom_type_grid_initialized(sz t) const { return m_grids[t].initialized(); }
    bool are_atom_types_grid_initialized(szv atom_types) const;
    void read(const std::string& str);
    void populate_from_data(const grid_dims& gd,
                            const std::vector<sz>& ad_types,
                            const std::vector<std::vector<double>>& aff_maps,
                            const std::vector<double>& elec_map,
                            const std::vector<double>& desolv_map);
    void write(const std::string& out_prefix, const szv& atom_types, const std::string& gpf_filename="NULL",
               const std::string& fld_filename="NULL", const std::string& receptor_filename="NULL");
    void set_reactive_state(const reactive_state& rs) { m_reactive = rs; }
    void clear_reactive_state()                        { m_reactive.clear(); }
    bool has_reactive_state() const                    { return m_reactive.ready(); }
    const reactive_state& get_reactive_state() const   { return m_reactive; }
    void set_metal_state(const ag4_metal_state& ms)    { m_metal_state = ms; }
    void clear_metal_state()                           { m_metal_state.clear(); }
    bool has_metal_state() const                       { return m_metal_state.enabled(); }
    const ag4_metal_state& get_metal_state() const     { return m_metal_state; }
    fl get_metal_rerank_terms(const model& m, fl& geometry_term, fl& water_term, fl& jt_term) const;
    bool get_reactive_terms(const model& m, fl& distance_energy, fl& angle_energy,
                             vec* ligand_minus_force = 0, vec* frame_minus_force = 0) const;
    bool get_reactive_geometry(const model& m, fl& dist_ang, fl& angle_deg) const;
    bool debug_check_reactive_gradient(const model& m, fl eps, std::ostream& out) const;
    void set_metal_soft_weight(fl w) { m_metal_soft_weight = w; }
    fl   get_metal_soft_weight() const { return m_metal_soft_weight; }
private:
    fl eval_metal_soft_grad(model& m, fl T = 5.0) const;
	grid_dims m_gd;
    reactive_state m_reactive;
    ag4_metal_state m_metal_state;
    fl m_metal_soft_weight = 0.0; // 0 = disabled; no effect on non-metal receptors
	fl m_slope; // does not get (de-)serialized
	std::vector<grid> m_grids;
};

#endif