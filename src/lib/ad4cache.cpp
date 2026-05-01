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

#include <cmath>
#include "ad4cache.h"


namespace fs = boost::filesystem;

bool ad4cache::get_reactive_terms(const model& m, fl& distance_energy, fl& angle_energy,
                                    vec* ligand_minus_force, vec* frame_minus_force) const {
	if (frame_minus_force) frame_minus_force->assign(0);
	distance_energy = 0;
	angle_energy = 0;
	if (ligand_minus_force) ligand_minus_force->assign(0);
	if (!m_reactive.ready()) return false;

	const vec& lxyz = m.coords[m_reactive.ligand_atom_index];
	const vec& rxyz = m_reactive.receptor_atom_xyz;
	fl dx = lxyz[0] - rxyz[0];
	fl dy = lxyz[1] - rxyz[1];
	fl dz = lxyz[2] - rxyz[2];
	fl r  = std::sqrt(dx*dx + dy*dy + dz*dz);
	if (r < (fl)1e-8) r = (fl)1e-8;
	fl dr    = r - (fl)m_reactive.bond_length;
	fl w     = (fl)m_reactive.attractor_width;
	fl A     = (fl)m_reactive.attractor_strength;
	fl gauss = std::exp(-(dr*dr) / (2.0 * w*w));
	distance_energy = -A * gauss;

	if (ligand_minus_force) {
		fl dEdr  = A * dr / (w*w) * gauss;
		fl inv_r = (fl)1.0 / r;
		(*ligand_minus_force)[0] += -dEdr * dx * inv_r;
		(*ligand_minus_force)[1] += -dEdr * dy * inv_r;
		(*ligand_minus_force)[2] += -dEdr * dz * inv_r;
	}

	// ── Receptor-frame angle term (angle at rec_atom between lig_atom and rec_frame) ──
	if (m_reactive.has_receptor_frame &&
	    m_reactive.geometry_mode == reactive_geometry_mode::angle) {
		const vec& fxyz = m_reactive.receptor_frame_xyz;
		fl ax = dx, ay = dy, az = dz;
		fl bx = fxyz[0] - rxyz[0], by = fxyz[1] - rxyz[1], bz = fxyz[2] - rxyz[2];
		fl a_len = r;
		fl b_len = std::sqrt(bx*bx + by*by + bz*bz);
		if (b_len < (fl)1e-8) b_len = (fl)1e-8;
		fl cos_th = (ax*bx + ay*by + az*bz) / (a_len*b_len);
		if (cos_th > (fl)1) cos_th = (fl)1;
		if (cos_th < (fl)-1) cos_th = (fl)-1;
		fl c_lo = (fl)m_reactive.cos_angle_lo;  // cos(target+width): lower cosine = upper angle
		fl c_hi = (fl)m_reactive.cos_angle_hi;  // cos(target-width): upper cosine = lower angle
		fl K    = (fl)m_reactive.angle_strength;
		// Flat-bottom in exact angle space: zero penalty within [target-width, target+width]
		fl pen;
		if      (cos_th < c_lo) pen = cos_th - c_lo;  // angle too large
		else if (cos_th > c_hi) pen = cos_th - c_hi;  // angle too small
		else                    pen = (fl)0;
		angle_energy += K * pen * pen;

		if (ligand_minus_force) {
			fl dE = 2.0 * K * pen;
			fl inv_a = (fl)1.0 / a_len;
			fl inv_b = (fl)1.0 / b_len;
			fl a_hat_x = ax * inv_a, a_hat_y = ay * inv_a, a_hat_z = az * inv_a;
			fl b_hat_x = bx * inv_b, b_hat_y = by * inv_b, b_hat_z = bz * inv_b;
			(*ligand_minus_force)[0] += -dE * (b_hat_x - a_hat_x * cos_th) * inv_a;
			(*ligand_minus_force)[1] += -dE * (b_hat_y - a_hat_y * cos_th) * inv_a;
			(*ligand_minus_force)[2] += -dE * (b_hat_z - a_hat_z * cos_th) * inv_a;
		}
	}

	// ── Ligand-frame approach-angle term (angle at lig_atom between rec_atom and lig_frame) ──
	if (m_reactive.has_ligand_frame &&
	    m_reactive.ligand_frame_atom_index < m.coords.size()) {
		const vec& lfxyz = m.coords[m_reactive.ligand_frame_atom_index];
		// u = rxyz - lxyz  (from lig_atom toward receptor anchor)
		fl ux = -dx, uy = -dy, uz = -dz;  // same as rxyz - lxyz
		// v = lfxyz - lxyz (from lig_atom toward lig_frame)
		fl vx = lfxyz[0] - lxyz[0], vy = lfxyz[1] - lxyz[1], vz = lfxyz[2] - lxyz[2];
		fl u_len = r;
		fl v_len = std::sqrt(vx*vx + vy*vy + vz*vz);
		if (v_len < (fl)1e-8) v_len = (fl)1e-8;
		fl cos_phi = (ux*vx + uy*vy + uz*vz) / (u_len*v_len);
		if (cos_phi > (fl)1)  cos_phi = (fl)1;
		if (cos_phi < (fl)-1) cos_phi = (fl)-1;
		fl c_lo_l = (fl)m_reactive.cos_angle_lo;
		fl c_hi_l = (fl)m_reactive.cos_angle_hi;
		fl K_l    = (fl)m_reactive.angle_strength;
		fl pen_l;
		if      (cos_phi < c_lo_l) pen_l = cos_phi - c_lo_l;
		else if (cos_phi > c_hi_l) pen_l = cos_phi - c_hi_l;
		else                       pen_l = (fl)0;
		angle_energy += K_l * pen_l * pen_l;

		if (pen_l != (fl)0 && (ligand_minus_force || frame_minus_force)) {
			fl dE_l  = 2.0f * K_l * pen_l;
			fl inv_u = (fl)1.0 / u_len;
			fl inv_v = (fl)1.0 / v_len;
			fl u_hat_x = ux*inv_u, u_hat_y = uy*inv_u, u_hat_z = uz*inv_u;
			fl v_hat_x = vx*inv_v, v_hat_y = vy*inv_v, v_hat_z = vz*inv_v;
			if (ligand_minus_force) {
				// Vertex L force: d(cos)/dL = -[(v_hat-u_hat*cos)/|u| + (u_hat-v_hat*cos)/|v|]
				(*ligand_minus_force)[0] += dE_l * ((v_hat_x - u_hat_x*cos_phi)*inv_u + (u_hat_x - v_hat_x*cos_phi)*inv_v);
				(*ligand_minus_force)[1] += dE_l * ((v_hat_y - u_hat_y*cos_phi)*inv_u + (u_hat_y - v_hat_y*cos_phi)*inv_v);
				(*ligand_minus_force)[2] += dE_l * ((v_hat_z - u_hat_z*cos_phi)*inv_u + (u_hat_z - v_hat_z*cos_phi)*inv_v);
			}
			if (frame_minus_force) {
				// Frame atom F force: d(cos)/dF = (u_hat - cos*v_hat)/|v|  → minus_force = -dE * d(cos)/dF
				(*frame_minus_force)[0] += -dE_l * (u_hat_x - v_hat_x*cos_phi) * inv_v;
				(*frame_minus_force)[1] += -dE_l * (u_hat_y - v_hat_y*cos_phi) * inv_v;
				(*frame_minus_force)[2] += -dE_l * (u_hat_z - v_hat_z*cos_phi) * inv_v;
			}
		}
	}

	return true;
}

bool ad4cache::get_reactive_geometry(const model& m, fl& dist_ang, fl& angle_deg) const {
	dist_ang  = 0;
	angle_deg = 0;
	if (!m_reactive.ready()) return false;

	const vec& lxyz = m.coords[m_reactive.ligand_atom_index];
	const vec& rxyz = m_reactive.receptor_atom_xyz;
	fl dx = lxyz[0]-rxyz[0], dy = lxyz[1]-rxyz[1], dz = lxyz[2]-rxyz[2];
	dist_ang = std::sqrt(dx*dx + dy*dy + dz*dz);

	if (m_reactive.has_receptor_frame &&
	    m_reactive.geometry_mode == reactive_geometry_mode::angle) {
		fl r = (dist_ang < (fl)1e-8 ? (fl)1e-8 : dist_ang);
		const vec& fxyz = m_reactive.receptor_frame_xyz;
		fl bx = fxyz[0]-rxyz[0], by = fxyz[1]-rxyz[1], bz = fxyz[2]-rxyz[2];
		fl b_len = std::sqrt(bx*bx+by*by+bz*bz);
		if (b_len < (fl)1e-8) b_len = (fl)1e-8;
		fl cos_th = (dx*bx + dy*by + dz*bz) / (r*b_len);
		if (cos_th > (fl)1)  cos_th = (fl)1;
		if (cos_th < (fl)-1) cos_th = (fl)-1;
		const fl pi_f = (fl)3.14159265358979;
		angle_deg = std::acos(cos_th) * (fl)180.0 / pi_f;
	}
	return true;
}

bool ad4cache::debug_check_reactive_gradient(const model& m, fl eps, std::ostream& out) const {
	if (!m_reactive.ready() || eps <= 0) return false;

	fl distance_energy = 0, angle_energy = 0;
	vec analytic_minus_force;
	analytic_minus_force.assign(0);
	if (!get_reactive_terms(m, distance_energy, angle_energy, &analytic_minus_force)) return false;

	std::ios::fmtflags old_flags = out.flags();
	std::streamsize old_precision = out.precision();
	out.setf(std::ios::scientific, std::ios::floatfield);
	out << std::setprecision(6);
	out << "Reactive gradient check (central difference, eps=" << eps << ")\n";
	fl max_abs_diff = 0;
	for (sz axis = 0; axis < 3; ++axis) {
		model plus = m;
		model minus = m;
		plus.coords[m_reactive.ligand_atom_index][axis] += eps;
		minus.coords[m_reactive.ligand_atom_index][axis] -= eps;

		fl eplus_d = 0, eplus_a = 0;
		fl eminus_d = 0, eminus_a = 0;
		get_reactive_terms(plus, eplus_d, eplus_a, 0);
		get_reactive_terms(minus, eminus_d, eminus_a, 0);

		fl numeric_grad = ((eplus_d + eplus_a) - (eminus_d + eminus_a)) / (2.0 * eps);
		fl numeric_minus_force = -numeric_grad;
		fl diff = std::abs(analytic_minus_force[axis] - numeric_minus_force);
		if (diff > max_abs_diff) max_abs_diff = diff;
		out << "  axis " << axis
		    << ": analytic_minus_force=" << std::scientific << std::setprecision(6) << analytic_minus_force[axis]
		    << "  numeric_minus_force=" << std::scientific << std::setprecision(6) << numeric_minus_force
		    << "  abs_diff=" << std::scientific << std::setprecision(6) << diff << "\n";
	}
	// ── Lig-frame atom gradient check (P2 term, if enabled) ─────────────────
	if (m_reactive.has_ligand_frame &&
	    m_reactive.ligand_frame_atom_index < m.coords.size()) {
		vec analytic_lfmf;
		analytic_lfmf.assign(0);
		fl dummy_de = 0, dummy_ae = 0;
		get_reactive_terms(m, dummy_de, dummy_ae, nullptr, &analytic_lfmf);
		out << "  [lig_frame atom " << m_reactive.ligand_frame_atom_index << "]\n";
		for (sz axis = 0; axis < 3; ++axis) {
			model plus = m;
			model minus = m;
			plus.coords [m_reactive.ligand_frame_atom_index][axis] += eps;
			minus.coords[m_reactive.ligand_frame_atom_index][axis] -= eps;
			fl eplus_d = 0, eplus_a = 0, eminus_d = 0, eminus_a = 0;
			get_reactive_terms(plus,  eplus_d,  eplus_a,  0);
			get_reactive_terms(minus, eminus_d, eminus_a, 0);
			fl numeric_mf = -((eplus_d + eplus_a) - (eminus_d + eminus_a)) / (2.0 * eps);
			fl diff = std::abs(analytic_lfmf[axis] - numeric_mf);
			if (diff > max_abs_diff) max_abs_diff = diff;
			out << "  axis " << axis
			    << ": analytic_frame_mf=" << std::scientific << std::setprecision(6) << analytic_lfmf[axis]
			    << "  numeric_frame_mf="  << std::scientific << std::setprecision(6) << numeric_mf
			    << "  abs_diff="          << std::scientific << std::setprecision(6) << diff << "\n";
		}
	}
	out << "  max_abs_diff=" << std::scientific << std::setprecision(6) << max_abs_diff << "\n";
	out.flags(old_flags);
	out.precision(old_precision);
	return true;
}

std::string get_adtype_str(sz& t) {
	switch(t) {
		case AD_TYPE_C : return "C";
		case AD_TYPE_A : return "A";
		case AD_TYPE_N : return "N";
		case AD_TYPE_O : return "O";
		case AD_TYPE_P : return "P";
		case AD_TYPE_S : return "S";
		case AD_TYPE_H : return "H";
		case AD_TYPE_F : return "F";
		case AD_TYPE_I : return "I";
		case AD_TYPE_NA: return "NA";
		case AD_TYPE_OA: return "OA";
		case AD_TYPE_SA: return "SA";
		case AD_TYPE_HD: return "HD";
		case AD_TYPE_Mg: return "Mg";
		case AD_TYPE_Mn: return "Mn";
		case AD_TYPE_Zn: return "Zn";
		case AD_TYPE_Ca: return "Ca";
		case AD_TYPE_Fe: return "Fe";
		case AD_TYPE_Cl: return "Cl";
		case AD_TYPE_Br: return "Br";
		case AD_TYPE_Si: return "Si";
		case AD_TYPE_At: return "At";
		case AD_TYPE_W : return "W";
		// --- Full-metal AD4 extension ---
		case AD_TYPE_He: return "He";
		case AD_TYPE_Li: return "Li";
		case AD_TYPE_Be: return "Be";
		case AD_TYPE_B : return "B";
		case AD_TYPE_Ne: return "Ne";
		case AD_TYPE_Na: return "Na";
		case AD_TYPE_Al: return "Al";
		case AD_TYPE_K : return "K";
		case AD_TYPE_Sc: return "Sc";
		case AD_TYPE_Ti: return "Ti";
		case AD_TYPE_V : return "V";
		case AD_TYPE_Co: return "Co";
		case AD_TYPE_Ni: return "Ni";
		case AD_TYPE_Cu: return "Cu";
		case AD_TYPE_Ga: return "Ga";
		case AD_TYPE_Ge: return "Ge";
		case AD_TYPE_As: return "As";
		case AD_TYPE_Kr: return "Kr";
		case AD_TYPE_Rb: return "Rb";
		case AD_TYPE_Sr: return "Sr";
		case AD_TYPE_Y : return "Y";
		case AD_TYPE_Zr: return "Zr";
		case AD_TYPE_Nb: return "Nb";
		case AD_TYPE_Mo: return "Mo";
		case AD_TYPE_Tc: return "Tc";
		case AD_TYPE_Ru: return "Ru";
		case AD_TYPE_Rh: return "Rh";
		case AD_TYPE_Pd: return "Pd";
		case AD_TYPE_Ag: return "Ag";
		case AD_TYPE_Cd: return "Cd";
		case AD_TYPE_In: return "In";
		case AD_TYPE_Sn: return "Sn";
		case AD_TYPE_Sb: return "Sb";
		case AD_TYPE_Te: return "Te";
		case AD_TYPE_Xe: return "Xe";
		case AD_TYPE_Cs: return "Cs";
		case AD_TYPE_Ba: return "Ba";
		case AD_TYPE_La: return "La";
		case AD_TYPE_Ce: return "Ce";
		case AD_TYPE_Pr: return "Pr";
		case AD_TYPE_Nd: return "Nd";
		case AD_TYPE_Pm: return "Pm";
		case AD_TYPE_Sm: return "Sm";
		case AD_TYPE_Eu: return "Eu";
		case AD_TYPE_Gd: return "Gd";
		case AD_TYPE_Tb: return "Tb";
		case AD_TYPE_Dy: return "Dy";
		case AD_TYPE_Ho: return "Ho";
		case AD_TYPE_Er: return "Er";
		case AD_TYPE_Tm: return "Tm";
		case AD_TYPE_Yb: return "Yb";
		case AD_TYPE_Lu: return "Lu";
		case AD_TYPE_Hf: return "Hf";
		case AD_TYPE_Ta: return "Ta";
		case AD_TYPE_Re: return "Re";
		case AD_TYPE_Os: return "Os";
		case AD_TYPE_Ir: return "Ir";
		case AD_TYPE_Pt: return "Pt";
		case AD_TYPE_Au: return "Au";
		case AD_TYPE_Hg: return "Hg";
		case AD_TYPE_Tl: return "Tl";
		case AD_TYPE_Pb: return "Pb";
		case AD_TYPE_Bi: return "Bi";
		case AD_TYPE_Po: return "Po";
		case AD_TYPE_Rn: return "Rn";
		case AD_TYPE_Fr: return "Fr";
		case AD_TYPE_Ra: return "Ra";
		case AD_TYPE_Ac: return "Ac";
		case AD_TYPE_Th: return "Th";
		case AD_TYPE_Pa: return "Pa";
		case AD_TYPE_U : return "U";
		case AD_TYPE_Np: return "Np";
		case AD_TYPE_Pu: return "Pu";
		case AD_TYPE_Am: return "Am";
		case AD_TYPE_Cm: return "Cm";
		case AD_TYPE_Bk: return "Bk";
		case AD_TYPE_Cf: return "Cf";
		case AD_TYPE_E : return "E";
		case AD_TYPE_Fm: return "Fm";
		case AD_TYPE_Cr1: return "Cr1";
		case AD_TYPE_Tg:  return "Tg";
		case AD_TYPE_Se:  return "Se";
		case AD_TYPE_TZ:  return "TZ";
		case AD_TYPE_SQ:  return "SQ";
		case AD_TYPE_MH:  return "MH";
		case AD_TYPE_JT:  return "JT";
		default: VINA_CHECK(false);
	}
	return std::string(); // placate compiler
}

static bool ad4_metal_probe_name(sz ad, std::string& name) {
	switch (ad) {
		case AD_TYPE_NA: name = "NA"; return true;
		case AD_TYPE_N : name = "N";  return true;
		case AD_TYPE_OA: name = "OA"; return true;
		case AD_TYPE_O : name = "O";  return true;
		case AD_TYPE_SA: name = "SA"; return true;
		case AD_TYPE_S : name = "S";  return true;
		default: return false;
	}
}

static fl ad4_gauss_fit(fl x, fl mu, fl sigma) {
	if (sigma <= 0) return 0;
	fl d = (x - mu) / sigma;
	return std::exp(-0.5 * d * d);
}

fl ad4cache::get_metal_rerank_terms(const model& m, fl& geometry_term, fl& water_term, fl& jt_term) const {
	geometry_term = 0;
	water_term = 0;
	jt_term = 0;
	if (!m_metal_state.enabled()) return 0;

	struct donor_atom_info {
		std::string probe;
		vec xyz;
	};
	std::vector<donor_atom_info> donors;
	donors.reserve(m.num_movable_atoms());
	VINA_FOR(i, m.num_movable_atoms()) {
		std::string probe;
		sz ad = m.get_atom(i).get(atom_type::AD);
		if (!ad4_metal_probe_name(ad, probe)) continue;
		donor_atom_info da;
		da.probe = probe;
		da.xyz = m.movable_coords(i);
		donors.push_back(da);
	}

	for (const auto& site : m_metal_state.sites) {
		fl best_geo = 0;
		fl best_axial = 0;
		fl best_equatorial = 0;
		for (const auto& donor : donors) {
			fl dx = donor.xyz[0] - site.metal_xyz[0];
			fl dy = donor.xyz[1] - site.metal_xyz[1];
			fl dz = donor.xyz[2] - site.metal_xyz[2];
			fl dist = std::sqrt(dx*dx + dy*dy + dz*dz);
			for (const auto& ov : site.direct_overrides) {
				if (ov.probe != donor.probe || ov.eps <= 0.0) continue;
				best_geo = std::max(best_geo, (fl)(ov.eps / 20.0) * ad4_gauss_fit(dist, (fl)ov.r_eq, 0.30));
			}
			if (site.jt_enabled && dist > 0.1) {
				fl dot = std::fabs((dx*site.jt_axis[0] + dy*site.jt_axis[1] + dz*site.jt_axis[2]) / dist);
				best_axial = std::max(best_axial, dot * ad4_gauss_fit(dist, (fl)site.jt_axial_dist, 0.35));
				best_equatorial = std::max(best_equatorial, (1.0f - dot) * ad4_gauss_fit(dist, (fl)site.coord_dist, 0.30));
			}
		}
		geometry_term -= 1.25 * best_geo;

		fl site_water = 0;
		for (const auto& ws : site.bridge_waters) {
			fl best_site = 0;
			for (const auto& donor : donors) {
				fl dx = donor.xyz[0] - ws.xyz[0];
				fl dy = donor.xyz[1] - ws.xyz[1];
				fl dz = donor.xyz[2] - ws.xyz[2];
				fl dist = std::sqrt(dx*dx + dy*dy + dz*dz);
				best_site = std::max(best_site, (fl)ws.weight * ad4_gauss_fit(dist, (fl)ws.target_dist, 0.40));
			}
			site_water += best_site;
		}
		water_term -= 0.90 * site_water;

		if (site.jt_enabled)
			jt_term -= 0.50 * (best_axial + 0.50 * best_equatorial);
	}

	return geometry_term + water_term + jt_term;
}

// Smooth metal-geometry soft constraint with analytical gradient.
// Uses the log-sum-exp identity to approximate max_j(v_j) ≈ (1/T)*log(Σ exp(T*v_j)),
// which is differentiable and converges to the hard max as T→∞.
// Weight w = m_metal_soft_weight scales the contribution during search.
fl ad4cache::eval_metal_soft_grad(model& m, fl T) const {
    if (!m_metal_state.enabled() || m_metal_soft_weight <= (fl)0) return (fl)0;

    fl e = (fl)0;
    const fl w = m_metal_soft_weight;

    VINA_FOR(i, m.num_movable_atoms()) {
        std::string probe;
        sz ad = m.get_atom(i).get(atom_type::AD);
        if (!ad4_metal_probe_name(ad, probe)) continue;
        const vec& xyz = m.movable_coords(i);

        for (const auto& site : m_metal_state.sites) {
            // ── geometry Gaussian soft-max ────────────────────────────────────
            fl sum_exp = (fl)0;
            fl sum_wt_dx = (fl)0, sum_wt_dy = (fl)0, sum_wt_dz = (fl)0;
            for (const auto& ov : site.direct_overrides) {
                if (ov.probe != probe || ov.eps <= 0.0) continue;
                fl dx = xyz[0] - site.metal_xyz[0];
                fl dy = xyz[1] - site.metal_xyz[1];
                fl dz = xyz[2] - site.metal_xyz[2];
                fl dist = std::sqrt(dx*dx + dy*dy + dz*dz);
                if (dist < (fl)1e-8) continue;
                fl sigma = (fl)0.30;
                fl dr   = dist - (fl)ov.r_eq;
                fl gval = (fl)(ov.eps / 20.0) * std::exp(-(fl)0.5 * dr*dr / (sigma*sigma));
                fl eT   = std::exp(T * gval);
                sum_exp += eT;
                // dG/d(atom_coord) = gval * (-dr/(sigma^2*dist)) * (coord_vec)
                fl dg_dd = gval * (-dr / (sigma*sigma * dist));
                sum_wt_dx += eT * dg_dd * dx;
                sum_wt_dy += eT * dg_dd * dy;
                sum_wt_dz += eT * dg_dd * dz;
            }
            if (sum_exp > (fl)0) {
                fl lse = (fl)1.0 / T * std::log(sum_exp);          // log-sum-exp value
                e += -1.25f * w * lse;                              // E_geo contribution
                fl inv = (fl)1.0 / sum_exp;                         // softmax denominator
                m.minus_forces[i][0] += -1.25f * w * sum_wt_dx * inv;
                m.minus_forces[i][1] += -1.25f * w * sum_wt_dy * inv;
                m.minus_forces[i][2] += -1.25f * w * sum_wt_dz * inv;
            }

            // ── bridge water Gaussian soft-max ───────────────────────────────
            for (const auto& ws : site.bridge_waters) {
                fl dx = xyz[0] - ws.xyz[0];
                fl dy = xyz[1] - ws.xyz[1];
                fl dz = xyz[2] - ws.xyz[2];
                fl dist = std::sqrt(dx*dx + dy*dy + dz*dz);
                if (dist < (fl)1e-8) continue;
                fl sigma = (fl)0.40;
                fl dr   = dist - (fl)ws.target_dist;
                fl gval = (fl)ws.weight * std::exp(-(fl)0.5 * dr*dr / (sigma*sigma));
                e += -(fl)0.90 * w * gval;
                fl dg_dd = gval * (-dr / (sigma*sigma * dist));
                m.minus_forces[i][0] += -(fl)0.90 * w * dg_dd * dx;
                m.minus_forces[i][1] += -(fl)0.90 * w * dg_dd * dy;
                m.minus_forces[i][2] += -(fl)0.90 * w * dg_dd * dz;
            }
        }
    }
    return e;
}

fl ad4cache::eval(const model& m, fl v) const {
	fl e = 0;
	sz nat = num_atom_types(atom_type::AD);

	VINA_FOR(i, m.num_movable_atoms()) {
		if(!m.is_atom_in_ligand(i)) continue; // we only want ligand interaction
		const atom& a = m.atoms[i];
		sz t = a.get(atom_type::AD);

		switch (t)
		{
			case AD_TYPE_G0:
			case AD_TYPE_G1:
			case AD_TYPE_G2:
			case AD_TYPE_G3:
				continue;
			case AD_TYPE_CG0:
			case AD_TYPE_CG1:
			case AD_TYPE_CG2:
			case AD_TYPE_CG3:
				t = AD_TYPE_C;
				break;
		}

		// hybrid mode: suppress or scale VdW+HB for the reactive ligand atom
		if (m_reactive.ready() &&
		    m_reactive.mode == reactive_mode::hybrid &&
		    i == m_reactive.ligand_atom_index) {
			fl vdw_scale = (fl)m_reactive.hybrid_vdw_scale;
			if (vdw_scale > 0.0f)
				e += m_grids[t].evaluate(m.coords[i], m_slope, v) * vdw_scale;
			e += m_grids[AD_TYPE_SIZE].evaluate(m.coords[i], m_slope, v) * a.charge;
			e += m_grids[AD_TYPE_SIZE + 1].evaluate(m.coords[i], m_slope, v) * std::abs(a.charge);
			continue;
		}

		// HB + vdW
		const grid& g = m_grids[t];
		e += g.evaluate(m.coords[i], m_slope, v);

		// elec
		const grid& ge = m_grids[AD_TYPE_SIZE];
		e += ge.evaluate(m.coords[i], m_slope, v) * a.charge;

		// desolv
		const grid& gd = m_grids[AD_TYPE_SIZE + 1];
		e += gd.evaluate(m.coords[i], m_slope, v) * std::abs(a.charge);
	}
	fl reactive_distance = 0, reactive_angle = 0;
	if (get_reactive_terms(m, reactive_distance, reactive_angle, 0))
		e += reactive_distance + reactive_angle;
	return e;
}

fl ad4cache::eval_intra(model& m, fl v) const {
	fl e = 0;
	sz nat = num_atom_types(atom_type::AD);

	VINA_FOR(i, m.num_movable_atoms()) {
		if(m.is_atom_in_ligand(i)) continue; // we only want flex-rigid interaction
		const atom& a = m.atoms[i];
		sz t = a.get(atom_type::AD);

		switch (t)
		{
			case AD_TYPE_G0:
			case AD_TYPE_G1:
			case AD_TYPE_G2:
			case AD_TYPE_G3:
				continue;
			case AD_TYPE_CG0:
			case AD_TYPE_CG1:
			case AD_TYPE_CG2:
			case AD_TYPE_CG3:
				t = AD_TYPE_C;
				break;
		}

		// HB + vdW
		const grid& g = m_grids[t];
		e += g.evaluate(m.coords[i], m_slope, v);

		// elec
		const grid& ge = m_grids[AD_TYPE_SIZE];
		e += ge.evaluate(m.coords[i], m_slope, v) * a.charge;

		// desolv
		const grid& gd = m_grids[AD_TYPE_SIZE + 1];
		e += gd.evaluate(m.coords[i], m_slope, v) * std::abs(a.charge);
	}
	return e;
}

fl ad4cache::eval_deriv(model& m, fl v) const { // sets m.minus_forces
	fl e = 0;
	sz nat = num_atom_types(atom_type::AD);

	VINA_FOR(i, m.num_movable_atoms()) {
		const atom& a = m.atoms[i];
		sz t = a.get(atom_type::AD);

		switch (t)
		{
			case AD_TYPE_G0:
			case AD_TYPE_G1:
			case AD_TYPE_G2:
			case AD_TYPE_G3:
				m.minus_forces[i].assign(0);
				continue;
			case AD_TYPE_CG0:
			case AD_TYPE_CG1:
			case AD_TYPE_CG2:
			case AD_TYPE_CG3:
				t = AD_TYPE_C;
				break;
		}

		// hybrid mode: suppress or scale VdW+HB for the reactive ligand atom
		if (m_reactive.ready() &&
		    m_reactive.mode == reactive_mode::hybrid &&
		    m.is_atom_in_ligand(i) &&
		    i == m_reactive.ligand_atom_index) {
			fl vdw_scale = (fl)m_reactive.hybrid_vdw_scale;
			vec d2;
			if (vdw_scale > 0.0f) {
				e += m_grids[t].evaluate(m.coords[i], m_slope, v, d2) * vdw_scale;
				d2 *= vdw_scale;
				m.minus_forces[i] = d2;
			} else {
				m.minus_forces[i].assign(0);
			}
			e += m_grids[AD_TYPE_SIZE].evaluate(m.coords[i], m_slope, v, d2) * a.charge;
			d2 *= a.charge; m.minus_forces[i] += d2;
			e += m_grids[AD_TYPE_SIZE+1].evaluate(m.coords[i], m_slope, v, d2) * std::abs(a.charge);
			d2 *= std::abs(a.charge); m.minus_forces[i] += d2;
			continue;
		}

		// HB + vdW
		vec deriv;
		const grid& g = m_grids[t];
		e += g.evaluate(m.coords[i], m_slope, v, deriv);
		m.minus_forces[i] = deriv;

		// elec
		const grid& ge = m_grids[AD_TYPE_SIZE];
		e += ge.evaluate(m.coords[i], m_slope, v, deriv) * a.charge;
		deriv *= a.charge;
		m.minus_forces[i] += deriv;

		// desolv
		const grid& gd = m_grids[AD_TYPE_SIZE + 1];
		e += gd.evaluate(m.coords[i], m_slope, v, deriv) * std::abs(a.charge);
		deriv *= std::abs(a.charge);
		m.minus_forces[i] += deriv;
	}
	fl reactive_distance = 0, reactive_angle = 0;
	vec reactive_minus_force, frame_minus_force;
	reactive_minus_force.assign(0);
	frame_minus_force.assign(0);
	if (get_reactive_terms(m, reactive_distance, reactive_angle, &reactive_minus_force, &frame_minus_force)) {
		e += reactive_distance + reactive_angle;
		m.minus_forces[m_reactive.ligand_atom_index] += reactive_minus_force;
		if (m_reactive.has_ligand_frame)
			m.minus_forces[m_reactive.ligand_frame_atom_index] += frame_minus_force;
	}
	// Metal geometry soft constraint (gradient-aware, for use during search).
	// Only active when m_metal_soft_weight > 0 (default 0 = off, backwards compatible).
	if (m_metal_soft_weight > (fl)0)
		e += eval_metal_soft_grad(m);
	return e;
}

bool ad4cache::is_in_grid(const model& m, fl margin) const {
	VINA_FOR(i, m.num_movable_atoms()) {
		if(m.atoms[i].is_hydrogen()) continue;

		const vec& a_coords = m.coords[i];
		VINA_FOR_IN(j, m_gd) {
			if(m_gd[j].n_voxels > 0)
				if(a_coords[j] < m_gd[j].begin - margin || a_coords[j] > m_gd[j].end + margin) 
					return false;
		}
	}
	return true;
}

bool ad4cache::are_atom_types_grid_initialized(szv atom_types) const {
	VINA_FOR_IN(i, atom_types) {
		sz t = atom_types[i];

		switch (t)
		{
			case AD_TYPE_G0:
			case AD_TYPE_G1:
			case AD_TYPE_G2:
			case AD_TYPE_G3:
				continue;
			case AD_TYPE_CG0:
			case AD_TYPE_CG1:
			case AD_TYPE_CG2:
			case AD_TYPE_CG3:
				t = AD_TYPE_C;
				break;
		}

		if (!is_atom_type_grid_initialized(t)) {
			std::cerr << "ERROR: Affinity map for atom type " << get_adtype_str(t) << " is not present.\n";
			return false;
		}
	}

	if (!is_atom_type_grid_initialized(AD_TYPE_SIZE)) {
		std::cerr << "ERROR: Electrostatic map is not present.\n";
		return false;
	}

	if (!is_atom_type_grid_initialized(AD_TYPE_SIZE + 1)) {
		std::cerr << "ERROR: Desolvation map is not present.\n";
		return false;
	}

	return true;
}

std::vector<std::string> split(std::string str) {
	std::vector<std::string> fields;
	std::string field;
	std::istringstream iss(str);
	while(std::getline(iss, field, ' '))
	{
		fields.push_back(field);
	};
	return fields;
}

void read_ad4_map(path& filename, std::vector<grid_dims>& gds, grid& g) {
	sz line_counter = 0;
	sz pt_counter = 0;
	sz x = 0;
	sz y = 0;
	sz z = 0;
	grid_dims gd;
	std::string line;
	fl spacing, center, halfspan;

	ifile in(filename);

	while(std::getline(in, line)) {
		line_counter++;
		if (line_counter == 4) {
			std::vector<std::string> fields = split(line);
			spacing = std::atof(fields[1].c_str());
		}
		if (line_counter == 5) {
			std::vector<std::string> fields = split(line);
			VINA_FOR(i, 3) {
				// n_voxels must be EVEN
				// because the number of sampled points in the grid is always ODD
				// (number of sampled points == n_voxels + 1)
				gd[i].n_voxels = std::atoi(fields[i + 1].c_str());
				if (gd[i].n_voxels % 2 == 1) {
					std::cerr << "ERROR: number of voxels (NELEMENTS) must be even\n";
					exit(EXIT_FAILURE);
				}
			}
		}
		if (line_counter == 6) {
			std::vector<std::string> fields = split(line);
			VINA_FOR(i, 3) {
				center = std::atof(fields[i+1].c_str());
				halfspan = (gd[i].n_voxels) * spacing / 2.0;
				gd[i].begin = center - halfspan;
				gd[i].end = center + halfspan;
				// std::cout << center << " " << halfspan << " " << gd[i].begin << " " << gd[i].end << "\n";
			}
			gds.push_back(gd);
			g.init(gd);
		}
		if (line_counter > 6) {
			// std::cout << pt_counter << " " << x << " " << y << " " << z << " " << std::atof(line.c_str()) << "\n";
			g.m_data(x, y, z) = std::atof(line.c_str());
			y += sz(x == (gd[0].n_voxels + 1));
			z += sz(y == (gd[1].n_voxels + 1));
			x = x % (gd[0].n_voxels + 1);
			y = y % (gd[1].n_voxels + 1);
			pt_counter++;
			x++;
		}
	} // line loop
}

void ad4cache::read(const std::string& map_prefix) {

	std::string type, filename;
	std::vector<grid_dims> gds; // to check all maps have same dims (TODO)

	bool got_C_already = false;

	VINA_FOR(atom_type, AD_TYPE_SIZE){
		sz t = atom_type;

		switch (t)
	{
		case AD_TYPE_G0:
		case AD_TYPE_G1:
		case AD_TYPE_G2:
		case AD_TYPE_G3:
			continue;
		case AD_TYPE_CG0:
		case AD_TYPE_CG1:
		case AD_TYPE_CG2:
		case AD_TYPE_CG3:
			if (got_C_already) continue;
			t = AD_TYPE_C;
			got_C_already = true;
			break;
	}

		type = get_adtype_str(t);
		filename = map_prefix + "." + type + ".map";
		path p(filename);
		if (fs::exists(p)) {
			read_ad4_map(p, gds, m_grids[t]);

		} // if file exists
	} // map loop

	//  elec map
	filename = map_prefix + ".e.map";
	path pe(filename);
	read_ad4_map(pe, gds, m_grids[AD_TYPE_SIZE]);

	//  dsolv map
	filename = map_prefix + ".d.map";
	path pd(filename);
	read_ad4_map(pd, gds, m_grids[AD_TYPE_SIZE + 1]);

	// Store in Ad4cache object
	m_gd = gds[0];
}

void ad4cache::populate_from_data(const grid_dims& gd,
                                   const std::vector<sz>& ad_types,
                                   const std::vector<std::vector<double>>& aff_maps,
                                   const std::vector<double>& elec_map,
                                   const std::vector<double>& desolv_map)
{
    m_gd = gd;
    int px = gd[0].n_voxels + 1;  // sample points per axis (n_voxels+1)
    int py = gd[1].n_voxels + 1;
    int pz = gd[2].n_voxels + 1;

    // ---- Affinity maps ----
    for (sz k = 0; k < ad_types.size(); k++) {
        sz t = ad_types[k];
        if (t >= (sz)m_grids.size()) continue;
        m_grids[t].init(gd);
        const std::vector<double>& flat = aff_maps[k];
        for (int iz = 0; iz < pz; iz++)
            for (int iy = 0; iy < py; iy++)
                for (int ix = 0; ix < px; ix++)
                    m_grids[t].m_data(ix, iy, iz) =
                        (fl)flat[iz * py * px + iy * px + ix];
    }

    // ---- Electrostatic map (index AD_TYPE_SIZE) ----
    {
        sz t = AD_TYPE_SIZE;
        m_grids[t].init(gd);
        for (int iz = 0; iz < pz; iz++)
            for (int iy = 0; iy < py; iy++)
                for (int ix = 0; ix < px; ix++)
                    m_grids[t].m_data(ix, iy, iz) =
                        (fl)elec_map[iz * py * px + iy * px + ix];
    }

    // ---- Desolvation map (index AD_TYPE_SIZE + 1) ----
    {
        sz t = AD_TYPE_SIZE + 1;
        m_grids[t].init(gd);
        for (int iz = 0; iz < pz; iz++)
            for (int iy = 0; iy < py; iy++)
                for (int ix = 0; ix < px; ix++)
                    m_grids[t].m_data(ix, iy, iz) =
                        (fl)desolv_map[iz * py * px + iy * px + ix];
    }
}

void ad4cache::write(const std::string& out_prefix, const szv& atom_types, const std::string& gpf_filename,
					      const std::string& fld_filename, const std::string& receptor_filename) {
	std::string atom_type;
	std::string filename;
	bool got_C_already = false;

	VINA_FOR_IN(i, atom_types) {
		sz t = atom_types[i];

		switch (t)
		{
			case AD_TYPE_G0:
			case AD_TYPE_G1:
			case AD_TYPE_G2:
			case AD_TYPE_G3:
				continue;
			case AD_TYPE_CG0:
			case AD_TYPE_CG1:
			case AD_TYPE_CG2:
			case AD_TYPE_CG3:
				t = AD_TYPE_C;
				break;
		}

		if (t == AD_TYPE_C && got_C_already)
			continue;
		if (t == AD_TYPE_C)
			got_C_already = true;

		if (m_grids[t].initialized()) {
			if (t < AD_TYPE_SIZE)
				atom_type = get_adtype_str(t);
			else if (t == AD_TYPE_SIZE)
				atom_type = "e";
			else if (t == AD_TYPE_SIZE + 1)
				atom_type = "d";

			filename = out_prefix + "." + atom_type + ".map";

			path p(filename);
			ofile out(p);

			// write header
			out << "GRID_PARAMETER_FILE " << gpf_filename << "\n";
			out << "GRID_DATA_FILE " << fld_filename << "\n";
			out << "MACROMOLECULE " << receptor_filename << "\n";

			// m_factor_inv is spacing
			// check that it's the same in every dimension (it must be)
			// check that == operator is OK
			if ((m_grids[t].m_factor_inv[0] != m_grids[t].m_factor_inv[1]) & (m_grids[t].m_factor_inv[0] != m_grids[t].m_factor_inv[2])) {
				printf("m_factor_inv x=%f, y=%f, z=%f\n", m_grids[t].m_factor_inv[0], m_grids[t].m_factor_inv[1], m_grids[t].m_factor_inv[2]);
				return;
			}

			out << "SPACING " << m_grids[t].m_factor_inv[0] << "\n";

			// The number of elements in the grid is an odd number. But NELEMENTS has to be an even number.
			int size_x = (m_grids[t].m_data.dim0() % 2 == 0) ? m_grids[t].m_data.dim0() : m_grids[t].m_data.dim0() - 1;
			int size_y = (m_grids[t].m_data.dim1() % 2 == 0) ? m_grids[t].m_data.dim1() : m_grids[t].m_data.dim1() - 1;
			int size_z = (m_grids[t].m_data.dim2() % 2 == 0) ? m_grids[t].m_data.dim2() : m_grids[t].m_data.dim2() - 1;
			out << "NELEMENTS " << size_x << " " << size_y  << " " << size_z << "\n";

			// center
			fl cx = m_grids[t].m_init[0] + m_grids[t].m_range[0] * 0.5;
			fl cy = m_grids[t].m_init[1] + m_grids[t].m_range[1] * 0.5;
			fl cz = m_grids[t].m_init[2] + m_grids[t].m_range[2] * 0.5;
			out << "CENTER " << cx << " " << cy << " " << cz << "\n";

			// write data
			VINA_FOR(z, m_grids[t].m_data.dim2()) {
				VINA_FOR(y, m_grids[t].m_data.dim1()) {
					VINA_FOR(x, m_grids[t].m_data.dim0()) {
						out << std::setprecision(4) << m_grids[t].m_data(x, y, z) << "\n"; // slow?
					} // x
				} // y
			} // z
		} // map initialized
	} // map atom type
} // cache::write
