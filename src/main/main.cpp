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

#include <iostream>
#include <set>
#include <string>
#include <vector> // ligand paths
#include <array>
#include <exception>
#include <cstdlib>
#include <algorithm>
#include <cstring>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <limits>
#include "parse_error.h"
#include <boost/program_options.hpp>
#include "vina.h"
#include "utils.h"
#include "scoring_function.h"
#include <unordered_map>
#include <boost/filesystem.hpp>

#ifdef _WIN32
#include <io.h>
#include <windows.h>
#else
#include <unistd.h>
#endif

struct ag4_detect_atom_rec {
	std::string type;
	double x, y, z;
};

static double ag4_bvs_r0(ag4_metal_mode mode, const std::string& donor_type) {
	if (mode == ag4_metal_mode::fe2) {
		if (donor_type == "OA" || donor_type == "O") return 1.76;
		if (donor_type == "NA" || donor_type == "N") return 1.79;
		if (donor_type == "SA" || donor_type == "S") return 2.05;
	} else if (mode == ag4_metal_mode::fe3) {
		if (donor_type == "OA" || donor_type == "O") return 1.73;
		if (donor_type == "NA" || donor_type == "N") return 1.76;
		if (donor_type == "SA" || donor_type == "S") return 1.98;
	} else if (mode == ag4_metal_mode::cu1) {
		if (donor_type == "OA" || donor_type == "O") return 1.72;
		if (donor_type == "NA" || donor_type == "N") return 1.74;
		if (donor_type == "SA" || donor_type == "S") return 1.96;
	} else if (mode == ag4_metal_mode::cu2 || mode == ag4_metal_mode::cu2_jt) {
		if (donor_type == "OA" || donor_type == "O") return 1.68;
		if (donor_type == "NA" || donor_type == "N") return 1.70;
		if (donor_type == "SA" || donor_type == "S") return 1.89;
	} else if (mode == ag4_metal_mode::mn2) {
		if (donor_type == "OA" || donor_type == "O") return 1.79;
		if (donor_type == "NA" || donor_type == "N") return 1.80;
	} else if (mode == ag4_metal_mode::mn3 || mode == ag4_metal_mode::mn3_jt) {
		if (donor_type == "OA" || donor_type == "O") return 1.76;
		if (donor_type == "NA" || donor_type == "N") return 1.77;
	} else if (mode == ag4_metal_mode::co2) {
		if (donor_type == "OA" || donor_type == "O") return 1.75;
		if (donor_type == "NA" || donor_type == "N") return 1.77;
	} else if (mode == ag4_metal_mode::co3) {
		if (donor_type == "OA" || donor_type == "O") return 1.70;
		if (donor_type == "NA" || donor_type == "N") return 1.72;
	} else if (mode == ag4_metal_mode::v4) {
		if (donor_type == "OA" || donor_type == "O") return 1.78;
		if (donor_type == "NA" || donor_type == "N") return 1.80;
	} else if (mode == ag4_metal_mode::v5) {
		if (donor_type == "OA" || donor_type == "O") return 1.80;
		if (donor_type == "NA" || donor_type == "N") return 1.82;
	} else if (mode == ag4_metal_mode::mo4) {
		if (donor_type == "OA" || donor_type == "O") return 1.90;
		if (donor_type == "SA" || donor_type == "S") return 2.12;
	} else if (mode == ag4_metal_mode::mo6) {
		if (donor_type == "OA" || donor_type == "O") return 1.86;
		if (donor_type == "NA" || donor_type == "N") return 1.92;
	// Ni BVS r0: Brese & O'Keeffe 1991; Ni2+-O=1.654, Ni2+-N=1.679, Ni2+-S=1.978
	} else if (mode == ag4_metal_mode::ni2) {
		if (donor_type == "OA" || donor_type == "O") return 1.654;
		if (donor_type == "NA" || donor_type == "N") return 1.679;
		if (donor_type == "SA" || donor_type == "S") return 1.978;
	// Ni3+-O r0 ~0.03 shorter; empirically from Ni-Fe hydrogenase CSD survey
	} else if (mode == ag4_metal_mode::ni3) {
		if (donor_type == "OA" || donor_type == "O") return 1.620;
		if (donor_type == "NA" || donor_type == "N") return 1.650;
		if (donor_type == "SA" || donor_type == "S") return 1.950;
	}
	return 0.0;
}

static double ag4_compute_bvs(ag4_metal_mode mode, const ag4_detect_atom_rec& metal,
	                          const std::vector<ag4_detect_atom_rec>& donors) {
	const double B = 0.37;
	double sum = 0.0;
	for (const auto& d : donors) {
		double dx = d.x - metal.x, dy = d.y - metal.y, dz = d.z - metal.z;
		double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
		if (dist > 3.2) continue;
		double r0 = ag4_bvs_r0(mode, d.type);
		if (r0 <= 0.0) continue;
		sum += std::exp((r0 - dist) / B);
	}
	return sum;
}

static ag4_metal_mode ag4_bvs_pick_mode(ag4_metal_mode base, const ag4_detect_atom_rec& metal,
	                                     const std::vector<ag4_detect_atom_rec>& donors) {
	std::vector<std::pair<ag4_metal_mode, double> > candidates;
	if (base == ag4_metal_mode::fe) candidates = {{ag4_metal_mode::fe2, 2.0}, {ag4_metal_mode::fe3, 3.0}};
	else if (base == ag4_metal_mode::cu) candidates = {{ag4_metal_mode::cu1, 1.0}, {ag4_metal_mode::cu2, 2.0}};
	else if (base == ag4_metal_mode::mn) candidates = {{ag4_metal_mode::mn2, 2.0}, {ag4_metal_mode::mn3, 3.0}};
	else if (base == ag4_metal_mode::co) candidates = {{ag4_metal_mode::co2, 2.0}, {ag4_metal_mode::co3, 3.0}};
	else if (base == ag4_metal_mode::v) candidates = {{ag4_metal_mode::v4, 4.0}, {ag4_metal_mode::v5, 5.0}};
	else if (base == ag4_metal_mode::mo) candidates = {{ag4_metal_mode::mo4, 4.0}, {ag4_metal_mode::mo6, 6.0}};
	else if (base == ag4_metal_mode::ni) candidates = {{ag4_metal_mode::ni2, 2.0}, {ag4_metal_mode::ni3, 3.0}};
	// Mo-Fe heteronuclear (nitrogenase FeMo-co) and polynuclear Mn (OEC) are not yet
	// handled here; their coupled BVS requires multi-centre optimisation (future work).
	if (candidates.empty()) return base;

	ag4_metal_mode best_mode = base;
	double best_delta = std::numeric_limits<double>::max();
	for (const auto& cand : candidates) {
		double bvs = ag4_compute_bvs(cand.first, metal, donors);
		double delta = std::fabs(bvs - cand.second);
		if (delta < best_delta) {
			best_delta = delta;
			best_mode = cand.first;
		}
	}
	if (best_delta < 1.25) return best_mode;
	return base;
}

// M5: Auto-detect metal atom type from receptor PDBQT (last token on ATOM/HETATM lines)
static ag4_metal_mode detect_metal_mode_from_pdbqt(const std::string& path) {
	static const std::pair<const char*, ag4_metal_mode> metal_map[] = {
		{"Zn",ag4_metal_mode::zn},
		{"Mg",ag4_metal_mode::mg},   {"Ca",ag4_metal_mode::ca},
		{"Mn",ag4_metal_mode::mn},   {"Fe",ag4_metal_mode::fe},
		{"Co",ag4_metal_mode::co},   {"Ni",ag4_metal_mode::ni},
		{"Cu",ag4_metal_mode::cu},   {"Pt",ag4_metal_mode::pt},
		{"Pd",ag4_metal_mode::pd},   {"Ru",ag4_metal_mode::ru},
		{"Ir",ag4_metal_mode::ir},   {"Au",ag4_metal_mode::au},
		{"Cd",ag4_metal_mode::cd},   {"Hg",ag4_metal_mode::hg},
		{"Na",ag4_metal_mode::na_ion},{"K", ag4_metal_mode::k_ion},
		// Early transition metals
		{"V", ag4_metal_mode::v},    {"Cr1",ag4_metal_mode::cr}, {"Cr",ag4_metal_mode::cr},
		{"Ti",ag4_metal_mode::ti},   {"Sc",ag4_metal_mode::sc},
		{"Y", ag4_metal_mode::y_ion},{"Zr",ag4_metal_mode::zr},
		{"Nb",ag4_metal_mode::nb},   {"Hf",ag4_metal_mode::hf},
		{"Ta",ag4_metal_mode::ta},   {"Tg",ag4_metal_mode::tg_mode},
		{"Mo",ag4_metal_mode::mo},
		// 2nd/3rd row TMs
		{"Rh",ag4_metal_mode::rh},   {"Ag",ag4_metal_mode::ag},
		{"Tc",ag4_metal_mode::tc},   {"Re",ag4_metal_mode::re},
		{"Os",ag4_metal_mode::os},
		// Post-transition / metalloids
		{"Ga",ag4_metal_mode::ga},   {"In",ag4_metal_mode::in_ion},
		{"Sn",ag4_metal_mode::sn},   {"Sb",ag4_metal_mode::sb},
		{"Bi",ag4_metal_mode::bi},   {"Tl",ag4_metal_mode::tl},
		{"Pb",ag4_metal_mode::pb},
		// s-block
		{"Li",ag4_metal_mode::li},   {"Al",ag4_metal_mode::al},
		{"Sr",ag4_metal_mode::sr},   {"Ba",ag4_metal_mode::ba},
		// Lanthanides
		{"La",ag4_metal_mode::la},   {"Ce",ag4_metal_mode::ce},
		{"Pr",ag4_metal_mode::pr},   {"Nd",ag4_metal_mode::nd},
		{"Sm",ag4_metal_mode::sm},   {"Eu",ag4_metal_mode::eu},
		{"Gd",ag4_metal_mode::gd},   {"Tb",ag4_metal_mode::tb},
		{"Dy",ag4_metal_mode::dy},   {"Ho",ag4_metal_mode::ho},
		{"Er",ag4_metal_mode::er},   {"Tm",ag4_metal_mode::tm},
		{"Yb",ag4_metal_mode::yb},   {"Lu",ag4_metal_mode::lu},
		// Metalloids
		{"Se",ag4_metal_mode::se},   {"As",ag4_metal_mode::as_met},
		{"Ge",ag4_metal_mode::ge},
		// Main-group
		{"B", ag4_metal_mode::b},    {"Be",ag4_metal_mode::be},
		{"Te",ag4_metal_mode::te},   {"Po",ag4_metal_mode::po},
		{"At",ag4_metal_mode::at},
		// Heavy alkali
		{"Rb",ag4_metal_mode::rb},   {"Cs",ag4_metal_mode::cs},
		// Pm + heavy alkaline-earth
		{"Pm",ag4_metal_mode::pm},   {"Ra",ag4_metal_mode::ra},
		// Actinides
		{"Ac",ag4_metal_mode::ac},   {"Th",ag4_metal_mode::th},
		{"Pa",ag4_metal_mode::pa},   {"U", ag4_metal_mode::u},
		{"Np",ag4_metal_mode::np},   {"Pu",ag4_metal_mode::pu},
		{"Am",ag4_metal_mode::am},   {"Cm",ag4_metal_mode::cm},
		{"Bk",ag4_metal_mode::bk},   {"Cf",ag4_metal_mode::cf},
		{"E", ag4_metal_mode::es},   {"Fm",ag4_metal_mode::fm},
	};
	std::ifstream f(path);
	if (!f.is_open()) return ag4_metal_mode::none;
	std::set<ag4_metal_mode> found;
	std::string line;
	while (std::getline(f, line)) {
		if (line.size() < 6) continue;
		std::string rec = line.substr(0,6);
		if (rec != "ATOM  " && rec != "HETATM") continue;
		std::istringstream iss(line);
		std::string tok, last;
		while (iss >> tok) last = tok;
		for (const auto& mp : metal_map)
			if (last == mp.first) { found.insert(mp.second); break; }
	}
	if (found.empty()) return ag4_metal_mode::none;
	if (found.size() > 1) return ag4_metal_mode::none; // multi-metal: cannot auto-infer single mode

	// O8: Re-parse to collect donor positions for oxidation-state inference
	ag4_metal_mode base = *found.begin();
	struct DonorCount { int oa=0, na=0, sa=0; };
	std::vector<ag4_detect_atom_rec> metal_pos;
	std::vector<ag4_detect_atom_rec> donor_pos;

	std::ifstream f2(path);
	if (f2.is_open()) {
		std::string line2;
		while (std::getline(f2, line2)) {
			if (line2.size() < 60) continue;
			std::string rec2 = line2.substr(0,6);
			if (rec2 != "ATOM  " && rec2 != "HETATM") continue;
			double x2=0,y2=0,z2=0;
			try { x2=std::stod(line2.substr(30,8)); y2=std::stod(line2.substr(38,8)); z2=std::stod(line2.substr(46,8)); }
			catch(...) { continue; }
			std::istringstream iss2(line2);
			std::string tok2, last2;
			while (iss2 >> tok2) last2=tok2;
			for (const auto& mp : metal_map)
				if (last2 == mp.first && mp.second == base) { metal_pos.push_back({last2,x2,y2,z2}); break; }
			if (last2=="OA" || last2=="O") donor_pos.push_back({last2,x2,y2,z2});
			else if (last2=="NA"||last2=="N") donor_pos.push_back({last2,x2,y2,z2});
			else if (last2=="SA"||last2=="S") donor_pos.push_back({last2,x2,y2,z2});
		}
	}

	if (metal_pos.empty()) return base;
	ag4_metal_mode bvs_mode = ag4_bvs_pick_mode(base, metal_pos[0], donor_pos);
	if (bvs_mode != base) return bvs_mode;

	// Count donors within 3.0 Å of the first metal atom
	DonorCount dc;
	const ag4_detect_atom_rec& M = metal_pos[0];
	for (const auto& d : donor_pos) {
		double dx=d.x-M.x, dy=d.y-M.y, dz=d.z-M.z;
		if (dx*dx+dy*dy+dz*dz < 9.0) { // 3.0Å cutoff
			if (d.type=="OA" || d.type=="O") dc.oa++;
			else if (d.type=="NA" || d.type=="N") dc.na++;
			else if (d.type=="SA" || d.type=="S") dc.sa++;
		}
	}

	// O8: Infer specific oxidation state from donor distribution
	if (base == ag4_metal_mode::fe) {
		if (dc.oa >= 3 && dc.oa > dc.na) return ag4_metal_mode::fe3;
		if (dc.na >= 2) return ag4_metal_mode::fe2;
	} else if (base == ag4_metal_mode::cu) {
		if (dc.sa >= 1) return ag4_metal_mode::cu1;
		return ag4_metal_mode::cu2;
	} else if (base == ag4_metal_mode::mn) {
		if (dc.oa >= 4) return ag4_metal_mode::mn3;
		return ag4_metal_mode::mn2;
	} else if (base == ag4_metal_mode::co) {
		if (dc.na >= 4) return ag4_metal_mode::co3;
		return ag4_metal_mode::co2;
	} else if (base == ag4_metal_mode::v) {
		if (dc.oa >= 3) return ag4_metal_mode::v5;
		return ag4_metal_mode::v4;
	} else if (base == ag4_metal_mode::mo) {
		if (dc.sa >= 2) return ag4_metal_mode::mo4;
		if (dc.oa >= 3) return ag4_metal_mode::mo6;
	}
	return base;
}
static const char* ag4_metal_mode_to_str(ag4_metal_mode m) {
	switch (m) {
		case ag4_metal_mode::zn: return "zn";
		case ag4_metal_mode::mg: return "mg"; case ag4_metal_mode::ca: return "ca";
		case ag4_metal_mode::mn: return "mn"; case ag4_metal_mode::fe: return "fe";
		case ag4_metal_mode::co: return "co"; case ag4_metal_mode::ni: return "ni";
		case ag4_metal_mode::cu: return "cu"; case ag4_metal_mode::pt: return "pt";
		case ag4_metal_mode::pd: return "pd"; case ag4_metal_mode::ru: return "ru";
		case ag4_metal_mode::ir: return "ir"; case ag4_metal_mode::au: return "au";
		case ag4_metal_mode::cd: return "cd"; case ag4_metal_mode::hg: return "hg";
		case ag4_metal_mode::na_ion: return "na"; case ag4_metal_mode::k_ion: return "k";
		case ag4_metal_mode::v:      return "v";   case ag4_metal_mode::cr:     return "cr";
		case ag4_metal_mode::ti:     return "ti";  case ag4_metal_mode::sc:     return "sc";
		case ag4_metal_mode::y_ion:  return "y";   case ag4_metal_mode::zr:     return "zr";
		case ag4_metal_mode::nb:     return "nb";  case ag4_metal_mode::hf:     return "hf";
		case ag4_metal_mode::ta:     return "ta";  case ag4_metal_mode::tg_mode:return "w";
		case ag4_metal_mode::mo:     return "mo";
		case ag4_metal_mode::rh:     return "rh";  case ag4_metal_mode::ag:     return "ag";
		case ag4_metal_mode::tc:     return "tc";  case ag4_metal_mode::re:     return "re";
		case ag4_metal_mode::os:     return "os";
		case ag4_metal_mode::ga:     return "ga";  case ag4_metal_mode::in_ion: return "in";
		case ag4_metal_mode::sn:     return "sn";  case ag4_metal_mode::sb:     return "sb";
		case ag4_metal_mode::bi:     return "bi";  case ag4_metal_mode::tl:     return "tl";
		case ag4_metal_mode::pb:     return "pb";
		case ag4_metal_mode::li:     return "li";  case ag4_metal_mode::al:     return "al";
		case ag4_metal_mode::sr:     return "sr";  case ag4_metal_mode::ba:     return "ba";
		case ag4_metal_mode::la:     return "la";  case ag4_metal_mode::ce:     return "ce";
		case ag4_metal_mode::pr:     return "pr";  case ag4_metal_mode::nd:     return "nd";
		case ag4_metal_mode::sm:     return "sm";  case ag4_metal_mode::eu:     return "eu";
		case ag4_metal_mode::gd:     return "gd";  case ag4_metal_mode::tb:     return "tb";
		case ag4_metal_mode::dy:     return "dy";  case ag4_metal_mode::ho:     return "ho";
		case ag4_metal_mode::er:     return "er";  case ag4_metal_mode::tm:     return "tm";
		case ag4_metal_mode::yb:     return "yb";  case ag4_metal_mode::lu:     return "lu";
		case ag4_metal_mode::se:     return "se";  case ag4_metal_mode::as_met: return "as";
		case ag4_metal_mode::ge:     return "ge";
		case ag4_metal_mode::b:      return "b";   case ag4_metal_mode::be:     return "be";
		case ag4_metal_mode::te:     return "te";  case ag4_metal_mode::po:     return "po";
		case ag4_metal_mode::at:     return "at";
		case ag4_metal_mode::rb:     return "rb";  case ag4_metal_mode::cs:     return "cs";
		case ag4_metal_mode::pm:     return "pm";  case ag4_metal_mode::ra:     return "ra";
		case ag4_metal_mode::ac:     return "ac";  case ag4_metal_mode::th:     return "th";
		case ag4_metal_mode::pa:     return "pa";  case ag4_metal_mode::u:      return "u";
		case ag4_metal_mode::np:     return "np";  case ag4_metal_mode::pu:     return "pu";
		case ag4_metal_mode::am:     return "am";  case ag4_metal_mode::cm:     return "cm";
		case ag4_metal_mode::bk:     return "bk";  case ag4_metal_mode::cf:     return "cf";
		case ag4_metal_mode::es:     return "es";  case ag4_metal_mode::fm:     return "fm";
		case ag4_metal_mode::fe2:    return "fe2"; case ag4_metal_mode::fe3:    return "fe3";
		case ag4_metal_mode::cu1:    return "cu1"; case ag4_metal_mode::cu2:    return "cu2";
		case ag4_metal_mode::cu2_jt: return "cu2_jt";
		case ag4_metal_mode::mn2:    return "mn2"; case ag4_metal_mode::mn3:    return "mn3";
		case ag4_metal_mode::mn3_jt: return "mn3_jt";
		case ag4_metal_mode::ni2:    return "ni2"; case ag4_metal_mode::ni3:    return "ni3";
		case ag4_metal_mode::co2:    return "co2"; case ag4_metal_mode::co3:    return "co3";
		case ag4_metal_mode::as3:    return "as3"; case ag4_metal_mode::as5:    return "as5";
		case ag4_metal_mode::sb3:    return "sb3"; case ag4_metal_mode::sb5:    return "sb5";
		case ag4_metal_mode::v4:     return "v4";  case ag4_metal_mode::v5:     return "v5";
		case ag4_metal_mode::mo4:    return "mo4"; case ag4_metal_mode::mo6:    return "mo6";
		case ag4_metal_mode::uo2:    return "uo2";
		case ag4_metal_mode::gd_dtpa:return "gd_dtpa";
		case ag4_metal_mode::gd_dota:return "gd_dota";
		case ag4_metal_mode::mg_aq:  return "mg_aq";
		case ag4_metal_mode::ca_aq:  return "ca_aq";
		case ag4_metal_mode::fe3_aq: return "fe3_aq";
		case ag4_metal_mode::mn2_aq: return "mn2_aq";
		case ag4_metal_mode::co2_aq: return "co2_aq";
		default: return "none";
	}
}

static ag4_metal_mode detect_pdbqt_metal_atom_type(const std::string& t) {
	static const std::pair<const char*, ag4_metal_mode> metal_map[] = {
		{"Zn",ag4_metal_mode::zn}, {"Mg",ag4_metal_mode::mg}, {"Ca",ag4_metal_mode::ca},
		{"Mn",ag4_metal_mode::mn}, {"Fe",ag4_metal_mode::fe}, {"Co",ag4_metal_mode::co},
		{"Ni",ag4_metal_mode::ni}, {"Cu",ag4_metal_mode::cu}, {"Pt",ag4_metal_mode::pt},
		{"Pd",ag4_metal_mode::pd}, {"Ru",ag4_metal_mode::ru}, {"Ir",ag4_metal_mode::ir},
		{"Au",ag4_metal_mode::au}, {"Cd",ag4_metal_mode::cd}, {"Hg",ag4_metal_mode::hg},
		{"Na",ag4_metal_mode::na_ion}, {"K",ag4_metal_mode::k_ion}, {"V",ag4_metal_mode::v},
		{"Cr1",ag4_metal_mode::cr}, {"Cr",ag4_metal_mode::cr}, {"Ti",ag4_metal_mode::ti},
		{"Sc",ag4_metal_mode::sc}, {"Y",ag4_metal_mode::y_ion}, {"Zr",ag4_metal_mode::zr},
		{"Nb",ag4_metal_mode::nb}, {"Hf",ag4_metal_mode::hf}, {"Ta",ag4_metal_mode::ta},
		{"Tg",ag4_metal_mode::tg_mode}, {"Mo",ag4_metal_mode::mo}, {"Rh",ag4_metal_mode::rh},
		{"Ag",ag4_metal_mode::ag}, {"Tc",ag4_metal_mode::tc}, {"Re",ag4_metal_mode::re},
		{"Os",ag4_metal_mode::os}, {"Ga",ag4_metal_mode::ga}, {"In",ag4_metal_mode::in_ion},
		{"Sn",ag4_metal_mode::sn}, {"Sb",ag4_metal_mode::sb}, {"Bi",ag4_metal_mode::bi},
		{"Tl",ag4_metal_mode::tl}, {"Pb",ag4_metal_mode::pb}, {"Li",ag4_metal_mode::li},
		{"Al",ag4_metal_mode::al}, {"Sr",ag4_metal_mode::sr}, {"Ba",ag4_metal_mode::ba},
		{"La",ag4_metal_mode::la}, {"Ce",ag4_metal_mode::ce}, {"Pr",ag4_metal_mode::pr},
		{"Nd",ag4_metal_mode::nd}, {"Sm",ag4_metal_mode::sm}, {"Eu",ag4_metal_mode::eu},
		{"Gd",ag4_metal_mode::gd}, {"Tb",ag4_metal_mode::tb}, {"Dy",ag4_metal_mode::dy},
		{"Ho",ag4_metal_mode::ho}, {"Er",ag4_metal_mode::er}, {"Tm",ag4_metal_mode::tm},
		{"Yb",ag4_metal_mode::yb}, {"Lu",ag4_metal_mode::lu}, {"Se",ag4_metal_mode::se},
		{"As",ag4_metal_mode::as_met}, {"Ge",ag4_metal_mode::ge}, {"B",ag4_metal_mode::b},
		{"Be",ag4_metal_mode::be}, {"Te",ag4_metal_mode::te}, {"Po",ag4_metal_mode::po},
		{"At",ag4_metal_mode::at}, {"Rb",ag4_metal_mode::rb}, {"Cs",ag4_metal_mode::cs},
		{"Pm",ag4_metal_mode::pm}, {"Ra",ag4_metal_mode::ra}, {"Ac",ag4_metal_mode::ac},
		{"Th",ag4_metal_mode::th}, {"Pa",ag4_metal_mode::pa}, {"U",ag4_metal_mode::u},
		{"Np",ag4_metal_mode::np}, {"Pu",ag4_metal_mode::pu}, {"Am",ag4_metal_mode::am},
		{"Cm",ag4_metal_mode::cm}, {"Bk",ag4_metal_mode::bk}, {"Cf",ag4_metal_mode::cf},
		{"E",ag4_metal_mode::es}, {"Fm",ag4_metal_mode::fm},
	};
	for (const auto& mp : metal_map)
		if (t == mp.first) return mp.second;
	return ag4_metal_mode::none;
}

static std::vector<ag4_metal_mode> detect_metal_modes_from_pdbqt(const std::string& path) {
	std::ifstream f(path);
	if (!f.is_open()) return {};
	std::set<ag4_metal_mode> found;
	std::string line;
	while (std::getline(f, line)) {
		if (line.size() < 6) continue;
		std::string rec = line.substr(0,6);
		if (rec != "ATOM  " && rec != "HETATM") continue;
		std::istringstream iss(line);
		std::string tok, last;
		while (iss >> tok) last = tok;
		ag4_metal_mode mode = detect_pdbqt_metal_atom_type(last);
		if (mode != ag4_metal_mode::none) found.insert(mode);
	}
	if (found.empty()) return {};
	if (found.size() == 1) {
		ag4_metal_mode refined = detect_metal_mode_from_pdbqt(path);
		if (refined != ag4_metal_mode::none) return {refined};
	}
	return std::vector<ag4_metal_mode>(found.begin(), found.end());
}

static std::string ag4_metal_modes_to_str(const std::vector<ag4_metal_mode>& modes) {
	std::ostringstream oss;
	bool first = true;
	for (ag4_metal_mode mode : modes) {
		if (mode == ag4_metal_mode::none) continue;
		if (!first) oss << ",";
		oss << ag4_metal_mode_to_str(mode);
		first = false;
	}
	return oss.str();
}

static void apply_auto_metal_modes(Vina& v, const std::vector<ag4_metal_mode>& modes) {
	v.clear_extra_metal_modes();
	v.set_zn_mode(false);
	bool first = true;
	for (ag4_metal_mode mode : modes) {
		if (mode == ag4_metal_mode::none) continue;
		if (mode == ag4_metal_mode::zn) v.set_zn_mode(true);
		if (first) {
			v.set_metal_mode(mode);
			first = false;
		} else {
			v.add_metal_mode(mode);
		}
	}
	if (first) v.set_metal_mode(ag4_metal_mode::none);
}

// O3: Max coordination number per metal mode
static int ag4_max_cn(ag4_metal_mode m) {
	switch (m) {
		case ag4_metal_mode::ca: case ag4_metal_mode::ca_aq: return 7;
		case ag4_metal_mode::k_ion: return 7;
		case ag4_metal_mode::au: case ag4_metal_mode::cu1: case ag4_metal_mode::ag:
		case ag4_metal_mode::hg: return 2;  // linear d10
		case ag4_metal_mode::pt: case ag4_metal_mode::pd: case ag4_metal_mode::ni: case ag4_metal_mode::ni2: case ag4_metal_mode::cu: case ag4_metal_mode::cu2: return 4;
		case ag4_metal_mode::cu2_jt: return 6;
		case ag4_metal_mode::mn3_jt: return 6;
		case ag4_metal_mode::gd: case ag4_metal_mode::gd_dtpa: case ag4_metal_mode::gd_dota:
		case ag4_metal_mode::la: case ag4_metal_mode::ce: case ag4_metal_mode::pr: case ag4_metal_mode::nd:
		case ag4_metal_mode::ba: return 9;
		case ag4_metal_mode::none: case ag4_metal_mode::zn: return 4;
		default: return 6;  // octahedral default
	}
}

// O6: Recommended weight_ad4_elec for high-charge-density metals
static const char* ag4_elec_weight_hint(ag4_metal_mode m) {
	switch (m) {
		case ag4_metal_mode::fe3: case ag4_metal_mode::mn3: case ag4_metal_mode::mn3_jt: return "0.0879";
		case ag4_metal_mode::v5:  case ag4_metal_mode::mo6: case ag4_metal_mode::as5: return "0.1172";
		case ag4_metal_mode::cu2: case ag4_metal_mode::cu2_jt: case ag4_metal_mode::fe2: return "0.0732";
		case ag4_metal_mode::mg:  case ag4_metal_mode::mg_aq: return "0.0732";
		case ag4_metal_mode::ti:  case ag4_metal_mode::zr: case ag4_metal_mode::hf: return "0.1172";
		case ag4_metal_mode::uo2: return "0.1465";
		default: return nullptr;
	}
}

// Dir4: Metal geometry check — print receptor metal-ligand distances vs expected r_eq
// Extended with O3 (CN residual), O7 (geoP score)
static void metal_geometry_check(const std::string& pdbqt_path, ag4_metal_mode mode) {
	std::ifstream f(pdbqt_path);
	if (!f.is_open()) {
		std::cerr << "[metal_geometry_check] Cannot open " << pdbqt_path << "\n";
		return;
	}
	static const std::set<std::string> metal_types = {
		"Zn","Mg","Ca","Mn","Fe","Co","Ni","Cu","Pt","Pd","Ru","Ir","Au","Cd","Hg",
		"Na","K","V","Cr1","Cr","Ti","Sc","Y","Zr","Nb","Hf","Ta","Tg","Mo",
		"Rh","Ag","Tc","Re","Os","Ga","In","Sn","Sb","Bi","Tl","Pb","Li","Al","Sr","Ba",
		"La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
		"Se","As","Ge","B","Be","Te","Po","At","Rb","Cs","Ra",
		"Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","E","Fm"
	};
	static const std::set<std::string> donor_types = {"OA","NA","N","SA","HD"};
	struct Atom { std::string atype, name; double x, y, z; };
	std::vector<Atom> metals, donors;
	std::string line;
	while (std::getline(f, line)) {
		if (line.size() < 60) continue;
		std::string rec = line.substr(0,6);
		if (rec != "ATOM  " && rec != "HETATM") continue;
		std::string name = (line.size()>=16) ? line.substr(12,4) : "    ";
		while (!name.empty() && name[0]==' ') name=name.substr(1);
		while (!name.empty() && name.back()==' ') name.pop_back();
		double x=0,y=0,z=0;
		try { x=std::stod(line.substr(30,8)); y=std::stod(line.substr(38,8)); z=std::stod(line.substr(46,8)); }
		catch(...) { continue; }
		std::string atype;
		if (line.size()>=79) { atype=line.substr(77,2); }
		else if (line.size()>=78) { atype=line.substr(77,1); }
		while (!atype.empty() && atype[0]==' ') atype=atype.substr(1);
		while (!atype.empty() && atype.back()==' ') atype.pop_back();
		if (metal_types.count(atype)) metals.push_back({atype,name,x,y,z});
		else if (donor_types.count(atype)) donors.push_back({atype,name,x,y,z});
	}
	if (metals.empty()) {
		std::cout << "[metal_geometry_check] No metal atoms found in " << pdbqt_path << "\n";
		return;
	}
	std::vector<ag4_nbp_override> overrides;
	ag4_apply_metal_mode(mode, overrides);
	std::cout << "\n[Metal Geometry Check]  " << pdbqt_path
	          << "  mode=" << ag4_metal_mode_to_str(mode) << "\n";
	std::cout << std::string(74, '-') << "\n";
	std::printf("%-6s %-4s  %-4s %-4s  %8s  %8s  %8s  %s\n",
	            "Metal","Type","Donor","Type","Dist(A)","Exp_r(A)","Delta","Quality");
	std::cout << std::string(74, '-') << "\n";
	const double GEOP_SIGMA = 0.20;  // O7: geoP Gaussian width (Å)
	for (const auto& m : metals) {
		double geoP_sum = 0.0;
		int    cn_actual = 0;
		for (const auto& d : donors) {
			double dx=m.x-d.x, dy=m.y-d.y, dz=m.z-d.z;
			double dist=std::sqrt(dx*dx+dy*dy+dz*dz);
			if (dist>3.5) continue;
			if (dist<2.8) cn_actual++;  // O3: count coordination partners
			double exp_r=-1.0;
			for (const auto& ov : overrides)
			if (ov.probe==d.atype && ov.receptor==m.atype) { exp_r=ov.r_eq; break; }
			if (exp_r>0) {
				double delta=dist-exp_r;
				const char* q = std::fabs(delta)<0.10 ? "ideal" :
				                std::fabs(delta)<0.20 ? "good"  :
				                std::fabs(delta)<0.35 ? "fair"  : "poor";
				std::printf("%-6s %-4s  %-4s %-4s  %8.3f  %8.3f  %+8.3f  %s\n",
				            m.name.c_str(),m.atype.c_str(),d.name.c_str(),d.atype.c_str(),
				            dist,exp_r,delta,q);
				// O7: geoP score contribution
				geoP_sum += std::exp(-(delta*delta)/(2.0*GEOP_SIGMA*GEOP_SIGMA));
			} else {
				std::printf("%-6s %-4s  %-4s %-4s  %8.3f  %8s  %8s  no_override\n",
				            m.name.c_str(),m.atype.c_str(),d.name.c_str(),d.atype.c_str(),
				            dist,"-","-");
			}
		}
		// O3: Coordination number residual vacancy
		int max_cn = ag4_max_cn(mode);
		int vacancies = (max_cn > cn_actual) ? (max_cn - cn_actual) : 0;
		std::printf("  >> %s (%s): CN_receptor=%d / max_CN=%d  vacancies=%d  geoP=%.3f\n",
		            m.name.c_str(), m.atype.c_str(),
		            cn_actual, max_cn, vacancies, geoP_sum);
	}
	std::cout << std::string(74, '-') << "\n";
}

// Parse a single metal mode token (e.g. "fe3", "gd_dtpa") → enum
static ag4_metal_mode parse_metal_mode_token(const std::string& tok) {
	if      (tok=="zn")       return ag4_metal_mode::zn;
	else if (tok=="mg")       return ag4_metal_mode::mg;
	else if (tok=="ca")       return ag4_metal_mode::ca;
	else if (tok=="mn")       return ag4_metal_mode::mn;
	else if (tok=="fe")       return ag4_metal_mode::fe;
	else if (tok=="co")       return ag4_metal_mode::co;
	else if (tok=="ni")       return ag4_metal_mode::ni;
	else if (tok=="cu")       return ag4_metal_mode::cu;
	else if (tok=="pt")       return ag4_metal_mode::pt;
	else if (tok=="pd")       return ag4_metal_mode::pd;
	else if (tok=="ru")       return ag4_metal_mode::ru;
	else if (tok=="ir")       return ag4_metal_mode::ir;
	else if (tok=="au")       return ag4_metal_mode::au;
	else if (tok=="cd")       return ag4_metal_mode::cd;
	else if (tok=="hg")       return ag4_metal_mode::hg;
	else if (tok=="na")       return ag4_metal_mode::na_ion;
	else if (tok=="k")        return ag4_metal_mode::k_ion;
	else if (tok=="v")        return ag4_metal_mode::v;
	else if (tok=="cr")       return ag4_metal_mode::cr;
	else if (tok=="ti")       return ag4_metal_mode::ti;
	else if (tok=="sc")       return ag4_metal_mode::sc;
	else if (tok=="y")        return ag4_metal_mode::y_ion;
	else if (tok=="zr")       return ag4_metal_mode::zr;
	else if (tok=="nb")       return ag4_metal_mode::nb;
	else if (tok=="hf")       return ag4_metal_mode::hf;
	else if (tok=="ta")       return ag4_metal_mode::ta;
	else if (tok=="w")        return ag4_metal_mode::tg_mode;
	else if (tok=="mo")       return ag4_metal_mode::mo;
	else if (tok=="rh")       return ag4_metal_mode::rh;
	else if (tok=="ag")       return ag4_metal_mode::ag;
	else if (tok=="tc")       return ag4_metal_mode::tc;
	else if (tok=="re")       return ag4_metal_mode::re;
	else if (tok=="os")       return ag4_metal_mode::os;
	else if (tok=="ga")       return ag4_metal_mode::ga;
	else if (tok=="in")       return ag4_metal_mode::in_ion;
	else if (tok=="sn")       return ag4_metal_mode::sn;
	else if (tok=="sb")       return ag4_metal_mode::sb;
	else if (tok=="bi")       return ag4_metal_mode::bi;
	else if (tok=="tl")       return ag4_metal_mode::tl;
	else if (tok=="pb")       return ag4_metal_mode::pb;
	else if (tok=="li")       return ag4_metal_mode::li;
	else if (tok=="al")       return ag4_metal_mode::al;
	else if (tok=="sr")       return ag4_metal_mode::sr;
	else if (tok=="ba")       return ag4_metal_mode::ba;
	else if (tok=="la")       return ag4_metal_mode::la;
	else if (tok=="ce")       return ag4_metal_mode::ce;
	else if (tok=="pr")       return ag4_metal_mode::pr;
	else if (tok=="nd")       return ag4_metal_mode::nd;
	else if (tok=="sm")       return ag4_metal_mode::sm;
	else if (tok=="eu")       return ag4_metal_mode::eu;
	else if (tok=="gd")       return ag4_metal_mode::gd;
	else if (tok=="tb")       return ag4_metal_mode::tb;
	else if (tok=="dy")       return ag4_metal_mode::dy;
	else if (tok=="ho")       return ag4_metal_mode::ho;
	else if (tok=="er")       return ag4_metal_mode::er;
	else if (tok=="tm")       return ag4_metal_mode::tm;
	else if (tok=="yb")       return ag4_metal_mode::yb;
	else if (tok=="lu")       return ag4_metal_mode::lu;
	else if (tok=="se")       return ag4_metal_mode::se;
	else if (tok=="as")       return ag4_metal_mode::as_met;
	else if (tok=="ge")       return ag4_metal_mode::ge;
	else if (tok=="b")        return ag4_metal_mode::b;
	else if (tok=="be")       return ag4_metal_mode::be;
	else if (tok=="te")       return ag4_metal_mode::te;
	else if (tok=="po")       return ag4_metal_mode::po;
	else if (tok=="at")       return ag4_metal_mode::at;
	else if (tok=="rb")       return ag4_metal_mode::rb;
	else if (tok=="cs")       return ag4_metal_mode::cs;
	else if (tok=="pm")       return ag4_metal_mode::pm;
	else if (tok=="ra")       return ag4_metal_mode::ra;
	else if (tok=="ac")       return ag4_metal_mode::ac;
	else if (tok=="th")       return ag4_metal_mode::th;
	else if (tok=="pa")       return ag4_metal_mode::pa;
	else if (tok=="u")        return ag4_metal_mode::u;
	else if (tok=="np")       return ag4_metal_mode::np;
	else if (tok=="pu")       return ag4_metal_mode::pu;
	else if (tok=="am")       return ag4_metal_mode::am;
	else if (tok=="cm")       return ag4_metal_mode::cm;
	else if (tok=="bk")       return ag4_metal_mode::bk;
	else if (tok=="cf")       return ag4_metal_mode::cf;
	else if (tok=="es")       return ag4_metal_mode::es;
	else if (tok=="fm")       return ag4_metal_mode::fm;
	// Oxidation-state variants
	else if (tok=="fe2")      return ag4_metal_mode::fe2;
	else if (tok=="fe3")      return ag4_metal_mode::fe3;
	else if (tok=="cu1")      return ag4_metal_mode::cu1;
	else if (tok=="cu2")      return ag4_metal_mode::cu2;
	else if (tok=="cu2_jt")   return ag4_metal_mode::cu2_jt;
	else if (tok=="mn2")      return ag4_metal_mode::mn2;
	else if (tok=="mn3")      return ag4_metal_mode::mn3;
	else if (tok=="mn3_jt")   return ag4_metal_mode::mn3_jt;
	else if (tok=="ni2")      return ag4_metal_mode::ni2;
	else if (tok=="ni3")      return ag4_metal_mode::ni3;
	else if (tok=="co2")      return ag4_metal_mode::co2;
	else if (tok=="co3")      return ag4_metal_mode::co3;
	else if (tok=="as3")      return ag4_metal_mode::as3;
	else if (tok=="as5")      return ag4_metal_mode::as5;
	else if (tok=="sb3")      return ag4_metal_mode::sb3;
	else if (tok=="sb5")      return ag4_metal_mode::sb5;
	else if (tok=="v4")       return ag4_metal_mode::v4;
	else if (tok=="v5")       return ag4_metal_mode::v5;
	else if (tok=="mo4")      return ag4_metal_mode::mo4;
	else if (tok=="mo6")      return ag4_metal_mode::mo6;
	else if (tok=="uo2")      return ag4_metal_mode::uo2;
	else if (tok=="gd_dtpa")  return ag4_metal_mode::gd_dtpa;
	else if (tok=="gd_dota")  return ag4_metal_mode::gd_dota;
	// Aqua variants (O4)
	else if (tok=="mg_aq")    return ag4_metal_mode::mg_aq;
	else if (tok=="ca_aq")    return ag4_metal_mode::ca_aq;
	else if (tok=="fe3_aq")   return ag4_metal_mode::fe3_aq;
	else if (tok=="mn2_aq")   return ag4_metal_mode::mn2_aq;
	else if (tok=="co2_aq")   return ag4_metal_mode::co2_aq;
	else {
		std::cerr << "ERROR: unknown --metal_mode token '" << tok
		          << "'. Examples: mg|ca|fe|fe2|fe3|cu|cu1|cu2|cu2_jt|mn|mn2|mn3|mn3_jt|co|co2|co3|ni|ni2|ni3|"
		             "pt|pd|ru|ir|au|zn|as3|as5|sb3|sb5|v4|v5|mo4|mo6|uo2|gd_dtpa|gd_dota\n"
		             "  Aqua variants: mg_aq|ca_aq|fe3_aq|mn2_aq|co2_aq\n"
		             "  Multi-metal: --metal_mode fe3,zn\n";
		exit(EXIT_FAILURE);
	}
}

struct usage_error : public std::runtime_error {
	usage_error(const std::string& message) : std::runtime_error(message) {}
};

struct options_occurrence {
	bool some;
	bool all;
	options_occurrence() : some(false), all(true) {} // convenience
	options_occurrence& operator+=(const options_occurrence& x) {
		some = some || x.some;
		all  = all  && x.all;
		return *this;
	}
};

options_occurrence get_occurrence(boost::program_options::variables_map& vm, boost::program_options::options_description& d) {
	options_occurrence tmp;
	VINA_FOR_IN(i, d.options())
		if(vm.count((*d.options()[i]).long_name()))
			tmp.some = true;
		else
			tmp.all = false;
	return tmp;
}

void check_occurrence(boost::program_options::variables_map& vm, boost::program_options::options_description& d) {
	VINA_FOR_IN(i, d.options()) {
		const std::string& str = (*d.options()[i]).long_name();
		if(!vm.count(str))
			std::cerr << "Required parameter --" << str << " is missing!\n";
	}
}

bool stdout_supports_color() {
	const char* no_color = std::getenv("NO_COLOR");
	if(no_color)
		return false;

	const char* clicolor_force = std::getenv("CLICOLOR_FORCE");
	if(clicolor_force && std::string(clicolor_force) != "0")
		return true;

#ifdef _WIN32
	if(!_isatty(_fileno(stdout)))
		return false;

	HANDLE h = GetStdHandle(STD_OUTPUT_HANDLE);
	if(h == INVALID_HANDLE_VALUE || h == NULL)
		return false;

	DWORD mode = 0;
	if(!GetConsoleMode(h, &mode)) {
		const char* ansicon = std::getenv("ANSICON");
		const char* wt_session = std::getenv("WT_SESSION");
		const char* conemu_ansi = std::getenv("ConEmuANSI");
		return ansicon || wt_session || (conemu_ansi && std::string(conemu_ansi) == "ON");
	}

	if(mode & ENABLE_VIRTUAL_TERMINAL_PROCESSING)
		return true;

	if(SetConsoleMode(h, mode | ENABLE_VIRTUAL_TERMINAL_PROCESSING))
		return true;

	const char* ansicon = std::getenv("ANSICON");
	const char* wt_session = std::getenv("WT_SESSION");
	const char* conemu_ansi = std::getenv("ConEmuANSI");
	return ansicon || wt_session || (conemu_ansi && std::string(conemu_ansi) == "ON");
#else
	if(!isatty(STDOUT_FILENO))
		return false;

	const char* term = std::getenv("TERM");
	if(!term)
		return false;
	if(std::string(term) == "dumb")
		return false;

	return true;
#endif
}

std::string make_lkina_ascii_banner(const std::string& version_string, bool use_color) {
	static const char* art[] = {
		"  _     _  ___ _             ",
		" | |   | |/ / (_)_ __   __ _ ",
		" | |   | ' /| | | '_ \\  / _` |",
		" | |___| . \\| | | | | | (_| |",
		" |_____|_|\\_\\_|_|_| |_|\\__,_|",
		nullptr
	};
	std::string result;
	if (!use_color) {
		for (int i = 0; art[i]; ++i)
			result += std::string(art[i]) + "\n";
		result += " " + version_string + " -- Metal-Enhanced AD4 Docking Engine\n";
		return result;
	}
	// Horizontal gradient: cyan-blue(91,200,240) -> purple(160,122,224) -> pink(232,120,180)
	int max_w = 0;
	for (int i = 0; art[i]; ++i) {
		int n = (int)std::strlen(art[i]);
		if (n > max_w) max_w = n;
	}
	for (int i = 0; art[i]; ++i) {
		std::string line(art[i]);
		int len = (int)line.size();
		for (int j = 0; j < len; ++j) {
			float t = (max_w > 1) ? (float)j / (max_w - 1) : 0.0f;
			int r, g, b;
			if (t < 0.5f) {
				float s = t * 2.0f;
				r = 91  + (int)((160 - 91 ) * s);
				g = 200 + (int)((122 - 200) * s);
				b = 240 + (int)((224 - 240) * s);
			} else {
				float s = (t - 0.5f) * 2.0f;
				r = 160 + (int)((232 - 160) * s);
				g = 122 + (int)((120 - 122) * s);
				b = 224 + (int)((180 - 224) * s);
			}
			r = std::max(0, std::min(255, r));
			g = std::max(0, std::min(255, g));
			b = std::max(0, std::min(255, b));
			result += "\033[38;2;" + std::to_string(r) + ";" +
			          std::to_string(g) + ";" + std::to_string(b) + "m";
			result += line[j];
		}
		result += "\033[0m\n";
	}
	result += "\033[1;97m " + version_string + " -- Metal-Enhanced AD4 Docking Engine\033[0m\n";
	return result;
}

int main(int argc, char* argv[]) {
	using namespace boost::program_options;
	const std::string git_version = VERSION;
	const std::string version_string = "LKina " + git_version;
	const bool use_color_banner = stdout_supports_color();
	const std::string error_message = "\n\n\
Please report bugs through the Issue Tracker on GitHub \n\
(https://github.com/ccsb-scripps/LKina/issues), so\n\
that this problem can be resolved. The reproducibility of the\n\
error may be vital, so please remember to include the following in\n\
your problem report:\n\
* the EXACT error message,\n\
* your version of the program,\n\
* the type of computer system you are running it on,\n\
* all command line options,\n\
* configuration file (if used),\n\
* ligand file as PDBQT,\n\
* receptor file as PDBQT,\n\
* flexible side chains file as PDBQT (if used),\n\
* output file as PDBQT (if any),\n\
* input (if possible),\n\
* random seed the program used (this is printed when the program starts).\n\
\n\
Thank you!\n";

	const std::string lkina_ascii_banner = make_lkina_ascii_banner(version_string, use_color_banner);

	const std::string lkina_metal_list =
		"\nLKina Supported AD_TYPE Tokens (117 total: 113 chemical elements + 4 pseudoatoms TZ/SQ/MH/JT, AD_TYPE_SIZE=117):\n"
		"  Standard AutoDock4/Vina types (index 0-30):\n"
		"    Organic : C   A   N   O   P   S   H   F   I   NA  OA  SA  HD\n"
		"    Metals  : Mg  Mn  Zn  Ca  Fe\n"
		"    Other   : Cl  Br  Si  At  G0  G1  G2  G3  CG0 CG1 CG2 CG3 W\n"
		"  LKina Extended elements (index 31-112):\n"
		"    He  Li  Be  B   Ne  Na  Al  K   Sc  Ti  V   Co  Ni  Cu  Ga\n"
		"    Ge  As  Kr  Rb  Sr  Y   Zr  Nb  Mo  Tc  Ru  Rh  Pd  Ag  Cd\n"
		"    In  Sn  Sb  Te  Xe  Cs  Ba  La  Ce  Pr  Nd  Pm  Sm  Eu  Gd\n"
		"    Tb  Dy  Ho  Er  Tm  Yb  Lu  Hf  Ta  Re  Os  Ir  Pt  Au  Hg\n"
		"    Tl  Pb  Bi  Po  Rn  Fr  Ra  Ac  Th  Pa  U   Np  Pu  Am  Cm\n"
		"    Bk  Cf  E   Fm  Cr1 Tg  Se\n\n";

	const std::string cite_message = "\
#################################################################\n\
#                    LKina Docking Engine                       #\n\
#      Metal-Enhanced AutoDock4 Force Field  (AD4 + inline      #\n\
#      AutoGrid, no external autogrid4 required)                #\n\
#                                                               #\n\
# Based on AutoDock Vina v1.2.7 (Eberhardt et al. 2021)        #\n\
# Extended with inline AG4 grid engine for full-metal docking   #\n\
# Supports 113 atom types incl. 80+ metals (AD4 force field)    #\n\
#                                                               #\n\
# Usage: LKina --scoring LKDock --generate_maps --receptor r.pdbqt #\n\
#              --center_x X --center_y Y --center_z Z           #\n\
#              --size_x N  --size_y N  --size_z N               #\n\
#              --ligand metal.pdbqt --out out.pdbqt             #\n\
#################################################################\n";

	try {
		std::string rigid_name;
		std::string flex_name;
		std::string config_name;
		std::string out_name;
		std::string out_dir;
		std::string out_maps;
		std::vector<std::string> ligand_names;
		std::vector<std::string> batch_ligand_names;
		std::string maps;
		std::string sf_name = "vina";
		double center_x;
		double center_y;
		double center_z;
		double size_x;
		double size_y;
		double size_z;
		int cpu = 0;
		int seed = 0;
		int exhaustiveness = 8;
		int max_evals = 0;
		int verbosity = 1;
		int num_modes = 9;
		double min_rmsd = 1.0;
		double energy_range = 3.0;
		double grid_spacing = 0.375;
		double buffer_size = 4;
		double unbound_energy = NAN;

		// autodock4.2 weights
		double weight_ad4_vdw   = 0.1662;
		double weight_ad4_hb    = 0.1209;
		double weight_ad4_elec  = 0.1406;
		double weight_ad4_dsolv = 0.1322;
		double weight_ad4_rot   = 0.2983;

		// vina weights
		double weight_gauss1      = -0.035579;
		double weight_gauss2      = -0.005156;
		double weight_repulsion   =  0.840245;
		double weight_hydrophobic = -0.035069;
		double weight_hydrogen    = -0.587439;
		double weight_rot         =  0.05846;

		// vinardo weights
		double weight_vinardo_gauss1 = -0.045;
		double weight_vinardo_repulsion = 0.8;
		double weight_vinardo_hydrophobic = -0.035;
		double weight_vinardo_hydrogen = -0.600;

		// macrocycle closure
		double weight_glue        = 50.000000; // linear attraction

		// reactive covalent docking
		std::string reactive_mode_str;
		std::string reactive_preset_name;
		std::string reactive_rec_atom;
		std::string reactive_lig_atom;
		double reactive_bond_length       = 0.0;
		double reactive_attractor_width   = 1.5;
		double reactive_attractor_strength = 8.0;
		double reactive_angle_strength     = 4.0;
		double reactive_gradcheck_eps      = 1e-4;
		std::string reactive_frame_atom;
		double reactive_hybrid_vdw_scale  = 0.0;
		double reactive_target_angle      = 180.0;
		double reactive_angle_width       = 0.0;
		std::string reactive_lig_frame_atom;
		bool   reactive_debug             = false;
		bool   reactive_debug_energy      = false;
		bool   reactive_gradcheck         = false;
		bool   reactive_programmatic_width = false;
		bool   reactive_programmatic_strength = false;

		bool score_only = false;
		bool local_only = false;
		bool no_refine = false;
		bool force_even_voxels = false;
		bool randomize_only = false;
		bool generate_maps = false;
		bool zn_mode = false;
		std::string metal_mode_str;
		bool metal_geometry_check_flag = false;
		bool metal_bias = false;
		double metal_bias_strength = 2.0;  // O5: Gaussian well depth (kcal/mol)
		double metal_bias_width    = 1.5;  // O5: Gaussian well sigma (Å)
		double metal_soft_weight   = 0.0;  // soft rerank gradient during search (0=off)
		double ligand_metal_geometry_weight = 0.0;
		bool no_auto_metal = false;  // M5 fix: opt-out of automatic metal-mode detection
		bool help = false;
		bool help_advanced = false;
		bool version = false; // FIXME
		bool autobox = false;
		variables_map vm;

		positional_options_description positional; // remains empty

		options_description inputs("Input");
		inputs.add_options()
			("receptor", value<std::string>(&rigid_name), "rigid part of the receptor (PDBQT)")
			("flex", value<std::string>(&flex_name), "flexible side chains, if any (PDBQT)")
			("ligand", value< std::vector<std::string> >(&ligand_names)->multitoken(), "ligand (PDBQT)")
			("batch", value< std::vector<std::string> >(&batch_ligand_names)->multitoken(), "batch directory or ligands (PDBQT)")
			("scoring", value<std::string>(&sf_name)->default_value(sf_name), "scoring function (LKDock/ad4, vina or vinardo)")
		;
		//options_description search_area("Search area (required, except with --score_only)");
		options_description search_area("Search space (required)");
		search_area.add_options()
			("maps", value<std::string>(&maps), "affinity maps for the autodock4.2 (ad4) or vina scoring function")
			("center_x", value<double>(&center_x), "X coordinate of the center (Angstrom)")
			("center_y", value<double>(&center_y), "Y coordinate of the center (Angstrom)")
			("center_z", value<double>(&center_z), "Z coordinate of the center (Angstrom)")
			("size_x", value<double>(&size_x), "size in the X dimension (Angstrom)")
			("size_y", value<double>(&size_y), "size in the Y dimension (Angstrom)")
			("size_z", value<double>(&size_z), "size in the Z dimension (Angstrom)")
			("autobox", bool_switch(&autobox), "set maps dimensions based on input ligand(s) (for --score_only and --local_only)")
		;
		//options_description outputs("Output prefixes (optional - by default, input names are stripped of .pdbqt\nare used as prefixes. _001.pdbqt, _002.pdbqt, etc. are appended to the prefixes to produce the output names");
		options_description outputs("Output (optional)");
		outputs.add_options()
			("out", value<std::string>(&out_name), "output models (PDBQT), the default is chosen based on the ligand file name")
			("dir", value<std::string>(&out_dir), "output directory for batch mode")
			("write_maps", value<std::string>(&out_maps), "output filename (directory + prefix name) for maps. Option --force_even_voxels may be needed to comply with .map format")
		;
		options_description advanced("Advanced options (see the manual)");
		advanced.add_options()
			("score_only",     bool_switch(&score_only),     "score only - search space can be omitted")
			("local_only",     bool_switch(&local_only),     "do local search only")
			("unbound_energy", value<double>(&unbound_energy)->default_value(unbound_energy), "Explicitly set the Unbound System's Energy for --score_only jobs")
			("no_refine", bool_switch(&no_refine),  "when --receptor is provided, do not use explicit receptor atoms (instead of precalculated grids) for: (1) local optimization and scoring after docking, (2) --local_only jobs, and (3) --score_only jobs")
			("force_even_voxels", bool_switch(&force_even_voxels),  "calculated grid maps will have an even number of voxels (intervals) in each dimension (odd number of grid points)")
			("randomize_only", bool_switch(&randomize_only), "randomize input, attempting to avoid clashes")
			("generate_maps", bool_switch(&generate_maps), "(AD4 only) compute affinity maps inline from --receptor + grid dimensions, no external autogrid4 needed")
			("zn_mode", bool_switch(&zn_mode), "(AD4 only) legacy alias for --metal_mode zn; enables AutoDock4Zn Zn coordination potentials")
			("metal_geometry_check", bool_switch(&metal_geometry_check_flag), "(AD4 only) print metal-ligand distances vs expected r_eq for the active metal_mode; useful for validating receptor preparation")
			("metal_bias", bool_switch(&metal_bias), "(AD4+generate_maps) auto-add a soft Gaussian attractor toward the receptor metal centre (MBD-style coordinate bias; O5)")
			("metal_bias_strength", value<double>(&metal_bias_strength)->default_value(2.0), "O5: metal bias Gaussian well depth (kcal/mol, default 2.0)")
			("metal_bias_width",    value<double>(&metal_bias_width)->default_value(1.5),  "O5: metal bias Gaussian sigma (Å, default 1.5)")
			("metal_soft_weight",   value<double>(&metal_soft_weight)->default_value(0.0), "(AD4+metal_mode) weight for metal geometry soft-constraint gradient during search (0=off; 0.1-0.5 recommended for metal-protein docking)")
			("no_auto_metal", bool_switch(&no_auto_metal), "(AD4) disable automatic metal-mode detection from receptor/ligand PDBQT atom types (keep pure AD4 scoring); metal_mode energies sit on a different scale than standard AD4/vina and should be enabled deliberately")
			("ligand_metal_geometry_weight", value<double>(&ligand_metal_geometry_weight)->default_value(0.0), "(AD4) optional rerank weight for ligand-side metallocomplex geometry QC: Pt/Pd square-planar and Ru/Os/Re octahedral (0=report only)")
			("metal_mode", value<std::string>(&metal_mode_str)->default_value(""), "(AD4 only) metal coordination mode — applies literature nbp_r_eps overrides; zn activates AutoDock4Zn. Biological: zn|mg|ca|mn|fe|co|ni|cu. Medicinal: pt|pd|ru|ir|au|rh|ag|tc|re|os. Toxicology: cd|hg|tl|pb|sb|bi|as. s-block: na|k|li|al|sr|ba. Early TM: v|cr|ti|sc|y|zr|nb|hf|ta|w|mo. Lanthanides: la|ce|pr|nd|sm|eu|gd|tb|dy|ho|er|tm|yb|lu. Metalloids: se|ge. Post-transition: ga|in|sn")
			("reactive_mode", value<std::string>(&reactive_mode_str)->default_value(""), "(AD4+generate_maps only) reactive covalent mode: distance or hybrid")
			("reactive_preset", value<std::string>(&reactive_preset_name)->default_value(""), "reaction-type preset: cys_michael|cys_sn2|ser_covalent|lys_targeting|boronic_acid|tyr_covalent (sets default bond/angle/strength params; individual flags override)")
			("reactive_rec_atom", value<std::string>(&reactive_rec_atom)->default_value(""), "receptor reactive atom: \"chain:resnum:name\", \"chain:resnum:icode:name\", or \"x,y,z\"")
			("reactive_lig_atom", value<std::string>(&reactive_lig_atom)->default_value(""), "ligand reactive atom: 1-based index, serial:N, name:ATOM, or atom name")
			("reactive_bond_length", value<double>(&reactive_bond_length)->default_value(0.0), "ideal covalent bond length (0=Gaussian centred at anchor)")
			("reactive_attractor_width", value<double>(&reactive_attractor_width)->default_value(1.5), "Gaussian well sigma (Angstrom)")
			("reactive_attractor_strength", value<double>(&reactive_attractor_strength)->default_value(8.0), "Gaussian well depth (kcal/mol)")
			("reactive_frame_atom", value<std::string>(&reactive_frame_atom)->default_value(""), "receptor frame atom for angular constraint: \"chain:resnum:name\" or \"x,y,z\"")
			("reactive_angle_strength", value<double>(&reactive_angle_strength)->default_value(4.0), "angular constraint K (kcal/mol), used with --reactive_frame_atom")
			("reactive_debug", bool_switch(&reactive_debug), "print reactive docking debug info")
			("reactive_two_step", value<bool>()->default_value(false)->implicit_value(true), "(AD4 only) C3 two-step strategy: phase-1 MC without reactive energy, phase-2 quasi-newton refinement on poses within --reactive_presample_dist of anchor")
			("reactive_presample_dist", value<double>()->default_value(10.0), "(AD4 only) two-step phase-1 distance filter: only refine poses within this Angstrom distance of the reactive anchor (default 10.0)")
			("reactive_hybrid_vdw_scale", value<double>(&reactive_hybrid_vdw_scale)->default_value(0.0), "hybrid mode VdW+HB grid scale (0=suppress, 1=full VdW+HB kept)")
			("reactive_target_angle", value<double>(&reactive_target_angle)->default_value(180.0), "ideal attack angle at reactive atom (degrees, default 180)")
			("reactive_angle_width", value<double>(&reactive_angle_width)->default_value(0.0), "flat-bottom half-width (degrees): 0=pure harmonic, >0=no penalty within N degrees of target")
			("reactive_lig_frame_atom", value<std::string>(&reactive_lig_frame_atom)->default_value(""), "ligand-side frame atom spec for approach-angle constraint at the reactive ligand atom")
			("reactive_debug_energy", bool_switch(&reactive_debug_energy), "print reactive distance/angle energy breakdown in score output")
			("reactive_gradcheck", bool_switch(&reactive_gradcheck), "run finite-difference gradient check for reactive term during score output")
			("reactive_gradcheck_eps", value<double>(&reactive_gradcheck_eps)->default_value(1e-4), "finite-difference step size (Angstrom) for --reactive_gradcheck")
			("reactive_weak_attractor", value<bool>()->default_value(false)->implicit_value(true), "(C3b) two-step phase 1 uses a broad/weak Gaussian bias toward anchor (σ×3, ε×0.15) instead of no reactive energy — closer to Goullieux 2023 attracting-cavities design")

			("weight_gauss1", value<double>(&weight_gauss1)->default_value(weight_gauss1),                "gauss_1 weight")
			("weight_gauss2", value<double>(&weight_gauss2)->default_value(weight_gauss2),                "gauss_2 weight")
			("weight_repulsion", value<double>(&weight_repulsion)->default_value(weight_repulsion),       "repulsion weight")
			("weight_hydrophobic", value<double>(&weight_hydrophobic)->default_value(weight_hydrophobic), "hydrophobic weight")
			("weight_hydrogen", value<double>(&weight_hydrogen)->default_value(weight_hydrogen),          "Hydrogen bond weight")
			("weight_rot", value<double>(&weight_rot)->default_value(weight_rot),                         "N_rot weight")

			("weight_vinardo_gauss1", value<double>(&weight_vinardo_gauss1)->default_value(weight_vinardo_gauss1), "Vinardo gauss_1 weight")
			("weight_vinardo_repulsion", value<double>(&weight_vinardo_repulsion)->default_value(weight_vinardo_repulsion), "Vinardo repulsion weight")
			("weight_vinardo_hydrophobic", value<double>(&weight_vinardo_hydrophobic)->default_value(weight_vinardo_hydrophobic), "Vinardo hydrophobic weight")
			("weight_vinardo_hydrogen", value<double>(&weight_vinardo_hydrogen)->default_value(weight_vinardo_hydrogen), "Vinardo Hydrogen bond weight")
			("weight_vinardo_rot", value<double>(&weight_rot)->default_value(weight_rot), "Vinardo N_rot weight")

			("weight_ad4_vdw", value<double>(&weight_ad4_vdw)->default_value(weight_ad4_vdw), "ad4_vdw weight")
			("weight_ad4_hb", value<double>(&weight_ad4_hb)->default_value(weight_ad4_hb), "ad4_hb weight")
			("weight_ad4_elec", value<double>(&weight_ad4_elec)->default_value(weight_ad4_elec), "ad4_elec weight")
			("weight_ad4_dsolv", value<double>(&weight_ad4_dsolv)->default_value(weight_ad4_dsolv), "ad4_dsolv weight")
			("weight_ad4_rot", value<double>(&weight_ad4_rot)->default_value(weight_ad4_rot), "ad4_rot weight")

			("weight_glue", value<double>(&weight_glue)->default_value(weight_glue),                      "macrocycle glue weight")
		;
		options_description misc("Misc (optional)");
		misc.add_options()
			("cpu", value<int>(&cpu)->default_value(0), "the number of CPUs to use (the default is to try to detect the number of CPUs or, failing that, use 1)")
			("seed", value<int>(&seed)->default_value(0), "explicit random seed")
			("exhaustiveness", value<int>(&exhaustiveness)->default_value(8), "exhaustiveness of the global search (roughly proportional to time): 1+")
			("max_evals", value<int>(&max_evals)->default_value(0), "number of evaluations in each MC run (if zero, which is the default, the number of MC steps is based on heuristics)")
			("num_modes", value<int>(&num_modes)->default_value(9), "maximum number of binding modes to generate")
			("min_rmsd", value<double>(&min_rmsd)->default_value(1.0), "minimum RMSD between output poses")
			("energy_range", value<double>(&energy_range)->default_value(3.0), "maximum energy difference between the best binding mode and the worst one displayed (kcal/mol)")
			("spacing", value<double>(&grid_spacing)->default_value(0.375), "grid spacing (Angstrom)")
			("verbosity", value<int>(&verbosity)->default_value(1), "verbosity (0=no output, 1=normal, 2=verbose)")
		;
		options_description config("Configuration file (optional)");
		config.add_options()
			("config", value<std::string>(&config_name), "the above options can be put here")
		;
		options_description info("Information (optional)");
		info.add_options()
			("help",          bool_switch(&help), "display usage summary")
			("help_advanced", bool_switch(&help_advanced), "display usage summary with advanced options")
			("version",       bool_switch(&version), "display program version")
		;
		options_description desc, desc_config, desc_simple;
		desc       .add(inputs).add(search_area).add(outputs).add(advanced).add(misc).add(config).add(info);
		desc_config.add(inputs).add(search_area).add(outputs).add(advanced).add(misc);
		desc_simple.add(inputs).add(search_area).add(outputs).add(misc).add(config).add(info);

		std::cout << lkina_ascii_banner << '\n';
		try {
			//store(parse_command_line(argc, argv, desc, command_line_style::default_style ^ command_line_style::allow_guessing), vm);
			store(command_line_parser(argc, argv)
				.options(desc)
				.style(command_line_style::default_style ^ command_line_style::allow_guessing)
				.positional(positional)
				.run(),
				vm);
			notify(vm);
		} catch(boost::program_options::error& e) {
			std::cerr << "Command line parse error: " << e.what() << '\n' << "\nCorrect usage:\n" << desc_simple << '\n';
			return 1;
		}

		if (vm.count("config")) {
			try {
				path name = make_path(config_name);
				ifile config_stream(name);
				store(parse_config_file(config_stream, desc_config), vm);
				notify(vm);
			}
			catch(boost::program_options::error& e) {
				std::cerr << "Configuration file parse error: " << e.what() << '\n' << "\nCorrect usage:\n" << desc_simple << '\n';
				return 1;
			}
		}

		auto option_explicit = [&](const char* name) -> bool {
			auto it = vm.find(name);
			return it != vm.end() && !it->second.defaulted();
		};

		if (help) {
			std::cout << desc_simple << '\n';
			std::cout << lkina_metal_list;
			return 0;
		}

		if (help_advanced) {
			std::cout << desc << '\n';
			std::cout << lkina_metal_list;
			return 0;
		}

		if (version) {
			return 0;
		}

		if (verbosity > 0) {
			std::cout << cite_message << '\n';
		}

		// C2: reactive docking with prebuilt maps needs --receptor for anchor resolution
		bool reactive_requested_early = (!reactive_mode_str.empty() && reactive_mode_str != "off")
		                              || !reactive_preset_name.empty();
		if (vm.count("receptor") && vm.count("maps")) {
			if (!reactive_requested_early) {
				std::cerr << "ERROR: Cannot specify both receptor and affinity maps at the same time, --flex argument is allowed with receptor or maps.\n";
				exit(EXIT_FAILURE);
			}
		}

		// Accept "LKDock" as an alias for the AD4 scoring function
		if (sf_name.compare("LKDock") == 0) sf_name = "ad4";
		if (reactive_requested_early && sf_name.compare("ad4") != 0) {
			std::cerr << "ERROR: --reactive_mode/--reactive_preset is only supported with --scoring ad4/LKDock.\n";
			exit(EXIT_FAILURE);
		}

		if (sf_name.compare("vina") == 0 || sf_name.compare("vinardo") == 0) {
			if (!vm.count("receptor") && !vm.count("maps")) {
				std::cerr << desc_simple << "ERROR: The receptor or affinity maps must be specified.\n";
				exit(EXIT_FAILURE);
			}
		} else if (sf_name.compare("ad4") == 0) {
			if (generate_maps) {
				// Inline map generation: receptor required, maps not required
				if (!vm.count("receptor")) {
					std::cerr << "ERROR: --generate_maps requires --receptor.\n";
					exit(EXIT_FAILURE);
				}
				if (!vm.count("center_x") || !vm.count("center_y") || !vm.count("center_z") ||
				    !vm.count("size_x")   || !vm.count("size_y")   || !vm.count("size_z")) {
					std::cerr << "ERROR: --generate_maps requires grid center and size (--center_x/y/z, --size_x/y/z).\n";
					exit(EXIT_FAILURE);
				}
			} else {
				if (vm.count("receptor") && !reactive_requested_early) {
					std::cerr << "ERROR: No receptor allowed, only --flex argument with the AD4 scoring function.\n";
					std::cerr << "       Use --generate_maps to compute maps inline from a receptor, or --reactive_preset/--reactive_mode with --receptor for reactive docking.\n";
					exit(EXIT_FAILURE);
				}
				if (!vm.count("maps")) {
					std::cerr << desc_simple << "\n\nERROR: Affinity maps are missing.\n";
					exit(EXIT_FAILURE);
				}
			}
		} else {
			std::cerr << desc_simple << "Scoring function " << sf_name << " unknown.\n";
			exit(EXIT_FAILURE);
		}

		if (!vm.count("ligand") && !vm.count("batch")) {
			std::cerr << desc_simple << "\n\nERROR: Missing ligand(s).\n";
			exit(EXIT_FAILURE);
		} else if (vm.count("ligand") && vm.count("batch")) {
			std::cerr << desc_simple << "\n\nERROR: Can't use both --ligand and --batch arguments simultaneously.\n";
			exit(EXIT_FAILURE);
		} else if (vm.count("batch") && !vm.count("dir")) {
			std::cerr << desc_simple << "\n\nERROR: Need to specify an output directory for batch mode.\n";
			exit(EXIT_FAILURE);
		} else if (vm.count("dir")) {
			if (!is_directory(out_dir)) {
				std::cerr << "ERROR: Directory " << out_dir << " does not exist.\n";
				exit(EXIT_FAILURE);
			}
		} else if (vm.count("ligand") && vm.count("dir")) {
			std::cerr << "WARNING: In ligand mode, --dir argument is ignored.\n";
		}

		if (vm.count("batch") && batch_ligand_names.size() == 1 && is_directory(batch_ligand_names[0])) {
			std::string in_dir = batch_ligand_names[0];
			batch_ligand_names.clear();
			for (const auto& entry : boost::filesystem::directory_iterator(in_dir)) {
				if (entry.path().extension() == ".pdbqt") {
					batch_ligand_names.push_back(entry.path().string());
				}
			}
		}

		if (!score_only) {
			if (!vm.count("out") && ligand_names.size() == 1) {
				out_name = default_output(ligand_names[0]);
				std::cout << "Output will be " << out_name << '\n';
			} else if (!vm.count("out") && ligand_names.size() >= 1) {
				std::cerr << desc_simple << "\n\nERROR: Output name must be defined when docking simultaneously multiple ligands.\n";
				exit(EXIT_FAILURE);
			}
		}
		if (ligand_metal_geometry_weight < 0.0) {
			std::cerr << "ERROR: --ligand_metal_geometry_weight must be >= 0.\n";
			exit(EXIT_FAILURE);
		}

		if (verbosity > 0) {
			std::cout << "Scoring function : " << (sf_name == "ad4" ? "LKDock (AD4)" : sf_name) << "\n";
			if (vm.count("receptor"))
				std::cout << "Rigid receptor: " << rigid_name << "\n";
			if (vm.count("flex"))
				std::cout << "Flex receptor: " << flex_name << "\n";
			if (ligand_names.size() == 1) {
				std::cout << "Ligand: " << ligand_names[0] << "\n";
			} else if (ligand_names.size() > 1) {
				std::cout << "Ligands:\n";
				VINA_RANGE(i, 0, ligand_names.size()) {
					std::cout << "  - " << ligand_names[i] << "\n";
				}
			} else if (batch_ligand_names.size() > 1) {
				std::cout << "Ligands (batch mode): " << batch_ligand_names.size() << " molecules\n";
			}
			if (!vm.count("maps") && !autobox) {
				std::cout << "Grid center: X " << center_x << " Y " << center_y << " Z " << center_z << "\n";
				std::cout << "Grid size  : X " << size_x << " Y " << size_y << " Z " << size_z << "\n";
				std::cout << "Grid space : " << grid_spacing << "\n";
			} else if (autobox) {
				std::cout << "Grid center: ligand center (autobox)\n";
				std::cout << "Grid size  : ligand size + " << buffer_size << " A in each dimension (autobox)\n";
				std::cout << "Grid space : " << grid_spacing << "\n";
			}
			std::cout << "Exhaustiveness: " << exhaustiveness << "\n";
			std::cout << "CPU: " << cpu << "\n";
			if (!vm.count("seed"))
				std::cout << "Seed: " << seed << "\n";
			std::cout << "Verbosity: " << verbosity << "\n";
			std::cout << "\n";
		}

		Vina v(sf_name, cpu, seed, verbosity, no_refine);
		v.set_ligand_metal_geometry_weight((float)ligand_metal_geometry_weight);
		bool ad4_metal_mode_configured = false;
		bool ad4_batch_ligand_metal_auto = false;

		// For AD4 + --maps: set_receptor errors on rigid receptor, but for reactive docking
		// we allow --receptor to be stored as the anchor-resolution path only (no model load).
		bool is_ad4 = (sf_name.compare("vina") != 0 && sf_name.compare("vinardo") != 0);
		if ((vm.count("receptor") || vm.count("flex")) && !generate_maps) {
			if (is_ad4 && vm.count("receptor") && !vm.count("flex")) {
				// AD4 + prebuilt maps: store receptor path for reactive anchor resolution only
				v.set_ad4_receptor_path(rigid_name);
			} else {
				v.set_receptor(rigid_name, flex_name);
			}
		}

		// Technically we don't have to initialize weights,
		// because they are initialized during the Vina object creation with the default weights
		// but we still do it in case the user decided to change them
		if (sf_name.compare("vina") == 0) {
			v.set_vina_weights(weight_gauss1, weight_gauss2, weight_repulsion,
							   weight_hydrophobic, weight_hydrogen, weight_glue, weight_rot);
		} else if (sf_name.compare("vinardo") == 0) {
			v.set_vinardo_weights(weight_vinardo_gauss1, weight_vinardo_repulsion,
								  weight_vinardo_hydrophobic, weight_vinardo_hydrogen, weight_glue, weight_rot);
		} else {
			v.set_ad4_weights(weight_ad4_vdw, weight_ad4_hb, weight_ad4_elec,
							  weight_ad4_dsolv, weight_glue, weight_ad4_rot);
			if (zn_mode) {
				v.set_zn_mode(true);
				if (metal_mode_str.empty()) {
					v.set_metal_mode(ag4_metal_mode::zn);
					ad4_metal_mode_configured = true;
				}
				if (verbosity > 0) std::cout << "Zn mode         : enabled (AD4Zn coordination potentials)\n";
			}
			// M5: auto-detect metal from receptor PDBQT when --metal_mode not given
			// (--no_auto_metal disables; warning printed to stderr at any verbosity so
			//  users are never silently moved onto the metal-mode energy scale)
			if (metal_mode_str.empty() && !ad4_metal_mode_configured && !no_auto_metal && !rigid_name.empty()) {
				ag4_metal_mode detected = detect_metal_mode_from_pdbqt(rigid_name);
				if (detected != ag4_metal_mode::none) {
					v.set_metal_mode(detected);
					if (detected == ag4_metal_mode::zn) v.set_zn_mode(true);
					ad4_metal_mode_configured = true;
					if (metal_soft_weight > 0.0)
						v.set_metal_soft_weight((float)metal_soft_weight);
					std::cerr << "[LKina] Auto-detected " << ag4_metal_mode_to_str(detected)
					          << " metal in receptor -> metal_mode applied automatically. "
					          << "Metal-mode absolute energies sit on a different scale than "
					          << "standard AD4/vina; use --no_auto_metal to force pure AD4 scoring.\n";
					if (verbosity > 0)
						std::cout << "Auto-detected   : " << ag4_metal_mode_to_str(detected)
						          << " metal in receptor → metal_mode applied automatically\n";
				}
			}
			if (metal_mode_str.empty() && !ad4_metal_mode_configured && !no_auto_metal && ligand_names.size() == 1) {
				std::vector<ag4_metal_mode> detected_modes = detect_metal_modes_from_pdbqt(ligand_names[0]);
				if (!detected_modes.empty()) {
					apply_auto_metal_modes(v, detected_modes);
					ad4_metal_mode_configured = true;
					if (metal_soft_weight > 0.0)
						v.set_metal_soft_weight((float)metal_soft_weight);
					std::cerr << "[LKina] Auto-detected " << ag4_metal_modes_to_str(detected_modes)
					          << " metal in ligand -> ligand-metal coordination maps enabled. "
					          << "Use --no_auto_metal to force pure AD4 scoring.\n";
					if (verbosity > 0)
						std::cout << "Auto-detected   : " << ag4_metal_modes_to_str(detected_modes)
						          << " metal in ligand → ligand-metal coordination maps enabled\n";
				}
			}
			if (metal_geometry_check_flag && !rigid_name.empty()) {
				ag4_metal_mode check_mode = ag4_metal_mode::none;
				if (!metal_mode_str.empty()) {
					std::istringstream tmp(metal_mode_str);
					std::string t;
					if (std::getline(tmp, t, ',')) {
						size_t s=t.find_first_not_of(' '), e=t.find_last_not_of(' ');
						if (s!=std::string::npos) check_mode=parse_metal_mode_token(t.substr(s,e-s+1));
					}
				} else {
					check_mode = detect_metal_mode_from_pdbqt(rigid_name);
				}
				metal_geometry_check(rigid_name, check_mode);
			}
			if (!metal_mode_str.empty()) {
				// Support comma-separated multi-metal: --metal_mode "fe3,zn"
				std::istringstream css(metal_mode_str);
				std::string tok;
				bool first_mode = true;
				bool parsed_contains_zn = false;
				while (std::getline(css, tok, ',')) {
					size_t s = tok.find_first_not_of(' ');
					size_t e = tok.find_last_not_of(' ');
					if (s == std::string::npos) continue;
					tok = tok.substr(s, e-s+1);
					ag4_metal_mode mm = parse_metal_mode_token(tok);
					if (mm == ag4_metal_mode::zn) {
						v.set_zn_mode(true);
						parsed_contains_zn = true;
					}
					if (first_mode) { v.set_metal_mode(mm); first_mode = false; }
					else v.add_metal_mode(mm);
				}
				if (zn_mode && !parsed_contains_zn)
					v.add_metal_mode(ag4_metal_mode::zn);
				ad4_metal_mode_configured = true;
				if (metal_soft_weight > 0.0) v.set_metal_soft_weight((float)metal_soft_weight);
			if (verbosity > 0) std::cout << "Metal mode      : " << metal_mode_str << " coordination potentials enabled\n";
			// O6: Print electrostatics weight hint for high-charge-density metals
			if (verbosity > 0) {
				// Check first mode token for elec hint
				std::istringstream tmp6(metal_mode_str);
				std::string t6;
				if (std::getline(tmp6, t6, ',')) {
					size_t s6=t6.find_first_not_of(' '), e6=t6.find_last_not_of(' ');
					if (s6!=std::string::npos) {
						ag4_metal_mode hm = parse_metal_mode_token(t6.substr(s6,e6-s6+1));
						const char* hint = ag4_elec_weight_hint(hm);
						if (hint)
							std::cout << "TIP (O6) High-charge metal: consider --weight_ad4_elec "
							          << hint << " (default 0.1406)\n";
					}
				}
			}
			}
			if (generate_maps && vm.count("batch") && metal_mode_str.empty() && !ad4_metal_mode_configured && !no_auto_metal) {
				for (const auto& batch_ligand_name : batch_ligand_names) {
					if (!detect_metal_modes_from_pdbqt(batch_ligand_name).empty()) {
						ad4_batch_ligand_metal_auto = true;
						if (verbosity > 0)
							std::cout << "Batch metal mode: per-ligand auto-detection enabled for ligand metals\n";
						break;
					}
				}
			}
			// O5: metal_bias — auto-add soft Gaussian attractor toward receptor metal centre (MBD-style)
			if (metal_bias && generate_maps && !rigid_name.empty()) {
				bool reactive_already = (!reactive_mode_str.empty() && reactive_mode_str != "off")
				                     || !reactive_preset_name.empty();
				if (reactive_already) {
					if (verbosity > 0)
						std::cout << "NOTE (O5): --metal_bias skipped because --reactive_mode is already active.\n";
				} else {
					// Parse receptor to find first metal atom position
					std::ifstream mbf(rigid_name);
					bool found_metal = false;
					double mbx=0,mby=0,mbz=0;
					static const std::set<std::string> met_types = {
						"Zn","Mg","Ca","Mn","Fe","Co","Ni","Cu","Pt","Pd","Ru","Ir","Au",
						"Cd","Hg","Na","K","V","Cr1","Cr","Ti","Sc","Mo","Rh","Re","Os"
					};
					if (mbf.is_open()) {
						std::string mbl;
						while (std::getline(mbf, mbl) && !found_metal) {
							if (mbl.size()<60) continue;
							std::string mrec=mbl.substr(0,6);
							if (mrec!="ATOM  " && mrec!="HETATM") continue;
							try { mbx=std::stod(mbl.substr(30,8)); mby=std::stod(mbl.substr(38,8)); mbz=std::stod(mbl.substr(46,8)); }
							catch(...) { continue; }
							std::istringstream miss(mbl);
							std::string mtok, mlast;
							while (miss >> mtok) mlast=mtok;
							if (met_types.count(mlast)) { found_metal=true; }
						}
					}
					if (found_metal) {
						// Inject as a distance-mode reactive attractor at the metal xyz
						reactive_mode_str = "distance";
						char mbuf[128];
						std::snprintf(mbuf,sizeof(mbuf),"%.4f,%.4f,%.4f",mbx,mby,mbz);
						reactive_rec_atom = mbuf;
						reactive_attractor_strength = metal_bias_strength;
						reactive_attractor_width    = metal_bias_width;
						reactive_programmatic_strength = true;
						reactive_programmatic_width = true;
						if (verbosity > 0)
							std::cout << "Metal bias (O5) : Gaussian attractor at (" << mbuf
							          << ")  depth=" << metal_bias_strength
							          << " kcal/mol  sigma=" << metal_bias_width << " A\n";
					} else if (verbosity > 0) {
						std::cout << "NOTE (O5): --metal_bias: no metal atom found in " << rigid_name << "\n";
					}
				}
			}
				// Set reactive options BEFORE map generation so the anchor is resolved
			bool reactive_requested = (!reactive_mode_str.empty() && reactive_mode_str != "off")
			                       || !reactive_preset_name.empty();
			if (reactive_requested) {
				if (!generate_maps && !vm.count("receptor")) {
					std::cerr << "ERROR: --reactive_mode/--reactive_preset with --maps requires --receptor for anchor resolution.\n";
					exit(EXIT_FAILURE);
				}
				if (generate_maps && !vm.count("receptor")) {
					std::cerr << "ERROR: --reactive_mode/--reactive_preset requires --receptor.\n";
					exit(EXIT_FAILURE);
				}
				ReactiveOptions ropts;
				// Step 1: apply preset defaults (fills mode, bond_length, attractor_*, angle_*)
				if (!reactive_preset_name.empty()) {
					ropts.preset_name = reactive_preset_name;
					if (!ropts.apply_preset()) {
						std::cerr << "ERROR: unknown --reactive_preset '" << reactive_preset_name
						          << "'. Valid: cys_michael, cys_sn2, ser_covalent, lys_targeting, boronic_acid, tyr_covalent\n";
						exit(EXIT_FAILURE);
					}
				}
				// Step 2: explicit --reactive_mode overrides preset mode
				if (!reactive_mode_str.empty()) {
					if (reactive_mode_str == "distance")    ropts.mode = reactive_mode::distance;
					else if (reactive_mode_str == "hybrid") ropts.mode = reactive_mode::hybrid;
					else {
						std::cerr << "ERROR: --reactive_mode must be 'distance' or 'hybrid'.\n";
						exit(EXIT_FAILURE);
					}
				}
				// Step 3: atom specs and debug flags (always passed through)
				ropts.receptor_atom_spec         = reactive_rec_atom;
				ropts.ligand_atom_spec           = reactive_lig_atom;
				ropts.receptor_frame_atom_spec   = reactive_frame_atom;
				ropts.ligand_frame_atom_spec     = reactive_lig_frame_atom;
				ropts.debug                      = reactive_debug;
				ropts.debug_energy               = reactive_debug_energy;
				ropts.gradient_check             = reactive_gradcheck;
				ropts.gradient_check_eps         = reactive_gradcheck_eps;
				if (vm.count("reactive_two_step") && vm["reactive_two_step"].as<bool>())
					ropts.two_step = true;
				if (option_explicit("reactive_presample_dist"))
					ropts.presample_max_dist = vm["reactive_presample_dist"].as<double>();
				if (vm.count("reactive_weak_attractor") && vm["reactive_weak_attractor"].as<bool>())
					ropts.weak_attractor = true;
				// Step 4: numeric CLI options override preset only when explicitly provided
				if (option_explicit("reactive_bond_length"))        ropts.bond_length        = reactive_bond_length;
				if (option_explicit("reactive_attractor_width") || reactive_programmatic_width)       ropts.attractor_width    = reactive_attractor_width;
				if (option_explicit("reactive_attractor_strength") || reactive_programmatic_strength) ropts.attractor_strength = reactive_attractor_strength;
				if (option_explicit("reactive_angle_strength"))     ropts.angle_strength     = reactive_angle_strength;
				if (option_explicit("reactive_hybrid_vdw_scale"))   ropts.hybrid_vdw_scale   = reactive_hybrid_vdw_scale;
				if (option_explicit("reactive_target_angle"))       ropts.target_angle_deg   = reactive_target_angle;
				if (option_explicit("reactive_angle_width"))        ropts.angle_width_deg    = reactive_angle_width;
				v.set_reactive_options(ropts);
				if (verbosity > 0) {
					const char* mstr = (ropts.mode == reactive_mode::hybrid) ? "hybrid" : "distance";
					std::cout << "Reactive mode   : " << mstr;
					if (!reactive_preset_name.empty()) std::cout << "  preset=" << reactive_preset_name;
					std::cout << "  rec=" << reactive_rec_atom
					          << "  lig=" << reactive_lig_atom << "\n";
					std::cout << "Reactive params : bond=" << ropts.bond_length
					          << " A  sigma=" << ropts.attractor_width
					          << " A  depth=" << ropts.attractor_strength
					          << " kcal/mol  angle=" << ropts.target_angle_deg
					          << " deg  width=" << ropts.angle_width_deg
					          << " deg  angle_k=" << ropts.angle_strength
					          << "  hybrid_vdw_scale=" << ropts.hybrid_vdw_scale << "\n";
				}
			}
				if (generate_maps && !ad4_batch_ligand_metal_auto) {
				// Inline computation — no external autogrid4 required
				v.compute_ad4_maps(rigid_name,
								   center_x, center_y, center_z,
								   size_x,   size_y,   size_z,
								   "", grid_spacing);
				if (vm.count("write_maps"))
					v.write_maps(out_maps);
			} else if (!generate_maps) {
				v.load_maps(maps);
				// It works, but why would you do this?!
				if (vm.count("write_maps"))
					v.write_maps(out_maps);
			}
		}

		if (vm.count("ligand")) {
			v.set_ligand_from_file(ligand_names);

			if (sf_name.compare("vina") == 0 || sf_name.compare("vinardo") == 0) {
				if (vm.count("maps")) {
					v.load_maps(maps);
				} else {
					// Will compute maps only for Vina atom types in the ligand(s)
					// In the case users ask for score and local only with the autobox arg, we compute the optimal box size for it/them.
					if ((score_only || local_only) && autobox) {
						std::vector<double> dim = v.grid_dimensions_from_ligand(buffer_size);
						v.compute_vina_maps(dim[0], dim[1], dim[2], dim[3], dim[4], dim[5], grid_spacing, force_even_voxels);
					} else {
						v.compute_vina_maps(center_x, center_y, center_z, size_x, size_y, size_z, grid_spacing, force_even_voxels);
					}

					if (vm.count("write_maps"))
						v.write_maps(out_maps);
				}
			}

			if (randomize_only) {
				v.randomize();
				v.write_pose(out_name);
			} else if (score_only) {
				std::vector<double> energies;
				if (std::isnan(unbound_energy)) {
					energies = v.score();
				} else {
					energies = v.score(unbound_energy);
				}
				v.show_score(energies);
			} else if (local_only) {
				std::vector<double> energies;
				energies = v.optimize();
				v.write_pose(out_name);
				v.show_score(energies);
			} else {
				v.global_search(exhaustiveness, num_modes, min_rmsd, max_evals);
				v.write_poses(out_name, num_modes, energy_range);
			}
		} else if (vm.count("batch")) {
			if (sf_name.compare("vina") == 0 || sf_name.compare("vinardo") == 0) {
				if (vm.count("maps")) {
					v.load_maps(maps);
				} else {
					// Will compute maps for all Vina atom types
					v.compute_vina_maps(center_x, center_y, center_z, size_x, size_y, size_z, grid_spacing);

					if (vm.count("write_maps"))
						v.write_maps(out_maps);
				}
			}

			std::set<std::string> repeated_names;
			std::set<std::string> raw_names;
			std::string name;
			VINA_RANGE(i, 0, batch_ligand_names.size()) {
				name = get_filename(batch_ligand_names[i]);
				if (raw_names.count(name)) {
					repeated_names.insert(name);
				}
				raw_names.insert(name);
			}
			std::unordered_map<std::string, int> instance_counter;

			std::size_t failed_ligand_parsing = 0;
			std::size_t failed_ligand_runtime = 0;
			VINA_RANGE(i, 0, batch_ligand_names.size()) {
				name = get_filename(batch_ligand_names[i]);
				if (repeated_names.count(name)) {
					if (instance_counter.count(name)) {
						instance_counter[name] += 1;
					} else {
						instance_counter[name] = 1;
					}
					out_name = default_output(name, out_dir, instance_counter[name]);
				} else {
					out_name = default_output(name, out_dir);
				}

				try {
					v.set_ligand_from_file(batch_ligand_names[i]);
					if (ad4_batch_ligand_metal_auto && !no_auto_metal) {
						std::vector<ag4_metal_mode> detected_modes = detect_metal_modes_from_pdbqt(batch_ligand_names[i]);
						apply_auto_metal_modes(v, detected_modes);
						if (metal_soft_weight > 0.0)
							v.set_metal_soft_weight((float)metal_soft_weight);
						if (verbosity > 0 && !detected_modes.empty())
							std::cout << "Auto-detected   : " << ag4_metal_modes_to_str(detected_modes)
							          << " metal in " << batch_ligand_names[i]
							          << " → ligand-metal coordination maps enabled\n";
						v.compute_ad4_maps(rigid_name,
										   center_x, center_y, center_z,
										   size_x,   size_y,   size_z,
										   "", grid_spacing);
					}

					if (randomize_only) {
						v.randomize();
						v.write_pose(out_name);
					} else if (score_only) {
						if (std::isnan(unbound_energy)) {
							v.score();
						} else {
							v.score(unbound_energy);
						}
					} else if (local_only) {
						v.optimize();
						v.write_pose(out_name);
					} else {
						v.global_search(exhaustiveness, num_modes, min_rmsd, max_evals);
						v.write_poses(out_name, num_modes, energy_range);
					}
				}
				catch(pdbqt_parse_error& e) {
					std::cerr << e.what();
					std::cout << "Failed parsing " << batch_ligand_names[i] << ". Skipping it.\n";
					failed_ligand_parsing++;
					continue;
				}
				catch(vina_runtime_error& e) {
					std::cerr << e.what();
					std::cout << "Failed processing " << batch_ligand_names[i] << ". Skipping it.\n";
					failed_ligand_runtime++;
					continue;
				}
			}
			if (repeated_names.size()) {
				std::cout << "Found " << repeated_names.size() << " repeated filenames in the input batch.\n";
				std::cout << "The corresponding output filenames are suffixed with _instance<n>_out.pdbqt\n";
			}
			if (failed_ligand_parsing) {
				std::cout << "Failed to parse " << failed_ligand_parsing << " ligands.\n";
			}
			if (failed_ligand_runtime) {
				std::cout << "Failed to process " << failed_ligand_runtime << " ligands.\n";
			}
		}
	}

	catch(pdbqt_parse_error& e) {
		std::cerr << e.what();
		return 1;
	}
	catch(file_error& e) {
		std::cerr << "\n\nError: could not open \"" << e.name.string() << "\" for " << (e.in ? "reading" : "writing") << ".\n";
		return 1;
	}
	catch(boost::filesystem::filesystem_error& e) {
		std::cerr << "\n\nFile system error: " << e.what() << '\n';
		return 1;
	}
	catch(usage_error& e) {
		std::cerr << "\n\nUsage error: " << e.what() << ".\n";
		return 1;
	}
	catch(std::bad_alloc&) {
		std::cerr << "\n\nError: insufficient memory!\n";
		return 1;
	}

	// Errors that shouldn't happen:
	catch(std::exception& e) {
		std::cerr << "\n\nAn error occurred: " << e.what() << ". " << error_message;
		return 1;
	}
	catch(internal_error& e) {
		std::cerr << "\n\nAn internal error occurred in " << e.file << "(" << e.line << "). " << error_message;
		return 1;
	}
	catch(...) {
		std::cerr << "\n\nAn unknown error occurred. " << error_message;
		return 1;
	}
}
