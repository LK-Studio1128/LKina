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

#ifndef VINA_ATOM_CONSTANTS_H
#define VINA_ATOM_CONSTANTS_H

#include "common.h"
#include "ad4_parameter_data.h"
#include <map>
#include <sstream>

// based on SY_TYPE_* but includes H
const sz EL_TYPE_H    =  0;
const sz EL_TYPE_C    =  1;
const sz EL_TYPE_N    =  2;
const sz EL_TYPE_O    =  3;
const sz EL_TYPE_S    =  4;
const sz EL_TYPE_P    =  5;
const sz EL_TYPE_F    =  6;
const sz EL_TYPE_Cl   =  7;
const sz EL_TYPE_Br   =  8;
const sz EL_TYPE_I    =  9;
const sz EL_TYPE_Si   = 10; // Silicon
const sz EL_TYPE_At   = 11; // Astatine
const sz EL_TYPE_Met  = 12;
const sz EL_TYPE_Dummy= 13;
const sz EL_TYPE_SIZE = 14;

// AutoDock4
const sz AD_TYPE_C    =  0;
const sz AD_TYPE_A    =  1;
const sz AD_TYPE_N    =  2;
const sz AD_TYPE_O    =  3;
const sz AD_TYPE_P    =  4;
const sz AD_TYPE_S    =  5;
const sz AD_TYPE_H    =  6; // non-polar hydrogen
const sz AD_TYPE_F    =  7;
const sz AD_TYPE_I    =  8;
const sz AD_TYPE_NA   =  9;
const sz AD_TYPE_OA   = 10;
const sz AD_TYPE_SA   = 11;
const sz AD_TYPE_HD   = 12;
const sz AD_TYPE_Mg   = 13;
const sz AD_TYPE_Mn   = 14;
const sz AD_TYPE_Zn   = 15;
const sz AD_TYPE_Ca   = 16;
const sz AD_TYPE_Fe   = 17;
const sz AD_TYPE_Cl   = 18;
const sz AD_TYPE_Br   = 19;
const sz AD_TYPE_Si   = 20; // Silicon
const sz AD_TYPE_At   = 21; // Astatine
const sz AD_TYPE_G0   = 22; // closure of cyclic molecules
const sz AD_TYPE_G1   = 23;
const sz AD_TYPE_G2   = 24;
const sz AD_TYPE_G3   = 25;
const sz AD_TYPE_CG0  = 26;
const sz AD_TYPE_CG1  = 27;
const sz AD_TYPE_CG2  = 28;
const sz AD_TYPE_CG3  = 29;
const sz AD_TYPE_W    = 30; // hydrated ligand
// --- Full-metal AD4 extension (auto-generated from AD4_parameters.dat; Phase 1 of FULL_METAL_AD4_PLAN.md) ---
const sz AD_TYPE_He   = 31;
const sz AD_TYPE_Li   = 32;
const sz AD_TYPE_Be   = 33;
const sz AD_TYPE_B    = 34;
const sz AD_TYPE_Ne   = 35;
const sz AD_TYPE_Na   = 36;
const sz AD_TYPE_Al   = 37;
const sz AD_TYPE_K    = 38;
const sz AD_TYPE_Sc   = 39;
const sz AD_TYPE_Ti   = 40;
const sz AD_TYPE_V    = 41;
const sz AD_TYPE_Co   = 42;
const sz AD_TYPE_Ni   = 43;
const sz AD_TYPE_Cu   = 44;
const sz AD_TYPE_Ga   = 45;
const sz AD_TYPE_Ge   = 46;
const sz AD_TYPE_As   = 47;
const sz AD_TYPE_Kr   = 48;
const sz AD_TYPE_Rb   = 49;
const sz AD_TYPE_Sr   = 50;
const sz AD_TYPE_Y    = 51;
const sz AD_TYPE_Zr   = 52;
const sz AD_TYPE_Nb   = 53;
const sz AD_TYPE_Mo   = 54;
const sz AD_TYPE_Tc   = 55;
const sz AD_TYPE_Ru   = 56;
const sz AD_TYPE_Rh   = 57;
const sz AD_TYPE_Pd   = 58;
const sz AD_TYPE_Ag   = 59;
const sz AD_TYPE_Cd   = 60;
const sz AD_TYPE_In   = 61;
const sz AD_TYPE_Sn   = 62;
const sz AD_TYPE_Sb   = 63;
const sz AD_TYPE_Te   = 64;
const sz AD_TYPE_Xe   = 65;
const sz AD_TYPE_Cs   = 66;
const sz AD_TYPE_Ba   = 67;
const sz AD_TYPE_La   = 68;
const sz AD_TYPE_Ce   = 69;
const sz AD_TYPE_Pr   = 70;
const sz AD_TYPE_Nd   = 71;
const sz AD_TYPE_Pm   = 72;
const sz AD_TYPE_Sm   = 73;
const sz AD_TYPE_Eu   = 74;
const sz AD_TYPE_Gd   = 75;
const sz AD_TYPE_Tb   = 76;
const sz AD_TYPE_Dy   = 77;
const sz AD_TYPE_Ho   = 78;
const sz AD_TYPE_Er   = 79;
const sz AD_TYPE_Tm   = 80;
const sz AD_TYPE_Yb   = 81;
const sz AD_TYPE_Lu   = 82;
const sz AD_TYPE_Hf   = 83;
const sz AD_TYPE_Ta   = 84;
const sz AD_TYPE_Re   = 85;
const sz AD_TYPE_Os   = 86;
const sz AD_TYPE_Ir   = 87;
const sz AD_TYPE_Pt   = 88;
const sz AD_TYPE_Au   = 89;
const sz AD_TYPE_Hg   = 90;
const sz AD_TYPE_Tl   = 91;
const sz AD_TYPE_Pb   = 92;
const sz AD_TYPE_Bi   = 93;
const sz AD_TYPE_Po   = 94;
const sz AD_TYPE_Rn   = 95;
const sz AD_TYPE_Fr   = 96;
const sz AD_TYPE_Ra   = 97;
const sz AD_TYPE_Ac   = 98;
const sz AD_TYPE_Th   = 99;
const sz AD_TYPE_Pa   = 100;
const sz AD_TYPE_U    = 101;
const sz AD_TYPE_Np   = 102;
const sz AD_TYPE_Pu   = 103;
const sz AD_TYPE_Am   = 104;
const sz AD_TYPE_Cm   = 105;
const sz AD_TYPE_Bk   = 106;
const sz AD_TYPE_Cf   = 107;
const sz AD_TYPE_E    = 108;
const sz AD_TYPE_Fm   = 109;
const sz AD_TYPE_Cr1  = 110;
const sz AD_TYPE_Tg   = 111; // Tungsten (Wolfram), PDBQT type "W" — distinct from AD_TYPE_W=30 (water probe)
const sz AD_TYPE_Se   = 112; // Selenium — promoted from atom_equivalence (Se→S) to own AD_TYPE with Se-specific parameters
const sz AD_TYPE_TZ   = 113; // Zinc tetrahedral pseudoatom (AutoDock4Zn) — placed at empty coordination sites around Zn
const sz AD_TYPE_SQ   = 114; // Square-planar coordination pseudoatom (M4) — placed at vacant sites around Pt/Pd/Ni/Cu (sq-planar) or Au (linear)
const sz AD_TYPE_MH   = 115; // Octahedral (Metal Hexacoordinate) pseudoatom — placed at vacant octahedral sites around Fe/Mg/Mn/Co/Ru/Ir etc.
const sz AD_TYPE_JT   = 116;
const sz AD_TYPE_SIZE = 117;

// X-Score
const sz XS_TYPE_C_H   =  0;
const sz XS_TYPE_C_P   =  1;
const sz XS_TYPE_N_P   =  2;
const sz XS_TYPE_N_D   =  3;
const sz XS_TYPE_N_A   =  4;
const sz XS_TYPE_N_DA  =  5;
const sz XS_TYPE_O_P   =  6;
const sz XS_TYPE_O_D   =  7;
const sz XS_TYPE_O_A   =  8;
const sz XS_TYPE_O_DA  =  9;
const sz XS_TYPE_S_P   = 10;
const sz XS_TYPE_P_P   = 11;
const sz XS_TYPE_F_H   = 12;
const sz XS_TYPE_Cl_H  = 13;
const sz XS_TYPE_Br_H  = 14;
const sz XS_TYPE_I_H   = 15;
const sz XS_TYPE_Si    = 16;
const sz XS_TYPE_At    = 17;
const sz XS_TYPE_Met_D = 18;
const sz XS_TYPE_C_H_CG0 = 19;
const sz XS_TYPE_C_P_CG0 = 20;
const sz XS_TYPE_G0      = 21;
const sz XS_TYPE_C_H_CG1 = 22;
const sz XS_TYPE_C_P_CG1 = 23;
const sz XS_TYPE_G1      = 24;
const sz XS_TYPE_C_H_CG2 = 25;
const sz XS_TYPE_C_P_CG2 = 26;
const sz XS_TYPE_G2      = 27;
const sz XS_TYPE_C_H_CG3 = 28;
const sz XS_TYPE_C_P_CG3 = 29;
const sz XS_TYPE_G3      = 30;
const sz XS_TYPE_W       = 31;
const sz XS_TYPE_SIZE    = 32;

// DrugScore-CSD
const sz SY_TYPE_C_3   = 0;
const sz SY_TYPE_C_2   = 1;
const sz SY_TYPE_C_ar  = 2;
const sz SY_TYPE_C_cat = 3;
const sz SY_TYPE_N_3   = 4;
const sz SY_TYPE_N_ar  = 5;
const sz SY_TYPE_N_am  = 6;
const sz SY_TYPE_N_pl3 = 7;
const sz SY_TYPE_O_3   = 8;
const sz SY_TYPE_O_2   = 9;
const sz SY_TYPE_O_co2 = 10;
const sz SY_TYPE_S     = 11;
const sz SY_TYPE_P     = 12;
const sz SY_TYPE_F     = 13;
const sz SY_TYPE_Cl    = 14;
const sz SY_TYPE_Br    = 15;
const sz SY_TYPE_I     = 16;
const sz SY_TYPE_Met   = 17;
const sz SY_TYPE_SIZE  = 18;

struct atom_kind {
	std::string name;
	fl radius;
	fl depth;
	fl hb_depth;
	fl hb_radius;
	fl solvation;
	fl volume;
	fl covalent_radius;
};

inline fl ad4_covalent_radius_from_name(const std::string& name) {
	if (name == "C" || name == "A" || name == "G0" || name == "G1" || name == "G2" || name == "G3" ||
	    name == "CG0" || name == "CG1" || name == "CG2" || name == "CG3") return 0.77;
	if (name == "N" || name == "NA") return 0.75;
	if (name == "O" || name == "OA") return 0.73;
	if (name == "P") return 1.06;
	if (name == "S" || name == "SA") return 1.02;
	if (name == "H" || name == "HD") return 0.37;
	if (name == "F") return 0.71;
	if (name == "Cl") return 0.99;
	if (name == "Br") return 1.14;
	if (name == "I") return 1.33;
	if (name == "Si") return 1.11;
	if (name == "At") return 1.44;
	if (name == "Mg") return 1.30;
	if (name == "Mn") return 1.39;
	if (name == "Zn") return 1.31;
	if (name == "Ca") return 1.74;
	if (name == "Fe") return 1.25;
	if (name == "Cr1") return 1.28;
	if (name == "Tg") return 1.41;
	if (name == "Se") return 1.20;
	if (name == "TZ" || name == "SQ" || name == "MH" || name == "JT") return 0.25;
	if (name == "W") return 0.00;
	return 1.75;
}

inline std::vector<atom_kind> make_atom_kind_data() {
	std::map<std::string, atom_kind> parsed;
	std::istringstream iss(ad4_param::builtin_ad4_parameter_text);
	std::string line;
	while (std::getline(iss, line)) {
		if (line.size() < 8 || line.compare(0, 8, "atom_par") != 0) continue;
		std::istringstream ls(line);
		std::string tag, name;
		double rii = 0.0, eps = 0.0, vol = 0.0, sol = 0.0;
		double hb_depth_raw = 0.0, hb_radius_raw = 0.0;
		int i1 = 0, i2 = 0, i3 = 0, i4 = 0;
		ls >> tag >> name >> rii >> eps >> vol >> sol >> hb_depth_raw >> hb_radius_raw >> i1 >> i2 >> i3 >> i4;
		atom_kind ak;
		ak.name = name;
		ak.radius = (fl)(rii * 0.5);
		ak.depth = (fl)eps;
		ak.hb_depth = 0.0;
		ak.hb_radius = 0.0;
		ak.solvation = (fl)sol;
		ak.volume = (fl)vol;
		ak.covalent_radius = ad4_covalent_radius_from_name(name);
		if (name == "NA") { ak.hb_depth = -5.0; ak.hb_radius = 1.9; }
		else if (name == "OA") { ak.hb_depth = -5.0; ak.hb_radius = 1.9; }
		else if (name == "SA") { ak.hb_depth = -1.0; ak.hb_radius = 2.5; }
		else if (name == "HD") { ak.hb_depth =  1.0; ak.hb_radius = 0.0; }
		parsed[name] = ak;
	}
	std::vector<atom_kind> out(AD_TYPE_SIZE);
	const auto set = [&](sz i, const std::string& name, const std::string& from) {
		std::map<std::string, atom_kind>::const_iterator found = parsed.find(from);
		VINA_CHECK(found != parsed.end());
		out[i] = found->second;
		out[i].name = name;
	};
	const auto set_same = [&](sz i, const std::string& name) {
		set(i, name, name);
	};
	const auto set_manual = [&](sz i, const std::string& name, fl radius, fl depth, fl hb_depth, fl hb_radius, fl solvation, fl volume, fl covalent_radius) {
		atom_kind ak;
		ak.name = name;
		ak.radius = radius;
		ak.depth = depth;
		ak.hb_depth = hb_depth;
		ak.hb_radius = hb_radius;
		ak.solvation = solvation;
		ak.volume = volume;
		ak.covalent_radius = covalent_radius;
		out[i] = ak;
	};
	set_same(AD_TYPE_C, "C");
	set_same(AD_TYPE_A, "A");
	set_same(AD_TYPE_N, "N");
	set(AD_TYPE_O, "O", "OA");
	out[AD_TYPE_O].hb_depth = 0.0;
	out[AD_TYPE_O].hb_radius = 0.0;
	set_same(AD_TYPE_P, "P");
	set_same(AD_TYPE_S, "S");
	set_same(AD_TYPE_H, "H");
	set_same(AD_TYPE_F, "F");
	set_same(AD_TYPE_I, "I");
	set_same(AD_TYPE_NA, "NA");
	set_same(AD_TYPE_OA, "OA");
	set_same(AD_TYPE_SA, "SA");
	set_same(AD_TYPE_HD, "HD");
	set_same(AD_TYPE_Mg, "Mg");
	set_same(AD_TYPE_Mn, "Mn");
	set_same(AD_TYPE_Zn, "Zn");
	set_same(AD_TYPE_Ca, "Ca");
	set_same(AD_TYPE_Fe, "Fe");
	set_same(AD_TYPE_Cl, "Cl");
	set_same(AD_TYPE_Br, "Br");
	set_same(AD_TYPE_Si, "Si");
	set_same(AD_TYPE_At, "At");
	set_manual(AD_TYPE_G0, "G0", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.77);
	set_manual(AD_TYPE_G1, "G1", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.77);
	set_manual(AD_TYPE_G2, "G2", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.77);
	set_manual(AD_TYPE_G3, "G3", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.77);
	set_manual(AD_TYPE_CG0, "CG0", 2.0, 0.15, 0.0, 0.0, -0.00143, 33.51030, 0.77);
	set_manual(AD_TYPE_CG1, "CG1", 2.0, 0.15, 0.0, 0.0, -0.00143, 33.51030, 0.77);
	set_manual(AD_TYPE_CG2, "CG2", 2.0, 0.15, 0.0, 0.0, -0.00143, 33.51030, 0.77);
	set_manual(AD_TYPE_CG3, "CG3", 2.0, 0.15, 0.0, 0.0, -0.00143, 33.51030, 0.77);
	set_manual(AD_TYPE_W, "W", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	set_same(AD_TYPE_He, "He");
	set_same(AD_TYPE_Li, "Li");
	set_same(AD_TYPE_Be, "Be");
	set_same(AD_TYPE_B, "B");
	set_same(AD_TYPE_Ne, "Ne");
	set_same(AD_TYPE_Na, "Na");
	set_same(AD_TYPE_Al, "Al");
	set_same(AD_TYPE_K, "K");
	set_same(AD_TYPE_Sc, "Sc");
	set_same(AD_TYPE_Ti, "Ti");
	set_same(AD_TYPE_V, "V");
	set_same(AD_TYPE_Co, "Co");
	set_same(AD_TYPE_Ni, "Ni");
	set_same(AD_TYPE_Cu, "Cu");
	set_same(AD_TYPE_Ga, "Ga");
	set_same(AD_TYPE_Ge, "Ge");
	set_same(AD_TYPE_As, "As");
	set_same(AD_TYPE_Kr, "Kr");
	set_same(AD_TYPE_Rb, "Rb");
	set_same(AD_TYPE_Sr, "Sr");
	set_same(AD_TYPE_Y, "Y");
	set_same(AD_TYPE_Zr, "Zr");
	set_same(AD_TYPE_Nb, "Nb");
	set_same(AD_TYPE_Mo, "Mo");
	set_same(AD_TYPE_Tc, "Tc");
	set_same(AD_TYPE_Ru, "Ru");
	set_same(AD_TYPE_Rh, "Rh");
	set_same(AD_TYPE_Pd, "Pd");
	set_same(AD_TYPE_Ag, "Ag");
	set_same(AD_TYPE_Cd, "Cd");
	set_same(AD_TYPE_In, "In");
	set_same(AD_TYPE_Sn, "Sn");
	set_same(AD_TYPE_Sb, "Sb");
	set_same(AD_TYPE_Te, "Te");
	set_same(AD_TYPE_Xe, "Xe");
	set_same(AD_TYPE_Cs, "Cs");
	set_same(AD_TYPE_Ba, "Ba");
	set_same(AD_TYPE_La, "La");
	set_same(AD_TYPE_Ce, "Ce");
	set_same(AD_TYPE_Pr, "Pr");
	set_same(AD_TYPE_Nd, "Nd");
	set_same(AD_TYPE_Pm, "Pm");
	set_same(AD_TYPE_Sm, "Sm");
	set_same(AD_TYPE_Eu, "Eu");
	set_same(AD_TYPE_Gd, "Gd");
	set_same(AD_TYPE_Tb, "Tb");
	set_same(AD_TYPE_Dy, "Dy");
	set_same(AD_TYPE_Ho, "Ho");
	set_same(AD_TYPE_Er, "Er");
	set_same(AD_TYPE_Tm, "Tm");
	set_same(AD_TYPE_Yb, "Yb");
	set_same(AD_TYPE_Lu, "Lu");
	set_same(AD_TYPE_Hf, "Hf");
	set_same(AD_TYPE_Ta, "Ta");
	set_same(AD_TYPE_Re, "Re");
	set_same(AD_TYPE_Os, "Os");
	set_same(AD_TYPE_Ir, "Ir");
	set_same(AD_TYPE_Pt, "Pt");
	set_same(AD_TYPE_Au, "Au");
	set_same(AD_TYPE_Hg, "Hg");
	set_same(AD_TYPE_Tl, "Tl");
	set_same(AD_TYPE_Pb, "Pb");
	set_same(AD_TYPE_Bi, "Bi");
	set_same(AD_TYPE_Po, "Po");
	set_same(AD_TYPE_Rn, "Rn");
	set_same(AD_TYPE_Fr, "Fr");
	set_same(AD_TYPE_Ra, "Ra");
	set_same(AD_TYPE_Ac, "Ac");
	set_same(AD_TYPE_Th, "Th");
	set_same(AD_TYPE_Pa, "Pa");
	set_same(AD_TYPE_U, "U");
	set_same(AD_TYPE_Np, "Np");
	set_same(AD_TYPE_Pu, "Pu");
	set_same(AD_TYPE_Am, "Am");
	set_same(AD_TYPE_Cm, "Cm");
	set_same(AD_TYPE_Bk, "Bk");
	set_same(AD_TYPE_Cf, "Cf");
	set_same(AD_TYPE_E, "E");
	set_same(AD_TYPE_Fm, "Fm");
	set_same(AD_TYPE_Cr1, "Cr1");
	set_same(AD_TYPE_Tg, "Tg");
	set_same(AD_TYPE_Se, "Se");
	set_same(AD_TYPE_TZ, "TZ");
	set_same(AD_TYPE_SQ, "SQ");
	set_same(AD_TYPE_MH, "MH");
	set_same(AD_TYPE_JT, "JT");
	return out;
}

inline const std::vector<atom_kind>& atom_kind_data_ref() {
	static const std::vector<atom_kind> data = make_atom_kind_data();
	return data;
}

const fl metal_solvation_parameter = -0.00110;

const fl metal_covalent_radius = 1.75; // for metals not on the list // FIXME this info should be moved to non_ad_metals

const sz atom_kinds_size = AD_TYPE_SIZE;

struct atom_equivalence {
	std::string name;
	std::string to;
};

const atom_equivalence atom_equivalence_data[] = {
	{"Cr",  "Cr1"} // Se promoted to own AD_TYPE; Cr still maps to Cr1 for PDBQT compatibility
};

const sz atom_equivalences_size = sizeof(atom_equivalence_data) / sizeof(const atom_equivalence);

struct acceptor_kind {
	sz ad_type;
	fl radius;
	fl depth;
};

const acceptor_kind acceptor_kind_data[] = { // ad_type, optimal length, depth
	{AD_TYPE_NA, 1.9, 5.0},
	{AD_TYPE_OA, 1.9, 5.0},
	{AD_TYPE_SA, 2.5, 1.0}
};

const sz acceptor_kinds_size = sizeof(acceptor_kind_data) / sizeof(acceptor_kind);

inline bool ad_is_hydrogen(sz ad) {
	return ad == AD_TYPE_H || ad == AD_TYPE_HD;
}

inline bool ad_is_heteroatom(sz ad) { // returns false for ad >= AD_TYPE_SIZE
	return ad != AD_TYPE_A && ad != AD_TYPE_C  && 
		   ad != AD_TYPE_H && ad != AD_TYPE_HD && 
		   ad < AD_TYPE_SIZE;
}

inline sz ad_type_to_el_type(sz t) {
	switch(t) {
		case AD_TYPE_C    : return EL_TYPE_C;
		case AD_TYPE_A    : return EL_TYPE_C;
		case AD_TYPE_N    : return EL_TYPE_N;
		case AD_TYPE_O    : return EL_TYPE_O;
		case AD_TYPE_P    : return EL_TYPE_P;
		case AD_TYPE_S    : return EL_TYPE_S;
		case AD_TYPE_H    : return EL_TYPE_H;
		case AD_TYPE_F    : return EL_TYPE_F;
		case AD_TYPE_I    : return EL_TYPE_I;
		case AD_TYPE_NA   : return EL_TYPE_N;
		case AD_TYPE_OA   : return EL_TYPE_O;
		case AD_TYPE_SA   : return EL_TYPE_S;
		case AD_TYPE_HD   : return EL_TYPE_H;
		case AD_TYPE_Mg   : return EL_TYPE_Met;
		case AD_TYPE_Mn   : return EL_TYPE_Met;
		case AD_TYPE_Zn   : return EL_TYPE_Met;
		case AD_TYPE_Ca   : return EL_TYPE_Met;
		case AD_TYPE_Fe   : return EL_TYPE_Met;
		case AD_TYPE_Cl   : return EL_TYPE_Cl;
		case AD_TYPE_Br   : return EL_TYPE_Br;
		case AD_TYPE_Si   : return EL_TYPE_Si;
		case AD_TYPE_At   : return EL_TYPE_At;
		case AD_TYPE_CG0  : return EL_TYPE_C;
		case AD_TYPE_CG1  : return EL_TYPE_C;
		case AD_TYPE_CG2  : return EL_TYPE_C;
		case AD_TYPE_CG3  : return EL_TYPE_C;
		case AD_TYPE_G0   : return EL_TYPE_Dummy;
		case AD_TYPE_G1   : return EL_TYPE_Dummy;
		case AD_TYPE_G2   : return EL_TYPE_Dummy;
		case AD_TYPE_G3   : return EL_TYPE_Dummy;
		case AD_TYPE_W    : return EL_TYPE_Dummy;
		// --- Full-metal AD4 extension: all extra elements map to EL_TYPE_Met ---
		case AD_TYPE_He   : return EL_TYPE_Met;
		case AD_TYPE_Li   : return EL_TYPE_Met;
		case AD_TYPE_Be   : return EL_TYPE_Met;
		case AD_TYPE_B    : return EL_TYPE_Met;
		case AD_TYPE_Ne   : return EL_TYPE_Met;
		case AD_TYPE_Na   : return EL_TYPE_Met;
		case AD_TYPE_Al   : return EL_TYPE_Met;
		case AD_TYPE_K    : return EL_TYPE_Met;
		case AD_TYPE_Sc   : return EL_TYPE_Met;
		case AD_TYPE_Ti   : return EL_TYPE_Met;
		case AD_TYPE_V    : return EL_TYPE_Met;
		case AD_TYPE_Co   : return EL_TYPE_Met;
		case AD_TYPE_Ni   : return EL_TYPE_Met;
		case AD_TYPE_Cu   : return EL_TYPE_Met;
		case AD_TYPE_Ga   : return EL_TYPE_Met;
		case AD_TYPE_Ge   : return EL_TYPE_Met;
		case AD_TYPE_As   : return EL_TYPE_Met;
		case AD_TYPE_Kr   : return EL_TYPE_Met;
		case AD_TYPE_Rb   : return EL_TYPE_Met;
		case AD_TYPE_Sr   : return EL_TYPE_Met;
		case AD_TYPE_Y    : return EL_TYPE_Met;
		case AD_TYPE_Zr   : return EL_TYPE_Met;
		case AD_TYPE_Nb   : return EL_TYPE_Met;
		case AD_TYPE_Mo   : return EL_TYPE_Met;
		case AD_TYPE_Tc   : return EL_TYPE_Met;
		case AD_TYPE_Ru   : return EL_TYPE_Met;
		case AD_TYPE_Rh   : return EL_TYPE_Met;
		case AD_TYPE_Pd   : return EL_TYPE_Met;
		case AD_TYPE_Ag   : return EL_TYPE_Met;
		case AD_TYPE_Cd   : return EL_TYPE_Met;
		case AD_TYPE_In   : return EL_TYPE_Met;
		case AD_TYPE_Sn   : return EL_TYPE_Met;
		case AD_TYPE_Sb   : return EL_TYPE_Met;
		case AD_TYPE_Te   : return EL_TYPE_Met;
		case AD_TYPE_Xe   : return EL_TYPE_Met;
		case AD_TYPE_Cs   : return EL_TYPE_Met;
		case AD_TYPE_Ba   : return EL_TYPE_Met;
		case AD_TYPE_La   : return EL_TYPE_Met;
		case AD_TYPE_Ce   : return EL_TYPE_Met;
		case AD_TYPE_Pr   : return EL_TYPE_Met;
		case AD_TYPE_Nd   : return EL_TYPE_Met;
		case AD_TYPE_Pm   : return EL_TYPE_Met;
		case AD_TYPE_Sm   : return EL_TYPE_Met;
		case AD_TYPE_Eu   : return EL_TYPE_Met;
		case AD_TYPE_Gd   : return EL_TYPE_Met;
		case AD_TYPE_Tb   : return EL_TYPE_Met;
		case AD_TYPE_Dy   : return EL_TYPE_Met;
		case AD_TYPE_Ho   : return EL_TYPE_Met;
		case AD_TYPE_Er   : return EL_TYPE_Met;
		case AD_TYPE_Tm   : return EL_TYPE_Met;
		case AD_TYPE_Yb   : return EL_TYPE_Met;
		case AD_TYPE_Lu   : return EL_TYPE_Met;
		case AD_TYPE_Hf   : return EL_TYPE_Met;
		case AD_TYPE_Ta   : return EL_TYPE_Met;
		case AD_TYPE_Re   : return EL_TYPE_Met;
		case AD_TYPE_Os   : return EL_TYPE_Met;
		case AD_TYPE_Ir   : return EL_TYPE_Met;
		case AD_TYPE_Pt   : return EL_TYPE_Met;
		case AD_TYPE_Au   : return EL_TYPE_Met;
		case AD_TYPE_Hg   : return EL_TYPE_Met;
		case AD_TYPE_Tl   : return EL_TYPE_Met;
		case AD_TYPE_Pb   : return EL_TYPE_Met;
		case AD_TYPE_Bi   : return EL_TYPE_Met;
		case AD_TYPE_Po   : return EL_TYPE_Met;
		case AD_TYPE_Rn   : return EL_TYPE_Met;
		case AD_TYPE_Fr   : return EL_TYPE_Met;
		case AD_TYPE_Ra   : return EL_TYPE_Met;
		case AD_TYPE_Ac   : return EL_TYPE_Met;
		case AD_TYPE_Th   : return EL_TYPE_Met;
		case AD_TYPE_Pa   : return EL_TYPE_Met;
		case AD_TYPE_U    : return EL_TYPE_Met;
		case AD_TYPE_Np   : return EL_TYPE_Met;
		case AD_TYPE_Pu   : return EL_TYPE_Met;
		case AD_TYPE_Am   : return EL_TYPE_Met;
		case AD_TYPE_Cm   : return EL_TYPE_Met;
		case AD_TYPE_Bk   : return EL_TYPE_Met;
		case AD_TYPE_Cf   : return EL_TYPE_Met;
		case AD_TYPE_E    : return EL_TYPE_Met;
		case AD_TYPE_Fm   : return EL_TYPE_Met;
		case AD_TYPE_Cr1  : return EL_TYPE_Met;
		case AD_TYPE_Tg   : return EL_TYPE_Met;
		case AD_TYPE_Se   : return EL_TYPE_Met;
		case AD_TYPE_TZ   : return EL_TYPE_Dummy; // TZ is a pseudoatom — no real element
		case AD_TYPE_SQ   : return EL_TYPE_Dummy; // SQ is a pseudoatom — no real element
		case AD_TYPE_MH   : return EL_TYPE_Dummy; // MH is a pseudoatom — no real element
		case AD_TYPE_JT   : return EL_TYPE_Dummy; // JT is a pseudoatom — no real element
		case AD_TYPE_SIZE : return EL_TYPE_SIZE;
		default: VINA_CHECK(false);
	}
	return EL_TYPE_SIZE; // to placate the compiler in case of warnings - it should never get here though
}

const fl xs_vdw_radii[] = {
	1.9, // C_H
	1.9, // C_P
	1.8, // N_P
	1.8, // N_D
	1.8, // N_A
	1.8, // N_DA
	1.7, // O_P
	1.7, // O_D
	1.7, // O_A
	1.7, // O_DA
	2.0, // S_P
	2.1, // P_P
	1.5, // F_H
	1.8, // Cl_H
	2.0, // Br_H
	2.2, // I_H
    2.2, // Si
    2.3, // At
	1.2, // Met_D
    1.9, // C_H_CG0
    1.9, // C_P_CG0
    1.9, // C_H_CG1
    1.9, // C_P_CG1
    1.9, // C_H_CG2
    1.9, // C_P_CG2
    1.9, // C_H_CG3
    1.9, // C_P_CG3
    0.0, // G0
    0.0, // G1
    0.0, // G2
    0.0, // G3
    0.0  // W
};

const fl xs_vinardo_vdw_radii[] = {
	2.0, // C_H
	2.0, // C_P
	1.7, // N_P
	1.7, // N_D
	1.7, // N_A
	1.7, // N_DA
	1.6, // O_P
	1.6, // O_D
	1.6, // O_A
	1.6, // O_DA
	2.0, // S_P
	2.1, // P_P
	1.5, // F_H
	1.8, // Cl_H
	2.0, // Br_H
	2.2, // I_H
	2.2, // Si
	2.3, // At
	1.2, // Met_D
	2.0, // C_H_CG0
	2.0, // C_P_CG0
	2.0, // C_H_CG1
	2.0, // C_P_CG1
	2.0, // C_H_CG2
	2.0, // C_P_CG2
	2.0, // C_H_CG3
	2.0, // C_P_CG3
	0.0, // G0
	0.0, // G1
	0.0, // G2
	0.0, // G3
	0.0	 // W
};

inline fl xs_radius(sz t) {
	const sz n = sizeof(xs_vdw_radii) / sizeof(const fl);
	assert(n == XS_TYPE_SIZE);
	assert(t < n);
	return xs_vdw_radii[t];
}

inline fl xs_vinardo_radius(sz t) {
	const sz n = sizeof(xs_vdw_radii) / sizeof(const fl);
	assert(n == XS_TYPE_SIZE);
	assert(t < n);
	return xs_vinardo_vdw_radii[t];
}

const std::string non_ad_metal_names[] = { // expand as necessary
	"Cu", "Fe", "Na", "K", "Hg", "Co", "U", "Cd", "Ni"
};

inline bool is_non_ad_metal_name(const std::string& name) {
	const sz s = sizeof(non_ad_metal_names) / sizeof(const std::string);
	VINA_FOR(i, s)
		if(non_ad_metal_names[i] == name)
			return true;
	return false;
}

inline bool xs_is_hydrophobic(sz xs) {
	return xs == XS_TYPE_C_H || 
		   xs == XS_TYPE_F_H ||
		   xs == XS_TYPE_Cl_H ||
		   xs == XS_TYPE_Br_H || 
		   xs == XS_TYPE_I_H;
}

inline bool xs_is_acceptor(sz xs) {
	return xs == XS_TYPE_N_A ||
		   xs == XS_TYPE_N_DA ||
		   xs == XS_TYPE_O_A ||
		   xs == XS_TYPE_O_DA;
}

inline bool xs_is_donor(sz xs) {
	return xs == XS_TYPE_N_D ||
		   xs == XS_TYPE_N_DA ||
		   xs == XS_TYPE_O_D ||
		   xs == XS_TYPE_O_DA ||
		   xs == XS_TYPE_Met_D;
}

inline bool xs_donor_acceptor(sz t1, sz t2) {
	return xs_is_donor(t1) && xs_is_acceptor(t2);
}

inline bool xs_h_bond_possible(sz t1, sz t2) {
	return xs_donor_acceptor(t1, t2) || xs_donor_acceptor(t2, t1);
}

inline const atom_kind& ad_type_property(sz i) {
	const std::vector<atom_kind>& data = atom_kind_data_ref();
	VINA_CHECK(data.size() == atom_kinds_size);
	VINA_CHECK(i < data.size());
	return data[i];
}

inline sz string_to_ad_type(const std::string& name) { // returns AD_TYPE_SIZE if not found (no exceptions thrown, because metals unknown to AD4 are not exceptional)
	const std::vector<atom_kind>& data = atom_kind_data_ref();
	VINA_CHECK(data.size() == atom_kinds_size);
	VINA_FOR(i, data.size())
		if(data[i].name == name)
			return i;
	VINA_FOR(i, atom_equivalences_size)
		if(atom_equivalence_data[i].name == name)
			return string_to_ad_type(atom_equivalence_data[i].to);
	return AD_TYPE_SIZE;
}

inline fl max_covalent_radius() {
	fl tmp = 0;
	const std::vector<atom_kind>& data = atom_kind_data_ref();
	VINA_CHECK(data.size() == atom_kinds_size);
	VINA_FOR(i, data.size())
		if(data[i].covalent_radius > tmp)
			tmp = data[i].covalent_radius;
	return tmp;
}

#endif
