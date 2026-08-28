/*
   LKina — AutoDock4 force field parameter data
   Copyright (C) 2025 LK-Studio1128 and LKina contributors
   SPDX-License-Identifier: GPL-3.0-or-later

   This file embeds the AD4_parameters.dat content from AutoDock4
   (GPL-2.0 or later), extended with additional metal atom types.
   It is therefore distributed under the GNU General Public License v3.

   Parameter reference:
     Morris et al., J. Comput. Chem. 19:1639-1662 (1998)
     Morris et al., J. Comput. Chem. 30:2785-2791 (2009)
   AutoDock4 upstream: https://github.com/ccsb-scripps/AutoDock4

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

#include "ad4_parameter_data.h"

namespace ad4_param {

const double FE_coeff_vdW    = 0.1560;
const double FE_coeff_hbond  = 0.0974;
const double FE_coeff_estat  = 0.1465;
const double FE_coeff_desolv = 0.1159;
const double FE_coeff_tors   = 0.2744;

const char builtin_ad4_parameter_text[] =
"FE_coeff_vdW    0.1560\n"
"FE_coeff_hbond  0.0974\n"
"FE_coeff_estat  0.1465\n"
"FE_coeff_desolv 0.1159\n"
"FE_coeff_tors   0.2744\n"
"atom_par C   4.00 0.150 33.5103 -0.00143 0.0 0.0 0 -1 -1 0\n"
"atom_par A   4.00 0.150 33.5103 -0.00052 0.0 0.0 0 -1 -1 0\n"
"atom_par N   3.50 0.160 22.4493 -0.00162 0.0 0.0 0 -1 -1 1\n"
"atom_par NA  3.50 0.160 22.4493 -0.00162 1.9 5.0 4 -1 -1 1\n"
"atom_par NS  3.50 0.160 22.4493 -0.00162 1.9 5.0 3 -1 -1 1\n"
"atom_par OA  3.20 0.200 17.1573 -0.00251 1.9 5.0 3 -1 -1 2\n"
"atom_par OS  3.20 0.200 17.1573 -0.00251 1.9 5.0 3 -1 -1 2\n"
"atom_par SA  4.00 0.200 33.5103 -0.00214 2.5 1.0 5 -1 -1 6\n"
"atom_par S   4.00 0.200 33.5103 -0.00214 0.0 0.0 0 -1 -1 6\n"
"atom_par H   2.00 0.020 0.0000  0.00051  0.0 0.0 0 -1 -1 3\n"
"atom_par HD  2.00 0.020 0.0000  0.00051  0.0 0.0 2 -1 -1 3\n"
"atom_par HS  2.00 0.020 0.0000  0.00051  0.0 0.0 1 -1 -1 3\n"
"atom_par P   4.20 0.200 38.7924 -0.00110 0.0 0.0 0 -1 -1 5\n"
"atom_par Br  4.33 0.389 42.5661 -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Ca  1.98 0.550 2.7700  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Cl  4.09 0.276 35.8235 -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par F   3.09 0.080 15.4480 -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Fe  1.30 0.010 1.8400  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par I   4.72 0.550 55.0585 -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Mg  1.30 0.875 1.5600  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Mn  1.30 0.875 2.1400  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Zn  1.48 0.550 1.7000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par He  2.36 0.056 15.240  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Li  2.45 0.025 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Be  2.76 0.085 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par B   4.08 0.180 12.052  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Ne  3.24 0.042 15.440  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Na  3.98 0.030 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Al  4.49 0.505 11.278  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Si  4.30 0.402 12.175  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par K   3.81 0.035 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Sc  3.30 0.019 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Ti  3.18 0.017 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par V   3.14 0.016 12.000  -0.00110 0.0 0.0 0 -1 -1 1\n"
"atom_par Cr1 2.75 0.010 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Tg  3.07 0.067 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Co  2.87 0.014 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Ni  2.83 0.015 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Cu  3.50 0.005 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Ga  4.38 0.415 11.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Ge  4.28 0.379 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par As  4.23 0.309 13.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Se  4.21 0.291 14.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Kr  4.14 0.220 16.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Rb  4.11 0.040 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Sr  3.64 0.235 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Y   3.35 0.072 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Zr  3.12 0.069 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Nb  3.17 0.059 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Mo  3.05 0.056 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Tc  3.00 0.048 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Ru  2.96 0.056 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Rh  2.93 0.053 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Pd  1.34 0.048 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Ag  3.15 0.036 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Cd  2.85 0.228 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par In  4.46 0.599 11.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Sn  4.39 0.567 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Sb  4.42 0.449 13.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Te  4.47 0.398 14.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Xe  4.40 0.332 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Cs  4.52 0.045 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Ba  3.70 0.364 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par La  3.52 0.017 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Ce  3.56 0.013 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Pr  3.61 0.010 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Nd  3.58 0.010 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Pm  3.55 0.009 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Sm  3.52 0.008 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Eu  3.49 0.008 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Gd  3.37 0.009 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Tb  3.45 0.007 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Dy  3.43 0.007 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Ho  3.41 0.007 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Er  3.39 0.007 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Tm  3.37 0.006 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Yb  3.36 0.228 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Lu  3.64 0.041 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Hf  3.41 0.072 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Ta  3.71 0.081 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par W   3.07 0.067 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Re  2.95 0.066 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Os  3.12 0.120 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Ir  2.84 0.073 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Pt  2.75 0.080 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Au  3.29 0.039 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Hg  2.71 0.385 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Tl  4.35 0.680 11.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Pb  4.30 0.663 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Bi  4.37 0.518 13.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Po  4.71 0.325 14.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par At  4.75 0.284 15.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Rn  4.77 0.248 16.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Fr  4.90 0.050 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Ra  3.68 0.404 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Ac  3.48 0.033 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Th  3.40 0.026 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Pa  3.42 0.022 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par U   3.40 0.022 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Np  3.42 0.019 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Pu  3.42 0.016 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Am  3.38 0.014 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Cm  3.33 0.014 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Bk  3.34 0.013 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Cf  3.31 0.013 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par E   3.30 0.012 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
"atom_par Fm  3.29 0.012 12.000  -0.00110 0.0 0.0 0 -1 -1 4\n"
// Coordination pseudoatom types (TZ tetrahedral from official AD4Zn.dat;
// SQ/MH/JT analogous).  In the official parameter files these carry NO vdW
// interaction (Rii=1.0, epsii=0.0) — all ligand contacts run through the
// explicit nbp_r_eps overrides, e.g. "nbp_r_eps 0.25 23.2135 12 6 NA TZ".
// Filling epsii with the unweighted override value (the previous bug) made
// every non-override probe see a spurious LJ well against the pseudoatoms.
"atom_par TZ  1.00 0.000 0.0000  0.00000 0.0 0.0 0 -1 -1 0\n"
"atom_par SQ  1.00 0.000 0.0000  0.00000 0.0 0.0 0 -1 -1 0\n"
"atom_par MH  1.00 0.000 0.0000  0.00000 0.0 0.0 0 -1 -1 0\n"
"atom_par JT  1.00 0.000 0.0000  0.00000 0.0 0.0 0 -1 -1 0\n";

} // namespace ad4_param
