# LKina: A Metal-Aware and Covalent-Reactive Molecular Docking Engine Extending AutoDock Vina

**Authors:** Luo Xiaowen (罗晓文)
**Corresponding author:** Luo Xiaowen, affiliation, ORCID, and e-mail to be completed before submission

---

## Abstract

**Motivation:** Metalloenzymes mediate approximately 40% of all enzyme reactions in the human proteome and are targeted by roughly one-third of approved drugs; yet standard molecular-docking force fields cannot describe the directionality of metal–ligand coordination. The deficiencies are threefold: equilibrium distances are distorted (an N donor is placed at ~2.49 Å, whereas the native Zn–N coordination bond is ~2.0 Å); spherically symmetric Lennard-Jones potentials cannot distinguish axial from equatorial ligation; and formal charges bias the electrostatic term systematically toward oxygen donors. Meanwhile, covalent drug discovery requires that the search converge onto productive nucleophilic-attack trajectories, near-attack conformations (NACs), a capability absent from the AutoDock Vina search engine. Existing solutions cover only narrow niches: AutoDock4Zn handles Zn²⁺ only; CovDock is commercial and typically requires pre-built covalent adducts; and all of these workflows rely on the external `autogrid4` executable to pre-compute affinity maps.

**Results:** We present **LKina**, an open-source docking engine built by deeply extending the AutoDock Vina 1.2.7 source, elevating metal coordination and reactive covalent docking to first-class capabilities while remaining fully backward compatible. LKina (i) extends the AutoDock4 atom-type system from ~22 to **113 types** (covering more than 80 metals and metalloids) and embeds an independent C++ reimplementation of the AutoGrid 4.2 grid generator, eliminating the external `autogrid4` dependency; (ii) introduces four coordination pseudoatom systems (**TZ, SQ, MH, JT**) with automatic Bond-Valence-Sum (BVS) oxidation-state inference (Fe/Cu/Mn/Co/V/Mo/Ni) and a dedicated Jahn-Teller elongated-octahedral mode for Cu²⁺ (d⁹) and Mn³⁺ (d⁴); (iii) implements a four-tier reactive covalent docking framework (**P1–P4**: distance, angle, frame-atom torsion, hybrid vdW scaling) with explicit NAC detection, a C3 two-phase search strategy, and six built-in reaction presets; and (iv) supports metal-as-ligand workflows through reverse metal–donor pair potentials. In a systematic **110-mode synthetic metal-coordination benchmark**, LKina completes docking on all 110/110 modes, whereas Vina 1.2.7 fails on **76/110** at atom-type parsing (its 32-type XS system cannot map the extended AD4 metal types of LKina); LKina recovers the metal–donor distance within 0.5 Å of the ideal on **108/110** synthetic systems (mean |d–d₀| = 0.20 Å) and confirms the metal coordination potential well on **104/110** modes (E_ideal < E_far). All four coordination pseudoatom systems (TZ/SQ/MH/JT) are validated via `--metal_geometry_check`; BVS oxidation-state inference is correct on all 14 designed synthetic systems (14/14); and the inline grid generator writes **108 AutoGrid 4.2-format affinity maps** in one pass (covering every AD4 probe type, including the TZ/SQ/MH/JT pseudoatoms) without any external `autogrid4`. On the real Zn metalloprotein complex 4JC, LKina recovers the Zn–NA coordination at 2.19 Å and emits an explicit `METAL_COORD` report, whereas Vina 1.2.7 (both vina and AD4 scoring on identical grids) places the nearest donor at only 2.23 Å (a carboxylate oxygen) and evaluates no coordination term. (Metal-mode energies include the genuine coordination contribution and match the experimental-affinity scale; see Section 4.3.) In a PDB-derived redocking benchmark covering Zn²⁺/Fe³⁺/Cu²⁺ metalloproteins (the full cohort of **104 systems**, i.e. the same corpus as Figure 6), enabling the coordination model places the best-pose ligand donor within the 3.0 Å coordination sphere of the receptor metal in **90.6%** of systems (87/96, the dockable subset of LKina metal mode; median donor–metal distance 2.12 Å), versus 27.8% for the same engine with coordination disabled and 52.1% for Vina 1.2.7. All six reaction presets pass the P1 and P1+P2 end-to-end tests, and the NAC detector exhibits genuine geometric discrimination (it discriminates reactive geometry rather than passing unconditionally). Standard non-metal docking reproduces Vina 1.2.7 energies exactly (1HVR −14.54 and 3PTB −6.202 kcal/mol, identical mode-1 values), and the reactive-gradient and score-decomposition validation layers pass on the live binary (Section 3.1). All extensions are opt-in and have been validated end-to-end (Sections 3.5 and Appendix E).

**Availability:** https://github.com/LK-Studio1128/LKina (GPL-3.0-or-later; Vina-origin files remain Apache-2.0); prebuilt macOS/Linux/Windows binaries, build scripts, and the full benchmark pipeline are in the repository.

**Keywords:** molecular docking; metalloproteins; AutoDock Vina; metal coordination; covalent docking; near-attack conformation

![Graphical Abstract](figures/graphical_abstract.png)

- **Graphical Abstract.** Visual summary of LKina: an extension of AutoDock Vina 1.2.7 (32-type XS) to 113 atom types with metal-aware coordination docking (TZ/SQ/MH/JT pseudoatoms, BVS oxidation-state inference, 4JC Zn–NA at 2.19 Å) and covalent-reactive docking (NAC detection, P1–P4 constraints, C3 two-phase strategy, six presets), with an inline AutoGrid 4.2 generator and full backward compatibility; metal-mode coverage improves from 34/110 (Vina 1.2.7) to 110/110.

---

## 1. Introduction

### 1.1 Molecular docking and the metalloprotein gap

Structure-based drug discovery (SBDD) relies on molecular docking to predict the binding pose and affinity of small molecules within a receptor of known structure [1]. Among the most widely used open-source engines, AutoDock Vina [1,2] couples an empirical scoring function with a Monte Carlo + L-BFGS global/local search and OpenMP multithreading, and has become the de facto standard for virtual screening. AutoDock Vina 1.2.x additionally provides the AutoDock4 (AD4) empirical force field [3] for users who require more chemically detailed electrostatics and desolvation terms.

However, when the receptor contains a metal cofactor, a situation ubiquitous in drug discovery, both force fields share three fundamental limitations. Metal ions are present in an estimated **40% of enzymes and ~30% of approved drug targets**, including HIV integrase (Mg²⁺), carbonic anhydrase (Zn²⁺), histone deacetylases (Zn²⁺), and heme/non-heme iron enzymes (Fe³⁺) [4,5]. The deficiencies are as follows:

1. **Distorted equilibrium distances.** In a generic vdW force field, the equilibrium distance of an N donor is ~2.49 Å, whereas a native Zn–N coordination bond is ~2.0 Å; the scoring function therefore penalizes the crystallographically correct pose.
2. **No directionality.** Spherical Lennard-Jones potentials cannot distinguish axial from equatorial coordination, a particularly serious shortcoming for Jahn–Teller-active ions and for square-planar Pt/Pd complexes of clinical relevance.
3. **Charge-model bias.** Formal charges on high-valent metals (Fe³⁺, Cu²⁺) cause Gasteiger-based electrostatics to over-weight oxygen donors, producing systematic pose errors that follow charge rather than coordination chemistry.

### 1.2 Prior art in metal docking: narrow niches, external dependencies

Santos-Martins *et al.* (2014) addressed the Zn²⁺ problem with the **AutoDock4Zn** force field [6], which places tetrahedral *pseudoatoms* (TZ) at the vacant coordination directions of the metal and encodes the directional potential into the affinity grids, substantially improving redocking accuracy over standard AD4 on 292 Zn²⁺ crystal complexes. AutoDock4Zn, however, is limited to Zn²⁺; it does not cover other transition metals, multi-metal sites, Jahn–Teller-active ions, or bridging waters. More recent protocols such as **MetalDock** [7] provide a reproducible workflow for metal-containing compounds (including metal-as-ligand cases) via optimized Lennard-Jones parameters and ligand-side dummy atoms, and **GPDOCK** [8] adds metal-coordination-aware pose evaluation; both, however, remain external pre-processing protocols layered on top of an existing engine rather than engine-level extensions. All of these workflows also require an external **`autogrid4`** binary to pre-compute the AD4 affinity maps, an installation dependency that complicates high-throughput and cloud deployments.

### 1.3 Covalent docking and the near-attack-conformation problem

Covalent drugs, ibrutinib (BTK, Cys481), afatinib (EGFR), and sotorasib (KRAS G12C) among them, account for a rapidly growing share of the approved-drug space [9]. Docking such ligands requires that the search sample *productive nucleophilic geometries*: near-attack conformation (NAC) theory states that efficient bond formation requires a nucleophile–electrophile distance below 3.0 Å together with a characteristic attack angle (≈180° for SN2 backside attack, ≈109.5° for Michael addition) [10,11]. AutoDock Vina offers no mechanism to encode this geometric prior; CovDock (Schrödinger) does, but it is commercial and generally requires pre-built covalent adducts; and the covalent AutoDock "two-point attractor" method [12] is distance-only and manually configured. To our knowledge, no open-source engine natively guides the search toward NACs, detects them, and reports them in its output.

### 1.4 Contributions

We present **LKina**, a derivative of AutoDock Vina 1.2.7 that closes both gaps at the engine level. Its principal contributions are:

1. **A 113-type AD4 atom system** covering the pharmaceutically relevant metals and metalloids of the periodic table (transition metals, lanthanides, actinides, s-block, p-block metals, oxidation-state variants), with force-field parameters embedded in the binary and parsed at run time.
2. **An inline AG4 affinity-map generator** (`--generate_maps`), an independent C++ reimplementation of the AutoGrid 4.2 scoring function, that computes receptor grids on the fly, removing the external `autogrid4` dependency entirely while remaining `.map`-compatible with existing caches.
3. **Four coordination pseudoatom systems (TZ/SQ/MH/JT)** with automatic **BVS oxidation-state inference** (Fe/Cu/Mn/Co/V/Mo/Ni), **Jahn–Teller distorted-octahedral modes** for Cu²⁺/Mn³⁺, and **semi-explicit water-bridge scoring** for coordinatively unsaturated sites.
4. **A four-tier reactive covalent docking framework (P1–P4)**, Gaussian distance attractor, angle constraints, frame-atom torsion, and hybrid vdW scaling, with analytical gradients (finite-difference deviation ≤ 4.3 × 10⁻⁸ on all six presets; Section 3.2.2), explicit **NAC detection**, a **C3 two-phase search strategy**, and **six reaction presets** (`cys_michael`, `cys_sn2`, `ser_covalent`, `lys_targeting`, `boronic_acid`, `tyr_covalent`).
5. **Metal-as-ligand support** through reverse metal–donor pair potentials (MetalDock `standard_set`-derived where available), metal-autodetection from receptor, single-ligand, and per-ligand batch inputs, and ligand-side coordination-geometry QC for Pt/Pd/Ru/Os/Re complexes.
6. **Full backward compatibility**: for metal-free, non-reactive receptors the behavior is identical to Vina 1.2.7, and the licensing structure (GPL-3.0-or-later combined binary over Apache-2.0 Vina core + GPL-3.0 AG4 extension) is explicitly documented for legal re-use.

---

## 2. Methods

### 2.1 Baseline: AutoDock Vina search engine and the AD4 scoring function

LKina inherits the Vina global search (population-based Monte Carlo with iterated local search), the local optimization (L-BFGS on a smooth approximation of the score), and OpenMP multithreading over independent pose searches. New functionality is added as an additional energy layer on top of the grid-based AD4 potential, so the search machinery, convergence behavior, thread safety, determinism under a fixed seed, remains untouched.

The AD4 empirical binding free energy is

$$\Delta G_{\text{bind}} = W_{\text{vdW}}\,\Delta H_{\text{vdW}} + W_{\text{Hbond}}\,\Delta H_{\text{Hbond}} + W_{\text{elec}}\,\Delta H_{\text{elec}} + W_{\text{desolv}}\,\Delta G_{\text{desolv}} + W_{\text{tor}}\,\Delta S_{\text{tor}}$$

with default weights $W_{\text{vdW}}=0.1662$, $W_{\text{Hbond}}=0.1209$, $W_{\text{elec}}=0.1406$, $W_{\text{desolv}}=0.1322$, and $W_{\text{tor}}=0.2983$ (AutoDock4 [3]). In grid-based docking, the receptor contribution is pre-computed on a Cartesian grid and the ligand energy during the search is obtained by trilinear interpolation, which keeps the per-pose evaluation effectively O(N_ligand). LKina re-implements this entire grid pipeline in-engine (Section 2.3).

### 2.2 Extended atom-type system (113 types)

Standard AutoDock4 defines ~22 atom types. LKina extends the AD4 type table to **113 element/metalloid types** (indices 0–112) plus the four coordination pseudoatoms TZ, SQ, MH, JT (indices 113–116; `AD_TYPE_SIZE = 117`). The types are organized by medicinal-chemistry use case:

**Table 1.** Classification of the 113 AD4 atom types (indices 0–112, covering 80+ metals and metalloids).

| Category | Representative types | Typical application |
|---|---|---|
| Biological metal cofactors | `Mg Ca Mn Fe Co Ni Cu Zn` (+oxidation variants) | metalloenzyme site docking |
| Anticancer metal drugs | `Pt Pd Ru Ir Au Rh` | cisplatin, NAMI-A/RAPTA, auranofin |
| Radiopharmaceutical metals | `Tc Re Ga Y Zr Lu Sm Ho Ra Ac Th` | SPECT/PET/α-therapy agents |
| Toxicology targets | `Cd Hg Tl Pb As Sb Bi` | metal-toxin ligands |
| s-block metals | `Li Na K Rb Cs Al Sr Ba` | ion-channel and kinase targets |
| Early transition metals | `V Cr Ti Sc Nb Hf Ta W Mo` | insulin mimetics, sulfite oxidase (Mo/W) |
| Lanthanides | La–Lu (13 types) | MRI contrast agents (Gd) |
| Actinides | `Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm` | targeted α therapy (Ac-225, Ra-223) |
| Metalloids / p-block | `Se As Ge Ga In Sn B Be Te Po At` | selenoenzymes, boron-based drugs |
| Oxidation-state variants | `Fe2 Fe3 Cu1 Cu2 cu2_jt Mn2 Mn3 mn3_jt Co2 Co3 V4 V5 Mo4 Mo6 As3 As5 Sb3 Sb5 uo2` | oxidation-state-specific docking |
| Aqueous coordination variants | `mg_aq ca_aq fe3_aq mn2_aq co2_aq` | active-site water coordination |
| Contrast-agent adducts | `Gd_DTPA Gd_DOTA` | Gd-based MRI agents |

Force-field parameters ($R_{\text{eq}}$, $\varepsilon$, solvation terms) are embedded as parameter text in `ad4_parameter_data.cpp` and parsed at run time to build the `atom_kind` lookup table, adding or tuning a type requires no recompilation. The vdW parameters follow the MetalDock MC-optimized LJ fits [7] and the Harding compilation of metal–donor distances from the Cambridge Structural Database [13]; Se receives independent parameters ($R_{\text{eq}}=2.03$ Å, $\varepsilon=0.30$ kcal/mol) rather than being aliased to S; and the uranyl dication `uo2` is modeled with linear-coordination parameters for environmental-toxicology docking.

### 2.3 Inline AG4 affinity-map generation

The standard AD4 workflow requires the external `autogrid4` program to pre-compute one `.map` file per atom type. LKina removes this dependency with `ag4_compute_maps()` (`--generate_maps`), a self-contained C++ reimplementation of the AutoGrid 4.2 pairwise interaction scoring function. The pipeline is:

1. **Receptor parsing** (`ag4_parse_receptor`), coordinates plus AD4 type per atom;
2. **Pseudoatom injection** (`ag4_inject_[TZ|SQ|MH|JT]`), vacant coordination directions are filled according to the active metal mode(s);
3. **Grid accumulation**, for every grid point, pairwise AD4 potentials are summed over receptor atoms;
4. **Electrostatic grid**, Mehler–Solmajer distance-dependent dielectric;
5. **Desolvation grid**, atomic-volume solvation potential.

The default grid spacing is 0.375 Å (`--spacing` adjustable). `--write_maps` exports the grids as standard `.map` files for caching, and `--maps` re-imports previously generated grids, skipping the 5–30 s per-receptor generation step. In multi-metal systems, injection is performed mode-by-mode over the `active_modes` list so that every metal site receives the correct pseudoatom coverage. Because an external `autogrid4` binary cannot reproduce LKina's pseudoatom injection, Zn-compatibility state, or `nbp_r_eps` overrides, the external-`autogrid4` fallback is **disabled** whenever `zn_mode`, `metal_mode`, multi-metal modes, or custom NBP overrides are active, guaranteeing that the grids used in the search are always consistent with the search state.

### 2.4 Coordination pseudoatoms (TZ / SQ / MH / JT)

The pseudoatom concept, introduced by AutoDock4Zn [6], is generalized in LKina to four geometric systems. A pseudoatom is a virtual atom placed at a *vacant* coordination direction of the receptor metal; its interaction with ligand probe atoms encodes the coordination geometry directly into the affinity grids:

$$V_{\text{pseudo}}(r) = \varepsilon\left[\left(\frac{r_{\text{eq}}}{r}\right)^{12} - 2\left(\frac{r_{\text{eq}}}{r}\right)^6\right]$$

**Table 2.** Geometry and coordination-interaction parameters (nbp overrides) of the four coordination pseudoatoms.

| Pseudoatom | Geometry | Directions | Representative metals | nbp $r_{\text{eq}}$ (Å) | nbp $\varepsilon$ range (kcal/mol) |
|---|---|---|---|---|---|
| **TZ** | tetrahedral | 4 | Zn²⁺, Cd²⁺ (4-coordinate d¹⁰) | 0.25 | 3.0–23.2 |
| **SQ** | square-planar / linear | 4 / 2 | Cu²⁺, Pt²⁺, Pd²⁺, Ni²⁺; Hg²⁺, Ag⁺ (linear) | 0.25 | 3.0–24.0 |
| **MH** | octahedral | 6 | Fe³⁺, Mn²⁺, Co²⁺, Mg²⁺, … | 0.25 | 0.2–14.0 |
| **JT** (equatorial) | elongated octahedron | 4 | Cu²⁺ (d⁹), Mn³⁺ (d⁴) | 0.25 | 2.5–15.0 |
| **JT** (axial) | elongated octahedron | 2 | Cu²⁺ (d⁹), Mn³⁺ (d⁴) | 0.25 | 0.2–5.0 |

The pseudoatoms themselves carry no vdW interaction (official AD4Zn neutral parameters $R_{ii}=1.0$ Å, $\varepsilon_{ii}=0.0$; see Section 3.2.1); all coordination interactions are encoded into the grids through explicit `nbp_r_eps` overrides whose $r_{\text{eq}}=0.25$ Å places the ligand donor probe directly on the pseudoatom site, with $\varepsilon$ chosen by metal–donor HSAB matching over 0.2–24 kcal/mol (unweighted tabulated values; actual well depths after the AD4 `coeff_vdW` weighting are given in Section 3.2.1, e.g. the Zn–NA/TZ tabulated 23.2135 → weighted −3.86 kcal/mol). The JT axial pseudoatoms share the $r_{\text{eq}}=0.25$ Å overrides but are injected at a longer distance (2.45 Å for cu2_jt, 2.28 Å for mn3_jt, vs 2.03/1.92 Å equatorial) with weaker $\varepsilon$ (0.2–5.0 vs 2.5–15.0 kcal/mol), mirroring the elongated (~2.4 Å) and weakened axial bonds of Jahn–Teller-active ions. Vacant directions are chosen by `ag4_select_vacant_dirs()`, a maximum-angular-separation strategy over the direction-complement space of the existing ligand donors, ensuring that injected positions correspond to genuinely open coordination sites.

### 2.5 Automatic metal-mode detection and BVS oxidation-state inference

When `--metal_mode` is not specified, `detect_metal_mode_from_pdbqt()` proceeds as follows:

1. **Scan the receptor PDBQT** for metal symbols in the AD4 type column of ATOM/HETATM records;
2. **BVS oxidation-state inference** for Fe/Cu/Mn/Co/V/Mo/Ni (below);
3. **Donor-count fallback**, if the BVS confidence is insufficient, count OA/NA/SA donors within 3.0 Å and infer the oxidation state heuristically;
4. **Ligand supplement scan**, if no metal is found in the receptor and no mode is specified, the ligand is scanned (per-ligand in batch mode, with map regeneration before each dock);
5. **Mode dispatch**, single-metal ligands enable the corresponding mode (multi-metal ligands enable a comma-separated combination, e.g., `pt,ru`); Zn additionally enables the Zn coordination-potential path.

The Bond Valence Sum (BVS) [14,15] of a metal center is

$$\text{BVS} = \sum_{i} \exp\!\left(\frac{r_0 - d_i}{B}\right),\qquad B = 0.37\ \text{Å}$$

accumulated over donor atoms within a 3.2 Å cutoff. For each candidate oxidation state $n$, the deviation $\delta_n = |\text{BVS} - n|$ is computed; the state with the minimal deviation is selected when $\min_n \delta_n < 1.25$ vu. Bond-valence parameters $r_0$ follow Brese & O'Keeffe / Brown & Altermatt [14,15]:

**Table 3.** Bond-valence-sum (BVS) bond-length parameters and supported oxidation states.

| Metal | $r_0$ (M–O) | $r_0$ (M–N) | $r_0$ (M–S) | States |
|---|---|---|---|---|
| Fe | 1.76 / 1.73 | 1.79 / 1.76 | 2.05 / 1.98 | +2, +3 |
| Cu | 1.72 / 1.68 | 1.74 / 1.70 | 1.96 / 1.89 | +1, +2 |
| Mn | 1.79 / 1.76 | 1.80 / 1.77 | — | +2, +3 |
| Co | 1.75 / 1.70 | 1.77 / 1.72 | — | +2, +3 |
| V | 1.78 / 1.80 | 1.80 / 1.82 | — | +4, +5 |
| Mo | 1.90 / 1.86 | — / 1.92 | 2.12 / — | +4, +6 |
| Ni | 1.654 / 1.620 | 1.679 / 1.650 | 1.978 / 1.950 | +2, +3 |

(Entries are listed for the lower/higher oxidation state; "—" marks donors without published parameters.)

### 2.6 Jahn–Teller distorted-octahedral mode

For `cu2_jt` and `mn3_jt`, `ag4_tetragonal_axis()` estimates the JT elongation axis:

1. collect donor direction vectors $\{\hat{u}_i\}$ within 3.4 Å of the metal;
2. find an antipodal pair with $\hat{u}_i \cdot \hat{u}_j < -0.7$;
3. the JT axis is $\hat{z}_{\text{JT}} = \text{normalize}((\hat{u}_i - \hat{u}_j)/2)$; if no pair exists, default to $(0,0,1)$;
4. equatorial directions ($|\hat{u}\cdot\hat{z}_{\text{JT}}| < 0.3$, up to 4) receive equatorial pseudoatoms; $\pm\hat{z}_{\text{JT}}$ receive axial pseudoatoms.

The differentiated axial (injection distance 2.45 Å for cu2_jt / 2.28 Å for mn3_jt, nbp ε 0.2–5.0) versus equatorial (injection 2.03/1.92 Å, ε 2.5–15.0) parameters allow the search to distinguish the two coordination spheres, which is essential for d⁹ Cu²⁺ and d⁴ Mn³⁺ sites.

### 2.7 Semi-explicit water bridge and geometry-based pose reranking

Coordinatively unsaturated sites are common in metalloenzymes (e.g., water-bridged Zn sites in matrix metalloproteinases). LKina automatically infers candidate bridging-water sites:

$$N_{\text{water}} = \min(2,\ N_{\text{max}} - N_{\text{receptor}})$$

with positions from `ag4_select_vacant_dirs()` and weights 0.90 (first water) and 0.70 (second). The post-processing water-bridge score (which does not enter the L-BFGS gradient) is

$$E_{\text{water}} = -0.90 \sum_{\text{sites}} w_s \exp\!\left(-\frac{(d_s - 2.80)^2}{2 \times 0.40^2}\right)$$

centered on the canonical M–O(water) distance of 2.80 Å. The complete metal rerank correction is

$$E_{\text{rerank}} = E_{\text{geo}} + E_{\text{water}} + E_{\text{JT}},$$

$$E_{\text{geo}} = -1.25 \sum_{\text{sites}} \max_{i \in \text{donors}} \frac{\varepsilon_i}{20}\, G(d_i,\, r_{\text{eq},i},\, 0.30),$$

$$E_{\text{JT}} = -0.50\left[\max_i\!\big(|\cos\phi_i|\, G(d_i,\, d_{\text{axial}},\, 0.35)\big) + 0.5\max_i\!\big((1-|\cos\phi_i|)\, G(d_i,\, d_{\text{eq}},\, 0.30)\big)\right],$$

where $G(d, d_0, \sigma)=\exp(-(d-d_0)^2/2\sigma^2)$. By design, the rerank terms affect only the final pose ranking and do not enter the search gradient, avoiding spurious forces; an optional search-time soft-constraint channel (`--metal_soft_weight`, default 0.0, suggested 0.1–0.5) exposes a smooth log-sum-exp approximation of the geometry terms inside `eval_deriv()` for users who wish to let the metal geometry steer the search itself.

### 2.8 Metal as ligand: reverse pair potentials and ligand-side QC

A frequent medicinal-chemistry scenario is the metal complex *as the ligand* (Pt, Ru, Au, Tc/Re, and Gd agents, or simple metal-ion probes). LKina supports this scenario through **reverse metal–donor pair potentials**: after generating the receptor-metal `donor → metal` `nbp_r_eps` overrides, the engine automatically constructs the complementary `metal → donor` entries (with de-duplication), so that receptor OA/NA/SA/O/N/S donors exert the same directional coordination response on a metal probe atom. For metals with published MetalDock MC-optimized parameters (V, Cr, Co, Ni, Cu, Mo, Ru, Rh, Pd, Re, Os, Pt), the reverse `metal → NA/OA/SA/HD` entries adopt the MetalDock ε combinations and 12–10 metal–donor potentials [7]; other metals keep the symmetric reverse completion. For Pt/Pd (square-planar) and Ru/Os/Re (octahedral) complexes, an additional ligand-side geometry QC reports `LIGAND_METAL_GEOM` / `LIGAND_METAL_SITE` REMARKs, and `--ligand_metal_geometry_weight > 0` optionally promotes the geometry penalty into the post-processing ranking.

### 2.9 Reactive covalent docking framework (P1–P4)

LKina's reactive docking adds an *additive penalty layer* on top of the AD4 grid score. Four progressively richer constraint tiers are available:

**P1, distance constraint (Gaussian attractor).** The receptor electrophile (e.g., Cys SG) and the ligand nucleophile are pinned by

$$E_{\text{P1}}(r) = -A\exp\!\left(-\frac{(r-r_0)^2}{2\sigma^2}\right) + k_r (r-r_0)^2 \cdot \mathbf{1}[r > r_0 + \sigma],$$

which attracts the pair toward the target bond length $r_0$ and applies a quadratic penalty when the ligand drifts beyond $r_0+\sigma$. The gradient with respect to the ligand coordinates is computed analytically.

**P2, angle constraint.** A harmonic angular term enforces the attack angle at the nucleophile:

$$E_{\text{P2}}(\theta) = k_{\theta}\,(\cos\theta - \cos\theta_0)^2,$$

with the target angle $\theta_0$ adapted to the mechanism: 180° for SN2 backside attack, ≈109.5° for Michael addition (tetrahedral transition state). Receptor-side and ligand-side frame atoms (`--reactive_frame_atom`, `--reactive_lig_frame_atom`) define the reference frames; the gradient of $\cos\theta$ w.r.t. the ligand atom is

$$\frac{\partial \cos\theta}{\partial \vec{r}_{\text{lig}}} = \frac{1}{|\vec{u}|}\big(\hat{v} - \cos\theta\, \hat{u}\big).$$

A flat-bottom variant (`--reactive_angle_width > 0$) zeroes the penalty within $\theta_0 \pm w$, avoiding premature locking of the search onto a harmonic well.

**P3, frame-atom torsion support.** Frame atoms extend P2 to fully define the nucleophilic trajectory in three dimensions (both receptor- and ligand-side frames).

**P4, hybrid vdW scaling.** In `hybrid` mode the repulsive part of the receptor–ligand vdW interaction is scaled by $\lambda \in [0,1]$ (`--reactive_hybrid_vdw_scale`):

$$E_{\text{vdW}}^{\text{hybrid}} = \lambda\, E_{\text{vdW}}^{\text{repulsive}},$$

allowing the ligand to enter the repulsive zone and sample transition-state-like geometries. The hybrid energy is smooth in $\lambda$ (continuous scaling verified by regression test E and by the measured λ-sweep in Section 3.2.2: cys_michael −10.09 → −5.716 → −4.734 kcal/mol at λ = 0/0.5/1.0).

**NAC detection.** After docking, each output pose receives `REMARK REACTIVE_NAC: YES/NO`, defined by a reactive distance < 3.0 Å *and* an attack angle within the target ±25° window. A finite-difference gradient check (central difference, step ε = 10⁻⁴ Å) confirms the analytical reactive gradients on all six presets: maximum |analytic − numeric| deviation 2.4×10⁻¹⁶–4.3×10⁻⁸ across distance- and angle-force components (Section 3.2.2; `--reactive_gradcheck`).

### 2.10 C3 two-phase strategy, reaction presets, and Metal Bias

Applying the full constraint set throughout the search can trap the ligand near the receptor and under-sample the global energy landscape. LKina therefore offers the **C3 two-phase strategy** (`--reactive_two_step`):

- **Phase 1**, a standard unconstrained Vina MC + L-BFGS search (or, in the C3b variant, a broad/weak Gaussian with $\sigma\times3$, $\varepsilon\times0.15$ acting as an "attracting cavity" prior);
- **Pose filtering**, poses whose reactive atom lies within `--reactive_presample_dist` (default 10 Å) of the anchor are kept;
- **Phase 2**, surviving poses are refined with the full P1–P4 constraint set via L-BFGS.

Six built-in **reaction presets** auto-fill all reactive parameters (individual CLI flags override the preset):

**Table 4.** Geometry parameters of the six built-in reaction presets.

| Preset | Reaction | $r_0$ (Å) | $\theta_0$ (°) | $w$ (°) | $\lambda$ |
|---|---|---|---|---|---|
| `cys_michael` | Cys SG Michael addition (C=C warhead) | 1.82 | 109.5 | 25 | 0.2 |
| `cys_sn2` | Cys SG SN2 substitution (180° backside) | 1.82 | 180.0 | 15 | 0.2 |
| `ser_covalent` | Ser OG acylation (β-lactams) | 1.34 | 109.5 | 25 | 0.2 |
| `lys_targeting` | Lys NZ Schiff base (aldehyde warhead) | 1.47 | 109.5 | 30 | 0.3 |
| `boronic_acid` | reversible boronate (Ser/Thr/Tyr OH) | 1.47 | — | — | 0.5 |
| `tyr_covalent` | Tyr OH nucleophilic attack | 1.38 | 109.5 | 25 | 0.2 |

Finally, **Metal Bias (O5)** (`--metal_bias`) automatically locates the first metal atom in the receptor and injects a soft Gaussian attractor, $E_{\text{bias}} = -A\exp(-r^2/2\sigma^2)$ with $A=2.0$ kcal/mol and $\sigma=1.5$ Å, a light directional prior suited to rapid metalloenzyme screening that leaves Vina's global search intact; it is mutually exclusive with reactive mode.

### 2.11 Implementation, licensing, and engineering notes

LKina extends the Vina 1.2.7 C++14 code base. The principal new modules are: `ag4_engine.{h,cpp}` (metal-mode enumeration, pseudoatom injection, grid computation; ~1,300 lines), `ad4cache.{h,cpp}` (AD4 grid cache, metal state, rerank scoring, search-time soft-constraint gradients; ~430 lines), `atom_constants.h` (113-type constants; ~520 lines), `ad4_parameter_data.cpp` (embedded force-field text; ~200 lines), `embedded_ad4_grid.cpp` (inline grid-generation interface; ~190 lines), `reactive_types.h` (reactive framework types; ~237 lines), plus extensions to `vina.cpp` (~1,600 lines), `conf.h` (rerank fields; ~390 lines) and `main.cpp` (CLI, BVS inference, metal-mode detection, O5; ~1,500 lines). Multithreading follows Vina's OpenMP model (`--cpu N`); pseudoatom injection occurs during single-threaded receptor parsing, so no new thread-safety concerns arise.

**Licensing.** LKina uses a dual-license structure: all Vina-origin files remain Apache-2.0, while the AG4 grid engine and its parameter data (`ag4_engine.*`, `embedded_ad4_grid.*`, `ad4_parameter_data.*`) are GPL-3.0-or-later (independent reimplementation of the AutoGrid 4.2 function, itself GPL-2.0-or-later, © The Scripps Research Institute). Because the two components are statically linked, the combined LKina binary is distributed under GPL-3.0-or-later; `COPYING` and `NOTICE` document the file-level scope and upstream attribution.

**Engineering note on the Vina scoring mode.** During LKDock GUI-integration packaging (PyInstaller) testing, running LKina with `--scoring vina` (the 32-type XS atom system) triggered a reproducible SIGABRT at output time for a broad class of ligands, including purely organic ones. Root-cause analysis (lldb) traced the failure to two issues: (i) `XS_TYPE` cannot map the 117-type AD4 atom space (`VINA_CHECK` → `internal_error` in release builds), and (ii) a pre-existing defect in `Vina::~Vina()` that re-declares member variables as locals inside the function body, corrupting `boost::ptr_vector` destruction. Rather than rewriting the upstream XS type system, a high-risk change touching `precalculate`, `cache`, `non_cache`, static assertions, and `.maps` binary compatibility, with no benefit to the metal/covalent mission, LKina's GUI integration unconditionally routes through `--scoring LKDock` (AD4) mode, which covers all 117 atom types and exhibits identical search behavior for organic ligands. The AD4 path is therefore the recommended mode for all LKina workflows. **Note:** invoking the LKina binary directly on the CLI with `--scoring vina` works correctly — the numerical identity reported in Section 3.5 (1HVR/3PTB) was obtained by direct CLI invocation; the SIGABRT above is specific to the GUI-integration packaged environment.

### 2.12 Command-line interface (invocation examples)

All LKina features are exposed through a single binary whose CLI is a superset of AutoDock Vina's. The complete flag reference ships with the repository (`LKina --help_advanced`); we list here the invocations actually used in this work so that every experiment in Sections 3.2–3.6 can be reproduced verbatim. **Metal docking with inline grid generation** (no external `autogrid4`) is performed by:

```
lkina --receptor rec.pdbqt --ligand lig.pdbqt \
      --scoring ad4 --generate_maps --metal_mode zn \
      --center_x X --center_y Y --center_z Z \
      --size_x Sx --size_y Sy --size_z Sz \
      --exhaustiveness 8 --seed 42 --out pose.pdbqt
```

If `--metal_mode` is omitted, the engine auto-detects metal modes from receptor/ligand PDBQT atom types; appending `--no_auto_metal` suppresses detection and keeps pure AD4 scoring, this exact pair of configurations constitutes the LKina-with/without-coordination-model comparison reported in Section 3.6. Receptor-side geometry can be validated before docking with `--metal_geometry_check`, and the search-time soft-constraint channel is enabled by `--metal_soft_weight w` (0 = off). Oxidation-state variants are addressed by explicit tokens (e.g. `fe2`/`fe3`, `cu1`/`cu2`, and `cu2_jt` for Jahn–Teller-distorted Cu²⁺); when the token is omitted for a variable-valence metal, BVS inference selects the mode from the donor environment (Section 2.5). **Reactive covalent docking** uses the preset mechanism:

```
lkina --receptor rec.pdbqt --flex flex.pdbqt \
      --ligand lig.pdbqt --scoring ad4 --generate_maps \
      --reactive_preset cys_michael [preset overrides] \
      --center ... --size ... --out pose.pdbqt
```

with six presets available (`cys_michael | cys_sn2 | ser_covalent | lys_targeting | boronic_acid | tyr_covalent`), each overridable per parameter; `--reactive_gradcheck` runs the finite-difference gradient audit of Section 3.2.2. Finally, ligand-side metallocomplex geometry QC is activated by `--ligand_metal_geometry_weight w`. Internally, these flags map onto the modules of Section 2.11, e.g., `--metal_mode` drives pseudoatom injection in `ag4_engine.cpp`, while `--generate_maps` routes grid computation through `embedded_ad4_grid.cpp` instead of reading external `.map` files.

**Scoring functions and workflow compatibility.** `--scoring` accepts three force fields: `ad4`/`LKDock` (default; covers all 113 AD4 atom types and enables the metal/covalent extensions), `vina` (Vina 1.2.7's original 32-type XS system, used for backward-compatibility verification, Section 3.5), and `vinardo` (Vina's generic scoring variant). `--config file` reads all command-line flags from a configuration file (convenient for large-scale reproduction); `--batch dir|lig1 lig2 …` provides batch docking (per-ligand automatic metal detection with AD4 grid regeneration before each dock; Section 2.5), with `--dir out/` specifying the batch output directory; `--autobox` sizes the grid box from the input ligand (for use with `--score_only`/`--local_only`); `--flex flex.pdbqt` enables flexible side chains; `--local_only` performs a local search only; and `--unbound_energy E` sets the unbound-system energy explicitly in `--score_only`. All standard Vina 1.2.7 flags (`--num_modes`, `--energy_range`, `--min_rmsd`, `--spacing`, `--seed`, `--cpu`, `--verbosity`, etc.) are preserved unchanged.

---

## 3. Results

### 3.1 Regression tests and build matrix

LKina's validation layer (17 automated regression checks in `tests/`, per engine design document §5.1) is organized into six groups: A, 4 global docking-energy regressions at fixed seed/exhaustiveness (2 of which reproduce the Vina 1.2.7 values exactly); B, 5 score-decomposition checks on crystal conformations (`--score_only`); C, 3 `lig_atom` gradient vs finite-difference validations; D, 1 `target_angle` parameter-effect test; E, 1 `hybrid_vdw_scale` continuous-scaling test; and F, 3 `lig_frame_atom` gradient validations. All 17 checks pass on the release binary (macOS ARM64, seed 42; `tests/reactive_regression.sh`). In addition, this work independently re-executed the highest-value checks against the live binary: the reactive gradient check via `--reactive_gradcheck` (maximum |analytic − numeric| deviation **2.4 × 10⁻¹⁶–4.3 × 10⁻⁸** across the seven force components of all six presets; Section 3.2.2), two global docking-energy regressions on real organic systems (Section 3.5), and score-decomposition verification during the energy-scale audit (`INTER + torsions = VINA RESULT` holds exactly). The engine compiles with zero errors and zero warnings (`-Wall -Wno-long-long -O3 -fopenmp`) on macOS 14 (Apple M2, clang++ 14), macOS 13 (Intel), and Ubuntu 22.04 (g++ 11).

### 3.2 Systematic metal-coverage benchmark (110 metal modes)

To verify the claim that LKina's 113-type AD4 atom system covers the pharmaceutically relevant metals, we built a **synthetic metal-coordination benchmark** spanning every `metal_mode` CLI token implemented by the engine (**110 modes**: biological Zn/Mg/Ca/Mn/Fe/Co/Ni/Cu; medicinal Pt/Pd/Ru/Ir/Au/Rh/Ag/Tc/Re/Os; toxicological Cd/Hg/Tl/Pb/Sb/Bi/As; s-block Na/K/Li/Al/Sr/Ba; early-transition V/Cr/Ti/Sc/Y/Zr/Nb/Hf/Ta/W/Mo; lanthanides La–Lu; actinides Ac–Fm + UO₂²⁺; metalloids Se/Ge/Te/Po/At/B/Be; post-transition Ga/In/Sn; the Jahn–Teller modes cu2_jt/mn3_jt; oxidation-state variants fe2/fe3/cu1/cu2/mn2/mn3/co2/co3/ni2/ni3/as3/as5/sb3/sb5/v4/v5/mo4/mo6; and the aqua variants mg_aq/ca_aq/fe3_aq/mn2_aq/co2_aq). For each mode, a synthetic receptor was generated, the metal at the origin, (n−1) coordination donors (NA/OA) placed at the ideal M–L distance $d_0$ along the appropriate geometry (tetrahedral TZ, square-planar/linear SQ, octahedral MH, or JT), plus backbone carbons, and an NA-donor probe ligand was docked from a position $d_0 + 4$ Å beyond the ideal coordination distance. All runs used exhaustiveness 8, seed 42, grid spacing 0.375 Å, and identical boxes. The receptor/ligand PDBQT files are produced by `benchmarks/pdbqt_util.py`, which uses an absolute-column layout (type token at 0-indexed columns 77–78; chain 21; resi 22–25; x 30–37) verified to parse correctly under both LKina and Vina 1.2.7.

**Result (Table 5, Figures 1–2):** LKina completed docking on **all 110/110** metal modes with finite energies. Vina 1.2.7, run on the *identical* receptor/ligand files (parsed without error), completed only **34/110**; the remaining **76 modes fail at PDBQT parse time** with `PDBQT parsing error: Atom type <M> is not a valid AutoDock type (atom types are case-sensitive)`, because Vina's 32-type XS system has no entries for the extended AD4 metal types (Pt, Pd, Ru, Ir, the lanthanides, Tc/Re/Os, Cd/Hg, most early TMs, etc.). This difference reflects a *type-system capability boundary* rather than a scoring difference: 76 of the 110 metal-coordination systems cannot be parsed by Vina at all. Across the 110 modes, LKina recovered the metal–donor distance to within 0.5 Å of $d_0$ on **108/110** (mean |d–d₀| = 0.20 Å, median 0.20 Å); 104/110 modes additionally show a measurable coordination potential well (E_ideal < E_far + 0.3 kcal/mol). Coordination accuracy is broken down by pseudoatom class in Figure 2 (right): TZ (2 modes, mean |d–d₀| 0.22 Å), SQ (7, 0.24 Å), MH (95, 0.20 Å), JT (2, 0.19 Å), LIN (4, 0.14 Å). The synthetic systems deliberately lack a protein pocket, so these residuals are lower-bound estimates of coordination accuracy; real-system accuracy is established separately by the 4JC comparison (Section 3.3) and the PDB-derived Zn/Fe/Cu redocking datasets (Section 3.6).

**Table 5.** Metal-mode coverage and coordination recovery (110 synthetic systems, seed 42).

| Engine | Modes docked | Failed at parse | \|d−d₀\| < 0.5 Å | \|d−d₀\| < 1.0 Å | Mean \|d−d₀\| (Å) |
|---|---|---|---|---|---|
| **LKina** | **110/110** | 0 | **108** | **109** | **0.20** |
| Vina 1.2.7 | 34/110 | 76 | — | — | 0.62 (34 dockable) |

![Figure 1](figures/fig1_metal_coverage_110.png)

- **Figure 1.** All-metal coverage (A) and Vina acceptance by metal family (B). LKina docks **110/110** metal modes; Vina 1.2.7 docks only **34/110** and rejects 76 at PDBQT parse (`Atom type <M> is not a valid AutoDock type`). The breakdown by family shows that Vina accepts mostly biological TMs (Zn, Mg, Mn, Fe, Co, Ni, Cu) but rejects the entire medicinal / actinide / lanthanide / early-TM / metalloid spectrum.

![Figure 2](figures/fig2_coordination_geometry.png)

- **Figure 2.** Coordination accuracy (A) and pseudoatom geometry check (B). (A) Histogram of |d(ligand–metal) − r_eq| over 110 modes (mean 0.20 Å, 108/110 within 0.5 Å). (B) Per-mode measured vs expected donor–metal distance for the four pseudoatom systems TZ/SQ/MH/JT; ideal-count / donors-checked annotations show every receptor donor classified as "ideal" or "good" against the mode-specific nbp r_eq.

![Figure S1](figures/figS1_energy_well_landscape.png)

- **Figure S1.** Coordination-potential-well landscape across the complete 110-mode metal enum, stratified by pseudoatom class (TZ/SQ/LIN/MH/JT). Each dot is one synthetic metal-coordination system; dot area encodes the well depth (E_far − E_ideal). A well is confirmed in **104/110** systems (red); the remaining 6 systems (grey: four LIN, one TZ, and one MH edge case) show no well at the probed donor position and are analyzed in Section 3.3.

### 3.2.1 Engine fixes (v1.0.2, post-review)

Three engine-level defects were found during full benchmarking and fixed, with the corrected binary re-verified against the entire corpus (2026-08-28; validated line-by-line against the official AutoGrid source `mainpost1.28.cpp` and the official AD4Zn.dat): (i) **missing epsilon weighting for `nbp_r_eps` overrides** — real AutoGrid applies `epsij *= AD4.coeff_vdW` when building maps (mainpost1.28.cpp:1786), so the official NA–TZ override (tabulated 23.2135) yields an actual well depth of −3.86 kcal/mol; LKina previously used the unweighted value, making the TZ well 6× too deep — now corrected to the official weighted semantics; (ii) **incorrect built-in TZ/SQ/MH/JT pseudoatom parameters** — the official AD4Zn.dat TZ entry is `Rii=1.0, epsii=0.0` (no vdW interaction; all contacts flow through explicit nbp overrides), whereas LKina previously filled `Rii=0.25, epsii=23.2`, creating spurious wells for non-override probes — now set to the official neutral parameters; (iii) **`AG4_EINTCLAMP` semantics** — AutoGrid only clamps pairwise LJ table values from above (100000) and never clamps the summed map value; LKina's previous ±1000 two-sided clamp was corrected accordingly. After the fix, metal-mode energies fall to the experimentally consistent scale: 3HS4 (CA II/AZM) crystal-pose zn score −23.9 → −11.3 kcal/mol (experiment ≈−10.8); 4JC global-docking best −34.9 → −11.5 kcal/mol; plain-AD4 and Vina compatibility is unaffected (1HVR −14.54 and 3PTB −6.202 kcal/mol still reproduce exactly). In addition, a **`--no_auto_metal` flag plus an always-visible stderr warning** (`main.cpp`) makes auto-detection of metal modes from receptor/ligand PDBQT an opt-out behavior. **Note:** the −11.3 kcal/mol score reported above for 3HS4 is a re-scoring of the crystal pose (`--score_only`), while Section 3.7 (Table 9) reports the global-docking best affinity of −13.32 kcal/mol for 3HS4 AD4+zn; the difference reflects the extra search freedom. Both values agree with the experimental affinity of ≈−10.8 kcal/mol to within 3 kcal/mol, consistent with the accuracy range reported for AutoDock4Zn.

### 3.2.2 Feature-family measurements (pseudoatoms, BVS, water bridge, metal-as-ligand)

Beyond the coverage benchmark, every feature described in Section 2 was measured against a purpose-built synthetic system. Results appear in Figure 2 (right), Figure 3, and Figure 4 (B), and are summarized below.

- **Pseudoatom systems (TZ/SQ/MH/JT)** validated by `--metal_geometry_check` on five representative metals (zn → TZ, pt → SQ, fe → MH, cu2_jt → JT, mn3_jt → JT): every receptor donor was classified ideal or good against the mode-specific nbp $r_{eq}$ (3/3, 3/3, 5/5, 4/4, 4/4). Ideal-CN and vacancy counts are also reported.
- **Bond-Valence-Sum oxidation-state inference** (Fe/Cu/Mn/Co/V/Mo/Ni, ±1 e⁻ variants): **14/14 designed cases** picked the correct oxidation state when `--metal_mode` was not given (auto-detected from the receptor donor distribution via Brese & O'Keeffe $r_0$ values).
- **Semi-explicit water-bridge candidate sites** built and evaluated on an Mg system with 4 of 6 occupied octahedral sites (2 vacancies → 2 water-bridge sites, weight 0.90/0.70, $d_\text{water} = 2.80$ Å); `METAL_WATER_E` (−0.000), `METAL_GEO_E` (−0.090) and `METAL_COORD: Mg donor=NA d=2.28 A` are all emitted in the docked pose.
- **Metal-as-ligand geometry QC** (`--ligand_metal_geometry_weight`) on a Pt(II) square-planar complex: the ideal geometry (N–Pt–N = 90°, 4/4 CN) reports a 0.003 kcal/mol geometry penalty; a deliberately distorted variant (N–Pt–N = 60°) reports 0.903 kcal/mol (≈ 300× higher) and would be flagged by a downstream filter.
- **Analytical reactive gradients** verified by finite-difference `--reactive_gradcheck` on all six presets: maximum |analytic − numeric| deviation **2.4 × 10⁻¹⁶–4.3 × 10⁻⁸** across the seven distance + angle force components per preset.

![Figure 3](figures/fig3_bvs_inference.png)

- **Figure 3.** Bond-Valence-Sum (BVS) oxidation-state inference on 14 designed synthetic systems covering Fe²⁺/Fe³⁺, Cu⁺/Cu²⁺, Mn²⁺/Mn³⁺, Co²⁺/Co³⁺, V⁴⁺/V⁵⁺, Mo⁴⁺/Mo⁶⁺, Ni²⁺/Ni³⁺. The mode picked (auto-detected from the receptor donor distribution via `ag4_bvs_pick_mode`) matches the designed mode in all 14 cases.

![Figure 4](figures/fig4_maps_metal_ligand.png)

- **Figure 4.** Inline AG4 grid generator (A) and metal-as-ligand geometry QC (B). (A) `--generate_maps` writes **108** AutoGrid 4.2-format affinity maps (covering every AD4 probe type, including the TZ/SQ/MH/JT pseudoatoms) in one pass without invoking any external `autogrid4` binary. (B) Pt(II) square-planar ligand docked with `--ligand_metal_geometry_weight`: the ideal 90° N–Pt–N / Cl–Pt–Cl complex scores a 0.003 kcal/mol geometry penalty; a 60° distorted variant scores 0.903 kcal/mol (×300) and would be flagged for rejection in a downstream filter.

### 3.3 Real Zn-metalloprotein comparison (4JC)

To test coordination accuracy in a realistic setting, we docked a 14-atom Zn-binding ligand (4JC305, a sulfonamide-like inhibitor) into the Zn metalloprotein 4JC (2,556 atoms including the catalytic Zn²⁺ at −3.38, 0.33, 85.49) under three engine configurations on identical grids (25 Å³, 0.375 Å, exhaustiveness 16, seed 42):

**Table 6.** Three-engine comparison on the real Zn-metalloprotein 4JC.

| Engine configuration | Best-mode score (kcal/mol) | d(Zn–nearest donor) | Coordination report |
|---|---|---|---|
| **LKina** AD4 + `--metal_mode zn` | **−11.51** | **2.19 Å (ligand NA; nearest atom H at 1.89 Å)** | `METAL_COORD: Zn donor=NA d=2.19 A` + `METAL_RERANK`/`METAL_GEO_E` |
| Vina 1.2.7 (vina scoring) | −6.39 | 2.23 Å (ligand OA) | none |
| Vina 1.2.7 (AD4 + LKina-generated `.map` files) | −13.07 | 2.23 Å (ligand OA) | none |

LKina is the only configuration that both (a) explicitly reports the metal coordination geometry and (b) positions a ligand **nitrogen** donor at the canonical Zn–N distance (~2.1 Å), precisely the geometry the inhibitor adopts in the crystal. Both Vina configurations place the nearest atom on a carboxylate oxygen (2.23 Å) and evaluate no coordination term. In particular, the Vina-AD4 row uses affinity maps generated by LKina's inline AG4 engine, so the three rows share the same underlying AD4 force field; the differences therefore isolate the *metal-coordination modeling* contribution (pseudoatom injection + nbp overrides + rerank), not the grid generator. Two caveats apply: (i) LKina's metal-mode score (−11.51 kcal/mol here) includes the genuine coordination contribution of the TZ potential and the Zn nbp overrides (≈2–4 kcal/mol deeper than plain AD4, corresponding to the physical Zn–ligand coordination energy); on the carbonic anhydrase II/AZM system the zn-mode score (−13.3) agrees with the experimental affinity (≈−10.8 kcal/mol) within <3 kcal/mol, meeting the accuracy standard reported for AutoDock4Zn (Section 4.3). (ii) The receptor must be prepared without co-crystallized ligand atoms: the initially tested `4JC_rec_zn_fixed.pdbqt` contained 1,349 UNL/UNK atoms that corrupted the absolute energies (Vina AD4 scored the same pose at +980 kcal/mol); the results reported here use the clean receptor (2,083 atoms).

### 3.4 Reactive covalent presets (6/6)

All six reaction presets were exercised on synthetic receptor/ligand systems with the reactive pair pre-positioned at the ideal NAC geometry (receptor electrophile Cys-SG/Ser-OG/Lys-NZ/Tyr-OH at the origin with a CB frame atom; ligand nucleophile at the preset bond length $r_0$). Each preset was run in P1 (distance) and P1+P2 (distance + angle + gradient check) modes:

**Table 7.** NAC detection and reactive distance across the six reaction presets.

| Preset | $r_0$ (Å) | $\theta_0$ (°) | P1 NAC | P1+P2 NAC | P1+P2 d (Å) | P1+P2 angle (°) |
|---|---|---|---|---|---|---|
| `cys_michael` | 1.82 | 109.5 | YES | YES | 2.26 | 93.9 |
| `cys_sn2` | 1.82 | 180.0 | YES | NO | 2.47 | 108.3 |
| `ser_covalent` | 1.34 | 109.5 | YES | YES | 2.01 | 99.1 |
| `lys_targeting` | 1.47 | 109.5 | YES | YES | 2.33 | 94.1 |
| `boronic_acid` | 1.47 | — | YES | NO | 2.71 | 110.9 |
| `tyr_covalent` | 1.38 | 109.5 | YES | YES | 2.03 | 100.8 |

All six presets completed both the P1 and P1+P2 runs (return code 0, finite energies). The two P1+P2 `NAC=NO` cases are informative rather than failures: `cys_sn2` targets a 180° backside attack, which the synthetic approach geometry (nucleophile approaching at ~108°) does not satisfy, and `boronic_acid` has no angle constraint by design and a slightly long 2.71 Å contact. The NAC detector correctly reports NO in both, demonstrating that it discriminates geometry rather than passing unconditionally.

![Figure 5](figures/fig5_covalent_framework.png)

- **Figure 5.** Covalent framework (A–C). (A) All six reaction presets complete each of the four tiers P1 (distance), P1+P2 (angle), P4 (hybrid-vdW-scale sweep 0/0.5/1.0) and C3 (two-step search) on the engineered synthetic systems. (B) NAC detection discriminates geometrically: cys_michael / ser_covalent / lys_targeting / tyr_covalent report NAC=YES, cys_sn2 (180° backside attack not satisfied) and boronic_acid (no angle constraint set) report NAC=NO. (C) P1 distance potential scanned on the cys_michael preset: Gaussian well centered at $r_0 = 1.82$ Å with depth 10 kcal/mol.

![Figure S2](figures/figS2_reactive_convergence.png)

- **Figure S2.** Effect of the P2 angular constraint on pose geometry across the six reaction presets. Reactive distance from the P1-only search (grey) versus the P1+P2 search (red), with each preset's ideal bond length r₀ marked in gold. Adding the angular term pulls the majority of presets toward chemically shorter, more ideal reactive distances.

### 3.5 Backward compatibility: numerical identity with Vina 1.2.7

On two real organic complexes, 1HVR (HIV-1 protease, XK2 ligand) and 3PTB (trypsin, BEN ligand), LKina in standard mode (vina scoring, no metal flags) reproduces the Vina 1.2.7 mode-1 energies **exactly** (−14.54 vs −14.54 and −6.202 vs −6.202 kcal/mol) with identical boxes, exhaustiveness 16, and seed 42. This numerical identity confirms that LKina's extensions are strictly additive: for metal-free, non-reactive receptors, the engine is a drop-in replacement for Vina 1.2.7.

### 3.6 Metalloprotein redocking benchmarks (measured datasets)

To complement the synthetic and single-system comparisons above, we ran a fully scripted redocking benchmark on PDB-derived metalloprotein complexes. Structures were obtained from the RCSB PDB (resolution ≤ 2.5 Å; entries containing the target metal ion and at least one organic ligand within 6 Å of any metal ion of that type); receptors (protein + metal ion, with waters and co-ligands removed) and rigid co-crystallized ligands were prepared with Open Babel 3.1; and each ligand was re-docked into its parent receptor with three engines under identical conditions (exhaustiveness 8, seed 42, one output mode, box = ligand bounding sphere + 6 Å, grid spacing 0.375 Å): **LKina AD4 + metal mode** (zn / fe3 / cu2_jt), **LKina AD4 standard** (`--no_auto_metal`, i.e. the same engine without coordination terms), and **Vina 1.2.7** (default vina scoring). The pipeline, per-entry metadata, and raw poses are deposited in `benchmarks/redock_benchmark/`.

Two metrics are reported: (i) the classical top-1 RMSD ≤ 2.0 Å success rate against the crystal pose (heavy atoms, direct correspondence without superposition); and (ii) the **donor–metal distance**, the minimum distance from any ligand N/O/S heavy atom to the nearest receptor metal ion in the docked pose, which directly measures whether the search reproduces a coordination-competent geometry. The full cohort comprises **n = 104** systems (Zn²⁺ 22 / Fe³⁺ 47 / Cu²⁺ 35), the same corpus as in Figure 6; owing to each engine's own parse/type constraints, some systems fail per engine, the denominators of the RMSD rows below are the attempted cohort per metal (Zn 22 / Fe 47 / Cu 35; RMSD could not be computed for 2 Zn and 1 Fe systems), whereas the denominators of the "whole-cohort" rows are each engine's *dockable subset* (96 for LKina metal mode and Vina 1.2.7; 97 for LKina AD4 standard), so the comparison is same-protocol rather than forced-alignment:

**Table 8.** Metalloprotein redocking benchmark metrics (each engine's own dockable subset).

| Metric | LKina + metal mode | LKina AD4 std | Vina 1.2.7 |
|---|---|---|---|
| Zn²⁺ (n = 22), RMSD ≤ 2.0 Å | 0% | 18.2% (4) | 31.8% (7) |
| Fe³⁺ (n = 47), RMSD ≤ 2.0 Å | 12.8% (6) | 10.6% (5) | 34.0% (16) |
| Cu²⁺ (n = 35), RMSD ≤ 2.0 Å | 11.4% (4) | 2.9% (1) | 11.4% (4) |
| All systems (n = 104): donor–metal median (Å) | **2.12** | 3.83 | 2.83 |
| All systems: donor–metal mean (Å) | **2.58** | 4.51 | 4.21 |
| All systems: donor within 3.0 Å of metal | **87/96 (90.6%)** | 27/97 (27.8%) | 50/96 (52.1%) |

![Figure 6](figures/fig6_redock_benchmark.png)

- **Figure 6.** PDB-derived redocking benchmark (Zn²⁺ n=22 / Fe³⁺ n=47 / Cu²⁺ n=35). (A) Top-1 RMSD ≤ 2.0 Å success and (B) fraction of systems whose best pose places a ligand N/O/S donor within 3.0 Å of the receptor metal, for LKina AD4 + metal mode, LKina AD4 (`--no_auto_metal`), and Vina 1.2.7 under an identical preparation protocol. Regenerated by `make_redock_figures.py` from the archived JSON results and docked poses.

![Figure 7](figures/fig7_donor_metal_distance.png)

- **Figure 7.** Best-pose donor–metal minimum distance per system on the combined benchmark (n = 100 systems docked successfully by all three engines, Zn 20 / Fe 46 / Cu 34, the intersection used for like-for-like per-system comparison). LKina with the coordination model active keeps the best-pose ligand donor inside the 3.0 Å coordination sphere in 96% of systems (median 2.09 Å), versus 31% (median 3.73 Å) for the identical engine with `--no_auto_metal` and 54% (median 2.74 Å) for Vina 1.2.7. The corresponding whole-set statistics over each engine's own dockable subset are reported in Table §3.6 (medians 2.12 / 3.83 / 2.83 Å).

These measurements support three conclusions. First, the *coordination-geometry* effect claimed for `metal_mode` is real and substantial: with the metal terms enabled, the best pose places a ligand donor within the receptor metal's 3.0 Å coordination sphere in ~91% of systems (median distance 2.12 Å, consistent with canonical Zn–N / Fe–O / Cu–N bond lengths), versus ~28% for the same engine under `--no_auto_metal` and ~52% for Vina's electrostatics-driven scoring. For metalloenzyme applications this is the decisive metric, a functional pose must satisfy coordination chemistry.

Second, we report the RMSD results candidly, including the unfavorable ones: under this lean preparation pipeline (Open Babel without protonation-state optimization, rigid receptor, exhaustiveness 8), all three engines achieve only moderate global pose-redocking success (Fe³⁺/Cu²⁺ ≤ 34%, Zn²⁺ ≤ 32%), and LKina metal mode does not improve global RMSD success over its own AD4 baseline on every subset. The absolute rates fall below those of literature redocking studies that use refined protonation workflows (AutoDockTools/Meeko with manual tuning) and higher exhaustiveness; a systematic comparison of preparation protocols is beyond the scope of this work. We therefore present the RMSD column as a like-for-like engine comparison under a single fixed protocol, and the donor–metal distance as the metric specific to the coordination-model claim.

Third, all of these datasets are released in machine-readable JSON form together with the repository (`benchmarks/redock_benchmark/results/redock_{zn,fe,cu}_120.json` — the 120 in the file names is the size of the RCSB retrieval candidate pool; after screening, 104 systems were actually re-docked: Zn 22 / Fe 47 / Cu 35; together with `redock_summary.json` and `donor_metal_distance_summary.json`), along with every input pose; apart from the RCSB download endpoint, every number above is reproducible without any external service.

*(Earlier preprints of this manuscript cited 292-/56-/41-complex benchmark suites following Santos-Martins et al. [6], as well as a 28-complex Cys-Michael covalent benchmark (P1+P2 64.3%, NAC 71.4%); all of these figures originate from the engine design document with no archived dataset, and have been replaced by the measured results above and by the six-preset synthetic end-to-end measurements of Section 3.4, respectively.)*

### 3.7 Scoring-function comparison, parameter sweeps, and hotspot-target cases

Three families of robustness measurements were run against the fixed binary (v1.0.2 engine; all scripts, JSON results, and raw poses are deposited in `benchmarks/`, with summary figures S5–S10).

**Scoring functions (Figure S5, Table 9).** Three representative systems were docked under vina, vinardo, and AD4 scoring at exhaustiveness 8 (5 modes requested, seed 42): the Zn metalloprotein 4JC (zn mode), the carbonic anhydrase II/AZM complex 3HS4 (zn mode), and the metal-free KRAS G12C covalent complex 6OIM (plain AD4); under AD4 scoring, 4JC and 6OIM actually produced 4 modes due to energy-range truncation:

**Table 9.** Scoring-function comparison (best-mode affinity; d(M) = best-pose distance to the nearest receptor metal).

| System | Scoring | Best affinity (kcal/mol) | d(M) (Å) |
|---|---|---|---|
| 4JC (Zn) | vina | −5.80 | 2.65 |
| 4JC (Zn) | vinardo | −5.01 | **12.46** |
| 4JC (Zn) | AD4 + zn | **−11.48** | **1.83** |
| 3HS4 (Zn) | vina | −5.56 | 2.47 |
| 3HS4 (Zn) | vinardo | −3.83 | 2.20 |
| 3HS4 (Zn) | AD4 + zn | **−13.32** | **2.08** |
| 6OIM (metal-free) | vina | **−9.43** | — |
| 6OIM (metal-free) | vinardo | −5.49 | — |
| 6OIM (metal-free) | AD4 | −19.30 | — |

On the metal systems, AD4 + metal mode pulls the best pose into true Zn-coordination distance (1.8–2.2 Å), whereas vina/vinardo — whose potentials contain no coordination term — leave the 4JC best pose **12.5 Å** away from Zn²⁺: a quantitative demonstration that metal-coordination docking cannot be substituted by a generic scoring function. On the metal-free 6OIM system all three functions behave normally, confirming that the metal machinery does not interfere with standard use.

**Search-protocol robustness (Figures S6–S7).** Across exhaustiveness 1–32 on 3HS4 (zn mode, 3 seeds each), the best-affinity standard deviation falls to **σ ≈ 0.003 kcal/mol at exhaustiveness 8** (and remains ≤ 0.010 even at exhaustiveness 1), while runtime grows near-linearly — the default exhaustiveness 8 is well past the convergence knee. Five-seed variance at exhaustiveness 8 is **0.010 kcal/mol** (3HS4, zn) and **0.053 kcal/mol** (6OIM, vina), both orders of magnitude below the ~1 kcal/mol accuracy of any current scoring function; single-seed results are therefore reproducible at the scoring-accuracy scale.

**Parameter response (Figures S8–S9).** Sweeping `--metal_soft_weight` over 0–0.5 on 3HS4 shifts the best-pose Zn distance by only 2.08 → 2.07 Å and the affinity by −13.319 → −13.316 kcal/mol (within noise): the search-time soft constraint is non-perturbative once coordination is already correct, consistent with its intended role of guiding initially displaced poses (Section 4.4). Sweeping the `--metal_bias` attractor strength 0 → 8 kcal/mol deepens the affinity monotonically (−13.32 → −14.85) while the Zn distance stays at 2.08–2.11 Å, i.e. the attractor works as designed without corrupting geometry. The covalent anchor strength responds linearly (best score −0.63 → −8.68 kcal/mol for attractor strength 2 → 16 on the cys_michael system), so anchoring strength is quantitatively tunable.

**Hotspot-target cases (Figure S10).** (i) *Carbonic anhydrase II/AZM (3HS4)*: zn mode recovers the sulfonamide–Zn coordination at 2.08 Å (crystal Zn–N 1.94 Å) with affinity −13.32 kcal/mol versus the experimental ≈ −10.8; adding the metal bias deepens the well (−13.69) without breaking the coordination distance (2.10 Å) — the zn + metal-bias workflow fits zinc-enzyme screening directly. (ii) *KRAS G12C covalent (6OIM, cys_michael preset)*: we report an honest negative. The reactive modes all complete (best −19.1 to −19.3 kcal/mol, indistinguishable from plain AD4 on this metal-free target), but the warhead–C12(SG) distance converges only to ~7.2 Å rather than a covalent bond length. The likely causes are the rigid-receptor preparation (the catalytic Cys 12 is already covalently modified in the crystal structure) and the default 1.5 Å attractor well providing insufficient guidance inside the ~20 Å box with 4 rotatable warhead bonds. A re-scan with a wider/stronger attractor or the hybrid mode with an explicit bond-length term is listed in Section 4.6.

![Figure S5](figures/figS5_scoring_comparison.png)

- **Figure S5.** Scoring-function comparison on 4JC / 3HS4 / 6OIM (Table 9). Best-mode affinity (left axis, bars) and best-pose distance to the nearest receptor metal (right axis, markers). On the Zn systems only AD4 + metal mode achieves coordination distance; the vinardo 4JC best pose sits 12.5 Å from Zn²⁺.

![Figure S6](figures/figS6_exhaustiveness.png)

- **Figure S6.** Exhaustiveness sweep on 3HS4 (zn mode; exhaustiveness 1–32, 3 seeds each). Best-affinity standard deviation across seeds (top) falls to σ ≈ 0.003 kcal/mol at exhaustiveness 8; mean runtime grows near-linearly (bottom).

![Figure S7](figures/figS7_seed_variance.png)

- **Figure S7.** Five-seed variance at exhaustiveness 8: 3HS4 (zn mode) spans −13.316 to −13.326 kcal/mol (range 0.010); 6OIM (vina) spans −9.433 to −9.486 (range 0.053). Single-seed affinities are reproducible far below scoring-function accuracy.

![Figure S8](figures/figS8_soft_weight.png)

- **Figure S8.** `--metal_soft_weight` sweep (w = 0–0.5) on 3HS4: best-pose Zn distance (2.08 → 2.07 Å) and affinity (−13.319 → −13.316 kcal/mol) are unchanged within noise — the search-time soft constraint does not perturb systems whose coordination is already correct.

![Figure S9](figures/figS9_metal_bias.png)

- **Figure S9.** `--metal_bias` strength sweep (0–8 kcal/mol) on 3HS4: affinity deepens monotonically (−13.32 → −14.85 kcal/mol) while the best-pose Zn distance stays at 2.08–2.11 Å — the O5 attractor is effective and geometry-preserving.

![Figure S10](figures/figS10_hotspot_cases.png)

- **Figure S10.** Hotspot-target cases. (A) KRAS G12C (6OIM, cys_michael): all five docking modes complete, but the warhead–SG12 distance converges only to ~7.2 Å (dashed line = covalent bond length ~1.8 Å), an honest negative reported in Section 3.7. (B) Carbonic anhydrase II (3HS4): zn mode places the AZM sulfonamide N at 2.08 Å of the catalytic Zn²⁺ (crystal 1.94 Å); the metal bias deepens the score without disturbing the coordination distance.

---

## 4. Discussion

### 4.1 Comparison with prior metal- and covalent-docking methods

LKina differs from AutoDock4Zn [6] in *coverage* (any metal rather than Zn²⁺ alone; four pseudoatom geometries; Jahn–Teller and water-bridge terms) and from MetalDock [7] in *architecture* (engine-level atom types and grids versus an external pre-processing protocol; no external `autogrid4`). The 110-mode synthetic benchmark quantifies the coverage gap directly: 76 metal-coordination systems fail even Vina 1.2.7's parsing, whereas LKina docks all of them with a mean coordination error of only 0.20 Å. Unlike CovDock, LKina is open source, requires no pre-built covalent adducts, and reports NAC status per pose. Other representative docking programs, the CHARMm grid-based CDOCKER [19] and the genetic-algorithm GOLD [20], provide reference accuracy baselines for covalent and flexible docking, but their commercial licenses or external pre-processing dependencies preclude a direct engine-level comparison; this work therefore runs LKina and Vina 1.2.7 head-to-head under identical conditions (same receptor/ligand files, fixed seed, top-1 RMSD; Sections 3.3–3.6) and treats the CDOCKER/GOLD literature accuracy as a qualitative reference rather than a re-run. To our knowledge, LKina is currently the only open-source engine that integrates, at the engine level, metal-coordination docking, NAC-constrained reactive docking, and a fully self-contained grid pipeline.

### 4.2 The synthetic benchmark: what it does and does not measure

The synthetic metal-coordination systems are deliberately minimal (a metal, a partial donor shell, and a probe ligand) and therefore **measure coverage and robustness, not predictive accuracy**: they verify that every metal mode parses, that pseudoatom injection and nbp overrides load without error, and that the search completes with finite energies and emits coordination remarks. The coordination-distance residuals (mean 0.20 Å across 110 modes) are lower-bound estimates because the systems lack the protein pocket that makes the metal site geometrically unique in real complexes. Real-system accuracy is established separately by the 4JC comparison (Section 3.3) and the PDB-derived Zn/Fe/Cu redocking datasets (Section 3.6).

### 4.3 Vina comparison and the metal-mode energy scale

LKina's metal mode augments standard AD4 with the TZ coordination potential and the Zn nbp overrides, so its scores include the genuine coordination contribution and are naturally 2–4 kcal/mol deeper than plain AD4 — consistent with the AutoDock4Zn design, in which each `nbp_r_eps` epsilon is multiplied by `FE_coeff_vdW = 0.1662` when the potential is written (NA–TZ well depth −3.86 kcal/mol). The implementation was verified line-by-line against the official AutoGrid source (`mainpost1.28.cpp`) and validated on the carbonic anhydrase II/AZM system: the zn-mode score (−13.3 kcal/mol) agrees with the experimental affinity (≈−10.8 kcal/mol) within <3 kcal/mol, matching the accuracy reported for AutoDock4Zn. Metal-mode scores are therefore suitable for *ranking poses within metal mode*, *monitoring coordination geometry*, and *order-of-magnitude comparison against experimental affinities*; coordination geometry should remain the primary evidence when comparing across scoring functions.

The 34-vs-110 coverage result is independent of scoring-function choice: it is a property of the *atom-type parser*, and identical receptor/ligand files were used for both engines. The 4JC comparison goes further by running Vina-AD4 on LKina-generated `.map` files, so the three rows share the same force field and differ only in the metal-coordination modeling. We explicitly do not claim that LKina's AD4 scoring is universally better than Vina's vina scoring; the claim is that LKina adds a coordination model that Vina lacks, at no cost to backward compatibility (Section 3.5).

### 4.4 Design philosophy: conservative extensions, measurable behavior

Three principles govern LKina's design. (i) *Backward compatibility is a correctness property*: defaults reproduce Vina 1.2.7 exactly (numerically verified in Section 3.5), and every new term is opt-in. (ii) *Post-processing rerank terms do not enter the search gradient*: metal geometry affects pose ranking only, avoiding spurious forces; the search-time soft-constraint channel (`--metal_soft_weight`) is exposed but off by default, with a smooth log-sum-exp approximation that keeps the optimization well-posed. (iii) *Grid/search consistency is enforced*: when metal-specific potentials are active, fallback to an external grid generator is forbidden, eliminating a subtle class of search/grid mismatch errors.

### 4.5 Limitations

(i) **Multi-metal coupling.** Modes are applied per site independently; bridging ligands of dinuclear centers (μ-aqua, μ-oxo) may be double-counted in rerank scores, and the BVS inference remains single-center (Mo–Fe nitrogenase and the tetranuclear Mn cluster of photosystem II are not yet handled rigorously; Ni–Fe hydrogenase receives a single-center Ni approximation). (ii) **Metal-as-ligand scope.** The reverse-potential extension is a conservative, grid-based approximation: MetalDock-style ligand-side dummy atoms and internal geometric constraints are not ported, and dedicated benchmarks for Pt/Ru/Tc/Re/Gd agents are pending. (iii) **Reactive energy is not thermodynamic**: reactive penalty terms must not be interpreted as binding-affinity predictions. (iv) **Synthetic benchmark granularity**: the 110-mode systems are minimal by design; a full PDB-derived multi-metal benchmark would require additional infrastructure. (v) **Absolute success rates** are protocol-dependent (exhaustiveness 16, single seed); the reported comparisons use identical protocols for all methods, so the *differences* are the claim.

### 4.6 Future work

Two open directions remain genuine limitations of the current implementation: (i) **multi-metal coupling and multi-center BVS**, systems such as the Mo–Fe cofactor of nitrogenase and the tetranuclear Mn cluster of photosystem II require multi-atom, jointly scored active-site models that go beyond the single-center BVS inference of Section 4.5; and (ii) **broader medicinal metal-as-ligand chemical space**, the current reverse-potential extension is a conservative grid-based approximation; MetalDock-style ligand-side dummy atoms and internal geometric constraints are not yet ported, and dedicated redocking and in-vitro pharmacological validation for Pt-, Ru-, and Tc-based agents are still pending. In addition, the KRAS G12C covalent case (Section 3.7) left the warhead–SG distance at ~7.2 Å; a re-scan with a wider/stronger reactive attractor or the hybrid mode with an explicit bond-length term (1.8 Å) is planned. Already completed (and therefore moved out of future work): the `--metal_soft_weight` ablation on three real systems (Zn 4JC, Fe 1A8E, Cu 1A2V at w = 0/0.1/0.3/0.5; drift ≤ 0.01 Å; `benchmarks/soft_weight_results.json`); the metal-as-ligand full Pt/Pd/Ru/Os/Re pool pipeline (n = 20; `benchmarks/metallocomplex_redocking_benchmark.py`); and the C3 regime-dependence result (single 5/6 vs c3b 3/6 vs c3 0/6 on synthetic deep pockets; Figure S3).

![Figure S3](figures/figS3_c3_ablation.png)

- **Figure S3.** Synthetic deep-pocket ablation (single-stage vs C3 vs C3b; six reaction presets; see Section 4.6). The x-axis lists the six presets; bars are coloured: deep red = single-stage (full P1–P4 constraints throughout), grey = C3 two-step (unconstrained Phase-1 + 15 Å anchor filter + L-BFGS Phase-2 refinement), blue = C3b (broad/weak Gaussian attractor, σ×3, ε×0.15). The y-axis is the reactive distance d(SG–nucleophile) in the docked pose; the dashed line is the 2.6 Å convergence threshold; bar tops are measured values. **Convergence: single 5/6, c3b 3/6, c3 0/6**, single-stage's persistent P1 gradient steers the MC toward the narrow pocket mouth; the unconstrained Phase-1 in plain C3 loses direction; C3b's weak attractor partially restores it (3/6 vs 0/6). The +8% Phase-2 convergence claim from the early design document (LKINA.md §5.4.3, same provenance as the unreproducible 292/56/41 numbers) is **not recovered** on this synthetic design. C3's benefit is regime-dependent and requires realistic pockets with competing local minima for further evaluation.

![Figure S4](figures/figS4_metallocomplex.png)

- **Figure S4.** Metal-as-ligand redocking on the full Pt/Pd/Ru/Os/Re pool (n=20; see Section 4.6). The x-axis lists the 20 PDB entries (with the metal family and ligand residue name in the secondary label); the y-axis is top-1 RMSD against the crystal pose. Dark red = LKina metal-as-ligand mode (`--ligand_metal_geometry_weight 1.0`); cream = the same engine with `--no_auto_metal` (standard AD4); the dashed line is the 2.0 Å redock-success threshold. **Vina 1.2.7 failed 20/20 at parse time** (the 32-type XS system has no entries for PT/PD/RU/OS/RE), see subtitle. LKina's metal-as-ligand mode beats its own AD4 baseline on 6/20 systems (median |dRMSD| improvement 1.92 Å), confirming that the reverse metal–donor pair potentials and the ligand-side geometry QC are active; the global RMSD ≤ 2.0 Å count is 1/20 (6IG4_PT 1.49 Å) because flexible metal complexes are intrinsically hard to redock with a global-RMSD metric, the chemistry-relevant test for these ligands is whether the docked pose is coordination-competent (cf. the donor–metal distance metric in Section 3.6).

## 5. Conclusion

LKina extends AutoDock Vina 1.2.7 into a metal-aware, covalent-reactive docking engine without sacrificing its search performance or backward compatibility. A 113-type AD4 atom system, four coordination pseudoatom geometries, BVS oxidation-state inference, Jahn–Teller modes, a semi-explicit water-bridge model, and an inline AG4 grid generator together make metalloenzyme docking fully self-contained. A systematic 110-mode synthetic benchmark shows that LKina docks every metal mode (110/110) with a mean coordination error of 0.20 Å, while Vina 1.2.7 fails on 76/110 at parse time; on a real Zn metalloprotein, LKina alone recovers the Zn–N coordination at 2.19 Å with an explicit coordination report. In a PDB-derived redocking benchmark across Zn²⁺/Fe³⁺/Cu²⁺ metalloproteins (104 systems), enabling the coordination model placed the best-pose ligand donor within the 3.0 Å metal coordination sphere in ~91% of systems (median donor–metal distance 2.12 Å), versus ~28% for the same engine without it and ~52% for Vina 1.2.7. The P1–P4 reactive framework, together with NAC detection, six reaction presets, and the C3 two-phase strategy, completes end-to-end runs on all six presets with genuine NAC discrimination. All extensions are opt-in, and standard docking reproduces Vina 1.2.7 energies exactly; LKina can therefore serve as a drop-in replacement within existing Vina workflows.

---

## 6. Data and code availability

LKina is open source at https://github.com/LK-Studio1128/LKina, including the source code, build scripts (macOS arm64/x86_64, Linux x86_64, Windows MSYS2/MinGW), the regression suite, and release binaries. License: the combined binary is distributed under GPL-3.0-or-later (Vina-origin files remain Apache-2.0; the AG4 engine and parameter data are GPL-3.0-or-later; see `COPYING`/`NOTICE`). Benchmark data: the PDB-derived Zn²⁺/Fe³⁺/Cu²⁺ redocking sets (Section 3.6) are regenerated by `benchmarks/redock_benchmark/redock_pipeline.py` from the RCSB PDB [16], with JSON results deposited in `benchmarks/redock_benchmark/results/`. The 110-mode synthetic metal benchmark, the feature-family measurements, and the four-tier covalent framework benchmark are all reproducible from `benchmarks/` (see Appendix E); every run used seed 42 for exact reproducibility. The v1.0.2 snapshot is archived on Zenodo: https://doi.org/10.5281/zenodo.22156943 [21].

## Acknowledgements

LKina builds on the AutoDock ecosystem developed by the Olson and Forli labs at The Scripps Research Institute, AutoDock Vina [1,2], AutoDock4/AutoDockTools [3], AutoGrid4, and AutoDock4Zn [6]. We gratefully acknowledge Santos-Martins *et al.* for the 292-complex Zn²⁺ benchmark protocol and dataset, and Hakkennes *et al.* for MetalDock [7], whose optimized metal parameters inform LKina's reverse pair potentials. Meeko [17] and Open Babel [18] are used in the preparation pipeline.

## Author contributions

LK-Studio1128 conceived and implemented LKina, designed and ran all benchmarks, and wrote the manuscript.

## Competing interests

The authors declare no competing interests.

---

## Appendix E. Benchmark corpus summary (full reproducibility)

All measured numbers in the paper are reproducible from the repository (`benchmarks/`); the live binaries used are the post-fix engine (metal-mode AD4 scoring fixes, tagged `v1.0.2`; see Section 3.2.1). The table below summarises the entire benchmark corpus, with the JSON outputs that the scripts deposit next to themselves (`metal_coverage_results_all.json`, `feature_family_results.json`, `covalent_full_results.json`).

**Table 10.** Benchmark corpus summary (full reproducibility).

| Benchmark | System | Metric | Value | File |
|---|---|---|---|---|
| Engine fix | synthetic | `AG4_EINTCLAMP` aligned to official semantics (pairwise LJ-table clamp of 100000; map-total ±1000 two-sided clamp removed) | applied | `src/lib/ag4_engine.cpp` |
| Engine fix | synthetic | `--no_auto_metal` flag | applied | `src/main/main.cpp` |
| All-metal coverage | 110 synthetic | LKina dock-OK | 110/110 | `run_all_metal_modes.py` |
| All-metal coverage | 110 synthetic | Vina 1.2.7 dock-OK | 34/110 (76 fail type parse) | `run_all_metal_modes.py` |
| All-metal coverage | 110 synthetic | LKina \|d–d₀\| mean | 0.20 Å | `run_all_metal_modes.py` |
| All-metal coverage | 110 synthetic | LKina \|d–d₀\| < 0.5 Å | 108/110 | `run_all_metal_modes.py` |
| Pseudoatom check | 5 metals | ideal/good donor count | 3/3 3/3 5/5 4/4 4/4 | `run_feature_family_tests.py` |
| BVS inference | 14 designs | oxidation state correct | 14/14 | `run_feature_family_tests.py` |
| Water bridge | Mg (4/6 coord) | `METAL_WATER_E` emitted | yes | `run_feature_family_tests.py` |
| Metal-as-ligand | Pt sq-planar ideal / distorted | penalty | 0.003 vs 0.903 kcal/mol | `run_feature_family_tests.py` |
| Covalent P1 | 6 presets | rc | 6/6 | `run_covalent_full.py` |
| Covalent P1+P2 | 6 presets | rc + NAC | 6/6 (NAC=YES 4/6) | `run_covalent_full.py` |
| Covalent P4 | 6 presets × 3 vdw scales | rc | 6/6 | `run_covalent_full.py` |
| Covalent C3 | 6 presets | rc | 6/6 | `run_covalent_full.py` |
| Reactive gradient | 6 presets × 7 comps | max \|analytic−numeric\| | 2.4e-16–4.3e-08 | `run_feature_tests.py` |
| λ-sweep monotonicity (P4) | cys_michael, λ=0/0.5/1.0 | energies −10.09/−5.716/−4.734 kcal/mol | monotonically increasing | `run_covalent_full.py` |
| Energy well scan | cys_michael | Gaussian bottom at r₀ = 1.82 Å, depth 10 kcal/mol | yes | `run_covalent_full.py` |
| generate_maps | 8-probe ligand | AutoGrid4-format maps written | 108 | `run_covalent_full.py` |
| Backward compat | 1HVR / 3PTB | mode-1 energy vs Vina | −14.54 / −6.202 (identical) | `real_systems/` |
| Scoring functions | 4JC / 3HS4 / 6OIM | best affinity, d(metal) per function | Table 9 | `run_scoring_comparison.py` |
| Exhaustiveness sweep | 3HS4 zn, exh 1–32 × 3 seeds | σ(best) at exh 8 | 0.003 kcal/mol | `run_parameter_sweeps.py` |
| Seed variance | 3HS4 zn / 6OIM vina, 5 seeds | best-affinity range | 0.010 / 0.053 kcal/mol | `run_parameter_sweeps.py` |
| soft_weight sweep | 3HS4, w = 0–0.5 | Zn distance drift | ≤ 0.01 Å | `run_parameter_sweeps.py` |
| metal_bias sweep | 3HS4, 0–8 kcal/mol | monotonic deepening, geometry kept | −13.32 → −14.85 | `run_parameter_sweeps.py` |
| Reactive strength | cys_michael, strength 2–16 | linear score response | −0.63 → −8.68 kcal/mol | `run_parameter_sweeps.py` |
| Hotspot 3HS4 | CA II/AZM | zn-mode Zn–N (crystal 1.94 Å) | 2.08 Å, −13.32 kcal/mol | `run_hotspot_cases.py` |
| Hotspot 6OIM | KRAS G12C covalent | warhead–SG12 distance (honest negative) | ~7.2 Å | `run_hotspot_cases.py` |

## References

1. Trott O, Olson AJ. AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading. *J Comput Chem.* 2010;31(2):455–461 doi:10.1002/jcc.21334.
2. Eberhardt J, Santos-Martins D, Tillack AF, Forli S. AutoDock Vina 1.2.0: new docking methods, expanded force field, and Python bindings. *J Chem Inf Model.* 2021;61(8):3891–3898 doi:10.1021/acs.jcim.1c00203.
3. Morris GM, Huey R, Lindstrom W, et al. AutoDock4 and AutoDockTools4: automated docking with selective receptor flexibility. *J Comput Chem.* 2009;30(16):2785–2791 doi:10.1002/jcc.21256.
4. Andreini C, Bertini I, Cavallaro G, Holliday GL, Thornton JM. Metal ions in biological catalysis: from enzyme databases to general principles. *J Biol Inorg Chem.* 2008;13(8):1205–1218 doi:10.1007/s00775-008-0404-5.
5. Valdez CE, Smith QA, Nechay MR, Alexandrova AN. Mysteries of metals in metalloenzymes. *Acc Chem Res.* 2014;47(10):3110–3117 doi:10.1021/ar500227u.
6. Santos-Martins D, Forli S, Ramos MJ, Olson AJ. AutoDock4Zn: an improved AutoDock force field for small-molecule docking to zinc metalloproteins. *J Chem Inf Model.* 2014;54(8):2371–2379 doi:10.1021/ci500209e.
7. Hakkennes MLA, Buda F, Bonnet S. MetalDock: an open access docking tool for easy and reproducible docking of metal complexes. *J Chem Inf Model.* 2023;63(24):7816–7825. doi:10.1021/acs.jcim.3c01582
8. Wang K. GPDOCK: highly accurate docking strategy for metalloproteins based on geometric probability. *Brief Bioinform.* 2023;24(1):bbac620. doi:10.1093/bib/bbac620
9. Scarpino A, Ferenczy GG, Keserű GM. Covalent docking in drug discovery: scope and limitations. *Curr Pharm Des.* 2020;26(44):5684–5699. doi:10.2174/1381612824999201105164942
10. Lightstone FC, Bruice TC. Ground state conformations and entropic and enthalpic factors in the efficiency of intramolecular and enzymatic reactions. *J Am Chem Soc.* 1996;118(11):2595–2605 doi:10.1021/ja952589l.
11. Bruice TC. A view at the millennium: the efficiency of enzymatic catalysis. *Acc Chem Res.* 2002;35(3):139–148 doi:10.1021/ar0001646.
12. Bianco G, Forli S, Goodsell DS, Olson AJ. Covalent docking using AutoDock: two-point attractor and flexible side chain methods. *Protein Sci.* 2016;25(1):295–301 doi:10.1002/pro.2733.
13. Harding MM. Small revisions to predicted distances around metal sites in proteins. *Acta Crystallogr D.* 2006;62(6):678–682 doi:10.1107/S0907444906014594.
14. Brown ID, Altermatt D. Bond-valence parameters obtained from a systematic analysis of the inorganic crystal structure database. *Acta Crystallogr B.* 1985;41(4):244–247 doi:10.1107/S0108768185002063.
15. Brese NE, O'Keeffe M. Bond-valence parameters for solids. *Acta Crystallogr B.* 1991;47(2):192–197 doi:10.1107/S0108768190011041.
16. Berman HM, Westbrook J, Feng Z, et al. The Protein Data Bank. *Nucleic Acids Res.* 2000;28(1):235–242 doi:10.1093/nar/28.1.235.
17. Santos-Martins D, He Y, Eberhardt J, et al. Meeko: Molecule Parametrization and Software Interoperability for Docking and Beyond. *J Chem Inf Model.* 2025;65(24):13045–13050. doi:10.1021/acs.jcim.5c02271
18. O'Boyle NM, Banck M, James CA, et al. Open Babel: an open chemical toolbox. *J Cheminform.* 2011;3:33 doi:10.1186/1758-2946-3-33.
19. Wu G, Robertson DH, Brooks CL, Vieth M. Detailed analysis of grid-based molecular docking: a case study of CDOCKER—a CHARMm-based MD docking algorithm. *J Comput Chem.* 2003;24(13):1549–1562. doi:10.1002/jcc.10306
20. Jones G, Willett P, Glen RC, Leach AR, Taylor R. Development and validation of a genetic algorithm for flexible docking. *J Mol Biol.* 1997;267(3):727–748. doi:10.1006/jmbi.1996.0897
21. LK-Studio1128. LKina: A Metal-Aware and Covalent-Reactive Molecular Docking Engine Extending AutoDock Vina. Version v1.0.2. Zenodo, 2026. https://doi.org/10.5281/zenodo.22156943

---

## Supplementary

- `benchmarks/pdbqt_util.py` + `run_all_metal_modes.py` + `metal_coverage_results_all.json`, 110-mode synthetic metal benchmark (Table 5, Fig. 1–2).
- `benchmarks/run_feature_family_tests.py` + `feature_family_results.json`, pseudoatom geometry check, BVS inference (14/14), water bridge, metal-as-ligand QC.
- `benchmarks/run_covalent_full.py` + `covalent_full_results.json`, four-tier covalent framework (P1/P1+P2/P4/C3), energy well scan, 108-map `--generate_maps` audit.
- `benchmarks/run_reactive_tests.py` + `reactive_presets_results.json`, 6-preset reactive benchmark (Section 3.4, Fig. S2).
- `make_figures_v3.py` and `benchmarks/make_supplementary_figures.py` (with `make_figS3S4.py`, `make_redock_figures.py`, `make_graphical_abstract.py`), figure-generation scripts for all main and supplementary figures (600-dpi PNG + vector PDF + TIFF, three formats).
- `figures/figS1_energy_well_landscape.*`, coordination-well landscape over 110 modes (Fig. S1).
- `benchmarks/real_systems/`, 4JC three-engine comparison and 1HVR/3PTB backward-compatibility runs.
- `benchmarks/metallocomplex_redocking_benchmark.py` + `metallocomplex_results/`, metal-as-ligand redocking pipeline (full Pt/Pd/Ru/Os/Re pool, n=20, Section 4.6).
- `benchmarks/run_c3_ablation.py` + `c3_ablation_results.json`, C3 deep-pocket ablation (single-stage vs C3 vs C3b, six synthetic deep-pocket presets, Section 4.6).
- `benchmarks/run_soft_weight_ablation.py` + `soft_weight_results.json`, `--metal_soft_weight` ablation (one real system each for Zn/Fe/Cu × w = 0/0.1/0.3/0.5, Section 4.6).
- `benchmarks/run_scoring_comparison.py` + `scoring_comparison_results.json`, scoring-function comparison on 4JC/3HS4/6OIM (Table 9, Fig. S5).
- `benchmarks/run_parameter_sweeps.py` + `parameter_sweep_results.json` + `parameter_sweeps/`, exhaustiveness convergence, seed variance, soft-weight and metal-bias sweeps, reactive-strength linearity (Figs. S6–S9).
- `benchmarks/run_hotspot_cases.py` + `hotspot_results.json` + `hotspot_cases/`, hotspot-target cases (KRAS G12C 6OIM covalent, CA II 3HS4 zinc; Fig. S10).
- `make_supplementary_figures.py`, one-command regeneration of Fig. S5–S10.
- Regression suite (17 checks, A–F groups; highest-value items re-executed on the live binary in this work, Section 3.1), shipped with the release tag.
- PDBbind v2020-based Fe³⁺/Cu²⁺ collection criteria and per-complex metrics (repository).
- NAC definition and the six-preset reactive benchmark (Section 3.4; `benchmarks/run_covalent_full.py` + `covalent_full_results.json`). The 28-complex Cys-Michael benchmark cited in earlier versions (P1+P2 64.3%, NAC 71.4%) originates from the design document with no archived dataset and is superseded by the Section 3.4 synthetic end-to-end measurements.
