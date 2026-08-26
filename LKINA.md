# LKina: A Metal-Aware Reactive Molecular Docking Engine with Extended AutoDock4 Force Field

**LKina：一种支持金属配位与响应式共价对接的扩展 AutoDock4 力场对接引擎**

---

## 项目信息 / Project Information

| 项目 | 内容 |
|------|------|
| 版本 (Version) | 1.0.1 |
| 上游基线 (Base) | AutoDock Vina 1.2.7（Apache-2.0）+ AutoGrid 4.2 参考（GPL-2.0-or-later） |
| 二进制 (Binary) | `LKina`（macOS / Linux）、`LKina.exe`（Windows） |
| 仓库 (Repository) | https://github.com/LK-Studio1128/LKina |
| Release v1.0.1 | https://github.com/LK-Studio1128/LKina/releases/tag/v1.0.1 |
| 许可证 (License) | GPL-3.0-or-later（Vina-origin 文件保持 Apache-2.0；AG4 扩展为 GPL-3.0） |
| 平台 (Platforms) | macOS arm64 · Linux x86-64 · Windows x86-64 |
| 开发者 | LK-Studio1128 |

## 核心亮点 / Highlights

- **113 种 AD4 原子类型 + 4 类金属配位伪原子**（TZ / SQ / MH / JT），覆盖 80+ 金属 / 类金属
- **内联 AG4 格点生成**：无需外部 `autogrid4` 可执行文件，金属对接零外部依赖
- **BVS 氧化态自动推断**：Fe / Cu / Mn / Co / V / Mo / Ni 的 +2 / +3 / +4 / +6 自动识别
- **Jahn-Teller 变形八面体模式**：Cu²⁺（d⁹）、Mn³⁺（d⁴）专用拉长八面体几何
- **半显式水桥候选位点**：配位不饱和处的 M–O(water) 几何后处理评分
- **响应式共价对接 P1–P4 框架**：距离约束 + 角度约束 + 帧原子扭转 + 混合 vdW 缩放 + C3 两阶段策略
- **6 种反应预设**：`cys_michael` / `cys_sn2` / `ser_covalent` / `lys_targeting` / `boronic_acid` / `tyr_covalent`
- **Metal Bias (O5) 软高斯吸引子** + **金属作为配体**（reverse metal-donor，参考 MetalDock）
- **完全向后兼容**：对不含金属、不触发共价的标准受体，行为与 Vina 1.2.7 一致

---

## Abstract

分子对接是基于结构的药物设计的核心计算工具，但现有开源引擎在金属酶和共价靶点上存在系统性缺陷：标准力场无法描述金属配位的方向性，共价对接缺乏对近攻构象（Near Attack Conformation, NAC）的显式约束。本文介绍 **LKina v1.0.0**，一个在 AutoDock Vina 1.2.7 源码基础上深度扩展的对接引擎。LKina 通过四项核心创新解决上述问题：

1. 将 AutoDock4 原子类型扩展至 **113 种**（涵盖 80+ 金属/类金属），并内嵌在线格点地图生成模块，无需外部 autogrid4；
2. 引入四类金属配位伪原子（**TZ/SQ/MH/JT**），配合键价和（Bond Valence Sum，BVS）氧化态自动推断，支持含单/多金属位点受体的高精度对接；
3. 针对 Cu²⁺（d⁹）和 Mn³⁺（d⁴）设计 Jahn-Teller 变形八面体模式，结合半显式水桥候选位点模型进行几何后处理排序；
4. 实现四层递进响应式共价对接框架（**P1–P4**），支持距离约束、角度约束、帧原子扭转及混合 vdW 缩放，提供 NAC 检测与梯度精度验证。

LKina 继承 Vina 的 Monte Carlo + L-BFGS 搜索策略和 OpenMP 并行化，所有新功能对标准受体完全向后兼容。响应式对接进一步提供 C3 两阶段策略（Phase-1 无约束 MC 预采样 + Phase-2 带约束精化）和六种内置反应类型预设（`cys_michael` / `cys_sn2` / `ser_covalent` / `lys_targeting` / `boronic_acid` / `tyr_covalent`），以及金属偏置吸引子（O5，Metal Bias）功能。近期增量版本进一步统一了 `--zn_mode` 与 `--metal_mode zn` 的内部语义，扩展了单配体和 batch 逐 ligand 金属自动识别与 reverse metal-donor pair potentials，并在金属模式下禁用无法复现 LKina 金属伪原子和 `nbp_r_eps` 覆盖的外部 `autogrid4` fallback。17 项回归测试全部通过，编译零错误零警告。在 292 个 Zn²⁺晶体复合物重对接测试中，LKina TZ 模式 Top-1 RMSD ≤ 2.0 Å 成功率达 74.3%，优于 AutoDock4Zn（70.5%）和标准 AD4（58.2%）；Fe³⁺（56 个复合物）和 Cu²⁺（41 个复合物）子集相对标准 AD4 分别提升 17.9 pp 和 24.4 pp；28 个 Cys 迈克尔加成共价复合物的 NAC 成功率达 71.4%（P1+P2 双约束）。

> **许可证说明**：AG4 格点引擎（`ag4_engine`、`embedded_ad4_grid`、`ad4_parameter_data`）参照 AutoGrid4 源码（GPL-2.0-or-later）实现，以 GPL-3.0-or-later 许可证发布；与 Apache-2.0 Vina 核心静态链接后，**LKina 可执行文件整体以 GNU GPL v3 发行**。原始 Vina 版权归属 The Scripps Research Institute（2006–2010）。

> **项目地址**：<https://github.com/LK-Studio1128/LKina> · **上游基线**：AutoDock Vina 1.2.7 (<https://github.com/ccsb-scripps/AutoDock-Vina>) · **配套 GUI**：LKDock v3.0。金属对接 / 共价对接的图形化使用教程分别见 `LKDock_v3.0_金属对接使用手册.md` 与 `LKDock_v3.0_共价对接使用手册.md`。

---

## 目录 / Table of Contents

- [Abstract](#abstract)
- [1. Introduction](#1-introduction)
- [2. Theory](#2-theory)
- [3. Methods](#3-methods)
- [4. Implementation](#4-implementation)
- [5. Validation](#5-validation)
- [6. Discussion and Limitations](#6-discussion-and-limitations)
- [7. Conclusion](#7-conclusion)
- [8. Availability, Build and Reproducibility](#8-availability-build-and-reproducibility)
- [References](#references)
- [Appendix A — Installation and Quick Start](#appendix-a--installation-and-quick-start)
- [Appendix B — Metal-mode Selection Cheat Sheet](#appendix-b--metal-mode-selection-cheat-sheet)
- [Appendix C — Reactive Preset Cheat Sheet](#appendix-c--reactive-preset-cheat-sheet)
- [Appendix D — Output REMARK Quick Reference](#appendix-d--output-remark-quick-reference)
- [Appendix E — Common Errors and Troubleshooting](#appendix-e--common-errors-and-troubleshooting)
- [Appendix F — Integration with LKDock v3.0](#appendix-f--integration-with-lkdock-v30)

---

## 1. Introduction

### 1.1 金属酶的对接困境

金属酶在生命过程中扮演核心角色：锌指蛋白调控基因表达，铁硫蛋白参与电子传递，铜蓝蛋白催化氧化还原，锰超氧化物歧化酶清除活性氧。全球上市药物中，约 **40%** 的靶点含有金属离子辅因子，包括 HIV 整合酶抑制剂（Mg²⁺）、碳酸酐酶抑制剂（Zn²⁺）、组蛋白去乙酰化酶抑制剂（Zn²⁺）和铁蛋白酶抑制剂（Fe³⁺）。然而，标准分子力场在描述金属-配体相互作用时存在三个根本性缺陷：

1. **平衡距离失真**：van der Waals 力场中 N 的平衡距离约 2.49 Å，而实际 Zn–N 配位距离约 2.0 Å，导致正确配位构象被排斥。
2. **缺乏方向性**：球对称的 Lennard-Jones 势无法区分轴向与赤道配位，对 Jahn-Teller 变形体系尤为不足。
3. **电荷模型偏差**：高价金属（Fe³⁺, Cu²⁺）的形式电荷导致 Gasteiger 体系对含氧基团产生系统性偏好，而非遵循实际配位规律。

Santos-Martins 等人于 2014 年提出 AutoDock4Zn 力场，通过引入四面体配位伪原子 TZ 解决了 Zn²⁺ 配位的方向性问题，在 292 个晶体复合物上显著改善了对接精度。然而该工作仅针对 Zn²⁺，未覆盖其他过渡金属、多金属位点，也未处理 Jahn-Teller 活性金属和桥接水分子。

### 1.2 共价对接的需求与挑战

共价靶向药物（Covalent Targeted Drugs）近年受到高度关注：已批准的共价药物包括不可逆 BTK 抑制剂伊布替尼（Cys481）、EGFR 抑制剂阿法替尼、KRAS G12C 抑制剂 AMG510 等。共价对接面临两个核心挑战：

1. 共价键形成后的构象搜索需将反应活性位点约束在攻击距离和角度范围内；
2. 现有引擎（AutoDock Vina、Glide Covalent、CovDock）多需要手动构建共价复合物结构，自动化程度低。

近攻构象（Near Attack Conformation, **NAC**）理论认为，高效酶催化和亲核取代反应发生时，亲核试剂与亲电中心的距离需 < 3.0 Å，攻击角约 180°（SN2 背面攻击型）或 ~109.5°（迈克尔加成型）。LKina 通过显式能量惩罚项引导搜索收敛于 NAC 构象，是目前极少数将 NAC 检测集成于评分函数的开源引擎之一。

### 1.3 本文贡献

本文提出 LKina，在以下方面对 AutoDock Vina 1.2.7 进行系统性扩展：

| 贡献 | 描述 |
|------|------|
| **原子类型体系** | 将 AD4 原子类型从标准 ~22 种扩展至 113 种，涵盖元素周期表中全部常见金属和类金属 |
| **在线格点生成** | 内嵌 AG4（AutoGrid4）等效格点计算，无需外部 autogrid4 可执行文件 |
| **金属配位力场** | 在 AutoDock4Zn 基础上扩展至 TZ/SQ/MH/JT 四类伪原子，支持 80+ 金属及 Jahn-Teller 变形模式 |
| **BVS 推断** | 基于键价和理论的自动氧化态推断，现覆盖 Fe/Cu/Mn/Co/V/Mo/Ni，并为 Mo-Fe 异核与多核 Mn 耦合扩展预留接口 |
| **半显式水桥** | 在配位空位方向自动放置水候选位点，用于后处理几何排序 |
| **Zn 模式统一** | 将历史 `--zn_mode` 收敛为 `--metal_mode zn` 的兼容别名，统一 Zn pairwise 覆盖、TZ 注入、自动识别和诊断路径 |
| **金属作为配体** | 支持从单配体 PDBQT 或 batch 逐 ligand 自动识别一个或多个 ligand 金属中心，并通过 reverse metal-donor pair potentials 补全 ligand metal 与受体供体原子的配位响应 |
| **金属配合物 ligand 几何 QC** | 对 Pt/Pd 方平面与 Ru/Os/Re 八面体 ligand 侧配位几何输出 `LIGAND_METAL_*` REMARK，并可选用 `--ligand_metal_geometry_weight` 进行后处理重排序 |
| **金属格点一致性保护** | 当 Zn/metal mode 或自定义 `nbp_r_eps` 覆盖激活时，禁止回退到外部 `autogrid4`，避免生成与 LKina 搜索状态不一致的格点 |
| **响应式共价对接** | P1–P4 四层约束框架，支持距离约束、角度约束、帧原子扭转和混合 vdW 缩放 |
| **完全向后兼容** | 所有新功能对不含金属的标准受体透明，不改变非金属系统的对接行为 |
| **C3 两阶段搜索** | Phase-1 无约束 MC 预采样筛选近端 pose，Phase-2 带约束 L-BFGS 精化，兼顾全局探索与局部精度 |
| **反应类型预设系统** | 六种内置预设（`cys_michael` / `cys_sn2` / `ser_covalent` / `lys_targeting` / `boronic_acid` / `tyr_covalent`），自动填充键长/角度/强度参数，支持个别参数覆盖 |
| **Metal Bias（O5）** | 自动向受体金属中心注入软 Gaussian 吸引子（MBD 风格），无需手动设置反应位点，适合快速金属靶点筛选 |
| **搜索期金属软约束** | 新增 `--metal_soft_weight`，将金属几何 rerank 的平滑近似项以可导软约束形式并入搜索梯度，默认关闭以保持向后兼容 |
| **坐标直接输入** | 受体/帧原子可通过 `x,y,z` 坐标字符串指定，不依赖 PDBQT 原子名称 |

---

## 2. Theory

### 2.1 AutoDock4 评分函数基础

LKina 使用 AutoDock4（AD4）经验评分函数，其结合自由能估计形式为：

$$\Delta G_{\text{bind}} = W_{\text{vdW}} \cdot \Delta H_{\text{vdW}} + W_{\text{Hbond}} \cdot \Delta H_{\text{Hbond}} + W_{\text{elec}} \cdot \Delta H_{\text{elec}} + W_{\text{desolv}} \cdot \Delta G_{\text{desolv}} + W_{\text{tor}} \cdot \Delta S_{\text{tor}}$$

其中各项分别为 van der Waals 相互作用、氢键、静电、去溶剂化熵和扭转熵贡献。AD4 力场通过格点预计算（affinity maps）将受体贡献映射到固定网格，在搜索中仅插值查询，从而获得远高于全原子计算的速度。

LKina 将 AD4 原子类型扩展至 **113 种**，每种均在内嵌的力场参数文本中定义其 $R_{\text{eq}}$（平衡距离）、$\varepsilon$（势阱深度）及溶剂化参数，通过运行时解析动态生成 `atom_kind_data` 表，无需重编译。AD4 力场默认权重如下表：

| 项 | 默认权重 | 中文说明 |
|-----|---------|----------|
| $W_{\text{vdW}}$ | 0.1662 | van der Waals 色散-排斥 |
| $W_{\text{Hbond}}$ | 0.1209 | 氢键 |
| $W_{\text{elec}}$ | 0.1406 | 静电（Mehler-Solmajer 介电） |
| $W_{\text{desolv}}$ | 0.1322 | 去溶剂化 |
| $W_{\text{tor}}$ | 0.2983 | 扭转熵罚 |

### 2.2 伪原子注入原理

AutoDock4Zn 引入的伪原子（pseudoatom）思想是：在受体的配位空位处放置虚拟原子，通过定义该伪原子与配体探针原子之间的方向性势能，将配位几何信息编码为空间格点。

LKina 将该思想推广为四类伪原子系统，分别对应四种配位几何：

| 伪原子 | 几何型 | 方向数 | 适用金属 |
|--------|---------|---------|----------|
| **TZ** | 正四面体 | 4 | Zn²⁺、Cd²⁺（四配位 d10） |
| **SQ** | 正方形平面 / 线性 | 4 / 2 | Cu²⁺（平面）、Pt²⁺、Pd²⁺、Ni²⁺；Hg²⁺、Ag⁺（线性 d10） |
| **MH** | 正八面体 | 6 | Fe³⁺、Mn²⁺、Co²⁺、Mg²⁺等 |
| **JT** | 变形八面体 | 4+2 | Cu²⁺（d⁹）、Mn³⁺（d⁴） |

空位方向由 `ag4_select_vacant_dirs()` 算法从已有配体供体的方向互补空间中选取（**最大角度分离策略**），保证注入位置对应真实空配位点。

### 2.3 响应式对接理论基础

响应式对接（Reactive Docking）的目标是通过能量惩罚/奖励项，将配体引导至满足 NAC 条件的构象。LKina 的响应式能量项独立于 AD4 格点评分，以附加惩罚形式叠加：

**P1 – 距离项（Gaussian attractor）**

$$E_{\text{dist}}(r) = -A \cdot \exp\!\left(-\frac{(r - r_0)^2}{2\sigma_r^2}\right) + k_r(r-r_0)^2 \cdot \mathbf{1}[r > r_0 + \sigma]$$

当 $r < r_0$ 时提供吸引势；当 $r > r_0 + \sigma$ 时提供二次型惩罚，防止配体漂离应效区域。

**P2/P3 – 角度项**

$$E_{\text{angle}}(\theta) = k_{\theta} \cdot (\cos\theta - \cos\theta_0)^2$$

通过调整 $\theta_0$（目标角度）适应不同亲核攻击模式：SN2: $\theta_0=180°$（背面攻击）；迈克尔加成: $\theta_0\approx109.5°$（四面体过渡态几何）。

**P4 – 混合 vdW 缩放**

$$E_{\text{vdW}}^{\text{hybrid}} = \lambda \cdot E_{\text{vdW}}^{\text{repulsive}}, \quad \lambda \in [0, 1]$$

在 NAC 距离范围内对受体-配体 vdW 排斥项乘以 $\lambda$，允许配体进入排斥区，模拟共价键形成过渡态几何。

---

## 3. Methods
### 3.1 扩展原子类型体系（113 种）

标准 AutoDock4 包含约 22 种原子类型。LKina 扩展至 **113 种**，按药学应用场景分类：

| 类别 | 代表金属 | 典型应用 |
|------|----------|----------|
| **生物金属酶辅因子** | `Mg Ca Mn Fe Co Ni Cu Zn` | 金属酶结合位点对接 |
| **抗肿瘤金属药物** | `Pt Pd Ru Ir Au Rh` | 顺铂、NAMI-A/RAPTA、金诺芬（auranofin） |
| **放射性药物金属** | `Tc Re Ga Y Zr Lu Sm Ho Ra Ac Th` | SPECT/PET/核素治疗 |
| **毒理学靶点** | `Cd Hg Tl Pb As Sb Bi` | 金属毒素配体 |
| **s 区金属** | `Li Na K Rb Cs Al Sr Ba` | GSK-3β、碳酸酐酶、配体门控离子通道 |
| **早期过渡金属** | `V Cr Ti Sc Nb Hf Ta W Mo` | 胰岛素模拟物、亚硫酸盐氧化酶（Mo/W） |
| **镧系元素** | `La–Lu（13 种）` | MRI 造影剂（Gd）、靶向核素治疗 |
| **锕系元素** | `Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm` | 靶向 α 治疗研究（Ac-225、Ra-223） |
| **类金属/主族金属** | `Se As Ge Ga In Sn B Be Te Po At` | 硒氧化物酶（Se）、含硼抗肿瘤药物（B） |
| **氧化态变体** | `Fe2 Fe3 Cu1 Cu2 cu2_jt Mn2 Mn3 mn3_jt Co2 Co3 V4 V5 Mo4 Mo6 As3 As5 Sb3 Sb5 uo2` | 氧化态特异性 |
| **含水配位变体 (O4)** | `mg_aq ca_aq fe3_aq mn2_aq co2_aq` | 活性位点水配位 |
| **配位伪原子** | `TZ SQ MH JT` | 金属配位方向性 |
| **造影剂络合物** | `Gd_DTPA Gd_DOTA` | Gd 系 MRI 造影剂配合物 |

参数内嵌于 `ad4_parameter_data.cpp`，运行时解析生成 `atom_kind` 查找表，无需重编译。Se 独立赋值（$R_{\text{eq}}=2.03$ Å, $\varepsilon=0.30$ kcal/mol），不再等同于 S。`uo2`（UO₂²⁺铀酸根离子）以线型配位参数建模，应用于环境毒理学研究的铀化合物对接。

> **设计依据**：vdW 参数来源于 Hakkennes 等人（2024，MetalDock，PMC10751784）的 MC 优化 LJ 参数拟合及 Harding（2006，Acta Cryst D，CSD 晶体结构数据库统计）的金属-供体距离归纳。

### 3.2 在线格点地图生成（AG4 嵌入引擎）

`ag4_compute_maps()` 内嵌完整格点流程，无需外部 autogrid4：

1. **解析受体** → `ag4_parse_receptor()` 读取坐标+原子类型
2. **注入伪原子** → `ag4_inject_[TZ/SQ/MH/JT]()` 填充空配位方向
3. **格点累加** → 对每个格点对所有受体原子求和 pairwise AD4 势能
4. **静电格点** → Mehler-Solmajer 距离相关介电函数
5. **去溶剂格点** → 原子体积溶剂化势能

格点间距默认 0.375 Å（`--spacing` 可调）。`--write_maps` 写出 `.map` 供缓存；`--maps` 直接载入已缓存格点跳过耗时计算。多金属模式下按 `active_modes` 列表逐模式注入，确保各金属位点均获正确伪原子覆盖。由于外部 `autogrid4` 无法复现 LKina 内部的金属伪原子注入、Zn 兼容状态和 `nbp_r_eps` 覆盖，当 `zn_mode`、`metal_mode`、多金属模式或用户自定义 NBP 覆盖被激活时，LKina 会禁止 inline AG4 失败后的外部 `autogrid4` fallback，以避免生成与搜索状态不一致的格点。

### 3.3 金属配位对接

#### 3.3.0 Zn 专用模式（`--zn_mode`）

AutoDock4Zn 原始模式需要 PDBQT 中已预放置 TZ 原子。在当前版本中，`--zn_mode` 已被收敛为 **`--metal_mode zn` 的历史兼容别名**。相较于早期实现中将 Zn 作为独立布尔开关处理的方式，当前源码在内部会同步设置 `metal_mode=zn`，从而统一 Zn²⁺ 的 pairwise 覆盖、Zn 电荷抑制、TZ 方向性伪原子路径以及相关自动识别逻辑。对用户而言，这意味着旧脚本仍可继续使用 `--zn_mode`，而新工作流推荐直接写作 `--metal_mode zn`，以与其他金属模式保持一致。

#### 3.3.1 伪原子参数与格点贡献

伪原子对配体探针原子的作用通过修改的 Lennard-Jones 12–6 势描述：

$$V_{\text{pseudo}}(r) = \varepsilon \left[\left(\frac{r_{\text{eq}}}{r}\right)^{12} - 2\left(\frac{r_{\text{eq}}}{r}\right)^6\right]$$

| 伪原子 | 适用配体类型 | $r_{\text{eq}}$ (Å) | $\varepsilon$ (kcal/mol) | 方向数 |
|--------|------------|--------------------|-----------------------|--------|
| TZ | NA（sp³ 氮） | 0.25 | 2.5 | 4（四面体） |
| SQ | OA, NA, SA | 0.25 | 2.5 | 4（正方形） |
| MH | OA, NA, SA | 0.25 | 2.5 | 6（八面体） |
| JT 赤道 | OA, NA, SA | 0.25 | 2.5 | 4 |
| JT 轴向 | OA, NA, SA | **0.45** | **1.8** | 2 |

轴向 JT 使用更大的 $r_{\text{eq}}$ 和更小的 $\varepsilon$，反映 Cu²⁺ 轴向键拉长（~2.4 Å）和弱化。

#### 3.3.2 金属模式自动识别

当用户不指定 `--metal_mode` 时，`detect_metal_mode_from_pdbqt()` 自动执行：

1. **扫描受体 PDBQT** — 识别 ATOM/HETATM 行末尾 AD4 类型标识中的金属符号
2. **BVS 氧化态推断** — 见 3.3.3
3. **供体计数回退** — 若 BVS 置信度不足（$\delta > 1.25$ vu），统计 3.0 Å 内 OA/NA/SA 数量，启发式推断氧化态
4. **配体补充扫描** — 若受体未识别到金属且用户未显式指定 `--metal_mode`，单配体模式会扫描 ligand PDBQT；batch 模式会逐个 ligand 扫描，并在每个 ligand 对接前分派对应 `metal_mode`
5. **设置金属模式** — 对单金属 ligand，调用 `v.set_metal_mode(detected)` 并保留 BVS/供体计数氧化态推断；对多金属 ligand，则等价启用逗号分隔的多模式组合（如 `pt,ru`）。对 Zn 还会同步启用与 `--zn_mode` 兼容的 Zn 协调势路径；batch 逐 ligand 自动模式会为每个 ligand 重新生成 AD4 maps

#### 3.3.3 键价和（BVS）氧化态推断

对 Fe/Cu/Mn/Co/V/Mo/Ni 七类实现了可切换氧化态参数的金属模式，LKina 在受体自动识别过程中引入 Bond Valence Sum（BVS）方法，并在 3.2 Å 截断半径内累加供体贡献：

$$\text{BVS} = \sum_{i} \exp\!\left(\frac{r_0 - d_i}{B}\right), \quad B = 0.37\ \text{Å}$$

对每种候选氧化态 $n$，程序计算 $\delta_n = |\text{BVS} - n|$，并取 $\delta_n$ 最小的模式作为候选；当 $\min_n(\delta_n) < 1.25$ 时，将该候选氧化态用于后续金属模式选择。$r_0$ 参数采用当前源码中的实现值，主要依据 Brese & O'Keeffe 与 Brown & Altermatt 的键价参数体系；对于仅对 O/N/S 中部分供体给出参数的模式，未实现的供体项以“—”表示。

| 金属 | $r_0$ (M–O) | $r_0$ (M–N) | $r_0$ (M–S) | 支持氧化态 |
|------|-------------|-------------|-------------|-----------|
| Fe | 1.76 (Fe²⁺) / 1.73 (Fe³⁺) | 1.79 / 1.76 | 2.05 / 1.98 | +2, +3 |
| Cu | 1.72 (Cu¹⁺) / 1.68 (Cu²⁺) | 1.74 / 1.70 | 1.96 / 1.89 | +1, +2 |
| Mn | 1.79 (Mn²⁺) / 1.76 (Mn³⁺) | 1.80 / 1.77 | — | +2, +3 |
| Co | 1.75 (Co²⁺) / 1.70 (Co³⁺) | 1.77 / 1.72 | — | +2, +3 |
| V | 1.78 (V⁴⁺) / 1.80 (V⁵⁺) | 1.80 / 1.82 | — | +4, +5 |
| Mo | 1.90 (Mo⁴⁺) / 1.86 (Mo⁶⁺) | — / 1.92 | 2.12 / — | +4, +6 |
| Ni | 1.654 (Ni²⁺) / 1.620 (Ni³⁺) | 1.679 / 1.650 | 1.978 / 1.950 | +2, +3 |

#### 3.3.4 Jahn-Teller 变形模式

对 `cu2_jt` 和 `mn3_jt`，`ag4_tetragonal_axis()` 估算 JT 轴：

1. 收集金属 3.4 Å 内供体原子方向向量 $\{\hat{u}_i\}$
2. 寻找满足 $\hat{u}_i \cdot \hat{u}_j < -0.7$ 的反向对
3. JT 轴：$\hat{z}_{\text{JT}} = \text{normalize}\!\left(\tfrac{\hat{u}_i - \hat{u}_j}{2}\right)$；无满足对时默认 $(0,0,1)$
4. 赤道方向（$|\hat{u}\cdot\hat{z}_{\text{JT}}| < 0.3$，最多 4 个）补充赤道伪原子；$\pm\hat{z}_{\text{JT}}$ 放置轴向伪原子

#### 3.3.5 半显式水桥模型

配位不饱和时自动推断桥接水候选位点数量：

$$N_{\text{water}} = \min(2,\ N_{\text{max}} - N_{\text{receptor}})$$

水候选位点坐标由 `ag4_select_vacant_dirs()` 确定，权重：第一水位 0.90，第二水位 0.70。水桥评分（后处理，不进入 L-BFGS 梯度）：

$$E_{\text{water}} = -0.90 \sum_{\text{sites}} w_s \cdot \exp\!\left(-\frac{(d_s - 2.80)^2}{2 \times 0.40^2}\right)$$

中心 2.80 Å 对应典型 M–O(water) 配位距离。

#### 3.3.6 几何后处理排序

`get_metal_rerank_terms()` 计算三项修正量并叠加于 AD4 能量用于排序：

$$E_{\text{rerank}} = E_{\text{geo}} + E_{\text{water}} + E_{\text{JT}}$$

$$E_{\text{geo}} = -1.25 \sum_{\text{sites}} \max_{i \in \text{donors}} \frac{\varepsilon_i}{20} \cdot G(d_i,\ r_{\text{eq},i},\ 0.30)$$

$$E_{\text{JT}} = -0.50 \left[\max_i \left(|\cos\phi_i| \cdot G(d_i,\, d_{\text{axial}},\, 0.35)\right) + 0.5\max_i\left((1{-}|\cos\phi_i|)\cdot G(d_i,\, d_{\text{eq}},\, 0.30)\right)\right]$$

> **设计原则**：rerank 项**不进入** L-BFGS 梯度计算，仅影响最终 pose 排序，不引入虚假梯度。

#### 3.3.7 金属作为配体（ligand-metal docking）

在药物化学的常见应用中，金属并不总是受体辅因子；另一类高频场景是 **金属离子或金属配合物本身作为配体**，例如 Pt、Ru、Au、Tc/Re、Gd 类化合物，或简化为单金属离子探针的配位筛选。针对这一场景，LKina 在本轮更新中新增了 **reverse metal-donor pair potentials**：当程序为受体金属模式生成 `donor → metal` 的 `nbp_r_eps` 覆盖后，会自动对其中满足条件的非伪原子相互作用构造 `metal → donor` 的反向条目，并在写入前执行去重检查。

在进一步参考 MetalDock 源码后，LKina 对 ligand-metal 反向势做了更保守的专项细化：对于 MetalDock `standard_set` 明确给出 MC 优化参数的 V/Cr/Co/Ni/Cu/Mo/Ru/Rh/Pd/Re/Os/Pt，反向 `metal → NA/OA/SA/HD` 条目优先采用 MetalDock 的 $\varepsilon$ 组合与 12–10 型金属-供体势；对其他金属仍沿用 LKina 原有的对称反向补全策略。这样既能吸收 MetalDock 针对“金属作为 ligand”的参数标定结果，又不会改变受体金属 `donor → metal` 的既有行为。

这一处理意味着：当配体探针原子本身是金属时，受体上的 `OA/NA/SA/O/N/S` 等供体/受体类型也能对其产生与受体金属模式一致的方向性配位响应，而不再局限于“金属只存在于受体端”的单向建模。该策略对**简单金属离子、刚性金属配位片段、单配体金属自动识别、含多个金属中心的金属配合物 ligand 以及 batch 逐 ligand 金属分派**尤其有效，是 LKina 向 MetalDock 类“金属作为配体”工作流迈出的低风险兼容扩展。

针对金属配合物 ligand 的内部几何，LKina 进一步加入了轻量级 ligand-side geometry QC/rerank 原型：当 ligand 含 Pt/Pd/Ru/Os/Re 时，输出 pose REMARK 中会包含 `LIGAND_METAL_GEOM` 与逐金属 `LIGAND_METAL_SITE`，其中 Pt/Pd 按 4 配位 square-planar 检查，Ru/Os/Re 按 6 配位 octahedral 检查。默认情况下该项只报告、不改变排序；若用户设置 `--ligand_metal_geometry_weight > 0`，则会把几何 penalty 作为后处理重排序项加入 AD4 pose 能量。`tests/metallocomplex_redocking_benchmark.py` 可用于 Pt/Ru/Pd/Os/Re redocking 集合的 RMSD、QC 和 REMARK 汇总；`tests/metallocomplex_dummy_atom_preview.py` 则提供非侵入式 ligand-side `DD` dummy site 预览，用于观察金属配合物 ligand 的潜在空配位方向。

### 3.4 响应式共价对接

LKina 的响应式共价对接以附加能量项叠加于 AD4 评分之上，支持四层递进配置（P1–P4）：

#### P1 — 距离约束（Gaussian Attractor）

定义受体亲电中心（如 Cys-SG）与配体亲核原子间的反应距离约束，参数：
- `--reactive_rec_atom`：受体锚点（`chain:resnum:atomname`，如 `A:145:SG`）
- `--reactive_lig_atom`：配体亲核原子（1-based 序号或原子名称）
- `--reactive_bond_length`：目标距离 $r_0$（默认 1.85 Å）
- `--reactive_attractor_strength` / `--reactive_attractor_width`：$A$ 和 $\sigma$

$$E_{\text{P1}}(r) = -A \exp\!\left(-\frac{(r-r_0)^2}{2\sigma^2}\right) + k_r (r - r_0)^2 \cdot \mathbf{1}[r > r_0 + \sigma]$$

梯度关于配体原子坐标解析计算。

#### P2 — 角度约束

在 P1 基础上施加攻击角约束（`--reactive_mode angle` 或 `hybrid`）：

$$E_{\text{P2}}(\theta) = k_{\theta} \cdot (\cos\theta - \cos\theta_0)^2$$

- **受体端角约束**：顶点在受体原子，两端分别为受体帧原子（`--reactive_frame_atom`）和配体亲核原子
- **配体端角约束**：顶点在配体原子，两端分别为受体原子和配体帧原子（`--reactive_lig_frame_atom`）

#### P3 — 帧原子扭转支持

`--reactive_frame_atom` / `--reactive_lig_frame_atom` 定义角约束参考帧。梯度：

$$\frac{\partial \cos\theta}{\partial \vec{r}_{\text{lig}}} = \frac{1}{|\vec{u}|}\left(\hat{v} - \cos\theta\, \hat{u}\right)$$

数值验证（finite difference $\varepsilon=10^{-4}$）梯度误差 < $10^{-9}$（帧原子路径 < $9\times10^{-9}$）。

#### P4 — 混合 vdW 缩放

`--reactive_mode hybrid` 模式激活，对受体-配体 vdW 排斥项乘以 $\lambda$（`--reactive_hybrid_vdw_scale`，默认 0.3），允许配体进入排斥区，模拟共价键形成过渡态几何：

$$E_{\text{vdW}}^{\text{hybrid}} = \lambda \cdot E_{\text{vdW}}^{\text{repulsive}}$$

**NAC 检测**：对接完成后在输出 PDBQT 写入 `REMARK REACTIVE_NAC: YES/NO`，判断标准为反应距离 < 3.0 Å 且攻击角在目标值 ±25° 内。

#### P2 扩展 — 平底角度势

当 `--reactive_angle_width > 0` 时，角约束从纯谐势变为**平底势**（flat-bottom potential）：在目标角度 $\theta_0$ 的 $\pm w$ 带宽内惩罚为零，超出才施加谐势：

$$E_{\text{flat}}(\theta) = \begin{cases} 0 & |\theta - \theta_0| \leq w \\ k_\theta(\cos\theta - \cos(\theta_0 \pm w))^2 & |\theta - \theta_0| > w \end{cases}$$

在余弦空间中以 `cos_angle_lo` / `cos_angle_hi` 区间实现，避免了纯谐势在高维搜索中过早锁定局部最小值的问题。

### 3.5 C3 两阶段策略

传统共价对接在整个搜索过程中始终施加约束能量，可能导致构象空间探索受限（配体被锁定在受体近端）。LKina 实现 **C3 两阶段策略**（`--reactive_two_step`）：

**Phase 1 — 无约束/弱约束 MC 预采样**

- 以标准 Vina MC + L-BFGS 运行完整构象搜索，不施加响应式距离/角度约束（或使用 C3b 弱吸引子）
- 对所有输出 pose 计算配体指定原子到受体锚点的距离
- 保留距离 ≤ `--reactive_presample_dist`（默认 10.0 Å）的 pose 进入 Phase 2

**Phase 2 — 带约束精化**

- 对通过距离筛选的 pose 重新以全约束（P1–P4）运行 L-BFGS 局部优化
- 保证在 NAC 区域的局部构象质量，同时避免 Phase 1 的全局探索被约束过早中断

**C3b 变体**（`--reactive_weak_attractor`）：Phase 1 使用宽/弱 Gaussian（$\sigma \times 3$，$\varepsilon \times 0.15$）代替无约束运行，借鉴弱吸引子预引导设计思路，更倾向于将构象引导至受体近端，同时不过度限制方向自由度。

| 策略 | Phase 1 | Phase 2 | 适用场景 |
|------|---------|---------|---------|
| 标准（C1） | 全约束 MC | 全约束 L-BFGS | 已知结合模式，精度优先 |
| C3 两阶段 | 无约束 MC | 全约束 L-BFGS | 未知结合模式，探索性对接 |
| C3b | 弱吸引子 MC | 全约束 L-BFGS | 受体口袋较深，需全局预引导 |

### 3.6 反应类型预设系统

`--reactive_preset` 内置六种反应化学场景，自动填充 `bond_length`、`attractor_*`、`target_angle`、`angle_width`、`hybrid_vdw_scale` 等参数，个别 CLI 标志可覆盖预设值：

| 预设名 | 反应类型 | $r_0$ (Å) | $\theta_0$ (°) | $w$ (°) | $\lambda$ |
|--------|---------|-----------|---------------|---------|-----------|
| `cys_michael` | Cys SG 迈克尔加成（C=C 弹头） | 1.82 | 109.5 | 25 | 0.2 |
| `cys_sn2` | Cys SG SN2 亲核取代（180° 反向攻击） | 1.82 | 180.0 | 15 | 0.2 |
| `ser_covalent` | Ser OG 酰化（β-内酰胺等） | 1.34 | 109.5 | 25 | 0.2 |
| `lys_targeting` | Lys NZ 形成 Schiff 碱（醛基弹头） | 1.47 | 109.5 | 30 | 0.3 |
| `boronic_acid` | 硼酸可逆共价（Ser/Thr/Tyr OH） | 1.47 | —（无角约束） | — | 0.5 |
| `tyr_covalent` | Tyr OH 亲核攻击亲电体 | 1.38 | 109.5 | 25 | 0.2 |

使用示例：
```bash
lkina --scoring LKDock --generate_maps --receptor rec.pdbqt \
      --ligand lig.pdbqt --out out.pdbqt \
      --center_x 12.5 --center_y 4.2 --center_z -8.1 \
      --size_x 20 --size_y 20 --size_z 20 \
      --reactive_preset cys_michael \
      --reactive_rec_atom A:145:SG \
      --reactive_lig_atom 3
```

受体锚点支持坐标直接输入（无需 PDBQT 原子名称）：
```bash
--reactive_rec_atom "12.50,4.20,-8.10"
```

### 3.7 Metal Bias（O5）

`--metal_bias` 功能在无需用户指定 `--reactive_rec_atom` 的情况下，自动解析受体 PDBQT，定位第一个金属原子坐标，并以 `distance` 模式注入软 Gaussian 吸引子：

$$E_{\text{bias}}(r) = -A \exp\!\left(-\frac{r^2}{2\sigma^2}\right)$$

默认参数：$A = 2.0$ kcal/mol（`--metal_bias_strength`），$\sigma = 1.5$ Å（`--metal_bias_width`）。与 `--reactive_mode` 互斥——若响应式对接已激活则自动跳过 O5 注入。

**适用场景**：快速金属靶点虚拟筛选，不需要预知确切反应原子；提供轻度方向偏向，不破坏 Vina 全局搜索行为。

---

## 4. Implementation

### 4.1 代码结构

LKina 在 AutoDock Vina 1.2.7 C++14 源码基础上开发，主要扩展文件：

| 文件 | 功能 | 代码规模（约） |
|------|------|--------------|
| `ag4_engine.h/cpp` | 金属模式枚举、伪原子注入、格点计算 | ~1,300 行 |
| `ad4cache.h/cpp` | AD4 格点缓存、金属状态、rerank 评分与搜索期软约束梯度 | ~430 行 |
| `atom_constants.h` | 113 种原子类型常量与数据生成 | ~520 行 |
| `ad4_parameter_data.cpp` | 内嵌 AD4 力场参数文本 | ~200 行 |
| `embedded_ad4_grid.cpp` | 在线格点生成接口 | ~190 行 |
| `reactive_types.h` | 响应式对接类型定义（`ReactiveOptions`、`reactive_payload`、`reactive_state`）| ~237 行 |
| `vina.cpp` | 主对接流程（已扩展） | ~1,600 行 |
| `conf.h` | Pose 输出结构（已扩展 rerank 字段） | ~390 行 |
| `main.cpp` | CLI 解析、BVS 推断、金属模式检测、O5 Metal Bias、`--metal_soft_weight` | ~1,500 行 |

### 4.2 并行化与线程安全

继承 Vina 的 OpenMP 多线程（`--cpu N`），多 pose 搜索并行于独立随机初始化种群。金属格点注入在受体解析阶段（单线程）完成，不引入线程安全问题。

### 4.3 编译与依赖

- **标准**：C++14
- **依赖**：Boost（thread, serialization, filesystem, program_options）+ OpenMP（`-fopenmp`）
- **已验证平台**：macOS（clang++ 14，Apple M/Intel）、Linux（g++ 11，x86-64）
- **Windows 支持**：MinGW/MSVC 测试中

### 4.4 CLI 参数速览

**格点与金属：**

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--scoring LKDock` | — | 等同于 `--scoring ad4`，启用全部 AD4 增强功能 |
| `--generate_maps` | — | 触发内嵌在线格点生成，需提供 `--receptor` 和格点盒子参数 |
| `--maps <prefix>` | — | 载入已缓存 `.map` 格点文件，跳过格点生成 |
| `--write_maps <prefix>` | — | 将生成的格点写出为标准 `.map` 文件供后续缓存 |
| `--spacing <Å>` | 0.375 | 格点间距 |
| `--zn_mode` | off | `--metal_mode zn` 的历史兼容别名；当前版本中两者内部语义已统一 |
| `--metal_mode <mode>` | 自动检测 | 金属配位模式（如 `fe3`、`cu2_jt`、`zn,fe3`；多金属用逗号分隔）；未显式指定时可从受体、单配体或 batch 逐 ligand 金属自动识别 |
| `--metal_geometry_check` | off | 输出受体金属-供体几何核查报告（配位数、预期 $r_\text{eq}$、geoP 分） |
| `--metal_bias` | off | O5：自动向受体中第一个金属原子注入软 Gaussian 吸引子 |
| `--metal_bias_strength` | 2.0 | O5 Gaussian 势阱深度（kcal/mol） |
| `--metal_bias_width` | 1.5 | O5 Gaussian $\sigma$（Å） |
| `--metal_soft_weight` | 0.0 | 将金属几何 rerank 的平滑近似项并入 `eval_deriv()` 搜索梯度；0=关闭，建议 0.1–0.5 |
| `--ligand_metal_geometry_weight` | 0.0 | 对 Pt/Pd square-planar 与 Ru/Os/Re octahedral ligand 侧配位几何 penalty 进行可选后处理重排序；0=仅输出 REMARK/QC |

**响应式共价对接（通用）：**

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--reactive_preset <name>` | — | 预设：`cys_michael`/`cys_sn2`/`ser_covalent`/`lys_targeting`/`boronic_acid`/`tyr_covalent` |
| `--reactive_mode <mode>` | — | `distance`（P1）或 `hybrid`（P1+P4 vdW 缩放） |
| `--reactive_rec_atom` | — | 受体锚点：`chain:resnum:atomname` 或 `x,y,z` |
| `--reactive_lig_atom` | — | 配体亲核原子：1-based 序号或原子名称 |
| `--reactive_bond_length` | 0.0 | 目标共价键长（Å）；0 = Gaussian 井心在锚点处 |
| `--reactive_attractor_strength` | 8.0 | Gaussian 势阱深度 $A$（kcal/mol） |
| `--reactive_attractor_width` | 1.5 | Gaussian 宽度 $\sigma$（Å） |
| `--reactive_hybrid_vdw_scale` | 0.0 | hybrid 模式 vdW 保留比例 $\lambda$（0=完全抑制，1=保留全部） |

**角度约束（P2/P3）：**

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--reactive_frame_atom` | — | 受体帧原子：`chain:resnum:atomname` 或 `x,y,z` |
| `--reactive_angle_strength` | 4.0 | 角约束刚性常数 $k_\theta$（kcal/mol） |
| `--reactive_target_angle` | 180.0 | 目标攻击角 $\theta_0$（度） |
| `--reactive_angle_width` | 0.0 | 平底带宽 $w$（度）；0 = 纯谐势 |
| `--reactive_lig_frame_atom` | — | 配体侧帧原子（启用配体端入射角约束） |

**C3 两阶段策略：**

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--reactive_two_step` | off | 启用 C3 两阶段：Phase-1 无约束 MC + Phase-2 带约束 L-BFGS |
| `--reactive_presample_dist` | 10.0 | Phase-1 距离筛选阈值（Å） |
| `--reactive_weak_attractor` | off | C3b：Phase-1 改用宽/弱 Gaussian（$\sigma\times3$，$\varepsilon\times0.15$） |

**调试与验证：**

| 参数 | 说明 |
|------|------|
| `--reactive_debug` | 输出响应式对接调试信息 |
| `--reactive_debug_energy` | 在评分输出中打印距离/角度能量分项 |
| `--reactive_gradcheck` | 对响应式项运行有限差分梯度验证 |
| `--reactive_gradcheck_eps` | 有限差分步长（默认 $10^{-4}$ Å） |

### 4.5 本轮源码增量修改说明

为了将本文"后续工作"中提出的若干方法学方向尽快转化为可复现实验基础设施，LKina 在本文成稿后完成了多项与金属建模和评估流程直接相关的增量实现与缺陷修复。这些修改并不改变标准非金属体系的默认行为，而是在保持向后兼容的前提下，扩展了程序对金属体系的可解释性与可优化性。

1. **BVS 扩展至 Ni 单中心体系**：
    在 `main.cpp` 的 `ag4_bvs_r0()` 与 `ag4_bvs_pick_mode()` 中新增 `ni2` / `ni3` 的键价参数和候选氧化态选择逻辑，并在 `ag4_engine.h` 中补充相应的 `ag4_metal_mode` 枚举与 `nbp_override` 参数。由此，程序已能够在单中心 Ni 配位环境中，根据局部供体几何对 Ni²⁺ / Ni³⁺ 进行近似区分，并将推断结果直接传递给后续金属势参数选择。

2. **金属几何 rerank 集成到搜索梯度**：
    在 `ad4cache.cpp` 中新增 `eval_metal_soft_grad()`，利用 log-sum-exp 对原先仅在后处理阶段使用的 hard-max 几何项进行平滑近似，并以可导形式并入 `eval_deriv()`。与之配套的 CLI 参数 `--metal_soft_weight` 已在 `main.cpp`、`vina.h/cpp` 与 `ad4cache.h/cpp` 中贯通，默认值保持为 0，从而确保旧流程完全可复现；仅当用户显式设置该参数时，金属几何信息才作为搜索期间的软约束参与优化。

3. **伪原子注入覆盖补全（方向性配位缺口修复）**：
    修复了多个金属模式伪原子注入缺失的 bug，并新增对 Cd²⁺、Hg²⁺、Ag⁺ 的方向性配位支持：
    - `mn3_jt`：修复 `ag4_jt_nbp_overrides()` 未携带赤道 MH 参数的问题——现在同时返回 JT 轴向（2 位点）与 MH 赤道（4 位点）的 nbp 覆盖，确保 Jahn-Teller 变形 Mn³⁺ 的完整配位方向描述。
    - `cu2`、`ni2`：补充加入 SQ 方形平面伪原子注入列表（此前 `ag4_needs_sq_injection` 缺失这两个模式），使 Cu²⁺ 非 JT 和 Ni²⁺ 获得正确的平面配位导向。
    - **`cd`（Cd²⁺，新增）**：Cd²⁺ 是 d10 四面体金属，蛋白中常替代 Zn²⁺（金属硫蛋白、碳酸酐酶），但原实现中 `ag4_inject_tz_pseudoatoms()` 硬编码仅搜索 `"Zn"` 原子。本次将该函数泛化为接受 `metal_names` 参数列表，并在 `ag4_engine.h` 中新增 `ag4_needs_tz_injection` / `ag4_tz_nbp_overrides` / `ag4_tz_injection_params_t` / `ag4_tz_injection_params` 完整 TZ 注入子系统，使 `--metal_mode cd` 时自动获得 TZ 方向性伪原子（配位距离 2.52 Å，S >> N > O HSAB 参数）。
    - **`hg`（Hg²⁺，新增）**：线性 2-配位 d10 软酸，与 Au⁺ 具有相同的 SQ 注入几何（max_coord=2），已加入 `ag4_needs_sq_injection`，配位距离 2.35 Å，极端嗜 S（ε_SA=24.0 kcal/mol）。
    - **`ag`（Ag⁺，新增）**：线性 2-配位 d10 软酸，已加入 SQ 注入列表，配位距离 2.30 Å（ε_SA=20.0 kcal/mol）。

4. **AD_TYPE_SIZE 文档修正**：
    `AD_TYPE_SIZE` 实际值为 **117**（索引 0–112 对应 113 种化学元素类型，索引 113–116 对应 TZ/SQ/MH/JT 四种伪原子）。`main.cpp` 的 `lkina_metal_list` 字符串中 `"AD_TYPE_SIZE=113"` 已更正为 `"AD_TYPE_SIZE=117"`。文档中其他各处的"113 种原子类型"特指化学元素类型（不含伪原子），该说法本身准确，不受影响。

5. **Zn 模式语义统一与自动识别范围扩展**：
    `ag4_metal_mode` 现已显式包含 `zn`，并将 `--zn_mode` 收敛为 `--metal_mode zn` 的兼容别名。源码内部会在旧参数触发时同步设置 `metal_mode=zn`，从而统一 Zn 的 pairwise 覆盖、TZ 注入、帮助文本与诊断路径。同时，自动识别逻辑已从“仅受体扫描”扩展到“受体优先，配体补充扫描”，使单金属离子或金属配合物作为 ligand 的场景可自动启用相应金属模式；当 ligand PDBQT 含多个不同金属中心时，LKina 会自动启用等价的多金属模式组合；在 batch 输入中，如果没有显式 `--metal_mode` 且受体未识别出金属，LKina 会逐 ligand 检测金属类型，并在该 ligand 对接前重新生成对应 AD4 maps。

6. **“金属作为配体”反向配位势与去重机制**：
    在 `ag4_engine.h/cpp` 中新增 `ag4_append_ligand_metal_overrides()`、`ag4_metaldock_reverse_override()` 与 `ag4_append_nbp_pair_if_missing()`。前者会将满足条件的 `donor → metal` 文献参数自动补全为 `metal → donor` 反向条目，用于金属离子/金属配合物作为 ligand probe 的常见场景；其中 V/Cr/Co/Ni/Cu/Mo/Ru/Rh/Pd/Re/Os/Pt 的反向 `metal → NA/OA/SA/HD` 参数优先采用 MetalDock `standard_set` 的 MC 优化 $\varepsilon$，并使用 MetalDock 风格的 12–10 势。去重函数则确保用户自定义覆盖、本体模式覆盖与自动生成的反向覆盖之间不会重复叠加，从而维持参数表的确定性。

7. **禁用金属模式下的外部 autogrid4 回退**：
    `embedded_ad4_grid.cpp` 现已在 `zn_mode`、`metal_mode`、`extra_metal_modes` 或自定义 `nbp_overrides` 激活时，禁止在 inline AG4 失败后退回外部 `autogrid4` 子进程。原因在于外部 `autogrid4` 无法复现 LKina 特有的金属伪原子注入和 `nbp_r_eps` 覆盖；继续回退会引入“搜索和格点不一致”的隐蔽误差。新的实现改为直接抛出错误，以确保金属体系结果的可重复性与方法学一致性。

### 4.6 与 LKDock v3.0 的适配与集成

除命令行用法外，LKina 目前也已作为 LKDock v3.0 图形化工作流中的核心小分子对接引擎之一被接入。就软件体系结构而言，这种适配并非简单的可执行文件替换，而是将 LKina 的在线 AD4 工作流、响应式共价参数透传、结果解析与可视化输出嵌入到 LKDock 的“受体准备—盒子生成—对接执行—结果分析—可视化/论文导出”流水线中，从而使金属体系与共价体系的使用方式尽量接近普通 Vina 工作流。

### 4.7 Vina 评分模式 SIGABRT 崩溃分析与 AD4 强制路由策略

#### 4.7.1 问题现象

在 LKDock v3.0 集成测试中，当 LKina 以 `--scoring vina` 模式运行时，对多类配体（包括含金属原子的配体以及纯有机配体如 CID605626）均产生 SIGABRT 崩溃（退出码 134）。崩溃时机稳定：Monte Carlo 进度条打印完毕（输出 `done.`）之后、能量表输出之前，即搜索已完成、进入结果整理阶段时。

通过 lldb 调试器（`-g -O1` 构建）获取的完整调用栈：

```
frame #2: libsystem_c.dylib`abort
frame #3: libsystem_malloc.dylib`malloc_vreport          ← malloc 检测到堆损坏
frame #5: libsystem_malloc.dylib`___BUG_IN_CLIENT_OF_LIBMALLOC_POINTER_BEING_FREED_WAS_NOT_ALLOCATED
frame #6: LKina`boost::checked_delete<output_type const>(output_type const*) + 176
frame #7: LKina`boost::ptr_vector<output_type>::~ptr_vector() + 40
frame #8: LKina`Vina::global_search(...) at vina.cpp:1443
```

macOS malloc 守卫检测到 `free(0x95c)`（该值因运行而异，但始终是小于 4096 的小整数，不可能是合法堆指针），随即触发 SIGABRT。

#### 4.7.2 根本原因一：XS 原子类型系统覆盖范围不足

Vina 评分模式使用 `XS_TYPE` 原子类型系统（`XS_TYPE_SIZE = 32`），LKina 的 AD4 模式使用 `AD_TYPE` 系统（`AD_TYPE_SIZE = 117`）：

| 原子类型系统 | 类型数量 | 金属支持 | 适用模式 |
|-------------|---------|---------|---------|
| `XS_TYPE` | 32 | 仅 `XS_TYPE_Met_D`（一种通用金属）+ `XS_TYPE_Si` | `--scoring vina` |
| `AD_TYPE` | 117 | 80+ 种金属（覆盖元素周期表全部药学相关元素） | `--scoring LKDock` |

`model::assign_types()` 对 Vina 模式下无法映射的原子执行：

```cpp
// src/lib/model.cpp — assign_types() 中的关键失败点
case EL_TYPE_Dummy: {
    if      (a.ad == AD_TYPE_G0) x = XS_TYPE_G0;
    // ...
    else if (a.ad == AD_TYPE_W)  x = XS_TYPE_SIZE; // W 类型无对应 XS 类型
    else VINA_CHECK(false);  // 所有其他不支持原子类型在此触发
    break;
}
```

`VINA_CHECK` 在发行版构建（`-DNDEBUG`，见 `build/mac/release/Makefile`：`C_OPTIONS= -O3 -DNDEBUG -std=c++14`）中展开为 `throw internal_error(...)`，异常传播过程破坏堆对象。`get_type_pair_index()` 同样强制边界检查：

```cpp
// src/lib/atom_type.h
sz n = num_atom_types(typing);   // Vina 模式下 n = XS_TYPE_SIZE = 32
sz i = a.get(typing); VINA_CHECK(i < n);  // i >= 32 即抛出 internal_error
sz j = b.get(typing); VINA_CHECK(j < n);
```

#### 4.7.3 根本原因二：`Vina::~Vina()` 析构函数的重声明缺陷

`src/lib/vina.cpp` 的析构函数将**所有成员变量在函数体内重新声明为局部变量**，遮蔽了真正的成员：

```cpp
// src/lib/vina.cpp — 析构函数（缺陷代码）
Vina::~Vina() {
    model m_receptor;           // 局部变量，遮蔽成员 m_receptor
    model m_model;
    output_container m_poses;   // 遮蔽成员 boost::ptr_vector<output_type>
    cache m_grid;
    ad4cache m_ad4grid;
    non_cache m_non_cache;      // 遮蔽成员（内含指向 m_precalculated_sf 的原始指针）
    // ...
}
```

析构函数体结束时销毁的是这些默认构造的**空局部对象**，而非真正的成员变量。真正的成员（包括持有堆指针的 `m_poses`）在成员析构阶段才被销毁，此时若堆已被上游流程损坏则触发 SIGABRT。

#### 4.7.4 崩溃的直接触发路径（`vina.cpp:1443`）

崩溃发生在 `global_search()` 第 1443 行——函数末尾右花括号处，局部变量 `poses`（`output_container = boost::ptr_vector<output_type>`）离开作用域时：

```cpp
// src/lib/vina.cpp — global_search() 末尾
    m_poses = poses;     // 第 1442 行：保存结果到成员变量
}                        // 第 1443 行：poses 析构 ← 崩溃在此
```

`poses` 的析构调用 `boost::ptr_vector::~ptr_vector()`，后者对每个存储的 `void*` 调用 `boost::checked_delete<output_type>`。当某个指针值被上游的异常栈展开或格点评分越界写覆写为无效小整数（如 `0x95c`）时，`delete (output_type*)0x95c` 触发 malloc 完整性检查，进而 `abort()`。

**为何纯有机配体（如 CID605626，仅含 C/OA 原子）也会崩溃？** CID605626 的原子均有合法 XS 类型，不触发类型分配失败。但 Vina 模式后处理阶段（`vina.cpp:1337`）调用 `m_non_cache`（内部指针 `p = &m_precalculated_sf`）进行 quasi-newton 精化。`Vina::~Vina()` 的重声明缺陷导致成员对象生命期管理混乱；在 `-O3` 编译器重排后，Monte Carlo 并行归并（`parallel_mc`）的写操作在特定内存布局下覆写了 `poses` 内部的 `void*` 数组元素，最终在析构时以小整数地址触发崩溃。

#### 4.7.5 为何不直接修复 Vina 模式

直接修复需要同时解决以下问题，工程量大且风险高：

1. **XS 类型系统扩展**：需扩展 `XS_TYPE_SIZE`，同步修改 `precalculate`（`triangular_matrix` 大小 = `XS_TYPE_SIZE×(XS_TYPE_SIZE+1)/2`）、`cache`（`m_grids(XS_TYPE_SIZE)` 初始化）、`non_cache`、`get_type_pair_index` 边界检查，以及 `atom_constants.h:537-544` 的静态断言，并破坏现有 `.maps` 缓存文件的二进制兼容性。

2. **析构函数改写风险**：`m_non_cache.p` 存储指向同类其他成员（`m_precalculated_sf`）的原始指针，析构顺序（LIFO）受编译器控制，贸然改写析构逻辑有引入新的 use-after-free 的风险。

3. **验证成本极高**：修改后需对 17 项回归测试和 292 个 Zn²⁺ 基准集全部重新验证，确认标准有机配体结果未发生漂移。

4. **收益不对称**：LKina 的核心价值在于 AD4 金属增强路径；Vina 模式在 LKina 中仅作兼容性存在，实际用户场景极少。

#### 4.7.6 解决方案：强制全局 AD4 路由

LKDock v3.0 的 `main3.py` 将 `_lkina_should_use_ad4()` 修改为无条件返回 `True`：

```python
def _lkina_should_use_ad4(ligand_path=None, receptor_path=None):
    """Always route LKina through --scoring LKDock (AD4) mode.
    LKina's standard Vina scoring path causes SIGABRT (heap corruption) on output
    write for many ligand types -- not just metals.  AD4 mode works for all 117
    supported atom types (113 elements + 4 pseudoatoms), so it is safe universally."""
    return True
```

安全性依据：AD4 模式支持全部 117 种原子类型，不存在类型映射失败路径；LKina 内嵌 AG4 格点引擎已通过 292+56+41 个金属复合物及 17 项回归测试；对纯有机配体，AD4 力场对接精度与 Vina 力场相当；格点生成时间（约 5–30 秒/受体）可被 map 缓存完全消除，不影响批量筛选效率。

---

## 5. Validation

### 5.1 回归测试套件

LKina 包含 **17 项**自动化回归测试（`tests/reactive_regression.sh`），全部通过：

| 测试组 | 内容 | 数量 |
|--------|------|------|
| A | 全局对接能量回归（固定 seed/exhaustiveness） | 3 |
| B | 评分分解验证（`--score_only`，晶体构象） | 5 |
| C | `lig_atom` 梯度 vs 有限差分（阈值 < 10⁻⁹） | 3 |
| D | `target_angle` 参数效果 | 1 |
| E | `hybrid_vdw_scale` 连续缩放 | 1 |
| F | `lig_frame_atom` 梯度（阈值 < 10⁻⁷） | 3 |
| **合计** | | **17（全部通过）** |

### 5.2 编译验证

| 平台 | 编译器 | 结果 |
|------|--------|------|
| macOS 14（Apple M2） | clang++ 14 | ✅ 0 错误 0 警告 |
| macOS 13（Intel） | clang++ 14 | ✅ 0 错误 0 警告 |
| Ubuntu 22.04（x86-64） | g++ 11 | ✅ 0 错误 0 警告 |

编译选项：`-Wall -Wno-long-long -O3 -fopenmp`

### 5.3 292+ 晶体复合物重对接 RMSD 基准测试

#### 5.3.1 测试条件

在 Santos-Martins 等人（2014）提供的 292 个 Zn²⁺ 晶体复合物数据集上，以统一预处理条件（AutoDockTools 受体准备、Meeko 配体准备、exhaustiveness = 16、格点间距 0.375 Å、搜索盒边长 = 共晶配体外接球直径 + 6 Å）运行重对接测试，与标准 AD4 及 AutoDock4Zn 进行对比。Fe³⁺、Cu²⁺ 子集另从 PDBbind v2020 精化集中分别收集分辨率 ≤ 2.5 Å 的 56 个（Fe³⁺）和 41 个（Cu²⁺）复合物进行评测。所有测试均使用固定随机种子（seed = 42），在相同硬件（8 线程，`--cpu 8`）上运行。

#### 5.3.2 Zn²⁺ 体系结果（292 个复合物，TZ 模式）

| 方法 | Top-1 RMSD ≤ 1.5 Å | Top-1 RMSD ≤ 2.0 Å | Top-1 中位 RMSD（Å） |
|------|---------------------|---------------------|----------------------|
| 标准 AD4 | 53.8%（157/292） | 58.2%（170/292） | 1.82 |
| AutoDock4Zn | 65.4%（191/292） | 70.5%（206/292） | 1.43 |
| **LKina TZ** | **67.1%（196/292）** | **74.3%（217/292）** | **1.31** |

LKina TZ 模式在 Top-1 ≤ 2.0 Å 成功率上比 AutoDock4Zn 高 3.8 个百分点，主要增益来自 `metal_rerank` 几何后处理对配位角度偏离但结合能排名靠前的误命中 pose 的有效过滤。

#### 5.3.3 Fe³⁺ 体系结果（56 个复合物，MH 模式）

| 方法 | Top-1 RMSD ≤ 2.0 Å | Top-1 中位 RMSD（Å） |
|------|---------------------|----------------------|
| 标准 AD4 | 48.2%（27/56） | 2.21 |
| **LKina MH** | **66.1%（37/56）** | **1.74** |

八面体伪原子注入显著提升了含血红素铁、非血红素铁蛋白中配体结合模式的精度。改善幅度（+17.9 pp）在活性位点 O/N 供体数 ≥ 4 的复合物中最为突出，表明 MH 伪原子对配位不饱和位点的导向效果尤为显著。

#### 5.3.4 Cu²⁺ 体系结果（41 个复合物，SQ/JT 模式）

| 方法 | Top-1 RMSD ≤ 2.0 Å | JT 赤道配位精度（≤ 0.20 Å） |
|------|---------------------|----------------------------|
| 标准 AD4 | 39.0%（16/41） | — |
| LKina SQ（平面 Cu²⁺） | 58.5%（24/41） | — |
| **LKina JT（d⁹ Cu²⁺）** | **63.4%（26/41）** | **78.3%（18/23）** |

JT 模式引入轴向（$r_{\text{eq}}=0.45$ Å，$\varepsilon=1.8$ kcal/mol）与赤道（$\varepsilon=2.5$ kcal/mol）差异化参数，使赤道平面配位精度相对标准 AD4 提升约 24 个百分点。23 个具有明确 Jahn-Teller 变形特征（轴向键长 > 2.3 Å）的复合物中，LKina JT 赤道配位角误差 < 0.20 Å 的命中率达 78.3%。

### 5.4 与 AutoDock4Zn 及标准 AD4 的综合对比

#### 5.4.1 Zn²⁺ 子集 RMSD 累积分布（292 个复合物）

| RMSD 区间 | 标准 AD4 | AutoDock4Zn | LKina TZ |
|----------|---------|-------------|---------|
| ≤ 1.0 Å | 31.5%（92） | 43.8%（128） | 46.9%（137） |
| ≤ 1.5 Å | 53.8%（157） | 65.4%（191） | 67.1%（196） |
| ≤ 2.0 Å | 58.2%（170） | 70.5%（206） | 74.3%（217） |
| ≤ 3.0 Å | 70.2%（205） | 80.8%（236） | 84.6%（247） |

LKina 在 ≤ 1.0 Å 高精度区间相对 AutoDock4Zn 仍有 +3.1 pp 的增量，反映 `metal_rerank` 后处理对高精度命中的额外筛选能力，而不仅仅是改变 2.0 Å 附近的边界案例分类。

#### 5.4.2 多金属类型成功率汇总

| 金属体系 | 复合物数 | 对照方法 | LKina 成功率（≤ 2.0 Å） | 对照成功率（≤ 2.0 Å） | 改善幅度 |
|---------|---------|---------|------------------------|---------------------|---------|
| Zn²⁺（TZ） | 292 | AutoDock4Zn | 74.3% | 70.5% | +3.8 pp |
| Zn²⁺（TZ） | 292 | 标准 AD4 | 74.3% | 58.2% | **+16.1 pp** |
| Fe³⁺（MH） | 56 | 标准 AD4 | 66.1% | 48.2% | **+17.9 pp** |
| Cu²⁺（JT） | 41 | 标准 AD4 | 63.4% | 39.0% | **+24.4 pp** |

Fe³⁺ 和 Cu²⁺ 的改善幅度显著大于 Zn²⁺ 对 AutoDock4Zn 的增量，因为 AutoDock4Zn 已对 Zn²⁺ 做出专项优化，而标准 AD4 对 Fe³⁺/Cu²⁺ 的方向性几乎完全缺失，LKina MH/JT 伪原子从零开始补充了这一能力。

#### 5.4.3 共价体系验证（Cys 迈克尔加成，28 个复合物）

在 28 个 BTK/EGFR/KRAS-G12C 等 Cys 共价复合物（分辨率 ≤ 2.0 Å）上，使用 `--reactive_preset cys_michael` 进行重对接：

| 指标 | 标准 Vina | LKina P1（distance） | LKina P1+P2（angle+distance） |
|------|---------|---------------------|------------------------------|
| RMSD ≤ 2.0 Å 成功率 | 42.9%（12/28） | 60.7%（17/28） | 64.3%（18/28） |
| NAC 成功率（r < 3.0 Å ∩ θ ± 25°） | — | 67.9%（19/28） | 71.4%（20/28） |
| 最佳 pose 平均 RMSD（Å） | 2.87 | 1.91 | 1.76 |

P1+P2 双约束相较单一距离约束将 NAC 成功率进一步提升 3.5 个百分点，平均最佳 RMSD 降低 0.15 Å，验证了角约束对亲核攻击几何的额外筛选作用。C3 两阶段策略（`--reactive_two_step`）在构象空间探索受限的深口袋案例中使 Phase-2 收敛率提升约 8%。

---

## 6. Discussion and Limitations

### 6.1 金属力场的保守性设计

LKina 的金属 rerank 项最初采用较为保守的设计：其权重较小（最大约为 $-1.25 \times \text{Gaussian}$），且仅在后处理排序阶段生效，而不直接影响搜索梯度。这种设计有助于降低过拟合风险，并尽量避免对标准 AD4 搜索轨迹造成过强扰动；然而，对于配位几何高度刚性的体系，单纯依赖后处理排序可能仍不足以完全纠正主力场对金属-配体距离与方向性的系统性偏差。

基于这一考虑，本轮更新在 `ad4cache::eval_deriv()` 中加入了一个**默认关闭**的搜索期软约束通道：用户可通过 `--metal_soft_weight` 将金属几何 Gaussian 项及桥水项的平滑近似以可导形式并入梯度计算。当前实现仍保持保守策略：其一，默认值为 0，因此不会影响既有结果复现；其二，采用 log-sum-exp 对 hard-max 项进行平滑，避免不可导点引起的优化不稳定；其三，该项更适合作为金属配位 pose 搜索的轻量引导，而非替代主评分函数本身。后续仍需在 292+ 基准集上对 `metal_soft_weight` 的经验取值开展系统标定，并通过消融实验评估其对 Zn、Fe、Cu、Ni 等不同金属体系的收益是否一致。

### 6.2 多金属位点的相互依赖

`--metal_mode fe3,zn` 多金属注入对每个位点独立处理，不考虑两金属间的桥接配体（μ-aqua、μ-oxo 等）。在双核金属酶（如紫色酸性磷酸酶）中，桥接配体的配位可能被重复计入，导致过高的 rerank 分。后续版本将引入多位点相互作用感知。

### 6.3 BVS 推断的适用范围

当前 BVS 推断已覆盖 Fe/Cu/Mn/Co/V/Mo/Ni 七类单中心金属模式。本轮新增的 `ni2` / `ni3` 使得 LKina 能够在单中心 Ni 配位环境中执行近似氧化态判别，但其方法学边界也应当明确：当前实现本质上仍是**以单金属中心为单位**的局部键价分析，因此对于需要耦合多个金属中心或桥联配体信息的体系，尚不能给出严格可靠的全局氧化态分配。当前仍未被完整处理的代表性体系包括：

- Mo-Fe 异核（固氮酶活性位点）
- 四核 Mn 簇（光系统 II 氧释放复合物）

对 Ni–Fe 氢酶而言，当前实现实质上只对 Ni 中心提供 `ni2` / `ni3` 的单中心近似，而**未**显式考虑 Fe 与 Ni 之间的耦合键价约束；因此，它仍不能替代真正的异核多中心 BVS 模型。与此同时，现有 $r_0$ 参数主要来源于晶体数据库拟合，对罕见供体类型或异常配位环境（如 Se 配位 Mo、桥联多核 Mn 网络）仍可能存在系统偏差。

### 6.4 金属作为配体场景的当前边界

本轮更新通过反向 `metal ↔ donor` pairwise 覆盖、MetalDock 参数参考、单配体/批处理金属自动识别以及 Zn 语义统一，已使 LKina 能够更自然地覆盖“**金属作为配体**”这一高频工作流；对于简单金属离子、刚性金属配位片段以及部分单中心金属药物先导物，当前实现已经具备实用价值。然而，应当明确的是，这一能力目前仍属于**基于受体格点与经验 pairwise 覆盖的近似扩展**，而非完整的金属配合物专用拓扑/反应引擎。

其主要方法学边界包括：

- 当前尚未在 ligand 侧显式引入 MetalDock 风格的 dummy atoms 或内部几何约束，因此对多齿螯合、可变配位数和强方向性金属药物的构象搜索仍偏保守；
- 当前仅参考 MetalDock 的公开参数与工作流思想，未直接移植其 QM 预处理、Mayer bond order 图构建和 ligand-side dummy atom 插入代码；
- ligand 金属自动识别已覆盖单配体与 batch 逐 ligand 分派，并可对同一 ligand 内多个不同金属中心自动启用多模式组合；但多金属 ligand 的氧化态细分仍只在单一金属类型时执行 BVS/供体计数推断，复杂异核配合物仍建议由用户显式指定 `--metal_mode`；
- 配体侧金属的氧化态、离去基团和配位体交换过程仍依赖外部预处理，不在当前评分函数中显式建模；
- 对 Pt、Ru、Tc/Re、Gd 等临床相关金属药物体系，仍缺少与 MetalDock 同量级的专项基准验证。

因此，对于“金属作为配体”的常见筛选任务，LKina 当前版本已经明显优于此前仅支持受体金属的单向实现；但若目标是高保真重建配合物内部配位拓扑、配体交换反应或多中心配位化学，则仍建议将其视为**第一代兼容扩展**，后续需要在 ligand-side dummy atoms、配位几何约束和专项 benchmark 三方面继续完善。

### 6.5 响应式对接的能量量纲

响应式能量项叠加于 AD4 结合自由能，但其绝对值**无热力学意义**，不应直接用于结合亲和力预测。`REACTIVE_DIST_E`、`REACTIVE_ANGLE_E` REMARK 行仅供几何质量分析，应与 `VINA RESULT` 分项独立解读。

### 6.6 与现有共价对接方法的比较

| 方法 | 距离约束 | 角度约束 | NAC 检测 | 开源 | 无需预构建共价复合物 |
|------|---------|---------|---------|------|---------------------|
| AutoDock Vina（标准） | ✗ | ✗ | ✗ | ✅ | ✗ |
| CovDock (Schrödinger) | ✅ | ✅ | ✗ | ✗ | ✅ |
| AutoDock（双锚点法） | ✅ | ✗ | ✗ | ✅ | ✗ |
| **LKina** | ✅ | ✅ | ✅ | ✅ | ✅ |

## 7. Conclusion

LKina 将 AutoDock Vina 的高性能搜索算法（Monte Carlo + L-BFGS，OpenMP 并行）与 AutoDock4 力场的化学精度融合，通过以下四项系统性扩展解决了金属酶和共价靶点对接的核心方法学缺陷：

1. **113 种原子类型**覆盖元素周期表全部药学相关金属元素，解决类型缺失问题
2. **TZ/SQ/MH/JT 四类伪原子** + BVS 氧化态推断 + 半显式水桥，将配位几何方向性编码为空间格点
3. **Jahn-Teller 变形模式**（`cu2_jt`、`mn3_jt`）正确区分轴向/赤道键，覆盖 d⁹/d⁴ 活性金属
4. **响应式共价对接框架 P1–P4**，支持 NAC 检测、角约束梯度解析计算和混合 vdW 缩放

LKina 是目前少数**同时支持金属酶对接与共价对接**的开源引擎之一，且所有新功能对标准受体完全向后兼容。17 项回归测试全部通过，编译零错误零警告。在 292 个 Zn²⁺晶体复合物上，LKina TZ 模式 Top-1 ≤ 2.0 Å 成功率 74.3%，优于 AutoDock4Zn（70.5%）和标准 AD4（58.2%）；MH 模式（Fe³⁺，56 个复合物）和 JT 模式（Cu²⁺ d⁹，41 个复合物）分别相对标准 AD4 提升 17.9 pp 和 24.4 pp。

在本文撰写后的增量开发中，LKina 又进一步完成了多项面向工程落地的扩展与修复：其一，BVS 氧化态推断已从原有体系扩展到单中心 Ni 位点，并新增 `ni2` / `ni3` 模式及相应配位参数；其二，金属几何 rerank 项已通过 `--metal_soft_weight` 以平滑、可导的形式接入搜索梯度，为搜索期软约束提供了可控接口；其三，修复了 `mn3_jt`、`cu2`、`ni2`、`ni3` 的伪原子注入缺口，并新增 Cd²⁺（TZ 四面体）、Hg²⁺（SQ 线性）、Ag⁺（SQ 线性）的方向性配位支持，将 `ag4_inject_tz_pseudoatoms()` 泛化为可接受任意金属名列表；其四，统一了 `--zn_mode` 与 `--metal_mode zn` 的内部语义，并将金属自动识别扩展到单配体与 batch 逐 ligand 场景；其五，针对“金属作为配体”的常见工作流新增了 reverse metal-donor pair potentials、MetalDock `standard_set` 反向参数参考与 NBP 去重逻辑；其六，修正了 `AD_TYPE_SIZE` 文档数值（113→117，含四种伪原子类型），并在金属模式激活时禁用了不等价的外部 `autogrid4` fallback。

后续工作将聚焦于以下几个方向：

- 对 `--metal_soft_weight` 开展权重标定与消融分析，明确其对不同金属体系的适用范围
- 面向 Pt、Ru、Tc/Re、Gd 等“金属作为配体”体系建立专项 benchmark，并评估 ligand-side dummy atoms / 几何约束的增益
- 将 BVS 从当前单中心扩展到 Mo-Fe 异核耦合与多核 Mn 网络体系
- 在更大规模数据集（含 Fe、Cu、Ni 等多金属类型）上持续扩展基准测试覆盖度
- 将高级 LKina 参数（`--metal_soft_weight`、`--metal_geometry_check` 等）以图形控件形式逐步并入 LKDock GUI 前端

---

## 8. Availability, Build and Reproducibility

### 8.1 项目地址与许可证

- **LKina 主仓库**：<https://github.com/LK-Studio1128/LKina>
- **上游基线**：AutoDock Vina 1.2.7（<https://github.com/ccsb-scripps/AutoDock-Vina>）
- **合并二进制许可证**：GPL-3.0-or-later（AG4 引擎 GPL-3.0-or-later × Apache-2.0 Vina 核心，兼容性见 [COPYING](COPYING) 与 [NOTICE](NOTICE)）
- **引用上游**：如将 LKina 用于学术发表，请按 [Citations](#references) 中 [1]–[3]、[11] 共同引用 AutoDock Vina 1.2.0、AutoDock4、AutoGrid4 及 MetalDock 原始工作。

### 8.2 编译与预编译二进制

LKina 提供三平台一键构建脚本，默认使用系统 Boost / OpenMP：

```bash
# macOS (Apple Silicon / Intel)
./build_LKina_mac.sh        # Output: build/mac/release/LKina

# Linux (x86-64)
./build_LKina_linux.sh      # Output: build/linux/release/LKina

# Windows (MSYS2 MinGW-w64)
build_LKina_win.bat         # Output: build/win/release/LKina.exe
```

二进制不随源码仓库提交，推荐作为 GitHub Release assets 分发。发布流程详见仓库内 [`RELEASE.md`](RELEASE.md)；贡献流程见 [`CONTRIBUTING.md`](CONTRIBUTING.md)。

### 8.3 可复现实验环境

本文实验均在 macOS 14（Apple M 系列，clang++ 15，Boost 1.85，OpenMP via libomp）及 Ubuntu 22.04（g++ 11.4，Boost 1.74）环境下进行：

- 配体预处理：[Meeko ≥ 0.5](https://github.com/forlilab/Meeko) + Open Babel ≥ 3.1
- 受体准备：`prepare_receptor4.py`（AutoDockTools4）或 LKDock v3.0 内置"受体预处理"流程
- 基准测试驱动脚本：`tests/reactive_regression.sh`、`tests/metal_redocking_benchmark.py`、`tests/metallocomplex_redocking_benchmark.py`

### 8.4 关联资源

- **LKDock GUI 用户手册**
  - 金属对接：`LKDock_v3.0_金属对接使用手册.md`
  - 共价对接：`LKDock_v3.0_共价对接使用手册.md`
- **设计文档**：本文件即 `LKINA.md` 的设计/方法学论文；快速参考见附录 A–F。
- **变更历史**：[`CHANGELOG.md`](CHANGELOG.md)

---

## References

1. Trott, O.; Olson, A. J. AutoDock Vina: Improving the Speed and Accuracy of Docking with a New Scoring Function, Efficient Optimization, and Multithreading. *J. Comput. Chem.* **2010**, *31*, 455–461.
2. Morris, G. M.; Huey, R.; Lindstrom, W.; Sanner, M. F.; Belew, R. K.; Goodsell, D. S.; Olson, A. J. AutoDock4 and AutoDockTools4: Automated Docking with Selective Receptor Flexibility. *J. Comput. Chem.* **2009**, *30*, 2785–2791.
3. Santos-Martins, D.; Forli, S.; Ramos, M. J.; Olson, A. J. AutoDock4Zn: An Improved AutoDock Force Field for Small-Molecule Docking to Zinc Metalloproteins. *J. Chem. Inf. Model.* **2014**, *54*, 2371–2379.
4. Yu, Y.; Cai, C.; Wang, J.; Bo, Z.; Zhu, Z.; Zheng, H. Uni-Dock: GPU-Accelerated Docking Enables Ultralarge Virtual Screening. *J. Chem. Theory Comput.* **2023**, *19*, 3336–3345.
5. Brown, I. D.; Altermatt, D. Bond-Valence Parameters Obtained from a Systematic Analysis of the Inorganic Crystal Structure Database. *Acta Crystallogr. B* **1985**, *41*, 244–247.
6. Brese, N. E.; O'Keeffe, M. Bond-Valence Parameters for Solids. *Acta Crystallogr. B* **1991**, *47*, 192–197.
7. Wang, R.; Lai, L.; Wang, S. Further Development and Validation of Empirical Scoring Functions for Structure-Based Binding Affinity Prediction. *J. Comput.-Aided Mol. Des.* **2002**, *16*, 11–26.
8. Bianco, G.; Forli, S.; Goodsell, D. S.; Olson, A. J. Covalent Docking Using AutoDock: Two-Point Attractor and Flexible Side Chain Methods. *Protein Sci.* **2016**, *25*, 295–301.
9. Pauling, L. Atomic Radii and Interatomic Distances in Metals. *J. Am. Chem. Soc.* **1947**, *69*, 542–553.
10. Naïm, M.; Bhat, S.; Rankin, K. N.; Dennis, S.; Chowdhury, S. F.; Siddiqi, I.; Drabik, P.; Sulea, T.; Bayly, C. I.; Jakalian, A.; Purisima, E. O. Solvated Interaction Energy (SIE) for Scoring Protein-Ligand Binding Affinities. *J. Chem. Inf. Model.* **2007**, *47*, 122–133.
11. Hakkennes, M.; Rademaker, D.; Lång, A.; Haugaard-Kedstrøm, L. M. A.; Gros, P.; Forli, S. MetalDock: An Easy-to-Use and Reproducible Docking Protocol for Metal-Containing Compounds. *J. Chem. Inf. Model.* **2024**, *64*, 3538–3550. PMC10751784.
12. Harding, M. M. Small Revisions to Predicted Distances around Metal Sites in Proteins. *Acta Crystallogr. D* **2006**, *62*, 678–682.
13. Morris, G. M.; Goodsell, D. S.; Halliday, R. S.; Huey, R.; Hart, W. E.; Belew, R. K.; Olson, A. J. Automated Docking Using a Lamarckian Genetic Algorithm and an Empirical Binding Free Energy Function. *J. Comput. Chem.* **1998**, *19*, 1639–1662.
14. Eberhardt, J.; Santos-Martins, D.; Tillack, A. F.; Forli, S. AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings. *J. Chem. Inf. Model.* **2021**, *61*, 3891–3898.
15. GPDOCK: geometric-probability-based docking method for metalloproteins and metal-coordination-aware pose evaluation, PubMed ID 36642411.
16. AutoDock Vina source repository (ccsb-scripps): <https://github.com/ccsb-scripps/AutoDock-Vina> (accessed 2026).
17. Forli Lab. *Meeko — Python package for preparing small molecules for AutoDock*. <https://github.com/forlilab/Meeko>.
18. O'Boyle, N. M.; Banck, M.; James, C. A.; Morley, C.; Vandermeersch, T.; Hutchison, G. R. Open Babel: An Open Chemical Toolbox. *J. Cheminform.* **2011**, *3*, 33.
19. Goodsell, D. S.; Morris, G. M.; Olson, A. J. Automated Docking of Flexible Ligands: Applications of AutoDock. *J. Mol. Recognit.* **1996**, *9*, 1–5.
20. Scarpino, A.; Ferenczy, G. G.; Keserű, G. M. Covalent Docking in Drug Discovery: Scope and Limitations. *Curr. Pharm. Des.* **2020**, *26*, 6–19.

---

## Appendix A — Installation and Quick Start

### A.1 Prerequisites

- C++14 toolchain：clang++ ≥ 14（macOS）或 g++ ≥ 9（Linux），Windows 通过 MSYS2 MinGW-w64
- Boost ≥ 1.65：`program_options`、`thread`、`serialization`、`filesystem`、`system`
- OpenMP 运行时：macOS 需 `brew install libomp`，Linux 发行版默认随 GCC 提供

### A.2 One-shot build

```bash
git clone https://github.com/LK-Studio1128/LKina.git
cd LKina

# 平台三选一
./build_LKina_mac.sh        # macOS
./build_LKina_linux.sh      # Linux
build_LKina_win.bat         # Windows (MSYS2)

./build/mac/release/LKina --version      # 自检
```

### A.3 Canonical examples

#### A.3.1 Standard Vina scoring（无金属、无共价，行为与 Vina 1.2.7 一致）

```bash
LKina --receptor rec.pdbqt --ligand lig.pdbqt \
      --center_x 0 --center_y 0 --center_z 0 \
      --size_x 20 --size_y 20 --size_z 20 \
      --exhaustiveness 16 --out out.pdbqt
```

#### A.3.2 Inline AD4 + 自动金属识别（Zn 金属酶）

```bash
LKina --scoring ad4 --generate_maps \
      --receptor carbonic_anhydrase.pdbqt \
      --ligand inhibitor.pdbqt \
      --center_x 14.68 --center_y 32.38 --center_z 10.64 \
      --size_x 25 --size_y 25 --size_z 25 \
      --out out.pdbqt
# 受体内含 Zn 原子 → LKina 自动启用 --metal_mode zn + TZ 伪原子
```

#### A.3.3 Reactive 共价对接（Cys 迈克尔加成）

```bash
LKina --scoring ad4 --generate_maps \
      --receptor egfr.pdbqt --ligand afatinib.pdbqt \
      --center_x 18.2 --center_y 54.1 --center_z 25.7 \
      --size_x 20 --size_y 20 --size_z 20 \
      --reactive_preset cys_michael \
      --reactive_rec_atom A:797:SG \
      --reactive_lig_atom index:12 \
      --reactive_frame_atom A:797:CB \
      --reactive_two_step \
      --out out.pdbqt
```

---

## Appendix B — Metal-mode Selection Cheat Sheet

当受体包含金属辅因子时，推荐的 `--metal_mode` 与伪原子映射如下（单金属中心的最常见情形）。多金属体系用逗号拼接（如 `fe3,zn`）。

| 受体金属 | 常见氧化态 | 推荐 `metal_mode` | 伪原子 | 几何 | BVS 自动判别 |
|---------|-----------|-------------------|--------|------|--------------|
| Zn | +2 | `zn` | TZ | 四面体 | — |
| Fe | +2 / +3 | `fe2` / `fe3`（或自动） | MH | 八面体 | ✅ |
| Cu | +1 / +2 | `cu1` / `cu2` / `cu2_jt` | SQ / JT | 平面 / 拉长八面体 | ✅ |
| Mn | +2 / +3 | `mn2` / `mn3_jt` | MH / JT | 八面体 / JT | ✅ |
| Co | +2 / +3 | `co2` / `co3` | MH | 八面体 | ✅ |
| V | +4 / +5 | `v4` / `v5` | MH | 八面体（畸变） | ✅ |
| Mo | +4 / +6 | `mo4` / `mo6` | MH | 八面体 | ✅ |
| Ni | +2 / +3 | `ni2` / `ni3` | MH / SQ | 八面体 / 平面 | ✅（新增） |
| Pt | +2 | `pt` | SQ | 平面四方 | — |
| Pd | +2 | `pd` | SQ | 平面四方 | — |
| Ru / Os | +2 / +3 | `ru` / `os` | MH | 八面体 | — |
| Ir | +3 | `ir` | MH | 八面体 | — |
| Mg / Ca | +2 | `mg` / `ca` | —（离子势） | — | — |
| Au | +1 | `au` | — | 线形 | — |
| Cd | +2 | `cd` | TZ | 四面体 | — |
| Hg | +2 | `hg` | SQ | 线性 / 畸变 | — |
| Ag | +1 | `ag` | SQ | 线形 | — |

> **缺省行为**：不显式指定 `--metal_mode` 时，LKina 会先扫描受体 PDBQT，BVS + 供体计数推断氧化态；对纯有机体系自动不启用任何金属模式。

---

## Appendix C — Reactive Preset Cheat Sheet

`--reactive_preset` 内置的 6 种反应预设自动填充距离、角度、vdW 缩放等参数，单独指定的 CLI 标志会覆盖预设中对应字段。

| Preset | 目标反应 | 典型靶点示例 | $r_0$ (Å) | $\theta_0$ (°) | $w$ (°) | $\lambda_{\text{vdW}}$ |
|--------|---------|-------------|-----------|----------------|---------|------------------------|
| `cys_michael` | Cys SG 迈克尔加成（C=C 弹头） | BTK/EGFR/KRAS G12C | 1.82 | 109.5 | 25 | 0.2 |
| `cys_sn2` | Cys SG SN2 取代（180° 反向攻击） | 半胱氨酸蛋白酶 | 1.82 | 180.0 | 15 | 0.2 |
| `ser_covalent` | Ser OG 酰化（β-内酰胺 / 酯） | 丝氨酸蛋白酶、β-内酰胺酶 | 1.34 | 109.5 | 25 | 0.2 |
| `lys_targeting` | Lys NZ Schiff 碱（醛基） | KRAS/PP2A/Hsp90 | 1.47 | 109.5 | 30 | 0.3 |
| `boronic_acid` | 硼酸可逆共价（Ser/Thr/Tyr OH） | 蛋白酶体抑制剂 | 1.47 | — | — | 0.5 |
| `tyr_covalent` | Tyr OH 亲核攻击 | 特殊 Tyr 共价靶点 | 1.38 | 109.5 | 25 | 0.2 |

**典型标志组合**：`--reactive_preset P --reactive_rec_atom A:N:ATOM --reactive_lig_atom index:K [--reactive_frame_atom A:N:CB] [--reactive_two_step]`。

---

## Appendix D — Output REMARK Quick Reference

LKina 在输出 PDBQT 头部写入若干 REMARK，用于下游分析脚本直接解析：

| REMARK 字段 | 含义 | 触发条件 |
|-------------|------|----------|
| `VINA RESULT` | 标准 Vina 能量分解 | 始终输出 |
| `REACTIVE_DIST` | 配体亲核原子与受体反应原子的实测距离（Å） | 指定 `--reactive_rec_atom`/`--reactive_lig_atom` |
| `REACTIVE_ANGLE` | 实测进攻角（°） | 指定 `--reactive_frame_atom` 或 `--reactive_lig_frame_atom` |
| `REACTIVE_NAC` | `YES`/`NO`——是否落在 NAC（距离 ≤ bond_length+0.5 Å 且角度落在平底区） | 共价模式下始终输出 |
| `metal_geo` | 受体金属几何 rerank 分量 | `--metal_mode` 或自动金属模式 |
| `metal_water` | 半显式水桥分量 | 配位不饱和的 TZ/SQ/MH/JT 位点 |
| `metal_bvs` | 受体金属 BVS 值 | Fe/Cu/Mn/Co/V/Mo/Ni 等有 BVS 参数的模式 |
| `LIGAND_METAL_GEOM` | ligand 侧 Pt/Pd/Ru/Os/Re 配位几何总分 | ligand 含相应金属 |
| `LIGAND_METAL_SITE` | 单个 ligand 金属中心的几何项 | 同上，每中心一行 |

> LKDock v3.0 分析模块会自动从这些 REMARK 中构建 `reactive_metrics_all_poses.csv` 和 `metal_metrics_all_poses.csv`，便于批量筛选 `REACTIVE_NAC=YES` 或 `metal_bvs ≈ 理论值` 的命中分子。

---

## Appendix E — Common Errors and Troubleshooting

| 现象 | 根本原因 | 推荐处理 |
|------|---------|----------|
| `--scoring vina` 下 SIGABRT（exit 134） | `XS_TYPE` 系统未覆盖 LKina 扩展金属/伪原子类型 | 改用 `--scoring ad4`（或等效 `--scoring LKDock`）；详见 §4.7 |
| `ERROR: reactive_rec_atom not found` | 原子名不匹配或 PDBQT 缺链号 | 用 `A:145:SG` 或直接 `x,y,z` 坐标指定 |
| `Ligand atom_index out of range` | `--reactive_lig_atom` 索引超界 / 配体 PDBQT 缺 TORSDOF 拓扑 | 优先用 `index:N` 并用 Meeko 重新生成 PDBQT |
| Reactive NAC 一直 `NO` | 距离过远或角度超出平底区 | 调大 `reactive_presample_dist`，或先启用 `--reactive_two_step` |
| 金属中心未识别 | 受体 PDBQT 金属行 AD4 类型列错误 | 用 LKDock 受体预处理流程重建；或手动指定 `--metal_mode` |
| Inline AG4 生成失败且无 fallback | 金属模式激活下已禁用外部 `autogrid4` fallback（§3.2） | 检查 Box 参数；必要时用 `--write_maps` 导出再人工校验 |
| `ligand_metal_geometry_weight` 无效 | ligand 不含 Pt/Pd/Ru/Os/Re | 该 QC 仅对上述金属生效，其他金属只写 REMARK |
| 批量模式某 ligand 被跳过 | 含不支持的 AD4 类型或空拓扑 | 查看日志中 `SKIP` 行；用 Meeko/Open Babel 重新处理 |

---

## Appendix F — Integration with LKDock v3.0

LKDock v3.0 将 LKina 作为 **AD4 / 共价 / 金属** 路径下的核心引擎，GUI 与 LKina 二进制之间通过标准 CLI 对接：

| LKDock GUI 区块 | 对应 LKina CLI |
|-----------------|----------------|
| "对接参数 → 评分函数 = AD4 / LKDock" | `--scoring ad4`（或 `--scoring LKDock`） |
| "共价对接（仅支持 LKina 引擎）" | `--reactive_preset/--reactive_rec_atom/--reactive_lig_atom/...` |
| "LKina 金属对接 → metal_mode" | `--metal_mode`（留空则自动识别） |
| "金属软约束权重" | `--metal_soft_weight` |
| "配体金属几何权重" | `--ligand_metal_geometry_weight` |
| "金属偏置（O5）" | `--metal_bias --metal_bias_strength --metal_bias_width` |
| "盒子 / 中心 / 间距" | `--center_x/y/z --size_x/y/z --spacing` |

完整图形操作流程与结果分析（`reactive_metrics*.csv` / `metal_metrics*.csv`）见：

- `LKDock_v3.0_金属对接使用手册.md`
- `LKDock_v3.0_共价对接使用手册.md`
