# LKina Zenodo DOI 存档材料包（已核验 · 2026-08-28）

> 目标：为 LKina 仓库（github.com/LK-Studio1128/LKina）的 **v1.0.1** release 在 Zenodo 上取得永久 DOI，回填论文 §6 数据可用性与参考文献 21，使"所有数字可复现"的声明有不可变锚点。

## 一、前置核验结果（已完成）

| 项目 | 状态 |
|---|---|
| 仓库存在且公开 | ✓ github.com/LK-Studio1128/LKina |
| **v1.0.1 release 已发布（2026-08-27）** | ✓ 含 ±1000 钳制、`--no_auto_metal`、三平台二进制（macOS ARM64 / Linux x64 / Win x64），5 个资产 |
| 论文附录 E 声明"修复后 v1.0.1" | ✓ 与 release 完全一致 |
| 三平台二进制 sha256 | ✓ linux df48c330… / mac ea269026… / win baf2f200… |
| v1.0.1 源码包已下载备存 | ✓ `/tmp/LKina_v1.0.1_src.zip`（162 文件，858 KB，commit ee52a51） |
| 论文引用复检（Crossref） | ✓ ref 7 MetalDock 2023;63(24):7816–7825、ref 8 GPDOCK 2023;24(1)、ref 17 Meeko 2025;65(24):13045–13050 全部正确 |

⚠️ 附注：GitHub v1.0.1 **release 说明**中的 MetalDock 引用写成了"2024, 64, 3538–3550"，与 Crossref 实际（2023, 63(24), 7816–7825）不符。论文内引用正确，**release 说明建议顺手更正**（可选）。

## 二、Zenodo 存档元数据（登录后可直接粘贴/核对）

| 字段 | 值 |
|---|---|
| 标题 | LKina v1.0.1: A metal-aware and covalent-reactive molecular docking engine extending AutoDock Vina |
| 版本 | v1.0.1 |
| 作者（Creators） | Luo Xiaowen（罗晓文） |
| 描述 | LKina extends AutoDock Vina 1.2.7 with a 113-type AD4 atom system (80+ metals), an inline AutoGrid 4.2 grid generator (no external autogrid4), four coordination pseudoatoms (TZ/SQ/MH/JT) with BVS oxidation-state inference and Jahn-Teller modes, a semi-explicit water bridge, a four-tier reactive covalent framework (P1–P4) with NAC detection, C3 two-step search and six reaction presets, and metal-as-ligand support. This snapshot (v1.0.1) contains the source code, three-platform release binaries and the full benchmark pipeline. Corresponding manuscript: "LKina: A Metal-Aware and Covalent-Reactive Molecular Docking Engine Extending AutoDock Vina". |
| 关键词 | molecular docking; metalloproteins; AutoDock Vina; metal coordination; covalent docking; near-attack conformation |
| License | GPL-3.0-or-later（合并二进制；Vina 来源文件 Apache-2.0，见 COPYING/NOTICE） |
| 关联文献（Related identifiers） | 预留：论文投稿后回填 DOI |
| 访问权限 | Public |

## 三、执行路径

### 路径 A：GitHub 联动（3 分钟，推荐，无需给我任何凭证）

1. 打开 https://zenodo.org ，点右上角 **Log in**，选择 **"Log in with GitHub"**（首次需授权 Zenodo 访问您的 GitHub 公开仓库）；
2. 登录后进入 **https://zenodo.org/account/settings/github/**；
3. 找到 **LKina** 仓库，把开关拨到 **ON**（Enable）；
4. 现有 release v1.0.1 不会被自动补归档（联动只作用于*之后*的 release）。因此补一个空触发的归档：
   - 方式一：在 GitHub 仓库 **New release** 页面直接给 v1.0.1 打勾的 release 重新点一次 **"Re-use release notes" 后 Publish**（即再发一个 v1.0.2/1.0.1.1 或直接创建 v1.0.1-archival tag）触发快照；
   - 方式二：等论文投稿前若有新修复版（如 v1.0.2）再一并归档（推荐，省一次操作）；
5. Zenodo 自动生成存档页，得到 DOI：形如 `10.5281/zenodo.XXXXXXXX`；
6. 把 DOI 发给我 → 我立即回填 §6 + ref 21 并重建 6 个载体。

### 路径 B：API Token（我可以全程代劳）

1. 登录 Zenodo → https://zenodo.org/account/settings/applications/tokens/ → 新建 **deposit** 权限的 token；
2. 把 token 发给我（仅本会话使用，用完您可立即吊销）；
3. 我会自动完成：创建 deposit → 上传 `/tmp/LKina_v1.0.1_src.zip`（含源码+README+LICENSE/COPYING/NOTICE）→ 填元数据 → publish → 取 DOI → 回填论文 §6 与 ref 21 → 重建 docx/HTML → 提交。

> 注：Zenodo 存档仅含**源码快照**（与 release 的 Source code (zip) 一致，162 文件），论文的 `benchmarks/` 全套数据已随 LKina 仓库 main 分支交付；如需把 benchmarks 数据也一并归档，可在描述中注明并以仓库 zip 整体上传（路径 B 时我可改用完整仓库 zip）。

## 四、拿到 DOI 后的回填计划

1. §6 数据可用性：`出版时将注册 Zenodo DOI [21]` → 改为实际 DOI 链接（https://doi.org/10.5281/zenodo.XXXXXXX）；
2. ref 21 更新为具体存档条目（含标题、版本、年份、DOI）；
3. GitHub README 加 Zenodo DOI 徽章（可选，我可代改）；
4. 重建中英 docx ×2 + HTML ×2，引用守恒校验，提交。

## 五、Checklist（2026-08-29 全部完成 ✅）

- [x] （您）Zenodo 登录 + GitHub 联动开启 LKina 仓库
- [x] （我）CITATION.cff 更新至 v1.0.2 + 创建 v1.0.2-archival release（tag 指向最新 main，含源码+benchmarks+paper/）
- [x] （我）核验 Zenodo 记录：**DOI 10.5281/zenodo.22156943**（concept DOI 10.5281/zenodo.22151067，status published，56.3 MB 仓库 zip，GPL-3.0-or-later）→ https://zenodo.org/records/22156943
- [x] （我）回填 §6 + ref 21（中英 md ×2、docx ×2、HTML ×2），引用守恒校验通过；本地 9bf9102 + GitHub 056414a 已推送
- [ ] （可选，您）更正 GitHub release 说明中 MetalDock 引用年份（2024→2023）
