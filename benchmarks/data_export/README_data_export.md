# LKina Benchmark Data Export

所有论文测试数据的独立可提取版本（UTF-8 BOM，Excel 可直接打开）。
源 JSON 归档于 GitHub 仓库 `benchmarks/`（byi/LKina）。

| 文件 | 行数 | 内容 | 论文对应 |
|---|---|---|---|
| `metal_modes_110.csv` | 110 | 110 种金属模式：模式/令牌/伪原子/d₀/实测距离/误差/E_ideal/E_far/势阱/Vina 状态 | §3.2、表 1、图 3 |
| `4JC_comparison.csv` | 3 | 4JC 三引擎对比：能量/Zn 距离/最近原子/配位报告 | §3.3 |
| `metallocomplex_pool.csv` | 20 | 金属作配体全池：PDB/金属/配体/LKina RMSD/AD4 RMSD/Vina 失败 | §4.6、图 S4 |
| `c3_ablation.csv` | 18 | C3 深口袋消融：6 预设 × 3 变体（single/c3/c3b）距离/角度/NAC/收敛 | §4.6、图 S3 |
| `soft_weight_ablation.csv` | 12 | --metal_soft_weight：3 体系 × 4 权重 供体-金属距离/几何能/总能量 | §4.6 |
| `covalent_presets.csv` | 6 | 六种反应预设：P1/P1+P2/P4 三档/C3 | §3.4、图 6 |
| `feature_family.csv` | 22 | 伪原子几何/BVS 推断/水桥/金属作配体 QC | §3.2.2、图 4-5 |
| `reactive_presets.csv` | 6 | 六预设 P1 vs P1+P2（NAC 判别） | §3.4、图 S2 |
| `redock_summary.csv` | 9 | 104 体系重对接按金属×引擎：成功率/均值 RMSD | §3.6 表 |
| `redock_donor_metal.csv` | 3 | 供体-金属距离汇总：均值/中位/3.0 Å 内计数 | §3.6 表 |

## 关键统计速查

- 110 金属模式：LKina 110/110 对接成功；Vina 1.2.7 仅 34/110（76 解析失败）
- 配位恢复：108/110 误差 < 0.5 Å（均值 |d–d₀| = 0.20 Å）；104/110 存在势阱
- 4JC：LKina Zn–NA 2.13 Å / −34.49；Vina vina −6.39、AD4+maps −13.07（均 2.23 Å OA）
- 重对接 104 体系：供体 3.0 Å 内 90.6%（87/96）vs AD4 27.8%（27/97）vs Vina 52.1%（50/96）
- 向后兼容：1HVR −14.54 / 3PTB −6.202（与 Vina 1.2.7 mode-1 完全一致）
- 六预设：P1 与 P1+P2 全部 rc=0；NAC 判别有效（cys_sn2/boronic_acid → NO）
- 金属作配体全池 n=20：LKina 20/20 完成 + 逐姿态 QC；Vina 20/20 解析失败
- C3 深口袋消融：single 5/6、c3b 3/6、c3 0/6 收敛
- soft_weight：w=0/0.1/0.3/0.5 供体-金属距离漂移 ≤ 0.01 Å（2.12-2.24 Å）

## 复现

所有 CSV 由 `../export_benchmark_data.py` 从 `benchmarks/` 归档 JSON 生成：
```
python3 ../export_benchmark_data.py
```
