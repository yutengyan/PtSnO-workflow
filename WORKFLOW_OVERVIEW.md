# PtSnO 纳米团簇分析工作流 - 完整概览

> **项目**: Pt-Sn-O纳米团簇MD模拟数据分析  
> **目标**: 扩散行为 + 结构演化 + 相态转变 + 热容计算  
> **脚本总数**: 21个Python脚本 + 3个文档  
> **创建日期**: 2025-11-08  
> **工作目录**: `C:\Users\11207\OneDrive\02_Code\work1-PtSnO\workflow\`

---

## 📋 目录

- [工作流架构](#工作流架构)
- [脚本分类与功能](#脚本分类与功能)
- [数据流程图](#数据流程图)
- [典型使用场景](#典型使用场景)
- [快速启动指南](#快速启动指南)
- [常见问题解答](#常见问题解答)

---

## 🏗 工作流架构

### 三大分析体系

```
┌─────────────────────────────────────────────────────────────────────────┐
│                         扩散分析体系 (Step 1-5)                          │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                         │
│  GMX MSD数据 → Step1 异常检测 → Step2 集合平均 → Step4 D值计算        │
│                                                                         │
│                                      ↓                                  │
│                                 Step3 MSD绘图                           │
│                                                                         │
│                                 Step4 D值计算 → Step5 Sn含量分析        │
│                                                                         │
└─────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────┐
│                         能量热容体系 (Step 6)                            │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                         │
│  LAMMPS能量 → Step6 热容计算 → subtract_support 载体扣除 →             │
│                                                                         │
│               → Step6.3 区域分析 + Step6.2 Cv系列分析                   │
│                                                                         │
└─────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────┐
│                         结构相态体系 (Step 7)                            │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                         │
│  Lindemann指数 ──┐                                                      │
│                  │                                                      │
│  LAMMPS能量 ─────┼──→ Step7.4 整合分析 → Step7.4.2 聚类相态分区        │
│                  │           ↓                                          │
│                  └──→ Step7 Lindemann相态判定                           │
│                                                                         │
│  配位数/Q6数据 → Step7.5 综合分析 + Step7.6 Q6对比                      │
│                                                                         │
└─────────────────────────────────────────────────────────────────────────┘
```

**流程说明**:
- **Step 1-5**: GMX MSD → 异常检测 → 集合平均 → D值计算 → Sn含量分析
- **Step 6**: LAMMPS能量 → 热容计算 → 载体扣除 → 区域/系列分析
- **Step 7**: Lindemann + 能量整合 → 相态判定/聚类 | 配位数/Q6 → 结构演化

### 数据源分布

| 分析体系 | 数据源 | 文件数量 | 数据格式 |
|---------|--------|---------|---------|
| Step 1-5 | GMX MSD | 9,659 | .xvg |
| Step 6 | LAMMPS能量 | 95 | .csv |
| Step 7.4 | LAMMPS能量 + Lindemann | 95 + 3,262 | .csv |
| Step 7.5/7.6 | 配位数/Q6时间序列 | ~1,000+ | .csv (v625/v626) |

---

## 📂 脚本分类与功能

### 🔵 扩散分析流程 (Step 1-5) - 5个脚本

#### 1. `step1_detect_outliers.py` - 异常检测

**功能**:
- 从9,659个GMX .xvg文件提取D值
- 使用3种方法检测异常（IQR/3σ/MAD）
- 生成异常run清单

**输入**:
```
d:/OneDrive/py/Cv/lin/MSD_Analysis_Collection/test-unwrap-new/file/
├── collected_gmx_msd/ (5,910个)
└── gmx_msd_results_*/ (3,749个)
```

**输出**:
- `results/large_D_outliers.csv` - 异常run清单
- `results/ensemble_comparison.csv` - 改进前后对比

**运行**:
```bash
python step1_detect_outliers.py
```

---

#### 2. `step2_ensemble_analysis.py` - 集合平均分析

**功能**:
- 舍弃异常run后重新集合平均
- Intercept + D值二次筛选
- 统计分析与质量报告

**输入**:
- GMX .xvg文件
- `results/large_D_outliers.csv`

**输出**:
- `results/ensemble_analysis_filtered.csv`
- `results/ensemble_analysis_results.csv`
- `results/run_quality_report.txt`

**运行**:
```bash
python step2_ensemble_analysis.py
```

**关键参数**:
```python
INTERCEPT_MAX = 20.0  # Intercept阈值 (Ų)
D_MAX_THRESHOLD = 0.1  # D值上限 (cm²/s)
```

---

#### 3. `step3_plot_msd.py` - MSD曲线绘制

**功能**:
- 绘制集合平均MSD曲线
- 叠加单次run（半透明）
- 标注异常run（红色虚线）

**输出**:
- `results/msd_curves/*.png` - 每个系统的MSD曲线

**运行**:
```bash
python step3_plot_msd.py
```

---

#### 4. `step4_calculate_ensemble_D.py` - 扩散系数计算

**功能**:
- 从集合平均MSD拟合D值
- 计算标准误差和置信区间
- 生成完整D值数据库

**输出**:
- `results/ensemble_D_analysis/ensemble_D_values.csv`
- `results/ensemble_D_analysis/D_calculation_report.txt`

**运行**:
```bash
python step4_calculate_ensemble_D.py
```

---

#### 5. `step5_analyze_sn_content.py` - Sn含量影响分析

**功能**:
- D vs Sn含量关系
- D vs 温度（不同Sn含量）
- Arrhenius分析 + 活化能计算

**输出**:
- `results/sn_content_analysis/activation_energies.csv`
- `results/sn_content_analysis/*.png` - 系列分析图

**运行**:
```bash
python step5_analyze_sn_content.py
```

---

### 🟢 能量热容分析 (Step 6) - 5个脚本

#### 6. `step6_energy_analysis_v2.py` - LAMMPS能量分析

**功能**:
- 分析LAMMPS总能量（团簇 + 载体）
- 计算总热容 Cv_total = dE/dT

**⚠️ 重要**:
- 输出包含载体（240个Al₂O₃原子）的贡献
- 载体热容 Cv_support 未知（估计值 ~38.2 meV/K）

**输入**:
- `data/lammps_energy/energy_master_20251016_121110.csv`

**输出**:
- `results/energy_analysis_v2/energy_per_system.csv`
- `results/energy_analysis_v2/heat_capacity_per_system.csv`

**运行**:
```bash
python step6_energy_analysis_v2.py
```

---

#### 7. `step5.9calculate_support_heat_capacity.py` - 载体热容计算

**功能**:
- 计算纯Al₂O₃载体热容
- 需要单独模拟数据

**⚠️ 注意**: 需要载体单独模拟数据（目前使用默认值）

**输出**:
- `results/support_heat_capacity_analysis/support_heat_capacity_results.csv`

**运行**:
```bash
python step5.9calculate_support_heat_capacity.py
```

---

#### 8. `subtract_support_v2.py` - 载体热容扣除工具

**功能**:
- 从总热容中扣除载体贡献
- 得到纯团簇热容 Cv_cluster

**公式**:
```
Cv_cluster = Cv_total - Cv_support
```

**运行**:
```bash
python subtract_support_v2.py
```

---

#### 9. `step6.2analyze_cv_series.py` - Cv系列专项分析

**功能**:
- Cv-1到Cv-5系列对比
- 时间演化分析
- 统计一致性检验

**运行**:
```bash
python step6.2analyze_cv_series.py
```

---

#### 10. `step6_3_adaptive_regional_heat_capacity.py` - 自适应区域热容

**功能**:
- 自动检测相变区间
- 三区域热容计算（固态/预熔化/液态）
- 熔化温度识别

**输出**:
- `results/adaptive_heat_capacity/adaptive_heat_capacity_summary.csv`

**运行**:
```bash
python step6_3_adaptive_regional_heat_capacity.py
```

---

### 🟡 结构相态分析 (Step 7) - 11个脚本

#### 11. `step7_lindemann_analysis.py` - Lindemann指数分析

**功能**:
- 分析Lindemann指数 δ（原子振动幅度）
- 相态判定: δ < 0.1 (固态), δ ≥ 0.1 (液态)
- 4级路径签名精确匹配

**输入**:
- `data/lindemann/lindemann_master_run_20251025_205545.csv`
- `results/large_D_outliers.csv` (可选)

**输出**:
- `results/step7_lindemann/*.png` - δ vs T曲线
- `results/step7_lindemann/lindemann_comparison_*.csv`

**运行**:
```bash
# 使用MSD过滤
python step7_lindemann_analysis.py

# 不过滤
python step7_lindemann_analysis.py --no-filter
```

---

#### 12. `step7_3_individual_runs_heat_capacity.py` - 单次运行热容分析

**功能**:
- 原型单系统分析（Cv系列）
- 三区域线性拟合

**输出**:
- `results/step7_3_individual_runs/step7_3_individual_runs_analysis.png`
- `results/step7_3_individual_runs/step7_3_merged_data.csv`

**运行**:
```bash
python step7_3_individual_runs_heat_capacity.py
```

---

#### 13. `step7_4_multi_system_heat_capacity.py` - 多体系热容分析

**功能**:
- 合并LAMMPS能量 + Lindemann指数
- 生成51个结构的完整数据集
- 三区域热容计算

**输入**:
- `data/lammps_energy/energy_master_20251016_121110.csv`
- `data/lindemann/lindemann_master_run_20251025_205545.csv`
- `results/large_D_outliers.csv` (可选)

**输出**:
- `results/step7_4_multi_system/step7_4_all_systems_data.csv` (2,692条记录)
- `results/step7_4_multi_system/*.png` - 51个结构分析图

**运行**:
```bash
# 使用MSD过滤
python step7_4_multi_system_heat_capacity.py --msd-filter

# 不过滤
python step7_4_multi_system_heat_capacity.py
```

**关键输出**: `step7_4_all_systems_data.csv` - **Step 7.4.2的必需输入**

---

#### 14. `step7_4_2_clustering_analysis.py` - 聚类相态分区

**功能**:
- 自动检测相边界（替代固定0.1/0.15阈值）
- K-means/层次聚类/DBSCAN
- 2D/3D/4D特征空间

**依赖**: 必须先运行 `step7_4_multi_system_heat_capacity.py`

**输入**:
- `results/step7_4_multi_system/step7_4_all_systems_data.csv` (**必需**)
- `results/ensemble_analysis_results.csv` (可选，用于D值)

**运行示例**:
```bash
# 基础2D聚类 (温度+δ)
python step7_4_2_clustering_analysis.py --structure pt6sn8 --n-partitions 3

# 3D聚类 (温度+δ+能量)
python step7_4_2_clustering_analysis.py --structure pt6sn8 --n-partitions 3 --use-energy

# 4D聚类 (温度+δ+能量+D值)
python step7_4_2_clustering_analysis.py --structure pt6sn8 --n-partitions 3 --use-energy --use-d-value

# 自动确定最优分区数
python step7_4_2_clustering_analysis.py --structure pt6sn8 --auto-partition --use-energy

# 批量分析所有结构
python step7_4_2_clustering_analysis.py --structure all --auto-partition --use-energy
```

**输出**:
- `results/step7_4_2_clustering/<structure>/*.png` - 聚类可视化
- `results/step7_4_2_clustering/<structure>/*.csv` - 聚类结果数据

---

#### 15. `v625_data_locator.py` - 数据定位工具

**功能**:
- 自动查找分散的v625/v626数据目录
- 支持多级目录递归搜索
- 统一数据访问接口

**用途**: **Step 7.5/7.6的基础工具类**

**使用**:
```python
from v625_data_locator import V625DataLocator

locator = V625DataLocator(base_path)
run_paths = locator.find_all_runs('Pt8')
```

---

#### 16. `step7-5-unified_multi_temp_v626_analysis.py` - 统一多温度分析

**功能**:
- 配位数 + Q6 + δ 综合分析
- 自动检测v625/v626格式
- 4-8次重复运行平均
- 键类型统计（Pt-Pt, Pt-Sn, Sn-Sn）

**数据源**: `data/coordination/coordination_time_series_results_sample_20251106_214943/`

**运行**:
```bash
# 分析Pt8Snx系列
python step7-5-unified_multi_temp_v626_analysis.py --series Pt8Snx

# 分析PtxSn8-x系列
python step7-5-unified_multi_temp_v626_analysis.py --series PtxSn8-x

# 分析Pt6Snx系列
python step7-5-unified_multi_temp_v626_analysis.py --series Pt6Snx

# 分析所有系列
python step7-5-unified_multi_temp_v626_analysis.py --all

# 启用MSD过滤
python step7-5-unified_multi_temp_v626_analysis.py --series Pt8Snx --enable-msd-filter
```

**输出**:
- `results/step7.5.unified/<series>_multi_temp_data.csv`
- `results/step7.5.unified/<series>_comprehensive_analysis.png`
- `results/step7.5.unified/<series>_heatmap.png`

---

#### 17. `step7-5-cv_pt6sn8o4_analysis.py` - Pt6Sn8O4氧化物体系分析 (**新增**)

**功能**:
- 专门分析Pt6Sn8O4氧化物体系
- 5次Cv模拟（Cv-1到Cv-5）
- 19个温度点（200K-1100K, 间隔50K）
- 完整键类型分析（Pt-Pt, Pt-Sn, Pt-O, Sn-O, O-O）

**体系信息**:
- 组成: Pt₆Sn₈O₄ (18原子)
- 路径: `dp-md/4090-ustc/o68/g-1535-Sn8Pt6O4/`
- 重复次数: 5次 (Cv-1, Cv-2, Cv-3, Cv-4, Cv-5)

**运行**:
```bash
# 基础分析
python step7-5-cv_pt6sn8o4_analysis.py

# 启用MSD过滤
python step7-5-cv_pt6sn8o4_analysis.py --enable-msd-filter
```

**输出**:
- `results/step7.5.cv_pt6sn8o4/pt6sn8o4_all_data.csv` - 完整数据表
- `results/step7.5.cv_pt6sn8o4/pt6sn8o4_comprehensive_analysis.png` - 9×19综合图
- `results/step7.5.cv_pt6sn8o4/pt6sn8o4_temperature_trends.png` - 温度趋势图
- `results/step7.5.cv_pt6sn8o4/pt6sn8o4_cv_comparison.png` - Cv运行对比图

**特色分析**:
- **氧键统计**: O-Pt键, O-Sn键分析
- **多运行平均**: 自动平均5次Cv模拟
- **温度依赖**: 200K-1100K完整温度扫描

---

#### 18. `step7-6-1_temp_side_by_side_comparison.py` - 温度并排对比

**功能**:
- 并排展示多个温度的Q6时间演化
- 适合温度效应对比（如300K vs 900K）

**运行**:
```powershell
# PowerShell (Windows) - 完整命令在一行
python step7-6-1_temp_side_by_side_comparison.py --series Pt8Snx --temps "300K,900K" --systems "pt8sn1-2-best,pt8sn2-1-best,pt8sn3-1-best,pt8sn4-1-best,pt8sn5-1-best,pt8sn6-1-best,pt8sn7-1-best,pt8sn8-1-best,pt8sn9-1-best,pt8sn10-2-best"
```

**输出**:
- `results/step7.6_q6_time/q6_time_comparison_300K_vs_900K_Pt8Snx_all.png`

---

#### 19. `step7-6-2_individual_system_temp_comparison.py` - 单系统多温度对比

**功能**:
- 为每个体系单独生成温度对比图
- 3行（cluster_q6, pt_q6, sn_q6）× 2列（300K, 900K）
- 统计信息盒子（均值、标准差、变异系数CV）

**运行**:
```powershell
# PowerShell (Windows)
python step7-6-2_individual_system_temp_comparison.py --series Pt8Snx --temps "300K,900K" --systems "pt8sn5-1-best,pt8sn6-1-best"
```

**输出**:
- `results/step7.6_individual_system/q6_comparison_<system>_300K_vs_900K.png`

---

#### 20. `step7-6-3_q6_stats_comparison.py` - Q6统计对比

**功能**:
- 对比不同体系在300K vs 900K的Q6均值和CV
- 柱状图 + 散点图
- 稳定性评分

**运行**:
```powershell
# PowerShell (Windows)
python step7-6-3_q6_stats_comparison.py --series Pt8 --temps "300K,900K" --systems "pt8sn1-2-best,pt8sn2-1-best,pt8sn3-1-best,pt8sn4-1-best,pt8sn5-1-best,pt8sn6-1-best,pt8sn7-1-best,pt8sn8-1-best,pt8sn9-1-best,pt8sn10-2-best"
```

**输出**:
- `results/step7.6_q6_stats/q6_stats_comparison_300K_vs_900K.png`
- `results/step7.6_q6_stats/q6_stats_comparison_300K_vs_900K.csv`

---

#### 21. `step7_3_individual_runs_heat_capacity.py` - 原型热容分析

**功能**:
- 单系统热容分析原型（Cv系列）
- 三区域线性拟合

**输出**:
- `results/step7_3_individual_runs/step7_3_individual_runs_analysis.png`

**运行**:
```bash
python step7_3_individual_runs_heat_capacity.py
```

---

## 🔄 数据流程图

### 完整流程概览

```
┌─────────────────────────────────────────────────────────────────┐
│                     原始数据 (3个独立来源)                      │
├──────────────────┬──────────────────┬───────────────────────────┤
│  GMX MSD (.xvg)  │ LAMMPS能量 (.csv)│ 配位数/Q6时间序列 (.csv) │
│   9,659 files    │     95 files     │      ~1,000+ files        │
└────────┬─────────┴─────────┬────────┴──────────────┬────────────┘
         │                   │                       │
         ▼                   ▼                       ▼
    ┌─────────┐         ┌─────────┐           ┌──────────────┐
    │ Step1-5 │         │ Step 6  │           │  Step 7.5/7.6│
    │ 扩散分析 │         │ 热容分析 │           │  结构演化   │
    └────┬────┘         └────┬────┘           └──────┬───────┘
         │                   │                        │
         │                   ▼                        │
         │         ┌──────────────────┐              │
         │         │    Step 7.4     │              │
         └────────►│ 能量+Lindemann  │◄─────────────┘
                   │   整合分析      │
                   └────────┬─────────┘
                            │
                            ▼
                   ┌──────────────────┐
                   │   Step 7.4.2    │
                   │   聚类相态分区  │
                   └──────────────────┘
```

### 详细依赖关系

```
Step 1 (异常检测)
    ↓ large_D_outliers.csv
Step 2 (集合平均) ← large_D_outliers.csv
    ↓ ensemble_analysis_results.csv
    ├─► Step 3 (MSD绘图)
    ├─► Step 4 (D值计算)
    │       ↓ ensemble_D_values.csv
    │   Step 5 (Sn含量分析)
    │
    └─► Step 7 Lindemann (可选过滤)
    └─► Step 7.4 (可选过滤)
    └─► Step 7.4.2 (可选D值特征)
    └─► Step 7.5 (可选过滤)

Step 6 (能量分析)
    ↓ energy_per_system.csv, heat_capacity_per_system.csv
Step 5.9 (载体热容) → subtract_support_v2 → Cv_cluster
    ↓
Step 6.2 (Cv系列)
Step 6.3 (区域热容)

Step 7.4 (多体系热容)
    ↓ step7_4_all_systems_data.csv (必需)
Step 7.4.2 (聚类分析)

v625_data_locator (工具类)
    ↓
Step 7.5 (统一多温度分析)
Step 7.5-cv (Pt6Sn8O4专项)
Step 7.6.1/7.6.2/7.6.3 (Q6对比)
```

---

## 🎯 典型使用场景

### 场景1: 完整扩散分析流程

**目标**: 从原始MSD数据到Sn含量效应分析

```bash
# Step 1: 检测异常run
python step1_detect_outliers.py

# Step 2: 集合平均（去除异常）
python step2_ensemble_analysis.py

# Step 3: 绘制MSD曲线（可视化）
python step3_plot_msd.py

# Step 4: 计算D值
python step4_calculate_ensemble_D.py

# Step 5: Sn含量影响分析
python step5_analyze_sn_content.py
```

**预计时间**: ~1-2小时（9,659个文件）

**输出**: 
- D vs T曲线
- Arrhenius图
- 活化能数据

---

### 场景2: 热容完整分析（含载体扣除）

**目标**: LAMMPS能量 → 总热容 → 团簇热容

```bash
# Step 1: 能量分析（含载体）
python step6_energy_analysis_v2.py

# Step 2: 载体热容计算（可选，需单独数据）
python step5.9calculate_support_heat_capacity.py

# Step 3: 扣除载体贡献
python subtract_support_v2.py

# Step 4: 区域热容分析
python step6_3_adaptive_regional_heat_capacity.py

# Step 5: Cv系列对比
python step6.2analyze_cv_series.py
```

**注意**: 载体热容为估计值（~38.2 meV/K），建议单独模拟验证

---

### 场景3: 相态分析（Lindemann + 聚类）

**目标**: 自动相态分区 + 相边界识别

```bash
# Step 1: Lindemann指数分析
python step7_lindemann_analysis.py

# Step 2: 生成Step7.4基础数据（必需）
python step7_4_multi_system_heat_capacity.py --msd-filter

# Step 3: 聚类相态分区
python step7_4_2_clustering_analysis.py --structure pt6sn8 --auto-partition --use-energy

# Step 4: 批量分析所有结构
python step7_4_2_clustering_analysis.py --structure all --auto-partition --use-energy
```

**输出**:
- δ vs T相图
- 自动聚类结果
- 相边界温度

---

### 场景4: 结构演化完整分析

**目标**: 配位数 + Q6 + 键类型统计

```powershell
# Step 1: 综合分析（所有系列）
python step7-5-unified_multi_temp_v626_analysis.py --all --enable-msd-filter

# Step 2: Pt6Sn8O4氧化物专项分析
python step7-5-cv_pt6sn8o4_analysis.py --enable-msd-filter

# Step 3: 温度对比（300K vs 900K）
python step7-6-1_temp_side_by_side_comparison.py --series Pt8Snx --temps "300K,900K" --systems "pt8sn1-2-best,pt8sn2-1-best,pt8sn3-1-best,pt8sn4-1-best,pt8sn5-1-best"

# Step 4: 单系统详细对比
python step7-6-2_individual_system_temp_comparison.py --series Pt8Snx --temps "300K,900K" --systems "pt8sn5-1-best"

# Step 5: Q6统计对比
python step7-6-3_q6_stats_comparison.py --series Pt8 --temps "300K,900K" --systems "pt8sn1-2-best,pt8sn2-1-best,pt8sn3-1-best"
```

**输出**:
- 配位数演化曲线
- Q6时间序列
- 键类型统计热图
- 温度依赖性分析

---

### 场景5: 快速诊断单个体系

**目标**: pt6sn8体系的扩散+结构+相态

```bash
# 1. 检查D值
python step4_calculate_ensemble_D.py  # 查看pt6sn8的D值数据

# 2. 聚类相态分区
python step7_4_2_clustering_analysis.py --structure pt6sn8 --auto-partition --use-energy

# 3. 结构演化分析
python step7-5-unified_multi_temp_v626_analysis.py --series Pt6Snx --enable-msd-filter
```

**预计时间**: ~10-15分钟

---

## 🚀 快速启动指南

### 环境配置

#### 1. Python环境

```bash
# 推荐使用conda
conda activate base  # 或您的环境名

# 检查Python版本
python --version  # 建议 ≥ 3.8
```

#### 2. 安装依赖

```bash
# 基础包
pip install pandas numpy scipy matplotlib

# Step 7.4.2额外依赖
pip install scikit-learn seaborn

# 可选（更好的进度条）
pip install tqdm
```

#### 3. 验证安装

```python
# 测试导入
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
print("✅ 所有依赖已安装")
```

---

### 目录结构检查

确保以下目录结构存在：

```
workflow/
├── data/
│   ├── coordination/
│   │   └── coordination_time_series_results_sample_20251106_214943/
│   ├── lammps_energy/
│   │   └── energy_master_20251016_121110.csv
│   └── lindemann/
│       └── lindemann_master_run_20251025_205545.csv
├── results/  # 自动创建
└── [所有脚本]
```

**检查命令**:
```powershell
# Windows PowerShell
Write-Host "=== 检查输入文件 ===" -ForegroundColor Green
Test-Path "data\lammps_energy\energy_master_20251016_121110.csv"
Test-Path "data\lindemann\lindemann_master_run_*.csv"
Test-Path "data\coordination\coordination_time_series_results_sample_20251106_214943"
```

---

### 首次运行推荐流程

#### 新手入门（3步）

```bash
# 1. 测试Step 7.4（最快，数据少）
python step7_4_multi_system_heat_capacity.py

# 2. 查看输出
# results/step7_4_multi_system/step7_4_all_systems_data.csv

# 3. 聚类分析（单个结构）
python step7_4_2_clustering_analysis.py --structure pt6sn8 --auto-partition
```

**预计时间**: ~3-5分钟

---

#### 完整体验（6步）

```bash
# 1. 扩散分析（数据多，时间长）
python step1_detect_outliers.py
python step2_ensemble_analysis.py

# 2. 热容分析
python step6_energy_analysis_v2.py

# 3. Lindemann分析
python step7_lindemann_analysis.py

# 4. Step7.4整合
python step7_4_multi_system_heat_capacity.py --msd-filter

# 5. 聚类分析
python step7_4_2_clustering_analysis.py --structure all --auto-partition --use-energy

# 6. 结构演化分析
python step7-5-unified_multi_temp_v626_analysis.py --series Pt8Snx --enable-msd-filter
```

**预计时间**: ~2-3小时

---

## ❓ 常见问题解答

### Q1: Step 7.4.2 提示找不到数据文件？

**错误信息**:
```
FileNotFoundError: results/step7_4_multi_system/step7_4_all_systems_data.csv
```

**原因**: 未运行Step 7.4生成基础数据

**解决方案**:
```bash
# 必须先运行Step 7.4
python step7_4_multi_system_heat_capacity.py

# 然后运行Step 7.4.2
python step7_4_2_clustering_analysis.py --structure pt6sn8
```

---

### Q2: Step 2 如何调整筛选阈值？

**场景**: 需要更严格的异常检测

**解决方案**: 编辑 `step2_ensemble_analysis.py`

```python
# 第30-35行（估计位置）
INTERCEPT_MAX = 20.0      # 改为 15.0（更严格）
D_MAX_THRESHOLD = 0.1     # 改为 0.05（更严格）
MIN_RUNS_REQUIRED = 2     # 改为 3（需要更多run）
```

**影响**: 更多run被筛除，数据质量提升，但可用数据减少

---

### Q3: 如何只分析特定系统？

**场景**: 只想分析pt8开头的系统

**解决方案**: 修改脚本中的 `SYSTEM_FILTER`

```python
# 在脚本开头添加
SYSTEM_FILTER = {
    'include_patterns': [r'^pt8'],  # 只包含pt8开头
    'exclude_patterns': [r'^[Oo]\d+']  # 排除含氧系统
}
```

**适用脚本**: Step 2, 4, 5, 7.4

---

### Q4: v625和v626数据有什么区别？

**区别**:

| 特征 | v625 | v626 |
|------|------|------|
| 格式 | 单次运行 | 多次运行 |
| 目录结构 | `300K/` | `T300.r3.gpu0/`, `T300.r4.gpu0/` |
| 数据量 | 少 | 多（4-8次重复） |
| 统计可靠性 | 低 | 高 |

**自动检测**: Step 7.5/7.6会自动识别格式，优先使用v626

---

### Q5: 如何确认数据已正确加载？

**方法**: 查看脚本输出的统计信息

**正常输出示例**:
```
[V625DataLocator] 初始化完成: coordination_time_series_results_sample_20251106_214943
[Pt8] 找到8个运行文件夹:
  1. 4090-ustc\more\run3\Pt8
  2. dp-md\4090-ustc\GPU-Pt8\Pt8
  ...
  
处理 Pt8Sn1 (pt8sn1-2-best, 8Pt+1Sn=9原子)...
  使用8个运行文件夹
  300K: 5次运行, Q6=0.245±0.012
  ✅ 数据加载成功
```

**异常输出示例**:
```
❌ 错误: 未找到Pt8的运行数据
ValueError: 数据根目录不存在
```

**解决**: 检查数据路径和目录结构

---

### Q6: 载体热容问题如何解决？

**问题**: Step 6输出的热容包含240个Al₂O₃原子

**临时方案**:
```bash
# 使用默认估计值（~38.2 meV/K）
python subtract_support_v2.py
```

**推荐方案**:
1. 单独模拟纯Al₂O₃体系（240原子）
2. 获取 Cv_support 实测值
3. 运行 `step5.9calculate_support_heat_capacity.py`
4. 使用实测值扣除

**公式**:
```
Cv_cluster = Cv_total - Cv_support
```

---

### Q7: 如何并行运行多个脚本？

**场景**: 加速分析流程

**可并行组合**:

```bash
# 终端1: 扩散分析
python step1_detect_outliers.py
python step2_ensemble_analysis.py

# 终端2: 热容分析（独立）
python step6_energy_analysis_v2.py

# 终端3: Lindemann分析（独立）
python step7_lindemann_analysis.py
```

**不可并行**:
- Step 7.4 必须在 Step 7.4.2 之前
- Step 1 必须在 Step 2-5 之前
- Step 2 必须在 Step 4-5 之前

---

### Q8: Windows编码问题如何解决？

**症状**: 中文乱码（锟斤拷）

**原因**: Windows GBK编码 vs UTF-8源文件

**已修复**: 所有脚本已自动处理编码

**手动修复**（如遇新脚本）:
```python
# 在脚本开头添加
import sys
import io
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8')
```

---

### Q9: 如何生成自定义温度对比？

**场景**: 对比200K, 500K, 800K, 1100K

**Step 7.6.2示例**:
```powershell
python step7-6-2_individual_system_temp_comparison.py --series Pt8Snx --temps "200K,500K,800K,1100K" --systems "pt8sn5-1-best"
```

**Step 7.6.3示例**:
```powershell
python step7-6-3_q6_stats_comparison.py --series Pt8 --temps "200K,1100K" --systems "pt8sn1-2-best,pt8sn5-1-best,pt8sn10-2-best"
```

---

### Q10: 如何查看中间数据文件？

**推荐工具**:
- **CSV文件**: Excel, Python pandas
- **图片**: 任意图片查看器
- **报告**: 文本编辑器

**快速查看CSV**:
```python
import pandas as pd

# 查看Step7.4输出
df = pd.read_csv('results/step7_4_multi_system/step7_4_all_systems_data.csv')
print(df.head(10))
print(df.describe())

# 筛选特定结构
df_pt6sn8 = df[df['structure'] == 'pt6sn8']
print(df_pt6sn8)
```

---

## 📚 相关文档

### 主要文档

- **README.md** - 工作流总体指南（基础版）
- **WORKFLOW_OVERVIEW.md** - 本文件（完整详细版）
- **SCRIPT_INDEX.md** - 脚本快速索引
- **STEP7_DATA_SOURCE_GUIDE.md** - Step 7系列详细数据源指南

### 脚本内文档

每个脚本开头都有详细的文档字符串（docstring），包括：
- 功能说明
- 输入输出
- 使用示例
- 注意事项

**查看方式**:
```bash
# 直接打开脚本文件查看前50行
head -n 50 step7_4_multi_system_heat_capacity.py

# 或使用Python
python -c "import step7_4_multi_system_heat_capacity; help(step7_4_multi_system_heat_capacity)"
```

---

## 📊 数据统计概览

### 输入数据规模

| 数据类型 | 文件数 | 总大小 | 平均大小 |
|---------|--------|--------|---------|
| GMX MSD (.xvg) | 9,659 | ~2.1 GB | ~220 KB |
| LAMMPS能量 (.csv) | 95 | ~12 MB | ~126 KB |
| Lindemann指数 (.csv) | 3,262 | ~45 MB | ~14 KB |
| 配位数/Q6 (.csv) | ~1,000+ | ~280 MB | ~280 KB |
| **总计** | **~13,000+** | **~2.4 GB** | - |

### 分析输出规模

| 输出类型 | 数量 | 总大小估计 |
|---------|------|-----------|
| CSV数据表 | ~50 | ~150 MB |
| PNG图片 | ~200+ | ~300 MB |
| TXT报告 | ~10 | ~2 MB |
| **总计** | **~260** | **~450 MB** |

---

## 🔄 版本历史

| 日期 | 版本 | 更新内容 |
|------|------|---------|
| 2025-11-08 | v2.0 | 新增Pt6Sn8O4专项分析 + 完整概览文档 |
| 2025-11-06 | v1.0 | 创建workflow文件夹 + 基础README |
| 2025-10-27 | - | Step 7.6系列脚本开发 |
| 2025-10-26 | - | v626数据格式支持 + v625_data_locator |
| 2025-10-22 | - | Step 7.4.2聚类分析功能 |
| 2025-10-16 | - | Step 1-5基础流程建立 |

---

## 👨‍💻 技术支持与贡献

### 脚本开发

- **主要开发**: GitHub Copilot
- **数据分析**: Pt-Sn-O纳米团簇MD模拟
- **工作流整合**: 2025-11-06 ~ 2025-11-08

### 反馈与改进

如有问题或建议，请：
1. 检查本文档的"常见问题解答"
2. 查看脚本内文档字符串
3. 检查 `results/` 目录的输出报告

---

## 📝 备注与注意事项

### 重要提示

1. **数据路径**: 所有脚本中的硬编码路径已统一为workflow目录结构
2. **并行运行**: Step 1-5, Step 6, Step 7可独立并行运行
3. **数据关联**: Step 7.4.2可选使用Step 2的D值进行聚类
4. **输出目录**: 首次运行会自动创建 `results/` 及子目录
5. **载体热容**: Step 6输出包含载体，需额外处理
6. **编码问题**: Windows环境已自动处理UTF-8编码

### 数据备份建议

```bash
# 定期备份results目录
xcopy /E /I results results_backup_20251108

# 或使用压缩
tar -czf results_20251108.tar.gz results/
```

---

## 🎓 学习路径推荐

### 初学者 (第1-2天)

1. 阅读本文档"快速启动指南"
2. 运行Step 7.4单个结构测试
3. 查看输出图片和CSV文件
4. 理解Step 7.4.2聚类结果

### 进阶使用 (第3-5天)

1. 运行完整Step 1-5扩散分析
2. 学习Step 6热容计算和载体扣除
3. 掌握Step 7.4.2多维聚类
4. 理解v625/v626数据结构

### 高级定制 (第6-7天)

1. 修改筛选阈值和参数
2. 自定义系统过滤规则
3. 扩展新的分析指标
4. 批量分析和结果汇总

---

## 🔗 快速链接

### 关键文件

- **Step 7.4基础数据**: `results/step7_4_multi_system/step7_4_all_systems_data.csv`
- **MSD异常清单**: `results/large_D_outliers.csv`
- **集合平均结果**: `results/ensemble_analysis_results.csv`
- **D值数据库**: `results/ensemble_D_analysis/ensemble_D_values.csv`

### 输出目录

- `results/msd_curves/` - MSD曲线图
- `results/step7_4_2_clustering/` - 聚类分析结果
- `results/step7.5.unified/` - 统一多温度分析
- `results/step7.5.cv_pt6sn8o4/` - Pt6Sn8O4专项分析
- `results/step7.6_q6_time/` - Q6时间序列对比

---

**文档版本**: v2.0  
**最后更新**: 2025-11-08  
**工作目录**: `C:\Users\11207\OneDrive\02_Code\work1-PtSnO\workflow\`  
**脚本总数**: 21个Python脚本 + 3个Markdown文档

---

**祝分析顺利！ 🎉**
