# 脚本快速索引

## 📁 文件清单（20个脚本 + 2个文档）

### 📘 文档
- `README.md` - 总体工作流指南
- `STEP7_DATA_SOURCE_GUIDE.md` - Step 7详细数据源说明

---

## 🔵 Step 1-5: MSD扩散分析（5个脚本）

| 脚本 | 功能 | 输入 | 输出 |
|------|------|------|------|
| `step1_detect_outliers.py` | 异常run检测 | GMX .xvg | large_D_outliers.csv |
| `step2_ensemble_analysis.py` | 集合平均 | .xvg + outliers | ensemble_analysis_results.csv |
| `step3_plot_msd.py` | MSD曲线 | .xvg + results | msd_curves/*.png |
| `step4_calculate_ensemble_D.py` | 计算D值 | .xvg + outliers | ensemble_D_values.csv |
| `step5_analyze_sn_content.py` | Sn含量分析 | ensemble_D_values | sn_content_analysis/* |

**运行顺序**: Step1 → Step2 → (Step3 + Step4) → Step5

---

## 🟢 Step 6: 能量与热容分析（5个脚本）

| 脚本 | 功能 | 关键特点 |
|------|------|----------|
| `step6_energy_analysis_v2.py` | LAMMPS能量分析 | ⚠️ 包含载体（240 Al₂O₃） |
| `step5.9calculate_support_heat_capacity.py` | 载体热容计算 | 需单独Al₂O₃模拟数据 |
| `subtract_support_v2.py` | 扣除载体热容 | Cv_cluster = Cv_total - Cv_support |
| `step6.2analyze_cv_series.py` | Cv系列对比 | Cv-1到Cv-5专项分析 |
| `step6_3_adaptive_regional_heat_capacity.py` | 自适应区域热容 | 自动检测相变区间 |

**核心问题**: 载体热容 Cv_support 未知（~18-21 meV/K估计值）

---

## 🟡 Step 7: Lindemann指数与结构分析（10个脚本）

### Step 7 - Lindemann分析（1个）

| 脚本 | 功能 | 判据 |
|------|------|------|
| `step7_lindemann_analysis.py` | Lindemann指数分析 | δ<0.1固态, δ≥0.1液态 |

### Step 7.4 - 热容与相态（2个）

| 脚本 | 功能 | 依赖 |
|------|------|------|
| `step7_4_multi_system_heat_capacity.py` | 基础数据生成 | LAMMPS能量 + Lindemann |
| `step7_4_2_clustering_analysis.py` | 聚类相态分区 | step7_4_all_systems_data.csv |

**运行顺序**: step7_4_multi_system → step7_4_2_clustering

### Step 7.5/7.6 - 结构演化（4个 + 1个工具）

| 脚本 | 功能 | 数据格式 |
|------|------|----------|
| `v625_data_locator.py` | **工具**: 数据定位 | 自动查找分散目录 |
| `step7-5-unified_multi_temp_v626_analysis.py` | CN/Q6综合分析 | v625/v626 |
| `step7-6-1_temp_side_by_side_comparison.py` | 温度并排对比 | v625/v626 |
| `step7-6-2_individual_system_temp_comparison.py` | 单系统多温度 | v625/v626 |
| `step7-6-3_q6_stats_comparison.py` | Q6统计对比 | v625/v626 |

**数据源**: `coordination_time_series_results_sample_*/`

---

## 🎯 典型使用场景

### 场景1: 完整扩散分析

```bash
python step1_detect_outliers.py
python step2_ensemble_analysis.py
python step3_plot_msd.py
python step4_calculate_ensemble_D.py
python step5_analyze_sn_content.py
```

**时间**: ~1-2小时（9659个文件）

---

### 场景2: 热容分析（含载体扣除）

```bash
# 1. 能量分析（含载体）
python step6_energy_analysis_v2.py

# 2. 扣除载体（估计值）
python subtract_support_v2.py

# 3. 区域热容
python step6_3_adaptive_regional_heat_capacity.py
```

**注意**: 载体热容为估计值，建议单独模拟验证

---

### 场景3: 相态分析（Lindemann + 聚类）

```bash
# 1. Lindemann分析
python step7_lindemann_analysis.py

# 2. 生成Step7.4基础数据
python step7_4_multi_system_heat_capacity.py

# 3. 聚类相态分区
python step7_4_2_clustering_analysis.py --structure all --auto-partition --use-energy
```

---

### 场景4: 结构演化分析

```bash
# 分析所有系列
python step7-5-unified_multi_temp_v626_analysis.py --all

# 温度对比（300K vs 900K）
python step7-6-1_temp_side_by_side_comparison.py --series Pt8Snx --temps 300K,900K
```

---

## 📊 数据流程图

```
【扩散流程】
GMX .xvg → Step1 → Step2 → Step3/4 → Step5
                      ↓
                large_D_outliers.csv
                      ↓
            (可选) Step7.4.2聚类

【能量流程】
LAMMPS能量 → Step6 → subtract_support → Cv_cluster
                ↓
          Step6.2/6.3 (Cv分析)

【结构流程】
Lindemann → Step7 → 熔化温度
能量+Lindemann → Step7.4 → step7_4_all_systems_data.csv → Step7.4.2聚类
v625/v626数据 → Step7.5/7.6 → CN/Q6分析
```

---

## 🔧 脚本大小统计

| 类别 | 脚本数 | 总大小 | 平均大小 |
|------|--------|--------|----------|
| Step 1-5 | 5 | ~111 KB | 22.2 KB |
| Step 6 | 5 | ~191 KB | 38.2 KB |
| Step 7 | 10 | ~376 KB | 37.6 KB |
| **总计** | **20** | **~678 KB** | **33.9 KB** |

**最大脚本**: `step7_4_2_clustering_analysis.py` (137.5 KB)  
**最小脚本**: `step7-6-1_temp_side_by_side_comparison.py` (8.8 KB)

---

## 📌 快速查找

### 按功能查找

- **异常检测**: step1, step7 (--no-filter)
- **集合平均**: step2, step4
- **可视化**: step3, step6, step7.4.2, step7-5, step7-6
- **统计分析**: step5, step6.2, step7-6-3
- **聚类分析**: step7.4.2
- **热容计算**: step6, step6.3, step7.4
- **数据工具**: v625_data_locator, subtract_support_v2

### 按数据源查找

- **GMX .xvg**: step1-5
- **LAMMPS能量**: step6, step7.4
- **Lindemann指数**: step7, step7.4
- **配位数/Q6**: step7-5, step7-6
- **混合数据**: step7.4.2 (能量+Lindemann+D值)

---

## ⚠️ 重要提示

1. **载体热容问题**: Step6输出的热容包含载体（240个Al₂O₃），需扣除
2. **路径硬编码**: 所有脚本中的数据路径需根据实际情况修改
3. **Step 7.4依赖**: step7_4_2必须先运行step7_4_multi_system生成基础数据
4. **v625/v626兼容**: step7-5/7-6会自动检测数据格式

---

**最后更新**: 2025-11-06  
**脚本总数**: 20个Python脚本 + 2个Markdown文档
