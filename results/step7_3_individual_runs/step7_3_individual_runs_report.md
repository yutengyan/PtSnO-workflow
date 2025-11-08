# Step 7.3: 单次模拟级别热容分析报告

Individual Run Level Heat Capacity Analysis Report

================================================================================

**报告生成时间 Report Date**: 2025-11-08 12:05:01

## 1. 分析方法 Analysis Method

### 1.1 核心改进 Core Improvements

相比于 Step 6.3 v2 (使用平均 Lindemann 指数 + 标准差)，本方法的核心改进：

Compared to Step 6.3 v2 (averaged Lindemann index with standard deviation):

- **数据粒度 Data Granularity**: 单次模拟级别 Individual run level (not averaged)
- **分类方法 Classification**: 硬阈值 Hard threshold (δ<0.1=solid, 0.1-0.15=premelting, >0.15=liquid)
- **样本量 Sample Size**: 从 19 个平均点 → **95 个单次模拟点** From 19 averaged → **95 individual runs**
- **优势 Advantage**: 避免标准差导致的相态混淆，保留物理真实性 Avoid phase confusion, preserve physical reality

### 1.2 为什么单次模拟更科学？Why Individual Runs are More Scientific?

**问题示例 Problem Example**: 600K 温度点

- **v2 方法 v2 Method**: δ = 0.098 ± 0.023 → δ_max = 0.121 → 分类为"固态/预熔化过渡" Classified as "solid/premelting transition"
- **实际情况 Reality**: 5次模拟可能是 3次固态(δ<0.1) + 2次预熔化(δ>0.1) Actually 3 solid + 2 premelting runs
- **v2 问题 v2 Problem**: 平均后反而混淆了相态信息！Averaging obscures phase information!

**7.3 方法 7.3 Method**: 每次模拟有明确的 δ 值 → 清晰的相态 → 能量与相态一一对应
Each run has clear δ value → Clear phase → Energy-phase correspondence

## 2. 数据统计 Data Statistics

- **总模拟点数 Total Points**: 95
- **期望点数 Expected**: 95 (5 structures × 19 temperatures)
- **匹配成功率 Match Rate**: 100.0%

### 2.1 相态分布 Phase Distribution

- **液态 Liquid**: 36 个点 points (37.9%)
- **预熔化 Premelting**: 20 个点 points (21.1%)
- **固态 Solid**: 39 个点 points (41.1%)

### 2.2 温度分布 Temperature Distribution

| 温度 Temp (K) | 固态 Solid | 预熔化 Premelting | 液态 Liquid | 总计 Total |
|---------------|-----------|------------------|------------|------------|
| 200 | 5 | 0 | 0 | 5 |
| 250 | 5 | 0 | 0 | 5 |
| 300 | 4 | 1 | 0 | 5 ⚠️ 固液共存 Coexistence |
| 350 | 5 | 0 | 0 | 5 |
| 400 | 2 | 3 | 0 | 5 ⚠️ 固液共存 Coexistence |
| 450 | 5 | 0 | 0 | 5 |
| 500 | 5 | 0 | 0 | 5 |
| 550 | 1 | 3 | 1 | 5 ⚠️ 固液共存 Coexistence |
| 600 | 3 | 2 | 0 | 5 ⚠️ 固液共存 Coexistence |
| 650 | 1 | 3 | 1 | 5 ⚠️ 固液共存 Coexistence |
| 700 | 2 | 2 | 1 | 5 ⚠️ 固液共存 Coexistence |
| 750 | 0 | 1 | 4 | 5 ⚠️ 液化转变 Liquefaction |
| 800 | 1 | 1 | 3 | 5 ⚠️ 固液共存 Coexistence |
| 850 | 0 | 1 | 4 | 5 ⚠️ 液化转变 Liquefaction |
| 900 | 0 | 2 | 3 | 5 ⚠️ 液化转变 Liquefaction |
| 950 | 0 | 0 | 5 | 5 |
| 1000 | 0 | 1 | 4 | 5 ⚠️ 液化转变 Liquefaction |
| 1050 | 0 | 0 | 5 | 5 |
| 1100 | 0 | 0 | 5 | 5 |

## 3. 热容计算结果 Heat Capacity Results

**Support 热容 Support Cv**: 38.2151 meV/K

### 3.1 三段分区热容汇总 Three-Region Heat Capacity Summary

| 区域 Region | 温度范围 Temp Range | 数据点 Points | Cv_cluster (meV/K) | R² | 评级 Grade |
|-------------|---------------------|--------------|-------------------|-----|------------|
| 固态 Solid | 200-800 K | 39 | 1.9140 ± 0.2188 | 0.998901 | 优秀 Excellent ✓✓ |
| 预熔化 Premelting | 300-1000 K | 20 | 3.3424 ± 0.2660 | 0.999263 | 优秀 Excellent ✓✓ |
| 液态 Liquid | 550-1100 K | 36 | 5.1666 ± 0.3082 | 0.998287 | 优秀 Excellent ✓✓ |

### 3.2 详细拟合结果 Detailed Fitting Results

#### 固态 Solid

- **温度范围 Temperature Range**: 200-800 K
- **数据点数 Data Points**: 39
- **线性拟合 Linear Fit**: E = 0.040129 × T + -90799.062
- **Cv_total**: 40.1291 ± 0.2188 meV/K
- **Cv_cluster**: 1.9140 meV/K
- **R²**: 0.998901 ✓✓ (优秀 Excellent)
- **p-value**: 2.35e-56

#### 预熔化 Premelting

- **温度范围 Temperature Range**: 300-1000 K
- **数据点数 Data Points**: 20
- **线性拟合 Linear Fit**: E = 0.041558 × T + -90799.553
- **Cv_total**: 41.5575 ± 0.2660 meV/K
- **Cv_cluster**: 3.3424 meV/K
- **R²**: 0.999263 ✓✓ (优秀 Excellent)
- **p-value**: 1.19e-29

#### 液态 Liquid

- **温度范围 Temperature Range**: 550-1100 K
- **数据点数 Data Points**: 36
- **线性拟合 Linear Fit**: E = 0.043382 × T + -90800.851
- **Cv_total**: 43.3817 ± 0.3082 meV/K
- **Cv_cluster**: 5.1666 meV/K
- **R²**: 0.998287 ✓✓ (优秀 Excellent)
- **p-value**: 1.28e-48

## 4. 与 Step 6.3 v2 的对比 Comparison with Step 6.3 v2

| 特性 Feature | Step 6.3 v2 | Step 7.3 | 改进 Improvement |
|--------------|-------------|----------|------------------|
| 数据粒度 Data Granularity | 平均后 Averaged (19点) | 单次模拟 Individual (95点) | +400% |
| δ 的性质 Delta Property | δ±σ (区间 interval) | δ (确定值 exact) | 明确 Clear |
| 分类依据 Classification | δ±σ vs threshold | δ vs threshold | 简化 Simplified |
| 相态混淆 Phase Confusion | 可能 Possible | 不会 No ✓ | 消除 Eliminated |
| 样本量 Sample Size | 19 | 95 | ×5 |
| R² (固态 Solid) | 0.909 | 0.999 | +10% |
| R² (预熔化 Premelting) | 0.937 | 0.999 | +6.6% |
| R² (液态 Liquid) | 0.979 | 0.998 | +1.9% |
| 物理解释 Physical Meaning | 宏观统计 Macro | 微观样本 Micro ✓ | 更真实 More realistic |


## 5. 数据样本 Data Samples (前20条 First 20)

| 结构 Structure | 温度 Temp | run | δ | 能量 Energy (eV) | 相态 Phase |
|----------------|-----------|-----|---|-----------------|------------|
| Cv-1 | 1000 | 16 | 0.2256 | -90757.670 | 液态 Liquid |
| Cv-1 | 1050 | 17 | 0.3044 | -90755.189 | 液态 Liquid |
| Cv-1 | 1100 | 18 | 0.3428 | -90752.760 | 液态 Liquid |
| Cv-1 | 250 | 1 | 0.0699 | -90788.970 | 固态 Solid |
| Cv-1 | 200 | 0 | 0.0669 | -90791.095 | 固态 Solid |
| Cv-1 | 300 | 2 | 0.1049 | -90786.872 | 预熔化 Premelting |
| Cv-1 | 350 | 3 | 0.0627 | -90785.345 | 固态 Solid |
| Cv-1 | 450 | 5 | 0.0518 | -90781.161 | 固态 Solid |
| Cv-1 | 400 | 4 | 0.1433 | -90782.644 | 预熔化 Premelting |
| Cv-1 | 500 | 6 | 0.0687 | -90779.230 | 固态 Solid |
| Cv-1 | 550 | 7 | 0.1380 | -90776.986 | 预熔化 Premelting |
| Cv-1 | 650 | 9 | 0.0665 | -90772.727 | 固态 Solid |
| Cv-1 | 600 | 8 | 0.0853 | -90774.671 | 固态 Solid |
| Cv-1 | 700 | 10 | 0.0849 | -90771.375 | 固态 Solid |
| Cv-1 | 800 | 12 | 0.1352 | -90765.966 | 预熔化 Premelting |
| Cv-1 | 850 | 13 | 0.2090 | -90763.852 | 液态 Liquid |
| Cv-1 | 750 | 11 | 0.1490 | -90768.630 | 预熔化 Premelting |
| Cv-1 | 900 | 14 | 0.1444 | -90762.104 | 预熔化 Premelting |
| Cv-1 | 950 | 15 | 0.2272 | -90759.761 | 液态 Liquid |
| Cv-2 | 1050 | 17 | 0.2587 | -90755.714 | 液态 Liquid |

## 6. 结论与建议 Conclusions and Recommendations

### 6.1 主要发现 Main Findings

1. **统计精度提升 Statistical Precision**: 样本量从 19 → 95，增加 400% Sample size increased by 400%
2. **拟合质量 Fitting Quality**: 
   - 固态 Solid R² = 0.998901 (0.909 → 0.999)
   - 预熔化 Premelting R² = 0.999263 (0.937 → 0.999)
   - 液态 Liquid R² = 0.998287 (0.979 → 0.998)
3. **相态分类 Phase Classification**: 清晰明确，无混淆 Clear and unambiguous
4. **物理意义 Physical Meaning**: 每次模拟是独立热力学样本 Each run is independent thermodynamic sample

### 6.2 建议 Recommendations

- ✅ **推荐使用 Step 7.3 方法 Recommend Step 7.3**用于最终分析和发表 for final analysis and publication
- ✅ 单次模拟级别保留了完整的物理信息 Individual runs preserve complete physical information
- ✅ 避免了平均后的标准差导致的相态混淆 Avoid phase confusion from averaging
- ✅ 更高的样本量提供了更可靠的统计结果 Higher sample size provides more reliable statistics
- ✅ 所有区域 R² > 0.998，接近完美线性关系 All regions R² > 0.998, nearly perfect linearity


---
**脚本版本 Script Version**: step7_3_individual_runs_heat_capacity.py v1.0
**基于 Based on**: step6_3_adaptive_regional_heat_capacity_3regions_v2.py
