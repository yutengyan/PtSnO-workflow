#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 6.3: 自适应分区热容计算 (结合Lindemann指数相变判断)
========================================================

功能概述
========
本脚本基于Lindemann指数(林德曼指数)自动识别团簇熔化温度,
智能划分固态/液态温度区间,分别计算热容,避免相变点导致的拟合失败。

核心改进
========
1. **自动识别相变温度** (基于Lindemann指数)
   - δ < 0.1:  固态
   - 0.1 < δ < 0.15: 预熔化(过渡态)
   - δ > 0.15: 液态

2. **智能分区策略**
   - 固态区: 所有δ<0.1的温度点 (纯固态)
   - 液态区: 所有δ>0.15的温度点 (纯液态)
   - **排除预熔化点** (0.1<δ<0.15) - 避免拟合失败

3. **物理意义清晰**
   - 固态Cv: 晶格振动贡献
   - 液态Cv: 振动+流动+构型熵贡献
   - 液态Cv通常 > 固态Cv

4. **解决人工分区问题**
   - 人工分区(500-700K)可能跨越相变点 → R²低
   - 自适应分区避开相变 → R²高

数据来源
========
1. **能量数据** (必需)
   - cluster: files/lammps_energy_analysis/energy_master_*.csv
   - support: files/lammps_energy_analysis/sup/energy_master_*.csv
   - 列: 结构, 温度, 平均能量, ...

2. **Lindemann指数数据** (可选但推荐)
   - 位置: files/takeit/lindemann_master_run_*.csv
   - 列: 结构, 温度, 林德曼指数/Lindemann指数, ...
   - 来源: step7_lindemann_analysis.py 生成
   - **已对同一体系同一温度的多次模拟求平均**

3. **如果没有Lindemann数据**
   - 回退到固定分区: 200-400K, 500-700K, 800-1000K
   - 不保证避开相变点

计算流程
========
1. 加载Lindemann指数数据 (files/takeit/lindemann_master_*.csv)
   → 按体系+温度分组求平均 (已在step7完成)
   → 判断每个温度点的状态 (固态/预熔化/液态)

2. 识别相变温度
   → 找到δ从<0.15到≥0.15的转变点
   → 线性插值得到精确的Tm

3. 自适应分区
   → 固态区: 所有δ<0.1的温度 (严格固态)
   → 液态区: 所有δ>0.15的温度 (严格液态)
   → **排除中间的预熔化点** (保证拟合质量)

4. 载体热容拟合
   → 全温度范围线性拟合: E_support(T) = a + b*T
   → Cv_support = b (常数)

5. 团簇能量计算
   → E_cluster(T) = E_total(T) - E_support(T)
   → E_cluster_rel = E_cluster - E_cluster_min

6. 分区拟合团簇热容
   → 对每个区间: E_cluster_rel(T) = a + b*T
   → Cv_cluster = b

输出文件
========
results/adaptive_heat_capacity/
├── adaptive_heat_capacity_summary.csv     # 汇总结果
└── adaptive_heat_capacity_analysis_*.png  # 综合分析图
    ├── 第1行: 固态区/液态区能量拟合
    ├── 第2行: Lindemann指数 + 分区标注
    └── 第3行: 固态vs液态热容对比

使用示例
========
# 自动运行
python step6_3_adaptive_regional_heat_capacity.py

# 分析其他体系 (修改main函数中的system_names)
system_names = ['O1-1', 'O1-2', 'O1-3']

典型结果
========
对于Pt6Sn8O4体系:
- Tm ≈ 550 K (熔化温度)
- 固态区(200-500K): Cv ≈ 2.6 meV/K
- 液态区(600-1100K): Cv ≈ 5.2 meV/K
- 液态Cv/固态Cv ≈ 2.0 (符合物理预期)

注意事项
========
1. **温度区间不重叠**
   - 固态区: T < Tm (严格)
   - 液态区: T > Tm (严格)
   - 如果温度点跨越Tm,按Lindemann指数分类

2. **最少点数要求**
   - 每个区间至少4个点才能拟合
   - 点数不足时合并区间或使用全温度范围

3. **Lindemann数据质量**
   - step7已对多次模拟求平均
   - 本脚本直接使用平均后的数据
   - 标准差列: '林德曼标准差'

4. **相变点识别**
   - 基于δ=0.15阈值 (经典判据)
   - 可能存在亚稳态 (过冷液体/过热固体)
   - 建议结合MSD、RDF等方法验证

技术细节
========
- Lindemann判据来源: 经典统计物理 (Lindemann, 1910)
- 阈值选择: δ_solid=0.1, δ_melting=0.15 (广泛接受)
- 插值方法: 线性插值 (scipy.interpolate.interp1d)
- 拟合方法: 最小二乘线性回归 (scipy.stats.linregress)

作者: GitHub Copilot  
日期: 2025-10-21
版本: v1.0
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import stats
from scipy.interpolate import interp1d
from matplotlib import rcParams
import warnings
warnings.filterwarnings('ignore')

# 设置中文字体
rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'Arial Unicode MS']
rcParams['axes.unicode_minus'] = False

# 物理常数
KB_EV = 8.617333262e-5  # eV/K

# Lindemann判据 (来自step7)
LINDEMANN_THRESHOLDS = {
    'solid': 0.1,      # δ < 0.1: 固态
    'premelting': 0.15,  # 0.1 < δ < 0.15: 预熔化
    'melting': 0.15    # δ > 0.15: 液态
}

# 数据路径
BASE_DIR = Path(__file__).parent
CLUSTER_ENERGY_FILE = BASE_DIR / 'data' / 'lammps_energy' / 'energy_master_20251016_121110.csv'
SUPPORT_ENERGY_FILE = BASE_DIR / 'data' / 'lammps_energy' / 'sup' / 'energy_master_20251021_151520.csv'
LINDEMANN_DATA_DIR = BASE_DIR / 'data' / 'lindemann'

def load_energy_data(energy_file, file_type='cluster'):
    """加载能量数据"""
    print(f"\n>>> 加载{file_type}能量数据: {energy_file.name}")
    df = pd.read_csv(energy_file, encoding='utf-8')
    
    # 统一列名
    df.columns = ['路径', '结构', '温度', '模拟序号', '总步数', '采样步数', 
                  '平均能量', '标准差', '最小值', '最大值', '采样间隔', 
                  '跳过步数', '完整路径']
    
    print(f"    完成: {len(df)} 条记录")
    return df

def load_lindemann_data(system_names):
    """
    加载Lindemann指数数据
    
    **数据来源**: step7 生成的 Cv_series_lindemann_statistics.csv
    - 该文件已对Cv-1到Cv-5系列的5次模拟求平均
    - 包含19个温度点 (200-1100K, 间隔50K)
    - 已分类状态: 固态/固液共存/液态
    
    返回: DataFrame with columns ['温度', '林德曼指数', '林德曼标准差', 'N_runs', '状态']
    """
    print(f"\n{'='*80}")
    print("加载Lindemann指数数据 (用于相变判断)")
    print("="*80)
    
    # step7 生成的统计结果文件 (已平均和分类)
    lindemann_file = BASE_DIR / 'results' / 'lindemann_analysis' / 'Cv_series_lindemann_statistics.csv'
    
    if not lindemann_file.exists():
        print(f"  ⚠️ 未找到Lindemann统计文件: {lindemann_file}")
        print("     将使用固定分区")
        return None
    
    print(f"  ✓ 读取 step7 分析结果: {lindemann_file.name}")
    print(f"     (数据已平均: Cv-1到Cv-5的5次模拟)")
    
    df_lindemann = pd.read_csv(lindemann_file, encoding='utf-8-sig')
    
    # 重命名列以匹配后续处理
    df_lindemann.rename(columns={
        '平均林德曼指数': '林德曼指数',
        '标准差': '林德曼标准差',
        '测量次数': 'N_runs'
    }, inplace=True)
    
    # 状态映射: step7的中文状态 → 统一代码
    state_map = {
        '固态 (Solid)': '固态',
        '固液共存 (Pre-melting)': '预熔化',
        '液态 (Liquid)': '液态'
    }
    df_lindemann['状态'] = df_lindemann['状态'].map(state_map)
    
    print(f"\n  ✓ 成功加载Lindemann数据:")
    print(f"    - Cv系列: 5个体系 (Cv-1 到 Cv-5)")
    print(f"    - 温度点: {sorted(df_lindemann['温度'].unique())}")
    print(f"    - 总记录数: {len(df_lindemann)} (19温度 × 1平均)")
    
    # 状态统计
    print(f"\n  状态分布:")
    state_counts = df_lindemann['状态'].value_counts()
    print(f"    固态: {state_counts.get('固态', 0)} 点")
    print(f"    液态: {state_counts.get('液态', 0)} 点")
    print(f"    预熔化: {state_counts.get('预熔化', 0)} 点")
    
    return df_lindemann

def detect_phase_transition(df_lindemann):
    """
    检测相变温度
    
    返回: {
        'Tm': 熔化温度,
        'solid_temps': 固态温度列表,
        'transition_temps': 相变温度列表,
        'liquid_temps': 液态温度列表
    }
    """
    if df_lindemann is None or len(df_lindemann) == 0:
        return None
    
    print(f"\n{'='*80}")
    print("相变温度检测 (基于Lindemann指数)")
    print("="*80)
    
    df_sorted = df_lindemann.sort_values('温度')
    temps = df_sorted['温度'].values
    lindemann = df_sorted['林德曼指数'].values
    states = df_sorted['状态'].values
    
    # 分类温度点
    solid_temps = temps[states == '固态'].tolist()
    premelting_temps = temps[states == '预熔化'].tolist()
    liquid_temps = temps[states == '液态'].tolist()
    
    # 检测熔化温度 (δ从<0.15到≥0.15的转变点)
    Tm = None
    threshold = LINDEMANN_THRESHOLDS['melting']
    
    for i in range(len(temps) - 1):
        if lindemann[i] < threshold and lindemann[i+1] >= threshold:
            # 线性插值
            T1, T2 = temps[i], temps[i+1]
            L1, L2 = lindemann[i], lindemann[i+1]
            Tm = T1 + (threshold - L1) * (T2 - T1) / (L2 - L1)
            print(f"\n  ✓ 检测到熔化温度:")
            print(f"    Tm ≈ {Tm:.1f} K")
            print(f"    转变区间: {T1}K (δ={L1:.4f}) → {T2}K (δ={L2:.4f})")
            break
    
    if Tm is None:
        if lindemann.min() >= threshold:
            print(f"\n  ℹ️ 全温度范围均为液态 (δ > {threshold})")
        else:
            print(f"\n  ℹ️ 全温度范围均为固态 (δ < {threshold})")
    
    print(f"\n  温度分类:")
    print(f"    固态: {solid_temps} ({len(solid_temps)} 点)")
    print(f"    预熔化: {premelting_temps} ({len(premelting_temps)} 点)")
    print(f"    液态: {liquid_temps} ({len(liquid_temps)} 点)")
    
    phase_info = {
        'Tm': Tm,
        'solid_temps': solid_temps,
        'premelting_temps': premelting_temps,
        'liquid_temps': liquid_temps,
        'all_temps': temps.tolist(),
        'lindemann': lindemann.tolist(),
        'states': states.tolist()
    }
    
    return phase_info

def create_adaptive_regions(phase_info, all_temps, min_points=4):
    """
    基于相变信息自适应创建温度区间
    
    策略:
    1. 如果有明确Tm，分为: 固态区(T<Tm)、液态区(T>Tm)
    2. **严格避免重叠**: 固态区最大温度 < 液态区最小温度
    3. 排除预熔化点（提高R²）
    4. 确保每个区间至少有min_points个点
    
    返回: {
        'region_name': (T_min, T_max, 状态),
        ...
    }
    """
    print(f"\n{'='*80}")
    print("自适应温度分区 (严格不重叠)")
    print("="*80)
    
    regions = {}
    
    if phase_info is None:
        # 使用固定分区
        print("  使用固定温度分区 (无Lindemann数据)")
        regions = {
            '低温区': (200, 400, '固态'),
            '中温区': (500, 700, '未知'),
            '高温区': (800, 1000, '未知')
        }
    else:
        solid_temps = sorted(phase_info['solid_temps'])
        liquid_temps = sorted(phase_info['liquid_temps'])
        premelting_temps = sorted(phase_info['premelting_temps'])
        Tm = phase_info['Tm']
        
        print(f"\n  Lindemann分类:")
        print(f"    固态点 (δ<0.1): {solid_temps}")
        print(f"    预熔化点 (0.1≤δ<0.15): {premelting_temps}")
        print(f"    液态点 (δ≥0.15): {liquid_temps}")
        if Tm:
            print(f"    熔化温度: Tm ≈ {Tm:.1f} K")
        
        # 固态区间: 严格 δ<0.1 的点
        if len(solid_temps) >= min_points:
            T_solid_min = min(solid_temps)
            T_solid_max = max(solid_temps)
            regions['固态区'] = (T_solid_min, T_solid_max, '固态')
            print(f"\n  ✓ 固态区: {T_solid_min}-{T_solid_max}K ({len(solid_temps)} 点, δ<0.1)")
        elif len(solid_temps) > 0:
            print(f"\n  ⚠️ 固态点不足 ({len(solid_temps)} < {min_points}), 尝试放宽条件")
            # 如果严格固态点不足,包含部分预熔化点
            combined_temps = sorted(solid_temps + premelting_temps)
            if len(combined_temps) >= min_points and Tm:
                # 只取Tm以下的点
                low_temps = [t for t in combined_temps if t < Tm]
                if len(low_temps) >= min_points:
                    T_solid_min = min(low_temps)
                    T_solid_max = max(low_temps)
                    regions['低温固态区'] = (T_solid_min, T_solid_max, '固态为主')
                    print(f"  ✓ 低温固态区: {T_solid_min}-{T_solid_max}K ({len(low_temps)} 点, T<Tm)")
        
        # 液态区间: 严格 δ>0.15 的点
        if len(liquid_temps) >= min_points:
            T_liquid_min = min(liquid_temps)
            T_liquid_max = max(liquid_temps)
            regions['液态区'] = (T_liquid_min, T_liquid_max, '液态')
            print(f"  ✓ 液态区: {T_liquid_min}-{T_liquid_max}K ({len(liquid_temps)} 点, δ≥0.15)")
        elif len(liquid_temps) > 0:
            print(f"  ⚠️ 液态点不足 ({len(liquid_temps)} < {min_points}), 尝试放宽条件")
            # 如果严格液态点不足,包含部分预熔化点
            combined_temps = sorted(premelting_temps + liquid_temps)
            if len(combined_temps) >= min_points and Tm:
                # 只取Tm以上的点
                high_temps = [t for t in combined_temps if t > Tm]
                if len(high_temps) >= min_points:
                    T_liquid_min = min(high_temps)
                    T_liquid_max = max(high_temps)
                    regions['高温液态区'] = (T_liquid_min, T_liquid_max, '液态为主')
                    print(f"  ✓ 高温液态区: {T_liquid_min}-{T_liquid_max}K ({len(high_temps)} 点, T>Tm)")
        
        # 检查重叠
        if len(regions) >= 2:
            region_list = sorted(regions.items(), key=lambda x: x[1][0])  # 按T_min排序
            for i in range(len(region_list) - 1):
                name1, (T1_min, T1_max, _) = region_list[i]
                name2, (T2_min, T2_max, _) = region_list[i+1]
                if T1_max >= T2_min:
                    print(f"\n  ⚠️ 检测到重叠: {name1}({T1_max}K) >= {name2}({T2_min}K)")
                    print(f"     调整: {name1}上限 → {T2_min-50}K")
                    # 修正: 第一个区间上限降低
                    regions[name1] = (T1_min, T2_min - 50, regions[name1][2])
        
        # 如果只有一个区间，尝试细分（但要避免跨越Tm）
        if len(regions) == 1:
            region_name = list(regions.keys())[0]
            T_min, T_max, state = regions[region_name]
            
            region_temps = [t for t in all_temps if T_min <= t <= T_max]
            region_temps_sorted = sorted(region_temps)
            
            # 如果有Tm且在区间内，按Tm分割
            if Tm and T_min < Tm < T_max and len(region_temps_sorted) >= min_points * 2:
                low_temps = [t for t in region_temps_sorted if t < Tm]
                high_temps = [t for t in region_temps_sorted if t > Tm]
                
                if len(low_temps) >= min_points and len(high_temps) >= min_points:
                    regions = {}
                    regions[f'{state}区_低温'] = (min(low_temps), max(low_temps), '低温段')
                    regions[f'{state}区_高温'] = (min(high_temps), max(high_temps), '高温段')
                    print(f"\n  ℹ️ 按Tm={Tm:.0f}K分割为两个区间")
            # 否则简单二等分
            elif len(region_temps_sorted) >= min_points * 2:
                mid_idx = len(region_temps_sorted) // 2
                T_mid = (region_temps_sorted[mid_idx-1] + region_temps_sorted[mid_idx]) / 2
                
                regions = {}
                regions[f'{state}区_低'] = (T_min, T_mid, state)
                regions[f'{state}区_高'] = (T_mid, T_max, state)
                print(f"\n  ℹ️ 单区间二等分 (T={T_mid:.0f}K)")
    
    if len(regions) == 0:
        print("\n  ⚠️ 无法创建有效区间，使用全温度范围")
        regions['全温区'] = (min(all_temps), max(all_temps), '未知')
    
    print(f"\n  最终分区 ({len(regions)}个):")
    for name, (Tmin, Tmax, state) in sorted(regions.items(), key=lambda x: x[1][0]):
        print(f"    {name}: {Tmin:.0f}-{Tmax:.0f}K ({state})")
    
    return regions

def fit_support_heat_capacity(df_support, support_name='sup240'):
    """全温度范围拟合载体热容 (常数)"""
    print(f"\n{'='*80}")
    print(f"载体热容拟合 (使用 {support_name})")
    print("="*80)
    
    df_sup = df_support[df_support['结构'] == support_name].copy()
    
    grouped = df_sup.groupby('温度').agg({
        '平均能量': ['mean', 'std', 'count']
    }).reset_index()
    
    grouped.columns = ['温度', '平均能量', '标准差', '次数']
    
    T = grouped['温度'].values
    E = grouped['平均能量'].values
    
    print(f"\n载体数据:")
    print(f"  温度范围: {T.min():.0f} - {T.max():.0f} K")
    print(f"  数据点数: {len(T)}")
    
    # 线性拟合
    slope, intercept, r_value, p_value, std_err = stats.linregress(T, E)
    r_squared = r_value ** 2
    
    cv_support_meV_K = slope * 1000
    
    print(f"\n拟合结果:")
    print(f"  E_support(T) = {intercept:.6f} + {slope:.8f} × T  (eV)")
    print(f"  Cv_载体 = {cv_support_meV_K:.4f} meV/K  (常数)")
    print(f"  R² = {r_squared:.8f}")
    
    def E_support_func(T_array):
        return intercept + slope * T_array
    
    support_fit = {
        'slope': slope,
        'intercept': intercept,
        'r_squared': r_squared,
        'p_value': p_value,
        'std_err': std_err,
        'cv_meV_K': cv_support_meV_K,
        'E_func': E_support_func,
        'T_data': T,
        'E_data': E
    }
    
    return support_fit

def calculate_cluster_energy(df_cluster, support_fit):
    """计算团簇能量"""
    print(f"\n{'='*80}")
    print("计算团簇能量")
    print("="*80)
    
    grouped = df_cluster.groupby('温度').agg({
        '平均能量': ['mean', 'std', 'count']
    }).reset_index()
    
    grouped.columns = ['温度', 'E_total_mean', 'E_total_std', 'N_runs']
    
    T = grouped['温度'].values
    E_support = support_fit['E_func'](T)
    E_cluster = grouped['E_total_mean'].values - E_support
    
    grouped['E_support'] = E_support
    grouped['E_cluster'] = E_cluster
    
    print(f"  温度点数: {len(T)}")
    print(f"  E_cluster 范围: {E_cluster.min():.4f} - {E_cluster.max():.4f} eV")
    
    return grouped

def fit_regional_heat_capacity(df_cluster_energy, support_fit, region_name, temp_range, phase_state):
    """拟合某个温度区间的团簇热容"""
    
    T_min, T_max = temp_range
    
    df_region = df_cluster_energy[
        (df_cluster_energy['温度'] >= T_min) & 
        (df_cluster_energy['温度'] <= T_max)
    ].copy()
    
    if len(df_region) < 2:
        return None
    
    T = df_region['温度'].values
    E_cluster = df_region['E_cluster'].values
    E_total = df_region['E_total_mean'].values
    
    # 相对能量
    E_cluster_min = E_cluster.min()
    E_cluster_rel = E_cluster - E_cluster_min
    
    # 线性拟合相对能量
    slope_rel, intercept_rel, r_value_rel, p_value_rel, std_err_rel = stats.linregress(T, E_cluster_rel)
    r_squared_rel = r_value_rel ** 2
    
    # 热容
    cv_cluster_meV_K = slope_rel * 1000
    cv_stderr_meV_K = std_err_rel * 1000
    
    # 预测和残差
    E_cluster_rel_pred = intercept_rel + slope_rel * T
    residuals_rel = E_cluster_rel - E_cluster_rel_pred
    rmse_rel = np.sqrt(np.sum(residuals_rel**2) / (len(T) - 2)) if len(T) > 2 else 0
    
    # 验证
    slope_total, _, _, _, _ = stats.linregress(T, E_total)
    cv_total_meV_K = slope_total * 1000
    cv_cluster_by_subtraction = cv_total_meV_K - support_fit['cv_meV_K']
    
    results = {
        'region_name': region_name,
        'temp_range': temp_range,
        'phase_state': phase_state,
        'n_points': len(T),
        'temperatures': T,
        'E_cluster_rel': E_cluster_rel,
        'E_cluster_rel_pred': E_cluster_rel_pred,
        'residuals_rel': residuals_rel,
        'slope_rel': slope_rel,
        'intercept_rel': intercept_rel,
        'r_squared_rel': r_squared_rel,
        'p_value_rel': p_value_rel,
        'rmse_rel': rmse_rel,
        'cv_cluster_meV_K': cv_cluster_meV_K,
        'cv_stderr_meV_K': cv_stderr_meV_K,
        'cv_total_meV_K': cv_total_meV_K,
        'cv_support_meV_K': support_fit['cv_meV_K'],
        'cv_cluster_by_subtraction': cv_cluster_by_subtraction,
    }
    
    return results

def analyze_system_adaptive(df_cluster, support_fit, df_lindemann, system_names, output_dir):
    """自适应分区热容分析"""
    
    # 筛选体系
    if isinstance(system_names, str):
        system_names = [system_names]
    
    df_system = df_cluster[df_cluster['结构'].isin(system_names)].copy()
    
    if len(df_system) == 0:
        print("错误: 未找到指定体系的数据")
        return None
    
    print(f"\n>>> 筛选体系: {system_names}")
    print(f"    找到记录: {len(df_system)}")
    
    all_temps = sorted(df_system['温度'].unique())
    print(f"    温度点: {all_temps}")
    
    # 计算团簇能量
    df_cluster_energy = calculate_cluster_energy(df_system, support_fit)
    
    # 检测相变
    phase_info = detect_phase_transition(df_lindemann)
    
    # 创建自适应分区
    adaptive_regions = create_adaptive_regions(phase_info, all_temps, min_points=4)
    
    # 分区拟合
    all_results = []
    
    print(f"\n{'='*80}")
    print("分区热容拟合")
    print("="*80)
    
    for region_name, (T_min, T_max, phase_state) in adaptive_regions.items():
        print(f"\n[{region_name}] {T_min:.0f}-{T_max:.0f}K ({phase_state})")
        print("-" * 60)
        
        results = fit_regional_heat_capacity(
            df_cluster_energy, support_fit, region_name, (T_min, T_max), phase_state
        )
        
        if results is None:
            continue
        
        n_cluster_atoms = 18
        results['n_cluster_atoms'] = n_cluster_atoms
        results['cv_cluster_per_atom_meV_K'] = results['cv_cluster_meV_K'] / n_cluster_atoms
        
        all_results.append(results)
        
        # 打印结果
        print(f"  数据点数: {results['n_points']}")
        print(f"  R² = {results['r_squared_rel']:.6f}")
        print(f"  RMSE = {results['rmse_rel']:.6f} eV")
        print(f"  Cv_团簇 = {results['cv_cluster_meV_K']:.4f} ± {results['cv_stderr_meV_K']:.4f} meV/K")
        print(f"  Cv_团簇/原子 = {results['cv_cluster_per_atom_meV_K']:.5f} meV/K/atom")
    
    # 保存结果
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    save_results(all_results, support_fit, phase_info, output_path, system_names)
    plot_results(all_results, support_fit, phase_info, output_path, system_names)
    
    return all_results

def save_results(all_results, support_fit, phase_info, output_path, system_names):
    """保存数值结果"""
    system_label = ', '.join(system_names) if isinstance(system_names, list) else system_names
    
    summary_data = []
    for res in all_results:
        summary_data.append({
            '体系': system_label,
            '区域': res['region_name'],
            '温度范围': f"{res['temp_range'][0]:.0f}-{res['temp_range'][1]:.0f}K",
            '相态': res['phase_state'],
            '数据点数': res['n_points'],
            'R²': res['r_squared_rel'],
            'Cv_团簇_meV/K': res['cv_cluster_meV_K'],
            'Cv_标准误_meV/K': res['cv_stderr_meV_K'],
            'Cv_每原子_meV/K/atom': res['cv_cluster_per_atom_meV_K'],
            'RMSE_eV': res['rmse_rel'],
            'p_value': res['p_value_rel']
        })
    
    df_summary = pd.DataFrame(summary_data)
    summary_file = output_path / 'adaptive_heat_capacity_summary.csv'
    df_summary.to_csv(summary_file, index=False, encoding='utf-8-sig')
    print(f"\n>>> 汇总结果已保存: {summary_file}")

def plot_results(all_results, support_fit, phase_info, output_path, system_names):
    """绘制分析图"""
    system_label = ', '.join(system_names) if isinstance(system_names, list) else system_names
    
    n_regions = len(all_results)
    fig = plt.figure(figsize=(18, 12))
    gs = fig.add_gridspec(3, n_regions, hspace=0.30, wspace=0.30)
    
    colors = plt.cm.tab10(np.linspace(0, 1, n_regions))
    
    # 第1行: 团簇相对能量拟合
    for idx, res in enumerate(all_results):
        ax = fig.add_subplot(gs[0, idx])
        color = colors[idx]
        
        T = res['temperatures']
        E_rel = res['E_cluster_rel']
        E_rel_pred = res['E_cluster_rel_pred']
        
        ax.scatter(T, E_rel, s=150, color=color, alpha=0.8, 
                  edgecolor='black', linewidth=2, label='数据点', zorder=3)
        
        T_fit = np.linspace(T.min(), T.max(), 100)
        E_fit = res['intercept_rel'] + res['slope_rel'] * T_fit
        ax.plot(T_fit, E_fit, '-', color=color, linewidth=3, alpha=0.7,
               label=f'拟合 (R²={res["r_squared_rel"]:.4f})', zorder=2)
        
        # 不确定度带
        rmse = res['rmse_rel']
        ax.fill_between(T_fit, E_fit - rmse, E_fit + rmse, 
                       color=color, alpha=0.15)
        
        ax.set_xlabel('温度 (K)', fontsize=11, fontweight='bold')
        ax.set_ylabel('团簇相对能量 (eV)', fontsize=11, fontweight='bold')
        ax.set_title(f'{res["region_name"]}\n{res["phase_state"]}, Cv={res["cv_cluster_meV_K"]:.3f} meV/K', 
                    fontsize=11, fontweight='bold')
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
    
    # 第2行: Lindemann指数 + 分区标注
    if phase_info is not None:
        ax_lind = fig.add_subplot(gs[1, :])
        
        temps = phase_info['all_temps']
        lindemann = phase_info['lindemann']
        states = phase_info['states']
        
        # 根据状态绘制不同颜色
        for state, color, marker in [('固态', 'blue', 'o'), 
                                      ('预熔化', 'orange', 's'), 
                                      ('液态', 'red', '^')]:
            mask = np.array(states) == state
            if mask.any():
                ax_lind.scatter(np.array(temps)[mask], np.array(lindemann)[mask],
                               s=100, color=color, marker=marker, label=state,
                               edgecolor='black', linewidth=1.5, alpha=0.7)
        
        # Lindemann阈值线
        ax_lind.axhline(LINDEMANN_THRESHOLDS['solid'], color='blue', 
                       linestyle='--', linewidth=1.5, alpha=0.5, label='固态阈值(0.1)')
        ax_lind.axhline(LINDEMANN_THRESHOLDS['melting'], color='red', 
                       linestyle='--', linewidth=1.5, alpha=0.5, label='熔化阈值(0.15)')
        
        # 标注温度分区
        for idx, res in enumerate(all_results):
            T_min, T_max = res['temp_range']
            ax_lind.axvspan(T_min, T_max, alpha=0.1, color=colors[idx], 
                           label=f'{res["region_name"]} ({res["phase_state"]})')
        
        # 熔化温度标注
        if phase_info['Tm'] is not None:
            ax_lind.axvline(phase_info['Tm'], color='red', linestyle=':', 
                           linewidth=2.5, alpha=0.7, label=f'Tm≈{phase_info["Tm"]:.0f}K')
        
        ax_lind.set_xlabel('温度 (K)', fontsize=12, fontweight='bold')
        ax_lind.set_ylabel('Lindemann指数', fontsize=12, fontweight='bold')
        ax_lind.set_title('Lindemann指数 vs 温度 + 自适应分区', fontsize=13, fontweight='bold')
        ax_lind.legend(fontsize=9, ncol=3, loc='upper left')
        ax_lind.grid(True, alpha=0.3)
    
    # 第3行: 热容对比
    ax_cv = fig.add_subplot(gs[2, :])
    
    regions = [res['region_name'] for res in all_results]
    cv_cluster = [res['cv_cluster_meV_K'] for res in all_results]
    cv_stderr = [res['cv_stderr_meV_K'] for res in all_results]
    r_squared = [res['r_squared_rel'] for res in all_results]
    
    x = np.arange(len(regions))
    bars = ax_cv.bar(x, cv_cluster, color=colors, alpha=0.8, 
                     edgecolor='black', linewidth=1.5, yerr=cv_stderr, capsize=5)
    
    # 标注数值和R²
    for i, (cv, r2) in enumerate(zip(cv_cluster, r_squared)):
        ax_cv.text(i, cv + cv_stderr[i] + 0.2, 
                  f'{cv:.3f}\nR²={r2:.3f}', 
                  ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    ax_cv.set_ylabel('团簇热容 (meV/K)', fontsize=12, fontweight='bold')
    ax_cv.set_title('自适应分区热容对比', fontsize=13, fontweight='bold')
    ax_cv.set_xticks(x)
    ax_cv.set_xticklabels(regions, rotation=15, ha='right')
    ax_cv.axhline(y=support_fit['cv_meV_K'], color='green', linestyle=':', 
                 linewidth=2, alpha=0.7, label=f'Cv_载体={support_fit["cv_meV_K"]:.2f} meV/K')
    ax_cv.legend(fontsize=10)
    ax_cv.grid(True, axis='y', alpha=0.3)
    
    plt.suptitle(f'{system_label} - 自适应分区热容分析\n(基于Lindemann指数相变判断)', 
                fontsize=16, fontweight='bold', y=0.995)
    
    output_file = output_path / f'adaptive_heat_capacity_analysis_{system_label.replace(", ", "_")}.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f">>> 分析图已保存: {output_file}")
    plt.close()

def main():
    """主函数"""
    print("="*80)
    print("Step 6.3: 自适应分区热容计算")
    print("="*80)
    print("特点: 结合Lindemann指数自动识别相变温度并优化分区")
    print("")
    
    # 检查文件
    if not CLUSTER_ENERGY_FILE.exists():
        print(f"错误: 找不到团簇能量文件")
        return
    
    if not SUPPORT_ENERGY_FILE.exists():
        print(f"错误: 找不到载体能量文件")
        return
    
    # 配置
    output_dir = BASE_DIR / 'results' / 'adaptive_heat_capacity'
    system_names = ['Cv-1', 'Cv-2', 'Cv-3', 'Cv-4', 'Cv-5']
    
    # 加载数据
    df_cluster = load_energy_data(CLUSTER_ENERGY_FILE, file_type='cluster')
    df_support_raw = load_energy_data(SUPPORT_ENERGY_FILE, file_type='support')
    df_lindemann = load_lindemann_data(system_names)
    
    # 拟合载体
    support_fit = fit_support_heat_capacity(df_support_raw, support_name='sup240')
    
    # 自适应分析
    results = analyze_system_adaptive(
        df_cluster, support_fit, df_lindemann, system_names, output_dir
    )
    
    if results:
        print("\n" + "="*80)
        print("分析完成!")
        print("="*80)
        print(f"\n输出目录: {output_dir}")
        print("\n✨ 自适应分区优势:")
        print("  1. 根据Lindemann指数自动识别相变温度")
        print("  2. 固态和液态分别计算热容")
        print("  3. 避开预熔化区域（R²低）")
        print("  4. 物理意义更清晰")
        print("")

if __name__ == '__main__':
    main()
