#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 6.1.1: 分区热容拟合图 - 论文出图专用

从 step6_1 提取的核心绘图功能，生成适合论文发表的分区热容拟合图：
- 整体拟合 vs 分区拟合对比
- 按温度平均的数据点（带误差棒）
- 多数投票规则避免温度交叉

输入: step6_1 生成的聚类结果 CSV 文件
输出: 高质量论文图 (PNG/PDF)

用法:
  python step6_1_1_partition_cv_plot.py --structure Pt8sn6
  python step6_1_1_partition_cv_plot.py --structure Air86 --format pdf
  python step6_1_1_partition_cv_plot.py --structure all --dpi 600
  python step6_1_1_partition_cv_plot.py --list

作者: AI Assistant
日期: 2025-11-30
"""

import os
import sys
import glob
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.stats import linregress
from pathlib import Path
from datetime import datetime

# 设置高质量论文图样式 - Arial (Nature/Science/ACS推荐)
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
plt.rcParams['mathtext.fontset'] = 'dejavusans'
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['xtick.major.width'] = 1.5
plt.rcParams['ytick.major.width'] = 1.5

# 字体大小常量
FONT_TICK = 28
FONT_LABEL = 34

# 相态颜色配置
PHASE_COLORS = {
    'Solid': '#3498db',      # 蓝色 - 固态
    'Pre-melting': '#f39c12', # 橙色 - 预熔化
    'Liquid': '#e74c3c',     # 红色 - 液态
    'Phase_1': '#3498db',
    'Phase_2': '#e74c3c',
    'Phase_3': '#f39c12',
}

# 载体热容 (meV/K)
CV_SUPPORT = 38.2151


def find_clustering_results(base_dir='results/step6_1_clustering', method='auto'):
    """
    查找所有可用的聚类结果
    
    Parameters:
    -----------
    base_dir : str
        聚类结果目录
    method : str
        聚类方法优先级: 'auto', 'lindemann-threshold', 'kmeans'
        - 'auto': 优先使用lindemann-threshold,如果不存在则使用kmeans
        - 'lindemann-threshold': 只查找lindemann阈值方法
        - 'kmeans': 只查找kmeans方法
    
    Returns:
    --------
    dict : {structure_name: file_path}
    """
    results = {}
    
    # 定义搜索模式优先级
    if method == 'auto':
        patterns = [
            '*_lindemann-threshold_n2_clustered_data.csv',
            '*_kmeans_n2_clustered_data.csv'
        ]
        suffixes = [
            '_lindemann-threshold_n2_clustered_data.csv',
            '_kmeans_n2_clustered_data.csv'
        ]
    elif method == 'lindemann-threshold':
        patterns = ['*_lindemann-threshold_n2_clustered_data.csv']
        suffixes = ['_lindemann-threshold_n2_clustered_data.csv']
    elif method == 'kmeans':
        patterns = ['*_kmeans_n2_clustered_data.csv']
        suffixes = ['_kmeans_n2_clustered_data.csv']
    else:
        raise ValueError(f"Unknown method: {method}")
    
    # 查找文件
    found_structures = set()
    for pattern, suffix in zip(patterns, suffixes):
        files = glob.glob(os.path.join(base_dir, pattern))
        for f in files:
            basename = os.path.basename(f)
            structure = basename.replace(suffix, '')
            
            # 优先级:只添加尚未找到的结构
            if structure not in found_structures:
                results[structure] = f
                found_structures.add(structure)
    
    return results


def classify_structure(name):
    """
    分类结构名称
    返回: dict with classification flags
    """
    import re
    name_lower = name.lower()
    
    pt_match = re.search(r'pt(\d+)', name_lower)
    sn_match = re.search(r'sn(\d+)', name_lower)
    o_match = re.search(r'o(\d+)', name_lower)
    
    pt_num = int(pt_match.group(1)) if pt_match else 0
    sn_num = int(sn_match.group(1)) if sn_match else 0
    o_num = int(o_match.group(1)) if o_match else 0
    
    result = {}
    
    # Air 系列
    if name_lower.startswith('air'):
        result['air'] = True
        return result
    
    # 含氧
    if o_num > 0:
        result['oxide'] = True
        return result
    
    # Pt8SnX 系列
    if pt_num == 8:
        result['pt8snx'] = True
    
    # Pt6SnX 系列
    if pt_num == 6:
        result['pt6snx'] = True
    
    # Sum8 系列 (Pt+Sn=8)
    if pt_num + sn_num == 8 and o_num == 0:
        result['sum8'] = True
    
    return result


def parse_exclude_points(exclude_args):
    """
    解析要排除的点
    
    参数:
        exclude_args: list of str, 例如 ["300K:1,2", "400K:0"]
    
    返回:
        dict: {temp: [indices]} 例如 {300: [1, 2], 400: [0]}
    """
    exclude_dict = {}
    if not exclude_args:
        return exclude_dict
    
    for arg in exclude_args:
        try:
            temp_str, indices_str = arg.split(':')
            temp = int(temp_str.replace('K', ''))
            indices = [int(i) for i in indices_str.split(',')]
            exclude_dict[temp] = indices
        except Exception as e:
            print(f"  ⚠️ 警告: 无法解析排除点 '{arg}': {e}")
    
    return exclude_dict


def filter_data_by_exclusion(df, exclude_dict, sort_by='delta'):
    """
    根据排除规则过滤数据
    
    参数:
        df: DataFrame 包含 'temp', 'delta', 'avg_energy' 列
        exclude_dict: dict {temp: [indices]}
        sort_by: 排序依据 ('delta' 或 'energy')
    
    返回:
        DataFrame: 过滤后的数据
    """
    if not exclude_dict:
        return df
    
    print(f"\n  [排除点过滤] 排序依据: {sort_by}")
    
    # 创建掩码
    mask = np.ones(len(df), dtype=bool)
    
    for temp, indices in exclude_dict.items():
        # 获取该温度的所有数据
        temp_mask = df['temp'] == temp
        temp_indices = np.where(temp_mask)[0]
        
        if len(temp_indices) == 0:
            print(f"    警告: 温度 {temp}K 没有数据点")
            continue
        
        # 按指定字段排序
        temp_df = df[temp_mask].copy()
        temp_df['original_idx'] = temp_indices
        
        if sort_by == 'delta':
            # 按 Lindemann 指数排序（从大到小）
            temp_df_sorted = temp_df.sort_values('delta', ascending=False)
            sort_label = 'Lindemann'
        else:  # sort_by == 'energy'
            # 按能量排序（从大到小）
            temp_df_sorted = temp_df.sort_values('avg_energy', ascending=False)
            sort_label = 'Energy'
        
        # 标记要排除的点
        for idx in indices:
            if idx < len(temp_df_sorted):
                original_idx = int(temp_df_sorted.iloc[idx]['original_idx'])
                mask[original_idx] = False
                
                if sort_by == 'delta':
                    value = temp_df_sorted.iloc[idx]['delta']
                    print(f"    排除: {temp}K 第{idx}个点 (delta={value:.4f})")
                else:
                    value = temp_df_sorted.iloc[idx]['avg_energy']
                    print(f"    排除: {temp}K 第{idx}个点 (energy={value:.4f} eV)")
            else:
                print(f"    警告: {temp}K 索引{idx}超出范围 (最大={len(temp_df_sorted)-1})")
    
    filtered_df = df[mask].copy()
    n_excluded = len(df) - len(filtered_df)
    print(f"  过滤结果: 原始{len(df)}条 → 保留{len(filtered_df)}条 (排除{n_excluded}条)\n")
    
    return filtered_df


def filter_structures_by_series(available, series_list):
    """
    根据系列筛选结构
    
    Args:
        available: dict of structure_name -> csv_path
        series_list: list of series names like ['pt8snx', 'pt6snx']
    
    Returns: filtered list of structure names
    """
    filtered = []
    for name in available.keys():
        classification = classify_structure(name)
        for series in series_list:
            if classification.get(series, False):
                filtered.append(name)
                break
    return filtered


def load_support_energy_data():
    """加载载体能量数据"""
    support_csv = 'data/lammps_energy/sup/energy_master_20251021_151520.csv'
    
    if not os.path.exists(support_csv):
        return None
    
    try:
        df_support = pd.read_csv(support_csv)
        if 'temp' in df_support.columns and 'avg_energy' in df_support.columns:
            T = df_support['temp'].values
            E = df_support['avg_energy'].values
            slope, intercept, r_value, _, _ = linregress(T, E)
            return slope, intercept, r_value**2
    except Exception as e:
        print(f"  警告: 读取载体能量数据失败: {e}")
    
    return None


def load_cluster_data(csv_path, exclude_dict=None, exclude_sort_by='delta'):
    """
    加载聚类结果数据并过滤
    
    参数:
        csv_path: CSV文件路径
        exclude_dict: 排除点字典 {temp: [indices]}
        exclude_sort_by: 排序依据 ('delta' 或 'energy')
    """
    try:
        df = pd.read_csv(csv_path)
        
        if exclude_dict:
            print(f"\n  加载数据: {csv_path}")
            print(f"  原始数据: {len(df)} 条")
            df = filter_data_by_exclusion(df, exclude_dict, sort_by=exclude_sort_by)
        
        return df
    except Exception as e:
        print(f"  错误: 无法读取 {csv_path}: {e}")
        return None


def determine_partitions_by_lindemann(df, threshold=0.1):
    """
    根据固定Lindemann阈值自动确定分区温度范围
    
    参数:
        df: DataFrame，必须包含 'temp' 和 'delta' 列
        threshold: Lindemann阈值（默认0.1）
        
    返回:
        tuple: (custom_partitions, temp_to_partition_dict)
            - custom_partitions: [(T_min1, T_max1), (T_min2, T_max2)]
            - temp_to_partition_dict: {temp: partition_id}
    
    逻辑:
        1. 对每个温度计算平均Lindemann指数
        2. δ_avg < threshold → 分区1 (固相)
        3. δ_avg >= threshold → 分区2 (液相)
        4. 找到分区边界温度
    """
    print(f"\n  [Lindemann阈值分区] 阈值 δ = {threshold:.3f}")
    
    # 按温度分组计算平均Lindemann指数
    temp_delta = df.groupby('temp')['delta'].agg(['mean', 'std', 'count']).reset_index()
    temp_delta = temp_delta.sort_values('temp')
    
    # 根据阈值分类每个温度
    temp_delta['partition'] = temp_delta['mean'].apply(
        lambda x: 1 if x < threshold else 2
    )
    
    print(f"\n  温度分区结果:")
    print(f"  {'温度(K)':<10} {'δ平均':<10} {'δ标准差':<10} {'点数':<6} {'分区':<6}")
    print(f"  {'-'*50}")
    for _, row in temp_delta.iterrows():
        print(f"  {row['temp']:<10.0f} {row['mean']:<10.4f} {row['std']:<10.4f} "
              f"{row['count']:<6.0f} {row['partition']:<6.0f}")
    
    # 创建温度到分区的映射
    temp_to_partition = dict(zip(temp_delta['temp'], temp_delta['partition']))
    
    # 确定分区温度范围
    partition1_temps = temp_delta[temp_delta['partition'] == 1]['temp'].values
    partition2_temps = temp_delta[temp_delta['partition'] == 2]['temp'].values
    
    custom_partitions = []
    
    if len(partition1_temps) > 0:
        T1_min = partition1_temps.min()
        T1_max = partition1_temps.max()
        custom_partitions.append((T1_min, T1_max))
        print(f"\n  分区1 (固相, δ<{threshold}): {T1_min:.0f}-{T1_max:.0f} K ({len(partition1_temps)} 个温度)")
    else:
        print(f"\n  警告: 没有温度点被分到分区1 (固相)")
    
    if len(partition2_temps) > 0:
        T2_min = partition2_temps.min()
        T2_max = partition2_temps.max()
        custom_partitions.append((T2_min, T2_max))
        print(f"  分区2 (液相, δ≥{threshold}): {T2_min:.0f}-{T2_max:.0f} K ({len(partition2_temps)} 个温度)")
    else:
        print(f"  警告: 没有温度点被分到分区2 (液相)")
    
    if len(custom_partitions) < 2:
        print(f"\n  警告: Lindemann阈值 {threshold} 未能产生两个分区，可能需要调整阈值")
        return None, None
    
    # 检查分区连续性
    if len(partition1_temps) > 0 and len(partition2_temps) > 0:
        T1_last = partition1_temps.max()
        T2_first = partition2_temps.min()
        if T2_first > T1_last:
            gap = T2_first - T1_last
            print(f"\n  分区边界: {T1_last:.0f} K (固相最高) → {T2_first:.0f} K (液相最低)")
            if gap > 0:
                print(f"  过渡区间: {gap:.0f} K")
        else:
            print(f"\n  警告: 分区温度交叉! 固相最高温 {T1_last:.0f} K ≥ 液相最低温 {T2_first:.0f} K")
            print(f"        建议调整阈值或检查数据")
    
    return custom_partitions, temp_to_partition


def plot_partition_cv(df, structure_name, output_dir, output_format='png', dpi=300, tick_params=None, custom_partitions=None, peak_method='fit', remove_outliers=False, outlier_iqr=1.5):
    """
    绘制分区热容拟合图（论文出图专用）
    
    核心逻辑：
    1. 按温度分组计算团簇能量平均值和标准差
    2. 使用多数投票规则将每个温度分配给唯一的相态（或使用自定义分区）
    3. 对每个相态的专属温度点进行线性拟合
    4. 绘制整体拟合线 vs 分区拟合线对比
    
    tick_params: 刻度参数字典
        - y_ticks_custom: 自定义Y轴刻度列表
        - cv_ticks_custom: 自定义Cv轴刻度列表
        - y_nticks: Y轴刻度数量
        - cv_nticks: Cv轴刻度数量
        - figsize: 图片尺寸 (宽, 高)
    
    custom_partitions: 自定义分区列表，格式为 [(T_min1, T_max1), (T_min2, T_max2), ...]
        例如: [(200, 700), (750, 950)] 表示第一分区200-700K，第二分区750-950K
        如果为 None，则使用聚类结果的多数投票规则
    
    peak_method: 热容峰计算方法
        - 'data': 数据点法 - 使用实际数据点的能量差计算过渡区热容
        - 'fit': 拟合线外推法 - 使用拟合线外推的能量差计算过渡区热容
    
    remove_outliers: 是否剔除离群点
    outlier_iqr: IQR倍数阈值，默认1.5
    """
    
    # 默认刻度参数
    if tick_params is None:
        tick_params = {}
    y_ticks_custom = tick_params.get('y_ticks_custom', None)
    cv_ticks_custom = tick_params.get('cv_ticks_custom', None)
    y_nticks = tick_params.get('y_nticks', 5)
    cv_nticks = tick_params.get('cv_nticks', 5)
    figsize = tick_params.get('figsize', (10, 10))
    
    print(f"\n>>> 绘制 {structure_name} 分区热容图...")
    
    # 检查必要列
    required_cols = ['temp', 'avg_energy', 'phase_clustered']
    if not all(col in df.columns for col in required_cols):
        print(f"  错误: 缺少必要列 {required_cols}")
        return None
    
    # 判断是否是 Air 系列（气相团簇）
    is_air_system = structure_name.startswith('Air') or structure_name in ['68', '86']
    
    # 加载载体能量数据
    if is_air_system:
        slope_support = 0.0
        intercept_support = 0.0
        print(f"  [Air系列] 气相纳米团簇，不扣除载体能量")
    else:
        support_fit = load_support_energy_data()
        if support_fit is not None:
            slope_support, intercept_support, R2_support = support_fit
            print(f"  [载体数据] Cv_support={slope_support*1000:.4f} meV/K, R²={R2_support:.6f}")
        else:
            slope_support = CV_SUPPORT / 1000  # meV/K -> eV/K
            T_min = df['temp'].min()
            E_total_min = df[df['temp'] == T_min]['avg_energy'].mean()
            intercept_support = E_total_min * 0.9 - slope_support * T_min
            print(f"  [警告] 使用默认Cv_support估算载体能量")
    
    # ========== 1. 按温度分组计算团簇能量（可选离群点剔除） ==========
    temp_groups = df.groupby('temp')
    temps_unique = []
    E_cluster_mean = []
    E_cluster_std = []
    outlier_stats = {'total_removed': 0, 'by_temp': {}, 'details': {}}
    
    # 解析离群点剔除参数
    outlier_method = 'iqr'  # 默认方法
    outlier_threshold = outlier_iqr
    outlier_iterations = 1  # 默认迭代次数
    
    if isinstance(remove_outliers, str):
        # 解析格式: "method:threshold:iterations" 或 "method:threshold" 或 "method"
        parts = remove_outliers.split(':')
        outlier_method = parts[0].lower()
        if len(parts) >= 2:
            outlier_threshold = float(parts[1])
        if len(parts) >= 3:
            outlier_iterations = int(parts[2])
    elif remove_outliers is True:
        outlier_method = 'iqr'
        outlier_threshold = outlier_iqr
    
    for temp, group in temp_groups:
        if is_air_system:
            E_cluster = group['avg_energy'].values.copy()
        else:
            E_support = slope_support * temp + intercept_support
            E_cluster = group['avg_energy'].values - E_support
        
        # 离群点剔除 - 仅当 remove_outliers 不为 False 时启用
        if remove_outliers:
            n_before = len(E_cluster)
            removed_values = []
            
            # 迭代剔除
            for iteration in range(outlier_iterations):
                if len(E_cluster) < 4:
                    break
                    
                if outlier_method == 'iqr':
                    # IQR 方法
                    Q1 = np.percentile(E_cluster, 25)
                    Q3 = np.percentile(E_cluster, 75)
                    IQR = Q3 - Q1
                    lower_bound = Q1 - outlier_threshold * IQR
                    upper_bound = Q3 + outlier_threshold * IQR
                    
                elif outlier_method == 'zscore':
                    # Z-score 方法
                    mean = np.mean(E_cluster)
                    std = np.std(E_cluster)
                    if std > 0:
                        lower_bound = mean - outlier_threshold * std
                        upper_bound = mean + outlier_threshold * std
                    else:
                        lower_bound = mean - 1e-10
                        upper_bound = mean + 1e-10
                        
                elif outlier_method == 'mad':
                    # MAD (Median Absolute Deviation) 方法 - 对异常值更鲁棒
                    median = np.median(E_cluster)
                    mad = np.median(np.abs(E_cluster - median))
                    if mad > 0:
                        lower_bound = median - outlier_threshold * 1.4826 * mad
                        upper_bound = median + outlier_threshold * 1.4826 * mad
                    else:
                        lower_bound = median - 1e-10
                        upper_bound = median + 1e-10
                        
                elif outlier_method == 'percentile':
                    # 百分位数方法 - 直接剔除最极端的点
                    lower_bound = np.percentile(E_cluster, outlier_threshold)
                    upper_bound = np.percentile(E_cluster, 100 - outlier_threshold)
                    
                else:
                    # 默认 IQR
                    Q1 = np.percentile(E_cluster, 25)
                    Q3 = np.percentile(E_cluster, 75)
                    IQR = Q3 - Q1
                    lower_bound = Q1 - outlier_threshold * IQR
                    upper_bound = Q3 + outlier_threshold * IQR
                
                mask = (E_cluster >= lower_bound) & (E_cluster <= upper_bound)
                outliers_this_iter = E_cluster[~mask]
                removed_values.extend(outliers_this_iter.tolist())
                E_cluster = E_cluster[mask]
                
                if len(outliers_this_iter) == 0:
                    break  # 没有更多离群点
            
            n_removed = n_before - len(E_cluster)
            if n_removed > 0:
                outlier_stats['total_removed'] += n_removed
                outlier_stats['by_temp'][temp] = n_removed
                outlier_stats['details'][temp] = removed_values
        
        temps_unique.append(temp)
        E_cluster_mean.append(np.mean(E_cluster))
        E_cluster_std.append(np.std(E_cluster) if len(E_cluster) > 1 else 0)
    
    # 输出离群点剔除统计
    if remove_outliers:
        method_names = {
            'iqr': f'IQR×{outlier_threshold}',
            'zscore': f'Z-score>{outlier_threshold}σ',
            'mad': f'MAD×{outlier_threshold}',
            'percentile': f'百分位{outlier_threshold}%-{100-outlier_threshold}%'
        }
        method_desc = method_names.get(outlier_method, outlier_method)
        iter_desc = f", {outlier_iterations}次迭代" if outlier_iterations > 1 else ""
        
        if outlier_stats['total_removed'] > 0:
            print(f"\n  [离群点剔除] 共剔除 {outlier_stats['total_removed']} 个点 ({method_desc}{iter_desc})")
            for temp, n in sorted(outlier_stats['by_temp'].items()):
                values = outlier_stats['details'].get(temp, [])
                values_str = ', '.join([f'{v:.4f}' for v in values[:3]])
                if len(values) > 3:
                    values_str += f'... (共{len(values)}个)'
                print(f"    T={temp:.0f}K: 剔除 {n} 个点 [{values_str}]")
        else:
            print(f"\n  [离群点剔除] 未发现离群点 ({method_desc}{iter_desc})")
    
    temps_unique = np.array(temps_unique)
    E_cluster_mean = np.array(E_cluster_mean)
    E_cluster_std = np.array(E_cluster_std)
    
    # 计算相对能量（相对于最低温度）
    E_cluster_ref = E_cluster_mean.min()
    E_cluster_mean_rel = E_cluster_mean - E_cluster_ref
    
    # ========== 2. 确定每个温度的分区 ==========
    temp_to_partition = {}
    
    if custom_partitions is not None:
        # 使用自定义分区
        print(f"\n  使用自定义分区:")
        for i, (T_min, T_max) in enumerate(custom_partitions):
            partition_name = f'partition{i+1}'
            print(f"    {partition_name}: {T_min}K - {T_max}K")
            for temp in temps_unique:
                if T_min <= temp <= T_max:
                    temp_to_partition[temp] = partition_name
        
        # 检查是否有未分配的温度
        unassigned = [t for t in temps_unique if t not in temp_to_partition]
        if unassigned:
            print(f"  警告: 以下温度未被分配到任何分区: {unassigned}")
            # 将未分配温度归入最近的分区
            for temp in unassigned:
                # 找到最近的分区边界
                min_dist = float('inf')
                nearest_partition = None
                for i, (T_min, T_max) in enumerate(custom_partitions):
                    partition_name = f'partition{i+1}'
                    dist = min(abs(temp - T_min), abs(temp - T_max))
                    if dist < min_dist:
                        min_dist = dist
                        nearest_partition = partition_name
                temp_to_partition[temp] = nearest_partition
                print(f"    T={temp}K → {nearest_partition} (按最近原则)")
    else:
        # 使用多数投票规则（原有逻辑）
        print(f"\n  多数投票温度分配:")
        for temp in temps_unique:
            df_temp = df[df['temp'] == temp]
            partition_counts = df_temp['phase_clustered'].value_counts()
            dominant_partition = partition_counts.idxmax()
            temp_to_partition[temp] = dominant_partition
            print(f"    T={temp:4.0f}K: {dict(partition_counts)} → {dominant_partition}")
    
    # ========== 3. 整体拟合 ==========
    if len(temps_unique) < 3:
        print(f"  错误: 温度点不足 ({len(temps_unique)} < 3)")
        return None
    
    slope_overall, intercept_overall, r_value_overall, _, std_err_overall = linregress(
        temps_unique, E_cluster_mean_rel)
    R2_overall = r_value_overall ** 2
    Cv_overall = slope_overall * 1000  # meV/K
    Cv_overall_err = std_err_overall * 1000
    
    print(f"\n  整体拟合: Cv={Cv_overall:.4f}±{Cv_overall_err:.4f} meV/K, R2={R2_overall:.4f}")
    
    # ========== 4. 分区拟合 ==========
    phases = df['phase_clustered'].unique()
    phase_fits = {}
    
    for phase in phases:
        phase_temps = [temp for temp, part in temp_to_partition.items() if part == phase]
        phase_temps = sorted(phase_temps)
        
        if len(phase_temps) >= 2:
            mask = np.isin(temps_unique, phase_temps)
            T_phase = temps_unique[mask]
            E_phase_rel = E_cluster_mean_rel[mask]
            E_phase_std = E_cluster_std[mask]
            
            slope_ph, intercept_ph, r_value_ph, _, std_err_ph = linregress(T_phase, E_phase_rel)
            R2_ph = r_value_ph ** 2
            Cv_ph = slope_ph * 1000
            Cv_ph_err = std_err_ph * 1000
            
            phase_fits[phase] = {
                'slope': slope_ph,
                'intercept': intercept_ph,
                'R2': R2_ph,
                'Cv': Cv_ph,
                'Cv_err': Cv_ph_err,
                'n_temps': len(T_phase),
                'T_range': (T_phase.min(), T_phase.max()),
                'T_data': T_phase,
                'E_data': E_phase_rel,
                'E_std': E_phase_std
            }
            
            print(f"  {phase}: Cv={Cv_ph:.4f}±{Cv_ph_err:.4f} meV/K, R²={R2_ph:.4f}, "
                  f"n={len(T_phase)}, T={T_phase.min():.0f}-{T_phase.max():.0f}K")
    
    # ========== 5. 绘制简洁的双Y轴图 ==========
    fig, ax1 = plt.subplots(figsize=figsize)
    
    # ----- 左Y轴: 能量-温度数据点（带误差棒）和拟合线 -----
    # 绘制数据点（带误差棒）
    ax1.errorbar(temps_unique, E_cluster_mean_rel, yerr=E_cluster_std,
                 fmt='o', markersize=10, color='black', 
                 ecolor='gray', elinewidth=2, capsize=4, capthick=2,
                 zorder=5, label='Data')
    
    # 绘制拟合线（黑色）
    phases_sorted = sorted(phase_fits.keys())
    for phase in phases_sorted:
        fit = phase_fits[phase]
        T_phase_fit = np.linspace(fit['T_range'][0], fit['T_range'][1], 50)
        E_phase_fit = fit['slope'] * T_phase_fit + fit['intercept']
        ax1.plot(T_phase_fit, E_phase_fit, '-', color='black', linewidth=2.5, zorder=4)
    
    # 连接两个分区之间的数据点（实线连接实际数据点，而非拟合线）
    if len(phases_sorted) >= 2:
        fit1 = phase_fits[phases_sorted[0]]
        fit2 = phase_fits[phases_sorted[1]]
        # 分区1的最后一个数据点
        T1_end = fit1['T_range'][1]
        idx1 = np.where(temps_unique == T1_end)[0]
        if len(idx1) > 0:
            E1_end = E_cluster_mean_rel[idx1[0]]
        else:
            E1_end = fit1['slope'] * T1_end + fit1['intercept']
        # 分区2的第一个数据点
        T2_start = fit2['T_range'][0]
        idx2 = np.where(temps_unique == T2_start)[0]
        if len(idx2) > 0:
            E2_start = E_cluster_mean_rel[idx2[0]]
        else:
            E2_start = fit2['slope'] * T2_start + fit2['intercept']
        # 用实线连接两个数据点
        ax1.plot([T1_end, T2_start], [E1_end, E2_start], '-', color='black', linewidth=2.5, zorder=4)
    
    ax1.set_xlabel('Temperature (K)', fontsize=FONT_LABEL)
    ax1.set_ylabel('Total Energy (eV)', fontsize=FONT_LABEL)
    ax1.tick_params(axis='both', labelsize=FONT_TICK)
    
    # 设置Y轴刻度
    E_ylim = ax1.get_ylim()
    if y_ticks_custom is not None:
        ax1.set_yticks(y_ticks_custom)
    else:
        y_ticks = np.linspace(E_ylim[0], E_ylim[1], y_nticks)
        y_ticks = np.round(y_ticks)
        ax1.set_yticks(y_ticks)
    
    # ----- 右Y轴: 热容曲线 -----
    ax2 = ax1.twinx()
    
    # 分界温度和热容值（用于导出）
    T_boundary = None
    Cv1 = None
    Cv2 = None
    Cv_peak = None
    
    if len(phases_sorted) >= 2:
        # 找到分区边界温度
        phase1_temps = [t for t, p in temp_to_partition.items() if p == phases_sorted[0]]
        phase2_temps = [t for t, p in temp_to_partition.items() if p == phases_sorted[1]]
        
        if phase1_temps and phase2_temps:
            T1_last = max(phase1_temps)   # 分区1最后一个温度
            T2_first = min(phase2_temps)  # 分区2第一个温度
            T_boundary = (T1_last + T2_first) / 2
            print(f"\n  分界温度: {T_boundary:.0f} K (过渡区: {T1_last:.0f}-{T2_first:.0f}K)")
            
            Cv1 = phase_fits[phases_sorted[0]]['Cv']
            Cv2 = phase_fits[phases_sorted[1]]['Cv']
            
            # ========== 方法1: 实际数据点的数值微分 ==========
            # 注意：这里使用的是按温度平均的能量（所有run的平均，不区分分区归属）
            idx1 = np.where(temps_unique == T1_last)[0]
            idx2 = np.where(temps_unique == T2_first)[0]
            if len(idx1) > 0 and len(idx2) > 0:
                E1_data = E_cluster_mean_rel[idx1[0]]
                E2_data = E_cluster_mean_rel[idx2[0]]
                Cv_transition_data = (E2_data - E1_data) / (T2_first - T1_last) * 1000  # meV/K
            else:
                Cv_transition_data = (Cv1 + Cv2) / 2
            
            # ========== 方法1b: 只用归属于该分区的点的能量 ==========
            # T1_last (如650K) 被分给分区1，只用分区1的点计算能量
            # T2_first (如700K) 被分给分区2，只用分区2的点计算能量
            df_T1 = df[df['temp'] == T1_last]
            df_T2 = df[df['temp'] == T2_first]
            
            # 获取该温度被分配到的分区
            partition_T1 = temp_to_partition[T1_last]
            partition_T2 = temp_to_partition[T2_first]
            
            # 筛选只属于该分区的点
            df_T1_filtered = df_T1[df_T1['phase_clustered'] == partition_T1]
            df_T2_filtered = df_T2[df_T2['phase_clustered'] == partition_T2]
            
            if len(df_T1_filtered) > 0 and len(df_T2_filtered) > 0:
                if is_air_system:
                    E1_partition = df_T1_filtered['avg_energy'].mean()
                    E2_partition = df_T2_filtered['avg_energy'].mean()
                else:
                    E_support_T1 = slope_support * T1_last + intercept_support
                    E_support_T2 = slope_support * T2_first + intercept_support
                    E1_partition = df_T1_filtered['avg_energy'].mean() - E_support_T1
                    E2_partition = df_T2_filtered['avg_energy'].mean() - E_support_T2
                
                # 转换为相对能量
                E1_partition_rel = E1_partition - E_cluster_ref - (E_cluster_mean[0] - E_cluster_ref)
                E2_partition_rel = E2_partition - E_cluster_ref - (E_cluster_mean[0] - E_cluster_ref)
                
                # 直接用绝对能量差计算
                Cv_transition_partition = (E2_partition - E1_partition) / (T2_first - T1_last) * 1000
                
                n_T1_total = len(df_T1)
                n_T1_used = len(df_T1_filtered)
                n_T2_total = len(df_T2)
                n_T2_used = len(df_T2_filtered)
                print(f"    分区点法: T1={T1_last}K 用{n_T1_used}/{n_T1_total}点({partition_T1}), "
                      f"T2={T2_first}K 用{n_T2_used}/{n_T2_total}点({partition_T2})")
                print(f"    分区点法: Cv_transition={Cv_transition_partition:.2f} meV/K")
            else:
                Cv_transition_partition = Cv_transition_data
                print(f"    分区点法: 数据不足，回退到全部数据点")
            
            # ========== 方法2: 拟合线外推的能量差 ==========
            # E1_fit: 用分区1拟合线代入T1_last
            # E2_fit: 用分区2拟合线代入T2_first
            fit1 = phase_fits[phases_sorted[0]]
            fit2 = phase_fits[phases_sorted[1]]
            E1_fit = fit1['slope'] * T1_last + fit1['intercept']
            E2_fit = fit2['slope'] * T2_first + fit2['intercept']
            Cv_transition_fit = (E2_fit - E1_fit) / (T2_first - T1_last) * 1000  # meV/K
            
            print(f"  热容峰计算方法: {peak_method}")
            print(f"    全点法(data): Cv={Cv_transition_data:.2f} meV/K (所有点平均)")
            print(f"    分区法(partition): Cv={Cv_transition_partition:.2f} meV/K (只用分区内点)")
            print(f"    拟合法(fit): Cv={Cv_transition_fit:.2f} meV/K (拟合线外推)")
            print(f"    Cv1={Cv1:.2f}, Cv2={Cv2:.2f} meV/K")
            
            # 根据选择的方法判断热容峰
            if peak_method == 'data':
                Cv_transition = Cv_transition_data
                method_name = "全点数据法"
            elif peak_method == 'partition':
                Cv_transition = Cv_transition_partition
                method_name = "分区点法"
            else:  # 'fit'
                Cv_transition = Cv_transition_fit
                method_name = "拟合线外推法"
            
            # 判断是否存在热容峰
            has_peak = Cv_transition > max(Cv1, Cv2)
            
            if has_peak:
                Cv_peak = Cv_transition
                print(f"  ★ 存在热容峰: Cv_peak={Cv_peak:.2f} meV/K ({method_name})")
                print(f"  热容: Cv1={Cv1:.2f}, Cv_peak={Cv_peak:.2f}, Cv2={Cv2:.2f} meV/K")
                
                # 绘制带平滑峰的热容曲线（使用高斯峰 + sigmoid过渡）
                T_plot = np.linspace(temps_unique.min(), temps_unique.max(), 500)
                Cv_plot = np.zeros_like(T_plot)
                
                # 峰的宽度参数
                sigma = (T2_first - T1_last) / 2  # 高斯宽度
                
                for i, T in enumerate(T_plot):
                    # 基线：sigmoid 从 Cv1 过渡到 Cv2
                    transition = 1 / (1 + np.exp(-(T - T_boundary) / (sigma * 0.5)))
                    baseline = Cv1 + (Cv2 - Cv1) * transition
                    
                    # 高斯峰叠加
                    gaussian = (Cv_peak - baseline) * np.exp(-0.5 * ((T - T_boundary) / sigma)**2)
                    Cv_plot[i] = baseline + gaussian
                
                ax2.plot(T_plot, Cv_plot, 'r-', linewidth=2.5, zorder=3)
                
                # 构建导出数据（关键点）
                T_cv = np.array([temps_unique.min(), T1_last, T_boundary, T2_first, temps_unique.max()])
                Cv_curve = np.array([Cv1, Cv1, Cv_peak, Cv2, Cv2])
            else:
                print(f"  热容: Cv1={Cv1:.2f} meV/K, Cv2={Cv2:.2f} meV/K (无峰)")
                
                # 绘制阶梯形热容曲线（无峰）
                ax2.plot([temps_unique.min(), T_boundary], [Cv1, Cv1], 'r-', linewidth=2.5, zorder=3)
                ax2.plot([T_boundary, T_boundary], [Cv1, Cv2], 'r--', linewidth=2, zorder=3)
                ax2.plot([T_boundary, temps_unique.max()], [Cv2, Cv2], 'r-', linewidth=2.5, zorder=3)
                
                T_cv = np.array([temps_unique.min(), T_boundary - 0.1, T_boundary, T_boundary + 0.1, temps_unique.max()])
                Cv_curve = np.array([Cv1, Cv1, (Cv1 + Cv2) / 2, Cv2, Cv2])
    else:
        Cv_single = list(phase_fits.values())[0]['Cv']
        T_cv = np.array([temps_unique.min(), temps_unique.max()])
        Cv_curve = np.array([Cv_single, Cv_single])
        ax2.plot(T_cv, Cv_curve, 'r-', linewidth=2.5, zorder=3)
        Cv1 = Cv_single
        Cv2 = Cv_single
    
    ax2.set_ylabel(r'$C_v$ (meV/K)', fontsize=FONT_LABEL, color='red')
    ax2.tick_params(axis='y', labelcolor='red', labelsize=FONT_TICK, color='red')
    ax2.spines['right'].set_color('red')
    
    # 设置Y轴范围（考虑峰值）
    cv_values = [Cv1, Cv2] if Cv1 and Cv2 else list(Cv_curve)
    if Cv_peak:
        cv_values.append(Cv_peak)
    cv_min = min(cv_values) * 0.85
    cv_max = max(cv_values) * 1.1
    ax2.set_ylim(cv_min, cv_max)
    
    # 设置Cv轴刻度
    Cv_ylim = ax2.get_ylim()
    if cv_ticks_custom is not None:
        ax2.set_yticks(cv_ticks_custom)
    else:
        cv_ticks = np.linspace(Cv_ylim[0], Cv_ylim[1], cv_nticks)
        cv_ticks = np.round(cv_ticks)
        ax2.set_yticks(cv_ticks)
    
    plt.tight_layout()
    
    # 保存图片
    output_file = Path(output_dir) / f'{structure_name}_partition_cv.{output_format}'
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close()
    
    print(f"\n  图已保存: {output_file}")
    
    # ========== 6. 导出数据供 Origin 使用 ==========
    # 导出能量数据
    df_energy = pd.DataFrame({
        'Temperature_K': temps_unique,
        'Energy_eV': E_cluster_mean_rel,
        'Energy_std_eV': E_cluster_std,
        'Partition': [temp_to_partition.get(t, 'unknown') for t in temps_unique]
    })
    energy_csv = Path(output_dir) / f'{structure_name}_energy_data.csv'
    df_energy.to_csv(energy_csv, index=False)
    print(f"  能量数据已导出: {energy_csv}")
    
    # 导出热容数据（阶梯函数关键点）
    df_cv = pd.DataFrame({
        'Temperature_K': T_cv,
        'Cv_meV_K': Cv_curve
    })
    cv_csv = Path(output_dir) / f'{structure_name}_cv_curve.csv'
    df_cv.to_csv(cv_csv, index=False)
    print(f"  热容曲线已导出: {cv_csv}")
    
    # 导出拟合参数汇总
    fit_summary = {
        'structure': structure_name,
        'T_boundary_K': T_boundary,
        'Cv_overall_meV_K': Cv_overall,
        'Cv_overall_err': Cv_overall_err,
        'R2_overall': R2_overall,
    }
    for i, (phase, fit) in enumerate(phase_fits.items()):
        fit_summary[f'phase_{i+1}_name'] = phase
        fit_summary[f'phase_{i+1}_Cv_meV_K'] = fit['Cv']
        fit_summary[f'phase_{i+1}_Cv_err'] = fit['Cv_err']
        fit_summary[f'phase_{i+1}_R2'] = fit['R2']
        fit_summary[f'phase_{i+1}_T_min_K'] = fit['T_range'][0]
        fit_summary[f'phase_{i+1}_T_max_K'] = fit['T_range'][1]
        fit_summary[f'phase_{i+1}_slope_eV_K'] = fit['slope']
        fit_summary[f'phase_{i+1}_intercept_eV'] = fit['intercept']
    
    fit_csv = Path(output_dir) / f'{structure_name}_fit_params.csv'
    pd.DataFrame([fit_summary]).to_csv(fit_csv, index=False)
    print(f"  拟合参数已导出: {fit_csv}")
    
    # 返回拟合结果
    return {
        'structure': structure_name,
        'overall': {'Cv': Cv_overall, 'Cv_err': Cv_overall_err, 'R2': R2_overall},
        'partitions': phase_fits
    }


def list_available_structures(base_dir='results/step6_1_clustering', method='auto'):
    """列出所有可用的结构"""
    results = find_clustering_results(base_dir, method=method)
    
    print("\n" + "=" * 60)
    print("可用结构列表")
    print("=" * 60)
    
    # 分类
    air_series = []
    pt6_series = []
    pt8_series = []
    oxide_series = []
    other = []
    
    for name in sorted(results.keys()):
        name_lower = name.lower()
        if 'air' in name_lower:
            air_series.append(name)
        elif name == 'Cv' or 'o' in name_lower:
            oxide_series.append(name)
        elif name_lower.startswith('pt6'):
            pt6_series.append(name)
        elif name_lower.startswith('pt8'):
            pt8_series.append(name)
        else:
            other.append(name)
    
    print(f"\n🔵 气相团簇 ({len(air_series)}): {', '.join(air_series) if air_series else '无'}")
    print(f"🟢 Pt6系列 ({len(pt6_series)}): {', '.join(sorted(pt6_series)) if pt6_series else '无'}")
    print(f"🟢 Pt8系列 ({len(pt8_series)}): {', '.join(sorted(pt8_series)) if pt8_series else '无'}")
    print(f"🟠 含氧团簇 ({len(oxide_series)}): {', '.join(sorted(oxide_series)) if oxide_series else '无'}")
    if other:
        print(f"⚪ 其他 ({len(other)}): {', '.join(sorted(other))}")
    
    print(f"\n总计: {len(results)} 个结构")
    print("=" * 60)
    
    return results


def parse_args():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description='分区热容拟合图 - 论文出图专用',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例:
  %(prog)s --structure Pt8sn6              # 单个结构
  %(prog)s --structure Air86 --format pdf  # 输出PDF
  %(prog)s --structure all --dpi 600       # 所有结构，高分辨率
  %(prog)s --only-series pt8snx            # 批量绘制 Pt8SnX 系列
  %(prog)s --only-series pt6snx,pt8snx     # 批量绘制多个系列
  %(prog)s --list                          # 列出可用结构
  %(prog)s --structure Pt8sn6 --y-ticks 0,2,4 --cv-ticks 3,4,5,6,7  # 自定义刻度
  
自定义分区:
  %(prog)s --structure o2pt7sn7 --partitions 200-700,750-950   # 手动指定分区
  %(prog)s --structure Pt6sn8 --partitions 200-550,600-950     # 第一分区200-550K，第二分区600-950K

Lindemann阈值自动分区（新功能）:
  %(prog)s --structure Pt8sn6 --use-lindemann-threshold 0.1    # δ<0.1固相，δ≥0.1液相
  %(prog)s --structure Pt6sn8 --use-lindemann-threshold 0.08   # 使用δ=0.08作为阈值
  %(prog)s --structure Air86 --use-lindemann-threshold 0.12    # 气相团簇，δ=0.12阈值
  %(prog)s --structure all --use-lindemann-threshold 0.1       # 批量处理，统一阈值

热容峰计算方法:
  %(prog)s --structure o2pt7sn7 --peak-method data      # 全点数据法（所有点平均）
  %(prog)s --structure o2pt7sn7 --peak-method partition # 分区点法（只用归属分区的点）★推荐
  %(prog)s --structure o2pt7sn7 --peak-method fit       # 拟合线外推法

离群点剔除（增强版）:
  %(prog)s --structure o2pt7sn7 --remove-outliers                    # 默认 IQR×1.5
  %(prog)s --structure o2pt7sn7 --remove-outliers --outlier-iqr 1.0  # 更严格 IQR×1.0
  %(prog)s --structure o2pt7sn7 --outlier-method zscore:2            # Z-score 2σ
  %(prog)s --structure o2pt7sn7 --outlier-method zscore:1.5:3        # Z-score 1.5σ, 3次迭代
  %(prog)s --structure o2pt7sn7 --outlier-method mad:2.5             # MAD方法（更鲁棒）
  %(prog)s --structure o2pt7sn7 --outlier-method percentile:5        # 剔除最极端5%%的点
        '''
    )
    
    parser.add_argument('--structure', '-s', type=str, default=None,
                        help='结构名称 (如 Pt8sn6, Air86) 或 "all" 处理所有')
    parser.add_argument('--only-series', type=str, default=None,
                        help='只处理指定系列（逗号分隔）: pt8snx, pt6snx, air, oxide, sum8')
    parser.add_argument('--list', '-l', action='store_true',
                        help='列出所有可用结构')
    parser.add_argument('--clustering-method', type=str, default='auto',
                        choices=['auto', 'lindemann-threshold', 'kmeans'],
                        help='聚类数据文件选择方法: '
                             'auto(默认,优先lindemann-threshold), '
                             'lindemann-threshold(仅Lindemann阈值数据), '
                             'kmeans(仅KMeans聚类数据)')
    parser.add_argument('--format', '-f', type=str, default='png',
                        choices=['png', 'pdf', 'svg', 'eps'],
                        help='输出格式 (默认: png)')
    parser.add_argument('--dpi', type=int, default=300,
                        help='输出分辨率 (默认: 300)')
    parser.add_argument('--output-dir', '-o', type=str, 
                        default='results/step6_1_1_partition_cv',
                        help='输出目录')
    parser.add_argument('--figsize', type=str, default='10x10',
                        help='图片尺寸，格式: 宽x高，例如 10x8 (默认: 10x10)')
    parser.add_argument('--y-ticks', type=str, default=None,
                        help='手动指定能量Y轴刻度，逗号分隔，例如: 0,2,4')
    parser.add_argument('--cv-ticks', type=str, default=None,
                        help='手动指定Cv轴刻度，逗号分隔，例如: 3,4,5,6,7')
    parser.add_argument('--y-nticks', type=int, default=5,
                        help='能量Y轴刻度数量 (默认: 5)，如果指定了 --y-ticks 则忽略')
    parser.add_argument('--cv-nticks', type=int, default=5,
                        help='Cv轴刻度数量 (默认: 5)，如果指定了 --cv-ticks 则忽略')
    parser.add_argument('--partitions', '-p', type=str, default=None,
                        help='手动指定分区温度范围，格式: T1_min-T1_max,T2_min-T2_max，'
                             '例如: 200-700,750-950 表示第一分区200-700K，第二分区750-950K')
    parser.add_argument('--peak-method', type=str, default='fit',
                        choices=['data', 'partition', 'fit'],
                        help='热容峰计算方法: data=全点数据法, partition=分区点法, fit=拟合线外推法 (默认: fit)')
    parser.add_argument('--remove-outliers', action='store_true',
                        help='启用离群点剔除（默认IQR法）')
    parser.add_argument('--outlier-iqr', type=float, default=1.5,
                        help='IQR倍数阈值 (默认: 1.5)，越小剔除越严格')
    parser.add_argument('--outlier-method', type=str, default=None,
                        help='离群点剔除方法，格式: method:threshold:iterations\n'
                             '  iqr:1.5      - IQR法 (默认)\n'
                             '  zscore:2     - Z-score法，2个标准差\n'
                             '  zscore:1.5:3 - Z-score法，1.5σ，迭代3次\n'
                             '  mad:2.5      - MAD法（对异常值更鲁棒）\n'
                             '  percentile:5 - 百分位法，剔除最极端5%%的点')
    parser.add_argument('--exclude', nargs='+', metavar='TEMP:INDICES',
                        help='手动排除特定温度的数据点，格式: "300K:0,1" "600K:0"\n'
                             '索引按Lindemann指数从大到小排序（默认），0表示最大值\n'
                             '也可使用 --exclude-sort-by energy 改为按能量排序')
    parser.add_argument('--exclude-sort-by', type=str, default='delta',
                        choices=['delta', 'energy'],
                        help='排除点的排序依据: delta=Lindemann指数(默认), energy=能量')
    parser.add_argument('--use-lindemann-threshold', type=float, default=None,
                        metavar='THRESHOLD',
                        help='使用固定Lindemann阈值自动分区（替代手动--partitions）\n'
                             '例如: --use-lindemann-threshold 0.1 表示 δ<0.1为固相，δ≥0.1为液相\n'
                             '注意: 此参数会覆盖 --partitions，基于Lindemann指数自动确定分区边界')
    
    return parser.parse_args()


def main():
    """主函数"""
    args = parse_args()
    
    print("=" * 70)
    print("Step 6.1.1: 分区热容拟合图 - 论文出图专用")
    print("=" * 70)
    
    # 解析figsize
    try:
        fig_w, fig_h = map(float, args.figsize.lower().split('x'))
        figsize = (fig_w, fig_h)
        print(f"  图片尺寸: {figsize[0]}x{figsize[1]}")
    except ValueError:
        print(f"警告: 无效的figsize格式 '{args.figsize}'，使用默认 10x10")
        figsize = (10, 10)
    
    # 解析自定义刻度
    y_ticks_custom = None
    cv_ticks_custom = None
    
    if args.y_ticks:
        try:
            y_ticks_custom = [float(x.strip()) for x in args.y_ticks.split(',')]
            print(f"  能量Y轴刻度: {y_ticks_custom}")
        except ValueError:
            print(f"警告: 无效的 --y-ticks 格式 '{args.y_ticks}'，将自动计算")
    
    if args.cv_ticks:
        try:
            cv_ticks_custom = [float(x.strip()) for x in args.cv_ticks.split(',')]
            print(f"  Cv轴刻度: {cv_ticks_custom}")
        except ValueError:
            print(f"警告: 无效的 --cv-ticks 格式 '{args.cv_ticks}'，将自动计算")
    
    # 解析自定义分区
    custom_partitions = None
    if args.partitions:
        try:
            custom_partitions = []
            for part in args.partitions.split(','):
                T_min, T_max = map(float, part.strip().split('-'))
                custom_partitions.append((T_min, T_max))
            print(f"  自定义分区: {custom_partitions}")
        except ValueError:
            print(f"警告: 无效的 --partitions 格式 '{args.partitions}'，将使用聚类结果")
            print(f"  正确格式: T1_min-T1_max,T2_min-T2_max，例如 200-700,750-950")
            custom_partitions = None
    
    # 解析排除点参数
    exclude_dict = parse_exclude_points(args.exclude)
    if exclude_dict:
        print(f"\n  手动排除点配置:")
        print(f"    排序依据: {args.exclude_sort_by}")
        print(f"    排除规则: {exclude_dict}")
    
    # 列出可用结构
    if args.list:
        list_available_structures(method=args.clustering_method)
        return
    
    if args.structure is None and args.only_series is None:
        print("错误: 请指定 --structure 或 --only-series，或使用 --list 查看可用结构")
        return
    
    # 创建输出目录
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 获取可用结构
    available = find_clustering_results(method=args.clustering_method)
    
    # 根据参数确定要处理的结构列表
    if args.only_series:
        # 按系列筛选
        series_list = [s.strip().lower() for s in args.only_series.split(',')]
        structures = filter_structures_by_series(available, series_list)
        if not structures:
            print(f"\n错误: 系列 {series_list} 中没有找到任何结构")
            return
        print(f"\n处理系列 {series_list}: 共 {len(structures)} 个结构")
        for s in sorted(structures):
            print(f"  - {s}")
    elif args.structure.lower() == 'all':
        structures = list(available.keys())
        print(f"\n处理所有 {len(structures)} 个结构...")
    else:
        structures = [args.structure]
    
    # 处理每个结构
    results = []
    success = 0
    failed = 0
    
    for structure in structures:
        # 查找结构（大小写不敏感）
        found_name = None
        for name in available.keys():
            if name.lower() == structure.lower():
                found_name = name
                break
        
        if found_name is None:
            print(f"\n警告: 未找到结构 '{structure}'")
            failed += 1
            continue
        
        csv_path = available[found_name]
        df = load_cluster_data(csv_path, exclude_dict=exclude_dict, exclude_sort_by=args.exclude_sort_by)
        
        if df is None:
            failed += 1
            continue
        
        # 如果使用Lindemann阈值自动分区
        lindemann_partitions = None
        if args.use_lindemann_threshold is not None:
            if 'delta' not in df.columns:
                print(f"\n  警告: 数据中没有'delta'列，无法使用Lindemann阈值分区")
                print(f"        将使用聚类结果或手动分区")
            else:
                lindemann_partitions, _ = determine_partitions_by_lindemann(
                    df, threshold=args.use_lindemann_threshold
                )
                if lindemann_partitions:
                    print(f"\n  使用Lindemann阈值自动分区: {lindemann_partitions}")
                    # 覆盖手动分区
                    custom_partitions = lindemann_partitions
                    
                    # 🔥 关键修复: 根据Lindemann阈值重新分类phase_clustered
                    print(f"\n  [重新分类] 根据δ阈值={args.use_lindemann_threshold}重新分配数据点...")
                    df['phase_clustered'] = df['delta'].apply(
                        lambda x: 'partition1' if x < args.use_lindemann_threshold else 'partition2'
                    )
                    
                    # 统计重新分类结果
                    reclassify_stats = df.groupby(['temp', 'phase_clustered']).size().unstack(fill_value=0)
                    print(f"\n  重新分类后的分区分布:")
                    print(f"  {'温度(K)':<10} {'P1(固相)':<12} {'P2(液相)':<12}")
                    print(f"  {'-'*35}")
                    for temp in sorted(df['temp'].unique()):
                        p1_count = reclassify_stats.loc[temp, 'partition1'] if 'partition1' in reclassify_stats.columns else 0
                        p2_count = reclassify_stats.loc[temp, 'partition2'] if 'partition2' in reclassify_stats.columns else 0
                        print(f"  {temp:<10.0f} {p1_count:<12.0f} {p2_count:<12.0f}")
                else:
                    print(f"\n  Lindemann阈值分区失败，将使用聚类结果或手动分区")
        
        # 构建刻度参数
        tick_params = {
            'y_ticks_custom': y_ticks_custom,
            'cv_ticks_custom': cv_ticks_custom,
            'y_nticks': args.y_nticks,
            'cv_nticks': args.cv_nticks,
            'figsize': figsize,
        }
        
        # 确定离群点剔除参数
        if args.outlier_method:
            remove_outliers_param = args.outlier_method
        elif args.remove_outliers:
            remove_outliers_param = True
        else:
            remove_outliers_param = False
        
        result = plot_partition_cv(df, found_name, output_dir, 
                                   args.format, args.dpi, tick_params, custom_partitions,
                                   args.peak_method, remove_outliers_param, args.outlier_iqr)
        
        if result:
            results.append(result)
            success += 1
        else:
            failed += 1
    
    # 汇总
    print("\n" + "=" * 70)
    print(f"处理完成: 成功 {success}, 失败 {failed}")
    print(f"输出目录: {output_dir}")
    print("=" * 70)
    
    # 生成汇总表格
    if results:
        summary_file = output_dir / 'partition_cv_summary.csv'
        rows = []
        for r in results:
            row = {
                'structure': r['structure'],
                'Cv_overall': r['overall']['Cv'],
                'Cv_overall_err': r['overall']['Cv_err'],
                'R2_overall': r['overall']['R2'],
            }
            for i, (phase, fit) in enumerate(r['partitions'].items()):
                row[f'phase_{i+1}'] = phase
                row[f'Cv_{i+1}'] = fit['Cv']
                row[f'Cv_{i+1}_err'] = fit['Cv_err']
                row[f'R2_{i+1}'] = fit['R2']
            rows.append(row)
        
        df_summary = pd.DataFrame(rows)
        df_summary.to_csv(summary_file, index=False)
        print(f"汇总表已保存: {summary_file}")


if __name__ == '__main__':
    main()
