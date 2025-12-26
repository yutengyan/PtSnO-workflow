#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 6.1.1.7: 双重残差分析优化工具

针对不同分区使用不同的残差评估标准：
1. 能量残差: E_actual - E_fit (eV)
2. Lindemann残差: δ_actual - δ_fit (无量纲)

支持为每个分区设置独立的排除阈值，更合理地识别异常点。

策略:
- 固相分区 (低温): 主要看能量残差 + Lindemann残差
- 液相分区 (高温): 主要看Lindemann残差 + 能量残差
- 相变区: 两者权重相当

作者: AI Assistant
日期: 2025-12-25
"""

import os
import sys
import glob
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import linregress
from collections import defaultdict

# 设置中文字体
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'Microsoft YaHei']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.dpi'] = 150
plt.rcParams['savefig.dpi'] = 300

CV_SUPPORT = 38.2151  # meV/K


def find_clustering_file(structure_name, base_dir='results/step6_1_clustering', method='auto'):
    """
    查找指定结构的聚类结果文件
    
    Parameters:
    -----------
    structure_name : str
        结构名称
    base_dir : str
        聚类结果目录
    method : str
        聚类方法: 'auto', 'lindemann-threshold', 'kmeans'
        - 'auto': 优先使用lindemann-threshold,如果不存在则使用kmeans
        - 'lindemann-threshold': 只使用lindemann阈值方法
        - 'kmeans': 只使用kmeans方法
    
    Returns:
    --------
    str or None : 找到的文件路径,如果不存在则返回None
    """
    # 定义搜索顺序
    if method == 'auto':
        search_patterns = [
            f'{structure_name}_lindemann-threshold_n2_clustered_data.csv',
            f'{structure_name}_kmeans_n2_clustered_data.csv'
        ]
    elif method == 'lindemann-threshold':
        search_patterns = [f'{structure_name}_lindemann-threshold_n2_clustered_data.csv']
    elif method == 'kmeans':
        search_patterns = [f'{structure_name}_kmeans_n2_clustered_data.csv']
    else:
        raise ValueError(f"Unknown method: {method}. Use 'auto', 'lindemann-threshold', or 'kmeans'")
    
    # 尝试每种模式
    for pattern in search_patterns:
        file_path = os.path.join(base_dir, pattern)
        if os.path.exists(file_path):
            return file_path
    
    # 大小写不敏感匹配
    all_files = glob.glob(os.path.join(base_dir, '*_n2_clustered_data.csv'))
    for pattern in search_patterns:
        pattern_lower = pattern.lower()
        for f in all_files:
            if os.path.basename(f).lower() == pattern_lower:
                return f
    
    return None


def parse_partitions(partition_str):
    """解析分区字符串"""
    partitions = []
    for part in partition_str.split(','):
        T_min, T_max = map(float, part.strip().split('-'))
        partitions.append((T_min, T_max))
    return partitions


def parse_partition_thresholds(threshold_dict_str):
    """
    解析分区专属阈值
    
    格式: "1:energy=2.0,delta=1.5;2:energy=2.5,delta=2.0"
    """
    partition_thresholds = {}
    
    if not threshold_dict_str:
        return partition_thresholds
    
    for part_config in threshold_dict_str.split(';'):
        part_config = part_config.strip()
        if ':' not in part_config:
            continue
        
        partition_id, thresholds_str = part_config.split(':', 1)
        partition_id = int(partition_id.strip())
        
        thresholds = {}
        for threshold_pair in thresholds_str.split(','):
            threshold_pair = threshold_pair.strip()
            if '=' not in threshold_pair:
                continue
            
            key, value = threshold_pair.split('=')
            thresholds[key.strip()] = float(value.strip())
        
        partition_thresholds[partition_id] = thresholds
    
    return partition_thresholds


def filter_by_lindemann_threshold(df, partitions, lindemann_threshold, mode='strict'):
    """
    根据Lindemann阈值筛选分区内的点
    
    逻辑:
        - 分区1 (固相): 保留 δ < threshold 的点
        - 分区2 (液相): 保留 δ ≥ threshold 的点
    
    参数:
        df: 原始数据
        partitions: 分区温度范围列表
        lindemann_threshold: Lindemann阈值 (如0.1)
        mode: 'strict'=剔除不符合点, 'soft'=标记但保留
    
    返回:
        tuple: (filtered_df, filter_stats)
            filtered_df: 筛选后的数据
            filter_stats: 筛选统计信息
    """
    print("\n" + "="*80)
    print(f"Lindemann阈值筛选 (δ阈值={lindemann_threshold}, 模式={mode})")
    print("="*80)
    
    filter_stats = {
        'total_before': len(df),
        'total_after': 0,
        'removed_count': 0,
        'by_partition': {},
        'by_temp': {}
    }
    
    # 为每个数据点分配分区
    df = df.copy()
    df['assigned_partition'] = None
    
    for idx, row in df.iterrows():
        temp = row['temp']
        for partition_id, (T_min, T_max) in enumerate(partitions):
            if T_min <= temp <= T_max:
                df.at[idx, 'assigned_partition'] = partition_id + 1
                break
    
    # 标记是否符合Lindemann阈值条件
    df['lindemann_valid'] = False
    
    for idx, row in df.iterrows():
        partition = row['assigned_partition']
        delta = row['delta']
        
        if partition is None:
            continue
        
        # 分区1 (固相): δ < threshold
        if partition == 1:
            df.at[idx, 'lindemann_valid'] = (delta < lindemann_threshold)
        # 分区2 (液相): δ >= threshold
        elif partition == 2:
            df.at[idx, 'lindemann_valid'] = (delta >= lindemann_threshold)
        else:
            # 更多分区,需要扩展逻辑
            df.at[idx, 'lindemann_valid'] = True
    
    # 统计筛选结果
    print("\n  分区筛选统计:")
    print(f"  {'分区':<8} {'温度范围':<15} {'原始点数':<10} {'符合点数':<10} {'剔除点数':<10} {'保留率':<10}")
    print(f"  {'-'*70}")
    
    for partition_id, (T_min, T_max) in enumerate(partitions):
        partition_num = partition_id + 1
        partition_data = df[df['assigned_partition'] == partition_num]
        
        n_total = len(partition_data)
        n_valid = partition_data['lindemann_valid'].sum()
        n_removed = n_total - n_valid
        retention_rate = (n_valid / n_total * 100) if n_total > 0 else 0
        
        filter_stats['by_partition'][partition_num] = {
            'total': n_total,
            'valid': n_valid,
            'removed': n_removed,
            'retention_rate': retention_rate
        }
        
        print(f"  {partition_num:<8} {T_min:.0f}-{T_max:.0f}K{'':<6} "
              f"{n_total:<10} {n_valid:<10} {n_removed:<10} {retention_rate:>6.1f}%")
    
    # 按温度统计
    print("\n  温度级别筛选统计:")
    print(f"  {'温度(K)':<10} {'分区':<8} {'原始':<8} {'符合':<8} {'剔除':<8} {'示例剔除点δ值':<30}")
    print(f"  {'-'*80}")
    
    # 记录每个温度被剔除的原始索引
    lindemann_excluded_indices = {}
    
    # 同时记录按能量排序的索引 (用于生成--exclude命令)
    lindemann_excluded_indices_energy = {}
    
    for temp in sorted(df['temp'].unique()):
        temp_data = df[df['temp'] == temp]
        partition = temp_data['assigned_partition'].iloc[0] if len(temp_data) > 0 else None
        
        n_total = len(temp_data)
        n_valid = temp_data['lindemann_valid'].sum()
        n_removed = n_total - n_valid
        
        # 获取被剔除点的δ值示例
        removed_deltas = temp_data[~temp_data['lindemann_valid']]['delta'].values
        example_str = ', '.join([f'{d:.4f}' for d in removed_deltas[:3]])
        if len(removed_deltas) > 3:
            example_str += f' ... ({len(removed_deltas)}个)'
        
        # 获取被剔除点在该温度下的索引 (按delta排序后的索引)
        if n_removed > 0:
            temp_data_sorted = temp_data.sort_values('delta').reset_index(drop=True)
            removed_indices = temp_data_sorted[~temp_data_sorted['lindemann_valid']].index.tolist()
            lindemann_excluded_indices[temp] = removed_indices
            
            # 🔥 同时计算按能量排序的索引 (用于--exclude命令)
            temp_data_energy_sorted = temp_data.sort_values('avg_energy', ascending=False).reset_index(drop=True)
            removed_indices_energy = temp_data_energy_sorted[~temp_data_energy_sorted['lindemann_valid']].index.tolist()
            lindemann_excluded_indices_energy[temp] = removed_indices_energy
        
        filter_stats['by_temp'][temp] = {
            'total': n_total,
            'valid': n_valid,
            'removed': n_removed,
            'removed_deltas': removed_deltas.tolist(),
            'removed_indices': lindemann_excluded_indices.get(temp, []),
            'removed_indices_energy': lindemann_excluded_indices_energy.get(temp, [])
        }
        
        print(f"  {temp:<10.0f} {partition if partition else 'N/A':<8} "
              f"{n_total:<8} {n_valid:<8} {n_removed:<8} {example_str:<30}")
    
    # 保存Lindemann剔除的索引信息
    filter_stats['lindemann_excluded_indices'] = lindemann_excluded_indices
    filter_stats['lindemann_excluded_indices_energy'] = lindemann_excluded_indices_energy  # 按能量排序的索引
    
    # 根据模式决定是否实际剔除
    if mode == 'strict':
        filtered_df = df[df['lindemann_valid']].copy()
        filter_stats['total_after'] = len(filtered_df)
        filter_stats['removed_count'] = filter_stats['total_before'] - filter_stats['total_after']
        
        print(f"\n  [严格模式] 剔除不符合点: {filter_stats['removed_count']} 个")
        print(f"  剩余数据点: {filter_stats['total_after']} 个 "
              f"({filter_stats['total_after']/filter_stats['total_before']*100:.1f}%)")
    else:
        filtered_df = df.copy()
        filter_stats['total_after'] = len(filtered_df)
        filter_stats['removed_count'] = 0
        
        print(f"\n  [软模式] 标记但保留所有点")
    
    return filtered_df, filter_stats


def compute_dual_residuals(df, partitions):
    """
    计算双重残差：能量残差 + Lindemann残差
    
    返回:
        pd.DataFrame: 包含能量和Lindemann残差的数据
    """
    residual_data = []
    
    # 按温度分组
    temp_groups = df.groupby('temp')
    
    # 为每个分区计算拟合
    partition_fits_energy = {}
    partition_fits_delta = {}
    
    for partition_id, (T_min, T_max) in enumerate(partitions):
        # 收集该分区的数据
        partition_temps = []
        partition_E_cluster = []
        partition_delta = []
        
        for temp, group in temp_groups:
            if not (T_min <= temp <= T_max):
                continue
            
            # 计算团簇能量
            E_support = CV_SUPPORT / 1000 * temp
            E_cluster = group['avg_energy'].values - E_support
            delta_values = group['delta'].values
            
            partition_temps.append(temp)
            partition_E_cluster.append(np.mean(E_cluster))
            partition_delta.append(np.mean(delta_values))
        
        if len(partition_temps) < 2:
            continue
        
        partition_temps = np.array(partition_temps)
        partition_E_cluster = np.array(partition_E_cluster)
        partition_delta = np.array(partition_delta)
        
        # 拟合能量
        slope_E, intercept_E, r_value_E, _, _ = linregress(partition_temps, partition_E_cluster)
        
        # 拟合Lindemann
        slope_delta, intercept_delta, r_value_delta, _, _ = linregress(partition_temps, partition_delta)
        
        partition_fits_energy[partition_id] = {
            'slope': slope_E,
            'intercept': intercept_E,
            'r2': r_value_E ** 2,
            'Cv_cluster': slope_E * 1000
        }
        
        partition_fits_delta[partition_id] = {
            'slope': slope_delta,
            'intercept': intercept_delta,
            'r2': r_value_delta ** 2
        }
    
    # 计算每个数据点的残差
    for temp, group in temp_groups:
        # 确定分区
        partition_id = None
        for pid, (T_min, T_max) in enumerate(partitions):
            if T_min <= temp <= T_max:
                partition_id = pid
                break
        
        if partition_id is None:
            continue
        
        # 计算团簇能量
        E_support = CV_SUPPORT / 1000 * temp
        E_cluster = group['avg_energy'].values - E_support
        delta_values = group['delta'].values
        
        # 拟合值
        E_fit = partition_fits_energy[partition_id]['slope'] * temp + partition_fits_energy[partition_id]['intercept']
        delta_fit = partition_fits_delta[partition_id]['slope'] * temp + partition_fits_delta[partition_id]['intercept']
        
        # 每个点的残差
        for idx, (E_c, delta) in enumerate(zip(E_cluster, delta_values)):
            residual_E = E_c - E_fit
            residual_delta = delta - delta_fit
            
            residual_data.append({
                'temp': temp,
                'partition': partition_id + 1,
                'point_idx': idx,
                'energy_cluster': E_c,
                'energy_fit': E_fit,
                'residual_energy': residual_E,
                'delta': delta,
                'delta_fit': delta_fit,
                'residual_delta': residual_delta
            })
    
    residual_df = pd.DataFrame(residual_data)
    
    return residual_df, partition_fits_energy, partition_fits_delta


def identify_outliers_dual(residual_df, partition_thresholds, global_threshold_energy=2.0, 
                           global_threshold_delta=2.0, max_exclude=3):
    """
    基于双重残差识别异常点（分区独立标准差）
    
    参数:
        residual_df: 残差数据
        partition_thresholds: 分区专属阈值字典
        global_threshold_energy: 全局能量阈值 (σ单位)
        global_threshold_delta: 全局Lindemann阈值 (σ单位)
        max_exclude: 每个温度最多排除点数
    
    返回:
        tuple: (exclude_recommendations, outlier_details)
            exclude_recommendations: {temp: [indices_to_exclude]}
            outlier_details: {temp: {'energy': [...], 'delta': [...], 'both': [...]}}
    """
    exclude_recommendations = defaultdict(list)
    outlier_details = defaultdict(lambda: {'energy': [], 'delta': [], 'both': []})
    
    print("\n" + "="*80)
    print("双重残差异常点识别（分区独立标准）")
    print("="*80)
    
    # 为每个分区单独计算标准差
    partition_stats = {}
    for partition_id in sorted(residual_df['partition'].unique()):
        partition_data = residual_df[residual_df['partition'] == partition_id]
        
        std_energy = partition_data['residual_energy'].std()
        std_delta = partition_data['residual_delta'].std()
        
        partition_stats[partition_id] = {
            'std_energy': std_energy,
            'std_delta': std_delta,
            'n_points': len(partition_data)
        }
        
        print(f"\n  分区{partition_id} 残差统计:")
        print(f"    样本数: {len(partition_data)}")
        print(f"    能量残差标准差: {std_energy*1000:.3f} meV")
        print(f"    Lindemann残差标准差: {std_delta:.4f}")
    
    print("\n" + "-"*80)
    print("异常点识别结果:")
    print("-"*80)
    
    # 按温度分组
    for temp in sorted(residual_df['temp'].unique()):
        temp_data = residual_df[residual_df['temp'] == temp].copy()
        partition_id = temp_data['partition'].iloc[0]
        
        # 获取该分区的标准差
        std_energy = partition_stats[partition_id]['std_energy']
        std_delta = partition_stats[partition_id]['std_delta']
        
        # 获取该分区的阈值
        if partition_id in partition_thresholds:
            threshold_energy = partition_thresholds[partition_id].get('energy', global_threshold_energy)
            threshold_delta = partition_thresholds[partition_id].get('delta', global_threshold_delta)
        else:
            threshold_energy = global_threshold_energy
            threshold_delta = global_threshold_delta
        
        # 计算该温度的标准化残差（使用分区标准差）
        temp_data['z_score_energy'] = np.abs(temp_data['residual_energy']) / std_energy
        temp_data['z_score_delta'] = np.abs(temp_data['residual_delta']) / std_delta
        
        # 综合评分 (能量和Lindemann的最大Z分数)
        temp_data['combined_score'] = temp_data[['z_score_energy', 'z_score_delta']].max(axis=1)
        
        # 识别超过阈值的点
        outliers_energy = temp_data[temp_data['z_score_energy'] > threshold_energy]
        outliers_delta = temp_data[temp_data['z_score_delta'] > threshold_delta]
        
        # 分类记录异常点
        energy_indices = set(outliers_energy['point_idx'].tolist())
        delta_indices = set(outliers_delta['point_idx'].tolist())
        both_indices = energy_indices & delta_indices  # 交集
        only_energy_indices = energy_indices - delta_indices
        only_delta_indices = delta_indices - energy_indices
        
        # 合并异常点（去重）
        outlier_indices = energy_indices | delta_indices
        
        if outlier_indices:
            # 按综合评分排序，选择最异常的点
            temp_data_sorted = temp_data.sort_values('combined_score', ascending=False)
            
            exclude_list = []
            for _, row in temp_data_sorted.iterrows():
                if row['point_idx'] in outlier_indices and len(exclude_list) < max_exclude:
                    idx = int(row['point_idx'])
                    exclude_list.append(idx)
                    
                    # 分类记录
                    if idx in both_indices:
                        outlier_details[temp]['both'].append(idx)
                    elif idx in only_energy_indices:
                        outlier_details[temp]['energy'].append(idx)
                    elif idx in only_delta_indices:
                        outlier_details[temp]['delta'].append(idx)
            
            if exclude_list:
                exclude_recommendations[temp] = exclude_list
                
                print(f"\n  {temp}K (分区{partition_id}): 识别到 {len(exclude_list)} 个异常点")
                print(f"    分区标准差: E={std_energy*1000:.3f} meV, δ={std_delta:.4f}")
                print(f"    判定阈值: 能量={threshold_energy:.1f}σ, Lindemann={threshold_delta:.1f}σ")
                
                for idx in exclude_list[:3]:  # 显示前3个
                    row = temp_data[temp_data['point_idx'] == idx].iloc[0]
                    reason = ''
                    if idx in both_indices:
                        reason = '(双重异常)'
                    elif idx in only_energy_indices:
                        reason = '(仅能量)'
                    elif idx in only_delta_indices:
                        reason = '(仅Lindemann)'
                    
                    print(f"      点{idx}: E_res={row['residual_energy']*1000:+.2f} meV ({row['z_score_energy']:.2f}σ), "
                          f"δ_res={row['residual_delta']:+.4f} ({row['z_score_delta']:.2f}σ) {reason}")
    
    return dict(exclude_recommendations), dict(outlier_details)


def plot_dual_residual_analysis(residual_df, partition_fits_energy, partition_fits_delta,
                                exclude_dict, outlier_details, structure_name, output_dir):
    """
    绘制双重残差分析图（仅两个散点分布图）
    
    布局:
    - 左: 能量残差分布（只标记能量异常点，红X）
    - 右: Lindemann残差分布（只标记Lindemann异常点，红X）
    """
    fig, axes = plt.subplots(1, 2, figsize=(18, 7))
    fig.suptitle(f'{structure_name} - 双重残差分析（分区独立标准）', fontsize=16, fontweight='bold')
    
    # 颜色方案
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    
    # 计算每个分区的标准差
    partition_stats = {}
    for partition_id in sorted(residual_df['partition'].unique()):
        partition_data = residual_df[residual_df['partition'] == partition_id]
        partition_stats[partition_id] = {
            'std_energy': partition_data['residual_energy'].std(),
            'std_delta': partition_data['residual_delta'].std()
        }
    
    # ============================================================
    # 左: 能量残差分布（只标记能量异常点）
    # ============================================================
    ax1 = axes[0]
    
    # 绘制所有数据点
    for partition_id in sorted(residual_df['partition'].unique()):
        partition_data = residual_df[residual_df['partition'] == partition_id]
        
        ax1.scatter(partition_data['temp'], partition_data['residual_energy'] * 1000,
                   color=colors[(partition_id-1) % len(colors)], 
                   alpha=0.6, s=80, label=f'分区{partition_id}', edgecolors='k', linewidths=0.5)
    
    # 只标记能量异常的点（红X）
    n_energy_outliers = 0
    for temp, details in outlier_details.items():
        temp_data = residual_df[residual_df['temp'] == temp]
        
        # 双重异常点 + 仅能量异常点
        energy_outlier_indices = details['both'] + details['energy']
        
        for idx in energy_outlier_indices:
            point = temp_data[temp_data['point_idx'] == idx]
            if not point.empty:
                ax1.scatter(point['temp'], point['residual_energy'] * 1000,
                           color='red', marker='x', s=300, linewidths=3.5, zorder=5)
                n_energy_outliers += 1
    
    # 零线
    ax1.axhline(y=0, color='k', linestyle='--', linewidth=2, zorder=1)
    
    # 为每个分区绘制独立的±阈值线
    for partition_id, stats in partition_stats.items():
        std_E = stats['std_energy'] * 1000
        partition_data = residual_df[residual_df['partition'] == partition_id]
        temp_range = [partition_data['temp'].min(), partition_data['temp'].max()]
        
        color = colors[(partition_id-1) % len(colors)]
        
        # ±2σ线
        ax1.plot(temp_range, [2*std_E, 2*std_E], 
                linestyle=':', color=color, linewidth=2.5, alpha=0.8,
                label=f'分区{partition_id}: ±2σ={2*std_E:.1f} meV')
        ax1.plot(temp_range, [-2*std_E, -2*std_E], 
                linestyle=':', color=color, linewidth=2.5, alpha=0.8)
    
    ax1.set_xlabel('温度 (K)', fontsize=14, fontweight='bold')
    ax1.set_ylabel('能量残差 (meV)', fontsize=14, fontweight='bold')
    ax1.set_title(f'能量残差分布 (标记{n_energy_outliers}个能量异常点)', fontsize=14, fontweight='bold')
    ax1.legend(loc='best', fontsize=10, framealpha=0.9)
    ax1.grid(True, alpha=0.3, linestyle='--')
    
    # ============================================================
    # 右: Lindemann残差分布（只标记Lindemann异常点）
    # ============================================================
    ax2 = axes[1]
    
    # 绘制所有数据点
    for partition_id in sorted(residual_df['partition'].unique()):
        partition_data = residual_df[residual_df['partition'] == partition_id]
        
        ax2.scatter(partition_data['temp'], partition_data['residual_delta'],
                   color=colors[(partition_id-1) % len(colors)], 
                   alpha=0.6, s=80, label=f'分区{partition_id}', edgecolors='k', linewidths=0.5)
    
    # 只标记Lindemann异常的点（红X）
    n_delta_outliers = 0
    for temp, details in outlier_details.items():
        temp_data = residual_df[residual_df['temp'] == temp]
        
        # 双重异常点 + 仅Lindemann异常点
        delta_outlier_indices = details['both'] + details['delta']
        
        for idx in delta_outlier_indices:
            point = temp_data[temp_data['point_idx'] == idx]
            if not point.empty:
                ax2.scatter(point['temp'], point['residual_delta'],
                           color='red', marker='x', s=300, linewidths=3.5, zorder=5)
                n_delta_outliers += 1
    
    # 零线
    ax2.axhline(y=0, color='k', linestyle='--', linewidth=2, zorder=1)
    
    # 为每个分区绘制独立的±阈值线
    for partition_id, stats in partition_stats.items():
        std_delta = stats['std_delta']
        partition_data = residual_df[residual_df['partition'] == partition_id]
        temp_range = [partition_data['temp'].min(), partition_data['temp'].max()]
        
        color = colors[(partition_id-1) % len(colors)]
        
        # ±2σ线
        ax2.plot(temp_range, [2*std_delta, 2*std_delta], 
                linestyle=':', color=color, linewidth=2.5, alpha=0.8,
                label=f'分区{partition_id}: ±2σ={2*std_delta:.4f}')
        ax2.plot(temp_range, [-2*std_delta, -2*std_delta], 
                linestyle=':', color=color, linewidth=2.5, alpha=0.8)
    
    ax2.set_xlabel('温度 (K)', fontsize=14, fontweight='bold')
    ax2.set_ylabel('Lindemann残差', fontsize=14, fontweight='bold')
    ax2.set_title(f'Lindemann残差分布 (标记{n_delta_outliers}个Lindemann异常点)', fontsize=14, fontweight='bold')
    ax2.legend(loc='best', fontsize=10, framealpha=0.9)
    ax2.grid(True, alpha=0.3, linestyle='--')
    
    plt.tight_layout()
    
    output_file = output_dir / f'{structure_name}_dual_residual_analysis.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n  双重残差分析图已保存: {output_file}")
    print(f"    能量异常点: {n_energy_outliers}")
    print(f"    Lindemann异常点: {n_delta_outliers}")
    
    plt.close()


def generate_optimized_command(structure_name, partitions, exclude_dict, 
                               sort_by='energy', additional_args='', 
                               use_lindemann_threshold=None):
    """生成优化后的命令行"""
    exclude_args = []
    
    for temp in sorted(exclude_dict.keys()):
        indices = exclude_dict[temp]
        indices_str = ','.join(map(str, indices))
        exclude_args.append(f'"{int(temp)}K:{indices_str}"')
    
    partition_str = ','.join([f'{int(T_min)}-{int(T_max)}' for T_min, T_max in partitions])
    
    command = f"python step6_1_1_partition_cv_plot.py \\\n"
    command += f"    --structure {structure_name} \\\n"
    
    # 如果使用Lindemann阈值,优先使用该参数
    if use_lindemann_threshold is not None:
        command += f"    --use-lindemann-threshold {use_lindemann_threshold} \\\n"
    else:
        command += f"    --partitions {partition_str} \\\n"
    
    if exclude_args:
        command += f"    --exclude {' '.join(exclude_args)} \\\n"
        command += f"    --exclude-sort-by {sort_by}"
    
    if additional_args:
        command += f" \\\n    {additional_args}"
    
    return command


def main():
    parser = argparse.ArgumentParser(
        description='双重残差分析优化工具（能量+Lindemann）',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例用法:

  # 基础分析（使用全局阈值）
  python step6_1_1_7_dual_residual_analysis.py \\
      --structure Pt8sn6 \\
      --partitions 200-500,600-1100
  
  # 为不同分区设置不同阈值
  python step6_1_1_7_dual_residual_analysis.py \\
      --structure Pt8sn6 \\
      --partitions 200-500,600-1100 \\
      --partition-thresholds "1:energy=2.0,delta=1.5;2:energy=2.5,delta=2.0"
  
  # 使用Lindemann阈值筛选分区内的点（推荐）
  python step6_1_1_7_dual_residual_analysis.py \\
      --structure Pt8sn6 \\
      --partitions 200-500,600-1100 \\
      --use-lindemann-threshold 0.1 \\
      --lindemann-filter-mode strict
  
  # 完整分析：Lindemann筛选 + 分区阈值 + 生成命令
  python step6_1_1_7_dual_residual_analysis.py \\
      --structure Pt8sn6 \\
      --partitions 200-500,600-1100 \\
      --use-lindemann-threshold 0.1 \\
      --partition-thresholds "1:energy=1.5,delta=1.2;2:energy=2.0,delta=1.5" \\
      --max-exclude 10 \\
      --generate-command \\
      --additional-args "--y-ticks 0,2,4 --cv-ticks 3,4,5,6 --figsize 10x8 --peak-method partition"
"""
    )
    
    parser.add_argument('--structure', '-s', type=str, required=True,
                       help='结构名称（例如: Pt8sn6）')
    parser.add_argument('--partitions', '-p', type=str, required=True,
                       help='分区定义，格式: 200-500,600-1100')
    parser.add_argument('--clustering-method', type=str, default='auto',
                       choices=['auto', 'lindemann-threshold', 'kmeans'],
                       help='聚类方法选择: auto=优先lindemann-threshold, lindemann-threshold=仅lindemann, kmeans=仅kmeans (默认: auto)')
    
    # 阈值设置
    parser.add_argument('--threshold-energy', type=float, default=2.0,
                       help='全局能量残差阈值 (σ单位，默认: 2.0)')
    parser.add_argument('--threshold-delta', type=float, default=2.0,
                       help='全局Lindemann残差阈值 (σ单位，默认: 2.0)')
    parser.add_argument('--partition-thresholds', type=str, default='',
                       help='分区专属阈值，格式: "1:energy=2.0,delta=1.5;2:energy=2.5,delta=2.0"')
    
    # 排除设置
    parser.add_argument('--max-exclude', type=int, default=3,
                       help='每个温度最多排除点数 (默认: 3)')
    parser.add_argument('--sort-by', type=str, default='energy',
                       choices=['energy', 'delta'],
                       help='排除点排序依据 (默认: energy)')
    
    # Lindemann阈值筛选
    parser.add_argument('--use-lindemann-threshold', type=float, default=None,
                       help='使用Lindemann阈值筛选分区内的点（例如: 0.1表示分区1保留δ<0.1的点，分区2保留δ≥0.1的点）')
    parser.add_argument('--lindemann-filter-mode', type=str, default='strict',
                       choices=['strict', 'soft'],
                       help='Lindemann筛选模式: strict=严格剔除不符合点, soft=标记但保留 (默认: strict)')
    parser.add_argument('--exclude-lindemann-points', action='store_true',
                       help='在--exclude参数中包含Lindemann跨界点(而不是使用--use-lindemann-threshold重新分类)')
    
    # 输出设置
    parser.add_argument('--output-dir', '-o', type=str,
                       default='results/step6_1_1_dual_residual',
                       help='输出目录')
    parser.add_argument('--generate-command', action='store_true',
                       help='生成优化后的命令行')
    parser.add_argument('--additional-args', type=str, default='',
                       help='传递给step6_1_1_partition_cv_plot.py的额外参数')
    
    args = parser.parse_args()
    
    print("="*80)
    print("Step 6.1.1.7: 双重残差分析优化工具")
    print("="*80)
    
    # 创建输出目录
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 查找聚类文件
    csv_path = find_clustering_file(args.structure, method=args.clustering_method)
    if csv_path is None:
        print(f"\n错误: 未找到结构 '{args.structure}' 的聚类结果 (方法: {args.clustering_method})")
        print(f"  请确保已运行:")
        if args.clustering_method == 'lindemann-threshold':
            print(f"  python step6_1_clustering_analysis.py --structure {args.structure} --method lindemann-threshold --lindemann-threshold 0.1 --n-partitions 2")
        elif args.clustering_method == 'kmeans':
            print(f"  python step6_1_clustering_analysis.py --structure {args.structure} --n-partitions 2")
        else:
            print(f"  python step6_1_clustering_analysis.py --structure {args.structure} --method lindemann-threshold --lindemann-threshold 0.1 --n-partitions 2")
            print(f"  或")
            print(f"  python step6_1_clustering_analysis.py --structure {args.structure} --n-partitions 2")
        return
    
    # 检测使用的方法
    used_method = 'lindemann-threshold' if 'lindemann-threshold' in csv_path else 'kmeans'
    
    print(f"\n  结构: {args.structure}")
    print(f"  聚类方法: {used_method}")
    print(f"  数据文件: {csv_path}")
    print(f"  分区: {args.partitions}")
    
    # 加载数据
    df = pd.read_csv(csv_path)
    print(f"  原始数据: {len(df)} 条")
    
    # 解析分区
    partitions = parse_partitions(args.partitions)
    print(f"  分区数: {len(partitions)}")
    
    # 🔥 Lindemann阈值筛选（如果启用）
    lindemann_filtered = False
    lindemann_filter_stats = None
    if args.use_lindemann_threshold is not None:
        df, filter_stats = filter_by_lindemann_threshold(
            df, partitions, args.use_lindemann_threshold, mode=args.lindemann_filter_mode
        )
        lindemann_filtered = True
        lindemann_filter_stats = filter_stats
        
        # 保存筛选统计
        filter_summary = {
            'threshold': args.use_lindemann_threshold,
            'mode': args.lindemann_filter_mode,
            'total_before': filter_stats['total_before'],
            'total_after': filter_stats['total_after'],
            'removed_count': filter_stats['removed_count']
        }
        
        # 如果是strict模式,同步更新数据计数
        if args.lindemann_filter_mode == 'strict':
            print(f"\n  [Lindemann筛选后] 数据: {len(df)} 条 "
                  f"(剔除 {filter_stats['removed_count']} 条)")
    
    # 解析分区阈值
    partition_thresholds = parse_partition_thresholds(args.partition_thresholds)
    
    if partition_thresholds:
        print(f"\n  分区专属阈值:")
        for pid, thresholds in partition_thresholds.items():
            print(f"    分区{pid}: 能量={thresholds.get('energy', args.threshold_energy):.1f}σ, "
                  f"Lindemann={thresholds.get('delta', args.threshold_delta):.1f}σ")
    else:
        print(f"\n  全局阈值: 能量={args.threshold_energy:.1f}σ, Lindemann={args.threshold_delta:.1f}σ")
    
    # 计算双重残差
    print(f"\n计算双重残差...")
    residual_df, partition_fits_energy, partition_fits_delta = compute_dual_residuals(df, partitions)
    
    # 打印拟合统计
    print(f"\n分区拟合统计:")
    print(f"  {'分区':>6} {'温度范围':>12} {'Cv':>10} {'R2(能量)':>12} {'R2(δ)':>10}")
    print(f"  {'-'*60}")
    for pid in sorted(partition_fits_energy.keys()):
        T_min, T_max = partitions[pid]
        fit_E = partition_fits_energy[pid]
        fit_delta = partition_fits_delta[pid]
        print(f"  {pid+1:>6} {T_min:>6.0f}-{T_max:<6.0f} "
              f"{fit_E['Cv_cluster']:>10.2f} {fit_E['r2']:>12.4f} {fit_delta['r2']:>10.4f}")
    
    # 识别异常点
    exclude_dict, outlier_details = identify_outliers_dual(
        residual_df,
        partition_thresholds,
        global_threshold_energy=args.threshold_energy,
        global_threshold_delta=args.threshold_delta,
        max_exclude=args.max_exclude
    )
    
    # 保存残差数据
    output_csv = output_dir / f'{args.structure}_dual_residual_report.csv'
    residual_df.to_csv(output_csv, index=False)
    print(f"\n  残差数据已保存: {output_csv}")
    
    # 绘制分析图
    print(f"\n绘制双重残差分析图...")
    plot_dual_residual_analysis(
        residual_df, partition_fits_energy, partition_fits_delta,
        exclude_dict, outlier_details, args.structure, output_dir
    )
    
    # 打印排除建议
    print(f"\n" + "="*80)
    print("排除建议汇总")
    print("="*80)
    
    # 如果有Lindemann预筛选,先显示
    if lindemann_filtered and lindemann_filter_stats:
        print(f"\n  【第1步】Lindemann阈值预筛选 (δ={args.use_lindemann_threshold})")
        print(f"  {'-'*76}")
        
        total_lindemann_removed = lindemann_filter_stats['removed_count']
        if total_lindemann_removed > 0:
            print(f"  预筛选剔除: {total_lindemann_removed} 个点")
            print(f"\n  按温度统计:")
            
            for temp in sorted(lindemann_filter_stats['by_temp'].keys()):
                temp_stats = lindemann_filter_stats['by_temp'][temp]
                n_removed = temp_stats['removed']
                
                if n_removed > 0:
                    removed_deltas = temp_stats['removed_deltas']
                    delta_str = ', '.join([f'{d:.4f}' for d in removed_deltas[:5]])
                    if len(removed_deltas) > 5:
                        delta_str += f' ... (共{len(removed_deltas)}个)'
                    
                    print(f"    {int(temp)}K: 剔除 {n_removed} 个点, δ值=[{delta_str}]")
        else:
            print(f"  无需预筛选 (所有点都符合分区定义)")
    
    # 显示残差异常点
    if exclude_dict:
        print(f"\n  【第2步】双重残差异常点识别")
        print(f"  {'-'*76}")
        
        total_excluded = sum(len(indices) for indices in exclude_dict.values())
        print(f"  残差异常点: {total_excluded} 个")
        print(f"  涉及温度数: {len(exclude_dict)}")
        
        print(f"\n  按温度列表 (点索引基于筛选后数据):")
        for temp in sorted(exclude_dict.keys()):
            indices = exclude_dict[temp]
            indices_str = ','.join(map(str, indices))
            print(f"    {int(temp)}K: [{indices_str}]")
        
        # 汇总统计 - 仅显示识别到的点数,不预测最终剔除数量
        print(f"\n  【汇总统计】")
        print(f"  {'-'*76}")
        if lindemann_filtered and lindemann_filter_stats:
            print(f"  Lindemann预筛选: {lindemann_filter_stats['removed_count']} 个点")
            print(f"  残差异常点:     {total_excluded} 个点")
            print(f"  合计识别:       {lindemann_filter_stats['removed_count'] + total_excluded} 个点")
            print(f"\n  [注] 实际剔除点数取决于后续选择的处理模式:")
            print(f"    - --use-lindemann-threshold: 仅剔除残差异常点 ({total_excluded}个)")
            print(f"    - --exclude-lindemann-points: 剔除全部识别点 (合并去重后)")
        else:
            print(f"  总计识别: {total_excluded} 个点")
    else:
        if not lindemann_filtered or lindemann_filter_stats['removed_count'] == 0:
            print(f"\n未识别到需要排除的异常点。")
        else:
            print(f"\n  【第2步】双重残差异常点识别")
            print(f"  {'-'*76}")
            print(f"  未识别到残差异常点 (Lindemann预筛选已充分清理数据)")
    
    # 生成命令
    if args.generate_command:
        print(f"\n" + "="*80)
        print("优化后的命令")
        print("="*80)
        
        # 整合排除索引
        combined_exclude_dict = {}
        use_threshold_in_command = False
        
        if lindemann_filtered and lindemann_filter_stats:
            if args.exclude_lindemann_points:
                # 🔥 模式1: 在--exclude中包含所有点(Lindemann+残差),不使用--use-lindemann-threshold
                # 直接合并Lindemann剔除索引和残差异常索引
                all_temps = set()
                if exclude_dict:
                    all_temps.update(exclude_dict.keys())
                if 'lindemann_excluded_indices_energy' in lindemann_filter_stats:
                    all_temps.update(lindemann_filter_stats['lindemann_excluded_indices_energy'].keys())
                
                for temp in all_temps:
                    # 🔥 使用按能量排序的Lindemann索引
                    lindemann_excluded = lindemann_filter_stats.get('lindemann_excluded_indices_energy', {}).get(temp, [])
                    residual_excluded = exclude_dict.get(temp, [])
                    
                    # 将残差索引映射回原始索引
                    # 注意:这里使用按δ排序的Lindemann索引来计算偏移
                    lindemann_excluded_delta_sorted = sorted(
                        lindemann_filter_stats.get('lindemann_excluded_indices', {}).get(temp, [])
                    )
                    original_residual_indices = []
                    for res_idx in residual_excluded:
                        n_excluded_before = sum(1 for ex_idx in lindemann_excluded_delta_sorted if ex_idx <= res_idx)
                        original_idx = res_idx + n_excluded_before
                        original_residual_indices.append(original_idx)
                    
                    # 合并并去重 (现在都是按能量排序的索引)
                    all_excluded = sorted(set(lindemann_excluded + original_residual_indices))
                    if all_excluded:
                        combined_exclude_dict[temp] = all_excluded
                
                print(f"\n  [说明]")
                total_combined = sum(len(v) for v in combined_exclude_dict.values())
                print(f"  --exclude参数包含全部{total_combined}个剔除点:")
                
                # 统计有多少Lindemann点和残差点
                lindemann_only = 0
                residual_only = 0
                overlap = 0
                
                for temp in sorted(set(list(lindemann_filter_stats.get('lindemann_excluded_indices_energy', {}).keys()) + 
                                      list(exclude_dict.keys()))):
                    lind_set = set(lindemann_filter_stats.get('lindemann_excluded_indices_energy', {}).get(temp, []))
                    resid_set = set(exclude_dict.get(temp, []))
                    
                    # 映射残差索引到原始索引
                    lindemann_excluded_delta_sorted = sorted(
                        lindemann_filter_stats.get('lindemann_excluded_indices', {}).get(temp, [])
                    )
                    resid_original = set()
                    for res_idx in resid_set:
                        n_excluded_before = sum(1 for ex_idx in lindemann_excluded_delta_sorted if ex_idx <= res_idx)
                        original_idx = res_idx + n_excluded_before
                        resid_original.add(original_idx)
                    
                    overlap += len(lind_set & resid_original)
                    lindemann_only += len(lind_set - resid_original)
                    residual_only += len(resid_original - lind_set)
                
                print(f"    - 仅Lindemann跨界点: {lindemann_only} 个")
                print(f"    - 仅残差异常点: {residual_only} 个")
                print(f"    - 两者重叠: {overlap} 个")
                print(f"    - 全部删除这些点,无需--use-lindemann-threshold")
                use_threshold_in_command = False
                
            else:
                # 🔥 模式2: 使用--use-lindemann-threshold重新分类,--exclude只包含残差异常点
                # step6_1_1_partition_cv_plot.py只重新分类,不删除数据
                # 所以--exclude应该只包含残差异常点,索引需要映射回原始数据
                
                if exclude_dict:
                    for temp, residual_indices in exclude_dict.items():
                        # 获取该温度被Lindemann剔除的索引(已排序)
                        lindemann_excluded_at_temp = sorted(lindemann_filter_stats.get('lindemann_excluded_indices', {}).get(temp, []))
                        
                        # 将残差异常点的索引(基于筛选后数据)映射回原始索引
                        original_indices = []
                        for res_idx in residual_indices:
                            # 计算原始索引:加上之前被Lindemann剔除的点数
                            n_excluded_before = sum(1 for ex_idx in lindemann_excluded_at_temp if ex_idx <= res_idx)
                            original_idx = res_idx + n_excluded_before
                            original_indices.append(original_idx)
                        
                        combined_exclude_dict[temp] = sorted(original_indices)
                
                print(f"\n  [说明]")
                print(f"  --use-lindemann-threshold {args.use_lindemann_threshold}:")
                print(f"    - 不删除数据,只重新分类phase_clustered标签")
                print(f"    - 确保分区边界严格按照δ阈值划分")
                print(f"  --exclude参数:")
                print(f"    - 仅包含{sum(len(v) for v in combined_exclude_dict.values())}个残差异常点")
                print(f"    - 索引基于原始100点数据")
                print(f"    - Lindemann预筛选的{lindemann_filter_stats.get('removed_count', 0)}个点通过重新分类自动处理")
                use_threshold_in_command = True
        else:
            # 没有Lindemann筛选,直接使用残差索引
            combined_exclude_dict = exclude_dict.copy() if exclude_dict else {}
            use_threshold_in_command = False
        
        command = generate_optimized_command(
            args.structure, partitions, combined_exclude_dict,
            sort_by=args.sort_by, additional_args=args.additional_args,
            use_lindemann_threshold=args.use_lindemann_threshold if use_threshold_in_command else None
        )
        
        print(f"\n{command}\n")
        
        # 打印排除点统计
        if combined_exclude_dict or (lindemann_filtered and lindemann_filter_stats):
            total_exclude = sum(len(indices) for indices in combined_exclude_dict.values())
            print(f"  [提示] 数据清理统计:")
            if lindemann_filtered and lindemann_filter_stats and lindemann_filter_stats['removed_count'] > 0:
                lindemann_count = lindemann_filter_stats['removed_count']
                
                # 判断是哪种模式
                if args.exclude_lindemann_points:
                    # --exclude-lindemann-points模式: total_exclude已包含Lindemann点
                    print(f"    - --exclude参数包含: {total_exclude} 个点")
                    print(f"    - 最终保留: {100 - total_exclude} 个点")
                else:
                    # --use-lindemann-threshold模式: Lindemann点通过重新分类处理
                    print(f"    - Lindemann预筛选: {lindemann_count} 个点 (通过重新分类处理)")
                    print(f"    - 残差异常点(--exclude): {total_exclude} 个点")
                    print(f"    - 总计清理: {lindemann_count + total_exclude} 个点")
                    print(f"    - 最终保留: {100 - lindemann_count - total_exclude} 个点")
            else:
                print(f"    - 残差异常点(--exclude): {total_exclude} 个点")
                print(f"    - 最终保留: {100 - total_exclude} 个点")
    
    print("\n完成!")


if __name__ == '__main__':
    main()
