#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
Step 5.1: 熔点分析 - 负载型纳米团簇熔点与组成关系
================================================================================

作者: GitHub Copilot
日期: 2025-12-02
版本: v1.0

功能概述
========
分析负载型 Pt-Sn 纳米团簇的熔点，探究熔点与以下因素的关系：
1. Pt/Sn 比例
2. 总原子数
3. 是否含氧

熔点定义方法
============
本脚本采用两种方法确定熔点，并进行对比：

**方法A: 林德曼指数跃变法 (Lindemann Index Jump)**
- 找到 δ 从 <0.1 跃变到 ≥0.1 的温度
- 使用线性插值确定精确熔点
- 物理意义：原子振动幅度达到临界值

**方法B: 聚类分区边界法 (Clustering Partition Boundary)**
- 从 K-means 聚类的 2 分区结果中提取分区边界
- 分区边界代表低温相→高温相的转变区域中心
- 需要结合林德曼阈值确定边界温度

数据来源
========
1. step6_0 主数据: results/step6_0_multi_system/step6_0_all_systems_data.csv
   - 包含温度、Lindemann-δ、能量等原始数据
   
2. step6_1 聚类数据: results/step6_1_clustering/*_kmeans_n2_clustered_data.csv
   - 包含聚类分区结果

输出
====
1. results/step5_1_melting_point/melting_point_summary.csv
   - 熔点汇总表（两种方法）
   
2. results/step5_1_melting_point/Tm_vs_composition.png
   - 熔点 vs 组成关系图
   
3. results/step5_1_melting_point/MELTING_POINT_ANALYSIS_REPORT.md
   - 分析报告

================================================================================
Usage (use vscode-1 environment):
================================================================================

# Run full analysis
& C:/Users/11207/.conda/envs/vscode-1/python.exe step5_1_melting_point_analysis.py

# Exclude specific structures (e.g., exclude Pt8sn0 and Pt6)
& C:/Users/11207/.conda/envs/vscode-1/python.exe step5_1_melting_point_analysis.py --exclude Pt8sn0 Pt6

================================================================================
"""

import os
import re
import glob
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import stats
from scipy.interpolate import interp1d
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# 配置
# ============================================================================

BASE_DIR = Path(__file__).parent

# 输入数据
DATA_FILE = BASE_DIR / 'results' / 'step6_0_multi_system' / 'step6_0_all_systems_data.csv'
CLUSTERING_DIR = BASE_DIR / 'results' / 'step6_1_clustering'

# 输出目录
OUTPUT_DIR = BASE_DIR / 'results' / 'step5_1_melting_point'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# 英文字体设置（适合期刊发表）
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.unicode_minus'] = False

print("=" * 70)
print("Step 5.1: Melting Point Analysis")
print("=" * 70)

# ============================================================================
# 命令行参数解析
# ============================================================================

def parse_args():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description='Melting Point Analysis for Pt-Sn Nanoclusters',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        '--exclude', '-e',
        nargs='+',
        default=[],
        help='List of structures to exclude from series analysis (e.g., --exclude Pt8sn0 Pt6)'
    )
    
    return parser.parse_args()

# 解析命令行参数
args = parse_args()


# ============================================================================
# 1. 组成解析函数
# ============================================================================

def parse_composition(name):
    """
    从结构名提取 (Pt数, Sn数, O数)
    
    支持格式：
    - Pt6sn8, pt8sn5, PT3SN5
    - Pt6sn8o4, O2pt4sn6, Sn6pt5o2
    - Air68, Air86
    - Cv, cv-1~5 -> Pt6Sn8O4
    """
    name_lower = name.lower().strip()
    
    # Cv 系列: 实际对应 Pt6Sn8O4
    if name_lower in ['cv', 'cv-1', 'cv-2', 'cv-3', 'cv-4', 'cv-5']:
        return (6, 8, 4)  # Pt6Sn8O4
    
    # Air 系列: air68 -> Pt6Sn8, air86 -> Pt8Sn6
    if name_lower.startswith('air'):
        digits = re.findall(r'\d+', name_lower)
        if digits:
            num = digits[0]
            if len(num) == 2:
                return (int(num[0]), int(num[1]), 0)
        return None
    
    # 通用解析: 提取 pt, sn, o 数量
    pt_match = re.search(r'pt\s*(\d+)', name_lower)
    sn_match = re.search(r'sn\s*(\d+)', name_lower)
    o_match = re.search(r'o\s*(\d+)', name_lower)
    
    pt = int(pt_match.group(1)) if pt_match else 0
    sn = int(sn_match.group(1)) if sn_match else 0
    o = int(o_match.group(1)) if o_match else 0
    
    if pt == 0 and sn == 0:
        return None
    
    return (pt, sn, o)


def classify_structure(name):
    """
    分类结构类型
    
    Returns:
        tuple: (type_key, type_label)
        - type_key: 'air', 'supported', 'oxide'
        - type_label: 显示标签
    """
    comp = parse_composition(name)
    if comp is None:
        return ('unknown', 'Unknown')
    
    pt, sn, o = comp
    name_lower = name.lower()
    
    # 气相 (Air)
    if name_lower.startswith('air'):
        return ('air', f'Gas-phase')
    
    # 含氧负载型
    if o > 0:
        return ('oxide', f'Supported Pt-Sn-O')
    
    # 负载型 (无氧)
    return ('supported', f'Supported Pt-Sn')


# ============================================================================
# 2. 林德曼指数熔点分析 (方法A)
# ============================================================================

def calculate_melting_point_lindemann(df_structure, threshold=0.1):
    """
    使用林德曼指数跃变法计算熔点
    
    Args:
        df_structure: 单个结构的数据
        threshold: 林德曼阈值 (默认 0.1)
    
    Returns:
        dict: {
            'Tm_lindemann': 熔点温度 (K),
            'Tm_lindemann_err': 熔点误差估计 (K),
            'method': 'lindemann',
            'delta_at_Tm': 熔点处的 δ 值
        }
    """
    # 按温度分组计算平均 Lindemann 指数
    df_avg = df_structure.groupby('temp').agg({
        'delta': ['mean', 'std', 'count']
    }).reset_index()
    df_avg.columns = ['temp', 'delta_mean', 'delta_std', 'count']
    df_avg = df_avg.sort_values('temp')
    
    temps = df_avg['temp'].values
    deltas = df_avg['delta_mean'].values
    
    if len(temps) < 2:
        return None
    
    # 找到 δ 从 < threshold 跃变到 >= threshold 的位置
    Tm = None
    for i in range(len(temps) - 1):
        if deltas[i] < threshold and deltas[i+1] >= threshold:
            # 线性插值
            T1, T2 = temps[i], temps[i+1]
            d1, d2 = deltas[i], deltas[i+1]
            
            if d2 != d1:
                Tm = T1 + (threshold - d1) * (T2 - T1) / (d2 - d1)
            else:
                Tm = (T1 + T2) / 2
            
            # 误差估计 (取温度步长的一半)
            Tm_err = (T2 - T1) / 2
            
            return {
                'Tm_lindemann': Tm,
                'Tm_lindemann_err': Tm_err,
                'method': 'lindemann',
                'delta_below': d1,
                'delta_above': d2,
                'T_below': T1,
                'T_above': T2
            }
    
    # 如果全程 δ < threshold (未熔化)
    if all(d < threshold for d in deltas):
        return {
            'Tm_lindemann': np.nan,
            'Tm_lindemann_err': np.nan,
            'method': 'lindemann',
            'note': 'Always solid (δ < 0.1)'
        }
    
    # 如果全程 δ >= threshold (低温已熔化)
    if all(d >= threshold for d in deltas):
        return {
            'Tm_lindemann': temps[0],  # 最低温度
            'Tm_lindemann_err': np.nan,
            'method': 'lindemann',
            'note': f'Melted at lowest T ({temps[0]} K)'
        }
    
    return None


# ============================================================================
# 2.5 林德曼指数跃变点分析 (方法A2: dδ/dT 最大值)
# ============================================================================

def calculate_melting_point_transition(df_structure):
    """
    使用林德曼指数跃变点法计算熔点 (dδ/dT 最大值)
    
    物理意义：找到 δ(T) 曲线斜率最大的点，即熔化转变最剧烈的温度
    
    Args:
        df_structure: 单个结构的数据
    
    Returns:
        dict: {
            'Tm_transition': 跃变点温度 (K),
            'Tm_transition_err': 温度误差估计 (K),
            'dDelta_dT_max': 最大斜率值,
            'delta_at_transition': 跃变点处的 δ 值
        }
    """
    from scipy.ndimage import gaussian_filter1d
    
    # 按温度分组计算平均 Lindemann 指数
    df_avg = df_structure.groupby('temp').agg({
        'delta': ['mean', 'std']
    }).reset_index()
    df_avg.columns = ['temp', 'delta_mean', 'delta_std']
    df_avg = df_avg.sort_values('temp')
    
    temps = df_avg['temp'].values
    deltas = df_avg['delta_mean'].values
    
    if len(temps) < 3:
        return None
    
    # 平滑处理（高斯滤波），避免噪声影响
    deltas_smooth = gaussian_filter1d(deltas, sigma=1)
    
    # 计算 dδ/dT (中心差分)
    dT = np.diff(temps)
    dDelta = np.diff(deltas_smooth)
    dDelta_dT = dDelta / dT
    
    # 对应的温度点（取中点）
    T_mid = (temps[:-1] + temps[1:]) / 2
    
    if len(dDelta_dT) == 0:
        return None
    
    # 找最大斜率点
    max_idx = np.argmax(dDelta_dT)
    Tm_transition = T_mid[max_idx]
    dDelta_dT_max = dDelta_dT[max_idx]
    
    # 估算误差（取相邻温度步长）
    Tm_err = dT[max_idx] / 2 if max_idx < len(dT) else dT[-1] / 2
    
    # 插值获取跃变点处的 δ 值
    f_interp = interp1d(temps, deltas_smooth, kind='linear', fill_value='extrapolate')
    delta_at_transition = float(f_interp(Tm_transition))
    
    # 检查是否有明显的跃变（斜率阈值）
    if dDelta_dT_max < 0.0001:  # 斜率太小，无明显跃变
        return {
            'Tm_transition': np.nan,
            'Tm_transition_err': np.nan,
            'dDelta_dT_max': dDelta_dT_max,
            'delta_at_transition': np.nan,
            'note': 'No significant transition detected'
        }
    
    return {
        'Tm_transition': Tm_transition,
        'Tm_transition_err': Tm_err,
        'dDelta_dT_max': dDelta_dT_max,
        'delta_at_transition': delta_at_transition,
        'method': 'transition'
    }


# ============================================================================
# 3. 聚类分区边界熔点分析 (方法B)
# ============================================================================

def calculate_melting_point_clustering(structure_name):
    """
    从聚类分区结果提取相变温度
    
    使用 K-means n=2 分区结果，找到分区边界
    
    Args:
        structure_name: 结构名称
    
    Returns:
        dict: {
            'Tm_clustering': 相变温度 (K),
            'Tm_clustering_err': 误差估计 (K),
            'method': 'clustering',
            'partition1_Tmax': 分区1最高温度,
            'partition2_Tmin': 分区2最低温度
        }
    """
    # 读取聚类数据
    clustered_file = CLUSTERING_DIR / f'{structure_name}_kmeans_n2_clustered_data.csv'
    
    if not clustered_file.exists():
        return None
    
    try:
        df = pd.read_csv(clustered_file)
    except Exception as e:
        print(f"  [WARNING] Failed to read {clustered_file}: {e}")
        return None
    
    if 'phase_clustered' not in df.columns:
        return None
    
    # 分区1 和 分区2 的温度范围
    df_p1 = df[df['phase_clustered'] == 'partition1']
    df_p2 = df[df['phase_clustered'] == 'partition2']
    
    if df_p1.empty or df_p2.empty:
        return None
    
    # 分区边界 = partition1 最高温度 和 partition2 最低温度 的中点
    p1_Tmax = df_p1['temp'].max()
    p2_Tmin = df_p2['temp'].min()
    
    # 相变温度 = 边界中点
    Tm = (p1_Tmax + p2_Tmin) / 2
    Tm_err = abs(p2_Tmin - p1_Tmax) / 2
    
    # 同时获取分区的 δ 均值
    delta_p1 = df_p1['delta'].mean()
    delta_p2 = df_p2['delta'].mean()
    
    return {
        'Tm_clustering': Tm,
        'Tm_clustering_err': Tm_err,
        'method': 'clustering',
        'partition1_Tmax': p1_Tmax,
        'partition2_Tmin': p2_Tmin,
        'delta_partition1': delta_p1,
        'delta_partition2': delta_p2
    }


# ============================================================================
# 4. 主分析函数
# ============================================================================

def analyze_all_melting_points():
    """
    分析所有结构的熔点
    """
    print(f"\n[1] Loading data...")
    
    if not DATA_FILE.exists():
        print(f"  [ERROR] Data file not found: {DATA_FILE}")
        return None
    
    df = pd.read_csv(DATA_FILE)
    print(f"  Loaded {len(df)} records from {len(df['structure'].unique())} structures")
    
    results = []
    
    structures = sorted(df['structure'].unique())
    
    print(f"\n[2] Analyzing melting points...")
    
    for structure in structures:
        df_struct = df[df['structure'] == structure]
        
        # 解析组成
        comp = parse_composition(structure)
        if comp is None:
            continue
        
        pt, sn, o = comp
        total_atoms = pt + sn + o
        pt_sn_ratio = pt / sn if sn > 0 else float('inf')
        # Sn fraction = Sn / (Pt + Sn + O)，考虑整个团簇的组成
        sn_fraction = sn / total_atoms if total_atoms > 0 else 0
        
        # 分类
        type_key, type_label = classify_structure(structure)
        
        # 方法A1: 林德曼指数阈值法 (δ = 0.1)
        result_lindemann = calculate_melting_point_lindemann(df_struct)
        
        # 方法A2: 林德曼指数跃变点法 (dδ/dT 最大值)
        result_transition = calculate_melting_point_transition(df_struct)
        
        # 方法B: 聚类分区
        result_clustering = calculate_melting_point_clustering(structure)
        
        # 汇总
        row = {
            'structure': structure,
            'type': type_key,
            'type_label': type_label,
            'Pt': pt,
            'Sn': sn,
            'O': o,
            'total_atoms': total_atoms,
            'Pt_Sn_ratio': pt_sn_ratio,
            'Sn_fraction': sn_fraction,
            'has_oxygen': o > 0,
        }
        
        # 林德曼法结果
        if result_lindemann:
            row['Tm_lindemann'] = result_lindemann.get('Tm_lindemann', np.nan)
            row['Tm_lindemann_err'] = result_lindemann.get('Tm_lindemann_err', np.nan)
            row['lindemann_note'] = result_lindemann.get('note', '')
        else:
            row['Tm_lindemann'] = np.nan
            row['Tm_lindemann_err'] = np.nan
            row['lindemann_note'] = 'No transition found'
        
        # 跃变点法结果 (dδ/dT 最大值)
        if result_transition:
            row['Tm_transition'] = result_transition.get('Tm_transition', np.nan)
            row['Tm_transition_err'] = result_transition.get('Tm_transition_err', np.nan)
            row['dDelta_dT_max'] = result_transition.get('dDelta_dT_max', np.nan)
            row['delta_at_transition'] = result_transition.get('delta_at_transition', np.nan)
            row['transition_note'] = result_transition.get('note', '')
        else:
            row['Tm_transition'] = np.nan
            row['Tm_transition_err'] = np.nan
            row['dDelta_dT_max'] = np.nan
            row['delta_at_transition'] = np.nan
            row['transition_note'] = 'No transition found'
        
        # 聚类法结果
        if result_clustering:
            row['Tm_clustering'] = result_clustering.get('Tm_clustering', np.nan)
            row['Tm_clustering_err'] = result_clustering.get('Tm_clustering_err', np.nan)
            row['delta_partition1'] = result_clustering.get('delta_partition1', np.nan)
            row['delta_partition2'] = result_clustering.get('delta_partition2', np.nan)
        else:
            row['Tm_clustering'] = np.nan
            row['Tm_clustering_err'] = np.nan
            row['delta_partition1'] = np.nan
            row['delta_partition2'] = np.nan
        
        # 方法间差异计算
        # Tm_lindemann vs Tm_clustering
        if not np.isnan(row['Tm_lindemann']) and not np.isnan(row['Tm_clustering']):
            row['Tm_diff_lind_clust'] = row['Tm_lindemann'] - row['Tm_clustering']
        else:
            row['Tm_diff_lind_clust'] = np.nan
        
        # Tm_lindemann vs Tm_transition (阈值法 vs 跃变点法)
        if not np.isnan(row['Tm_lindemann']) and not np.isnan(row['Tm_transition']):
            row['Tm_diff_lind_trans'] = row['Tm_lindemann'] - row['Tm_transition']
        else:
            row['Tm_diff_lind_trans'] = np.nan
        
        results.append(row)
    
    df_results = pd.DataFrame(results)
    
    # 保存结果
    output_csv = OUTPUT_DIR / 'melting_point_summary.csv'
    df_results.to_csv(output_csv, index=False)
    print(f"\n  [OK] Saved: {output_csv}")
    
    return df_results
    
    return df_results


# ============================================================================
# 5. 可视化
# ============================================================================

def plot_melting_point_analysis(df_results):
    """
    生成熔点分析可视化
    """
    print(f"\n[3] Generating visualizations...")
    
    # 过滤有效数据
    df_valid = df_results[~df_results['Tm_lindemann'].isna()].copy()
    
    if df_valid.empty:
        print("  [WARNING] No valid melting point data")
        return
    
    # ========================================================================
    # 图1: 熔点 vs Sn 含量 (按类型分组)
    # ========================================================================
    fig, axes = plt.subplots(3, 2, figsize=(14, 16))
    
    # 获取屏蔽列表
    exclude_list = [s.lower() for s in args.exclude]
    
    # 样式设置 - 使用中英文混合标签
    colors = {
        'air': '#1E90FF',      # 蓝色 - 气相
        'supported': '#2E8B57', # 绿色 - 负载型
        'oxide': '#FF6347'      # 红色 - 含氧
    }
    markers = {
        'air': 'o',
        'supported': 's', 
        'oxide': '^'
    }
    labels = {
        'air': 'Gas-phase Alloy',
        'supported': 'Supported Alloy',
        'oxide': 'Supported Oxide'
    }
    
    # --- (a) Tm (Lindemann) vs Sn fraction - 仅负载型无氧合金 ---
    ax = axes[0, 0]
    df_supported = df_valid[
        (df_valid['type'] == 'supported') & 
        (df_valid['has_oxygen'] == False) &
        (~df_valid['structure'].str.lower().isin(exclude_list))  # 应用过滤
    ]
    
    if not df_supported.empty:
        ax.scatter(df_supported['Sn_fraction'], df_supported['Tm_lindemann'],
                  c='#2E8B57', marker='s', s=80, alpha=0.8,
                  label='Supported Alloy')
        
        # 线性拟合
        if len(df_supported) > 2:
            slope, intercept, r, p, se = stats.linregress(
                df_supported['Sn_fraction'], df_supported['Tm_lindemann']
            )
            x_fit = np.linspace(df_supported['Sn_fraction'].min(), 
                               df_supported['Sn_fraction'].max(), 100)
            y_fit = slope * x_fit + intercept
            ax.plot(x_fit, y_fit, '--', color='gray', linewidth=2,
                   label=f'Fit: {slope:.0f}x + {intercept:.0f}\n$R^2$={r**2:.2f}, p={p:.4f}')
    
    ax.set_xlabel('Sn Fraction = Sn/(Pt+Sn+O)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Melting Point $T_m$ (K)', fontsize=11, fontweight='bold')
    ax.set_title('(a) $T_m$ vs Sn Fraction (Supported Alloy)', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    
    # --- (b) Tm (Lindemann) vs total atoms - 仅负载型无氧合金 ---
    ax = axes[0, 1]
    df_supported_b = df_valid[
        (df_valid['type'] == 'supported') & 
        (df_valid['has_oxygen'] == False) &
        (~df_valid['structure'].str.lower().isin(exclude_list))  # 应用过滤
    ]
    
    if not df_supported_b.empty:
        ax.scatter(df_supported_b['total_atoms'], df_supported_b['Tm_lindemann'],
                  c='#2E8B57', marker='s', s=80, alpha=0.8,
                  label='Supported Alloy')
        
        # 线性拟合
        if len(df_supported_b) > 2:
            slope, intercept, r, p, se = stats.linregress(
                df_supported_b['total_atoms'], df_supported_b['Tm_lindemann']
            )
            x_fit = np.linspace(df_supported_b['total_atoms'].min(), 
                               df_supported_b['total_atoms'].max(), 100)
            y_fit = slope * x_fit + intercept
            ax.plot(x_fit, y_fit, '--', color='gray', linewidth=2,
                   label=f'Fit: {slope:.1f}x + {intercept:.0f}\n$R^2$={r**2:.2f}, p={p:.4f}')
    
    ax.set_xlabel('Total Atoms (Pt + Sn)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Melting Point $T_m$ (K)', fontsize=11, fontweight='bold')
    ax.set_title('(b) $T_m$ vs Total Atoms (Supported Alloy)', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    
    # --- (c) Tm (Lindemann) vs Tm (Clustering) comparison ---
    ax = axes[1, 0]
    df_both = df_valid[~df_valid['Tm_clustering'].isna()]
    
    for type_key in ['air', 'supported', 'oxide']:
        df_type = df_both[df_both['type'] == type_key]
        if not df_type.empty:
            ax.scatter(df_type['Tm_lindemann'], df_type['Tm_clustering'],
                      c=colors[type_key], marker=markers[type_key],
                      label=labels[type_key], s=80, alpha=0.8)
    
    # 对角线
    if not df_both.empty:
        lim_min = min(df_both['Tm_lindemann'].min(), df_both['Tm_clustering'].min()) - 50
        lim_max = max(df_both['Tm_lindemann'].max(), df_both['Tm_clustering'].max()) + 50
        ax.plot([lim_min, lim_max], [lim_min, lim_max], 'k--', alpha=0.5, label='y=x')
        ax.set_xlim(lim_min, lim_max)
        ax.set_ylim(lim_min, lim_max)
    
    ax.set_xlabel('$T_m$ (Lindemann) [K]', fontsize=11, fontweight='bold')
    ax.set_ylabel('$T_m$ (Clustering) [K]', fontsize=11, fontweight='bold')
    ax.set_title('(c) Comparison: Lindemann vs Clustering', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')
    
    # --- (d) Tm difference histogram ---
    ax = axes[1, 1]
    df_diff = df_both[~df_both['Tm_diff_lind_clust'].isna()]
    
    if not df_diff.empty:
        for type_key in ['supported', 'oxide']:
            df_type = df_diff[df_diff['type'] == type_key]
            if not df_type.empty:
                ax.hist(df_type['Tm_diff_lind_clust'], bins=15, alpha=0.6,
                       color=colors[type_key], label=labels[type_key],
                       edgecolor='black')
        
        # 统计信息
        mean_diff = df_diff['Tm_diff_lind_clust'].mean()
        std_diff = df_diff['Tm_diff_lind_clust'].std()
        ax.axvline(mean_diff, color='red', linestyle='--', linewidth=2,
                   label=f'Mean = {mean_diff:.1f} K')
        ax.axvline(0, color='black', linestyle='-', linewidth=1, alpha=0.5)
    
    ax.set_xlabel('$\\Delta T_m$ = $T_m$(Lindemann) - $T_m$(Clustering) [K]', fontsize=11, fontweight='bold')
    ax.set_ylabel('Count', fontsize=11, fontweight='bold')
    ax.set_title('(d) $T_m$ Difference Distribution', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, axis='y')
    
    # --- (e) Tm vs Sn fraction - 所有类型 (支撑材料，应用筛选) ---
    ax = axes[2, 0]
    for type_key in ['air', 'supported', 'oxide']:
        df_type = df_valid[
            (df_valid['type'] == type_key) &
            (~df_valid['structure'].str.lower().isin(exclude_list))  # 应用筛选
        ]
        if not df_type.empty:
            ax.scatter(df_type['Sn_fraction'], df_type['Tm_lindemann'],
                      c=colors[type_key], marker=markers[type_key],
                      label=labels[type_key], s=80, alpha=0.8)
    
    ax.set_xlabel('Sn Fraction = Sn/(Pt+Sn+O)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Melting Point $T_m$ (K)', fontsize=11, fontweight='bold')
    ax.set_title('(e) $T_m$ vs Sn Fraction (All Types)', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    
    # --- (f) Tm vs total atoms - 所有类型 (支撑材料，应用筛选) ---
    ax = axes[2, 1]
    for type_key in ['air', 'supported', 'oxide']:
        df_type = df_valid[
            (df_valid['type'] == type_key) &
            (~df_valid['structure'].str.lower().isin(exclude_list))  # 应用筛选
        ]
        if not df_type.empty:
            ax.scatter(df_type['total_atoms'], df_type['Tm_lindemann'],
                      c=colors[type_key], marker=markers[type_key],
                      label=labels[type_key], s=80, alpha=0.8)
    
    ax.set_xlabel('Total Atoms (Pt + Sn + O)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Melting Point $T_m$ (K)', fontsize=11, fontweight='bold')
    ax.set_title('(f) $T_m$ vs Total Atoms (All Types)', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    output_file = OUTPUT_DIR / 'Tm_vs_composition.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  [OK] Saved: {output_file}")
    plt.close()
    
    # ========================================================================
    # 图2: 含氧 vs 无氧对比
    # ========================================================================
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # --- 按是否含氧分组 ---
    df_no_o = df_valid[df_valid['has_oxygen'] == False]
    df_with_o = df_valid[df_valid['has_oxygen'] == True]
    
    # (a) 箱线图
    ax = axes[0]
    data_box = []
    labels_box = []
    
    if not df_no_o.empty:
        data_box.append(df_no_o['Tm_lindemann'].dropna().values)
        labels_box.append(f'No Oxygen\n(n={len(df_no_o)})')
    if not df_with_o.empty:
        data_box.append(df_with_o['Tm_lindemann'].dropna().values)
        labels_box.append(f'With Oxygen\n(n={len(df_with_o)})')
    
    if data_box:
        bp = ax.boxplot(data_box, labels=labels_box, patch_artist=True)
        colors_box = ['#2E8B57', '#FF6347']
        for patch, color in zip(bp['boxes'], colors_box):
            patch.set_facecolor(color)
            patch.set_alpha(0.6)
    
    ax.set_ylabel('Melting Point $T_m$ (K)', fontsize=11, fontweight='bold')
    ax.set_title('(a) $T_m$ by Oxygen Content', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')
    
    # 统计检验
    if len(df_no_o) > 1 and len(df_with_o) > 1:
        t_stat, p_value = stats.ttest_ind(
            df_no_o['Tm_lindemann'].dropna(),
            df_with_o['Tm_lindemann'].dropna()
        )
        ax.text(0.5, 0.95, f't-test p = {p_value:.4f}',
               transform=ax.transAxes, fontsize=10,
               verticalalignment='top', horizontalalignment='center',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # (b) Sn fraction 对比
    ax = axes[1]
    if not df_no_o.empty:
        ax.scatter(df_no_o['Sn_fraction'], df_no_o['Tm_lindemann'],
                  c='#2E8B57', marker='s', s=80, label='No Oxygen', alpha=0.7)
    if not df_with_o.empty:
        ax.scatter(df_with_o['Sn_fraction'], df_with_o['Tm_lindemann'],
                  c='#FF6347', marker='^', s=80, label='With Oxygen', alpha=0.7)
    
    ax.set_xlabel('Sn Fraction = Sn/(Pt+Sn+O)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Melting Point $T_m$ (K)', fontsize=11, fontweight='bold')
    ax.set_title('(b) $T_m$ vs Sn Fraction (by Oxygen)', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    output_file = OUTPUT_DIR / 'Tm_oxygen_comparison.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  [OK] Saved: {output_file}")
    plt.close()
    
    # ========================================================================
    # 图3: Pt8SnX 和 Pt6SnX 系列分析 (仅负载型无氧合金)
    # ========================================================================
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    
    # 获取屏蔽列表
    exclude_list = [s.lower() for s in args.exclude]
    
    # Pt8SnX 系列 - 仅负载型无氧合金 (排除 Air, 排除氧化物, 排除屏蔽的结构)
    ax = axes[0]
    df_pt8 = df_valid[
        (df_valid['Pt'] == 8) & 
        (df_valid['O'] == 0) & 
        (df_valid['type'] == 'supported') &  # 只要负载型
        (~df_valid['structure'].str.lower().isin(exclude_list))  # 排除屏蔽的结构
    ].copy()
    
    if not df_pt8.empty:
        df_pt8_sorted = df_pt8.sort_values('Sn')
        ax.plot(df_pt8_sorted['Sn'], df_pt8_sorted['Tm_lindemann'],
               'o-', color='black', markersize=10,
               linewidth=2, label=r'Pt$_8$Sn$_x$')
        
        # 线性拟合
        if len(df_pt8_sorted) > 2:
            slope, intercept, r, p, se = stats.linregress(
                df_pt8_sorted['Sn'], df_pt8_sorted['Tm_lindemann']
            )
            x_fit = np.linspace(df_pt8_sorted['Sn'].min(), df_pt8_sorted['Sn'].max(), 100)
            y_fit = slope * x_fit + intercept
            ax.plot(x_fit, y_fit, '--', color='gray', alpha=0.7,
                   label=f'Fit: slope={slope:.1f} K/Sn, R²={r**2:.3f}')
    
    ax.set_xlabel('Sn Atoms (x)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Melting Point $T_m$ (K)', fontsize=11, fontweight='bold')
    ax.set_title(r'(a) Pt$_8$Sn$_x$ Series', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9, loc='best')
    ax.grid(False)
    
    # Pt6SnX 系列 - 仅负载型无氧合金
    ax = axes[1]
    df_pt6 = df_valid[
        (df_valid['Pt'] == 6) & 
        (df_valid['O'] == 0) & 
        (df_valid['type'] == 'supported') &  # 只要负载型
        (~df_valid['structure'].str.lower().isin(exclude_list))  # 排除屏蔽的结构
    ].copy()
    
    if not df_pt6.empty:
        df_pt6_sorted = df_pt6.sort_values('Sn')
        ax.plot(df_pt6_sorted['Sn'], df_pt6_sorted['Tm_lindemann'],
               'o-', color='black', markersize=10,
               linewidth=2, label=r'Pt$_6$Sn$_x$')
        
        # 线性拟合
        if len(df_pt6_sorted) > 2:
            slope, intercept, r, p, se = stats.linregress(
                df_pt6_sorted['Sn'], df_pt6_sorted['Tm_lindemann']
            )
            x_fit = np.linspace(df_pt6_sorted['Sn'].min(), df_pt6_sorted['Sn'].max(), 100)
            y_fit = slope * x_fit + intercept
            ax.plot(x_fit, y_fit, '--', color='gray', alpha=0.7,
                   label=f'Fit: slope={slope:.1f} K/Sn, R²={r**2:.3f}')
    
    ax.set_xlabel('Sn Atoms (x)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Melting Point $T_m$ (K)', fontsize=11, fontweight='bold')
    ax.set_title(r'(b) Pt$_6$Sn$_x$ Series', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9, loc='best')
    ax.grid(False)
    
    # 总原子数为8的系列 (Pt+Sn=8, 无氧, 负载型)
    ax = axes[2]
    df_total8 = df_valid[
        ((df_valid['Pt'] + df_valid['Sn']) == 8) & 
        (df_valid['O'] == 0) & 
        (df_valid['type'] == 'supported') &  # 只要负载型
        (~df_valid['structure'].str.lower().isin(exclude_list))  # 排除屏蔽的结构
    ].copy()
    
    if not df_total8.empty:
        df_total8_sorted = df_total8.sort_values('Sn')
        ax.plot(df_total8_sorted['Sn'], df_total8_sorted['Tm_lindemann'],
               'o-', color='black', markersize=10,
               linewidth=2, label=r'Pt$_{8-x}$Sn$_x$')
        
        # 线性拟合
        if len(df_total8_sorted) > 2:
            slope, intercept, r, p, se = stats.linregress(
                df_total8_sorted['Sn'], df_total8_sorted['Tm_lindemann']
            )
            x_fit = np.linspace(df_total8_sorted['Sn'].min(), df_total8_sorted['Sn'].max(), 100)
            y_fit = slope * x_fit + intercept
            ax.plot(x_fit, y_fit, '--', color='gray', alpha=0.7,
                   label=f'Fit: slope={slope:.1f} K/Sn, R²={r**2:.3f}')
    
    ax.set_xlabel('Sn Atoms (x)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Melting Point $T_m$ (K)', fontsize=11, fontweight='bold')
    ax.set_title(r'(c) Pt$_{8-x}$Sn$_x$ Series (Total=8)', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9, loc='best')
    ax.grid(False)
    
    plt.tight_layout()
    
    output_file = OUTPUT_DIR / 'Tm_series_analysis.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  [OK] Saved: {output_file}")
    plt.close()


# ============================================================================
# 6. 生成报告
# ============================================================================

def generate_report(df_results):
    """
    生成分析报告
    """
    print(f"\n[4] Generating report...")
    
    df_valid = df_results[~df_results['Tm_lindemann'].isna()]
    df_both = df_valid[~df_valid['Tm_clustering'].isna()]
    
    # 统计
    n_total = len(df_results)
    n_valid = len(df_valid)
    
    # 按类型统计
    type_stats = df_valid.groupby('type').agg({
        'Tm_lindemann': ['mean', 'std', 'min', 'max', 'count']
    }).round(1)
    
    # 两种方法对比
    if len(df_both) > 0:
        mean_diff = df_both['Tm_diff_lind_clust'].mean()
        std_diff = df_both['Tm_diff_lind_clust'].std()
        corr = df_both['Tm_lindemann'].corr(df_both['Tm_clustering'])
    else:
        mean_diff = std_diff = corr = np.nan
    
    # 含氧 vs 无氧
    df_no_o = df_valid[df_valid['has_oxygen'] == False]
    df_with_o = df_valid[df_valid['has_oxygen'] == True]
    
    report = f"""# Melting Point Analysis Report

**Generated**: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}

---

## 1. Overview

- **Total structures analyzed**: {n_total}
- **Valid melting points (Lindemann)**: {n_valid}
- **Both methods available**: {len(df_both)}

---

## 2. Melting Point Methods Comparison

### Method A: Lindemann Index Jump
- Threshold: δ = 0.1
- Identifies temperature where δ crosses from <0.1 to ≥0.1
- Uses linear interpolation for precision

### Method B: Clustering Partition Boundary
- Uses K-means n=2 clustering results
- Boundary = midpoint between partition1 max T and partition2 min T

### Comparison Statistics
| Metric | Value |
|--------|-------|
| Mean difference (Lindemann - Clustering) | {mean_diff:.1f} K |
| Std. deviation of difference | {std_diff:.1f} K |
| Correlation coefficient | {corr:.4f} |

---

## 3. Melting Point by Structure Type

| Type | Mean Tm (K) | Std (K) | Min (K) | Max (K) | Count |
|------|-------------|---------|---------|---------|-------|
"""
    
    for type_key in ['air', 'supported', 'oxide']:
        df_type = df_valid[df_valid['type'] == type_key]
        if not df_type.empty:
            report += f"| {type_key.capitalize()} | {df_type['Tm_lindemann'].mean():.1f} | {df_type['Tm_lindemann'].std():.1f} | {df_type['Tm_lindemann'].min():.1f} | {df_type['Tm_lindemann'].max():.1f} | {len(df_type)} |\n"
    
    report += f"""
---

## 4. Oxygen Effect Analysis

| Group | Mean Tm (K) | Std (K) | Count |
|-------|-------------|---------|-------|
| Without Oxygen | {df_no_o['Tm_lindemann'].mean():.1f} | {df_no_o['Tm_lindemann'].std():.1f} | {len(df_no_o)} |
| With Oxygen | {df_with_o['Tm_lindemann'].mean():.1f} | {df_with_o['Tm_lindemann'].std():.1f} | {len(df_with_o)} |

"""
    
    # t-test
    if len(df_no_o) > 1 and len(df_with_o) > 1:
        t_stat, p_value = stats.ttest_ind(
            df_no_o['Tm_lindemann'].dropna(),
            df_with_o['Tm_lindemann'].dropna()
        )
        report += f"**t-test**: t = {t_stat:.3f}, p = {p_value:.4f}\n\n"
        if p_value < 0.05:
            report += "> ⚠️ **Significant difference** between oxygen-containing and oxygen-free structures (p < 0.05)\n"
        else:
            report += "> ✅ No significant difference (p ≥ 0.05)\n"
    
    report += f"""
---

## 5. Series Analysis

### Pt₈Snₓ Series
"""
    df_pt8 = df_valid[(df_valid['Pt'] == 8) & (df_valid['O'] == 0)]
    if not df_pt8.empty and len(df_pt8) > 2:
        slope, intercept, r, p, se = stats.linregress(df_pt8['Sn'], df_pt8['Tm_lindemann'])
        report += f"- Linear fit: Tm = {slope:.1f} × Sn + {intercept:.1f}\n"
        report += f"- R² = {r**2:.4f}\n"
        report += f"- **Effect of Sn**: {slope:.1f} K per Sn atom\n"
    
    report += f"""
### Pt₆Snₓ Series
"""
    df_pt6 = df_valid[(df_valid['Pt'] == 6) & (df_valid['O'] == 0)]
    if not df_pt6.empty and len(df_pt6) > 2:
        slope, intercept, r, p, se = stats.linregress(df_pt6['Sn'], df_pt6['Tm_lindemann'])
        report += f"- Linear fit: Tm = {slope:.1f} × Sn + {intercept:.1f}\n"
        report += f"- R² = {r**2:.4f}\n"
        report += f"- **Effect of Sn**: {slope:.1f} K per Sn atom\n"
    
    report += f"""
---

## 6. Full Results Table

| Structure | Type | Pt | Sn | O | Tm (Lindemann) | Tm (Transition) | Tm (Clustering) | Δ(L-T) |
|-----------|------|----|----|---|----------------|-----------------|-----------------|--------|
"""
    
    for _, row in df_results.iterrows():
        Tm_L = f"{row['Tm_lindemann']:.0f}" if not np.isnan(row['Tm_lindemann']) else "N/A"
        Tm_T = f"{row['Tm_transition']:.0f}" if not np.isnan(row.get('Tm_transition', np.nan)) else "N/A"
        Tm_C = f"{row['Tm_clustering']:.0f}" if not np.isnan(row['Tm_clustering']) else "N/A"
        diff = f"{row['Tm_diff_lind_trans']:.0f}" if not np.isnan(row.get('Tm_diff_lind_trans', np.nan)) else "N/A"
        report += f"| {row['structure']} | {row['type']} | {row['Pt']} | {row['Sn']} | {row['O']} | {Tm_L} | {Tm_T} | {Tm_C} | {diff} |\n"
    
    report += """
---

*Report generated by step5_1_melting_point_analysis.py*
"""
    
    output_file = OUTPUT_DIR / 'MELTING_POINT_ANALYSIS_REPORT.md'
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(report)
    print(f"  [OK] Saved: {output_file}")


# ============================================================================
# 主函数
# ============================================================================

def main():
    # 1. 分析熔点
    df_results = analyze_all_melting_points()
    
    if df_results is None or df_results.empty:
        print("[ERROR] No results to analyze")
        return
    
    # 2. 可视化
    plot_melting_point_analysis(df_results)
    
    # 3. 生成报告
    generate_report(df_results)
    
    # 4. 打印摘要
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    
    df_valid = df_results[~df_results['Tm_lindemann'].isna()]
    
    print(f"\nMelting Point Statistics (Lindemann method):")
    print(f"  Total structures with Tm: {len(df_valid)}")
    print(f"  Mean Tm: {df_valid['Tm_lindemann'].mean():.1f} ± {df_valid['Tm_lindemann'].std():.1f} K")
    print(f"  Range: {df_valid['Tm_lindemann'].min():.0f} - {df_valid['Tm_lindemann'].max():.0f} K")
    
    # 按类型
    print(f"\nBy structure type:")
    for type_key in ['air', 'supported', 'oxide']:
        df_type = df_valid[df_valid['type'] == type_key]
        if not df_type.empty:
            print(f"  {type_key.capitalize():12s}: {df_type['Tm_lindemann'].mean():.1f} ± {df_type['Tm_lindemann'].std():.1f} K (n={len(df_type)})")
    
    print("\n" + "=" * 70)
    print("Done!")
    print("=" * 70)


if __name__ == '__main__':
    main()
