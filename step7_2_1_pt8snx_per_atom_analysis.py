#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
Step 7.2.1: Pt8SnX 系列 Per-Atom 动力学分析
================================================================================

专注于 Pt8SnX (x=0,1,2,...,10) 系列的每原子分析:
- 每个原子的林德曼指数 δ (分 Pt/Sn 元素)
- 每个原子的扩散系数 D (分 Pt/Sn 元素)
- 不同 Sn 含量对动力学的影响

数据来源:
- Per-atom Lindemann: data/lindemann/per-atoms/per_atom_master_run_*.csv
- Per-atom D: data/gmx_msd/per-atom/collected_gmx_per_atom_msd/*.csv

输出:
- 热力图: δ(Pt), δ(Sn), D(Pt), D(Sn) vs (Sn含量, 温度)
- 对比图: Pt vs Sn 的动力学差异
- 统计表: 分元素、分温度的统计
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import re
import argparse
from scipy import stats

# ============================================================================
# 配置
# ============================================================================

# 使用 Times New Roman 字体
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['font.size'] = 12

# 路径
BASE_DIR = Path(__file__).parent
LINDEMANN_DIR = BASE_DIR / 'data' / 'lindemann' / 'per-atoms'
MSD_DIR = BASE_DIR / 'data' / 'gmx_msd' / 'per-atom' / 'collected_gmx_per_atom_msd'
OUTPUT_DIR = BASE_DIR / 'results' / 'pt8snx_per_atom_dynamics'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# 阈值
HIGH_LINDEMANN = 0.10  # 高振动阈值
HIGH_D = 0.01  # 高扩散阈值 (10^-5 cm²/s)

# 颜色
ELEMENT_COLORS = {'Pt': '#1f77b4', 'Sn': '#ff7f0e', 'O': '#2ca02c'}


# ============================================================================
# 数据加载
# ============================================================================

def extract_sn_count(structure_name):
    """从结构名提取 Sn 原子数量"""
    # pt8sn0-2-best -> 0
    # pt8sn10-2-best -> 10
    match = re.search(r'pt8sn(\d+)', str(structure_name).lower())
    if match:
        return int(match.group(1))
    return None


def extract_full_run_key(path):
    """从完整路径提取唯一的 run key (使用最后4层目录)
    
    完整路径: /home/.../GPU-Pt8/Pt8-2/pt8sn0-2-best/T1000.r7.gpu0
    提取:     GPU-Pt8/Pt8-2/pt8sn0-2-best/T1000.r7.gpu0
    
    这4层包含:
    - GPU-Pt8: GPU类型
    - Pt8-2: 批次 (Pt8, Pt8-2, Pt8-3)
    - pt8sn0-2-best: 结构名
    - T1000.r7.gpu0: 温度.run号.gpu
    """
    if pd.isna(path):
        return None
    parts = str(path).split('/')
    # 取最后4层作为唯一标识
    if len(parts) >= 4:
        return '/'.join(parts[-4:])
    return str(path)


def load_lindemann_data():
    """加载 Per-atom Lindemann 数据"""
    print("\n[1] 加载 Per-atom Lindemann 数据...")
    
    files = sorted(LINDEMANN_DIR.glob('per_atom_master_run_*.csv'))
    if not files:
        print("  [ERROR] 未找到 Lindemann 文件")
        return None
    
    df_list = []
    for f in files:
        df = pd.read_csv(f)
        df_list.append(df)
        print(f"  - {f.name}: {len(df)} records")
    
    df = pd.concat(df_list, ignore_index=True)
    
    # 使用完整目录作为唯一标识
    df['full_run_key'] = df['目录'].apply(extract_full_run_key)
    
    # 从目录提取结构名 (倒数第2层)
    def extract_structure(path):
        if pd.isna(path):
            return None
        parts = str(path).split('/')
        if len(parts) >= 2:
            # 倒数第2层是结构名
            return parts[-2].lower()
        return None
    
    df['structure'] = df['目录'].apply(extract_structure)
    
    # 提取温度 (从倒数第1层)
    def extract_temp(path):
        if pd.isna(path):
            return None
        parts = str(path).split('/')
        if len(parts) >= 1:
            match = re.search(r'T(\d+)', parts[-1])
            return int(match.group(1)) if match else None
        return None
    
    df['temp'] = df['目录'].apply(extract_temp)
    
    # 筛选 Pt8SnX 系列
    df = df[df['structure'].str.contains('pt8sn', case=False, na=False)]
    df['sn_count'] = df['structure'].apply(extract_sn_count)
    
    # 重命名列
    df = df.rename(columns={
        'lindemann_index': 'delta',
        'element': 'element'
    })
    
    # 验证唯一性
    dup_check = df.groupby(['full_run_key', 'atom_id']).size()
    print(f"  [OK] Pt8SnX Lindemann 数据: {len(df)} records")
    print(f"       Sn 含量: {sorted(df['sn_count'].unique())}")
    print(f"       元素: {df['element'].unique()}")
    print(f"       唯一 full_run_key 数: {df['full_run_key'].nunique()}")
    print(f"       每个 (full_run_key, atom_id) 记录数: {dup_check.value_counts().to_dict()}")
    
    return df[['structure', 'temp', 'atom_id', 'element', 'delta', 'sn_count', 'full_run_key']]


def load_msd_data():
    """加载 Per-atom MSD/D 数据"""
    print("\n[2] 加载 Per-atom MSD 数据...")
    
    files = sorted(MSD_DIR.glob('per_atom_diffusion_coefficients_*.csv'))
    if not files:
        print("  [ERROR] 未找到 MSD 文件")
        return None
    
    df_list = []
    for f in files:
        df = pd.read_csv(f)
        df_list.append(df)
        print(f"  - {f.name}: {len(df)} records")
    
    df = pd.concat(df_list, ignore_index=True)
    
    # 使用完整目录作为唯一标识
    df['full_run_key'] = df['完整目录路径'].apply(extract_full_run_key)
    
    # 筛选 Pt8SnX 系列
    df = df[df['结构'].str.contains('pt8sn', case=False, na=False)]
    df['sn_count'] = df['结构'].apply(extract_sn_count)
    
    # 重命名列
    df = df.rename(columns={
        '结构': 'structure',
        '温度(K)': 'temp',
        '元素': 'element',
        'D(1e-5 cm²/s)': 'D'
    })
    
    # 验证唯一性
    dup_check = df.groupby(['full_run_key', 'atom_id']).size()
    print(f"  [OK] Pt8SnX MSD 数据: {len(df)} records")
    print(f"       Sn 含量: {sorted(df['sn_count'].unique())}")
    print(f"       唯一 full_run_key 数: {df['full_run_key'].nunique()}")
    print(f"       每个 (full_run_key, atom_id) 记录数: {dup_check.value_counts().to_dict()}")
    
    return df[['structure', 'temp', 'atom_id', 'element', 'D', 'sn_count', 'full_run_key']]


def merge_data(df_lindemann, df_msd):
    """合并 Lindemann 和 MSD 数据 - 使用完整目录路径精确匹配"""
    print("\n[3] 合并数据 (按 full_run_key + atom_id 精确匹配)...")
    
    # 检查 full_run_key 重叠情况
    lin_keys = set(df_lindemann['full_run_key'].dropna().unique())
    msd_keys = set(df_msd['full_run_key'].dropna().unique())
    common_keys = lin_keys & msd_keys
    
    print(f"  Lindemann full_run_keys: {len(lin_keys)}")
    print(f"  MSD full_run_keys: {len(msd_keys)}")
    print(f"  共同 full_run_keys: {len(common_keys)}")
    
    if len(common_keys) == 0:
        print("  [WARNING] 没有共同的 full_run_key!")
        print("  Lindemann 示例:", list(lin_keys)[:3])
        print("  MSD 示例:", list(msd_keys)[:3])
        return None
    
    # 按 full_run_key + atom_id 精确匹配
    df = pd.merge(
        df_lindemann,
        df_msd[['full_run_key', 'atom_id', 'D']],
        on=['full_run_key', 'atom_id'],
        how='inner'
    )
    
    # 验证：每个 (full_run_key, atom_id) 应该只有一条记录
    dup_check = df.groupby(['full_run_key', 'atom_id']).size()
    if dup_check.max() > 1:
        print(f"  [ERROR] 存在重复匹配! 最大重复数: {dup_check.max()}")
        print(f"  重复分布: {dup_check.value_counts().to_dict()}")
    else:
        print(f"  [OK] 每个 (full_run_key, atom_id) 只有一条记录 ✓")
    
    print(f"  [OK] 合并后: {len(df)} records")
    print(f"       结构数: {df['structure'].nunique()}")
    print(f"       温度: {sorted(df['temp'].unique())}")
    print(f"       运行数: {df['full_run_key'].nunique()}")
    
    return df


# ============================================================================
# 统计计算
# ============================================================================

def calculate_statistics(df):
    """计算分元素、分Sn含量、分温度的统计"""
    print("\n[4] 计算统计...")
    
    stats_list = []
    
    for (sn_count, temp, element), group in df.groupby(['sn_count', 'temp', 'element']):
        stats = {
            'sn_count': sn_count,
            'temp': temp,
            'element': element,
            'n_atoms': len(group),
            'delta_mean': group['delta'].mean(),
            'delta_std': group['delta'].std(),
            'delta_median': group['delta'].median(),
            'D_mean': group['D'].mean(),
            'D_std': group['D'].std(),
            'D_median': group['D'].median(),
            'high_delta_ratio': (group['delta'] > HIGH_LINDEMANN).mean(),
            'high_D_ratio': (group['D'] > HIGH_D).mean(),
        }
        stats_list.append(stats)
    
    df_stats = pd.DataFrame(stats_list)
    print(f"  [OK] 统计记录数: {len(df_stats)}")
    
    return df_stats


# ============================================================================
# 可视化 - 核心热力图
# ============================================================================

def plot_heatmaps(df_stats):
    """绘制核心热力图: δ(Pt), δ(Sn), D(Pt), D(Sn) vs (Sn含量, 温度)"""
    print("\n[5] 绘制热力图...")
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # 准备数据
    elements = ['Pt', 'Sn']
    metrics = [('delta_mean', 'Lindemann Index (δ)'), ('D_mean', 'D (10⁻⁵ cm²/s)')]
    
    for row_idx, (metric, metric_label) in enumerate(metrics):
        for col_idx, element in enumerate(elements):
            ax = axes[row_idx, col_idx]
            
            # 筛选元素数据
            elem_data = df_stats[df_stats['element'] == element]
            
            # 创建透视表
            pivot = elem_data.pivot_table(
                index='temp', 
                columns='sn_count', 
                values=metric, 
                aggfunc='mean'
            )
            
            # 对 D 值使用对数刻度
            if metric == 'D_mean':
                pivot_plot = np.log10(pivot + 1e-10)
                fmt = '.2f'
                cbar_label = f'log₁₀({metric_label})'
                cmap = 'YlOrRd'
            else:
                pivot_plot = pivot
                fmt = '.3f'
                cbar_label = metric_label
                cmap = 'RdYlGn_r'  # 红色=高振动
            
            # 绘制热力图
            sns.heatmap(pivot_plot, annot=True, fmt=fmt, cmap=cmap, ax=ax,
                       cbar_kws={'label': cbar_label}, linewidths=0.5)
            
            ax.set_xlabel('Sn Count in Pt8SnX', fontsize=12)
            ax.set_ylabel('Temperature (K)', fontsize=12)
            ax.set_title(f'{element}: {metric_label}', fontsize=13, fontweight='bold')
    
    plt.suptitle('Pt8SnX Series: Per-Atom Dynamics Heatmap\n'
                 '(Left: Pt atoms, Right: Sn atoms)', 
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    
    output_path = OUTPUT_DIR / 'pt8snx_heatmap_delta_D.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"  [SAVED] {output_path}")
    
    return fig


def plot_high_ratio_heatmaps(df_stats):
    """绘制高振动/高扩散比例热力图"""
    print("\n[6] 绘制高活跃比例热力图...")
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    elements = ['Pt', 'Sn']
    metrics = [('high_delta_ratio', 'High δ Ratio (δ>0.10)'), 
               ('high_D_ratio', 'High D Ratio (D>0.01)')]
    
    for row_idx, (metric, metric_label) in enumerate(metrics):
        for col_idx, element in enumerate(elements):
            ax = axes[row_idx, col_idx]
            
            elem_data = df_stats[df_stats['element'] == element]
            pivot = elem_data.pivot_table(
                index='temp', 
                columns='sn_count', 
                values=metric, 
                aggfunc='mean'
            ) * 100  # 转换为百分比
            
            sns.heatmap(pivot, annot=True, fmt='.1f', cmap='YlOrRd', ax=ax,
                       cbar_kws={'label': f'{metric_label} (%)'}, 
                       vmin=0, vmax=100, linewidths=0.5)
            
            ax.set_xlabel('Sn Count in Pt8SnX', fontsize=12)
            ax.set_ylabel('Temperature (K)', fontsize=12)
            ax.set_title(f'{element}: {metric_label}', fontsize=13, fontweight='bold')
    
    plt.suptitle('Pt8SnX Series: High Activity Ratio Heatmap\n'
                 '(Percentage of atoms exceeding threshold)', 
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    
    output_path = OUTPUT_DIR / 'pt8snx_heatmap_high_ratio.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"  [SAVED] {output_path}")
    
    return fig


def plot_element_comparison(df_stats):
    """绘制 Pt vs Sn 元素对比图"""
    print("\n[7] 绘制元素对比图...")
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    
    # 特征温度
    feature_temps = [300, 600, 900]
    
    # 第一行: δ vs Sn含量
    for col_idx, temp in enumerate(feature_temps):
        ax = axes[0, col_idx]
        
        temp_data = df_stats[df_stats['temp'] == temp]
        if len(temp_data) == 0:
            # 找最接近的温度
            available = df_stats['temp'].unique()
            temp = min(available, key=lambda x: abs(x - temp))
            temp_data = df_stats[df_stats['temp'] == temp]
        
        for element in ['Pt', 'Sn']:
            elem_data = temp_data[temp_data['element'] == element].sort_values('sn_count')
            if len(elem_data) > 0:
                ax.errorbar(elem_data['sn_count'], elem_data['delta_mean'],
                           yerr=elem_data['delta_std'], 
                           fmt='o-', color=ELEMENT_COLORS[element],
                           linewidth=2, markersize=8, capsize=4,
                           label=element)
        
        ax.axhline(HIGH_LINDEMANN, color='red', linestyle='--', alpha=0.5)
        ax.set_xlabel('Sn Count (x in Pt8SnX)', fontsize=12)
        ax.set_ylabel('Lindemann Index (δ)', fontsize=12)
        ax.set_title(f'T = {int(temp)} K', fontsize=13, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_xticks(range(0, 11))
    
    # 第二行: D vs Sn含量
    for col_idx, temp in enumerate(feature_temps):
        ax = axes[1, col_idx]
        
        temp_data = df_stats[df_stats['temp'] == temp]
        if len(temp_data) == 0:
            available = df_stats['temp'].unique()
            temp = min(available, key=lambda x: abs(x - temp))
            temp_data = df_stats[df_stats['temp'] == temp]
        
        for element in ['Pt', 'Sn']:
            elem_data = temp_data[temp_data['element'] == element].sort_values('sn_count')
            if len(elem_data) > 0:
                ax.errorbar(elem_data['sn_count'], elem_data['D_mean'],
                           yerr=elem_data['D_std'], 
                           fmt='s-', color=ELEMENT_COLORS[element],
                           linewidth=2, markersize=8, capsize=4,
                           label=element)
        
        ax.axhline(HIGH_D, color='blue', linestyle='--', alpha=0.5)
        ax.set_xlabel('Sn Count (x in Pt8SnX)', fontsize=12)
        ax.set_ylabel('D (10⁻⁵ cm²/s)', fontsize=12)
        ax.set_title(f'T = {int(temp)} K', fontsize=13, fontweight='bold')
        ax.legend()
        ax.set_yscale('log')
        ax.grid(True, alpha=0.3)
        ax.set_xticks(range(0, 11))
    
    plt.suptitle('Pt8SnX Series: Element Comparison\n'
                 '(How does Sn content affect Pt and Sn atom dynamics?)', 
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    
    output_path = OUTPUT_DIR / 'pt8snx_element_comparison.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"  [SAVED] {output_path}")
    
    return fig


def plot_temperature_curves(df_stats):
    """绘制 δ 和 D vs 温度曲线 (按 Sn 含量分组)"""
    print("\n[8] 绘制温度曲线...")
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # 选择几个代表性的 Sn 含量
    sn_counts = [0, 2, 4, 6, 8, 10]
    colors = plt.cm.viridis(np.linspace(0, 1, len(sn_counts)))
    
    # (1) Pt: δ vs T
    ax = axes[0, 0]
    for sn, color in zip(sn_counts, colors):
        data = df_stats[(df_stats['sn_count'] == sn) & (df_stats['element'] == 'Pt')]
        data = data.sort_values('temp')
        if len(data) > 0:
            ax.plot(data['temp'], data['delta_mean'], 'o-', color=color,
                   linewidth=2, markersize=6, label=f'Pt8Sn{sn}')
    ax.axhline(HIGH_LINDEMANN, color='red', linestyle='--', alpha=0.5)
    ax.set_xlabel('Temperature (K)', fontsize=12)
    ax.set_ylabel('Lindemann Index (δ)', fontsize=12)
    ax.set_title('Pt atoms: δ vs Temperature', fontsize=13, fontweight='bold')
    ax.legend(ncol=2, fontsize=9)
    ax.grid(True, alpha=0.3)
    
    # (2) Sn: δ vs T
    ax = axes[0, 1]
    for sn, color in zip(sn_counts, colors):
        if sn == 0:
            continue  # Pt8Sn0 没有 Sn 原子
        data = df_stats[(df_stats['sn_count'] == sn) & (df_stats['element'] == 'Sn')]
        data = data.sort_values('temp')
        if len(data) > 0:
            ax.plot(data['temp'], data['delta_mean'], 's-', color=color,
                   linewidth=2, markersize=6, label=f'Pt8Sn{sn}')
    ax.axhline(HIGH_LINDEMANN, color='red', linestyle='--', alpha=0.5)
    ax.set_xlabel('Temperature (K)', fontsize=12)
    ax.set_ylabel('Lindemann Index (δ)', fontsize=12)
    ax.set_title('Sn atoms: δ vs Temperature', fontsize=13, fontweight='bold')
    ax.legend(ncol=2, fontsize=9)
    ax.grid(True, alpha=0.3)
    
    # (3) Pt: D vs T
    ax = axes[1, 0]
    for sn, color in zip(sn_counts, colors):
        data = df_stats[(df_stats['sn_count'] == sn) & (df_stats['element'] == 'Pt')]
        data = data.sort_values('temp')
        if len(data) > 0:
            ax.plot(data['temp'], data['D_mean'], 'o-', color=color,
                   linewidth=2, markersize=6, label=f'Pt8Sn{sn}')
    ax.axhline(HIGH_D, color='blue', linestyle='--', alpha=0.5)
    ax.set_xlabel('Temperature (K)', fontsize=12)
    ax.set_ylabel('D (10⁻⁵ cm²/s)', fontsize=12)
    ax.set_title('Pt atoms: D vs Temperature', fontsize=13, fontweight='bold')
    ax.legend(ncol=2, fontsize=9)
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3)
    
    # (4) Sn: D vs T
    ax = axes[1, 1]
    for sn, color in zip(sn_counts, colors):
        if sn == 0:
            continue
        data = df_stats[(df_stats['sn_count'] == sn) & (df_stats['element'] == 'Sn')]
        data = data.sort_values('temp')
        if len(data) > 0:
            ax.plot(data['temp'], data['D_mean'], 's-', color=color,
                   linewidth=2, markersize=6, label=f'Pt8Sn{sn}')
    ax.axhline(HIGH_D, color='blue', linestyle='--', alpha=0.5)
    ax.set_xlabel('Temperature (K)', fontsize=12)
    ax.set_ylabel('D (10⁻⁵ cm²/s)', fontsize=12)
    ax.set_title('Sn atoms: D vs Temperature', fontsize=13, fontweight='bold')
    ax.legend(ncol=2, fontsize=9)
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3)
    
    plt.suptitle('Pt8SnX Series: Temperature Dependence\n'
                 '(δ and D vs Temperature for different Sn contents)', 
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    
    output_path = OUTPUT_DIR / 'pt8snx_temperature_curves.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"  [SAVED] {output_path}")
    
    return fig


def plot_delta_vs_D_scatter(df):
    """绘制 δ vs D 散点图 (按 Sn 含量着色)"""
    print("\n[9] 绘制 δ vs D 散点图...")
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # 选择高温数据 (900K)
    temp = 900
    available = df['temp'].unique()
    if temp not in available:
        temp = max(available)
    
    df_plot = df[df['temp'] == temp]
    
    # (1) Pt 原子
    ax = axes[0]
    pt_data = df_plot[df_plot['element'] == 'Pt']
    
    scatter = ax.scatter(pt_data['delta'], pt_data['D'], 
                        c=pt_data['sn_count'], cmap='viridis',
                        s=30, alpha=0.6)
    plt.colorbar(scatter, ax=ax, label='Sn Count')
    
    ax.axvline(HIGH_LINDEMANN, color='red', linestyle='--', alpha=0.5)
    ax.axhline(HIGH_D, color='blue', linestyle='--', alpha=0.5)
    ax.set_xlabel('Lindemann Index (δ)', fontsize=12)
    ax.set_ylabel('D (10⁻⁵ cm²/s)', fontsize=12)
    ax.set_title(f'Pt atoms at {int(temp)}K', fontsize=13, fontweight='bold')
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3)
    
    # (2) Sn 原子
    ax = axes[1]
    sn_data = df_plot[df_plot['element'] == 'Sn']
    
    scatter = ax.scatter(sn_data['delta'], sn_data['D'], 
                        c=sn_data['sn_count'], cmap='viridis',
                        s=30, alpha=0.6)
    plt.colorbar(scatter, ax=ax, label='Sn Count')
    
    ax.axvline(HIGH_LINDEMANN, color='red', linestyle='--', alpha=0.5)
    ax.axhline(HIGH_D, color='blue', linestyle='--', alpha=0.5)
    ax.set_xlabel('Lindemann Index (δ)', fontsize=12)
    ax.set_ylabel('D (10⁻⁵ cm²/s)', fontsize=12)
    ax.set_title(f'Sn atoms at {int(temp)}K', fontsize=13, fontweight='bold')
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3)
    
    plt.suptitle(f'Pt8SnX Series: δ vs D Correlation at {int(temp)}K\n'
                 '(Color = Sn content in Pt8SnX)', 
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    
    output_path = OUTPUT_DIR / f'pt8snx_delta_vs_D_{int(temp)}K.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"  [SAVED] {output_path}")
    
    return fig


def plot_delta_D_correlation_analysis(df):
    """绘制 δ vs D 相关性深度分析图"""
    print("\n[10] 绘制 δ-D 相关性分析图...")
    
    from scipy import stats
    
    # ===== 图1: 分温度相关性散点图 =====
    fig1, axes1 = plt.subplots(2, 3, figsize=(15, 10))
    
    feature_temps = [300, 600, 900]
    
    for col, temp in enumerate(feature_temps):
        td = df[df['temp'] == temp]
        
        # 上行: Pt
        ax = axes1[0, col]
        pt = td[td['element'] == 'Pt']
        if len(pt) > 0:
            ax.scatter(pt['delta'], pt['D'], c=ELEMENT_COLORS['Pt'], alpha=0.4, s=20)
            r, p = stats.pearsonr(pt['delta'], pt['D'])
            ax.set_xlabel('Lindemann Index (δ)', fontsize=11)
            ax.set_ylabel('D (10⁻⁵ cm²/s)', fontsize=11)
            ax.set_title(f'Pt @ {temp}K (r={r:.3f})', fontsize=12, fontweight='bold')
            ax.set_yscale('log')
            ax.axvline(HIGH_LINDEMANN, color='red', linestyle='--', alpha=0.5)
            ax.axhline(HIGH_D, color='blue', linestyle='--', alpha=0.5)
            ax.grid(True, alpha=0.3)
        
        # 下行: Sn
        ax = axes1[1, col]
        sn = td[td['element'] == 'Sn']
        if len(sn) > 0:
            ax.scatter(sn['delta'], sn['D'], c=ELEMENT_COLORS['Sn'], alpha=0.4, s=20)
            r, p = stats.pearsonr(sn['delta'], sn['D'])
            ax.set_xlabel('Lindemann Index (δ)', fontsize=11)
            ax.set_ylabel('D (10⁻⁵ cm²/s)', fontsize=11)
            ax.set_title(f'Sn @ {temp}K (r={r:.3f})', fontsize=12, fontweight='bold')
            ax.set_yscale('log')
            ax.axvline(HIGH_LINDEMANN, color='red', linestyle='--', alpha=0.5)
            ax.axhline(HIGH_D, color='blue', linestyle='--', alpha=0.5)
            ax.grid(True, alpha=0.3)
    
    plt.suptitle('Per-Atom δ vs D Correlation by Temperature\n'
                 '(Each point = one atom in one simulation)', 
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    
    output_path = OUTPUT_DIR / 'delta_D_correlation_by_temp.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"  [SAVED] {output_path}")
    
    # ===== 图2: 相关系数 vs 温度 =====
    fig2, axes2 = plt.subplots(1, 2, figsize=(12, 5))
    
    temps = sorted(df['temp'].unique())
    
    # 分元素计算相关系数
    corr_data = {'temp': [], 'element': [], 'r': [], 'n': []}
    
    for temp in temps:
        td = df[df['temp'] == temp]
        for elem in ['Pt', 'Sn']:
            ed = td[td['element'] == elem]
            if len(ed) > 10:
                r, p = stats.pearsonr(ed['delta'], ed['D'])
                corr_data['temp'].append(temp)
                corr_data['element'].append(elem)
                corr_data['r'].append(r)
                corr_data['n'].append(len(ed))
    
    df_corr = pd.DataFrame(corr_data)
    
    # (1) r vs T 曲线
    ax = axes2[0]
    for elem in ['Pt', 'Sn']:
        ed = df_corr[df_corr['element'] == elem]
        ax.plot(ed['temp'], ed['r'], 'o-', color=ELEMENT_COLORS[elem], 
               linewidth=2, markersize=8, label=elem)
    
    ax.axhline(0, color='black', linestyle='-', alpha=0.3)
    ax.axhline(0.3, color='gray', linestyle='--', alpha=0.5, label='Weak/Medium threshold')
    ax.set_xlabel('Temperature (K)', fontsize=12)
    ax.set_ylabel('Pearson Correlation (r)', fontsize=12)
    ax.set_title('δ-D Correlation vs Temperature', fontsize=13, fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-0.1, 0.6)
    
    # (2) 相关性强度热力图
    ax = axes2[1]
    
    # 按温度和Sn含量计算相关性
    corr_matrix = []
    sn_counts = sorted(df['sn_count'].unique())
    
    for temp in temps:
        row = []
        for sn in sn_counts:
            subset = df[(df['temp'] == temp) & (df['sn_count'] == sn)]
            if len(subset) > 5:
                r, p = stats.pearsonr(subset['delta'], subset['D'])
                row.append(r)
            else:
                row.append(np.nan)
        corr_matrix.append(row)
    
    corr_df = pd.DataFrame(corr_matrix, index=temps, columns=sn_counts)
    
    sns.heatmap(corr_df, annot=True, fmt='.2f', cmap='RdYlGn', ax=ax,
               center=0, vmin=-0.5, vmax=0.7, linewidths=0.5,
               cbar_kws={'label': 'Pearson r'})
    ax.set_xlabel('Sn Count', fontsize=12)
    ax.set_ylabel('Temperature (K)', fontsize=12)
    ax.set_title('δ-D Correlation Heatmap\n(by Sn content and Temperature)', 
                fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    
    output_path = OUTPUT_DIR / 'delta_D_correlation_summary.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"  [SAVED] {output_path}")
    
    # ===== 图3: 四象限分类图 =====
    fig3, axes3 = plt.subplots(1, 3, figsize=(15, 5))
    
    for col, temp in enumerate([300, 600, 900]):
        ax = axes3[col]
        td = df[df['temp'] == temp]
        
        # 分类
        td = td.copy()
        td['category'] = 'Stable'
        td.loc[(td['delta'] > HIGH_LINDEMANN) & (td['D'] <= HIGH_D), 'category'] = 'Vibrating'
        td.loc[(td['delta'] <= HIGH_LINDEMANN) & (td['D'] > HIGH_D), 'category'] = 'Diffusing'
        td.loc[(td['delta'] > HIGH_LINDEMANN) & (td['D'] > HIGH_D), 'category'] = 'Active'
        
        # 统计
        cat_counts = td['category'].value_counts()
        cat_pct = cat_counts / len(td) * 100
        
        colors_map = {'Stable': '#2ca02c', 'Vibrating': '#ff7f0e', 
                     'Diffusing': '#1f77b4', 'Active': '#d62728'}
        
        # 绘制散点
        for cat in ['Stable', 'Vibrating', 'Diffusing', 'Active']:
            cat_data = td[td['category'] == cat]
            if len(cat_data) > 0:
                pct = cat_pct.get(cat, 0)
                ax.scatter(cat_data['delta'], cat_data['D'], 
                          c=colors_map[cat], alpha=0.4, s=15,
                          label=f'{cat} ({pct:.1f}%)')
        
        ax.axvline(HIGH_LINDEMANN, color='red', linestyle='--', alpha=0.7)
        ax.axhline(HIGH_D, color='blue', linestyle='--', alpha=0.7)
        ax.set_xlabel('Lindemann Index (δ)', fontsize=11)
        ax.set_ylabel('D (10⁻⁵ cm²/s)', fontsize=11)
        ax.set_title(f'{temp}K', fontsize=13, fontweight='bold')
        ax.set_yscale('log')
        ax.legend(loc='lower right', fontsize=9)
        ax.grid(True, alpha=0.3)
        
        # 添加象限标签
        ax.text(0.02, 0.98, 'Stable\n(Low δ, Low D)', transform=ax.transAxes, 
               fontsize=9, va='top', ha='left', color='green')
        ax.text(0.98, 0.98, 'Vibrating\n(High δ, Low D)', transform=ax.transAxes, 
               fontsize=9, va='top', ha='right', color='orange')
        ax.text(0.02, 0.02, 'Diffusing\n(Low δ, High D)', transform=ax.transAxes, 
               fontsize=9, va='bottom', ha='left', color='blue')
        ax.text(0.98, 0.02, 'Active\n(High δ, High D)', transform=ax.transAxes, 
               fontsize=9, va='bottom', ha='right', color='red')
    
    plt.suptitle('Atom Mobility Classification: Four Quadrants\n'
                 f'(Thresholds: δ={HIGH_LINDEMANN}, D={HIGH_D})', 
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    
    output_path = OUTPUT_DIR / 'delta_D_quadrant_analysis.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"  [SAVED] {output_path}")
    
    # 保存相关性统计
    df_corr.to_csv(OUTPUT_DIR / 'delta_D_correlation_stats.csv', index=False, encoding='utf-8-sig')
    print(f"  [SAVED] {OUTPUT_DIR / 'delta_D_correlation_stats.csv'}")
    
    return fig1, fig2, fig3


# ============================================================================
# 报告生成
# ============================================================================

def generate_report(df, df_stats):
    """生成分析报告"""
    print("\n[10] 生成报告...")
    
    report = []
    report.append("=" * 80)
    report.append("Pt8SnX 系列 Per-Atom 动力学分析报告")
    report.append("=" * 80)
    report.append("")
    
    report.append("1. 数据概况")
    report.append("-" * 40)
    report.append(f"   总原子记录数: {len(df)}")
    report.append(f"   Sn 含量范围: Pt8Sn{df['sn_count'].min()} - Pt8Sn{df['sn_count'].max()}")
    report.append(f"   温度范围: {df['temp'].min()}K - {df['temp'].max()}K")
    report.append(f"   元素: {', '.join(df['element'].unique())}")
    report.append("")
    
    report.append("2. 分元素统计 (所有温度平均)")
    report.append("-" * 40)
    for element in ['Pt', 'Sn']:
        elem_data = df[df['element'] == element]
        if len(elem_data) > 0:
            report.append(f"   {element}:")
            report.append(f"      原子数: {len(elem_data)}")
            report.append(f"      δ = {elem_data['delta'].mean():.4f} ± {elem_data['delta'].std():.4f}")
            report.append(f"      D = {elem_data['D'].mean():.4f} ± {elem_data['D'].std():.4f}")
            report.append(f"      高振动比例 (δ>{HIGH_LINDEMANN}): {(elem_data['delta'] > HIGH_LINDEMANN).mean()*100:.1f}%")
            report.append(f"      高扩散比例 (D>{HIGH_D}): {(elem_data['D'] > HIGH_D).mean()*100:.1f}%")
    report.append("")
    
    report.append("3. Sn含量影响 (高温 ~900K)")
    report.append("-" * 40)
    
    high_temp = df_stats['temp'].max()
    high_temp_data = df_stats[df_stats['temp'] == high_temp]
    
    for element in ['Pt', 'Sn']:
        elem_data = high_temp_data[high_temp_data['element'] == element].sort_values('sn_count')
        if len(elem_data) > 0:
            report.append(f"   {element} @ {int(high_temp)}K:")
            for _, row in elem_data.iterrows():
                report.append(f"      Pt8Sn{int(row['sn_count'])}: δ={row['delta_mean']:.4f}, D={row['D_mean']:.4f}")
    report.append("")
    
    report.append("4. 林德曼指数-扩散系数相关性分析")
    report.append("-" * 40)
    
    # 全局相关性
    valid_mask = df['delta'].notna() & df['D'].notna() & (df['D'] > 0)
    valid_df = df[valid_mask]
    if len(valid_df) > 10:
        r_global, p_global = stats.pearsonr(valid_df['delta'], valid_df['D'])
        report.append(f"   全局相关性: r = {r_global:.4f}, p = {p_global:.2e}")
    
    # 分元素相关性
    report.append("   分元素相关性:")
    for element in ['Pt', 'Sn']:
        elem_data = valid_df[valid_df['element'] == element]
        if len(elem_data) > 10:
            r, p = stats.pearsonr(elem_data['delta'], elem_data['D'])
            report.append(f"      {element}: r = {r:.4f}, p = {p:.2e}, n = {len(elem_data)}")
    
    # 分温度相关性
    report.append("   分温度相关性:")
    temps = sorted(valid_df['temp'].unique())
    for temp in temps:
        temp_data = valid_df[valid_df['temp'] == temp]
        if len(temp_data) > 10:
            r, p = stats.pearsonr(temp_data['delta'], temp_data['D'])
            report.append(f"      {int(temp)}K: r = {r:.4f}")
    
    report.append("")
    
    report.append("5. 关键发现")
    report.append("-" * 40)
    
    # 找哪个元素先活跃
    pt_first_high = df_stats[(df_stats['element'] == 'Pt') & (df_stats['high_delta_ratio'] > 0.5)]
    sn_first_high = df_stats[(df_stats['element'] == 'Sn') & (df_stats['high_delta_ratio'] > 0.5)]
    
    if len(pt_first_high) > 0 and len(sn_first_high) > 0:
        pt_temp = pt_first_high['temp'].min()
        sn_temp = sn_first_high['temp'].min()
        if pt_temp < sn_temp:
            report.append(f"   Pt 原子先达到高振动 (50%): {int(pt_temp)}K")
            report.append(f"   Sn 原子后达到高振动 (50%): {int(sn_temp)}K")
        else:
            report.append(f"   Sn 原子先达到高振动 (50%): {int(sn_temp)}K")
            report.append(f"   Pt 原子后达到高振动 (50%): {int(pt_temp)}K")
    
    # 添加相关性解释
    report.append("")
    report.append("   相关性解释:")
    if len(valid_df) > 10:
        r_global, _ = stats.pearsonr(valid_df['delta'], valid_df['D'])
        if r_global > 0.6:
            report.append(f"   δ 与 D 显著正相关 (r={r_global:.2f}): 高振动原子通常扩散快")
        elif r_global > 0.3:
            report.append(f"   δ 与 D 中度正相关 (r={r_global:.2f}): 振动与扩散有一定关联")
        elif r_global > 0:
            report.append(f"   δ 与 D 弱正相关 (r={r_global:.2f}): 振动与扩散关联较弱")
        else:
            report.append(f"   δ 与 D 无明显相关 (r={r_global:.2f})")
        
        # 温度依赖分析
        low_temp_r = []
        high_temp_r = []
        for temp in temps:
            temp_data = valid_df[valid_df['temp'] == temp]
            if len(temp_data) > 10:
                r, _ = stats.pearsonr(temp_data['delta'], temp_data['D'])
                if temp <= 500:
                    low_temp_r.append(r)
                else:
                    high_temp_r.append(r)
        
        if low_temp_r and high_temp_r:
            avg_low = np.mean(low_temp_r)
            avg_high = np.mean(high_temp_r)
            if avg_low > avg_high + 0.1:
                report.append(f"   低温相关性较强 (r≈{avg_low:.2f}), 高温相关性减弱 (r≈{avg_high:.2f})")
            elif avg_high > avg_low + 0.1:
                report.append(f"   高温相关性较强 (r≈{avg_high:.2f}), 低温相关性较弱 (r≈{avg_low:.2f})")
            else:
                report.append(f"   相关性随温度变化不大 (低温r≈{avg_low:.2f}, 高温r≈{avg_high:.2f})")
    
    report.append("")
    report.append("=" * 80)
    
    report_path = OUTPUT_DIR / 'analysis_report.txt'
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(report))
    print(f"  [SAVED] {report_path}")
    
    # 打印到终端
    print("\n" + "\n".join(report))


# ============================================================================
# 主函数
# ============================================================================

def main():
    parser = argparse.ArgumentParser(description='Pt8SnX Per-Atom 动力学分析')
    parser.add_argument('--interactive', action='store_true', help='显示图形')
    args = parser.parse_args()
    
    print("=" * 80)
    print("Step 7.2.1: Pt8SnX 系列 Per-Atom 动力学分析")
    print("=" * 80)
    
    # 加载数据
    df_lindemann = load_lindemann_data()
    df_msd = load_msd_data()
    
    if df_lindemann is None or df_msd is None:
        print("[ERROR] 数据加载失败")
        return
    
    # 合并数据
    df = merge_data(df_lindemann, df_msd)
    
    if len(df) == 0:
        print("[ERROR] 合并后无数据")
        return
    
    # 保存合并数据
    df.to_csv(OUTPUT_DIR / 'merged_per_atom_data.csv', index=False, encoding='utf-8-sig')
    print(f"[SAVED] {OUTPUT_DIR / 'merged_per_atom_data.csv'}")
    
    # 计算统计
    df_stats = calculate_statistics(df)
    df_stats.to_csv(OUTPUT_DIR / 'statistics.csv', index=False, encoding='utf-8-sig')
    print(f"[SAVED] {OUTPUT_DIR / 'statistics.csv'}")
    
    # 绘图
    plot_heatmaps(df_stats)
    plot_high_ratio_heatmaps(df_stats)
    plot_element_comparison(df_stats)
    plot_temperature_curves(df_stats)
    plot_delta_vs_D_scatter(df)
    
    # 相关性分析
    plot_delta_D_correlation_analysis(df)
    
    # 生成报告
    generate_report(df, df_stats)
    
    print("\n" + "=" * 80)
    print(f"分析完成！输出目录: {OUTPUT_DIR}")
    print("=" * 80)
    
    if args.interactive:
        plt.show()


if __name__ == '__main__':
    main()
