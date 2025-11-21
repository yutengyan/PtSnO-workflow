#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
Step 7.8.6: Per-Atom 林德曼指数分析
================================================================================

作者: GitHub Copilot
日期: 2025-01-12
版本: v1.0

功能概述
========
本脚本用于分析per-atom林德曼指数数据,识别高振动原子、表面熔化现象,
以及分析不同元素的空间异质性。

Per-Atom vs 整体林德曼指数
==========================
- **整体方法** (step7): 仅计算同质原子对 (Pt-Pt, Sn-Sn)
  δ_Pt = mean([δ_ij for i,j in Pt-Pt pairs])
  
- **Per-Atom方法** (step7.8.6): 每个原子对所有其他原子
  δ_i = mean([δ_ij for j in all atoms])
  
关键区别:
- Per-atom包含异质相互作用 (Pt-Sn, Pt-O)
- 可识别单个原子的振动状态
- 揭示空间异质性 (表面 vs 核心)

核心功能
========
1. **数据加载与处理**
   - 读取per_atom_master CSV文件
   - 按体系-温度-元素分组统计
   - 计算高振动原子比例

2. **元素对比分析**
   - Pt vs Sn vs O的林德曼指数分布
   - 各元素的熔化趋势
   - 元素异质性分析

3. **空间异质性识别**
   - 识别高振动原子 (δ>0.10)
   - 稳定vs活跃原子分类
   - 表面熔化检测

4. **温度依赖性分析**
   - δ vs T曲线 (按元素分组)
   - 元素熔化温度对比
   - 熔化序列识别 (哪个元素先熔化)

5. **数据可视化**
   - Per-atom δ分布直方图
   - 元素对比箱线图
   - δ vs T曲线 (按元素)
   - 高振动原子比例热力图

输入文件
========
per_atom_master_run_*.csv - Per-atom林德曼指数数据
列: 目录, 结构, 温度(K), atom_id, element, lindemann_index, 高振动原子数, 总原子数, 时间戳

输出文件
========
results/per_atom_lindemann_analysis/
├── per_atom_statistics.csv              # 统计汇总
├── element_comparison.csv               # 元素对比
├── high_vibration_atoms.csv             # 高振动原子列表
├── per_atom_analysis_report.txt         # 分析报告
├── per_atom_delta_distribution.png      # δ分布直方图
├── element_comparison_boxplot.png       # 元素对比箱线图
├── delta_vs_temperature_*.png           # δ vs T曲线
└── high_vibration_ratio_heatmap.png     # 高振动比例热力图

使用示例
========
# 分析单个per-atom文件
python step7_8_6_per_atom_lindemann_analysis.py

# 分析特定体系
python step7_8_6_per_atom_lindemann_analysis.py --filter "Sn5O2Pt4"

依赖库
======
- pandas: 数据处理
- numpy: 数值计算
- matplotlib: 绘图
- seaborn: 高级绘图

相关脚本
========
- step7_lindemann_analysis.py: 整体林德曼指数分析
- lindemann_per_atom_analysis.py: Per-atom计算工具
- run_per_atom_lindemann_batch.sh: 批处理脚本
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import re
from collections import defaultdict
import argparse

# ============================================================================
# 配置部分
# ============================================================================

# 中文显示
plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# 基础路径
BASE_DIR = Path(__file__).parent

# 输入文件 - 从服务器批处理结果目录
# 修改为实际的per-atom批处理输出目录
DATA_DIR = BASE_DIR / 'data' / 'per_atom_lindemann'

# 输出目录
OUTPUT_DIR = BASE_DIR / 'results' / 'per_atom_lindemann_analysis'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# 阈值设置
HIGH_VIBRATION_THRESHOLD = 0.10  # δ>0.10为高振动

# 颜色方案
ELEMENT_COLORS = {
    'Pt': '#1f77b4',  # 蓝色
    'Sn': '#ff7f0e',  # 橙色
    'O': '#2ca02c',   # 绿色
    'Al': '#d62728',  # 红色
}


# ============================================================================
# 辅助函数
# ============================================================================

def classify_system_series(structure_name):
    """分类体系所属系列 (与step7保持一致)"""
    if not structure_name or pd.isna(structure_name):
        return 'Other'
    
    structure_lower = str(structure_name).lower()
    
    # 优先检查含氧系列
    if re.search(r'[Oo]1(?:[^0-9]|$)', structure_lower):
        return 'O1'
    if re.search(r'[Oo]2(?:[^0-9]|$)', structure_lower):
        return 'O2'
    if re.search(r'[Oo]3(?:[^0-9]|$)', structure_lower):
        return 'O3'
    if re.search(r'[Oo]4(?:[^0-9]|$)', structure_lower):
        return 'O4'
    
    # Cv系列
    if structure_lower.startswith('cv-'):
        return 'Cv'
    
    # 提取Pt和Sn原子数
    pt_match = re.search(r'pt(\d+)', structure_lower)
    sn_match = re.search(r'sn(\d+)', structure_lower)
    
    if pt_match and sn_match:
        n_pt = int(pt_match.group(1))
        n_sn = int(sn_match.group(1))
        
        if n_pt == 8:
            return 'Pt8SnX'
        if n_pt == 6:
            return 'Pt6SnX'
        if n_pt + n_sn == 8:
            return 'Pt(8-x)SnX'
    
    if pt_match and not sn_match:
        n_pt = int(pt_match.group(1))
        if n_pt == 8:
            return 'Pt8SnX'
        if n_pt == 6:
            return 'Pt6SnX'
    
    return 'Other'


# ============================================================================
# 主分析函数
# ============================================================================

def load_per_atom_data(data_dir, structure_filter=None):
    """
    加载per-atom林德曼指数数据
    
    Args:
        data_dir: 数据目录路径
        structure_filter: 结构名筛选 (可选)
    
    Returns:
        df: DataFrame包含所有per-atom数据
    """
    print("=" * 80)
    print("Step 7.8.6: Per-Atom 林德曼指数分析")
    print("=" * 80)
    
    print("\n[*] Loading per-atom lindemann data...")
    
    # 查找所有per_atom_master文件
    per_atom_files = sorted(Path(data_dir).glob('per_atom_master_*.csv'))
    
    if len(per_atom_files) == 0:
        raise FileNotFoundError(f"No per-atom files found in {data_dir}")
    
    print(f"  [INFO] Found {len(per_atom_files)} per-atom files")
    
    # 加载并合并
    dfs = []
    for file in per_atom_files:
        df = pd.read_csv(file)
        dfs.append(df)
        print(f"    - {file.name}: {len(df)} atom records")
    
    df_all = pd.concat(dfs, ignore_index=True)
    print(f"  [OK] Merged: {len(df_all)} atom records")
    
    # 如果结构列为空，从目录中提取结构名
    if df_all['结构'].isna().all() or (df_all['结构'] == '').all():
        print("  [INFO] '结构' column is empty, extracting from '目录'...")
        def extract_structure_from_path(path):
            """从路径中提取结构名，例如 pt8sn0-2-best"""
            if pd.isna(path):
                return 'Unknown'
            # 找到最后一个包含pt/sn的文件夹
            parts = str(path).split('/')
            for part in reversed(parts):
                if 'pt' in part.lower() and ('sn' in part.lower() or part.lower().startswith('pt')):
                    return part.split('.')[0]  # 移除.rX.gpuX等后缀
            return 'Unknown'
        
        df_all['结构'] = df_all['目录'].apply(extract_structure_from_path)
        print(f"    - Extracted {df_all['结构'].nunique()} unique structures")
    
    # 添加系列分类
    df_all['系列'] = df_all['结构'].apply(classify_system_series)
    
    # 应用结构筛选（在提取结构名之后）
    if structure_filter:
        original_count = len(df_all)
        df_all = df_all[df_all['结构'].str.contains(structure_filter, case=False, na=False)]
        print(f"  [FILTER] Structure filter '{structure_filter}': {original_count} → {len(df_all)} records")
    
    # 统计信息
    print(f"\n[*] Data Summary:")
    print(f"  - Total atoms: {len(df_all)}")
    print(f"  - Unique structures: {df_all['结构'].nunique()}")
    print(f"  - Temperature range: {df_all['温度(K)'].min():.0f} - {df_all['温度(K)'].max():.0f} K")
    print(f"  - Elements: {sorted(df_all['element'].unique())}")
    print(f"  - δ range: {df_all['lindemann_index'].min():.6f} - {df_all['lindemann_index'].max():.6f}")
    print(f"  - High vibration atoms (δ>0.10): {(df_all['lindemann_index'] > HIGH_VIBRATION_THRESHOLD).sum()} ({(df_all['lindemann_index'] > HIGH_VIBRATION_THRESHOLD).mean()*100:.1f}%)")
    
    # 元素统计
    print(f"\n[*] Element Statistics:")
    elem_stats = df_all.groupby('element').agg({
        'lindemann_index': ['count', 'mean', 'std', 'min', 'max']
    })
    print(elem_stats)
    
    # 系列统计
    print(f"\n[*] Series Statistics:")
    series_counts = df_all.groupby('系列').size().sort_values(ascending=False)
    for series, count in series_counts.items():
        print(f"    - {series}: {count} atoms")
    
    return df_all


def calculate_statistics(df):
    """
    计算统计数据
    
    按体系-温度-元素分组,计算:
    - 平均δ
    - 标准差
    - 高振动原子比例
    - 原子数
    """
    print("\n[*] Calculating statistics...")
    
    stats_list = []
    
    # 按体系-温度-元素分组
    for (struct, temp, elem), group in df.groupby(['结构', '温度(K)', 'element']):
        n_atoms = len(group)
        mean_delta = group['lindemann_index'].mean()
        std_delta = group['lindemann_index'].std()
        min_delta = group['lindemann_index'].min()
        max_delta = group['lindemann_index'].max()
        
        n_high = (group['lindemann_index'] > HIGH_VIBRATION_THRESHOLD).sum()
        high_ratio = n_high / n_atoms * 100
        
        # 获取系列信息
        series = group['系列'].iloc[0]
        
        stats_list.append({
            '结构': struct,
            '系列': series,
            '温度(K)': temp,
            '元素': elem,
            '原子数': n_atoms,
            '平均δ': mean_delta,
            '标准差': std_delta,
            '最小δ': min_delta,
            '最大δ': max_delta,
            '高振动原子数': n_high,
            '高振动比例(%)': high_ratio
        })
    
    df_stats = pd.DataFrame(stats_list)
    
    # 保存
    output_file = OUTPUT_DIR / 'per_atom_statistics.csv'
    df_stats.to_csv(output_file, index=False, encoding='utf-8-sig')
    print(f"  [OK] Saved: {output_file.name}")
    
    return df_stats


def analyze_element_comparison(df):
    """
    元素对比分析
    
    对比Pt, Sn, O的林德曼指数分布
    """
    print("\n[*] Analyzing element comparison...")
    
    comparison_list = []
    
    # 按体系-温度分组,对比各元素
    for (struct, temp), group in df.groupby(['结构', '温度(K)']):
        series = group['系列'].iloc[0]
        
        # 各元素统计
        elem_data = {}
        for elem in ['Pt', 'Sn', 'O']:
            elem_group = group[group['element'] == elem]
            if len(elem_group) > 0:
                elem_data[f'{elem}_mean'] = elem_group['lindemann_index'].mean()
                elem_data[f'{elem}_count'] = len(elem_group)
                elem_data[f'{elem}_high_ratio'] = (elem_group['lindemann_index'] > HIGH_VIBRATION_THRESHOLD).mean() * 100
            else:
                elem_data[f'{elem}_mean'] = np.nan
                elem_data[f'{elem}_count'] = 0
                elem_data[f'{elem}_high_ratio'] = 0
        
        comparison_list.append({
            '结构': struct,
            '系列': series,
            '温度(K)': temp,
            **elem_data
        })
    
    df_comp = pd.DataFrame(comparison_list)
    
    # 保存
    output_file = OUTPUT_DIR / 'element_comparison.csv'
    df_comp.to_csv(output_file, index=False, encoding='utf-8-sig')
    print(f"  [OK] Saved: {output_file.name}")
    
    return df_comp


def identify_high_vibration_atoms(df):
    """
    识别并记录所有高振动原子
    """
    print("\n[*] Identifying high vibration atoms...")
    
    df_high = df[df['lindemann_index'] > HIGH_VIBRATION_THRESHOLD].copy()
    
    print(f"  [INFO] High vibration atoms: {len(df_high)} / {len(df)} ({len(df_high)/len(df)*100:.1f}%)")
    
    # 按元素统计
    elem_high = df_high.groupby('element').size()
    print(f"  [INFO] By element:")
    for elem, count in elem_high.items():
        total = len(df[df['element'] == elem])
        print(f"    - {elem}: {count} / {total} ({count/total*100:.1f}%)")
    
    # 保存
    output_file = OUTPUT_DIR / 'high_vibration_atoms.csv'
    df_high.to_csv(output_file, index=False, encoding='utf-8-sig')
    print(f"  [OK] Saved: {output_file.name}")
    
    return df_high


def plot_delta_distribution(df):
    """
    绘制per-atom δ分布直方图
    """
    print("\n[*] Plotting delta distribution...")
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Per-Atom 林德曼指数分布', fontsize=16, fontweight='bold')
    
    # 子图1: 总体分布 (按元素)
    ax1 = axes[0, 0]
    elements = sorted(df['element'].unique())
    for elem in elements:
        df_elem = df[df['element'] == elem]
        ax1.hist(df_elem['lindemann_index'], bins=50, alpha=0.5, 
                label=f'{elem} (n={len(df_elem)})',
                color=ELEMENT_COLORS.get(elem, 'gray'))
    
    ax1.axvline(HIGH_VIBRATION_THRESHOLD, color='red', linestyle='--', 
               linewidth=2, label=f'阈值 (δ={HIGH_VIBRATION_THRESHOLD})')
    ax1.set_xlabel('林德曼指数 δ', fontsize=12)
    ax1.set_ylabel('原子数', fontsize=12)
    ax1.set_title('总体分布 (按元素)', fontsize=13, fontweight='bold')
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    
    # 子图2: 累积分布函数 (CDF)
    ax2 = axes[0, 1]
    for elem in elements:
        df_elem = df[df['element'] == elem]
        sorted_delta = np.sort(df_elem['lindemann_index'])
        cdf = np.arange(1, len(sorted_delta) + 1) / len(sorted_delta)
        ax2.plot(sorted_delta, cdf, label=elem, linewidth=2,
                color=ELEMENT_COLORS.get(elem, 'gray'))
    
    ax2.axvline(HIGH_VIBRATION_THRESHOLD, color='red', linestyle='--', 
               linewidth=2, alpha=0.7)
    ax2.set_xlabel('林德曼指数 δ', fontsize=12)
    ax2.set_ylabel('累积概率', fontsize=12)
    ax2.set_title('累积分布函数 (CDF)', fontsize=13, fontweight='bold')
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)
    
    # 子图3: 核密度估计 (KDE)
    ax3 = axes[1, 0]
    for elem in elements:
        df_elem = df[df['element'] == elem]
        df_elem['lindemann_index'].plot.kde(ax=ax3, label=elem, linewidth=2,
                                           color=ELEMENT_COLORS.get(elem, 'gray'))
    
    ax3.axvline(HIGH_VIBRATION_THRESHOLD, color='red', linestyle='--', 
               linewidth=2, alpha=0.7, label=f'阈值 (δ={HIGH_VIBRATION_THRESHOLD})')
    ax3.set_xlabel('林德曼指数 δ', fontsize=12)
    ax3.set_ylabel('密度', fontsize=12)
    ax3.set_title('核密度估计 (KDE)', fontsize=13, fontweight='bold')
    ax3.legend(fontsize=10)
    ax3.grid(True, alpha=0.3)
    ax3.set_xlim(0, df['lindemann_index'].max() * 1.1)
    
    # 子图4: 统计汇总表
    ax4 = axes[1, 1]
    ax4.axis('off')
    
    # 创建统计表
    stats_data = []
    for elem in elements:
        df_elem = df[df['element'] == elem]
        stats_data.append([
            elem,
            len(df_elem),
            f"{df_elem['lindemann_index'].mean():.4f}",
            f"{df_elem['lindemann_index'].std():.4f}",
            f"{(df_elem['lindemann_index'] > HIGH_VIBRATION_THRESHOLD).sum()}",
            f"{(df_elem['lindemann_index'] > HIGH_VIBRATION_THRESHOLD).mean()*100:.1f}%"
        ])
    
    table = ax4.table(cellText=stats_data,
                     colLabels=['元素', '原子数', '平均δ', '标准差', '高振动数', '高振动比例'],
                     cellLoc='center',
                     loc='center',
                     bbox=[0, 0.3, 1, 0.6])
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 2)
    
    # 设置表头样式
    for i in range(6):
        table[(0, i)].set_facecolor('#4472C4')
        table[(0, i)].set_text_props(weight='bold', color='white')
    
    # 设置行颜色
    for i in range(1, len(stats_data) + 1):
        for j in range(6):
            if i % 2 == 0:
                table[(i, j)].set_facecolor('#E7E6E6')
    
    ax4.set_title('统计汇总', fontsize=13, fontweight='bold', pad=20)
    
    plt.tight_layout()
    filename = OUTPUT_DIR / 'per_atom_delta_distribution.png'
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  [OK] Saved: {filename.name}")


def plot_element_comparison_boxplot(df):
    """
    绘制元素对比箱线图
    """
    print("\n[*] Plotting element comparison boxplot...")
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    fig.suptitle('元素林德曼指数对比', fontsize=16, fontweight='bold')
    
    # 子图1: 总体箱线图
    ax1 = axes[0]
    elements = sorted(df['element'].unique())
    data_list = [df[df['element'] == elem]['lindemann_index'].values for elem in elements]
    
    bp = ax1.boxplot(data_list, labels=elements, patch_artist=True,
                    widths=0.6, showmeans=True,
                    meanprops=dict(marker='D', markerfacecolor='red', markersize=8))
    
    # 设置颜色
    for patch, elem in zip(bp['boxes'], elements):
        patch.set_facecolor(ELEMENT_COLORS.get(elem, 'gray'))
        patch.set_alpha(0.6)
    
    ax1.axhline(HIGH_VIBRATION_THRESHOLD, color='red', linestyle='--', 
               linewidth=2, alpha=0.7, label=f'阈值 (δ={HIGH_VIBRATION_THRESHOLD})')
    ax1.set_ylabel('林德曼指数 δ', fontsize=12)
    ax1.set_title('元素对比箱线图', fontsize=13, fontweight='bold')
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3, axis='y')
    
    # 子图2: 小提琴图
    ax2 = axes[1]
    df_plot = df[['element', 'lindemann_index']].copy()
    sns.violinplot(data=df_plot, x='element', y='lindemann_index', ax=ax2,
                  palette=ELEMENT_COLORS, inner='box')
    
    ax2.axhline(HIGH_VIBRATION_THRESHOLD, color='red', linestyle='--', 
               linewidth=2, alpha=0.7)
    ax2.set_xlabel('元素', fontsize=12)
    ax2.set_ylabel('林德曼指数 δ', fontsize=12)
    ax2.set_title('元素分布小提琴图', fontsize=13, fontweight='bold')
    ax2.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    filename = OUTPUT_DIR / 'element_comparison_boxplot.png'
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  [OK] Saved: {filename.name}")


def plot_delta_vs_temperature(df, df_stats):
    """
    绘制δ vs T曲线 (按元素和系列)
    """
    print("\n[*] Plotting delta vs temperature...")
    
    # 按系列分组绘制
    series_list = sorted(df['系列'].unique())
    
    for series in series_list:
        df_series = df[df['系列'] == series].copy()
        
        if len(df_series) == 0:
            continue
        
        structures = sorted(df_series['结构'].unique())
        
        if len(structures) == 0:
            continue
        
        # 为每个结构创建子图
        n_structures = len(structures)
        n_cols = min(3, n_structures)
        n_rows = (n_structures + n_cols - 1) // n_cols
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(6*n_cols, 5*n_rows))
        if n_structures == 1:
            axes = np.array([axes])
        axes = axes.flatten()
        
        fig.suptitle(f'{series}系列: Per-Atom 林德曼指数 vs 温度', 
                    fontsize=16, fontweight='bold')
        
        for idx, struct in enumerate(structures):
            ax = axes[idx]
            df_struct = df_series[df_series['结构'] == struct]
            
            # 按元素绘制
            elements = sorted(df_struct['element'].unique())
            
            for elem in elements:
                df_elem = df_struct[df_struct['element'] == elem]
                
                # 计算每个温度点的统计量
                temp_stats = df_elem.groupby('温度(K)').agg({
                    'lindemann_index': ['mean', 'std', 'count']
                }).reset_index()
                temp_stats.columns = ['温度', 'mean', 'std', 'count']
                temp_stats = temp_stats.sort_values('温度')
                
                # 绘制均值曲线
                ax.plot(temp_stats['温度'], temp_stats['mean'],
                       marker='o', linewidth=2, markersize=6,
                       color=ELEMENT_COLORS.get(elem, 'gray'),
                       label=f'{elem} (n={temp_stats["count"].iloc[0]:.0f})')
                
                # 添加误差范围
                ax.fill_between(temp_stats['温度'],
                               temp_stats['mean'] - temp_stats['std'],
                               temp_stats['mean'] + temp_stats['std'],
                               alpha=0.2,
                               color=ELEMENT_COLORS.get(elem, 'gray'))
            
            # 添加阈值线
            ax.axhline(HIGH_VIBRATION_THRESHOLD, color='red', linestyle='--',
                      linewidth=2, alpha=0.5, label=f'阈值 (δ={HIGH_VIBRATION_THRESHOLD})')
            
            ax.set_xlabel('温度 (K)', fontsize=11)
            ax.set_ylabel('林德曼指数 δ', fontsize=11)
            ax.set_title(struct, fontsize=12, fontweight='bold')
            ax.legend(fontsize=9, loc='best')
            ax.grid(True, alpha=0.3)
        
        # 隐藏多余的子图
        for idx in range(n_structures, len(axes)):
            axes[idx].axis('off')
        
        plt.tight_layout()
        filename = OUTPUT_DIR / f'delta_vs_temperature_{series}.png'
        plt.savefig(filename, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  [OK] Saved: {filename.name}")


def plot_high_vibration_ratio_heatmap(df_stats):
    """
    绘制高振动原子比例热力图
    """
    print("\n[*] Plotting high vibration ratio heatmap...")
    
    # 按系列分组
    series_list = sorted(df_stats['系列'].unique())
    
    for series in series_list:
        df_series = df_stats[df_stats['系列'] == series].copy()
        
        if len(df_series) == 0:
            continue
        
        # 为每个元素创建热力图
        elements = sorted(df_series['元素'].unique())
        
        n_elements = len(elements)
        fig, axes = plt.subplots(1, n_elements, figsize=(6*n_elements, 8))
        if n_elements == 1:
            axes = [axes]
        
        fig.suptitle(f'{series}系列: 高振动原子比例热力图', 
                    fontsize=16, fontweight='bold')
        
        for idx, elem in enumerate(elements):
            ax = axes[idx]
            df_elem = df_series[df_series['元素'] == elem]
            
            # 创建透视表
            pivot_data = df_elem.pivot_table(
                index='温度(K)',
                columns='结构',
                values='高振动比例(%)',
                aggfunc='mean'
            )
            
            if pivot_data.empty:
                ax.text(0.5, 0.5, f'No data for {elem}',
                       ha='center', va='center', transform=ax.transAxes)
                continue
            
            # 绘制热力图
            im = ax.imshow(pivot_data.values, aspect='auto', cmap='RdYlBu_r',
                          interpolation='nearest', origin='lower',
                          vmin=0, vmax=100)
            
            # 设置坐标轴
            ax.set_xticks(np.arange(len(pivot_data.columns)))
            ax.set_yticks(np.arange(len(pivot_data.index)))
            ax.set_xticklabels(pivot_data.columns, fontsize=9)
            ax.set_yticklabels([f'{int(t)}K' for t in pivot_data.index], fontsize=9)
            
            plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
            
            # 标注数值
            for i in range(len(pivot_data.index)):
                for j in range(len(pivot_data.columns)):
                    value = pivot_data.values[i, j]
                    if not np.isnan(value):
                        text_color = 'white' if value > 50 else 'black'
                        ax.text(j, i, f'{value:.0f}%',
                               ha="center", va="center", color=text_color,
                               fontsize=8)
            
            ax.set_title(f'{elem} 元素', fontsize=12, fontweight='bold')
            ax.set_xlabel('结构', fontsize=11)
            ax.set_ylabel('温度', fontsize=11)
            
            # 添加颜色条
            cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
            cbar.set_label('高振动原子比例 (%)', fontsize=10)
        
        plt.tight_layout()
        filename = OUTPUT_DIR / f'high_vibration_ratio_heatmap_{series}.png'
        plt.savefig(filename, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  [OK] Saved: {filename.name}")


def generate_report(df, df_stats, df_comp, df_high):
    """
    生成分析报告
    """
    print("\n[*] Generating analysis report...")
    
    report_lines = []
    report_lines.append("=" * 80)
    report_lines.append("Per-Atom 林德曼指数分析报告")
    report_lines.append("=" * 80)
    report_lines.append(f"生成时间: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report_lines.append("")
    
    # 1. 数据概览
    report_lines.append("[1] 数据概览")
    report_lines.append("-" * 80)
    report_lines.append(f"总原子数: {len(df)}")
    report_lines.append(f"体系数量: {df['结构'].nunique()}")
    report_lines.append(f"温度范围: {df['温度(K)'].min():.0f} - {df['温度(K)'].max():.0f} K")
    report_lines.append(f"δ范围: {df['lindemann_index'].min():.6f} - {df['lindemann_index'].max():.6f}")
    report_lines.append(f"高振动原子: {len(df_high)} / {len(df)} ({len(df_high)/len(df)*100:.1f}%)")
    report_lines.append("")
    
    # 2. 元素统计
    report_lines.append("[2] 元素统计")
    report_lines.append("-" * 80)
    elem_stats = df.groupby('element').agg({
        'lindemann_index': ['count', 'mean', 'std', 'min', 'max']
    })
    for elem in sorted(df['element'].unique()):
        stats = elem_stats.loc[elem]
        count = int(stats[('lindemann_index', 'count')])
        mean = stats[('lindemann_index', 'mean')]
        std = stats[('lindemann_index', 'std')]
        min_val = stats[('lindemann_index', 'min')]
        max_val = stats[('lindemann_index', 'max')]
        
        high_count = len(df_high[df_high['element'] == elem])
        high_ratio = high_count / count * 100
        
        report_lines.append(f"{elem}:")
        report_lines.append(f"  - 原子数: {count}")
        report_lines.append(f"  - δ平均: {mean:.6f} ± {std:.6f}")
        report_lines.append(f"  - δ范围: [{min_val:.6f}, {max_val:.6f}]")
        report_lines.append(f"  - 高振动: {high_count} ({high_ratio:.1f}%)")
    report_lines.append("")
    
    # 3. 系列统计
    report_lines.append("[3] 系列统计")
    report_lines.append("-" * 80)
    series_stats = df.groupby('系列').agg({
        '结构': 'nunique',
        'lindemann_index': ['mean', 'std']
    })
    for series in sorted(df['系列'].unique()):
        stats = series_stats.loc[series]
        n_struct = int(stats[('结构', 'nunique')])
        mean = stats[('lindemann_index', 'mean')]
        std = stats[('lindemann_index', 'std')]
        
        report_lines.append(f"{series}: {n_struct} 结构, δ平均={mean:.6f}±{std:.6f}")
    report_lines.append("")
    
    # 4. 高振动原子TOP10
    report_lines.append("[4] 高振动原子 TOP 10")
    report_lines.append("-" * 80)
    df_top10 = df_high.nlargest(10, 'lindemann_index')
    report_lines.append(f"{'结构':<18} {'温度(K)':<8} {'原子ID':<8} {'元素':<6} {'δ'}")
    report_lines.append("-" * 80)
    for _, row in df_top10.iterrows():
        report_lines.append(f"{row['结构']:<18} {row['温度(K)']:<8.0f} "
                          f"{row['atom_id']:<8} {row['element']:<6} {row['lindemann_index']:.6f}")
    report_lines.append("")
    
    # 5. 元素对比总结
    report_lines.append("[5] 元素熔化趋势对比")
    report_lines.append("-" * 80)
    
    # 按温度统计各元素的平均δ
    temp_elem_stats = df.groupby(['温度(K)', 'element'])['lindemann_index'].mean().unstack()
    
    report_lines.append("各温度点元素平均δ:")
    report_lines.append(f"{'温度(K)':<10} " + "  ".join([f"{elem:<12}" for elem in temp_elem_stats.columns]))
    report_lines.append("-" * 80)
    for temp in sorted(temp_elem_stats.index):
        line = f"{temp:<10.0f} "
        for elem in temp_elem_stats.columns:
            val = temp_elem_stats.loc[temp, elem]
            if pd.notna(val):
                marker = " ⚠️" if val > HIGH_VIBRATION_THRESHOLD else ""
                line += f"{val:.6f}{marker:<6} "
            else:
                line += f"{'N/A':<12} "
        report_lines.append(line)
    
    report_lines.append("")
    report_lines.append("=" * 80)
    
    # 保存报告
    report_file = OUTPUT_DIR / 'per_atom_analysis_report.txt'
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write('\n'.join(report_lines))
    
    print(f"  [OK] Saved: {report_file.name}")
    
    # 打印到控制台
    print("\n" + '\n'.join(report_lines))


def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description='Step 7.8.6: Per-Atom 林德曼指数分析'
    )
    parser.add_argument(
        '--data-dir',
        type=str,
        default=str(DATA_DIR),
        help='Per-atom数据目录路径'
    )
    parser.add_argument(
        '--filter',
        type=str,
        default=None,
        help='结构名筛选 (支持正则表达式)'
    )
    args = parser.parse_args()
    
    # 1. 加载数据
    df = load_per_atom_data(args.data_dir, args.filter)
    
    # 2. 计算统计
    df_stats = calculate_statistics(df)
    
    # 3. 元素对比分析
    df_comp = analyze_element_comparison(df)
    
    # 4. 识别高振动原子
    df_high = identify_high_vibration_atoms(df)
    
    # 5. 绘图
    plot_delta_distribution(df)
    plot_element_comparison_boxplot(df)
    plot_delta_vs_temperature(df, df_stats)
    plot_high_vibration_ratio_heatmap(df_stats)
    
    # 6. 生成报告
    generate_report(df, df_stats, df_comp, df_high)
    
    print("\n" + "=" * 80)
    print(f"[DONE] Analysis complete! Results saved to: {OUTPUT_DIR}")
    print("=" * 80)


if __name__ == '__main__':
    main()
