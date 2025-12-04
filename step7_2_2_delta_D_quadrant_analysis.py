#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 7.2.2: δ-D 四象限分析图
============================

分析 Lindemann Index (δ) 与 Diffusion Coefficient (D) 的关系，
按四象限分类原子动力学行为。

输出:
- delta_D_quadrant_300K.png  (独立图, 20:15)
- delta_D_quadrant_600K.png  (独立图, 20:15)
- delta_D_quadrant_900K.png  (独立图, 20:15)
- delta_D_quadrant_combined.png (合并图, 共享坐标轴)

分区定义 (标准象限划分, 右上角开始逆时针):
- I:   高δ高D (Active) - 右上
- II:  低δ高D (Diffusing) - 左上
- III: 低δ低D (Stable) - 左下
- IV:  高δ低D (Vibrating) - 右下

================================================================================
命令行参数
================================================================================
--only-series    : 只绘制指定系列（逗号分隔）
                   可选系列:
                   - pt8snx   : Pt=8 的系列 (Pt8Sn0~Pt8Sn10，无氧)
                   - pt6snx   : Pt=6 的系列 (Pt6Sn1~Pt6Sn9，无氧)
                   - sum8     : Pt+Sn=8 (如Pt7Sn1,Pt6Sn2,...，无氧)
                   例: --only-series sum8
                   例: --only-series pt8snx,pt6snx

--exclude, -e    : 排除指定组分，格式 "(pt,sn)" 或 "(pt1,sn1);(pt2,sn2)"
                   例: --exclude "(8,0)" 排除 Pt8Sn0
                   例: --exclude "(8,0);(6,0)" 排除 Pt8Sn0 和 Pt6Sn0

--temps, -t      : 指定温度列表（逗号分隔）
                   例: --temps 300,600,900
                   默认: 300,600,900

================================================================================
使用示例
================================================================================

# 默认: pt8snx，不排除
python step7_2_2_delta_D_quadrant_analysis.py

# pt8snx 系列，排除 Pt8Sn0
python step7_2_2_delta_D_quadrant_analysis.py --only-series pt8snx --exclude "(8,0)"

# sum8 + pt8snx + pt6snx，排除 Pt8Sn0 和 Pt6Sn0
python step7_2_2_delta_D_quadrant_analysis.py --only-series sum8,pt8snx,pt6snx --exclude "(8,0);(6,0)"

# 指定温度
python step7_2_2_delta_D_quadrant_analysis.py --temps 200,500,800,1100

================================================================================

Author: AI Assistant
Date: 2024-12
"""

import argparse
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ============== 配置 ==============
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'

# 原始数据路径
LINDEMANN_DIR = Path('data/lindemann/per-atoms')
MSD_DIR = Path('data/gmx_msd/per-atom/collected_gmx_per_atom_msd')
OUTPUT_DIR = Path('results/per_atom_quadrant_analysis')

# 单位转换: cm²/s → Å²/fs
# 1 cm² = 1e16 Å², 1 s = 1e15 fs
# 因此 1 cm²/s = 1e16/1e15 = 10 Å²/fs
CM2_S_TO_A2_FS = 10.0

# 阈值 (Å²/fs 单位)
HIGH_LINDEMANN = 0.10
HIGH_D = 1e-6  # Å²/fs (对应 1e-7 cm²/s)

# 简洁配色 (灰色系) - 标准象限划分 (右上角开始逆时针)
COLORS_MAP = {
    'I': '#000000',      # 黑色 (Active: 高δ高D) - 右上
    'II': '#A0A0A0',     # 浅灰 (Diffusing: 低δ高D) - 左上
    'III': '#808080',    # 灰色 (Stable: 低δ低D) - 左下
    'IV': '#404040'      # 深灰 (Vibrating: 高δ低D) - 右下
}


def parse_structure_name(structure_name):
    """
    解析结构名称，提取 Pt, Sn, O 数量
    
    支持格式:
    - pt8sn0-2-best -> (8, 0, 0)
    - pt8sn5-1-best -> (8, 5, 0)
    - pt6sn8o4 -> (6, 8, 4)
    
    Returns: (pt_count, sn_count, o_count) or None
    """
    name_lower = structure_name.lower()
    
    pt_match = re.search(r'pt(\d+)', name_lower)
    sn_match = re.search(r'sn(\d+)', name_lower)
    o_match = re.search(r'o(\d+)', name_lower)
    
    pt_num = int(pt_match.group(1)) if pt_match else 0
    sn_num = int(sn_match.group(1)) if sn_match else 0
    o_num = int(o_match.group(1)) if o_match else 0
    
    if pt_num > 0 or sn_num > 0:
        return (pt_num, sn_num, o_num)
    
    return None


def classify_structure(structure_name):
    """
    分类结构
    
    Returns: dict with classification flags
    - pt8snx: Pt=8 系列 (无氧)
    - pt6snx: Pt=6 系列 (无氧)
    - sum8: Pt+Sn=8 (无氧)
    - has_oxide: 含氧
    
    注意: 含氧体系不归入 sum8/pt8snx/pt6snx
    """
    comp = parse_structure_name(structure_name)
    if not comp:
        return {'primary': 'other'}
    
    pt, sn, o = comp
    result = {}
    
    # 含氧体系单独分类，不归入其他系列
    if o > 0:
        result['has_oxide'] = True
        result[f'o{o}'] = True
        result['primary'] = f'o{o}'
        return result
    
    result['has_oxide'] = False
    
    # 无氧体系
    # 按 Pt 原子数分类
    if pt == 8:
        result['pt8snx'] = True
    if pt == 6:
        result['pt6snx'] = True
    
    # 按 Pt+Sn 总数分类 (只有无氧体系)
    total = pt + sn
    if total == 8:
        result['sum8'] = True
    
    # 确定主分类
    if total == 8:
        result['primary'] = 'sum8'
    elif pt == 8:
        result['primary'] = 'pt8snx'
    elif pt == 6:
        result['primary'] = 'pt6snx'
    else:
        result['primary'] = 'other'
    
    return result


def match_series(classification, target_series_list):
    """
    检查分类是否匹配目标系列列表中的任一个
    
    Args:
        classification: classify_structure 返回的分类字典
        target_series_list: 目标系列列表，如 ['sum8', 'pt8snx']
    
    Returns: True if matches any
    """
    for series in target_series_list:
        if classification.get(series, False):
            return True
    return False


def extract_full_run_key(path):
    """从完整路径提取唯一的 run key"""
    if pd.isna(path):
        return None
    parts = str(path).split('/')
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
    
    # 重命名列
    df = df.rename(columns={
        'lindemann_index': 'delta',
        'element': 'element'
    })
    
    print(f"  [OK] Lindemann 数据: {len(df)} records")
    print(f"       结构数: {df['structure'].nunique()}")
    
    return df[['structure', 'temp', 'atom_id', 'element', 'delta', 'full_run_key']]


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
    
    # 重命名列
    df = df.rename(columns={
        '结构': 'structure',
        '温度(K)': 'temp',
        '元素': 'element',
        'D(1e-5 cm²/s)': 'D'
    })
    
    # D 值换算：原单位是 10⁻⁵ cm²/s，转换为真实值 cm²/s
    df['D'] = df['D'] * 1e-5
    
    # 结构名转小写
    df['structure'] = df['structure'].str.lower()
    
    print(f"  [OK] MSD 数据: {len(df)} records")
    print(f"       结构数: {df['structure'].nunique()}")
    
    return df[['structure', 'temp', 'atom_id', 'element', 'D', 'full_run_key']]


def load_data():
    """加载并合并 Lindemann 和 MSD 数据"""
    print("=" * 60)
    print("加载原始数据...")
    print("=" * 60)
    
    # 加载两个数据源
    df_lind = load_lindemann_data()
    df_msd = load_msd_data()
    
    if df_lind is None or df_msd is None:
        raise FileNotFoundError("无法加载数据文件")
    
    # 合并数据
    print("\n[3] 合并 Lindemann 和 MSD 数据...")
    df = pd.merge(
        df_lind,
        df_msd[['full_run_key', 'atom_id', 'D']],
        on=['full_run_key', 'atom_id'],
        how='inner'
    )
    
    print(f"  [OK] 合并后数据: {len(df)} records")
    print(f"       结构数: {df['structure'].nunique()}")
    
    # 添加组分信息
    df['composition'] = df['structure'].apply(parse_structure_name)
    df['pt_count'] = df['composition'].apply(lambda x: x[0] if x else 0)
    df['sn_count'] = df['composition'].apply(lambda x: x[1] if x else 0)
    df['o_count'] = df['composition'].apply(lambda x: x[2] if x else 0)
    
    # 单位转换: cm²/s → Å²/fs
    df['D'] = df['D'] * CM2_S_TO_A2_FS
    
    print(f"\n原始数据: {len(df)} 条记录")
    
    return df


def filter_data(df, only_series=None, exclude_compositions=None):
    """
    根据参数筛选数据
    
    Args:
        df: 原始数据
        only_series: 只包含的系列列表，如 ['sum8', 'pt8snx']
        exclude_compositions: 排除的组分列表，如 [(8, 0), (6, 0)]
    
    Returns: 筛选后的数据
    """
    original_count = len(df)
    df = df.copy()
    
    # 添加分类信息 (缓存以提高性能)
    structure_classification = {}
    for structure in df['structure'].unique():
        structure_classification[structure] = classify_structure(structure)
    
    df['classification'] = df['structure'].map(structure_classification)
    
    # 按系列筛选
    if only_series:
        mask = df['classification'].apply(lambda c: match_series(c, only_series))
        df = df[mask]
        print(f"筛选系列 {only_series}: {original_count} -> {len(df)} 条")
    
    # 排除组分 (按 pt,sn 组合)
    if exclude_compositions:
        before = len(df)
        for pt, sn in exclude_compositions:
            mask = ~((df['pt_count'] == pt) & (df['sn_count'] == sn))
            df = df[mask]
        excluded_str = ', '.join([f'({pt},{sn})' for pt, sn in exclude_compositions])
        print(f"排除组分 {excluded_str}: {before} -> {len(df)} 条")
    
    # 显示最终包含的组分
    structures = df[['structure', 'sn_count']].drop_duplicates().sort_values('sn_count')
    print(f"\n最终包含的结构:")
    for _, row in structures.iterrows():
        comp = parse_structure_name(row['structure'])
        if comp:
            pt, sn, o = comp
            label = f"Pt{pt}Sn{sn}" + (f"O{o}" if o > 0 else "")
            print(f"  - {row['structure']}: {label}")
    
    print(f"\n最终数据: {len(df)} 条记录 (D单位: Å²/fs)")
    
    return df


def classify_atoms(df):
    """按四象限分类原子 (标准象限: 右上角开始逆时针)"""
    df = df.copy()
    # I: 右上 (高δ高D) - Active
    # II: 左上 (低δ高D) - Diffusing
    # III: 左下 (低δ低D) - Stable
    # IV: 右下 (高δ低D) - Vibrating
    df['category'] = 'III'  # 默认: Stable (低δ低D) - 左下
    df.loc[(df['delta'] > HIGH_LINDEMANN) & (df['D'] > HIGH_D), 'category'] = 'I'     # Active - 右上
    df.loc[(df['delta'] <= HIGH_LINDEMANN) & (df['D'] > HIGH_D), 'category'] = 'II'   # Diffusing - 左上
    df.loc[(df['delta'] > HIGH_LINDEMANN) & (df['D'] <= HIGH_D), 'category'] = 'IV'   # Vibrating - 右下
    return df


def plot_single_quadrant(ax, df, temp, show_ylabel=True):
    """
    在指定axes上绘制单个温度的四象限图
    """
    td = df[df['temp'] == temp].copy()
    td = classify_atoms(td)
    
    # 统计
    cat_counts = td['category'].value_counts()
    cat_pct = cat_counts / len(td) * 100
    
    # 绘制散点
    for cat in ['I', 'II', 'III', 'IV']:
        cat_data = td[td['category'] == cat]
        if len(cat_data) > 0:
            pct = cat_pct.get(cat, 0)
            ax.scatter(cat_data['delta'], cat_data['D'], 
                      c=COLORS_MAP[cat], alpha=0.5, s=40,
                      label=f'{cat} ({pct:.1f}%)')
    
    # 阈值线
    ax.axvline(HIGH_LINDEMANN, color='black', linestyle='--', linewidth=1.5, alpha=0.8)
    ax.axhline(HIGH_D, color='black', linestyle='--', linewidth=1.5, alpha=0.8)
    
    # 坐标轴标签 (34号字体)
    ax.set_xlabel('Lindemann Index', fontsize=34)
    if show_ylabel:
        ax.set_ylabel('D (Å²/fs)', fontsize=34)
    
    # 坐标轴刻度数字 (28号字体)
    ax.tick_params(axis='both', labelsize=28)
    
    # Y轴对数刻度
    ax.set_yscale('log')
    
    # 设置X轴刻度 (5个刻度)
    ax.set_xlim(0, 0.45)
    ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4])
    
    # 设置Y轴刻度 (Å²/fs 单位, 6个刻度)
    ax.set_ylim(1e-8, 1e-3)
    ax.set_yticks([1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3])
    
    # 图注 (22号字体, 无边框)
    ax.legend(loc='lower right', fontsize=22, frameon=False)
    
    return cat_pct


def generate_individual_plots(df, temperatures):
    """生成独立的四象限图 (每个20:15比例)"""
    print("\n生成独立四象限图...")
    
    for temp in temperatures:
        if temp not in df['temp'].unique():
            print(f"  [SKIP] {temp}K 无数据")
            continue
        
        # 每张图独立，20:15 比例 (宽:高)
        fig, ax = plt.subplots(figsize=(10, 7.5))
        
        cat_pct = plot_single_quadrant(ax, df, temp, show_ylabel=True)
        
        plt.tight_layout()
        
        # 保存
        output_path = OUTPUT_DIR / f'delta_D_quadrant_{temp}K.png'
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"  [SAVED] {output_path}")
        
        # 打印统计
        print(f"    {temp}K 分区: ", end='')
        for cat in ['I', 'II', 'III', 'IV']:
            pct = cat_pct.get(cat, 0)
            print(f"{cat}={pct:.1f}% ", end='')
        print()


def generate_combined_plot(df, temperatures):
    """生成合并图 (子图共享坐标轴)"""
    print("\n生成合并四象限图...")
    
    # 筛选有数据的温度
    valid_temps = [t for t in temperatures if t in df['temp'].unique()]
    if not valid_temps:
        print("  [SKIP] 无有效温度数据")
        return
    
    n_temps = len(valid_temps)
    
    # 创建共享Y轴的子图
    fig, axes = plt.subplots(1, n_temps, figsize=(8*n_temps, 9), sharey=True)
    if n_temps == 1:
        axes = [axes]
    
    for idx, temp in enumerate(valid_temps):
        ax = axes[idx]
        show_ylabel = (idx == 0)
        
        cat_pct = plot_single_quadrant(ax, df, temp, show_ylabel=show_ylabel)
        
        # 添加温度标注在图上方
        ax.text(0.5, 1.02, f'{temp} K', transform=ax.transAxes, 
               fontsize=34, fontweight='bold', ha='center', va='bottom')
        
        # 非第一个子图隐藏Y轴刻度标签
        if idx > 0:
            ax.tick_params(axis='y', labelleft=False)
    
    plt.tight_layout()
    
    # 保存
    output_path = OUTPUT_DIR / 'delta_D_quadrant_combined.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  [SAVED] {output_path}")


def generate_comparison_plot(df, temp_low=300, temp_high=900):
    """
    生成两个温度的对比图 (默认 300K vs 900K)
    
    两个子图并排，共享坐标轴，20:15 比例
    """
    print(f"\n生成 {temp_low}K vs {temp_high}K 对比图...")
    
    # 检查数据
    temps_available = df['temp'].unique()
    if temp_low not in temps_available or temp_high not in temps_available:
        print(f"  [SKIP] 需要 {temp_low}K 和 {temp_high}K 数据")
        return
    
    # 创建 2 个子图，共享Y轴，每个子图接近 20:15 比例
    fig, axes = plt.subplots(1, 2, figsize=(16, 7.5), sharey=True)
    
    for idx, temp in enumerate([temp_low, temp_high]):
        ax = axes[idx]
        show_ylabel = (idx == 0)
        
        cat_pct = plot_single_quadrant(ax, df, temp, show_ylabel=show_ylabel)
        
        # 添加温度标注在图上方
        ax.text(0.5, 1.02, f'{temp} K', transform=ax.transAxes, 
               fontsize=34, fontweight='bold', ha='center', va='bottom')
        
        # 非第一个子图隐藏Y轴刻度标签
        if idx > 0:
            ax.tick_params(axis='y', labelleft=False)
    
    plt.tight_layout()
    
    # 保存
    output_path = OUTPUT_DIR / f'delta_D_quadrant_{temp_low}K_vs_{temp_high}K.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  [SAVED] {output_path}")


def generate_statistics_report(df):
    """生成统计报告"""
    print("\n===== δ-D 四象限统计报告 =====")
    print(f"阈值: δ={HIGH_LINDEMANN}, D={HIGH_D:.0e} Å²/fs")
    print()
    
    df = classify_atoms(df)
    
    # 总体统计
    print("总体统计:")
    total_counts = df['category'].value_counts()
    total_pct = total_counts / len(df) * 100
    for cat in ['I', 'II', 'III', 'IV']:
        n = total_counts.get(cat, 0)
        pct = total_pct.get(cat, 0)
        print(f"  {cat}: {n:5d} ({pct:5.1f}%)")
    print()
    
    # 按温度统计
    print("按温度统计:")
    print("-" * 50)
    print(f"{'Temp':>6s} | {'I':>10s} | {'II':>10s} | {'III':>10s} | {'IV':>10s}")
    print("-" * 50)
    
    for temp in sorted(df['temp'].unique()):
        td = df[df['temp'] == temp]
        counts = td['category'].value_counts()
        pcts = counts / len(td) * 100
        row = f"{temp:5d}K |"
        for cat in ['I', 'II', 'III', 'IV']:
            pct = pcts.get(cat, 0)
            row += f" {pct:8.1f}% |"
        print(row)
    print("-" * 50)


def parse_args():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description='δ-D 四象限分析图',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例:
  # 默认: sum8，排除 sn=0
  python step7_2_2_delta_D_quadrant_analysis.py

  # pt8snx 系列，排除 sn=0
  python step7_2_2_delta_D_quadrant_analysis.py --only-series pt8snx --exclude "0"

  # sum8 系列，不排除任何组分
  python step7_2_2_delta_D_quadrant_analysis.py --only-series sum8 --exclude ""

  # 指定温度
  python step7_2_2_delta_D_quadrant_analysis.py --temps 200,500,800,1100

  # 只生成 300K vs 900K 对比图
  python step7_2_2_delta_D_quadrant_analysis.py --compare 300,900
'''
    )
    
    parser.add_argument('--only-series', '-s', type=str, default='pt8snx',
                        help='只绘制指定系列（逗号分隔）: sum8, pt8snx, pt6snx')
    parser.add_argument('--exclude', '-e', type=str, default='',
                        help='排除的组分，格式 "(pt,sn)" 或 "(pt1,sn1);(pt2,sn2)"')
    parser.add_argument('--temps', '-t', type=str, default='300,600,900',
                        help='温度列表（逗号分隔），如 "300,600,900"')
    parser.add_argument('--compare', '-c', type=str, default='',
                        help='生成两温度对比图，格式 "低温,高温"，如 "300,900"')
    
    return parser.parse_args()


def parse_exclude_arg(exclude_str):
    """
    解析排除参数
    
    支持格式:
    - "(8,0)" -> [(8, 0)]
    - "(8,0);(6,0)" -> [(8, 0), (6, 0)]
    
    Returns: list of (pt, sn) tuples
    """
    if not exclude_str or not exclude_str.strip():
        return []
    
    exclude_list = []
    # 按分号分割
    parts = exclude_str.split(';')
    for part in parts:
        part = part.strip()
        if not part:
            continue
        # 解析 (pt,sn) 格式
        match = re.match(r'\((\d+),(\d+)\)', part)
        if match:
            pt = int(match.group(1))
            sn = int(match.group(2))
            exclude_list.append((pt, sn))
        else:
            print(f"  [WARNING] 无法解析排除参数: {part}")
    
    return exclude_list


def main():
    """主函数"""
    args = parse_args()
    
    print("=" * 60)
    print("Step 7.2.2: δ-D 四象限分析")
    print("=" * 60)
    
    # 解析参数
    only_series = [s.strip() for s in args.only_series.split(',')] if args.only_series else None
    exclude_compositions = parse_exclude_arg(args.exclude)
    temperatures = [int(t.strip()) for t in args.temps.split(',')]
    
    # 解析对比参数
    compare_temps = None
    if args.compare:
        try:
            parts = args.compare.split(',')
            if len(parts) == 2:
                compare_temps = (int(parts[0].strip()), int(parts[1].strip()))
        except:
            print(f"  [WARNING] 无法解析 --compare 参数: {args.compare}")
    
    print(f"\n参数:")
    print(f"  --only-series: {only_series}")
    print(f"  --exclude: {exclude_compositions}")
    print(f"  --temps: {temperatures}")
    if compare_temps:
        print(f"  --compare: {compare_temps[0]}K vs {compare_temps[1]}K")
    
    # 确保输出目录存在
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # 加载数据
    df = load_data()
    
    # 筛选数据
    df = filter_data(df, only_series=only_series, exclude_compositions=exclude_compositions)
    
    if len(df) == 0:
        print("\n[ERROR] 筛选后无数据!")
        return
    
    # 生成独立图
    generate_individual_plots(df, temperatures)
    
    # 生成合并图
    generate_combined_plot(df, temperatures)
    
    # 生成300K vs 900K对比图
    if compare_temps:
        generate_comparison_plot(df, temp_low=compare_temps[0], temp_high=compare_temps[1])
    
    # 生成统计报告
    generate_statistics_report(df)
    
    print("\n" + "=" * 60)
    print("完成!")
    print("=" * 60)


if __name__ == '__main__':
    main()
