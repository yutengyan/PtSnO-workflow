#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Arrhenius 分析：从每原子扩散系数中提取迁移能（Ea）和指前因子（D0）。

用法示例:
    # 默认运行（所有结构）
    python step7_2_2_arrhenius_analysis.py

    # 只分析 sum8 和 pt8snx 系列，排除 Pt8Sn0
    python step7_2_2_arrhenius_analysis.py --only-series sum8,pt8snx --exclude "(8,0)"

    # 只分析 Pt 和 Sn 元素（排除 all）
    python step7_2_2_arrhenius_analysis.py --only-series sum8,pt8snx --exclude "(8,0)" --elements Pt,Sn

输出:
    results/arrhenius/arrhenius_per_structure_element.csv
    results/arrhenius/arrhenius_per_structure.csv
    results/arrhenius/plots/<structure>_<element>_arrhenius.png

策略:
  - 从 data/gmx_msd/per-atom/collected_gmx_per_atom_msd/ 中读取 per-atom diffusion CSV
  - 将 D（原文件单位为 1e-5 cm^2/s）转换为 Å^2/fs（仓库中使用转换常量）
  - 按 --only-series / --exclude 筛选结构
  - 对每个 structure (和 structure+element) 计算温度序列的平均 D
  - 对 ln(D) vs 1/T 做线性拟合, 得到 Ea (eV) 和 D0

作者: 自动生成
日期: 2025
"""

import argparse
from pathlib import Path
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math

# 常量
K_B_EV_PER_K = 8.617333262145e-5  # eV/K
CM2_S_TO_A2_FS = 10.0  # 仓库中约定的转换常量 (1 cm^2/s -> 10 Å^2/fs)

# 数据路径
MSD_DIR = Path('data/gmx_msd/per-atom/collected_gmx_per_atom_msd')
CLUSTERING_DIR = Path('results/step6_1_clustering')
OUTPUT_DIR = Path('results/arrhenius')
PLOTS_DIR = OUTPUT_DIR / 'plots'

# ============== 结构分类函数 (复用自 step7_2_2_delta_D_quadrant_analysis.py) ==============

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
    """
    comp = parse_structure_name(structure_name)
    if not comp:
        return {'primary': 'other'}
    
    pt, sn, o = comp
    result = {}
    
    # 含氧体系单独分类
    if o > 0:
        result['has_oxide'] = True
        result[f'o{o}'] = True
        result['primary'] = f'o{o}'
        return result
    
    result['has_oxide'] = False
    
    # 无氧体系
    if pt == 8:
        result['pt8snx'] = True
    if pt == 6:
        result['pt6snx'] = True
    
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
    """检查分类是否匹配目标系列列表中的任一个"""
    for series in target_series_list:
        if classification.get(series, False):
            return True
    return False


def filter_data(df, only_series=None, exclude_compositions=None, exclude_structures=None):
    """
    根据参数筛选数据
    
    Args:
        df: 原始数据 (必须有 'structure' 列)
        only_series: 只包含的系列列表，如 ['sum8', 'pt8snx']
        exclude_compositions: 排除的组分列表，如 [(8, 0), (6, 0)]
        exclude_structures: 排除的结构名称列表，如 ['pt8sn5-1-best']
    
    Returns: 筛选后的数据
    """
    original_count = len(df)
    df = df.copy()
    
    # 解析组分
    df['composition'] = df['structure'].apply(parse_structure_name)
    df['pt_count'] = df['composition'].apply(lambda x: x[0] if x else 0)
    df['sn_count'] = df['composition'].apply(lambda x: x[1] if x else 0)
    df['o_count'] = df['composition'].apply(lambda x: x[2] if x else 0)
    
    # 添加分类信息
    structure_classification = {}
    for structure in df['structure'].unique():
        structure_classification[structure] = classify_structure(structure)
    df['classification'] = df['structure'].map(structure_classification)
    
    # 按系列筛选
    if only_series:
        mask = df['classification'].apply(lambda c: match_series(c, only_series))
        df = df[mask]
        print(f"筛选系列 {only_series}: {original_count} -> {len(df)} 条")
    
    # 排除组分
    if exclude_compositions:
        before = len(df)
        for pt, sn in exclude_compositions:
            mask = ~((df['pt_count'] == pt) & (df['sn_count'] == sn))
            df = df[mask]
        excluded_str = ', '.join([f'({pt},{sn})' for pt, sn in exclude_compositions])
        print(f"排除组分 {excluded_str}: {before} -> {len(df)} 条")

    # 排除特定结构
    if exclude_structures:
        before = len(df)
        df = df[~df['structure'].isin(exclude_structures)]
        print(f"排除结构 {exclude_structures}: {before} -> {len(df)} 条")
    
    # 显示最终包含的结构
    structures = df[['structure', 'sn_count']].drop_duplicates().sort_values('sn_count')
    print(f"\n最终包含的结构:")
    for _, row in structures.iterrows():
        comp = parse_structure_name(row['structure'])
        if comp:
            pt, sn, o = comp
            label = f"Pt{pt}Sn{sn}" + (f"O{o}" if o > 0 else "")
            print(f"  - {row['structure']}: {label}")
    
    print(f"\n最终数据: {len(df)} 条记录")
    
    return df


def parse_exclude_arg(exclude_str):
    """解析 --exclude 参数，如 '(8,0);(6,0)' -> [(8,0), (6,0)]"""
    if not exclude_str:
        return None
    result = []
    parts = exclude_str.replace(' ', '').split(';')
    for p in parts:
        match = re.match(r'\((\d+),(\d+)\)', p)
        if match:
            result.append((int(match.group(1)), int(match.group(2))))
    return result if result else None


def get_auto_partitions(structure_name):
    """
    从 step6_1 的聚类结果中自动获取分区
    
    返回: dict {phase_name: (T_min, T_max)}
    """
    # 提取基础名称，例如 pt8sn6-1-best -> pt8sn6
    base_name = structure_name.split('-')[0].lower()
    
    # 尝试匹配文件名
    matching_files = list(CLUSTERING_DIR.glob("*_kmeans_n2_clustered_data.csv"))
    
    target_file = None
    for f in matching_files:
        f_base = f.name.split('_')[0].lower()
        if f_base == base_name:
            target_file = f
            break
            
    if not target_file:
        # 尝试更宽松的匹配
        for f in matching_files:
            if base_name in f.name.lower():
                target_file = f
                break
                
    if not target_file:
        return None
        
    try:
        df_cluster = pd.read_csv(target_file)
        if 'temp' not in df_cluster.columns or 'phase_clustered' not in df_cluster.columns:
            return None
            
        # 按 phase_clustered 分组获取温度范围
        partitions = {}
        # 统计每个温度最频繁出现的 phase
        temp_phases = df_cluster.groupby('temp')['phase_clustered'].agg(lambda x: x.value_counts().idxmax())
        
        # 确保温度是连续的，或者至少按顺序处理
        unique_phases = []
        for p in temp_phases.values:
            if not unique_phases or p != unique_phases[-1]:
                unique_phases.append(p)
        
        for phase in unique_phases:
            phase_temps = temp_phases[temp_phases == phase].index
            if len(phase_temps) > 0:
                partitions[phase] = (min(phase_temps), max(phase_temps))
        
        return partitions
    except Exception as e:
        print(f"  ⚠️ 自动分区失败 ({structure_name}): {e}")
        return None

# ==============================================================================


def find_msd_file(msd_dir=MSD_DIR):
    files = sorted(msd_dir.glob('per_atom_diffusion_coefficients_*.csv'))
    if not files:
        files = sorted(msd_dir.glob('*diffusion*.csv'))
    return files[0] if files else None


def load_per_atom_D(msd_file):
    print(f"读取 MSD 文件: {msd_file}")
    df = pd.read_csv(msd_file)

    # 兼容多种列名
    col_map = {}
    if '结构' in df.columns:
        col_map['结构'] = 'structure'
    if 'structure' in df.columns:
        col_map['structure'] = 'structure'
    if '温度(K)' in df.columns:
        col_map['温度(K)'] = 'temp'
    if 'temp' in df.columns:
        col_map['temp'] = 'temp'
    if '元素' in df.columns:
        col_map['元素'] = 'element'
    if 'element' in df.columns:
        col_map['element'] = 'element'
    # D 可能以多种形式出现
    if 'D(1e-5 cm²/s)' in df.columns:
        col_map['D(1e-5 cm²/s)'] = 'D_input'
    elif 'D(1e-5 cm2/s)' in df.columns:
        col_map['D(1e-5 cm2/s)'] = 'D_input'
    elif 'D' in df.columns:
        col_map['D'] = 'D_input'

    df = df.rename(columns=col_map)

    if 'structure' not in df.columns or 'temp' not in df.columns or 'D_input' not in df.columns:
        raise ValueError('无法识别 MSD 文件中的必要列 (structure, temp, D)')

    # 规范化
    df['structure'] = df['structure'].astype(str).str.lower()
    df['temp'] = df['temp'].astype(float)

    # 将 D_input 从 1e-5 cm^2/s 转为真实 cm^2/s，如果值看起来很小/很大我们仍做乘法保护
    # 如果原列名表示 1e-5 cm^2/s，我们按乘以 1e-5
    # 不能完全可靠地识别；这里按仓库原代码做法: 输入为 1e-5 cm^2/s
    df['D_cm2_per_s'] = df['D_input'] * 1e-5

    # 转换到 Å^2/fs
    df['D'] = df['D_cm2_per_s'] * CM2_S_TO_A2_FS

    # 保留必要列
    if 'atom_id' in df.columns:
        keep = ['structure', 'temp', 'atom_id', 'element', 'D']
    else:
        # 如果没有 atom_id，就按结构/温度/元素汇总
        keep = ['structure', 'temp', 'element', 'D']

    return df[keep]


def summarize_mean_D(df):
    """按 structure,temp,element 计算平均 D，并同时计算每个 structure,temp 的整体平均（忽略元素区分）"""
    # 平均按 atom 聚合（如果存在 atom_id）
    group_cols = ['structure', 'temp', 'element'] if 'element' in df.columns else ['structure', 'temp']
    mean_by_element = df.groupby(group_cols, as_index=False)['D'].mean().rename(columns={'D':'D_mean'})

    mean_overall = df.groupby(['structure', 'temp'], as_index=False)['D'].mean().rename(columns={'D':'D_mean'})
    mean_overall['element'] = 'all'

    # 将 element 汇总结果与 overall 合并为单个表格方便统一拟合
    mean_by_element_union = pd.concat([
        mean_by_element.assign(element=lambda x: x['element'].astype(str)),
        mean_overall[['structure','temp','element','D_mean']]
    ], ignore_index=True, sort=False)

    # 对温度排序
    mean_by_element_union = mean_by_element_union.sort_values(['structure','element','temp'])
    return mean_by_element_union


def fit_arrhenius(temp_list, D_list, temp_range_name='all'):
    """对给定温度列表和对应平均 D 进行 Arrhenius 拟合。
    返回: dict with keys: n_points, slope, intercept, Ea_eV, D0, r2, range
    """
    T = np.array(temp_list, dtype=float)
    D = np.array(D_list, dtype=float)

    # 仅保留 D>0 值
    mask = np.isfinite(D) & (D > 0) & np.isfinite(T) & (T > 0)
    T = T[mask]
    D = D[mask]

    n = len(D)
    if n < 2:
        return None

    x = 1.0 / T
    y = np.log(D)

    # 线性拟合
    coeffs = np.polyfit(x, y, 1)
    slope, intercept = coeffs[0], coeffs[1]

    # 计算 R^2
    y_pred = slope * x + intercept
    ss_res = np.sum((y - y_pred)**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot != 0 else np.nan

    Ea_eV = -slope * K_B_EV_PER_K
    D0 = math.exp(intercept)

    return {
        'range': temp_range_name,
        'n_points': int(n),
        'slope': float(slope),
        'intercept': float(intercept),
        'Ea_eV': float(Ea_eV),
        'D0': float(D0),
        'r2': float(r2)
    }


def run_analysis(msd_df, min_points=2, do_plot=True, output_dir=OUTPUT_DIR, elements_filter=None, partitions=None, auto_partition=False):
    output_dir.mkdir(parents=True, exist_ok=True)
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    mean_df = summarize_mean_D(msd_df)

    results = []

    # 遍历每个 structure 与 element
    for (structure, element), group in mean_df.groupby(['structure','element']):
        temps_all = group['temp'].values
        Ds_all = group['D_mean'].values
        
        # 准备要拟合的任务列表：(temps, Ds, range_name)
        fit_tasks = [(temps_all, Ds_all, 'overall')]
        
        # 1. 自动分区
        if auto_partition:
            auto_parts = get_auto_partitions(structure)
            if auto_parts:
                print(f"  [Auto] {structure}: 发现分区 {list(auto_parts.keys())}")
                for phase_name, (p_min, p_max) in auto_parts.items():
                    mask = (temps_all >= p_min) & (temps_all <= p_max)
                    if np.sum(mask) >= min_points:
                        fit_tasks.append((temps_all[mask], Ds_all[mask], f"{phase_name}({int(p_min)}-{int(p_max)}K)"))
            else:
                print(f"  [Auto] {structure}: 未找到聚类结果，跳过自动分区")

        # 2. 手动分区
        if partitions:
            for p_min, p_max in partitions:
                mask = (temps_all >= p_min) & (temps_all <= p_max)
                if np.sum(mask) >= min_points:
                    fit_tasks.append((temps_all[mask], Ds_all[mask], f'{int(p_min)}-{int(p_max)}K'))

        # 执行拟合
        current_fits = []
        for temps, Ds, range_name in fit_tasks:
            fit = fit_arrhenius(temps, Ds, range_name)
            if fit:
                fit_res = {
                    'structure': structure,
                    'element': element,
                    'temp_range': range_name,
                    **fit
                }
                results.append(fit_res)
                current_fits.append(fit_res)

        # 绘图
        if do_plot and current_fits:
            try:
                fig, ax = plt.subplots(figsize=(8, 6))
                x_all = 1.0 / temps_all
                y_all = np.log(Ds_all)
                ax.scatter(x_all, y_all, color='black', label='Data', zorder=5)

                colors = ['C1', 'C2', 'C3', 'C4']
                for i, fit in enumerate(current_fits):
                    # 整体拟合用虚线，分区拟合用实线
                    is_overall = fit['temp_range'] == 'overall'
                    ls = '--' if is_overall else '-'
                    color = 'gray' if is_overall else colors[i % len(colors)]
                    alpha = 0.5 if is_overall else 1.0
                    
                    # 确定拟合线的范围
                    if is_overall:
                        xi = np.linspace(min(x_all), max(x_all), 100)
                    else:
                        # 解析范围字符串，例如 "partition1(200-700K)" 或 "300-400K"
                        range_str = fit['temp_range']
                        if '(' in range_str:
                            # 处理 "partition1(200-700K)" 格式
                            inner = range_str.split('(')[1].split(')')[0]
                        else:
                            # 处理 "300-400K" 格式
                            inner = range_str
                        
                        t_min, t_max = map(float, inner.replace('K','').split('-'))
                        xi = np.linspace(1.0/t_max, 1.0/t_min, 50)
                        
                    yi = fit['slope'] * xi + fit['intercept']
                    ax.plot(xi, yi, color=color, linestyle=ls, alpha=alpha, linewidth=2,
                            label=f'{fit["temp_range"]}: Ea={fit["Ea_eV"]:.3f} eV')

                ax.set_xlabel('1/T (1/K)', fontsize=12)
                ax.set_ylabel('ln(D (Å²/fs))', fontsize=12)
                ax.legend(fontsize=10)
                ax.grid(True, linestyle=':', alpha=0.6)
                ax.set_title(f'Arrhenius Plot: {structure} ({element})', fontsize=14)

                plot_path = PLOTS_DIR / f'{structure}_{element}_arrhenius.png'
                plt.tight_layout()
                plt.savefig(plot_path, dpi=150)
                plt.close(fig)
            except Exception as e:
                print(f"绘图失败: {structure} {element}: {e}")

    res_df = pd.DataFrame(results)
    # 保存CSV
    res_df.to_csv(output_dir / 'arrhenius_per_structure_element.csv', index=False)

    # 另外保存按 structure 的 overall (element=='all')
    overall = res_df[res_df['element']=='all'].copy()
    overall.to_csv(output_dir / 'arrhenius_per_structure.csv', index=False)

    print(f"保存拟合结果: {output_dir / 'arrhenius_per_structure_element.csv'}")
    
    # 如果指定了元素筛选，返回筛选后的结果
    if elements_filter:
        elements_lower = [e.lower() for e in elements_filter]
        res_df = res_df[res_df['element'].str.lower().isin(elements_lower)]

    return res_df


def plot_ea_trends(res_df, output_dir=OUTPUT_DIR, hide_low_t=False):
    """分析 Ea 随 Sn 含量的变化趋势并绘图 (高定制化版本)"""
    if res_df.empty:
        return

    import matplotlib.ticker as ticker
    # 设置全局字体
    try:
        plt.rcParams['font.family'] = 'Arial'
    except:
        pass

    # 提取 Sn 含量
    res_df = res_df.copy()
    res_df['sn_count'] = res_df['structure'].apply(lambda x: parse_structure_name(x)[1] if parse_structure_name(x) else 0)
    
    # 只看 element='all' 的结果
    df_all = res_df[res_df['element'] == 'all'].copy()
    if df_all.empty:
        return

    # 识别分区类型
    def categorize_range(r):
        if 'partition1' in r: return 'Low-T Phase'
        if 'partition2' in r: return 'High-T Phase'
        return None
    
    df_all['range_cat'] = df_all['temp_range'].apply(categorize_range)
    df_all = df_all.dropna(subset=['range_cat'])
    
    # 如果隐藏低温区
    if hide_low_t:
        df_all = df_all[df_all['range_cat'] != 'Low-T Phase']
    
    # 绘图设置: 10x8 英寸
    fig, ax = plt.subplots(figsize=(10, 8))
    
    markers = {'Low-T Phase': 'o', 'High-T Phase': 's'}
    colors = {'Low-T Phase': '#1f77b4', 'High-T Phase': '#d62728'}
    
    # 数据点大小 50 (markersize = sqrt(50))
    msize = np.sqrt(50)

    categories = ['High-T Phase'] if hide_low_t else ['Low-T Phase', 'High-T Phase']
    for cat in categories:
        subset = df_all[df_all['range_cat'] == cat].sort_values('sn_count')
        if not subset.empty:
            ax.plot(subset['sn_count'], subset['Ea_eV'], marker=markers[cat], 
                     color=colors[cat], label=cat, linewidth=3, 
                     markersize=msize, markeredgecolor='black', markeredgewidth=1)

    # 坐标轴标签: 34号字，不加粗，使用 LaTeX 下标
    ax.set_xlabel('Sn Count', fontsize=34, fontweight='normal', labelpad=15)
    ax.set_ylabel('$E_a$ (eV)', fontsize=34, fontweight='normal', labelpad=15)
    
    # 坐标轴数字: 28号字，刻度朝外
    ax.tick_params(axis='both', which='major', labelsize=28, direction='out', length=10, width=2)
    
    # 刻度设置: 4-7个，对称，整数化
    ax.set_xticks([0, 2, 4, 6, 8, 10])
    
    # 根据数据动态调整 Y 轴刻度
    if not df_all.empty:
        y_max = df_all['Ea_eV'].max()
        if y_max < 0.2:
            ax.set_yticks([0.0, 0.05, 0.1, 0.15, 0.2])
            ax.set_ylim(-0.01, 0.21)
        else:
            ax.set_yticks([0.0, 0.3, 0.6, 0.9, 1.2])
            ax.set_ylim(-0.05, 1.3)

    ax.set_xlim(-0.5, 10.5)

    # 辅助线: 不要
    ax.grid(False)
    
    # 四个框: 都要，加粗
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(2)

    # 图注: 28号字，不要边框
    ax.legend(fontsize=28, frameon=False, loc='upper right')
    
    # 透明背景保存
    suffix = '_high_t_only' if hide_low_t else ''
    trend_plot_path = output_dir / f'ea_vs_sn_trends{suffix}.png'
    plt.tight_layout()
    plt.savefig(trend_plot_path, dpi=300, transparent=True)
    print(f"\n定制化趋势分析图已保存至: {trend_plot_path}")
    plt.close()


def parse_args():
    p = argparse.ArgumentParser(description='从每原子扩散系数中提取迁移能 (Arrhenius)')
    p.add_argument('--msd-file', type=str, default=None, help='指定 per-atom diffusion CSV 文件路径')
    p.add_argument('--min-points', type=int, default=2, help='拟合所需最少温度点数 (默认2)')
    p.add_argument('--no-plot', action='store_true', help='禁用绘图')
    
    # 新增筛选参数
    p.add_argument('--only-series', type=str, default=None,
                   help='只分析指定系列（逗号分隔）: pt8snx, pt6snx, sum8')
    p.add_argument('--exclude', '-e', type=str, default=None,
                   help='排除组分，格式 "(pt,sn)" 或 "(pt1,sn1);(pt2,sn2)"')
    p.add_argument('--exclude-structures', type=str, default=None,
                   help='排除特定结构名称（逗号分隔），如 "pt8sn5-1-best"')
    p.add_argument('--elements', type=str, default=None,
                   help='只输出指定元素的结果（逗号分隔），如 "Pt,Sn"')
    p.add_argument('--partitions', '-p', type=str, default=None,
                   help='手动指定拟合温度区间，格式: T1_min-T1_max,T2_min-T2_max，'
                        '例如: 300-500,800-1100')
    p.add_argument('--auto-partition', action='store_true',
                   help='根据 step6_1 的聚类结果自动进行分区拟合')
    p.add_argument('--hide-low-t-ea', action='store_true',
                   help='在趋势图中隐藏低温区的 Ea 数据')
    
    return p.parse_args()


def main():
    args = parse_args()

    msd_file = Path(args.msd_file) if args.msd_file else find_msd_file()
    if msd_file is None or not msd_file.exists():
        print('未找到 per-atom diffusion CSV 文件，请检查 data 路径')
        return

    msd_df = load_per_atom_D(msd_file)
    
    # 解析筛选参数
    only_series = args.only_series.split(',') if args.only_series else None
    exclude_compositions = parse_exclude_arg(args.exclude)
    exclude_structures = args.exclude_structures.split(',') if args.exclude_structures else None
    
    # 解析分区参数
    partitions = []
    if args.partitions:
        try:
            for part in args.partitions.split(','):
                t_min, t_max = map(float, part.strip().split('-'))
                partitions.append((t_min, t_max))
        except Exception as e:
            print(f"解析分区参数失败: {e}")
    
    # 筛选数据
    if only_series or exclude_compositions or exclude_structures:
        msd_df = filter_data(msd_df, only_series=only_series, 
                             exclude_compositions=exclude_compositions,
                             exclude_structures=exclude_structures)
    
    # 指定元素筛选
    elements_filter = None
    if args.elements:
        elements_filter = [e.strip() for e in args.elements.split(',')]

    res = run_analysis(msd_df, min_points=args.min_points, do_plot=(not args.no_plot),
                       elements_filter=elements_filter, partitions=partitions,
                       auto_partition=args.auto_partition)
    
    # 趋势分析绘图
    if args.only_series and 'pt8snx' in args.only_series:
        plot_ea_trends(res, hide_low_t=args.hide_low_t_ea)
    
    # 打印结果摘要
    # ...existing code...


if __name__ == '__main__':
    main()
