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


def filter_data(df, only_series=None, exclude_compositions=None):
    """
    根据参数筛选数据
    
    Args:
        df: 原始数据 (必须有 'structure' 列)
        only_series: 只包含的系列列表，如 ['sum8', 'pt8snx']
        exclude_compositions: 排除的组分列表，如 [(8, 0), (6, 0)]
    
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


def fit_arrhenius(temp_list, D_list):
    """对给定温度列表和对应平均 D 进行 Arrhenius 拟合。
    返回: dict with keys: n_points, slope, intercept, Ea_eV, D0, r2
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
        'n_points': int(n),
        'slope': float(slope),
        'intercept': float(intercept),
        'Ea_eV': float(Ea_eV),
        'D0': float(D0),
        'r2': float(r2)
    }


def run_analysis(msd_df, min_points=3, do_plot=True, output_dir=OUTPUT_DIR, elements_filter=None):
    output_dir.mkdir(parents=True, exist_ok=True)
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    mean_df = summarize_mean_D(msd_df)

    results = []

    # 遍历每个 structure 与 element
    for (structure, element), group in mean_df.groupby(['structure','element']):
        temps = group['temp'].values
        Ds = group['D_mean'].values
        if len(temps) < min_points:
            # 记录不足样本
            results.append({
                'structure': structure,
                'element': element,
                'n_points': int(len(temps)),
                'slope': np.nan,
                'intercept': np.nan,
                'Ea_eV': np.nan,
                'D0': np.nan,
                'r2': np.nan
            })
            continue

        fit = fit_arrhenius(temps, Ds)
        if fit is None:
            results.append({
                'structure': structure,
                'element': element,
                'n_points': int(len(temps)),
                'slope': np.nan,
                'intercept': np.nan,
                'Ea_eV': np.nan,
                'D0': np.nan,
                'r2': np.nan
            })
            continue

        results.append({
            'structure': structure,
            'element': element,
            'n_points': fit['n_points'],
            'slope': fit['slope'],
            'intercept': fit['intercept'],
            'Ea_eV': fit['Ea_eV'],
            'D0': fit['D0'],
            'r2': fit['r2']
        })

        # 绘图
        if do_plot:
            try:
                fig, ax = plt.subplots(figsize=(6,5))
                x = 1.0 / temps
                y = np.log(Ds)
                ax.scatter(x, y, label=f'{element} data')

                # 画拟合线
                xi = np.linspace(min(x), max(x), 100)
                yi = fit['slope'] * xi + fit['intercept']
                ax.plot(xi, yi, color='C1', linestyle='--', label=f'fit: Ea={fit["Ea_eV"]:.3f} eV')

                ax.set_xlabel('1/T (1/K)')
                ax.set_ylabel('ln(D (Å²/fs))')
                ax.legend()
                ax.grid(True)
                ax.set_title(f'{structure} - {element}')

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
    print(f"保存总体结果: {output_dir / 'arrhenius_per_structure.csv'}")
    
    # 如果指定了元素筛选，返回筛选后的结果
    if elements_filter:
        elements_lower = [e.lower() for e in elements_filter]
        res_df = res_df[res_df['element'].str.lower().isin(elements_lower)]

    return res_df


def parse_args():
    p = argparse.ArgumentParser(description='从每原子扩散系数中提取迁移能 (Arrhenius)')
    p.add_argument('--msd-file', type=str, default=None, help='指定 per-atom diffusion CSV 文件路径')
    p.add_argument('--min-points', type=int, default=3, help='拟合所需最少温度点数 (默认3)')
    p.add_argument('--no-plot', action='store_true', help='禁用绘图')
    
    # 新增筛选参数
    p.add_argument('--only-series', type=str, default=None,
                   help='只分析指定系列（逗号分隔）: pt8snx, pt6snx, sum8')
    p.add_argument('--exclude', '-e', type=str, default=None,
                   help='排除组分，格式 "(pt,sn)" 或 "(pt1,sn1);(pt2,sn2)"')
    p.add_argument('--elements', type=str, default=None,
                   help='只输出指定元素的结果（逗号分隔），如 "Pt,Sn"')
    
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
    
    # 筛选数据
    if only_series or exclude_compositions:
        msd_df = filter_data(msd_df, only_series=only_series, exclude_compositions=exclude_compositions)
    
    # 指定元素筛选
    elements_filter = None
    if args.elements:
        elements_filter = [e.strip() for e in args.elements.split(',')]

    res = run_analysis(msd_df, min_points=args.min_points, do_plot=(not args.no_plot),
                       elements_filter=elements_filter)
    
    # 打印 Pt 和 Sn 的结果摘要
    print("\n" + "="*70)
    print("迁移能垒 (Ea) 结果摘要")
    print("="*70)
    
    for elem in ['Pt', 'Sn', 'all']:
        elem_df = res[res['element'].str.lower() == elem.lower()]
        if len(elem_df) > 0:
            print(f"\n【{elem}】")
            for _, row in elem_df.iterrows():
                if pd.notna(row['Ea_eV']):
                    print(f"  {row['structure']:20s}  Ea = {row['Ea_eV']:8.4f} eV  "
                          f"D0 = {row['D0']:.2e}  R² = {row['r2']:.4f}  (n={row['n_points']})")
                else:
                    print(f"  {row['structure']:20s}  拟合失败 (n={row['n_points']})")

    print('\nDone')


if __name__ == '__main__':
    main()
