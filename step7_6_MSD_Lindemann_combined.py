#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
Step 7.6: MSD 曲线 + Lindemann Index 合并图 (双子图，Publication风格)
================================================================================

功能概述
========
将两个独立图合并为一张竖向双子图：

  上子图 (A): Ensemble MSD 曲线 vs 时间 (ps)
    - 逻辑来自 step7_5_plot_ensemble_MSD_by_element.py
    - 900K, Pt / Sn 分元素，蓝橙配色
    - 可选误差带 (--errorbar)

  下子图 (B): Lindemann Index δ vs 温度 (K)
    - 逻辑来自 plot_lindemann_filtered.py
    - 读取 per-atom CSV 进行三级平均
    - 蓝橙配色（与上子图一致）

两个子图 x 轴完全独立（时间 vs 温度），共享颜色图例。

使用示例
========
  # 基本用法 (默认 pt8sn6 @ 900K, 含误差带)
  python step7_6_MSD_Lindemann_combined.py

  # 指定结构
  python step7_6_MSD_Lindemann_combined.py --structure pt8sn6

  # 关闭误差带
  python step7_6_MSD_Lindemann_combined.py --no-errorbar

  # 指定 per-atom 文件
  python step7_6_MSD_Lindemann_combined.py \\
    --per-atom-file data/lindemann/per-atoms/sup86-50k/per_atom_master_run_20260311_173823.csv

  # 调整图片尺寸
  python step7_6_MSD_Lindemann_combined.py --figsize 10x14

  # 指定拟合范围
  python step7_6_MSD_Lindemann_combined.py --fit-range 20-140

输出
====
  results/combined_MSD_Lindemann/
    └── {structure}_MSD_Lindemann_combined.png
"""

import re
import argparse
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

warnings.filterwarnings('ignore')

# ============================================================================
# 全局配置
# ============================================================================

BASE_DIR = Path(__file__).parent

# ── 数据路径 ──
DATA_PATHS = {
    'standard': BASE_DIR / 'data' / 'gmx_msd' / 'unwrap' / 'gmx_msd_results_20251118_152614',
    'air':      BASE_DIR / 'data' / 'gmx_msd' / 'unwrap' / 'air' / 'gmx_msd_results_20251124_170114',
}

DEFAULT_PER_ATOM_FILE = (BASE_DIR / 'data' / 'lindemann' / 'per-atoms' /
                         'sup86-50k' / 'per_atom_master_run_20260311_173823.csv')

STRUCTURE_LEVEL_FILE = (BASE_DIR / 'results' /
                        'step7_8_2_alloy_series_sup86_50k' /
                        'structure_level_lindemann.csv')

OUTPUT_DIR = BASE_DIR / 'results' / 'combined_MSD_Lindemann'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ── 配色（与 step7_5 / plot_lindemann_filtered 保持一致）──
ELEMENT_COLORS = {
    'Pt':   '#1f77b4',  # 蓝色
    'Sn':   '#ff7f0e',  # 橙色
    'PtSn': '#2ca02c',  # 绿色（MSD 整体曲线备用）
}

ELEMENT_MARKERS = {
    'Pt': 'o',
    'Sn': 's',
}


# ============================================================================
# ── 部分 1: MSD 数据加载（来自 step7_5）──
# ============================================================================

def read_gmx_msd_xvg(filepath):
    """读取 GMX MSD .xvg 文件，返回 (time_ps, msd_A2)。"""
    time_data, msd_data = [], []
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#') or line.startswith('@'):
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        t = float(parts[0])
                        msd_a2 = float(parts[1]) * 100  # nm² → Å²
                        time_data.append(t)
                        msd_data.append(msd_a2)
                    except ValueError:
                        continue
    except Exception:
        return None, None
    if not time_data:
        return None, None
    return np.array(time_data), np.array(msd_data)


def build_file_index(data_path, structure, temperatures=('900K',)):
    """构建 {(temp, element): [xvg_path, ...]} 索引。"""
    file_index = defaultdict(list)
    for xvg_file in data_path.rglob('*_msd_*.xvg'):
        try:
            parts = xvg_file.parts
            filename = xvg_file.stem
            if '_msd_' not in filename:
                continue
            element = filename.split('_msd_')[-1]
            if element not in ('Pt', 'Sn', 'PtSn'):
                continue
            temperature = composition = None
            for i in range(len(parts) - 1, 0, -1):
                if parts[i].endswith('K'):
                    temperature = parts[i]
                    composition = parts[i - 1]
                    break
            if not temperature or not composition:
                continue
            if not composition.lower().startswith(structure.lower()):
                continue
            if temperature not in temperatures:
                continue
            file_index[(temperature, element)].append(xvg_file)
        except Exception:
            continue
    return file_index


def load_msd_data(file_index, temperatures=('900K',)):
    """加载所有 MSD 时间序列，返回 msd_cache。"""
    msd_cache = {}
    for temp in temperatures:
        for element in ('Pt', 'Sn', 'PtSn'):
            key = (temp, element)
            files = file_index.get(key, [])
            msd_list = []
            for fp in files:
                t, m = read_gmx_msd_xvg(fp)
                if t is not None:
                    msd_list.append((t, m))
            if msd_list:
                msd_cache[key] = msd_list
    return msd_cache


def _calc_D_single(time, msd, fit_range):
    """单窗口线性拟合，返回 (D_cm2_s, r2)。"""
    from scipy import stats
    mask = (time >= fit_range[0]) & (time <= fit_range[1])
    if mask.sum() < 10:
        n = len(time)
        mask = np.zeros(n, dtype=bool)
        mask[int(n * 0.2):int(n * 0.8)] = True
    slope, _, r, _, _ = stats.linregress(time[mask], msd[mask])
    return slope / 6.0 * 1e-4, r ** 2


# ============================================================================
# ── 部分 2: Lindemann 数据加载（来自 plot_lindemann_filtered）──
# ============================================================================

def load_per_atom_as_structure_level(per_atom_file, structure_filter=None):
    """
    读取 per-atom CSV，三级聚合 (原子 × run → 元素 → 温度)，
    返回与 structure_level_lindemann.csv 兼容的宽格式 DataFrame。
    """
    per_atom_path = Path(per_atom_file)
    if not per_atom_path.is_absolute():
        per_atom_path = BASE_DIR / per_atom_path
    if not per_atom_path.exists():
        print(f"  [ERROR] per-atom 文件不存在: {per_atom_path}")
        return None

    print(f"  [INFO] 读取 per-atom 文件: {per_atom_path.name}")
    df = pd.read_csv(per_atom_path)
    df = df.rename(columns={'温度(K)': 'temperature', '目录': 'path',
                             '结构': 'structure_col', 'lindemann_index': 'delta'})

    def _struct_from_path(p):
        m = re.search(r'/(pt\d+sn\d+[^/]*)/T\d+', str(p), re.IGNORECASE)
        return m.group(1).lower() if m else 'unknown'

    df['structure'] = df['path'].apply(_struct_from_path)

    if structure_filter:
        df = df[df['structure'].str.contains(structure_filter, case=False, na=False)]
        if len(df) == 0:
            print(f"  [ERROR] 过滤后无数据: structure_filter='{structure_filter}'")
            return None

    print(f"  [INFO] {len(df)} 原子行，结构: {sorted(df['structure'].unique())}")

    agg = (df.groupby(['structure', 'temperature', 'element'])['delta']
             .agg(mean='mean', std='std').reset_index())

    mean_wide = agg.pivot_table(index=['structure', 'temperature'],
                                columns='element', values='mean').reset_index()
    std_wide  = agg.pivot_table(index=['structure', 'temperature'],
                                columns='element', values='std').reset_index()
    mean_wide.columns.name = std_wide.columns.name = None

    rename_m = {e: f'delta_{e}' for e in mean_wide.columns
                if e not in ('structure', 'temperature')}
    rename_s = {e: f'delta_{e}_std' for e in std_wide.columns
                if e not in ('structure', 'temperature')}
    mean_wide = mean_wide.rename(columns=rename_m)
    std_wide  = std_wide.rename(columns=rename_s)
    result = mean_wide.merge(std_wide, on=['structure', 'temperature'], how='left')
    result['anomaly_flag'] = ''
    result = result.sort_values(['structure', 'temperature']).reset_index(drop=True)
    return result


def load_structure_level_lindemann(filepath=None):
    """读取 structure_level_lindemann.csv（step7_8_2 输出）。"""
    path = Path(filepath) if filepath else STRUCTURE_LEVEL_FILE
    if not path.is_absolute():
        path = BASE_DIR / path
    if not path.exists():
        print(f"  [ERROR] 文件不存在: {path}")
        return None
    df = pd.read_csv(path, encoding='utf-8-sig')
    print(f"  [INFO] 加载 structure_level 数据: {len(df)} 条")
    return df


def extract_lindemann_stats(df, structure):
    """
    从宽格式 DataFrame 提取某结构的长格式统计 (temp, element, delta_mean, delta_std)。
    若存在多个变体（如 pt8sn6-1-best + pt8sn6-1-best-2），按温度取 ensemble 平均。
    """
    df_s = df[df['structure'].str.contains(structure, case=False)].copy()
    if len(df_s) == 0:
        print(f"  [ERROR] 未找到结构 {structure}")
        return None

    # 多变体 → ensemble 平均
    if df_s['structure'].nunique() > 1:
        df_s = df_s.groupby('temperature').agg(
            delta_Pt=('delta_Pt', 'mean'),
            delta_Sn=('delta_Sn', 'mean'),
            delta_Pt_std=('delta_Pt', 'std'),
            delta_Sn_std=('delta_Sn', 'std'),
        ).reset_index()

    rows = []
    for _, row in df_s.iterrows():
        temp = row['temperature']
        for elem in ('Pt', 'Sn'):
            rows.append({
                'temp':       temp,
                'element':    elem,
                'delta_mean': row[f'delta_{elem}'],
                'delta_std':  row[f'delta_{elem}_std'],
            })
    df_stats = pd.DataFrame(rows).sort_values(['element', 'temp'])
    return df_stats


# ============================================================================
# ── 部分 3.5: CSV 数据导出 ──
# ============================================================================

def _export_csvs(structure, temperature, msd_export, df_lind_stats, output_dir):
    """
    导出两份 CSV 供数据审查：
      1. {structure}_{temperature}_MSD_data.csv   —— MSD 时间序列（均值 ± std）+ D
      2. {structure}_Lindemann_summary.csv         —— Lindemann δ vs T（宽格式）
    """
    # ── 1. MSD 时间序列 CSV ──────────────────────────────────────────────────
    if msd_export:
        # 以 Pt 的时间轴为基准
        ref_elem = next(iter(msd_export))
        t_ref = msd_export[ref_elem][0]
        df_msd = pd.DataFrame({'time_ps': t_ref})

        D_summary_rows = []
        for elem, (t, mean, std, D, r2, n_runs) in msd_export.items():
            # 对齐到 t_ref 长度（若不同取最短）
            n = min(len(t_ref), len(mean))
            df_msd[f'MSD_{elem}_mean_A2'] = np.nan
            df_msd[f'MSD_{elem}_std_A2']  = np.nan
            df_msd.loc[:n-1, f'MSD_{elem}_mean_A2'] = mean[:n]
            df_msd.loc[:n-1, f'MSD_{elem}_std_A2']  = std[:n]
            D_summary_rows.append({
                'element':    elem,
                'temperature': temperature,
                'D_cm2_s':    D,
                'D_1e5_cm2_s': round(D * 1e5, 6),
                'r2':          round(r2, 6),
                'n_runs':      n_runs,
            })

        msd_csv = output_dir / f'{structure}_{temperature}_MSD_data.csv'
        df_msd.to_csv(msd_csv, index=False, float_format='%.6f')
        print(f"[SAVED CSV] {msd_csv}  ({len(df_msd)} 行, {df_msd.shape[1]} 列)")

        # D 值汇总打印
        print('\n── MSD 扩散系数汇总 ──────────────────────────────────')
        df_D = pd.DataFrame(D_summary_rows)
        print(df_D.to_string(index=False))
        D_csv = output_dir / f'{structure}_{temperature}_D_values.csv'
        df_D.to_csv(D_csv, index=False, float_format='%.6f')
        print(f"[SAVED CSV] {D_csv}")

    # ── 2. Lindemann 汇总 CSV ────────────────────────────────────────────────
    if df_lind_stats is not None and len(df_lind_stats) > 0:
        df_pivot = df_lind_stats.pivot(index='temp', columns='element',
                                       values=['delta_mean', 'delta_std'])
        df_pivot.columns = [f'delta_{e}' if stat == 'delta_mean'
                            else f'delta_{e}_std'
                            for stat, e in df_pivot.columns]
        df_pivot = df_pivot.reset_index().rename(columns={'temp': 'temperature_K'})
        # 比值
        if 'delta_Pt' in df_pivot.columns and 'delta_Sn' in df_pivot.columns:
            df_pivot['ratio_Pt_Sn'] = (df_pivot['delta_Pt'] /
                                        df_pivot['delta_Sn'].replace(0, float('nan'))).round(4)
        col_order = ['temperature_K', 'delta_Pt', 'delta_Pt_std',
                     'delta_Sn', 'delta_Sn_std', 'ratio_Pt_Sn']
        df_pivot = df_pivot[[c for c in col_order if c in df_pivot.columns]]
        df_pivot = df_pivot.sort_values('temperature_K').reset_index(drop=True)

        lind_csv = output_dir / f'{structure}_Lindemann_summary.csv'
        df_pivot.to_csv(lind_csv, index=False, float_format='%.6f')
        print(f"\n── Lindemann Index 汇总 ───────────────────────────────")
        print(df_pivot.to_string(index=False))
        print(f"[SAVED CSV] {lind_csv}")


# ============================================================================
# ── 部分 4: 绘图 ──
# ============================================================================

def _apply_pub_style(ax, xlabel, ylabel,
                     tick_size=24, label_size=30, spine_lw=1.5,
                     tick_direction='out', tick_length=6):
    """统一应用 Publication 样式到单个 Axes。"""
    ax.set_xlabel(xlabel, fontsize=label_size)
    ax.set_ylabel(ylabel, fontsize=label_size)
    ax.tick_params(axis='both', which='major',
                   labelsize=tick_size, direction=tick_direction,
                   length=tick_length, width=spine_lw)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(spine_lw)
    ax.grid(False)


def plot_combined(msd_cache, df_lind_stats, structure,
                  temperature='900K',
                  fit_range=(20, 140),
                  errorbar=True,
                  figsize=(10, 14),
                  output_dir=OUTPUT_DIR,
                  panel_labels=('(a)', '(b)'),
                  ylabel_b=None):
    """
    绘制双子图（上: MSD 曲线；下: Lindemann δ vs T）并保存。
    """
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.sans-serif'] = ['Arial']
    plt.rcParams['mathtext.default'] = 'regular'

    # 是否需要为面板标签留出顶部空间
    need_panel_label_space = bool(panel_labels[0] or panel_labels[1])

    fig = plt.figure(figsize=figsize)
    gs  = gridspec.GridSpec(2, 1, figure=fig)
    ax_msd  = fig.add_subplot(gs[0])
    ax_lind = fig.add_subplot(gs[1])

    # ── 上子图: MSD ──────────────────────────────────────────────────────────
    printed_D = {}
    msd_export = {}   # 用于 CSV 导出: {element: (t_common, msd_mean, msd_std, D, r2)}
    for element in ('Pt', 'Sn'):
        key = (temperature, element)
        msd_list = msd_cache.get(key)
        if not msd_list:
            print(f"  [WARNING] 无 {temperature} {element} MSD 数据")
            continue

        min_len = min(len(m) for _, m in msd_list)
        msd_arr = np.array([m[:min_len] for _, m in msd_list])
        t_common = msd_list[0][0][:min_len]

        msd_mean = msd_arr.mean(axis=0)
        msd_std  = (msd_arr.std(axis=0, ddof=1)
                    if len(msd_arr) > 1 else np.zeros_like(msd_mean))

        color = ELEMENT_COLORS[element]

        if errorbar:
            ax_msd.fill_between(t_common,
                                msd_mean - msd_std,
                                msd_mean + msd_std,
                                color=color, alpha=0.25, zorder=1)

        ax_msd.plot(t_common, msd_mean,
                    color=color, linewidth=3.5,
                    label=element, zorder=3)

        # 计算并打印 D
        D, r2 = _calc_D_single(t_common, msd_mean, fit_range)
        printed_D[element] = D
        msd_export[element] = (t_common, msd_mean, msd_std, D, r2, len(msd_list))
        print(f"  MSD {temperature} {element}: D = {D*1e5:.4f} x10^-5 cm2/s  (R2 = {r2:.4f})")

    # x 轴刻度
    t_max = max((msd_cache[(temperature, e)][0][0][-1]
                 for e in ('Pt', 'Sn')
                 if (temperature, e) in msd_cache), default=175)
    ax_msd.set_xlim(0, t_max)
    ax_msd.set_ylim(bottom=0)
    if t_max <= 200:
        ax_msd.set_xticks(np.arange(0, t_max + 1, 50))
    elif t_max <= 500:
        ax_msd.set_xticks(np.arange(0, t_max + 1, 100))
    else:
        ax_msd.set_xticks(np.arange(0, t_max + 1, 200))

    _apply_pub_style(ax_msd,
                     xlabel='Time (ps)',
                     ylabel=r'MSD ($\AA^2$)')
    ax_msd.legend(fontsize=24, loc='upper left', frameon=False)

    # 面板标签 —— clip_on=False 防止被裁剪
    if panel_labels[0]:
        ax_msd.text(-0.12, 1.06, panel_labels[0],
                    transform=ax_msd.transAxes,
                    fontsize=28, fontweight='bold', va='bottom', ha='left',
                    clip_on=False)

    # ── 下子图: Lindemann ──────────────────────────────────────────────────────
    if df_lind_stats is not None and len(df_lind_stats) > 0:
        for element in ('Pt', 'Sn'):
            edata = df_lind_stats[df_lind_stats['element'] == element].sort_values('temp')
            if len(edata) == 0:
                continue

            temps      = edata['temp'].values
            delta_mean = edata['delta_mean'].values
            delta_std  = edata['delta_std'].values

            color  = ELEMENT_COLORS[element]
            marker = ELEMENT_MARKERS[element]

            ax_lind.fill_between(temps,
                                 delta_mean - delta_std,
                                 delta_mean + delta_std,
                                 color=color, alpha=0.2, zorder=1)
            ax_lind.plot(temps, delta_mean,
                         marker=marker, color=color,
                         linewidth=3.5, markersize=9,
                         markeredgewidth=1.5, markeredgecolor='white',
                         label=element, zorder=3)

        ax_lind.set_xlim(
            df_lind_stats['temp'].min() - 20,
            df_lind_stats['temp'].max() + 20
        )
        ax_lind.set_ylim(bottom=0)
        ax_lind.locator_params(axis='x', nbins=6)
        ax_lind.locator_params(axis='y', nbins=6)
    else:
        ax_lind.text(0.5, 0.5, 'No Lindemann data',
                     ha='center', va='center',
                     transform=ax_lind.transAxes, fontsize=20)

    _apply_pub_style(ax_lind,
                     xlabel='Temperature (K)',
                     ylabel=ylabel_b if ylabel_b is not None else 'Lindemann Index δ')
    ax_lind.legend(fontsize=24, loc='upper left', frameon=False)

    if panel_labels[1]:
        ax_lind.text(-0.12, 1.06, panel_labels[1],
                     transform=ax_lind.transAxes,
                     fontsize=28, fontweight='bold', va='bottom', ha='left',
                     clip_on=False)

    # ── CSV 导出 ──────────────────────────────────────────────────────────────
    _export_csvs(structure, temperature, msd_export, df_lind_stats, output_dir)

    # ── 保存图片 ──────────────────────────────────────────────────────────────
    # tight_layout 自动调整内部间距，top 留出面板标签空间
    top_pad = 0.90 if need_panel_label_space else 0.96
    plt.tight_layout(rect=[0, 0, 1, top_pad])
    # bbox_inches=None → 严格按 figsize 输出，不再二次裁剪
    out_path = output_dir / f'{structure}_MSD_Lindemann_combined.png'
    plt.savefig(out_path, dpi=300, bbox_inches=None,
                transparent=True, facecolor='none')
    print(f"\n[SAVED] {out_path}")
    plt.show()
    plt.close()
    return out_path


# ============================================================================
# ── 主函数 ──
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Step 7.6: MSD + Lindemann 双子图 (Publication风格)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument('--structure', type=str, default='pt8sn6',
                        help='结构名称 (默认: pt8sn6)')
    parser.add_argument('--temp', type=str, default='900K',
                        help='MSD 子图使用的温度 (默认: 900K)')
    parser.add_argument('--data-type', type=str, default='standard',
                        choices=['standard', 'air'],
                        help='MSD 数据类型 (默认: standard)')
    parser.add_argument('--fit-range', type=str, default='20-140',
                        help='MSD 拟合范围 ps, 格式: start-end (默认: 20-140)')
    parser.add_argument('--no-errorbar', action='store_true',
                        help='关闭 MSD 误差带')
    parser.add_argument('--per-atom-file', type=str, default=None,
                        help=('Lindemann per-atom CSV 路径。'
                              '默认使用 data/lindemann/per-atoms/sup86-50k/...'))
    parser.add_argument('--lindemann-file', type=str, default=None,
                        help=('直接指定 structure_level_lindemann.csv。'
                              '若同时提供 --per-atom-file，则 per-atom 优先。'))
    parser.add_argument('--figsize', type=str, default='10x14',
                        help='图片尺寸 宽x高 (默认: 10x14)')
    parser.add_argument('--no-panel-labels', action='store_true',
                        help='不显示 (a)(b) 面板标签')
    parser.add_argument('--ylabel-b', type=str, default=None,
                        nargs='?', const='',
                        help=('下子图 (Lindemann) 的 Y 轴标签。'
                              '默认: "Lindemann Index delta"。'
                              '使用 --ylabel-b 不加值可隐藏标签，'
                              '或 --ylabel-b "delta" 自定义'))
    parser.add_argument('--output-dir', type=str, default=None,
                        help='输出目录 (默认: results/combined_MSD_Lindemann/)')

    args = parser.parse_args()

    # ── 解析参数 ──
    try:
        r_start, r_end = [float(x) for x in args.fit_range.split('-')]
        fit_range = (r_start, r_end)
    except Exception:
        fit_range = (20, 140)
        print(f"  [WARNING] --fit-range 格式错误，使用默认 20-140")

    try:
        fw, fh = [float(x) for x in args.figsize.lower().split('x')]
        figsize = (fw, fh)
    except Exception:
        figsize = (10, 14)

    temperature = args.temp if args.temp.endswith('K') else f'{args.temp}K'

    output_dir = Path(args.output_dir) if args.output_dir else OUTPUT_DIR
    if not output_dir.is_absolute():
        output_dir = BASE_DIR / output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    panel_labels = ('', '') if args.no_panel_labels else ('(a)', '(b)')

    print('=' * 72)
    print('Step 7.6: MSD + Lindemann 双子图')
    print(f'  结构: {args.structure}')
    print(f'  MSD 温度: {temperature}')
    print(f'  拟合范围: {fit_range[0]:.0f}-{fit_range[1]:.0f} ps')
    print(f'  误差带: {"关" if args.no_errorbar else "开"}')
    print(f'  图片尺寸: {figsize[0]:.0f}x{figsize[1]:.0f} 英寸')
    print('=' * 72)

    # ── 加载 MSD 数据 ──
    data_path = DATA_PATHS[args.data_type]
    if not data_path.exists():
        print(f"[ERROR] MSD 数据路径不存在: {data_path}")
        return

    print('\n[1/2] 加载 MSD 数据...')
    file_index = build_file_index(data_path, args.structure, (temperature,))
    if not file_index:
        print(f"  [ERROR] 未找到 {args.structure} @ {temperature} 的 xvg 文件")
        print(f"         数据路径: {data_path}")
        msd_cache = {}
    else:
        msd_cache = load_msd_data(file_index, (temperature,))
        total = sum(len(v) for v in msd_cache.values())
        print(f"  [OK] 加载 {total} 条 MSD 曲线")

    # ── 加载 Lindemann 数据 ──
    print('\n[2/2] 加载 Lindemann 数据...')
    df_lind_stats = None

    per_atom_file = args.per_atom_file or (
        str(DEFAULT_PER_ATOM_FILE) if DEFAULT_PER_ATOM_FILE.exists() else None
    )

    if per_atom_file:
        df_wide = load_per_atom_as_structure_level(per_atom_file,
                                                   structure_filter=args.structure)
        if df_wide is not None:
            df_lind_stats = extract_lindemann_stats(df_wide, args.structure)
    elif args.lindemann_file or STRUCTURE_LEVEL_FILE.exists():
        lpath = args.lindemann_file or str(STRUCTURE_LEVEL_FILE)
        df_wide = load_structure_level_lindemann(lpath)
        if df_wide is not None:
            df_lind_stats = extract_lindemann_stats(df_wide, args.structure)
    else:
        print('  [WARNING] 未找到 Lindemann 数据，下子图将留空')

    if df_lind_stats is not None:
        print(f"  [OK] Lindemann 数据: {df_lind_stats['temp'].nunique()} 个温度点")

    # ── 绘图 ──
    print('\n绘图...')
    plot_combined(
        msd_cache      = msd_cache,
        df_lind_stats  = df_lind_stats,
        structure      = args.structure,
        temperature    = temperature,
        fit_range      = fit_range,
        errorbar       = not args.no_errorbar,
        figsize        = figsize,
        output_dir     = output_dir,
        panel_labels   = panel_labels,
        ylabel_b       = args.ylabel_b,
    )

    print('\n' + '=' * 72)
    print(f'完成！输出目录: {output_dir}')
    print('=' * 72)


if __name__ == '__main__':
    main()
