#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
Step 5.1.1: 系列热图与系列图分析
================================================================================

功能概述
========
本脚本绘制指定系列的：
1. Lindemann 指数热图：横轴为组分变量，纵轴温度，颜色表示 δ 值
2. 熔点 vs 组分系列图：展示熔点随组分变化的趋势

支持的系列
==========
--series 参数可选值：
  pt8snx    : Pt8Snx 系列（负载型，O=0）- 默认，x轴=Sn数量
  pt8-xsnx  : Pt8-xSnx 系列（总金属=8，O=0），x轴=Sn数量
  pt6snx    : Pt6Snx 系列（负载型，O=0），x轴=Sn数量
  o1        : O1 系列（含1个O的氧化物），x轴=金属原子数(Pt+Sn)
  o2        : O2 系列（含2个O的氧化物），x轴=金属原子数(Pt+Sn)
  o3        : O3 系列（含3个O的氧化物），x轴=金属原子数(Pt+Sn)
  o4        : O4 系列（含4个O的氧化物），x轴=金属原子数(Pt+Sn)
  oxide     : 所有氧化物（O>0），x轴=总原子数(Pt+Sn+O)
  all       : 所有负载型结构，x轴=总原子数

排序规则
========
1. 首先按 x 轴变量（Sn数量/金属原子数/总原子数）升序排列
2. 当 x 值相同时，按 Pt 比例降序排列（Pt 多的排前面）

数据来源
========
1. Lindemann 原始数据: results/step6_0_multi_system/step6_0_all_systems_data.csv
2. 熔点汇总数据: results/step5_1_melting_point/melting_point_summary.csv
   (需先运行 step5_1_melting_point_analysis.py 生成)

输出文件
========
results/step5_1_1_pt8snx/
├── {series}_delta_heatmap.png   # Lindemann 指数热图
└── {series}_tm_vs_comp.png      # 熔点 vs 组分系列图

命令行参数
==========
--series, -s     : 选择分析的系列（默认 pt8snx）
                   可选: pt8snx, pt6snx, o1, o2, o3, o4, oxide, all
--exclude, -e    : 屏蔽指定结构，不参与绑图（可指定多个）
                   例: --exclude Pt8sn0 Pt8sn1
--decimals, -d   : 热图格子中数值的小数位数，可选 2 或 3（默认 3）
                   例: --decimals 2
--fontscale, -f  : 字体缩放比例（默认 1.0，1.2 表示放大 20%）
                   例: --fontscale 1.2
--no-values      : 不在热图格子中显示数值
--no-title       : 不显示图片标题
--no-fit         : 不显示系列图的线性拟合线
--no-show        : 只保存图片，不弹出交互式窗口

使用示例
========
# Pt8Snx 系列（默认），屏蔽 Pt8sn0
& python step5_1_1_pt8snx_plot.py -e Pt8sn0 --no-title --no-show

# O1 含氧系列
& python step5_1_1_pt8snx_plot.py -s o1 --no-title --no-show

# O2 含氧系列，放大字体 1.2 倍
& python step5_1_1_pt8snx_plot.py -s o2 -f 1.2 --no-title --no-show

# O3 含氧系列
& python step5_1_1_pt8snx_plot.py -s o3 --no-title --no-show

# O4 含氧系列
& python step5_1_1_pt8snx_plot.py -s o4 --no-title --no-show

# 所有氧化物系列
& python step5_1_1_pt8snx_plot.py -s oxide --no-title --no-show

# Pt6Snx 系列
& python step5_1_1_pt8snx_plot.py -s pt6snx --no-title --no-show

# 2位小数，不显示数值
& python step5_1_1_pt8snx_plot.py -s o1 -d 2 --no-values --no-title --no-show

注意事项
========
1. 结构名称大小写不敏感（Pt8sn0 和 PT8SN0 等效）
2. 热图使用 RdYlBu_r 色图，δ=0.1 处有黑色虚线标记熔点阈值
3. 系列图使用黑色实心圆点+实线，灰色虚线表示线性拟合
4. 含氧系列：x 轴相同的多个结构会被平均后显示在同一列
5. 热图等高线标记 δ=0.1 熔点阈值

================================================================================
"""
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy import stats
from scipy.interpolate import interp1d

BASE_DIR = Path(__file__).parent
DATA_FILE = BASE_DIR / 'results' / 'step6_0_multi_system' / 'step6_0_all_systems_data.csv'
MP_SUMMARY = BASE_DIR / 'results' / 'step5_1_melting_point' / 'melting_point_summary.csv'
OUTPUT_DIR = BASE_DIR / 'results' / 'step5_1_1_pt8snx'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# 英文字体设置（适合期刊发表）
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.unicode_minus'] = False


def format_structure_label(name):
    """
    将结构名格式化为 LaTeX 下标形式，统一顺序为 Pt_x Sn_y O_z
    例如：Sn5o4pt2 -> Pt$_2$Sn$_5$O$_4$
    """
    import re
    # 提取各元素的数量
    elements = {}
    pattern = r'([A-Za-z]+)(\d+)'
    for match in re.finditer(pattern, name, re.IGNORECASE):
        elem = match.group(1).lower()
        num = int(match.group(2))
        elements[elem] = num
    
    # 按 Pt -> Sn -> O 顺序输出
    result = ''
    for elem in ['pt', 'sn', 'o']:
        if elem in elements and elements[elem] > 0:
            elem_name = elem.capitalize() if elem != 'pt' else 'Pt'
            result += f'{elem_name}$_{{{elements[elem]}}}$'
    
    return result if result else name


def calculate_tm_from_threshold(df_all, structure, threshold=0.1):
    """
    根据自定义阈值计算熔点
    
    Args:
        df_all: 包含所有结构数据的 DataFrame
        structure: 结构名
        threshold: Lindemann 阈值 (默认 0.1)
    
    Returns:
        tuple: (Tm, Tm_err) 或 (nan, nan)
    """
    df_s = df_all[df_all['structure'] == structure]
    df_avg = df_s.groupby('temp')['delta'].mean().reset_index()
    df_avg = df_avg.sort_values('temp')
    
    temps = df_avg['temp'].values
    deltas = df_avg['delta'].values
    
    if len(temps) < 2:
        return np.nan, np.nan
    
    for i in range(len(temps) - 1):
        if deltas[i] < threshold and deltas[i+1] >= threshold:
            T1, T2 = temps[i], temps[i+1]
            d1, d2 = deltas[i], deltas[i+1]
            if d2 != d1:
                Tm = T1 + (threshold - d1) * (T2 - T1) / (d2 - d1)
            else:
                Tm = (T1 + T2) / 2
            Tm_err = (T2 - T1) / 2
            return Tm, Tm_err
    
    return np.nan, np.nan


# 系列配置
SERIES_CONFIG = {
    'pt8snx': {
        'filter': lambda df: (df['Pt'] == 8) & (df['O'] == 0) & (df['type'] == 'supported'),
        'x_col': 'Sn',
        'x_label': 'Number of Sn atoms',
        'title': r'Pt$_8$Sn$_x$',
        'legend': r'Pt$_8$Sn$_x$',
    },
    'pt8-xsnx': {
        'filter': lambda df: (df['Pt'] + df['Sn'] == 8) & (df['O'] == 0) & (df['type'] == 'supported'),
        'x_col': 'Sn',
        'x_label': 'Number of Sn atoms (x)',
        'title': r'Pt$_{8-x}$Sn$_x$',
        'legend': r'Pt$_{8-x}$Sn$_x$',
    },
    'pt6snx': {
        'filter': lambda df: (df['Pt'] == 6) & (df['O'] == 0) & (df['type'] == 'supported'),
        'x_col': 'Sn',
        'x_label': 'Number of Sn atoms',
        'title': r'Pt$_6$Sn$_x$',
        'legend': r'Pt$_6$Sn$_x$',
    },
    'o1': {
        'filter': lambda df: (df['O'] == 1) & (df['type'] == 'oxide'),
        'x_col': 'metal_total',  # 用于排序
        'x_label': 'Structure',
        'title': r'Pt$_x$Sn$_y$O$_1$ (oxide)',
        'legend': r'Pt$_x$Sn$_y$O$_1$',
        'use_structure_labels': True,  # 使用结构名作为 x 轴
    },
    'o2': {
        'filter': lambda df: (df['O'] == 2) & (df['type'] == 'oxide'),
        'x_col': 'metal_total',
        'x_label': 'Structure',
        'title': r'Pt$_x$Sn$_y$O$_2$ (oxide)',
        'legend': r'Pt$_x$Sn$_y$O$_2$',
        'use_structure_labels': True,
    },
    'o3': {
        'filter': lambda df: (df['O'] == 3) & (df['type'] == 'oxide'),
        'x_col': 'metal_total',
        'x_label': 'Structure',
        'title': r'Pt$_x$Sn$_y$O$_3$ (oxide)',
        'legend': r'Pt$_x$Sn$_y$O$_3$',
        'use_structure_labels': True,
    },
    'o4': {
        'filter': lambda df: (df['O'] == 4) & (df['type'] == 'oxide'),
        'x_col': 'metal_total',
        'x_label': 'Structure',
        'title': r'Pt$_x$Sn$_y$O$_4$ (oxide)',
        'legend': r'Pt$_x$Sn$_y$O$_4$',
        'use_structure_labels': True,
    },
    'oxide': {
        'filter': lambda df: (df['O'] > 0) & (df['type'] == 'oxide'),
        'x_col': 'total_atoms',  # 用于排序
        'x_label': 'Structure',
        'title': 'All oxides',
        'legend': 'Oxide',
        'use_structure_labels': True,
    },
    'all': {
        'filter': lambda df: df['type'] == 'supported',
        'x_col': 'total_atoms',
        'x_label': 'Total atoms',
        'title': 'All supported structures',
        'legend': 'Supported',
    },
}


def parse_args():
    p = argparse.ArgumentParser(description='Plot series Lindemann heatmap and Tm vs composition')
    p.add_argument('--series', '-s', default='pt8snx', 
                   choices=list(SERIES_CONFIG.keys()),
                   help='选择分析的系列 (默认 pt8snx)')
    p.add_argument('--threshold', '-t', type=float, default=0.10,
                   help='Lindemann 指数熔点阈值 (默认 0.10，可设为 0.15 等)')
    p.add_argument('--exclude', '-e', nargs='+', default=[], help='屏蔽指定结构（可多个，大小写不敏感）')
    p.add_argument('--decimals', '-d', type=int, default=3, choices=[2, 3], help='热图数值小数位数 (2 或 3，默认 3)')
    p.add_argument('--fontscale', '-f', type=float, default=1.0, help='字体缩放比例 (默认 1.0)')
    p.add_argument('--no-values', action='store_true', help='不在热图格子中显示数值')
    p.add_argument('--no-title', action='store_true', help='不显示图片标题')
    p.add_argument('--no-fit', action='store_true', help='不显示系列图的线性拟合线')
    p.add_argument('--no-show', action='store_true', help='只保存图片，不弹出交互式窗口')
    return p.parse_args()


def main():
    args = parse_args()

    if not DATA_FILE.exists():
        print(f"[ERROR] Data file not found: {DATA_FILE}")
        return
    if not MP_SUMMARY.exists():
        print(f"[ERROR] Melting point summary not found: {MP_SUMMARY}")
        return

    # 获取系列配置
    series_name = args.series
    config = SERIES_CONFIG[series_name]
    print(f"[INFO] Analyzing series: {series_name} ({config['title']})")

    df_all = pd.read_csv(DATA_FILE)
    df_mp = pd.read_csv(MP_SUMMARY)
    
    # 添加辅助列
    df_mp['metal_total'] = df_mp['Pt'] + df_mp['Sn']

    # 应用系列筛选
    mask = config['filter'](df_mp)
    df_series = df_mp[mask].copy()
    if df_series.empty:
        print(f'[WARN] No entries found for series: {series_name}')
        return

    # 应用屏蔽列表
    exclude_l = [s.lower() for s in args.exclude]
    if exclude_l:
        excluded = df_series[df_series['structure'].str.lower().isin(exclude_l)]['structure'].tolist()
        if excluded:
            print(f"[INFO] Excluding structures: {excluded}")
        df_series = df_series[~df_series['structure'].str.lower().isin(exclude_l)].copy()

    # 排序：按 x_col 排序，同 x_col 时按 Pt 比例降序（Pt 多的排前面）
    df_series['Pt_ratio'] = df_series['Pt'] / (df_series['Pt'] + df_series['Sn'] + 0.001)
    df_series = df_series.sort_values([config['x_col'], 'Pt_ratio'], ascending=[True, False])
    
    structures = list(df_series['structure'].unique())  # 保持排序顺序
    x_col = config['x_col']
    use_structure_labels = config.get('use_structure_labels', False)
    x_values = list(df_series[x_col].unique())
    print(f"Found {len(structures)} structures for {series_name}: {x_col} values = {list(x_values)}")
    
    # 显示排序后的结构列表
    print(f"  Sorted structures: {structures}")

    # Now extract delta per temp from main data for those structures
    df_sel = df_all[df_all['structure'].isin(structures)].copy()
    if df_sel.empty:
        print('[WARN] No raw lindemann records found for selected structures')
        return

    # group by structure and temp to compute mean delta
    df_mean = df_sel.groupby(['structure', 'temp']).agg(delta_mean=('delta', 'mean')).reset_index()

    # 判断是否使用结构名作为 x 轴（含氧系列）
    if use_structure_labels:
        # 使用结构名作为 x 轴，每个结构独立一列
        struct_to_idx = {s: i for i, s in enumerate(structures)}
        df_mean['x_val'] = df_mean['structure'].map(struct_to_idx)
        pivot_table = df_mean.pivot(index='temp', columns='x_val', values='delta_mean')
        pivot_table = pivot_table.sort_index(ascending=True)
        # x 轴标签使用格式化的结构名
        x_labels = [format_structure_label(s) for s in structures]
    else:
        # 原逻辑：按 x_col 值分组（可能会平均）
        struct_to_x = dict(zip(df_series['structure'], df_series[x_col]))
        df_mean['x_val'] = df_mean['structure'].map(struct_to_x)
        pivot = df_mean.groupby(['temp', 'x_val']).agg(delta_mean=('delta_mean', 'mean')).reset_index()
        pivot_table = pivot.pivot(index='temp', columns='x_val', values='delta_mean')
        pivot_table = pivot_table.sort_index(ascending=True)
        x_labels = [f'{int(c)}' for c in pivot_table.columns.values]

    # Prepare heatmap grid
    temps = pivot_table.index.values
    cols = pivot_table.columns.values  # x values
    data = pivot_table.values

    # 字体缩放
    fs = args.fontscale

    # ========================================================================
    # 绘制热图 (参考 step5_0 风格: RdYlBu_r, bilinear, 标注数值, 阈值线)
    # ========================================================================
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # 使用 imshow, bilinear 插值, RdYlBu_r 色图
    im = ax.imshow(data, aspect='auto', origin='lower', cmap=plt.cm.RdYlBu_r,
                   interpolation='bilinear', vmin=0, vmax=0.3)

    # 设置坐标轴 ticks
    ax.set_xticks(np.arange(len(cols)))
    ax.set_xticklabels(x_labels, fontsize=int(12*fs), rotation=45 if use_structure_labels else 0, ha='right' if use_structure_labels else 'center')
    ax.set_yticks(np.arange(len(temps)))
    ax.set_yticklabels([f'{int(t)}K' for t in temps], fontsize=int(12*fs))
    
    ax.set_xlabel(config['x_label'], fontsize=int(14*fs), fontweight='bold')
    ax.set_ylabel('Temperature (K)', fontsize=int(14*fs), fontweight='bold')
    if not args.no_title:
        ax.set_title(f"{config['title']}: Lindemann Index Heatmap", fontsize=int(16*fs), fontweight='bold')
    
    # 添加 colorbar (增大 pad 防止标签重叠)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.08)
    cbar.set_label('Lindemann Index δ', fontsize=int(14*fs))
    cbar.ax.tick_params(labelsize=int(12*fs))
    # 在 colorbar 上标记阈值线 (只画线，不加文字标注避免重叠)
    cbar.ax.axhline(0.1, color='black', linestyle='--', linewidth=2)

    # 在每个格子中标注数值 (可通过 --no-values 关闭)
    if not args.no_values:
        # 根据 RdYlBu_r 色图: 低值(蓝色)用白字, 中间值(黄/白色)用黑字, 高值(红色)用白字
        for i in range(len(temps)):
            for j in range(len(cols)):
                value = data[i, j]
                if not np.isnan(value):
                    # RdYlBu_r: 0→蓝, 0.1→浅蓝/白, 0.15→黄/白, 0.2→橙, 0.3→红
                    # 中间区域 (0.08-0.18) 是浅色/白色区域，用黑字
                    if 0.08 < value < 0.18:
                        text_color = 'black'
                    else:
                        text_color = 'white'
                    ax.text(j, i, f'{value:.{args.decimals}f}',
                           ha="center", va="center", color=text_color,
                           fontsize=int(10*fs), fontweight='bold')

    # ========================================================================
    # 在热图上画等高线 (表示熔点阈值)
    # ========================================================================
    threshold = args.threshold
    X, Y = np.meshgrid(np.arange(len(cols)), np.arange(len(temps)))
    # 绘制等高线 (contour) at δ=threshold
    cs = ax.contour(X, Y, data, levels=[threshold], colors=['black'], linewidths=2, linestyles='--')
    
    # 在 colorbar 上标记阈值
    cbar.ax.axhline(threshold, color='black', linestyle='--', linewidth=2)
    
    # 添加阈值信息到标题
    if not args.no_title:
        ax.set_title(f"{config['title']} - Lindemann Index (δ threshold = {threshold})",
                    fontsize=int(14*fs), fontweight='bold')

    plt.tight_layout()
    # 输出文件名包含阈值信息（仅当非默认值0.10时）
    thresh_suffix = f'_d{threshold:.2f}'.replace('.', '') if threshold != 0.10 else ''
    out1 = OUTPUT_DIR / f'{series_name}_delta_heatmap{thresh_suffix}.png'
    fig.savefig(out1, dpi=200, bbox_inches='tight')
    print(f"Saved heatmap: {out1}")
    if not args.no_show:
        plt.show()
    plt.close(fig)

    # ========================================================================
    # 系列图: Tm vs x_col (或结构名)
    # 使用自定义阈值重新计算熔点
    # ========================================================================
    print(f"[INFO] Calculating Tm with threshold δ={threshold}")
    tm_data = []
    for struct in structures:
        Tm, Tm_err = calculate_tm_from_threshold(df_all, struct, threshold)
        x_val = df_series[df_series['structure'] == struct][x_col].values[0]
        tm_data.append({'structure': struct, x_col: x_val, 'Tm': Tm, 'Tm_err': Tm_err})
    
    df_tm = pd.DataFrame(tm_data)
    df_tm = df_tm.dropna(subset=['Tm'])
    # 保持与热图一致的排序
    df_tm['_sort_idx'] = df_tm['structure'].map({s: i for i, s in enumerate(structures)})
    df_tm = df_tm.sort_values('_sort_idx')
    
    print(f"  Valid Tm data: {len(df_tm)} structures")

    fig, ax = plt.subplots(figsize=(8 if use_structure_labels else 7, 5))
    
    if use_structure_labels:
        # 含氧系列：x轴使用结构索引，标签使用格式化结构名
        x_vals = range(len(df_tm))
        ax.plot(x_vals, df_tm['Tm'], 'o-', 
                color='black', markersize=10, linewidth=2, 
                label=config['legend'])
        ax.set_xticks(x_vals)
        ax.set_xticklabels([format_structure_label(s) for s in df_tm['structure']], 
                          rotation=45, ha='right', fontsize=int(12*fs))
        # 含氧系列不做线性拟合（x轴是分类变量）
    else:
        # 非含氧系列：x轴使用数值
        ax.plot(df_tm[x_col], df_tm['Tm'], 'o-', 
                color='black', markersize=10, linewidth=2, 
                label=config['legend'])
        # 线性拟合 (灰色虚线) - 可通过 --no-fit 关闭
        if not args.no_fit and len(df_tm) > 2:
            slope, intercept, r, p, se = stats.linregress(df_tm[x_col], df_tm['Tm'])
            x_fit = np.linspace(df_tm[x_col].min(), df_tm[x_col].max(), 100)
            y_fit = slope * x_fit + intercept
            ax.plot(x_fit, y_fit, '--', color='gray', alpha=0.7, linewidth=1.5,
                   label=f'Fit: slope={slope:.1f} K/atom, R²={r**2:.3f}')

    ax.set_xlabel(config['x_label'], fontsize=int(14*fs), fontweight='bold')
    ax.set_ylabel('Melting Point $T_m$ (K)', fontsize=int(14*fs), fontweight='bold')
    ax.tick_params(axis='both', labelsize=int(12*fs))
    if not args.no_title:
        ax.set_title(f"{config['title']}: $T_m$ (δ={threshold})", fontsize=int(16*fs), fontweight='bold')
    ax.legend(fontsize=int(11*fs), loc='best')
    ax.grid(False)

    plt.tight_layout()
    out2 = OUTPUT_DIR / f'{series_name}_tm_vs_comp{thresh_suffix}.png'
    fig.savefig(out2, dpi=300, bbox_inches='tight')
    print(f"Saved series plot: {out2}")
    if not args.no_show:
        plt.show()
    plt.close(fig)

    # Print summary
    print('\nSummary:')
    print(f'  Series: {series_name}')
    print(f'  Threshold: δ={threshold}')
    print(f'  Structures plotted: {len(structures)}')
    print(f'  {x_col} values: {sorted(list(cols))}')
    if len(df_tm) > 0:
        print(f"  Mean Tm: {df_tm['Tm'].mean():.1f} K ± {df_tm['Tm'].std():.1f} K")


if __name__ == '__main__':
    main()