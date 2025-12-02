#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
Step 5.1.1: Pt8Snx 系列热图与系列图 (负载型, O=0)
================================================================================

功能概述
========
本脚本专门绘制 Pt8Snx 系列（负载型、不含氧）的：
1. Lindemann 指数热图：横轴 Sn 数量，纵轴温度，颜色表示 δ 值
2. 熔点 vs Sn 系列图：展示熔点随 Sn 含量的变化趋势

数据来源
========
1. Lindemann 原始数据: results/step6_0_multi_system/step6_0_all_systems_data.csv
2. 熔点汇总数据: results/step5_1_melting_point/melting_point_summary.csv
   (需先运行 step5_1_melting_point_analysis.py 生成)

输出文件
========
results/step5_1_1_pt8snx/
├── pt8snx_delta_heatmap.png   # Lindemann 指数热图
└── pt8snx_tm_vs_sn.png        # 熔点 vs Sn 系列图

命令行参数
==========
--exclude, -e    : 屏蔽指定结构，不参与绑图（可指定多个）
                   例: --exclude Pt8sn0 Pt8sn1
--decimals, -d   : 热图格子中数值的小数位数，可选 2 或 3（默认 3）
                   例: --decimals 2
--fontscale, -f  : 字体缩放比例（默认 1.0，1.2 表示放大 20%）
                   例: --fontscale 1.2
--no-values      : 不在热图格子中显示数值
--no-title       : 不显示图片标题
--no-show        : 只保存图片，不弹出交互式窗口

使用示例
========
# 基本使用（显示所有结构，3位小数，弹出窗口）
& C:/Users/11207/.conda/envs/vscode-1/python.exe step5_1_1_pt8snx_plot.py

# 屏蔽 Pt8sn0，只保存图片
& C:/Users/11207/.conda/envs/vscode-1/python.exe step5_1_1_pt8snx_plot.py --exclude Pt8sn0 --no-show

# 屏蔽 Pt8sn0，2位小数，不显示数值，不显示标题
& C:/Users/11207/.conda/envs/vscode-1/python.exe step5_1_1_pt8snx_plot.py -e Pt8sn0 -d 2 --no-values --no-title --no-show

# 放大字体 1.3 倍，不显示标题
& C:/Users/11207/.conda/envs/vscode-1/python.exe step5_1_1_pt8snx_plot.py -f 1.3 --no-title --no-show

# 屏蔽多个结构
& C:/Users/11207/.conda/envs/vscode-1/python.exe step5_1_1_pt8snx_plot.py --exclude Pt8sn0 Pt8sn1 Pt8sn2

注意事项
========
1. 结构名称大小写不敏感（Pt8sn0 和 PT8SN0 等效）
2. 热图使用 RdYlBu_r 色图，δ=0.1 处有黑色虚线标记熔点阈值
3. 系列图使用黑色实心圆点+实线，灰色虚线表示线性拟合

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


def parse_args():
    p = argparse.ArgumentParser(description='Plot Pt8Snx Lindemann heatmap and Tm vs Sn')
    p.add_argument('--exclude', '-e', nargs='+', default=[], help='屏蔽指定结构（可多个，大小写不敏感）')
    p.add_argument('--decimals', '-d', type=int, default=3, choices=[2, 3], help='热图数值小数位数 (2 或 3，默认 3)')
    p.add_argument('--fontscale', '-f', type=float, default=1.0, help='字体缩放比例 (默认 1.0)')
    p.add_argument('--no-values', action='store_true', help='不在热图格子中显示数值')
    p.add_argument('--no-title', action='store_true', help='不显示图片标题')
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

    df_all = pd.read_csv(DATA_FILE)
    df_mp = pd.read_csv(MP_SUMMARY)

    # Filter summary for Pt=8, O=0, supported
    df_pt8 = df_mp[(df_mp['Pt'] == 8) & (df_mp['O'] == 0) & (df_mp['type'] == 'supported')].copy()
    if df_pt8.empty:
        print('[WARN] No Pt8 supported entries found in melting point summary')
        return

    # 应用屏蔽列表
    exclude_l = [s.lower() for s in args.exclude]
    if exclude_l:
        excluded = df_pt8[df_pt8['structure'].str.lower().isin(exclude_l)]['structure'].tolist()
        if excluded:
            print(f"[INFO] Excluding structures: {excluded}")
        df_pt8 = df_pt8[~df_pt8['structure'].str.lower().isin(exclude_l)].copy()

    structures = sorted(df_pt8['structure'].unique())
    sn_values = sorted(df_pt8['Sn'].unique())
    print(f"Found {len(structures)} Pt8 supported structures: Sn values = {sn_values}")

    # Now extract delta per temp from main data for those structures
    df_sel = df_all[df_all['structure'].isin(structures)].copy()
    if df_sel.empty:
        print('[WARN] No raw lindemann records found for selected structures')
        return

    # group by structure and temp to compute mean delta
    df_mean = df_sel.groupby(['structure', 'temp']).agg(delta_mean=('delta', 'mean')).reset_index()

    # map structure -> Sn
    struct_to_sn = dict(zip(df_pt8['structure'], df_pt8['Sn']))
    df_mean['Sn'] = df_mean['structure'].map(struct_to_sn)

    # pivot: index temp, columns Sn, values delta_mean (average across structures with same Sn)
    pivot = df_mean.groupby(['temp', 'Sn']).agg(delta_mean=('delta_mean', 'mean')).reset_index()
    pivot_table = pivot.pivot(index='temp', columns='Sn', values='delta_mean')
    pivot_table = pivot_table.sort_index(ascending=True)

    # Prepare heatmap grid
    temps = pivot_table.index.values
    cols = pivot_table.columns.values  # Sn values
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
    ax.set_xticklabels([f'{int(c)}' for c in cols], fontsize=int(12*fs))
    ax.set_yticks(np.arange(len(temps)))
    ax.set_yticklabels([f'{int(t)}K' for t in temps], fontsize=int(12*fs))
    
    ax.set_xlabel('Number of Sn atoms', fontsize=int(14*fs), fontweight='bold')
    ax.set_ylabel('Temperature (K)', fontsize=int(14*fs), fontweight='bold')
    if not args.no_title:
        ax.set_title(r'Pt$_8$Sn$_x$: Lindemann Index Heatmap (supported, O=0)', fontsize=int(16*fs), fontweight='bold')
    
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
    # 在热图上画 δ=0.1 等高线 (表示熔点)
    # ========================================================================
    # 需要在连续坐标上插值来画等高线
    X, Y = np.meshgrid(np.arange(len(cols)), np.arange(len(temps)))
    # 绘制等高线 (contour) at δ=0.1
    cs = ax.contour(X, Y, data, levels=[0.1], colors=['black'], linewidths=2, linestyles='--')
    ax.clabel(cs, inline=True, fontsize=int(12*fs), fmt='δ=0.1')

    plt.tight_layout()
    out1 = OUTPUT_DIR / 'pt8snx_delta_heatmap.png'
    fig.savefig(out1, dpi=200, bbox_inches='tight')
    print(f"Saved heatmap: {out1}")
    if not args.no_show:
        plt.show()
    plt.close(fig)

    # ========================================================================
    # 系列图: Tm vs Sn (参考 Tm_series_analysis.png (a) 风格)
    # ========================================================================
    df_tm = df_pt8[['structure', 'Sn', 'Tm_lindemann', 'Tm_lindemann_err']].copy()
    df_tm = df_tm.dropna(subset=['Tm_lindemann']).sort_values('Sn')

    fig, ax = plt.subplots(figsize=(7, 5))
    
    # 黑色实心圆点+实线
    ax.plot(df_tm['Sn'], df_tm['Tm_lindemann'], 'o-', 
            color='black', markersize=10, linewidth=2, 
            label=r'Pt$_8$Sn$_x$')

    # 线性拟合 (灰色虚线)
    if len(df_tm) > 2:
        slope, intercept, r, p, se = stats.linregress(df_tm['Sn'], df_tm['Tm_lindemann'])
        x_fit = np.linspace(df_tm['Sn'].min(), df_tm['Sn'].max(), 100)
        y_fit = slope * x_fit + intercept
        ax.plot(x_fit, y_fit, '--', color='gray', alpha=0.7, linewidth=1.5,
               label=f'Fit: slope={slope:.1f} K/Sn, R²={r**2:.3f}')

    ax.set_xlabel('Number of Sn atoms', fontsize=int(14*fs), fontweight='bold')
    ax.set_ylabel('Melting Point $T_m$ (K)', fontsize=int(14*fs), fontweight='bold')
    ax.tick_params(axis='both', labelsize=int(12*fs))
    if not args.no_title:
        ax.set_title(r'Pt$_8$Sn$_x$: $T_m$ (Lindemann) vs Sn (supported, O=0)', fontsize=int(16*fs), fontweight='bold')
    ax.legend(fontsize=int(11*fs), loc='best')
    ax.grid(False)

    plt.tight_layout()
    out2 = OUTPUT_DIR / 'pt8snx_tm_vs_sn.png'
    fig.savefig(out2, dpi=300, bbox_inches='tight')
    print(f"Saved series plot: {out2}")
    if not args.no_show:
        plt.show()
    plt.close(fig)

    # Print summary
    print('\nSummary:')
    print(f'  Structures plotted: {len(structures)}')
    print(f'  Sn values: {sorted(list(cols))}')
    if len(df_tm) > 0:
        print(f"  Mean Tm: {df_tm['Tm_lindemann'].mean():.1f} K ± {df_tm['Tm_lindemann'].std():.1f} K")


if __name__ == '__main__':
    main()
