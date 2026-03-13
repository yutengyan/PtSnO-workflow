#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
Step 5.1.1.1: 多系列合并热图（带系列分隔线）
================================================================================

功能概述
========
本脚本绘制多系列合并的 Lindemann 指数热图：
- 支持 pt8snx, pt6snx, sum8 等系列的自由组合
- 系列之间用分隔线区分
- 保留热图中的数值显示和等高线标记
- 支持交互式调整图片尺寸、colorbar 位置等参数

支持的系列
==========
--only-series 参数可选值（逗号分隔）：
  pt8snx    : Pt8Snx 系列（Pt=8，无氧）
  pt6snx    : Pt6Snx 系列（Pt=6，无氧）
  sum8      : 总金属原子=8（Pt+Sn=8，无氧）
  air       : 气相合金（Air68, Air86）

输出文件
========
results/step5_1_1_1_multi_series/
└── multi_series_heatmap.png

命令行参数
==========
--only-series, -s : 只绘制指定系列（逗号分隔）
                    例: --only-series pt8snx,pt6snx,sum8
--exclude, -e     : 排除指定结构，支持简写格式
                    例: -e 80 60  (排除 Pt8Sn0, Pt6Sn0)
--separator       : 分隔线颜色（默认 black）
                    可选: black, white, none
--fontscale, -f   : 字体缩放比例（默认 1.0）
--decimals, -d    : 热图格子中数值的小数位数（默认 2）
--no-values       : 不在热图格子中显示数值
--no-title        : 不显示图片标题
--no-show         : 只保存图片，不弹出交互式窗口
--threshold, -t   : Lindemann 阈值（默认 0.10）
--interactive     : 开启交互式调整模式（可调整尺寸、colorbar位置等）
--figsize         : 图片尺寸，格式 WxH（默认 16x7.5）
--cbar-pad        : colorbar 与图片的间距（默认 0.02）

使用示例
========
# 绘制 pt8snx + pt6snx + sum8 三个系列
python step5_1_1_1_multi_series_heatmap.py --only-series pt8snx,pt6snx,sum8 -e 80 60 --no-title --no-show

# 只绘制 sum8 系列
python step5_1_1_1_multi_series_heatmap.py --only-series sum8 --no-title --no-show

# 交互式调整模式
python step5_1_1_1_multi_series_heatmap.py --only-series pt8snx,sum8 -e 80 --interactive

# 自定义图片尺寸
python step5_1_1_1_multi_series_heatmap.py --only-series pt8snx,sum8 --figsize 20x8

================================================================================
"""
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, TextBox
from scipy.interpolate import interp1d
import re

BASE_DIR = Path(__file__).parent
DATA_FILE = BASE_DIR / 'results' / 'step6_0_multi_system' / 'step6_0_all_systems_data.csv'
# 与 step6_1_3_lindemann_only 保持一致的 Pt8Sn6 数据源
DATA_FILE_50K = BASE_DIR / 'results' / 'step6_1_clustering' / 'Pt8Sn6_lindemann-threshold_n2_clustered_data.csv'
MP_SUMMARY = BASE_DIR / 'results' / 'step5_1_melting_point' / 'melting_point_summary.csv'
OUTPUT_DIR = BASE_DIR / 'results' / 'step5_1_1_1_multi_series'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# 英文字体设置（适合期刊发表）
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.unicode_minus'] = False


# 系列配置
SERIES_CONFIG = {
    'pt8snx': {
        'name': 'Pt$_8$Sn$_x$',
        'filter': lambda row: row['Pt'] == 8 and row['O'] == 0 and row['type'] != 'air',
        'sort_key': lambda row: (row['Sn'], -row['Pt']),
        'color': '#1f77b4'
    },
    'pt6snx': {
        'name': 'Pt$_6$Sn$_x$',
        'filter': lambda row: row['Pt'] == 6 and row['Sn'] > 0 and row['O'] == 0 and row['type'] != 'air' and (row['Pt'] + row['Sn']) != 8,
        'sort_key': lambda row: (row['Sn'], -row['Pt']),
        'color': '#ff7f0e'
    },
    'sum8': {
        'name': 'Pt$_{8-x}$Sn$_x$',
        'filter': lambda row: (row['Pt'] + row['Sn']) == 8 and row['Sn'] > 0 and row['O'] == 0 and row['type'] != 'air',
        'sort_key': lambda row: (row['Sn'], -row['Pt']),
        'color': '#2ca02c'
    },
    'air': {
        'name': 'Air',
        'filter': lambda row: row['type'] == 'air',
        'sort_key': lambda row: (row['Sn'], -row['Pt']),
        'color': '#d62728'
    }
}


def parse_composition(name):
    """解析结构名中的 Pt、Sn、O 原子数"""
    name_lower = name.lower()
    pt_match = re.search(r'pt(\d+)', name_lower)
    sn_match = re.search(r'sn(\d+)', name_lower)
    o_match = re.search(r'o(\d+)', name_lower)
    
    pt = int(pt_match.group(1)) if pt_match else 0
    sn = int(sn_match.group(1)) if sn_match else 0
    o = int(o_match.group(1)) if o_match else 0
    
    return pt, sn, o


def expand_structure_name(name):
    """扩展简写结构名，如 80 -> Pt8Sn0, 62 -> Pt6Sn2"""
    name = name.strip()
    if name.isdigit():
        if len(name) == 3:
            return f"Pt{name[0]}Sn{name[1]}O{name[2]}"
        elif len(name) == 2:
            return f"Pt{name[0]}Sn{name[1]}"
    return name


def should_exclude(struct_name, exclude_list):
    """检查结构是否应该被排除（基于 Pt、Sn、O 原子数匹配）"""
    if not exclude_list:
        return False
    struct_comp = parse_composition(struct_name)
    for excl in exclude_list:
        excl_expanded = expand_structure_name(excl)
        excl_comp = parse_composition(excl_expanded)
        # 只有当所有成分都匹配时才排除
        if excl_comp[0] > 0 or excl_comp[1] > 0:  # 至少有 Pt 或 Sn
            if excl_comp == struct_comp:
                return True
    return False


def format_structure_label(name):
    """格式化结构名为 (x,y) 或 (x,y,z) 形式"""
    pt, sn, o = parse_composition(name)
    if o > 0:
        return f'({pt},{sn},{o})'
    else:
        return f'({pt},{sn})'


def classify_structure(row):
    """
    对结构进行分类
    返回匹配的系列列表（一个结构可能属于多个系列）
    """
    matched = []
    for series_key, config in SERIES_CONFIG.items():
        if config['filter'](row):
            matched.append(series_key)
    return matched


def parse_args():
    """解析命令行参数"""
    p = argparse.ArgumentParser(
        description='多系列合并热图（带分隔线）',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例:
  python step5_1_1_1_multi_series_heatmap.py --only-series pt8snx,pt6snx,sum8 -e 80 60 --no-title --no-show
  python step5_1_1_1_multi_series_heatmap.py --only-series sum8 --no-title --no-show
''')
    
    p.add_argument('--only-series', '-s', type=str, default='pt8snx',
                   help='只绘制指定系列（逗号分隔）: pt8snx, pt6snx, sum8, air')
    p.add_argument('--exclude', '-e', nargs='+', default=[],
                   help='排除指定结构，如 80 60 (排除 Pt8Sn0, Pt6Sn0)')
    p.add_argument('--separator', type=str, default='black', choices=['black', 'white', 'none'],
                   help='系列分隔线颜色（默认 black）')
    p.add_argument('--fontscale', '-f', type=float, default=1.0,
                   help='字体缩放比例（默认 1.0）')
    p.add_argument('--decimals', '-d', type=int, default=2, choices=[2, 3],
                   help='热图数值小数位数（默认 2）')
    p.add_argument('--threshold', '-t', type=float, default=0.10,
                   help='Lindemann 阈值（默认 0.10）')
    p.add_argument('--no-values', action='store_true',
                   help='不显示热图格子中的数值')
    p.add_argument('--no-title', action='store_true',
                   help='不显示图片标题')
    p.add_argument('--no-show', action='store_true',
                   help='只保存图片，不弹出窗口')
    p.add_argument('--interactive', action='store_true',
                   help='开启交互式调整模式（可调整尺寸、colorbar位置等）')
    p.add_argument('--figsize', type=str, default='16x7.5',
                   help='图片尺寸，格式 WxH（默认 16x7.5）')
    p.add_argument('--cbar-pad', type=float, default=0.02,
                   help='colorbar 与图片的间距（默认 0.02）')
    
    # 字体大小精细调整
    p.add_argument('--tick-fontsize', type=int, default=None,
                   help='坐标轴刻度字体大小（默认 28*fontscale），同时控制 x/y 刻度')
    p.add_argument('--xtick-fontsize', type=int, default=None,
                   help='x 轴结构名刻度字体大小（默认与 --tick-fontsize 一致）')
    p.add_argument('--xtick-rotation', type=int, default=45,
                   help='x 轴刻度标签旋转角度，0=正着显示，45=倾斜（默认 45）')
    p.add_argument('--label-fontsize', type=int, default=None,
                   help='坐标轴标签字体大小（默认 34*fontscale）')
    p.add_argument('--value-fontsize', type=int, default=None,
                   help='热图数值字体大小（默认 9*fontscale）')
    p.add_argument('--cbar-fontsize', type=int, default=None,
                   help='colorbar 字体大小（默认 28*fontscale）')
    
    return p.parse_args()


def interactive_mode(params, out_path):
    """
    交互式调整模式
    允许用户调整图片尺寸、colorbar 位置等参数，并导出图片
    """
    fig = params['fig']
    ax = params['ax']
    im = params['im']
    cbar = params['cbar']
    
    # 当前参数
    current_params = {
        'fig_width': params['fig_w'],
        'fig_height': params['fig_h'],
        'cbar_pad': params['args'].cbar_pad,
        'dpi': 200
    }
    
    print("\n" + "=" * 60)
    print("交互式调整模式")
    print("=" * 60)
    print("使用滑块调整参数，点击 'Save' 保存图片")
    print("关闭窗口退出")
    print("=" * 60)
    
    # 调整图的布局，为控件留出空间
    plt.subplots_adjust(bottom=0.25)
    
    # 创建滑块
    ax_width = plt.axes([0.15, 0.15, 0.25, 0.03])
    ax_height = plt.axes([0.15, 0.10, 0.25, 0.03])
    ax_cbar_pad = plt.axes([0.55, 0.15, 0.25, 0.03])
    ax_dpi = plt.axes([0.55, 0.10, 0.25, 0.03])
    
    slider_width = Slider(ax_width, 'Width', 8, 30, valinit=current_params['fig_width'], valstep=0.5)
    slider_height = Slider(ax_height, 'Height', 4, 15, valinit=current_params['fig_height'], valstep=0.5)
    slider_cbar_pad = Slider(ax_cbar_pad, 'CBar Pad', 0.01, 0.15, valinit=current_params['cbar_pad'], valstep=0.01)
    slider_dpi = Slider(ax_dpi, 'DPI', 100, 400, valinit=current_params['dpi'], valstep=50)
    
    # 保存按钮
    ax_save = plt.axes([0.4, 0.02, 0.1, 0.04])
    btn_save = Button(ax_save, 'Save')
    
    # 重置按钮
    ax_reset = plt.axes([0.52, 0.02, 0.1, 0.04])
    btn_reset = Button(ax_reset, 'Reset')
    
    # 状态显示
    ax_status = plt.axes([0.15, 0.02, 0.2, 0.04])
    ax_status.set_xticks([])
    ax_status.set_yticks([])
    status_text = ax_status.text(0.5, 0.5, '', transform=ax_status.transAxes, 
                                  ha='center', va='center', fontsize=10)
    
    def update(val):
        """更新图片尺寸"""
        new_width = slider_width.val
        new_height = slider_height.val
        fig.set_size_inches(new_width, new_height)
        current_params['fig_width'] = new_width
        current_params['fig_height'] = new_height
        current_params['cbar_pad'] = slider_cbar_pad.val
        current_params['dpi'] = int(slider_dpi.val)
        fig.canvas.draw_idle()
    
    def save_figure(event):
        """保存图片 - 只保存热图部分，不包含控件"""
        dpi = int(slider_dpi.val)
        w = current_params['fig_width']
        h = current_params['fig_height']
        
        # 创建新图，只包含热图和 colorbar
        fig_save, ax_save_new = plt.subplots(figsize=(w, h))
        
        # 复制热图
        im_save = ax_save_new.imshow(params['data'], aspect='auto', origin='lower', 
                                      cmap=plt.cm.RdYlBu_r, interpolation='bilinear', 
                                      vmin=0, vmax=0.3)
        
        # 设置坐标轴
        ax_save_new.set_xticks(np.arange(len(params['cols'])))
        _rot = params['args'].xtick_rotation
        _ha  = 'center' if _rot == 0 else 'right'
        ax_save_new.set_xticklabels(params['x_labels'], fontsize=params['XTICK_FONTSIZE'],
                                     rotation=_rot, ha=_ha)
        ax_save_new.set_yticks(np.arange(len(params['temps'])))
        ax_save_new.set_yticklabels([f'{int(t)}' for t in params['temps']], 
                                     fontsize=params['TICK_FONTSIZE'])
        ax_save_new.set_xlabel('Pt$_x$Sn$_y$ (x,y)', fontsize=params['LABEL_FONTSIZE'])
        ax_save_new.set_ylabel('Temperature (K)', fontsize=params['LABEL_FONTSIZE'])
        
        # 添加数值标注
        if not params['args'].no_values:
            for i in range(len(params['temps'])):
                for j in range(len(params['cols'])):
                    value = params['data'][i, j]
                    if not np.isnan(value):
                        text_color = 'black' if 0.08 < value < 0.18 else 'white'
                        ax_save_new.text(j, i, f'{value:.{params["args"].decimals}f}',
                                        ha="center", va="center", color=text_color,
                                        fontsize=params['VALUE_FONTSIZE'], fontweight='bold')
        
        # 绘制等高线（分段）
        threshold = params['threshold']
        for seg_idx in range(len(params['all_boundaries']) - 1):
            start_col = params['all_boundaries'][seg_idx]
            end_col = params['all_boundaries'][seg_idx + 1]
            if end_col <= start_col:
                continue
            seg_data = params['contour_data'][:, start_col:end_col]
            seg_cols = np.arange(start_col, end_col)
            X_seg, Y_seg = np.meshgrid(seg_cols, params['contour_y'])
            ax_save_new.contour(X_seg, Y_seg, seg_data, levels=[threshold], 
                               colors=['black'], linewidths=2, linestyles='--')
        
        # 绘制分隔线
        if params['args'].separator != 'none' and params['series_boundaries']:
            sep_color = params['args'].separator
            for boundary in params['series_boundaries']:
                ax_save_new.axvline(x=boundary - 0.5, color=sep_color, linewidth=2, linestyle='-')
        
        # Colorbar
        cbar_save = fig_save.colorbar(im_save, ax=ax_save_new, fraction=0.046, 
                                       pad=current_params['cbar_pad'])
        cbar_save.set_label('Lindemann Index δ', fontsize=params['LABEL_FONTSIZE'])
        cbar_save.set_ticks([0, 0.1, 0.2, 0.3])
        cbar_save.ax.tick_params(labelsize=params['CBAR_FONTSIZE'])
        cbar_save.ax.axhline(threshold, color='black', linestyle='--', linewidth=2)
        
        fig_save.tight_layout()
        
        # 保存
        save_path = out_path.parent / f"{out_path.stem}_{w}x{h}_dpi{dpi}.png"
        fig_save.savefig(save_path, dpi=dpi, bbox_inches='tight')
        plt.close(fig_save)
        
        status_text.set_text(f'Saved: {save_path.name}')
        print(f"[SAVED] {save_path}")
        fig.canvas.draw_idle()
    
    def reset(event):
        """重置参数"""
        slider_width.reset()
        slider_height.reset()
        slider_cbar_pad.reset()
        slider_dpi.reset()
        status_text.set_text('Reset')
        fig.canvas.draw_idle()
    
    slider_width.on_changed(update)
    slider_height.on_changed(update)
    slider_cbar_pad.on_changed(update)
    slider_dpi.on_changed(update)
    btn_save.on_clicked(save_figure)
    btn_reset.on_clicked(reset)
    
    # 打印当前参数到控制台的函数
    def on_close(event):
        print(f"\n最终参数:")
        print(f"  --figsize {current_params['fig_width']}x{current_params['fig_height']}")
        print(f"  --cbar-pad {current_params['cbar_pad']}")
        print(f"  DPI: {current_params['dpi']}")
    
    fig.canvas.mpl_connect('close_event', on_close)
    
    plt.show()


def main():
    args = parse_args()
    
    print("=" * 60)
    print("Step 5.1.1.1: 多系列合并热图")
    print("=" * 60)
    
    # 解析系列参数
    target_series = [s.strip().lower() for s in args.only_series.split(',')]
    print(f"目标系列: {target_series}")
    print(f"排除结构: {args.exclude}")
    
    # 加载数据
    if not DATA_FILE.exists():
        print(f"[ERROR] 数据文件不存在: {DATA_FILE}")
        return
    if not MP_SUMMARY.exists():
        print(f"[ERROR] 熔点汇总文件不存在: {MP_SUMMARY}")
        return
    
    df_all = pd.read_csv(DATA_FILE)
    df_mp = pd.read_csv(MP_SUMMARY)    
    # 检查是否有与 step6_1_3 一致的 Pt8Sn6 50K 聚类数据
    df_pt8sn6_50k = None
    if DATA_FILE_50K.exists():
        print(f"  [OK] 发现 50K 步长数据：{DATA_FILE_50K.name}")
        df_50k_raw = pd.read_csv(DATA_FILE_50K)
        if not {'temp', 'delta'}.issubset(df_50k_raw.columns):
            print("  [WARN] 50K 文件缺少 temp/delta 列，将回退到 100K 数据")
            df_50k_raw = None
        if df_50k_raw is not None:
            # 该文件已是 Pt8Sn6 聚类结果，按 structure 字段再做一次稳健筛选
            df_pt8sn6_50k = df_50k_raw[df_50k_raw['structure'].astype(str).str.lower() == 'pt8sn6'].copy()
            if df_pt8sn6_50k.empty:
                df_pt8sn6_50k = df_50k_raw.copy()
            print(f"    - 筛选 Pt8Sn6 相关结构：{len(df_pt8sn6_50k)} 行数据")
            print(f"    - 映射后列名：{df_pt8sn6_50k.columns.tolist()}")
    else:
        print(f"  [INFO] 未找到 50K 步长数据，将使用 100K 插值")    
    # 热图显示使用 100K 网格（视觉更简洁）
    unified_temps = np.arange(200, 1200, 100)
    # 等值线计算使用更细温度网格，保证 Pt8Sn6 的 0.1 交点不因显示网格粗化而偏移
    contour_temps = np.arange(200, 1101, 10)
    
    fs = args.fontscale
    
    # 收集各系列的结构
    all_structures = []  # [(series_key, structure_name, row_data), ...]
    series_boundaries = []  # 系列边界索引
    
    for series_key in target_series:
        if series_key not in SERIES_CONFIG:
            print(f"[WARN] 未知系列: {series_key}")
            continue
        
        config = SERIES_CONFIG[series_key]
        
        # 筛选该系列的结构
        series_structures = []
        for _, row in df_mp.iterrows():
            if config['filter'](row):
                struct_name = row['structure']
                # 检查排除
                if should_exclude(struct_name, args.exclude):
                    print(f"  排除: {struct_name}")
                    continue
                series_structures.append((series_key, struct_name, row))
        
        # 按系列内排序规则排序
        series_structures.sort(key=lambda x: config['sort_key'](x[2]))
        
        if series_structures:
            if all_structures:
                series_boundaries.append(len(all_structures))
            all_structures.extend(series_structures)
            print(f"  {series_key}: {len(series_structures)} 个结构")
    
    if not all_structures:
        print("[ERROR] 没有找到符合条件的结构")
        return
    
    structures = [s[1] for s in all_structures]
    struct_series = {s[1]: s[0] for s in all_structures}
    
    print(f"\n总计: {len(structures)} 个结构")
    print(f"系列边界索引: {series_boundaries}")
    
    # 提取 Lindemann 数据并重采样
    df_sel = df_all[df_all['structure'].isin(structures)].copy()
    df_mean = df_sel.groupby(['structure', 'temp']).agg(delta_mean=('delta', 'mean')).reset_index()
    
    resampled_data = []
    struct_curves = {}
    for struct in structures:
        # Pt8Sn6 仅使用与 step6_1_3 一致的 50K 聚类数据
        if df_pt8sn6_50k is not None and 'pt8sn6' in struct.lower():
            df_struct = df_pt8sn6_50k[['temp', 'delta']].copy()
            df_struct = df_struct.groupby('temp').agg(delta_mean=('delta', 'mean')).reset_index()
            df_struct = df_struct.sort_values('temp')
            print(f"  [INFO] {struct}: 使用 step6_1_3 同源 50K 聚类数据 ({len(df_struct)} 温度点)")
        else:
            # 其他结构使用 100K 数据插值
            df_struct = df_mean[df_mean['structure'] == struct].sort_values('temp')

        if len(df_struct) >= 2:
            struct_curves[struct] = df_struct[['temp', 'delta_mean']].copy()
        
        if len(df_struct) >= 2:
            f = interp1d(df_struct['temp'], df_struct['delta_mean'],
                        kind='linear', bounds_error=False, fill_value=np.nan)
            for t in unified_temps:
                resampled_data.append({
                    'structure': struct,
                    'temp': t,
                    'delta_mean': f(t)
                })
    
    df_resampled = pd.DataFrame(resampled_data)
    struct_to_idx = {s: i for i, s in enumerate(structures)}
    df_resampled['x_val'] = df_resampled['structure'].map(struct_to_idx)
    pivot_table = df_resampled.pivot(index='temp', columns='x_val', values='delta_mean')
    pivot_table = pivot_table.sort_index(ascending=True)
    
    temps = pivot_table.index.values
    cols = pivot_table.columns.values
    data = pivot_table.values.astype(float)

    # 构建用于等值线的细网格数据（x: 结构列；y: 细温度网格）
    contour_data = np.full((len(contour_temps), len(cols)), np.nan)
    for struct in structures:
        if struct not in struct_curves:
            continue
        curve = struct_curves[struct]
        if len(curve) < 2:
            continue
        f_curve = interp1d(curve['temp'], curve['delta_mean'],
                           kind='linear', bounds_error=False, fill_value=np.nan)
        x_idx = struct_to_idx[struct]
        contour_data[:, x_idx] = f_curve(contour_temps)

    y_step = float(temps[1] - temps[0]) if len(temps) >= 2 else 100.0
    contour_y = (contour_temps - float(temps[0])) / y_step
    
    # 字体大小设置（与 5.1.2 一致的默认值）
    fs = args.fontscale
    # 解析图片尺寸
    try:
        fig_w, fig_h = map(float, args.figsize.lower().split('x'))
    except ValueError:
        print(f"[WARN] 无法解析图片尺寸 '{args.figsize}'，使用默认 16x7.5")
        fig_w, fig_h = 16, 7.5

    # 字体缩放：若用户未手动指定任何字号，则按 figsize 相对基准 16x7.5
    # 的几何平均比例自动缩放（保持字号与图面积视觉一致）
    BASE_W, BASE_H = 16.0, 7.5
    auto_scale = ((fig_w / BASE_W) * (fig_h / BASE_H)) ** 0.5  # 面积平方根比
    fs_effective = fs * auto_scale  # 在 --fontscale 基础上叠加尺寸自动缩放

    # 若用户手动指定了字号，则直接用；否则用自动缩放后的 fs_effective
    TICK_FONTSIZE  = args.tick_fontsize  if args.tick_fontsize  else int(28 * fs_effective)
    LABEL_FONTSIZE = args.label_fontsize if args.label_fontsize else int(34 * fs_effective)
    VALUE_FONTSIZE = args.value_fontsize if args.value_fontsize else int(9  * fs_effective)
    CBAR_FONTSIZE  = args.cbar_fontsize  if args.cbar_fontsize  else int(28 * fs_effective)
    # x 轴结构名刻度：默认与 TICK_FONTSIZE 一致，可单独用 --xtick-fontsize 覆盖
    XTICK_FONTSIZE = args.xtick_fontsize if args.xtick_fontsize else TICK_FONTSIZE

    print(f"\n字体大小:")
    print(f"  y刻度: {TICK_FONTSIZE}pt, x刻度: {XTICK_FONTSIZE}pt, 标签: {LABEL_FONTSIZE}pt, 数值: {VALUE_FONTSIZE}pt, colorbar: {CBAR_FONTSIZE}pt")
    print(f"  (自动缩放比: {auto_scale:.3f}, figsize: {fig_w}x{fig_h})")
    print(f"  colorbar 间距: {args.cbar_pad}")

    # 创建图
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    
    # 绘制热图
    im = ax.imshow(data, aspect='auto', origin='lower', cmap=plt.cm.RdYlBu_r,
                   interpolation='bilinear', vmin=0, vmax=0.3)
    
    # x 轴标签 - 使用与 5.1.2 一致的 (x,y) 格式
    x_labels = [format_structure_label(s) for s in structures]
    ax.set_xticks(np.arange(len(cols)))
    xtick_rot = args.xtick_rotation
    xtick_ha = 'center' if xtick_rot == 0 else 'right'
    ax.set_xticklabels(x_labels, fontsize=XTICK_FONTSIZE, rotation=xtick_rot, ha=xtick_ha)
    
    # y 轴（去掉K单位，与 5.1.2 一致）
    ax.set_yticks(np.arange(len(temps)))
    ax.set_yticklabels([f'{int(t)}' for t in temps], fontsize=TICK_FONTSIZE)
    
    # 轴标签
    ax.set_xlabel('Pt$_x$Sn$_y$ (x,y)', fontsize=LABEL_FONTSIZE)
    ax.set_ylabel('Temperature (K)', fontsize=LABEL_FONTSIZE)
    
    if not args.no_title:
        title_parts = [SERIES_CONFIG[s]['name'] for s in target_series if s in SERIES_CONFIG]
        ax.set_title(f"Lindemann Index: {' + '.join(title_parts)}", 
                    fontsize=int(34*fs), fontweight='bold')
    
    # 在每个格子中标注数值
    if not args.no_values:
        for i in range(len(temps)):
            for j in range(len(cols)):
                value = data[i, j]
                if not np.isnan(value):
                    # 根据颜色深浅选择文字颜色
                    if 0.08 < value < 0.18:
                        text_color = 'black'
                    else:
                        text_color = 'white'
                    ax.text(j, i, f'{value:.{args.decimals}f}',
                           ha="center", va="center", color=text_color,
                           fontsize=VALUE_FONTSIZE, fontweight='bold')
    
    # 绘制等高线（熔点阈值）- 分段绘制，不跨越系列边界
    threshold = args.threshold
    
    # 确定每个系列的列范围
    all_boundaries = [0] + series_boundaries + [len(cols)]
    
    for seg_idx in range(len(all_boundaries) - 1):
        start_col = all_boundaries[seg_idx]
        end_col = all_boundaries[seg_idx + 1]
        
        if end_col <= start_col:
            continue
        
        # 提取该系列用于等值线的细网格数据
        seg_data = contour_data[:, start_col:end_col]
        seg_cols = np.arange(start_col, end_col)
        
        X_seg, Y_seg = np.meshgrid(seg_cols, contour_y)
        
        # 绘制该段的等高线
        ax.contour(X_seg, Y_seg, seg_data, levels=[threshold], 
                  colors=['black'], linewidths=2, linestyles='--')
    
    # 绘制系列分隔线
    if args.separator != 'none' and series_boundaries:
        sep_color = args.separator
        for boundary in series_boundaries:
            # 在 boundary-0.5 处画线（介于两个格子之间）
            ax.axvline(x=boundary - 0.5, color=sep_color, linewidth=2, linestyle='-')
    
    # Colorbar - 靠近热图
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=args.cbar_pad)
    cbar.set_label('Lindemann Index δ', fontsize=LABEL_FONTSIZE)
    cbar.set_ticks([0, 0.1, 0.2, 0.3])
    cbar.ax.tick_params(labelsize=CBAR_FONTSIZE)
    cbar.ax.axhline(threshold, color='black', linestyle='--', linewidth=2)
    
    plt.tight_layout()
    
    # 保存参数（用于交互式模式）
    plot_params = {
        'fig': fig, 'ax': ax, 'im': im, 'cbar': cbar,
        'data': data, 'temps': temps, 'cols': cols,
        'structures': structures, 'x_labels': x_labels,
        'series_boundaries': series_boundaries,
        'threshold': threshold, 'args': args,
        'TICK_FONTSIZE': TICK_FONTSIZE, 'LABEL_FONTSIZE': LABEL_FONTSIZE,
        'VALUE_FONTSIZE': VALUE_FONTSIZE, 'CBAR_FONTSIZE': CBAR_FONTSIZE,
        'XTICK_FONTSIZE': XTICK_FONTSIZE,
        'target_series': target_series, 'all_boundaries': all_boundaries,
        'contour_data': contour_data, 'contour_y': contour_y,
        'fig_w': fig_w, 'fig_h': fig_h
    }
    
    # 保存
    series_suffix = '_'.join(target_series)
    out_path = OUTPUT_DIR / f'multi_series_heatmap_{series_suffix}.png'
    fig.savefig(out_path, dpi=200, bbox_inches='tight')
    print(f"\n[SAVED] {out_path}")
    
    # 交互式模式
    if args.interactive:
        interactive_mode(plot_params, out_path)
    elif not args.no_show:
        plt.show()
    
    plt.close(fig)
    
    print("\n" + "=" * 60)
    print("完成!")
    print("=" * 60)


if __name__ == '__main__':
    main()
