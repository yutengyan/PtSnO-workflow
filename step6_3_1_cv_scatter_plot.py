#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 6.3.1: Heat Capacity Scatter Plot - Cv1 vs Cv2

Creates a scatter plot comparing low-T and high-T heat capacity.
Labels use short format: (x,y) or (x,y,z) for Pt_x Sn_y O_z

Author: AI Assistant
Date: 2025-12-01

================================================================================
功能概述
================================================================================
本脚本绑制热容散点图，对比低温区 Cv1 和高温区 Cv2。
支持三种分类模式：
  1. 简单分类 (--classify simple)：气相、负载、负载含氧（默认）
  2. 系列分类 (--classify series)：按 Pt8Snx、Pt6Snx、其他负载、O1-O4 等分组
  3. 详细分类 (--classify detailed)：每个系列独立颜色和图例

体系分类说明
============
【气相合金】
  - Air系列: Air68, Air86 (Pt₆Sn₈, Pt₈Sn₆)

【负载Pt-Sn（无氧）】
  - Pt8Snx系列: Pt8Sn0~Pt8Sn10 (固定Pt=8，变Sn)
  - Pt6Snx系列: Pt6Sn1~Pt6Sn9 (固定Pt=6，变Sn)
  - 其他负载: Pt3Sn5, Pt5Sn3, Pt4Sn4, Pt7Sn1 等

【负载含氧】
  - O1系列: 1个氧原子 (Pt2Sn2O1, Pt3Sn2O1, ...)
  - O2系列: 2个氧原子 (Pt3Sn3O2, Pt4Sn6O2, ...)
  - O3系列: 3个氧原子 (Pt2Sn3O3, Pt5Sn7O3, ...)
  - O4系列: 4个氧原子 (Pt3Sn6O4, Pt6Sn8O4/Cv, ...)

================================================================================
命令行参数
================================================================================
--classify, -c   : 分类模式
                   - simple   : 简单3分类（气相/负载/含氧），默认
                   - series   : 系列分类（灵活组合）
                   - detailed : 详细分类（每系列独立）

--no-air         : 不绘制气相合金
--no-supported   : 不绘制负载Pt-Sn（无氧）
--no-oxide       : 不绘制负载含氧

--only-series    : 只绘制指定系列（逗号分隔，可自由组合）
                   可选系列:
                   - air       : 气相合金
                   - sum8      : 总原子数=8 (Pt+Sn=8，如Pt6Sn2,Pt4Sn4)
                   - pt8snx    : Pt=8 的系列 (Pt8Sn0,Pt8Sn2,...)
                   - pt6snx    : Pt=6 的系列 (Pt6Sn2,Pt6Sn8,...)
                   - o1,o2,o3,o4 : 按氧原子数分类
                   例: --only-series sum8,air
                   例: --only-series pt8snx,pt6snx,air
                   例: --only-series o1,o2,o3,o4

--exclude, -e    : 排除指定组分，格式 "x,y,z;x,y"
                   例: --exclude "8,0" 排除 Pt8Sn0
                   例: --exclude "8,0;6,0" 排除 Pt8Sn0 和 Pt6Sn0

--no-errorbars   : 不显示误差棒
--no-labels      : 不显示数据点标签
--fontscale, -f  : 字体缩放比例（默认1.0）
--markerscale, -m: 标记缩放比例（默认1.0）
-o, --output     : 输出文件路径

================================================================================
使用示例
================================================================================

# ★★★ 推荐命令 ★★★

# 1. 简单分类（默认）- 适合初步浏览
& C:/Users/11207/.conda/envs/vscode-1/python.exe step6_3_1_cv_scatter_plot.py --no-labels --no-errorbars --markerscale 2 --fontscale 1.5

# 2. 系列分类 - 区分 Pt8Snx、Pt6Snx、O1-O4
& C:/Users/11207/.conda/envs/vscode-1/python.exe step6_3_1_cv_scatter_plot.py --classify series --no-labels --no-errorbars --markerscale 2

# 3. 详细分类 - 每个子系列独立颜色
& C:/Users/11207/.conda/envs/vscode-1/python.exe step6_3_1_cv_scatter_plot.py --classify detailed --no-labels --no-errorbars --markerscale 2

# ★★★ 只绘制特定系列 ★★★

# 只绘制 Pt8Snx 和 Pt6Snx（对比固定Pt系列）
& C:/Users/11207/.conda/envs/vscode-1/python.exe step6_3_1_cv_scatter_plot.py --only-series pt8snx,pt6snx --no-labels --markerscale 2

# 只绘制含氧系列 O1-O4
& C:/Users/11207/.conda/envs/vscode-1/python.exe step6_3_1_cv_scatter_plot.py --only-series o1,o2,o3,o4 --no-labels --markerscale 2

# 只绘制气相合金
& C:/Users/11207/.conda/envs/vscode-1/python.exe step6_3_1_cv_scatter_plot.py --only-series air --no-labels --markerscale 2

# ★★★ 排除特定类型 ★★★

# 不绘制气相
& C:/Users/11207/.conda/envs/vscode-1/python.exe step6_3_1_cv_scatter_plot.py --no-air --no-labels --markerscale 2

# 只绘制负载无氧
& C:/Users/11207/.conda/envs/vscode-1/python.exe step6_3_1_cv_scatter_plot.py --no-air --no-oxide --no-labels --markerscale 2

# 只绘制负载含氧
& C:/Users/11207/.conda/envs/vscode-1/python.exe step6_3_1_cv_scatter_plot.py --no-air --no-supported --no-labels --markerscale 2

# ★★★ 排除特定组分 ★★★

# 排除 (3,5,3) 和 (6,8)
& C:/Users/11207/.conda/envs/vscode-1/python.exe step6_3_1_cv_scatter_plot.py --exclude "3,5,3;6,8" --no-labels --markerscale 2

# ★★★ 自定义输出 ★★★

# 保存到指定路径
& C:/Users/11207/.conda/envs/vscode-1/python.exe step6_3_1_cv_scatter_plot.py -o "results/cv_series_plot.png" --classify series

================================================================================
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
try:
    from adjustText import adjust_text
    HAS_ADJUSTTEXT = True
except ImportError:
    HAS_ADJUSTTEXT = False

# Font settings for academic journals - Arial (Nature/Science/ACS推荐)
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'sans-serif']
plt.rcParams['mathtext.fontset'] = 'dejavusans'  # 数学字体与 Arial 风格一致
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['font.size'] = 10


def parse_composition(name):
    """
    Parse structure name to extract Pt, Sn, O numbers
    
    Handles various formats:
    - Pt6sn8, Pt6sn8o4
    - O2pt4sn6, O2sn8pt7
    - Sn6pt5o2, Sn3o2pt2
    - Air68, Air86
    - Cv
    
    Returns: (pt_num, sn_num, o_num) or None
    """
    import re
    name_lower = name.lower()
    
    # Special case: Cv = Pt6Sn8O4
    if name == 'Cv':
        return (6, 8, 4)
    
    # Air series: Air68 = Pt6Sn8
    if 'air' in name_lower:
        match = re.search(r'air(\d+)(\d+)', name_lower)
        if match:
            return (int(match.group(1)), int(match.group(2)), 0)
        return None
    
    # Extract all components using flexible regex
    pt_match = re.search(r'pt(\d+)', name_lower)
    sn_match = re.search(r'sn(\d+)', name_lower)
    o_match = re.search(r'o(\d+)', name_lower)
    
    pt_num = int(pt_match.group(1)) if pt_match else 0
    sn_num = int(sn_match.group(1)) if sn_match else 0
    o_num = int(o_match.group(1)) if o_match else 0
    
    if pt_num > 0 or sn_num > 0:
        return (pt_num, sn_num, o_num)
    
    return None


def format_label_short(name):
    """
    Format structure name in short format for plot labels
    
    Rules:
    - 只显示 (x,y) 或 (x,y,z) 格式，不带系列前缀
    - 因为颜色/形状已经区分了系列，标签只需显示组成
    - Gas-phase: (x,y) e.g., (6,8)
    - Supported Pt-Sn: (x,y) e.g., (6,8)
    - Supported Pt-Sn-O: (x,y,z) e.g., (6,8,4)
    """
    comp = parse_composition(name)
    if comp is None:
        return name
    
    pt, sn, o = comp
    
    # 统一格式，不区分系列（颜色已区分）
    if o > 0:
        return f'({pt},{sn},{o})'
    else:
        return f'({pt},{sn})'


def classify_structure(name):
    """
    Classify structure type
    
    Returns: (type_key, display_name)
    """
    name_lower = name.lower()
    
    if 'air' in name_lower:
        return 'air', 'Gas-phase'
    
    # Use composition parser
    comp = parse_composition(name)
    if comp:
        pt, sn, o = comp
        if o > 0:
            return 'oxide', 'Supported Pt-Sn-O'
        else:
            return 'supported', 'Supported Pt-Sn'
    
    return 'other', 'Other'


def classify_structure_detailed(name):
    """
    Classify structure into detailed series
    
    Returns: dict with multiple classification keys
    
    Classification keys:
    - series_sum: 按总原子数分类 (sum8, sum6, other)
    - series_pt: 按Pt原子数分类 (pt8snx, pt6snx, other)
    - series_o: 按氧原子数分类 (o1, o2, o3, o4, no_oxide)
    - air: 是否气相
    
    使用时根据 only_series 参数选择匹配哪种分类
    """
    name_lower = name.lower()
    result = {}
    
    # 气相合金
    if 'air' in name_lower:
        result['air'] = True
        result['primary'] = 'air'
        return result
    
    result['air'] = False
    
    # 解析组分
    comp = parse_composition(name)
    if not comp:
        result['primary'] = 'other'
        return result
    
    pt, sn, o = comp
    
    # 按氧原子数分类 (o1, o2, o3, o4)
    # 含氧体系只归类到 o1-o4，不归到 sum8/pt6snx 等
    if o > 0:
        result[f'o{o}'] = True
        result['has_oxide'] = True
        result['primary'] = f'o{o}'
        return result  # 含氧体系直接返回，不再归入其他系列
    
    result['has_oxide'] = False
    
    # 以下只针对无氧体系
    # 按总原子数分类 (sum8)
    total = pt + sn
    if total == 8:
        result['sum8'] = True
    
    # 按Pt原子数分类 (pt8snx, pt6snx)
    # 注意：如果已经归入 sum8，则不再归入 pt6snx/pt8snx，避免重复显示
    if pt == 8 and total != 8:
        result['pt8snx'] = True
    if pt == 6 and total != 8:
        result['pt6snx'] = True
    
    # 确定主分类（用于默认显示）
    if total == 8:
        result['primary'] = 'sum8'
    elif pt == 8:
        result['primary'] = 'pt8snx'
    elif pt == 6:
        result['primary'] = 'pt6snx'
    else:
        result['primary'] = 'other_supported'
    
    return result


def match_series(classification, target_series):
    """
    检查分类结果是否匹配目标系列
    
    Args:
        classification: classify_structure_detailed 返回的分类字典
        target_series: 目标系列名 (如 'sum8', 'pt8snx', 'o1', 'air')
    
    Returns: True if matches
    """
    if target_series == 'air':
        return classification.get('air', False)
    
    # 直接检查分类字典中是否有该标记
    return classification.get(target_series, False)


# 系列样式配置
SERIES_STYLES = {
    # 气相: air-Pt_xSn_y
    'air': {'color': '#000000', 'marker': 'o', 'label': 'air-Pt$_x$Sn$_y$'},  # 黑色圆形
    
    # 按总原子数分类 (Pt+Sn=N)
    'sum8': {'color': '#000000', 'marker': 's', 'label': 'sup-Pt$_{8-x}$Sn$_x$'},  # 黑色方形
    
    # 按Pt原子数分类 (Pt=N)
    'pt8snx': {'color': '#000000', 'marker': '^', 'label': 'sup-Pt$_8$Sn$_x$'},  # 黑色上三角
    'pt6snx': {'color': '#000000', 'marker': 'v', 'label': 'sup-Pt$_6$Sn$_x$'},  # 黑色下三角
    
    'other_supported': {'color': '#808080', 'marker': 'p', 'label': 'sup-Pt$_x$Sn$_y$'},
    
    # 负载含氧: sup-Pt_xSn_yO_z (x,y,z) - 分开，使用不同标记和学术配色
    'o1': {'color': '#000000', 'marker': 'o', 'label': 'sup-Pt$_x$Sn$_y$O$_1$'},  # 黑色圆形
    'o2': {'color': '#000000', 'marker': 's', 'label': 'sup-Pt$_x$Sn$_y$O$_2$'},  # 黑色方形
    'o3': {'color': '#000000', 'marker': '^', 'label': 'sup-Pt$_x$Sn$_y$O$_3$'},  # 黑色上三角
    'o4': {'color': '#000000', 'marker': 'D', 'label': 'sup-Pt$_x$Sn$_y$O$_4$'},  # 黑色菱形
    
    # 负载含氧: sup-Pt_xSn_yO_z (x,y,z) - 合并为一个系列
    'oxide_all': {'color': '#E74C3C', 'marker': 'o', 'label': 'sup-Pt$_x$Sn$_y$O$_z$'},
    
    # 简单分类样式
    'supported': {'color': '#2E8B57', 'marker': 's', 'label': 'sup-Pt$_x$Sn$_y$'},
    'oxide': {'color': '#E74C3C', 'marker': 'o', 'label': 'sup-Pt$_x$Sn$_y$O$_z$'},
}


def load_partition_data(structure, base_dir='results/step6_1_clustering'):
    """Load Cv1 and Cv2 from quality metrics file"""
    filepath = os.path.join(base_dir, f'{structure}_auto2_quality_metrics.csv')
    
    if not os.path.exists(filepath):
        return None
    
    try:
        df = pd.read_csv(filepath)
        
        if len(df) < 2:
            return None
        
        return {
            'cv1': df['Cv_cluster'].iloc[0],
            'err1': df['Cv_cluster_err'].iloc[0],
            'cv2': df['Cv_cluster'].iloc[1],
            'err2': df['Cv_cluster_err'].iloc[1],
        }
    except Exception as e:
        print(f"  Warning: Failed to read {structure}: {e}")
        return None


def find_all_structures(base_dir='results/step6_1_clustering'):
    """Find all available structures"""
    pattern = os.path.join(base_dir, '*_auto2_quality_metrics.csv')
    files = glob.glob(pattern)
    
    structures = []
    for f in files:
        basename = os.path.basename(f)
        structure = basename.replace('_auto2_quality_metrics.csv', '')
        structures.append(structure)
    
    return sorted(structures)


def should_exclude(name, exclude_list):
    """
    Check if a structure should be excluded based on composition
    
    Args:
        name: structure name
        exclude_list: list of tuples like [(3,5,3), (3,4,1)] to exclude
    
    Returns: True if should be excluded
    """
    if not exclude_list:
        return False
    
    comp = parse_composition(name)
    if comp is None:
        return False
    
    pt, sn, o = comp
    
    for exc in exclude_list:
        if len(exc) == 2:
            # (x, y) format - match Pt and Sn
            if pt == exc[0] and sn == exc[1]:
                return True
        elif len(exc) == 3:
            # (x, y, z) format - match Pt, Sn, O
            if pt == exc[0] and sn == exc[1] and o == exc[2]:
                return True
    
    return False


def create_cv_scatter_plot(output_path='results/step6_1_clustering/cv1_vs_cv2_scatter.png',
                           show_air=True, show_supported=True, show_oxide=True,
                           exclude_list=None, show_errorbars=True, show_labels=True,
                           fontscale=1.0, markerscale=1.0,
                           classify_mode='simple', only_series=None, merge_oxide=False,
                           interactive=False, no_stroke=False, figsize=(8, 8)):
    """
    Create Cv1 vs Cv2 scatter plot with non-overlapping labels
    
    Args:
        output_path: output file path
        show_air: whether to show gas-phase data
        show_supported: whether to show supported Pt-Sn (no oxygen) data
        show_oxide: whether to show supported Pt-Sn-O data
        exclude_list: list of compositions to exclude, e.g., [(3,5,3), (3,4,1)]
        show_errorbars: whether to show error bars (default True)
        show_labels: whether to show data point labels (default True)
        fontscale: scale factor for all fonts (default 1.0)
        markerscale: scale factor for marker sizes (default 1.0)
        classify_mode: 'simple', 'series', or 'detailed'
        only_series: list of series keys to plot (e.g., ['pt8snx', 'pt6snx'])
        merge_oxide: whether to merge O1-O4 into one series (default False)
        interactive: whether to enable interactive label dragging (default False)
        no_stroke: whether to disable white stroke around labels (default False)
        figsize: figure size as tuple (width, height), default (8, 8)
    """
    
    # Base font sizes - 大字体配合大图片 (20:15 比例)
    FONT_LABEL = 26 * fontscale      # data point labels (26pt)
    FONT_AXIS = 34 * fontscale       # axis labels (坐标轴标签)
    FONT_TITLE = 36 * fontscale      # title
    FONT_LEGEND = 26 * fontscale     # legend (图注)
    FONT_TICK = 28 * fontscale       # tick labels (坐标轴数字)
    STROKE_WIDTH = 5 * fontscale     # text stroke width (增大白边避免遮挡)
    
    # Base marker size (will be scaled)
    MARKER_SIZE = 200 * markerscale
    
    base_dir = 'results/step6_1_clustering'
    structures = find_all_structures(base_dir)
    
    print(f"Found {len(structures)} systems")
    
    # Print filter settings
    print(f"\nFilter settings:")
    print(f"  Classify mode: {classify_mode}")
    print(f"  Merge O1-O4: {merge_oxide}")
    print(f"  Show Gas-phase: {show_air}")
    print(f"  Show Supported Pt-Sn: {show_supported}")
    print(f"  Show Supported Pt-Sn-O: {show_oxide}")
    print(f"  Show Error bars: {show_errorbars}")
    print(f"  Show Labels: {show_labels}")
    print(f"  Font scale: {fontscale}")
    print(f"  Marker scale: {markerscale}")
    if only_series:
        print(f"  Only series: {only_series}")
    if exclude_list:
        print(f"  Excluded compositions: {exclude_list}")
    
    # 根据分类模式选择分类函数和样式
    if classify_mode == 'simple':
        # 简单3分类
        all_series_keys = ['air', 'supported', 'oxide']
    else:
        # 详细分类 - 根据 only_series 动态确定要显示的系列
        # 支持: air, sum8, sum6, pt8snx, pt6snx, o1, o2, o3, o4, oxide_all
        if only_series:
            # 使用用户指定的系列
            if merge_oxide:
                # 把 o1-o4 映射到 oxide_all
                all_series_keys = []
                for s in only_series:
                    if s in ['o1', 'o2', 'o3', 'o4']:
                        if 'oxide_all' not in all_series_keys:
                            all_series_keys.append('oxide_all')
                    else:
                        if s not in all_series_keys:
                            all_series_keys.append(s)
            else:
                all_series_keys = list(only_series)
        else:
            # 默认显示所有系列
            if merge_oxide:
                all_series_keys = ['air', 'sum8', 'pt8snx', 'pt6snx', 'other_supported', 'oxide_all']
            else:
                all_series_keys = ['air', 'sum8', 'pt8snx', 'pt6snx', 'other_supported', 'o1', 'o2', 'o3', 'o4']
    
    # 加载并分类数据
    data_by_series = {k: [] for k in all_series_keys}
    excluded_count = 0
    
    for struct in structures:
        # 检查排除列表
        if should_exclude(struct, exclude_list):
            comp = parse_composition(struct)
            print(f"    Excluding: {struct} -> {comp}")
            excluded_count += 1
            continue
        
        d = load_partition_data(struct, base_dir)
        if d is None:
            continue
        
        # 获取分类信息
        if classify_mode == 'simple':
            series_key, _ = classify_structure(struct)
            matched_series = [series_key]
        else:
            classification = classify_structure_detailed(struct)
            # 找到该结构匹配的所有系列
            matched_series = []
            for target in all_series_keys:
                if target == 'oxide_all':
                    # oxide_all 匹配任何含氧体系
                    if classification.get('has_oxide', False):
                        matched_series.append('oxide_all')
                elif match_series(classification, target):
                    matched_series.append(target)
        
        # 检查类型过滤
        if classify_mode != 'simple':
            if classification.get('air', False) and not show_air:
                continue
            if not classification.get('air', False) and not classification.get('has_oxide', False) and not show_supported:
                continue
            if classification.get('has_oxide', False) and not show_oxide:
                continue
        
        short_name = format_label_short(struct)
        
        # 添加到所有匹配的系列中
        for series_key in matched_series:
            if series_key in data_by_series:
                data_by_series[series_key].append({
                    'name': struct,
                    'display': short_name,
                    'cv1': d['cv1'],
                    'cv2': d['cv2'],
                    'err1': d['err1'],
                    'err2': d['err2'],
                })
    
    if excluded_count > 0:
        print(f"  Total excluded: {excluded_count}")
    
    # 打印统计
    print(f"\nClassification statistics (after filtering):")
    for key in all_series_keys:
        if data_by_series[key]:
            print(f"  {SERIES_STYLES.get(key, {}).get('label', key)}: {len(data_by_series[key])}")
    
    # 确定要绘制的系列
    series_to_plot = [k for k in all_series_keys if data_by_series[k]]
    
    if not series_to_plot:
        print("No data to plot!")
        return None
    
    # Create figure - 正方形数据配正方形图，减少空白
    # 由于 Cv1 vs Cv2 是对角线对称图，使用接近正方形的尺寸
    fig, ax = plt.subplots(figsize=figsize)
    
    # Increase tick label font size (scaled)
    ax.tick_params(axis='both', labelsize=FONT_TICK)
    
    all_points = []
    
    # Plot data points by series
    for series_key in series_to_plot:
        data = data_by_series[series_key]
        if not data:
            continue
        
        s = SERIES_STYLES.get(series_key, {'color': 'gray', 'marker': 'o', 'label': series_key})
        
        cv1_vals = [d['cv1'] for d in data]
        cv2_vals = [d['cv2'] for d in data]
        err1_vals = [d['err1'] for d in data]
        err2_vals = [d['err2'] for d in data]
        
        # Error bars (optional)
        if show_errorbars:
            ax.errorbar(cv1_vals, cv2_vals, xerr=err1_vals, yerr=err2_vals,
                        fmt='none', ecolor=s['color'], alpha=0.3, capsize=2, zorder=5)
        
        # Scatter points - 增加高亮边框（黑色边框让点更突出）
        ax.scatter(cv1_vals, cv2_vals, c=s['color'], marker=s['marker'],
                   s=MARKER_SIZE, alpha=0.95, label=s['label'], 
                   edgecolors='black', linewidths=1.5, zorder=15)
        
        # Collect label data - 使用唯一ID (series_key + name) 来标识每个点
        for d in data:
            unique_id = f"{series_key}_{d['name']}"  # 唯一标识符
            all_points.append((d['cv1'], d['cv2'], d['display'], s['color'], series_key, unique_id, d['name']))
    
    # Add diagonal line (Cv1 = Cv2)
    lims = [
        min(ax.get_xlim()[0], ax.get_ylim()[0]),
        max(ax.get_xlim()[1], ax.get_ylim()[1])
    ]
    ax.plot(lims, lims, 'k--', alpha=0.4, linewidth=1.5, label=r'$C_{v,1} = C_{v,2}$', zorder=1)
    
    # 固定刻度为 0, 2, 4, 6, 8 - 简洁对称，从0开始无空白
    ticks = np.array([0, 2, 4, 6, 8])
    ax_min = 0
    ax_max = 8
    
    ax.set_xlim(ax_min, ax_max)
    ax.set_ylim(ax_min, ax_max)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    
    # 处理原点0重叠：隐藏Y轴的0刻度标签
    yticklabels = [str(int(t)) if t != 0 else '' for t in ticks]
    ax.set_yticklabels(yticklabels)
    
    # 标签位置文件路径
    label_positions_file = os.path.splitext(output_path)[0] + '_label_positions.json'
    
    # 尝试加载已保存的标签位置
    saved_positions = {}
    if os.path.exists(label_positions_file):
        try:
            import json
            with open(label_positions_file, 'r', encoding='utf-8') as f:
                saved_positions = json.load(f)
            print(f"📂 已加载保存的标签位置: {label_positions_file}")
        except Exception as e:
            print(f"⚠️ 无法加载标签位置文件: {e}")
    
    # Add labels with adjustText (optional)
    texts = []
    text_to_uid = {}   # text对象到唯一ID的映射
    uid_to_text = {}   # 唯一ID到text对象的映射
    if show_labels:
        for x, y, label, color, series_key, unique_id, name in all_points:
            # 检查是否有保存的位置（用唯一ID查找）
            if unique_id in saved_positions:
                pos_x, pos_y = saved_positions[unique_id]
            else:
                pos_x, pos_y = x, y
            
            # 根据 no_stroke 参数决定是否加白边
            if no_stroke:
                txt = ax.text(pos_x, pos_y, label, fontsize=FONT_LABEL, color='black', 
                              ha='center', va='center', zorder=20, picker=True)
            else:
                txt = ax.text(pos_x, pos_y, label, fontsize=FONT_LABEL, color='black', 
                              ha='center', va='center',
                              path_effects=[path_effects.withStroke(linewidth=STROKE_WIDTH, foreground='white')],
                              zorder=20, picker=True)
            texts.append(txt)
            text_to_uid[txt] = unique_id
            uid_to_text[unique_id] = txt
        
        # Get scatter point positions for adjustText to avoid
        x_points = [p[0] for p in all_points]
        y_points = [p[1] for p in all_points]
        
        # 使用 adjustText 自动优化标签位置（如果没有已保存的位置）
        # 交互模式下也会先自动优化，然后再允许手动微调
        if HAS_ADJUSTTEXT and not saved_positions:
            print("🔧 正在使用 adjustText 自动优化标签位置...")
            # 注意：arrowprops=None 禁止 adjustText 画线，我们统一管理连接线
            adjust_text(texts, 
                        x=x_points, y=y_points,
                        ax=ax,
                        arrowprops=None,            # 禁止自动画线，我们统一管理
                        expand_points=(6.0, 6.0),   # 更大的点周围空间
                        expand_text=(4.0, 4.0),     # 更大的文字间距
                        force_text=(3.0, 3.0),      # 更强文字排斥力
                        force_points=(3.0, 3.0),    # 更强点排斥力
                        lim=2000,                   # 更多迭代确保收敛
                        only_move={'points': 'xy', 'texts': 'xy'})  # 允许任意方向移动
            print("✅ adjustText 优化完成")
    
    # Labels (带单位，不要标题)
    ax.set_xlabel(r'$C_{v,1}$ (meV/K)', fontsize=FONT_AXIS, fontweight='bold')
    ax.set_ylabel(r'$C_{v,2}$ (meV/K)', fontsize=FONT_AXIS, fontweight='bold')
    # 不设置标题
    
    # Legend - 无框框
    ncol = 1 if len(series_to_plot) <= 4 else 2
    ax.legend(loc='lower right', fontsize=FONT_LEGEND, frameon=False,
              handlelength=1.5, handletextpad=0.5, borderpad=0.4, ncol=ncol)
    
    ax.set_aspect('equal', adjustable='box')
    
    # 减少边距
    plt.tight_layout(pad=0.5)
    
    # 建立唯一ID到数据点坐标的映射
    uid_to_point = {unique_id: (x, y) for x, y, label, color, series_key, unique_id, name in all_points}
    
    # 绘制连接线（从数据点到标签）- 无论交互还是非交互模式都需要
    lines = {}  # uid -> line对象
    if show_labels:
        for txt in texts:
            uid = text_to_uid.get(txt)
            if uid and uid in uid_to_point:
                px, py = uid_to_point[uid]
                tx, ty = txt.get_position()
                # 判断标签是否移动了（阈值设小一点，让更多标签有连接线）
                dist = ((tx - px)**2 + (ty - py)**2)**0.5
                if dist > 0.02:  # 非常小的阈值，几乎所有移动都会显示连接线
                    line, = ax.plot([px, tx], [py, ty], 'gray', alpha=0.5, lw=0.8, zorder=5)
                    lines[uid] = line
    
    # 交互模式：允许拖动标签
    if interactive and show_labels:
        import json
        
        print("\n" + "=" * 60)
        print("交互模式 - 可拖动标签调整位置")
        print("=" * 60)
        print("操作说明:")
        print("  - 鼠标左键拖动标签到新位置")
        print("  - 按 S 键保存图片和标签位置")
        print("  - 按 R 键重置所有标签到数据点位置")
        print("  - 按 Q 键或关闭窗口退出")
        if saved_positions:
            print(f"  📂 已加载 {len(saved_positions)} 个保存的标签位置")
        print("=" * 60)
        
        # 拖动状态
        drag_state = {'text': None, 'offset': (0, 0), 'uid': None}
        
        def on_pick(event):
            if event.artist in texts:
                drag_state['text'] = event.artist
                drag_state['uid'] = text_to_uid.get(event.artist)
                # 计算鼠标到文本中心的偏移
                x0, y0 = drag_state['text'].get_position()
                drag_state['offset'] = (x0 - event.mouseevent.xdata, y0 - event.mouseevent.ydata)
        
        def on_motion(event):
            if drag_state['text'] is not None and event.xdata is not None:
                new_x = event.xdata + drag_state['offset'][0]
                new_y = event.ydata + drag_state['offset'][1]
                drag_state['text'].set_position((new_x, new_y))
                
                # 更新连接线
                uid = drag_state['uid']
                if uid and uid in uid_to_point:
                    px, py = uid_to_point[uid]
                    if uid in lines:
                        lines[uid].set_data([px, new_x], [py, new_y])
                    else:
                        # 创建新连接线
                        line, = ax.plot([px, new_x], [py, new_y], 'gray', alpha=0.5, lw=0.8, zorder=5)
                        lines[uid] = line
                
                fig.canvas.draw_idle()
        
        def on_release(event):
            drag_state['text'] = None
            drag_state['uid'] = None
        
        def save_positions():
            """保存所有标签位置到JSON文件（使用唯一ID作为key）"""
            positions = {}
            for txt in texts:
                uid = text_to_uid.get(txt)
                if uid:
                    x, y = txt.get_position()
                    positions[uid] = [x, y]
            
            with open(label_positions_file, 'w', encoding='utf-8') as f:
                json.dump(positions, f, indent=2, ensure_ascii=False)
            print(f"💾 标签位置已保存: {label_positions_file}")
        
        def on_key(event):
            if event.key.lower() == 's':
                # 保存图片
                plt.savefig(output_path, dpi=300, bbox_inches='tight')
                print(f"\n✅ Figure saved: {output_path}")
                # 保存标签位置
                save_positions()
            elif event.key.lower() == 'r':
                # 重置所有标签到数据点位置
                for txt in texts:
                    uid = text_to_uid.get(txt)
                    if uid and uid in uid_to_point:
                        px, py = uid_to_point[uid]
                        txt.set_position((px, py))
                        # 移除连接线
                        if uid in lines:
                            lines[uid].remove()
                            del lines[uid]
                fig.canvas.draw_idle()
                print("🔄 已重置所有标签位置")
            elif event.key.lower() == 'q':
                plt.close()
        
        fig.canvas.mpl_connect('pick_event', on_pick)
        fig.canvas.mpl_connect('motion_notify_event', on_motion)
        fig.canvas.mpl_connect('button_release_event', on_release)
        fig.canvas.mpl_connect('key_press_event', on_key)
        
        plt.show()
    else:
        # 非交互模式：已经在上面绘制了连接线，这里不需要重复
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"\n✅ Figure saved: {output_path}")
    
    # 导出数据供 Origin 使用
    export_data_for_origin(all_points, data_by_series, output_path)
    
    return output_path


def export_data_for_origin(all_points, data_by_series, output_path):
    """
    导出数据为 CSV 格式，方便在 Origin 中绑图
    
    输出文件:
    1. {output_path}_data.csv - 所有数据点
    2. {output_path}_by_series.csv - 按系列分组的数据
    """
    import os
    
    base_path = os.path.splitext(output_path)[0]
    
    # 1. 导出所有数据点
    # all_points 格式: (cv1, cv2, display, color, series_key, unique_id, name)
    data_rows = []
    for point in all_points:
        cv1, cv2, display, color, series_key = point[0], point[1], point[2], point[3], point[4]
        data_rows.append({
            'Label': display,
            'Cv1': cv1,
            'Cv2': cv2,
            'Series': series_key,
            'Color': color
        })
    
    df_all = pd.DataFrame(data_rows)
    csv_path = f"{base_path}_origin_data.csv"
    df_all.to_csv(csv_path, index=False)
    print(f"✅ Origin data exported: {csv_path}")
    
    # 2. 按系列分组导出（更适合 Origin 绑图）
    # Origin 喜欢每列一个系列的格式
    series_data = {}
    max_len = 0
    
    for series_key, data in data_by_series.items():
        if data:
            series_data[series_key] = data
            max_len = max(max_len, len(data))
    
    # 创建宽格式数据框
    wide_data = {}
    for series_key, data in series_data.items():
        style = SERIES_STYLES.get(series_key, {})
        label = style.get('label', series_key).replace('$', '').replace('_', '').replace('{', '').replace('}', '')
        
        cv1_col = f'{series_key}_Cv1'
        cv2_col = f'{series_key}_Cv2'
        name_col = f'{series_key}_Name'
        
        wide_data[cv1_col] = [d['cv1'] for d in data] + [None] * (max_len - len(data))
        wide_data[cv2_col] = [d['cv2'] for d in data] + [None] * (max_len - len(data))
        wide_data[name_col] = [d['display'] for d in data] + [None] * (max_len - len(data))
    
    df_wide = pd.DataFrame(wide_data)
    csv_wide_path = f"{base_path}_origin_wide.csv"
    df_wide.to_csv(csv_wide_path, index=False)
    print(f"✅ Origin wide format exported: {csv_wide_path}")
    
    # 3. 导出 Origin 绑图设置说明
    settings_path = f"{base_path}_origin_settings.txt"
    with open(settings_path, 'w', encoding='utf-8') as f:
        f.write("=" * 60 + "\n")
        f.write("Origin 绑图设置说明\n")
        f.write("=" * 60 + "\n\n")
        
        f.write("【数据文件】\n")
        f.write(f"  长格式: {os.path.basename(csv_path)}\n")
        f.write(f"  宽格式: {os.path.basename(csv_wide_path)}\n\n")
        
        f.write("【图形类型】\n")
        f.write("  Scatter Plot (XY散点图)\n\n")
        
        f.write("【坐标轴设置】\n")
        f.write("  X轴: Cv1 (Low-T)\n")
        f.write("  Y轴: Cv2 (High-T)\n")
        f.write("  坐标轴标签字体: Arial Bold, 34pt\n")
        f.write("  刻度数字字体: Arial, 28pt\n")
        f.write("  刻度: 4-7个，对称\n\n")
        
        f.write("【系列颜色和符号】\n")
        for series_key, style in SERIES_STYLES.items():
            if series_key in series_data:
                f.write(f"  {series_key}:\n")
                f.write(f"    颜色: {style['color']}\n")
                f.write(f"    符号: {style['marker']}\n")
                f.write(f"    标签: {style['label']}\n\n")
        
        f.write("【图例设置】\n")
        f.write("  位置: 右下角\n")
        f.write("  字体: Arial, 26pt\n")
        f.write("  无边框\n\n")
        
        f.write("【对角线】\n")
        f.write("  y = x 虚线, 黑色, alpha=0.4\n\n")
        
        f.write("【符号对应表】\n")
        f.write("  o = 圆形 (Circle)\n")
        f.write("  s = 方形 (Square)\n")
        f.write("  D = 菱形 (Diamond)\n")
        f.write("  p = 五角形 (Pentagon)\n")
        f.write("  ^ = 上三角 (Triangle Up)\n")
        f.write("  v = 下三角 (Triangle Down)\n")
        f.write("  < = 左三角 (Triangle Left)\n")
        f.write("  > = 右三角 (Triangle Right)\n")
    
    print(f"✅ Origin settings exported: {settings_path}")


def parse_exclude_arg(exclude_str):
    """
    Parse exclude argument string to list of tuples
    
    Examples:
        "3,5,3" -> [(3,5,3)]
        "3,5,3;3,4,1" -> [(3,5,3), (3,4,1)]
        "6,8" -> [(6,8)]
    """
    if not exclude_str:
        return None
    
    exclude_list = []
    for item in exclude_str.split(';'):
        parts = item.strip().split(',')
        try:
            nums = tuple(int(p.strip()) for p in parts)
            if len(nums) in [2, 3]:
                exclude_list.append(nums)
        except ValueError:
            print(f"  Warning: Could not parse '{item}'")
    
    return exclude_list if exclude_list else None


def main():
    """Main function with command-line arguments"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Step 6.3.1: Heat Capacity Scatter Plot - Cv1 vs Cv2',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Plot all data with simple classification (default)
  python step6_3_1_cv_scatter_plot.py
  
  # Plot with series classification (Pt8Snx, Pt6Snx, O1-O4)
  python step6_3_1_cv_scatter_plot.py --classify series
  
  # Only plot Pt8Snx and Pt6Snx series
  python step6_3_1_cv_scatter_plot.py --only-series pt8snx,pt6snx
  
  # Only plot oxide series O1-O4
  python step6_3_1_cv_scatter_plot.py --only-series o1,o2,o3,o4
  
  # Plot without error bars (cleaner)
  python step6_3_1_cv_scatter_plot.py --no-errorbars
  
  # Exclude specific compositions: (3,5,3) and (3,4,1)
  python step6_3_1_cv_scatter_plot.py --exclude "3,5,3;3,4,1"
        """)
    
    parser.add_argument('--classify', '-c', type=str, default='simple',
                        choices=['simple', 'series', 'detailed'],
                        help='Classification mode: simple (3 types), series (by series), detailed')
    parser.add_argument('--only-series', type=str, default=None,
                        help='Only plot these series (comma-separated): air,pt8snx,pt6snx,other_supported,o1,o2,o3,o4')
    parser.add_argument('--no-air', action='store_true',
                        help='Do not plot gas-phase data')
    parser.add_argument('--no-supported', action='store_true',
                        help='Do not plot supported Pt-Sn (no oxygen) data')
    parser.add_argument('--no-oxide', action='store_true',
                        help='Do not plot supported Pt-Sn-O data')
    parser.add_argument('--merge-oxide', action='store_true',
                        help='Merge O1-O4 into one oxide series (no distinction)')
    parser.add_argument('--no-errorbars', action='store_true',
                        help='Do not show error bars (cleaner plot)')
    parser.add_argument('--no-labels', action='store_true',
                        help='Do not show data point labels (scatter only)')
    parser.add_argument('--interactive', '-i', action='store_true',
                        help='Interactive mode: drag labels with mouse, press S to save')
    parser.add_argument('--no-stroke', action='store_true',
                        help='No white stroke around labels')
    parser.add_argument('--exclude', '-e', type=str, default=None,
                        help='Compositions to exclude, e.g., "3,5,3;3,4,1" or "6,8"')
    parser.add_argument('--fontscale', '-f', type=float, default=1.0,
                        help='Scale factor for all fonts (default 1.0, try 1.5 for larger)')
    parser.add_argument('--markerscale', '-m', type=float, default=1.0,
                        help='Scale factor for marker sizes (default 1.0, try 2 for larger)')
    parser.add_argument('--figsize', type=str, default='8x8',
                        help='Figure size in inches, format WxH (default 8x8)')
    parser.add_argument('-o', '--output', type=str, 
                        default='results/step6_1_clustering/cv1_vs_cv2_scatter.png',
                        help='Output file path')
    
    args = parser.parse_args()
    
    print("=" * 60)
    print("Step 6.3.1: Heat Capacity Scatter Plot - Cv1 vs Cv2")
    print("=" * 60)
    
    # Parse exclude list
    exclude_list = parse_exclude_arg(args.exclude)
    
    # Parse only_series
    only_series = None
    if args.only_series:
        only_series = [s.strip().lower() for s in args.only_series.split(',')]
        print(f"Only plotting series: {only_series}")
    
    # 如果指定了 only_series，自动使用 series 分类模式
    classify_mode = args.classify
    if only_series and classify_mode == 'simple':
        classify_mode = 'series'
        print(f"Auto-switching to series classification mode")
    
    output_path = create_cv_scatter_plot(
        output_path=args.output,
        show_air=not args.no_air,
        show_supported=not args.no_supported,
        show_oxide=not args.no_oxide,
        exclude_list=exclude_list,
        show_errorbars=not args.no_errorbars,
        show_labels=not args.no_labels,
        fontscale=args.fontscale,
        markerscale=args.markerscale,
        classify_mode=classify_mode,
        only_series=only_series,
        merge_oxide=args.merge_oxide,
        interactive=args.interactive,
        no_stroke=args.no_stroke,
        figsize=tuple(map(float, args.figsize.lower().split('x')))
    )
    
    print("\n" + "=" * 60)
    print("Done!")
    print("=" * 60)


if __name__ == '__main__':
    main()
