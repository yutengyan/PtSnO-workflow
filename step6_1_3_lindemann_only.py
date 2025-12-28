#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 6.1.3: Air68 vs Air86 Lindemann指数对比图（简化版）

仅绘制 Lindemann 指数分布图（a1, a2），字体大小按照 6.1.1 标准。
支持按温度排序后排除指定数据点。

用法:
  python step6_1_3_lindemann_only.py
  python step6_1_3_lindemann_only.py --figsize 16x8 --dpi 300
  python step6_1_3_lindemann_only.py --exclude-68 "300K:1,2" "400K:0"
  python step6_1_3_lindemann_only.py --exclude-86 "500K:3,4,5"

作者: AI Assistant
日期: 2025-12-10
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.collections import LineCollection
from pathlib import Path

# 设置高质量论文图样式 - 与 6.1.1 保持一致
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
plt.rcParams['mathtext.fontset'] = 'dejavusans'
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['xtick.major.width'] = 1.5
plt.rcParams['ytick.major.width'] = 1.5

# 字体大小常量 - 与 6.1.1 保持一致
FONT_TICK = 28
FONT_LABEL = 34
FONT_LEGEND = 26

# 分区颜色（科学绘图标准配色）
PARTITION_COLORS = {
    1: '#0173B2',  # 蓝色 (科学绘图标准蓝)
    2: '#DE4343',  # 红色 (科学绘图标准红)
}


def parse_exclude_points(exclude_args):
    """
    解析要排除的点
    
    参数:
        exclude_args: list of str, 例如 ["300K:1,2", "400K:0"]
    
    返回:
        dict: {temp: [indices]} 例如 {300: [1, 2], 400: [0]}
    """
    exclude_dict = {}
    if not exclude_args:
        return exclude_dict
    
    for arg in exclude_args:
        try:
            temp_str, indices_str = arg.split(':')
            temp = int(temp_str.replace('K', ''))
            indices = [int(i) for i in indices_str.split(',')]
            exclude_dict[temp] = indices
        except Exception as e:
            print(f"  ⚠️ 警告: 无法解析排除点 '{arg}': {e}")
    
    return exclude_dict


def filter_data(df, exclude_dict):
    """
    根据排除规则过滤数据
    
    参数:
        df: DataFrame 包含 'temp', 'delta', 'phase_clustered' 列
        exclude_dict: dict {temp: [indices]}
    
    返回:
        DataFrame: 过滤后的数据
    """
    if not exclude_dict:
        return df
    
    # 创建掩码
    mask = np.ones(len(df), dtype=bool)
    
    for temp, indices in exclude_dict.items():
        # 获取该温度的所有数据
        temp_mask = df['temp'] == temp
        temp_indices = np.where(temp_mask)[0]
        
        # 按 Lindemann 指数排序（从大到小）
        temp_df = df[temp_mask].copy()
        temp_df['original_idx'] = temp_indices
        temp_df_sorted = temp_df.sort_values('delta', ascending=False)
        
        # 标记要排除的点
        for idx in indices:
            if idx < len(temp_df_sorted):
                original_idx = temp_df_sorted.iloc[idx]['original_idx']
                mask[original_idx] = False
                lindemann_val = temp_df_sorted.iloc[idx]['delta']
                print(f"    排除: {temp}K 第{idx}个点 (delta={lindemann_val:.4f})")
    
    return df[mask]


def load_clustering_data(csv_path, exclude_dict=None):
    """加载聚类结果数据并过滤"""
    try:
        df = pd.read_csv(csv_path)
        print(f"  ✓ 加载: {os.path.basename(csv_path)}")
        print(f"    原始数据: {len(df)} 条")
        
        if exclude_dict:
            df = filter_data(df, exclude_dict)
            print(f"    过滤后: {len(df)} 条")
        
        return df
    except Exception as e:
        print(f"  ✗ 错误: 无法加载 {csv_path}: {e}")
        return None


def plot_lindemann_single_fig(data, title, figsize=(10, 8), dpi=300, 
                              y_ticks=None, x_ticks=None, x_nticks=None, custom_partitions=None,
                              show_average_line=False, show_error_bars=False, y_lim=None, hide_x_label=False):
    """
    绘制单个系统的 Lindemann 指数图
    
    参数:
        data: 数据 (DataFrame)
        title: 标题
        figsize: 图片尺寸
        dpi: 分辨率
        y_ticks: 自定义Y轴刻度 (list)
        x_ticks: 自定义X轴刻度 (list)
        x_nticks: X轴刻度数量
        custom_partitions: 自定义分区 [(T_min1, T_max1), (T_min2, T_max2)]
        show_average_line: 是否显示平均值连线
        show_error_bars: 是否显示误差棒
        y_lim: Y轴范围 (tuple): (y_min, y_max)
        hide_x_label: 是否隐藏X轴标签和刻度标签
    """
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    
    plot_lindemann_single(ax, data, title, y_ticks=y_ticks, x_ticks=x_ticks,
                         x_nticks=x_nticks, custom_partitions=custom_partitions,
                         show_average_line=show_average_line, show_error_bars=show_error_bars,
                         y_lim=y_lim, hide_x_label=hide_x_label)
    
    # 如果隐藏X轴，使用固定的subplot参数保持占位
    if hide_x_label:
        # 固定边距，确保X轴区域保留
        plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.15)
    else:
        plt.tight_layout()
    
    return fig


def plot_lindemann_single(ax, df, title, y_ticks=None, x_ticks=None, x_nticks=None, custom_partitions=None,
                         show_average_line=False, show_error_bars=False, y_lim=None, hide_x_label=False):
    """绘制单个系统的 Lindemann 指数分布"""
    temps = sorted(df['temp'].unique())
    
    # 如果指定了自定义分区，重新分配分区
    if custom_partitions is not None:
        df = df.copy()
        for i, (T_min, T_max) in enumerate(custom_partitions):
            mask = (df['temp'] >= T_min) & (df['temp'] <= T_max)
            df.loc[mask, 'phase_clustered'] = f'partition{i+1}'
    
    # 按分区绘制散点（只在不显示误差棒时绘制）
    if not show_error_bars:
        for phase in sorted(df['phase_clustered'].unique()):
            df_phase = df[df['phase_clustered'] == phase]
            
            # 提取分区编号（处理 'partition1' 或 1 两种格式）
            if isinstance(phase, str):
                phase_num = int(phase.replace('partition', ''))
            else:
                phase_num = int(phase)
            
            color = PARTITION_COLORS.get(phase_num, '#999999')
            label = f'Partition {phase_num}'
            
            ax.scatter(df_phase['temp'], df_phase['delta'],
                      c=color, s=80, alpha=0.7, edgecolors='black', linewidth=0.5,
                      label=label, zorder=3)
    
    # 绘制误差棒（如果启用）
    if show_error_bars:
        for temp in temps:
            df_temp = df[df['temp'] == temp]
            mean_delta = df_temp['delta'].mean()
            std_delta = df_temp['delta'].std()
            
            # 获取该温度点的分区信息（用于配色）
            phase = df_temp['phase_clustered'].iloc[0]
            if isinstance(phase, str):
                phase_num = int(phase.replace('partition', ''))
            else:
                phase_num = int(phase)
            color = PARTITION_COLORS.get(phase_num, '#999999')
            
            # 绘制误差棒（使用对应分区的颜色）
            ax.errorbar(temp, mean_delta, yerr=std_delta,
                       fmt='none', ecolor=color, elinewidth=2.5,
                       capsize=6, capthick=2.5, alpha=0.5, zorder=1)
    
    # 按分区分段绘制平均值连线（如果启用或显示误差棒时自动启用）
    if show_average_line or show_error_bars:
        all_phase_data = {}  # 存储每个分区的数据用于跨分区连接
        
        for phase in sorted(df['phase_clustered'].unique()):
            df_phase = df[df['phase_clustered'] == phase]
            temps_phase = sorted(df_phase['temp'].unique())
            
            # 计算该分区内每个温度的平均值
            temp_means = []
            temp_values = []
            for temp in temps_phase:
                df_temp = df_phase[df_phase['temp'] == temp]
                mean_delta = df_temp['delta'].mean()
                temp_means.append(mean_delta)
                temp_values.append(temp)
            
            # 提取分区编号
            if isinstance(phase, str):
                phase_num = int(phase.replace('partition', ''))
            else:
                phase_num = int(phase)
            
            # 保存数据用于跨分区连接
            all_phase_data[phase_num] = (temp_values, temp_means)
            
            color = PARTITION_COLORS.get(phase_num, '#999999')
            label = f'Partition {phase_num}' if show_error_bars else None
            
            # 绘制分区内的平均值连线（用分区颜色，实线，实心标记点）
            ax.plot(temp_values, temp_means, 
                    color=color, linestyle='-', linewidth=2.5, 
                    alpha=0.8, zorder=3,
                    marker='o', markersize=14.14, markerfacecolor=color,
                    markeredgecolor=color, markeredgewidth=0,
                    label=label)
        
        # 绘制分区之间的连接线（蓝色到红色渐变）
        if len(all_phase_data) == 2:
            # 获取两个分区的最后一个点和第一个点
            phase1_temps, phase1_means = all_phase_data[1]
            phase2_temps, phase2_means = all_phase_data[2]
            
            # 创建渐变色连接线
            x_start, x_end = phase1_temps[-1], phase2_temps[0]
            y_start, y_end = phase1_means[-1], phase2_means[0]
            
            # 创建插值点以实现渐变效果
            n_points = 100
            x_interp = np.linspace(x_start, x_end, n_points)
            y_interp = np.linspace(y_start, y_end, n_points)
            
            # 创建线段集合
            points = np.array([x_interp, y_interp]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            
            # 创建颜色渐变（从蓝色到红色）
            from matplotlib.colors import LinearSegmentedColormap
            cmap = LinearSegmentedColormap.from_list('blue_red', 
                                                     [PARTITION_COLORS[1], PARTITION_COLORS[2]])
            
            # 绘制渐变线
            lc = LineCollection(segments, cmap=cmap, linewidth=2.5, alpha=0.8, zorder=3)
            lc.set_array(np.linspace(0, 1, len(segments)))
            ax.add_collection(lc)
    
    # 设置标签（不加粗，参考6.1.1.3）
    if not hide_x_label:
        ax.set_xlabel('Temperature (K)', fontsize=FONT_LABEL)
    else:
        # 保持X轴标签占位，但设为空白
        ax.set_xlabel('', fontsize=FONT_LABEL)
    ax.set_ylabel('Lindemann Index', fontsize=FONT_LABEL)
    
    # 设置刻度
    if hide_x_label:
        # 完全隐藏X轴刻度线和标签
        ax.tick_params(axis='x', which='both', bottom=False, top=False, 
                      labelbottom=False, labelsize=FONT_TICK)
        ax.tick_params(axis='y', labelsize=FONT_TICK, width=1.5, length=6)
    else:
        ax.tick_params(axis='both', labelsize=FONT_TICK, width=1.5, length=6)
    
    # 设置X轴范围和刻度
    ax.set_xlim(150, 1150)
    if x_ticks is not None:
        # 使用自定义X轴刻度
        ax.set_xticks(x_ticks)
    elif x_nticks is not None:
        # 自定义X轴刻度数量
        x_ticks_auto = np.linspace(min(temps), max(temps), x_nticks)
        ax.set_xticks(x_ticks_auto)
    else:
        ax.set_xticks(temps)
    
    # 设置Y轴范围和刻度
    if y_lim is not None:
        ax.set_ylim(y_lim)
    
    if y_ticks is not None:
        ax.set_yticks(y_ticks)
    
    # 图例（无边框，放在右下方）
    ax.legend(loc='lower right', fontsize=FONT_LEGEND, frameon=False)


def parse_partitions(partition_str):
    """解析分区字符串，例如 '200-400,500-1100' -> [(200, 400), (500, 1100)]"""
    if not partition_str:
        return None
    try:
        partitions = []
        for part in partition_str.split(','):
            T_min, T_max = map(int, part.split('-'))
            partitions.append((T_min, T_max))
        return partitions
    except Exception as e:
        print(f"  ⚠️ 警告: 无法解析分区 '{partition_str}': {e}")
        return None


def parse_y_ticks(y_ticks_str):
    """解析Y轴刻度，例如 '0,0.1,0.2,0.3' -> [0, 0.1, 0.2, 0.3]"""
    if not y_ticks_str:
        return None
    try:
        return [float(x) for x in y_ticks_str.split(',')]
    except Exception as e:
        print(f"  ⚠️ 警告: 无法解析Y轴刻度 '{y_ticks_str}': {e}")
        return None


def parse_x_ticks(x_ticks_str):
    """解析X轴刻度，例如 '200,400,600,800,1000' -> [200, 400, 600, 800, 1000]"""
    if not x_ticks_str:
        return None
    try:
        return [float(x) for x in x_ticks_str.split(',')]
    except Exception as e:
        print(f"  ⚠️ 警告: 无法解析X轴刻度 '{x_ticks_str}': {e}")
        return None


def main():
    parser = argparse.ArgumentParser(
        description='绘制 Air68 和 Air86 的 Lindemann 指数图（分开绘制）',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # 基本用法（分开绘制）
  python step6_1_3_lindemann_only.py
  
  # 自定义尺寸和刻度（参考6.1.1.3）
  python step6_1_3_lindemann_only.py --figsize 10x8 --x-nticks 5
  
  # 自定义分区和Y轴刻度
  python step6_1_3_lindemann_only.py --partitions-68 200-400,500-1100 --y-ticks 0,0.1,0.2,0.3
  
  # 排除异常点（序号从大到小，0表示最大）
  python step6_1_3_lindemann_only.py --exclude-68 "300K:0,1" "400K:0"
  
  # 完整示例
  python step6_1_3_lindemann_only.py --figsize 10x8 --partitions-68 200-400,500-1100 --y-ticks 0,0.1,0.2,0.3 --x-nticks 5

注意:
  - 点的序号从0开始，按每个温度的Lindemann指数从大到小排序
  - 例如 "300K:0,1" 表示排除300K温度下Lindemann指数最大和第2大的点
  - 两个系统的图会分别保存
"""
    )
    
    parser.add_argument('--figsize', type=str, default='10x8',
                       help='图片尺寸 (默认: 10x8)')
    parser.add_argument('--dpi', type=int, default=300,
                       help='图片分辨率 (默认: 300)')
    parser.add_argument('--partitions-68', type=str, metavar='T1-T2,T3-T4',
                       help='Air68 自定义分区，格式: 200-400,500-1100')
    parser.add_argument('--partitions-86', type=str, metavar='T1-T2,T3-T4',
                       help='Air86 自定义分区，格式: 200-700,800-1100')
    parser.add_argument('--y-ticks', type=str, metavar='Y1,Y2,Y3',
                       help='自定义Y轴刻度，格式: 0,0.1,0.2,0.3')
    parser.add_argument('--x-ticks', type=str, metavar='X1,X2,X3',
                       help='自定义X轴刻度，格式: 200,400,600,800,1000')
    parser.add_argument('--x-nticks', type=int, metavar='N',
                       help='X轴刻度数量，例如: 5')
    parser.add_argument('--exclude-68', nargs='+', metavar='TEMP:INDICES',
                       help='Air68 要排除的点，格式: "300K:0,1" "400K:0"')
    parser.add_argument('--exclude-86', nargs='+', metavar='TEMP:INDICES',
                       help='Air86 要排除的点，格式: "500K:0,1,2"')
    parser.add_argument('--show-average-line', action='store_true',
                       help='显示按分区的平均值连线')
    parser.add_argument('--show-error-bars', action='store_true',
                       help='显示每个温度点的误差棒（标准差）')
    parser.add_argument('--align-y-axis', action='store_true',
                       help='对齐两个系统的Y轴范围，便于直接对比')
    parser.add_argument('--hide-x-label', action='store_true',
                       help='隐藏X轴标签和刻度标签（用于组合图）')
    parser.add_argument('--transparent', action='store_true',
                       help='保存为透明背景图片')
    parser.add_argument('--add-sup86', action='store_true',
                       help='添加 sup86 (负载型 Pt8Sn6) 数据对比')
    parser.add_argument('--partitions-sup86', type=str, metavar='T1-T2,T3-T4',
                       help='sup86 自定义分区，格式: 200-400,500-1100')
    parser.add_argument('--exclude-sup86', nargs='+', metavar='TEMP:INDICES',
                       help='sup86 要排除的点，格式: "500K:0,1,2"')
    
    args = parser.parse_args()
    
    # 解析图片尺寸
    try:
        width, height = map(float, args.figsize.lower().split('x'))
        figsize = (width, height)
    except:
        print(f"⚠️ 警告: 无法解析 figsize '{args.figsize}'，使用默认值 10x8")
        figsize = (10, 8)
    
    # 解析参数
    exclude_68 = parse_exclude_points(args.exclude_68)
    exclude_86 = parse_exclude_points(args.exclude_86)
    exclude_sup86 = parse_exclude_points(args.exclude_sup86) if args.add_sup86 else {}
    partitions_68 = parse_partitions(args.partitions_68)
    partitions_86 = parse_partitions(args.partitions_86)
    partitions_sup86 = parse_partitions(args.partitions_sup86) if args.add_sup86 else None
    y_ticks = parse_y_ticks(args.y_ticks)
    x_ticks = parse_x_ticks(args.x_ticks)
    x_nticks = args.x_nticks
    
    system_names = "Air68、Air86"
    if args.add_sup86:
        system_names += "、sup86"
    print(f"\nStep 6.1.3: {system_names} Lindemann指数图（分开绘制）")
    print(f"图片尺寸: {figsize[0]}x{figsize[1]}, 分辨率: {args.dpi} dpi")
    
    # 数据路径
    base_dir = Path('results/step6_1_clustering')
    csv_68 = base_dir / 'Air68_kmeans_n2_clustered_data.csv'
    csv_86 = base_dir / 'Air86_kmeans_n2_clustered_data.csv'
    csv_sup86 = base_dir / 'Pt8sn6_kmeans_n2_clustered_data.csv' if args.add_sup86 else None
    
    # 检查文件是否存在
    if not csv_68.exists():
        print(f"\n✗ 错误: 找不到 Air68 数据文件: {csv_68}")
        print("  请先运行: python step6_1_clustering_analysis.py")
        return 1
    
    if not csv_86.exists():
        print(f"\n✗ 错误: 找不到 Air86 数据文件: {csv_86}")
        print("  请先运行: python step6_1_clustering_analysis.py")
        return 1
    
    if args.add_sup86 and not csv_sup86.exists():
        print(f"\n✗ 错误: 找不到 sup86 数据文件: {csv_sup86}")
        print("  请先运行: python step6_1_clustering_analysis.py --structure Pt8sn6")
        return 1
    
    output_dir = Path('results/step6_1_partition_cv')
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # ==================== Air68 ====================
    print(f"\n>>> 处理 Air68 (Pt6Sn8)...")
    if exclude_68:
        print(f"  排除点:")
        for temp, indices in exclude_68.items():
            print(f"    {temp}K: {indices}")
    if partitions_68:
        print(f"  自定义分区: {partitions_68}")
    
    data_68 = load_clustering_data(csv_68, exclude_68)
    if data_68 is None:
        return 1
    
    # ==================== Air86 ====================
    print(f"\n>>> 处理 Air86 (Pt8Sn6)...")
    if exclude_86:
        print(f"  排除点:")
        for temp, indices in exclude_86.items():
            print(f"    {temp}K: {indices}")
    if partitions_86:
        print(f"  自定义分区: {partitions_86}")
    
    data_86 = load_clustering_data(csv_86, exclude_86)
    if data_86 is None:
        return 1
    
    # ==================== sup86 (可选) ====================
    data_sup86 = None
    if args.add_sup86:
        print(f"\n>>> 处理 sup86 (负载型 Pt8Sn6)...")
        if exclude_sup86:
            print(f"  排除点:")
            for temp, indices in exclude_sup86.items():
                print(f"    {temp}K: {indices}")
        if partitions_sup86:
            print(f"  自定义分区: {partitions_sup86}")
        
        data_sup86 = load_clustering_data(csv_sup86, exclude_sup86)
        if data_sup86 is None:
            return 1
    
    # 计算统一的Y轴范围（如果启用对齐）
    y_lim = None
    if args.align_y_axis:
        all_data = [data_68, data_86]
        if data_sup86 is not None:
            all_data.append(data_sup86)
        
        y_min = min(d['delta'].min() for d in all_data)
        y_max = max(d['delta'].max() for d in all_data)
        # 添加5%的边距
        y_margin = (y_max - y_min) * 0.05
        y_lim = (y_min - y_margin, y_max + y_margin)
        print(f"\n✓ Y轴对齐: [{y_lim[0]:.4f}, {y_lim[1]:.4f}]")
    
    # ==================== 绘制 Air68 ====================
    print(f"\n>>> 绘制 Air68 图片...")
    fig_68 = plot_lindemann_single_fig(data_68, "Air68 (Pt$_6$Sn$_8$)", 
                                        figsize=figsize, dpi=args.dpi,
                                        y_ticks=y_ticks, x_ticks=x_ticks, x_nticks=x_nticks,
                                        custom_partitions=partitions_68,
                                        show_average_line=args.show_average_line,
                                        show_error_bars=args.show_error_bars,
                                        y_lim=y_lim,
                                        hide_x_label=args.hide_x_label)
    
    output_file_68 = output_dir / 'Air68_lindemann.png'
    # 隐藏X轴时不使用tight，保持原始比例
    if args.hide_x_label:
        fig_68.savefig(output_file_68, dpi=args.dpi, transparent=args.transparent)
    else:
        fig_68.savefig(output_file_68, dpi=args.dpi, bbox_inches='tight', transparent=args.transparent)
    plt.close(fig_68)
    print(f"  ✓ 已保存: {output_file_68}")
    
    # ==================== 绘制 Air86 ====================
    print(f"\n>>> 绘制 Air86 图片...")
    fig_86 = plot_lindemann_single_fig(data_86, "Air86 (Pt$_8$Sn$_6$)", 
                                        figsize=figsize, dpi=args.dpi,
                                        y_ticks=y_ticks, x_ticks=x_ticks, x_nticks=x_nticks,
                                        custom_partitions=partitions_86,
                                        show_average_line=args.show_average_line,
                                        show_error_bars=args.show_error_bars,
                                        y_lim=y_lim,
                                        hide_x_label=args.hide_x_label)
    
    output_file_86 = output_dir / 'Air86_lindemann.png'
    # 隐藏X轴时不使用tight，保持原始比例
    if args.hide_x_label:
        fig_86.savefig(output_file_86, dpi=args.dpi, transparent=args.transparent)
    else:
        fig_86.savefig(output_file_86, dpi=args.dpi, bbox_inches='tight', transparent=args.transparent)
    plt.close(fig_86)
    print(f"  ✓ 已保存: {output_file_86}")
    
    # ==================== 绘制 sup86 (可选) ====================
    if args.add_sup86 and data_sup86 is not None:
        print(f"\n>>> 绘制 sup86 图片...")
        fig_sup86 = plot_lindemann_single_fig(data_sup86, "sup86 (Pt$_8$Sn$_6$/support)", 
                                            figsize=figsize, dpi=args.dpi,
                                            y_ticks=y_ticks, x_ticks=x_ticks, x_nticks=x_nticks,
                                            custom_partitions=partitions_sup86,
                                            show_average_line=args.show_average_line,
                                            show_error_bars=args.show_error_bars,
                                            y_lim=y_lim,
                                            hide_x_label=args.hide_x_label)
        
        output_file_sup86 = output_dir / 'sup86_lindemann.png'
        # 隐藏X轴时不使用tight，保持原始比例
        if args.hide_x_label:
            fig_sup86.savefig(output_file_sup86, dpi=args.dpi, transparent=args.transparent)
        else:
            fig_sup86.savefig(output_file_sup86, dpi=args.dpi, bbox_inches='tight', transparent=args.transparent)
        plt.close(fig_sup86)
        print(f"  ✓ 已保存: {output_file_sup86}")
    
    num_systems = 2 if not args.add_sup86 else 3
    print(f"\n✅ 完成! 已生成 {num_systems} 个独立的图片文件。")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
