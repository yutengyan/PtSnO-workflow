#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
Step 5.1.2: 氧化物系列综合对比分析
================================================================================

功能概述
========
本脚本将 O1、O2、O3、O4 系列绑制在一起，提供多种可视化方式：

1. 热图子图模式 (--mode heatmap)
   - 2x2 子图排列，每个子图显示一个氧含量系列的 Lindemann 热图
   
2. 熔点对比模式 (--mode tm)
   - 单图叠加所有系列的 Tm vs 组分曲线，不同颜色区分
   
3. 综合模式 (--mode both) - 默认
   - 同时生成热图子图和熔点对比图
   
4. 合并热图模式 (--mode combined)
   - 将所有系列合并到一张热图中，支持按原子数或氧系列排序

数据来源
========
1. Lindemann 原始数据: results/step6_0_multi_system/step6_0_all_systems_data.csv
2. 熔点汇总数据: results/step5_1_melting_point/melting_point_summary.csv

输出文件
========
results/step5_1_2_oxide_series/
├── oxide_heatmap_grid.png       # 2x2 热图子图
├── oxide_combined_heatmap.png   # 合并热图（combined模式）
├── oxide_tm_comparison.png      # Tm 对比曲线图
└── oxide_tm_vs_metal.png        # Tm vs 金属原子数（分氧含量着色）

命令行参数
==========
必选参数：无

可选参数：
  --mode, -m       : 绘图模式，可选值：
                     - heatmap  : 2x2热图子图，每个子图一个氧系列
                     - tm       : Tm对比曲线图
                     - both     : 同时生成heatmap和tm（默认）
                     - combined : 合并热图，所有系列在一张图中
                     
  --sort, -s       : 排序方式，仅对combined模式有效，可选值：
                     - atoms    : 按原子总数排序（默认）
                     - oxygen   : 按O1→O2→O3→O4氧系列分组，组内按原子数排序
                     
  --tm-method, -t  : 熔点判据选择，可选值：
                     - threshold  : δ=0.1阈值法（默认），黑色等高线
                     - transition : dδ/dT最大值跃变点法，绿色等高线
                     - both       : 同时显示两种判据
                     
  --exclude, -e    : 排除指定结构，支持多个，支持简写格式：
                     - 794       : 表示 Pt7Sn9O4
                     - 68o4      : 表示 Pt6Sn8O4
                     - Pt7Sn9O4  : 完整名称
                     
  --fontscale, -f  : 字体缩放比例，默认1.0
  
  --no-title       : 不显示图片标题
  
  --no-show        : 只保存图片，不弹出交互式窗口
  
  --separator      : 氧系列分隔线颜色（仅--sort oxygen时有效）
                     - white : 白色分隔线（默认）
                     - black : 黑色分隔线
                     - none  : 不绘制分隔线

使用示例
========
# 1. 基础用法 - 综合模式（默认生成 heatmap + tm 两张图）
python step5_1_2_oxide_series_plot.py --no-show

# 2. 只生成2x2热图子图
python step5_1_2_oxide_series_plot.py -m heatmap --no-title --no-show

# 3. 只生成 Tm 对比图，放大字体
python step5_1_2_oxide_series_plot.py -m tm --no-title --no-show -f 1.2

# 4. 合并热图模式 - 按原子数排序（默认）
python step5_1_2_oxide_series_plot.py --mode combined --no-show

# 5. 合并热图模式 - 按O1→O2→O3→O4系列分组排序（白线分隔）
python step5_1_2_oxide_series_plot.py --mode combined --sort oxygen --no-show

# 6. 按氧系列排序，使用黑色分隔线
python step5_1_2_oxide_series_plot.py --mode combined --sort oxygen --separator black --no-show

# 7. 排除特定结构（如排除 Pt7Sn9O4）
python step5_1_2_oxide_series_plot.py --mode combined -s oxygen -e 794 --no-show

# 8. 同时排除多个结构
python step5_1_2_oxide_series_plot.py --mode combined -e 794 68o4 --no-show

# 9. 使用跃变点法显示熔点等高线
python step5_1_2_oxide_series_plot.py --mode combined --tm-method transition --no-show

# 10. 同时显示阈值法和跃变点法的熔点等高线
python step5_1_2_oxide_series_plot.py --mode combined --tm-method both --no-show

# 11. 使用热容分区边界线
python step5_1_2_oxide_series_plot.py --mode combined --tm-method cv-partition --no-show

# 12. 完整参数示例：按氧系列排序，排除794，使用热容分区，不显示标题
python step5_1_2_oxide_series_plot.py -m combined -s oxygen -e 794 -t cv-partition --no-title --no-show

熔点判据说明
============
1. 阈值法 (threshold)：
   - 定义：Lindemann指数 δ 首次达到0.1的温度
   - 物理含义：表面熔化开始的温度
   - 可视化：黑色虚线等高线
   
2. 跃变点法 (transition)：
   - 定义：dδ/dT（Lindemann指数随温度变化率）达到最大值的温度
   - 物理含义：熔化过程最剧烈的温度，接近体相熔化
   - 可视化：绿色虚线等高线
   
3. 热容分区法 (cv-partition)：
   - 定义：step6_1聚类分析得到的热容分区边界温度
   - 物理含义：基于能量-温度关系的相变温度
   - 可视化：品红色虚线
   - 数据来源：results/step6_1_clustering/*_kmeans_n2_clustered_data.csv

排序方式说明
============
1. atoms 模式：
   - 所有结构按总原子数(Pt+Sn+O)从小到大排列
   - 适合查看尺寸效应
   
2. oxygen 模式：
   - 首先按氧原子数分组：O1 → O2 → O3 → O4
   - 组内按总原子数排序
   - 自动绘制白线/黑线分隔不同氧系列
   - 适合对比不同氧含量系列的规律

================================================================================
"""
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

BASE_DIR = Path(__file__).parent
DATA_FILE = BASE_DIR / 'results' / 'step6_0_multi_system' / 'step6_0_all_systems_data.csv'
MP_SUMMARY = BASE_DIR / 'results' / 'step5_1_melting_point' / 'melting_point_summary.csv'
OUTPUT_DIR = BASE_DIR / 'results' / 'step5_1_2_oxide_series'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# 英文字体设置
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.unicode_minus'] = False

# 氧化物系列配置
OXIDE_SERIES = {
    'O1': {'O': 1, 'color': '#1f77b4', 'marker': 'o'},  # 蓝色
    'O2': {'O': 2, 'color': '#ff7f0e', 'marker': 's'},  # 橙色
    'O3': {'O': 3, 'color': '#2ca02c', 'marker': '^'},  # 绿色
    'O4': {'O': 4, 'color': '#d62728', 'marker': 'D'},  # 红色
}


# 特殊结构名映射（如 Cv -> Pt6Sn8O4）
STRUCTURE_NAME_MAP = {
    'cv': 'Pt6Sn8O4',
}


def format_structure_label(name):
    """
    将结构名格式化为 (x,y,z) 形式，代表 Pt_x Sn_y O_z
    例如：Sn5o4pt2 -> (2,5,4), Cv -> (6,8,4)
    """
    import re
    
    # 检查特殊名称映射
    name_lower = name.lower()
    if name_lower in STRUCTURE_NAME_MAP:
        name = STRUCTURE_NAME_MAP[name_lower]
    
    # 提取各元素的数量
    elements = {}
    pattern = r'([A-Za-z]+)(\d+)'
    for match in re.finditer(pattern, name, re.IGNORECASE):
        elem = match.group(1).lower()
        num = int(match.group(2))
        elements[elem] = num
    
    # 返回 (Pt, Sn, O) 格式
    pt = elements.get('pt', 0)
    sn = elements.get('sn', 0)
    o = elements.get('o', 0)
    
    if pt > 0 or sn > 0 or o > 0:
        return f'({pt},{sn},{o})'
    return name


def format_structure_label_full(name):
    """
    将结构名格式化为 LaTeX 下标形式，统一顺序为 Pt_x Sn_y O_z
    例如：Sn5o4pt2 -> Pt$_2$Sn$_5$O$_4$, Cv -> Pt$_6$Sn$_8$O$_4$
    """
    import re
    
    # 检查特殊名称映射
    name_lower = name.lower()
    if name_lower in STRUCTURE_NAME_MAP:
        name = STRUCTURE_NAME_MAP[name_lower]
    
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


def parse_args():
    p = argparse.ArgumentParser(description='Plot oxide series (O1-O4) comparison')
    p.add_argument('--mode', '-m', default='both', choices=['heatmap', 'tm', 'both', 'combined'],
                   help='绘图模式: heatmap=2x2热图, tm=Tm对比, both=两者, combined=合并热图')
    p.add_argument('--sort', '-s', default='atoms', choices=['atoms', 'oxygen'],
                   help='排序方式: atoms=按原子总数, oxygen=按O1-O4氧系列分组')
    p.add_argument('--tm-method', '-t', default='threshold', choices=['threshold', 'transition', 'both', 'cv-partition'],
                   help='熔点判据: threshold=δ=0.1阈值法, transition=dδ/dT跃变点法, both=两者都显示, cv-partition=热容分区边界')
    p.add_argument('--exclude', '-e', nargs='+', default=[], 
                   help='排除指定结构，支持简写如 794 表示 Pt7Sn9O4')
    p.add_argument('--fontscale', '-f', type=float, default=1.0, help='字体缩放比例')
    p.add_argument('--no-title', action='store_true', help='不显示图片标题')
    p.add_argument('--no-show', action='store_true', help='只保存图片，不弹出交互式窗口')
    p.add_argument('--separator', default='white', choices=['white', 'black', 'none'],
                   help='氧系列分隔线颜色 (仅--sort oxygen时有效): white=白线, black=黑线, none=不绘制')
    return p.parse_args()


def expand_structure_name(name):
    """
    扩展简写结构名，如 794 -> Pt7Sn9O4
    支持格式：794, 68o4, pt7sn9o4 等
    """
    name = name.strip()
    # 如果是纯数字如 794，解析为 Pt7Sn9O4
    if name.isdigit() and len(name) == 3:
        return f"Pt{name[0]}Sn{name[1]}O{name[2]}"
    if name.isdigit() and len(name) == 2:
        return f"Pt{name[0]}Sn{name[1]}"
    # 否则返回原名（会做大小写不敏感匹配）
    return name


def load_data():
    """加载数据"""
    df_all = pd.read_csv(DATA_FILE)
    df_mp = pd.read_csv(MP_SUMMARY)
    df_mp['metal_total'] = df_mp['Pt'] + df_mp['Sn']
    return df_all, df_mp


def load_cv_partition_boundary(structure):
    """
    加载热容分区边界温度
    返回分区边界温度（partition1和partition2之间的边界）
    """
    import glob
    clustering_dir = BASE_DIR / 'results' / 'step6_1_clustering'
    
    # 尝试匹配结构名（大小写不敏感）
    pattern = str(clustering_dir / '*_kmeans_n2_clustered_data.csv')
    files = glob.glob(pattern)
    
    for f in files:
        fname = Path(f).name
        struct_in_file = fname.replace('_kmeans_n2_clustered_data.csv', '')
        if struct_in_file.lower() == structure.lower():
            try:
                df = pd.read_csv(f)
                # 找出每个分区的温度范围
                p1_temps = df[df['phase_clustered'] == 'partition1']['temp'].unique()
                p2_temps = df[df['phase_clustered'] == 'partition2']['temp'].unique()
                
                if len(p1_temps) > 0 and len(p2_temps) > 0:
                    # 边界是 partition1 最高温度和 partition2 最低温度的中间
                    boundary = (max(p1_temps) + min(p2_temps)) / 2
                    return boundary
            except Exception as e:
                pass
    return None


def plot_heatmap_grid(df_all, df_mp, args):
    """绑制 2x2 热图子图"""
    fs = args.fontscale
    
    # 统一的温度网格 (100K 间隔)
    unified_temps = np.arange(200, 1150, 100)
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    axes = axes.flatten()
    
    for idx, (series_name, config) in enumerate(OXIDE_SERIES.items()):
        ax = axes[idx]
        o_val = config['O']
        
        # 筛选该氧含量的结构
        df_series = df_mp[(df_mp['O'] == o_val) & (df_mp['type'] == 'oxide')].copy()
        if df_series.empty:
            ax.text(0.5, 0.5, f'No data for {series_name}', ha='center', va='center', 
                   transform=ax.transAxes, fontsize=14)
            ax.set_title(f'Pt$_x$Sn$_y$O$_{{{o_val}}}$', fontsize=int(14*fs), fontweight='bold')
            continue
        
        # 排序
        df_series['Pt_ratio'] = df_series['Pt'] / (df_series['Pt'] + df_series['Sn'] + 0.001)
        df_series = df_series.sort_values(['metal_total', 'Pt_ratio'], ascending=[True, False])
        structures = list(df_series['structure'].unique())
        
        # 提取 Lindemann 数据
        df_sel = df_all[df_all['structure'].isin(structures)].copy()
        if df_sel.empty:
            continue
            
        df_mean = df_sel.groupby(['structure', 'temp']).agg(delta_mean=('delta', 'mean')).reset_index()
        
        # 重采样到统一温度网格
        resampled_data = []
        for struct in structures:
            df_struct = df_mean[df_mean['structure'] == struct].sort_values('temp')
            if len(df_struct) >= 2:
                # 线性插值到统一温度网格
                from scipy.interpolate import interp1d
                f = interp1d(df_struct['temp'], df_struct['delta_mean'], 
                           kind='linear', bounds_error=False, fill_value=np.nan)
                for t in unified_temps:
                    resampled_data.append({
                        'structure': struct,
                        'temp': t,
                        'delta_mean': f(t)
                    })
        
        df_resampled = pd.DataFrame(resampled_data)
        
        # 使用结构索引作为 x 轴
        struct_to_idx = {s: i for i, s in enumerate(structures)}
        df_resampled['x_val'] = df_resampled['structure'].map(struct_to_idx)
        pivot_table = df_resampled.pivot(index='temp', columns='x_val', values='delta_mean')
        pivot_table = pivot_table.sort_index(ascending=True)
        
        temps = pivot_table.index.values
        cols = pivot_table.columns.values
        data = pivot_table.values.astype(float)  # 确保是 float 类型
        x_labels = [format_structure_label(s) for s in structures]
        
        # 绘制热图
        im = ax.imshow(data, aspect='auto', origin='lower', cmap=plt.cm.RdYlBu_r,
                      interpolation='bilinear', vmin=0, vmax=0.3)
        
        # 坐标轴
        ax.set_xticks(np.arange(len(cols)))
        ax.set_xticklabels(x_labels, fontsize=int(9*fs), rotation=45, ha='right')
        ax.set_yticks(np.arange(0, len(temps), max(1, len(temps)//6)))
        ax.set_yticklabels([f'{int(temps[i])}K' for i in range(0, len(temps), max(1, len(temps)//6))], 
                          fontsize=int(9*fs))
        
        ax.set_xlabel('Structure', fontsize=int(11*fs), fontweight='bold')
        ax.set_ylabel('Temperature (K)', fontsize=int(11*fs), fontweight='bold')
        ax.set_title(f'Pt$_x$Sn$_y$O$_{{{o_val}}}$', fontsize=int(14*fs), fontweight='bold')
        
        # 等值线
        X, Y = np.meshgrid(np.arange(len(cols)), np.arange(len(temps)))
        ax.contour(X, Y, data, levels=[0.1], colors=['black'], linewidths=1.5, linestyles='--')
    
    # 添加共享 colorbar
    fig.subplots_adjust(right=0.88, hspace=0.35, wspace=0.25)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label('Lindemann Index δ', fontsize=int(12*fs))
    cbar.ax.tick_params(labelsize=int(10*fs))
    cbar.ax.axhline(0.1, color='black', linestyle='--', linewidth=2)
    
    if not args.no_title:
        fig.suptitle('Oxide Series: Lindemann Index Heatmaps', fontsize=int(16*fs), fontweight='bold', y=0.98)
    
    out = OUTPUT_DIR / 'oxide_heatmap_grid.png'
    fig.savefig(out, dpi=200, bbox_inches='tight')
    print(f"Saved: {out}")
    if not args.no_show:
        plt.show()
    plt.close(fig)


def plot_combined_heatmap(df_all, df_mp, args):
    """绑制 O1-O4 合并热图（横向排列，按总原子数排序）"""
    fs = args.fontscale
    from scipy.interpolate import interp1d
    
    # 处理排除列表：扩展简写并转为小写
    exclude_list = [expand_structure_name(e).lower() for e in args.exclude]
    
    # 统一的温度网格 (100K 间隔)
    unified_temps = np.arange(200, 1150, 100)
    
    # 收集所有氧化物结构
    all_structures = []
    for series_name, config in OXIDE_SERIES.items():
        o_val = config['O']
        df_series = df_mp[(df_mp['O'] == o_val) & (df_mp['type'] == 'oxide')].copy()
        # 排除指定结构（大小写不敏感）
        df_series = df_series[~df_series['structure'].str.lower().isin(exclude_list)]
        if not df_series.empty:
            df_series['series'] = series_name
            df_series['color'] = config['color']
            all_structures.append(df_series)
    
    if exclude_list:
        print(f"  Excluding structures: {args.exclude} -> {exclude_list}")
    
    if not all_structures:
        print("[WARN] No oxide structures found")
        return
    
    df_all_struct = pd.concat(all_structures, ignore_index=True)
    
    # 计算辅助排序列
    df_all_struct['total_atoms'] = df_all_struct['Pt'] + df_all_struct['Sn'] + df_all_struct['O']
    df_all_struct['Pt_ratio'] = df_all_struct['Pt'] / (df_all_struct['Pt'] + df_all_struct['Sn'] + 0.001)
    
    # 根据排序模式选择排序方式
    if args.sort == 'oxygen':
        # 按氧系列分组 (O1 -> O2 -> O3 -> O4)，组内按原子数和Pt比例排序
        df_all_struct = df_all_struct.sort_values(['O', 'total_atoms', 'Pt_ratio'], ascending=[True, True, False])
        sort_desc = "sorted by O series (O1→O4)"
    else:
        # 默认：按总原子数排序，同总原子数按 O 数量排序，同 O 数量按 Pt 比例排序
        df_all_struct = df_all_struct.sort_values(['total_atoms', 'O', 'Pt_ratio'], ascending=[True, True, False])
        sort_desc = "sorted by total atoms"
    
    structures = list(df_all_struct['structure'].unique())
    struct_info = {row['structure']: {'O': row['O'], 'series': row['series'], 'color': row['color']} 
                   for _, row in df_all_struct.iterrows()}
    
    print(f"  Combined heatmap: {len(structures)} structures, {sort_desc}")
    
    # 提取 Lindemann 数据并重采样
    df_sel = df_all[df_all['structure'].isin(structures)].copy()
    df_mean = df_sel.groupby(['structure', 'temp']).agg(delta_mean=('delta', 'mean')).reset_index()
    
    resampled_data = []
    for struct in structures:
        df_struct = df_mean[df_mean['structure'] == struct].sort_values('temp')
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
    
    # 创建图（增大尺寸以适应更大的字体）
    fig_width = max(24, len(structures) * 1.2)
    fig, ax = plt.subplots(figsize=(fig_width, 10))
    
    # 绘制热图
    im = ax.imshow(data, aspect='auto', origin='lower', cmap=plt.cm.RdYlBu_r,
                  interpolation='bilinear', vmin=0, vmax=0.3)
    
    # x 轴标签：(x,y,z) 格式
    x_labels = [format_structure_label(s) for s in structures]
    ax.set_xticks(np.arange(len(cols)))
    ax.set_xticklabels(x_labels, fontsize=int(28*fs), rotation=45, ha='right')
    
    # y 轴（去掉K单位）
    ax.set_yticks(np.arange(len(temps)))
    ax.set_yticklabels([f'{int(t)}' for t in temps], fontsize=int(28*fs))
    
    ax.set_xlabel('Pt$_x$Sn$_y$O$_z$ (x,y,z)', fontsize=int(34*fs), fontweight='bold')
    ax.set_ylabel('Temperature', fontsize=int(34*fs), fontweight='bold')
    
    # 等值线 - 根据判据选择
    X, Y = np.meshgrid(np.arange(len(cols)), np.arange(len(temps)))
    
    if args.tm_method in ['threshold', 'both']:
        # 阈值法: δ=0.1 等值线 - 对每个结构取最低温度点
        threshold_points = []
        for col_idx in range(len(cols)):
            col_data = data[:, col_idx]
            if np.isnan(col_data).all():
                threshold_points.append(np.nan)
                continue
            # 找到第一个（最低温度）达到0.1的点
            above_threshold = np.where(col_data >= 0.1)[0]
            if len(above_threshold) > 0:
                first_idx = above_threshold[0]
                # 线性插值找到精确位置
                if first_idx > 0 and col_data[first_idx-1] < 0.1:
                    # 在 first_idx-1 和 first_idx 之间插值
                    y0, y1 = col_data[first_idx-1], col_data[first_idx]
                    t_interp = first_idx - 1 + (0.1 - y0) / (y1 - y0)
                    threshold_points.append(t_interp)
                else:
                    threshold_points.append(first_idx)
            else:
                threshold_points.append(np.nan)
        
        # 绘制阈值点连线（不跨越氧系列边界）
        # 找出氧系列边界
        series_boundaries = [0]  # 每个系列的起始索引
        current_o = struct_info[structures[0]]['O']
        for i, struct in enumerate(structures):
            o_val = struct_info[struct]['O']
            if o_val != current_o:
                series_boundaries.append(i)
                current_o = o_val
        series_boundaries.append(len(structures))
        
        # 分段绘制（按氧系列）
        for seg_idx in range(len(series_boundaries) - 1):
            start, end = series_boundaries[seg_idx], series_boundaries[seg_idx + 1]
            seg_x = list(range(start, end))
            seg_y = [threshold_points[i] for i in seg_x]
            valid_pairs = [(x, y) for x, y in zip(seg_x, seg_y) if not np.isnan(y)]
            if len(valid_pairs) > 1:
                vx, vy = zip(*valid_pairs)
                ax.plot(vx, vy, 'w-', linewidth=2.5, alpha=0.8)
                ax.plot(vx, vy, 'k--', linewidth=2)
    
    # 辅助函数：获取氧系列边界（用于分段绘制）
    def get_series_boundaries():
        boundaries = [0]
        current_o = struct_info[structures[0]]['O']
        for i, struct in enumerate(structures):
            o_val = struct_info[struct]['O']
            if o_val != current_o:
                boundaries.append(i)
                current_o = o_val
        boundaries.append(len(structures))
        return boundaries
    
    if args.tm_method in ['transition', 'both']:
        # 跃变点法: 绘制 dδ/dT 最大值对应的线
        # 计算每个结构的 dδ/dT，找最大值点
        from scipy.ndimage import gaussian_filter1d
        transition_points = []
        for col_idx in range(len(cols)):
            col_data = data[:, col_idx]
            if np.isnan(col_data).all():
                transition_points.append(np.nan)
                continue
            # 平滑后求导
            smooth_data = gaussian_filter1d(col_data, sigma=1)
            dDelta = np.diff(smooth_data)
            if len(dDelta) > 0:
                max_idx = np.argmax(dDelta)
                transition_points.append(max_idx + 0.5)  # 中点位置
            else:
                transition_points.append(np.nan)
        
        # 分段绘制跃变点连线（按氧系列）
        series_boundaries = get_series_boundaries()
        for seg_idx in range(len(series_boundaries) - 1):
            start, end = series_boundaries[seg_idx], series_boundaries[seg_idx + 1]
            seg_x = list(range(start, end))
            seg_y = [transition_points[i] for i in seg_x]
            valid_pairs = [(x, y) for x, y in zip(seg_x, seg_y) if not np.isnan(y)]
            if len(valid_pairs) > 1:
                vx, vy = zip(*valid_pairs)
                ax.plot(vx, vy, 'w-', linewidth=2.5, alpha=0.8)
                ax.plot(vx, vy, 'g--', linewidth=2)
    
    # 热容分区边界线 (cv-partition)
    if args.tm_method == 'cv-partition':
        cv_boundary_points = []
        for col_idx, struct in enumerate(structures):
            boundary_temp = load_cv_partition_boundary(struct)
            if boundary_temp is not None:
                # 将温度转换为 y 轴索引
                y_idx = (boundary_temp - temps[0]) / (temps[1] - temps[0])
                cv_boundary_points.append(y_idx)
            else:
                cv_boundary_points.append(np.nan)
        
        # 分段绘制热容分区边界线（按氧系列）
        series_boundaries = get_series_boundaries()
        for seg_idx in range(len(series_boundaries) - 1):
            start, end = series_boundaries[seg_idx], series_boundaries[seg_idx + 1]
            seg_x = list(range(start, end))
            seg_y = [cv_boundary_points[i] for i in seg_x]
            valid_pairs = [(x, y) for x, y in zip(seg_x, seg_y) if not np.isnan(y)]
            if len(valid_pairs) > 1:
                vx, vy = zip(*valid_pairs)
                ax.plot(vx, vy, 'w-', linewidth=2.5, alpha=0.8)
                ax.plot(vx, vy, 'm--', linewidth=2)
    
    # 绘制氧系列分隔线（仅当 --sort oxygen 时）
    if args.sort == 'oxygen' and args.separator != 'none':
        sep_color = args.separator
        # 找到每个氧系列的边界位置
        current_o = None
        sep_positions = []
        for i, struct in enumerate(structures):
            o_val = struct_info[struct]['O']
            if current_o is not None and o_val != current_o:
                sep_positions.append(i - 0.5)  # 在两个系列之间
            current_o = o_val
        
        # 绘制分隔线
        for pos in sep_positions:
            ax.axvline(x=pos, color=sep_color, linewidth=2, linestyle='-', alpha=0.8)
    
    # Colorbar
    cbar = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.02)
    cbar.set_label('Lindemann Index δ', fontsize=int(34*fs))
    cbar.ax.tick_params(labelsize=int(28*fs))
    if args.tm_method in ['threshold', 'both']:
        cbar.ax.axhline(0.1, color='black', linestyle='--', linewidth=2)
    
    # 图例（简化版 - 只显示判据类型，不显示O1-O4颜色，无框框）
    from matplotlib.lines import Line2D
    legend_elements = []
    # 添加判据图例
    if args.tm_method in ['threshold', 'both']:
        legend_elements.append(Line2D([0], [0], color='black', linestyle='--', lw=2, 
                                      label='δ=0.1'))
    if args.tm_method in ['transition', 'both']:
        legend_elements.append(Line2D([0], [0], color='green', linestyle='--', lw=2, 
                                      label='max dδ/dT'))
    if args.tm_method == 'cv-partition':
        legend_elements.append(Line2D([0], [0], color='magenta', linestyle='--', lw=2, 
                                      label='Cv Partition'))
    # 只有当有图例元素时才显示
    if legend_elements:
        ax.legend(handles=legend_elements, loc='upper left', fontsize=int(26*fs), 
                 frameon=False)  # frameon=False 去掉框框
    
    # 标题显示判据类型
    method_label = {'threshold': '(δ=0.1)', 'transition': '(max dδ/dT)', 'both': '(Both Methods)', 
                    'cv-partition': '(Cv Partition)'}
    if not args.no_title:
        ax.set_title(f'Oxide Series: Lindemann Index (O1-O4 Combined) {method_label[args.tm_method]}', 
                    fontsize=int(14*fs), fontweight='bold')
    
    plt.tight_layout()
    # 输出文件名包含判据类型
    suffix = {'threshold': '', 'transition': '_transition', 'both': '_both', 'cv-partition': '_cv_partition'}
    out = OUTPUT_DIR / f'oxide_heatmap_combined{suffix[args.tm_method]}.png'
    fig.savefig(out, dpi=200, bbox_inches='tight')
    print(f"Saved: {out}")
    if not args.no_show:
        plt.show()
    plt.close(fig)


def plot_tm_comparison(df_mp, args):
    """绘制 Tm 对比曲线图 - 按结构排列，支持双判据"""
    fs = args.fontscale
    
    # cv-partition 模式不绘制 Tm 对比图
    if args.tm_method == 'cv-partition':
        print("  [INFO] cv-partition mode: skipping Tm comparison plot")
        return
    
    # 确定使用哪个 Tm 列
    tm_columns = []
    if args.tm_method in ['threshold', 'both']:
        tm_columns.append(('Tm_lindemann', 'Threshold (δ=0.1)', 'o', 1.0))
    if args.tm_method in ['transition', 'both']:
        tm_columns.append(('Tm_transition', 'Transition (max dδ/dT)', '^', 0.6))
    
    fig, ax = plt.subplots(figsize=(14, 6))
    
    for tm_col, tm_label, marker_style, alpha in tm_columns:
        all_data = []
        
        for series_name, config in OXIDE_SERIES.items():
            o_val = config['O']
            df_series = df_mp[(df_mp['O'] == o_val) & (df_mp['type'] == 'oxide')].copy()
            df_series = df_series.dropna(subset=[tm_col])
            
            if df_series.empty:
                continue
            
            # 排序
            df_series['Pt_ratio'] = df_series['Pt'] / (df_series['Pt'] + df_series['Sn'] + 0.001)
            df_series = df_series.sort_values(['metal_total', 'Pt_ratio'], ascending=[True, False])
            
            for _, row in df_series.iterrows():
                all_data.append({
                    'structure': row['structure'],
                    'series': series_name,
                    'O': o_val,
                    'metal_total': row['metal_total'],
                    'Pt_ratio': row['Pt_ratio'],
                    'Tm': row[tm_col],
                    'color': config['color'],
                    'marker': config['marker'],
                })
        
        if not all_data:
            continue
            
        df_plot = pd.DataFrame(all_data)
        
        # 按 metal_total 和 Pt_ratio 排序
        df_plot = df_plot.sort_values(['metal_total', 'Pt_ratio'], ascending=[True, False])
        df_plot['x_idx'] = range(len(df_plot))
        
        # 分系列绘制
        for series_name, config in OXIDE_SERIES.items():
            df_s = df_plot[df_plot['series'] == series_name]
            if not df_s.empty:
                label = f'{series_name} ({tm_label})' if len(tm_columns) > 1 else f'Pt$_x$Sn$_y$O$_{{{config["O"]}}}$'
                ax.scatter(df_s['x_idx'], df_s['Tm'], 
                          c=config['color'], marker=marker_style, s=100, alpha=alpha,
                          label=label, zorder=3, edgecolors='black' if alpha < 1 else 'none')
    
    # x 轴标签（使用最后一个 df_plot）
    if 'df_plot' in dir() and not df_plot.empty:
        ax.set_xticks(df_plot['x_idx'])
        ax.set_xticklabels([format_structure_label(s) for s in df_plot['structure']], 
                           rotation=45, ha='right', fontsize=int(9*fs))
    
    ax.set_xlabel('Structure', fontsize=int(12*fs), fontweight='bold')
    ax.set_ylabel('Melting Point $T_m$ (K)', fontsize=int(12*fs), fontweight='bold')
    ax.tick_params(axis='y', labelsize=int(11*fs))
    
    method_label = {'threshold': '(Threshold)', 'transition': '(Transition)', 'both': '(Both Methods)',
                    'cv-partition': '(Cv Partition)'}
    if not args.no_title:
        ax.set_title(f'Oxide Series: Melting Point Comparison {method_label[args.tm_method]}', 
                    fontsize=int(14*fs), fontweight='bold')
    
    ax.legend(fontsize=int(9*fs), loc='best', ncol=2 if len(tm_columns) > 1 else 1)
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    suffix = {'threshold': '', 'transition': '_transition', 'both': '_both', 'cv-partition': '_cv_partition'}
    out = OUTPUT_DIR / f'oxide_tm_comparison{suffix[args.tm_method]}.png'
    fig.savefig(out, dpi=300, bbox_inches='tight')
    print(f"Saved: {out}")
    if not args.no_show:
        plt.show()
    plt.close(fig)


def plot_tm_vs_metal(df_mp, args):
    """绘制 Tm vs 金属原子数，按氧含量着色"""
    fs = args.fontscale
    
    fig, ax = plt.subplots(figsize=(9, 6))
    
    for series_name, config in OXIDE_SERIES.items():
        o_val = config['O']
        df_series = df_mp[(df_mp['O'] == o_val) & (df_mp['type'] == 'oxide')].copy()
        df_series = df_series.dropna(subset=['Tm_lindemann'])
        
        if df_series.empty:
            continue
        
        # 绘制散点
        ax.scatter(df_series['metal_total'], df_series['Tm_lindemann'],
                  c=config['color'], marker=config['marker'], s=120,
                  label=f'Pt$_x$Sn$_y$O$_{{{o_val}}}$', zorder=3, alpha=0.8)
        
        # 添加趋势线（如果点数足够）
        if len(df_series) >= 3:
            z = np.polyfit(df_series['metal_total'], df_series['Tm_lindemann'], 1)
            p = np.poly1d(z)
            x_line = np.linspace(df_series['metal_total'].min(), df_series['metal_total'].max(), 50)
            ax.plot(x_line, p(x_line), '--', color=config['color'], alpha=0.5, linewidth=1.5)
    
    ax.set_xlabel('Number of Metal Atoms (Pt+Sn)', fontsize=int(12*fs), fontweight='bold')
    ax.set_ylabel('Melting Point $T_m$ (K)', fontsize=int(12*fs), fontweight='bold')
    ax.tick_params(axis='both', labelsize=int(11*fs))
    
    if not args.no_title:
        ax.set_title('Oxide Series: $T_m$ vs Metal Content', fontsize=int(14*fs), fontweight='bold')
    
    ax.legend(fontsize=int(11*fs), loc='best')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    out = OUTPUT_DIR / 'oxide_tm_vs_metal.png'
    fig.savefig(out, dpi=300, bbox_inches='tight')
    print(f"Saved: {out}")
    if not args.no_show:
        plt.show()
    plt.close(fig)


def main():
    args = parse_args()
    
    if not DATA_FILE.exists():
        print(f"[ERROR] Data file not found: {DATA_FILE}")
        return
    if not MP_SUMMARY.exists():
        print(f"[ERROR] Melting point summary not found: {MP_SUMMARY}")
        return
    
    print(f"[INFO] Loading data...")
    df_all, df_mp = load_data()
    
    # 统计各系列结构数
    for series_name, config in OXIDE_SERIES.items():
        count = len(df_mp[(df_mp['O'] == config['O']) & (df_mp['type'] == 'oxide')])
        print(f"  {series_name}: {count} structures")
    
    if args.mode in ['heatmap', 'both']:
        print(f"\n[INFO] Plotting heatmap grid (2x2)...")
        plot_heatmap_grid(df_all, df_mp, args)
    
    if args.mode in ['combined']:
        print(f"\n[INFO] Plotting combined heatmap (O1-O4 in one figure)...")
        plot_combined_heatmap(df_all, df_mp, args)
    
    if args.mode in ['tm', 'both', 'combined']:
        print(f"\n[INFO] Plotting Tm comparison...")
        plot_tm_comparison(df_mp, args)
        plot_tm_vs_metal(df_mp, args)
    
    print("\n[DONE] All plots saved.")


if __name__ == '__main__':
    main()
