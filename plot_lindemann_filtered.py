#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
单独绘制 Lindemann Index vs Temperature 图 (Publication风格) - 带智能筛选
================================================================================

数据源:
  • structure_level_lindemann.csv (Step 7.8.5已完成的Ensemble级数据)
    - 每个结构在每个温度的Ensemble平均已完成
    - 包含delta_Pt, delta_Sn及其标准差
    - 数据已经过质量控制

功能: 
1. 直接使用Step 7.8.5的Ensemble级数据
2. 提取特定结构(如pt8sn6-1-best)或所有相关结构
3. 对异常温度点进行可选筛选
4. 绘制Publication风格的图

优势:
  • 与Step 7.8.5完全一致的数据和计算方法
  • 无需重新计算，避免精度损失
  • 数据质量有保证
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import argparse

# 全局配置
BASE_DIR = Path(__file__).parent
OUTPUT_DIR = BASE_DIR / 'results' / 'per_atom_dynamics'

# Step 7.8.2的数据文件 (默认: 50K步长新结果；旧路径: results/step7_8_2_alloy_series/...)
LINDEMANN_FILE = BASE_DIR / 'results' / 'step7_8_2_alloy_series_sup86_50k' / 'structure_level_lindemann.csv'

# 元素颜色和标记
ELEMENT_COLORS = {
    'Pt': '#1f77b4',  # 蓝色
    'Sn': '#ff7f0e',  # 橙色
}

ELEMENT_MARKERS = {
    'Pt': 'o',  # 圆形
    'Sn': 's',  # 方形
}

# Lindemann阈值
HIGH_LINDEMANN_THRESHOLD = 0.10


def load_per_atom_as_structure_level(per_atom_file, structure_filter=None):
    """
    读取 per-atom Lindemann CSV，按 原子→元素→run→温度 三级平均，
    生成与 structure_level_lindemann.csv 兼容的 DataFrame。

    聚合层级（更科学）：
      1. 每个 run 内：不做原子内平均（保留原子粒度）
      2. 按 (结构, 温度, element) 分组：对全部 run × 全部原子取 mean/std
      这样 δ_Pt = mean(全部 Pt 原子 × 全部 run) ，体现真实 ensemble 统计。

    返回的 DataFrame 列：
      structure, composition, temperature, delta_Pt, delta_Pt_std,
      delta_Sn, delta_Sn_std, n_pt, n_sn, anomaly_flag
    """
    import re as _re

    per_atom_path = Path(per_atom_file)
    if not per_atom_path.is_absolute():
        per_atom_path = BASE_DIR / per_atom_path
    if not per_atom_path.exists():
        print(f"  [ERROR] per-atom 文件不存在: {per_atom_path}")
        return None

    print(f"  [INFO] 读取 per-atom 文件: {per_atom_path.name}")
    df = pd.read_csv(per_atom_path)

    # 标准化列名
    df = df.rename(columns={'温度(K)': 'temperature', '目录': 'path',
                             '结构': 'structure_col', 'lindemann_index': 'delta'})

    # 从路径推断结构名（结构列通常为 NaN）
    def _struct_from_path(p):
        m = _re.search(r'/(pt\d+sn\d+[^/]*)/T\d+', str(p), _re.IGNORECASE)
        return m.group(1).lower() if m else 'unknown'

    df['structure'] = df['path'].apply(_struct_from_path)

    if structure_filter:
        df = df[df['structure'].str.contains(structure_filter, case=False, na=False)]
        if len(df) == 0:
            print(f"  [ERROR] 过滤后无数据: structure_filter='{structure_filter}'")
            return None

    print(f"  [INFO] 共 {len(df)} 原子行，结构: {sorted(df['structure'].unique())}")
    print(f"  [INFO] 温度点: {sorted(df['temperature'].unique())}")

    # ── 三级聚合：结构 × 温度 × 元素 → mean/std（跨全部run×全部原子）──
    agg = (df.groupby(['structure', 'temperature', 'element'])['delta']
             .agg(mean='mean', std='std')
             .reset_index())

    # 透视成宽格式（每温度一行）
    mean_wide = agg.pivot_table(index=['structure', 'temperature'],
                                columns='element', values='mean').reset_index()
    std_wide  = agg.pivot_table(index=['structure', 'temperature'],
                                columns='element', values='std').reset_index()
    mean_wide.columns.name = None
    std_wide.columns.name  = None

    # 重命名列
    rename_m = {e: f'delta_{e}' for e in mean_wide.columns if e not in ('structure','temperature')}
    rename_s = {e: f'delta_{e}_std' for e in std_wide.columns if e not in ('structure','temperature')}
    mean_wide = mean_wide.rename(columns=rename_m)
    std_wide  = std_wide.rename(columns=rename_s)
    result = mean_wide.merge(std_wide, on=['structure', 'temperature'], how='left')

    # 补充 composition / n_pt / n_sn
    def _comp(s):
        m = _re.search(r'pt(\d+)sn(\d+)', str(s).lower())
        return (f"pt{m.group(1)}sn{m.group(2)}", int(m.group(1)), int(m.group(2))) if m else (s, None, None)
    result[['composition', 'n_pt', 'n_sn']] = result['structure'].apply(
        lambda x: pd.Series(_comp(x)))
    result['anomaly_flag'] = ''
    result = result.sort_values(['structure', 'temperature']).reset_index(drop=True)

    n_runs = df.groupby(['structure','temperature'])['path'].nunique().mean()
    n_atoms = df.groupby(['structure','temperature','element']).size().mean()
    print(f"  [INFO] 平均 {n_runs:.0f} runs/温度，平均 {n_atoms:.0f} 原子/run/元素")
    print(f"  [INFO] 生成结构级数据: {len(result)} 行（{result['structure'].nunique()} 结构，{result['temperature'].nunique()} 温度）")

    return result


def load_structure_level_lindemann():
    """
    加载Step 7.8.5的structure-level Lindemann数据
    这是已经完成Ensemble平均的高质量数据
    """
    if not LINDEMANN_FILE.exists():
        print(f"  [ERROR] 文件不存在: {LINDEMANN_FILE}")
        return None
    
    try:
        df = pd.read_csv(LINDEMANN_FILE, encoding='utf-8-sig')
        print(f"  [INFO] 加载Step 7.8.5数据: {len(df)} 条记录")
        print(f"  [INFO] 包含结构: {df['structure'].nunique()} 个")
        print(f"  [INFO] 温度范围: {sorted(df['temperature'].unique())}")
        
        # 检查是否有异常标记列
        if 'anomaly_flag' in df.columns:
            anomalies = df[df['anomaly_flag'].notna() & (df['anomaly_flag'] != '') & (df['anomaly_flag'] != 'nan')]
            if len(anomalies) > 0:
                print(f"  [WARNING] 检测到 {len(anomalies)} 个异常标记:")
                for _, row in anomalies.iterrows():
                    print(f"    - {row['structure']} @ {row['temperature']:.0f}K: {row['anomaly_flag']}")
        
        return df
    except Exception as e:
        print(f"  [ERROR] 读取文件失败: {e}")
        return None




def extract_data_for_structure(df, structure):
    """
    从structure_level_lindemann.csv提取某个结构的数据
    
    参数:
        df: structure_level_lindemann数据
        structure: 结构名称 (如 'pt8sn6')
    
    返回:
        DataFrame with columns: temp, element, delta_mean, delta_std
    """
    # 筛选包含该结构的所有记录
    df_struct = df[df['structure'].str.contains(structure, case=False)].copy()
    
    if len(df_struct) == 0:
        print(f"  [ERROR] 未找到结构 {structure} 的数据")
        return None
    
    print(f"\n[*] 提取 {structure} 数据...")
    print(f"  找到 {len(df_struct)} 个温度点")
    
    # 自动筛除有异常标记的数据
    if 'anomaly_flag' in df_struct.columns:
        anomalies = df_struct[df_struct['anomaly_flag'].notna() & (df_struct['anomaly_flag'] != '') & (df_struct['anomaly_flag'] != 'nan')]
        if len(anomalies) > 0:
            print(f"  [*] 自动筛除 {len(anomalies)} 个标记为异常的温度点:")
            for _, row in anomalies.iterrows():
                print(f"    - {row['temperature']:.0f}K: {row['anomaly_flag']}")
            df_struct = df_struct[(df_struct['anomaly_flag'].isna()) | (df_struct['anomaly_flag'] == '') | (df_struct['anomaly_flag'] == 'nan')]

    # ── 若同一组分有多个结构变体（如 pt8sn6-1-best 和 pt8sn6-1-best-2），
    #    按温度取 ensemble 平均，合并为每温度一行 ──
    n_variants = df_struct['structure'].nunique()
    if n_variants > 1:
        print(f"  [*] 发现 {n_variants} 个结构变体: {df_struct['structure'].unique().tolist()}")
        print(f"  [*] 按温度取 ensemble 平均后合并为单曲线")
        df_struct = df_struct.groupby('temperature').agg(
            delta_Pt=('delta_Pt', 'mean'),
            delta_Sn=('delta_Sn', 'mean'),
            delta_Pt_std=('delta_Pt', 'std'),
            delta_Sn_std=('delta_Sn', 'std'),
        ).reset_index()
    else:
        df_struct = df_struct.rename(columns={
            'delta_Pt_std': 'delta_Pt_std',
            'delta_Sn_std': 'delta_Sn_std'
        })

    # 转换为长格式 (便于绘图)
    stats_list = []
    
    for _, row in df_struct.iterrows():
        temp = row['temperature']
        
        # Pt数据
        stats_list.append({
            'temp': temp,
            'element': 'Pt',
            'delta_mean': row['delta_Pt'],
            'delta_std': row['delta_Pt_std']
        })
        
        # Sn数据
        stats_list.append({
            'temp': temp,
            'element': 'Sn',
            'delta_mean': row['delta_Sn'],
            'delta_std': row['delta_Sn_std']
        })
    
    df_stats = pd.DataFrame(stats_list)
    df_stats = df_stats.sort_values(['element', 'temp'])
    
    print(f"  温度范围: {sorted(df_stats['temp'].unique())}")
    print(f"  元素: {df_stats['element'].unique().tolist()}")
    
    return df_stats


def filter_anomalous_temperatures(df_stats, filter_method='200K_vs_300K'):
    """
    筛选异常温度点
    
    参数:
        df_stats: 统计数据
        filter_method: 筛选方法
            - '200K_vs_300K': 如果δ(200K) > δ(300K)，移除200K
            - 'none': 不筛选
    
    返回:
        筛选后的df_stats
    """
    if filter_method == 'none':
        print("\n[*] 不进行温度筛选")
        return df_stats
    
    if filter_method == '200K_vs_300K':
        print("\n[*] 检查200K异常...")
        
        df_filtered_list = []
        removed_temps = set()
        
        for element in df_stats['element'].unique():
            df_elem = df_stats[df_stats['element'] == element].copy()
            
            # 检查200K vs 300K
            delta_200 = df_elem[df_elem['temp'] == 200]['delta_mean']
            delta_300 = df_elem[df_elem['temp'] == 300]['delta_mean']
            
            if len(delta_200) > 0 and len(delta_300) > 0:
                if delta_200.values[0] > delta_300.values[0]:
                    print(f"  {element}: δ(200K)={delta_200.values[0]:.5f} > δ(300K)={delta_300.values[0]:.5f}")
                    print(f"         → 移除200K点")
                    df_elem = df_elem[df_elem['temp'] != 200]
                    removed_temps.add(200)
            
            df_filtered_list.append(df_elem)
        
        df_filtered = pd.concat(df_filtered_list, ignore_index=True)
        
        if len(removed_temps) > 0:
            print(f"  已移除温度点: {sorted(removed_temps)}")
        else:
            print(f"  ✓ 无异常温度点")
        
        return df_filtered
    
    return df_stats
    
    print(f"\n  [OK] 统计完成")
    print(f"  温度点数: {len(df_stats['temp'].unique())}")
    print(f"  元素数: {len(df_stats['element'].unique())}")
    
    return df_stats


def plot_lindemann_publication(df_stats, structure, output_dir, markersize=12,
                               figsize=(10, 8), show_threshold=True):
    """
    绘制Publication风格的Lindemann Index vs Temperature图
    """
    print("\n" + "=" * 80)
    print("绘制 Lindemann Index vs Temperature (Publication风格)")
    print("=" * 80)
    
    if df_stats is None or len(df_stats) == 0:
        print("  [ERROR] 无数据可绘制")
        return
    
    # 设置全局字体为Arial
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.sans-serif'] = ['Arial']
    plt.rcParams['mathtext.default'] = 'regular'
    
    # 创建图表
    fig, ax = plt.subplots(figsize=figsize)
    
    # 设置透明背景
    fig.patch.set_alpha(0.0)
    ax.patch.set_alpha(0.0)
    
    print("\n[*] 绘制数据曲线...")
    
    # 绘制每个元素
    for element in ['Pt', 'Sn']:
        elem_data = df_stats[df_stats['element'] == element].sort_values('temp')
        
        if len(elem_data) == 0:
            print(f"  [WARNING] 无 {element} 数据")
            continue
        
        temps = elem_data['temp'].values
        delta_mean = elem_data['delta_mean'].values
        delta_std = elem_data['delta_std'].values
        
        color = ELEMENT_COLORS.get(element, 'gray')
        marker = ELEMENT_MARKERS.get(element, 'o')
        
        # 绘制误差带 (半透明)
        ax.fill_between(temps,
                        delta_mean - delta_std,
                        delta_mean + delta_std,
                        color=color, alpha=0.2, zorder=1)
        
        # 绘制主线 (线宽4)
        ax.plot(temps, delta_mean,
                marker=marker, color=color, linewidth=4, markersize=markersize,
                label=element, zorder=3, markeredgewidth=1.5, markeredgecolor='white')
        
        print(f"    {element}: {len(temps)} 个温度点")
    
    # 添加Lindemann阈值线
    if show_threshold:
        ax.axhline(HIGH_LINDEMANN_THRESHOLD, color='red', linestyle='--', 
                   linewidth=2.5, alpha=0.6, zorder=2)
    
    # 坐标轴标签 (字体34, 不加粗, Arial)
    ax.set_xlabel('Temperature (K)', fontsize=34)
    ax.set_ylabel('Lindemann Index δ', fontsize=34)
    
    # 图例 (无边框, 字体28, 左上角)
    legend = ax.legend(fontsize=28, loc='upper left', frameon=False)
    
    # 坐标轴刻度字体 (28号)
    ax.tick_params(axis='both', which='major', labelsize=28, 
                   direction='out', length=8, width=1.5)
    
    # 设置刻度线朝外
    ax.tick_params(axis='x', direction='out')
    ax.tick_params(axis='y', direction='out')
    
    # 自动调整刻度数量 (4-7个, 尽量对称整数)
    ax.locator_params(axis='x', nbins=6)
    ax.locator_params(axis='y', nbins=6)
    
    # 移除网格线
    ax.grid(False)
    
    # 显示4个边框
    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(True)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(True)
    
    # 设置边框线宽
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
    
    # 调整布局
    plt.tight_layout()
    
    # 保存图片 (透明背景)
    fig_path = output_dir / f'{structure}_lindemann_filtered.png'
    
    plt.savefig(fig_path, dpi=300, bbox_inches='tight', 
                transparent=True, facecolor='none')
    
    print(f"\n[SAVED] {fig_path}")
    print(f"  尺寸: 10×8 英寸")
    print(f"  DPI: 300")
    print(f"  背景: 透明")
    print(f"  字体: Arial")
    print(f"  数据点大小: {markersize}")
    print(f"  线宽: 4")
    
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description='绘制Publication风格的Lindemann Index图 (基于Step 7.8.5数据)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # 绘制pt8sn6 (自动使用Step 7.8.2标记的异常)
  python %(prog)s --structure pt8sn6
  
  # 额外检查200K vs 300K
  python %(prog)s --structure pt8sn6 --check-200K
  
  # 不进行任何筛选 (保留所有数据)
  python %(prog)s --structure pt8sn6 --no-filter
  
  # 自定义数据点大小
  python %(prog)s --structure pt8sn6 --markersize 15
        """
    )
    
    parser.add_argument('--structure', type=str, default='pt8sn6',
                       help='结构名称 (如: pt8sn6, pt8sn1, etc.)')
    parser.add_argument('--no-filter', action='store_true',
                       help='不筛选任何异常 (包括Step 7.8.2标记的)')
    parser.add_argument('--check-200K', action='store_true',
                       help='额外检查并移除 δ(200K) > δ(300K) 的点')
    parser.add_argument('--markersize', type=float, default=12,
                       help='数据点大小 (默认: 12)')
    parser.add_argument('--lindemann-file', type=str, default=None,
                       help=('指定 structure_level_lindemann.csv 路径。'
                             '默认: results/step7_8_2_alloy_series_sup86_50k/structure_level_lindemann.csv；'
                             '旧数据: results/step7_8_2_alloy_series/structure_level_lindemann.csv'))
    parser.add_argument('--no-threshold', action='store_true',
                       help='隐去 δ=0.1 的阈值参考线')
    parser.add_argument('--figsize', type=str, default='10x8',
                       help='图片尺寸，格式: 宽x高 (默认: 10x8，例: 12x9)')
    parser.add_argument('--per-atom-file', type=str, default=None,
                       help=('直接读取 per-atom Lindemann CSV，在脚本内完成原子→元素→温度三级平均。'
                             '比 structure_level_lindemann.csv 更科学（保留原子粒度统计）。'
                             '例: data/lindemann/per-atoms/sup86-50k/per_atom_master_run_20260311_173823.csv'))
    
    args = parser.parse_args()

    # 解析 figsize
    try:
        fw, fh = [float(x) for x in args.figsize.lower().split('x')]
        figsize = (fw, fh)
    except Exception:
        print(f"  [WARNING] --figsize 格式错误 '{args.figsize}'，使用默认 10x8")
        figsize = (10, 8)

    # 若用户通过命令行指定了 --lindemann-file，覆盖全局路径
    global LINDEMANN_FILE
    if args.lindemann_file:
        p = Path(args.lindemann_file)
        LINDEMANN_FILE = p if p.is_absolute() else BASE_DIR / p
    
    print("=" * 80)
    print("Lindemann Index Publication Plot (基于Step 7.8.2数据)")
    print(f"  结构: {args.structure}")
    print(f"  数据源: {LINDEMANN_FILE}")
    
    if args.no_filter:
        print(f"  筛选: 不筛选 (保留所有数据)")
    else:
        print(f"  筛选: 自动筛除Step 7.8.2标记的异常")
        if args.check_200K:
            print(f"         + 额外检查200K vs 300K")
    
    print(f"  数据点大小: {args.markersize}")
    print(f"  图片尺寸: {figsize[0]:.0f}x{figsize[1]:.0f}")
    print(f"  阈值线(δ=0.1): {'隐藏' if args.no_threshold else '显示'}")
    if args.per_atom_file:
        print(f"  数据模式: per-atom（原子→元素→温度三级平均）")
        print(f"  per-atom文件: {args.per_atom_file}")
    else:
        print(f"  数据模式: structure-level（来自 step7_8_2 输出）")
    print("=" * 80)

    # ── 加载数据：per-atom 优先，否则读 structure_level_lindemann.csv ──
    if args.per_atom_file:
        df = load_per_atom_as_structure_level(args.per_atom_file,
                                              structure_filter=args.structure)
    else:
        df = load_structure_level_lindemann()
    
    if df is None:
        print("\n[ERROR] 无法加载数据")
        return
    
    # 如果不筛选，移除anomaly_flag列
    if args.no_filter and 'anomaly_flag' in df.columns:
        print("\n[*] --no-filter 模式: 忽略所有异常标记")
        df['anomaly_flag'] = ''
    
    # 提取指定结构的数据 (自动筛除有anomaly_flag的)
    df_stats = extract_data_for_structure(df, args.structure)
    
    if df_stats is None:
        return
    
    # 额外的200K vs 300K检查
    if args.check_200K and not args.no_filter:
        df_stats = filter_anomalous_temperatures(df_stats, filter_method='200K_vs_300K')
    
    # 确保输出目录存在
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # ── 导出审查用宽格式 CSV：每温度一行，含 δ_Pt、δ_Sn、比值 ──
    df_wide = df_stats.pivot(index='temp', columns='element', values=['delta_mean', 'delta_std'])
    df_wide.columns = [f'delta_{elem}' if stat == 'delta_mean' else f'delta_{elem}_std'
                       for stat, elem in df_wide.columns]
    df_wide = df_wide.reset_index().rename(columns={'temp': 'temperature'})
    # 比值：δ_Pt / δ_Sn（Sn 为 0 时置 NaN）
    if 'delta_Pt' in df_wide.columns and 'delta_Sn' in df_wide.columns:
        df_wide['ratio_Pt_Sn'] = df_wide['delta_Pt'] / df_wide['delta_Sn'].replace(0, float('nan'))
        df_wide['ratio_Pt_Sn'] = df_wide['ratio_Pt_Sn'].round(4)
    # 整理列顺序
    col_order = ['temperature', 'delta_Pt', 'delta_Pt_std', 'delta_Sn', 'delta_Sn_std', 'ratio_Pt_Sn']
    df_wide = df_wide[[c for c in col_order if c in df_wide.columns]]
    df_wide = df_wide.sort_values('temperature')

    stats_path = OUTPUT_DIR / f'{args.structure}_lindemann_summary.csv'
    df_wide.to_csv(stats_path, index=False, encoding='utf-8-sig', float_format='%.6f')
    print(f"\n[SAVED] {stats_path}")
    print(df_wide.to_string(index=False))
    
    # 绘制
    plot_lindemann_publication(df_stats, args.structure, OUTPUT_DIR,
                               markersize=args.markersize,
                               figsize=figsize,
                               show_threshold=not args.no_threshold)
    
    print("\n" + "=" * 80)
    print("完成！")
    print("=" * 80)


if __name__ == '__main__':
    main()
