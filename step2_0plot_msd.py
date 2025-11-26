#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step3: 绘制MSD曲线 - 高速版 (文件索引优化)
========================================
创建时间: 2025-10-16
最后更新: 2025-11-26
作者: GitHub Copilot

功能:
    1. 直接扫描 GMX MSD 数据目录构建文件索引
    2. 自动发现所有系统和温度
    3. 按系统和温度绘制MSD曲线
    4. 支持异常值筛选 (可选, 使用 step1_2 的 large_D_outliers.csv)
    5. 同时显示有效runs和被筛选runs (不同样式)

工作流程:
    - 可独立运行,无需先运行其他脚本
    - 如需筛选异常runs,需要 step1_2 生成的 large_D_outliers.csv

依赖关系:
    - 可选: step1_2 的 large_D_outliers.csv 用于筛选异常runs
    - 无必需依赖 (系统列表从文件直接发现)

性能优化:
    1. 一次性构建文件索引 (避免重复rglob)
    2. 缓存文件路径映射
    3. 预计速度提升: 10-20倍

输入:
    - GMX MSD原始数据 (.xvg文件,从GMX_DATA_DIRS读取)
    - [可选] step1_2的异常清单 (large_D_outliers.csv)

输出: results/msd_curves/
    ├── {system}_all_temps_GMX.png      (各系统的MSD曲线图)
    └── filtering_statistics.txt        (筛选统计报告)

使用方法:
--------
python step3_plot_msd.py              # 默认启用异常值筛选
python step3_plot_msd.py --nofilter   # 关闭筛选，绘制所有曲线
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from datetime import datetime
import re
import warnings
from collections import defaultdict
from tqdm import tqdm
import argparse
warnings.filterwarnings('ignore')

# ===== 全局配置 =====
BASE_DIR = Path(__file__).parent  # workflow目录

# 数据源: step1_2 的输出结果
# - ensemble_analysis_results.csv: 有效runs的集成平均 (用于获取系统列表)
# - large_D_outliers.csv: 异常runs清单 (用于筛选绘图)
RESULTS_CSV = BASE_DIR / 'results' / 'ensemble_analysis_results.csv'
OUTLIERS_CSV = BASE_DIR / 'results' / 'large_D_outliers.csv'

GMX_DATA_DIRS = [
    # BASE_DIR / 'data' / 'gmx_msd' / 'collected_gmx_msd',
    # BASE_DIR / 'data' / 'gmx_msd' / 'gmx_msd_results_20251015_184626_collected',
    # 新版unwrap per-atom MSD数据 (2025-11-18)
    BASE_DIR / 'data' / 'gmx_msd' / 'unwrap' / 'gmx_msd_results_20251118_152614',
    # BASE_DIR / 'data' / 'gmx_msd' / 'unwrap' / 'air' / 'gmx_msd_results_20251124_170114'  # 🌬️ 气象数据
]

OUTPUT_DIR = BASE_DIR / 'results' / 'msd_curves'
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

# ===== 异常值筛选开关（可通过命令行参数覆盖） =====
# 默认值：True = 启用筛选，False = 绘制所有曲线
# 命令行：--nofilter 将覆盖此设置
ENABLE_OUTLIER_FILTERING = True  # 默认启用筛选

# ===== 系统过滤配置 (可选) =====
# 如果为空列表,则绘制所有系统
# 如果指定模式,则只绘制匹配的系统
SYSTEM_FILTER = {
    'include_patterns': [
        # 示例: 只绘制Pt8相关系统
        r'^pt8',           # pt8开头的所有系统
        # r'pt8sn\d+',       # pt8sn0, pt8sn1, ..., pt8sn10
        # r'pt\d+sn\d+',     # 所有ptXsnY格式
        # r'^\d+$',          # 气象数据 (纯数字命名，如 68)
    ],
    'exclude_patterns': [
        # 示例: 排除含氧系统
        # r'^[Oo]\d+',       # O1, O2, o3, o4等开头
        # r'[Oo]\d+Pt',      # O2Pt4Sn6等
        # r'Pt\d+Sn\d+O',    # Pt2Sn2O1等
    ]
}

COLORS = {
    'Pt': '#E74C3C',
    'Sn': '#3498DB',
    'PtSn': '#2ECC71'
}

plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False


def load_large_D_outliers(enable_filtering=None):
    """
    加载大D值异常run清单
    
    Parameters:
    -----------
    enable_filtering : bool, optional
        是否启用筛选。如果为 None，使用全局配置
    
    Returns:
    --------
    outliers : set
        异常文件路径集合
    """
    # 使用参数或全局配置
    if enable_filtering is None:
        enable_filtering = ENABLE_OUTLIER_FILTERING
    
    # 检查是否启用筛选
    if not enable_filtering:
        print(f"   [!] Outlier filtering is DISABLED (--nofilter)")
        print(f"   [!] Will plot ALL runs (including outliers)")
        return set()
    
    try:
        df_outliers = pd.read_csv(OUTLIERS_CSV)
        outlier_files = set(df_outliers['filepath'].values)
        print(f"   [OK] Loaded outlier list: {len(outlier_files)} files")
        print(f"        (These runs will be excluded from plotting)")
        return outlier_files
    except FileNotFoundError:
        print(f"   [!] Outlier file not found: {OUTLIERS_CSV}")
        print(f"   [!] Will plot all runs without filtering")
        return set()


def build_file_index(outlier_files=None):
    """
    构建文件索引 - 一次性扫描所有文件,分别索引有效runs和被筛选runs
    
    Parameters:
    -----------
    outlier_files : set, optional
        异常文件路径集合
    
    Returns:
    --------
    file_index : dict
        {(composition, temperature, element): [file_path1, file_path2, ...]} - 有效runs
    file_index_filtered : dict
        {(composition, temperature, element): [file_path1, file_path2, ...]} - 被筛选runs
    filter_stats : dict
        {(composition, temperature, element): {'total': int, 'kept': int, 'filtered': int}}
    """
    if outlier_files is None:
        outlier_files = set()
    
    print("\n[*] Building file index...")
    file_index = defaultdict(list)
    file_index_filtered = defaultdict(list)  # 新增: 被筛选runs的索引
    filter_stats = defaultdict(lambda: {'total': 0, 'kept': 0, 'filtered': 0})
    
    total_files = 0
    filtered_files = 0
    seen_files = set()  # 用于去重
    
    for gmx_dir in GMX_DATA_DIRS:
        print(f"  Scanning: {gmx_dir.name}...")
        
        # 查找所有MSD .xvg文件
        for xvg_file in gmx_dir.rglob("*_msd_*.xvg"):
            # 去重检查
            try:
                normalized_path = xvg_file.resolve()
                if normalized_path in seen_files:
                    continue  # 跳过重复文件
                seen_files.add(normalized_path)
            except:
                pass  # 如果resolve()失败,继续处理
            try:
                # 解析路径: 从父目录层级提取信息
                parts = xvg_file.parts
                
                # 提取element (从文件名)
                filename = xvg_file.stem  # T1000.r7.gpu0_msd_Pt
                if '_msd_' in filename:
                    element = filename.split('_msd_')[-1]  # Pt, Sn, PtSn
                else:
                    continue
                
                # 提取temperature和composition (从路径向上查找)
                temperature = None
                composition = None
                for i in range(len(parts)-1, 0, -1):
                    if parts[i].endswith('K'):
                        temperature = parts[i]
                        composition = parts[i-1]
                        break
                
                if not temperature or not composition:
                    continue
                
                key = (composition, temperature, element)
                filter_stats[key]['total'] += 1
                total_files += 1
                
                # 检查是否是异常run
                if str(xvg_file) in outlier_files:
                    filter_stats[key]['filtered'] += 1
                    file_index_filtered[key].append(xvg_file)  # 添加到筛选索引
                    filtered_files += 1
                else:
                    filter_stats[key]['kept'] += 1
                    file_index[key].append(xvg_file)
                
            except Exception as e:
                continue
    
    print(f"   [OK] Indexed {total_files - filtered_files} valid files")
    print(f"   [OK] Indexed {filtered_files} filtered files")
    print(f"   [OK] Filtered out {filtered_files} outlier runs")
    
    return file_index, file_index_filtered, filter_stats


def extract_base_system(composition_name):
    """提取基础体系名称"""
    # 转换为字符串（防止整数等类型）
    composition_name = str(composition_name)
    match = re.match(r'^(Cv)-\d+$', composition_name)
    if match:
        return match.group(1)
    return composition_name


def group_compositions_by_system(df):
    """按体系分组"""
    df['base_system'] = df['composition'].apply(extract_base_system)
    
    system_groups = {}
    for base_sys in df['base_system'].unique():
        comps = df[df['base_system'] == base_sys]['composition'].unique()
        system_groups[base_sys] = sorted(comps)
    
    return system_groups


def filter_systems(system_groups, include_patterns=None, exclude_patterns=None):
    """
    根据正则表达式模式过滤系统
    
    Parameters:
    -----------
    system_groups : dict
        {base_system: [compositions]}
    include_patterns : list of str
        包含模式列表 (正则表达式)
        如果为空,则包含所有系统
    exclude_patterns : list of str
        排除模式列表 (正则表达式)
    
    Returns:
    --------
    filtered_groups : dict
        过滤后的系统字典
    """
    if include_patterns is None:
        include_patterns = []
    if exclude_patterns is None:
        exclude_patterns = []
    
    filtered_groups = {}
    
    for base_sys, comps in system_groups.items():
        # 检查是否匹配include模式
        if include_patterns:
            included = any(re.search(pattern, base_sys, re.IGNORECASE) 
                          for pattern in include_patterns)
            if not included:
                continue
        
        # 检查是否匹配exclude模式
        if exclude_patterns:
            excluded = any(re.search(pattern, base_sys, re.IGNORECASE) 
                          for pattern in exclude_patterns)
            if excluded:
                continue
        
        # 通过过滤
        filtered_groups[base_sys] = comps
    
    return filtered_groups


def read_gmx_msd_xvg(filepath):
    """读取GMX MSD .xvg文件"""
    time_data = []
    msd_data = []
    
    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#') or line.startswith('@'):
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        t = float(parts[0])
                        msd_nm2 = float(parts[1])
                        msd_a2 = msd_nm2 * 100
                        time_data.append(t)
                        msd_data.append(msd_a2)
                    except ValueError:
                        continue
    except:
        return None, None
    
    if len(time_data) == 0:
        return None, None
    
    return np.array(time_data), np.array(msd_data)


def plot_system_all_temps_fast(base_system, compositions, df_all, file_index, file_index_filtered, filter_stats, max_temps=None, use_file_temps=False):
    """
    快速绘制 - 使用文件索引,同时绘制有效runs和被筛选runs
    
    Parameters:
    -----------
    base_system : str
        基础体系名称
    compositions : list
        组成列表
    df_all : DataFrame
        所有数据
    file_index : dict
        有效runs的文件索引
    file_index_filtered : dict
        被筛选runs的文件索引
    filter_stats : dict
        筛选统计
    max_temps : int, optional
        最大温度数
    use_file_temps : bool, optional
        是否从文件索引获取温度列表（用于 --nofilter 模式）
    """
    
    print(f"\n{'='*80}")
    print(f"System: {base_system}")
    if len(compositions) > 1:
        print(f"Compositions: {', '.join(compositions)}")
    print(f"{'='*80}")
    
    df_sys = df_all[df_all['composition'].isin(compositions)]
    
    # 获取温度列表
    if use_file_temps:
        # 从文件索引获取所有可用温度（用于 --nofilter）
        temp_set = set()
        for key in list(file_index.keys()) + list(file_index_filtered.keys()):
            comp_str, temp, elem = key
            # 检查组成是否匹配
            if any(str(c) == comp_str for c in compositions):
                temp_set.add(temp)
        
        if not temp_set:
            print(f"[X] No files found in file index")
            return
        
        temperatures = sorted(temp_set, key=lambda x: int(x.replace('K', '')))
        print(f"Temperatures (from files): {len(temperatures)} - {temperatures}")
    else:
        # 从 DataFrame 获取温度（默认行为）
        if len(df_sys) == 0:
            print(f"[X] No data in DataFrame")
            return
        
        temperatures = sorted(df_sys['temperature'].unique(),
                             key=lambda x: int(x.replace('K', '')))
    
    print(f"Data points: {len(df_sys)}")
    print(f"Temperatures: {len(temperatures)} - {temperatures}")
    
    if max_temps and len(temperatures) > max_temps:
        indices = np.linspace(0, len(temperatures)-1, max_temps, dtype=int)
        temperatures = [temperatures[i] for i in indices]
    
    n_temps = len(temperatures)
    n_cols = min(3, n_temps)
    n_rows = (n_temps + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(7*n_cols, 5*n_rows))
    if n_temps == 1:
        axes = [axes]
    else:
        axes = axes.flatten()
    
    # 一次性读取所有数据 + 找最大值
    print("\n[*] Loading MSD data...")
    
    msd_cache = {}  # {(temp, element): [(time, msd), ...]} - 有效runs
    msd_cache_filtered = {}  # {(temp, element): [(time, msd), ...]} - 被筛选runs
    global_max_msd = 0
    
    # 统计每个温度的筛选情况
    temp_filter_stats = {}  # {temp: {'total': int, 'kept': int, 'filtered': int}}
    
    for temp in temperatures:
        temp_stats = {'total': 0, 'kept': 0, 'filtered': 0}
        
        for element in ['Pt', 'Sn', 'PtSn']:
            msd_list = []
            msd_list_filtered = []
            
            # 从索引查找文件
            for comp in compositions:
                # 确保 composition 是字符串（文件索引的key是字符串）
                key = (str(comp), temp, element)
                
                # 累加统计
                if key in filter_stats:
                    temp_stats['total'] += filter_stats[key]['total']
                    temp_stats['kept'] += filter_stats[key]['kept']
                    temp_stats['filtered'] += filter_stats[key]['filtered']
                
                # 读取有效runs
                files = file_index.get(key, [])
                for filepath in files:
                    time, msd = read_gmx_msd_xvg(filepath)
                    if time is not None:
                        msd_list.append((time, msd))
                        
                        # 更新全局最大值
                        max_val = np.max(msd)
                        if max_val > global_max_msd:
                            global_max_msd = max_val
                
                # 读取被筛选runs
                files_filtered = file_index_filtered.get(key, [])
                for filepath in files_filtered:
                    time, msd = read_gmx_msd_xvg(filepath)
                    if time is not None:
                        msd_list_filtered.append((time, msd))
                        
                        # 更新全局最大值
                        max_val = np.max(msd)
                        if max_val > global_max_msd:
                            global_max_msd = max_val
            
            if msd_list:
                msd_cache[(temp, element)] = msd_list
            if msd_list_filtered:
                msd_cache_filtered[(temp, element)] = msd_list_filtered
        
        temp_filter_stats[temp] = temp_stats
    
    total_curves = sum(len(v) for v in msd_cache.values())
    total_curves_filt = sum(len(v) for v in msd_cache_filtered.values())
    print(f"  [OK] Loaded {total_curves} valid curves + {total_curves_filt} filtered curves")
    print(f"  [OK] Global max MSD: {global_max_msd:.2f} A^2")
    
    # 绘制
    print("\n[*] Plotting...")
    
    for idx, temp in enumerate(temperatures):
        ax = axes[idx]
        has_data = False
        
        # 获取该温度的筛选统计
        stats = temp_filter_stats.get(temp, {'total': 0, 'kept': 0, 'filtered': 0})
        
        # === 绘制有效runs (实线,饱和颜色) ===
        for element in ['Pt', 'Sn', 'PtSn']:
            msd_list = msd_cache.get((temp, element), [])
            
            if not msd_list:
                continue
            
            # 对齐时间轴
            all_times = [t for t, m in msd_list]
            min_time = max(t.min() for t in all_times)
            max_time = min(t.max() for t in all_times)
            common_time = np.linspace(min_time, max_time, 500)
            
            # 插值
            msd_interp_list = []
            for time, msd in msd_list:
                try:
                    msd_interp = np.interp(common_time, time, msd)
                    msd_interp_list.append(msd_interp)
                    
                    # 透明单次 (有效runs)
                    ax.plot(common_time, msd_interp, alpha=0.2, linewidth=1,
                           color=COLORS[element], zorder=1)
                except:
                    continue
            
            if not msd_interp_list:
                continue
            
            # 集合平均 (有效runs)
            if len(msd_interp_list) >= 2:
                msd_array = np.array(msd_interp_list)
                msd_mean = np.mean(msd_array, axis=0)
                msd_std = np.std(msd_array, axis=0, ddof=1)
                
                ax.plot(common_time, msd_mean, '-', linewidth=3,
                       color=COLORS[element],
                       label=f'{element} Valid (n={len(msd_interp_list)})', zorder=10)
                
                ax.fill_between(common_time, msd_mean - msd_std, msd_mean + msd_std,
                               alpha=0.2, color=COLORS[element], zorder=5)
            else:
                ax.plot(common_time, msd_interp_list[0], '-', linewidth=3,
                       color=COLORS[element],
                       label=f'{element} Valid (n=1)', zorder=10)
            
            has_data = True
        
        # === 绘制被筛选runs (虚线,灰色) ===
        for element in ['Pt', 'Sn', 'PtSn']:
            msd_list_filt = msd_cache_filtered.get((temp, element), [])
            
            if not msd_list_filt:
                continue
            
            # 对齐时间轴
            all_times_filt = [t for t, m in msd_list_filt]
            min_time_filt = max(t.min() for t in all_times_filt)
            max_time_filt = min(t.max() for t in all_times_filt)
            common_time_filt = np.linspace(min_time_filt, max_time_filt, 500)
            
            # 插值
            msd_interp_list_filt = []
            for time, msd in msd_list_filt:
                try:
                    msd_interp = np.interp(common_time_filt, time, msd)
                    msd_interp_list_filt.append(msd_interp)
                    
                    # 透明单次 (被筛选runs - 灰色虚线)
                    ax.plot(common_time_filt, msd_interp, '--', alpha=0.15, linewidth=1,
                           color='gray', zorder=0.5)
                except:
                    continue
            
            if not msd_interp_list_filt:
                continue
            
            # 集合平均 (被筛选runs - 深灰色虚线)
            if len(msd_interp_list_filt) >= 1:
                msd_array_filt = np.array(msd_interp_list_filt)
                msd_mean_filt = np.mean(msd_array_filt, axis=0)
                
                ax.plot(common_time_filt, msd_mean_filt, '--', linewidth=2.5,
                       color='dimgray', alpha=0.6,
                       label=f'{element} Filtered (n={len(msd_interp_list_filt)})', zorder=8)
            
            has_data = True
        
        if not has_data:
            ax.text(0.5, 0.5, 'No Data', ha='center', va='center',
                   transform=ax.transAxes, fontsize=14)
        
        # 添加筛选信息到标题
        if stats['total'] > 0:
            filter_text = f" [{stats['kept']}/{stats['total']} runs]"
            if stats['filtered'] > 0:
                title_text = f"{temp}{filter_text}"
                title_color = 'darkred' if stats['filtered'] > stats['total'] * 0.3 else 'black'
            else:
                title_text = f"{temp}{filter_text}"
                title_color = 'black'
        else:
            title_text = f"{temp}"
            title_color = 'black'
        
        ax.set_xlabel('Time (ps)', fontsize=11, fontweight='bold')
        ax.set_ylabel('MSD (A^2)', fontsize=11, fontweight='bold')
        ax.set_title(title_text, fontsize=12, fontweight='bold', color=title_color)
        ax.legend(loc='upper left', fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        
        if global_max_msd > 0:
            ax.set_ylim(0, global_max_msd * 1.1)
    
    for idx in range(n_temps, len(axes)):
        axes[idx].axis('off')
    
    title_text = f'{base_system} - MSD Curves\n'
    if len(compositions) > 1:
        title_text += f'Combined: {", ".join(compositions)}\n'
    title_text += 'Thick: Average | Shaded: ±1σ | Thin: Individual | Y-axis unified'
    
    fig.suptitle(title_text, fontsize=14, fontweight='bold', y=0.995)
    
    # 使用try-except避免tight_layout错误
    try:
        plt.tight_layout(rect=[0, 0, 1, 0.99])
    except Exception as e:
        print(f"  [!] Warning: tight_layout failed ({e}), continuing...")
    
    output_file = OUTPUT_DIR / f'{base_system}_all_temps_GMX.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    print(f"  [OK] Saved: {output_file.name}")
    print(f"{'='*80}\n")


def main(enable_filtering=None):
    """
    主函数
    
    Parameters:
    -----------
    enable_filtering : bool, optional
        是否启用异常值筛选。None = 使用全局配置
    """
    print("\n" + "="*80)
    print("MSD Curves - Fast Version (with file indexing + outlier filtering)")
    print("="*80)
    print(f"Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*80)
    
    # 0. 加载异常清单
    print("\n[0/5] Loading outlier list...")
    outlier_files = load_large_D_outliers(enable_filtering)
    
    # 1. 构建文件索引 (一次性,分别索引有效runs和被筛选runs)
    file_index, file_index_filtered, filter_stats = build_file_index(outlier_files)
    
    # 2. 从文件索引直接发现系统列表 (不再依赖 ensemble_analysis_results.csv)
    print("\n[2/5] Discovering systems from file index...")
    
    # 从文件索引提取所有 (composition, temperature) 组合
    all_compositions = set()
    for key in list(file_index.keys()) + list(file_index_filtered.keys()):
        comp, temp, elem = key
        all_compositions.add(comp)
    
    if not all_compositions:
        print(f"\n[X] Error: No files found in file index")
        return
    
    print(f"   [OK] Found {len(all_compositions)} compositions from files")
    
    # 创建一个虚拟 DataFrame 用于 group_compositions_by_system
    df = pd.DataFrame({'composition': list(all_compositions)})
    
    # 3. 按体系分组
    system_groups = group_compositions_by_system(df)
    
    print(f"\n[*] Systems identified: {len(system_groups)}")
    
    # 4. 应用系统过滤
    include_patterns = SYSTEM_FILTER.get('include_patterns', [])
    exclude_patterns = SYSTEM_FILTER.get('exclude_patterns', [])
    
    if include_patterns or exclude_patterns:
        print(f"\n[*] Applying system filters...")
        if include_patterns:
            print(f"    Include patterns: {include_patterns}")
        if exclude_patterns:
            print(f"    Exclude patterns: {exclude_patterns}")
        
        system_groups = filter_systems(system_groups, include_patterns, exclude_patterns)
        print(f"    Filtered: {len(system_groups)} systems selected")
    
    print(f"[*] Output dir: {OUTPUT_DIR}")
    
    print("\nSystem list:")
    for i, (base_sys, comps) in enumerate(sorted(system_groups.items()), 1):
        # 统计该系统的文件数
        n_files = 0
        for comp in comps:
            for key in list(file_index.keys()) + list(file_index_filtered.keys()):
                if key[0] == str(comp):
                    n_files += len(file_index.get(key, [])) + len(file_index_filtered.get(key, []))
        
        if len(comps) > 1:
            print(f"  {i:2d}. {base_sys:30s} ({n_files:3d} files) <- {len(comps)} comps")
        else:
            print(f"  {i:2d}. {base_sys:30s} ({n_files:3d} files)")
    
    print("\n" + "="*80)
    print("Start plotting...")
    print("="*80)
    
    success = 0
    failed = []
    
    # 始终从文件索引获取温度列表 (不再依赖 CSV)
    use_file_temps = True
    
    for idx, (base_sys, comps) in enumerate(sorted(system_groups.items()), 1):
        print(f"\n[{idx}/{len(system_groups)}] {base_sys}...")
        
        try:
            plot_system_all_temps_fast(base_sys, comps, df, file_index, file_index_filtered, filter_stats, 
                                      max_temps=None, use_file_temps=use_file_temps)
            success += 1
        except Exception as e:
            print(f"[X] Failed: {e}")
            import traceback
            traceback.print_exc()
            failed.append((base_sys, str(e)))
    
    print("\n" + "="*80)
    print("DONE!")
    print("="*80)
    print(f"Expected: {len(system_groups)}")
    print(f"Success:  {success}")
    print(f"Failed:   {len(failed)}")
    
    if failed:
        print("\nFailed systems:")
        for sys, err in failed:
            print(f"  - {sys}: {err}")
    
    print(f"\nOutput: {OUTPUT_DIR}")
    print("="*80)
    
    # 生成筛选统计报告
    print("\n[*] Generating filtering statistics report...")
    generate_filter_report(filter_stats, OUTPUT_DIR)


def generate_filter_report(filter_stats, output_dir):
    """生成筛选统计报告"""
    report_file = output_dir / 'filtering_statistics.txt'
    
    # 按体系和温度分组统计
    system_stats = defaultdict(lambda: defaultdict(lambda: {'total': 0, 'kept': 0, 'filtered': 0}))
    
    for (comp, temp, elem), stats in filter_stats.items():
        # 提取基础体系名称
        base_sys = extract_base_system(comp)
        system_stats[base_sys][temp]['total'] += stats['total']
        system_stats[base_sys][temp]['kept'] += stats['kept']
        system_stats[base_sys][temp]['filtered'] += stats['filtered']
    
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write("="*80 + "\n")
        f.write("Step1 Filtering Statistics - By System and Temperature\n")
        f.write("="*80 + "\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("="*80 + "\n\n")
        
        # 总体统计
        total_runs = sum(s['total'] for stats in filter_stats.values() for s in [stats])
        kept_runs = sum(s['kept'] for stats in filter_stats.values() for s in [stats])
        filtered_runs = sum(s['filtered'] for stats in filter_stats.values() for s in [stats])
        
        f.write("Overall Statistics:\n")
        f.write(f"  Total runs:    {total_runs}\n")
        if total_runs > 0:
            f.write(f"  Kept runs:     {kept_runs} ({kept_runs/total_runs*100:.1f}%)\n")
            f.write(f"  Filtered runs: {filtered_runs} ({filtered_runs/total_runs*100:.1f}%)\n")
        else:
            f.write(f"  Kept runs:     {kept_runs}\n")
            f.write(f"  Filtered runs: {filtered_runs}\n")
        f.write("\n" + "="*80 + "\n\n")
        
        # 按体系详细统计
        for base_sys in sorted(system_stats.keys()):
            temps = system_stats[base_sys]
            
            f.write(f"System: {base_sys}\n")
            f.write("-"*80 + "\n")
            f.write(f"{'Temperature':<15} {'Total':<10} {'Kept':<10} {'Filtered':<12} {'Rate':<10}\n")
            f.write("-"*80 + "\n")
            
            for temp in sorted(temps.keys(), key=lambda x: int(x.replace('K', ''))):
                stats = temps[temp]
                if stats['total'] > 0:
                    rate = stats['kept'] / stats['total'] * 100
                    f.write(f"{temp:<15} {stats['total']:<10} {stats['kept']:<10} "
                           f"{stats['filtered']:<12} {rate:>6.1f}%\n")
            
            # 体系总计
            sys_total = sum(s['total'] for s in temps.values())
            sys_kept = sum(s['kept'] for s in temps.values())
            sys_filtered = sum(s['filtered'] for s in temps.values())
            sys_rate = sys_kept / sys_total * 100 if sys_total > 0 else 0
            
            f.write("-"*80 + "\n")
            f.write(f"{'Total':<15} {sys_total:<10} {sys_kept:<10} "
                   f"{sys_filtered:<12} {sys_rate:>6.1f}%\n")
            f.write("\n")
        
        f.write("="*80 + "\n")
        f.write("Note: [Kept/Total runs] shown in each subplot title\n")
        f.write("      Red title = >30% filtered\n")
        f.write("="*80 + "\n")
    
    print(f"  [OK] Saved: {report_file.name}")


if __name__ == '__main__':
    # 解析命令行参数
    parser = argparse.ArgumentParser(
        description='绘制MSD曲线 (支持异常值筛选)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  python step3_plot_msd.py              # 默认启用筛选
  python step3_plot_msd.py --nofilter   # 关闭筛选，绘制所有曲线
        """
    )
    parser.add_argument(
        '--nofilter',
        action='store_true',
        help='关闭异常值筛选，绘制所有曲线（包括异常值）'
    )
    
    args = parser.parse_args()
    
    # 根据命令行参数决定是否筛选
    if args.nofilter:
        enable_filtering = False
    else:
        enable_filtering = None  # 使用全局配置
    
    main(enable_filtering)
