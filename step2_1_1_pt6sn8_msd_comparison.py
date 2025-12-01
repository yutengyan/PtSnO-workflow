#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 2.1.1: 负载型 Pt6Sn8 MSD 对比图 (300K vs 900K)

只绘制 PtSn 整体的 MSD 曲线:
- 300K: 蓝色
- 900K: 红色
- 多次模拟平均 + 误差带
- 导出绘图数据到 CSV
- 支持排除 Step1 检测的异常 run

Author: AI Assistant
Date: 2025-12-01
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from datetime import datetime
import re
import warnings
from collections import defaultdict
from scipy.signal import savgol_filter

warnings.filterwarnings('ignore')

# ===== 数据集选择 =====
# 可选: 'pt6sn8' (负载型Pt6Sn8) 或 'cv' (Pt6Sn8O4)
DATASET = 'pt6sn8'  # 🔧 修改此处切换数据集

# ===== 异常排除开关 =====
EXCLUDE_OUTLIERS = False  # 🔧 是否排除 Step1 检测的异常 run (基于 large_D_outliers.csv)
OUTLIER_CSV = Path(__file__).parent / 'results' / 'large_D_outliers.csv'

# ===== 内置 IQR 异常过滤 =====
APPLY_IQR_FILTER = False   # 🔧 是否在计算系综平均时应用 IQR 过滤 (推荐开启)
IQR_MULTIPLIER = 1.0      # IQR 倍数 (1.5 = 标准, 1.0 = 严格)

# 数据集配置
DATASET_CONFIGS = {
    'pt6sn8': {
        'name': '负载型 Pt6Sn8',
        'system_pattern': r'^pt6sn8',
        'output_subdir': 'step2.1.1_pt6sn8_msd',
        'output_prefix': 'pt6sn8_loaded',
        'title': r'Pt$_6$Sn$_8$ (Loaded)'
    },
    'cv': {
        'name': 'Pt6Sn8O4 (Cv)',
        'system_pattern': r'^Cv',
        'output_subdir': 'step2.1.1_cv_msd',
        'output_prefix': 'pt6sn8o4_cv',
        'title': r'Pt$_6$Sn$_8$O$_4$ (Cv)'
    }
}

# ===== 配置 =====
BASE_DIR = Path(__file__).parent
GMX_DATA_DIR = BASE_DIR / 'data' / 'gmx_msd' / 'unwrap' / 'gmx_msd_results_20251118_152614'

# 根据选择设置配置
current_config = DATASET_CONFIGS[DATASET]
OUTPUT_DIR = BASE_DIR / 'results' / current_config['output_subdir']
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

TARGET_SYSTEM = current_config['system_pattern']
TARGET_TEMPS = ['300K', '900K']
TARGET_ELEMENT = 'PtSn'  # 只绘制整体 MSD

# 绘图配置
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['font.size'] = 10

COLORS = {
    '300K': '#1E90FF',  # 蓝色
    '900K': '#DC143C'   # 红色
}


def load_outlier_list(csv_path, element='PtSn'):
    """
    加载异常 run 清单
    
    Parameters:
    -----------
    csv_path : Path
        large_D_outliers.csv 路径
    element : str
        目标元素 (Pt/Sn/PtSn)
    
    Returns:
    --------
    set : 异常文件路径集合 (使用标准化路径)
    """
    if not csv_path.exists():
        print(f"    ⚠️ 异常清单不存在: {csv_path}")
        return set()
    
    try:
        df = pd.read_csv(csv_path)
        # 筛选目标元素
        df_elem = df[df['element'] == element]
        # 标准化路径 (统一使用正斜杠)
        outlier_paths = set()
        for fp in df_elem['filepath']:
            # 标准化路径用于匹配
            normalized = Path(fp).as_posix()
            outlier_paths.add(normalized)
        return outlier_paths
    except Exception as e:
        print(f"    ⚠️ 加载异常清单失败: {e}")
        return set()


def scan_msd_files(data_dir, system_pattern, temps, element, exclude_outliers=False):
    """
    扫描 MSD 数据文件
    
    Returns:
    --------
    data_dict : dict
        {temperature: [file_path1, file_path2, ...]}
    """
    data_dir = Path(data_dir)
    data_dict = {temp: [] for temp in temps}
    
    # 加载异常清单
    outlier_paths = set()
    if exclude_outliers and OUTLIER_CSV.exists():
        outlier_paths = load_outlier_list(OUTLIER_CSV, element)
        print(f"\n>>> 加载异常清单: {len(outlier_paths)} 个异常 run")
    
    print(f"\n>>> 扫描 MSD 数据...")
    print(f"    目录: {data_dir.name}")
    print(f"    体系: {system_pattern}")
    print(f"    元素: {element}")
    if exclude_outliers:
        print(f"    排除异常: 已启用")
    
    excluded_count = 0
    
    for xvg_file in data_dir.rglob(f"*_msd_{element}.xvg"):
        try:
            # 检查是否在异常清单中
            if exclude_outliers and outlier_paths:
                file_posix = xvg_file.as_posix()
                # 检查路径是否匹配异常清单
                is_outlier = False
                for outlier_path in outlier_paths:
                    # 使用路径后缀匹配 (因为绝对路径可能不同)
                    if outlier_path.endswith(file_posix.split('gmx_msd_results_')[-1]) or \
                       file_posix.endswith(outlier_path.split('gmx_msd_results_')[-1]):
                        is_outlier = True
                        break
                if is_outlier:
                    excluded_count += 1
                    continue
            
            parts = xvg_file.parts
            
            # 提取 temperature 和 composition
            temperature = None
            composition = None
            for i in range(len(parts)-1, 0, -1):
                if parts[i].endswith('K'):
                    temperature = parts[i]
                    composition = parts[i-1]
                    break
            
            if not temperature or not composition:
                continue
            
            # 检查是否匹配目标体系
            if not re.match(system_pattern, composition, re.IGNORECASE):
                continue
            
            # 检查是否是目标温度
            if temperature not in temps:
                continue
            
            data_dict[temperature].append(xvg_file)
            
        except Exception as e:
            continue
    
    for temp in temps:
        print(f"    {temp}: {len(data_dict[temp])} 个文件")
    
    if exclude_outliers and excluded_count > 0:
        print(f"    已排除: {excluded_count} 个异常文件")
    
    return data_dict


def read_gmx_msd_xvg(filepath):
    """读取 GMX MSD .xvg 文件"""
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
                        t = float(parts[0])  # ps
                        msd_nm2 = float(parts[1])
                        msd_a2 = msd_nm2 * 100  # nm^2 -> Å^2
                        time_data.append(t)
                        msd_data.append(msd_a2)
                    except ValueError:
                        continue
    except:
        return None, None
    
    if len(time_data) == 0:
        return None, None
    
    return np.array(time_data), np.array(msd_data)


def detect_outliers_iqr(values, multiplier=1.5):
    """
    使用 IQR 方法检测异常值
    
    Parameters:
    -----------
    values : array-like
        数据值
    multiplier : float
        IQR 倍数 (默认 1.5)
    
    Returns:
    --------
    mask : ndarray
        True = 正常, False = 异常
    """
    values = np.array(values)
    Q1 = np.percentile(values, 25)
    Q3 = np.percentile(values, 75)
    IQR = Q3 - Q1
    lower = Q1 - multiplier * IQR
    upper = Q3 + multiplier * IQR
    return (values >= lower) & (values <= upper)


def compute_ensemble_msd(file_list, apply_iqr_filter=True, iqr_multiplier=1.5):
    """
    计算系综平均 MSD
    
    Parameters:
    -----------
    file_list : list
        MSD 文件路径列表
    apply_iqr_filter : bool
        是否应用 IQR 异常过滤
    iqr_multiplier : float
        IQR 倍数
    
    Returns:
    --------
    time : ndarray
        公共时间轴
    msd_mean : ndarray
        平均 MSD
    msd_std : ndarray
        标准差
    n_runs : int
        有效轨迹数
    """
    all_msd = []
    all_time = []
    file_paths = []
    
    for filepath in file_list:
        time, msd = read_gmx_msd_xvg(filepath)
        if time is not None:
            all_time.append(time)
            all_msd.append(msd)
            file_paths.append(filepath)
    
    if not all_msd:
        return None, None, None, 0
    
    # 找公共时间范围
    min_time = max(t.min() for t in all_time)
    max_time = min(t.max() for t in all_time)
    
    # 创建公共时间轴
    n_points = min(len(t) for t in all_time)
    common_time = np.linspace(min_time, max_time, n_points)
    
    # 插值到公共时间轴
    msd_interp = []
    for time, msd in zip(all_time, all_msd):
        msd_i = np.interp(common_time, time, msd)
        msd_interp.append(msd_i)
    
    msd_array = np.array(msd_interp)
    
    # ===== IQR 异常过滤 =====
    if apply_iqr_filter and len(msd_array) >= 3:
        # 使用 MSD 终值进行异常检测
        final_msd_values = msd_array[:, -1]
        normal_mask = detect_outliers_iqr(final_msd_values, iqr_multiplier)
        
        n_outliers = (~normal_mask).sum()
        if n_outliers > 0:
            # 显示被排除的文件
            outlier_indices = np.where(~normal_mask)[0]
            print(f"    [IQR过滤] 排除 {n_outliers} 条异常轨迹:")
            for idx in outlier_indices:
                fp = file_paths[idx]
                # 提取 Cv-x 目录名
                cv_match = re.search(r'(Cv-\d+)', str(fp))
                cv_name = cv_match.group(1) if cv_match else Path(fp).parent.parent.name
                print(f"      - {cv_name}: MSD终值 = {final_msd_values[idx]:.2f} Å²")
            
            # 只保留正常数据
            msd_array = msd_array[normal_mask]
    
    msd_mean = np.mean(msd_array, axis=0)
    msd_std = np.std(msd_array, axis=0, ddof=1) if len(msd_array) > 1 else np.zeros_like(msd_mean)
    
    return common_time, msd_mean, msd_std, len(msd_array)


def plot_msd_comparison(data_300k, data_900k, output_dir):
    """
    绘制 300K vs 900K MSD 对比图
    """
    output_dir = Path(output_dir)
    
    # 计算系综平均 (应用 IQR 过滤)
    print(f"\n>>> 计算系综平均...")
    if APPLY_IQR_FILTER:
        print(f"    IQR过滤: 已启用 (multiplier={IQR_MULTIPLIER})")
    
    time_300, msd_300, std_300, n_300 = compute_ensemble_msd(
        data_300k, apply_iqr_filter=APPLY_IQR_FILTER, iqr_multiplier=IQR_MULTIPLIER)
    time_900, msd_900, std_900, n_900 = compute_ensemble_msd(
        data_900k, apply_iqr_filter=APPLY_IQR_FILTER, iqr_multiplier=IQR_MULTIPLIER)
    
    if time_300 is None or time_900 is None:
        print("  ⚠️ 数据不足")
        return
    
    print(f"    300K: {n_300} 条轨迹, MSD范围 [{msd_300.min():.2f}, {msd_300.max():.2f}] Å²")
    print(f"    900K: {n_900} 条轨迹, MSD范围 [{msd_900.min():.2f}, {msd_900.max():.2f}] Å²")
    
    # 创建图表
    fig, ax = plt.subplots(figsize=(10/2.54, 8/2.54))
    
    # 绘制 300K (蓝色)
    ax.fill_between(time_300, msd_300 - std_300, msd_300 + std_300,
                    alpha=0.3, color=COLORS['300K'])
    ax.plot(time_300, msd_300, color=COLORS['300K'], linewidth=1.5,
            label='300K')
    
    # 绘制 900K (红色)
    ax.fill_between(time_900, msd_900 - std_900, msd_900 + std_900,
                    alpha=0.3, color=COLORS['900K'])
    ax.plot(time_900, msd_900, color=COLORS['900K'], linewidth=1.5,
            label='900K')
    
    # 标签和格式
    ax.set_xlabel('Time (ps)', fontsize=10, fontweight='bold')
    ax.set_ylabel(r'MSD ($\AA^2$)', fontsize=10, fontweight='bold')
    
    # 图例 - 小而紧凑
    ax.legend(loc='upper left', fontsize=8, framealpha=0.8,
              handletextpad=0.3, borderpad=0.3, labelspacing=0.2)
    ax.tick_params(axis='both', which='major', labelsize=9)
    ax.set_xlim(0, max(time_300.max(), time_900.max()))
    ax.set_ylim(0, None)
    
    plt.tight_layout()
    
    # 保存图片
    prefix = current_config['output_prefix']
    output_file = output_dir / f'{prefix}_msd_300K_vs_900K.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"\n📊 图片已保存: {output_file}")
    
    # 导出 CSV 数据
    export_msd_csv(time_300, msd_300, std_300, n_300,
                   time_900, msd_900, std_900, n_900, output_dir, prefix)
    
    return output_file


def export_msd_csv(time_300, msd_300, std_300, n_300,
                   time_900, msd_900, std_900, n_900, output_dir, prefix):
    """导出 MSD 数据到 CSV"""
    output_dir = Path(output_dir)
    
    # 导出 300K 数据
    df_300 = pd.DataFrame({
        'Time_ps': time_300,
        'MSD_mean_A2': msd_300,
        'MSD_std_A2': std_300,
        'N_runs': n_300
    })
    csv_300 = output_dir / f'{prefix}_msd_300K.csv'
    df_300.to_csv(csv_300, index=False)
    print(f"    CSV: {csv_300}")
    
    # 导出 900K 数据
    df_900 = pd.DataFrame({
        'Time_ps': time_900,
        'MSD_mean_A2': msd_900,
        'MSD_std_A2': std_900,
        'N_runs': n_900
    })
    csv_900 = output_dir / f'{prefix}_msd_900K.csv'
    df_900.to_csv(csv_900, index=False)
    print(f"    CSV: {csv_900}")
    
    # 导出统计汇总
    summary = pd.DataFrame({
        'Temperature': ['300K', '900K'],
        'N_runs': [n_300, n_900],
        'MSD_final_mean_A2': [msd_300[-1], msd_900[-1]],
        'MSD_final_std_A2': [std_300[-1], std_900[-1]],
        'MSD_max_A2': [msd_300.max(), msd_900.max()]
    })
    csv_summary = output_dir / f'{prefix}_msd_summary.csv'
    summary.to_csv(csv_summary, index=False)
    print(f"    汇总: {csv_summary}")
    
    # 打印统计
    print(f"\n>>> MSD 统计汇总:")
    print(f"    {'Temperature':<12} {'N_runs':<8} {'MSD_final (Å²)':<20} {'MSD_max (Å²)':<15}")
    print(f"    {'-'*55}")
    print(f"    {'300K':<12} {n_300:<8} {msd_300[-1]:.2f} ± {std_300[-1]:.2f} {'':<5} {msd_300.max():.2f}")
    print(f"    {'900K':<12} {n_900:<8} {msd_900[-1]:.2f} ± {std_900[-1]:.2f} {'':<5} {msd_900.max():.2f}")
    print(f"    {'-'*55}")
    if msd_300[-1] > 0:
        print(f"    900K/300K 比值: {msd_900[-1]/msd_300[-1]:.1f}x")


def main():
    print(f"\n{'='*60}")
    dataset_name = current_config['name']
    print(f"Step 2.1.1: {dataset_name} MSD 对比 (300K vs 900K)")
    print(f"当前数据集: {DATASET}")
    print(f"预排除异常(Step1): {'是' if EXCLUDE_OUTLIERS else '否'}")
    print(f"IQR实时过滤: {'是 (x' + str(IQR_MULTIPLIER) + ')' if APPLY_IQR_FILTER else '否'}")
    print(f"{'='*60}")
    
    # 扫描数据 (传入排除开关)
    data_dict = scan_msd_files(GMX_DATA_DIR, TARGET_SYSTEM, TARGET_TEMPS, TARGET_ELEMENT, 
                                exclude_outliers=EXCLUDE_OUTLIERS)
    
    if not data_dict['300K'] or not data_dict['900K']:
        print(f"\n[X] 错误: 数据不足")
        return
    
    # 绘图
    print(f"\n>>> 绘制 MSD 对比图...")
    plot_msd_comparison(data_dict['300K'], data_dict['900K'], OUTPUT_DIR)
    
    print(f"\n{'='*60}")
    print(f"✅ 完成!")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
