#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 6.1.1.3: Air68 vs Air86 分区热容组合图

将 Air68 和 Air86 的分区热容图并排显示，统一Y轴范围便于对比

用法:
  python step6_1_1_3_air_cv_combined.py
  python step6_1_1_3_air_cv_combined.py --figsize 10x8

作者: AI Assistant
日期: 2025-12-01
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress
from pathlib import Path

# 设置高质量论文图样式 - Arial (Nature/Science/ACS推荐)
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
plt.rcParams['mathtext.fontset'] = 'dejavusans'
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['xtick.major.width'] = 1.5
plt.rcParams['ytick.major.width'] = 1.5

# 字体大小常量
FONT_TICK = 28
FONT_LABEL = 34


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


def load_cluster_data(csv_path, exclude_dict=None):
    """加载聚类结果数据并过滤"""
    try:
        df = pd.read_csv(csv_path)
        if exclude_dict:
            print(f"  原始数据: {len(df)} 条")
            df = filter_data(df, exclude_dict)
            print(f"  过滤后: {len(df)} 条")
        return df
    except Exception as e:
        print(f"  错误: 无法读取 {csv_path}: {e}")
        return None


def compute_partition_data(df, structure_name, custom_partitions=None, peak_method='fit'):
    """计算分区热容数据
    
    参数:
        df: 聚类结果数据
        structure_name: 结构名称
        custom_partitions: 自定义分区列表，格式为 [(T_min1, T_max1), (T_min2, T_max2), ...]
        peak_method: 热容峰计算方法 ('data', 'partition', 'fit')
    """
    
    # 按温度分组计算能量
    temp_groups = df.groupby('temp')
    temps_unique = []
    E_mean = []
    E_std = []
    
    for temp, group in temp_groups:
        E_cluster = group['avg_energy'].values
        temps_unique.append(temp)
        E_mean.append(np.mean(E_cluster))
        E_std.append(np.std(E_cluster))
    
    temps_unique = np.array(temps_unique)
    E_mean = np.array(E_mean)
    E_std = np.array(E_std)
    
    # 相对能量
    E_ref = E_mean.min()
    E_rel = E_mean - E_ref
    
    # 多数投票确定每个温度的相态
    temp_to_partition = {}
    for temp in temps_unique:
        df_temp = df[df['temp'] == temp]
        partition_counts = df_temp['phase_clustered'].value_counts()
        temp_to_partition[temp] = partition_counts.idxmax()
    
    # 如果指定了自定义分区，使用自定义分区；否则使用聚类结果
    if custom_partitions is not None:
        # 使用自定义分区
        for i, (T_min, T_max) in enumerate(custom_partitions):
            for temp in temps_unique:
                if T_min <= temp <= T_max:
                    temp_to_partition[temp] = i
        phases = list(range(len(custom_partitions)))
    else:
        phases = sorted(df['phase_clustered'].unique())
    
    # 分区拟合
    phase_fits = {}
    
    for phase in phases:
        phase_temps = [t for t, p in temp_to_partition.items() if p == phase]
        phase_temps = sorted(phase_temps)
        
        if len(phase_temps) >= 2:
            mask = np.isin(temps_unique, phase_temps)
            T_phase = temps_unique[mask]
            E_phase = E_rel[mask]
            E_phase_std = E_std[mask]
            
            slope, intercept, r_value, _, std_err = linregress(T_phase, E_phase)
            
            phase_fits[phase] = {
                'slope': slope,
                'intercept': intercept,
                'Cv': slope * 1000,  # meV/K
                'Cv_err': std_err * 1000,
                'R2': r_value ** 2,
                'T_range': (T_phase.min(), T_phase.max()),
                'T_data': T_phase,
                'E_data': E_phase,
                'E_std': E_phase_std
            }
    
    # 分界温度和热容峰计算
    T_boundary = None
    Cv_peak = None
    peak_method_used = None
    
    if len(phases) >= 2:
        phases_sorted = sorted(phases)
        phase1_temps = [t for t, p in temp_to_partition.items() if p == phases_sorted[0]]
        phase2_temps = [t for t, p in temp_to_partition.items() if p == phases_sorted[1]]
        
        if phase1_temps and phase2_temps and phases_sorted[0] in phase_fits and phases_sorted[1] in phase_fits:
            T1_last = max(phase1_temps)
            T2_first = min(phase2_temps)
            T_boundary = (T1_last + T2_first) / 2
            
            fit1 = phase_fits[phases_sorted[0]]
            fit2 = phase_fits[phases_sorted[1]]
            Cv1 = fit1['Cv']
            Cv2 = fit2['Cv']
            
            # ========== 方法1a: 所有点平均 ==========
            idx1 = np.where(temps_unique == T1_last)[0]
            idx2 = np.where(temps_unique == T2_first)[0]
            
            if len(idx1) > 0 and len(idx2) > 0:
                E1_data = E_rel[idx1[0]]
                E2_data = E_rel[idx2[0]]
                Cv_transition_data = (E2_data - E1_data) / (T2_first - T1_last) * 1000  # meV/K
            else:
                Cv_transition_data = (Cv1 + Cv2) / 2
            
            # ========== 方法1b: 只用归属于该分区的点 ==========
            df_T1 = df[df['temp'] == T1_last]
            df_T2 = df[df['temp'] == T2_first]
            
            partition_T1 = temp_to_partition[T1_last]
            partition_T2 = temp_to_partition[T2_first]
            
            df_T1_filtered = df_T1[df_T1['phase_clustered'] == partition_T1]
            df_T2_filtered = df_T2[df_T2['phase_clustered'] == partition_T2]
            
            if len(df_T1_filtered) > 0 and len(df_T2_filtered) > 0:
                E1_partition = df_T1_filtered['avg_energy'].mean()
                E2_partition = df_T2_filtered['avg_energy'].mean()
                
                Cv_transition_partition = (E2_partition - E1_partition) / (T2_first - T1_last) * 1000
                
                n_T1_total = len(df_T1)
                n_T1_used = len(df_T1_filtered)
                n_T2_total = len(df_T2)
                n_T2_used = len(df_T2_filtered)
                print(f"  {structure_name} 分区点法: T1={T1_last}K 用{n_T1_used}/{n_T1_total}点({partition_T1}), "
                      f"T2={T2_first}K 用{n_T2_used}/{n_T2_total}点({partition_T2}), "
                      f"Cv={Cv_transition_partition:.2f} meV/K")
            else:
                Cv_transition_partition = Cv_transition_data
            
            # ========== 方法2: 拟合线外推 ==========
            E1_fit = fit1['slope'] * T1_last + fit1['intercept']
            E2_fit = fit2['slope'] * T2_first + fit2['intercept']
            Cv_transition_fit = (E2_fit - E1_fit) / (T2_first - T1_last) * 1000  # meV/K
            
            # 根据选择的方法选择热容峰
            if peak_method == 'data':
                Cv_peak = Cv_transition_data
                peak_method_used = "全点数据法"
            elif peak_method == 'partition':
                Cv_peak = Cv_transition_partition
                peak_method_used = "分区点法"
            else:  # 'fit'
                Cv_peak = Cv_transition_fit
                peak_method_used = "拟合线外推法"
    
    return {
        'temps': temps_unique,
        'E_rel': E_rel,
        'E_std': E_std,
        'temp_to_partition': temp_to_partition,
        'phase_fits': phase_fits,
        'T_boundary': T_boundary,
        'phases': phases,
        'Cv_peak': Cv_peak,
        'peak_method_used': peak_method_used
    }


def compute_unified_ylims(data_68, data_86):
    """计算统一的Y轴范围"""
    # 能量Y轴
    all_E = np.concatenate([data_68['E_rel'], data_86['E_rel']])
    all_E_std = np.concatenate([data_68['E_std'], data_86['E_std']])
    E_min = (all_E - all_E_std).min()
    E_max = (all_E + all_E_std).max()
    E_margin = (E_max - E_min) * 0.1
    E_ylim = (E_min - E_margin, E_max + E_margin)
    
    # 热容Y轴
    all_Cv = []
    for data in [data_68, data_86]:
        for fit in data['phase_fits'].values():
            all_Cv.append(fit['Cv'])
        # 检查是否有峰
        if data['Cv_peak'] is not None:
            all_Cv.append(data['Cv_peak'])
    
    if all_Cv:
        Cv_min = min(all_Cv) * 0.85
        Cv_max = max(all_Cv) * 1.1
    else:
        Cv_min, Cv_max = 0, 10
    
    Cv_ylim = (Cv_min, Cv_max)
    
    return E_ylim, Cv_ylim


def plot_single_partition(data, title, output_path, E_ylim, Cv_ylim, figsize=(10, 8)):
    """绘制单个分区热容图（独立子图）"""
    fig, ax1 = plt.subplots(figsize=figsize)
    
    temps = data['temps']
    E_rel = data['E_rel']
    E_std = data['E_std']
    phase_fits = data['phase_fits']
    phases = sorted(data['phases'])
    T_boundary = data['T_boundary']
    Cv_peak = data.get('Cv_peak')
    
    # 左Y轴: 能量数据点
    ax1.errorbar(temps, E_rel, yerr=E_std,
                 fmt='o', markersize=10, color='black',
                 ecolor='gray', elinewidth=2, capsize=4, capthick=2,
                 zorder=5)
    
    # 拟合线
    for phase in phases:
        if phase in phase_fits:
            fit = phase_fits[phase]
            T_fit = np.linspace(fit['T_range'][0], fit['T_range'][1], 50)
            E_fit = fit['slope'] * T_fit + fit['intercept']
            ax1.plot(T_fit, E_fit, '-', color='black', linewidth=2.5, zorder=4)
    
    # 连接分区
    if len(phases) >= 2 and phases[0] in phase_fits and phases[1] in phase_fits:
        fit1 = phase_fits[phases[0]]
        fit2 = phase_fits[phases[1]]
        T1_end = fit1['T_range'][1]
        T2_start = fit2['T_range'][0]
        idx1 = np.where(temps == T1_end)[0]
        idx2 = np.where(temps == T2_start)[0]
        if len(idx1) > 0 and len(idx2) > 0:
            E1 = E_rel[idx1[0]]
            E2 = E_rel[idx2[0]]
            ax1.plot([T1_end, T2_start], [E1, E2], '-', color='black', linewidth=2.5, zorder=4)
    
    ax1.set_xlabel('Temperature (K)', fontsize=FONT_LABEL)
    ax1.set_ylabel('Total Energy (eV)', fontsize=FONT_LABEL)
    ax1.set_ylim(E_ylim)
    ax1.tick_params(axis='both', labelsize=FONT_TICK)
    ax1.yaxis.set_major_locator(plt.MaxNLocator(5))  # 限制Y轴刻度数量
    
    # 右Y轴: 热容
    ax2 = ax1.twinx()
    
    if len(phases) >= 2 and phases[0] in phase_fits and phases[1] in phase_fits:
        Cv1 = phase_fits[phases[0]]['Cv']
        Cv2 = phase_fits[phases[1]]['Cv']
        
        # 使用从 compute_partition_data 计算的 Cv_peak
        if Cv_peak is not None:
            Cv_transition = Cv_peak
        else:
            # 备用: 如果 Cv_peak 为空，计算平均值
            Cv_transition = (Cv1 + Cv2) / 2
        
        has_peak = Cv_transition > max(Cv1, Cv2)
        
        if has_peak:
            # 带峰的热容曲线
            T_plot = np.linspace(temps.min(), temps.max(), 500)
            sigma = (T_boundary - temps.min() + temps.max() - T_boundary) / 4  # 自适应宽度
            Cv_plot = np.zeros_like(T_plot)
            
            for i, T in enumerate(T_plot):
                transition = 1 / (1 + np.exp(-(T - T_boundary) / (sigma * 0.5)))
                baseline = Cv1 + (Cv2 - Cv1) * transition
                gaussian = (Cv_transition - baseline) * np.exp(-0.5 * ((T - T_boundary) / sigma)**2)
                Cv_plot[i] = baseline + gaussian
            
            ax2.plot(T_plot, Cv_plot, 'r-', linewidth=2.5, zorder=3)
        else:
            # 阶梯形热容曲线
            ax2.plot([temps.min(), T_boundary], [Cv1, Cv1], 'r-', linewidth=2.5, zorder=3)
            ax2.plot([T_boundary, T_boundary], [Cv1, Cv2], 'r--', linewidth=2, zorder=3)
            ax2.plot([T_boundary, temps.max()], [Cv2, Cv2], 'r-', linewidth=2.5, zorder=3)
    else:
        if phase_fits:
            Cv_single = list(phase_fits.values())[0]['Cv']
            ax2.axhline(y=Cv_single, color='red', linewidth=2.5, zorder=3)
    
    ax2.set_ylabel(r'$C_v$ (meV/K)', fontsize=FONT_LABEL, color='red')
    ax2.tick_params(axis='y', labelcolor='red', labelsize=FONT_TICK, color='red')
    ax2.spines['right'].set_color('red')
    ax2.set_ylim(Cv_ylim)
    ax2.yaxis.set_major_locator(plt.MaxNLocator(nbins=5))
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    print(f"📊 已保存: {output_path}")


def plot_single_linear_fit(data, title, output_path, E_ylim, Cv_ylim, figsize=(10, 8)):
    """绘制单一线性拟合图（不分区，独立子图）"""
    fig, ax1 = plt.subplots(figsize=figsize)
    
    temps = data['temps']
    E_rel = data['E_rel']
    E_std = data['E_std']
    
    # 左Y轴: 能量数据点
    ax1.errorbar(temps, E_rel, yerr=E_std,
                 fmt='o', markersize=10, color='black',
                 ecolor='gray', elinewidth=2, capsize=4, capthick=2,
                 zorder=5)
    
    # 整体线性拟合（单一拟合线）
    slope, intercept, r_value, _, std_err = linregress(temps, E_rel)
    Cv_overall = slope * 1000  # meV/K
    Cv_err = std_err * 1000
    R2 = r_value ** 2
    
    # 绘制拟合线
    T_fit = np.linspace(temps.min(), temps.max(), 100)
    E_fit = slope * T_fit + intercept
    ax1.plot(T_fit, E_fit, '-', color='black', linewidth=2.5, zorder=4)
    
    ax1.set_xlabel('Temperature (K)', fontsize=FONT_LABEL)
    ax1.set_ylabel('Total Energy (eV)', fontsize=FONT_LABEL)
    ax1.set_ylim(E_ylim)
    ax1.tick_params(axis='both', labelsize=FONT_TICK)
    ax1.yaxis.set_major_locator(plt.MaxNLocator(nbins=5))
    
    # 右Y轴: 热容（单一水平线）
    ax2 = ax1.twinx()
    ax2.axhline(y=Cv_overall, color='red', linewidth=2.5, zorder=3)
    
    ax2.set_ylabel(r'$C_v$ (meV/K)', fontsize=FONT_LABEL, color='red')
    ax2.tick_params(axis='y', labelcolor='red', labelsize=FONT_TICK, color='red')
    ax2.spines['right'].set_color('red')
    ax2.set_ylim(Cv_ylim)
    ax2.yaxis.set_major_locator(plt.MaxNLocator(nbins=5))
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    print(f"📊 已保存: {output_path}")
    print(f"    单一拟合: Cv={Cv_overall:.2f}±{Cv_err:.2f} meV/K, R²={R2:.4f}")
    
    return Cv_overall


def plot_combined_cv(data_68, data_86, output_dir, figsize=(10, 8)):
    """绘制独立子图（替代原来的组合图）"""
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # ========== 计算统一Y轴范围 ==========
    E_ylim, Cv_ylim = compute_unified_ylims(data_68, data_86)
    print(f"\n  统一Y轴范围:")
    print(f"    能量: {E_ylim[0]:.3f} ~ {E_ylim[1]:.3f} eV")
    print(f"    热容: {Cv_ylim[0]:.1f} ~ {Cv_ylim[1]:.1f} meV/K")
    
    # ========== 绘制三个独立子图 ==========
    # 子图1: Pt8Sn6 (Air86) - 分区拟合
    print(f"\n>>> 绘制 Pt8Sn6 分区拟合图...")
    plot_single_partition(data_86, r'Pt$_8$Sn$_6$', 
                          output_dir / 'Air86_Pt8Sn6_partition_cv.png',
                          E_ylim, Cv_ylim, figsize)
    
    # 子图2: Pt6Sn8 (Air68) - 分区拟合
    print(f"\n>>> 绘制 Pt6Sn8 分区拟合图...")
    plot_single_partition(data_68, r'Pt$_6$Sn$_8$ (partition)', 
                          output_dir / 'Air68_Pt6Sn8_partition_cv.png',
                          E_ylim, Cv_ylim, figsize)
    
    # 子图3: Pt6Sn8 (Air68) - 单一拟合
    print(f"\n>>> 绘制 Pt6Sn8 单一拟合图...")
    Cv_68_single = plot_single_linear_fit(data_68, r'Pt$_6$Sn$_8$ (single fit)', 
                                          output_dir / 'Air68_Pt6Sn8_single_fit_cv.png',
                                          E_ylim, Cv_ylim, figsize)
    
    return [
        output_dir / 'Air86_Pt8Sn6_partition_cv.png',
        output_dir / 'Air68_Pt6Sn8_partition_cv.png',
        output_dir / 'Air68_Pt6Sn8_single_fit_cv.png'
    ]


def interactive_adjust_plot(data_68, data_86, output_dir, figsize):
    """交互式调整图片参数"""
    from matplotlib.ticker import MaxNLocator, MultipleLocator
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 计算统一Y轴范围
    E_ylim, Cv_ylim = compute_unified_ylims(data_68, data_86)
    
    # 默认参数
    params = {
        'figsize': figsize,
        'y_nticks': 5,           # Y轴刻度数量
        'y_integer': True,       # Y轴使用整数
        'cv_nticks': 5,          # Cv轴刻度数量
        'cv_integer': True,      # Cv轴使用整数
    }
    
    plt.ion()  # 开启交互模式
    
    def create_preview_figure():
        """创建预览图"""
        fig, axes = plt.subplots(1, 3, figsize=(params['figsize'][0]*3, params['figsize'][1]))
        
        # 绘制三个子图预览
        for idx, (data, title) in enumerate([
            (data_86, r'Pt$_8$Sn$_6$'),
            (data_68, r'Pt$_6$Sn$_8$ (partition)'),
            (data_68, r'Pt$_6$Sn$_8$ (single fit)')
        ]):
            ax1 = axes[idx]
            temps = data['temps']
            E_rel = data['E_rel']
            E_std = data['E_std']
            
            ax1.errorbar(temps, E_rel, yerr=E_std, fmt='o', markersize=8, color='black',
                        ecolor='gray', elinewidth=1.5, capsize=3, capthick=1.5, zorder=5)
            
            ax1.set_xlabel('Temperature (K)', fontsize=FONT_LABEL)
            ax1.set_ylabel('Total Energy (eV)', fontsize=FONT_LABEL)
            ax1.set_ylim(E_ylim)
            ax1.tick_params(axis='both', labelsize=FONT_TICK)
            ax1.set_title(title, fontsize=FONT_LABEL)
            
            # 应用Y轴刻度设置
            if params['y_integer']:
                ax1.yaxis.set_major_locator(MaxNLocator(nbins=params['y_nticks'], integer=True))
            else:
                ax1.yaxis.set_major_locator(MaxNLocator(nbins=params['y_nticks']))
            
            # 右Y轴
            ax2 = ax1.twinx()
            ax2.axhline(y=50, color='red', linewidth=2, zorder=3)  # 示意线
            ax2.set_ylabel(r'$C_v$ (meV/K)', fontsize=FONT_LABEL, color='red')
            ax2.tick_params(axis='y', labelcolor='red', labelsize=FONT_TICK, color='red')
            ax2.spines['right'].set_color('red')
            ax2.set_ylim(Cv_ylim)
            
            # 应用Cv轴刻度设置
            if params['cv_integer']:
                ax2.yaxis.set_major_locator(MaxNLocator(nbins=params['cv_nticks'], integer=True))
            else:
                ax2.yaxis.set_major_locator(MaxNLocator(nbins=params['cv_nticks']))
        
        plt.tight_layout()
        return fig
    
    fig = create_preview_figure()
    plt.show(block=False)
    
    print("\n" + "="*60)
    print("🎨 交互式调整模式")
    print("="*60)
    
    while True:
        print(f"\n当前参数:")
        print(f"  [1] figsize: {params['figsize'][0]}x{params['figsize'][1]}")
        print(f"  [2] 能量Y轴刻度数: {params['y_nticks']} (整数: {'是' if params['y_integer'] else '否'})")
        print(f"  [3] Cv轴刻度数: {params['cv_nticks']} (整数: {'是' if params['cv_integer'] else '否'})")
        print(f"\n命令: 输入数字修改参数, 'r'刷新预览, 's'保存并退出, 'q'不保存退出")
        
        cmd = input(">>> ").strip().lower()
        
        if cmd == 'q':
            plt.close(fig)
            print("已取消")
            return None
        
        elif cmd == 's':
            plt.close(fig)
            print("\n正在保存最终图片...")
            # 使用调整后的参数保存
            plot_combined_cv_with_params(data_68, data_86, output_dir, params)
            return params
        
        elif cmd == 'r':
            plt.close(fig)
            fig = create_preview_figure()
            plt.show(block=False)
            print("✅ 预览已刷新")
        
        elif cmd == '1':
            val = input("  输入新figsize (格式 宽x高, 如 10x8): ").strip()
            try:
                w, h = map(float, val.lower().split('x'))
                params['figsize'] = (w, h)
                print(f"  ✅ figsize 设为 {w}x{h}")
            except:
                print("  ❌ 格式错误")
        
        elif cmd == '2':
            val = input("  输入能量Y轴刻度数 (如 5): ").strip()
            try:
                params['y_nticks'] = int(val)
                print(f"  ✅ 能量Y轴刻度数 设为 {params['y_nticks']}")
            except:
                print("  ❌ 格式错误")
            
            val2 = input("  使用整数刻度? (y/n, 默认y): ").strip().lower()
            params['y_integer'] = val2 != 'n'
            print(f"  ✅ 整数刻度: {'是' if params['y_integer'] else '否'}")
        
        elif cmd == '3':
            val = input("  输入Cv轴刻度数 (如 5): ").strip()
            try:
                params['cv_nticks'] = int(val)
                print(f"  ✅ Cv轴刻度数 设为 {params['cv_nticks']}")
            except:
                print("  ❌ 格式错误")
            
            val2 = input("  使用整数刻度? (y/n, 默认y): ").strip().lower()
            params['cv_integer'] = val2 != 'n'
            print(f"  ✅ 整数刻度: {'是' if params['cv_integer'] else '否'}")
        
        else:
            print("  未知命令，请重试")


def plot_combined_cv_with_params(data_68, data_86, output_dir, params):
    """使用交互参数绘制并保存图片"""
    from matplotlib.ticker import MaxNLocator
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    E_ylim, Cv_ylim = compute_unified_ylims(data_68, data_86)
    
    figsize = params['figsize']
    y_nticks = params['y_nticks']
    y_integer = params['y_integer']
    cv_nticks = params['cv_nticks']
    cv_integer = params['cv_integer']
    y_ticks_custom = params.get('y_ticks_custom', None)
    cv_ticks_custom = params.get('cv_ticks_custom', None)
    show_error_bars = params.get('show_error_bars', False)
    transparent = params.get('transparent', False)
    
    # 子图1: Pt8Sn6 (Air86)
    plot_single_partition_with_params(data_86, r'Pt$_8$Sn$_6$', 
                                      output_dir / 'Air86_Pt8Sn6_partition_cv.png',
                                      E_ylim, Cv_ylim, figsize, y_nticks, y_integer, 
                                      cv_nticks, cv_integer, y_ticks_custom, cv_ticks_custom,
                                      show_error_bars, transparent)
    
    # 子图2: Pt6Sn8 (Air68) - 分区
    plot_single_partition_with_params(data_68, r'Pt$_6$Sn$_8$ (partition)', 
                                      output_dir / 'Air68_Pt6Sn8_partition_cv.png',
                                      E_ylim, Cv_ylim, figsize, y_nticks, y_integer, 
                                      cv_nticks, cv_integer, y_ticks_custom, cv_ticks_custom,
                                      show_error_bars, transparent)
    
    # 子图3: Pt6Sn8 (Air68) - 单一拟合
    plot_single_linear_fit_with_params(data_68, r'Pt$_6$Sn$_8$ (single fit)', 
                                       output_dir / 'Air68_Pt6Sn8_single_fit_cv.png',
                                       E_ylim, Cv_ylim, figsize, y_nticks, y_integer, 
                                       cv_nticks, cv_integer, y_ticks_custom, cv_ticks_custom,
                                       show_error_bars, transparent)


def plot_combined_cv_with_params_three_systems(data_68, data_86, data_sup86, output_dir, params):
    """三系统版本：使用交互参数绘制并保存图片"""
    from matplotlib.ticker import MaxNLocator
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 计算统一Y轴范围（包含三个系统）
    all_E = np.concatenate([data_68['E_rel'], data_86['E_rel'], data_sup86['E_rel']])
    all_E_std = np.concatenate([data_68['E_std'], data_86['E_std'], data_sup86['E_std']])
    E_min = (all_E - all_E_std).min()
    E_max = (all_E + all_E_std).max()
    E_margin = (E_max - E_min) * 0.1
    E_ylim = (E_min - E_margin, E_max + E_margin)
    
    all_Cv = []
    for data in [data_68, data_86, data_sup86]:
        for phase, fit in data['phase_fits'].items():
            all_Cv.append(fit['Cv'])
    
    if all_Cv:
        Cv_min = min(all_Cv) * 0.9
        Cv_max = max(all_Cv) * 1.1
    else:
        Cv_min, Cv_max = 0, 10
    
    Cv_ylim = (Cv_min, Cv_max)
    
    figsize = params['figsize']
    y_nticks = params['y_nticks']
    y_integer = params['y_integer']
    cv_nticks = params['cv_nticks']
    cv_integer = params['cv_integer']
    y_ticks_custom = params.get('y_ticks_custom', None)
    cv_ticks_custom = params.get('cv_ticks_custom', None)
    show_error_bars = params.get('show_error_bars', False)
    transparent = params.get('transparent', False)
    
    # 子图1: Pt8Sn6 (Air86)
    plot_single_partition_with_params(data_86, r'Pt$_8$Sn$_6$ (Air)', 
                                      output_dir / 'Air86_Pt8Sn6_partition_cv.png',
                                      E_ylim, Cv_ylim, figsize, y_nticks, y_integer, 
                                      cv_nticks, cv_integer, y_ticks_custom, cv_ticks_custom,
                                      show_error_bars, transparent)
    
    # 子图2: Pt6Sn8 (Air68)
    plot_single_partition_with_params(data_68, r'Pt$_6$Sn$_8$ (Air)', 
                                      output_dir / 'Air68_Pt6Sn8_partition_cv.png',
                                      E_ylim, Cv_ylim, figsize, y_nticks, y_integer, 
                                      cv_nticks, cv_integer, y_ticks_custom, cv_ticks_custom,
                                      show_error_bars, transparent)
    
    # 子图3: Pt8Sn6 (sup86) - 负载型
    plot_single_partition_with_params(data_sup86, r'Pt$_8$Sn$_6$ (support)', 
                                      output_dir / 'sup86_Pt8Sn6_partition_cv.png',
                                      E_ylim, Cv_ylim, figsize, y_nticks, y_integer, 
                                      cv_nticks, cv_integer, y_ticks_custom, cv_ticks_custom,
                                      show_error_bars, transparent)


def plot_single_partition_with_params(data, title, output_path, E_ylim, Cv_ylim, 
                                      figsize, y_nticks, y_integer, cv_nticks, cv_integer,
                                      y_ticks_custom=None, cv_ticks_custom=None, show_error_bars=False, transparent=False):
    """带参数的分区热容图绘制"""
    from matplotlib.ticker import MaxNLocator, LinearLocator, MultipleLocator
    
    fig, ax1 = plt.subplots(figsize=figsize)
    
    temps = data['temps']
    E_rel = data['E_rel']
    E_std = data['E_std']
    phase_fits = data['phase_fits']
    phases = sorted(data['phases'])
    T_boundary = data['T_boundary']
    
    # 根据 show_error_bars 参数选择绘制方式
    if show_error_bars:
        # 显示误差棒模式：黑色实心大点 + 误差棒，不显示散点
        ax1.errorbar(temps, E_rel, yerr=E_std, fmt='o', markersize=14, 
                     markerfacecolor='black', markeredgecolor='black',
                     ecolor='black', elinewidth=2.5, capsize=6, capthick=2.5, zorder=5)
    else:
        # 默认模式：空心点 + 灰色误差棒
        ax1.errorbar(temps, E_rel, yerr=E_std, fmt='o', markersize=10, color='black',
                     ecolor='gray', elinewidth=2, capsize=4, capthick=2, zorder=5)
    
    for phase in phases:
        if phase in phase_fits:
            fit = phase_fits[phase]
            T_min, T_max = fit['T_range']
            mask = (temps >= T_min) & (temps <= T_max)
            T_fit = np.linspace(T_min, T_max, 100)
            E_fit = fit['intercept'] + (fit['Cv']/1000) * T_fit
            ax1.plot(T_fit, E_fit, '-', color='black', linewidth=2.5, zorder=4)
    
    if len(phases) >= 2 and phases[0] in phase_fits and phases[1] in phase_fits:
        fit1 = phase_fits[phases[0]]
        fit2 = phase_fits[phases[1]]
        T1_end = fit1['T_range'][1]
        T2_start = fit2['T_range'][0]
        idx1 = np.where(temps == T1_end)[0]
        idx2 = np.where(temps == T2_start)[0]
        if len(idx1) > 0 and len(idx2) > 0:
            E1 = E_rel[idx1[0]]
            E2 = E_rel[idx2[0]]
            ax1.plot([T1_end, T2_start], [E1, E2], '-', color='black', linewidth=2.5, zorder=4)
    
    ax1.set_xlabel('Temperature (K)', fontsize=FONT_LABEL)
    ax1.set_ylabel('Total Energy (eV)', fontsize=FONT_LABEL)
    ax1.set_ylim(E_ylim)
    ax1.tick_params(axis='both', labelsize=FONT_TICK)
    
    # 设置Y轴刻度
    if y_ticks_custom is not None:
        # 使用自定义刻度
        ax1.set_yticks(y_ticks_custom)
    elif y_integer:
        # 整数刻度：生成整数刻度
        y_ticks = np.linspace(E_ylim[0], E_ylim[1], y_nticks)
        y_ticks = np.round(y_ticks)
        ax1.set_yticks(y_ticks)
    else:
        y_ticks = np.linspace(E_ylim[0], E_ylim[1], y_nticks)
        ax1.set_yticks(y_ticks)
    
    # 右Y轴
    ax2 = ax1.twinx()
    
    if len(phases) >= 2 and phases[0] in phase_fits and phases[1] in phase_fits:
        Cv1 = phase_fits[phases[0]]['Cv']
        Cv2 = phase_fits[phases[1]]['Cv']
        fit1 = phase_fits[phases[0]]
        fit2 = phase_fits[phases[1]]
        T1_last = fit1['T_range'][1]
        T2_first = fit2['T_range'][0]
        idx1 = np.where(temps == T1_last)[0]
        idx2 = np.where(temps == T2_first)[0]
        
        if len(idx1) > 0 and len(idx2) > 0:
            E1 = E_rel[idx1[0]]
            E2 = E_rel[idx2[0]]
            Cv_transition = (E2 - E1) / (T2_first - T1_last) * 1000
        else:
            Cv_transition = (Cv1 + Cv2) / 2
        
        has_peak = Cv_transition > max(Cv1, Cv2)
        
        if has_peak:
            T_plot = np.linspace(temps.min(), temps.max(), 500)
            sigma = (T2_first - T1_last) / 2
            Cv_plot = np.zeros_like(T_plot)
            for i, T in enumerate(T_plot):
                transition = 1 / (1 + np.exp(-(T - T_boundary) / (sigma * 0.5)))
                baseline = Cv1 + (Cv2 - Cv1) * transition
                gaussian = (Cv_transition - baseline) * np.exp(-0.5 * ((T - T_boundary) / sigma)**2)
                Cv_plot[i] = baseline + gaussian
            ax2.plot(T_plot, Cv_plot, 'r-', linewidth=2.5, zorder=3)
        else:
            ax2.plot([temps.min(), T_boundary], [Cv1, Cv1], 'r-', linewidth=2.5, zorder=3)
            ax2.plot([T_boundary, T_boundary], [Cv1, Cv2], 'r--', linewidth=2, zorder=3)
            ax2.plot([T_boundary, temps.max()], [Cv2, Cv2], 'r-', linewidth=2.5, zorder=3)
    else:
        Cv_single = list(phase_fits.values())[0]['Cv']
        ax2.axhline(y=Cv_single, color='red', linewidth=2.5, zorder=3)
    
    ax2.set_ylabel(r'$C_v$ (meV/K)', fontsize=FONT_LABEL, color='red')
    ax2.tick_params(axis='y', labelcolor='red', labelsize=FONT_TICK, color='red')
    ax2.spines['right'].set_color('red')
    ax2.set_ylim(Cv_ylim)
    
    # 设置Cv轴刻度
    if cv_ticks_custom is not None:
        # 使用自定义刻度
        ax2.set_yticks(cv_ticks_custom)
    elif cv_integer:
        cv_ticks = np.linspace(Cv_ylim[0], Cv_ylim[1], cv_nticks)
        cv_ticks = np.round(cv_ticks)
        ax2.set_yticks(cv_ticks)
    else:
        cv_ticks = np.linspace(Cv_ylim[0], Cv_ylim[1], cv_nticks)
        ax2.set_yticks(cv_ticks)
    
    plt.tight_layout()
    if transparent:
        plt.savefig(output_path, dpi=300, bbox_inches='tight', transparent=True)
    else:
        plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"📊 已保存: {output_path}")


def plot_single_linear_fit_with_params(data, title, output_path, E_ylim, Cv_ylim,
                                       figsize, y_nticks, y_integer, cv_nticks, cv_integer,
                                       y_ticks_custom=None, cv_ticks_custom=None, show_error_bars=False, transparent=False):
    """带参数的单一拟合图绘制"""
    from matplotlib.ticker import MaxNLocator, LinearLocator, MultipleLocator
    
    fig, ax1 = plt.subplots(figsize=figsize)
    
    temps = data['temps']
    E_rel = data['E_rel']
    E_std = data['E_std']
    
    # 根据 show_error_bars 参数选择绘制方式
    if show_error_bars:
        # 显示误差棒模式：黑色实心大点 + 误差棒
        ax1.errorbar(temps, E_rel, yerr=E_std, fmt='o', markersize=14,
                     markerfacecolor='black', markeredgecolor='black',
                     ecolor='black', elinewidth=2.5, capsize=6, capthick=2.5, zorder=5)
    else:
        # 默认模式：空心点 + 灰色误差棒
        ax1.errorbar(temps, E_rel, yerr=E_std, fmt='o', markersize=10, color='black',
                     ecolor='gray', elinewidth=2, capsize=4, capthick=2, zorder=5)
    
    slope, intercept, r_value, _, std_err = linregress(temps, E_rel)
    Cv_overall = slope * 1000
    Cv_err = std_err * 1000
    
    T_fit = np.linspace(temps.min(), temps.max(), 100)
    E_fit = slope * T_fit + intercept
    ax1.plot(T_fit, E_fit, '-', color='black', linewidth=2.5, zorder=4)
    
    ax1.set_xlabel('Temperature (K)', fontsize=FONT_LABEL)
    ax1.set_ylabel('Total Energy (eV)', fontsize=FONT_LABEL)
    ax1.set_ylim(E_ylim)
    ax1.tick_params(axis='both', labelsize=FONT_TICK)
    
    # 设置Y轴刻度
    if y_ticks_custom is not None:
        ax1.set_yticks(y_ticks_custom)
    elif y_integer:
        y_ticks = np.linspace(E_ylim[0], E_ylim[1], y_nticks)
        y_ticks = np.round(y_ticks)
        ax1.set_yticks(y_ticks)
    else:
        y_ticks = np.linspace(E_ylim[0], E_ylim[1], y_nticks)
        ax1.set_yticks(y_ticks)
    
    ax2 = ax1.twinx()
    ax2.axhline(y=Cv_overall, color='red', linewidth=2.5, zorder=3)
    ax2.set_ylabel(r'$C_v$ (meV/K)', fontsize=FONT_LABEL, color='red')
    ax2.tick_params(axis='y', labelcolor='red', labelsize=FONT_TICK, color='red')
    ax2.spines['right'].set_color('red')
    ax2.set_ylim(Cv_ylim)
    
    # 设置Cv轴刻度
    if cv_ticks_custom is not None:
        ax2.set_yticks(cv_ticks_custom)
    elif cv_integer:
        cv_ticks = np.linspace(Cv_ylim[0], Cv_ylim[1], cv_nticks)
        cv_ticks = np.round(cv_ticks)
        ax2.set_yticks(cv_ticks)
    else:
        cv_ticks = np.linspace(Cv_ylim[0], Cv_ylim[1], cv_nticks)
        ax2.set_yticks(cv_ticks)
    
    plt.tight_layout()
    if transparent:
        plt.savefig(output_path, dpi=300, bbox_inches='tight', transparent=True)
    else:
        plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"📊 已保存: {output_path}")


def main():
    # 解析命令行参数
    parser = argparse.ArgumentParser(
        description='Air68 vs Air86 分区热容独立子图',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例:
  python step6_1_1_3_air_cv_combined.py --y-nticks 3 --cv-nticks 4
  python step6_1_1_3_air_cv_combined.py --y-ticks 0,2,4 --cv-ticks 3,4,5,6,7
  python step6_1_1_3_air_cv_combined.py --figsize 12x10 --y-ticks 0,2,4
  python step6_1_1_3_air_cv_combined.py --peak-method partition --partitions-68 200-500,600-1100 --partitions-86 200-700,800-1100
  python step6_1_1_3_air_cv_combined.py --peak-method partition --partitions-both 200-400,500-1100
'''
    )
    parser.add_argument('--figsize', type=str, default='10x8',
                       help='图片尺寸，格式: 宽x高，例如 10x8 (默认: 10x8)')
    parser.add_argument('--interactive', '-i', action='store_true',
                       help='开启交互模式，可以实时调整Y轴刻度等参数')
    parser.add_argument('--y-nticks', type=int, default=5,
                       help='能量Y轴刻度数量 (默认: 5)，如果指定了 --y-ticks 则忽略')
    parser.add_argument('--cv-nticks', type=int, default=5,
                       help='Cv轴刻度数量 (默认: 5)，如果指定了 --cv-ticks 则忽略')
    parser.add_argument('--y-ticks', type=str, default=None,
                       help='手动指定能量Y轴刻度，逗号分隔，例如: 0,2,4')
    parser.add_argument('--cv-ticks', type=str, default=None,
                       help='手动指定Cv轴刻度，逗号分隔，例如: 3,4,5,6,7')
    parser.add_argument('--no-integer', action='store_true',
                       help='Y轴不使用整数刻度（仅在未指定 --y-ticks 时有效）')
    parser.add_argument('--partitions-68', '-p68', type=str, default=None,
                       help='Air68 分区，格式: T_min1-T_max1,T_min2-T_max2，例如: 200-500,600-1100')
    parser.add_argument('--partitions-86', '-p86', type=str, default=None,
                       help='Air86 分区，格式: T_min1-T_max1,T_min2-T_max2，例如: 200-700,800-1100')
    parser.add_argument('--partitions-both', '-pb', type=str, default=None,
                       help='同时应用到两个系统的分区，格式: T_min1-T_max1,T_min2-T_max2，例如: 200-400,500-1100')
    parser.add_argument('--peak-method', type=str, default='fit',
                       choices=['data', 'partition', 'fit'],
                       help='热容峰计算方法: data(全点法), partition(分区点法★推荐), fit(拟合法) (默认: fit)')
    parser.add_argument('--exclude-68', nargs='+', metavar='TEMP:INDICES',
                       help='Air68 要排除的点，格式: "300K:0,1" "400K:0"')
    parser.add_argument('--exclude-86', nargs='+', metavar='TEMP:INDICES',
                       help='Air86 要排除的点，格式: "500K:0,1" "600K:0"')
    parser.add_argument('--add-sup86', action='store_true',
                       help='添加 sup86 (负载型 Pt8Sn6) 数据对比')
    parser.add_argument('--partitions-sup86', '-psup', type=str, default=None,
                       help='sup86 分区，格式: T_min1-T_max1,T_min2-T_max2')
    parser.add_argument('--exclude-sup86', nargs='+', metavar='TEMP:INDICES',
                       help='sup86 要排除的点，格式: "500K:0,1" "600K:0"')
    parser.add_argument('--show-error-bars', action='store_true',
                       help='显示误差棒（不显示散点，用黑色实心大点+误差棒替代）')
    parser.add_argument('--transparent', action='store_true',
                       help='保存为透明背景图片')
    args = parser.parse_args()
    
    # 解析figsize
    try:
        fig_w, fig_h = map(float, args.figsize.lower().split('x'))
        figsize = (fig_w, fig_h)
    except ValueError:
        print(f"警告: 无效的figsize格式 '{args.figsize}'，使用默认 10x8")
        figsize = (10, 8)
    
    # 解析分区参数：支持 --partitions-68, --partitions-86, 或 --partitions-both
    custom_partitions_68 = None
    custom_partitions_86 = None
    
    # 优先级: --partitions-68/86 > --partitions-both
    if args.partitions_68 or args.partitions_86:
        # 分别指定两个系统
        if args.partitions_68:
            try:
                custom_partitions_68 = []
                for part in args.partitions_68.split(','):
                    T_min, T_max = map(int, part.split('-'))
                    custom_partitions_68.append((T_min, T_max))
                print(f"  Air68 分区: {custom_partitions_68}")
            except ValueError:
                print(f"警告: 无效的 --partitions-68 格式 '{args.partitions_68}'，将使用聚类结果")
                custom_partitions_68 = None
        
        if args.partitions_86:
            try:
                custom_partitions_86 = []
                for part in args.partitions_86.split(','):
                    T_min, T_max = map(int, part.split('-'))
                    custom_partitions_86.append((T_min, T_max))
                print(f"  Air86 分区: {custom_partitions_86}")
            except ValueError:
                print(f"警告: 无效的 --partitions-86 格式 '{args.partitions_86}'，将使用聚类结果")
                custom_partitions_86 = None
    
    elif args.partitions_both:
        # 同时应用到两个系统
        try:
            partitions_both = []
            for part in args.partitions_both.split(','):
                T_min, T_max = map(int, part.split('-'))
                partitions_both.append((T_min, T_max))
            custom_partitions_68 = partitions_both
            custom_partitions_86 = partitions_both
            print(f"  两个系统分区: {partitions_both}")
        except ValueError:
            print(f"警告: 无效的 --partitions-both 格式 '{args.partitions_both}'，将使用聚类结果")
            custom_partitions_68 = None
            custom_partitions_86 = None
    
    # 解析排除点参数
    exclude_68 = parse_exclude_points(args.exclude_68)
    exclude_86 = parse_exclude_points(args.exclude_86)
    
    if exclude_68:
        print(f"  Air68 排除点: {exclude_68}")
    if exclude_86:
        print(f"  Air86 排除点: {exclude_86}")
    
    # 解析自定义刻度
    y_ticks_custom = None
    cv_ticks_custom = None
    
    if args.y_ticks:
        try:
            y_ticks_custom = [float(x.strip()) for x in args.y_ticks.split(',')]
            print(f"  能量Y轴刻度: {y_ticks_custom}")
        except ValueError:
            print(f"警告: 无效的 --y-ticks 格式 '{args.y_ticks}'，将自动计算")
    
    if args.cv_ticks:
        try:
            cv_ticks_custom = [float(x.strip()) for x in args.cv_ticks.split(',')]
            print(f"  Cv轴刻度: {cv_ticks_custom}")
        except ValueError:
            print(f"警告: 无效的 --cv-ticks 格式 '{args.cv_ticks}'，将自动计算")
    
    system_names = "Air68 vs Air86"
    if args.add_sup86:
        system_names = "Air68 vs Air86 vs sup86"
    
    print("=" * 60)
    print(f"Step 6.1.1.3: {system_names} 分区热容独立子图")
    print("=" * 60)
    print(f"  图片尺寸: {figsize[0]}x{figsize[1]}")
    print(f"  热容峰方法: {args.peak_method}")
    
    # 加载数据
    base_dir = Path('results/step6_1_clustering')
    
    csv_68 = base_dir / 'Air68_kmeans_n2_clustered_data.csv'
    csv_86 = base_dir / 'Air86_kmeans_n2_clustered_data.csv'
    csv_sup86 = base_dir / 'Pt8sn6_kmeans_n2_clustered_data.csv' if args.add_sup86 else None
    
    if not csv_68.exists():
        print(f"错误: 找不到 {csv_68}")
        return
    if not csv_86.exists():
        print(f"错误: 找不到 {csv_86}")
        return
    if args.add_sup86 and not csv_sup86.exists():
        print(f"错误: 找不到 {csv_sup86}")
        print("  请先运行: python step6_1_clustering_analysis.py --structure Pt8sn6")
        return
    
    print(f"\n>>> 加载数据...")
    print("  Air68:")
    df_68 = load_cluster_data(csv_68, exclude_68)
    print("  Air86:")
    df_86 = load_cluster_data(csv_86, exclude_86)
    
    df_sup86 = None
    if args.add_sup86:
        print("  sup86:")
        exclude_sup86 = parse_exclude_points(args.exclude_sup86)
        custom_partitions_sup86 = parse_partitions(args.partitions_sup86) if args.partitions_sup86 else None
        if args.partitions_both and custom_partitions_sup86 is None:
            custom_partitions_sup86 = parse_partitions(args.partitions_both)
        df_sup86 = load_cluster_data(csv_sup86, exclude_sup86)
    
    if df_68 is None or df_86 is None or (args.add_sup86 and df_sup86 is None):
        return
    
    print(f"\n  最终数据:")
    print(f"    Air68: {len(df_68)} 条记录")
    print(f"    Air86: {len(df_86)} 条记录")
    if args.add_sup86:
        print(f"    sup86: {len(df_sup86)} 条记录")
    
    # 计算分区数据
    print(f"\n>>> 计算分区热容...")
    data_68 = compute_partition_data(df_68, 'Air68', custom_partitions=custom_partitions_68, peak_method=args.peak_method)
    data_86 = compute_partition_data(df_86, 'Air86', custom_partitions=custom_partitions_86, peak_method=args.peak_method)
    
    data_sup86 = None
    if args.add_sup86:
        data_sup86 = compute_partition_data(df_sup86, 'sup86', custom_partitions=custom_partitions_sup86, peak_method=args.peak_method)
    
    # 打印热容信息
    systems_list = [('Air68 (Pt6Sn8)', data_68), ('Air86 (Pt8Sn6)', data_86)]
    if args.add_sup86:
        systems_list.append(('sup86 (Pt8Sn6/support)', data_sup86))
    
    for name, data in systems_list:
        print(f"\n  {name}:")
        for phase, fit in data['phase_fits'].items():
            print(f"    {phase}: Cv={fit['Cv']:.2f}±{fit['Cv_err']:.2f} meV/K, "
                  f"T={fit['T_range'][0]:.0f}-{fit['T_range'][1]:.0f}K")
        if data['T_boundary']:
            print(f"    分界温度: {data['T_boundary']:.0f} K")
        if data['Cv_peak'] is not None:
            print(f"    热容峰({data['peak_method_used']}): Cv_peak={data['Cv_peak']:.2f} meV/K")
    
    # 绘制图片
    output_dir = Path('results/step6_1_1_partition_cv')
    
    if args.interactive:
        # 交互模式（暂不支持三系统）
        if args.add_sup86:
            print(f"\n⚠️ 警告: 交互模式暂不支持三系统对比，将使用非交互模式")
        else:
            print(f"\n>>> 进入交互模式...")
            interactive_adjust_plot(data_68, data_86, output_dir, figsize)
            return
    
    # 非交互模式，使用命令行参数
    print(f"\n>>> 绘制图片...")
    y_integer = not args.no_integer
    params = {
        'figsize': figsize,
        'y_nticks': args.y_nticks,
        'y_integer': y_integer,
        'cv_nticks': args.cv_nticks,
        'cv_integer': y_integer,
        'y_ticks_custom': y_ticks_custom,      # 自定义Y轴刻度
        'cv_ticks_custom': cv_ticks_custom,    # 自定义Cv轴刻度
        'show_error_bars': args.show_error_bars,  # 是否显示误差棒
        'transparent': args.transparent,        # 是否透明背景
    }
    
    if args.add_sup86:
        plot_combined_cv_with_params_three_systems(data_68, data_86, data_sup86, output_dir, params)
    else:
        plot_combined_cv_with_params(data_68, data_86, output_dir, params)
    
    # 导出CSV数据
    print(f"\n>>> 导出CSV数据...")
    if args.add_sup86:
        export_cv_data_to_csv_three_systems(data_68, data_86, data_sup86, output_dir)
    else:
        export_cv_data_to_csv(data_68, data_86, output_dir)
    
    print(f"\n{'='*60}")
    print("✅ 完成!")
    print("=" * 60)


def export_cv_data_to_csv(data_68, data_86, output_dir):
    """导出热容数据到CSV文件"""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. 导出原始数据点 (用于Origin绘图)
    for name, data in [('Air68_Pt6Sn8', data_68), ('Air86_Pt8Sn6', data_86)]:
        df_raw = pd.DataFrame({
            'Temperature_K': data['temps'],
            'Energy_meV': data['E_rel'],
            'Energy_std_meV': data['E_std']
        })
        csv_path = output_dir / f'{name}_energy_data.csv'
        df_raw.to_csv(csv_path, index=False)
        print(f"    原始数据: {csv_path}")
    
    # 2. 导出拟合参数汇总
    summary_rows = []
    
    # Air68 分区拟合
    for phase, fit in data_68['phase_fits'].items():
        summary_rows.append({
            'System': 'Air68_Pt6Sn8',
            'Fit_Type': 'partition',
            'Phase': phase,
            'T_min_K': fit['T_range'][0],
            'T_max_K': fit['T_range'][1],
            'Cv_meV_K': fit['Cv'],
            'Cv_err_meV_K': fit['Cv_err'],
            'Intercept_meV': fit['intercept'],
            'R_squared': fit['R2']
        })
    
    # Air68 单一拟合 (计算)
    temps = np.array(data_68['temps'])
    E_rel = np.array(data_68['E_rel'])
    slope, intercept, r_value, _, std_err = linregress(temps, E_rel)
    summary_rows.append({
        'System': 'Air68_Pt6Sn8',
        'Fit_Type': 'single_linear',
        'Phase': 'all',
        'T_min_K': temps.min(),
        'T_max_K': temps.max(),
        'Cv_meV_K': slope * 1000,
        'Cv_err_meV_K': std_err * 1000,
        'Intercept_meV': intercept,
        'R_squared': r_value**2
    })
    
    # Air86 分区拟合
    for phase, fit in data_86['phase_fits'].items():
        summary_rows.append({
            'System': 'Air86_Pt8Sn6',
            'Fit_Type': 'partition',
            'Phase': phase,
            'T_min_K': fit['T_range'][0],
            'T_max_K': fit['T_range'][1],
            'Cv_meV_K': fit['Cv'],
            'Cv_err_meV_K': fit['Cv_err'],
            'Intercept_meV': fit['intercept'],
            'R_squared': fit['R2']
        })
    
    df_summary = pd.DataFrame(summary_rows)
    csv_summary = output_dir / 'Air_cv_fitting_summary.csv'
    df_summary.to_csv(csv_summary, index=False)
    print(f"    拟合汇总: {csv_summary}")
    
    # 3. 导出拟合线数据 (用于Origin精确绘制拟合线)
    fit_lines = []
    
    # Air68 partition fits
    for phase, fit in data_68['phase_fits'].items():
        T_range = np.linspace(fit['T_range'][0], fit['T_range'][1], 50)
        E_fit = fit['intercept'] + (fit['Cv']/1000) * T_range
        for t, e in zip(T_range, E_fit):
            fit_lines.append({
                'System': 'Air68_Pt6Sn8',
                'Fit_Type': 'partition',
                'Phase': phase,
                'Temperature_K': t,
                'Energy_fit_meV': e
            })
    
    # Air68 single linear fit
    T_full = np.linspace(temps.min(), temps.max(), 100)
    E_single = intercept + slope * T_full
    for t, e in zip(T_full, E_single):
        fit_lines.append({
            'System': 'Air68_Pt6Sn8',
            'Fit_Type': 'single_linear',
            'Phase': 'all',
            'Temperature_K': t,
            'Energy_fit_meV': e
        })
    
    # Air86 partition fits
    for phase, fit in data_86['phase_fits'].items():
        T_range = np.linspace(fit['T_range'][0], fit['T_range'][1], 50)
        E_fit = fit['intercept'] + (fit['Cv']/1000) * T_range
        for t, e in zip(T_range, E_fit):
            fit_lines.append({
                'System': 'Air86_Pt8Sn6',
                'Fit_Type': 'partition',
                'Phase': phase,
                'Temperature_K': t,
                'Energy_fit_meV': e
            })
    
    df_fit_lines = pd.DataFrame(fit_lines)
    csv_fit = output_dir / 'Air_cv_fitting_lines.csv'
    df_fit_lines.to_csv(csv_fit, index=False)
    print(f"    拟合线: {csv_fit}")


def export_cv_data_to_csv_three_systems(data_68, data_86, data_sup86, output_dir):
    """三系统版本：导出热容数据到CSV文件"""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. 导出原始数据点 (用于Origin绘图)
    for name, data in [('Air68_Pt6Sn8', data_68), ('Air86_Pt8Sn6', data_86), ('sup86_Pt8Sn6', data_sup86)]:
        df_raw = pd.DataFrame({
            'Temperature_K': data['temps'],
            'Energy_meV': data['E_rel'],
            'Energy_std_meV': data['E_std']
        })
        csv_path = output_dir / f'{name}_energy_data.csv'
        df_raw.to_csv(csv_path, index=False)
        print(f"    原始数据: {csv_path}")
    
    # 2. 导出拟合参数汇总
    summary_rows = []
    
    # 为每个系统添加拟合数据
    for system_name, data in [('Air68_Pt6Sn8', data_68), ('Air86_Pt8Sn6', data_86), ('sup86_Pt8Sn6', data_sup86)]:
        # 分区拟合
        for phase, fit in data['phase_fits'].items():
            summary_rows.append({
                'System': system_name,
                'Fit_Type': 'partition',
                'Phase': phase,
                'T_min_K': fit['T_range'][0],
                'T_max_K': fit['T_range'][1],
                'Cv_meV_K': fit['Cv'],
                'Cv_err_meV_K': fit['Cv_err'],
                'Intercept_meV': fit['intercept'],
                'R_squared': fit['R2']
            })
    
    df_summary = pd.DataFrame(summary_rows)
    csv_summary = output_dir / 'Three_systems_cv_fitting_summary.csv'
    df_summary.to_csv(csv_summary, index=False)
    print(f"    拟合汇总: {csv_summary}")
    
    # 3. 导出拟合线数据 (用于Origin精确绘制拟合线)
    fit_lines = []
    
    for system_name, data in [('Air68_Pt6Sn8', data_68), ('Air86_Pt8Sn6', data_86), ('sup86_Pt8Sn6', data_sup86)]:
        for phase, fit in data['phase_fits'].items():
            T_range = np.linspace(fit['T_range'][0], fit['T_range'][1], 50)
            E_fit = fit['intercept'] + (fit['Cv']/1000) * T_range
            for t, e in zip(T_range, E_fit):
                fit_lines.append({
                    'System': system_name,
                    'Fit_Type': 'partition',
                    'Phase': phase,
                    'Temperature_K': t,
                    'Energy_fit_meV': e
                })
    
    df_fit_lines = pd.DataFrame(fit_lines)
    csv_fit = output_dir / 'Three_systems_cv_fitting_lines.csv'
    df_fit_lines.to_csv(csv_fit, index=False)
    print(f"    拟合线: {csv_fit}")


if __name__ == '__main__':
    main()
