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


def load_cluster_data(csv_path):
    """加载聚类结果数据"""
    try:
        df = pd.read_csv(csv_path)
        return df
    except Exception as e:
        print(f"  错误: 无法读取 {csv_path}: {e}")
        return None


def compute_partition_data(df, structure_name):
    """计算分区热容数据"""
    
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
    
    # 分区拟合
    phases = sorted(df['phase_clustered'].unique())
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
    
    # 分界温度
    T_boundary = None
    if len(phases) >= 2:
        phase1_temps = [t for t, p in temp_to_partition.items() if p == phases[0]]
        phase2_temps = [t for t, p in temp_to_partition.items() if p == phases[1]]
        if phase1_temps and phase2_temps:
            T1_last = max(phase1_temps)
            T2_first = min(phase2_temps)
            T_boundary = (T1_last + T2_first) / 2
    
    return {
        'temps': temps_unique,
        'E_rel': E_rel,
        'E_std': E_std,
        'temp_to_partition': temp_to_partition,
        'phase_fits': phase_fits,
        'T_boundary': T_boundary,
        'phases': phases
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
        if len(data['phases']) >= 2:
            phases = sorted(data['phases'])
            fit1 = data['phase_fits'].get(phases[0])
            fit2 = data['phase_fits'].get(phases[1])
            if fit1 and fit2:
                T1_last = fit1['T_range'][1]
                T2_first = fit2['T_range'][0]
                idx1 = np.where(data['temps'] == T1_last)[0]
                idx2 = np.where(data['temps'] == T2_first)[0]
                if len(idx1) > 0 and len(idx2) > 0:
                    E1 = data['E_rel'][idx1[0]]
                    E2 = data['E_rel'][idx2[0]]
                    Cv_trans = (E2 - E1) / (T2_first - T1_last) * 1000
                    all_Cv.append(Cv_trans)
    
    Cv_min = min(all_Cv) * 0.85
    Cv_max = max(all_Cv) * 1.1
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
        
        # 计算过渡区热容
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
            # 带峰的热容曲线
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
            # 阶梯形热容曲线
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
    
    # 子图1: Pt8Sn6 (Air86)
    plot_single_partition_with_params(data_86, r'Pt$_8$Sn$_6$', 
                                      output_dir / 'Air86_Pt8Sn6_partition_cv.png',
                                      E_ylim, Cv_ylim, figsize, y_nticks, y_integer, 
                                      cv_nticks, cv_integer, y_ticks_custom, cv_ticks_custom)
    
    # 子图2: Pt6Sn8 (Air68) - 分区
    plot_single_partition_with_params(data_68, r'Pt$_6$Sn$_8$ (partition)', 
                                      output_dir / 'Air68_Pt6Sn8_partition_cv.png',
                                      E_ylim, Cv_ylim, figsize, y_nticks, y_integer, 
                                      cv_nticks, cv_integer, y_ticks_custom, cv_ticks_custom)
    
    # 子图3: Pt6Sn8 (Air68) - 单一拟合
    plot_single_linear_fit_with_params(data_68, r'Pt$_6$Sn$_8$ (single fit)', 
                                       output_dir / 'Air68_Pt6Sn8_single_fit_cv.png',
                                       E_ylim, Cv_ylim, figsize, y_nticks, y_integer, 
                                       cv_nticks, cv_integer, y_ticks_custom, cv_ticks_custom)


def plot_single_partition_with_params(data, title, output_path, E_ylim, Cv_ylim, 
                                      figsize, y_nticks, y_integer, cv_nticks, cv_integer,
                                      y_ticks_custom=None, cv_ticks_custom=None):
    """带参数的分区热容图绘制"""
    from matplotlib.ticker import MaxNLocator, LinearLocator, MultipleLocator
    
    fig, ax1 = plt.subplots(figsize=figsize)
    
    temps = data['temps']
    E_rel = data['E_rel']
    E_std = data['E_std']
    phase_fits = data['phase_fits']
    phases = sorted(data['phases'])
    T_boundary = data['T_boundary']
    
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
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"📊 已保存: {output_path}")


def plot_single_linear_fit_with_params(data, title, output_path, E_ylim, Cv_ylim,
                                       figsize, y_nticks, y_integer, cv_nticks, cv_integer,
                                       y_ticks_custom=None, cv_ticks_custom=None):
    """带参数的单一拟合图绘制"""
    from matplotlib.ticker import MaxNLocator, LinearLocator, MultipleLocator
    
    fig, ax1 = plt.subplots(figsize=figsize)
    
    temps = data['temps']
    E_rel = data['E_rel']
    E_std = data['E_std']
    
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
    args = parser.parse_args()
    
    # 解析figsize
    try:
        fig_w, fig_h = map(float, args.figsize.lower().split('x'))
        figsize = (fig_w, fig_h)
    except ValueError:
        print(f"警告: 无效的figsize格式 '{args.figsize}'，使用默认 10x8")
        figsize = (10, 8)
    
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
    
    print("=" * 60)
    print("Step 6.1.1.3: Air68 vs Air86 分区热容独立子图")
    print("=" * 60)
    print(f"  图片尺寸: {figsize[0]}x{figsize[1]}")
    
    # 加载数据
    base_dir = Path('results/step6_1_clustering')
    
    csv_68 = base_dir / 'Air68_kmeans_n2_clustered_data.csv'
    csv_86 = base_dir / 'Air86_kmeans_n2_clustered_data.csv'
    
    if not csv_68.exists():
        print(f"错误: 找不到 {csv_68}")
        return
    if not csv_86.exists():
        print(f"错误: 找不到 {csv_86}")
        return
    
    print(f"\n>>> 加载数据...")
    df_68 = load_cluster_data(csv_68)
    df_86 = load_cluster_data(csv_86)
    
    if df_68 is None or df_86 is None:
        return
    
    print(f"    Air68: {len(df_68)} 条记录")
    print(f"    Air86: {len(df_86)} 条记录")
    
    # 计算分区数据
    print(f"\n>>> 计算分区热容...")
    data_68 = compute_partition_data(df_68, 'Air68')
    data_86 = compute_partition_data(df_86, 'Air86')
    
    # 打印热容信息
    for name, data in [('Air68 (Pt6Sn8)', data_68), ('Air86 (Pt8Sn6)', data_86)]:
        print(f"\n  {name}:")
        for phase, fit in data['phase_fits'].items():
            print(f"    {phase}: Cv={fit['Cv']:.2f}±{fit['Cv_err']:.2f} meV/K, "
                  f"T={fit['T_range'][0]:.0f}-{fit['T_range'][1]:.0f}K")
        if data['T_boundary']:
            print(f"    分界温度: {data['T_boundary']:.0f} K")
    
    # 绘制图片
    output_dir = Path('results/step6_1_1_partition_cv')
    
    if args.interactive:
        # 交互模式
        print(f"\n>>> 进入交互模式...")
        interactive_adjust_plot(data_68, data_86, output_dir, figsize)
    else:
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
        }
        plot_combined_cv_with_params(data_68, data_86, output_dir, params)
    
    # 导出CSV数据
    print(f"\n>>> 导出CSV数据...")
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


if __name__ == '__main__':
    main()
