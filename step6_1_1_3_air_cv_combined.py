#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 6.1.1.3: Air68 vs Air86 分区热容组合图

将 Air68 和 Air86 的分区热容图并排显示，统一Y轴范围便于对比

用法:
  python step6_1_1_3_air_cv_combined.py

作者: AI Assistant
日期: 2025-12-01
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress
from pathlib import Path

# 设置高质量论文图样式（与原图一致）
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'SimHei']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['font.size'] = 11
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['xtick.major.width'] = 1.2
plt.rcParams['ytick.major.width'] = 1.2


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


def plot_combined_cv(data_68, data_86, output_dir):
    """绘制组合图"""
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 创建图表 - 三列，保持原图高宽比 (8:6)
    fig, (ax1_left, ax1_mid, ax1_right) = plt.subplots(1, 3, figsize=(24, 6))
    
    # ========== 统一Y轴范围 ==========
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
    
    # ========== 绘制函数 ==========
    def plot_single(ax1, data, title):
        """绘制单个分区热容图"""
        temps = data['temps']
        E_rel = data['E_rel']
        E_std = data['E_std']
        phase_fits = data['phase_fits']
        phases = sorted(data['phases'])
        T_boundary = data['T_boundary']
        
        # 左Y轴: 能量数据点
        ax1.errorbar(temps, E_rel, yerr=E_std,
                     fmt='o', markersize=7, color='black',
                     ecolor='gray', elinewidth=1.5, capsize=3, capthick=1.5,
                     zorder=5)
        
        # 拟合线
        for phase in phases:
            if phase in phase_fits:
                fit = phase_fits[phase]
                T_fit = np.linspace(fit['T_range'][0], fit['T_range'][1], 50)
                E_fit = fit['slope'] * T_fit + fit['intercept']
                ax1.plot(T_fit, E_fit, '-', color='black', linewidth=2, zorder=4)
        
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
                ax1.plot([T1_end, T2_start], [E1, E2], '-', color='black', linewidth=2, zorder=4)
        
        ax1.set_xlabel('Temperature (K)', fontsize=13, fontweight='bold')
        ax1.set_ylabel('Total Energy (eV)', fontsize=13, fontweight='bold')
        ax1.set_ylim(E_ylim)
        ax1.tick_params(axis='both', labelsize=11)
        ax1.set_title(title, fontsize=14, fontweight='bold', pad=10)
        
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
                
                ax2.plot(T_plot, Cv_plot, 'r-', linewidth=2, zorder=3)
            else:
                # 阶梯形热容曲线
                ax2.plot([temps.min(), T_boundary], [Cv1, Cv1], 'r-', linewidth=2, zorder=3)
                ax2.plot([T_boundary, T_boundary], [Cv1, Cv2], 'r--', linewidth=1.5, zorder=3)
                ax2.plot([T_boundary, temps.max()], [Cv2, Cv2], 'r-', linewidth=2, zorder=3)
        else:
            Cv_single = list(phase_fits.values())[0]['Cv']
            ax2.axhline(y=Cv_single, color='red', linewidth=2, zorder=3)
        
        ax2.set_ylabel(r'$C_v$ (meV/K)', fontsize=13, fontweight='bold', color='red')
        ax2.tick_params(axis='y', labelcolor='red', labelsize=11, color='red')
        ax2.spines['right'].set_color('red')
        ax2.set_ylim(Cv_ylim)
        
        return ax2
    
    def plot_single_fit(ax1, data, title):
        """绘制单一线性拟合图（不分区）"""
        temps = data['temps']
        E_rel = data['E_rel']
        E_std = data['E_std']
        
        # 左Y轴: 能量数据点
        ax1.errorbar(temps, E_rel, yerr=E_std,
                     fmt='o', markersize=7, color='black',
                     ecolor='gray', elinewidth=1.5, capsize=3, capthick=1.5,
                     zorder=5)
        
        # 整体线性拟合（单一拟合线）
        slope, intercept, r_value, _, std_err = linregress(temps, E_rel)
        Cv_overall = slope * 1000  # meV/K
        Cv_err = std_err * 1000
        R2 = r_value ** 2
        
        # 绘制拟合线
        T_fit = np.linspace(temps.min(), temps.max(), 100)
        E_fit = slope * T_fit + intercept
        ax1.plot(T_fit, E_fit, '-', color='black', linewidth=2, zorder=4)
        
        ax1.set_xlabel('Temperature (K)', fontsize=13, fontweight='bold')
        ax1.set_ylabel('Total Energy (eV)', fontsize=13, fontweight='bold')
        ax1.set_ylim(E_ylim)
        ax1.tick_params(axis='both', labelsize=11)
        ax1.set_title(title, fontsize=14, fontweight='bold', pad=10)
        
        # 右Y轴: 热容（单一水平线）
        ax2 = ax1.twinx()
        ax2.axhline(y=Cv_overall, color='red', linewidth=2, zorder=3)
        
        ax2.set_ylabel(r'$C_v$ (meV/K)', fontsize=13, fontweight='bold', color='red')
        ax2.tick_params(axis='y', labelcolor='red', labelsize=11, color='red')
        ax2.spines['right'].set_color('red')
        ax2.set_ylim(Cv_ylim)
        
        print(f"    单一拟合: Cv={Cv_overall:.2f}±{Cv_err:.2f} meV/K, R²={R2:.4f}")
        
        return ax2, Cv_overall
    
    # ========== 绘制三个子图 ==========
    # 左图: Pt8Sn6 (Air86) - 分区拟合
    ax2_left = plot_single(ax1_left, data_86, r'Pt$_8$Sn$_6$')
    
    # 中图: Pt6Sn8 (Air68) - 分区拟合
    ax2_mid = plot_single(ax1_mid, data_68, r'Pt$_6$Sn$_8$ (partition)')
    
    # 右图: Pt6Sn8 (Air68) - 单一拟合
    print(f"\n  Pt6Sn8 单一拟合:")
    ax2_right, Cv_68_single = plot_single_fit(ax1_right, data_68, r'Pt$_6$Sn$_8$ (single fit)')
    
    plt.tight_layout()
    
    # 保存
    output_file = output_dir / 'Air68_Air86_cv_combined.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    print(f"📊 组合图已保存: {output_file}")
    
    return output_file


def main():
    print("=" * 60)
    print("Step 6.1.1.3: Air68 vs Air86 分区热容组合图")
    print("=" * 60)
    
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
    
    # 绘制组合图
    print(f"\n>>> 绘制组合图...")
    output_dir = Path('results/step6_1_1_partition_cv')
    plot_combined_cv(data_68, data_86, output_dir)
    
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
