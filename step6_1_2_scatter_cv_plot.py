#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 6.1.2: 分区散点图 + 热容双Y轴图 + f图（分区分布）

================================================================================
                              计算逻辑详解
================================================================================

1. 数据来源
-----------
   输入文件: step6_1 生成的聚类结果 CSV 文件
   路径: results/step6_1_clustering/{structure}_kmeans_n2_clustered_data.csv
   
   CSV 文件包含以下关键列:
   - temp: 温度 (K)
   - avg_energy: 平均能量 (eV)
   - phase_clustered: K-means 聚类结果 (partition1/partition2)
   
   每个温度点有多个 run（重复模拟），每个 run 被独立聚类。

2. 能量处理
-----------
   对于载体支撑的团簇（非 Air 系列）:
   - 加载载体能量数据: data/lammps_energy/sup/energy_master_*.csv
   - 拟合载体能量-温度关系: E_support = slope_support * T + intercept_support
   - 团簇能量 = 总能量 - 载体能量
   
   对于气相团簇（Air 系列）:
   - 直接使用总能量作为团簇能量（不扣除载体）
   
   相对能量计算:
   - E_ref = 最低温度点的平均团簇能量
   - E_cluster_rel = E_cluster - E_ref

3. 多数投票机制（Majority Voting）
----------------------------------
   目的: 将每个温度点分配到唯一的分区（用于热容拟合）
   
   方法: 对于每个温度 T:
   1. 统计该温度所有 run 的分区归属: {partition1: n1, partition2: n2}
   2. 选择数量最多的分区作为该温度的分区
   
   示例:
   - T=700K: {partition2: 9, partition1: 1} → 分配给 partition2
   - T=600K: {partition1: 10} → 分配给 partition1
   
   注意: 这与 f 图不同，f 图显示的是原始分布，不进行投票

4. 分区热容拟合
---------------
   对每个分区内的温度点进行线性拟合:
   
   E_mean(T) = slope * T + intercept
   
   热容计算:
   Cv = slope × 1000  (单位: meV/K)
   
   输出:
   - Cv: 热容值 (meV/K)
   - R²: 拟合优度
   - T_range: 温度范围

5. 热容峰检测
-------------
   分界温度 T_boundary = (T1_last + T2_first) / 2
   其中:
   - T1_last: partition1 的最高温度
   - T2_first: partition2 的最低温度
   
   数值微分计算过渡区热容:
   Cv_transition = (E2 - E1) / (T2_first - T1_last) × 1000
   
   判断热容峰:
   - 如果 Cv_transition > max(Cv1, Cv2): 存在热容峰
   - 否则: 无热容峰（连续相变）

6. 热容曲线绘制
---------------
   有热容峰时 (使用 sigmoid + Gaussian):
   - baseline(T) = Cv1 + (Cv2 - Cv1) × sigmoid(T - T_boundary)
   - gaussian(T) = (Cv_peak - baseline) × exp(-0.5 × ((T - T_boundary)/σ)²)
   - Cv(T) = baseline + gaussian
   
   无热容峰时:
   - 使用阶梯函数: T < T_boundary → Cv1, 否则 → Cv2

7. f 图：分区分布堆叠柱状图
---------------------------
   显示每个温度点各分区的 run 数量分布
   
   与多数投票的区别:
   - 多数投票: 用于热容拟合，每个温度只属于一个分区
   - f 图: 显示原始聚类结果，可以看到分区共存现象
   
   例如 T=700K 可能有 9 个 run 在 partition2，1 个 run 在 partition1

================================================================================

输出文件:
---------
1. {structure}_scatter_cv.png/pdf  - 主图（散点+Cv图 + f图）
2. {structure}_scatter_data.csv    - 散点数据（每个run的数据）
3. {structure}_mean_data.csv       - 平均能量数据（按温度）
4. {structure}_cv_curve.csv        - 热容曲线数据
5. {structure}_phase_distribution.csv - 分区分布数据（f图数据）
6. {structure}_fit_params.csv      - 拟合参数汇总

用法:
-----
  python step6_1_2_scatter_cv_plot.py --structure Pt8sn6
  python step6_1_2_scatter_cv_plot.py --structure Air86 --format pdf
  python step6_1_2_scatter_cv_plot.py --structure all
  python step6_1_2_scatter_cv_plot.py --structure Pt8sn6 --no-f  # 不绘制f图
  python step6_1_2_scatter_cv_plot.py --list

作者: AI Assistant
日期: 2025-11-30
"""

import os
import sys
import glob
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress
from pathlib import Path
from datetime import datetime

# 设置高质量论文图样式
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'SimHei']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['font.size'] = 11
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['xtick.major.width'] = 1.2
plt.rcParams['ytick.major.width'] = 1.2

# 分区颜色配置（与 step6_1 h 图一致）
PARTITION_COLORS = {
    'partition1': '#3498db',  # 蓝色 - 分区1（低温相/固态）
    'partition2': '#f39c12',  # 橙色 - 分区2（高温相/液态）
    'partition3': '#e74c3c',  # 红色 - 分区3（预熔化）
    'Phase_1': '#3498db',
    'Phase_2': '#f39c12',
    'Solid': '#3498db',
    'Liquid': '#f39c12',
}

# 载体热容 (meV/K)
CV_SUPPORT = 38.2151


def find_clustering_results(base_dir='results/step6_1_clustering'):
    """查找所有可用的聚类结果"""
    results = {}
    pattern = os.path.join(base_dir, '*_kmeans_n2_clustered_data.csv')
    files = glob.glob(pattern)
    
    for f in files:
        basename = os.path.basename(f)
        structure = basename.replace('_kmeans_n2_clustered_data.csv', '')
        results[structure] = f
    
    return results


def load_support_energy_data():
    """加载载体能量数据"""
    support_csv = 'data/lammps_energy/sup/energy_master_20251021_151520.csv'
    
    if not os.path.exists(support_csv):
        return None
    
    try:
        df_support = pd.read_csv(support_csv)
        if 'temp' in df_support.columns and 'avg_energy' in df_support.columns:
            T = df_support['temp'].values
            E = df_support['avg_energy'].values
            slope, intercept, r_value, _, _ = linregress(T, E)
            return slope, intercept, r_value**2
    except Exception as e:
        print(f"  警告: 读取载体能量数据失败: {e}")
    
    return None


def plot_scatter_cv(df, structure_name, output_dir, output_format='png', dpi=300, add_f_plot=True):
    """
    绘制分区散点图 + 热容双Y轴图 + f图（堆叠柱状图）
    
    特点：
    1. 散点图显示每个 run 的数据点（按分区着色）
    2. 左Y轴：能量
    3. 右Y轴：热容曲线
    4. 热容峰：使用分界点的数值微分
    5. f图：各温度点的相态分布堆叠柱状图
    """
    
    print(f"\n>>> 绘制 {structure_name} 分区散点+热容图...")
    
    # 检查必要列
    required_cols = ['temp', 'avg_energy', 'phase_clustered']
    if not all(col in df.columns for col in required_cols):
        print(f"  错误: 缺少必要列 {required_cols}")
        return None
    
    # 判断是否是 Air 系列
    is_air_system = structure_name.startswith('Air') or structure_name in ['68', '86']
    
    # 加载载体能量数据
    if is_air_system:
        slope_support = 0.0
        intercept_support = 0.0
        print(f"  [Air系列] 气相纳米团簇，不扣除载体能量")
    else:
        support_fit = load_support_energy_data()
        if support_fit is not None:
            slope_support, intercept_support, R2_support = support_fit
        else:
            slope_support = CV_SUPPORT / 1000
            T_min = df['temp'].min()
            E_total_min = df[df['temp'] == T_min]['avg_energy'].mean()
            intercept_support = E_total_min * 0.9 - slope_support * T_min
            print(f"  [警告] 使用默认Cv_support估算载体能量")
    
    # ========== 1. 计算团簇能量 ==========
    if is_air_system:
        df['E_cluster'] = df['avg_energy']
    else:
        df['E_cluster'] = df['avg_energy'] - (slope_support * df['temp'] + intercept_support)
    
    # 计算相对能量（相对于最低温度的平均值）
    E_ref = df[df['temp'] == df['temp'].min()]['E_cluster'].mean()
    df['E_cluster_rel'] = df['E_cluster'] - E_ref
    
    # ========== 2. 按温度分组统计 ==========
    temp_groups = df.groupby('temp')
    temps_unique = np.array(sorted(df['temp'].unique()))
    E_mean = np.array([temp_groups.get_group(t)['E_cluster_rel'].mean() for t in temps_unique])
    E_std = np.array([temp_groups.get_group(t)['E_cluster_rel'].std() for t in temps_unique])
    
    # ========== 3. 多数投票确定每个温度的分区 ==========
    temp_to_partition = {}
    print(f"\n  多数投票温度分配:")
    
    for temp in temps_unique:
        df_temp = df[df['temp'] == temp]
        partition_counts = df_temp['phase_clustered'].value_counts()
        dominant_partition = partition_counts.idxmax()
        temp_to_partition[temp] = dominant_partition
        print(f"    T={temp:4.0f}K: {dict(partition_counts)} -> {dominant_partition}")
    
    # ========== 4. 分区拟合 ==========
    phases = sorted(df['phase_clustered'].unique())
    phase_fits = {}
    
    for phase in phases:
        phase_temps = [temp for temp, part in temp_to_partition.items() if part == phase]
        phase_temps = sorted(phase_temps)
        
        if len(phase_temps) >= 2:
            mask = np.isin(temps_unique, phase_temps)
            T_phase = temps_unique[mask]
            E_phase = E_mean[mask]
            
            slope_ph, intercept_ph, r_value_ph, _, std_err_ph = linregress(T_phase, E_phase)
            
            phase_fits[phase] = {
                'slope': slope_ph,
                'intercept': intercept_ph,
                'R2': r_value_ph ** 2,
                'Cv': slope_ph * 1000,  # meV/K
                'Cv_err': std_err_ph * 1000,
                'T_range': (T_phase.min(), T_phase.max()),
            }
            
            print(f"  {phase}: Cv={slope_ph*1000:.4f} meV/K, R2={r_value_ph**2:.4f}, T={T_phase.min():.0f}-{T_phase.max():.0f}K")
    
    # ========== 5. 计算热容峰 ==========
    phases_sorted = sorted(phase_fits.keys())
    T_boundary = None
    Cv_peak = None
    Cv1 = None
    Cv2 = None
    
    if len(phases_sorted) >= 2:
        phase1_temps = [t for t, p in temp_to_partition.items() if p == phases_sorted[0]]
        phase2_temps = [t for t, p in temp_to_partition.items() if p == phases_sorted[1]]
        
        if phase1_temps and phase2_temps:
            T1_last = max(phase1_temps)
            T2_first = min(phase2_temps)
            T_boundary = (T1_last + T2_first) / 2
            
            Cv1 = phase_fits[phases_sorted[0]]['Cv']
            Cv2 = phase_fits[phases_sorted[1]]['Cv']
            
            # 用分界点两侧的平均能量求数值微分
            idx1 = np.where(temps_unique == T1_last)[0]
            idx2 = np.where(temps_unique == T2_first)[0]
            
            if len(idx1) > 0 and len(idx2) > 0:
                E1 = E_mean[idx1[0]]
                E2 = E_mean[idx2[0]]
                Cv_transition = (E2 - E1) / (T2_first - T1_last) * 1000  # meV/K
                
                # 判断是否存在热容峰
                has_peak = Cv_transition > max(Cv1, Cv2)
                if has_peak:
                    Cv_peak = Cv_transition
                    print(f"\n  分界温度: {T_boundary:.0f} K")
                    print(f"  * 存在热容峰: Cv_peak={Cv_peak:.2f} meV/K")
                else:
                    print(f"\n  分界温度: {T_boundary:.0f} K (无峰)")
            
            print(f"  热容: Cv1={Cv1:.2f}, Cv2={Cv2:.2f} meV/K")
    
    # ========== 6. 计算每个温度的分区分布（用于f图）==========
    temp_phase_counts = df.groupby(['temp', 'phase_clustered']).size().unstack(fill_value=0)
    temp_sorted = sorted(df['temp'].unique())
    
    # ========== 7. 绘制双Y轴散点图（如果add_f_plot则创建双图布局）==========
    if add_f_plot:
        fig, (ax1, ax3) = plt.subplots(1, 2, figsize=(16, 7), gridspec_kw={'width_ratios': [1.2, 1]})
    else:
        fig, ax1 = plt.subplots(figsize=(10, 7))
    
    # ----- 左Y轴: 能量散点图（按分区着色）-----
    for phase in phases_sorted:
        df_phase = df[df['phase_clustered'] == phase]
        color = PARTITION_COLORS.get(phase, '#95a5a6')
        
        ax1.scatter(df_phase['temp'], df_phase['E_cluster_rel'],
                   c=color, s=50, alpha=0.6,
                   edgecolors='white', linewidth=0.5,
                   label=f'{phase} (n={len(df_phase)})',
                   zorder=3)
    
    # 绘制分区拟合线
    for phase in phases_sorted:
        fit = phase_fits[phase]
        color = PARTITION_COLORS.get(phase, '#95a5a6')
        T_fit = np.linspace(fit['T_range'][0], fit['T_range'][1], 50)
        E_fit = fit['slope'] * T_fit + fit['intercept']
        ax1.plot(T_fit, E_fit, '-', color=color, linewidth=2.5, zorder=4)
    
    # 连接两个分区的数据点
    if len(phases_sorted) >= 2:
        fit1 = phase_fits[phases_sorted[0]]
        fit2 = phase_fits[phases_sorted[1]]
        T1_end = fit1['T_range'][1]
        T2_start = fit2['T_range'][0]
        
        idx1 = np.where(temps_unique == T1_end)[0]
        idx2 = np.where(temps_unique == T2_start)[0]
        if len(idx1) > 0 and len(idx2) > 0:
            E1_end = E_mean[idx1[0]]
            E2_start = E_mean[idx2[0]]
            ax1.plot([T1_end, T2_start], [E1_end, E2_start], 
                    '--', color='gray', linewidth=1.5, zorder=2)
    
    ax1.set_xlabel('Temperature (K)', fontsize=13, fontweight='bold')
    ax1.set_ylabel('Total Energy (eV)', fontsize=13, fontweight='bold')
    ax1.tick_params(axis='both', labelsize=11)
    ax1.legend(loc='lower right', fontsize=10, framealpha=0.9)
    
    # ----- 右Y轴: 热容曲线 -----
    ax2 = ax1.twinx()
    
    if len(phases_sorted) >= 2 and Cv1 is not None and Cv2 is not None:
        T_plot = np.linspace(temps_unique.min(), temps_unique.max(), 500)
        Cv_plot = np.zeros_like(T_plot)
        
        if Cv_peak is not None:
            # 有热容峰：平滑高斯峰
            sigma = (T2_first - T1_last) / 2
            
            for i, T in enumerate(T_plot):
                transition = 1 / (1 + np.exp(-(T - T_boundary) / (sigma * 0.5)))
                baseline = Cv1 + (Cv2 - Cv1) * transition
                gaussian = (Cv_peak - baseline) * np.exp(-0.5 * ((T - T_boundary) / sigma)**2)
                Cv_plot[i] = baseline + gaussian
        else:
            # 无热容峰：阶梯函数
            for i, T in enumerate(T_plot):
                if T < T_boundary:
                    Cv_plot[i] = Cv1
                else:
                    Cv_plot[i] = Cv2
        
        ax2.plot(T_plot, Cv_plot, 'r-', linewidth=2.5, zorder=5, label='Cv')
        
        # 设置Y轴范围
        cv_values = [Cv1, Cv2]
        if Cv_peak:
            cv_values.append(Cv_peak)
        cv_min = min(cv_values) * 0.8
        cv_max = max(cv_values) * 1.15
        ax2.set_ylim(cv_min, cv_max)
    else:
        # 单分区
        Cv_single = list(phase_fits.values())[0]['Cv'] if phase_fits else 0
        ax2.axhline(y=Cv_single, color='red', linewidth=2, zorder=5)
        Cv1 = Cv_single
        Cv2 = Cv_single
    
    ax2.set_ylabel('Cv (meV/K)', fontsize=13, fontweight='bold', color='red')
    ax2.tick_params(axis='y', labelcolor='red', labelsize=11, color='red')
    ax2.spines['right'].set_color('red')
    
    # 标题
    ax1.set_title(f'(a) {structure_name}: Energy & Cv', fontsize=14, fontweight='bold', pad=10)
    
    # ========== 8. 绘制 f 图：各温度点的相态分布堆叠柱状图 ==========
    if add_f_plot:
        # 准备堆叠数据（使用 phases_sorted 中的分区顺序）
        x_pos = np.arange(len(temp_sorted))
        bottom = np.zeros(len(temp_sorted))
        
        for phase in phases_sorted:
            if phase in temp_phase_counts.columns:
                counts = [temp_phase_counts.loc[t, phase] if t in temp_phase_counts.index else 0 
                          for t in temp_sorted]
            else:
                counts = [0] * len(temp_sorted)
            
            color = PARTITION_COLORS.get(phase, '#95a5a6')
            ax3.bar(x_pos, counts, bottom=bottom, label=phase, 
                   color=color, alpha=0.8, edgecolor='black', linewidth=0.5)
            bottom = bottom + np.array(counts)
        
        ax3.set_xlabel('Temperature (K)', fontsize=13, fontweight='bold')
        ax3.set_ylabel('Number of Runs', fontsize=13, fontweight='bold')
        ax3.set_title(f'(b) Phase Distribution', fontsize=14, fontweight='bold', pad=10)
        
        # 只显示部分温度标签（避免拥挤）
        step = max(1, len(temp_sorted) // 10)
        ax3.set_xticks(x_pos[::step])
        ax3.set_xticklabels([f'{int(t)}' for t in temp_sorted[::step]], rotation=45, ha='right')
        ax3.legend(fontsize=10, loc='upper right', framealpha=0.9)
        ax3.grid(True, alpha=0.3, linestyle=':', axis='y', linewidth=0.8)
    
    plt.tight_layout()
    
    # 保存图片
    output_file = Path(output_dir) / f'{structure_name}_scatter_cv.{output_format}'
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close()
    
    print(f"\n  图已保存: {output_file}")
    
    # ========== 9. 导出数据供 Origin 使用 ==========
    # 导出散点数据
    df_scatter = df[['temp', 'E_cluster_rel', 'phase_clustered']].copy()
    df_scatter.columns = ['Temperature_K', 'Energy_eV', 'Partition']
    scatter_csv = Path(output_dir) / f'{structure_name}_scatter_data.csv'
    df_scatter.to_csv(scatter_csv, index=False)
    print(f"  散点数据已导出: {scatter_csv}")
    
    # 导出平均能量数据
    df_mean = pd.DataFrame({
        'Temperature_K': temps_unique,
        'Energy_mean_eV': E_mean,
        'Energy_std_eV': E_std,
        'Partition': [temp_to_partition.get(t, 'unknown') for t in temps_unique]
    })
    mean_csv = Path(output_dir) / f'{structure_name}_mean_data.csv'
    df_mean.to_csv(mean_csv, index=False)
    print(f"  平均能量数据已导出: {mean_csv}")
    
    # 导出热容曲线
    if 'T_plot' in dir() and 'Cv_plot' in dir():
        df_cv = pd.DataFrame({
            'Temperature_K': T_plot,
            'Cv_meV_K': Cv_plot
        })
        cv_csv = Path(output_dir) / f'{structure_name}_cv_curve.csv'
        df_cv.to_csv(cv_csv, index=False)
        print(f"  热容曲线已导出: {cv_csv}")
    
    # 导出 f 图数据：各温度的分区分布
    if add_f_plot:
        f_data = {'Temperature_K': temp_sorted}
        for phase in phases_sorted:
            if phase in temp_phase_counts.columns:
                f_data[phase] = [temp_phase_counts.loc[t, phase] if t in temp_phase_counts.index else 0 
                                 for t in temp_sorted]
            else:
                f_data[phase] = [0] * len(temp_sorted)
        df_f = pd.DataFrame(f_data)
        f_csv = Path(output_dir) / f'{structure_name}_phase_distribution.csv'
        df_f.to_csv(f_csv, index=False)
        print(f"  分区分布数据已导出: {f_csv}")
    
    # 导出拟合参数
    fit_summary = {
        'structure': structure_name,
        'T_boundary_K': T_boundary,
        'Cv_peak_meV_K': Cv_peak,
    }
    for i, (phase, fit) in enumerate(phase_fits.items()):
        fit_summary[f'phase_{i+1}_name'] = phase
        fit_summary[f'phase_{i+1}_Cv_meV_K'] = fit['Cv']
        fit_summary[f'phase_{i+1}_R2'] = fit['R2']
        fit_summary[f'phase_{i+1}_T_min_K'] = fit['T_range'][0]
        fit_summary[f'phase_{i+1}_T_max_K'] = fit['T_range'][1]
    
    fit_csv = Path(output_dir) / f'{structure_name}_fit_params.csv'
    pd.DataFrame([fit_summary]).to_csv(fit_csv, index=False)
    print(f"  拟合参数已导出: {fit_csv}")
    
    return {
        'structure': structure_name,
        'T_boundary': T_boundary,
        'Cv_peak': Cv_peak,
        'partitions': phase_fits
    }


def list_available_structures(base_dir='results/step6_1_clustering'):
    """列出所有可用的结构"""
    results = find_clustering_results(base_dir)
    
    print("\n" + "=" * 60)
    print("可用结构列表")
    print("=" * 60)
    
    # 分类
    air_series = []
    pt6_series = []
    pt8_series = []
    oxide_series = []
    other = []
    
    for name in sorted(results.keys()):
        name_lower = name.lower()
        if 'air' in name_lower:
            air_series.append(name)
        elif name == 'Cv' or 'o' in name_lower:
            oxide_series.append(name)
        elif name_lower.startswith('pt6'):
            pt6_series.append(name)
        elif name_lower.startswith('pt8'):
            pt8_series.append(name)
        else:
            other.append(name)
    
    print(f"\n[Air] 气相团簇 ({len(air_series)}): {', '.join(air_series) if air_series else 'None'}")
    print(f"[Pt6] Pt6系列 ({len(pt6_series)}): {', '.join(sorted(pt6_series)) if pt6_series else 'None'}")
    print(f"[Pt8] Pt8系列 ({len(pt8_series)}): {', '.join(sorted(pt8_series)) if pt8_series else 'None'}")
    print(f"[O] 含氧团簇 ({len(oxide_series)}): {', '.join(sorted(oxide_series)) if oxide_series else 'None'}")
    if other:
        print(f"[Other] ({len(other)}): {', '.join(sorted(other))}")
    
    print(f"\nTotal: {len(results)} structures")
    print("=" * 60)
    
    return results


def parse_args():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description='Step 6.1.2: 分区散点图 + 热容双Y轴图 + f图',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  %(prog)s --structure Pt8sn6
  %(prog)s --structure Air86 --format pdf
  %(prog)s --structure all --dpi 600
  %(prog)s --structure Pt8sn6 --no-f  # 不绘制f图
  %(prog)s --list
        '''
    )
    
    parser.add_argument('--structure', '-s', type=str, default=None,
                        help='Structure name (e.g., Pt8sn6, Air86) or "all"')
    parser.add_argument('--list', '-l', action='store_true',
                        help='List all available structures')
    parser.add_argument('--format', '-f', type=str, default='png',
                        choices=['png', 'pdf', 'svg', 'eps'],
                        help='Output format (default: png)')
    parser.add_argument('--dpi', type=int, default=300,
                        help='Output resolution (default: 300)')
    parser.add_argument('--output-dir', '-o', type=str, 
                        default='results/step6_1_2_scatter_cv',
                        help='Output directory')
    parser.add_argument('--no-f', action='store_true',
                        help='Do not include f-plot (stacked bar chart)')
    
    return parser.parse_args()


def main():
    """主函数"""
    args = parse_args()
    
    print("=" * 70)
    print("Step 6.1.2: 分区散点图 + 热容双Y轴图")
    print("=" * 70)
    
    if args.list:
        list_available_structures()
        return
    
    if args.structure is None:
        print("Error: Please specify --structure or use --list")
        return
    
    # 创建输出目录
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 获取可用结构
    available = find_clustering_results()
    
    if args.structure.lower() == 'all':
        structures = list(available.keys())
        print(f"\nProcessing all {len(structures)} structures...")
    else:
        structures = [args.structure]
    
    # 处理每个结构
    results = []
    success = 0
    failed = 0
    
    for structure in structures:
        # 查找结构（大小写不敏感）
        found_name = None
        for name in available.keys():
            if name.lower() == structure.lower():
                found_name = name
                break
        
        if found_name is None:
            print(f"\nWarning: Structure '{structure}' not found")
            failed += 1
            continue
        
        csv_path = available[found_name]
        try:
            df = pd.read_csv(csv_path)
        except Exception as e:
            print(f"  Error reading {csv_path}: {e}")
            failed += 1
            continue
        
        result = plot_scatter_cv(df, found_name, output_dir, args.format, args.dpi, 
                                  add_f_plot=not args.no_f)
        
        if result:
            results.append(result)
            success += 1
        else:
            failed += 1
    
    # 汇总
    print("\n" + "=" * 70)
    print(f"Completed: {success} success, {failed} failed")
    print(f"Output directory: {output_dir}")
    print("=" * 70)
    
    # 生成汇总表格
    if results:
        summary_file = output_dir / 'scatter_cv_summary.csv'
        rows = []
        for r in results:
            row = {
                'structure': r['structure'],
                'T_boundary_K': r['T_boundary'],
                'Cv_peak_meV_K': r['Cv_peak'],
            }
            for i, (phase, fit) in enumerate(r['partitions'].items()):
                row[f'phase_{i+1}'] = phase
                row[f'Cv_{i+1}_meV_K'] = fit['Cv']
                row[f'R2_{i+1}'] = fit['R2']
            rows.append(row)
        
        df_summary = pd.DataFrame(rows)
        df_summary.to_csv(summary_file, index=False)
        print(f"Summary saved: {summary_file}")


if __name__ == '__main__':
    main()
