#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 6.1.1: 分区热容拟合图 - 论文出图专用

从 step6_1 提取的核心绘图功能，生成适合论文发表的分区热容拟合图：
- 整体拟合 vs 分区拟合对比
- 按温度平均的数据点（带误差棒）
- 多数投票规则避免温度交叉

输入: step6_1 生成的聚类结果 CSV 文件
输出: 高质量论文图 (PNG/PDF)

用法:
  python step6_1_1_partition_cv_plot.py --structure Pt8sn6
  python step6_1_1_partition_cv_plot.py --structure Air86 --format pdf
  python step6_1_1_partition_cv_plot.py --structure all --dpi 600
  python step6_1_1_partition_cv_plot.py --list

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
from matplotlib.gridspec import GridSpec
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


# 相态颜色配置
PHASE_COLORS = {
    'Solid': '#3498db',      # 蓝色 - 固态
    'Pre-melting': '#f39c12', # 橙色 - 预熔化
    'Liquid': '#e74c3c',     # 红色 - 液态
    'Phase_1': '#3498db',
    'Phase_2': '#e74c3c',
    'Phase_3': '#f39c12',
}

# 载体热容 (meV/K)
CV_SUPPORT = 38.2151


def find_clustering_results(base_dir='results/step6_1_clustering'):
    """查找所有可用的聚类结果"""
    results = {}
    # 使用实际的文件命名模式
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


def load_cluster_data(csv_path):
    """加载聚类结果数据"""
    try:
        df = pd.read_csv(csv_path)
        return df
    except Exception as e:
        print(f"  错误: 无法读取 {csv_path}: {e}")
        return None


def plot_partition_cv(df, structure_name, output_dir, output_format='png', dpi=300):
    """
    绘制分区热容拟合图（论文出图专用）
    
    核心逻辑：
    1. 按温度分组计算团簇能量平均值和标准差
    2. 使用多数投票规则将每个温度分配给唯一的相态
    3. 对每个相态的专属温度点进行线性拟合
    4. 绘制整体拟合线 vs 分区拟合线对比
    """
    
    print(f"\n>>> 绘制 {structure_name} 分区热容图...")
    
    # 检查必要列
    required_cols = ['temp', 'avg_energy', 'phase_clustered']
    if not all(col in df.columns for col in required_cols):
        print(f"  错误: 缺少必要列 {required_cols}")
        return None
    
    # 判断是否是 Air 系列（气相团簇）
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
            print(f"  [载体数据] Cv_support={slope_support*1000:.4f} meV/K, R²={R2_support:.6f}")
        else:
            slope_support = CV_SUPPORT / 1000  # meV/K -> eV/K
            T_min = df['temp'].min()
            E_total_min = df[df['temp'] == T_min]['avg_energy'].mean()
            intercept_support = E_total_min * 0.9 - slope_support * T_min
            print(f"  [警告] 使用默认Cv_support估算载体能量")
    
    # ========== 1. 按温度分组计算团簇能量 ==========
    temp_groups = df.groupby('temp')
    temps_unique = []
    E_cluster_mean = []
    E_cluster_std = []
    
    for temp, group in temp_groups:
        if is_air_system:
            E_cluster = group['avg_energy'].values
        else:
            E_support = slope_support * temp + intercept_support
            E_cluster = group['avg_energy'].values - E_support
        
        temps_unique.append(temp)
        E_cluster_mean.append(np.mean(E_cluster))
        E_cluster_std.append(np.std(E_cluster))
    
    temps_unique = np.array(temps_unique)
    E_cluster_mean = np.array(E_cluster_mean)
    E_cluster_std = np.array(E_cluster_std)
    
    # 计算相对能量（相对于最低温度）
    E_cluster_ref = E_cluster_mean.min()
    E_cluster_mean_rel = E_cluster_mean - E_cluster_ref
    
    # ========== 2. 多数投票确定每个温度的专属相态 ==========
    temp_to_partition = {}
    print(f"\n  多数投票温度分配:")
    
    for temp in temps_unique:
        df_temp = df[df['temp'] == temp]
        partition_counts = df_temp['phase_clustered'].value_counts()
        dominant_partition = partition_counts.idxmax()
        temp_to_partition[temp] = dominant_partition
        print(f"    T={temp:4.0f}K: {dict(partition_counts)} → {dominant_partition}")
    
    # ========== 3. 整体拟合 ==========
    if len(temps_unique) < 3:
        print(f"  错误: 温度点不足 ({len(temps_unique)} < 3)")
        return None
    
    slope_overall, intercept_overall, r_value_overall, _, std_err_overall = linregress(
        temps_unique, E_cluster_mean_rel)
    R2_overall = r_value_overall ** 2
    Cv_overall = slope_overall * 1000  # meV/K
    Cv_overall_err = std_err_overall * 1000
    
    print(f"\n  整体拟合: Cv={Cv_overall:.4f}±{Cv_overall_err:.4f} meV/K, R²={R2_overall:.4f}")
    
    # ========== 4. 分区拟合 ==========
    phases = df['phase_clustered'].unique()
    phase_fits = {}
    
    for phase in phases:
        phase_temps = [temp for temp, part in temp_to_partition.items() if part == phase]
        phase_temps = sorted(phase_temps)
        
        if len(phase_temps) >= 2:
            mask = np.isin(temps_unique, phase_temps)
            T_phase = temps_unique[mask]
            E_phase_rel = E_cluster_mean_rel[mask]
            E_phase_std = E_cluster_std[mask]
            
            slope_ph, intercept_ph, r_value_ph, _, std_err_ph = linregress(T_phase, E_phase_rel)
            R2_ph = r_value_ph ** 2
            Cv_ph = slope_ph * 1000
            Cv_ph_err = std_err_ph * 1000
            
            phase_fits[phase] = {
                'slope': slope_ph,
                'intercept': intercept_ph,
                'R2': R2_ph,
                'Cv': Cv_ph,
                'Cv_err': Cv_ph_err,
                'n_temps': len(T_phase),
                'T_range': (T_phase.min(), T_phase.max()),
                'T_data': T_phase,
                'E_data': E_phase_rel,
                'E_std': E_phase_std
            }
            
            print(f"  {phase}: Cv={Cv_ph:.4f}±{Cv_ph_err:.4f} meV/K, R²={R2_ph:.4f}, "
                  f"n={len(T_phase)}, T={T_phase.min():.0f}-{T_phase.max():.0f}K")
    
    # ========== 5. 绘制简洁的双Y轴图 ==========
    fig, ax1 = plt.subplots(figsize=(8, 6))
    
    # ----- 左Y轴: 能量-温度数据点（带误差棒）和拟合线 -----
    # 绘制数据点（带误差棒）
    ax1.errorbar(temps_unique, E_cluster_mean_rel, yerr=E_cluster_std,
                 fmt='o', markersize=7, color='black', 
                 ecolor='gray', elinewidth=1.5, capsize=3, capthick=1.5,
                 zorder=5, label='Data')
    
    # 绘制拟合线（黑色）
    phases_sorted = sorted(phase_fits.keys())
    for phase in phases_sorted:
        fit = phase_fits[phase]
        T_phase_fit = np.linspace(fit['T_range'][0], fit['T_range'][1], 50)
        E_phase_fit = fit['slope'] * T_phase_fit + fit['intercept']
        ax1.plot(T_phase_fit, E_phase_fit, '-', color='black', linewidth=2, zorder=4)
    
    # 连接两个分区之间的数据点（实线连接实际数据点，而非拟合线）
    if len(phases_sorted) >= 2:
        fit1 = phase_fits[phases_sorted[0]]
        fit2 = phase_fits[phases_sorted[1]]
        # 分区1的最后一个数据点
        T1_end = fit1['T_range'][1]
        idx1 = np.where(temps_unique == T1_end)[0]
        if len(idx1) > 0:
            E1_end = E_cluster_mean_rel[idx1[0]]
        else:
            E1_end = fit1['slope'] * T1_end + fit1['intercept']
        # 分区2的第一个数据点
        T2_start = fit2['T_range'][0]
        idx2 = np.where(temps_unique == T2_start)[0]
        if len(idx2) > 0:
            E2_start = E_cluster_mean_rel[idx2[0]]
        else:
            E2_start = fit2['slope'] * T2_start + fit2['intercept']
        # 用实线连接两个数据点
        ax1.plot([T1_end, T2_start], [E1_end, E2_start], '-', color='black', linewidth=2, zorder=4)
    
    ax1.set_xlabel('Temperature (K)', fontsize=13, fontweight='bold')
    ax1.set_ylabel('Total Energy (eV)', fontsize=13, fontweight='bold')
    ax1.tick_params(axis='both', labelsize=11)
    
    # ----- 右Y轴: 热容曲线 -----
    ax2 = ax1.twinx()
    
    # 分界温度和热容值（用于导出）
    T_boundary = None
    Cv1 = None
    Cv2 = None
    Cv_peak = None
    
    if len(phases_sorted) >= 2:
        # 找到分区边界温度
        phase1_temps = [t for t, p in temp_to_partition.items() if p == phases_sorted[0]]
        phase2_temps = [t for t, p in temp_to_partition.items() if p == phases_sorted[1]]
        
        if phase1_temps and phase2_temps:
            T1_last = max(phase1_temps)   # 分区1最后一个温度
            T2_first = min(phase2_temps)  # 分区2第一个温度
            T_boundary = (T1_last + T2_first) / 2
            print(f"\n  分界温度: {T_boundary:.0f} K (过渡区: {T1_last:.0f}-{T2_first:.0f}K)")
            
            Cv1 = phase_fits[phases_sorted[0]]['Cv']
            Cv2 = phase_fits[phases_sorted[1]]['Cv']
            
            # 计算过渡区热容（数值微分）
            idx1 = np.where(temps_unique == T1_last)[0]
            idx2 = np.where(temps_unique == T2_first)[0]
            if len(idx1) > 0 and len(idx2) > 0:
                E1 = E_cluster_mean_rel[idx1[0]]
                E2 = E_cluster_mean_rel[idx2[0]]
                Cv_transition = (E2 - E1) / (T2_first - T1_last) * 1000  # meV/K
            else:
                Cv_transition = (Cv1 + Cv2) / 2
            
            # 判断是否存在热容峰
            has_peak = Cv_transition > max(Cv1, Cv2)
            
            if has_peak:
                Cv_peak = Cv_transition
                print(f"  ★ 存在热容峰: Cv_peak={Cv_peak:.2f} meV/K (过渡区)")
                print(f"  热容: Cv1={Cv1:.2f}, Cv_peak={Cv_peak:.2f}, Cv2={Cv2:.2f} meV/K")
                
                # 绘制带平滑峰的热容曲线（使用高斯峰 + sigmoid过渡）
                T_plot = np.linspace(temps_unique.min(), temps_unique.max(), 500)
                Cv_plot = np.zeros_like(T_plot)
                
                # 峰的宽度参数
                sigma = (T2_first - T1_last) / 2  # 高斯宽度
                
                for i, T in enumerate(T_plot):
                    # 基线：sigmoid 从 Cv1 过渡到 Cv2
                    transition = 1 / (1 + np.exp(-(T - T_boundary) / (sigma * 0.5)))
                    baseline = Cv1 + (Cv2 - Cv1) * transition
                    
                    # 高斯峰叠加
                    gaussian = (Cv_peak - baseline) * np.exp(-0.5 * ((T - T_boundary) / sigma)**2)
                    Cv_plot[i] = baseline + gaussian
                
                ax2.plot(T_plot, Cv_plot, 'r-', linewidth=2, zorder=3)
                
                # 构建导出数据（关键点）
                T_cv = np.array([temps_unique.min(), T1_last, T_boundary, T2_first, temps_unique.max()])
                Cv_curve = np.array([Cv1, Cv1, Cv_peak, Cv2, Cv2])
            else:
                print(f"  热容: Cv1={Cv1:.2f} meV/K, Cv2={Cv2:.2f} meV/K (无峰)")
                
                # 绘制阶梯形热容曲线（无峰）
                ax2.plot([temps_unique.min(), T_boundary], [Cv1, Cv1], 'r-', linewidth=2, zorder=3)
                ax2.plot([T_boundary, T_boundary], [Cv1, Cv2], 'r--', linewidth=1.5, zorder=3)
                ax2.plot([T_boundary, temps_unique.max()], [Cv2, Cv2], 'r-', linewidth=2, zorder=3)
                
                T_cv = np.array([temps_unique.min(), T_boundary - 0.1, T_boundary, T_boundary + 0.1, temps_unique.max()])
                Cv_curve = np.array([Cv1, Cv1, (Cv1 + Cv2) / 2, Cv2, Cv2])
    else:
        Cv_single = list(phase_fits.values())[0]['Cv']
        T_cv = np.array([temps_unique.min(), temps_unique.max()])
        Cv_curve = np.array([Cv_single, Cv_single])
        ax2.plot(T_cv, Cv_curve, 'r-', linewidth=2, zorder=3)
        Cv1 = Cv_single
        Cv2 = Cv_single
    
    ax2.set_ylabel('Cv (meV/K)', fontsize=13, fontweight='bold', color='red')
    ax2.tick_params(axis='y', labelcolor='red', labelsize=11, color='red')
    ax2.spines['right'].set_color('red')
    
    # 设置Y轴范围（考虑峰值）
    cv_values = [Cv1, Cv2] if Cv1 and Cv2 else list(Cv_curve)
    if Cv_peak:
        cv_values.append(Cv_peak)
    cv_min = min(cv_values) * 0.85
    cv_max = max(cv_values) * 1.1
    ax2.set_ylim(cv_min, cv_max)
    
    ax1.set_title(f'{structure_name}', fontsize=14, fontweight='bold', pad=10)
    
    plt.tight_layout()
    
    # 保存图片
    output_file = Path(output_dir) / f'{structure_name}_partition_cv.{output_format}'
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close()
    
    print(f"\n  图已保存: {output_file}")
    
    # ========== 6. 导出数据供 Origin 使用 ==========
    # 导出能量数据
    df_energy = pd.DataFrame({
        'Temperature_K': temps_unique,
        'Energy_eV': E_cluster_mean_rel,
        'Energy_std_eV': E_cluster_std,
        'Partition': [temp_to_partition.get(t, 'unknown') for t in temps_unique]
    })
    energy_csv = Path(output_dir) / f'{structure_name}_energy_data.csv'
    df_energy.to_csv(energy_csv, index=False)
    print(f"  能量数据已导出: {energy_csv}")
    
    # 导出热容数据（阶梯函数关键点）
    df_cv = pd.DataFrame({
        'Temperature_K': T_cv,
        'Cv_meV_K': Cv_curve
    })
    cv_csv = Path(output_dir) / f'{structure_name}_cv_curve.csv'
    df_cv.to_csv(cv_csv, index=False)
    print(f"  热容曲线已导出: {cv_csv}")
    
    # 导出拟合参数汇总
    fit_summary = {
        'structure': structure_name,
        'T_boundary_K': T_boundary,
        'Cv_overall_meV_K': Cv_overall,
        'Cv_overall_err': Cv_overall_err,
        'R2_overall': R2_overall,
    }
    for i, (phase, fit) in enumerate(phase_fits.items()):
        fit_summary[f'phase_{i+1}_name'] = phase
        fit_summary[f'phase_{i+1}_Cv_meV_K'] = fit['Cv']
        fit_summary[f'phase_{i+1}_Cv_err'] = fit['Cv_err']
        fit_summary[f'phase_{i+1}_R2'] = fit['R2']
        fit_summary[f'phase_{i+1}_T_min_K'] = fit['T_range'][0]
        fit_summary[f'phase_{i+1}_T_max_K'] = fit['T_range'][1]
        fit_summary[f'phase_{i+1}_slope_eV_K'] = fit['slope']
        fit_summary[f'phase_{i+1}_intercept_eV'] = fit['intercept']
    
    fit_csv = Path(output_dir) / f'{structure_name}_fit_params.csv'
    pd.DataFrame([fit_summary]).to_csv(fit_csv, index=False)
    print(f"  拟合参数已导出: {fit_csv}")
    
    # 返回拟合结果
    return {
        'structure': structure_name,
        'overall': {'Cv': Cv_overall, 'Cv_err': Cv_overall_err, 'R2': R2_overall},
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
    
    print(f"\n🔵 气相团簇 ({len(air_series)}): {', '.join(air_series) if air_series else '无'}")
    print(f"🟢 Pt6系列 ({len(pt6_series)}): {', '.join(sorted(pt6_series)) if pt6_series else '无'}")
    print(f"🟢 Pt8系列 ({len(pt8_series)}): {', '.join(sorted(pt8_series)) if pt8_series else '无'}")
    print(f"🟠 含氧团簇 ({len(oxide_series)}): {', '.join(sorted(oxide_series)) if oxide_series else '无'}")
    if other:
        print(f"⚪ 其他 ({len(other)}): {', '.join(sorted(other))}")
    
    print(f"\n总计: {len(results)} 个结构")
    print("=" * 60)
    
    return results


def parse_args():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description='分区热容拟合图 - 论文出图专用',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例:
  %(prog)s --structure Pt8sn6              # 单个结构
  %(prog)s --structure Air86 --format pdf  # 输出PDF
  %(prog)s --structure all --dpi 600       # 所有结构，高分辨率
  %(prog)s --list                          # 列出可用结构
        '''
    )
    
    parser.add_argument('--structure', '-s', type=str, default=None,
                        help='结构名称 (如 Pt8sn6, Air86) 或 "all" 处理所有')
    parser.add_argument('--list', '-l', action='store_true',
                        help='列出所有可用结构')
    parser.add_argument('--format', '-f', type=str, default='png',
                        choices=['png', 'pdf', 'svg', 'eps'],
                        help='输出格式 (默认: png)')
    parser.add_argument('--dpi', type=int, default=300,
                        help='输出分辨率 (默认: 300)')
    parser.add_argument('--output-dir', '-o', type=str, 
                        default='results/step6_1_1_partition_cv',
                        help='输出目录')
    
    return parser.parse_args()


def main():
    """主函数"""
    args = parse_args()
    
    print("=" * 70)
    print("Step 6.1.1: 分区热容拟合图 - 论文出图专用")
    print("=" * 70)
    
    # 列出可用结构
    if args.list:
        list_available_structures()
        return
    
    if args.structure is None:
        print("错误: 请指定 --structure 或使用 --list 查看可用结构")
        return
    
    # 创建输出目录
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 获取可用结构
    available = find_clustering_results()
    
    if args.structure.lower() == 'all':
        structures = list(available.keys())
        print(f"\n处理所有 {len(structures)} 个结构...")
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
            print(f"\n警告: 未找到结构 '{structure}'")
            failed += 1
            continue
        
        csv_path = available[found_name]
        df = load_cluster_data(csv_path)
        
        if df is None:
            failed += 1
            continue
        
        result = plot_partition_cv(df, found_name, output_dir, 
                                   args.format, args.dpi)
        
        if result:
            results.append(result)
            success += 1
        else:
            failed += 1
    
    # 汇总
    print("\n" + "=" * 70)
    print(f"处理完成: 成功 {success}, 失败 {failed}")
    print(f"输出目录: {output_dir}")
    print("=" * 70)
    
    # 生成汇总表格
    if results:
        summary_file = output_dir / 'partition_cv_summary.csv'
        rows = []
        for r in results:
            row = {
                'structure': r['structure'],
                'Cv_overall': r['overall']['Cv'],
                'Cv_overall_err': r['overall']['Cv_err'],
                'R2_overall': r['overall']['R2'],
            }
            for i, (phase, fit) in enumerate(r['partitions'].items()):
                row[f'phase_{i+1}'] = phase
                row[f'Cv_{i+1}'] = fit['Cv']
                row[f'Cv_{i+1}_err'] = fit['Cv_err']
                row[f'R2_{i+1}'] = fit['R2']
            rows.append(row)
        
        df_summary = pd.DataFrame(rows)
        df_summary.to_csv(summary_file, index=False)
        print(f"汇总表已保存: {summary_file}")


if __name__ == '__main__':
    main()
