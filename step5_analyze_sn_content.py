#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
分析Sn含量对扩散行为的影响 (Pt8SnX系列)
===========================================

分析内容:
1. D vs Sn含量 (不同温度)
2. D vs 温度 (不同Sn含量)
3. Arrhenius图 (ln(D) vs 1/T)
4. 活化能计算
5. 统计对比表
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import stats
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# ===== 閰嶇疆 =====
BASE_DIR = Path(__file__).parent
DATA_FILE = BASE_DIR / 'results' / 'ensemble_D_analysis' / 'ensemble_D_values.csv'
OUTPUT_DIR = BASE_DIR / 'results' / 'sn_content_analysis'
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

# 颜色方案
COLORS = {
    'Pt': '#E74C3C',
    'Sn': '#3498DB',
    'PtSn': '#2ECC71'
}

# 温度列表
TEMPS = [200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100]

plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False


def extract_sn_content(composition):
    """
    从组成名称提取Sn含量
    
    Examples:
    ---------
    pt8sn0-2-best -> 0
    pt8sn5-1-best -> 5
    pt8sn10-2-best -> 10
    """
    import re
    match = re.search(r'pt8sn(\d+)', composition, re.IGNORECASE)
    if match:
        return int(match.group(1))
    return None


def load_and_prepare_data():
    """加载并预处理数据"""
    print("\n[*] Loading data...")
    df = pd.read_csv(DATA_FILE)
    
    # 添加Sn含量列
    df['sn_content'] = df['composition'].apply(extract_sn_content)
    
    # 过滤掉无法识别Sn含量的数据
    df = df[df['sn_content'].notna()]
    
    # 添加倒数温度 (用于Arrhenius图)
    df['inv_T_1000K'] = 1000.0 / df['temp_K']  # 1000/T (K^-1)
    
    # 添加ln(D) (只对正D值)
    df['ln_D'] = np.nan
    positive_D = df['D_cm2_s'] > 0
    df.loc[positive_D, 'ln_D'] = np.log(df.loc[positive_D, 'D_cm2_s'])
    
    print(f"  [OK] Loaded {len(df)} data points")
    print(f"  Sn content range: {df['sn_content'].min():.0f} - {df['sn_content'].max():.0f}")
    print(f"  Temperature range: {df['temp_K'].min():.0f} - {df['temp_K'].max():.0f} K")
    print(f"  Elements: {', '.join(df['element'].unique())}")
    
    return df


def plot_D_vs_sn_content(df, output_dir):
    """
    绘制 D vs Sn含量 (不同温度)
    每个温度一个子图,显示3个元素
    """
    print("\n[*] Plotting D vs Sn content...")
    
    # 只选择部分代表性温度
    selected_temps = [300, 500, 700, 800, 900, 1000, 1100]
    
    n_temps = len(selected_temps)
    n_cols = 3
    n_rows = (n_temps + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 5*n_rows))
    axes = axes.flatten()
    
    for idx, temp in enumerate(selected_temps):
        ax = axes[idx]
        df_temp = df[df['temp_K'] == temp]
        
        if len(df_temp) == 0:
            ax.text(0.5, 0.5, 'No Data', ha='center', va='center',
                   transform=ax.transAxes, fontsize=14)
            ax.set_title(f'{temp}K', fontsize=12, fontweight='bold')
            continue
        
        # 按元素绘制
        for element in ['Pt', 'Sn', 'PtSn']:
            df_elem = df_temp[df_temp['element'] == element]
            
            if len(df_elem) == 0:
                continue
            
            # 按Sn含量排序
            df_elem = df_elem.sort_values('sn_content')
            
            # 绘制数据点和曲线
            ax.plot(df_elem['sn_content'], df_elem['D_cm2_s'], 
                   'o-', label=element, color=COLORS[element],
                   markersize=8, linewidth=2, alpha=0.7)
        
        ax.set_xlabel('Sn Content (x in Pt8Snx)', fontsize=11, fontweight='bold')
        ax.set_ylabel('D (cm²/s)', fontsize=11, fontweight='bold')
        ax.set_title(f'{temp}K', fontsize=12, fontweight='bold')
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        ax.set_xlim(-0.5, 10.5)
        ax.set_xticks(range(0, 11))
    
    # 隐藏多余子图
    for idx in range(n_temps, len(axes)):
        axes[idx].axis('off')
    
    plt.suptitle('Diffusion Coefficient vs Sn Content\n'
                'Pt8SnX Series at Different Temperatures',
                fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    
    output_file = output_dir / 'D_vs_SnContent_by_temperature.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    print(f"  [OK] Saved: {output_file.name}")


def plot_D_vs_temperature_comparison(df, output_dir):
    """
    绘制 D vs 温度 (不同Sn含量的对比)
    3个子图分别对应Pt, Sn, PtSn
    """
    print("\n[*] Plotting D vs Temperature (Sn content comparison)...")
    
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))
    
    # 定义颜色映射 (Sn含量: 0-10)
    cmap = plt.cm.viridis
    sn_contents = sorted(df['sn_content'].unique())
    colors_sn = {sn: cmap(i/len(sn_contents)) for i, sn in enumerate(sn_contents)}
    
    for idx, element in enumerate(['Pt', 'Sn', 'PtSn']):
        ax = axes[idx]
        df_elem = df[df['element'] == element]
        
        if len(df_elem) == 0:
            ax.text(0.5, 0.5, 'No Data', ha='center', va='center',
                   transform=ax.transAxes, fontsize=14)
            ax.set_title(f'{element}', fontsize=12, fontweight='bold')
            continue
        
        # 按Sn含量分组绘制
        for sn in sn_contents:
            df_sn = df_elem[df_elem['sn_content'] == sn]
            df_sn = df_sn.sort_values('temp_K')
            
            if len(df_sn) == 0:
                continue
            
            ax.plot(df_sn['temp_K'], df_sn['D_cm2_s'],
                   'o-', label=f'Sn{int(sn)}', 
                   color=colors_sn[sn],
                   markersize=6, linewidth=2, alpha=0.7)
        
        ax.set_xlabel('Temperature (K)', fontsize=11, fontweight='bold')
        ax.set_ylabel('D (cm²/s)', fontsize=11, fontweight='bold')
        ax.set_title(f'{element}', fontsize=12, fontweight='bold', color=COLORS[element])
        ax.legend(fontsize=8, ncol=2, loc='upper left')
        ax.grid(True, alpha=0.3)
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        ax.set_xlim(150, 1150)
    
    plt.suptitle('Diffusion Coefficient vs Temperature\n'
                'Comparison of Different Sn Contents (Pt8SnX)',
                fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    
    output_file = output_dir / 'D_vs_Temperature_SnContent_comparison.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    print(f"  [OK] Saved: {output_file.name}")


def plot_arrhenius(df, output_dir):
    """
    绘制Arrhenius图: ln(D) vs 1000/T
    并计算活化能
    """
    print("\n[*] Plotting Arrhenius diagram...")
    
    # 只使用高温数据 (≥700K) 和正D值
    df_high_temp = df[(df['temp_K'] >= 700) & (df['D_cm2_s'] > 0)].copy()
    
    if len(df_high_temp) == 0:
        print("  [!] No high-temperature data available")
        return None
    
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))
    
    # 定义颜色映射
    cmap = plt.cm.viridis
    sn_contents = sorted(df_high_temp['sn_content'].unique())
    colors_sn = {sn: cmap(i/len(sn_contents)) for i, sn in enumerate(sn_contents)}
    
    # 存储活化能结果
    activation_energies = []
    
    for idx, element in enumerate(['Pt', 'Sn', 'PtSn']):
        ax = axes[idx]
        df_elem = df_high_temp[df_high_temp['element'] == element]
        
        if len(df_elem) == 0:
            ax.text(0.5, 0.5, 'No Data', ha='center', va='center',
                   transform=ax.transAxes, fontsize=14)
            ax.set_title(f'{element}', fontsize=12, fontweight='bold')
            continue
        
        # 按Sn含量分组
        for sn in sn_contents:
            df_sn = df_elem[df_elem['sn_content'] == sn]
            df_sn = df_sn.sort_values('temp_K')
            
            if len(df_sn) < 3:  # 至少需要3个点才能拟合
                continue
            
            # 绘制数据点
            ax.plot(df_sn['inv_T_1000K'], df_sn['ln_D'],
                   'o', label=f'Sn{int(sn)}', 
                   color=colors_sn[sn],
                   markersize=8, alpha=0.7)
            
            # 线性拟合 (Arrhenius方程: ln(D) = ln(D0) - Ea/(R*T))
            try:
                slope, intercept, r_value, p_value, std_err = stats.linregress(
                    df_sn['inv_T_1000K'], df_sn['ln_D']
                )
                
                # 绘制拟合线
                x_fit = np.array([df_sn['inv_T_1000K'].min(), 
                                 df_sn['inv_T_1000K'].max()])
                y_fit = slope * x_fit + intercept
                ax.plot(x_fit, y_fit, '--', color=colors_sn[sn], 
                       linewidth=1.5, alpha=0.5)
                
                # 计算活化能 (Ea = -slope * R * 1000)
                # R = 8.314 J/(mol·K)
                Ea_kJ_mol = -slope * 8.314  # kJ/mol
                
                activation_energies.append({
                    'element': element,
                    'sn_content': int(sn),
                    'Ea_kJ_mol': Ea_kJ_mol,
                    'r2': r_value**2,
                    'n_points': len(df_sn),
                    'T_range': f"{df_sn['temp_K'].min()}-{df_sn['temp_K'].max()}K"
                })
                
            except Exception as e:
                continue
        
        ax.set_xlabel('1000/T (K⁻¹)', fontsize=11, fontweight='bold')
        ax.set_ylabel('ln(D) [D in cm²/s]', fontsize=11, fontweight='bold')
        ax.set_title(f'{element}', fontsize=12, fontweight='bold', color=COLORS[element])
        ax.legend(fontsize=8, ncol=2, loc='best')
        ax.grid(True, alpha=0.3)
        ax.invert_xaxis()  # 高温在右
    
    plt.suptitle('Arrhenius Plot: ln(D) vs 1000/T\n'
                'High Temperature Region (≥700K) | Dashed lines: Linear fits',
                fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    
    output_file = output_dir / 'Arrhenius_plot_SnContent_comparison.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    print(f"  [OK] Saved: {output_file.name}")
    
    return activation_energies


def plot_activation_energy(activation_energies, output_dir):
    """绘制活化能 vs Sn含量"""
    if not activation_energies:
        print("  [!] No activation energy data available")
        return
    
    print("\n[*] Plotting activation energy vs Sn content...")
    
    df_ea = pd.DataFrame(activation_energies)
    
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    
    for element in ['Pt', 'Sn', 'PtSn']:
        df_elem = df_ea[df_ea['element'] == element]
        
        if len(df_elem) == 0:
            continue
        
        df_elem = df_elem.sort_values('sn_content')
        
        ax.plot(df_elem['sn_content'], df_elem['Ea_kJ_mol'],
               'o-', label=element, color=COLORS[element],
               markersize=10, linewidth=2.5, alpha=0.7)
        
        # 标注R²值
        for _, row in df_elem.iterrows():
            if row['r2'] >= 0.8:  # 只标注高质量拟合
                ax.annotate(f"R²={row['r2']:.2f}", 
                          xy=(row['sn_content'], row['Ea_kJ_mol']),
                          xytext=(5, 5), textcoords='offset points',
                          fontsize=7, alpha=0.6)
    
    ax.set_xlabel('Sn Content (x in Pt8Snx)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Activation Energy (kJ/mol)', fontsize=12, fontweight='bold')
    ax.set_title('Activation Energy vs Sn Content\n'
                'Calculated from Arrhenius Plot (T≥700K)',
                fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(-0.5, 10.5)
    ax.set_xticks(range(0, 11))
    
    plt.tight_layout()
    
    output_file = output_dir / 'ActivationEnergy_vs_SnContent.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    print(f"  [OK] Saved: {output_file.name}")
    
    return df_ea


def generate_statistics_table(df, activation_energies, output_dir):
    """生成统计对比表"""
    print("\n[*] Generating statistics table...")
    
    report_file = output_dir / 'SnContent_effect_report.txt'
    
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write("="*80 + "\n")
        f.write("Sn Content Effect on Diffusion Behavior - Statistical Analysis\n")
        f.write("="*80 + "\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"System: Pt8SnX (X = 0-10)\n")
        f.write("="*80 + "\n\n")
        
        # 1. 按Sn含量统计平均D值
        f.write("1. Average Diffusion Coefficient by Sn Content\n")
        f.write("-"*80 + "\n")
        f.write(f"{'Sn':<5} {'Element':<8} {'<D> (cm²/s)':<15} {'Std':<15} {'N':<5} {'T_range':<15}\n")
        f.write("-"*80 + "\n")
        
        for sn in sorted(df['sn_content'].unique()):
            for element in ['Pt', 'Sn', 'PtSn']:
                df_subset = df[(df['sn_content'] == sn) & 
                              (df['element'] == element) & 
                              (df['D_cm2_s'] > 0)]
                
                if len(df_subset) == 0:
                    continue
                
                mean_D = df_subset['D_cm2_s'].mean()
                std_D = df_subset['D_cm2_s'].std()
                n = len(df_subset)
                T_range = f"{df_subset['temp_K'].min():.0f}-{df_subset['temp_K'].max():.0f}K"
                
                f.write(f"{int(sn):<5} {element:<8} {mean_D:>14.2e} {std_D:>14.2e} {n:<5} {T_range:<15}\n")
        
        f.write("\n" + "="*80 + "\n\n")
        
        # 2. 按温度统计 (选择代表性温度)
        f.write("2. D Value Comparison at Selected Temperatures\n")
        f.write("-"*80 + "\n")
        
        for temp in [700, 800, 900, 1000, 1100]:
            f.write(f"\nTemperature: {temp}K\n")
            f.write(f"{'Sn':<5} {'Pt (cm²/s)':<15} {'Sn (cm²/s)':<15} {'PtSn (cm²/s)':<15}\n")
            f.write("-"*80 + "\n")
            
            df_temp = df[(df['temp_K'] == temp) & (df['D_cm2_s'] > 0)]
            
            for sn in sorted(df_temp['sn_content'].unique()):
                df_sn = df_temp[df_temp['sn_content'] == sn]
                
                D_pt = df_sn[df_sn['element'] == 'Pt']['D_cm2_s'].values
                D_sn = df_sn[df_sn['element'] == 'Sn']['D_cm2_s'].values
                D_ptsn = df_sn[df_sn['element'] == 'PtSn']['D_cm2_s'].values
                
                D_pt_str = f"{D_pt[0]:.2e}" if len(D_pt) > 0 else "N/A"
                D_sn_str = f"{D_sn[0]:.2e}" if len(D_sn) > 0 else "N/A"
                D_ptsn_str = f"{D_ptsn[0]:.2e}" if len(D_ptsn) > 0 else "N/A"
                
                f.write(f"{int(sn):<5} {D_pt_str:<15} {D_sn_str:<15} {D_ptsn_str:<15}\n")
        
        f.write("\n" + "="*80 + "\n\n")
        
        # 3. 活化能统计
        if activation_energies:
            f.write("3. Activation Energy (High Temperature: T≥700K)\n")
            f.write("-"*80 + "\n")
            f.write(f"{'Sn':<5} {'Element':<8} {'Ea (kJ/mol)':<15} {'R²':<8} {'N':<5} {'T_range':<15}\n")
            f.write("-"*80 + "\n")
            
            for ea in activation_energies:
                f.write(f"{ea['sn_content']:<5} {ea['element']:<8} "
                       f"{ea['Ea_kJ_mol']:>14.2f} {ea['r2']:>7.3f} "
                       f"{ea['n_points']:<5} {ea['T_range']:<15}\n")
            
            f.write("\n" + "="*80 + "\n\n")
        
        # 4. 关键发现
        f.write("4. Key Findings\n")
        f.write("-"*80 + "\n")
        
        # 找出D值最大和最小的Sn含量
        df_high_temp = df[(df['temp_K'] >= 800) & (df['D_cm2_s'] > 0)]
        
        for element in ['Pt', 'Sn', 'PtSn']:
            df_elem = df_high_temp[df_high_temp['element'] == element]
            
            if len(df_elem) == 0:
                continue
            
            # 按Sn含量分组计算平均D
            D_by_sn = df_elem.groupby('sn_content')['D_cm2_s'].mean()
            
            if len(D_by_sn) > 0:
                max_sn = D_by_sn.idxmax()
                min_sn = D_by_sn.idxmin()
                
                f.write(f"\n{element}:\n")
                f.write(f"  Highest D at high T (≥800K): Sn{int(max_sn)} "
                       f"(D_avg = {D_by_sn[max_sn]:.2e} cm²/s)\n")
                f.write(f"  Lowest D at high T (≥800K):  Sn{int(min_sn)} "
                       f"(D_avg = {D_by_sn[min_sn]:.2e} cm²/s)\n")
        
        f.write("\n" + "="*80 + "\n")
        f.write("End of Report\n")
        f.write("="*80 + "\n")
    
    print(f"  [OK] Saved: {report_file.name}")


def main():
    print("\n" + "="*80)
    print("Analyzing Sn Content Effect on Diffusion Behavior")
    print("="*80)
    print(f"Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*80)
    
    # 1. 加载数据
    df = load_and_prepare_data()
    
    # 2. 绘图分析
    print("\n" + "="*80)
    print("Generating Visualizations")
    print("="*80)
    
    # 2.1 D vs Sn含量 (不同温度)
    plot_D_vs_sn_content(df, OUTPUT_DIR)
    
    # 2.2 D vs 温度 (不同Sn含量)
    plot_D_vs_temperature_comparison(df, OUTPUT_DIR)
    
    # 2.3 Arrhenius图
    activation_energies = plot_arrhenius(df, OUTPUT_DIR)
    
    # 2.4 活化能 vs Sn含量
    if activation_energies:
        df_ea = plot_activation_energy(activation_energies, OUTPUT_DIR)
        
        # 保存活化能数据
        ea_file = OUTPUT_DIR / 'activation_energies.csv'
        df_ea.to_csv(ea_file, index=False, float_format='%.4f')
        print(f"  [OK] Saved: {ea_file.name}")
    
    # 3. 生成统计报告
    print("\n" + "="*80)
    print("Generating Statistical Report")
    print("="*80)
    
    generate_statistics_table(df, activation_energies, OUTPUT_DIR)
    
    # 完成
    print("\n" + "="*80)
    print("DONE!")
    print("="*80)
    print(f"Output directory: {OUTPUT_DIR}")
    print("\nGenerated files:")
    print("  - D_vs_SnContent_by_temperature.png")
    print("  - D_vs_Temperature_SnContent_comparison.png")
    print("  - Arrhenius_plot_SnContent_comparison.png")
    print("  - ActivationEnergy_vs_SnContent.png")
    print("  - activation_energies.csv")
    print("  - SnContent_effect_report.txt")
    print("="*80)


if __name__ == '__main__':
    main()
