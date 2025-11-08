#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
专门分析Cv系列(Cv-1到Cv-5)的数据
该系列为Sn8Pt6O4组分,跑了5遍,温度间隔50K
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from pathlib import Path

# 设置中文字体
matplotlib.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei']
matplotlib.rcParams['axes.unicode_minus'] = False

# 文件路径
base_dir = Path(r'd:\OneDrive\py\Cv\lin\MSD_Analysis_Collection\v3_simplified_workflow')
result_dir = base_dir / 'results' / 'energy_analysis_v2_no_filter'
heat_capacity_csv = result_dir / 'heat_capacity_per_system.csv'
energy_csv = result_dir / 'energy_per_system.csv'

# 创建输出目录
output_dir = result_dir / 'Cv_series_analysis'
output_dir.mkdir(exist_ok=True)

print("="*80)
print("Cv系列(Sn8Pt6O4)专项分析")
print("="*80)

# ========== 读取数据 ==========
df_cv = pd.read_csv(heat_capacity_csv)
df_cv = df_cv[df_cv['结构'].str.contains('Cv-', na=False)].copy()

df_energy = pd.read_csv(energy_csv)
df_energy = df_energy[df_energy['结构'].str.contains('Cv-', na=False)].copy()

print(f"\n[1] 数据概况")
print("-"*80)
print(f"热容数据: {len(df_cv)} 条记录")
print(f"能量数据: {len(df_energy)} 条记录")
print(f"结构: {sorted(df_cv['结构'].unique())}")
print(f"化学式: {df_cv['化学式'].unique()[0]}")
print(f"原子数: N = {df_cv['N_cluster'].unique()[0]}")
print(f"温度范围: {df_cv['温度'].min():.0f} - {df_cv['温度'].max():.0f} K")
print(f"温度点数: {df_cv['温度'].nunique()}")

# 各温度的重复次数
print(f"\n各温度的模拟次数:")
temp_counts = df_cv.groupby('温度').size()
for temp, count in temp_counts.items():
    print(f"  {temp:4.0f}K: {count}次")

# 各结构的温度覆盖
print(f"\n各次模拟的温度覆盖:")
for struct in sorted(df_cv['结构'].unique()):
    temps = sorted(df_cv[df_cv['结构']==struct]['温度'].unique())
    print(f"  {struct}: {len(temps)}个温度点")
    print(f"    温度: {temps}")

# ========== 统计分析 ==========
print(f"\n[2] 热容统计分析(按温度)")
print("-"*80)

stats_by_temp = df_cv.groupby('温度').agg({
    'Cv_total_meV_K': ['mean', 'std', 'min', 'max'],
    'Cv_per_atom_meV_K': ['mean', 'std', 'min', 'max']
}).round(3)

print("\n总热容 Cv_total [meV/K]:")
print(stats_by_temp['Cv_total_meV_K'])

print("\n每原子热容 Cv/N [meV/K/atom]:")
print(stats_by_temp['Cv_per_atom_meV_K'])

# 保存统计结果
stats_csv = output_dir / 'Cv_series_statistics_by_temperature.csv'
stats_by_temp.to_csv(stats_csv)
print(f"\n✓ 已保存统计结果: {stats_csv}")

# ========== 能量统计 ==========
print(f"\n[3] 能量统计分析(按温度)")
print("-"*80)

energy_stats = df_energy.groupby('温度').agg({
    '每原子能量': ['mean', 'std', 'min', 'max'],
    '相对能量': ['mean', 'std', 'min', 'max']
}).round(3)

print("\n每原子能量 [eV/atom]:")
print(energy_stats['每原子能量'])

print("\n相对能量 [eV]:")
print(energy_stats['相对能量'])

# ========== 绘图1: 热容vs温度(5次重复) ==========
print(f"\n[4] 绘制热容-温度曲线")
print("-"*80)

fig, axes = plt.subplots(2, 2, figsize=(16, 12))
fig.suptitle('Cv系列(Sn8Pt6O4) 热容分析 - 5次独立模拟', fontsize=16, fontweight='bold')

# 颜色映射
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
markers = ['o', 's', '^', 'D', 'v']

# 子图1: Cv_total vs T (各次模拟)
ax1 = axes[0, 0]
for i, struct in enumerate(sorted(df_cv['结构'].unique())):
    df_struct = df_cv[df_cv['结构'] == struct].sort_values('温度')
    ax1.plot(df_struct['温度'], df_struct['Cv_total_meV_K'], 
            marker=markers[i], color=colors[i], label=struct,
            linewidth=2, markersize=8, alpha=0.7)

# 添加平均值曲线
df_mean = df_cv.groupby('温度')['Cv_total_meV_K'].mean().reset_index()
ax1.plot(df_mean['温度'], df_mean['Cv_total_meV_K'],
        'k--', linewidth=3, label='平均值', zorder=10)

ax1.set_xlabel('温度 [K]', fontsize=12, fontweight='bold')
ax1.set_ylabel('Cv_total [meV/K]', fontsize=12, fontweight='bold')
ax1.set_title('总热容(含载体)', fontsize=13, fontweight='bold')
ax1.legend(loc='best', fontsize=10)
ax1.grid(True, alpha=0.3)

# 子图2: Cv/N vs T (各次模拟)
ax2 = axes[0, 1]
for i, struct in enumerate(sorted(df_cv['结构'].unique())):
    df_struct = df_cv[df_cv['结构'] == struct].sort_values('温度')
    ax2.plot(df_struct['温度'], df_struct['Cv_per_atom_meV_K'], 
            marker=markers[i], color=colors[i], label=struct,
            linewidth=2, markersize=8, alpha=0.7)

df_mean = df_cv.groupby('温度')['Cv_per_atom_meV_K'].mean().reset_index()
ax2.plot(df_mean['温度'], df_mean['Cv_per_atom_meV_K'],
        'k--', linewidth=3, label='平均值', zorder=10)

# 添加经典极限参考线
classical_limit = 3 * 0.0862  # 3kB per atom
ax2.axhline(classical_limit, color='red', linestyle=':', linewidth=2, 
           label=f'经典极限 ({classical_limit:.3f})')

ax2.set_xlabel('温度 [K]', fontsize=12, fontweight='bold')
ax2.set_ylabel('Cv/N [meV/K/atom]', fontsize=12, fontweight='bold')
ax2.set_title('每原子热容(表观值)', fontsize=13, fontweight='bold')
ax2.legend(loc='best', fontsize=10)
ax2.grid(True, alpha=0.3)

# 子图3: 平均值±标准差
ax3 = axes[1, 0]
df_stats = df_cv.groupby('温度').agg({
    'Cv_total_meV_K': ['mean', 'std']
}).reset_index()
df_stats.columns = ['温度', 'mean', 'std']

ax3.errorbar(df_stats['温度'], df_stats['mean'], yerr=df_stats['std'],
            fmt='o-', color='darkblue', linewidth=2, markersize=8,
            capsize=5, capthick=2, label='Mean ± Std')
ax3.fill_between(df_stats['温度'], 
                df_stats['mean']-df_stats['std'],
                df_stats['mean']+df_stats['std'],
                alpha=0.2, color='blue')

ax3.set_xlabel('温度 [K]', fontsize=12, fontweight='bold')
ax3.set_ylabel('Cv_total [meV/K]', fontsize=12, fontweight='bold')
ax3.set_title('总热容(平均值±标准差)', fontsize=13, fontweight='bold')
ax3.legend(loc='best', fontsize=10)
ax3.grid(True, alpha=0.3)

# 子图4: 相对标准差(RSD%)
ax4 = axes[1, 1]
df_rsd = df_cv.groupby('温度').agg({
    'Cv_total_meV_K': lambda x: (x.std() / x.mean() * 100) if x.mean() != 0 else 0
}).reset_index()
df_rsd.columns = ['温度', 'RSD']

ax4.bar(df_rsd['温度'], df_rsd['RSD'], width=40, 
       color='coral', edgecolor='darkred', linewidth=2, alpha=0.7)
ax4.axhline(5, color='green', linestyle='--', linewidth=2, label='5% 参考线')
ax4.axhline(10, color='orange', linestyle='--', linewidth=2, label='10% 参考线')

ax4.set_xlabel('温度 [K]', fontsize=12, fontweight='bold')
ax4.set_ylabel('相对标准差 RSD [%]', fontsize=12, fontweight='bold')
ax4.set_title('热容测量的重复性', fontsize=13, fontweight='bold')
ax4.legend(loc='best', fontsize=10)
ax4.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
fig_path = output_dir / 'Cv_series_heat_capacity_analysis.png'
plt.savefig(fig_path, dpi=300, bbox_inches='tight')
print(f"✓ 已保存图片: {fig_path}")
plt.close()

# ========== 绘图2: 能量vs温度 ==========
print(f"\n[5] 绘制能量-温度曲线")
print("-"*80)

fig, axes = plt.subplots(1, 2, figsize=(16, 6))
fig.suptitle('Cv系列(Sn8Pt6O4) 能量分析 - 5次独立模拟', fontsize=16, fontweight='bold')

# 子图1: 每原子能量
ax1 = axes[0]
for i, struct in enumerate(sorted(df_energy['结构'].unique())):
    df_struct = df_energy[df_energy['结构'] == struct].sort_values('温度')
    ax1.plot(df_struct['温度'], df_struct['每原子能量'], 
            marker=markers[i], color=colors[i], label=struct,
            linewidth=2, markersize=8, alpha=0.7)

df_mean = df_energy.groupby('温度')['每原子能量'].mean().reset_index()
ax1.plot(df_mean['温度'], df_mean['每原子能量'],
        'k--', linewidth=3, label='平均值', zorder=10)

ax1.set_xlabel('温度 [K]', fontsize=12, fontweight='bold')
ax1.set_ylabel('每原子能量 [eV/atom]', fontsize=12, fontweight='bold')
ax1.set_title('每原子能量 vs 温度', fontsize=13, fontweight='bold')
ax1.legend(loc='best', fontsize=10)
ax1.grid(True, alpha=0.3)

# 子图2: 相对能量
ax2 = axes[1]
for i, struct in enumerate(sorted(df_energy['结构'].unique())):
    df_struct = df_energy[df_energy['结构'] == struct].sort_values('温度')
    ax2.plot(df_struct['温度'], df_struct['相对能量'], 
            marker=markers[i], color=colors[i], label=struct,
            linewidth=2, markersize=8, alpha=0.7)

df_mean = df_energy.groupby('温度')['相对能量'].mean().reset_index()
ax2.plot(df_mean['温度'], df_mean['相对能量'],
        'k--', linewidth=3, label='平均值', zorder=10)

ax2.set_xlabel('温度 [K]', fontsize=12, fontweight='bold')
ax2.set_ylabel('相对能量 [eV]', fontsize=12, fontweight='bold')
ax2.set_title('相对能量 vs 温度', fontsize=13, fontweight='bold')
ax2.legend(loc='best', fontsize=10)
ax2.grid(True, alpha=0.3)

plt.tight_layout()
fig_path = output_dir / 'Cv_series_energy_analysis.png'
plt.savefig(fig_path, dpi=300, bbox_inches='tight')
print(f"✓ 已保存图片: {fig_path}")
plt.close()

# ========== 输出详细数据表 ==========
print(f"\n[6] 生成详细数据表")
print("-"*80)

# 热容数据表(宽格式)
cv_pivot = df_cv.pivot(index='温度', columns='结构', values='Cv_total_meV_K')
cv_pivot['平均值'] = df_cv.groupby('温度')['Cv_total_meV_K'].mean()
cv_pivot['标准差'] = df_cv.groupby('温度')['Cv_total_meV_K'].std()
cv_pivot['RSD%'] = (cv_pivot['标准差'] / cv_pivot['平均值'] * 100).round(2)

csv_path = output_dir / 'Cv_series_heat_capacity_table.csv'
cv_pivot.to_csv(csv_path)
print(f"✓ 已保存热容数据表: {csv_path}")

print("\n热容数据汇总 [meV/K]:")
print(cv_pivot.round(3))

# ========== 团簇热容分析(扣除载体) ==========
print(f"\n[7] 团簇热容分析(扣除载体38.3 meV/K)")
print("-"*80)

SUPPORT_CV = 38.3
df_cv['Cv_cluster'] = df_cv['Cv_total_meV_K'] - SUPPORT_CV
df_cv['Cv_cluster_per_atom'] = df_cv['Cv_cluster'] / df_cv['N_cluster']

print("\n团簇热容统计:")
cluster_stats = df_cv.groupby('温度').agg({
    'Cv_cluster': ['mean', 'std', 'min', 'max'],
    'Cv_cluster_per_atom': ['mean', 'std', 'min', 'max']
}).round(3)
print(cluster_stats)

# 绘制团簇热容
fig, axes = plt.subplots(1, 2, figsize=(16, 6))
fig.suptitle(f'Cv系列(Sn8Pt6O4) 团簇热容分析 (扣除载体{SUPPORT_CV} meV/K)', 
            fontsize=16, fontweight='bold')

# 子图1: 团簇总热容
ax1 = axes[0]
for i, struct in enumerate(sorted(df_cv['结构'].unique())):
    df_struct = df_cv[df_cv['结构'] == struct].sort_values('温度')
    # 只绘制正值
    df_positive = df_struct[df_struct['Cv_cluster'] > 0]
    if len(df_positive) > 0:
        ax1.plot(df_positive['温度'], df_positive['Cv_cluster'], 
                marker=markers[i], color=colors[i], label=struct,
                linewidth=2, markersize=8, alpha=0.7)

df_mean = df_cv[df_cv['Cv_cluster'] > 0].groupby('温度')['Cv_cluster'].mean().reset_index()
if len(df_mean) > 0:
    ax1.plot(df_mean['温度'], df_mean['Cv_cluster'],
            'k--', linewidth=3, label='平均值', zorder=10)

ax1.axhline(0, color='red', linestyle=':', linewidth=2, alpha=0.5)
ax1.set_xlabel('温度 [K]', fontsize=12, fontweight='bold')
ax1.set_ylabel('Cv_cluster [meV/K]', fontsize=12, fontweight='bold')
ax1.set_title('团簇热容(已扣除载体)', fontsize=13, fontweight='bold')
ax1.legend(loc='best', fontsize=10)
ax1.grid(True, alpha=0.3)

# 子图2: 团簇每原子热容
ax2 = axes[1]
for i, struct in enumerate(sorted(df_cv['结构'].unique())):
    df_struct = df_cv[df_cv['结构'] == struct].sort_values('温度')
    df_positive = df_struct[df_struct['Cv_cluster_per_atom'] > 0]
    if len(df_positive) > 0:
        ax2.plot(df_positive['温度'], df_positive['Cv_cluster_per_atom'], 
                marker=markers[i], color=colors[i], label=struct,
                linewidth=2, markersize=8, alpha=0.7)

df_mean = df_cv[df_cv['Cv_cluster_per_atom'] > 0].groupby('温度')['Cv_cluster_per_atom'].mean().reset_index()
if len(df_mean) > 0:
    ax2.plot(df_mean['温度'], df_mean['Cv_cluster_per_atom'],
            'k--', linewidth=3, label='平均值', zorder=10)

ax2.axhline(classical_limit, color='red', linestyle=':', linewidth=2, 
           label=f'经典极限 ({classical_limit:.3f})')
ax2.axhline(0, color='gray', linestyle=':', linewidth=1, alpha=0.5)

ax2.set_xlabel('温度 [K]', fontsize=12, fontweight='bold')
ax2.set_ylabel('Cv_cluster/N [meV/K/atom]', fontsize=12, fontweight='bold')
ax2.set_title('团簇每原子热容', fontsize=13, fontweight='bold')
ax2.legend(loc='best', fontsize=10)
ax2.grid(True, alpha=0.3)

plt.tight_layout()
fig_path = output_dir / 'Cv_series_cluster_heat_capacity.png'
plt.savefig(fig_path, dpi=300, bbox_inches='tight')
print(f"✓ 已保存图片: {fig_path}")
plt.close()

# 保存团簇热容数据
cluster_cv_pivot = df_cv.pivot(index='温度', columns='结构', values='Cv_cluster')
cluster_cv_pivot['平均值'] = df_cv.groupby('温度')['Cv_cluster'].mean()
cluster_cv_pivot['标准差'] = df_cv.groupby('温度')['Cv_cluster'].std()

csv_path = output_dir / 'Cv_series_cluster_heat_capacity_table.csv'
cluster_cv_pivot.to_csv(csv_path)
print(f"✓ 已保存团簇热容数据表: {csv_path}")

print("\n" + "="*80)
print("分析完成!")
print(f"结果保存在: {output_dir}")
print("="*80)
