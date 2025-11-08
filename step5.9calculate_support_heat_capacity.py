#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
计算Al2O3载体的热容并拟合分析
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
from scipy import stats
from pathlib import Path

# 设置中文字体
matplotlib.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei']
matplotlib.rcParams['axes.unicode_minus'] = False

print("="*80)
print("Al2O3载体热容计算与拟合分析")
print("="*80)

# 载体能量数据 (来自载体计算)
# 数据来源: /home/scms/jychen/tools/cp2k/md/msd/nnmD/dp-md/20251009/ener/run1
print("\n[数据来源] CP2K计算的Al2O3载体能量")
print("路径: /home/scms/jychen/tools/cp2k/md/msd/nnmD/dp-md/20251009/ener/run1")
print("-"*80)

# 从你提供的数据提取温度和平均总能
data_raw = """
200,0,-68659.91151634,0.22964220
300,1,-68656.10964780,0.34494389
400,2,-68652.30600862,0.46054243
500,3,-68648.49661528,0.57602641
600,4,-68644.68103945,0.69292374
700,5,-68640.86380918,0.80979003
800,6,-68637.04735694,0.92379395
900,7,-68633.22271738,1.04191388
1000,8,-68629.38836185,1.16122292
1100,9,-68625.55944410,1.28091887
"""

# 解析数据
data_lines = [line.strip() for line in data_raw.strip().split('\n') if line.strip()]
temperatures = []
energies = []
std_devs = []

for line in data_lines:
    parts = line.split(',')
    temp = float(parts[0])
    energy = float(parts[2])
    std = float(parts[3])
    temperatures.append(temp)
    energies.append(energy)
    std_devs.append(std)

temperatures = np.array(temperatures)
energies = np.array(energies)
std_devs = np.array(std_devs)

# 数据展示
print("\n[1] 输入数据")
print("-"*80)
df_input = pd.DataFrame({
    '温度 [K]': temperatures,
    '能量 [eV]': energies,
    '标准差 [eV]': std_devs
})
print(df_input.to_string(index=False))

# 计算热容 Cv = dE/dT
print("\n[2] 热容计算 (Cv = dE/dT)")
print("-"*80)

# 使用中心差分法计算导数
heat_capacities_eV = np.gradient(energies, temperatures)  # eV/K
heat_capacities_meV = heat_capacities_eV * 1000  # meV/K

# 计算平均热容
Cv_avg_eV = np.mean(heat_capacities_eV)
Cv_std_eV = np.std(heat_capacities_eV)
Cv_avg_meV = Cv_avg_eV * 1000
Cv_std_meV = Cv_std_eV * 1000

print(f"\n各温度点的热容:")
df_cv = pd.DataFrame({
    '温度 [K]': temperatures,
    '能量 [eV]': energies,
    '能量标准差 [eV]': std_devs,
    'Cv [eV/K]': heat_capacities_eV,
    'Cv [meV/K]': heat_capacities_meV
})
print(df_cv.to_string(index=False))

print(f"\n平均热容: Cv_avg = {Cv_avg_meV:.3f} ± {Cv_std_meV:.3f} meV/K")
print(f"           Cv_avg = {Cv_avg_eV:.6f} ± {Cv_std_eV:.6f} eV/K")

# 线性拟合 E(T) = a + b*T
print("\n[3] 线性拟合 E(T) = a + b*T")
print("-"*80)

slope, intercept, r_value, p_value, std_err = stats.linregress(temperatures, energies)
r_squared = r_value ** 2

print(f"\n拟合结果:")
print(f"  斜率 (b) = {slope:.6f} eV/K = {slope*1000:.3f} meV/K")
print(f"  截距 (a) = {intercept:.3f} eV")
print(f"  R² = {r_squared:.6f}")
print(f"  p-value = {p_value:.3e}")
print(f"  标准误差 = {std_err:.6f} eV/K")

print(f"\n热容 Cv = dE/dT = {slope*1000:.3f} meV/K (来自线性拟合)")

# 计算拟合值和残差
fitted_energies = intercept + slope * temperatures
residuals = energies - fitted_energies
rmse = np.sqrt(np.mean(residuals**2))

print(f"\nRMSE = {rmse:.6f} eV")
print(f"最大残差 = {np.max(np.abs(residuals)):.6f} eV")

# 如果R²不够高,尝试多项式拟合
if r_squared < 0.99:
    print("\n[4] 多项式拟合尝试 (R²<0.99)")
    print("-"*80)
    
    # 二次拟合 E(T) = a + b*T + c*T²
    coeffs_2 = np.polyfit(temperatures, energies, 2)
    poly_2 = np.poly1d(coeffs_2)
    fitted_2 = poly_2(temperatures)
    r_squared_2 = 1 - np.sum((energies - fitted_2)**2) / np.sum((energies - np.mean(energies))**2)
    
    # Cv = dE/dT = b + 2*c*T (温度依赖)
    Cv_func_2 = lambda T: (coeffs_2[1] + 2*coeffs_2[0]*T) * 1000  # meV/K
    
    print(f"\n二次拟合: E(T) = {coeffs_2[2]:.3f} + {coeffs_2[1]:.6f}*T + {coeffs_2[0]:.9f}*T²")
    print(f"  R² = {r_squared_2:.6f}")
    print(f"  Cv(T) = {coeffs_2[1]*1000:.3f} + {2*coeffs_2[0]*1000:.6f}*T meV/K")
    
    # 计算不同温度下的Cv
    print(f"\n温度依赖的热容:")
    for T in [200, 400, 600, 800, 1000]:
        Cv_T = Cv_func_2(T)
        print(f"  Cv({T}K) = {Cv_T:.3f} meV/K")
    
    # 三次拟合
    coeffs_3 = np.polyfit(temperatures, energies, 3)
    poly_3 = np.poly1d(coeffs_3)
    fitted_3 = poly_3(temperatures)
    r_squared_3 = 1 - np.sum((energies - fitted_3)**2) / np.sum((energies - np.mean(energies))**2)
    
    print(f"\n三次拟合: R² = {r_squared_3:.6f}")

# 绘图
print("\n[5] 生成可视化图表")
print("-"*80)

fig, axes = plt.subplots(2, 2, figsize=(16, 12))
fig.suptitle('Al₂O₃载体热容分析', fontsize=16, fontweight='bold')

# 子图1: 能量 vs 温度 + 线性拟合
ax1 = axes[0, 0]
ax1.errorbar(temperatures, energies, yerr=std_devs, fmt='o', markersize=10, 
            color='darkblue', ecolor='lightblue', capsize=5, capthick=2,
            label='实际数据 (含误差棒)', zorder=3)
ax1.plot(temperatures, fitted_energies, '--', linewidth=2, color='red',
        label=f'线性拟合 (R²={r_squared:.6f})', zorder=2)

# 如果有多项式拟合
if r_squared < 0.99:
    T_fine = np.linspace(temperatures.min(), temperatures.max(), 100)
    ax1.plot(T_fine, poly_2(T_fine), ':', linewidth=2, color='green',
            label=f'二次拟合 (R²={r_squared_2:.4f})', zorder=1)

ax1.set_xlabel('温度 [K]', fontsize=12, fontweight='bold')
ax1.set_ylabel('能量 [eV]', fontsize=12, fontweight='bold')
ax1.set_title('能量 vs 温度', fontsize=13, fontweight='bold')
ax1.legend(loc='best', fontsize=10)
ax1.grid(True, alpha=0.3)

# 子图2: 残差图
ax2 = axes[0, 1]
ax2.plot(temperatures, residuals * 1000, 'o-', markersize=8, color='purple',
        linewidth=2)
ax2.axhline(0, color='red', linestyle='--', linewidth=2, alpha=0.7)
ax2.fill_between(temperatures, 
                residuals * 1000 - rmse*1000, 
                residuals * 1000 + rmse*1000,
                alpha=0.2, color='purple')

ax2.set_xlabel('温度 [K]', fontsize=12, fontweight='bold')
ax2.set_ylabel('残差 [meV]', fontsize=12, fontweight='bold')
ax2.set_title(f'拟合残差 (RMSE={rmse*1000:.3f} meV)', fontsize=13, fontweight='bold')
ax2.grid(True, alpha=0.3)

# 子图3: 热容 vs 温度
ax3 = axes[1, 0]
ax3.plot(temperatures, heat_capacities_meV, 'o-', markersize=10, 
        color='darkred', linewidth=2, label='数值微分 Cv')
ax3.axhline(Cv_avg_meV, color='blue', linestyle='--', linewidth=2,
           label=f'平均值 = {Cv_avg_meV:.3f} meV/K')
ax3.axhline(slope*1000, color='green', linestyle=':', linewidth=2,
           label=f'线性拟合斜率 = {slope*1000:.3f} meV/K')

# 如果有温度依赖
if r_squared < 0.99:
    T_fine = np.linspace(temperatures.min(), temperatures.max(), 100)
    Cv_fine = Cv_func_2(T_fine)
    ax3.plot(T_fine, Cv_fine, '--', linewidth=2, color='orange',
            label='Cv(T) 二次拟合')

ax3.set_xlabel('温度 [K]', fontsize=12, fontweight='bold')
ax3.set_ylabel('热容 Cv [meV/K]', fontsize=12, fontweight='bold')
ax3.set_title('热容 vs 温度', fontsize=13, fontweight='bold')
ax3.legend(loc='best', fontsize=10)
ax3.grid(True, alpha=0.3)

# 子图4: Q-Q图 (检验正态性)
ax4 = axes[1, 1]
from scipy import stats as sp_stats
sp_stats.probplot(residuals, dist="norm", plot=ax4)
ax4.set_title('Q-Q图 (残差正态性检验)', fontsize=13, fontweight='bold')
ax4.grid(True, alpha=0.3)

plt.tight_layout()

# 保存
output_dir = Path('results/support_heat_capacity_analysis')
output_dir.mkdir(parents=True, exist_ok=True)

fig_path = output_dir / 'support_heat_capacity_fitting.png'
plt.savefig(fig_path, dpi=300, bbox_inches='tight')
print(f"✓ 已保存图片: {fig_path}")
plt.close()

# 保存结果
csv_path = output_dir / 'support_heat_capacity_results.csv'
df_results = pd.DataFrame({
    '温度 [K]': temperatures,
    '能量 [eV]': energies,
    '能量标准差 [eV]': std_devs,
    '拟合能量 [eV]': fitted_energies,
    '残差 [meV]': residuals * 1000,
    'Cv数值微分 [meV/K]': heat_capacities_meV
})
df_results.to_csv(csv_path, index=False, encoding='utf-8-sig')
print(f"✓ 已保存数据: {csv_path}")

# 保存摘要
summary_path = output_dir / 'support_heat_capacity_summary.txt'
with open(summary_path, 'w', encoding='utf-8') as f:
    f.write("="*80 + "\n")
    f.write("Al2O3载体热容分析摘要\n")
    f.write("="*80 + "\n\n")
    
    f.write("[1] 数据信息\n")
    f.write("-"*80 + "\n")
    f.write(f"温度范围: {temperatures.min():.0f} - {temperatures.max():.0f} K\n")
    f.write(f"数据点数: {len(temperatures)}\n")
    f.write(f"能量范围: {energies.min():.3f} - {energies.max():.3f} eV\n\n")
    
    f.write("[2] 线性拟合结果\n")
    f.write("-"*80 + "\n")
    f.write(f"拟合方程: E(T) = {intercept:.3f} + {slope:.6f}*T eV\n")
    f.write(f"R² = {r_squared:.6f}\n")
    f.write(f"p-value = {p_value:.3e}\n")
    f.write(f"RMSE = {rmse*1000:.3f} meV\n\n")
    
    f.write("[3] 热容结果\n")
    f.write("-"*80 + "\n")
    f.write(f"线性拟合热容: Cv = {slope*1000:.3f} meV/K\n")
    f.write(f"数值微分平均: Cv = {Cv_avg_meV:.3f} ± {Cv_std_meV:.3f} meV/K\n\n")
    
    if r_squared < 0.99:
        f.write("[4] 高阶拟合\n")
        f.write("-"*80 + "\n")
        f.write(f"二次拟合 R² = {r_squared_2:.6f}\n")
        f.write(f"三次拟合 R² = {r_squared_3:.6f}\n")
        f.write(f"Cv(T) = {coeffs_2[1]*1000:.3f} + {2*coeffs_2[0]*1000:.6f}*T meV/K\n\n")
    
    f.write("[5] 推荐值\n")
    f.write("-"*80 + "\n")
    if r_squared >= 0.99:
        f.write(f"✓ 线性拟合良好,推荐使用: Cv = {slope*1000:.3f} meV/K\n")
    else:
        f.write(f"! 线性拟合R²={r_squared:.4f}<0.99,建议使用温度依赖的Cv(T)\n")
        f.write(f"  或在特定温度范围内使用平均值: {Cv_avg_meV:.3f} meV/K\n")

print(f"✓ 已保存摘要: {summary_path}")

print("\n" + "="*80)
print("分析完成!")
print("="*80)
print(f"\n推荐载体热容值:")
print(f"  Cv_support = {slope*1000:.3f} meV/K (线性拟合, R²={r_squared:.4f})")
print(f"  Cv_support = {Cv_avg_meV:.3f} ± {Cv_std_meV:.3f} meV/K (数值平均)")
print(f"\n结果文件保存在: {output_dir}")
print("="*80)
