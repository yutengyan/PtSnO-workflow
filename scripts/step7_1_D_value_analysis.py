#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
Step7.1: 扩散系数(D值)分析
================================================================================

作者: GitHub Copilot
日期: 2025-10-29
版本: v2.0

功能概述
========
本脚本用于分析LAMMPS模拟的扩散系数(Diffusion Coefficient, D值)数据,
判断团簇的流动性状态,支持多体系分析和异常筛选。

扩散系数判据
============
- D < 1e-8 cm²/s:  低流动性 (Low Mobility)
- 1e-8 < D < 1e-7: 中等流动性 (Medium Mobility)  
- D > 1e-7:        高流动性 (High Mobility)

D值来源说明
============
D值来自 ensemble_analysis_results.csv 文件中的 D_ensemble 字段。

计算方法:
1. 对同一体系、同一温度、同一元素的多个独立模拟run进行**集合平均**
2. 从集合平均后的MSD曲线重新进行线性拟合: MSD = 6Dt + intercept
3. 提取斜率得到 D_ensemble (扩散系数,单位: cm²/s)

与原始D值的区别:
- 原始D值: 每个run单独拟合得到的D值
- D_ensemble: 多个run的MSD先平均,再拟合得到的D值
- D_ensemble更稳定,能更好地反映系综行为

数据字段:
- composition: 结构名称 (例如: Pt8Sn6, g-1-O1Sn4Pt3)
- temperature: 温度字符串 (例如: 300K, 1000K)
- temp_value: 温度数值 (例如: 300, 1000)
- element: 元素类型 (Pt, Sn, PtSn)
- D_ensemble: 集合平均扩散系数 (cm²/s)
- R2_ensemble: 拟合优度
- n_runs: 有效run数量
- n_runs_original: 原始run数量
- n_runs_filtered: 被异常筛选过滤的run数量

核心功能
========
1. **智能数据加载**
   - 加载ensemble_analysis_results.csv (包含D_ensemble)
   - 支持结构排除过滤 (--exclude参数)
   - 自动识别系列分类 (O0, O1-3, O4+)

2. **扩散系数分析**
   - 计算各体系在不同温度下的D值
   - 识别流动性转变温度
   - 统计分析

3. **可视化** (类似林德曼指数分析)
   - D值热力图 (按系列分组)
   - D值-温度曲线 (所有结构汇总在一张图)
   - 流动性转变温度分析

输出文件
========
results/D_value_analysis/
├── D_value_heatmap_*.png           # D值热力图 (按系列)
├── D_vs_temperature_*.png          # D值-温度曲线 (按系列)
├── transition_temperature_*.png    # 流动性转变温度分析
├── D_value_analysis_report.txt     # 分析报告
├── D_value_statistics.csv          # 统计数据
└── mobility_transition_temperatures.csv  # 转变温度数据
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import argparse
import warnings
import re
warnings.filterwarnings('ignore')

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# ===== 全局配置 =====
BASE_DIR = Path(__file__).parent.parent
RESULTS_DIR = BASE_DIR / 'results'
ENSEMBLE_CSV = RESULTS_DIR / 'ensemble_analysis_results.csv'
OUTPUT_DIR = RESULTS_DIR / 'D_value_analysis'
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

# D值分类阈值 (cm²/s)
D_LOW_THRESHOLD = 1e-8
D_HIGH_THRESHOLD = 1e-7


def load_D_value_data(exclude_structures=None):
    """
    加载D值数据
    
    Parameters:
    -----------
    exclude_structures : list, optional
        要排除的结构列表
    
    Returns:
    --------
    df : DataFrame
        包含D值的数据框
    """
    print("\n" + "="*80)
    print("加载D值数据")
    print("="*80)
    
    if not ENSEMBLE_CSV.exists():
        print(f"[X] 错误: 未找到数据文件 {ENSEMBLE_CSV}")
        return None
    
    # 加载数据
    df = pd.read_csv(ENSEMBLE_CSV)
    print(f"  [√] 加载了 {len(df)} 条记录")
    
    # 重命名列
    df = df.rename(columns={
        'composition': '结构',
        'temperature': '温度',
        'temp_value': '温度值',
        'element': '元素',
        'D_ensemble': 'D值'
    })
    
    # 提取原子数信息
    def extract_atom_counts(composition):
        """从composition提取Pt, Sn, O原子数"""
        import re
        pt = sn = o = 0
        
        # 匹配 Pt 数量
        pt_match = re.search(r'Pt(\d+)', composition, re.IGNORECASE)
        if pt_match:
            pt = int(pt_match.group(1))
        
        # 匹配 Sn 数量
        sn_match = re.search(r'Sn(\d+)', composition, re.IGNORECASE)
        if sn_match:
            sn = int(sn_match.group(1))
        
        # 匹配 O 数量
        o_match = re.search(r'O(\d+)', composition, re.IGNORECASE)
        if o_match:
            o = int(o_match.group(1))
        
        return pt, sn, o
    
    # 应用提取函数
    atom_counts = df['结构'].apply(extract_atom_counts)
    df['Pt原子数'] = atom_counts.apply(lambda x: x[0])
    df['Sn原子数'] = atom_counts.apply(lambda x: x[1])
    df['O原子数'] = atom_counts.apply(lambda x: x[2])
    
    # 确定系列 (使用与林德曼分析相同的分类方法)
    def classify_series(structure_name):
        """
        分类体系所属系列 (与step7_lindemann_analysis.py保持一致)
        
        Returns:
        --------
        series : str
            'O1', 'O2', 'O3', 'O4', 'Pt6SnX', 'Pt8SnX', 'Pt(8-x)SnX', 'Cv', 'Other'
        """
        if not structure_name:
            return 'Other'
        
        structure_lower = structure_name.lower()
        
        # 优先检查含氧系列
        if re.search(r'[Oo]1(?:[^0-9]|$)', structure_lower):
            return 'O1'
        if re.search(r'[Oo]2(?:[^0-9]|$)', structure_lower):
            return 'O2'
        if re.search(r'[Oo]3(?:[^0-9]|$)', structure_lower):
            return 'O3'
        if re.search(r'[Oo]4(?:[^0-9]|$)', structure_lower):
            return 'O4'
        
        # Cv系列
        if structure_lower.startswith('cv-'):
            return 'Cv'
        
        # 提取Pt和Sn原子数
        pt_match = re.search(r'pt(\d+)', structure_lower)
        sn_match = re.search(r'sn(\d+)', structure_lower)
        
        if pt_match and sn_match:
            n_pt = int(pt_match.group(1))
            n_sn = int(sn_match.group(1))
            
            # Pt8SnX系列 (固定8个Pt原子)
            if n_pt == 8:
                return 'Pt8SnX'
            
            # Pt6SnX系列 (固定6个Pt原子)
            if n_pt == 6:
                return 'Pt6SnX'
            
            # Pt(8-x)SnX系列: 总原子数=8
            if n_pt + n_sn == 8:
                return 'Pt(8-x)SnX'
        
        # 纯Pt系列 (pt8, pt6等)
        if pt_match and not sn_match:
            n_pt = int(pt_match.group(1))
            if n_pt == 8:
                return 'Pt8SnX'  # 归入Pt8SnX系列
            if n_pt == 6:
                return 'Pt6SnX'  # 归入Pt6SnX系列
        
        return 'Other'
    
    df['系列'] = df['结构'].apply(classify_series)
    
    print(f"  系列分布:")
    for series, count in df.groupby('系列').size().items():
        print(f"    {series}: {count} 条记录")
    
    # 应用结构排除过滤
    if exclude_structures:
        original_count = len(df)
        df = df[~df['结构'].isin(exclude_structures)]
        filtered_count = original_count - len(df)
        
        print(f"\n  应用结构过滤:")
        print(f"    排除的结构: {', '.join(exclude_structures)}")
        print(f"    过滤前: {original_count} 条")
        print(f"    过滤后: {len(df)} 条")
        print(f"    已删除: {filtered_count} 条")
    
    return df


def classify_mobility(D_value):
    """
    根据D值分类流动性
    
    Parameters:
    -----------
    D_value : float
        扩散系数 (cm²/s)
    
    Returns:
    --------
    category : str
        流动性分类
    """
    if pd.isna(D_value):
        return '未知'
    elif D_value < D_LOW_THRESHOLD:
        return '低流动性'
    elif D_value < D_HIGH_THRESHOLD:
        return '中等流动性'
    else:
        return '高流动性'


def analyze_D_values(df):
    """
    分析D值统计
    
    Parameters:
    -----------
    df : DataFrame
        D值数据
    
    Returns:
    --------
    df_stats : DataFrame
        统计结果
    """
    print("\n" + "="*80)
    print("D值统计分析")
    print("="*80)
    
    # 按结构和温度分组
    stats_list = []
    
    for struct in sorted(df['结构'].unique()):
        df_struct = df[df['结构'] == struct]
        
        # 获取元素信息
        row_sample = df_struct.iloc[0]
        series = row_sample['系列']
        pt = row_sample['Pt原子数']
        sn = row_sample['Sn原子数']
        o = row_sample['O原子数']
        
        # 按温度统计
        for temp in sorted(df_struct['温度值'].unique()):
            df_temp = df_struct[df_struct['温度值'] == temp]
            
            # 计算各元素的D值
            D_Pt = df_temp[df_temp['元素'] == 'Pt']['D值'].mean()
            D_Sn = df_temp[df_temp['元素'] == 'Sn']['D值'].mean()
            D_PtSn = df_temp[df_temp['元素'] == 'PtSn']['D值'].mean()
            
            # 平均D值
            D_mean = df_temp['D值'].mean()
            
            stats_list.append({
                '结构': struct,
                '系列': series,
                'Pt原子数': pt,
                'Sn原子数': sn,
                'O原子数': o,
                '温度': temp,
                'D_Pt': D_Pt,
                'D_Sn': D_Sn,
                'D_PtSn': D_PtSn,
                'D_mean': D_mean,
                '流动性': classify_mobility(D_mean)
            })
    
    df_stats = pd.DataFrame(stats_list)
    
    # 保存统计结果
    stats_file = OUTPUT_DIR / 'D_value_statistics.csv'
    df_stats.to_csv(stats_file, index=False, encoding='utf-8-sig')
    print(f"  [√] 统计结果已保存: {stats_file.name}")
    
    return df_stats


def plot_D_value_heatmap(df_stats):
    """
    绘制D值热力图
    
    Parameters:
    -----------
    df_stats : DataFrame
        统计数据
    """
    print("\n[*] 绘制D值热力图...")
    
    # 按系列分组
    series_list = sorted(df_stats['系列'].unique())
    
    for series in series_list:
        df_series = df_stats[df_stats['系列'] == series].copy()
        
        if len(df_series) == 0:
            continue
        
        # 提取结构信息并排序
        struct_info = []
        for struct in df_series['结构'].unique():
            row = df_series[df_series['结构'] == struct].iloc[0]
            pt = int(row['Pt原子数'])
            sn = int(row['Sn原子数'])
            o = int(row['O原子数'])
            total = pt + sn + o
            pt_ratio = pt / (pt + sn) if (pt + sn) > 0 else 0
            struct_info.append((struct, pt, sn, o, total, pt_ratio))
        
        # 排序: 按总原子数, 然后按Pt比例
        struct_info_sorted = sorted(struct_info, key=lambda x: (x[4], x[5]))
        sorted_structures = [s[0] for s in struct_info_sorted]
        
        # 准备数据透视表 (使用D_mean)
        pivot_data = df_series.groupby(['温度', '结构'])['D_mean'].mean().reset_index()
        heatmap_data = pivot_data.pivot(index='温度', columns='结构', values='D_mean')
        
        # 重新排列列顺序
        heatmap_data = heatmap_data[sorted_structures]
        
        # 转换为log10标度
        heatmap_data_log = np.log10(heatmap_data.values + 1e-15)  # 避免log(0)
        
        # 创建图形
        fig, ax = plt.subplots(figsize=(14, 8))
        
        # 绘制热力图 (使用log标度)
        im = ax.imshow(heatmap_data_log, aspect='auto', cmap=plt.cm.viridis,
                      interpolation='bilinear', origin='lower',
                      vmin=-15, vmax=-6)
        
        # 设置坐标轴
        ax.set_xticks(np.arange(len(heatmap_data.columns)))
        ax.set_yticks(np.arange(len(heatmap_data.index)))
        
        # X轴标签
        ax.set_xticklabels(heatmap_data.columns, rotation=45, ha='right', fontsize=10)
        
        # Y轴标签
        ax.set_yticklabels([f'{int(t)}K' for t in heatmap_data.index], fontsize=10)
        
        ax.set_xlabel('结构 (Structure)', fontsize=12, fontweight='bold')
        ax.set_ylabel('温度 (Temperature)', fontsize=12, fontweight='bold')
        ax.set_title(f'扩散系数热力图 - {series}系列\nDiffusion Coefficient Heatmap (log₁₀ scale, cm²/s)',
                    fontsize=14, fontweight='bold', pad=15)
        
        # 添加colorbar
        cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label('log₁₀(D) [cm²/s]', fontsize=11, fontweight='bold')
        
        # 添加参考线 (流动性阈值)
        # 计算对应的log10值
        low_line = np.log10(D_LOW_THRESHOLD)
        high_line = np.log10(D_HIGH_THRESHOLD)
        
        # 在colorbar上添加参考线
        cbar.ax.axhline(y=low_line, color='cyan', linestyle='--', linewidth=2, 
                       label=f'低流动性: {D_LOW_THRESHOLD:.0e}')
        cbar.ax.axhline(y=high_line, color='red', linestyle='--', linewidth=2,
                       label=f'高流动性: {D_HIGH_THRESHOLD:.0e}')
        
        plt.tight_layout()
        
        # 保存图片
        filename = OUTPUT_DIR / f'D_value_heatmap_{series}.png'
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"  [√] 已保存: {filename.name}")
        plt.close()


def plot_D_vs_temperature(df_stats):
    """
    绘制D值-温度曲线 (类似林德曼分析,所有结构在一张图上)
    
    Parameters:
    -----------
    df_stats : DataFrame
        统计数据
    """
    print("\n[*] 绘制D值-温度曲线...")
    
    # 按系列分组
    series_list = sorted(df_stats['系列'].unique())
    
    for series in series_list:
        df_series = df_stats[df_stats['系列'] == series].copy()
        
        if len(df_series) == 0:
            continue
        
        # 提取唯一结构
        structures = sorted(df_series['结构'].unique())
        
        # 创建单图,所有结构在一起
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # 为每个结构绘制曲线
        import matplotlib.cm as cm
        colors = cm.tab20(np.linspace(0, 1, len(structures)))
        
        for idx, struct in enumerate(structures):
            df_struct = df_series[df_series['结构'] == struct].sort_values('温度')
            temps = df_struct['温度'].values
            D_mean = df_struct['D_mean'].values
            
            # 只绘制平均D值
            ax.semilogy(temps, D_mean, 'o-', label=struct,
                       color=colors[idx], linewidth=2, markersize=6, alpha=0.7)
        
        # 添加流动性阈值线
        ax.axhline(y=D_LOW_THRESHOLD, color='cyan', linestyle='--', 
                  linewidth=2, alpha=0.7, label=f'低流动性阈值: {D_LOW_THRESHOLD:.0e}')
        ax.axhline(y=D_HIGH_THRESHOLD, color='orange', linestyle='--',
                  linewidth=2, alpha=0.7, label=f'高流动性阈值: {D_HIGH_THRESHOLD:.0e}')
        
        ax.set_xlabel('温度 (K)', fontsize=13, fontweight='bold')
        ax.set_ylabel('D (cm²/s)', fontsize=13, fontweight='bold')
        ax.set_title(f'扩散系数-温度关系 - {series}系列\nDiffusion Coefficient vs Temperature',
                    fontsize=15, fontweight='bold')
        ax.legend(fontsize=9, loc='best', ncol=2)
        ax.grid(True, alpha=0.3, which='both')
        
        plt.tight_layout()
        
        # 保存图片
        filename = OUTPUT_DIR / f'D_vs_temperature_{series}.png'
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"  [√] 已保存: {filename.name}")
        plt.close()


def plot_individual_structure_analysis(df_stats):
    """
    为每个结构绘制独立的详细分析图
    
    Parameters:
    -----------
    df_stats : DataFrame
        统计数据
    """
    print("\n[*] 绘制各结构的详细分析图...")
    
    structures = sorted(df_stats['结构'].unique())
    
    for struct in structures:
        df_struct = df_stats[df_stats['结构'] == struct].copy()
        df_struct = df_struct.sort_values('温度')
        
        temps = df_struct['温度'].values
        D_Pt = df_struct['D_Pt'].values
        D_Sn = df_struct['D_Sn'].values
        D_PtSn = df_struct['D_PtSn'].values
        D_mean = df_struct['D_mean'].values
        
        # 创建4子图布局
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        fig.suptitle(f'{struct} - 扩散系数详细分析', fontsize=15, fontweight='bold')
        
        # 子图1: D值随温度变化 (log scale)
        ax1 = axes[0, 0]
        if not np.all(np.isnan(D_Pt)):
            ax1.semilogy(temps, D_Pt, 'o-', label='Pt', color='#E74C3C', 
                        linewidth=2, markersize=8)
        if not np.all(np.isnan(D_Sn)):
            ax1.semilogy(temps, D_Sn, 's-', label='Sn', color='#3498DB',
                        linewidth=2, markersize=8)
        if not np.all(np.isnan(D_PtSn)):
            ax1.semilogy(temps, D_PtSn, '^-', label='PtSn', color='#2ECC71',
                        linewidth=2, markersize=8)
        
        ax1.axhline(D_LOW_THRESHOLD, color='cyan', linestyle='--', 
                   linewidth=1.5, alpha=0.7, label=f'低流动性 ({D_LOW_THRESHOLD:.0e})')
        ax1.axhline(D_HIGH_THRESHOLD, color='orange', linestyle='--',
                   linewidth=1.5, alpha=0.7, label=f'高流动性 ({D_HIGH_THRESHOLD:.0e})')
        
        ax1.set_xlabel('温度 (K)', fontsize=11, fontweight='bold')
        ax1.set_ylabel('D (cm²/s, log scale)', fontsize=11, fontweight='bold')
        ax1.set_title('扩散系数 vs 温度 (对数坐标)', fontsize=12, fontweight='bold')
        ax1.legend(fontsize=9, loc='best')
        ax1.grid(True, alpha=0.3, which='both')
        
        # 子图2: D值随温度变化 (linear scale)
        ax2 = axes[0, 1]
        if not np.all(np.isnan(D_Pt)):
            ax2.plot(temps, D_Pt * 1e8, 'o-', label='Pt', color='#E74C3C',
                    linewidth=2, markersize=8)
        if not np.all(np.isnan(D_Sn)):
            ax2.plot(temps, D_Sn * 1e8, 's-', label='Sn', color='#3498DB',
                    linewidth=2, markersize=8)
        if not np.all(np.isnan(D_PtSn)):
            ax2.plot(temps, D_PtSn * 1e8, '^-', label='PtSn', color='#2ECC71',
                    linewidth=2, markersize=8)
        
        ax2.axhline(D_LOW_THRESHOLD * 1e8, color='cyan', linestyle='--',
                   linewidth=1.5, alpha=0.7)
        ax2.axhline(D_HIGH_THRESHOLD * 1e8, color='orange', linestyle='--',
                   linewidth=1.5, alpha=0.7)
        
        ax2.set_xlabel('温度 (K)', fontsize=11, fontweight='bold')
        ax2.set_ylabel('D (×10⁻⁸ cm²/s)', fontsize=11, fontweight='bold')
        ax2.set_title('扩散系数 vs 温度 (线性坐标)', fontsize=12, fontweight='bold')
        ax2.legend(fontsize=9, loc='best')
        ax2.grid(True, alpha=0.3)
        
        # 子图3: Pt/Sn D值比值
        ax3 = axes[1, 0]
        if not np.all(np.isnan(D_Pt)) and not np.all(np.isnan(D_Sn)):
            ratio = D_Pt / (D_Sn + 1e-15)  # 避免除零
            ax3.plot(temps, ratio, 'D-', color='purple', linewidth=2, markersize=8)
            ax3.axhline(1.0, color='gray', linestyle='--', linewidth=1.5, alpha=0.7,
                       label='D_Pt = D_Sn')
        
        ax3.set_xlabel('温度 (K)', fontsize=11, fontweight='bold')
        ax3.set_ylabel('D_Pt / D_Sn', fontsize=11, fontweight='bold')
        ax3.set_title('Pt/Sn扩散系数比值', fontsize=12, fontweight='bold')
        ax3.legend(fontsize=9, loc='best')
        ax3.grid(True, alpha=0.3)
        
        # 子图4: 流动性状态分布
        ax4 = axes[1, 1]
        mobility_counts = df_struct['流动性'].value_counts()
        colors_mob = {'低流动性': '#3498DB', '中等流动性': '#F39C12', '高流动性': '#E74C3C'}
        
        bars = ax4.bar(range(len(mobility_counts)), mobility_counts.values,
                      color=[colors_mob.get(m, 'gray') for m in mobility_counts.index],
                      edgecolor='black', linewidth=2, alpha=0.7)
        ax4.set_xticks(range(len(mobility_counts)))
        ax4.set_xticklabels(mobility_counts.index, rotation=15, ha='right')
        ax4.set_ylabel('温度点数', fontsize=11, fontweight='bold')
        ax4.set_title('流动性状态分布', fontsize=12, fontweight='bold')
        ax4.grid(True, alpha=0.3, axis='y')
        
        # 添加数值标签
        for bar in bars:
            height = bar.get_height()
            ax4.text(bar.get_x() + bar.get_width()/2., height,
                    f'{int(height)}',
                    ha='center', va='bottom', fontsize=10, fontweight='bold')
        
        plt.tight_layout()
        
        # 保存图片
        filename = OUTPUT_DIR / f'D_analysis_{struct}.png'
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"  [√] 已保存: {filename.name}")
        plt.close()


def detect_mobility_transitions(df_stats):
    """
    检测流动性转变温度
    
    Parameters:
    -----------
    df_stats : DataFrame
        统计数据
    
    Returns:
    --------
    df_transitions : DataFrame
        转变温度数据
    """
    print("\n[*] 检测流动性转变温度...")
    
    from scipy.interpolate import interp1d
    
    transitions = []
    
    for struct in sorted(df_stats['结构'].unique()):
        df_struct = df_stats[df_stats['结构'] == struct].copy()
        df_struct = df_struct.sort_values('温度')
        
        if len(df_struct) < 3:
            continue
        
        temps = df_struct['温度'].values
        D_mean = df_struct['D_mean'].values
        
        # 检测从低流动性到高流动性的转变
        # 转变1: D从 < D_LOW 到 > D_LOW (进入中等流动性)
        low_temps = temps[D_mean < D_LOW_THRESHOLD]
        mid_temps = temps[(D_mean >= D_LOW_THRESHOLD) & (D_mean < D_HIGH_THRESHOLD)]
        
        T_transition_low = None
        if len(low_temps) > 0 and len(mid_temps) > 0:
            T_below = low_temps.max()
            T_above = mid_temps.min()
            
            if T_above - T_below <= 100:
                idx_below = np.where(temps == T_below)[0][0]
                idx_above = np.where(temps == T_above)[0][0]
                
                D1 = D_mean[idx_below]
                D2 = D_mean[idx_above]
                
                if D2 > D1:
                    try:
                        f = interp1d([D1, D2], [T_below, T_above], kind='linear')
                        T_transition_low = float(f(D_LOW_THRESHOLD))
                    except:
                        T_transition_low = (T_below + T_above) / 2
        
        # 转变2: D从 < D_HIGH 到 > D_HIGH (进入高流动性)
        mid_high_temps = temps[D_mean < D_HIGH_THRESHOLD]
        high_temps = temps[D_mean >= D_HIGH_THRESHOLD]
        
        T_transition_high = None
        if len(mid_high_temps) > 0 and len(high_temps) > 0:
            T_below = mid_high_temps.max()
            T_above = high_temps.min()
            
            if T_above - T_below <= 100:
                idx_below = np.where(temps == T_below)[0][0]
                idx_above = np.where(temps == T_above)[0][0]
                
                D1 = D_mean[idx_below]
                D2 = D_mean[idx_above]
                
                if D2 > D1:
                    try:
                        f = interp1d([D1, D2], [T_below, T_above], kind='linear')
                        T_transition_high = float(f(D_HIGH_THRESHOLD))
                    except:
                        T_transition_high = (T_below + T_above) / 2
        
        # 记录转变温度
        series = df_struct.iloc[0]['系列']
        transitions.append({
            '结构': struct,
            '系列': series,
            'T_低到中流动性': T_transition_low,
            'T_中到高流动性': T_transition_high
        })
        
        if T_transition_low is not None:
            print(f"  {struct}: 低→中流动性转变 ≈ {T_transition_low:.1f} K")
        if T_transition_high is not None:
            print(f"  {struct}: 中→高流动性转变 ≈ {T_transition_high:.1f} K")
    
    df_transitions = pd.DataFrame(transitions)
    
    # 保存转变温度数据
    trans_file = OUTPUT_DIR / 'mobility_transition_temperatures.csv'
    df_transitions.to_csv(trans_file, index=False, encoding='utf-8-sig')
    print(f"  [√] 转变温度数据已保存: {trans_file.name}")
    
    return df_transitions


def plot_transition_temperature_analysis(df_transitions):
    """
    绘制流动性转变温度分析图
    
    Parameters:
    -----------
    df_transitions : DataFrame
        转变温度数据
    """
    print("\n[*] 绘制流动性转变温度分析图...")
    
    # 过滤有效数据
    df_valid = df_transitions[df_transitions['T_低到中流动性'].notna() | 
                             df_transitions['T_中到高流动性'].notna()].copy()
    
    if len(df_valid) == 0:
        print("  [!] 没有检测到有效的流动性转变温度")
        return
    
    # 提取Sn含量
    def extract_sn_content(struct):
        import re
        match = re.search(r'Sn(\d+)', struct, re.IGNORECASE)
        if match:
            return int(match.group(1))
        return 0
    
    df_valid['Sn含量'] = df_valid['结构'].apply(extract_sn_content)
    df_valid = df_valid.sort_values('Sn含量')
    
    # 按系列分组绘图
    series_list = sorted(df_valid['系列'].unique())
    
    for series in series_list:
        df_series = df_valid[df_valid['系列'] == series].copy()
        
        if len(df_series) == 0:
            continue
        
        fig, ax = plt.subplots(figsize=(12, 7))
        
        sn_content = df_series['Sn含量'].values
        structures = df_series['结构'].values
        T_low_to_mid = df_series['T_低到中流动性'].values
        T_mid_to_high = df_series['T_中到高流动性'].values
        
        # 绘制转变温度
        has_low_to_mid = ~np.isnan(T_low_to_mid)
        has_mid_to_high = ~np.isnan(T_mid_to_high)
        
        if np.any(has_low_to_mid):
            ax.plot(sn_content[has_low_to_mid], T_low_to_mid[has_low_to_mid],
                   'o-', label='低→中流动性转变', color='#3498DB',
                   linewidth=2, markersize=10)
        
        if np.any(has_mid_to_high):
            ax.plot(sn_content[has_mid_to_high], T_mid_to_high[has_mid_to_high],
                   's-', label='中→高流动性转变', color='#E74C3C',
                   linewidth=2, markersize=10)
        
        # 添加结构标签
        for i, struct in enumerate(structures):
            if not np.isnan(T_low_to_mid[i]):
                ax.text(sn_content[i], T_low_to_mid[i] - 20, struct,
                       ha='center', va='top', fontsize=8, rotation=45, alpha=0.7)
            elif not np.isnan(T_mid_to_high[i]):
                ax.text(sn_content[i], T_mid_to_high[i] + 20, struct,
                       ha='center', va='bottom', fontsize=8, rotation=45, alpha=0.7)
        
        ax.set_xlabel('Sn含量 (原子数)', fontsize=12, fontweight='bold')
        ax.set_ylabel('转变温度 (K)', fontsize=12, fontweight='bold')
        ax.set_title(f'流动性转变温度 vs Sn含量 - {series}系列\nMobility Transition Temperature vs Sn Content',
                    fontsize=14, fontweight='bold')
        ax.legend(fontsize=11, loc='best')
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        filename = OUTPUT_DIR / f'transition_temperature_{series}.png'
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"  [√] 已保存: {filename.name}")
        plt.close()


def generate_report(df_stats, df_transitions):
    """
    生成分析报告
    
    Parameters:
    -----------
    df_stats : DataFrame
        统计数据
    df_transitions : DataFrame
        转变温度数据
    """
    print("\n[*] 生成分析报告...")
    
    report_lines = []
    report_lines.append("=" * 80)
    report_lines.append("扩散系数(D值)分析报告")
    report_lines.append("Diffusion Coefficient (D-value) Analysis Report")
    report_lines.append("=" * 80)
    report_lines.append("")
    
    # 总体统计
    report_lines.append("【总体统计】")
    report_lines.append(f"  分析结构数: {df_stats['结构'].nunique()}")
    report_lines.append(f"  温度点数: {df_stats['温度'].nunique()}")
    report_lines.append(f"  数据点总数: {len(df_stats)}")
    report_lines.append("")
    
    # D值范围
    report_lines.append("【D值范围】")
    report_lines.append(f"  最小D值: {df_stats['D_mean'].min():.3e} cm²/s")
    report_lines.append(f"  最大D值: {df_stats['D_mean'].max():.3e} cm²/s")
    report_lines.append(f"  平均D值: {df_stats['D_mean'].mean():.3e} cm²/s")
    report_lines.append("")
    
    # 流动性分类统计
    report_lines.append("【流动性分类统计】")
    mobility_counts = df_stats['流动性'].value_counts()
    for mobility, count in mobility_counts.items():
        percentage = count / len(df_stats) * 100
        report_lines.append(f"  {mobility}: {count} ({percentage:.1f}%)")
    report_lines.append("")
    
    # 按系列统计
    report_lines.append("【按系列统计】")
    report_lines.append(f"{'系列':<10} {'结构数':<10} {'平均D值':<15} {'温度范围':<15}")
    report_lines.append("-" * 80)
    
    for series in sorted(df_stats['系列'].unique()):
        df_series = df_stats[df_stats['系列'] == series]
        n_structs = df_series['结构'].nunique()
        mean_D = df_series['D_mean'].mean()
        temp_range = f"{df_series['温度'].min():.0f}-{df_series['温度'].max():.0f}K"
        
        report_lines.append(f"{series:<10} {n_structs:<10} {mean_D:<15.3e} {temp_range:<15}")
    
    report_lines.append("")
    
    # 流动性转变温度统计
    if df_transitions is not None and len(df_transitions) > 0:
        report_lines.append("【流动性转变温度】")
        report_lines.append(f"{'结构':<15} {'低→中流动性':<20} {'中→高流动性':<20}")
        report_lines.append("-" * 80)
        
        for _, row in df_transitions.iterrows():
            struct = row['结构']
            T_low = f"{row['T_低到中流动性']:.1f}K" if not pd.isna(row['T_低到中流动性']) else "N/A"
            T_high = f"{row['T_中到高流动性']:.1f}K" if not pd.isna(row['T_中到高流动性']) else "N/A"
            report_lines.append(f"{struct:<15} {T_low:<20} {T_high:<20}")
        
        report_lines.append("")
    
    report_lines.append("=" * 80)
    
    # 保存报告
    report_file = OUTPUT_DIR / 'D_value_analysis_report.txt'
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write('\n'.join(report_lines))
    
    print(f"  [√] 报告已保存: {report_file.name}")
    
    # 打印到控制台
    print("\n" + '\n'.join(report_lines))


def main():
    """
    主函数
    
    命令行参数:
    -----------
    --exclude : 排除指定的结构,支持多个,例如: --exclude Pt8Sn0 Pt8Sn1
    
    示例:
    python step7_1_D_value_analysis.py                    # 分析所有结构
    python step7_1_D_value_analysis.py --exclude Pt8Sn0   # 排除Pt8Sn0
    python step7_1_D_value_analysis.py --exclude Pt8Sn0 Pt8Sn1  # 排除多个
    """
    # 解析命令行参数
    parser = argparse.ArgumentParser(
        description='Step7.1: 扩散系数(D值)分析'
    )
    parser.add_argument(
        '--exclude',
        nargs='+',
        default=[],
        help='排除指定的结构,支持多个结构名称'
    )
    parser.add_argument(
        '--no-filter',
        action='store_true',
        help='不进行MSD数据过滤,分析所有结构'
    )
    args = parser.parse_args()
    
    print("\n" + "="*80)
    print("Step7.1: 扩散系数(D值)分析")
    print("="*80)
    
    # 1. 加载数据
    exclude_list = [] if args.no_filter else args.exclude
    df = load_D_value_data(exclude_structures=exclude_list)
    
    if df is None or len(df) == 0:
        print("\n[X] 错误: 没有可分析的数据!")
        return
    
    # 2. 统计分析
    df_stats = analyze_D_values(df)
    
    # 3. 绘制热力图
    plot_D_value_heatmap(df_stats)
    
    # 4. 绘制D值-温度曲线 (类似林德曼分析)
    plot_D_vs_temperature(df_stats)
    
    # 5. 检测流动性转变温度
    df_transitions = detect_mobility_transitions(df_stats)
    
    # 6. 绘制转变温度分析图
    plot_transition_temperature_analysis(df_transitions)
    
    # 7. 生成报告
    generate_report(df_stats, df_transitions)
    
    print("\n" + "="*80)
    print("分析完成!")
    print(f"输出目录: {OUTPUT_DIR}")
    print("="*80)


if __name__ == '__main__':
    main()
