#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
Step4.1: 扩散系数(D值)综合分析
================================================================================

作者: GitHub Copilot  
日期: 2025-11-25
版本: v2.0 (简化版 - 仅CSV模式)

功能概述
========
从 step2 的 CSV 结果读取 D 值,进行统计分析和可视化:
- D-T 关系图 (按体系分组)
- D值热力图 (带R²标注)
- R²热力图 (拟合质量检查)
- 流动性分析 (低/中/高流动性)
- 转变温度检测

扩散系数判据
============
- D < 1e-8 cm²/s:  低流动性 (Low Mobility)
- 1e-8 < D < 1e-7: 中等流动性 (Medium Mobility)  
- D > 1e-7:        高流动性 (High Mobility)

使用方法
========
基本用法:
  python step4_1_D_value_comprehensive_analysis.py

高级选项:
  python step4_1_D_value_comprehensive_analysis.py --nofilter  # 使用完整数据(含异常)
  python step4_1_D_value_comprehensive_analysis.py --exclude Pt8Sn0  # 排除特定结构

输入数据
========
- ensemble_analysis_results.csv (默认,step2筛选后数据)
- ensemble_analysis_all.csv (--nofilter时使用,含异常数据)

输出文件
========
results/D_value_comprehensive/
├── ensemble_D_values.csv               # D值汇总表
├── D_value_statistics.csv              # 统计数据 (含R²字段)
├── D_vs_T_*.png                        # D-T关系图 (按体系)
├── D_value_heatmap_*.png               # D值热力图 (带R²标注)
├── R2_heatmap_*.png                    # R²热力图 (拟合质量检查)
├── D_vs_temperature_*.png              # D-T曲线 (所有结构)
├── transition_temperature_*.png        # 转变温度分析
├── mobility_transition_temperatures.csv # 转变温度数据
└── D_analysis_report.txt               # 综合报告
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from datetime import datetime
import re
import warnings
from scipy.interpolate import interp1d
import argparse
from tqdm import tqdm
warnings.filterwarnings('ignore')

# ===== 全局配置 =====
BASE_DIR = Path(__file__).parent
RESULTS_DIR = BASE_DIR / 'results'
ENSEMBLE_CSV = RESULTS_DIR / 'ensemble_analysis_results.csv'
ENSEMBLE_CSV_ALL = RESULTS_DIR / 'ensemble_analysis_all.csv'  # 完整数据 (含异常)

OUTPUT_DIR = RESULTS_DIR / 'D_value_comprehensive'
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

# ===== D值分类阈值 =====
D_LOW_THRESHOLD = 1e-8
D_HIGH_THRESHOLD = 1e-7

# ===== 配色方案 =====
COLORS = {
    'Pt': '#E74C3C',
    'Sn': '#3498DB',
    'PtSn': '#2ECC71'
}

plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False


# ==================== 数据加载模块 ====================

def load_D_from_csv(exclude_structures=None, enable_filtering=True):
    """从CSV读取D值
    
    参数:
        exclude_structures: 要排除的结构列表
        enable_filtering: True=使用筛选后数据, False=使用完整数据(含异常)
    """
    print(f"\n[1/6] 从 CSV 文件读取 D 值...")
    
    # 根据 enable_filtering 选择数据源
    if enable_filtering:
        csv_file = ENSEMBLE_CSV
        data_type = "筛选后数据"
    else:
        csv_file = ENSEMBLE_CSV_ALL
        data_type = "完整数据 (含异常)"
    
    if not csv_file.exists():
        print(f"[X] 错误: 未找到 {csv_file}")
        if not enable_filtering:
            print(f"[!] 提示: 请先运行 step2 生成完整数据文件")
            print(f"[!] 或使用默认模式 (不带 --nofilter)")
        return None
    
    df = pd.read_csv(csv_file)
    print(f"  [OK] 加载了 {len(df)} 条记录 ({data_type})")
    
    if exclude_structures:
        original_count = len(df)
        df = df[~df['composition'].isin(exclude_structures)]
        print(f"  [OK] 排除 {', '.join(exclude_structures)}: {original_count} -> {len(df)} 条")
    
    return df


# ==================== 数据处理模块 ====================

def extract_base_system(composition_name):
    """提取基础体系名称"""
    match = re.match(r'^(Cv)-\d+$', composition_name)
    if match:
        return match.group(1)
    return composition_name


def classify_series(structure_name):
    """分类体系所属系列"""
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
        
        if n_pt == 8:
            return 'Pt8SnX'
        if n_pt == 6:
            return 'Pt6SnX'
        if n_pt + n_sn == 8:
            return 'Pt(8-x)SnX'
    
    if pt_match and not sn_match:
        n_pt = int(pt_match.group(1))
        if n_pt == 8:
            return 'Pt8SnX'
        if n_pt == 6:
            return 'Pt6SnX'
    
    return 'Other'


def extract_atom_counts(composition):
    """从composition提取Pt, Sn, O原子数"""
    pt = sn = o = 0
    
    pt_match = re.search(r'Pt(\d+)', composition, re.IGNORECASE)
    if pt_match:
        pt = int(pt_match.group(1))
    
    sn_match = re.search(r'Sn(\d+)', composition, re.IGNORECASE)
    if sn_match:
        sn = int(sn_match.group(1))
    
    o_match = re.search(r'O(\d+)', composition, re.IGNORECASE)
    if o_match:
        o = int(o_match.group(1))
    
    return pt, sn, o


def enrich_dataframe(df):
    """为DataFrame添加额外的分析列"""
    # 确保 composition 是字符串类型
    df['composition'] = df['composition'].astype(str)
    
    # 提取原子数
    atom_counts = df['composition'].apply(extract_atom_counts)
    df['Pt原子数'] = atom_counts.apply(lambda x: x[0])
    df['Sn原子数'] = atom_counts.apply(lambda x: x[1])
    df['O原子数'] = atom_counts.apply(lambda x: x[2])
    
    # 分类系列
    df['系列'] = df['composition'].apply(classify_series)
    
    # 基础体系
    df['base_system'] = df['composition'].apply(extract_base_system)
    
    return df


def classify_mobility(D_value):
    """根据D值分类流动性"""
    if pd.isna(D_value):
        return '未知'
    elif D_value < D_LOW_THRESHOLD:
        return '低流动性'
    elif D_value < D_HIGH_THRESHOLD:
        return '中等流动性'
    else:
        return '高流动性'


def analyze_D_statistics(df):
    """分析D值统计"""
    print("\n[*] D值统计分析...")
    
    stats_list = []
    
    for struct in sorted(df['composition'].unique()):
        df_struct = df[df['composition'] == struct]
        
        row_sample = df_struct.iloc[0]
        series = row_sample['系列']
        pt = row_sample['Pt原子数']
        sn = row_sample['Sn原子数']
        o = row_sample['O原子数']
        
        for temp in sorted(df_struct['temp_value'].unique()):
            df_temp = df_struct[df_struct['temp_value'] == temp]
            
            D_Pt = df_temp[df_temp['element'] == 'Pt']['D_ensemble'].mean()
            D_Sn = df_temp[df_temp['element'] == 'Sn']['D_ensemble'].mean()
            D_PtSn = df_temp[df_temp['element'] == 'PtSn']['D_ensemble'].mean()
            D_mean = df_temp['D_ensemble'].mean()
            
            # 添加R²统计
            R2_Pt = df_temp[df_temp['element'] == 'Pt']['R2_ensemble'].mean()
            R2_Sn = df_temp[df_temp['element'] == 'Sn']['R2_ensemble'].mean()
            R2_PtSn = df_temp[df_temp['element'] == 'PtSn']['R2_ensemble'].mean()
            R2_mean = df_temp['R2_ensemble'].mean()
            
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
                'R2_Pt': R2_Pt,
                'R2_Sn': R2_Sn,
                'R2_PtSn': R2_PtSn,
                'R2_mean': R2_mean,
                '流动性': classify_mobility(D_mean)
            })
    
    df_stats = pd.DataFrame(stats_list)
    
    stats_file = OUTPUT_DIR / 'D_value_statistics.csv'
    df_stats.to_csv(stats_file, index=False, encoding='utf-8-sig')
    print(f"  [√] 统计结果已保存: {stats_file.name}")
    
    return df_stats


# ==================== 可视化模块 ====================

def plot_D_vs_T_by_system(df, output_dir):
    """绘制D-T关系图 (按体系分组)"""
    print("\n[*] 绘制 D-T 关系图 (按体系)...")
    
    systems = df['base_system'].unique()
    
    for base_sys in tqdm(systems, desc="绘图中"):
        df_sys = df[df['base_system'] == base_sys]
        
        if len(df_sys) == 0:
            continue
        
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))
        
        for idx, element in enumerate(['Pt', 'Sn', 'PtSn']):
            ax = axes[idx]
            df_elem = df_sys[df_sys['element'] == element]
            
            if len(df_elem) == 0:
                ax.text(0.5, 0.5, 'No Data', ha='center', va='center',
                       transform=ax.transAxes, fontsize=14)
                ax.set_title(f'{element}', fontsize=12, fontweight='bold')
                continue
            
            for comp in df_elem['composition'].unique():
                df_comp = df_elem[df_elem['composition'] == comp]
                df_comp = df_comp.sort_values('temp_value')
                
                ax.plot(df_comp['temp_value'], df_comp['D_ensemble'], 
                       'o-', label=comp, alpha=0.7, markersize=6)
            
            ax.set_xlabel('Temperature (K)', fontsize=11, fontweight='bold')
            ax.set_ylabel(r'D (cm$^2$/s)', fontsize=11, fontweight='bold')
            ax.set_title(f'{element}', fontsize=12, fontweight='bold', color=COLORS[element])
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.3)
            ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        
        plt.suptitle(f'{base_sys} - Diffusion Coefficient vs Temperature',
                    fontsize=14, fontweight='bold')
        plt.tight_layout()
        
        output_file = output_dir / f'D_vs_T_{base_sys}.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close(fig)
    
    print(f"  [√] 生成 {len(systems)} 张 D-T 图")


def plot_D_value_heatmap(df_stats, output_dir, filter_low_quality=False):
    """绘制D值热力图 (带R²标注)
    
    参数:
        filter_low_quality: 如果为True，R²<0.1或D<0的数据显示为空缺
    """
    if filter_low_quality:
        print("\n[*] 绘制D值热力图 (过滤低质量数据: R²<0.1 或 D<0)...")
        sub_dir = output_dir / 'filtered_heatmaps'
        sub_dir.mkdir(exist_ok=True)
        save_dir = sub_dir
    else:
        print("\n[*] 绘制D值热力图 (带R2标注)...")
        save_dir = output_dir
    
    series_list = sorted(df_stats['系列'].unique())
    
    for series in series_list:
        df_series = df_stats[df_stats['系列'] == series].copy()
        
        if len(df_series) == 0:
            continue
        
        # 排序结构
        struct_info = []
        for struct in df_series['结构'].unique():
            row = df_series[df_series['结构'] == struct].iloc[0]
            pt = int(row['Pt原子数'])
            sn = int(row['Sn原子数'])
            o = int(row['O原子数'])
            total = pt + sn + o
            pt_ratio = pt / (pt + sn) if (pt + sn) > 0 else 0
            struct_info.append((struct, pt, sn, o, total, pt_ratio))
        
        struct_info_sorted = sorted(struct_info, key=lambda x: (x[4], x[5]))
        sorted_structures = [s[0] for s in struct_info_sorted]
        
        # 准备D值热力图数据
        pivot_data = df_series.groupby(['温度', '结构'])['D_mean'].mean().reset_index()
        heatmap_data = pivot_data.pivot(index='温度', columns='结构', values='D_mean')
        heatmap_data = heatmap_data[sorted_structures]
        
        # 准备R²数据 (用于标注和过滤)
        pivot_r2 = df_series.groupby(['温度', '结构'])['R2_mean'].mean().reset_index()
        r2_data = pivot_r2.pivot(index='温度', columns='结构', values='R2_mean')
        r2_data = r2_data[sorted_structures]
        
        # 处理热力图数据
        heatmap_values = heatmap_data.values.copy()
        r2_values = r2_data.values.copy()
        
        if filter_low_quality:
            # 过滤模式：R²<0.1 或 D<0 的设为 NaN (显示为空缺)
            bad_mask = (heatmap_values <= 0) | (r2_values < 0.1) | np.isnan(r2_values)
            heatmap_values = np.where(bad_mask, np.nan, heatmap_values)
            heatmap_data_log = np.log10(np.where(np.isnan(heatmap_values), 1, heatmap_values))
            heatmap_data_log = np.where(bad_mask, np.nan, heatmap_data_log)
        else:
            # 原始模式：负D值设为极小值，保留所有数据
            negative_mask = heatmap_values < 0
            heatmap_values = np.where(heatmap_values > 0, heatmap_values, 1e-15)
            heatmap_data_log = np.log10(heatmap_values)
        
        # 绘图
        fig, ax = plt.subplots(figsize=(16, 9))
        
        # 使用带有空缺显示的 colormap
        cmap = plt.cm.viridis.copy()
        cmap.set_bad(color='lightgray')  # 空缺值显示为浅灰色
        
        im = ax.imshow(heatmap_data_log, aspect='auto', cmap=cmap,
                      interpolation='nearest' if filter_low_quality else 'bilinear', 
                      origin='lower', vmin=-15, vmax=-6)
        
        # 添加R²标注
        for i in range(len(heatmap_data.index)):
            for j in range(len(heatmap_data.columns)):
                r2_val = r2_data.iloc[i, j]
                d_val = heatmap_data.iloc[i, j]
                
                if filter_low_quality:
                    # 过滤模式：空缺格子不标注
                    if pd.isna(d_val) or d_val <= 0 or pd.isna(r2_val) or r2_val < 0.1:
                        continue
                    # 只标注有效数据的 R²
                    text_color = 'white' if r2_val >= 0.5 else 'orange'
                    text = f'{r2_val:.2f}'
                    ax.text(j, i, text, ha='center', va='center',
                           color=text_color, fontsize=8, fontweight='normal',
                           bbox=dict(boxstyle='round,pad=0.3', facecolor='black', 
                                   alpha=0.3, edgecolor='none'))
                else:
                    # 原始模式：显示所有标注
                    is_negative = d_val < 0 if pd.notna(d_val) else False
                    
                    if is_negative:
                        ax.text(j, i, '✗', ha='center', va='center',
                               color='red', fontsize=12, fontweight='bold')
                        continue
                    
                    if pd.notna(r2_val):
                        if r2_val < 0.1:
                            text_color = 'red'
                            text = f'R$^2$={r2_val:.2f}' if r2_val >= 0.01 else f'R$^2$={r2_val:.3f}'
                            fontweight = 'bold'
                        elif r2_val < 0.5:
                            text_color = 'orange'
                            text = f'{r2_val:.2f}'
                            fontweight = 'normal'
                        else:
                            text_color = 'white'
                            text = f'{r2_val:.2f}'
                            fontweight = 'normal'
                        
                        ax.text(j, i, text, ha='center', va='center',
                               color=text_color, fontsize=8, fontweight=fontweight,
                               bbox=dict(boxstyle='round,pad=0.3', facecolor='black', 
                                       alpha=0.3, edgecolor='none'))
        
        ax.set_xticks(np.arange(len(heatmap_data.columns)))
        ax.set_yticks(np.arange(len(heatmap_data.index)))
        ax.set_xticklabels(heatmap_data.columns, rotation=45, ha='right', fontsize=10)
        ax.set_yticklabels([f'{int(t)}K' for t in heatmap_data.index], fontsize=10)
        
        ax.set_xlabel('结构 (Structure)', fontsize=12, fontweight='bold')
        ax.set_ylabel('温度 (Temperature)', fontsize=12, fontweight='bold')
        
        if filter_low_quality:
            ax.set_title(f'扩散系数热力图 - {series}系列 (已过滤: R$^2$<0.1 或 D≤0)\nFiltered Diffusion Coefficient Heatmap',
                        fontsize=14, fontweight='bold', pad=15)
        else:
            ax.set_title(f'扩散系数热力图 - {series}系列 (标注R$^2$值)\nDiffusion Coefficient Heatmap with R$^2$ labels',
                        fontsize=14, fontweight='bold', pad=15)
        
        cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label(r'$\log_{10}$(D) [cm$^2$/s]', fontsize=11, fontweight='bold')
        
        # 添加图例说明
        from matplotlib.patches import Patch
        if filter_low_quality:
            legend_elements = [
                Patch(facecolor='lightgray', edgecolor='black', label='空缺 (R$^2$<0.1 或 D≤0)'),
                Patch(facecolor='orange', label=r'0.1 ≤ R$^2$ < 0.5'),
                Patch(facecolor='white', edgecolor='black', label=r'R$^2$ ≥ 0.5')
            ]
        else:
            legend_elements = [
                Patch(facecolor='red', label=r'R$^2$ < 0.1 (不可靠)'),
                Patch(facecolor='orange', label=r'0.1 ≤ R$^2$ < 0.5 (一般)'),
                Patch(facecolor='white', label=r'R$^2$ ≥ 0.5 (可靠)')
            ]
        ax.legend(handles=legend_elements, loc='upper left', fontsize=9,
                 framealpha=0.9, edgecolor='black')
        
        plt.tight_layout()
        
        if filter_low_quality:
            filename = save_dir / f'D_value_heatmap_filtered_{series}.png'
        else:
            filename = save_dir / f'D_value_heatmap_{series}.png'
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"  [√] 已保存: {filename.name}")
        plt.close()


def plot_D_vs_temperature_series(df_stats, output_dir):
    """绘制D值-温度曲线 (所有结构在一张图)"""
    print("\n[*] 绘制 D-T 曲线 (按系列汇总)...")
    
    series_list = sorted(df_stats['系列'].unique())
    
    for series in series_list:
        df_series = df_stats[df_stats['系列'] == series].copy()
        
        if len(df_series) == 0:
            continue
        
        structures = sorted(df_series['结构'].unique())
        
        fig, ax = plt.subplots(figsize=(12, 8))
        
        import matplotlib.cm as cm
        colors = cm.tab20(np.linspace(0, 1, len(structures)))
        
        for idx, struct in enumerate(structures):
            df_struct = df_series[df_series['结构'] == struct].sort_values('温度')
            temps = df_struct['温度'].values
            D_mean = df_struct['D_mean'].values
            
            ax.semilogy(temps, D_mean, 'o-', label=struct,
                       color=colors[idx], linewidth=2, markersize=6, alpha=0.7)
        
        ax.axhline(y=D_LOW_THRESHOLD, color='cyan', linestyle='--', 
                  linewidth=2, alpha=0.7, label=f'低流动性阈值: {D_LOW_THRESHOLD:.0e}')
        ax.axhline(y=D_HIGH_THRESHOLD, color='orange', linestyle='--',
                  linewidth=2, alpha=0.7, label=f'高流动性阈值: {D_HIGH_THRESHOLD:.0e}')
        
        ax.set_xlabel('温度 (K)', fontsize=13, fontweight='bold')
        ax.set_ylabel(r'D (cm$^2$/s)', fontsize=13, fontweight='bold')
        ax.set_title(f'扩散系数-温度关系 - {series}系列',
                    fontsize=15, fontweight='bold')
        ax.legend(fontsize=9, loc='best', ncol=2)
        ax.grid(True, alpha=0.3, which='both')
        
        plt.tight_layout()
        
        filename = output_dir / f'D_vs_temperature_{series}.png'
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"  [√] 已保存: {filename.name}")
        plt.close()


def plot_R2_heatmap(df_stats, output_dir):
    """绘制R²热力图 - 拟合质量检查"""
    print("\n[*] 绘制R2热力图 (拟合质量检查)...")
    
    series_list = sorted(df_stats['系列'].unique())
    
    for series in series_list:
        df_series = df_stats[df_stats['系列'] == series].copy()
        
        if len(df_series) == 0:
            continue
        
        # 排序结构
        struct_info = []
        for struct in df_series['结构'].unique():
            row = df_series[df_series['结构'] == struct].iloc[0]
            pt = int(row['Pt原子数'])
            sn = int(row['Sn原子数'])
            o = int(row['O原子数'])
            total = pt + sn + o
            pt_ratio = pt / (pt + sn) if (pt + sn) > 0 else 0
            struct_info.append((struct, pt, sn, o, total, pt_ratio))
        
        struct_info_sorted = sorted(struct_info, key=lambda x: (x[4], x[5]))
        sorted_structures = [s[0] for s in struct_info_sorted]
        
        # 准备R²热力图数据
        pivot_r2 = df_series.groupby(['温度', '结构'])['R2_mean'].mean().reset_index()
        r2_heatmap = pivot_r2.pivot(index='温度', columns='结构', values='R2_mean')
        r2_heatmap = r2_heatmap[sorted_structures]
        
        # 绘图
        fig, ax = plt.subplots(figsize=(16, 9))
        
        # 使用 RdYlGn 配色 (红-黄-绿)
        # R²低=红色(差), R²高=绿色(好)
        im = ax.imshow(r2_heatmap.values, aspect='auto', cmap='RdYlGn',
                      interpolation='bilinear', origin='lower',
                      vmin=0.0, vmax=1.0)
        
        # 添加数值标注
        for i in range(len(r2_heatmap.index)):
            for j in range(len(r2_heatmap.columns)):
                r2_val = r2_heatmap.iloc[i, j]
                
                if pd.notna(r2_val):
                    # 根据背景色选择文字颜色
                    text_color = 'white' if r2_val < 0.5 else 'black'
                    
                    # 格式化显示
                    if r2_val < 0.01:
                        text = f'{r2_val:.3f}'
                    else:
                        text = f'{r2_val:.2f}'
                    
                    # 标注
                    ax.text(j, i, text, ha='center', va='center',
                           color=text_color, fontsize=9, fontweight='bold')
        
        ax.set_xticks(np.arange(len(r2_heatmap.columns)))
        ax.set_yticks(np.arange(len(r2_heatmap.index)))
        ax.set_xticklabels(r2_heatmap.columns, rotation=45, ha='right', fontsize=10)
        ax.set_yticklabels([f'{int(t)}K' for t in r2_heatmap.index], fontsize=10)
        
        ax.set_xlabel('结构 (Structure)', fontsize=12, fontweight='bold')
        ax.set_ylabel('温度 (Temperature)', fontsize=12, fontweight='bold')
        ax.set_title(f'拟合质量热力图 (R$^2$) - {series}系列\nFit Quality Heatmap (R$^2$ values)',
                    fontsize=14, fontweight='bold', pad=15)
        
        cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label(r'R$^2$ (拟合优度)', fontsize=11, fontweight='bold')
        
        # 添加质量区间标注
        cbar.ax.axhline(y=0.1, color='red', linestyle='--', linewidth=2, alpha=0.7)
        cbar.ax.axhline(y=0.5, color='orange', linestyle='--', linewidth=2, alpha=0.7)
        cbar.ax.axhline(y=0.95, color='green', linestyle='--', linewidth=2, alpha=0.7)
        
        # 添加文字标注在colorbar旁边
        cbar.ax.text(1.5, 0.05, '差', fontsize=9, color='red', fontweight='bold')
        cbar.ax.text(1.5, 0.3, '一般', fontsize=9, color='orange', fontweight='bold')
        cbar.ax.text(1.5, 0.7, '良好', fontsize=9, color='green', fontweight='bold')
        cbar.ax.text(1.5, 0.97, '优秀', fontsize=9, color='darkgreen', fontweight='bold')
        
        plt.tight_layout()
        
        filename = output_dir / f'R2_heatmap_{series}.png'
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"  [√] 已保存: {filename.name}")
        plt.close()



def detect_mobility_transitions(df_stats, output_dir):
    """检测流动性转变温度"""
    print("\n[*] 检测流动性转变温度...")
    
    transitions = []
    
    for struct in sorted(df_stats['结构'].unique()):
        df_struct = df_stats[df_stats['结构'] == struct].copy()
        df_struct = df_struct.sort_values('温度')
        
        if len(df_struct) < 3:
            continue
        
        temps = df_struct['温度'].values
        D_mean = df_struct['D_mean'].values
        
        # 检测转变1: 低→中流动性
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
        
        # 检测转变2: 中→高流动性
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
        
        series = df_struct.iloc[0]['系列']
        transitions.append({
            '结构': struct,
            '系列': series,
            'T_低到中流动性': T_transition_low,
            'T_中到高流动性': T_transition_high
        })
        
        if T_transition_low is not None:
            print(f"  {struct}: 低->中 ~= {T_transition_low:.1f} K")
        if T_transition_high is not None:
            print(f"  {struct}: 中->高 ~= {T_transition_high:.1f} K")
    
    df_transitions = pd.DataFrame(transitions)
    
    trans_file = output_dir / 'mobility_transition_temperatures.csv'
    df_transitions.to_csv(trans_file, index=False, encoding='utf-8-sig')
    print(f"  [√] 转变温度数据已保存: {trans_file.name}")
    
    return df_transitions


def plot_transition_analysis(df_transitions, output_dir):
    """绘制流动性转变温度分析图"""
    print("\n[*] 绘制转变温度分析图...")
    
    df_valid = df_transitions[df_transitions['T_低到中流动性'].notna() | 
                             df_transitions['T_中到高流动性'].notna()].copy()
    
    if len(df_valid) == 0:
        print("  [!] 未检测到有效转变温度")
        return
    
    def extract_sn_content(struct):
        match = re.search(r'Sn(\d+)', struct, re.IGNORECASE)
        if match:
            return int(match.group(1))
        return 0
    
    df_valid['Sn含量'] = df_valid['结构'].apply(extract_sn_content)
    df_valid = df_valid.sort_values('Sn含量')
    
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
        
        has_low_to_mid = ~np.isnan(T_low_to_mid)
        has_mid_to_high = ~np.isnan(T_mid_to_high)
        
        if np.any(has_low_to_mid):
            ax.plot(sn_content[has_low_to_mid], T_low_to_mid[has_low_to_mid],
                   'o-', label='低→中流动性', color='#3498DB',
                   linewidth=2, markersize=10)
        
        if np.any(has_mid_to_high):
            ax.plot(sn_content[has_mid_to_high], T_mid_to_high[has_mid_to_high],
                   's-', label='中→高流动性', color='#E74C3C',
                   linewidth=2, markersize=10)
        
        ax.set_xlabel('Sn含量 (原子数)', fontsize=12, fontweight='bold')
        ax.set_ylabel('转变温度 (K)', fontsize=12, fontweight='bold')
        ax.set_title(f'流动性转变温度 vs Sn含量 - {series}系列',
                    fontsize=14, fontweight='bold')
        ax.legend(fontsize=11, loc='best')
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        filename = output_dir / f'transition_temperature_{series}.png'
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"  [√] 已保存: {filename.name}")
        plt.close()


# ==================== 报告生成 ====================

def generate_comprehensive_report(df, df_stats, df_transitions, output_dir):
    """生成综合分析报告"""
    print("  生成综合报告...")
    
    report_file = output_dir / 'D_analysis_report.txt'
    
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write("="*80 + "\n")
        f.write("扩散系数(D值)综合分析报告\n")
        f.write("Diffusion Coefficient Comprehensive Analysis Report\n")
        f.write("="*80 + "\n")
        f.write(f"生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"数据来源: step2 集成平均结果 (ensemble_analysis_results.csv)\n")
        f.write("="*80 + "\n\n")
        
        # 总体统计
        f.write("总体统计:\n")
        f.write(f"  数据组数: {len(df)}\n")
        f.write(f"  体系数: {df['base_system'].nunique()}\n")
        f.write(f"  组成数: {df['composition'].nunique()}\n")
        f.write(f"  温度点数: {df['temp_value'].nunique()}\n")
        f.write(f"  元素类型: {', '.join(df['element'].unique())}\n\n")
        
        # D值范围
        f.write("D值范围:\n")
        for element in ['Pt', 'Sn', 'PtSn']:
            df_elem = df[df['element'] == element]
            if len(df_elem) > 0:
                f.write(f"  {element}:\n")
                f.write(f"    最小: {df_elem['D_ensemble'].min():.2e} cm²/s\n")
                f.write(f"    最大: {df_elem['D_ensemble'].max():.2e} cm²/s\n")
                f.write(f"    平均: {df_elem['D_ensemble'].mean():.2e} cm²/s\n")
        f.write("\n")
        
        # 流动性分类
        f.write("流动性分类统计:\n")
        f.write(f"  阈值: 低 < {D_LOW_THRESHOLD:.0e} < 中 < {D_HIGH_THRESHOLD:.0e} < 高\n")
        mobility_counts = df_stats['流动性'].value_counts()
        for mobility, count in mobility_counts.items():
            percentage = count / len(df_stats) * 100
            f.write(f"  {mobility}: {count} ({percentage:.1f}%)\n")
        f.write("\n")
        
        # 系列统计
        f.write("按系列统计:\n")
        f.write(f"{'系列':<15} {'结构数':<10} {'平均D值':<15}\n")
        f.write("-"*80 + "\n")
        
        for series in sorted(df_stats['系列'].unique()):
            df_series = df_stats[df_stats['系列'] == series]
            n_structs = df_series['结构'].nunique()
            mean_D = df_series['D_mean'].mean()
            f.write(f"{series:<15} {n_structs:<10} {mean_D:<15.3e}\n")
        
        f.write("\n")
        
        # 转变温度
        if df_transitions is not None and len(df_transitions) > 0:
            f.write("流动性转变温度:\n")
            f.write(f"{'结构':<20} {'低→中':<15} {'中→高':<15}\n")
            f.write("-"*80 + "\n")
            
            for _, row in df_transitions.iterrows():
                struct = row['结构']
                T_low = f"{row['T_低到中流动性']:.1f}K" if not pd.isna(row['T_低到中流动性']) else "N/A"
                T_high = f"{row['T_中到高流动性']:.1f}K" if not pd.isna(row['T_中到高流动性']) else "N/A"
                f.write(f"{struct:<20} {T_low:<15} {T_high:<15}\n")
        
        f.write("\n" + "="*80 + "\n")
        f.write("输出文件:\n")
        f.write(f"  - ensemble_D_values.csv\n")
        f.write(f"  - D_value_statistics.csv\n")
        f.write(f"  - D_vs_T_*.png\n")
        f.write(f"  - D_value_heatmap_*.png\n")
        f.write(f"  - D_vs_temperature_*.png\n")
        f.write(f"  - transition_temperature_*.png\n")
        f.write("="*80 + "\n")
    
    print(f"  [√] 报告已保存: {report_file.name}")


# ==================== 主函数 ====================

def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description='Step4.1: 扩散系数(D值)综合分析',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  # 基本用法 (使用筛选后数据)
  python step4_1_D_value_comprehensive_analysis.py
  
  # 使用完整数据 (含异常run)
  python step4_1_D_value_comprehensive_analysis.py --nofilter
  
  # 排除特定结构
  python step4_1_D_value_comprehensive_analysis.py --exclude Pt8Sn0
  python step4_1_D_value_comprehensive_analysis.py --exclude Pt8Sn0 Pt6Sn0
        """
    )
    parser.add_argument(
        '--nofilter',
        action='store_true',
        help='使用完整数据 (含异常run),默认使用筛选后数据'
    )
    parser.add_argument(
        '--exclude',
        nargs='+',
        default=[],
        help='排除指定结构'
    )
    
    args = parser.parse_args()
    
    print("\n" + "="*80)
    print("Step4.1: 扩散系数(D值)综合分析")
    print("="*80)
    print(f"时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*80)
    
    # 加载数据
    df = load_D_from_csv(
        exclude_structures=args.exclude,
        enable_filtering=not args.nofilter
    )
    
    if df is None or len(df) == 0:
        print("\n[X] 错误: 没有可分析的数据!")
        return
    
    # 数据增强
    print("\n[2/6] 数据处理...")
    df = enrich_dataframe(df)
    
    # 保存D值表
    d_values_file = OUTPUT_DIR / 'ensemble_D_values.csv'
    df.to_csv(d_values_file, index=False, float_format='%.6e')
    print(f"  [OK] D值表已保存: {d_values_file.name} ({len(df)} 组)")
    
    # 统计分析
    print("\n[3/6] 统计分析...")
    df_stats = analyze_D_statistics(df)
    
    # 可视化
    print("\n[4/6] 生成可视化...")
    plot_D_vs_T_by_system(df, OUTPUT_DIR)
    plot_D_value_heatmap(df_stats, OUTPUT_DIR)  # 原始热图（含所有数据）
    plot_D_value_heatmap(df_stats, OUTPUT_DIR, filter_low_quality=True)  # 过滤后热图
    plot_R2_heatmap(df_stats, OUTPUT_DIR)
    plot_D_vs_temperature_series(df_stats, OUTPUT_DIR)
    
    # 转变温度分析
    print("\n[5/6] 转变温度分析...")
    df_transitions = detect_mobility_transitions(df_stats, OUTPUT_DIR)
    plot_transition_analysis(df_transitions, OUTPUT_DIR)
    
    # 生成报告
    print("\n[6/6] 生成综合报告...")
    generate_comprehensive_report(df, df_stats, df_transitions, OUTPUT_DIR)
    
    print("\n" + "="*80)
    print("分析完成!")
    print("="*80)
    print(f"输出目录: {OUTPUT_DIR}")
    print(f"  - ensemble_D_values.csv ({len(df)} 组)")
    print(f"  - D_value_statistics.csv ({len(df_stats)} 条)")
    print(f"  - D_vs_T_*.png (D-T关系图, 按体系)")
    print(f"  - D_value_heatmap_*.png (D值热力图, 带R2标注)")
    print(f"  - R2_heatmap_*.png (R2热力图, 拟合质量检查)")
    print(f"  - D_vs_temperature_*.png (D-T汇总曲线)")
    print(f"  - transition_temperature_*.png (流动性转变温度)")
    print(f"  - D_analysis_report.txt (综合报告)")
    print("="*80)


if __name__ == '__main__':
    main()
