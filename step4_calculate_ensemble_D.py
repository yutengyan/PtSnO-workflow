#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step4: 从集合平均MSD计算扩散系数D
==========================================

功能:
1. 读取step3生成的集合平均MSD数据
2. 对每个(composition, temperature, element)组合拟合D值
3. 支持系统过滤 (与step3一致的配置)
4. 生成详细的D值统计报告和可视化

输入:
- GMX MSD原始数据 (从file_index读取)
- step1的异常清单 (large_D_outliers.csv)

输出:
- ensemble_D_values.csv: D值汇总表
- D_vs_temperature_*.png: D-T关系图
- D_comparison_report.txt: 统计报告
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from datetime import datetime
import re
import warnings
from collections import defaultdict
from scipy import stats
from tqdm import tqdm
warnings.filterwarnings('ignore')

# ===== 全局配置 =====
BASE_DIR = Path(__file__).parent  # workflow目录
OUTLIERS_CSV = BASE_DIR / 'results' / 'large_D_outliers.csv'

GMX_DATA_DIRS = [
    BASE_DIR / 'data' / 'gmx_msd' / 'collected_gmx_msd',
    BASE_DIR / 'data' / 'gmx_msd' / 'gmx_msd_results_20251015_184626_collected',
    # 新版unwrap per-atom MSD数据 (2025-11-18)
    BASE_DIR / 'data' / 'gmx_msd' / 'unwrap' / 'gmx_msd_results_20251118_152614'
]

OUTPUT_DIR = BASE_DIR / 'results' / 'ensemble_D_analysis'
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

# ===== 系统过滤配置 (与step3保持一致) =====
SYSTEM_FILTER = {
    'include_patterns': [
        # 示例: 只处理Pt8相关系统
        r'^pt8',           # pt8开头的所有系统
        # r'pt8sn\d+',       # pt8sn0, pt8sn1, ..., pt8sn10
        # r'pt\d+sn\d+',     # 所有ptXsnY格式
    ],
    'exclude_patterns': [
        # 示例: 排除含氧系统
        # r'^[Oo]\d+',       # O1, O2, o3, o4等开头
        # r'[Oo]\d+Pt',      # O2Pt4Sn6等
        # r'Pt\d+Sn\d+O',    # Pt2Sn2O1等
    ]
}

# ===== 拟合参数 =====
FIT_CONFIG = {
    'fit_range': (50.0, 500.0),  # 拟合时间范围 (ps)
    'min_points': 10,             # 最小数据点数
    'r2_threshold': 0.95,         # R²阈值 (低于此值会标记)
}

COLORS = {
    'Pt': '#E74C3C',
    'Sn': '#3498DB',
    'PtSn': '#2ECC71'
}

plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False


def load_large_D_outliers():
    """加载大D值异常run清单"""
    try:
        df_outliers = pd.read_csv(OUTLIERS_CSV)
        outlier_files = set(df_outliers['filepath'].values)
        print(f"   [OK] Loaded outlier list: {len(outlier_files)} files")
        return outlier_files
    except FileNotFoundError:
        print(f"   [!] Outlier file not found, will use all runs")
        return set()


def build_file_index(outlier_files=None):
    """
    构建文件索引 - 只索引有效runs (排除异常runs)
    
    Returns:
    --------
    file_index : dict
        {(composition, temperature, element): [file_path1, file_path2, ...]}
    """
    if outlier_files is None:
        outlier_files = set()
    
    print("\n[*] Building file index (valid runs only)...")
    file_index = defaultdict(list)
    
    total_files = 0
    filtered_files = 0
    seen_files = set()  # 用于去重
    
    for gmx_dir in GMX_DATA_DIRS:
        print(f"  Scanning: {gmx_dir.name}...")
        
        for xvg_file in gmx_dir.rglob("*_msd_*.xvg"):
            # 去重检查
            try:
                normalized_path = xvg_file.resolve()
                if normalized_path in seen_files:
                    continue  # 跳过重复文件
                seen_files.add(normalized_path)
            except:
                pass  # 如果resolve()失败,继续处理
            try:
                parts = xvg_file.parts
                
                # 提取element
                filename = xvg_file.stem
                if '_msd_' in filename:
                    element = filename.split('_msd_')[-1]
                else:
                    continue
                
                # 提取temperature和composition
                temperature = None
                composition = None
                for i in range(len(parts)-1, 0, -1):
                    if parts[i].endswith('K'):
                        temperature = parts[i]
                        composition = parts[i-1]
                        break
                
                if not temperature or not composition:
                    continue
                
                total_files += 1
                
                # 只保留有效runs
                if str(xvg_file) not in outlier_files:
                    key = (composition, temperature, element)
                    file_index[key].append(xvg_file)
                else:
                    filtered_files += 1
                
            except Exception as e:
                continue
    
    print(f"   [OK] Indexed {total_files - filtered_files} valid files")
    print(f"   [OK] Filtered out {filtered_files} outlier runs")
    
    return file_index


def read_gmx_msd_xvg(filepath):
    """读取GMX MSD .xvg文件"""
    time_data = []
    msd_data = []
    
    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#') or line.startswith('@'):
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        t = float(parts[0])
                        msd_nm2 = float(parts[1])
                        msd_a2 = msd_nm2 * 100
                        time_data.append(t)
                        msd_data.append(msd_a2)
                    except ValueError:
                        continue
    except:
        return None, None
    
    if len(time_data) == 0:
        return None, None
    
    return np.array(time_data), np.array(msd_data)


def ensemble_average_msd(file_list):
    """
    计算集合平均MSD
    
    Parameters:
    -----------
    file_list : list
        文件路径列表
    
    Returns:
    --------
    time_common : ndarray
        公共时间轴
    msd_mean : ndarray
        集合平均MSD
    msd_std : ndarray
        标准差
    n_runs : int
        有效run数
    """
    if not file_list:
        return None, None, None, 0
    
    # 读取所有MSD数据
    msd_list = []
    time_list = []
    
    for filepath in file_list:
        time, msd = read_gmx_msd_xvg(filepath)
        if time is not None:
            time_list.append(time)
            msd_list.append(msd)
    
    if not msd_list:
        return None, None, None, 0
    
    # 对齐时间轴
    min_time = max(t.min() for t in time_list)
    max_time = min(t.max() for t in time_list)
    time_common = np.linspace(min_time, max_time, 500)
    
    # 插值
    msd_interp_list = []
    for time, msd in zip(time_list, msd_list):
        try:
            msd_interp = np.interp(time_common, time, msd)
            msd_interp_list.append(msd_interp)
        except:
            continue
    
    if not msd_interp_list:
        return None, None, None, 0
    
    # 计算统计量
    msd_array = np.array(msd_interp_list)
    msd_mean = np.mean(msd_array, axis=0)
    msd_std = np.std(msd_array, axis=0, ddof=1) if len(msd_interp_list) > 1 else np.zeros_like(msd_mean)
    
    return time_common, msd_mean, msd_std, len(msd_interp_list)


def fit_msd_to_diffusion(time, msd, fit_range=None):
    """
    拟合MSD到扩散方程: MSD = 6Dt + b
    
    Parameters:
    -----------
    time : ndarray
        时间 (ps)
    msd : ndarray
        MSD (A^2)
    fit_range : tuple, optional
        拟合时间范围 (t_min, t_max) in ps
    
    Returns:
    --------
    D : float
        扩散系数 (cm^2/s)
    r2 : float
        拟合优度
    intercept : float
        截距 (A^2)
    fit_info : dict
        拟合详细信息
    """
    if fit_range is None:
        fit_range = FIT_CONFIG['fit_range']
    
    # 选择拟合范围
    mask = (time >= fit_range[0]) & (time <= fit_range[1])
    
    if np.sum(mask) < FIT_CONFIG['min_points']:
        return np.nan, np.nan, np.nan, {'error': 'insufficient_points'}
    
    t_fit = time[mask]
    msd_fit = msd[mask]
    
    # 线性拟合
    try:
        slope, intercept, r_value, p_value, std_err = stats.linregress(t_fit, msd_fit)
        
        # 转换单位: MSD [A^2] = 6D [cm^2/s] * t [ps] * 10^-4
        # slope [A^2/ps] = 6D * 10^-4
        # D = slope / 6.0 * 1e-4 [cm^2/s]
        D = slope / 6.0 * 1e-4
        
        r2 = r_value ** 2
        
        fit_info = {
            'slope': slope,
            'intercept': intercept,
            'r2': r2,
            'p_value': p_value,
            'std_err': std_err,
            'n_points': len(t_fit),
            'fit_range': fit_range
        }
        
        return D, r2, intercept, fit_info
        
    except Exception as e:
        return np.nan, np.nan, np.nan, {'error': str(e)}


def extract_base_system(composition_name):
    """提取基础体系名称"""
    match = re.match(r'^(Cv)-\d+$', composition_name)
    if match:
        return match.group(1)
    return composition_name


def filter_systems(compositions, include_patterns=None, exclude_patterns=None):
    """
    根据正则表达式模式过滤系统
    
    Parameters:
    -----------
    compositions : list
        组成列表
    include_patterns : list of str
        包含模式列表
    exclude_patterns : list of str
        排除模式列表
    
    Returns:
    --------
    filtered_comps : list
        过滤后的组成列表
    """
    if include_patterns is None:
        include_patterns = []
    if exclude_patterns is None:
        exclude_patterns = []
    
    filtered_comps = []
    
    for comp in compositions:
        base_sys = extract_base_system(comp)
        
        # 检查include
        if include_patterns:
            included = any(re.search(pattern, base_sys, re.IGNORECASE) 
                          for pattern in include_patterns)
            if not included:
                continue
        
        # 检查exclude
        if exclude_patterns:
            excluded = any(re.search(pattern, base_sys, re.IGNORECASE) 
                          for pattern in exclude_patterns)
            if excluded:
                continue
        
        filtered_comps.append(comp)
    
    return filtered_comps


def calculate_ensemble_D_values(file_index, compositions):
    """
    计算所有(composition, temperature, element)组合的集合平均D值
    
    Parameters:
    -----------
    file_index : dict
        文件索引
    compositions : list
        要处理的组成列表
    
    Returns:
    --------
    results : list of dict
        结果列表
    """
    results = []
    
    # 获取所有唯一的keys
    keys = [k for k in file_index.keys() if k[0] in compositions]
    
    print(f"\n[*] Calculating ensemble D values for {len(keys)} groups...")
    
    for comp, temp, element in tqdm(keys, desc="Processing"):
        file_list = file_index[(comp, temp, element)]
        
        # 计算集合平均MSD
        time, msd_mean, msd_std, n_runs = ensemble_average_msd(file_list)
        
        if time is None or n_runs == 0:
            continue
        
        # 拟合D值
        D, r2, intercept, fit_info = fit_msd_to_diffusion(time, msd_mean)
        
        # 保存结果
        result = {
            'composition': comp,
            'base_system': extract_base_system(comp),
            'temperature': temp,
            'temp_K': int(temp.replace('K', '')),
            'element': element,
            'D_cm2_s': D,
            'r2': r2,
            'intercept_A2': intercept,
            'n_runs': n_runs,
            'slope': fit_info.get('slope', np.nan),
            'p_value': fit_info.get('p_value', np.nan),
            'std_err': fit_info.get('std_err', np.nan),
            'n_fit_points': fit_info.get('n_points', 0),
            'fit_range_ps': f"{FIT_CONFIG['fit_range'][0]}-{FIT_CONFIG['fit_range'][1]}"
        }
        
        results.append(result)
    
    return results


def plot_D_vs_temperature(df_results, output_dir):
    """
    绘制D-T关系图 (按体系分组)
    
    Parameters:
    -----------
    df_results : DataFrame
        结果数据
    output_dir : Path
        输出目录
    """
    print("\n[*] Plotting D vs Temperature...")
    
    # 按base_system分组
    systems = df_results['base_system'].unique()
    
    for base_sys in tqdm(systems, desc="Plotting"):
        df_sys = df_results[df_results['base_system'] == base_sys]
        
        if len(df_sys) == 0:
            continue
        
        # 创建图形
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))
        
        for idx, element in enumerate(['Pt', 'Sn', 'PtSn']):
            ax = axes[idx]
            
            df_elem = df_sys[df_sys['element'] == element]
            
            if len(df_elem) == 0:
                ax.text(0.5, 0.5, 'No Data', ha='center', va='center',
                       transform=ax.transAxes, fontsize=14)
                ax.set_title(f'{element}', fontsize=12, fontweight='bold')
                continue
            
            # 按composition分组绘制
            for comp in df_elem['composition'].unique():
                df_comp = df_elem[df_elem['composition'] == comp]
                df_comp = df_comp.sort_values('temp_K')
                
                # 绘制D值
                ax.plot(df_comp['temp_K'], df_comp['D_cm2_s'], 
                       'o-', label=comp, alpha=0.7, markersize=6)
            
            ax.set_xlabel('Temperature (K)', fontsize=11, fontweight='bold')
            ax.set_ylabel('D (cm²/s)', fontsize=11, fontweight='bold')
            ax.set_title(f'{element}', fontsize=12, fontweight='bold', color=COLORS[element])
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.3)
            ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        
        plt.suptitle(f'{base_sys} - Diffusion Coefficient vs Temperature\n'
                    f'Ensemble Average (Valid Runs Only)',
                    fontsize=14, fontweight='bold')
        
        plt.tight_layout()
        
        output_file = output_dir / f'D_vs_T_{base_sys}.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close(fig)
    
    print(f"  [OK] Generated {len(systems)} D-T plots")


def generate_report(df_results, output_dir):
    """生成统计报告"""
    report_file = output_dir / 'D_calculation_report.txt'
    
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write("="*80 + "\n")
        f.write("Step4: Ensemble Average Diffusion Coefficient Calculation Report\n")
        f.write("="*80 + "\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("="*80 + "\n\n")
        
        # 拟合参数
        f.write("Fitting Configuration:\n")
        f.write(f"  Fit range: {FIT_CONFIG['fit_range'][0]}-{FIT_CONFIG['fit_range'][1]} ps\n")
        f.write(f"  Min points: {FIT_CONFIG['min_points']}\n")
        f.write(f"  R² threshold: {FIT_CONFIG['r2_threshold']}\n\n")
        
        # 总体统计
        f.write("Overall Statistics:\n")
        f.write(f"  Total groups: {len(df_results)}\n")
        f.write(f"  Systems: {df_results['base_system'].nunique()}\n")
        f.write(f"  Compositions: {df_results['composition'].nunique()}\n")
        f.write(f"  Temperatures: {df_results['temperature'].nunique()}\n")
        f.write(f"  Elements: {', '.join(df_results['element'].unique())}\n\n")
        
        # R²统计
        low_r2 = df_results[df_results['r2'] < FIT_CONFIG['r2_threshold']]
        f.write(f"Fit Quality:\n")
        f.write(f"  R² ≥ {FIT_CONFIG['r2_threshold']}: {len(df_results) - len(low_r2)} groups ({(len(df_results)-len(low_r2))/len(df_results)*100:.1f}%)\n")
        f.write(f"  R² < {FIT_CONFIG['r2_threshold']}: {len(low_r2)} groups ({len(low_r2)/len(df_results)*100:.1f}%)\n")
        if len(low_r2) > 0:
            f.write(f"  Mean R² (low quality): {low_r2['r2'].mean():.4f}\n")
        f.write(f"  Mean R² (all): {df_results['r2'].mean():.4f}\n\n")
        
        # D值范围
        f.write("D Value Range:\n")
        for element in ['Pt', 'Sn', 'PtSn']:
            df_elem = df_results[df_results['element'] == element]
            if len(df_elem) > 0:
                f.write(f"  {element}:\n")
                f.write(f"    Min: {df_elem['D_cm2_s'].min():.2e} cm²/s\n")
                f.write(f"    Max: {df_elem['D_cm2_s'].max():.2e} cm²/s\n")
                f.write(f"    Mean: {df_elem['D_cm2_s'].mean():.2e} cm²/s\n")
                f.write(f"    Median: {df_elem['D_cm2_s'].median():.2e} cm²/s\n")
        
        f.write("\n" + "="*80 + "\n")
        f.write("Detailed Results saved to: ensemble_D_values.csv\n")
        f.write("="*80 + "\n")
    
    print(f"  [OK] Saved: {report_file.name}")


def main():
    print("\n" + "="*80)
    print("Step4: Calculate Ensemble Average Diffusion Coefficients")
    print("="*80)
    print(f"Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*80)
    
    # 1. 加载异常清单
    print("\n[1/5] Loading outlier list...")
    outlier_files = load_large_D_outliers()
    
    # 2. 构建文件索引
    file_index = build_file_index(outlier_files)
    
    # 3. 获取所有组成并应用过滤
    all_compositions = set(k[0] for k in file_index.keys())
    print(f"\n[2/5] Total compositions: {len(all_compositions)}")
    
    include_patterns = SYSTEM_FILTER.get('include_patterns', [])
    exclude_patterns = SYSTEM_FILTER.get('exclude_patterns', [])
    
    if include_patterns or exclude_patterns:
        print(f"\n[*] Applying system filters...")
        if include_patterns:
            print(f"    Include patterns: {include_patterns}")
        if exclude_patterns:
            print(f"    Exclude patterns: {exclude_patterns}")
        
        compositions = filter_systems(list(all_compositions), include_patterns, exclude_patterns)
        print(f"    Filtered: {len(compositions)} compositions selected")
    else:
        compositions = list(all_compositions)
    
    systems = set(extract_base_system(c) for c in compositions)
    print(f"    Systems to process: {len(systems)}")
    
    # 4. 计算集合平均D值
    print(f"\n[3/5] Calculating ensemble D values...")
    results = calculate_ensemble_D_values(file_index, compositions)
    
    if not results:
        print("[X] No results generated!")
        return
    
    # 5. 保存结果
    print(f"\n[4/5] Saving results...")
    df_results = pd.DataFrame(results)
    
    # 按体系、温度、元素排序
    df_results = df_results.sort_values(['base_system', 'temp_K', 'element'])
    
    output_csv = OUTPUT_DIR / 'ensemble_D_values.csv'
    df_results.to_csv(output_csv, index=False, float_format='%.6e')
    print(f"  [OK] Saved: {output_csv}")
    print(f"       Total groups: {len(df_results)}")
    
    # 6. 生成可视化和报告
    print(f"\n[5/5] Generating visualizations and report...")
    plot_D_vs_temperature(df_results, OUTPUT_DIR)
    generate_report(df_results, OUTPUT_DIR)
    
    print("\n" + "="*80)
    print("DONE!")
    print("="*80)
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"  - ensemble_D_values.csv: D值汇总表 ({len(df_results)} groups)")
    print(f"  - D_vs_T_*.png: D-T关系图 ({len(systems)} systems)")
    print(f"  - D_calculation_report.txt: 统计报告")
    print("="*80)


if __name__ == '__main__':
    main()
