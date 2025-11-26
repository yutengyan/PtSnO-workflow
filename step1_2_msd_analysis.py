#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
=============================================================================
Step1.2: GMX MSD综合分析 - 异常检测 + 集合平均
=============================================================================
创建时间: 2025-11-26
作者: GitHub Copilot

功能说明:
    整合原 step1_detect_outliers.py 和 step2_ensemble_analysis.py 的功能:
    
    [阶段1] Run级别质量筛选 (原step1)
        - 扫描所有 .xvg 文件
        - 提取 GMX D值 (从注释行)
        - IQR/3sigma/MAD 异常检测
        - Intercept 阈值筛选
        - D值上限筛选
        - 同模拟联动筛选 (一个元素异常则该模拟全部元素标记)
    
    [阶段2] 集合平均分析 (原step2)
        - 对筛选后的 run 进行 MSD 曲线平均
        - 拟合得到扩散系数 D
        - 分别保存有效数据和异常数据

输出文件: results/
    ├── large_D_outliers.csv           (异常run清单)
    ├── ensemble_analysis_results.csv  (筛选后的集合平均结果)
    ├── ensemble_analysis_filtered.csv (被筛选runs的集合平均)
    ├── ensemble_analysis_all.csv      (完整数据,用于--nofilter)
    ├── ensemble_comparison.csv        (筛选前后对比)
    └── run_quality_report.txt         (质量报告)

使用方法:
    python step1_2_msd_analysis.py      # 完整流程
=============================================================================
"""

import pandas as pd
import numpy as np
from scipy.stats import linregress
from pathlib import Path
import os
import re
from collections import defaultdict
import sys

# 设置UTF-8输出(Windows兼容)
if sys.platform == 'win32':
    try:
        sys.stdout.reconfigure(encoding='utf-8')
    except:
        import os
        os.environ['PYTHONIOENCODING'] = 'utf-8'

# ==================== 可调参数 ====================

# 基础目录
BASE_DIR = Path(__file__).parent

# 数据源目录
GMX_DATA_DIRS = [
    BASE_DIR / 'data' / 'gmx_msd' / 'unwrap' / 'air' / 'gmx_msd_results_20251124_170114'
]

# 输出目录
OUTPUT_DIR = BASE_DIR / 'results'

# ==================== 异常检测参数 ====================

# 异常检测方法: 'IQR', '3sigma', 'MAD'
OUTLIER_METHOD = 'IQR'
IQR_MULTIPLIER = 1.5
SIGMA_MULTIPLIER = 3.0
MAD_MULTIPLIER = 3.0

# 最小run数量(少于此数量不做异常检测)
MIN_RUNS_FOR_OUTLIER = 3

# Intercept阈值 (A^2) - 筛除结构松弛异常
INTERCEPT_MAX = 10.0

# D值上限 (cm^2/s) - 筛除团簇漂移异常
D_MAX_THRESHOLD = 5e-5

# ==================== 拟合参数 ====================

FIT_START_RATIO = 0.1
FIT_END_RATIO = 0.9

# ==================== 核心函数 ====================

def extract_gmx_D_value(filepath):
    """
    从.xvg文件注释中提取GMX计算的D值
    
    查找格式: # D[element] = 0.2585 (+/- 0.7552) (1e-5 cm^2/s)
    
    返回: D值(cm^2/s),如果未找到返回None
    """
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            for line in f:
                if not (line.startswith('#') or line.startswith('@')):
                    break
                
                if 'D[' in line and '=' in line and 'cm' in line:
                    match = re.search(r'D\[\s*(\w+)\s*\]\s*=\s*([+-]?[\d.]+)', line)
                    if match:
                        D_value = float(match.group(2))
                        
                        unit_match = re.search(r'\(1e([+-]?\d+)\s*cm', line)
                        if unit_match:
                            exponent = int(unit_match.group(1))
                            D_value *= 10**exponent
                        else:
                            D_value *= 1e-5
                        
                        return D_value
        return None
    except Exception:
        return None


def read_xvg(filepath):
    """
    读取GMX .xvg文件
    
    单位转换: MSD (nm^2 -> A^2), Time (ps)
    """
    try:
        data = []
        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith('#') or line.startswith('@'):
                    continue
                parts = line.strip().split()
                if len(parts) >= 2:
                    try:
                        t = float(parts[0])
                        msd = float(parts[1]) * 100  # nm^2 -> A^2
                        data.append([t, msd])
                    except ValueError:
                        continue
        
        if len(data) == 0:
            return None, None
        
        data = np.array(data)
        return data[:, 0], data[:, 1]
    except Exception:
        return None, None


def parse_gmx_filename(filepath):
    """
    解析GMX MSD文件名
    
    路径示例: .../68/1000K/T1000.r24.gpu0_msd_Pt.xvg
    
    返回: dict with composition, temperature, element, run_id, filepath
    """
    path = Path(filepath)
    filename = path.name
    
    # 提取元素
    if filename.endswith('_Pt.xvg'):
        element = 'Pt'
    elif filename.endswith('_Sn.xvg'):
        element = 'Sn'
    elif filename.endswith('_PtSn.xvg'):
        element = 'PtSn'
    else:
        return None
    
    # 提取温度和run_id: T1000.r24.gpu0_msd_Pt.xvg
    match = re.match(r'T(\d+)\.r(\d+)\.gpu\d+_msd_\w+\.xvg', filename)
    if not match:
        return None
    
    temp_value = int(match.group(1))
    run_id = int(match.group(2))
    temperature = f"{temp_value}K"
    
    # 提取组成 (温度目录的父目录)
    parts = path.parts
    temp_dir_idx = None
    for i, part in enumerate(parts):
        if part == temperature:
            temp_dir_idx = i
            break
    
    if temp_dir_idx is None or temp_dir_idx == 0:
        return None
    
    composition = parts[temp_dir_idx - 1]
    
    return {
        'composition': composition,
        'temperature': temperature,
        'temp_value': temp_value,
        'element': element,
        'run_id': run_id,
        'filepath': path,
        'full_key': f"{composition}_{temperature}"
    }


def fit_msd_to_diffusion(time, msd, fit_start_ratio=0.1, fit_end_ratio=0.9):
    """
    线性拟合MSD -> 扩散系数
    
    输入: time (ps), msd (A^2)
    输出: D (cm^2/s), r2, intercept (A^2)
    """
    if len(time) < 10:
        return None, None, None
    
    start_idx = int(len(time) * fit_start_ratio)
    end_idx = int(len(time) * fit_end_ratio)
    
    if end_idx <= start_idx:
        return None, None, None
    
    t_fit = time[start_idx:end_idx]
    msd_fit = msd[start_idx:end_idx]
    
    try:
        slope, intercept, r_value, p_value, std_err = linregress(t_fit, msd_fit)
        D = slope / 6.0 * 1e-4  # A^2/ps -> cm^2/s
        r2 = r_value ** 2
        return D, r2, intercept
    except Exception:
        return None, None, None


def detect_outliers_iqr(values):
    """使用IQR方法检测异常值"""
    values = np.array(values)
    Q1 = np.percentile(values, 25)
    Q3 = np.percentile(values, 75)
    IQR = Q3 - Q1
    
    lower_bound = Q1 - IQR_MULTIPLIER * IQR
    upper_bound = Q3 + IQR_MULTIPLIER * IQR
    
    outlier_indices = np.where((values < lower_bound) | (values > upper_bound))[0]
    
    stats = {
        'Q1': Q1, 'Q3': Q3, 'IQR': IQR,
        'lower_bound': lower_bound, 'upper_bound': upper_bound,
        'median': np.median(values), 'mean': np.mean(values)
    }
    
    return outlier_indices.tolist(), stats


def detect_outliers_3sigma(values):
    """使用3sigma方法检测异常值"""
    values = np.array(values)
    mean = np.mean(values)
    std = np.std(values)
    
    lower_bound = mean - SIGMA_MULTIPLIER * std
    upper_bound = mean + SIGMA_MULTIPLIER * std
    
    outlier_indices = np.where((values < lower_bound) | (values > upper_bound))[0]
    
    stats = {
        'mean': mean, 'std': std,
        'lower_bound': lower_bound, 'upper_bound': upper_bound
    }
    
    return outlier_indices.tolist(), stats


def detect_outliers_mad(values):
    """使用MAD方法检测异常值"""
    values = np.array(values)
    median = np.median(values)
    mad = np.median(np.abs(values - median))
    
    if mad == 0:
        return [], {'median': median, 'mad': 0}
    
    lower_bound = median - MAD_MULTIPLIER * mad
    upper_bound = median + MAD_MULTIPLIER * mad
    
    outlier_indices = np.where((values < lower_bound) | (values > upper_bound))[0]
    
    stats = {
        'median': median, 'mad': mad,
        'lower_bound': lower_bound, 'upper_bound': upper_bound
    }
    
    return outlier_indices.tolist(), stats


def ensemble_average_msd(file_list):
    """集合平均多个run的MSD"""
    all_times = []
    all_msds = []
    
    for filepath in file_list:
        time, msd = read_xvg(filepath)
        if time is not None and msd is not None:
            all_times.append(time)
            all_msds.append(msd)
    
    if len(all_times) == 0:
        return None
    
    # 对齐到最短长度
    min_length = min(len(t) for t in all_times)
    time_avg = all_times[0][:min_length]
    
    msd_arrays = [msd[:min_length] for msd in all_msds]
    msd_avg = np.mean(msd_arrays, axis=0)
    msd_std = np.std(msd_arrays, axis=0, ddof=1) if len(msd_arrays) > 1 else np.zeros_like(msd_avg)
    
    return {
        'time': time_avg,
        'msd_mean': msd_avg,
        'msd_std': msd_std,
        'n_runs': len(all_times)
    }


# ==================== 主程序 ====================

def main():
    """主程序 - 整合异常检测和集合平均"""
    
    print("=" * 80)
    print("GMX MSD综合分析 - Step1.2 (异常检测 + 集合平均)".center(80))
    print("=" * 80)
    print()
    
    # 创建输出目录
    Path(OUTPUT_DIR).mkdir(exist_ok=True, parents=True)
    
    # ========================================================================
    # 阶段1: 扫描和解析文件
    # ========================================================================
    print("[1/7] 扫描.xvg文件...")
    
    all_files = []
    for data_dir in GMX_DATA_DIRS:
        if not Path(data_dir).exists():
            print(f"  [!] 目录不存在: {data_dir}")
            continue
        
        files = list(Path(data_dir).rglob('*_msd_*.xvg'))
        all_files.extend(files)
        print(f"  [OK] {Path(data_dir).name}: {len(files)} 个文件")
    
    print(f"\n  总计: {len(all_files)} 个.xvg文件")
    
    if len(all_files) == 0:
        print("\n[ERROR] 未找到任何.xvg文件!")
        return
    
    # ========================================================================
    # 阶段2: 提取GMX D值和run信息
    # ========================================================================
    print("\n[2/7] 提取GMX D值和run信息...")
    
    run_data = []
    fail_stats = {'parse': 0, 'gmx_D': 0, 'msd': 0, 'fit': 0}
    
    for i, filepath in enumerate(all_files):
        if (i + 1) % 100 == 0:
            print(f"  处理中: {i+1}/{len(all_files)}", end='\r')
        
        info = parse_gmx_filename(filepath)
        if not info:
            fail_stats['parse'] += 1
            continue
        
        gmx_D = extract_gmx_D_value(filepath)
        if gmx_D is None:
            fail_stats['gmx_D'] += 1
            continue
        
        time, msd = read_xvg(filepath)
        if time is None:
            fail_stats['msd'] += 1
            continue
        
        our_D, our_r2, intercept = fit_msd_to_diffusion(time, msd, FIT_START_RATIO, FIT_END_RATIO)
        if our_D is None:
            fail_stats['fit'] += 1
            continue
        
        run_data.append({
            'composition': info['composition'],
            'temperature': info['temperature'],
            'temp_value': info['temp_value'],
            'element': info['element'],
            'run_id': info['run_id'],
            'filepath': str(info['filepath'].resolve()),  # 使用绝对路径
            'gmx_D': gmx_D,
            'our_D': our_D,
            'our_r2': our_r2,
            'intercept': intercept,
            'time': time,
            'msd': msd
        })
    
    print(f"\n  [OK] 成功提取: {len(run_data)} 个run")
    print(f"  [INFO] 失败统计: parse={fail_stats['parse']}, gmx_D={fail_stats['gmx_D']}, "
          f"msd={fail_stats['msd']}, fit={fail_stats['fit']}")
    
    if len(run_data) == 0:
        print("\n[ERROR] 没有提取到任何run数据!")
        return
    
    df_runs = pd.DataFrame(run_data)
    
    # ========================================================================
    # 阶段3: 按组成-温度-元素分组
    # ========================================================================
    print("\n[3/7] 按组成-温度-元素分组...")
    
    df_runs['group_key'] = (df_runs['composition'] + '_' + 
                            df_runs['temperature'] + '_' + 
                            df_runs['element'])
    
    grouped = df_runs.groupby('group_key')
    print(f"  [OK] 共 {len(grouped)} 个唯一组合")
    
    # ========================================================================
    # 阶段4: run级别质量筛选
    # ========================================================================
    print(f"\n[4/7] Run级别质量筛选...")
    print(f"  方法1: {OUTLIER_METHOD}异常检测")
    print(f"  方法2: Intercept <= {INTERCEPT_MAX} A^2")
    print(f"  方法3: D < {D_MAX_THRESHOLD} cm^2/s")
    
    outlier_records = []
    iqr_filtered_count = 0
    intercept_filtered_count = 0
    d_filtered_count = 0
    total_runs_checked = 0
    
    for group_key, group_df in grouped:
        if len(group_df) < MIN_RUNS_FOR_OUTLIER:
            continue
        
        gmx_D_values = group_df['gmx_D'].values
        
        # 方法1: 统计异常检测
        if OUTLIER_METHOD == 'IQR':
            outlier_indices, stats = detect_outliers_iqr(gmx_D_values)
        elif OUTLIER_METHOD == '3sigma':
            outlier_indices, stats = detect_outliers_3sigma(gmx_D_values)
        elif OUTLIER_METHOD == 'MAD':
            outlier_indices, stats = detect_outliers_mad(gmx_D_values)
        else:
            outlier_indices, stats = detect_outliers_iqr(gmx_D_values)
        
        # 记录统计异常
        for idx in outlier_indices:
            row = group_df.iloc[idx]
            outlier_records.append({
                'group_key': group_key,
                'composition': row['composition'],
                'temperature': row['temperature'],
                'element': row['element'],
                'run_id': row['run_id'],
                'gmx_D': row['gmx_D'],
                'filepath': row['filepath'],
                'reason': f"{OUTLIER_METHOD}_outlier"
            })
            iqr_filtered_count += 1
        
        # 方法2+3: Intercept和D值筛选
        for idx in range(len(group_df)):
            total_runs_checked += 1
            
            if idx in outlier_indices:
                continue
            
            row = group_df.iloc[idx]
            intercept_run = row['intercept'] if 'intercept' in row else None
            D_run = row['our_D']
            
            # Intercept筛选
            if intercept_run is not None and intercept_run > INTERCEPT_MAX:
                outlier_records.append({
                    'group_key': group_key,
                    'composition': row['composition'],
                    'temperature': row['temperature'],
                    'element': row['element'],
                    'run_id': row['run_id'],
                    'gmx_D': row['gmx_D'],
                    'filepath': row['filepath'],
                    'reason': f'Intercept>{INTERCEPT_MAX}A^2'
                })
                intercept_filtered_count += 1
                continue
            
            # D值上限筛选
            if D_run is not None and D_run >= D_MAX_THRESHOLD:
                outlier_records.append({
                    'group_key': group_key,
                    'composition': row['composition'],
                    'temperature': row['temperature'],
                    'element': row['element'],
                    'run_id': row['run_id'],
                    'gmx_D': row['gmx_D'],
                    'filepath': row['filepath'],
                    'reason': f'D>={D_MAX_THRESHOLD}cm^2/s'
                })
                d_filtered_count += 1
    
    print(f"\n  [OK] 初步筛除: {len(outlier_records)} run / {total_runs_checked} run")
    print(f"    - {OUTLIER_METHOD}异常检测: {iqr_filtered_count} run")
    print(f"    - Intercept>{INTERCEPT_MAX}A^2: {intercept_filtered_count} run")
    print(f"    - D>={D_MAX_THRESHOLD}cm^2/s: {d_filtered_count} run")
    
    # ========================================================================
    # 阶段4.5: 同模拟联动筛选
    # ========================================================================
    print(f"\n[4.5/7] 同模拟联动筛选...")
    
    bad_simulation_set = set()
    for rec in outlier_records:
        bad_simulation_set.add((rec['composition'], rec['temperature'], rec['run_id']))
    
    print(f"  [INFO] 初步标记坏模拟: {len(bad_simulation_set)} 个")
    
    additional_outliers = []
    for _, row in df_runs.iterrows():
        sim_key = (row['composition'], row['temperature'], row['run_id'])
        
        if sim_key in bad_simulation_set:
            already_marked = any(
                rec['composition'] == row['composition'] and
                rec['temperature'] == row['temperature'] and
                rec['run_id'] == row['run_id'] and
                rec['element'] == row['element']
                for rec in outlier_records
            )
            
            if not already_marked:
                additional_outliers.append({
                    'group_key': f"{row['composition']}_{row['temperature']}_{row['element']}",
                    'composition': row['composition'],
                    'temperature': row['temperature'],
                    'element': row['element'],
                    'run_id': row['run_id'],
                    'gmx_D': row['gmx_D'],
                    'filepath': row['filepath'],
                    'reason': 'linked_bad_simulation'
                })
    
    n_additional = len(additional_outliers)
    outlier_records.extend(additional_outliers)
    
    print(f"  [OK] 联动标记额外 {n_additional} 个元素run")
    print(f"  [OK] 最终筛除: {len(outlier_records)} run (初步{len(outlier_records)-n_additional} + 联动{n_additional})")
    
    # 保存异常清单 (供其他脚本如step3/step4使用)
    if len(outlier_records) > 0:
        df_outliers = pd.DataFrame(outlier_records)
        outlier_file = Path(OUTPUT_DIR) / 'large_D_outliers.csv'
        df_outliers.to_csv(outlier_file, index=False, encoding='utf-8-sig')
        print(f"  [OK] 异常run清单: {outlier_file} (供step3/step4等脚本使用)")
    
    # 构建异常文件路径集合 (小写,用于匹配)
    outlier_filepath_set = set(rec['filepath'].lower() for rec in outlier_records)
    
    # ========================================================================
    # 阶段5: 集合平均分析
    # ========================================================================
    print(f"\n[5/7] 集合平均分析...")
    print(f"  方法: MSD曲线平均 -> 拟合 -> D值")
    
    # 重新按文件分组
    file_grouped = defaultdict(list)
    for _, row in df_runs.iterrows():
        key = f"{row['composition']}_{row['temperature']}_{row['element']}"
        file_grouped[key].append({
            'filepath': Path(row['filepath']),
            'run_id': row['run_id']
        })
    
    results = []
    results_filtered = []
    results_all = []  # 全部runs的集合平均 (不筛选)
    comparison_results = []
    
    for key, file_info_list in file_grouped.items():
        parts = key.rsplit('_', 2)
        if len(parts) != 3:
            continue
        
        composition = parts[0]
        temperature = parts[1]
        element = parts[2]
        temp_value = int(temperature.replace('K', ''))
        
        all_files_in_group = [f['filepath'] for f in file_info_list]
        original_n = len(all_files_in_group)
        
        # 分离有效和异常文件
        valid_files = [f for f in all_files_in_group 
                       if str(f.resolve()).lower() not in outlier_filepath_set]
        filtered_files = [f for f in all_files_in_group 
                          if str(f.resolve()).lower() in outlier_filepath_set]
        
        # === 对比分析 (全部run vs 筛选后) ===
        # 全部run集合平均
        avg_all = ensemble_average_msd(all_files_in_group)
        D_all, r2_all, intercept_all = (None, None, None)
        if avg_all:
            D_all, r2_all, intercept_all = fit_msd_to_diffusion(
                avg_all['time'], avg_all['msd_mean'], FIT_START_RATIO, FIT_END_RATIO
            )
        
        # === 保存全部runs结果 (不筛选, 用于--nofilter) ===
        if D_all is not None:
            results_all.append({
                'composition': composition,
                'temperature': temperature,
                'temp_value': temp_value,
                'element': element,
                'n_runs': original_n,
                'n_runs_original': original_n,
                'n_runs_filtered': 0,  # 不筛选
                'D_ensemble': D_all,
                'R2_ensemble': r2_all,
                'intercept': intercept_all
            })
        
        # 筛选后集合平均
        D_clean, r2_clean = D_all, r2_all
        if len(valid_files) >= 2:
            avg_clean = ensemble_average_msd(valid_files)
            if avg_clean:
                D_clean, r2_clean, _ = fit_msd_to_diffusion(
                    avg_clean['time'], avg_clean['msd_mean'], FIT_START_RATIO, FIT_END_RATIO
                )
        
        # 记录对比结果
        if D_all is not None and D_clean is not None:
            delta_r2 = (r2_clean - r2_all) if r2_all else 0
            comparison_results.append({
                'composition': composition,
                'temperature': temperature,
                'element': element,
                'n_runs_total': original_n,
                'n_runs_outlier': len(filtered_files),
                'n_runs_clean': len(valid_files),
                'D_all': D_all,
                'r2_all': r2_all,
                'D_clean': D_clean,
                'r2_clean': r2_clean,
                'delta_r2': delta_r2,
                'improved': delta_r2 > 0
            })
        
        # === 保存有效runs结果 ===
        if len(valid_files) > 0:
            avg_result = ensemble_average_msd(valid_files)
            if avg_result:
                D, r2, intercept = fit_msd_to_diffusion(
                    avg_result['time'], avg_result['msd_mean'],
                    FIT_START_RATIO, FIT_END_RATIO
                )
                if D is not None:
                    results.append({
                        'composition': composition,
                        'temperature': temperature,
                        'temp_value': temp_value,
                        'element': element,
                        'n_runs': len(valid_files),
                        'n_runs_original': original_n,
                        'n_runs_filtered': len(filtered_files),
                        'D_ensemble': D,
                        'R2_ensemble': r2,
                        'intercept': intercept
                    })
        
        # === 保存异常runs结果 ===
        if len(filtered_files) > 0:
            avg_filt = ensemble_average_msd(filtered_files)
            if avg_filt:
                D_filt, r2_filt, intercept_filt = fit_msd_to_diffusion(
                    avg_filt['time'], avg_filt['msd_mean'],
                    FIT_START_RATIO, FIT_END_RATIO
                )
                if D_filt is not None:
                    results_filtered.append({
                        'composition': composition,
                        'temperature': temperature,
                        'temp_value': temp_value,
                        'element': element,
                        'n_runs': len(filtered_files),
                        'D_ensemble': D_filt,
                        'R2_ensemble': r2_filt,
                        'intercept': intercept_filt
                    })
    
    print(f"  [OK] 有效组合: {len(results)} 组")
    print(f"  [OK] 异常组合: {len(results_filtered)} 组")
    
    # ========================================================================
    # 阶段6: 保存结果
    # ========================================================================
    print(f"\n[6/7] 保存结果...")
    
    # 有效runs结果
    df_results = pd.DataFrame(results)
    output_file = Path(OUTPUT_DIR) / 'ensemble_analysis_results.csv'
    df_results.to_csv(output_file, index=False, encoding='utf-8-sig')
    total_valid_runs = df_results['n_runs'].sum() if len(results) > 0 else 0
    print(f"  [OK] 筛选后结果: {output_file}")
    print(f"       ({total_valid_runs} runs -> {len(results)} 组)")
    
    # 异常runs结果
    if len(results_filtered) > 0:
        df_filtered = pd.DataFrame(results_filtered)
        output_file_filt = Path(OUTPUT_DIR) / 'ensemble_analysis_filtered.csv'
        df_filtered.to_csv(output_file_filt, index=False, encoding='utf-8-sig')
        total_filtered_runs = df_filtered['n_runs'].sum()
        print(f"  [OK] 异常runs结果: {output_file_filt}")
        print(f"       ({total_filtered_runs} runs -> {len(results_filtered)} 组)")
    
    # 完整数据 (用于 --nofilter) - 使用全部runs的集合平均
    if len(results_all) > 0:
        df_all = pd.DataFrame(results_all)
        output_file_all = Path(OUTPUT_DIR) / 'ensemble_analysis_all.csv'
        df_all.to_csv(output_file_all, index=False, encoding='utf-8-sig')
        total_all_runs = df_all['n_runs'].sum()
        print(f"  [OK] 完整数据(不筛选): {output_file_all}")
        print(f"       ({total_all_runs} runs -> {len(results_all)} 组)")
    
    # 对比结果
    if len(comparison_results) > 0:
        df_comparison = pd.DataFrame(comparison_results)
        comparison_file = Path(OUTPUT_DIR) / 'ensemble_comparison.csv'
        df_comparison.to_csv(comparison_file, index=False, encoding='utf-8-sig')
        print(f"  [OK] 筛选对比: {comparison_file}")
    
    # ========================================================================
    # 阶段7: 生成报告
    # ========================================================================
    print(f"\n[7/7] 生成报告...")
    
    # 统计
    total_runs = len(df_runs)
    n_outliers = len(outlier_records)
    n_improved = sum(1 for r in comparison_results if r.get('improved', False))
    
    # 文本报告
    report_file = Path(OUTPUT_DIR) / 'run_quality_report.txt'
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write("=" * 80 + "\n")
        f.write("GMX MSD综合分析报告 (Step1.2)\n")
        f.write("=" * 80 + "\n\n")
        
        f.write(f"异常检测方法: {OUTLIER_METHOD}\n")
        f.write(f"Intercept阈值: {INTERCEPT_MAX} A^2\n")
        f.write(f"D值上限: {D_MAX_THRESHOLD} cm^2/s\n\n")
        
        f.write(f"总run数: {total_runs}\n")
        f.write(f"异常run数: {n_outliers}\n")
        f.write(f"异常比例: {n_outliers/total_runs*100:.2f}%\n\n")
        
        f.write(f"有效组合: {len(results)}\n")
        f.write(f"异常组合: {len(results_filtered)}\n")
        f.write(f"改进组合: {n_improved}\n\n")
        
        f.write("异常原因统计:\n")
        reason_counts = defaultdict(int)
        for rec in outlier_records:
            reason_counts[rec['reason']] += 1
        for reason, count in sorted(reason_counts.items(), key=lambda x: -x[1]):
            f.write(f"  {reason}: {count}\n")
    
    print(f"  [OK] 质量报告: {report_file}")
    
    # ========================================================================
    # 统计汇总
    # ========================================================================
    print("\n" + "=" * 80)
    print("筛选策略与统计汇总")
    print("=" * 80)
    
    # 筛选策略详细说明
    print(f"\n  [筛选策略] 3级run质量筛选:")
    print(f"    ---------------------------------------------------------------")
    print(f"    | 规则 | 阈值 | 物理含义 | 筛除runs |")
    print(f"    ---------------------------------------------------------------")
    print(f"    | 1. {OUTLIER_METHOD}异常检测 | Q1-1.5*IQR ~ Q3+1.5*IQR | D值统计离群 | {iqr_filtered_count} runs |")
    print(f"    | 2. Intercept筛选 | > {INTERCEPT_MAX} A^2 | 结构松弛异常 | {intercept_filtered_count} runs |")
    print(f"    | 3. D值上限筛选 | >= {D_MAX_THRESHOLD:.0e} cm^2/s | 团簇漂移异常 | {d_filtered_count} runs |")
    print(f"    | 4. 联动筛选 | 同模拟任一元素异常 | 整个模拟标记 | {n_additional} runs |")
    print(f"    ---------------------------------------------------------------")
    print(f"    | 合计 | - | - | {n_outliers} runs ({n_outliers/total_runs*100:.1f}%) |")
    print(f"    ---------------------------------------------------------------")
    
    # 计算涉及的组数
    affected_groups = set()
    for rec in outlier_records:
        affected_groups.add(rec['group_key'])
    
    print(f"\n  [数据流] runs -> 组的转换:")
    print(f"    - 原始: {total_runs} runs -> {len(file_grouped)} 组 (每组6个run)")
    print(f"    - 异常runs: {n_outliers} runs (涉及 {len(affected_groups)} 个组)")
    print(f"    - 有效runs: {total_runs - n_outliers} runs -> {len(results)} 组 (部分组无有效run)")
    print(f"    - 异常runs单独分析: {n_outliers} runs -> {len(results_filtered)} 组")
    
    print(f"\n  [输出文件说明]:")
    print(f"    - large_D_outliers.csv: {n_outliers} 条记录 (供step3绘图/step4分析使用)")
    print(f"    - ensemble_analysis_results.csv: {len(results)} 组 (筛选后,step4.1默认读取)")
    print(f"    - ensemble_analysis_filtered.csv: {len(results_filtered)} 组 (异常组)")
    print(f"    - ensemble_analysis_all.csv: {len(results_all)} 组 (step4.1 --nofilter读取)")
    
    print(f"\n  [筛选效果] 改进统计:")
    print(f"    - 对比组合: {len(comparison_results)} 组")
    print(f"    - R2提升组: {n_improved} 组 ({n_improved/len(comparison_results)*100:.1f}%)" if len(comparison_results) > 0 else "    - 无对比数据")
    
    if len(results) > 0:
        print(f"\n  [数据质量] 有效组统计:")
        print(f"    - R2范围: {df_results['R2_ensemble'].min():.4f} ~ {df_results['R2_ensemble'].max():.4f}")
        print(f"    - R2平均: {df_results['R2_ensemble'].mean():.4f}")
        print(f"    - R2 < 0.5: {len(df_results[df_results['R2_ensemble'] < 0.5])} 组")
    
    print("\n" + "=" * 80)
    print("分析完成!")
    print("=" * 80)

if __name__ == '__main__':
    main()
