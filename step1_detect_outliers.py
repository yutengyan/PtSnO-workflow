#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
=============================================================================
Step1: GMX Run质量分析 - 基于GMX自带D值的异常检测
=============================================================================
创建时间: 2025-10-16 09:00
最后更新: 2025-11-25
作者: GitHub Copilot

核心思想:
    在集合平均前,检查每个run的GMX D值(从.xvg文件注释中提取)
    标记异常run,供后续步骤(step2/3/4)选择性使用,提高数据质量

算法流程:
    1. 扫描 GMX_DATA_DIRS 中的所有 .xvg 文件
    2. 提取GMX D值(从注释行 "# D[element] = XXX")
    3. 对同一(组成-温度-元素)组的多个run计算统计量
    4. 使用多种方法检测异常值(IQR/3σ/MAD)
    5. 额外检查: Intercept过大、关联坏模拟
    6. 生成异常run清单,供后续步骤使用

异常检测方法:
    方法1: IQR原则(推荐,鲁棒性好)
        - 计算四分位数Q1, Q3
        - IQR = Q3 - Q1
        - 标记 D < Q1-1.5×IQR 或 D > Q3+1.5×IQR 的run
    
    方法2: 3σ原则
        - 计算D值的均值μ和标准差σ
        - 标记 |D - μ| > 3σ 的run
    
    方法3: MAD原则(最鲁棒)
        - 计算中位数M和中位数绝对偏差MAD
        - 标记 |D - M| > 3×MAD 的run
    
    方法4: Intercept检查
        - 标记 Intercept > 20.0 A² 的run
    
    方法5: 关联检测
        - 如果某run的多个元素都异常,标记为"坏模拟"

输入:
    - GMX MSD原始数据 (.xvg文件,从GMX_DATA_DIRS读取)

输出: results/
    ├── large_D_outliers.csv            (异常run清单,供step2/3/4使用)
    └── run_quality_analysis_results/   (详细分析结果,可选)
        ├── run_quality_report.txt
        └── ensemble_comparison.csv

使用方法:
    python step1_detect_outliers.py     # 生成异常值清单
    ├── quality_improvement_summary.png (可视化对比)
    └── gmx_D_distribution.png          (GMX D值分布)
=============================================================================
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from pathlib import Path

# 基础目录
BASE_DIR = Path(__file__).parent

# 基础目录
BASE_DIR = Path(__file__).parent

import re
from collections import defaultdict
import sys

# 设置UTF-8输出(Windows兼容)
if sys.platform == 'win32':
    try:
        sys.stdout.reconfigure(encoding='utf-8')
    except:
        # 如果reconfigure失败,使用环境变量
        import os
        os.environ['PYTHONIOENCODING'] = 'utf-8'

# 设置中文显示
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['font.size'] = 10

# ==================== 配置参数 ====================


# GMX数据目录
GMX_DATA_DIRS = [
    # BASE_DIR / 'data' / 'gmx_msd' / 'collected_gmx_msd',
    # BASE_DIR / 'data' / 'gmx_msd' / 'gmx_msd_results_20251015_184626_collected',
    # BASE_DIR / 'data' / 'gmx_msd' / 'unwrap' / 'gmx_msd_results_20251118_152614',  # ✨ 新版 unwrap per-atom
    BASE_DIR / 'data' / 'gmx_msd' / 'unwrap' / 'air' / 'gmx_msd_results_20251124_170114'  # 🌬️ 气象数据
]

# 异常检测参数
OUTLIER_METHOD = 'IQR'  # 可选: '3sigma', 'IQR', 'MAD'
IQR_MULTIPLIER = 1.5    # IQR方法的倍数
SIGMA_MULTIPLIER = 3.0  # 3σ方法的倍数
MAD_MULTIPLIER = 3.0    # MAD方法的倍数

# 最小run数量(少于此数量不做异常检测)
MIN_RUNS_FOR_OUTLIER = 3

# ==================== run级别筛选阈值 (新增) ====================

# Intercept阈值 (A²) - 筛除结构松弛异常的run
INTERCEPT_MAX = 10.0

# D值上限 (cm²/s) - 筛除团簇漂移异常的run
# 物理意义: 纳米团簇固态/液态扩散系数量级通常 1e-5 ~ 1e-2 cm²/s
#          超过0.01视为团簇整体漂移或计算异常
D_MAX_THRESHOLD = 5e-5   # 0.01 cm²/s (原0.1太大) # 但其实GMX 单次run的D值都 < 0.01 cm²/s*

# 拟合参数
FIT_START_RATIO = 0.1
FIT_END_RATIO = 0.9

# 输出目录
OUTPUT_DIR = BASE_DIR / 'results'

# ==================== 核心函数 ====================

def extract_gmx_D_value(filepath):
    """
    从.xvg文件注释中提取GMX计算的D值
    
    查找格式: # D[element] = 0.2585 (+/- 0.7552) (1e-5 cm^2/s)
    
    返回: D值(cm²/s),如果未找到返回None
    """
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            for line in f:
                # 只读注释行(# 或 @开头)
                if not (line.startswith('#') or line.startswith('@')):
                    break
                
                # 查找D值行: # D[        Pt] = 0.2585 (+/- 0.7552) (1e-5 cm^2/s)
                # 注意:D值可能是负数(低温扩散极慢,拟合误差导致)
                if 'D[' in line and '=' in line and 'cm' in line:
                    # 提取D值 (允许元素名前后有空格,允许负号)
                    match = re.search(r'D\[\s*(\w+)\s*\]\s*=\s*([+-]?[\d.]+)', line)
                    if match:
                        element = match.group(1)
                        D_value = float(match.group(2))
                        
                        # 提取单位(默认1e-5 cm²/s)
                        unit_match = re.search(r'\(1e([+-]?\d+)\s*cm', line)
                        if unit_match:
                            exponent = int(unit_match.group(1))
                            D_value *= 10**exponent  # 转换为cm²/s
                        else:
                            # 默认单位: 1e-5 cm²/s
                            D_value *= 1e-5
                        
                        return D_value
        
        return None
    
    except Exception as e:
        print(f"⚠️  提取D值失败: {filepath}")
        print(f"   错误: {e}")
        return None


def read_gmx_msd_xvg(filepath):
    """读取GMX .xvg文件的MSD数据"""
    try:
        data = []
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#') or line.startswith('@'):
                    continue
                if not line:
                    continue
                
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        time = float(parts[0])  # ps
                        msd = float(parts[1])   # nm²
                        data.append([time, msd])
                    except ValueError:
                        continue
        
        if len(data) == 0:
            return None, None
        
        data = np.array(data)
        time = data[:, 0]
        msd = data[:, 1] * 100.0  # nm² → Ų
        
        return time, msd
    
    except Exception as e:
        return None, None


def parse_gmx_filename(filepath):
    """解析文件名提取元数据"""
    path = Path(filepath)
    filename = path.name
    
    # T1000.r24.gpu0_msd_Pt.xvg
    match = re.match(r'T(\d+)\.r(\d+)\.gpu\d+_msd_(\w+)\.xvg', filename)
    
    if not match:
        return None
    
    temp_value = int(match.group(1))
    run_id = int(match.group(2))
    element = match.group(3)
    
    # 提取组成
    parts = path.parts
    composition = None
    if len(parts) >= 3:
        composition = parts[-3]
    
    # 排除分类目录
    exclude_patterns = [
        r'^Pt8$', r'^PtxSn8-x$', r'^o\d+(-\d+)?$',
        r'^\d+K$', r'^run\d+$', r'^more$', r'^o68$'
    ]
    
    if composition:
        for pattern in exclude_patterns:
            if re.match(pattern, composition, re.IGNORECASE):
                for part in reversed(parts[:-2]):
                    is_valid = True
                    for exc_pattern in exclude_patterns:
                        if re.match(exc_pattern, part, re.IGNORECASE):
                            is_valid = False
                            break
                    if is_valid and any(x in part.lower() for x in ['pt', 'sn', 'o']):
                        composition = part
                        break
                break
    
    if not composition:
        return None
    
    temperature = f"{temp_value}K"
    
    return {
        'composition': composition,
        'temperature': temperature,
        'temp_value': temp_value,
        'element': element,
        'run_id': run_id,
        'filepath': str(filepath)
    }


def detect_outliers_iqr(values):
    """
    使用IQR方法检测异常值
    
    返回: (outlier_indices, stats_dict)
    """
    values = np.array(values)
    
    Q1 = np.percentile(values, 25)
    Q3 = np.percentile(values, 75)
    IQR = Q3 - Q1
    
    lower_bound = Q1 - IQR_MULTIPLIER * IQR
    upper_bound = Q3 + IQR_MULTIPLIER * IQR
    
    outlier_mask = (values < lower_bound) | (values > upper_bound)
    outlier_indices = np.where(outlier_mask)[0]
    
    stats = {
        'method': 'IQR',
        'Q1': Q1,
        'Q3': Q3,
        'IQR': IQR,
        'lower_bound': lower_bound,
        'upper_bound': upper_bound,
        'n_outliers': len(outlier_indices),
        'n_total': len(values)
    }
    
    return outlier_indices, stats


def detect_outliers_3sigma(values):
    """使用3σ方法检测异常值"""
    values = np.array(values)
    
    mean = np.mean(values)
    std = np.std(values)
    
    lower_bound = mean - SIGMA_MULTIPLIER * std
    upper_bound = mean + SIGMA_MULTIPLIER * std
    
    outlier_mask = (values < lower_bound) | (values > upper_bound)
    outlier_indices = np.where(outlier_mask)[0]
    
    stats = {
        'method': '3sigma',
        'mean': mean,
        'std': std,
        'lower_bound': lower_bound,
        'upper_bound': upper_bound,
        'n_outliers': len(outlier_indices),
        'n_total': len(values)
    }
    
    return outlier_indices, stats


def detect_outliers_mad(values):
    """使用中位数绝对偏差(MAD)检测异常值"""
    values = np.array(values)
    
    median = np.median(values)
    mad = np.median(np.abs(values - median))
    
    # MAD = 0时特殊处理
    if mad == 0:
        mad = np.std(values) * 0.6745  # 使用标准差估计
    
    lower_bound = median - MAD_MULTIPLIER * mad
    upper_bound = median + MAD_MULTIPLIER * mad
    
    outlier_mask = (values < lower_bound) | (values > upper_bound)
    outlier_indices = np.where(outlier_mask)[0]
    
    stats = {
        'method': 'MAD',
        'median': median,
        'mad': mad,
        'lower_bound': lower_bound,
        'upper_bound': upper_bound,
        'n_outliers': len(outlier_indices),
        'n_total': len(values)
    }
    
    return outlier_indices, stats


def fit_msd_to_diffusion(time, msd):
    """拟合MSD ~ t得到D"""
    start_idx = int(len(time) * FIT_START_RATIO)
    end_idx = int(len(time) * FIT_END_RATIO)
    
    time_fit = time[start_idx:end_idx]
    msd_fit = msd[start_idx:end_idx]
    
    if len(time_fit) < 3:
        return None, None, None
    
    slope, intercept, r_value, p_value, std_err = linregress(time_fit, msd_fit)
    
    # MSD = 6Dt, time单位ps, msd单位Ų
    # D = slope/6, 单位: Ų/ps = 1e-16 cm²/s / 1e-12 s = 1e-4 cm²/s
    D = slope / 6.0 * 1e-4  # cm²/s
    r2 = r_value ** 2
    
    return D, r2, intercept


def ensemble_average_msd(time_arrays, msd_arrays):
    """集合平均MSD"""
    min_length = min(len(t) for t in time_arrays)
    
    time_avg = time_arrays[0][:min_length]
    msd_arrays_trimmed = [msd[:min_length] for msd in msd_arrays]
    msd_avg = np.mean(msd_arrays_trimmed, axis=0)
    msd_std = np.std(msd_arrays_trimmed, axis=0)
    
    return time_avg, msd_avg, msd_std


# ==================== 主分析流程 ====================

def main():
    """主程序"""
    
    print("=" * 80)
    print("GMX Run质量分析 - 基于GMX自带D值的异常检测".center(80))
    print("=" * 80)
    print()
    
    # 创建输出目录
    Path(OUTPUT_DIR).mkdir(exist_ok=True, parents=True)
    
    # ========== [1/6] 扫描所有.xvg文件 ==========
    print("[1/6] 扫描.xvg文件...")
    all_files = []
    for data_dir in GMX_DATA_DIRS:
        if not Path(data_dir).exists():
            print(f"  ⚠️  目录不存在: {data_dir}")
            continue
        
        files = list(Path(data_dir).rglob('*_msd_*.xvg'))
        all_files.extend(files)
        print(f"  ✓ {Path(data_dir).name}: {len(files)}个文件")
    
    print(f"\n  总计: {len(all_files)}个.xvg文件")
    
    # ========== [2/6] 提取GMX D值和run信息 ==========
    print("\n[2/6] 提取GMX D值和run信息...")
    
    run_data = []  # 每个run的完整信息
    
    fail_stats = {'parse': 0, 'gmx_D': 0, 'msd': 0, 'fit': 0}
    
    for i, filepath in enumerate(all_files):
        if (i + 1) % 1000 == 0:
            print(f"  处理中: {i+1}/{len(all_files)}", end='\r')
        
        # 解析文件名
        info = parse_gmx_filename(filepath)
        if not info:
            fail_stats['parse'] += 1
            continue
        
        # 提取GMX D值
        gmx_D = extract_gmx_D_value(filepath)
        if gmx_D is None:
            fail_stats['gmx_D'] += 1
            if i < 10:  # 前10个失败打印详情
                print(f"\n  [DEBUG] GMX D提取失败: {filepath}")
            continue
        
        # 读取MSD数据
        time, msd = read_gmx_msd_xvg(filepath)
        if time is None:
            fail_stats['msd'] += 1
            continue
        
        # 重新拟合得到我们的D值
        our_D, our_r2, intercept = fit_msd_to_diffusion(time, msd)
        if our_D is None:
            fail_stats['fit'] += 1
            continue
        
        run_data.append({
            'composition': info['composition'],
            'temperature': info['temperature'],
            'temp_value': info['temp_value'],
            'element': info['element'],
            'run_id': info['run_id'],
            'filepath': info['filepath'],
            'gmx_D': gmx_D,          # GMX计算的D值
            'our_D': our_D,          # 我们拟合的D值
            'our_r2': our_r2,        # 我们拟合的R²
            'time': time,
            'msd': msd
        })
    
    print(f"\n  ✓ 成功提取: {len(run_data)}个run")
    print(f"  ✓ 失败统计: parse={fail_stats['parse']}, gmx_D={fail_stats['gmx_D']}, "
          f"msd={fail_stats['msd']}, fit={fail_stats['fit']}")
    
    if len(run_data) == 0:
        print("\n❌ 没有提取到任何run数据,请检查:")
        print("  1. .xvg文件是否包含GMX D值注释行")
        print("  2. 文件路径和格式是否正确")
        print("  3. 查看上方[DEBUG]输出")
        return
    
    df_runs = pd.DataFrame(run_data)
    
    # ========== [3/6] 按组成-温度-元素分组 ==========
    print("\n[3/6] 按组成-温度-元素分组...")
    
    df_runs['group_key'] = (df_runs['composition'] + '_' + 
                            df_runs['temperature'] + '_' + 
                            df_runs['element'])
    
    grouped = df_runs.groupby('group_key')
    print(f"  ✓ 共{len(grouped)}个唯一组合")
    
    # ========== [4/6] run级别质量筛选 ==========
    print(f"\n[4/6] run级别质量筛选...")
    print(f"  方法1: {OUTLIER_METHOD}异常检测")
    print(f"  方法2: Intercept≤{INTERCEPT_MAX} A²")
    print(f"  方法3: D<{D_MAX_THRESHOLD} cm²/s")
    
    outlier_records = []
    quality_stats = []
    
    # 统计三种筛选方式
    iqr_filtered_count = 0
    intercept_filtered_count = 0
    d_filtered_count = 0
    total_runs_checked = 0
    
    for group_key, group_df in grouped:
        if len(group_df) < MIN_RUNS_FOR_OUTLIER:
            # run数量太少,不做异常检测
            continue
        
        gmx_D_values = group_df['gmx_D'].values
        
        # 方法1: IQR/3sigma/MAD异常检测
        if OUTLIER_METHOD == 'IQR':
            outlier_indices, stats = detect_outliers_iqr(gmx_D_values)
        elif OUTLIER_METHOD == '3sigma':
            outlier_indices, stats = detect_outliers_3sigma(gmx_D_values)
        elif OUTLIER_METHOD == 'MAD':
            outlier_indices, stats = detect_outliers_mad(gmx_D_values)
        else:
            outlier_indices, stats = detect_outliers_iqr(gmx_D_values)
        
        # 记录IQR检测的异常run
        if len(outlier_indices) > 0:
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
        
        # 方法2+3: 对每个run进行Intercept和D值筛选
        for idx in range(len(group_df)):
            total_runs_checked += 1
            
            # 跳过已经被IQR标记的run
            if idx in outlier_indices:
                continue
            
            row = group_df.iloc[idx]
            filepath = Path(row['filepath'])
            time, msd = read_gmx_msd_xvg(filepath)
            
            if time is None or msd is None:
                continue
            
            # 拟合单个run
            D_run, r2_run, intercept_run = fit_msd_to_diffusion(time, msd)
            
            if D_run is None or r2_run is None or intercept_run is None:
                continue
            
            # Intercept筛选
            if intercept_run > INTERCEPT_MAX:
                outlier_records.append({
                    'group_key': group_key,
                    'composition': row['composition'],
                    'temperature': row['temperature'],
                    'element': row['element'],
                    'run_id': row['run_id'],
                    'gmx_D': row['gmx_D'],
                    'filepath': row['filepath'],
                    'reason': f'Intercept>{INTERCEPT_MAX}A²'
                })
                intercept_filtered_count += 1
                continue
            
            # D值上限筛选
            if D_run >= D_MAX_THRESHOLD:
                outlier_records.append({
                    'group_key': group_key,
                    'composition': row['composition'],
                    'temperature': row['temperature'],
                    'element': row['element'],
                    'run_id': row['run_id'],
                    'gmx_D': row['gmx_D'],
                    'filepath': row['filepath'],
                    'reason': f'D≥{D_MAX_THRESHOLD}cm²/s'
                })
                d_filtered_count += 1
        
        # 记录统计信息
        quality_stats.append({
            'group_key': group_key,
            'n_runs': len(group_df),
            'n_outliers': len(outlier_indices),
            **stats
        })
    
    print(f"\n  ✓ 初步筛除: {len(outlier_records)} run / {total_runs_checked} run ({len(outlier_records)/total_runs_checked*100:.1f}%)")
    print(f"    ├─ {OUTLIER_METHOD}异常检测: {iqr_filtered_count} run")
    print(f"    ├─ Intercept>{INTERCEPT_MAX}A²: {intercept_filtered_count} run")
    print(f"    └─ D≥{D_MAX_THRESHOLD}cm²/s: {d_filtered_count} run")
    print(f"  ✓ 涉及 {len(quality_stats)} 个组合")
    
    # ========== [4.5/6] 同模拟run全元素联动筛选 ==========
    print(f"\n[4.5/6] 同模拟run全元素联动筛选...")
    print(f"  规则: 若某模拟的任意元素被标记为outlier,则该模拟的所有元素均标记为坏结构")
    
    # 构建已标记的(composition, temperature, run_id)集合
    bad_simulation_set = set()
    for rec in outlier_records:
        bad_simulation_set.add((rec['composition'], rec['temperature'], rec['run_id']))
    
    print(f"  ✓ 初步标记坏模拟: {len(bad_simulation_set)} 个 (composition, temperature, run_id)")
    
    # 对所有数据进行联动筛选,找出同一模拟的其他元素
    additional_outliers = []
    for _, row in df_runs.iterrows():
        sim_key = (row['composition'], row['temperature'], row['run_id'])
        
        # 如果该模拟已被标记为坏,但该元素尚未被标记
        if sim_key in bad_simulation_set:
            # 检查是否已经在outlier_records中
            already_marked = False
            for rec in outlier_records:
                if (rec['composition'] == row['composition'] and 
                    rec['temperature'] == row['temperature'] and 
                    rec['run_id'] == row['run_id'] and 
                    rec['element'] == row['element']):
                    already_marked = True
                    break
            
            if not already_marked:
                additional_outliers.append({
                    'group_key': f"{row['composition']}_{row['temperature']}_{row['element']}",
                    'composition': row['composition'],
                    'temperature': row['temperature'],
                    'element': row['element'],
                    'run_id': row['run_id'],
                    'gmx_D': row['gmx_D'],
                    'filepath': row['filepath'],
                    'reason': 'linked_bad_simulation'  # 联动标记
                })
    
    # 合并额外标记的outliers
    n_additional = len(additional_outliers)
    outlier_records.extend(additional_outliers)
    
    print(f"  ✓ 联动标记额外 {n_additional} 个元素run")
    print(f"  ✓ 最终筛除: {len(outlier_records)} run (初步{len(outlier_records)-n_additional} + 联动{n_additional})")
    
    # 保存异常run清单
    if len(outlier_records) > 0:
        df_outliers = pd.DataFrame(outlier_records)
        outlier_file = Path(OUTPUT_DIR) / 'large_D_outliers.csv'
        df_outliers.to_csv(outlier_file, index=False)
        print(f"  ✓ 保存异常run清单: {outlier_file}")
    
    # ========== [5/6] 对比集合平均(舍弃前 vs 舍弃后) ==========
    print("\n[5/6] 对比集合平均效果...")
    
    comparison_results = []
    
    for group_key, group_df in grouped:
        if len(group_df) < 2:
            continue
        
        composition = group_df.iloc[0]['composition']
        temperature = group_df.iloc[0]['temperature']
        element = group_df.iloc[0]['element']
        
        # ===== 方法1: 全部run集合平均 =====
        time_all = [row['time'] for _, row in group_df.iterrows()]
        msd_all = [row['msd'] for _, row in group_df.iterrows()]
        
        time_avg_all, msd_avg_all, msd_std_all = ensemble_average_msd(time_all, msd_all)
        D_all, r2_all, _ = fit_msd_to_diffusion(time_avg_all, msd_avg_all)
        
        # ===== 方法2: 舍弃异常run后集合平均 =====
        outlier_run_ids = set()
        for rec in outlier_records:
            if rec['group_key'] == group_key:
                outlier_run_ids.add(rec['run_id'])
        
        if len(outlier_run_ids) > 0 and len(group_df) - len(outlier_run_ids) >= 2:
            # 有异常且剩余run≥2
            clean_df = group_df[~group_df['run_id'].isin(outlier_run_ids)]
            
            if len(clean_df) >= 2:  # 确保有足够run
                time_clean = [row['time'] for _, row in clean_df.iterrows()]
                msd_clean = [row['msd'] for _, row in clean_df.iterrows()]
                
                if len(time_clean) > 0 and len(msd_clean) > 0:  # 再次检查
                    time_avg_clean, msd_avg_clean, msd_std_clean = ensemble_average_msd(time_clean, msd_clean)
                    D_clean, r2_clean, _ = fit_msd_to_diffusion(time_avg_clean, msd_avg_clean)
                else:
                    D_clean = D_all
                    r2_clean = r2_all
            else:
                D_clean = D_all
                r2_clean = r2_all
        else:
            # 无异常或舍弃后不足2个
            D_clean = D_all
            r2_clean = r2_all
        
        # 计算改进
        if D_all is not None and D_clean is not None:
            delta_r2 = r2_clean - r2_all
            delta_D = D_clean - D_all
            delta_D_pct = (delta_D / D_all * 100) if D_all != 0 else 0
            
            comparison_results.append({
                'composition': composition,
                'temperature': temperature,
                'element': element,
                'n_runs_total': len(group_df),
                'n_runs_outlier': len(outlier_run_ids),
                'n_runs_clean': len(group_df) - len(outlier_run_ids),
                'D_all': D_all,
                'r2_all': r2_all,
                'D_clean': D_clean,
                'r2_clean': r2_clean,
                'delta_r2': delta_r2,
                'delta_D': delta_D,
                'delta_D_pct': delta_D_pct,
                'improved': delta_r2 > 0
            })
    
    df_comparison = pd.DataFrame(comparison_results)
    comparison_file = Path(OUTPUT_DIR) / 'ensemble_comparison.csv'
    df_comparison.to_csv(comparison_file, index=False)
    
    print(f"  ✓ 对比结果: {len(comparison_results)}个组合")
    print(f"  ✓ 保存至: {comparison_file}")
    
    # ========== [6/6] 生成报告和可视化 ==========
    print("\n[6/6] 生成报告和可视化...")
    
    # 统计汇总
    n_improved = df_comparison['improved'].sum() if len(df_comparison) > 0 else 0
    avg_delta_r2 = df_comparison['delta_r2'].mean() if len(df_comparison) > 0 else 0
    
    # 生成文本报告
    report_file = Path(OUTPUT_DIR) / 'run_quality_report.txt'
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write("=" * 80 + "\n")
        f.write("GMX Run质量分析报告\n")
        f.write("=" * 80 + "\n\n")
        
        f.write(f"异常检测方法: {OUTLIER_METHOD}\n")
        f.write(f"最小run数量: {MIN_RUNS_FOR_OUTLIER}\n\n")
        
        f.write(f"总run数: {len(df_runs)}\n")
        f.write(f"异常run数: {len(outlier_records)}\n")
        f.write(f"异常比例: {len(outlier_records)/len(df_runs)*100:.2f}%\n\n")
        
        f.write(f"对比组合数: {len(df_comparison)}\n")
        f.write(f"改进组合数: {n_improved}\n")
        f.write(f"平均ΔR²: {avg_delta_r2:.6f}\n\n")
        
        if len(df_comparison) > 0:
            f.write("改进最显著的前10个组合:\n")
            top10 = df_comparison.nlargest(10, 'delta_r2')
            for _, row in top10.iterrows():
                f.write(f"  {row['composition']}_{row['temperature']}_{row['element']}: "
                       f"ΔR²={row['delta_r2']:.6f}, "
                       f"舍弃{row['n_runs_outlier']}/{row['n_runs_total']}个run\n")
    
    print(f"  ✓ 文本报告: {report_file}")
    
    # 可视化: GMX D值分布
    if len(df_runs) > 0:
        fig, axes = plt.subplots(1, 3, figsize=(15, 4))
        
        for i, element in enumerate(['Pt', 'Sn', 'PtSn']):
            df_elem = df_runs[df_runs['element'] == element]
            if len(df_elem) > 0:
                axes[i].hist(df_elem['gmx_D'] * 1e4, bins=50, alpha=0.7, edgecolor='black')
                axes[i].set_xlabel('D (GMX) [10⁻⁴ cm²/s]')
                axes[i].set_ylabel('频数')
                axes[i].set_title(f'{element} (n={len(df_elem)})')
                axes[i].grid(alpha=0.3)
        
        plt.tight_layout()
        fig_file = Path(OUTPUT_DIR) / 'gmx_D_distribution.png'
        plt.savefig(fig_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"  ✓ D值分布图: {fig_file}")
    
    # 可视化: 改进效果
    if len(df_comparison) > 0:
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # (1) R²改进散点图
        axes[0, 0].scatter(df_comparison['r2_all'], df_comparison['r2_clean'], 
                          alpha=0.5, s=20)
        axes[0, 0].plot([0, 1], [0, 1], 'r--', label='y=x')
        axes[0, 0].set_xlabel('R² (全部run)')
        axes[0, 0].set_ylabel('R² (舍弃异常后)')
        axes[0, 0].set_title('R²改进效果')
        axes[0, 0].legend()
        axes[0, 0].grid(alpha=0.3)
        
        # (2) ΔR²直方图
        axes[0, 1].hist(df_comparison['delta_r2'], bins=50, alpha=0.7, edgecolor='black')
        axes[0, 1].axvline(0, color='red', linestyle='--', label='ΔR²=0')
        axes[0, 1].set_xlabel('ΔR² (舍弃后 - 全部)')
        axes[0, 1].set_ylabel('频数')
        axes[0, 1].set_title(f'R²改进分布 (改进率: {n_improved/len(df_comparison)*100:.1f}%)')
        axes[0, 1].legend()
        axes[0, 1].grid(alpha=0.3)
        
        # (3) D值改进散点图
        axes[1, 0].scatter(df_comparison['D_all'] * 1e4, 
                          df_comparison['D_clean'] * 1e4, 
                          alpha=0.5, s=20)
        max_D = max(df_comparison['D_all'].max(), df_comparison['D_clean'].max()) * 1e4
        axes[1, 0].plot([0, max_D], [0, max_D], 'r--', label='y=x')
        axes[1, 0].set_xlabel('D (全部run) [10⁻⁴ cm²/s]')
        axes[1, 0].set_ylabel('D (舍弃异常后) [10⁻⁴ cm²/s]')
        axes[1, 0].set_title('D值改进效果')
        axes[1, 0].legend()
        axes[1, 0].grid(alpha=0.3)
        
        # (4) 异常run数量分布
        axes[1, 1].hist(df_comparison['n_runs_outlier'], 
                       bins=range(0, df_comparison['n_runs_outlier'].max()+2), 
                       alpha=0.7, edgecolor='black')
        axes[1, 1].set_xlabel('异常run数量')
        axes[1, 1].set_ylabel('频数')
        axes[1, 1].set_title('异常run分布')
        axes[1, 1].grid(alpha=0.3)
        
        plt.tight_layout()
        fig_file = Path(OUTPUT_DIR) / 'quality_improvement_summary.png'
        plt.savefig(fig_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"  ✓ 改进效果图: {fig_file}")
    
    print("\n" + "=" * 80)
    print("✅ 分析完成!".center(80))
    print("=" * 80)
    print(f"\n输出目录: {OUTPUT_DIR}")
    print("\n关键发现:")
    print(f"  • 异常run: {len(outlier_records)}/{len(df_runs)} ({len(outlier_records)/len(df_runs)*100:.2f}%)")
    print(f"  • 改进组合: {n_improved}/{len(df_comparison)} ({n_improved/len(df_comparison)*100:.1f}%)")
    print(f"  • 平均ΔR²: {avg_delta_r2:.6f}")
    print("\n推荐查看:")
    print("  1. large_D_outliers.csv - 异常run清单")
    print("  2. ensemble_comparison.csv - 改进效果对比")
    print("  3. quality_improvement_summary.png - 可视化总结")
    print()


if __name__ == '__main__':
    main()

