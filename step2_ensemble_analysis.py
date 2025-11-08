#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
=============================================================================
GMX MSD集合平均分析工具 (v3.0 简化版)
=============================================================================
创建时间: 2025-10-16
作者: GitHub Copilot

工作流程:
    1. [可选] 运行 step1_detect_outliers.py 生成异常run清单
    2. 运行本脚本进行集合平均和筛选
    3. 运行 step3_plot_msd.py 绘制MSD曲线

筛选策略 (v3.0):
    第0步: 大D值预筛选 (单次run级别,可选)
        - 如果存在 large_D_outliers.csv,则舍弃异常run
        - 如果不存在,跳过此步骤
    
    第1步: Intercept筛选
        - 阈值: INTERCEPT_MAX (默认20.0 A²,可调整)
        - 物理意义: 过大表示结构松弛异常
    
    第2步: D值物理合理性筛选
        - 阈值: D < D_MAX_THRESHOLD (默认0.1 cm²/s)
        - 去除团簇整体漂移

参数说明:
    - INTERCEPT_MAX: Intercept上限,默认20.0 A²
      * 15.0: 严格模式,适合寻找高质量数据
      * 20.0: 宽松模式,适合MSD曲线绘图
      * 可根据需求调整
    
    - D_MAX_THRESHOLD: D值上限,默认0.1 cm²/s
      * 筛除异常大的扩散系数(团簇漂移)

运行方法:
    conda activate base
    cd v3_simplified_workflow
    python scripts/step2_ensemble_analysis.py

输出:
    results/
    ├── ensemble_analysis_results.csv  (筛选后数据)
    ├── statistics_by_element.csv      (元素统计)
    └── filtering_report.txt           (筛选报告)
=============================================================================
"""

import pandas as pd
import numpy as np
from scipy.stats import linregress
from pathlib import Path
import os
import re
from collections import defaultdict

# ==================== 可调参数 ====================

# 基础目录 (自动检测脚本所在的v3_simplified_workflow/)
BASE_DIR = Path(__file__).parent

# 数据源目录
GMX_DATA_DIRS = [
    BASE_DIR / 'data' / 'gmx_msd' / 'collected_gmx_msd',
    BASE_DIR / 'data' / 'gmx_msd' / 'gmx_msd_results_20251015_184626_collected'
]

# 输出目录
OUTPUT_DIR = BASE_DIR / 'results'

# 异常run清单 (可选,如果不存在则跳过预筛选)
OUTLIERS_FILE = BASE_DIR / 'results' / 'large_D_outliers.csv'

# ==================== 筛选阈值 (可调整) ====================

# Intercept阈值 (A²)
# 建议值:
#   15.0 - 严格模式 (高质量数据)
#   20.0 - 宽松模式 (适合MSD绘图,默认)
#   25.0 - 极宽松 (几乎不筛除)
INTERCEPT_MAX = 20.0

# D值上限 (cm²/s)
# 超过此值认为是团簇整体漂移
D_MAX_THRESHOLD = 0.1

# 拟合参数
FIT_START_RATIO = 0.1  # 拟合起始比例
FIT_END_RATIO = 0.9    # 拟合结束比例

# ==================== 核心函数 ====================

def load_large_D_outliers():
    """
    加载大D值异常run清单 (可选)
    
    返回:
        set: 异常文件路径集合,如果文件不存在则返回空集
    """
    outlier_file = Path(OUTLIERS_FILE)
    try:
        df_outliers = pd.read_csv(outlier_file)
        outlier_files = set(df_outliers['filepath'].values)
        print(f"   [OK] 加载异常run清单: {len(outlier_files)} 个")
        return outlier_files
    except FileNotFoundError:
        print(f"   [!] 未找到异常清单: {outlier_file}")
        print(f"   [!] 将不进行预筛选,使用所有run")
        return set()


def parse_gmx_filename(filepath):
    """
    解析GMX MSD文件名
    
    路径示例:
        .../g-1-O1Sn4Pt3/1000K/T1000.r24.gpu0_msd_Pt.xvg
        .../pt8sn2-1-best/300K/T300.r0.gpu0_msd_Sn.xvg
    
    返回:
        dict: {'composition', 'temperature', 'element', 'filepath'}
    """
    filename = filepath.stem
    
    # 提取元素
    if filename.endswith('_Pt'):
        element = 'Pt'
    elif filename.endswith('_Sn'):
        element = 'Sn'
    elif filename.endswith('_PtSn'):
        element = 'PtSn'
    else:
        return None
    
    # 提取温度 (从文件名)
    temp_match = re.match(r'T(\d+)', filename)
    if not temp_match:
        return None
    temperature = f"{temp_match.group(1)}K"
    
    # 提取组成 (从路径,温度目录的父目录)
    parts = filepath.parts
    if len(parts) < 3:
        return None
    
    # 找温度目录
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
        'element': element,
        'filepath': filepath,
        'full_key': f"{composition}_{temperature}"
    }


def read_xvg(filepath):
    """
    读取GMX .xvg文件
    
    注意单位转换:
        - GMX输出: MSD (nm²), Time (ps)
        - 本函数输出: MSD (A²), Time (ps)
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
                        t = float(parts[0])           # ps
                        msd = float(parts[1]) * 100   # nm² → A² (×100)
                        data.append([t, msd])
                    except ValueError:
                        continue
        
        if len(data) == 0:
            return None, None
        
        data = np.array(data)
        return data[:, 0], data[:, 1]
    
    except Exception as e:
        return None, None


def fit_msd_to_diffusion(time, msd, fit_start_ratio=0.1, fit_end_ratio=0.9):
    """
    线性拟合MSD → 扩散系数
    
    输入:
        time: ps (皮秒)
        msd: A² (埃²)
    
    输出:
        D: cm²/s
        r2: 拟合优度
        intercept: A²
    
    公式推导:
        MSD (A²) = 6D (cm²/s) × t (ps) × (10^-16 cm²/A²) / (10^-12 s/ps)
        MSD = 6D × t × 10^-4
        slope (A²/ps) = 6D × 10^-4
        D (cm²/s) = slope / (6 × 10^-4)
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
        # ✅ 修复单位转换: slope (A²/ps) → D (cm²/s)
        D = slope / 6.0 * 1e-4  # 原来缺少 ×10^-4
        r2 = r_value ** 2
        return D, r2, intercept
    except Exception as e:
        return None, None, None


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


def main():
    """主分析流程"""
    
    print("\n" + "="*80)
    print("GMX MSD集合平均分析工具 (v3.0 简化版)")
    print("="*80)
    
    # 创建输出目录
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # 步骤0: 加载异常清单(可选)
    print("\n[0/5] 加载大D值异常清单...")
    large_D_outliers = load_large_D_outliers()
    
    # 步骤1: 收集文件
    print("\n[1/5] 收集GMX MSD文件...")
    all_files = []
    for data_dir in GMX_DATA_DIRS:
        gmx_dir = Path(data_dir)
        if not gmx_dir.exists():
            print(f"  ⚠️ 目录不存在: {data_dir}")
            continue
        
        files_in_dir = list(gmx_dir.rglob('*_msd_*.xvg'))
        print(f"  {gmx_dir.name}: {len(files_in_dir)} 文件")
        all_files.extend(files_in_dir)
    
    print(f"\n   总计: {len(all_files)} 文件")
    
    # 步骤2: 解析和分组
    print("\n[2/5] 解析文件信息...")
    file_info = []
    for filepath in all_files:
        info = parse_gmx_filename(filepath)
        if info:
            file_info.append(info)
    
    print(f"   成功解析: {len(file_info)} 文件")
    
    # 按体系-温度-元素分组
    print("\n[3/5] 分组数据...")
    grouped = defaultdict(list)
    for info in file_info:
        key = f"{info['full_key']}_{info['element']}"
        grouped[key].append(info['filepath'])
    
    print(f"   独立组合: {len(grouped)} 组 (每组包含多次模拟run)")
    
    # 统计run数分布
    run_counts = [len(files) for files in grouped.values()]
    avg_runs = sum(run_counts) / len(run_counts)
    min_runs = min(run_counts)
    max_runs = max(run_counts)
    print(f"   每组run数: 最少{min_runs}, 最多{max_runs}, 平均{avg_runs:.1f}")
    
    # 预筛选统计
    total_files_before = sum(len(files) for files in grouped.values())
    files_to_remove = 0
    for files in grouped.values():
        for filepath in files:
            if str(filepath) in large_D_outliers:
                files_to_remove += 1
    
    print(f"\n   预筛选(step1异常检测):")
    print(f"     总run数: {total_files_before}")
    print(f"     已舍弃异常run: {files_to_remove} ({files_to_remove/total_files_before*100:.1f}%)")
    print(f"     保留有效run: {total_files_before - files_to_remove}")
    
    # 步骤4: 集合平均分析 (使用step1筛选后的run)
    print(f"\n[4/5] 对 {len(grouped)} 组进行集合平均分析...")
    print(f"   方法: MSD曲线平均 → 拟合 → 求D值")
    print(f"   (非直接平均D值,而是先平均MSD时间序列再拟合)")
    results = []
    results_filtered = []  # 新增: 保存被筛选掉的runs的集合平均
    processed = 0
    processed_filtered = 0
    step0_filtered = 0  # step1舍弃的run数
    total_runs_before_filter = 0
    
    for key, files in grouped.items():
        # 解析key
        parts = key.rsplit('_', 1)
        comp_temp = parts[0]
        element = parts[1]
        
        comp_temp_parts = comp_temp.rsplit('_', 1)
        composition = comp_temp_parts[0]
        temperature = comp_temp_parts[1]
        temp_value = int(temperature.replace('K', ''))
        
        # 统计原始run数
        original_n_files = len(files)
        total_runs_before_filter += original_n_files
        
        # 使用step1筛选后的run (舍弃large_D_outliers中的run)
        valid_files = [f for f in files if str(f) not in large_D_outliers]
        filtered_files = [f for f in files if str(f) in large_D_outliers]  # 新增: 被筛选掉的runs
        n_filtered = original_n_files - len(valid_files)
        step0_filtered += n_filtered
        
        # === 处理有效runs (原有逻辑) ===
        if len(valid_files) > 0:
            # 集合平均(只用通过筛选的run)
            avg_result = ensemble_average_msd(valid_files)
            if avg_result is not None:
                # 拟合集合平均的MSD
                D, r2, intercept = fit_msd_to_diffusion(
                    avg_result['time'],
                    avg_result['msd_mean'],
                    FIT_START_RATIO,
                    FIT_END_RATIO
                )
                
                if D is not None and r2 is not None:
                    processed += 1
                    results.append({
                        'composition': composition,
                        'temperature': temperature,
                        'temp_value': temp_value,
                        'element': element,
                        'n_runs': len(valid_files),
                        'n_runs_original': original_n_files,
                        'n_runs_filtered': n_filtered,
                        'D_ensemble': D,
                        'R2_ensemble': r2,
                        'intercept': intercept
                    })
        
        # === 新增: 处理被筛选掉的runs ===
        if len(filtered_files) > 0:
            avg_result_filt = ensemble_average_msd(filtered_files)
            if avg_result_filt is not None:
                D_filt, r2_filt, intercept_filt = fit_msd_to_diffusion(
                    avg_result_filt['time'],
                    avg_result_filt['msd_mean'],
                    FIT_START_RATIO,
                    FIT_END_RATIO
                )
                
                if D_filt is not None and r2_filt is not None:
                    processed_filtered += 1
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
        
        if (processed + processed_filtered) % 100 == 0:
            print(f"   [+] 已处理 {processed} 个有效组 + {processed_filtered} 个异常组...")
    
    # 步骤5: 保存结果
    print(f"\n[5/5] 保存结果...")
    
    # 保存有效runs的集合平均
    df_results = pd.DataFrame(results)
    output_file = Path(OUTPUT_DIR) / 'ensemble_analysis_results.csv'
    df_results.to_csv(output_file, index=False, encoding='utf-8-sig')
    print(f"   [OK] 有效runs结果: {output_file}")
    
    # 保存被筛选掉runs的集合平均
    if len(results_filtered) > 0:
        df_filtered = pd.DataFrame(results_filtered)
        output_file_filt = Path(OUTPUT_DIR) / 'ensemble_analysis_filtered.csv'
        df_filtered.to_csv(output_file_filt, index=False, encoding='utf-8-sig')
        print(f"   [OK] 被筛选runs结果: {output_file_filt}")
    else:
        print(f"   ⚠️ 无被筛选掉的runs可做集合平均")
    
    # 统计报告
    print("\n" + "="*80)
    print("集合平均统计")
    print("="*80)
    print(f"  原始数据:")
    print(f"    • 总run数: {total_runs_before_filter} 个run")
    print(f"    • 独立组合: {len(grouped)} 组 (体系×温度×元素)")
    print(f"    • 平均每组: {total_runs_before_filter/len(grouped):.1f} run")
    print(f"\n  Step1质量筛选:")
    print(f"    • 舍弃run: {step0_filtered} run ({step0_filtered/total_runs_before_filter*100:.1f}%)")
    print(f"    • 保留run: {total_runs_before_filter - step0_filtered} run")
    print(f"\n  集合平均结果:")
    print(f"    • 方法: MSD曲线平均 → 线性拟合 → D值")
    print(f"    • 有效组合: {len(results)} 组")
    print(f"    • 异常组合: {len(results_filtered)} 组 (被step1筛选掉的runs)")
    print(f"    • Step2筛选: 无 (保留所有可处理组合)")
    print(f"    • 组合可用率: {len(results)/len(grouped)*100:.1f}% ({len(results)}/{len(grouped)})")

    
    # 数据统计
    if len(results) > 0:
        print(f"\n数据统计:")
        print(f"  R² 范围: {df_results['R2_ensemble'].min():.4f} ~ {df_results['R2_ensemble'].max():.4f}")
        print(f"  R² 平均: {df_results['R2_ensemble'].mean():.4f}")
        print(f"  R² < 0.5: {len(df_results[df_results['R2_ensemble'] < 0.5])} 组 ({len(df_results[df_results['R2_ensemble'] < 0.5])/len(df_results)*100:.1f}%) - 拟合质量差")
        print(f"  Intercept 范围: {df_results['intercept'].min():.2f} ~ {df_results['intercept'].max():.2f} A²")
        print(f"  Intercept >10: {len(df_results[df_results['intercept'] > 10])}")
        print(f"  Intercept >15: {len(df_results[df_results['intercept'] > 15])}")
        print(f"  n_runs = 1: {len(df_results[df_results['n_runs'] == 1])} 组 ({len(df_results[df_results['n_runs'] == 1])/len(df_results)*100:.1f}%) - 无法集合平均")
        print(f"  n_runs ≥ 3: {len(df_results[df_results['n_runs'] >= 3])} 组 ({len(df_results[df_results['n_runs'] >= 3])/len(df_results)*100:.1f}%)")

    
    print("\n" + "="*80)
    print("分析完成!")
    print("="*80)


if __name__ == '__main__':
    main()
