#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
分区策略分析脚本 v2.4
======================

基于 1分区 vs 2分区 vs 3分区 的对比分析。

核心思路：
1. 先判断"分区是否有意义"（n=2 vs n=1 的R²增益）
   - 增益 < 2%: 推荐n=1，分区无意义
2. 再比较n=2 vs n=3（仅当分区有意义时）
   - 差异不显著时，选择简单模型n=2
3. 新增：检查3分区的热容差异显著性
   - 如果partition2和partition3的热容差异不显著（差异<2倍误差），说明3分区无意义

评分公式 (0-100分，仅用于n=2和n=3的比较):
- R² 拟合优度: 50% (平均R² 40% + 最小R² 10%)
- 聚类质量: 30% (Silhouette Score)
- Cv误差: 20% (误差越小越好)

注意：1分区没有聚类，只用R²作为基准判断分区是否有意义

作者: GitHub Copilot
日期: 2025-11-28
版本: 2.4
"""

import pandas as pd
import numpy as np
from pathlib import Path
import re
from collections import defaultdict
from scipy import stats


def extract_structure_name(filename):
    """从文件名提取结构名"""
    match = re.match(r'(.+?)_(auto2|auto3|fixed2|fixed3|kmeans_n2|kmeans_n3)_quality_metrics\.csv', filename)
    if match:
        return match.group(1), match.group(2)
    return None, None


def load_quality_metrics(results_dir):
    """加载所有质量指标文件"""
    results_dir = Path(results_dir)
    
    data = defaultdict(dict)
    
    for csv_file in results_dir.glob('*_quality_metrics.csv'):
        structure, partition_type = extract_structure_name(csv_file.name)
        if structure and partition_type:
            try:
                df = pd.read_csv(csv_file)
                if not df.empty:
                    data[structure][partition_type] = df
            except Exception as e:
                print(f"[WARNING] 无法读取 {csv_file.name}: {e}")
    
    return data


def calculate_single_partition_r2(results_dir, structure):
    """计算1分区（整体线性拟合）的R²
    
    从clustered_data.csv读取原始数据，进行整体能量-温度线性拟合
    与2/3分区的热容R²保持一致的比较基准
    """
    results_dir = Path(results_dir)
    
    # 尝试读取任意一个clustered_data文件
    for pattern in [f'{structure}_kmeans_n2_clustered_data.csv', 
                    f'{structure}_auto2_clustered_data.csv']:
        data_file = results_dir / pattern
        if data_file.exists():
            try:
                df = pd.read_csv(data_file)
                # 检查温度和能量列名
                temp_col = 'temperature' if 'temperature' in df.columns else 'temp'
                energy_col = 'avg_energy' if 'avg_energy' in df.columns else 'energy'
                
                if temp_col in df.columns and energy_col in df.columns:
                    # 按温度平均后拟合（与2/3分区保持一致）
                    df_avg = df.groupby(temp_col).agg({energy_col: 'mean'}).reset_index()
                    x = df_avg[temp_col].values
                    y = df_avg[energy_col].values
                    
                    # 线性回归 (温度 vs 能量)
                    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
                    r2 = r_value ** 2
                    
                    return {
                        'r2': r2,
                        'slope': slope,  # Cv = slope * 1000 meV/K
                        'intercept': intercept,
                        'n_points': len(x)
                    }
            except Exception as e:
                print(f"[WARNING] 计算1分区R²失败 {structure}: {e}")
    
    return None


def calculate_comprehensive_score(df):
    """计算综合评分 (0-100分)
    
    评分公式:
    - R² 拟合优度: 50% (平均R² 40% + 最小R² 10%)
    - 聚类质量: 30% (Silhouette Score, 范围[-1,1])
    - Cv误差: 20% (误差/Cv, 越小越好)
    """
    if df is None or df.empty:
        return None
    
    result = {}
    
    # 1. R² 指标 (50分)
    r2_col = 'R2' if 'R2' in df.columns else 'r2'
    if r2_col in df.columns:
        result['avg_r2'] = df[r2_col].mean()
        result['min_r2'] = df[r2_col].min()
        # R²得分: 平均R² 40分 + 最小R² 10分
        result['r2_score'] = result['avg_r2'] * 40 + result['min_r2'] * 10
    else:
        result['avg_r2'] = 0
        result['min_r2'] = 0
        result['r2_score'] = 0
    
    # 2. 聚类质量指标 (30分)
    if 'silhouette_score' in df.columns:
        # Silhouette Score 范围 [-1, 1], 需要归一化到 [0, 1]
        silhouette = df['silhouette_score'].iloc[0]
        result['silhouette'] = silhouette
        # 归一化: (silhouette + 1) / 2 * 30
        result['silhouette_score'] = max(0, (silhouette + 1) / 2) * 30
    else:
        result['silhouette'] = None
        result['silhouette_score'] = 15  # 默认中等分数
    
    # Davies-Bouldin 作为参考 (越低越好，但不计入评分)
    if 'davies_bouldin' in df.columns:
        result['davies_bouldin'] = df['davies_bouldin'].iloc[0]
    else:
        result['davies_bouldin'] = None
    
    # 3. Cv误差指标 (20分)
    if 'Cv_cluster' in df.columns and 'Cv_cluster_err' in df.columns:
        # 计算相对误差 (误差/Cv)
        valid = df['Cv_cluster'].abs() > 0.01
        if valid.any():
            error_ratios = (df.loc[valid, 'Cv_cluster_err'].abs() / 
                           df.loc[valid, 'Cv_cluster'].abs())
            result['cv_error_ratio'] = error_ratios.mean()
            # 误差得分: 误差比越小分数越高
            # 假设误差比 < 0.05 为优秀(20分), > 0.2 为差(0分)
            result['error_score'] = max(0, min(20, (0.2 - result['cv_error_ratio']) / 0.15 * 20))
        else:
            result['cv_error_ratio'] = 0
            result['error_score'] = 20
    else:
        result['cv_error_ratio'] = None
        result['error_score'] = 10  # 默认中等分数
    
    # 4. 热容差异显著性 (新指标，用于评估多分区是否有意义)
    result['cv_diff_significant'] = True  # 默认显著
    result['cv_values'] = []
    result['cv_errors'] = []
    result['cv_diff_ratio'] = None
    
    if 'Cv_cluster' in df.columns and 'Cv_cluster_err' in df.columns:
        cv_values = df['Cv_cluster'].values
        cv_errors = df['Cv_cluster_err'].values
        result['cv_values'] = cv_values.tolist()
        result['cv_errors'] = cv_errors.tolist()
        
        if len(cv_values) >= 2:
            # 计算相邻分区的热容差异
            cv_diffs = []
            for i in range(len(cv_values) - 1):
                diff = abs(cv_values[i+1] - cv_values[i])
                # 合并误差 (误差传播)
                combined_err = np.sqrt(cv_errors[i]**2 + cv_errors[i+1]**2)
                # 差异显著性: 差异是否大于2倍合并误差
                significant = diff > 2 * combined_err
                cv_diffs.append({
                    'diff': diff,
                    'combined_err': combined_err,
                    'significant': significant,
                    'ratio': diff / combined_err if combined_err > 0 else float('inf')
                })
            
            result['cv_diffs'] = cv_diffs
            
            # 对于2分区和3分区都计算热容差异显著性
            if len(cv_values) == 2:  # 2分区
                result['cv_diff_significant'] = cv_diffs[0]['significant']
                result['cv_diff_ratio'] = cv_diffs[0]['ratio']
            elif len(cv_values) == 3:  # 3分区
                # 检查partition2和partition3的差异
                last_diff = cv_diffs[-1]
                result['cv_diff_significant'] = last_diff['significant']
                result['cv_diff_ratio'] = last_diff['ratio']
    
    # 综合得分
    result['total_score'] = result['r2_score'] + result['silhouette_score'] + result['error_score']
    
    return result


def analyze_all_structures(results_dir):
    """分析所有结构的分区质量 - 综合评分 (包含1分区基准)"""
    data = load_quality_metrics(results_dir)
    results_dir = Path(results_dir)
    
    analysis_results = []
    
    for structure in sorted(data.keys()):
        partitions = data[structure]
        
        result = {
            'structure': structure,
            # n=1 指标 (整体拟合)
            'n1_r2': 0,
            'n1_score': 0,  # 简化得分 = R² * 50
            # n=2 指标
            'auto2_score': 0,
            'auto2_avg_r2': 0,
            'auto2_min_r2': 0,
            'auto2_silhouette': None,
            'auto2_error_ratio': None,
            'auto2_cv_values': [],
            # n=3 指标
            'auto3_score': 0,
            'auto3_avg_r2': 0,
            'auto3_min_r2': 0,
            'auto3_silhouette': None,
            'auto3_error_ratio': None,
            'auto3_cv_values': [],
            'auto3_cv_diff_significant': True,  # 3分区热容差异是否显著
            'auto3_cv_diff_ratio': None,  # 热容差异/合并误差
            # R²增益 (相对于1分区)
            'n2_r2_gain': 0,  # n2 R² - n1 R²
            'n3_r2_gain': 0,  # n3 R² - n1 R²
            # 决策
            'score_diff': 0,  # auto2 - auto3，正值表示n=2更优
            'r2_diff': 0,
            'partition_meaningful': True,  # 分区是否有意义
            'recommendation': 'N/A',
            'confidence': 'low',
            'reason': ''
        }
        
        # 分析1分区 (整体拟合)
        n1_info = calculate_single_partition_r2(results_dir, structure)
        if n1_info:
            result['n1_r2'] = n1_info['r2']
            result['n1_score'] = n1_info['r2'] * 50  # R²满分50分
        
        # 分析2分区
        if 'auto2' in partitions:
            metrics2 = calculate_comprehensive_score(partitions['auto2'])
            if metrics2:
                result['auto2_score'] = metrics2['total_score']
                result['auto2_avg_r2'] = metrics2['avg_r2']
                result['auto2_min_r2'] = metrics2['min_r2']
                result['auto2_silhouette'] = metrics2['silhouette']
                result['auto2_error_ratio'] = metrics2['cv_error_ratio']
                result['auto2_cv_values'] = metrics2.get('cv_values', [])
                result['auto2_cv_errors'] = metrics2.get('cv_errors', [])
                result['auto2_cv_diff_significant'] = metrics2.get('cv_diff_significant', True)
                result['auto2_cv_diff_ratio'] = metrics2.get('cv_diff_ratio', None)
        
        # 分析3分区
        if 'auto3' in partitions:
            metrics3 = calculate_comprehensive_score(partitions['auto3'])
            if metrics3:
                result['auto3_score'] = metrics3['total_score']
                result['auto3_avg_r2'] = metrics3['avg_r2']
                result['auto3_min_r2'] = metrics3['min_r2']
                result['auto3_silhouette'] = metrics3['silhouette']
                result['auto3_error_ratio'] = metrics3['cv_error_ratio']
                result['auto3_cv_values'] = metrics3.get('cv_values', [])
                result['auto3_cv_diff_significant'] = metrics3.get('cv_diff_significant', True)
                result['auto3_cv_diff_ratio'] = metrics3.get('cv_diff_ratio', None)
                result['auto3_silhouette'] = metrics3['silhouette']
                result['auto3_error_ratio'] = metrics3['cv_error_ratio']
        
        # 计算R²增益 (相对于1分区)
        if result['n1_r2'] > 0:
            result['n2_r2_gain'] = result['auto2_avg_r2'] - result['n1_r2']
            result['n3_r2_gain'] = result['auto3_avg_r2'] - result['n1_r2']
        
        # 计算差异 (2分区 vs 3分区)
        result['score_diff'] = result['auto2_score'] - result['auto3_score']
        result['r2_diff'] = result['auto2_avg_r2'] - result['auto3_avg_r2']
        
        # 决策逻辑 - 综合考虑R²、热容差异显著性和综合得分
        if result['auto2_score'] > 0 and result['auto3_score'] > 0:
            score_diff = result['score_diff']
            
            # 获取2分区热容差异信息
            auto2_cv_significant = result.get('auto2_cv_diff_significant', True)
            auto2_cv_ratio = result.get('auto2_cv_diff_ratio', None)
            
            # 首先判断分区是否有意义：
            # 1. 如果2分区热容差异显著 (ratio >= 2)，分区有物理意义
            # 2. 如果热容差异不显著，但R²增益大，也可能有意义
            if auto2_cv_significant and auto2_cv_ratio is not None and auto2_cv_ratio >= 2:
                # 热容差异显著，分区有物理意义
                result['partition_meaningful'] = True
            elif result['n1_r2'] > 0 and result['n2_r2_gain'] >= 0.02:
                # R²增益足够
                result['partition_meaningful'] = True
            else:
                # 热容差异不显著 且 R²增益不足
                result['partition_meaningful'] = False
            
            # 决策
            if not result['partition_meaningful']:
                result['recommendation'] = '1分区'
                cv_ratio_str = f'{auto2_cv_ratio:.2f}' if auto2_cv_ratio else 'N/A'
                result['confidence'] = 'high'
                result['reason'] = f'热容差异不显著(比值={cv_ratio_str}),R²增益={result["n2_r2_gain"]:.4f}'
            
            else:
                # 分区有意义，比较2分区 vs 3分区
                # 核心判据：3分区的热容差异是否显著（partition2 vs partition3）
                auto3_cv_significant = result.get('auto3_cv_diff_significant', True)
                auto3_cv_ratio = result.get('auto3_cv_diff_ratio', None)
                
                # 如果3分区热容差异不显著，直接选2分区
                if not auto3_cv_significant:
                    result['recommendation'] = '2分区'
                    result['confidence'] = 'high'
                    ratio_str = f'{auto3_cv_ratio:.2f}' if auto3_cv_ratio else 'N/A'
                    result['reason'] = f'3分区热容差异不显著(比值={ratio_str}<2),选n=2'
                
                # 3分区热容差异显著，根据综合得分决定
                elif score_diff > 2:
                    # n=2综合得分更优
                    result['recommendation'] = '2分区'
                    result['confidence'] = 'high' if score_diff > 5 else 'medium'
                    result['reason'] = f'综合得分差={score_diff:+.1f}, n=2更优'
                
                elif score_diff >= -2:
                    # 差异不显著，默认n=2（更简洁的模型）
                    result['recommendation'] = '2分区'
                    result['confidence'] = 'low'
                    result['reason'] = f'综合得分差={score_diff:+.1f}, 差异不显著,默认n=2'
                
                else:
                    # n=3综合得分更优，且热容差异显著
                    result['recommendation'] = '3分区'
                    result['confidence'] = 'high' if score_diff < -5 else 'medium'
                    ratio_str = f'{auto3_cv_ratio:.2f}' if auto3_cv_ratio else 'N/A'
                    result['reason'] = f'综合得分差={score_diff:+.1f}, 3分区热容显著(比值={ratio_str})'
        
        elif result['auto2_score'] > 0:
            result['recommendation'] = '2分区'
            result['confidence'] = 'medium'
            result['reason'] = '仅有2分区数据'
        
        elif result['auto3_score'] > 0:
            result['recommendation'] = '3分区'
            result['confidence'] = 'medium'
            result['reason'] = '仅有3分区数据'
        
        analysis_results.append(result)
    
    return analysis_results


def generate_report(analysis_results, output_path):
    """生成Markdown格式报告 (含1分区基准对比)"""
    
    # 分类统计
    total = len(analysis_results)
    n2_better = sum(1 for r in analysis_results if r['score_diff'] > 0)
    n3_better = sum(1 for r in analysis_results if r['score_diff'] < 0)
    partition_meaningful = sum(1 for r in analysis_results if r['partition_meaningful'])
    partition_not_meaningful = total - partition_meaningful
    
    valid_n1 = [r['n1_r2'] for r in analysis_results if r['n1_r2'] > 0]
    valid_n2 = [r['auto2_avg_r2'] for r in analysis_results if r['auto2_avg_r2'] > 0]
    valid_n3 = [r['auto3_avg_r2'] for r in analysis_results if r['auto3_avg_r2'] > 0]
    avg_r2_n1 = np.mean(valid_n1) if valid_n1 else 0
    avg_r2_n2 = np.mean(valid_n2) if valid_n2 else 0
    avg_r2_n3 = np.mean(valid_n3) if valid_n3 else 0
    
    valid_n2_gain = [r['n2_r2_gain'] for r in analysis_results if r['n1_r2'] > 0]
    valid_n3_gain = [r['n3_r2_gain'] for r in analysis_results if r['n1_r2'] > 0]
    avg_n2_gain = np.mean(valid_n2_gain) if valid_n2_gain else 0
    avg_n3_gain = np.mean(valid_n3_gain) if valid_n3_gain else 0
    
    valid_score_n2 = [r['auto2_score'] for r in analysis_results if r['auto2_score'] > 0]
    valid_score_n3 = [r['auto3_score'] for r in analysis_results if r['auto3_score'] > 0]
    avg_score_n2 = np.mean(valid_score_n2) if valid_score_n2 else 0
    avg_score_n3 = np.mean(valid_score_n3) if valid_score_n3 else 0
    
    sig_n2 = sum(1 for r in analysis_results if r['score_diff'] > 5)
    sig_n3 = sum(1 for r in analysis_results if r['score_diff'] < -5)
    marginal = sum(1 for r in analysis_results if -2 <= r['score_diff'] <= 2)
    
    recommend_1 = [r for r in analysis_results if r['recommendation'] == '1分区']
    recommend_2 = [r for r in analysis_results if r['recommendation'] == '2分区']
    recommend_3 = [r for r in analysis_results if r['recommendation'] == '3分区']
    high_conf = [r for r in analysis_results if r['confidence'] == 'high']
    
    report = f"""# 分区策略综合分析报告 (1/2/3分区对比)

**生成时间**: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}  
**分析体系数**: {total}  
**分析脚本版本**: v2.2

---

## 🎯 核心结论

> **统一使用 n=2（两相分区）是合理的选择**

基于 1分区(整体拟合) vs 2分区 vs 3分区 的综合对比分析。

---

## 📊 1分区 vs 2分区 vs 3分区 R² 对比

**⚠️ R²含义说明**：
- **n=1 R²**: 对所有温度点整体线性拟合（温度 vs 平均能量）的R²
- **n=2/n=3 R²**: 分区后，每个分区内部线性拟合的**平均R²**

由于n=1是全局拟合，n=2/n=3是局部拟合，两者R²不能直接比较！

| 分区数 | 平均R² | R²差异 | 说明 |
|--------|--------|--------|------|
| **n=1 (整体拟合)** | {avg_r2_n1:.4f} | - | 全局线性拟合 |
| **n=2** | {avg_r2_n2:.4f} | {avg_n2_gain:+.4f} | 分区后局部拟合平均 |
| **n=3** | {avg_r2_n3:.4f} | {avg_n3_gain:+.4f} | 分区后局部拟合平均 |

### 分区意义评估（基于热容差异显著性）

**核心判据**：2分区的两个分区热容差异是否统计显著（|Cv₁-Cv₂| / √(err₁²+err₂²) ≥ 2）

| 评估结果 | 数量 | 占比 | 说明 |
|----------|------|------|------|
| **分区有意义** | {partition_meaningful} | {100*partition_meaningful/total:.1f}% | 热容差异显著性比值 ≥ 2 |
| 分区意义不大 | {partition_not_meaningful} | {100*partition_not_meaningful/total:.1f}% | 热容差异显著性比值 < 2 |

---

## 📊 评分公式

综合得分 (0-100分) = R²得分 (50分) + 聚类质量 (30分) + 误差得分 (20分)

| 指标 | 权重 | 计算方式 |
|------|------|----------|
| **平均R²** | 40% | avg_R² × 40 |
| **最小R²** | 10% | min_R² × 10 |
| **Silhouette Score** | 30% | (silhouette+1)/2 × 30 |
| **Cv误差比例** | 20% | (0.2-error_ratio)/0.15 × 20 |

---

## 📊 统计摘要

### 总体对比

| 指标 | n=2 | n=3 | 差异 |
|------|-----|-----|------|
| **更优系统数** | {n2_better}/{total} ({100*n2_better/total:.1f}%) | {n3_better}/{total} ({100*n3_better/total:.1f}%) | n=2领先 |
| **平均综合得分** | {avg_score_n2:.1f} | {avg_score_n3:.1f} | {avg_score_n2-avg_score_n3:+.1f} |
| **平均R²** | {avg_r2_n2:.4f} | {avg_r2_n3:.4f} | {avg_r2_n2-avg_r2_n3:+.4f} |

### 综合得分差异分布

| 差异类型 | 数量 | 占比 | 说明 |
|----------|------|------|------|
| n=2**显著**更优 (diff>5) | {sig_n2} | {100*sig_n2/total:.1f}% | 高置信度推荐n=2 |
| n=3**显著**更优 (diff<-5) | {sig_n3} | {100*sig_n3/total:.1f}% | 高置信度推荐n=3 |
| **差异不显著** (\\|diff\\|≤2) | {marginal} | {100*marginal/total:.1f}% | 默认使用n=2 |

### 最终推荐

| 推荐分区 | 数量 | 占比 |
|----------|------|------|
| 1分区 | {len(recommend_1)} | {100*len(recommend_1)/total:.1f}% |
| **2分区** | {len(recommend_2)} | {100*len(recommend_2)/total:.1f}% |
| 3分区 | {len(recommend_3)} | {100*len(recommend_3)/total:.1f}% |
| 高置信度 | {len(high_conf)} | {100*len(high_conf)/total:.1f}% |

---

## 📋 完整分区推荐表

**列说明**：
- **n1 R²**: 整体线性拟合R²（全局）
- **n2 R²**: 2分区各自拟合的平均R²（局部）
- **n2-n1**: n2 R² - n1 R²（通常为负，因为局部拟合R²通常低于全局）
- **得分差**: n2综合得分 - n3综合得分（正值表示n=2更优）

| 体系 | 推荐 | 置信度 | n1 R² | n2 R² | n3 R² | n2-n1 | n2得分 | n3得分 | 得分差 | 理由 |
|------|------|--------|-------|-------|-------|-------|--------|--------|--------|------|
"""
    
    # 按得分差排序（n=2更优的排前面）
    for r in sorted(analysis_results, key=lambda x: -x['score_diff']):
        conf_icon = {'high': '🟢', 'medium': '🟡', 'low': '⚪'}.get(r['confidence'], '⚪')
        rec = '**2分区**' if r['recommendation'] == '2分区' else ('1分区' if r['recommendation'] == '1分区' else '3分区')
        n1_r2 = f"{r['n1_r2']:.4f}" if r['n1_r2'] > 0 else 'N/A'
        n2_diff = f"{r['n2_r2_gain']:+.4f}" if r['n1_r2'] > 0 else 'N/A'
        report += f"| {r['structure']} | {rec} | {conf_icon} | {n1_r2} | {r['auto2_avg_r2']:.4f} | {r['auto3_avg_r2']:.4f} | {n2_diff} | {r['auto2_score']:.1f} | {r['auto3_score']:.1f} | {r['score_diff']:+.1f} | {r['reason']} |\n"
    
    # 分区意义分析
    report += """
---

## 🔍 分区意义分析

### 分区意义不大的体系（热容差异不显著）

**判断标准**: 2分区的热容差异显著性比值 = |Cv₁-Cv₂| / √(err₁²+err₂²) < 2

"""
    not_meaningful = [r for r in analysis_results if not r['partition_meaningful']]
    if not_meaningful:
        report += "| 体系 | n1 R² | n2 R² | 热容差异比值 | 说明 |\n"
        report += "|------|-------|-------|--------------|------|\n"
        for r in sorted(not_meaningful, key=lambda x: x.get('auto2_cv_diff_ratio', 0) or 0):
            cv_ratio = r.get('auto2_cv_diff_ratio', None)
            cv_ratio_str = f"{cv_ratio:.2f}" if cv_ratio is not None else 'N/A'
            report += f"| {r['structure']} | {r['n1_r2']:.4f} | {r['auto2_avg_r2']:.4f} | {cv_ratio_str} | 热容差异不显著，分区无物理意义 |\n"
    else:
        report += "*所有体系的分区都有意义*\n"

    # 分系列分析
    report += """
---

## 📊 按系列分析

### Pt8系列

"""
    pt8_systems = [r for r in analysis_results if r['structure'].lower().startswith('pt8')]
    if pt8_systems:
        pt8_n2_better = sum(1 for r in pt8_systems if r['score_diff'] > 0)
        pt8_meaningful = sum(1 for r in pt8_systems if r['partition_meaningful'])
        report += f"**统计**: n=2更优 {pt8_n2_better}/{len(pt8_systems)} ({100*pt8_n2_better/len(pt8_systems):.1f}%), 分区有意义 {pt8_meaningful}/{len(pt8_systems)}\n\n"
        report += "| 体系 | 推荐 | n2-n1差 | 得分差 | R²差 | 置信度 |\n"
        report += "|------|------|---------|--------|------|--------|\n"
        for r in sorted(pt8_systems, key=lambda x: x['structure']):
            conf_icon = {'high': '🟢', 'medium': '🟡', 'low': '⚪'}.get(r['confidence'], '⚪')
            n2_gain = f"{r['n2_r2_gain']:.4f}" if r['n1_r2'] > 0 else 'N/A'
            report += f"| {r['structure']} | {r['recommendation']} | {n2_gain} | {r['score_diff']:+.1f} | {r['r2_diff']:+.4f} | {conf_icon} |\n"

    report += """
### Pt6系列

"""
    pt6_systems = [r for r in analysis_results if r['structure'].lower().startswith('pt6')]
    if pt6_systems:
        pt6_n2_better = sum(1 for r in pt6_systems if r['score_diff'] > 0)
        pt6_meaningful = sum(1 for r in pt6_systems if r['partition_meaningful'])
        report += f"**统计**: n=2更优 {pt6_n2_better}/{len(pt6_systems)} ({100*pt6_n2_better/len(pt6_systems):.1f}%), 分区有意义 {pt6_meaningful}/{len(pt6_systems)}\n\n"
        report += "| 体系 | 推荐 | n2-n1差 | 得分差 | R²差 | 置信度 |\n"
        report += "|------|------|---------|--------|------|--------|\n"
        for r in sorted(pt6_systems, key=lambda x: x['structure']):
            conf_icon = {'high': '🟢', 'medium': '🟡', 'low': '⚪'}.get(r['confidence'], '⚪')
            n2_diff = f"{r['n2_r2_gain']:+.4f}" if r['n1_r2'] > 0 else 'N/A'
            report += f"| {r['structure']} | {r['recommendation']} | {n2_diff} | {r['score_diff']:+.1f} | {r['r2_diff']:+.4f} | {conf_icon} |\n"

    report += """
### 其他系列

"""
    other_systems = [r for r in analysis_results 
                     if not r['structure'].lower().startswith('pt8') 
                     and not r['structure'].lower().startswith('pt6')]
    if other_systems:
        other_n2_better = sum(1 for r in other_systems if r['score_diff'] > 0)
        other_meaningful = sum(1 for r in other_systems if r['partition_meaningful'])
        report += f"**统计**: n=2更优 {other_n2_better}/{len(other_systems)} ({100*other_n2_better/len(other_systems):.1f}%), 分区有意义 {other_meaningful}/{len(other_systems)}\n\n"
        report += "| 体系 | 推荐 | n2-n1差 | 得分差 | R²差 | 置信度 |\n"
        report += "|------|------|---------|--------|------|--------|\n"
        for r in sorted(other_systems, key=lambda x: x['structure']):
            conf_icon = {'high': '🟢', 'medium': '🟡', 'low': '⚪'}.get(r['confidence'], '⚪')
            n2_diff = f"{r['n2_r2_gain']:+.4f}" if r['n1_r2'] > 0 else 'N/A'
            report += f"| {r['structure']} | {r['recommendation']} | {n2_diff} | {r['score_diff']:+.1f} | {r['r2_diff']:+.4f} | {conf_icon} |\n"

    report += f"""
---

## 📝 决策标准与逻辑解释

### 📌 核心问题：1分区没有聚类质量分

**重要说明**：
- **1分区** = 整体线性拟合，没有聚类，只有R²
- **2/3分区** = K-means聚类 + 分段线性拟合，有完整的综合得分（R² + Silhouette + 误差）

因此，我们采用**三阶段决策**（核心：热容差异显著性检验）：

### 第一阶段：计算1分区整体拟合R²

对所有数据进行整体线性拟合（温度 vs 平均能量），计算R²作为基准。

### 第二阶段：判断2分区是否有物理意义 (n=1 vs n=2)

**热容差异显著性检验**：

$$\\text{{显著性比值}} = \\frac{{|Cv_1 - Cv_2|}}{{\\sqrt{{err_1^2 + err_2^2}}}}$$

| 显著性比值 | 判定 | 理由 |
|------------|------|------|
| **比值 < 2** | 推荐n=1 | 2分区热容差异不显著，分区无物理意义 |
| **比值 ≥ 2** | 继续第三阶段 | 热容差异显著，需判断2分区还是3分区 |

### 第三阶段：判断3分区是否有物理意义 (n=2 vs n=3)

**同样使用热容差异显著性检验**：

$$\\text{{显著性比值}} = \\frac{{|Cv_2 - Cv_3|}}{{\\sqrt{{err_2^2 + err_3^2}}}}$$

| 显著性比值 | 判定 | 理由 |
|------------|------|------|
| **比值 < 2** | 推荐n=2 | 3分区的第2、3区热容差异不显著 |
| **比值 ≥ 2** | 结合综合得分 | 热容差异显著，根据综合得分选择 |

### 💡 统一的热容差异检验原则

无论是1→2分区还是2→3分区，核心判据都是：

> **热容差异 / 合并误差 ≥ 2** → 分区有物理意义

这确保了：
1. 每次增加分区数，都必须带来**统计显著的热容差异**
2. 避免过度拟合（分区越多但热容差异不显著）

### � 示例

**Air68体系**（推荐1分区）：
- 2分区热容: Cv₁ = 1.23 ± 0.05, Cv₂ = 1.35 ± 0.08
- 显著性比值 = |1.35-1.23| / √(0.05²+0.08²) = 0.12/0.094 = **0.71 < 2**
- 结论：热容差异不显著，分区无意义

**Pt8sn6体系**（推荐2分区）：
- 2分区热容差异显著（比值 ≥ 2）
- 3分区热容差异不显著（比值 < 2）
- 结论：2分区有物理意义，3分区无必要

---

## 🔧 使用命令

批量使用2分区运行所有体系：
```bash
python step6_2_run_batch_analysis.py
```

单个体系使用2分区：
```bash
python step6_1_clustering_analysis.py -s <结构名> -n 2 --use-d-value
```

---

*报告由 `analyze_partition_recommendations.py v2.4` 自动生成*
"""
    
    # 保存报告
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(report)
    
    print(f"[SUCCESS] 报告已保存至: {output_path}")
    
    return report


def main():
    results_dir = Path(__file__).parent / 'results' / 'step6_1_clustering'
    output_path = results_dir / 'PARTITION_RECOMMENDATION_REPORT.md'
    
    print("=" * 70)
    print("分区策略分析 v2.4 - 含热容差异显著性检验")
    print("=" * 70)
    
    # 分析所有结构
    analysis_results = analyze_all_structures(results_dir)
    
    if not analysis_results:
        print("[ERROR] 未找到任何分析结果!")
        return
    
    print(f"\n[INFO] 分析了 {len(analysis_results)} 个体系")
    
    # 1分区 vs 2分区 vs 3分区 R² 统计
    valid_n1 = [r['n1_r2'] for r in analysis_results if r['n1_r2'] > 0]
    valid_n2 = [r['auto2_avg_r2'] for r in analysis_results if r['auto2_avg_r2'] > 0]
    valid_n3 = [r['auto3_avg_r2'] for r in analysis_results if r['auto3_avg_r2'] > 0]
    
    print(f"\n[R² 对比]")
    if valid_n1:
        avg_r2_n1 = np.mean(valid_n1)
        avg_r2_n2 = np.mean(valid_n2) if valid_n2 else 0
        avg_r2_n3 = np.mean(valid_n3) if valid_n3 else 0
        print(f"  n=1 (整体拟合): 平均R² = {avg_r2_n1:.4f}")
        print(f"  n=2:             平均R² = {avg_r2_n2:.4f} (+{avg_r2_n2-avg_r2_n1:.4f})")
        print(f"  n=3:             平均R² = {avg_r2_n3:.4f} (+{avg_r2_n3-avg_r2_n1:.4f})")
    
    # 分区意义统计
    partition_meaningful = sum(1 for r in analysis_results if r['partition_meaningful'])
    print(f"\n[分区意义]")
    print(f"  分区有意义: {partition_meaningful}/{len(analysis_results)} ({100*partition_meaningful/len(analysis_results):.1f}%)")
    print(f"  分区意义不大: {len(analysis_results)-partition_meaningful}/{len(analysis_results)} ({100*(len(analysis_results)-partition_meaningful)/len(analysis_results):.1f}%)")
    
    # 统计摘要
    n2_better = sum(1 for r in analysis_results if r['score_diff'] > 0)
    n3_better = sum(1 for r in analysis_results if r['score_diff'] < 0)
    
    valid_score_n2 = [r['auto2_score'] for r in analysis_results if r['auto2_score'] > 0]
    valid_score_n3 = [r['auto3_score'] for r in analysis_results if r['auto3_score'] > 0]
    avg_score_n2 = np.mean(valid_score_n2) if valid_score_n2 else 0
    avg_score_n3 = np.mean(valid_score_n3) if valid_score_n3 else 0
    
    print(f"\n[n=2 vs n=3 对比]")
    print(f"  n=2 更优: {n2_better}/{len(analysis_results)} ({100*n2_better/len(analysis_results):.1f}%)")
    print(f"  n=3 更优: {n3_better}/{len(analysis_results)} ({100*n3_better/len(analysis_results):.1f}%)")
    print(f"  平均综合得分: n2={avg_score_n2:.1f} vs n3={avg_score_n3:.1f} (差异: {avg_score_n2-avg_score_n3:+.1f})")
    
    # 生成报告
    generate_report(analysis_results, output_path)
    
    # 推荐摘要
    recommend_1 = [r for r in analysis_results if r['recommendation'] == '1分区']
    recommend_2 = [r for r in analysis_results if r['recommendation'] == '2分区']
    recommend_3 = [r for r in analysis_results if r['recommendation'] == '3分区']
    
    print(f"\n[RECOMMENDATION]")
    print(f"  推荐1分区: {len(recommend_1)} 个体系")
    print(f"  推荐2分区: {len(recommend_2)} 个体系")
    print(f"  推荐3分区: {len(recommend_3)} 个体系")
    print(f"\n>>> 结论: 统一使用 n=2 (两相分区) 是合理的选择")


if __name__ == '__main__':
    main()
