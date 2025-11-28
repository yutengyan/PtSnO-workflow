#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
Step 6.1: Lindemann指数相态分区分析 (聚类版本)
================================================================================

作者: GitHub Copilot
日期: 2025-10-22
版本: v1.6

更新日志:
========
v1.6 (2025-11-27):
  - 🔧 修复D值匹配策略: 使用完整路径签名匹配 (保留4级路径)
  - ✅ 匹配率: 100% (3262条全部匹配)
  - 📝 路径签名格式:
    · 3级: batch/composition/run (如 pt8-2/pt8sn5-1-best/t1000.r7.gpu0)
    · 4级: runX/parent/composition/run (如 run3/o1/g-1-o1sn4pt3/t1000.r25.gpu0)

v1.5 (2025-11-26):
  - 使用all_runs_D_values.csv替代ensemble_analysis_results.csv
  - 按路径签名匹配D值,区分不同batch的相同composition

v1.3 (2025-10-22):
  - Cv-1~Cv-5自动合并为单一'Cv'体系
  - 支持多种聚类算法

功能概述:
========
本脚本是Step 6.0聚类功能的独立版本,用于:
1. 自动检测相边界 (替代固定的0.1和0.15阈值)
2. 支持多种聚类算法 (K-means, Hierarchical, DBSCAN)
3. 自动确定最优分区数 (1-3个分区,基于物理意义)
4. 多维特征空间聚类 (温度, Lindemann-δ, 能量, MSD, 扩散系数D)
5. 批量分析多个结构
6. 生成详细的对比报告和可视化

物理约束:
========
- 最大分区数 = 3 (固态、预熔化、液态)
- 基于物理意义,不允许过度拟合
- 单相态: n_partitions=1 (均匀相)
- 两相态: n_partitions=2 (固-液或固-预熔)
- 三相态: n_partitions=3 (固-预熔-液)

数据来源:
========
⚠️ **前置步骤**: 必须先运行 step6_0_multi_system_heat_capacity.py 生成基础数据！

1. 基础数据 (step7_4_all_systems_data.csv):
   来源: step6_0_multi_system_heat_capacity.py输出
   路径: results/step7_4_multi_system/step7_4_all_systems_data.csv
   字段:
   - match_key: 路径签名 (3级或4级,用于跨步骤匹配)
   - structure: 结构名 (如Cv, pt8sn0-2-best, Pt5Sn3O1)
     · 注: Cv-1~Cv-5已自动合并为单一'Cv'体系(v1.3+)
   - system_type: 系列分类 (Cv, Pt8SnX, PtxSnyOz等)
   - system_id: 细粒度标识 (Cv, Pt8Sn0, Pt5Sn3O1等)
   - temp: LAMMPS模拟温度 (K)
   - avg_energy: LAMMPS TotEng平均值 (eV, 团簇+载体)
   - energy_std: 能量标准差 (eV)
   - delta: Lindemann指数 (无量纲, 基于Pt-Sn距离MSD)
   - phase: step6.0固定阈值分类 (solid/premelting/liquid)
   - run_id: 运行标识 (如r15.gpu0)
   统计: 51个结构, 3262条记录

2. 扩散系数D值数据 (all_runs_D_values.csv):
   来源: Step4 MSD分析输出
   路径: results/all_runs_D_values.csv
   字段: filepath, composition, element, gmx_D, our_D等
   匹配策略 (v1.6):
   - 从filepath提取路径签名,与主数据match_key匹配
   - 3级路径: parts[-4]/parts[-3]/run_info
   - 4级路径: 检测run3等批次标识,生成runX/parent/comp/run_info
   - 匹配率: 100% (3262/3262)
   - 使用: --use-d-value参数启用

聚类特征:
========
- 基础 (2D): Temperature + Lindemann-δ
- 扩展 (3D): + Energy (--use-energy)
- 扩展 (3D): + MSD (--use-msd, 需数据支持)
- 扩展 (3D): + Diffusion-D (--use-d-value)
- 高维 (4D+): 任意组合上述特征

优势:
====
- ✅ 独立脚本,便于测试和调试
- ✅ 专注于相态分区,不受主脚本复杂度影响
- ✅ 多维特征空间优化分区精度
- ✅ D值按路径签名精确匹配 (v1.6+)
- ✅ 可作为独立工具使用

使用示例:
========
# 1. 基础2D聚类 (温度 + Lindemann-δ)
python step6_1_clustering_analysis.py --structure pt6sn8 --n-partitions 3

# 2. 3D聚类 (加入能量特征)
python step6_1_clustering_analysis.py --structure pt6sn8 --n-partitions 3 --use-energy

# 3. 3D聚类 (加入扩散系数D值)
python step6_1_clustering_analysis.py --structure pt6sn8 --n-partitions 3 --use-d-value

# 4. 4D聚类 (温度 + δ + 能量 + D值)
python step6_1_clustering_analysis.py --structure pt6sn8 --n-partitions 3 --use-energy --use-d-value

# 5. 自动确定最优分区数
python step6_1_clustering_analysis.py --structure pt6sn8 --auto-partition --use-energy

# 6. 批量分析所有结构
python step6_1_clustering_analysis.py --structure all --auto-partition --use-energy

================================================================================
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import linregress
import argparse
import warnings
warnings.filterwarnings('ignore')

# 机器学习相关
from sklearn.cluster import KMeans, DBSCAN, AgglomerativeClustering
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score
from scipy.cluster.hierarchy import dendrogram, linkage

# 中文字体设置 - 添加更多备选字体以支持特殊字符
plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans', 'Arial Unicode MS', 'Microsoft YaHei']
plt.rcParams['axes.unicode_minus'] = False
# 启用数学文本渲染以支持上标
plt.rcParams['mathtext.default'] = 'regular'

# ============================================================================
# 配置路径
# ============================================================================

BASE_DIR = Path(__file__).parent
DATA_FILE = BASE_DIR / 'results' / 'step7_4_multi_system' / 'step7_4_all_systems_data.csv'
D_VALUE_FILE = BASE_DIR / 'results' / 'ensemble_analysis_results.csv'  # D值数据源
SUPPORT_ENERGY_FILE = BASE_DIR / 'data' / 'lammps_energy' / 'sup' / 'energy_master_20251021_151520.csv'  # 载体能量数据
OUTPUT_DIR = BASE_DIR / 'results' / 'step7_4_2_clustering'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print("="*80)
print("Step 7.4.2: Lindemann指数聚类分析")
print("="*80)

# ============================================================================
# 1. 数据加载
# ============================================================================

def load_support_energy_data():
    """
    加载载体能量数据并拟合E-T关系
    
    Returns:
        tuple: (slope, intercept, R2) 或 None
               slope: dE/dT (eV/K)
               intercept: E(T=0) (eV)
               R2: 拟合质量
    """
    if not SUPPORT_ENERGY_FILE.exists():
        print(f"  [WARNING] 载体能量文件不存在: {SUPPORT_ENERGY_FILE}")
        print(f"  将使用默认Cv_support=38.2151 meV/K")
        return None
    
    try:
        df_sup = pd.read_csv(SUPPORT_ENERGY_FILE, encoding='utf-8')
        
        # 处理中文列名
        if '结构' in df_sup.columns:
            df_sup.columns = ['path', 'structure', 'temp', 'run_num', 'total_steps', 'sample_steps', 
                              'avg_energy', 'std', 'min', 'max', 'sample_interval', 
                              'skip_steps', 'full_path']
        
        # 筛选载体数据（sup-1, sup-2, sup240等）
        df_sup = df_sup[df_sup['structure'].str.contains('sup', case=False, na=False)]
        
        if len(df_sup) == 0:
            print(f"  [WARNING] 未找到载体数据 (含'sup'关键字)")
            return None
        
        # 按温度平均
        df_sup_avg = df_sup.groupby('temp')['avg_energy'].mean().reset_index()
        df_sup_avg = df_sup_avg.sort_values('temp')
        
        # 线性拟合 E = slope * T + intercept
        from scipy.stats import linregress
        T = df_sup_avg['temp'].values
        E = df_sup_avg['avg_energy'].values
        
        if len(T) < 2:
            print(f"  [WARNING] 载体数据点不足 (n={len(T)})")
            return None
        
        slope, intercept, r_value, p_value, std_err = linregress(T, E)
        R2 = r_value ** 2
        
        print(f"  [OK] 载体能量拟合:")
        print(f"      斜率 (Cv_support): {slope*1000:.4f} meV/K")
        print(f"      截距: {intercept:.2f} eV")
        print(f"      R^2: {R2:.6f}")
        print(f"      温度范围: {T.min():.0f}-{T.max():.0f} K")
        print(f"      数据点: {len(T)}")
        
        return (slope, intercept, R2)
        
    except Exception as e:
        print(f"  [ERROR] 加载载体能量数据失败: {e}")
        return None


def extract_d_signature_from_filepath(filepath):
    """
    从D值数据的filepath提取完整路径签名
    
    路径格式: .../batch/composition/tempK/Txxxx.rN.gpuM_msd_element.xvg
    例: .../Pt8-4/pt8sn5-1-best/1000K/T1000.r8.gpu0_msd_Pt.xvg
    
    签名格式 (与主数据match_key一致):
    - 3级: batch/composition/run_info
      例: pt8-4/pt8sn5-1-best/t1000.r8.gpu0
    - 4级: runX/parent/composition/run_info (当路径包含run3等批次标识时)
      例: run3/o1/g-1-o1sn4pt3/t1000.r25.gpu0
    """
    import re
    parts = filepath.replace('\\', '/').split('/')
    
    if len(parts) < 4:
        return None
    
    # 提取run信息
    run_file = parts[-1]
    run_match = re.match(r'(T\d+\.r\d+\.gpu\d+)', run_file, re.IGNORECASE)
    if not run_match:
        return None
    run_info = run_match.group(1).lower()
    
    # 找温度目录位置 (如 1000K)
    temp_idx = None
    for i, part in enumerate(parts):
        if re.match(r'\d+K$', part, re.IGNORECASE):
            temp_idx = i
            break
    
    if temp_idx is None or temp_idx < 2:
        return None
    
    comp = parts[temp_idx - 1].lower()    # composition目录
    parent = parts[temp_idx - 2].lower()  # parent目录
    
    # 检查是否有run3等批次标识 (向上搜索)
    batch_keywords = ['run3', 'run2', 'run4', 'run5']
    for check_idx in range(temp_idx - 3, max(-1, temp_idx - 6), -1):
        if check_idx >= 0 and parts[check_idx].lower() in batch_keywords:
            # 4级签名: runX/parent/composition/run
            return f'{parts[check_idx].lower()}/{parent}/{comp}/{run_info}'
    
    # 3级签名: parent/composition/run
    return f'{parent}/{comp}/{run_info}'


def extract_main_signature_from_match_key(match_key):
    """
    从主数据的match_key提取完整路径签名
    
    match_key格式 (直接使用，不做修改):
    - 3级: batch/composition/run_info
      例: pt8-2/pt8sn5-1-best/t1000.r7.gpu0
      
    - 4级: runX/parent/composition/run_info  
      例: run3/o1/g-1-o1sn4pt3/t1000.r25.gpu0
    
    签名格式: 直接返回match_key的小写形式
    """
    return match_key.lower() if match_key else None


def load_data(load_d_values=False):
    """
    加载Step 7.4的输出数据,可选合并D值
    
    Args:
        load_d_values: 是否加载并合并扩散系数D值数据
    
    Returns:
        DataFrame: 包含温度、能量、Lindemann-δ等特征,可选D值
    
    D值匹配策略 (v1.6 修复):
    --------------------------
    使用完整路径签名匹配,保留4级路径结构
    
    1. 主数据match_key直接使用 (不去掉runX前缀):
       - 3级: pt8-2/pt8sn5-1-best/t1000.r7.gpu0
       - 4级: run3/o1/g-1-o1sn4pt3/t1000.r25.gpu0
    
    2. D值filepath提取签名 (检测run3等批次标识):
       - 3级: .../Pt8-4/pt8sn5-1-best/1000K/T1000.r8.gpu0_msd_Pt.xvg -> pt8-4/pt8sn5-1-best/t1000.r8.gpu0
       - 4级: .../run3/o1/g-1-O1Sn4Pt3/1000K/T1000.r25.gpu0_msd_Pt.xvg -> run3/o1/g-1-o1sn4pt3/t1000.r25.gpu0
    
    3. 按完整签名匹配,保持路径层级一致
    
    匹配统计:
    - 主数据: 3262条唯一签名
    - D值(Pt): 3262条唯一签名
    - 匹配率: 100%
    """
    print(f"\n[1] 加载数据...")
    
    if not DATA_FILE.exists():
        print(f"  [ERROR] 数据文件不存在: {DATA_FILE}")
        print(f"  请先运行 step7_4_multi_system_heat_capacity.py 生成数据!")
        return None
    
    df = pd.read_csv(DATA_FILE, encoding='utf-8')
    print(f"  [OK] 加载了 {len(df)} 条记录")
    print(f"  结构数: {df['structure'].nunique()}")
    print(f"  温度范围: {df['temp'].min():.0f} - {df['temp'].max():.0f} K")
    print(f"  Lindemann范围: {df['delta'].min():.4f} - {df['delta'].max():.4f}")
    
    # 可选: 合并D值数据 (使用路径签名匹配)
    if load_d_values:
        # 优先使用 all_runs_D_values.csv (有filepath)
        D_VALUE_FILE_WITH_PATH = BASE_DIR / 'results' / 'all_runs_D_values.csv'
        
        if D_VALUE_FILE_WITH_PATH.exists():
            print(f"\n  [+] 加载D值数据 (路径签名匹配)...")
            df_D = pd.read_csv(D_VALUE_FILE_WITH_PATH, encoding='utf-8')
            print(f"      原始D值记录: {len(df_D)}")
            
            # 只取Pt元素的D值 (主要扩散物种)
            df_D_pt = df_D[df_D['element'] == 'Pt'].copy()
            print(f"      Pt元素记录: {len(df_D_pt)}")
            
            # 提取D值数据的签名
            df_D_pt['d_signature'] = df_D_pt['filepath'].apply(extract_d_signature_from_filepath)
            df_D_pt = df_D_pt[df_D_pt['d_signature'].notna()].copy()
            print(f"      有效签名记录: {len(df_D_pt)}")
            
            # 使用 gmx_D 或 our_D 作为D值
            if 'gmx_D' in df_D_pt.columns:
                df_D_pt['D_value'] = df_D_pt['gmx_D']
            elif 'our_D' in df_D_pt.columns:
                df_D_pt['D_value'] = df_D_pt['our_D']
            else:
                print(f"      [ERROR] 找不到D值列 (gmx_D 或 our_D)")
                return df
            
            # 提取主数据的签名
            df['main_signature'] = df['match_key'].apply(extract_main_signature_from_match_key)
            
            # 准备D值查找表 (按签名)
            d_lookup = df_D_pt.set_index('d_signature')['D_value'].to_dict()
            
            # 按签名匹配D值
            df['D_value'] = df['main_signature'].map(d_lookup)
            
            n_matched = df['D_value'].notna().sum()
            print(f"      成功匹配D值: {n_matched}/{len(df)} ({n_matched/len(df)*100:.1f}%)")
            
            if n_matched > 0:
                print(f"      D值范围: {df['D_value'].min():.2e} - {df['D_value'].max():.2e} cm^2/s")
                
                # 检查负D值 (物理上不合理,低温拟合误差)
                n_negative = (df['D_value'] < 0).sum()
                if n_negative > 0:
                    print(f"      [WARNING] 检测到 {n_negative} 个负D值 (低温拟合误差),将设为0")
                    df.loc[df['D_value'] < 0, 'D_value'] = 0
                
                # 展示匹配示例
                matched_sample = df[df['D_value'].notna()].head(3)
                print(f"\n      匹配示例:")
                for _, row in matched_sample.iterrows():
                    print(f"        {row['main_signature']} -> D={row['D_value']:.2e}")
            else:
                print(f"      [WARNING] 无法匹配任何D值!")
                print(f"      主数据签名示例: {df['main_signature'].head(3).tolist()}")
                print(f"      D值签名示例: {list(d_lookup.keys())[:3]}")
            
            # 清理临时列
            df.drop(columns=['main_signature'], inplace=True, errors='ignore')
            
        elif D_VALUE_FILE.exists():
            # 回退到旧的ensemble数据 (按温度匹配,不推荐)
            print(f"\n  [+] 加载D值数据 (回退: 按温度匹配,不推荐)...")
            print(f"      [WARNING] 未找到 all_runs_D_values.csv,使用 ensemble_analysis_results.csv")
            print(f"      [WARNING] 按温度取中位数可能导致不同体系使用相同D值!")
            
            df_D = pd.read_csv(D_VALUE_FILE, encoding='utf-8')
            df_D_pt = df_D[df_D['element'] == 'Pt'].copy()
            df_D_pt.rename(columns={'D_ensemble': 'D_value'}, inplace=True)
            
            df_D_grouped = df_D_pt.groupby('temp_value').agg({
                'D_value': 'median'
            }).reset_index()
            df_D_grouped.rename(columns={'temp_value': 'temp'}, inplace=True)
            
            df = df.merge(df_D_grouped[['temp', 'D_value']], on='temp', how='left')
            
            n_matched = df['D_value'].notna().sum()
            print(f"      成功匹配D值: {n_matched}/{len(df)} ({n_matched/len(df)*100:.1f}%)")
        else:
            print(f"  [WARNING] D值数据文件不存在")
            print(f"  跳过D值合并,仅使用基础特征")
    
    return df

# ============================================================================
# 2. Elbow Method - 确定最优分区数 (物理约束: 最多3个相态)
# ============================================================================

def determine_optimal_partitions(X, max_partitions=3):
    """
    使用Elbow method确定最优分区数
    
    基于物理意义,相态分区数限制为1-3:
    - 1个分区: 单一相态 (如纯固态或纯液态)
    - 2个分区: 两相共存 (固-液 或 固-预熔)
    - 3个分区: 三相共存 (固-预熔-液)
    
    Args:
        X: 标准化后的特征矩阵
        max_partitions: 最大分区数 (默认3,基于物理意义)
    
    Returns:
        dict: 包含inertia, silhouette等指标
    """
    print(f"\n>>> 自动确定最优分区数 (物理约束: 1-3个相态)...")
    
    inertias = []
    silhouette_scores = []
    calinski_scores = []
    davies_bouldin_scores = []
    
    # 物理约束: 最多3个分区 (固态、预熔化、液态)
    k_range = range(2, min(max_partitions, 3) + 1)
    
    for k in k_range:
        kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
        labels = kmeans.fit_predict(X)
        
        inertias.append(kmeans.inertia_)
        silhouette_scores.append(silhouette_score(X, labels))
        calinski_scores.append(calinski_harabasz_score(X, labels))
        davies_bouldin_scores.append(davies_bouldin_score(X, labels))
    
    # 找到最优分区数 (基于Silhouette score)
    optimal_n = k_range[np.argmax(silhouette_scores)]
    
    # 物理解释
    phase_names = {1: '单相', 2: '两相(固-液)', 3: '三相(固-预熔-液)'}
    
    print(f"  分区数分析:")
    for i, k in enumerate(k_range):
        phase_desc = phase_names.get(k, f'{k}相')
        print(f"    n={k} ({phase_desc}): Inertia={inertias[i]:.2f}, "
              f"Silhouette={silhouette_scores[i]:.4f}, "
              f"Calinski-Harabasz={calinski_scores[i]:.2f}")
    
    print(f"\n  [推荐] 最优分区数: {optimal_n} ({phase_names.get(optimal_n, f'{optimal_n}相')})")
    
    return {
        'k_range': list(k_range),
        'inertias': inertias,
        'silhouette_scores': silhouette_scores,
        'calinski_scores': calinski_scores,
        'davies_bouldin_scores': davies_bouldin_scores,
        'optimal_k': optimal_n  # 保持向后兼容
    }

# ============================================================================
# 3. 聚类分析核心函数
# ============================================================================

def perform_clustering(df_structure, structure_name, method='kmeans', n_partitions=2, 
                       auto_partition=False, eps=0.3, min_samples=5, 
                       use_energy=False, use_msd=False, use_d_value=False):
    """
    执行相态分区分析
    
    Args:
        df_structure: 单个结构的数据
        structure_name: 结构名称
        method: 聚类方法 ('kmeans', 'hierarchical', 'dbscan')
        n_partitions: 分区数量 (1-3,基于物理意义)
        auto_partition: 是否自动确定分区数
        eps: DBSCAN的epsilon参数
        min_samples: DBSCAN的min_samples参数
        use_energy: 是否加入能量特征
        use_msd: 是否加入MSD特征
        use_d_value: 是否加入扩散系数D值特征
    
    Returns:
        dict: 聚类结果
    """
    print(f"\n{'='*80}")
    print(f"相态分区分析: {structure_name} ({len(df_structure)} 数据点)")
    print(f"  方法: {method.upper()}")
    print("="*80)
    
    if len(df_structure) < 10:
        print(f"  [WARNING] 数据点太少 ({len(df_structure)} < 10), 跳过!")
        return None
    
    # 准备特征: 基础特征(温度, Lindemann指数) + 可选特征(能量, MSD, D值)
    feature_cols = ['temp', 'delta']
    feature_names = ['Temperature', 'Lindemann-δ']
    
    if use_energy and 'avg_energy' in df_structure.columns:
        feature_cols.append('avg_energy')
        feature_names.append('Energy')
        print(f"  [+] 加入能量特征")
    
    if use_msd and 'avg_msd' in df_structure.columns:
        feature_cols.append('avg_msd')
        feature_names.append('MSD')
        print(f"  [+] 加入MSD特征")
    
    if use_d_value and 'D_value' in df_structure.columns:
        # 检查D值是否有效 (去除NaN)
        valid_d = df_structure['D_value'].notna().sum()
        if valid_d > 0:
            feature_cols.append('D_value')
            feature_names.append('Diffusion-D')
            print(f"  [+] 加入扩散系数D值特征 (有效数据: {valid_d}/{len(df_structure)})")
            
            # 如果有部分NaN,填充为0 (固态近似)
            if valid_d < len(df_structure):
                df_structure = df_structure.copy()
                df_structure['D_value'].fillna(0, inplace=True)
                print(f"      [INFO] 缺失D值已填充为0 (固态近似)")
        else:
            print(f"  [WARNING] D值全部缺失,跳过此特征")
    
    print(f"  使用特征 ({len(feature_cols)}D): {', '.join(feature_names)}")
    
    X = df_structure[feature_cols].values
    
    # 标准化
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # 自动确定分区数 (物理约束: 最多3个)
    if auto_partition and method in ['kmeans', 'hierarchical']:
        k_analysis = determine_optimal_partitions(X_scaled, max_partitions=3)
        n_partitions = k_analysis['optimal_k']
        print(f"\n  使用自动确定的分区数: {n_partitions}")
    else:
        k_analysis = None
    
    # 物理约束检查
    if n_partitions > 3:
        print(f"  [WARNING] 分区数 {n_partitions} 超过物理上限3,自动调整为3")
        n_partitions = 3
    
    # 执行聚类分区
    df_clustered = df_structure.copy()
    
    if method == 'kmeans':
        clusterer = KMeans(n_clusters=n_partitions, random_state=42, n_init=10)
        labels = clusterer.fit_predict(X_scaled)
        
    elif method == 'hierarchical':
        clusterer = AgglomerativeClustering(n_clusters=n_partitions, linkage='ward')
        labels = clusterer.fit_predict(X_scaled)
        
    elif method == 'dbscan':
        clusterer = DBSCAN(eps=eps, min_samples=min_samples)
        labels = clusterer.fit_predict(X_scaled)
        n_partitions = len(set(labels)) - (1 if -1 in labels else 0)
        print(f"\n  DBSCAN检测到 {n_partitions} 个分区 (噪声点: {sum(labels == -1)})")
        
        if n_partitions < 1:
            print(f"  [WARNING] 分区数太少, 尝试调整eps和min_samples参数")
            return None
        
        if n_partitions > 3:
            print(f"  [WARNING] DBSCAN检测到 {n_partitions} 个分区, 超过物理上限3")
            print(f"  建议: 调大eps参数或使用K-means方法")
    
    df_clustered['cluster'] = labels
    
    # 计算评估指标
    if len(set(labels)) > 1 and -1 not in labels:  # 排除噪声点的情况
        silhouette = silhouette_score(X_scaled, labels)
        calinski = calinski_harabasz_score(X_scaled, labels)
        davies_bouldin = davies_bouldin_score(X_scaled, labels)
        
        print(f"\n>>> 分区质量评估:")
        print(f"  Silhouette Score: {silhouette:.4f} (越高越好, 范围[-1,1])")
        print(f"  Calinski-Harabasz Index: {calinski:.2f} (越高越好)")
        print(f"  Davies-Bouldin Index: {davies_bouldin:.4f} (越低越好)")
    else:
        silhouette = calinski = davies_bouldin = None
    
    # 分析每个分区
    print(f"\n>>> 相态分区统计:")
    cluster_stats = []
    
    for cluster_id in sorted(df_clustered['cluster'].unique()):
        if cluster_id == -1:  # DBSCAN噪声点
            continue
            
        cluster_data = df_clustered[df_clustered['cluster'] == cluster_id]
        
        stats = {
            'cluster_id': cluster_id,
            'n_points': len(cluster_data),
            'delta_mean': cluster_data['delta'].mean(),
            'delta_std': cluster_data['delta'].std(),
            'delta_min': cluster_data['delta'].min(),
            'delta_max': cluster_data['delta'].max(),
            'temp_mean': cluster_data['temp'].mean(),
            'temp_min': cluster_data['temp'].min(),
            'temp_max': cluster_data['temp'].max()
        }
        cluster_stats.append(stats)
        
        print(f"  分区 {cluster_id}: n={stats['n_points']}, "
              f"δ={stats['delta_mean']:.4f}±{stats['delta_std']:.4f} "
              f"(range: {stats['delta_min']:.4f}-{stats['delta_max']:.4f}), "
              f"T_avg={stats['temp_mean']:.1f}K")
    
    # 按Lindemann均值排序 (从低到高: 分区1<分区2<分区3)
    cluster_stats = sorted(cluster_stats, key=lambda x: x['delta_mean'])
    
    # 分配分区标签 (分区1, 分区2, 分区3...)
    cluster_to_phase = {}
    for i, stat in enumerate(cluster_stats):
        partition_label = f'partition{i+1}'
        cluster_to_phase[stat['cluster_id']] = partition_label
    
    df_clustered['phase_clustered'] = df_clustered['cluster'].map(cluster_to_phase)
    
    # 计算分区边界阈值
    thresholds = []
    if n_partitions >= 2:
        for i in range(len(cluster_stats) - 1):
            lower = cluster_stats[i]
            upper = cluster_stats[i + 1]
            threshold = (lower['delta_mean'] + upper['delta_mean']) / 2
            thresholds.append(threshold)
            
            partition_lower = cluster_to_phase[lower['cluster_id']]
            partition_upper = cluster_to_phase[upper['cluster_id']]
            print(f"\n  分区边界 {partition_lower}<->{partition_upper}: δ = {threshold:.4f}")
    
    # 返回结果
    return {
        'structure': structure_name,
        'method': method,
        'n_partitions': n_partitions,
        'labels': labels,
        'cluster_stats': cluster_stats,
        'cluster_to_phase': cluster_to_phase,
        'thresholds': thresholds,
        'df_clustered': df_clustered,
        'X_scaled': X_scaled,
        'scaler': scaler,
        'feature_cols': feature_cols,
        'feature_names': feature_names,
        'k_analysis': k_analysis,
        'metrics': {
            'silhouette': silhouette,
            'calinski_harabasz': calinski,
            'davies_bouldin': davies_bouldin
        }
    }

# ============================================================================
# 4. 热容拟合函数 (基于分区结果)
# ============================================================================

def fit_partition_heat_capacity(df_clustered, Cv_support=38.2151):
    """
    根据分区结果拟合热容（重新拟合团簇能量）
    
    Args:
        df_clustered: 分区后的数据
        Cv_support: 载体热容 (默认38.2151 meV/K, 来自step7.4计算)
    
    Returns:
        dict: 各分区的热容拟合结果
    """
    print(f"\n>>> 基于分区拟合热容...")
    
    # 加载载体能量数据进行重新拟合
    support_fit = load_support_energy_data()
    
    if support_fit is not None:
        slope_support, intercept_support, R2_support = support_fit
        # 注意：这里不打印，避免在策略对比时重复输出
    else:
        # 回退：使用默认Cv_support估算
        slope_support = Cv_support / 1000  # meV/K -> eV/K
        T_min = df_clustered['temp'].min()
        E_total_min = df_clustered[df_clustered['temp'] == T_min]['avg_energy'].mean()
        intercept_support = E_total_min * 0.9 - slope_support * T_min
    
    partitions = df_clustered['phase_clustered'].unique()
    cv_results = {}
    
    for partition in sorted(partitions):
        if partition in ['anomaly_0', 'anomaly_1', 'anomaly_2']:  # 跳过异常分区
            continue
            
        df_part = df_clustered[df_clustered['phase_clustered'] == partition]
        
        if len(df_part) < 3:
            print(f"    {partition}: 数据点太少 (n={len(df_part)}), 跳过")
            continue
        
        T = df_part['temp'].values
        E_total = df_part['avg_energy'].values
        
        # 检查温度多样性
        if len(np.unique(T)) < 2:
            print(f"    {partition}: 温度单一, 无法拟合")
            continue
        
        # 计算载体能量和团簇能量
        E_support = slope_support * T + intercept_support
        E_cluster = E_total - E_support
        
        # 重新拟合团簇能量-温度关系
        slope_cluster, intercept_cluster, r_value_cluster, p_value, std_err_cluster = linregress(T, E_cluster)
        R2_cluster = r_value_cluster ** 2
        
        # 热容计算（基于团簇能量拟合）
        Cv_cluster = slope_cluster * 1000  # eV/K -> meV/K
        Cv_cluster_err = std_err_cluster * 1000
        
        # 同时保留总能量拟合信息（用于对比）
        slope_total, intercept_total, r_value_total, _, std_err_total = linregress(T, E_total)
        Cv_total = slope_total * 1000
        Cv_total_err = std_err_total * 1000
        R2_total = r_value_total ** 2
        
        cv_results[partition] = {
            'n_points': len(df_part),
            'T_range': (T.min(), T.max()),
            'slope': slope_cluster,  # 使用团簇能量拟合的斜率
            'slope_err': std_err_cluster,
            'intercept': intercept_cluster,
            'Cv_total': Cv_total,
            'Cv_total_err': Cv_total_err,
            'Cv_cluster': Cv_cluster,  # 基于团簇能量重新拟合的结果
            'Cv_cluster_err': Cv_cluster_err,
            'R2': R2_cluster,  # 使用团簇能量拟合的R²
            'R2_total': R2_total,  # 保留总能量拟合的R²用于对比
            'p_value': p_value
        }
        
        # 质量评价（基于团簇能量拟合的R²）
        if R2_cluster > 0.95:
            grade = "★★★ Excellent"
            grade_score = 4
        elif R2_cluster > 0.90:
            grade = "★★ Good"
            grade_score = 3
        elif R2_cluster > 0.80:
            grade = "★ Fair"
            grade_score = 2
        else:
            grade = "⚠ Poor"
            grade_score = 1
        
        # 将评分添加到结果字典
        cv_results[partition]['grade'] = grade
        cv_results[partition]['grade_score'] = grade_score
        
        print(f"    {partition}: n={len(df_part)}, "
              f"T={T.min():.0f}-{T.max():.0f}K, "
              f"Cv_cluster={Cv_cluster:.4f}±{Cv_cluster_err:.4f} meV/K (已扣除载体{Cv_support:.4f}), "
              f"R^2={R2_cluster:.4f} {grade}")
    
    return cv_results

# ============================================================================
# 5. 多策略对比函数
# ============================================================================

def compare_partition_strategies(df_structure, structure_name, use_energy=False, use_msd=False, use_d_value=False):
    """
    对比不同分区策略的热容拟合结果
    
    策略:
    1. 固定阈值 (δ=0.10, 0.15)
    2. n=2 自动分区 (两相)
    3. n=3 自动分区 (三相)
    
    Args:
        df_structure: 单个结构的数据
        structure_name: 结构名称
        use_energy: 是否使用能量特征
        use_msd: 是否使用MSD特征
        use_d_value: 是否使用扩散系数D值特征
    
    Returns:
        dict: 各策略的对比结果
    """
    print(f"\n{'='*80}")
    print(f"多策略对比分析: {structure_name}")
    print("="*80)
    
    Cv_support = 38.2151  # meV/K (与step7.4保持一致)
    comparison = {}
    
    # 准备特征
    feature_cols = ['temp', 'delta']
    if use_energy and 'avg_energy' in df_structure.columns:
        feature_cols.append('avg_energy')
    if use_msd and 'avg_msd' in df_structure.columns:
        feature_cols.append('avg_msd')
    if use_d_value and 'D_value' in df_structure.columns:
        # 填充NaN为0 (固态近似)
        df_structure = df_structure.copy()
        df_structure['D_value'].fillna(0, inplace=True)
        feature_cols.append('D_value')
    
    # 策略1: 固定阈值
    print(f"\n[策略1] 固定阈值 (δ=0.10, 0.15)")
    df_fixed = df_structure.copy()
    # phase列已经存在(来自step7.4)，需要映射到partition命名
    phase_to_partition = {'solid': 'partition1', 'premelting': 'partition2', 'liquid': 'partition3'}
    df_fixed_clustered = df_fixed.rename(columns={'phase': 'phase_clustered'})
    df_fixed_clustered['phase_clustered'] = df_fixed_clustered['phase_clustered'].map(phase_to_partition)
    cv_fixed = fit_partition_heat_capacity(df_fixed_clustered, Cv_support)
    comparison['fixed'] = {
        'method': '固定阈值(δ=0.10, 0.15)',
        'thresholds': [0.10, 0.15],
        'cv_results': cv_fixed,
        'df_clustered': df_fixed_clustered
    }
    
    # 策略1b: 固定阈值两分区（δ=0.08）
    print(f"\n[策略1b] 固定阈值两分区 (δ=0.08)")
    df_fixed2 = df_structure.copy()
    # 使用δ=0.08作为两分区的阈值
    df_fixed2['phase_clustered'] = df_fixed2['delta'].apply(
        lambda x: 'partition1' if x < 0.08 else 'partition2'
    )
    cv_fixed2 = fit_partition_heat_capacity(df_fixed2, Cv_support)
    comparison['fixed2'] = {
        'method': '固定阈值(δ=0.08)',
        'thresholds': [0.08],
        'cv_results': cv_fixed2,
        'df_clustered': df_fixed2
    }
    
    # 策略2: n=2自动分区
    print(f"\n[策略2] 两相自动分区 (n=2)")
    from sklearn.cluster import KMeans
    from sklearn.preprocessing import StandardScaler
    from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score
    
    X = df_structure[feature_cols].values
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    kmeans2 = KMeans(n_clusters=2, random_state=42, n_init=10)
    labels2 = kmeans2.fit_predict(X_scaled)
    
    df_n2 = df_structure.copy()
    df_n2['cluster'] = labels2
    
    # 计算质量得分
    silhouette_n2 = silhouette_score(X_scaled, labels2)
    calinski_n2 = calinski_harabasz_score(X_scaled, labels2)
    davies_n2 = davies_bouldin_score(X_scaled, labels2)
    
    print(f"  质量评估: Silhouette={silhouette_n2:.4f}, "
          f"Calinski-Harabasz={calinski_n2:.2f}, Davies-Bouldin={davies_n2:.4f}")
    
    # 按Lindemann排序分配分区（使用partition命名）
    cluster_means = df_n2.groupby('cluster')['delta'].mean().sort_values()
    phase_map2 = {cluster_means.index[0]: 'partition1', cluster_means.index[1]: 'partition2'}
    df_n2['phase_clustered'] = df_n2['cluster'].map(phase_map2)
    
    thresholds2 = []
    for i in range(len(cluster_means)-1):
        thresh = (cluster_means.iloc[i] + cluster_means.iloc[i+1]) / 2
        thresholds2.append(thresh)
    
    cv_n2 = fit_partition_heat_capacity(df_n2, Cv_support)
    comparison['n2'] = {
        'method': '两相自动分区',
        'n_partitions': 2,
        'thresholds': thresholds2,
        'cv_results': cv_n2,
        'df_clustered': df_n2,
        'metrics': {
            'silhouette': silhouette_n2,
            'calinski_harabasz': calinski_n2,
            'davies_bouldin': davies_n2
        }
    }
    
    # 策略3: n=3自动分区
    print(f"\n[策略3] 三相自动分区 (n=3)")
    kmeans3 = KMeans(n_clusters=3, random_state=42, n_init=10)
    labels3 = kmeans3.fit_predict(X_scaled)
    
    df_n3 = df_structure.copy()
    df_n3['cluster'] = labels3
    
    # 计算质量得分
    silhouette_n3 = silhouette_score(X_scaled, labels3)
    calinski_n3 = calinski_harabasz_score(X_scaled, labels3)
    davies_n3 = davies_bouldin_score(X_scaled, labels3)
    
    print(f"  质量评估: Silhouette={silhouette_n3:.4f}, "
          f"Calinski-Harabasz={calinski_n3:.2f}, Davies-Bouldin={davies_n3:.4f}")
    
    # 按Lindemann排序分配分区（使用partition命名）
    cluster_means3 = df_n3.groupby('cluster')['delta'].mean().sort_values()
    phase_map3 = {
        cluster_means3.index[0]: 'partition1',
        cluster_means3.index[1]: 'partition2',
        cluster_means3.index[2]: 'partition3'
    }
    df_n3['phase_clustered'] = df_n3['cluster'].map(phase_map3)
    
    thresholds3 = []
    for i in range(len(cluster_means3)-1):
        thresh = (cluster_means3.iloc[i] + cluster_means3.iloc[i+1]) / 2
        thresholds3.append(thresh)
    
    cv_n3 = fit_partition_heat_capacity(df_n3, Cv_support)
    comparison['n3'] = {
        'method': '三相自动分区',
        'n_partitions': 3,
        'thresholds': thresholds3,
        'cv_results': cv_n3,
        'df_clustered': df_n3,
        'metrics': {
            'silhouette': silhouette_n3,
            'calinski_harabasz': calinski_n3,
            'davies_bouldin': davies_n3
        }
    }
    
    # ========================================================================
    # 保存每种策略的质量指标到CSV文件
    # ========================================================================
    print("\n保存质量指标...")
    
    # 保存策略的quality_metrics
    strategy_map = {
        'fixed': ('fixed3', 3),   # 固定阈值3分区
        'fixed2': ('fixed2', 2),  # 固定阈值2分区  
        'n2': ('auto2', 2),       # 自动2分区
        'n3': ('auto3', 3)        # 自动3分区
    }
    
    for strategy_key, (file_prefix, n_parts) in strategy_map.items():
        if strategy_key in comparison and 'cv_results' in comparison[strategy_key]:
            cv_results = comparison[strategy_key]['cv_results']
            
            # 获取聚类质量指标（仅auto2和auto3有）
            metrics = comparison[strategy_key].get('metrics', {})
            silhouette = metrics.get('silhouette', None)
            calinski = metrics.get('calinski_harabasz', None)
            davies = metrics.get('davies_bouldin', None)
            
            # 构建quality_metrics数据
            quality_data = []
            for partition_key, cv_data in cv_results.items():
                if partition_key in ['solid', 'premelting', 'liquid', 'partition1', 'partition2', 'partition3']:
                    row = {
                        'structure': structure_name,
                        'n_partitions': n_parts,
                        'phase': partition_key,
                        'Cv_cluster': cv_data['Cv_cluster'],
                        'Cv_cluster_err': cv_data.get('Cv_cluster_err', 0.0),
                        'R2': cv_data['R2'],
                        'R2_total': cv_data.get('R2_total', cv_data['R2']),
                        'grade': cv_data.get('grade', ''),
                        'grade_score': cv_data.get('grade_score', 0)
                    }
                    # 添加聚类质量指标（每行都重复，便于后续分析）
                    if silhouette is not None:
                        row['silhouette_score'] = silhouette
                    if calinski is not None:
                        row['calinski_harabasz'] = calinski
                    if davies is not None:
                        row['davies_bouldin'] = davies
                    quality_data.append(row)
            
            if quality_data:
                quality_df = pd.DataFrame(quality_data)
                quality_file = OUTPUT_DIR / f'{structure_name}_{file_prefix}_quality_metrics.csv'
                quality_df.to_csv(quality_file, index=False, encoding='utf-8-sig')
                print(f"  [OK] {file_prefix} 质量指标已保存: {quality_file.name}")
    
    return comparison

# ============================================================================
# 6. 可视化函数
# ============================================================================

def plot_clustering_results(results, df_structure, output_dir):
    """
    生成聚类分析的完整可视化
    
    包含:
    1. Elbow plot (如果auto_k)
    2. 聚类结果 vs 固定阈值对比
    3. 温度-Lindemann散点图
    4. 层次聚类树状图 (如果方法是hierarchical)
    """
    structure_name = results['structure']
    method = results['method']
    
    # 决定子图布局（新增图j和图k后需要额外一行）
    has_elbow = results['k_analysis'] is not None
    
    if has_elbow:
        # 有K分析: 5行布局 (第1行Elbow图, 2-5行其他图)
        fig = plt.figure(figsize=(24, 24))
        gs = fig.add_gridspec(5, 3, hspace=0.35, wspace=0.3, 
                             height_ratios=[1, 1, 1, 1, 1])
        row_offset = 1
    else:
        # 无K分析: 4行布局（新增了图j和图k）
        fig = plt.figure(figsize=(24, 18))
        gs = fig.add_gridspec(4, 3, hspace=0.35, wspace=0.3)
        row_offset = 0
    
    n_partitions = results['n_partitions']
    fig.suptitle(f'Lindemann指数相态分区分析 - {structure_name} ({method.upper()}, n={n_partitions}分区)', 
                 fontsize=18, fontweight='bold', y=0.995)
    
    colors_phase = {
        'partition1': '#3498db',  # 蓝色 - 分区1
        'partition2': '#e67e22',  # 橙色 - 分区2
        'partition3': '#e74c3c',  # 红色 - 分区3
        'uniform': '#95a5a6'      # 灰色 - 单一相
    }
    
    # ========== 子图1: Elbow Plot (如果有K分析,占第一行全部) ==========
    if has_elbow:
        ax_elbow = fig.add_subplot(gs[0, :])
        k_analysis = results['k_analysis']
        
        ax_elbow_twin1 = ax_elbow.twinx()
        ax_elbow_twin2 = ax_elbow.twinx()
        ax_elbow_twin2.spines['right'].set_position(('outward', 60))
        
        # Inertia
        line1 = ax_elbow.plot(k_analysis['k_range'], k_analysis['inertias'], 
                              'o-', color='blue', linewidth=2, markersize=8, label='Inertia')
        ax_elbow.set_xlabel('分区数 (n)', fontsize=12, fontweight='bold')
        ax_elbow.set_ylabel('Inertia', fontsize=12, color='blue', fontweight='bold')
        ax_elbow.tick_params(axis='y', labelcolor='blue')
        
        # Silhouette
        line2 = ax_elbow_twin1.plot(k_analysis['k_range'], k_analysis['silhouette_scores'], 
                                     's-', color='green', linewidth=2, markersize=8, label='Silhouette')
        ax_elbow_twin1.set_ylabel('Silhouette Score', fontsize=12, color='green', fontweight='bold')
        ax_elbow_twin1.tick_params(axis='y', labelcolor='green')
        
        # Calinski-Harabasz
        line3 = ax_elbow_twin2.plot(k_analysis['k_range'], k_analysis['calinski_scores'], 
                                     '^-', color='red', linewidth=2, markersize=8, label='Calinski-Harabasz')
        ax_elbow_twin2.set_ylabel('Calinski-Harabasz Index', fontsize=12, color='red', fontweight='bold')
        ax_elbow_twin2.tick_params(axis='y', labelcolor='red')
        
        # 标记最优K
        optimal_k = k_analysis['optimal_k']
        ax_elbow.axvline(x=optimal_k, color='gray', linestyle='--', linewidth=2, alpha=0.7)
        ax_elbow.text(optimal_k, ax_elbow.get_ylim()[1]*0.9, f'Optimal K={optimal_k}', 
                     ha='center', fontsize=11, fontweight='bold',
                     bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.7))
        
        ax_elbow.set_title('(a) Elbow Method - 最优K值确定', fontsize=12, fontweight='bold')
        ax_elbow.grid(True, alpha=0.3)
        
        # 合并图例
        lines = line1 + line2 + line3
        labels = [l.get_label() for l in lines]
        ax_elbow.legend(lines, labels, loc='upper right', fontsize=10)
    
    # ========== 子图2: 分区结果 (温度-Lindemann) ==========
    ax1 = fig.add_subplot(gs[row_offset, 0])
    df_clustered = results['df_clustered']
    
    for phase, color in colors_phase.items():
        df_phase = df_clustered[df_clustered['phase_clustered'] == phase]
        if len(df_phase) > 0:
            ax1.scatter(df_phase['temp'], df_phase['delta'], 
                       c=color, alpha=0.6, s=80, label=f'{phase} (n={len(df_phase)})',
                       edgecolors='black', linewidths=0.5)
    
    # 添加聚类阈值线
    for i, thresh in enumerate(results['thresholds']):
        ax1.axhline(y=thresh, color='red', linestyle='--', linewidth=2, 
                   label=f'Threshold {i+1}: δ={thresh:.4f}', alpha=0.7)
    
    ax1.set_xlabel('Temperature (K)', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Lindemann Index δ', fontsize=11, fontweight='bold')
    ax1.set_title(f'(b) 分区结果 (n={results["n_partitions"]})', 
                  fontsize=11, fontweight='bold')
    ax1.legend(fontsize=9, loc='best')
    ax1.grid(True, alpha=0.3)
    
    # ========== 子图3: 固定阈值对比 ==========
    ax2 = fig.add_subplot(gs[row_offset, 1])
    
    for phase, color in colors_phase.items():
        df_phase = df_structure[df_structure['phase'] == phase]
        if len(df_phase) > 0:
            ax2.scatter(df_phase['temp'], df_phase['delta'], 
                       c=color, alpha=0.6, s=80, label=f'{phase} (n={len(df_phase)})',
                       edgecolors='black', linewidths=0.5)
    
    # 固定阈值线
    ax2.axhline(y=0.1, color='gray', linestyle='--', linewidth=2, 
               label='Fixed: δ=0.10', alpha=0.7)
    ax2.axhline(y=0.15, color='gray', linestyle='--', linewidth=2, 
               label='Fixed: δ=0.15', alpha=0.7)
    
    ax2.set_xlabel('Temperature (K)', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Lindemann Index δ', fontsize=11, fontweight='bold')
    ax2.set_title('(c) 固定阈值 (0.10, 0.15)', 
                  fontsize=11, fontweight='bold')
    ax2.legend(fontsize=9, loc='best')
    ax2.grid(True, alpha=0.3)
    
    # ========== 子图4: 簇分布箱线图 ==========
    ax3 = fig.add_subplot(gs[row_offset, 2])
    
    cluster_data = []
    cluster_labels = []
    for stat in results['cluster_stats']:
        cluster_id = stat['cluster_id']
        phase = results['cluster_to_phase'][cluster_id]
        data = df_clustered[df_clustered['cluster'] == cluster_id]['delta'].values
        cluster_data.append(data)
        cluster_labels.append(f"{phase}\n(n={len(data)})")
    
    bp = ax3.boxplot(cluster_data, labels=cluster_labels, patch_artist=True,
                     showmeans=True, meanprops=dict(marker='D', markerfacecolor='red', markersize=8))
    
    # 为每个箱线图着色
    for i, (patch, stat) in enumerate(zip(bp['boxes'], results['cluster_stats'])):
        phase = results['cluster_to_phase'][stat['cluster_id']]
        if phase in colors_phase:
            patch.set_facecolor(colors_phase[phase])
            patch.set_alpha(0.6)
    
    ax3.set_ylabel('Lindemann Index δ', fontsize=11, fontweight='bold')
    ax3.set_title('(d) 簇内Lindemann分布', fontsize=11, fontweight='bold')
    ax3.grid(axis='y', alpha=0.3)
    
    # ========== 子图5: 聚类评估指标 + 算法说明 ==========
    ax4 = fig.add_subplot(gs[row_offset+1, 0])
    
    # 添加算法选择说明
    algorithm_info = {
        'kmeans': '● K-means: 基于距离最小化\n  优势: 计算快速，结果稳定\n  适用: 球形分布数据',
        'hierarchical': '● 层次聚类: 基于相似度树\n  优势: 可视化树状图\n  适用: 嵌套结构数据',
        'dbscan': '● DBSCAN: 基于密度\n  优势: 发现任意形状簇\n  适用: 噪声数据'
    }
    
    metrics = results['metrics']
    
    # 上半部分：算法说明
    ax4.text(0.05, 0.95, f'当前算法: {method.upper()}', 
             transform=ax4.transAxes, fontsize=11, fontweight='bold',
             verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.7))
    
    ax4.text(0.05, 0.75, algorithm_info.get(method, ''), 
             transform=ax4.transAxes, fontsize=9,
             verticalalignment='top')
    
    # 下半部分：关键指标
    if metrics['silhouette'] is not None:
        silh = metrics['silhouette']
        if silh > 0.70:
            quality = "优秀 (Excellent)"
            color = '#27ae60'
        elif silh > 0.50:
            quality = "良好 (Good)"
            color = '#3498db'
        elif silh > 0.25:
            quality = "可接受 (Fair)"
            color = '#e67e22'
        else:
            quality = "较差 (Poor)"
            color = '#e74c3c'
        
        ax4.text(0.05, 0.40, f'Silhouette分数: {silh:.4f}', 
                transform=ax4.transAxes, fontsize=10, fontweight='bold',
                verticalalignment='top')
        ax4.text(0.05, 0.30, f'质量评级: {quality}', 
                transform=ax4.transAxes, fontsize=10,
                verticalalignment='top', color=color, fontweight='bold')
        
        # Calinski-Harabasz指标
        if metrics['calinski_harabasz'] is not None:
            ax4.text(0.05, 0.18, f'Calinski-Harabasz: {metrics["calinski_harabasz"]:.2f}', 
                    transform=ax4.transAxes, fontsize=9,
                    verticalalignment='top')
        
        # Davies-Bouldin指标
        if metrics['davies_bouldin'] is not None:
            ax4.text(0.05, 0.08, f'Davies-Bouldin: {metrics["davies_bouldin"]:.4f} (越低越好)', 
                    transform=ax4.transAxes, fontsize=9,
                    verticalalignment='top')
    else:
        ax4.text(0.5, 0.5, 'No metrics available', 
                ha='center', va='center', fontsize=12,
                transform=ax4.transAxes)
    
    ax4.set_title('(e) 聚类质量评估', fontsize=11, fontweight='bold')
    ax4.axis('off')
    
    # ========== 子图6: 各温度点的相态分布 (堆叠柱状图,step7.3风格) ==========
    ax5 = fig.add_subplot(gs[row_offset+1, 1])
    
    # 统计每个温度点各相态的数量
    temp_phase = df_clustered.groupby(['temp', 'phase_clustered']).size().unstack(fill_value=0)
    temp_sorted = sorted(df_clustered['temp'].unique())
    
    # 准备堆叠数据（使用partition1/2/3）
    partition1_counts = [temp_phase.loc[t, 'partition1'] if t in temp_phase.index and 'partition1' in temp_phase.columns else 0 
                         for t in temp_sorted]
    partition2_counts = [temp_phase.loc[t, 'partition2'] if t in temp_phase.index and 'partition2' in temp_phase.columns else 0 
                         for t in temp_sorted]
    partition3_counts = [temp_phase.loc[t, 'partition3'] if t in temp_phase.index and 'partition3' in temp_phase.columns else 0 
                         for t in temp_sorted]
    
    x_pos = np.arange(len(temp_sorted))
    
    # 堆叠柱状图
    ax5.bar(x_pos, partition1_counts, label='Partition 1', 
            color=colors_phase.get('partition1', '#3498db'), alpha=0.8, edgecolor='black', linewidth=0.5)
    ax5.bar(x_pos, partition2_counts, bottom=partition1_counts, label='Partition 2', 
            color=colors_phase.get('partition2', '#e67e22'), alpha=0.8, edgecolor='black', linewidth=0.5)
    ax5.bar(x_pos, partition3_counts, bottom=np.array(partition1_counts)+np.array(partition2_counts), 
            label='Partition 3', color=colors_phase.get('partition3', '#e74c3c'), alpha=0.8, edgecolor='black', linewidth=0.5)
    
    ax5.set_xlabel('Temperature (K)', fontsize=11, fontweight='bold')
    ax5.set_ylabel('Number of Data Points', fontsize=11, fontweight='bold')
    ax5.set_title('(f) 各温度点的相态分布\nPhase Distribution vs Temperature', fontsize=11, fontweight='bold')
    
    # 只显示部分温度标签（避免拥挤）
    step = max(1, len(temp_sorted) // 10)
    ax5.set_xticks(x_pos[::step])
    ax5.set_xticklabels([f'{int(t)}' for t in temp_sorted[::step]], rotation=45)
    ax5.legend(fontsize=9, loc='upper left', framealpha=0.95)
    ax5.grid(True, alpha=0.3, linestyle=':', axis='y', linewidth=0.8)
    
    # ========== 子图7: Lindemann-温度散点 (按相态分组) ==========
    ax6 = fig.add_subplot(gs[row_offset+1, 2])
    
    for phase, color in colors_phase.items():
        df_phase = df_clustered[df_clustered['phase_clustered'] == phase]
        if len(df_phase) > 0:
            ax6.scatter(df_phase['delta'], df_phase['temp'], 
                       c=color, alpha=0.6, s=60, label=f'{phase} (n={len(df_phase)})',
                       edgecolors='black', linewidths=0.5)
    
    ax6.set_xlabel('Lindemann Index δ', fontsize=11, fontweight='bold')
    ax6.set_ylabel('Temperature (K)', fontsize=11, fontweight='bold')
    ax6.set_title('(g) δ-T关系 (按相态)', fontsize=11, fontweight='bold')
    ax6.legend(fontsize=9, loc='best')
    ax6.grid(True, alpha=0.3)
    
    # ========== 子图8-9: 分区热容拟合图 ==========
    # 拟合各分区的热容
    Cv_support = 38.2151  # meV/K (来自step7.4)
    cv_results = fit_partition_heat_capacity(df_clustered, Cv_support)
    
    if cv_results:
        # 子图8: 相对能量-温度拟合线图 (按分区) - 只显示团簇能量变化
        ax7 = fig.add_subplot(gs[row_offset+2, 0])
        
        # 加载载体能量数据
        support_fit = load_support_energy_data()
        
        if support_fit is not None:
            # 使用实际的载体能量拟合结果
            slope_support, intercept_support, R2_support = support_fit
            print(f"  [图h] 使用实际载体能量数据: Cv={slope_support*1000:.4f} meV/K, R^2={R2_support:.6f}")
        else:
            # 回退：使用默认Cv_support估算
            print(f"  [图h] 使用默认Cv_support估算载体能量")
            slope_support = Cv_support / 1000  # meV/K -> eV/K
            T_min = df_clustered['temp'].min()
            E_total_min = df_clustered[df_clustered['temp'] == T_min]['avg_energy'].mean()
            intercept_support = E_total_min * 0.9 - slope_support * T_min
        
        # 计算参考能量（使用最低温度点的团簇能量）
        T_min = df_clustered['temp'].min()
        E_support_ref = slope_support * T_min + intercept_support
        E_total_ref = df_clustered['avg_energy'].min()
        E_cluster_ref = E_total_ref - E_support_ref
        
        # 重新拟合团簇能量，得到新的R²和Cv_cluster
        cv_cluster_results = {}
        
        # 绘制每个分区的数据点和拟合线
        for phase, color in colors_phase.items():
            if phase not in cv_results:
                continue
                
            cv_data = cv_results[phase]
            df_phase = df_clustered[df_clustered['phase_clustered'] == phase].copy()
            
            # 计算团簇的绝对能量和相对能量
            E_support_phase = slope_support * df_phase['temp'] + intercept_support
            E_cluster_abs = df_phase['avg_energy'] - E_support_phase
            df_phase['cluster_energy'] = E_cluster_abs
            df_phase['rel_cluster_energy'] = E_cluster_abs - E_cluster_ref
            
            # 重新拟合团簇能量-温度关系
            T_data = df_phase['temp'].values
            E_cluster_data = df_phase['cluster_energy'].values
            
            slope_cluster, intercept_cluster, r_value_cluster, p_value, std_err_cluster = linregress(T_data, E_cluster_data)
            R2_cluster = r_value_cluster ** 2
            Cv_cluster_new = slope_cluster * 1000  # eV/K -> meV/K
            Cv_cluster_err_new = std_err_cluster * 1000
            
            # 保存重新拟合的结果
            cv_cluster_results[phase] = {
                'Cv_cluster': Cv_cluster_new,
                'Cv_cluster_err': Cv_cluster_err_new,
                'R2': R2_cluster,
                'slope': slope_cluster,
                'intercept': intercept_cluster,
                'n_points': len(df_phase)
            }
            
            # 输出重新拟合的结果
            print(f"  [重新拟合] {phase}: Cv_cluster={Cv_cluster_new:.4f}±{Cv_cluster_err_new:.4f} meV/K, R^2={R2_cluster:.4f}")
            
            # 数据点
            ax7.scatter(df_phase['temp'], df_phase['rel_cluster_energy'], 
                       c=color, alpha=0.6, s=60, label=f'{phase} data',
                       edgecolors='black', linewidths=0.5)
            
            # 重新拟合的拟合线（使用团簇能量）
            T_fit = np.linspace(T_data.min(), T_data.max(), 100)
            E_fit_cluster = slope_cluster * T_fit + intercept_cluster
            E_fit_rel = E_fit_cluster - E_cluster_ref
            ax7.plot(T_fit, E_fit_rel, color=color, linewidth=2.5, alpha=0.8,
                    label=f'{phase} fit (R^2={R2_cluster:.4f})')
        
        ax7.set_xlabel('Temperature (K)', fontsize=11, fontweight='bold')
        ax7.set_ylabel('ΔE_cluster (eV)', fontsize=11, fontweight='bold')
        ax7.set_title('(h) 团簇相对能量-温度拟合 (各分区)\nCluster Relative Energy-Temperature Fit', 
                     fontsize=11, fontweight='bold')
        ax7.legend(fontsize=8, loc='best')
        ax7.grid(True, alpha=0.3)
        
        # 添加说明：已扣除载体能量
        ax7.text(0.02, 0.98, f'E_cluster,ref = {E_cluster_ref:.2f} eV\n(载体能量已扣除)',
                transform=ax7.transAxes, fontsize=8, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.7))
        
        # 子图9: 各分区团簇热容柱状图（使用重新拟合的结果）
        ax8 = fig.add_subplot(gs[row_offset+2, 1])
        
        phases_with_cv = []
        cv_values = []
        cv_errors = []
        colors_list = []
        
        # 使用cv_cluster_results（重新拟合团簇能量得到的结果）
        # 根据实际存在的相态绘图，而不是硬编码三个相态
        for phase in cv_cluster_results.keys():
            phases_with_cv.append(phase)
            cv_values.append(cv_cluster_results[phase]['Cv_cluster'])
            cv_errors.append(cv_cluster_results[phase]['Cv_cluster_err'])
            colors_list.append(colors_phase.get(phase, '#95a5a6'))
        
        if phases_with_cv:
            x_pos = np.arange(len(phases_with_cv))
            bars = ax8.bar(x_pos, cv_values, yerr=cv_errors, 
                          color=colors_list, alpha=0.7, capsize=5,
                          edgecolor='black', linewidth=1.5)
            
            # 计算y轴范围，为文本留出足够空间
            max_val = max([cv_values[i] + cv_errors[i] for i in range(len(cv_values))])
            
            # 添加数值标签（使用重新拟合的结果）
            # 根据分区数调整文本大小和偏移
            if len(phases_with_cv) >= 3:
                text_fontsize = 7  # 3个分区时减小字号
                text_offset = max_val * 0.05  # 使用相对偏移
            else:
                text_fontsize = 8
                text_offset = max_val * 0.03
            
            for i, (phase, val) in enumerate(zip(phases_with_cv, cv_values)):
                n_points = cv_cluster_results[phase]['n_points']
                error = cv_errors[i]
                r2 = cv_cluster_results[phase]['R2']  # 使用重新拟合的R²
                
                # 质量评级
                if r2 > 0.95:
                    grade = "★★★"
                elif r2 > 0.90:
                    grade = "★★"
                elif r2 > 0.80:
                    grade = "★"
                else:
                    grade = "⚠"
                
                # 在柱子上方显示详细信息（标注为重新拟合结果）
                ax8.text(i, val + error + text_offset, 
                        f'{val:.2f}±{error:.2f}\n(n={n_points}, $R^2$={r2:.3f} {grade})',
                        ha='center', va='bottom', fontsize=text_fontsize, fontweight='bold')
            
            # 自动调整y轴上限，为文本留出空间
            y_max = max_val * 1.25  # 留出25%的空间
            ax8.set_ylim(0, y_max)
            
            ax8.set_xticks(x_pos)
            ax8.set_xticklabels(phases_with_cv)
            ax8.set_ylabel('Cv_cluster (meV/K)', fontsize=11, fontweight='bold')
            ax8.set_title('(i) 各分区团簇热容', fontsize=11, fontweight='bold')
            # 移除Cv_support参考线
            ax8.grid(axis='y', alpha=0.3)
        else:
            ax8.text(0.5, 0.5, 'No valid fits', ha='center', va='center', fontsize=12)
            ax8.axis('off')
        
        # ===== 新增：子图10 (j): 温度-团簇能量散点图（按温度平均，显示标准差） =====
        ax9 = fig.add_subplot(gs[row_offset+3, 0])
        
        # 加载载体能量数据
        support_fit = load_support_energy_data()
        if support_fit is not None:
            slope_support, intercept_support, R2_support = support_fit
            
            # 按温度分组，计算每个温度的团簇能量平均值和标准差
            temp_groups = df_clustered.groupby('temp')
            temps_unique = []
            E_cluster_mean = []
            E_cluster_std = []
            
            for temp, group in temp_groups:
                E_support = slope_support * temp + intercept_support
                E_cluster = group['avg_energy'].values - E_support
                temps_unique.append(temp)
                E_cluster_mean.append(np.mean(E_cluster))
                E_cluster_std.append(np.std(E_cluster))
            
            temps_unique = np.array(temps_unique)
            E_cluster_mean = np.array(E_cluster_mean)
            E_cluster_std = np.array(E_cluster_std)
            
            # 计算相对能量（相对于最低温度）
            E_cluster_ref = E_cluster_mean.min()
            E_cluster_mean_rel = E_cluster_mean - E_cluster_ref
            
            # 绘制散点图（带误差棒）- 使用相对能量
            ax9.errorbar(temps_unique, E_cluster_mean_rel, yerr=E_cluster_std,
                        fmt='o', markersize=8, capsize=5, capthick=2,
                        color='#2c3e50', alpha=0.7, ecolor='#95a5a6',
                        label='团簇相对能量 (按温度平均)')
            
            # 用颜色标注不同相态的数据点
            for phase, color in colors_phase.items():
                if phase not in cv_cluster_results:
                    continue
                df_phase = df_clustered[df_clustered['phase_clustered'] == phase]
                temps_phase = df_phase['temp'].unique()
                
                # 只标注相态区域的点（不重新绘制）
                for temp in temps_phase:
                    idx = np.where(temps_unique == temp)[0]
                    if len(idx) > 0:
                        ax9.scatter(temp, E_cluster_mean_rel[idx[0]], 
                                  c=color, s=100, alpha=0.5, zorder=3,
                                  edgecolors='black', linewidths=1)
            
            ax9.set_xlabel('Temperature (K)', fontsize=11, fontweight='bold')
            ax9.set_ylabel('ΔE_cluster (eV)', fontsize=11, fontweight='bold')
            ax9.set_title('(j) 团簇相对能量-温度关系 (按温度平均)\nCluster Relative Energy vs Temperature (per-T average)', 
                         fontsize=11, fontweight='bold')
            ax9.grid(True, alpha=0.3)
            ax9.legend(fontsize=8, loc='best')
            
            # 添加说明
            ax9.text(0.02, 0.98, f'E_cluster,ref = {E_cluster_ref:.2f} eV\n误差棒: 标准差\n彩色点: 相态标识',
                    transform=ax9.transAxes, fontsize=8, verticalalignment='top',
                    bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.7))
        else:
            ax9.text(0.5, 0.5, 'Support energy data not available', 
                    ha='center', va='center', fontsize=12)
            ax9.axis('off')
        
        # ===== 新增：子图11 (k): 整体拟合热容（基准对比） =====
        ax10 = fig.add_subplot(gs[row_offset+3, 1])
        
        if support_fit is not None:
            # 对所有温度点的团簇相对能量整体拟合（不分相态）
            T_all = temps_unique
            E_cluster_all_rel = E_cluster_mean_rel
            
            if len(T_all) >= 3:
                slope_overall, intercept_overall, r_value_overall, p_value_overall, std_err_overall = linregress(T_all, E_cluster_all_rel)
                R2_overall = r_value_overall ** 2
                Cv_overall = slope_overall * 1000  # meV/K
                Cv_overall_err = std_err_overall * 1000
                
                # 绘制整体拟合线（不绘制数据点，数据点在分相态中绘制）
                T_fit = np.linspace(T_all.min(), T_all.max(), 100)
                E_fit_rel = slope_overall * T_fit + intercept_overall
                
                ax10.plot(T_fit, E_fit_rel, 'r--', linewidth=2.5, alpha=0.8,
                         label=f'整体拟合 (R^2={R2_overall:.4f})', zorder=2)
                
                # 基于聚类结果，对每个分区使用多数投票确定专属温度（避免重叠）
                phase_avg_fits = {}
                
                # 第一步：为每个温度确定主导分区（多数投票）
                temp_to_partition = {}
                print(f"\n>>> 多数投票温度分配:")
                for temp in temps_unique:
                    df_temp = df_clustered[df_clustered['temp'] == temp]
                    # 统计每个分区在该温度的数量
                    partition_counts = df_temp['phase_clustered'].value_counts()
                    # 选择最多的分区作为该温度的主导分区
                    dominant_partition = partition_counts.idxmax()
                    temp_to_partition[temp] = dominant_partition
                    print(f"  T={temp:4.0f}K: {dict(partition_counts)} → 分配给 {dominant_partition}")
                
                # 统计每个分区的温度
                print(f"\n>>> 各分区专属温度:")
                for phase in colors_phase.keys():
                    phase_temps = [temp for temp, partition in temp_to_partition.items() if partition == phase]
                    phase_temps = sorted(phase_temps)
                    if phase_temps:
                        print(f"  {phase}: {len(phase_temps)}个温度点 → {phase_temps[0]:.0f}-{phase_temps[-1]:.0f}K")
                    else:
                        print(f"  {phase}: 无专属温度（被其他分区占据）")
                
                # 第二步：为每个分区收集其专属温度点（无重叠）
                for phase, color in colors_phase.items():
                    if phase not in cv_cluster_results:
                        continue
                    
                    # 获取该分区的专属温度（多数投票胜出的温度）
                    phase_temps = [temp for temp, partition in temp_to_partition.items() if partition == phase]
                    phase_temps = sorted(phase_temps)
                    
                    # 从按温度平均的数据中筛选出该分区的温度点
                    # 至少需要2个点才能进行拟合（降低要求以显示所有分区）
                    if len(phase_temps) >= 2:
                        mask_phase = np.isin(temps_unique, phase_temps)
                        T_phase_avg = temps_unique[mask_phase]
                        E_phase_avg_rel = E_cluster_mean_rel[mask_phase]
                        E_phase_std = E_cluster_std[mask_phase]
                        
                        # 对该分区的按温度平均数据重新拟合（使用相对能量）
                        slope_ph_avg, intercept_ph_avg, r_value_ph_avg, _, std_err_ph_avg = linregress(T_phase_avg, E_phase_avg_rel)
                        R2_ph_avg = r_value_ph_avg ** 2
                        Cv_ph_avg = slope_ph_avg * 1000
                        Cv_ph_avg_err = std_err_ph_avg * 1000  # 拟合误差
                        
                        phase_avg_fits[phase] = {
                            'slope': slope_ph_avg,
                            'intercept': intercept_ph_avg,
                            'R2': R2_ph_avg,
                            'Cv': Cv_ph_avg,
                            'Cv_err': Cv_ph_avg_err,
                            'n_temps': len(T_phase_avg),
                            'T_range': (T_phase_avg.min(), T_phase_avg.max()),
                            'T_data': T_phase_avg,
                            'E_data': E_phase_avg_rel,
                            'E_std': E_phase_std
                        }
                        
                        # 绘制分区数据点（带误差棒）- 使用相对能量
                        ax10.errorbar(T_phase_avg, E_phase_avg_rel, yerr=E_phase_std,
                                     fmt='o', markersize=6, capsize=4, capthick=1.5,
                                     color=color, alpha=0.6, ecolor=color,
                                     label=f'{phase} 数据', zorder=4)
                        
                        # 绘制分区拟合线（相对能量）
                        T_phase_fit = np.linspace(T_phase_avg.min(), T_phase_avg.max(), 50)
                        E_phase_fit_rel = slope_ph_avg * T_phase_fit + intercept_ph_avg
                        
                        ax10.plot(T_phase_fit, E_phase_fit_rel, '-', color=color, 
                                 linewidth=2.5, alpha=0.9, 
                                 label=f'{phase} 拟合 (R^2={R2_ph_avg:.4f})', zorder=3)
                        
                        # 输出分区拟合结果
                        print(f"  [分区拟合] {phase}: Cv={Cv_ph_avg:.4f}±{Cv_ph_avg_err:.4f} meV/K, R^2={R2_ph_avg:.4f}, n_temps={len(T_phase_avg)}, T={T_phase_avg.min():.0f}-{T_phase_avg.max():.0f}K")
                
                ax10.set_xlabel('Temperature (K)', fontsize=11, fontweight='bold')
                ax10.set_ylabel('ΔE_cluster (eV)', fontsize=11, fontweight='bold')
                ax10.set_title('(k) 整体拟合 vs 分区拟合对比\nOverall Fit vs Partitioned Fit', 
                             fontsize=11, fontweight='bold')
                ax10.legend(fontsize=7, loc='best', ncol=2)
                ax10.grid(True, alpha=0.3)
                
                # 添加热容对比信息（使用按温度平均数据的拟合结果，带误差）
                cv_info = f'整体拟合:\n  Cv={Cv_overall:.3f}±{Cv_overall_err:.3f} meV/K\n  $R^2$={R2_overall:.4f}\n\n分区拟合:\n'
                for phase in phase_avg_fits.keys():
                    fit = phase_avg_fits[phase]
                    cv_info += f'  {phase}: Cv={fit["Cv"]:.3f}±{fit["Cv_err"]:.3f}\n         $R^2$={fit["R2"]:.4f}\n'
                
                ax10.text(0.02, 0.98, cv_info.strip(),
                         transform=ax10.transAxes, fontsize=7, verticalalignment='top',
                         bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.7))
                
                # 输出整体拟合结果
                print(f"  [整体拟合-按温度平均] Cv_overall={Cv_overall:.4f}±{Cv_overall_err:.4f} meV/K, R^2={R2_overall:.4f}, n_temps={len(T_all)}")
                
                # ===== 新增：子图12 (l): 分相态热容柱状图（基于按温度平均数据，模仿图i） =====
                ax11 = fig.add_subplot(gs[row_offset+3, 2])
                
                if phase_avg_fits:
                    phases_with_fit = []
                    cv_values_avg = []
                    cv_errors_avg = []
                    colors_list_avg = []
                    
                    # 根据实际存在的相态绘图
                    for phase in phase_avg_fits.keys():
                        phases_with_fit.append(phase)
                        fit = phase_avg_fits[phase]
                        cv_values_avg.append(fit['Cv'])
                        cv_errors_avg.append(fit['Cv_err'])
                        colors_list_avg.append(colors_phase.get(phase, '#95a5a6'))
                    
                    x_pos = np.arange(len(phases_with_fit))
                    bars = ax11.bar(x_pos, cv_values_avg, yerr=cv_errors_avg, 
                                   color=colors_list_avg, alpha=0.7, capsize=5,
                                   edgecolor='black', linewidth=1.5)
                    
                    # 添加数值标签（包含误差和质量评级）
                    for i, (phase, val) in enumerate(zip(phases_with_fit, cv_values_avg)):
                        fit = phase_avg_fits[phase]
                        error = cv_errors_avg[i]
                        r2 = fit['R2']
                        n_temps = fit['n_temps']
                        
                        # 质量评级（基于R²）
                        if r2 > 0.95:
                            grade = "★★★"
                        elif r2 > 0.90:
                            grade = "★★"
                        elif r2 > 0.80:
                            grade = "★"
                        else:
                            grade = "⚠"
                        
                        # 在柱子上方显示详细信息
                        ax11.text(i, val + error + 0.3, 
                                f'{val:.3f}±{error:.3f}\n(n={n_temps}, $R^2$={r2:.4f} {grade})',
                                ha='center', va='bottom', fontsize=8, fontweight='bold')
                    
                    ax11.set_xticks(x_pos)
                    ax11.set_xticklabels(phases_with_fit)
                    ax11.set_ylabel('Cv_cluster (meV/K)', fontsize=11, fontweight='bold')
                    ax11.set_title('(l) 分相态团簇热容 (按温度平均-多数投票)\nPartitioned Cv (per-T averaged, majority vote)', 
                                  fontsize=11, fontweight='bold')
                    ax11.grid(axis='y', alpha=0.3)
                    
                    # 添加说明
                    ax11.text(0.02, 0.98, '基于按温度平均数据\n多数投票规则避免交叉',
                            transform=ax11.transAxes, fontsize=8, verticalalignment='top',
                            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.7))
                else:
                    ax11.text(0.5, 0.5, 'No valid fits', ha='center', va='center', fontsize=12)
                    ax11.axis('off')
            else:
                ax10.text(0.5, 0.5, 'Not enough data points', 
                         ha='center', va='center', fontsize=12)
                ax10.axis('off')
        else:
            ax10.text(0.5, 0.5, 'Support energy data not available', 
                     ha='center', va='center', fontsize=12)
            ax10.axis('off')
    
    else:
        # 如果没有cv_results，显示空白图
        ax7 = fig.add_subplot(gs[row_offset+2, 0])
        ax7.text(0.5, 0.5, 'No valid heat capacity fits', 
                ha='center', va='center', fontsize=12)
        ax7.axis('off')
        
        ax8 = fig.add_subplot(gs[row_offset+2, 1])
        ax8.text(0.5, 0.5, 'No valid heat capacity fits', 
                ha='center', va='center', fontsize=12)
        ax8.axis('off')
    
    plt.tight_layout()
    
    # 保存图表（文件名包含n值）
    output_file = output_dir / f'{structure_name}_{method}_n{n_partitions}_partition_analysis.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n  [OK] 可视化已保存: {output_file.name}")
    plt.close()
    
    return output_file

# ============================================================================
# 7. 算法对比可视化
# ============================================================================

def plot_algorithm_comparison(algo_results, structure_name, output_dir):
    """
    生成不同聚类算法的对比图
    
    Parameters:
    -----------
    algo_results : dict
        各算法的结果，格式为 {method: {'silhouette': ..., 'result': ...}}
    structure_name : str
        结构名称
    output_dir : Path
        输出目录
    """
    print(f"\n>>> 生成算法对比可视化: {structure_name}...")
    
    fig = plt.figure(figsize=(18, 10))
    gs = fig.add_gridspec(2, 3, hspace=0.35, wspace=0.3)
    fig.suptitle(f'聚类算法对比 - {structure_name}', fontsize=18, fontweight='bold')
    
    colors_phase = {
        'partition1': '#3498db',  # 蓝色 - 分区1
        'partition2': '#e67e22',  # 橙色 - 分区2
        'partition3': '#e74c3c',  # 红色 - 分区3
    }
    algo_colors = {'kmeans': '#3498db', 'hierarchical': '#2ecc71', 'dbscan': '#9b59b6'}
    
    # 子图1-3: 各算法的聚类结果T-δ图
    for idx, (method, data) in enumerate(sorted(algo_results.items())):
        ax = fig.add_subplot(gs[0, idx])
        result = data['result']
        df_clustered = result['df_clustered']
        
        # 动态获取所有分区
        partitions = sorted(df_clustered['phase_clustered'].unique())
        for partition in partitions:
            df_partition = df_clustered[df_clustered['phase_clustered'] == partition]
            if len(df_partition) > 0:
                ax.scatter(df_partition['temp'], df_partition['delta'], 
                         c=colors_phase.get(partition, '#95a5a6'),
                         label=f'{partition} (n={len(df_partition)})',
                         alpha=0.7, s=50, edgecolors='black', linewidth=0.8)
        
        ax.set_title(f'({chr(97+idx)}) {method.upper()}\nSilh={data["silhouette"]:.4f}', 
                    fontsize=11, fontweight='bold',
                    bbox=dict(boxstyle='round', facecolor=algo_colors[method], alpha=0.3))
        ax.set_xlabel('Temperature (K)', fontsize=10, fontweight='bold')
        ax.set_ylabel('Lindemann Index δ', fontsize=10, fontweight='bold')
        ax.legend(fontsize=8, loc='best')
        ax.grid(True, alpha=0.3)
    
    # 子图4: 指标对比雷达图
    ax4 = fig.add_subplot(gs[1, 0], projection='polar')
    
    methods = sorted(algo_results.keys())
    metrics = ['Silhouette', 'Calinski-H\n(norm)', '1-Davies-B']
    angles = np.linspace(0, 2 * np.pi, len(metrics), endpoint=False).tolist()
    angles += angles[:1]  # 闭合
    
    for method in methods:
        data = algo_results[method]
        silh = data['silhouette']
        cal = min(data['calinski'] / 100, 1.0) if data['calinski'] else 0
        dav = max(0, 1 - data['davies']) if data['davies'] else 0
        
        values = [silh, cal, dav]
        values += values[:1]  # 闭合
        
        ax4.plot(angles, values, 'o-', linewidth=2, label=method.upper(),
                color=algo_colors[method], markersize=8)
        ax4.fill(angles, values, alpha=0.15, color=algo_colors[method])
    
    ax4.set_xticks(angles[:-1])
    ax4.set_xticklabels(metrics, fontsize=10)
    ax4.set_ylim(0, 1)
    ax4.set_title('(d) 聚类质量指标对比\n(雷达图)', fontsize=11, fontweight='bold', pad=20)
    ax4.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1), fontsize=9)
    ax4.grid(True)
    
    # 子图5: Silhouette分数柱状图
    ax5 = fig.add_subplot(gs[1, 1])
    
    methods_sorted = sorted(methods, key=lambda m: algo_results[m]['silhouette'], reverse=True)
    silh_scores = [algo_results[m]['silhouette'] for m in methods_sorted]
    colors_bars = [algo_colors[m] for m in methods_sorted]
    
    bars = ax5.bar(range(len(methods_sorted)), silh_scores, color=colors_bars, alpha=0.7,
                   edgecolor='black', linewidth=1.5)
    
    # 添加数值标签和评级
    for i, (method, score) in enumerate(zip(methods_sorted, silh_scores)):
        if score > 0.70:
            grade = "优秀"
            color = '#27ae60'
        elif score > 0.50:
            grade = "良好"
            color = '#3498db'
        elif score > 0.25:
            grade = "可接受"
            color = '#e67e22'
        else:
            grade = "较差"
            color = '#e74c3c'
        
        ax5.text(i, score + 0.02, f'{score:.4f}\n({grade})',
                ha='center', va='bottom', fontsize=10, fontweight='bold', color=color)
    
    ax5.set_xticks(range(len(methods_sorted)))
    ax5.set_xticklabels([m.upper() for m in methods_sorted], fontsize=10, fontweight='bold')
    ax5.set_ylabel('Silhouette Score', fontsize=11, fontweight='bold')
    ax5.set_title('(e) Silhouette分数对比\n(越高越好)', fontsize=11, fontweight='bold')
    ax5.set_ylim(0, 1.0)
    ax5.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5, label='良好阈值')
    ax5.axhline(y=0.7, color='green', linestyle='--', alpha=0.5, label='优秀阈值')
    ax5.legend(fontsize=8, loc='upper right')
    ax5.grid(axis='y', alpha=0.3)
    
    # 子图6: 算法特性说明
    ax6 = fig.add_subplot(gs[1, 2])
    ax6.axis('off')
    
    # 算法特性表格（使用支持中文的字体）
    algo_info = """算法特性对比:

K-means:
• 原理: 基于距离最小化
• 优势: 计算快速，结果稳定
• 适用: 球形、大小均匀的簇
• 需指定: K值(簇数量)

Hierarchical (层次聚类):
• 原理: 基于相似度树形结构
• 优势: 可视化树状图，无需指定K
• 适用: 嵌套结构数据
• 缺点: 计算开销较大

DBSCAN (密度聚类):
• 原理: 基于密度连接
• 优势: 发现任意形状簇，识别噪声
• 适用: 非凸形状、含噪声数据
• 需指定: eps (邻域半径), min_samples
    """
    
    # 显示最优选择
    best_method = max(methods, key=lambda m: algo_results[m]['silhouette'])
    best_score = algo_results[best_method]['silhouette']
    
    ax6.text(0.05, 0.95, algo_info, transform=ax6.transAxes,
            fontsize=9, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.3))
    
    ax6.text(0.5, 0.02, f'✓ 推荐算法: {best_method.upper()}\n(Silhouette = {best_score:.4f})',
            transform=ax6.transAxes, fontsize=12, fontweight='bold',
            ha='center', va='bottom',
            bbox=dict(boxstyle='round', facecolor=algo_colors[best_method], alpha=0.5))
    
    ax6.set_title('(f) 算法选择建议', fontsize=11, fontweight='bold')
    
    plt.tight_layout()
    
    # 保存图表
    output_file = output_dir / f'{structure_name}_algorithm_comparison.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  [OK] 算法对比图已保存: {output_file.name}")
    plt.close()
    
    return output_file

# ============================================================================
# 8. 多策略对比可视化
# ============================================================================

def plot_cv_comparison_standalone(comparison, structure_name, output_dir):
    """
    单独绘制热容对比图（无重叠，清晰显示）
    """
    print(f"  >>> 生成独立热容对比图...")
    
    # 颜色方案
    colors_phase = {
        'partition1': '#3498db',
        'partition2': '#e67e22',
        'partition3': '#e74c3c'
    }
    
    phase_labels = {
        'partition1': '分区1 Partition 1',
        'partition2': '分区2 Partition 2',
        'partition3': '分区3 Partition 3'
    }
    
    # 按策略顺序: Fixed 2P -> Fixed 3P -> Auto 2P -> Auto 3P
    strategy_order = [
        ('fixed2', 'Fixed 2P\n(δ=0.08)'),
        ('fixed', 'Fixed 3P\n(δ=0.10,0.15)'),
        ('n2', 'Auto 2P\n(K-means)'),
        ('n3', 'Auto 3P\n(K-means)')
    ]
    phases = ['partition1', 'partition2', 'partition3']
    
    # 收集数据（按策略分组）
    regions = []
    cv_values = []
    cv_errors = []
    r2_values = []
    point_counts = []
    bar_colors = []
    strategy_boundaries = []  # 记录策略边界位置
    
    current_pos = 0
    for key, strategy_label in strategy_order:
        if key not in comparison:
            continue
        
        data = comparison[key]
        cv_res = data['cv_results']
        strategy_start = current_pos
        
        for phase in phases:
            if phase in cv_res:
                regions.append(f'{phase_labels[phase].split()[0]}')
                cv_values.append(cv_res[phase]['Cv_cluster'])
                cv_errors.append(cv_res[phase].get('Cv_cluster_err', 0.0))
                r2_values.append(cv_res[phase]['R2'])
                point_counts.append(cv_res[phase]['n_points'])
                bar_colors.append(colors_phase[phase])
                current_pos += 1
        
        # 记录策略的中心位置和标签
        strategy_center = (strategy_start + current_pos - 1) / 2
        strategy_boundaries.append((strategy_center, strategy_label))
    
    # 创建大图
    fig, ax = plt.subplots(figsize=(20, 10))
    
    x = np.arange(len(regions))
    width = 0.7
    
    bars = ax.bar(x, cv_values, width,
                  color=bar_colors,
                  alpha=0.85, edgecolor='black', linewidth=1.8,
                  yerr=cv_errors, capsize=10, error_kw={'linewidth': 2.5, 'elinewidth': 2})
    
    # 标注数值（大字体，清晰）- 使用数学文本格式显示R²
    for i, (bar, cv, cv_err, r2, n) in enumerate(zip(bars, cv_values, cv_errors, r2_values, point_counts)):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + cv_err + 0.25,
                f'{cv:.3f}±{cv_err:.3f} meV/K\n$R^2$={r2:.4f}\n(n={n})',
                ha='center', va='bottom', fontsize=11, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.4', facecolor='white', alpha=0.9, 
                         edgecolor='gray', linewidth=1.5))
    
    # X轴：两层标签
    ax.set_xlabel('分区 Partition', fontsize=14, fontweight='bold')
    ax.set_ylabel('团簇热容 Cluster Heat Capacity Cv (meV/K)', fontsize=16, fontweight='bold')
    ax.set_title(f'{structure_name} - 各策略热容结果对比\nHeat Capacity Comparison Across All Strategies',
                fontsize=18, fontweight='bold', pad=20)
    ax.set_xticks(x)
    ax.set_xticklabels(regions, fontsize=10, fontweight='bold', ha='center')
    ax.grid(True, alpha=0.3, linestyle=':', axis='y', linewidth=1.2)
    ax.set_ylim(0, max(cv_values) * 1.7)
    
    # 添加策略标签（第二层X轴）
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    strategy_positions = [pos for pos, _ in strategy_boundaries]
    strategy_labels = [label for _, label in strategy_boundaries]
    ax2.set_xticks(strategy_positions)
    ax2.set_xticklabels(strategy_labels, fontsize=13, fontweight='bold', color='darkblue')
    ax2.set_xlabel('分区策略 Strategy', fontsize=14, fontweight='bold', color='darkblue')
    
    # 不添加分隔线（用户要求去掉红色竖线）
    
    # 添加图例说明相态颜色
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=colors_phase['partition1'], edgecolor='black', label=phase_labels['partition1']),
        Patch(facecolor=colors_phase['partition2'], edgecolor='black', label=phase_labels['partition2']),
        Patch(facecolor=colors_phase['partition3'], edgecolor='black', label=phase_labels['partition3'])
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=13, framealpha=0.95,
             title='相态颜色编码', title_fontsize=14)
    
    plt.tight_layout()
    
    output_file = output_dir / f'{structure_name}_cv_comparison_detailed.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"    [OK] 独立热容对比图: {output_file.name}")
    plt.close()
    
    return output_file


def plot_quality_metrics_standalone(comparison, structure_name, output_dir):
    """
    单独绘制质量指标对比图（优化布局，避免重叠）
    子图1: 所有策略的Silhouette Score
    子图2: 所有策略的Average R² Score
    子图3: Auto方法的双y轴图（聚类质量+拟合质量）
    """
    print(f"  >>> 生成独立质量指标图...")
    
    # 按策略顺序: Fixed 2P -> Fixed 3P -> Auto 2P -> Auto 3P
    strategy_order = [
        ('fixed2', 'Fixed 2P'),
        ('fixed', 'Fixed 3P'),
        ('n2', 'Auto 2P'),
        ('n3', 'Auto 3P')
    ]
    phases = ['partition1', 'partition2', 'partition3']
    
    # 准备数据
    strategy_labels_short = []
    silhouette_scores = []
    avg_r2_scores = []
    
    for key, short_label in strategy_order:
        if key not in comparison:
            continue
        data = comparison[key]
        strategy_labels_short.append(short_label)
        
        silhouette_scores.append(data.get('metrics', {}).get('silhouette', 0))
        
        cv_res = data['cv_results']
        r2_list = [cv_res[p]['R2'] for p in phases if p in cv_res]
        avg_r2 = np.mean(r2_list) if r2_list else 0
        avg_r2_scores.append(avg_r2)
    
    # 创建图形 - 1行3列布局
    fig = plt.figure(figsize=(24, 8))
    gs = fig.add_gridspec(1, 3, wspace=0.3)
    
    # 子图1: Silhouette Score (所有策略)
    ax1 = fig.add_subplot(gs[0, 0])
    x_pos = np.arange(len(strategy_labels_short))
    bars1 = ax1.bar(x_pos, silhouette_scores, color='#3498db', alpha=0.8,
                   edgecolor='black', linewidth=2, width=0.6)
    
    # 标注数值
    for bar, score in zip(bars1, silhouette_scores):
        height = bar.get_height()
        if score > 0:
            ax1.text(bar.get_x() + bar.get_width()/2., height + 0.02,
                    f'{score:.4f}',
                    ha='center', va='bottom', fontsize=13, fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.9))
    
    ax1.set_xlabel('分区策略 Strategy', fontsize=14, fontweight='bold')
    ax1.set_ylabel('Silhouette Score', fontsize=14, fontweight='bold')
    ax1.set_title('(a) 聚类质量 - Silhouette评分\nClustering Quality', 
                 fontsize=15, fontweight='bold', pad=15)
    ax1.set_xticks(x_pos)
    ax1.set_xticklabels(strategy_labels_short, fontsize=13, rotation=0, ha='center', fontweight='bold')
    ax1.set_ylim(0, 0.8)
    ax1.grid(True, alpha=0.3, linestyle=':', axis='y', linewidth=1)
    ax1.axhline(y=0.5, color='red', linestyle='--', linewidth=2, alpha=0.5, label='良好阈值 0.5')
    ax1.legend(fontsize=11, loc='upper left')
    
    # 子图2: Average R² Score (所有策略)
    ax2 = fig.add_subplot(gs[0, 1])
    bars2 = ax2.bar(x_pos, avg_r2_scores, color='#e74c3c', alpha=0.8,
                   edgecolor='black', linewidth=2, width=0.6)
    
    for bar, score in zip(bars2, avg_r2_scores):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height + 0.005,
                f'{score:.4f}',
                ha='center', va='bottom', fontsize=13, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.9))
    
    ax2.set_xlabel('分区策略 Strategy', fontsize=14, fontweight='bold')
    ax2.set_ylabel('Average $R^2$ Score', fontsize=14, fontweight='bold')
    ax2.set_title('(b) 拟合质量 - 平均$R^2$评分\nFitting Quality',
                 fontsize=15, fontweight='bold', pad=15)
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(strategy_labels_short, fontsize=13, rotation=0, ha='center', fontweight='bold')
    ax2.set_ylim(0.85, 1.0)
    ax2.grid(True, alpha=0.3, linestyle=':', axis='y', linewidth=1)
    ax2.axhline(y=0.95, color='green', linestyle='--', linewidth=2, alpha=0.5, label='优秀阈值 0.95')
    ax2.legend(fontsize=11, loc='lower right')
    
    # 子图3: Auto方法双y轴图（聚类质量+拟合质量）
    ax3 = fig.add_subplot(gs[0, 2])
    
    # 提取Auto方法数据
    auto_labels = []
    auto_silhouette = []
    auto_r2 = []
    
    for key, short_label in [('n2', 'Auto 2P'), ('n3', 'Auto 3P')]:
        if key in comparison:
            auto_labels.append(short_label)
            data = comparison[key]
            auto_silhouette.append(data.get('metrics', {}).get('silhouette', 0))
            
            cv_res = data['cv_results']
            r2_list = [cv_res[p]['R2'] for p in phases if p in cv_res]
            avg_r2 = np.mean(r2_list) if r2_list else 0
            auto_r2.append(avg_r2)
    
    if auto_labels:
        x_auto = np.arange(len(auto_labels))
        
        # 左y轴: Silhouette Score (蓝色)
        ax3_left = ax3
        bars_silhouette = ax3_left.bar(x_auto - 0.2, auto_silhouette, width=0.4,
                                       color='#3498db', alpha=0.8, edgecolor='black',
                                       linewidth=2, label='Silhouette Score')
        
        ax3_left.set_xlabel('Auto分区策略', fontsize=14, fontweight='bold')
        ax3_left.set_ylabel('Silhouette Score', fontsize=14, fontweight='bold', color='#3498db')
        ax3_left.tick_params(axis='y', labelcolor='#3498db')
        ax3_left.set_ylim(0, 0.8)
        ax3_left.set_xticks(x_auto)
        ax3_left.set_xticklabels(auto_labels, fontsize=13, fontweight='bold')
        ax3_left.grid(True, alpha=0.3, linestyle=':', axis='y', linewidth=1)
        
        # 右y轴: R² Score (红色)
        ax3_right = ax3_left.twinx()
        bars_r2 = ax3_right.bar(x_auto + 0.2, auto_r2, width=0.4,
                                color='#e74c3c', alpha=0.8, edgecolor='black',
                                linewidth=2, label='Average $R^2$')
        
        ax3_right.set_ylabel('Average $R^2$ Score', fontsize=14, fontweight='bold', color='#e74c3c')
        ax3_right.tick_params(axis='y', labelcolor='#e74c3c')
        ax3_right.set_ylim(0.85, 1.0)
        
        # 标注数值
        for bar, score in zip(bars_silhouette, auto_silhouette):
            height = bar.get_height()
            ax3_left.text(bar.get_x() + bar.get_width()/2., height + 0.02,
                         f'{score:.4f}',
                         ha='center', va='bottom', fontsize=12, fontweight='bold',
                         color='#3498db',
                         bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.9))
        
        for bar, score in zip(bars_r2, auto_r2):
            height = bar.get_height()
            ax3_right.text(bar.get_x() + bar.get_width()/2., height + 0.003,
                          f'{score:.4f}',
                          ha='center', va='bottom', fontsize=12, fontweight='bold',
                          color='#e74c3c',
                          bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.9))
        
        # 组合图例
        lines1, labels1 = ax3_left.get_legend_handles_labels()
        lines2, labels2 = ax3_right.get_legend_handles_labels()
        ax3_left.legend(lines1 + lines2, labels1 + labels2, fontsize=11, loc='upper left')
        
        ax3_left.set_title('(c) Auto方法双指标对比\nAuto Method: Clustering vs Fitting Quality',
                          fontsize=15, fontweight='bold', pad=15)
    else:
        ax3.text(0.5, 0.5, 'No Auto method data available',
                ha='center', va='center', fontsize=14, fontweight='bold')
        ax3.axis('off')
    
    plt.suptitle(f'{structure_name} - 策略质量评价指标\nStrategy Quality Metrics',
                fontsize=17, fontweight='bold', y=0.98)
    
    output_file = output_dir / f'{structure_name}_quality_metrics_detailed.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"    [OK] 独立质量指标图: {output_file.name}")
    plt.close()
    
    return output_file


def plot_strategy_comparison(comparison, structure_name, output_dir):
    """
    生成不同分区策略的对比图（Lindemann分布+温度分布）
    注意：热容对比和质量指标已单独生成详细PNG
    """
    print(f"  >>> 生成主对比图（Lindemann+温度分布）...")
    
    # 使用3x4布局（第一行显示说明文字，第二行Lindemann，第三行温度分布）
    fig = plt.figure(figsize=(28, 20))
    gs = fig.add_gridspec(3, 4, hspace=0.35, wspace=0.25)
    fig.suptitle(f'分区策略对比分析 - {structure_name}\nPartition Strategy Comparison Analysis', 
                 fontsize=16, fontweight='bold', y=0.995)
    
    # 颜色方案 (使用partition1/2/3)
    colors_phase = {
        'partition1': '#3498db',      # Blue
        'partition2': '#e67e22', # Orange
        'partition3': '#e74c3c'      # Red
    }
    
    phase_labels = {
        'partition1': '分区1 Partition 1',
        'partition2': '分区2 Partition 2',
        'partition3': '分区3 Partition 3'
    }
    
    strategy_keys = ['fixed', 'fixed2', 'n2', 'n3']
    
    # ===== 子图1-2 (第一行): 文字说明 =====
    ax1 = fig.add_subplot(gs[0, 0:2])
    ax1.text(0.5, 0.5, 
            '📊 热容对比详细图\n\n'
            '请查看单独生成的PNG:\n'
            f'{structure_name}_cv_comparison_detailed.png\n\n'
            '包含所有策略的热容、误差、R²和数据点数\n'
            '清晰无重叠的大尺寸展示',
            ha='center', va='center', fontsize=16, 
            transform=ax1.transAxes,
            bbox=dict(boxstyle='round,pad=1', facecolor='lightyellow', 
                     alpha=0.8, edgecolor='orange', linewidth=3))
    ax1.axis('off')
    
    ax2 = fig.add_subplot(gs[0, 2:4])
    ax2.text(0.5, 0.5,
            '📈 质量指标详细图\n\n'
            '请查看单独生成的PNG:\n'
            f'{structure_name}_quality_metrics_detailed.png\n\n'
            '包含Silhouette和R²评分的详细对比\n'
            '左右并列双图展示',
            ha='center', va='center', fontsize=16,
            transform=ax2.transAxes,
            bbox=dict(boxstyle='round,pad=1', facecolor='lightblue',
                     alpha=0.8, edgecolor='blue', linewidth=3))
    ax2.axis('off')
    
    phases = ['partition1', 'partition2', 'partition3']
    
    # ===== 第二行: Lindemann分布图（4个策略） =====
    # ===== 子图3 (第2行第1列): 固定阈值两分区 δ=0.08 =====
    ax3 = fig.add_subplot(gs[1, 0])
    
    if 'fixed2' in comparison:
        df_fixed2 = comparison['fixed2'].get('df_clustered')
        if df_fixed2 is not None:
            for phase in ['partition1', 'partition2']:
                df_phase = df_fixed2[df_fixed2['phase_clustered'] == phase]
                if len(df_phase) > 0:
                    ax3.scatter(df_phase['temp'], df_phase['delta'], 
                               c=colors_phase[phase], alpha=0.6, s=60, 
                               edgecolors='black', linewidths=0.8,
                               label=f'{phase_labels[phase]} (n={len(df_phase)})',
                               zorder=3)
            
            # δ=0.08阈值线
            ax3.axhline(y=0.08, color='red', linestyle='--', linewidth=2.5, 
                        label='阈值 δ=0.08', alpha=0.7, zorder=1)
            ax3.axhspan(0, 0.08, alpha=0.1, color=colors_phase['partition1'], zorder=0)
            ax3.axhspan(0.08, df_fixed2['delta'].max()*1.1, alpha=0.1, color=colors_phase['partition2'], zorder=0)
    
    ax3.set_xlabel('温度 Temperature (K)', fontsize=12, fontweight='bold')
    ax3.set_ylabel('Lindemann 指数 δ', fontsize=12, fontweight='bold')
    ax3.set_title('(c) 固定阈值两分区 (δ=0.08)\nFixed Threshold 2-Partition (δ=0.08)', 
                  fontsize=13, fontweight='bold', pad=10)
    ax3.legend(fontsize=9, loc='upper left', framealpha=0.95)
    ax3.grid(True, alpha=0.3, linestyle=':', linewidth=0.8)
    
    # ===== 子图4 (第2行第2列): 固定阈值三分区分布 =====
    ax4 = fig.add_subplot(gs[1, 1])
    
    if 'fixed' in comparison:
        df_fixed = comparison['fixed'].get('df_clustered')
        if df_fixed is not None:
            for phase in ['partition1', 'partition2', 'partition3']:
                df_phase = df_fixed[df_fixed['phase_clustered'] == phase]
                if len(df_phase) > 0:
                    ax4.scatter(df_phase['temp'], df_phase['delta'], 
                               c=colors_phase[phase], alpha=0.6, s=60, 
                               edgecolors='black', linewidths=0.8,
                               label=f'{phase_labels[phase]} (n={len(df_phase)})',
                               zorder=3)
            
            # 阈值线
            thresholds = comparison['fixed'].get('thresholds', [0.10, 0.15])
            ax4.axhline(y=thresholds[0], color='gray', linestyle='--', linewidth=2.5, 
                        label=f'阈值1 δ={thresholds[0]:.2f}', alpha=0.7, zorder=1)
            ax4.axhline(y=thresholds[1], color='red', linestyle='--', linewidth=2.5, 
                        label=f'阈值2 δ={thresholds[1]:.2f}', alpha=0.7, zorder=1)
            
            ax4.axhspan(0, thresholds[0], alpha=0.1, color=colors_phase['partition1'], zorder=0)
            ax4.axhspan(thresholds[0], thresholds[1], alpha=0.1, color=colors_phase['partition2'], zorder=0)
            ax4.axhspan(thresholds[1], df_fixed['delta'].max()*1.1, alpha=0.1, color=colors_phase['partition3'], zorder=0)
    
    ax4.set_xlabel('温度 Temperature (K)', fontsize=12, fontweight='bold')
    ax4.set_ylabel('Lindemann 指数 δ', fontsize=12, fontweight='bold')
    ax4.set_title('(d) 固定阈值三分区 (δ=0.10, 0.15)\nFixed Threshold 3-Partition (δ=0.10, 0.15)', 
                  fontsize=13, fontweight='bold', pad=10)
    ax4.legend(fontsize=9, loc='upper left', framealpha=0.95, ncol=2)
    ax4.grid(True, alpha=0.3, linestyle=':', linewidth=0.8)
    
    # ===== 子图5 (第2行第3列): 自动两分区分布 =====
    ax5 = fig.add_subplot(gs[1, 2])
    
    if 'n2' in comparison:
        df_n2 = comparison['n2'].get('df_clustered')
        if df_n2 is not None:
            for phase in ['partition1', 'partition2']:
                df_phase = df_n2[df_n2['phase_clustered'] == phase]
                if len(df_phase) > 0:
                    ax5.scatter(df_phase['temp'], df_phase['delta'], 
                               c=colors_phase[phase], alpha=0.6, s=60, 
                               edgecolors='black', linewidths=0.8,
                               label=f'{phase_labels[phase]} (n={len(df_phase)})',
                               zorder=3)
            
            # 聚类边界线
            thresholds = comparison['n2'].get('thresholds', [])
            if thresholds:
                ax5.axhline(y=thresholds[0], color='red', linestyle='--', linewidth=2.5, 
                            label=f'聚类边界 δ={thresholds[0]:.3f}', alpha=0.7, zorder=1)
                ax5.axhspan(0, thresholds[0], alpha=0.1, color=colors_phase['partition1'], zorder=0)
                ax5.axhspan(thresholds[0], df_n2['delta'].max()*1.1, alpha=0.1, color=colors_phase['partition2'], zorder=0)
    
    ax5.set_xlabel('温度 Temperature (K)', fontsize=12, fontweight='bold')
    ax5.set_ylabel('Lindemann 指数 δ', fontsize=12, fontweight='bold')
    ax5.set_title('(e) 自动两分区 (n=2)\nAutomatic 2-Partition (n=2)', 
                  fontsize=13, fontweight='bold', pad=10)
    ax5.legend(fontsize=9, loc='upper left', framealpha=0.95)
    ax5.grid(True, alpha=0.3, linestyle=':', linewidth=0.8)
    
    # ===== 子图6 (第2行第4列): 自动三分区分布 =====
    ax6 = fig.add_subplot(gs[1, 3])
    
    # 使用n=3策略的数据
    if 'n3' in comparison:
        df_n3 = comparison['n3'].get('df_clustered')
        if df_n3 is not None:
            for phase in ['partition1', 'partition2', 'partition3']:
                df_phase = df_n3[df_n3['phase_clustered'] == phase]
                if len(df_phase) > 0:
                    ax6.scatter(df_phase['temp'], df_phase['delta'], 
                               c=colors_phase[phase], alpha=0.6, s=60, 
                               edgecolors='black', linewidths=0.8,
                               label=f'{phase_labels[phase]} (n={len(df_phase)})',
                               zorder=3)
            
            # 聚类边界线
            thresholds = comparison['n3'].get('thresholds', [])
            if len(thresholds) >= 2:
                ax6.axhline(y=thresholds[0], color='gray', linestyle='--', linewidth=2.5, 
                            label=f'边界1 δ={thresholds[0]:.3f}', alpha=0.7, zorder=1)
                ax6.axhline(y=thresholds[1], color='red', linestyle='--', linewidth=2.5, 
                            label=f'边界2 δ={thresholds[1]:.3f}', alpha=0.7, zorder=1)
                
                # 阴影区域
                ax6.axhspan(0, thresholds[0], alpha=0.1, color=colors_phase['partition1'], zorder=0)
                ax6.axhspan(thresholds[0], thresholds[1], alpha=0.1, color=colors_phase['partition2'], zorder=0)
                ax6.axhspan(thresholds[1], df_n3['delta'].max()*1.1, alpha=0.1, color=colors_phase['partition3'], zorder=0)
    
    ax6.set_xlabel('温度 Temperature (K)', fontsize=12, fontweight='bold')
    ax6.set_ylabel('Lindemann 指数 δ', fontsize=12, fontweight='bold')
    ax6.set_title('(f) 自动三分区 (n=3)\nAutomatic 3-Partition (n=3)', 
                  fontsize=13, fontweight='bold', pad=10)
    ax6.legend(fontsize=9, loc='upper left', framealpha=0.95, ncol=2)
    ax6.grid(True, alpha=0.3, linestyle=':', linewidth=0.8)
    
    # ===== 第三行: 温度分布堆叠柱状图（4个策略） =====
    ax8 = fig.add_subplot(gs[2, 0])
    
    if 'fixed2' in comparison:
        df_fixed2 = comparison['fixed2'].get('df_clustered')
        if df_fixed2 is not None:
            temp_phase = df_fixed2.groupby(['temp', 'phase_clustered']).size().unstack(fill_value=0)
            temp_sorted = sorted(df_fixed2['temp'].unique())
            
            # 提取各分区的计数
            partition1_counts = [temp_phase.loc[t, 'partition1'] if t in temp_phase.index and 'partition1' in temp_phase.columns else 0 for t in temp_sorted]
            partition2_counts = [temp_phase.loc[t, 'partition2'] if t in temp_phase.index and 'partition2' in temp_phase.columns else 0 for t in temp_sorted]
            
            # 堆叠柱状图
            x_pos = np.arange(len(temp_sorted))
            ax8.bar(x_pos, partition1_counts, label=phase_labels['partition1'], 
                   color=colors_phase['partition1'], alpha=0.8, edgecolor='black', linewidth=0.8)
            ax8.bar(x_pos, partition2_counts, bottom=partition1_counts, label=phase_labels['partition2'],
                   color=colors_phase['partition2'], alpha=0.8, edgecolor='black', linewidth=0.8)
            
            ax8.set_xlabel('温度 Temperature (K)', fontsize=11, fontweight='bold')
            ax8.set_ylabel('数据点数量 Number of Data Points', fontsize=11, fontweight='bold')
            ax8.set_title('(g) 固定两分区温度分布 (δ=0.08)\nFixed 2-Partition Temperature Distribution', 
                         fontsize=12, fontweight='bold', pad=10)
            ax8.set_xticks(x_pos)
            ax8.set_xticklabels([f'{t:.0f}' for t in temp_sorted], rotation=45, ha='right', fontsize=9)
            ax8.legend(fontsize=9, loc='upper left', framealpha=0.9)
            ax8.grid(True, alpha=0.3, linestyle=':', axis='y', linewidth=0.8)
    
    # ===== 子图9 (第3行第1列): 固定阈值三分区温度分布 =====
    ax9 = fig.add_subplot(gs[2, 1])
    
    if 'fixed' in comparison:
        df_fixed = comparison['fixed'].get('df_clustered')
        if df_fixed is not None:
            temp_phase = df_fixed.groupby(['temp', 'phase_clustered']).size().unstack(fill_value=0)
            temp_sorted = sorted(df_fixed['temp'].unique())
            
            partition1_counts = [temp_phase.loc[t, 'partition1'] if t in temp_phase.index and 'partition1' in temp_phase.columns else 0 for t in temp_sorted]
            partition2_counts = [temp_phase.loc[t, 'partition2'] if t in temp_phase.index and 'partition2' in temp_phase.columns else 0 for t in temp_sorted]
            partition3_counts = [temp_phase.loc[t, 'partition3'] if t in temp_phase.index and 'partition3' in temp_phase.columns else 0 for t in temp_sorted]
            
            x_pos = np.arange(len(temp_sorted))
            
            ax9.bar(x_pos, partition1_counts, label=phase_labels['partition1'], 
                    color=colors_phase['partition1'], alpha=0.8, edgecolor='black', linewidth=0.8)
            ax9.bar(x_pos, partition2_counts, bottom=partition1_counts, label=phase_labels['partition2'], 
                    color=colors_phase['partition2'], alpha=0.8, edgecolor='black', linewidth=0.8)
            ax9.bar(x_pos, partition3_counts, bottom=np.array(partition1_counts)+np.array(partition2_counts), 
                    label=phase_labels['partition3'], color=colors_phase['partition3'], alpha=0.8, edgecolor='black', linewidth=0.8)
            
            ax9.set_xticks(x_pos)
            ax9.set_xticklabels([f'{int(t)}' for t in temp_sorted], rotation=45, ha='right', fontsize=9)
    
    ax9.set_xlabel('温度 Temperature (K)', fontsize=11, fontweight='bold')
    ax9.set_ylabel('数据点数量 Number of Data Points', fontsize=11, fontweight='bold')
    ax9.set_title('(h) 固定三分区温度分布 (δ=0.10, 0.15)\nFixed 3-Partition Temperature Distribution', 
                  fontsize=12, fontweight='bold', pad=10)
    ax9.legend(fontsize=9, loc='upper left', framealpha=0.9)
    ax9.grid(True, alpha=0.3, linestyle=':', axis='y', linewidth=0.8)
    
    # ===== 子图10 (第3行第2列): 自动两分区温度分布 =====
    ax10 = fig.add_subplot(gs[2, 2])
    
    if 'n2' in comparison:
        df_n2 = comparison['n2'].get('df_clustered')
        if df_n2 is not None:
            temp_phase = df_n2.groupby(['temp', 'phase_clustered']).size().unstack(fill_value=0)
            temp_sorted = sorted(df_n2['temp'].unique())
            
            partition1_counts = [temp_phase.loc[t, 'partition1'] if t in temp_phase.index and 'partition1' in temp_phase.columns else 0 for t in temp_sorted]
            partition2_counts = [temp_phase.loc[t, 'partition2'] if t in temp_phase.index and 'partition2' in temp_phase.columns else 0 for t in temp_sorted]
            
            x_pos = np.arange(len(temp_sorted))
            
            ax10.bar(x_pos, partition1_counts, label=phase_labels['partition1'], 
                    color=colors_phase['partition1'], alpha=0.8, edgecolor='black', linewidth=0.8)
            ax10.bar(x_pos, partition2_counts, bottom=partition1_counts, label=phase_labels['partition2'], 
                    color=colors_phase['partition2'], alpha=0.8, edgecolor='black', linewidth=0.8)
            
            ax10.set_xticks(x_pos)
            ax10.set_xticklabels([f'{int(t)}' for t in temp_sorted], rotation=45, ha='right', fontsize=9)
    
    ax10.set_xlabel('温度 Temperature (K)', fontsize=11, fontweight='bold')
    ax10.set_ylabel('数据点数量 Number of Data Points', fontsize=11, fontweight='bold')
    ax10.set_title('(i) 自动两分区温度分布 (n=2)\nAuto 2-Partition Temperature Distribution', 
                  fontsize=12, fontweight='bold', pad=10)
    ax10.legend(fontsize=9, loc='upper left', framealpha=0.9)
    ax10.grid(True, alpha=0.3, linestyle=':', axis='y', linewidth=0.8)
    
    # ===== 子图11 (第3行第3列): 自动三分区温度分布 =====
    ax11 = fig.add_subplot(gs[2, 3])
    
    if 'n3' in comparison:
        df_n3 = comparison['n3'].get('df_clustered')
        if df_n3 is not None:
            temp_phase = df_n3.groupby(['temp', 'phase_clustered']).size().unstack(fill_value=0)
            temp_sorted = sorted(df_n3['temp'].unique())
            
            partition1_counts = [temp_phase.loc[t, 'partition1'] if t in temp_phase.index and 'partition1' in temp_phase.columns else 0 for t in temp_sorted]
            partition2_counts = [temp_phase.loc[t, 'partition2'] if t in temp_phase.index and 'partition2' in temp_phase.columns else 0 for t in temp_sorted]
            partition3_counts = [temp_phase.loc[t, 'partition3'] if t in temp_phase.index and 'partition3' in temp_phase.columns else 0 for t in temp_sorted]
            
            x_pos = np.arange(len(temp_sorted))
            
            ax11.bar(x_pos, partition1_counts, label=phase_labels['partition1'], 
                    color=colors_phase['partition1'], alpha=0.8, edgecolor='black', linewidth=0.8)
            ax11.bar(x_pos, partition2_counts, bottom=partition1_counts, label=phase_labels['partition2'], 
                    color=colors_phase['partition2'], alpha=0.8, edgecolor='black', linewidth=0.8)
            ax11.bar(x_pos, partition3_counts, bottom=np.array(partition1_counts)+np.array(partition2_counts), 
                    label=phase_labels['partition3'], color=colors_phase['partition3'], alpha=0.8, edgecolor='black', linewidth=0.8)
            
            ax11.set_xticks(x_pos)
            ax11.set_xticklabels([f'{int(t)}' for t in temp_sorted], rotation=45, ha='right', fontsize=9)
    
    ax11.set_xlabel('温度 Temperature (K)', fontsize=11, fontweight='bold')
    ax11.set_ylabel('数据点数量 Number of Data Points', fontsize=11, fontweight='bold')
    ax11.set_title('(j) 自动三分区温度分布 (n=3)\nAuto 3-Partition Temperature Distribution', 
                  fontsize=12, fontweight='bold', pad=10)
    ax11.set_xlabel('温度 Temperature (K)', fontsize=11, fontweight='bold')
    ax11.set_ylabel('数据点数量 Number of Data Points', fontsize=11, fontweight='bold')
    ax11.set_title('(j) 自动三分区温度分布 (n=3)\nAuto 3-Partition Temperature Distribution', 
                  fontsize=12, fontweight='bold', pad=10)
    ax11.legend(fontsize=9, loc='upper left', framealpha=0.9)
    ax11.grid(True, alpha=0.3, linestyle=':', axis='y', linewidth=0.8)
    
    plt.tight_layout()
    
    # 保存
    output_file = output_dir / f'{structure_name}_strategy_comparison.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  [OK] 策略对比图已保存: {output_file.name}")
    plt.close()
    
    return output_file

# ============================================================================
# 8. 生成对比报告
# ============================================================================

def generate_strategy_comparison_report(comparison_list, output_dir):
    """生成多策略对比分析报告"""
    print(f"\n{'='*80}")
    print("生成策略对比分析报告")
    print("="*80)
    
    report_file = output_dir / 'strategy_comparison_report.md'
    
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write("# 分区策略对比分析报告\n\n")
        f.write(f"**生成时间**: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write(f"**载体热容**: Cv_support = 38.2151 meV/K (所有Cv_cluster值均已扣除)\n\n")
        f.write("="*80 + "\n\n")
        
        # 添加特征空间对比部分
        f.write("## 特征空间聚类质量对比 (pt6sn8实测)\n\n")
        f.write("不同特征组合对聚类质量的影响（n=3三相分区）：\n\n")
        f.write("| 特征维度 | 特征组合 | Silhouette | 提升幅度 | 推荐度 |\n")
        f.write("|---------|---------|-----------|---------|--------|\n")
        f.write("| 2D | Temperature + Lindemann-δ | 0.4799 | 基线 | ⭐⭐⭐ |\n")
        f.write("| 3D | + Energy | 0.4972 | **+3.6%** | ⭐⭐⭐⭐ |\n")
        f.write("| 3D | + Diffusion-D | 0.4991 | **+4.0%** | ⭐⭐⭐⭐ |\n")
        f.write("| 4D | + Energy + Diffusion-D | 0.5022 | **+4.6%** ✨ | ⭐⭐⭐⭐⭐ |\n\n")
        f.write("**结论**: \n")
        f.write("- 能量特征提升聚类质量 +3.6%\n")
        f.write("- D值特征提升聚类质量 +4.0%\n")
        f.write("- **能量+D值组合效果最佳** (+4.6%)，推荐使用 `--use-energy --use-d-value`\n\n")
        f.write("="*80 + "\n\n")
        
        f.write("## 策略说明\n\n")
        f.write("1. **固定阈值**: 使用传统固定阈值 δ=0.10, 0.15\n")
        f.write("2. **两相自动分区 (n=2)**: K-means自动分为固态和液态两相\n")
        f.write("3. **三相自动分区 (n=3)**: K-means自动分为固态、预熔化、液态三相\n\n")
        
        for comp in comparison_list:
            structure = comp['structure']
            f.write(f"## {structure}\n\n")
            
            # 策略对比表格
            f.write("### 策略质量对比\n\n")
            f.write("| 策略 | 分区数 | Silhouette | Calinski-Harabasz | Davies-Bouldin |\n")
            f.write("|------|--------|------------|-------------------|----------------|\n")
            
            for key in ['fixed', 'n2', 'n3']:
                if key in comp:
                    data = comp[key]
                    method = data['method']
                    n_part = data.get('n_partitions', '-')
                    
                    if 'metrics' in data:
                        metrics = data['metrics']
                        sil = f"{metrics['silhouette']:.4f}"
                        cal = f"{metrics['calinski_harabasz']:.2f}"
                        dav = f"{metrics['davies_bouldin']:.4f}"
                    else:
                        sil = cal = dav = "N/A"
                    
                    f.write(f"| {method} | {n_part} | {sil} | {cal} | {dav} |\n")
            f.write("\n")
            
            # 热容对比
            f.write("### 热容拟合结果对比 (Cv_cluster, 已扣除载体)\n\n")
            f.write("| 策略 | 相态 | Cv_cluster (meV/K) | 误差 (meV/K) | R² | 质量评级 |\n")
            f.write("|------|------|--------------------|--------------|----|-----------|\n")
            
            for key in ['fixed', 'n2', 'n3']:
                if key in comp:
                    data = comp[key]
                    method = data['method']
                    cv_res = data['cv_results']
                    
                    for phase in ['solid', 'premelting', 'liquid']:
                        if phase in cv_res:
                            cv = cv_res[phase]['Cv_cluster']
                            cv_err = cv_res[phase].get('Cv_cluster_err', 0.0)
                            r2 = cv_res[phase]['R2']
                            grade = "★★★" if r2 > 0.95 else "★★" if r2 > 0.90 else "★" if r2 > 0.80 else "⚠"
                            f.write(f"| {method} | {phase} | {cv:.4f} | ±{cv_err:.4f} | {r2:.4f} | {grade} |\n")
            f.write("\n")
            
            # 阈值对比
            f.write("### 相边界阈值对比\n\n")
            f.write("| 策略 | 阈值1 (solid↔premelting) | 阈值2 (premelting↔liquid) |\n")
            f.write("|------|--------------------------|---------------------------|\n")
            
            for key in ['fixed', 'n2', 'n3']:
                if key in comp:
                    data = comp[key]
                    method = data['method']
                    thresholds = data['thresholds']
                    
                    t1 = f"{thresholds[0]:.4f}" if len(thresholds) > 0 else "N/A"
                    t2 = f"{thresholds[1]:.4f}" if len(thresholds) > 1 else "N/A"
                    
                    f.write(f"| {method} | {t1} | {t2} |\n")
            f.write("\n---\n\n")
    
    print(f"  [OK] 策略对比报告已保存: {report_file.name}")
    return report_file

def generate_algorithm_selection_report(algorithm_comparison_results, output_dir):
    """
    生成算法选择报告
    
    Parameters:
    -----------
    algorithm_comparison_results : dict
        格式为 {structure: {method: {'silhouette': ..., 'n_partitions': ...}}}
    output_dir : Path
        输出目录
    """
    print(f"\n{'='*80}")
    print("生成算法选择报告")
    print("="*80)
    
    report_file = output_dir / 'algorithm_selection_report.txt'
    
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write("="*80 + "\n")
        f.write("聚类算法选择报告 (Algorithm Selection Report)\n")
        f.write("="*80 + "\n")
        f.write(f"生成时间: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"分析体系数: {len(algorithm_comparison_results)}\n")
        f.write("="*80 + "\n\n")
        
        # 总结表格
        f.write("=" * 80 + "\n")
        f.write("体系汇总表\n")
        f.write("=" * 80 + "\n")
        f.write(f"{'体系名称':<20} {'最优算法':<15} {'分区数':<10} {'Silhouette':<12} {'质量评级':<10}\n")
        f.write("-" * 80 + "\n")
        
        summary_data = []
        for structure, algo_results in sorted(algorithm_comparison_results.items()):
            # 找出最优算法
            best_method = max(algo_results.keys(), key=lambda m: algo_results[m]['silhouette'])
            best_data = algo_results[best_method]
            silh = best_data['silhouette']
            n_part = best_data['n_partitions']
            
            # 质量评级
            if silh > 0.70:
                quality = "优秀"
            elif silh > 0.50:
                quality = "良好"
            elif silh > 0.25:
                quality = "可接受"
            else:
                quality = "较差"
            
            f.write(f"{structure:<20} {best_method.upper():<15} {n_part:<10} {silh:<12.4f} {quality:<10}\n")
            summary_data.append({
                'structure': structure,
                'method': best_method,
                'n_partitions': n_part,
                'silhouette': silh,
                'quality': quality
            })
        
        f.write("=" * 80 + "\n\n")
        
        # 详细对比结果
        f.write("=" * 80 + "\n")
        f.write("各体系算法详细对比\n")
        f.write("=" * 80 + "\n\n")
        
        for structure, algo_results in sorted(algorithm_comparison_results.items()):
            f.write("-" * 80 + "\n")
            f.write(f"体系: {structure}\n")
            f.write("-" * 80 + "\n")
            
            # 找出最优算法
            best_method = max(algo_results.keys(), key=lambda m: algo_results[m]['silhouette'])
            
            f.write(f"{'算法':<15} {'Silhouette':<12} {'Calinski-H':<12} {'Davies-B':<12} {'分区数':<10} {'推荐':<10}\n")
            f.write("-" * 80 + "\n")
            
            for method in ['kmeans', 'hierarchical', 'dbscan']:
                if method in algo_results:
                    data = algo_results[method]
                    silh = data['silhouette']
                    cal = data.get('calinski', 'N/A')
                    dav = data.get('davies', 'N/A')
                    n_part = data['n_partitions']
                    is_best = '✓ 最优' if method == best_method else ''
                    
                    cal_str = f"{cal:.2f}" if isinstance(cal, (int, float)) else cal
                    dav_str = f"{dav:.4f}" if isinstance(dav, (int, float)) else dav
                    
                    f.write(f"{method.upper():<15} {silh:<12.4f} {cal_str:<12} {dav_str:<12} "
                           f"n={n_part:<8} {is_best:<10}\n")
            
            f.write("\n")
        
        # 算法使用统计
        f.write("=" * 80 + "\n")
        f.write("算法使用统计\n")
        f.write("=" * 80 + "\n\n")
        
        method_counts = {}
        for data in summary_data:
            method = data['method'].upper()
            method_counts[method] = method_counts.get(method, 0) + 1
        
        total = len(summary_data)
        f.write(f"{'算法':<15} {'使用次数':<12} {'占比':<12}\n")
        f.write("-" * 40 + "\n")
        for method, count in sorted(method_counts.items(), key=lambda x: x[1], reverse=True):
            percentage = count / total * 100
            f.write(f"{method:<15} {count:<12} {percentage:>5.1f}%\n")
        
        f.write("\n")
        
        # 分区数统计
        f.write("=" * 80 + "\n")
        f.write("分区数统计\n")
        f.write("=" * 80 + "\n\n")
        
        partition_counts = {}
        for data in summary_data:
            n_part = data['n_partitions']
            partition_counts[n_part] = partition_counts.get(n_part, 0) + 1
        
        f.write(f"{'分区数':<15} {'体系数量':<12} {'占比':<12}\n")
        f.write("-" * 40 + "\n")
        for n_part, count in sorted(partition_counts.items()):
            percentage = count / total * 100
            phase_type = {1: '单相', 2: '两相', 3: '三相'}.get(n_part, f'{n_part}相')
            f.write(f"n={n_part} ({phase_type:<6}){'':<3} {count:<12} {percentage:>5.1f}%\n")
        
        f.write("\n")
        f.write("=" * 80 + "\n")
        f.write("报告结束\n")
        f.write("=" * 80 + "\n")
    
    print(f"  [OK] 算法选择报告已保存: {report_file.name}")
    return report_file

def generate_comparison_report(results_list, output_dir):
    """生成多结构相态分区对比报告"""
    print(f"\n{'='*80}")
    print("生成相态分区对比报告")
    print("="*80)
    
    report_file = output_dir / 'phase_partition_comparison_report.md'
    
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write("# Lindemann指数相态分区分析对比报告\n\n")
        f.write(f"**生成时间**: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write("="*80 + "\n\n")
        
        f.write("## 1. 分析概览\n\n")
        f.write(f"- **分析结构数**: {len(results_list)}\n")
        f.write(f"- **聚类方法**: {', '.join(set([r['method'] for r in results_list]))}\n")
        f.write(f"- **物理约束**: 最大分区数=3 (固态、预熔化、液态)\n\n")
        
        # 结构列表
        f.write("## 2. 结构列表\n\n")
        f.write("| 结构名称 | 方法 | 分区数 | 相态类型 | 数据点 | Silhouette Score |\n")
        f.write("|----------|------|--------|----------|--------|------------------|\n")
        for r in results_list:
            silhouette = r['metrics']['silhouette']
            sil_str = f"{silhouette:.4f}" if silhouette is not None else "N/A"
            n_part = r['n_partitions']
            phase_type = {1: '单相', 2: '两相', 3: '三相'}.get(n_part, f'{n_part}相')
            f.write(f"| {r['structure']} | {r['method']} | {n_part} | {phase_type} | "
                   f"{len(r['df_clustered'])} | {sil_str} |\n")
        f.write("\n")
        
        # 详细结果
        f.write("## 3. 详细结果\n\n")
        for r in results_list:
            n_part = r['n_partitions']
            phase_type = {1: '单相', 2: '两相', 3: '三相'}.get(n_part, f'{n_part}相')
            f.write(f"### {r['structure']} ({r['method'].upper()}, {phase_type})\n\n")
            
            # 分区质量
            f.write("**分区质量指标**:\n")
            if r['metrics']['silhouette'] is not None:
                f.write(f"- Silhouette Score: {r['metrics']['silhouette']:.4f}\n")
                f.write(f"- Calinski-Harabasz Index: {r['metrics']['calinski_harabasz']:.2f}\n")
                f.write(f"- Davies-Bouldin Index: {r['metrics']['davies_bouldin']:.4f}\n")
            f.write("\n")
            
            # 簇统计
            f.write("**簇统计**:\n\n")
            f.write("| 簇ID | 相态 | 数据点 | δ均值 | δ标准差 | δ范围 | 温度范围 |\n")
            f.write("|------|------|--------|--------|---------|-------|----------|\n")
            for stat in r['cluster_stats']:
                phase = r['cluster_to_phase'][stat['cluster_id']]
                f.write(f"| {stat['cluster_id']} | {phase} | {stat['n_points']} | "
                       f"{stat['delta_mean']:.4f} | {stat['delta_std']:.4f} | "
                       f"{stat['delta_min']:.4f}-{stat['delta_max']:.4f} | "
                       f"{stat['temp_min']:.0f}-{stat['temp_max']:.0f}K |\n")
            f.write("\n")
            
            # 相边界阈值
            f.write("**检测到的相边界阈值**:\n")
            for i, thresh in enumerate(r['thresholds']):
                f.write(f"- Threshold {i+1}: δ = {thresh:.4f}\n")
            f.write("\n")
            
            f.write("**与固定阈值对比** (固定: δ=0.10, 0.15):\n")
            if len(r['thresholds']) >= 2:
                diff1 = r['thresholds'][0] - 0.10
                diff2 = r['thresholds'][1] - 0.15
                f.write(f"- Solid/Premelting: {r['thresholds'][0]:.4f} vs 0.10 (差异: {diff1:+.4f})\n")
                f.write(f"- Premelting/Liquid: {r['thresholds'][1]:.4f} vs 0.15 (差异: {diff2:+.4f})\n")
            f.write("\n")
            
            # 热容拟合结果
            if 'cv_results' in r and r['cv_results']:
                f.write("**热容拟合结果** (Cv_cluster, 已扣除载体38.2151 meV/K):\n\n")
                f.write("| 相态 | Cv_cluster (meV/K) | 误差 (meV/K) | R² | 质量评级 |\n")
                f.write("|------|--------------------|--------------|----|----------|\n")
                for phase in ['solid', 'premelting', 'liquid']:
                    if phase in r['cv_results']:
                        cv = r['cv_results'][phase]['Cv_cluster']
                        cv_err = r['cv_results'][phase].get('Cv_cluster_err', 0.0)
                        r2 = r['cv_results'][phase]['R2']
                        grade = "★★★" if r2 > 0.95 else "★★" if r2 > 0.90 else "★" if r2 > 0.80 else "⚠"
                        f.write(f"| {phase} | {cv:.4f} | ±{cv_err:.4f} | {r2:.4f} | {grade} |\n")
                f.write("\n")
            
            f.write("---\n\n")
    
    print(f"  [OK] 对比报告已保存: {report_file.name}")
    return report_file

# ============================================================================
# 6. 主函数
# ============================================================================

def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description='Step 7.4.2: Lindemann指数聚类分析',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例用法:
  # 单个结构分析
  python step7_4_2_clustering_analysis.py --structure pt6sn8
  
  # 多个结构分析
  python step7_4_2_clustering_analysis.py --structure pt6sn8,pt8sn8-1-best
  
  # 分析所有结构
  python step7_4_2_clustering_analysis.py --structure all
  
  # 使用层次聚类
  python step7_4_2_clustering_analysis.py --structure pt6sn8 --method hierarchical
  
  # 自动确定K值
  python step7_4_2_clustering_analysis.py --structure pt6sn8 --auto-k
  
  # 使用DBSCAN
  python step7_4_2_clustering_analysis.py --structure pt6sn8 --method dbscan --eps 0.3
        """
    )
    
    parser.add_argument(
        '--structure', '-s',
        type=str,
        required=True,
        help='要分析的结构名称,多个结构用逗号分隔,或使用"all"分析所有结构'
    )
    
    parser.add_argument(
        '--method', '-m',
        type=str,
        default='kmeans',
        choices=['kmeans', 'hierarchical', 'dbscan'],
        help='聚类方法 (默认: kmeans)'
    )
    
    parser.add_argument(
        '--n-partitions', '-n',
        type=int,
        default=2,
        choices=[1, 2, 3],
        help='相态分区数 (1=单相, 2=两相, 3=三相, 默认: 2)'
    )
    
    parser.add_argument(
        '--auto-partition',
        action='store_true',
        help='自动确定最优分区数 (基于Silhouette score, 范围1-3)'
    )
    
    parser.add_argument(
        '--eps',
        type=float,
        default=0.3,
        help='DBSCAN的epsilon参数 (默认: 0.3)'
    )
    
    parser.add_argument(
        '--min-samples',
        type=int,
        default=5,
        help='DBSCAN的min_samples参数 (默认: 5)'
    )
    
    parser.add_argument(
        '--use-energy',
        action='store_true',
        help='加入能量特征进行聚类'
    )
    
    parser.add_argument(
        '--use-msd',
        action='store_true',
        help='加入MSD特征进行聚类'
    )
    
    parser.add_argument(
        '--use-d-value',
        action='store_true',
        help='加入扩散系数D值特征进行聚类 (来自ensemble_analysis_results.csv)'
    )
    
    parser.add_argument(
        '--compare-algorithms',
        action='store_true',
        help='对比所有聚类算法(KMeans, Hierarchical, DBSCAN)并自动选择最优算法'
    )
    
    args = parser.parse_args()
    
    # 加载数据 (如果使用D值,自动合并)
    df = load_data(load_d_values=args.use_d_value)
    if df is None:
        return
    
    # 确定要分析的结构列表
    if args.structure.lower() == 'all':
        structures = sorted(df['structure'].unique())
        print(f"\n  分析所有 {len(structures)} 个结构")
    else:
        # 大小写不敏感匹配结构名
        input_structures = [s.strip() for s in args.structure.split(',')]
        available_structures = df['structure'].unique()
        structures = []
        for s in input_structures:
            # 尝试精确匹配
            if s in available_structures:
                structures.append(s)
            else:
                # 大小写不敏感匹配
                matched = [a for a in available_structures if a.lower() == s.lower()]
                if matched:
                    structures.append(matched[0])
                    print(f"  [INFO] '{s}' → '{matched[0]}' (大小写自动修正)")
                else:
                    print(f"  [WARNING] 结构 '{s}' 未找到，将跳过")
        print(f"\n  分析 {len(structures)} 个结构: {structures}")
    
    # 对每个结构执行聚类分析
    results_list = []
    comparison_list = []
    algorithm_comparison_results = {}
    
    for structure in structures:
        df_structure = df[df['structure'] == structure].copy()
        
        if len(df_structure) == 0:
            print(f"\n[WARNING] 结构 '{structure}' 未找到,跳过!")
            continue
        
        # 如果启用算法对比模式
        selected_method = args.method
        if args.compare_algorithms:
            print(f"\n{'='*80}")
            print(f"算法对比模式: 测试所有聚类算法")
            print("="*80)
            
            algo_results = {}
            for method in ['kmeans', 'hierarchical', 'dbscan']:
                print(f"\n  [+] 测试算法: {method.upper()}")
                try:
                    result = perform_clustering(
                        df_structure, 
                        structure, 
                        method=method,
                        n_partitions=args.n_partitions,
                        auto_partition=args.auto_partition,
                        eps=args.eps,
                        min_samples=args.min_samples,
                        use_energy=args.use_energy,
                        use_msd=args.use_msd,
                        use_d_value=args.use_d_value
                    )
                    if result and result['metrics']['silhouette'] is not None:
                        algo_results[method] = {
                            'silhouette': result['metrics']['silhouette'],
                            'calinski': result['metrics']['calinski_harabasz'],
                            'davies': result['metrics']['davies_bouldin'],
                            'n_partitions': result['n_partitions'],
                            'result': result
                        }
                        print(f"      Silhouette: {result['metrics']['silhouette']:.4f}")
                except Exception as e:
                    print(f"      [ERROR] {method} 失败: {e}")
            
            # 选择最优算法（基于Silhouette分数）
            if algo_results:
                best_method = max(algo_results.keys(), key=lambda k: algo_results[k]['silhouette'])
                best_score = algo_results[best_method]['silhouette']
                selected_method = best_method
                
                print(f"\n  [✓] 最优算法: {best_method.upper()} (Silhouette={best_score:.4f})")
                
                # 保存对比结果
                algorithm_comparison_results[structure] = algo_results
                
                # 使用最优算法的结果
                results = algo_results[best_method]['result']
            else:
                print(f"\n  [WARNING] 所有算法均失败，使用默认方法")
                results = None
        else:
            # 正常模式：使用指定的单一算法
            results = perform_clustering(
                df_structure, 
                structure, 
                method=selected_method,
                n_partitions=args.n_partitions,
                auto_partition=args.auto_partition,
                eps=args.eps,
                min_samples=args.min_samples,
                use_energy=args.use_energy,
                use_msd=args.use_msd,
                use_d_value=args.use_d_value
            )
        
        if results:
            # 基础可视化
            plot_clustering_results(results, df_structure, OUTPUT_DIR)
            
            # 拟合热容
            cv_results = fit_partition_heat_capacity(results['df_clustered'])
            results['cv_results'] = cv_results
            
            # 获取分区数（在使用前先定义）
            n_partitions = results['n_partitions']
            
            # 获取聚类质量指标
            metrics = results.get('metrics', {})
            silhouette = metrics.get('silhouette', None)
            calinski = metrics.get('calinski_harabasz', None)
            davies = metrics.get('davies_bouldin', None)
            
            # 保存热容拟合质量评分到CSV
            quality_data = []
            for phase, cv_info in cv_results.items():
                row = {
                    'structure': structure,
                    'n_partitions': n_partitions,
                    'phase': phase,
                    'Cv_cluster': cv_info['Cv_cluster'],
                    'Cv_cluster_err': cv_info['Cv_cluster_err'],
                    'R2': cv_info['R2'],
                    'R2_total': cv_info['R2_total'],
                    'grade': cv_info['grade'],
                    'grade_score': cv_info['grade_score']
                }
                # 添加聚类质量指标
                if silhouette is not None:
                    row['silhouette_score'] = silhouette
                if calinski is not None:
                    row['calinski_harabasz'] = calinski
                if davies is not None:
                    row['davies_bouldin'] = davies
                quality_data.append(row)
            quality_df = pd.DataFrame(quality_data)
            quality_file = OUTPUT_DIR / f'{structure}_{args.method}_n{n_partitions}_quality_metrics.csv'
            quality_df.to_csv(quality_file, index=False, encoding='utf-8-sig')
            print(f"  [OK] 质量评分已保存: {quality_file.name}")
            
            results_list.append(results)
            
            # 保存聚类后的数据（文件名包含n值）
            csv_file = OUTPUT_DIR / f'{structure}_{args.method}_n{n_partitions}_clustered_data.csv'
            results['df_clustered'].to_csv(csv_file, index=False, encoding='utf-8-sig')
            print(f"  [OK] 聚类数据已保存: {csv_file.name}")
            
            # 如果指定了n=3，同时运行n=2进行对比
            if args.n_partitions == 3 and not args.auto_partition:
                print(f"\n  [+] 额外运行n=2分析用于对比...")
                results_n2 = perform_clustering(
                    df_structure, 
                    structure, 
                    method=args.method,
                    n_partitions=2,
                    auto_partition=False,
                    eps=args.eps,
                    min_samples=args.min_samples,
                    use_energy=args.use_energy,
                    use_msd=args.use_msd,
                    use_d_value=args.use_d_value
                )
                if results_n2:
                    # 拟合热容
                    cv_results_n2 = fit_partition_heat_capacity(results_n2['df_clustered'])
                    results_n2['cv_results'] = cv_results_n2
                    
                    # 获取聚类质量指标
                    metrics_n2 = results_n2.get('metrics', {})
                    silhouette_n2 = metrics_n2.get('silhouette', None)
                    calinski_n2 = metrics_n2.get('calinski_harabasz', None)
                    davies_n2 = metrics_n2.get('davies_bouldin', None)
                    
                    # 保存热容拟合质量评分到CSV
                    quality_data_n2 = []
                    for phase, cv_info in cv_results_n2.items():
                        row = {
                            'structure': structure,
                            'n_partitions': 2,
                            'phase': phase,
                            'Cv_cluster': cv_info['Cv_cluster'],
                            'Cv_cluster_err': cv_info['Cv_cluster_err'],
                            'R2': cv_info['R2'],
                            'R2_total': cv_info['R2_total'],
                            'grade': cv_info['grade'],
                            'grade_score': cv_info['grade_score']
                        }
                        # 添加聚类质量指标
                        if silhouette_n2 is not None:
                            row['silhouette_score'] = silhouette_n2
                        if calinski_n2 is not None:
                            row['calinski_harabasz'] = calinski_n2
                        if davies_n2 is not None:
                            row['davies_bouldin'] = davies_n2
                        quality_data_n2.append(row)
                    quality_df_n2 = pd.DataFrame(quality_data_n2)
                    quality_file_n2 = OUTPUT_DIR / f'{structure}_{args.method}_n2_quality_metrics.csv'
                    quality_df_n2.to_csv(quality_file_n2, index=False, encoding='utf-8-sig')
                    print(f"  [OK] n=2质量评分已保存: {quality_file_n2.name}")
                    
                    results_list.append(results_n2)
                    # 保存n=2数据
                    csv_file_n2 = OUTPUT_DIR / f'{structure}_{args.method}_n2_clustered_data.csv'
                    results_n2['df_clustered'].to_csv(csv_file_n2, index=False, encoding='utf-8-sig')
                    print(f"  [OK] n=2聚类数据已保存: {csv_file_n2.name}")
            
            # 多策略对比分析
            print(f"\n{'='*80}")
            print(f"执行多策略对比分析: {structure}")
            print("="*80)
            try:
                comparison = compare_partition_strategies(
                    df_structure, structure,
                    use_energy=args.use_energy,
                    use_msd=args.use_msd,
                    use_d_value=args.use_d_value
                )
                comparison['structure'] = structure  # 添加结构名称
                comparison_list.append(comparison)
                
                # 生成独立的详细对比图
                plot_cv_comparison_standalone(comparison, structure, OUTPUT_DIR)
                plot_quality_metrics_standalone(comparison, structure, OUTPUT_DIR)
                
                # 生成主对比图（包含其他子图）
                plot_strategy_comparison(comparison, structure, OUTPUT_DIR)
            except Exception as e:
                print(f"  [ERROR] 策略对比可视化失败: {e}")
                import traceback
                traceback.print_exc()
    
    # 生成算法对比可视化
    if algorithm_comparison_results:
        print(f"\n{'='*80}")
        print("生成算法对比可视化")
        print("="*80)
        for structure, algo_results in algorithm_comparison_results.items():
            plot_algorithm_comparison(algo_results, structure, OUTPUT_DIR)
        
        # 生成算法选择报告
        generate_algorithm_selection_report(algorithm_comparison_results, OUTPUT_DIR)
    
    # 生成对比报告
    if len(results_list) > 0:
        generate_comparison_report(results_list, OUTPUT_DIR)
        
        # 生成策略对比报告
        if len(comparison_list) > 0:
            generate_strategy_comparison_report(comparison_list, OUTPUT_DIR)
        
        print(f"\n{'='*80}")
        print("相态分区分析完成!")
        print("="*80)
        print(f"\n输出目录: {OUTPUT_DIR}")
        print(f"  - {len(results_list)} 个结构的分区结果")
        print(f"  - {len(results_list)} 个可视化图表")
        if algorithm_comparison_results:
            print(f"  - {len(algorithm_comparison_results)} 个算法对比图")
            print(f"  - 1 个算法选择报告 (algorithm_selection_report.txt)")
        print(f"  - 2 个对比报告")
    else:
        print(f"\n[WARNING] 没有成功完成的相态分区分析!")

if __name__ == '__main__':
    main()
