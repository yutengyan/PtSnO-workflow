#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
Step7: 林德曼指数分析 (改进版) + 4级路径签名算法
================================================================================

作者: GitHub Copilot
日期: 2025-10-18
版本: v2.1

功能概述
========
本脚本用于分析LAMMPS模拟的林德曼指数(Lindemann Index)数据,
判断团簇熔化状态,支持多体系分析、MSD异常筛选和熔化温度识别。

林德曼指数判据
==============
- δ < 0.1:  固态 (Solid)
- δ ≥ 0.1:  液态/熔化 (Liquid/Melted)

核心功能
========
1. **智能数据加载与筛选**
   - 合并多个lindemann_master文件 (覆盖r0-r29)
   - 去除完整路径重复
   - 可选MSD异常筛选 (与step1交叉验证)
   - 4级路径签名算法,精确匹配批次信息 (run3, run2等)

2. **林德曼指数分析**
   - 对同一体系同一温度的多个run求平均
   - 熔化温度识别 (δ从<0.1到≥0.1的转变温度)
   - 按体系系列分类 (O1, O2, Pt8SnX, Pt6SnX等)
   - 支持正则表达式体系筛选

3. **熔化温度检测**
   - 基于林德曼指数阈值 (δ=0.1)
   - 线性插值精确定位熔化温度
   - 统计熔化温度范围和趋势

4. **数据可视化**
   - 林德曼指数-温度曲线 (按系列分组)
   - 林德曼指数热力图 (结构 vs 温度)
   - 熔化温度分析图 (Tm vs Sn含量)
   - Sn含量-温度相图

5. **4级路径签名算法** (核心改进)
   - 解决批次重复问题 (run3, run2, run4, run5)
   - MSD路径示例: run3/ptxsn8-x/pt3sn5/t1100.r27.gpu0
   - Lindemann路径示例: run3/ptxsn8-x/pt3sn5/t1100.r27.gpu0
   - 自动检测批次标识符,无批次时使用3级签名
   - 精确1:1匹配,避免误判

命令行参数
==========
--no-filter : 不使用MSD异常筛选,分析所有Lindemann数据
              (适用于某些体系数据全被筛选时,如Pt3Sn5 @ 1100K)

使用示例
========
# 默认模式: 使用MSD异常筛选
python step7_lindemann_analysis.py

# 不筛选模式: 使用所有数据
python step7_lindemann_analysis.py --no-filter

配置方式
========
在脚本中修改 SYSTEM_FILTER 字典:
- include_patterns: 要包含的体系正则表达式列表
- exclude_patterns: 要排除的体系正则表达式列表

示例配置:
SYSTEM_FILTER = {
    'include_patterns': [
        r'[Oo]1',      # O1系列
        r'^pt8sn',     # Pt8SnX体系
    ],
    'exclude_patterns': []
}

输入文件
========
1. lindemann_master_run_*.csv - 林德曼指数数据文件
   列: 结构, 温度(K)/温度, Lindemann指数/林德曼指数, 目录

2. large_D_outliers.csv - MSD异常记录 (可选)
   列: composition, filepath, element, temperature, D_value, ...

3. convergence_master_run_*.csv - 收敛性数据 (可选)
4. lindemann_comparison_run_*.csv - 林德曼对比数据 (可选)

输出文件
========
results/lindemann_analysis/
├── melting_temperatures.csv           # 熔化温度汇总
├── lindemann_analysis_report.txt      # 分析报告
├── Lindemann_vs_T_*.png              # 林德曼指数-温度曲线
├── Lindemann_heatmap_*.png           # 林德曼指数热力图
├── MeltingTemp_vs_SnContent.png      # 熔化温度 vs Sn含量
└── MeltingTemp_PhaseDiagram.png      # 熔化温度相图

关键改进 (v2.1)
===============
1. **4级路径签名算法**
   - 修改前: 890签名 → 900记录 (批次重复,10个run被重复计数)
   - 修改后: 892签名 → 892记录 (完美1:1匹配)
   - 最终结果: 3262 - 892 = 2370 有效记录 (精确过滤)

2. **通用路径提取函数**
   - extract_path_signature(filepath, is_msd_path)
   - 支持MSD路径(有温度目录)和Lindemann路径(无温度目录)
   - 自动检测批次标识符(run3, run2, run4, run5)
   - 向上搜索最多3级目录

3. **命令行参数支持**
   - --no-filter: 跳过MSD异常筛选
   - 解决Pt3Sn5 @ 1100K等体系全被筛选的问题

4. **多次模拟求平均**
   - 自动检测同一体系同一温度的多个run
   - 计算平均值和标准差
   - 统计单次和多次模拟数量

技术细节
========
路径签名算法逻辑:
1. 从文件路径提取run信息 (T1000.r24.gpu0)
2. 定位关键目录索引:
   - MSD路径: 温度目录 (1000K)
   - Lindemann路径: run所在位置
3. 提取composition和parent目录
4. 向上搜索批次标识符 (最多3级)
5. 构建签名:
   - 有批次: batch/parent/composition/run (4级)
   - 无批次: parent/composition/run (3级)

批次标识符: ['run3', 'run2', 'run4', 'run5']

熔化温度检测方法:
1. 对每个体系按温度排序
2. 找到δ从<0.1到≥0.1的转变点
3. 线性插值计算精确的Tm
4. 记录Tm前后的δ值

数据去重策略:
1. 第一步: 去除完整路径重复 (真重复)
2. 第二步: 应用MSD异常筛选 (可选)
3. 第三步: 对同一结构+温度的多个run求平均

注意事项
========
1. 路径签名算法要求MSD和Lindemann数据使用相同的目录结构
2. 支持合并多个lindemann_master文件以覆盖完整的run范围
3. 去重时保留最后一条记录 (keep='last')
4. 熔化温度检测需要δ跨越0.1阈值

依赖库
======
- pandas: 数据处理
- numpy: 数值计算
- matplotlib: 绘图
- scipy: 插值
- argparse: 命令行参数解析

相关脚本
========
- step6_energy_analysis_v2.py: 能量分析 (使用相同的路径签名算法)
- step1_*.py: MSD异常检测 (生成large_D_outliers.csv)

更新日志
========
v2.1 (2025-10-18):
- 新增4级路径签名算法
- 新增--no-filter命令行参数
- 提取通用路径签名函数
- 修复批次重复问题 (892→892精确匹配)

v2.0:
- 支持多体系分析
- 合并多个lindemann文件
- 多次模拟求平均
- 系列分类功能

v1.0:
- 基础林德曼指数分析
- 熔化温度检测
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import re
import argparse
from scipy.interpolate import interp1d


# ===== 通用路径签名算法 =====
def extract_path_signature(filepath, is_msd_path=True):
    """
    从文件路径提取路径签名 (支持批次标识符如run3)
    
    算法策略:
    - 4级签名: batch/parent/composition/run (当检测到run3等批次标识符)
    - 3级签名: parent/composition/run (无批次标识符时)
    
    Args:
        filepath: 完整文件路径
        is_msd_path: True=MSD路径(有温度目录), False=能量路径(无温度目录)
    
    Returns:
        path_signature: 路径签名字符串,如 "run3/o2/o2pt4sn6/t1000.r24.gpu0"
    
    Examples:
        MSD路径:
        >>> extract_path_signature(
        ...     "D:/data/more/run3/o2/O2Pt4Sn6/1000K/T1000.r24.gpu0_msd_Pt.xvg",
        ...     is_msd_path=True
        ... )
        'run3/o2/o2pt4sn6/t1000.r24.gpu0'
        
        Lindemann路径:
        >>> extract_path_signature(
        ...     "/home/data/run3/o2/O2Pt7Sn7/T200.r0.gpu0",
        ...     is_msd_path=False
        ... )
        'run3/o2/o2pt7sn7/t200.r0.gpu0'
    """
    if not filepath:
        return None
    
    # 1. 提取run信息 (T1000.r24.gpu0)
    run_match = re.search(r'(T\d+\.r\d+\.gpu\d+)', filepath, re.IGNORECASE)
    if not run_match:
        return None
    run_info = run_match.group(1).lower()
    
    # 2. 分割路径
    parts = re.split(r'[\\/]', filepath)
    
    # 3. 找到关键目录的索引
    if is_msd_path:
        # MSD路径: 找温度目录 (1000K)
        key_idx = None
        for i, part in enumerate(parts):
            if re.match(r'\d+K$', part, re.IGNORECASE):
                key_idx = i
                break
    else:
        # Lindemann路径: 找run所在位置
        key_idx = None
        for i, part in enumerate(parts):
            if re.search(r'T\d+\.r\d+\.gpu\d+', part, re.IGNORECASE):
                key_idx = i
                break
    
    if key_idx is None or key_idx < 2:
        # 无法提取足够的层级,返回简化签名
        return run_info
    
    # 4. 提取目录层级
    composition_dir = parts[key_idx - 1].lower()  # O2Pt4Sn6 或 pt8sn5-1-best
    parent_dir = parts[key_idx - 2].lower()       # o2 或 Pt8
    
    # 5. 检查批次标识符 (run3, run2, run4, run5)
    batch_keywords = ['run3', 'run2', 'run4', 'run5']
    path_signature = f"{parent_dir}/{composition_dir}/{run_info}"
    
    # 向上搜索批次标识符 (最多向上3级)
    if key_idx >= 3:
        for check_idx in range(key_idx - 3, max(-1, key_idx - 6), -1):
            if check_idx < 0 or check_idx >= len(parts):
                break
            check_dir = parts[check_idx].lower()
            if check_dir in batch_keywords:
                # 找到批次标识,构建4级签名
                path_signature = f"{check_dir}/{parent_dir}/{composition_dir}/{run_info}"
                break
    
    return path_signature


# ============================================================================
# 配置部分
# ============================================================================

# 中文显示
plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# 基础路径
BASE_DIR = Path(__file__).parent

# 输入文件
DATA_DIR = BASE_DIR / 'data' / 'lindemann'
# ⚠️ 修改: 加载所有lindemann_master文件以获得完整覆盖(r0-r29)
LINDEMANN_FILES = sorted(DATA_DIR.glob('lindemann_master_run_*.csv'))
CONVERGENCE_FILE = DATA_DIR / 'convergence_master_run_20251015_170013.csv'
COMPARISON_FILE = DATA_DIR / 'lindemann_comparison_run_20251015_170013.csv'

# 异常值文件 (与step6保持一致)
OUTLIERS_FILE = BASE_DIR / 'results' / 'large_D_outliers.csv'

# 输出目录
OUTPUT_DIR = BASE_DIR / 'results' / 'lindemann_analysis'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# 系统筛选配置 (与step6保持一致)
SYSTEM_FILTER = {
    'include_patterns': [
        # ===== 含氧体系 (Pt-Sn-O三元系统) =====
        r'[Oo]1',          # O1系列
        r'[Oo]2',          # O2系列
        r'[Oo]3',          # O3系列
        r'[Oo]4',          # O4系列
        
        # ===== Pt(8-x)SnX体系 (总原子数=8,变Pt/Sn比) =====
        r'^pt7sn1',        # Pt7Sn1: 7+1=8
        r'^pt6sn2',        # Pt6Sn2: 6+2=8
        r'^pt5sn3',        # Pt5Sn3: 5+3=8
        r'^pt4sn4',        # Pt4Sn4: 4+4=8
        r'^pt3sn5',        # Pt3Sn5: 3+5=8
        r'^pt2sn6',        # Pt2Sn6: 2+6=8
        r'^pt8$',          # Pt8: 纯Pt (用于对比)
        
        # ===== Pt6SnX体系 (固定6个Pt原子,变Sn含量) =====
        r'^pt6',           # Pt6SnX体系: pt6, pt6sn1, ..., pt6sn9 (用于与Pt8对比)
        
        # ===== Pt8SnX体系 (固定8个Pt原子,变Sn含量) =====
        r'^pt8sn',         # Pt8SnX体系: pt8sn0, pt8sn1, ..., pt8sn10
        
        # ===== Cv系列 (特殊命名规则) =====
        r'^Cv-',           # Cv系列: Cv-1, Cv-2, Cv-3, Cv-4, Cv-5
    ],
    'exclude_patterns': []
}

# 林德曼指数阈值
LINDEMANN_THRESHOLDS = {
    'solid': 0.1,      # δ < 0.1: 固态
    'melting': 0.1     # δ ≥ 0.1: 液态/融化
}


# ============================================================================
# 辅助函数
# ============================================================================

def build_path_filter_set(outliers_df):
    """
    从outliers构建路径过滤集合 (使用通用路径签名算法)
    
    MSD路径示例:
    1. Pt8SnX系列: d:\...\GPU-Pt8\Pt8\pt8sn0-2-best\700K\T700.r4.gpu0_msd_Pt.xvg
    2. 含氧系列: d:\...\o68\o1-2\g-1-O1Sn4Pt3\1000K\T1000.r24.gpu0_msd_Pt.xvg
    3. run3批次: d:\...\more\run3\o2\O2Pt4Sn6\1000K\T1000.r24.gpu0_msd_Pt.xvg
    
    ✅ 使用通用extract_path_signature函数
       - 4级签名: run3/o2/o2pt4sn6/t1000.r24.gpu0
       - 3级签名: o1-2/g-1-o1sn4pt3/t1100.r9.gpu0
    """
    filter_path_signature_set = set()
    
    for _, row in outliers_df.iterrows():
        filepath = row.get('filepath', '')
        if not filepath:
            continue
        
        # 使用通用路径签名函数 (MSD路径有温度目录)
        path_signature = extract_path_signature(filepath, is_msd_path=True)
        if path_signature:
            filter_path_signature_set.add(path_signature)
    
    print(f"  [INFO] Built path signature filter set (Complete Algorithm with Batch Info):")
    print(f"    - Unique path signatures: {len(filter_path_signature_set)}")
    print(f"    - Expected: ~890 unique paths (2227 MSD anomalies)")
    if len(filter_path_signature_set) > 0:
        print(f"    - Sample signatures:")
        for idx, sig in enumerate(sorted(list(filter_path_signature_set))[:5]):
            print(f"      {idx+1}. {sig}")
    
    return filter_path_signature_set


def should_filter_lindemann_record(row, filter_path_signature_set):
    """
    判断林德曼记录是否应该被筛选掉 (使用通用路径签名算法)
    
    ✅ 使用通用extract_path_signature函数
    
    林德曼记录路径示例 (Linux格式,无温度目录):
    1. Pt8SnX: /home/.../GPU-Pt8/Pt8/pt8sn0-2-best/T400.r1.gpu0
    2. 含氧系列: /home/.../o68/o1-2/g-1-O1Sn4Pt3/T1000.r24.gpu0
    3. run3批次: /home/.../more/run3/o2/O2Pt7Sn7/T200.r0.gpu0
    """
    full_path = row.get('目录', '')  # 林德曼数据的路径列名是'目录'
    if not full_path:
        return False
    
    # 使用通用路径签名函数 (Lindemann路径无温度目录)
    path_signature = extract_path_signature(full_path, is_msd_path=False)
    if not path_signature:
        return False
    
    return path_signature in filter_path_signature_set


def extract_chemical_formula(structure_name):
    """
    从结构名提取化学式
    
    Examples:
    - pt8sn5-2-best -> pt8sn5
    - g-1535-Sn8Pt6O4 -> Sn8Pt6O4
    - O2Pt4Sn6 -> O2Pt4Sn6
    """
    if not structure_name:
        return None
    
    # 如果是g-开头的格式
    if structure_name.startswith('g-'):
        match = re.search(r'g-\d+-(.+)', structure_name, re.IGNORECASE)
        if match:
            return match.group(1)
    
    # 移除后缀 (-1-best, -2-best等)
    clean_name = re.sub(r'-\d+-best$', '', structure_name, flags=re.IGNORECASE)
    
    return clean_name


def classify_system_series(structure_name):
    """
    分类体系所属系列
    
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


def extract_pt_sn_o_atoms(structure_name, filepath=None):
    """
    从结构名提取Pt、Sn、O原子数
    对于Cv系列,从路径中的父目录提取化学式
    
    Parameters:
    -----------
    structure_name : str
        结构名称 (如 "Cv-1", "pt8sn5", "O2Pt4Sn6")
    filepath : str, optional
        完整路径 (用于Cv系列从父目录提取化学式)
    
    Returns:
    --------
    (n_pt, n_sn, n_o) : tuple
    """
    # 对于Cv系列,尝试从路径中提取化学式
    if structure_name and structure_name.startswith('Cv-') and filepath:
        parts = re.split(r'[\\/]', filepath)
        # Cv系列路径示例: /home/.../g-1535-Sn8Pt6O4/Cv-1/T1000.r16.gpu0
        # 找到Cv-X所在位置,向前一级提取g-XXXX-ChemicalFormula
        for i, part in enumerate(parts):
            if part.startswith('Cv-') and i > 0:
                parent_dir = parts[i-1]
                if parent_dir.startswith('g-'):
                    # 从 g-1535-Sn8Pt6O4 中提取 Sn8Pt6O4
                    match = re.search(r'g-\d+-(.+)', parent_dir, re.IGNORECASE)
                    if match:
                        chemical_formula = match.group(1).lower()
                        # 从化学式中提取原子数
                        structure_lower = chemical_formula
                        break
        else:
            structure_lower = structure_name.lower()
    else:
        structure_lower = structure_name.lower()
    
    n_o = 0
    n_pt = 0
    n_sn = 0
    
    # 提取O原子数
    o_matches = list(re.finditer(r'o(\d+)', structure_lower))
    if o_matches:
        n_o = sum(int(m.group(1)) for m in o_matches)
    
    # 提取Pt原子数
    pt_match = re.search(r'pt(\d+)', structure_lower)
    if pt_match:
        n_pt = int(pt_match.group(1))
    
    # 提取Sn原子数
    sn_match = re.search(r'sn(\d+)', structure_lower)
    if sn_match:
        n_sn = int(sn_match.group(1))
    
    if n_pt == 0 and n_sn == 0 and n_o == 0:
        return (None, None, None)
    
    return (n_pt, n_sn, n_o)


def filter_structures(structures, include_patterns=None, exclude_patterns=None):
    """
    根据正则表达式筛选结构
    """
    if not include_patterns and not exclude_patterns:
        return structures
    
    filtered = []
    
    for structure in structures:
        # 检查include模式
        if include_patterns:
            included = any(re.search(pattern, structure, re.IGNORECASE) 
                          for pattern in include_patterns)
            if not included:
                continue
        
        # 检查exclude模式
        if exclude_patterns:
            excluded = any(re.search(pattern, structure, re.IGNORECASE) 
                          for pattern in exclude_patterns)
            if excluded:
                continue
        
        filtered.append(structure)
    
    return filtered


def estimate_melting_temperature(df_struct, method='inflection'):
    """
    估算熔化温度
    
    Parameters:
    -----------
    df_struct : DataFrame
        单个结构的林德曼指数数据
    method : str
        'threshold': δ=0.1时的温度
        'inflection': 曲线拐点
        'derivative': 导数最大值
    
    Returns:
    --------
    Tm : float or None
        熔化温度(K)
    """
    if len(df_struct) < 3:
        return None
    
    df_sorted = df_struct.sort_values('温度')
    temps = df_sorted['温度'].values
    lindemann = df_sorted['林德曼指数'].values
    
    if method == 'threshold':
        # 找到δ=0.1时的温度
        threshold = LINDEMANN_THRESHOLDS['melting']
        
        # 如果最大值都小于阈值,说明未熔化
        if lindemann.max() < threshold:
            return None
        
        # 如果最小值都大于阈值,说明全部熔化
        if lindemann.min() > threshold:
            return temps[0]
        
        # 线性插值
        try:
            f = interp1d(lindemann, temps, kind='linear', fill_value='extrapolate')
            Tm = float(f(threshold))
            return Tm
        except:
            return None
    
    elif method == 'derivative':
        # 计算导数最大值对应的温度
        if len(temps) < 4:
            return None
        
        # 数值导数
        dL_dT = np.gradient(lindemann, temps)
        max_idx = np.argmax(dL_dT)
        
        return temps[max_idx]
    
    elif method == 'inflection':
        # 拐点法: 二阶导数为0
        if len(temps) < 5:
            return None
        
        # 二阶导数
        d2L_dT2 = np.gradient(np.gradient(lindemann, temps), temps)
        
        # 找到二阶导数最接近0的点
        zero_cross_idx = np.argmin(np.abs(d2L_dT2))
        
        return temps[zero_cross_idx]
    
    return None


# ============================================================================
# 主分析函数
# ============================================================================

def load_and_filter_data(no_filter=False):
    """
    加载并筛选数据
    
    Args:
        no_filter: bool, 如果为True,跳过MSD异常筛选,使用所有数据
    
    改进:
    1. 筛选掉step1标记的异常run (可选)
    2. 对同一体系同一温度的多个run进行平均
    """
    print("=" * 80)
    print("Step7: 林德曼指数分析 (改进版)")
    if no_filter:
        print("[!] 模式: 不使用MSD异常筛选 (--no-filter)")
    else:
        print("[+] 模式: 使用MSD异常筛选 (默认)")
    print("=" * 80)
    
    # 1. 加载林德曼数据 (合并所有文件以获得完整覆盖)
    print("\n[*] Loading lindemann data...")
    print(f"  [INFO] Found {len(LINDEMANN_FILES)} lindemann files")
    
    lindemann_dfs = []
    for file in LINDEMANN_FILES:
        df = pd.read_csv(file)
        df.rename(columns={'温度(K)': '温度', 'Lindemann指数': '林德曼指数'}, inplace=True)
        lindemann_dfs.append(df)
        print(f"    - {file.name}: {len(df)} records")
    
    # 合并所有数据
    df_lindemann = pd.concat(lindemann_dfs, ignore_index=True)
    print(f"  [OK] Merged: {len(df_lindemann)} records")
    
    # Step 1: 去除完整路径重复(真重复)
    before_path_dedup = len(df_lindemann)
    df_lindemann = df_lindemann.drop_duplicates(subset=['目录'], keep='last')
    after_path_dedup = len(df_lindemann)
    
    if before_path_dedup > after_path_dedup:
        print(f"  [DEDUP] 去除完整路径重复: -{before_path_dedup - after_path_dedup} 条 → {after_path_dedup} 条")
    
    # Step 2: 统计多次模拟(同一结构+温度,不同路径)
    multi_sim_count = df_lindemann.groupby(['结构', '温度']).size()
    multi_sim_groups = (multi_sim_count > 1).sum()
    if multi_sim_groups > 0:
        print(f"  [INFO] 检测到 {multi_sim_groups} 组体系有多次模拟(将在后续求平均)")
    
    print(f"  [OK] 去重后数据: {after_path_dedup} 条记录")
    
    # 统计run范围
    import re
    runs = []
    for path in df_lindemann['目录']:
        run_match = re.search(r'\.r(\d+)\.gpu', path, re.IGNORECASE)
        if run_match:
            runs.append(int(run_match.group(1)))
    if runs:
        print(f"  [INFO] Run range: r{min(runs)} - r{max(runs)} ({len(set(runs))} unique runs)")
    
    # 2. 加载并应用outliers筛选 (可选)
    if no_filter:
        print("\n[*] [!] Skipping MSD outlier filtering (--no-filter mode)")
        print(f"  [INFO] Using all {len(df_lindemann)} records without filtering")
    elif OUTLIERS_FILE.exists():
        print("\n[*] Loading outliers for filtering...")
        df_outliers = pd.read_csv(OUTLIERS_FILE)
        print(f"  [OK] Loaded {len(df_outliers)} outlier runs")
        
        filter_path_signature_set = build_path_filter_set(df_outliers)
        
        df_lindemann['is_outlier'] = df_lindemann.apply(
            lambda row: should_filter_lindemann_record(row, filter_path_signature_set),
            axis=1
        )
        
        n_outliers = df_lindemann['is_outlier'].sum()
        print(f"  [FILTER] Identified {n_outliers} outlier records in lindemann data")
        print(f"    - Using precise path signature matching (避免误判)")
        
        # 筛选掉outliers
        df_lindemann = df_lindemann[~df_lindemann['is_outlier']].copy()
        print(f"  [OK] After filtering: {len(df_lindemann)} valid records")
    else:
        print(f"\n[WARNING] Outliers file not found: {OUTLIERS_FILE}")
        print(f"  [INFO] Proceeding without outlier filtering")
    
    # 3. 应用系统筛选
    include = SYSTEM_FILTER.get('include_patterns', [])
    exclude = SYSTEM_FILTER.get('exclude_patterns', [])
    
    if include or exclude:
        print("\n[*] Applying system filters...")
        print(f"    Include patterns: {include}")
        print(f"    Exclude patterns: {exclude}")
        
        all_structures = df_lindemann['结构'].unique()
        filtered_structures = filter_structures(all_structures, include, exclude)
        
        df_lindemann = df_lindemann[df_lindemann['结构'].isin(filtered_structures)].copy()
        
        print(f"    Filtered: {len(filtered_structures)}/{len(all_structures)} structures selected")
        print(f"    Records: {len(df_lindemann)}")
    
    # 4. 对同一体系同一温度的多个run进行平均(包含多次模拟)
    print("\n[*] Averaging multiple runs per structure-temperature...")
    print(f"  [INFO] Before averaging: {len(df_lindemann)} records")
    
    # 按结构和温度分组,计算平均值和标准差
    df_averaged = df_lindemann.groupby(['结构', '温度']).agg({
        '林德曼指数': ['mean', 'std', 'count'],  # 平均值、标准差、样本数
        '目录': 'first',      # 保留第一个路径(用于记录)
        '方法': 'first',      # 保留第一个方法
    }).reset_index()
    
    # 展平列名
    df_averaged.columns = ['结构', '温度', '林德曼指数', 'δ标准差', 'run_count', '目录', '方法']
    
    # 统计信息
    multi_runs = df_averaged[df_averaged['run_count'] > 1]
    print(f"  [OK] After averaging: {len(df_averaged)} unique structure-temperature points")
    print(f"  [INFO] 单次模拟: {len(df_averaged[df_averaged['run_count'] == 1])} points")
    print(f"  [INFO] 多次模拟: {len(multi_runs)} points (包含多次模拟的体系)")
    if len(multi_runs) > 0:
        print(f"    - 平均模拟次数: {multi_runs['run_count'].mean():.1f}")
        print(f"    - 最多模拟次数: {multi_runs['run_count'].max()}")
    
    # 5. 添加系列分类和原子数信息
    print("\n[*] Classifying systems...")
    df_averaged['系列'] = df_averaged['结构'].apply(classify_system_series)
    
    # 提取原子数 (传入路径信息以支持Cv系列)
    atom_info = df_averaged.apply(
        lambda row: extract_pt_sn_o_atoms(row['结构'], row['目录']),
        axis=1
    )
    df_averaged['Pt原子数'] = atom_info.apply(lambda x: x[0] if x[0] is not None else 0)
    df_averaged['Sn原子数'] = atom_info.apply(lambda x: x[1] if x[1] is not None else 0)
    df_averaged['O原子数'] = atom_info.apply(lambda x: x[2] if x[2] is not None else 0)
    df_averaged['总原子数'] = df_averaged['Pt原子数'] + df_averaged['Sn原子数'] + df_averaged['O原子数']
    
    # 系列统计
    series_counts = df_averaged.groupby('系列').size()
    for series, count in series_counts.items():
        print(f"    - {series}: {count} records")
    
    return df_averaged


def analyze_melting_temperatures(df_lindemann):
    """
    分析熔化温度
    """
    print("\n[*] Analyzing melting temperatures...")
    
    melting_data = []
    
    for structure in df_lindemann['结构'].unique():
        df_struct = df_lindemann[df_lindemann['结构'] == structure].copy()
        
        # 使用三种方法估算Tm
        Tm_threshold = estimate_melting_temperature(df_struct, method='threshold')
        Tm_derivative = estimate_melting_temperature(df_struct, method='derivative')
        Tm_inflection = estimate_melting_temperature(df_struct, method='inflection')
        
        # 获取其他信息
        series = df_struct['系列'].iloc[0]
        pt = df_struct['Pt原子数'].iloc[0]
        sn = df_struct['Sn原子数'].iloc[0]
        o = df_struct['O原子数'].iloc[0]
        total = df_struct['总原子数'].iloc[0]
        
        melting_data.append({
            '结构': structure,
            '系列': series,
            'Pt原子数': int(pt),
            'Sn原子数': int(sn),
            'O原子数': int(o),
            '总原子数': int(total),
            'Tm_threshold': Tm_threshold,
            'Tm_derivative': Tm_derivative,
            'Tm_inflection': Tm_inflection
        })
    
    df_melting = pd.DataFrame(melting_data)
    
    # 保存
    output_file = OUTPUT_DIR / 'melting_temperatures.csv'
    df_melting.to_csv(output_file, index=False, encoding='utf-8-sig')
    print(f"  [OK] Saved: {output_file.name}")
    
    return df_melting


def analyze_cv_series_lindemann(df_lindemann):
    """
    专门分析Cv系列(Cv-1到Cv-5)的林德曼指数
    该系列为Sn8Pt6O4组分,跑了5遍,温度间隔50K
    
    Parameters:
    -----------
    df_lindemann : DataFrame
        完整的林德曼指数数据
        
    Returns:
    --------
    df_cv_stats : DataFrame
        Cv系列的统计分析结果
    """
    print("\n" + "="*80)
    print("[*] Step 7.2: 分析Cv系列(Sn8Pt6O4)的林德曼指数")
    print("="*80)
    
    # 筛选Cv系列数据
    df_cv = df_lindemann[df_lindemann['结构'].str.contains('Cv-', na=False)].copy()
    
    if len(df_cv) == 0:
        print("  [WARNING] No Cv series data found!")
        return None
    
    print(f"\n[1] 数据概况")
    print("-"*80)
    print(f"  记录数: {len(df_cv)}")
    print(f"  结构: {sorted(df_cv['结构'].unique())}")
    print(f"  温度范围: {df_cv['温度'].min():.0f} - {df_cv['温度'].max():.0f} K")
    print(f"  温度点数: {df_cv['温度'].nunique()}")
    
    # 统计分析
    print(f"\n[2] 林德曼指数统计(按温度)")
    print("-"*80)
    
    stats_by_temp = df_cv.groupby('温度').agg({
        '林德曼指数': ['mean', 'std', 'min', 'max', 'count']
    }).round(4)
    
    print("\n林德曼指数统计:")
    print(stats_by_temp)
    
    # 判断熔化状态
    print(f"\n[3] 熔化状态分析")
    print("-"*80)
    
    # 按温度计算平均林德曼指数并判断状态
    temp_analysis = []
    for temp in sorted(df_cv['温度'].unique()):
        df_temp = df_cv[df_cv['温度'] == temp]
        mean_lind = df_temp['林德曼指数'].mean()
        std_lind = df_temp['林德曼指数'].std()
        count = len(df_temp)
        
        # 判断状态
        if mean_lind < 0.1:
            state = '固态 (Solid)'
        else:
            state = '液态 (Liquid)'
        
        temp_analysis.append({
            '温度': temp,
            '平均林德曼指数': mean_lind,
            '标准差': std_lind,
            '测量次数': count,
            '状态': state
        })
        
        print(f"  {temp:4.0f}K: δ = {mean_lind:.4f} ± {std_lind:.4f} ({count}次) -> {state}")
    
    df_cv_stats = pd.DataFrame(temp_analysis)
    
    # 估算熔化温度
    print(f"\n[4] 熔化温度估算")
    print("-"*80)
    
    # 找到δ从<0.1到≥0.1的转变温度
    solid_temps = df_cv_stats[df_cv_stats['平均林德曼指数'] < 0.1]['温度'].values
    liquid_temps = df_cv_stats[df_cv_stats['平均林德曼指数'] >= 0.1]['温度'].values
    
    if len(solid_temps) > 0 and len(liquid_temps) > 0:
        T_solid_max = solid_temps.max()
        T_liquid_min = liquid_temps.min()
        
        # 线性插值估算Tm
        if T_liquid_min - T_solid_max <= 100:  # 温度间隔合理
            df_below = df_cv_stats[df_cv_stats['温度'] == T_solid_max].iloc[0]
            df_above = df_cv_stats[df_cv_stats['温度'] == T_liquid_min].iloc[0]
            
            # 线性插值: Tm = T1 + (0.1 - δ1)/(δ2 - δ1) * (T2 - T1)
            delta1 = df_below['平均林德曼指数']
            delta2 = df_above['平均林德曼指数']
            T1 = T_solid_max
            T2 = T_liquid_min
            
            if delta2 > delta1:
                Tm = T1 + (0.1 - delta1) / (delta2 - delta1) * (T2 - T1)
                print(f"  熔化温度 Tm ≈ {Tm:.1f} K")
                print(f"  (基于 {T_solid_max:.0f}K: δ={delta1:.4f} 和 {T_liquid_min:.0f}K: δ={delta2:.4f})")
            else:
                print(f"  无法估算 (林德曼指数未单调增加)")
        else:
            print(f"  熔化温度范围: {T_solid_max:.0f} - {T_liquid_min:.0f} K")
    elif len(solid_temps) == 0:
        print(f"  所有温度均为液态 (Tm < {df_cv_stats['温度'].min():.0f}K)")
    else:
        print(f"  所有温度均为固态 (Tm > {df_cv_stats['温度'].max():.0f}K)")
    
    # 绘图
    print(f"\n[5] 生成可视化图表")
    print("-"*80)
    
    import matplotlib.pyplot as plt
    import matplotlib
    matplotlib.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei']
    matplotlib.rcParams['axes.unicode_minus'] = False
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Cv系列(Sn8Pt6O4) 林德曼指数分析 - 5次独立模拟', 
                fontsize=16, fontweight='bold')
    
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    markers = ['o', 's', '^', 'D', 'v']
    
    # 子图1: 各次模拟的林德曼指数
    ax1 = axes[0, 0]
    for i, struct in enumerate(sorted(df_cv['结构'].unique())):
        df_struct = df_cv[df_cv['结构'] == struct].sort_values('温度')
        ax1.plot(df_struct['温度'], df_struct['林德曼指数'], 
                marker=markers[i], color=colors[i], label=struct,
                linewidth=2, markersize=8, alpha=0.7)
    
    # 添加平均值曲线
    df_mean = df_cv.groupby('温度')['林德曼指数'].mean().reset_index()
    ax1.plot(df_mean['温度'], df_mean['林德曼指数'],
            'k--', linewidth=3, label='平均值', zorder=10)
    
    # 添加阈值线
    ax1.axhline(0.10, color='green', linestyle=':', linewidth=2, 
               label='熔化阈值 (δ=0.10)', alpha=0.7)
    
    ax1.set_xlabel('温度 [K]', fontsize=12, fontweight='bold')
    ax1.set_ylabel('林德曼指数 δ', fontsize=12, fontweight='bold')
    ax1.set_title('林德曼指数 vs 温度', fontsize=13, fontweight='bold')
    ax1.legend(loc='best', fontsize=9, ncol=2)
    ax1.grid(True, alpha=0.3)
    
    # 子图2: 平均值±标准差
    ax2 = axes[0, 1]
    df_stats_plot = df_cv.groupby('温度').agg({
        '林德曼指数': ['mean', 'std']
    }).reset_index()
    df_stats_plot.columns = ['温度', 'mean', 'std']
    
    ax2.errorbar(df_stats_plot['温度'], df_stats_plot['mean'], 
                yerr=df_stats_plot['std'],
                fmt='o-', color='darkblue', linewidth=2, markersize=8,
                capsize=5, capthick=2, label='Mean ± Std')
    ax2.fill_between(df_stats_plot['温度'], 
                    df_stats_plot['mean']-df_stats_plot['std'],
                    df_stats_plot['mean']+df_stats_plot['std'],
                    alpha=0.2, color='blue')
    
    ax2.axhline(0.10, color='green', linestyle=':', linewidth=2, alpha=0.7, label='熔化阈值 (δ=0.10)')
    
    ax2.set_xlabel('温度 [K]', fontsize=12, fontweight='bold')
    ax2.set_ylabel('林德曼指数 δ', fontsize=12, fontweight='bold')
    ax2.set_title('平均林德曼指数(含误差棒)', fontsize=13, fontweight='bold')
    ax2.legend(loc='best', fontsize=10)
    ax2.grid(True, alpha=0.3)
    
    # 子图3: 相对标准差(RSD%)
    ax3 = axes[1, 0]
    df_rsd = df_cv.groupby('温度').agg({
        '林德曼指数': lambda x: (x.std() / x.mean() * 100) if x.mean() != 0 else 0
    }).reset_index()
    df_rsd.columns = ['温度', 'RSD']
    
    ax3.bar(df_rsd['温度'], df_rsd['RSD'], width=40, 
           color='coral', edgecolor='darkred', linewidth=2, alpha=0.7)
    ax3.axhline(10, color='orange', linestyle='--', linewidth=2, 
               label='10% 参考线')
    ax3.axhline(20, color='red', linestyle='--', linewidth=2, 
               label='20% 参考线')
    
    ax3.set_xlabel('温度 [K]', fontsize=12, fontweight='bold')
    ax3.set_ylabel('相对标准差 RSD [%]', fontsize=12, fontweight='bold')
    ax3.set_title('测量重复性', fontsize=13, fontweight='bold')
    ax3.legend(loc='best', fontsize=10)
    ax3.grid(True, alpha=0.3, axis='y')
    
    # 子图4: 状态分布
    ax4 = axes[1, 1]
    
    # 统计各状态的温度点数
    state_counts = df_cv_stats['状态'].value_counts()
    colors_state = {'固态 (Solid)': 'blue', 
                   '液态 (Liquid)': 'red'}
    
    bars = ax4.bar(range(len(state_counts)), state_counts.values, 
                  color=[colors_state.get(s, 'gray') for s in state_counts.index],
                  edgecolor='black', linewidth=2, alpha=0.7)
    ax4.set_xticks(range(len(state_counts)))
    ax4.set_xticklabels(state_counts.index, rotation=15, ha='right')
    ax4.set_ylabel('温度点数', fontsize=12, fontweight='bold')
    ax4.set_title('熔化状态分布', fontsize=13, fontweight='bold')
    ax4.grid(True, alpha=0.3, axis='y')
    
    # 添加数值标签
    for i, bar in enumerate(bars):
        height = bar.get_height()
        ax4.text(bar.get_x() + bar.get_width()/2., height,
                f'{int(height)}',
                ha='center', va='bottom', fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    
    # 保存图片
    fig_path = OUTPUT_DIR / 'Cv_series_lindemann_analysis.png'
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    print(f"  [OK] 已保存图片: {fig_path.name}")
    plt.close()
    
    # 保存统计数据
    csv_path = OUTPUT_DIR / 'Cv_series_lindemann_statistics.csv'
    df_cv_stats.to_csv(csv_path, index=False, encoding='utf-8-sig')
    print(f"  [OK] 已保存统计数据: {csv_path.name}")
    
    print("\n" + "="*80)
    print("[DONE] Cv系列林德曼指数分析完成!")
    print("="*80)
    
    return df_cv_stats


def plot_lindemann_vs_temperature(df_lindemann, df_melting):
    """
    绘制林德曼指数 vs 温度曲线
    """
    print("\n[*] Plotting Lindemann index vs temperature...")
    
    # 按系列分组绘制
    series_list = sorted(df_lindemann['系列'].unique())
    
    for series in series_list:
        df_series = df_lindemann[df_lindemann['系列'] == series].copy()
        
        if len(df_series) == 0:
            continue
        
        # 获取该系列的所有结构
        structures = sorted(df_series['结构'].unique())
        
        if len(structures) == 0:
            continue
        
        # 创建图形
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # 颜色映射
        colors = plt.cm.viridis(np.linspace(0, 1, len(structures)))
        
        for idx, structure in enumerate(structures):
            df_struct = df_series[df_series['结构'] == structure].sort_values('温度')
            
            # 绘制曲线
            ax.plot(df_struct['温度'], df_struct['林德曼指数'],
                   marker='o', color=colors[idx], linewidth=2,
                   markersize=6, alpha=0.7, label=structure)
            
            # 标记熔化温度
            if df_melting is not None:
                tm_row = df_melting[df_melting['结构'] == structure]
                if len(tm_row) > 0 and pd.notna(tm_row['Tm_threshold'].values[0]):
                    Tm = tm_row['Tm_threshold'].values[0]
                    ax.axvline(Tm, color=colors[idx], linestyle='--', alpha=0.3)
                    ax.text(Tm, ax.get_ylim()[1]*0.95, f'Tm',
                           fontsize=8, color=colors[idx], rotation=90,
                           va='top', ha='right')
        
        # 添加阈值线
        ax.axhline(LINDEMANN_THRESHOLDS['solid'], color='blue',
                  linestyle='--', linewidth=2, alpha=0.5,
                  label=f'固态阈值 (δ={LINDEMANN_THRESHOLDS["solid"]})')
        ax.axhline(LINDEMANN_THRESHOLDS['melting'], color='red',
                  linestyle='--', linewidth=2, alpha=0.5,
                  label=f'熔化阈值 (δ={LINDEMANN_THRESHOLDS["melting"]})')
        
        # 标注区域
        ax.fill_between([200, 1100], 0, LINDEMANN_THRESHOLDS['solid'],
                       color='blue', alpha=0.1, label='固态区')
        ax.fill_between([200, 1100],
                       LINDEMANN_THRESHOLDS['solid'],
                       LINDEMANN_THRESHOLDS['melting'],
                       color='yellow', alpha=0.1, label='预熔区')
        ax.fill_between([200, 1100], LINDEMANN_THRESHOLDS['melting'], 0.4,
                       color='red', alpha=0.1, label='液态区')
        
        ax.set_xlabel('温度 (K)', fontsize=14)
        ax.set_ylabel('林德曼指数 δ', fontsize=14)
        ax.set_title(f'{series}系列: 林德曼指数 vs 温度', fontsize=16, fontweight='bold')
        ax.legend(fontsize=9, loc='upper left', ncol=2)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(150, 1150)
        ax.set_ylim(0, max(0.4, df_series['林德曼指数'].max() * 1.1))
        
        plt.tight_layout()
        filename = OUTPUT_DIR / f'Lindemann_vs_T_{series}.png'
        plt.savefig(filename, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  [OK] Saved: {filename.name}")


def plot_lindemann_heatmap(df_lindemann):
    """
    绘制林德曼指数热熔图
    """
    print("\n[*] Plotting Lindemann index heatmap...")
    
    # 按系列分组
    series_list = sorted(df_lindemann['系列'].unique())
    
    for series in series_list:
        df_series = df_lindemann[df_lindemann['系列'] == series].copy()
        
        if len(df_series) == 0:
            continue
        
        # 提取结构信息并排序(与step6保持一致)
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
        
        # 准备数据透视表
        pivot_data = df_series.groupby(['温度', '结构'])['林德曼指数'].mean().reset_index()
        heatmap_data = pivot_data.pivot(index='温度', columns='结构', values='林德曼指数')
        
        # 重新排列列顺序
        heatmap_data = heatmap_data[sorted_structures]
        
        # 创建图形
        fig, ax = plt.subplots(figsize=(14, 8))
        
        # 绘制热熔图
        im = ax.imshow(heatmap_data.values, aspect='auto', cmap=plt.cm.RdYlBu_r,
                      interpolation='bilinear', origin='lower',
                      vmin=0, vmax=0.3)
        
        # 设置坐标轴
        ax.set_xticks(np.arange(len(heatmap_data.columns)))
        ax.set_yticks(np.arange(len(heatmap_data.index)))
        
        # 生成x轴标签(同分异构体用*标记)
        x_labels = []
        formula_count = {}
        formula_index = {}
        
        for struct_name in heatmap_data.columns:
            row = df_series[df_series['结构'] == struct_name].iloc[0]
            pt_num = int(row['Pt原子数'])
            sn_num = int(row['Sn原子数'])
            formula = f'Pt{pt_num}Sn{sn_num}'
            formula_count[formula] = formula_count.get(formula, 0) + 1
        
        for struct_name in heatmap_data.columns:
            row = df_series[df_series['结构'] == struct_name].iloc[0]
            pt_num = int(row['Pt原子数'])
            sn_num = int(row['Sn原子数'])
            formula = f'Pt{pt_num}Sn{sn_num}'
            
            if formula_count[formula] > 1:
                idx = formula_index.get(formula, 0)
                formula_index[formula] = idx + 1
                markers = ['', '*', '**', '***', '****', '*****']  # 扩展标记列表
                # 确保idx不超出范围
                marker = markers[min(idx, len(markers)-1)]
                x_labels.append(f'{formula}{marker}')
            else:
                x_labels.append(formula)
        
        ax.set_xticklabels(x_labels, fontsize=10)
        ax.set_yticklabels([f'{int(temp)}K' for temp in heatmap_data.index], fontsize=10)
        
        # 旋转x轴标签
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
        
        # 在每个格子中标注数值
        for i in range(len(heatmap_data.index)):
            for j in range(len(heatmap_data.columns)):
                value = heatmap_data.values[i, j]
                if not np.isnan(value):
                    text_color = 'white' if value < 0.15 else 'black'
                    ax.text(j, i, f'{value:.3f}',
                           ha="center", va="center", color=text_color,
                           fontsize=8, fontweight='bold')
        
        # 设置标题和颜色条
        ax.set_xlabel('组分配比', fontsize=12)
        ax.set_ylabel('温度', fontsize=12)
        ax.set_title(f'{series}系列: 林德曼指数热熔图', fontsize=14, fontweight='bold')
        
        cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label('林德曼指数 δ', fontsize=12)
        
        # 添加阈值线
        threshold_y = np.where(heatmap_data.index >= 200)[0]
        if len(threshold_y) > 0:
            for threshold in [LINDEMANN_THRESHOLDS['solid'], LINDEMANN_THRESHOLDS['melting']]:
                # 在colorbar上标记阈值
                cbar.ax.axhline(threshold, color='black', linestyle='--', linewidth=2)
        
        plt.tight_layout()
        filename = OUTPUT_DIR / f'Lindemann_heatmap_{series}.png'
        plt.savefig(filename, dpi=200, bbox_inches='tight')
        plt.close()
        print(f"  [OK] Saved: {filename.name}")


def plot_melting_temperature_analysis(df_melting):
    """
    绘制熔化温度分析图
    """
    print("\n[*] Plotting melting temperature analysis...")
    
    # 筛选有Tm数据的结构
    df_valid = df_melting[df_melting['Tm_threshold'].notna()].copy()
    
    if len(df_valid) == 0:
        print("  [WARNING] No valid melting temperatures found")
        return
    
    # 创建图形: 2x2布局
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # 图1: Tm vs 总原子数
    ax1 = axes[0, 0]
    for series in sorted(df_valid['系列'].unique()):
        df_s = df_valid[df_valid['系列'] == series]
        ax1.scatter(df_s['总原子数'], df_s['Tm_threshold'],
                   s=100, alpha=0.7, label=series)
    
    ax1.set_xlabel('总原子数', fontsize=12)
    ax1.set_ylabel('熔化温度 Tm (K)', fontsize=12)
    ax1.set_title('熔化温度 vs 总原子数', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    
    # 图2: Tm vs Pt比例
    ax2 = axes[0, 1]
    df_valid['Pt比例'] = df_valid['Pt原子数'] / (df_valid['Pt原子数'] + df_valid['Sn原子数'])
    
    for series in sorted(df_valid['系列'].unique()):
        df_s = df_valid[df_valid['系列'] == series]
        ax2.scatter(df_s['Pt比例'], df_s['Tm_threshold'],
                   s=100, alpha=0.7, label=series)
    
    ax2.set_xlabel('Pt/(Pt+Sn) 比例', fontsize=12)
    ax2.set_ylabel('熔化温度 Tm (K)', fontsize=12)
    ax2.set_title('熔化温度 vs Pt比例', fontsize=14, fontweight='bold')
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)
    
    # 图3: Tm vs Sn含量
    ax3 = axes[1, 0]
    for series in sorted(df_valid['系列'].unique()):
        df_s = df_valid[df_valid['系列'] == series]
        ax3.scatter(df_s['Sn原子数'], df_s['Tm_threshold'],
                   s=100, alpha=0.7, label=series)
    
    ax3.set_xlabel('Sn原子数', fontsize=12)
    ax3.set_ylabel('熔化温度 Tm (K)', fontsize=12)
    ax3.set_title('熔化温度 vs Sn含量', fontsize=14, fontweight='bold')
    ax3.legend(fontsize=10)
    ax3.grid(True, alpha=0.3)
    
    # 图4: Tm vs O含量
    ax4 = axes[1, 1]
    for series in sorted(df_valid['系列'].unique()):
        df_s = df_valid[df_valid['系列'] == series]
        ax4.scatter(df_s['O原子数'], df_s['Tm_threshold'],
                   s=100, alpha=0.7, label=series)
    
    ax4.set_xlabel('O原子数', fontsize=12)
    ax4.set_ylabel('熔化温度 Tm (K)', fontsize=12)
    ax4.set_title('熔化温度 vs O含量', fontsize=14, fontweight='bold')
    ax4.legend(fontsize=10)
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    filename = OUTPUT_DIR / 'Melting_Temperature_Analysis.png'
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  [OK] Saved: {filename.name}")


def generate_report(df_lindemann, df_melting):
    """
    生成分析报告
    """
    print("\n[*] Generating analysis report...")
    
    report_lines = []
    report_lines.append("=" * 80)
    report_lines.append("林德曼指数分析报告")
    report_lines.append("=" * 80)
    report_lines.append(f"生成时间: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report_lines.append("")
    
    # 1. 数据概览
    report_lines.append("[1] 数据概览")
    report_lines.append("-" * 80)
    report_lines.append(f"总记录数: {len(df_lindemann)}")
    report_lines.append(f"体系数量: {df_lindemann['结构'].nunique()}")
    report_lines.append(f"温度范围: {df_lindemann['温度'].min():.0f} - {df_lindemann['温度'].max():.0f} K")
    report_lines.append(f"林德曼指数范围: {df_lindemann['林德曼指数'].min():.4f} - {df_lindemann['林德曼指数'].max():.4f}")
    report_lines.append("")
    
    # 2. 系列统计
    report_lines.append("[2] 系列统计")
    report_lines.append("-" * 80)
    series_stats = df_lindemann.groupby('系列').agg({
        '结构': 'nunique',
        '林德曼指数': ['mean', 'min', 'max']
    })
    for series in sorted(df_lindemann['系列'].unique()):
        stats = series_stats.loc[series]
        n_struct = int(stats[('结构', 'nunique')])
        mean_l = stats[('林德曼指数', 'mean')]
        min_l = stats[('林德曼指数', 'min')]
        max_l = stats[('林德曼指数', 'max')]
        report_lines.append(f"{series}: {n_struct} 结构, δ平均={mean_l:.4f}, 范围=[{min_l:.4f}, {max_l:.4f}]")
    report_lines.append("")
    
    # 3. 熔化温度统计
    report_lines.append("[3] 熔化温度统计")
    report_lines.append("-" * 80)
    df_valid_tm = df_melting[df_melting['Tm_threshold'].notna()]
    report_lines.append(f"有效熔化温度: {len(df_valid_tm)}/{len(df_melting)} 结构")
    report_lines.append(f"平均Tm: {df_valid_tm['Tm_threshold'].mean():.1f} K")
    report_lines.append(f"Tm范围: {df_valid_tm['Tm_threshold'].min():.1f} - {df_valid_tm['Tm_threshold'].max():.1f} K")
    report_lines.append("")
    
    # 各系列Tm
    for series in sorted(df_valid_tm['系列'].unique()):
        df_s = df_valid_tm[df_valid_tm['系列'] == series]
        report_lines.append(f"  {series}: Tm平均={df_s['Tm_threshold'].mean():.1f} K, "
                          f"范围=[{df_s['Tm_threshold'].min():.1f}, {df_s['Tm_threshold'].max():.1f}] K, "
                          f"结构数={len(df_s)}")
    report_lines.append("")
    
    # 4. 各结构详细信息
    report_lines.append("[4] 各结构林德曼指数统计")
    report_lines.append("-" * 80)
    report_lines.append(f"{'结构':<18} {'系列':<8} {'总原子数':<8} {'Tm(K)':<10} {'δ平均':<10} {'δ范围'}")
    report_lines.append("-" * 80)
    
    for _, row in df_melting.sort_values(['系列', '总原子数']).iterrows():
        struct = row['结构']
        series = row['系列']
        total = int(row['总原子数'])
        Tm = row['Tm_threshold']
        
        df_struct = df_lindemann[df_lindemann['结构'] == struct]
        mean_l = df_struct['林德曼指数'].mean()
        min_l = df_struct['林德曼指数'].min()
        max_l = df_struct['林德曼指数'].max()
        
        Tm_str = f"{Tm:.1f}" if pd.notna(Tm) else "未熔化"
        
        report_lines.append(f"{struct:<18} {series:<8} {total:<8} {Tm_str:<10} {mean_l:<10.4f} [{min_l:.4f}, {max_l:.4f}]")
    
    report_lines.append("")
    report_lines.append("=" * 80)
    
    # 保存报告
    report_file = OUTPUT_DIR / 'lindemann_analysis_report.txt'
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write('\n'.join(report_lines))
    
    print(f"  [OK] Saved: {report_file.name}")
    
    # 打印到控制台
    print("\n" + '\n'.join(report_lines))


def main():
    """
    主函数
    
    命令行参数:
    -----------
    --no-filter : 不使用MSD异常筛选,分析所有Lindemann数据
                  (适用于当某些体系的数据全被筛选掉时)
    --exclude : 排除指定的结构,支持多个,例如: --exclude Pt8Sn0 Pt8Sn1
                (适用于绘制时不想包含某些特定结构)
    
    示例:
    python step7_lindemann_analysis.py                    # 默认:使用异常筛选
    python step7_lindemann_analysis.py --no-filter        # 不筛选,使用全部数据
    python step7_lindemann_analysis.py --exclude Pt8Sn0   # 排除Pt8Sn0
    python step7_lindemann_analysis.py --exclude Pt8Sn0 Pt8Sn1  # 排除多个结构
    """
    # 解析命令行参数
    parser = argparse.ArgumentParser(
        description='Step7: 林德曼指数分析 (改进版)'
    )
    parser.add_argument(
        '--no-filter',
        action='store_true',
        help='不使用MSD异常筛选,分析所有Lindemann数据'
    )
    parser.add_argument(
        '--exclude',
        nargs='+',
        default=[],
        help='排除指定的结构,支持多个结构名称'
    )
    args = parser.parse_args()
    
    # 1. 加载和筛选数据 (传入no_filter参数)
    df_lindemann = load_and_filter_data(no_filter=args.no_filter)
    
    # 1.5. 应用结构排除过滤
    if args.exclude:
        excluded_structures = args.exclude
        original_count = len(df_lindemann)
        df_lindemann = df_lindemann[~df_lindemann['结构'].isin(excluded_structures)]
        filtered_count = original_count - len(df_lindemann)
        
        print(f"\n{'='*80}")
        print(f"应用结构过滤")
        print(f"{'='*80}")
        print(f"  排除的结构: {', '.join(excluded_structures)}")
        print(f"  过滤前数据点: {original_count}")
        print(f"  过滤后数据点: {len(df_lindemann)}")
        print(f"  已删除数据点: {filtered_count}")
        
        if len(df_lindemann) == 0:
            print("\n[X] 错误: 过滤后没有数据剩余!")
            return
    
    # 2. 分析熔化温度
    df_melting = analyze_melting_temperatures(df_lindemann)
    
    # 2.5. Cv系列专项分析
    df_cv_stats = analyze_cv_series_lindemann(df_lindemann)
    
    # 3. 绘制林德曼指数 vs 温度
    plot_lindemann_vs_temperature(df_lindemann, df_melting)
    
    # 4. 绘制林德曼指数热熔图
    plot_lindemann_heatmap(df_lindemann)
    
    # 5. 绘制熔化温度分析
    plot_melting_temperature_analysis(df_melting)
    
    # 6. 生成报告
    generate_report(df_lindemann, df_melting)
    
    print("\n" + "=" * 80)
    print(f"[DONE] Analysis complete! Results saved to: {OUTPUT_DIR}")
    print("=" * 80)


if __name__ == '__main__':
    main()
