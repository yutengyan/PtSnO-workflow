#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
Step 7.4: 多体系单次运行热容分析 (Multi-System Heat Capacity Analysis)
================================================================================

作者: GitHub Copilot
日期: 2025-10-22
版本: v1.2

最近更新:
---------
v1.3 (2025-10-22):
  - 🔄 **重要修复**: Cv体系数据合并
    · 问题: Cv-1, Cv-2, Cv-3, Cv-4, Cv-5被错误地识别为5个独立结构
    · 实际: 这些是同一个体系(Sn8Pt6O4)的5次重复模拟
    · 修复: detect_system_type()函数统一返回('Cv', 'Cv')
    · 结果: Cv数据从5×19=95个点合并为单一体系(95个点,19个温度×5次重复)
    · 影响: step7_4输出中Cv现为1个结构(之前5个), step7_4_2聚类分析使用Cv而非Cv-1~5

v1.2 (2025-10-22):
  - 🔧 修复: 统一match_key生成算法,使用4-level路径签名 (与MSD筛选一致)
  - 📊 修复: 数据量从错误的8030→5296更正为正确的3262→2370
  - ✅ 验证: 筛选率从错误的34.05%更正为正确的27.35%
  - 📝 增强: 添加详细的数据来源和去重说明

v1.1 (2025-10-21):
  - ✨ 新增: Lindemann数据去重功能 (4012→3262条)
  - 📈 新增: 筛选统计报告自动生成
  - 🎯 新增: K-means聚类分析 (可选)
  - 🐛 修复: Unicode编码问题 (Windows兼容性)

================================================================================
功能概述
================================================================================
本脚本用于分析多个纳米团簇体系的热容,基于LAMMPS能量输出和Lindemann指数,
支持单次运行级别的详细分析,并提供多种数据质量控制方法。

**核心特性**:
1. ✅ **多体系支持**: Cv、Pt8SnX、Pt6SnX、Pt8-xSnx、PtxSny、PtxSnyOz (6大系列, 51个结构)
   - 注: Cv-1~Cv-5自动合并为单一Cv体系(5次重复模拟,95个数据点)
2. ✅ **智能体系识别**: 自动分类和命名标准化 (system_type + system_id)
   - 特殊处理: Cv-X → ('Cv', 'Cv') 统一命名,避免重复模拟被拆分
3. ✅ **路径签名筛选**: 继承Step 6的4级路径签名算法 (与MSD筛选完全一致)
4. ✅ **双重数据质量控制**:
   - Method 1 (opt-in): Step 1 MSD异常筛选 (892个路径签名, 27.35%移除率)
   - Method 2 (opt-in): IQR统计异常筛选 (Lindemann + Energy, ~3.5%额外移除)
5. ✅ **单结构细粒度分析**: 每个结构独立分析 (system_id级别), 生成55个独立图表
6. ✅ **三区域热容拟合**: Solid / Premelting / Liquid区域独立线性拟合
7. ✅ **聚类分析** (可选): K-means自动检测相边界 (针对特定体系如Pt6Sn8)
8. ✅ **筛选透明度**: 自动生成详细的筛选统计报告

================================================================================
输入文件
================================================================================
**必需输入** (自动加载):

1. **能量数据** (LAMMPS总能量)
   路径: data/lammps_energy/energy_master_20251016_121110.csv
   来源: LAMMPS MD模拟的etotal输出 (团簇 + 载体总能量)
   列: path, structure, temp, run_num, total_steps, sample_steps,
       avg_energy, std, min, max, sample_interval, skip_steps, full_path
   统计: 3262条记录 → 2370条 (应用--msd-filter后)
   路径示例: /home/scms/jychen/.../pt8sn8-1-best/T900.r15.gpu0

2. **Lindemann指数数据** (相态分类)
   路径: data/lindemann/lindemann_master_run_*.csv (3个文件)
   来源: MSD轨迹计算的Lindemann指数 (基于Pt-Sn距离MSD)
   列: 目录, 结构, 温度(K), Lindemann指数, 方法, 耗时(s), 时间戳
   统计: 
     - 原始合并: 4012条 (3个CSV文件)
     - 去重后: 3262条 (移除750条重复, 18.7%)
     - 筛选后: 2370条 (应用--msd-filter后)
   路径示例: /home/scms/jychen/.../pt8sn8-1-best/T900.r15.gpu0
   注意: 去重基于(目录, 结构, 温度)三元组, 保留最新记录

3. **Step 1 MSD异常过滤数据** (数据质量控制)
   路径: results/large_D_outliers.csv
   来源: Step 1的Pt和Sn元素MSD异常检测
   列: group_key, composition, temperature, element, run_id, 
       gmx_D, filepath, reason
   统计: 
     - 2227条异常记录
     - 892个唯一路径签名 (4-level)
   异常原因分布:
     - linked_bad_simulation: 804条 (36.1%)
     - IQR_outlier: 717条 (32.2%)
     - Intercept>10.0A²: 706条 (31.7%)
   关键: 基于Pt/Sn **元素MSD**, 非Pt-Sn距离MSD

4. **载体热容数据** (可选,用于计算纯团簇热容)
   路径: data/lammps_energy/sup/energy_master_20251021_151520.csv
   来源: 纯Al₂O₃载体的LAMMPS模拟
   用途: 从总热容中扣除载体贡献
   默认值: 38.2151 meV/K (如果文件不存在)

================================================================================
路径签名算法 (v1.2 关键修复)
================================================================================
**关键修复 (v1.2)**: `normalize_path`和`extract_path_signature`现在使用**相同**的算法!

**问题背景** (v1.1及之前):
- `extract_path_signature`: 4-level签名 → 用于MSD筛选
- `normalize_path`: 简单的结构+温度+run → 用于Energy-Lindemann合并
- **后果**: 导致数据量虚高 (8030条) 和筛选率错误 (34.05%)

**修复方案 (v1.2)**:
```python
def normalize_path(path, ...):
    # 现在直接调用extract_path_signature!
    is_msd_path = 'msd_' in path.lower() or re.search(r'\\d+K[/\\\\]', path)
    return extract_path_signature(path, is_msd_path=is_msd_path)
```

**4级路径签名算法** (统一使用):
提取格式: batch/parent/composition/run
  MSD路径: run3/o2/o2pt4sn6/t1000.r24.gpu0
  能量路径: run3/o2/o2pt7sn7/t200.r0.gpu0
  
签名构建步骤:
1. 提取run信息: T\\d+\\.r\\d+\\.gpu\\d+ → 小写化 → t1000.r24.gpu0
2. 向上提取composition目录: O2Pt4Sn6 → 小写化 → o2pt4sn6
3. 再向上提取parent目录: o2 → 小写化
4. 检测batch标识符: run3, run2, run4, run5 (向上最多3级)
5. 构建4级或3级签名 (取决于是否检测到批次)

**路径兼容性**:
- Windows路径: 支持反斜杠 (\\)
- Linux路径: 支持正斜杠 (/)
- MSD路径: 识别温度目录 (1000K/)
- Energy/Lindemann路径: 直接提取run信息

**筛选效果** (v1.2修复后):
- 能量数据: 3262 → 2370 (筛除892条, 27.3%) ✅ 正确
- Lindemann数据: 3262 → 2370 (筛除892条, 27.3%) ✅ 正确
- 最终合并: 3262条 (1:1匹配) → 2370条 ✅ 正确

================================================================================
输出文件
================================================================================
输出目录: results/step7_4_multi_system/

**1. 主数据文件**:
   - step7_4_all_systems_data.csv
     内容: 所有匹配的能量-Lindemann数据
     列: match_key, structure, system_type, system_id, temp, 
         avg_energy, energy_std, delta, phase, run_id
     统计: 2370条 (--msd-filter) 或 2287条 (--msd-filter --iqr-filter)

**2. 分析报告**:
   - step7_4_multi_system_comparison.md
     内容: 55个结构的热容分析汇总
     包含: 固态/预熔/液态三区域热容、R²值、温度范围

**3. 综合对比图**:
   - step7_4_multi_system_comparison.png
     内容: 6个体系类型的热容对比 (4×3子图布局)
     子图: 
       - 固态Cv vs Sn含量
       - 预熔Cv vs Sn含量  
       - 液态Cv vs Sn含量
       - 熔化温度 vs Sn含量

**4. 单结构详细分析图** (55个PNG文件):
   目录: individual_structure_plots/
   命名: {structure_name}_individual_runs_analysis.png
   示例: 
     - Cv-1_individual_runs_analysis.png
     - Pt8Sn0_individual_runs_analysis.png
     - Pt6Sn3_individual_runs_analysis.png
     - O2Pt4Sn6_individual_runs_analysis.png
   
   图表布局 (2×2):
     (a) 相对能量 vs 温度
         - 散点: 单次运行数据点 (按相态着色)
         - 线: 固态/预熔/液态三区域拟合
         - 标注: 斜率(热容)、R²、温度范围
     
     (b) 热容柱状图
         - 三柱: 固态、预熔、液态热容
         - 载体线: Cv_support参考值
         - 标注: 数值、误差、R²值
     
     (c) 温度-相态分布图 (堆叠柱状图)
         - Y轴: 温度 (200-1100K)
         - 颜色: 固态(蓝)、预熔(橙)、液态(红)
         - 显示: 每个温度的相态分布
     
     (d) Lindemann散点图
         - 散点: Lindemann指数 vs 温度
         - 阈值线: 0.1 (固态上限)、0.15 (液态下限)
         - 区域: 固态区(蓝)、预熔区(橙)、液态区(红)

================================================================================
分析算法
================================================================================
**三区域线性拟合** (固态 → 预熔 → 液态):

1. 相态分类 (基于Lindemann指数):
   - 固态: δ < 0.1
   - 预熔: 0.1 ≤ δ < 0.15
   - 液态: δ ≥ 0.15

2. 温度区间识别 (每个相态):
   - 自动检测连续温度范围
   - 最小数据点: 5个 (保证统计可靠性)

3. 线性拟合 (每个区域):
   E_total(T) = slope × T + intercept
   Cv_total = slope × 1000 (eV/K → meV/K)
   Cv_cluster = Cv_total - Cv_support

4. 统计评估:
   - R²值: 拟合优度 (>0.995为优秀)
   - p值: 统计显著性 (<0.05为显著)
   - 标准误差: 参数不确定性

================================================================================
体系分类策略
================================================================================
自动识别6大类55个结构:

1. **Cv系列** (5个结构):
   - 匹配模式: ^Cv-\\d+
   - 示例: Cv-1, Cv-2, Cv-3, Cv-4, Cv-5
   - 特点: 参考体系,固定成分

2. **Pt6系列** (1个结构):
   - 匹配模式: ^pt6$
   - 示例: Pt6 (纯Pt)

3. **Pt6SnX系列** (9个结构):
   - 匹配模式: ^pt6sn\\d+
   - 示例: Pt6Sn1, Pt6Sn2, ..., Pt6Sn9
   - 特点: 固定6个Pt,变Sn含量

4. **Pt8SnX系列** (11个结构):
   - 匹配模式: ^pt8sn\\d+
   - 示例: Pt8Sn0, Pt8Sn1, ..., Pt8Sn10
   - 特点: 固定8个Pt,变Sn含量

5. **PtxSny系列** (4个结构):
   - 匹配模式: ^pt\\d+sn\\d+$ (无O)
   - 示例: Pt3Sn5, Pt4Sn4, Pt5Sn3, Pt7Sn1
   - 特点: 变Pt/Sn比,无氧

6. **PtxSnyOz系列** (25个结构):
   - 匹配模式: 包含O的三元系统
   - 示例: O2Pt4Sn6, Pt7Sn5O1, Sn10Pt7O4等
   - 特点: Pt-Sn-O三元合金

================================================================================
使用说明
================================================================================
**命令行参数**:
    --msd-filter         : 开启Step 1的MSD异常值筛选(路径签名匹配)
    --iqr-filter         : 开启IQR统计筛选(Lindemann指数 + 能量值)
    --iqr-factor FLOAT   : IQR倍数阈值(默认3.0,更严格则用更大值)

**筛选策略**:
    方法1 (可选): Step 1 MSD异常值筛选
      - 基于路径签名匹配 (4级签名: batch/parent/composition/run)
      - 892个异常模拟路径签名
      - 典型筛除27-28%数据
      - 来源: large_D_outliers.csv
      - 使用方法: --msd-filter

    方法2 (可选): IQR统计筛选 (Lindemann + Energy)
      - 按(结构,温度)分组进行统计分析
      - **Lindemann IQR**: 超出 [Q1-k*IQR, Q3+k*IQR] 范围的δ值
      - **Energy IQR**: 超出 [Q1-k*IQR, Q3+k*IQR] 范围的能量值
      - 可调阈值k (默认3.0*IQR, 严于标准1.5*IQR)
      - 要求每组 ≥5 个数据点才进行检测
      - 使用方法: --iqr-filter

**运行示例**:
    # 示例1: 默认运行 (无筛选,使用全部原始数据)
    python step7_4_multi_system_heat_capacity.py

    # 示例2: 仅方法1筛选
    python step7_4_multi_system_heat_capacity.py --msd-filter

    # 示例3: 仅方法2统计筛选 (Lindemann + Energy)
    python step7_4_multi_system_heat_capacity.py --iqr-filter

    # 示例4: 双重筛选 (方法1 + 方法2)
    python step7_4_multi_system_heat_capacity.py --msd-filter --iqr-filter

    # 示例5: 自定义IQR阈值 (更宽松的筛选)
    python step7_4_multi_system_heat_capacity.py --iqr-filter --iqr-factor 5.0

**输出确认**:
    - 终端: 显示55个结构的分析进度和统计
    - 文件: 检查 results/step7_4_multi_system/ 目录

**典型数据量** (默认方法1筛选):
    - 筛选前: 10355条匹配记录 (未过滤)
    - 筛选后: 6814条匹配记录 (应用Step 1过滤)
    - 质量提升: 34.2% (移除异常模拟)
    - 不同筛选策略会产生不同数据量

**注意事项**:
    1. 热容值为总热容 (团簇 + 载体)
    2. 若需纯团簇热容,需扣除Cv_support
    3. 载体热容默认38.2151 meV/K (如无载体数据)
    4. R² > 0.995 视为优秀拟合
    5. 数据点 < 5 的相态区域会跳过拟合
    6. 筛选策略选择建议:
       - 默认方法1: 适用于大多数情况
       - 双重筛选: 最严格的质量控制
       - 无筛选: 用于对比验证筛选效果
       - 仅方法2: 纯统计方法,独立于Step 1

================================================================================
依赖关系
================================================================================
上游依赖:
    - Step 1: large_D_outliers.csv (异常值筛选)
    - LAMMPS: energy_master_*.csv (能量数据)
    - MSD分析: lindemann_master_run_*.csv (Lindemann指数)

下游应用:
    - Step 7.5: 热容系统性分析
    - 科学发现: Sn含量对热容的影响
    - 相变研究: 熔化温度随成分变化

================================================================================
版本历史
================================================================================
v1.1 (2025-10-21):
    - ✅ 添加Step 1路径签名筛选 (继承Step 6算法)
    - ✅ 筛除892个异常模拟,数据质量提升34.2%
    - ✅ 详细文档化输入输出和算法原理

v1.0 (2025-10-20):
    - ✅ 初始版本,支持55个结构分析
    - ✅ 三区域热容拟合算法
    - ✅ 生成个体分析图和综合报告

================================================================================
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import linregress, iqr
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from matplotlib import rcParams
import warnings
import re
import argparse
from datetime import datetime
warnings.filterwarnings('ignore')

# Chinese font settings
rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'Arial Unicode MS']
rcParams['axes.unicode_minus'] = False

# Lindemann thresholds
LINDEMANN_THRESHOLDS = {
    'solid': 0.1,
    'melting': 0.15
}

# File paths
BASE_DIR = Path(__file__).parent
CLUSTER_ENERGY_FILE = BASE_DIR / 'data' / 'lammps_energy' / 'energy_master_20251016_121110.csv'
SUPPORT_ENERGY_FILE = BASE_DIR / 'data' / 'lammps_energy' / 'sup' / 'energy_master_20251021_151520.csv'
LINDEMANN_DIR = BASE_DIR / 'data' / 'lindemann'
RESULTS_DIR = BASE_DIR / 'results' / 'step7_4_multi_system'

# Step 1 filtering results
OUTLIERS_FILE = BASE_DIR / 'results' / 'large_D_outliers.csv'


def detect_system_type(structure_name):
    """
    Detect system type from structure name
    
    Returns:
        tuple: (system_type, system_id)
        
    Examples:
        'Cv-1' -> ('Cv', 'Cv')  # All Cv-X runs merged into single 'Cv' system
        'pt8sn0-2-best' -> ('Pt8SnX', 'Pt8Sn0')
        'pt6sn3' -> ('Pt6SnX', 'Pt6Sn3')
        'pt3sn5' -> ('PtxSny', 'Pt3Sn5')
        'Pt5Sn3O1' -> ('PtxSnyOz', 'Pt5Sn3O1')
    """
    name = structure_name.strip()
    
    # Cv series - MERGE ALL Cv-1, Cv-2, ..., Cv-5 into single 'Cv' system
    # These are repeat simulations of the same structure (Sn8Pt6O4)
    if re.match(r'^Cv-\d+', name, re.IGNORECASE):
        return ('Cv', 'Cv')  # Use 'Cv' as unified system_id
    
    # Pt8SnX series (pt8sn0-2-best, pt8sn1-2-best, etc.)
    if re.match(r'^pt8sn\d+', name, re.IGNORECASE):
        match = re.match(r'^(pt8sn\d+)', name, re.IGNORECASE)
        base_name = match.group(1)
        return ('Pt8SnX', base_name.capitalize())
    
    # Pt6SnX series (pt6sn1, pt6sn2, etc.)
    if re.match(r'^pt6sn\d+', name, re.IGNORECASE):
        match = re.match(r'^(pt6sn\d+)', name, re.IGNORECASE)
        base_name = match.group(1)
        return ('Pt6SnX', base_name.capitalize())
    
    # Pt6 alone
    if re.match(r'^pt6$', name, re.IGNORECASE):
        return ('Pt6', 'Pt6')
    
    # PtxSny series (pt3sn5, pt4sn4, pt5sn3, pt7sn1, etc.)
    # Match ptXsnY where X and Y are numbers
    if re.match(r'^pt\d+sn\d+(?!o)', name, re.IGNORECASE):
        match = re.match(r'^(pt\d+sn\d+)', name, re.IGNORECASE)
        base_name = match.group(1)
        return ('PtxSny', base_name.capitalize())
    
    # PtxSnyOz series (with oxygen)
    if re.search(r'o\d*', name, re.IGNORECASE):
        # Extract composition
        match = re.match(r'^([a-z0-9]+)', name, re.IGNORECASE)
        base_name = match.group(1) if match else name
        return ('PtxSnyOz', base_name.capitalize())
    
    # Other types
    return ('Other', name)


def extract_path_signature(filepath, is_msd_path=True):
    """
    从文件路径提取路径签名 (与Step 6保持一致)
    
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
        # 能量/Lindemann路径: 找run所在位置
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


def load_outliers():
    """
    Load Step 1 outliers and build path signature filter set
    
    Returns:
        set: Set of path signatures to exclude
    """
    if not OUTLIERS_FILE.exists():
        print(f"\n>>> Warning: Outliers file not found: {OUTLIERS_FILE}")
        print("    Continuing without Step 1 filtering")
        return set()
    
    print(f"\n>>> Loading Step 1 outliers: {OUTLIERS_FILE.name}")
    df_outliers = pd.read_csv(OUTLIERS_FILE, encoding='utf-8')
    
    print(f"    Loaded: {len(df_outliers)} outlier records")
    print(f"    Reasons breakdown:")
    reason_counts = df_outliers['reason'].value_counts()
    for reason, count in reason_counts.items():
        # Use ASCII-safe printing to avoid encoding issues on Windows
        try:
            print(f"      - {reason}: {count}")
        except UnicodeEncodeError:
            print(f"      - [encoding issue]: {count}")
    
    # Build path signature filter set (same as Step 6)
    filter_signatures = set()
    for _, row in df_outliers.iterrows():
        filepath = row.get('filepath', '')
        if not filepath:
            continue
        
        # MSD paths have temperature directories (e.g., .../1000K/T1000.r24.gpu0_msd_Pt.xvg)
        path_signature = extract_path_signature(filepath, is_msd_path=True)
        if path_signature:
            filter_signatures.add(path_signature)
    
    print(f"    Built path signature filter set:")
    print(f"      - Unique signatures: {len(filter_signatures)}")
    if len(filter_signatures) > 0:
        print(f"      - Sample signatures:")
        for idx, sig in enumerate(sorted(list(filter_signatures))[:3]):
            print(f"        {idx+1}. {sig}")
    
    return filter_signatures


def detect_lindemann_outliers_iqr(df_merged):
    """
    基于 Lindemann 指数的 IQR 异常值检测
    
    算法原理:
    1. 按 (composition, temperature) 分组
    2. 计算每组的 Lindemann 指数 IQR (四分位距)
    3. 识别异常值: Q1 - 3*IQR 或 Q3 + 3*IQR 之外的点
    
    Args:
        df_merged: 合并后的数据 (包含 structure, temp, delta 列)
    
    Returns:
        set: 异常记录的 match_key 集合
    """
    print(f"\n>>> Detecting Lindemann IQR outliers")
    
    outlier_keys = set()
    total_groups = 0
    groups_with_outliers = 0
    total_outliers = 0
    
    # 按 (structure, temp) 分组分析
    for (structure, temp), group in df_merged.groupby(['structure', 'temp']):
        total_groups += 1
        
        # 至少需要5个数据点才能进行IQR分析
        if len(group) < 5:
            continue
        
        delta_values = group['delta'].values
        
        # 计算 IQR
        Q1 = np.percentile(delta_values, 25)
        Q3 = np.percentile(delta_values, 75)
        IQR = Q3 - Q1
        
        # 定义异常值边界 (使用3倍IQR,比1.5倍更严格)
        lower_bound = Q1 - 3.0 * IQR
        upper_bound = Q3 + 3.0 * IQR
        
        # 识别异常值
        outlier_mask = (delta_values < lower_bound) | (delta_values > upper_bound)
        
        if outlier_mask.any():
            groups_with_outliers += 1
            group_outliers = group[outlier_mask]
            total_outliers += len(group_outliers)
            
            # 添加到异常集合
            outlier_keys.update(group_outliers['match_key'].values)
            
            # 详细日志 (仅显示前5个异常组)
            if groups_with_outliers <= 5:
                print(f"    {structure} @ {temp}K: {len(group_outliers)} outliers")
                print(f"      IQR range: [{lower_bound:.4f}, {upper_bound:.4f}]")
                print(f"      Outlier δ: {group_outliers['delta'].values}")
    
    print(f"\n    Summary:")
    print(f"      Total groups analyzed: {total_groups}")
    print(f"      Groups with outliers: {groups_with_outliers}")
    print(f"      Total outlier records: {total_outliers}")
    print(f"      Unique match_keys: {len(outlier_keys)}")
    
    return outlier_keys


def detect_energy_outliers_iqr(df_merged, iqr_factor=3.0):
    """
    基于能量值的 IQR 异常值检测
    
    算法原理:
    1. 按 (composition, temperature) 分组
    2. 计算每组的能量值 IQR (四分位距)
    3. 识别异常值: Q1 - k*IQR 或 Q3 + k*IQR 之外的点
    
    Args:
        df_merged: 合并后的数据 (包含 structure, temp, energy_cluster 列)
        iqr_factor: IQR倍数 (默认3.0, 更严格则用更大值)
    
    Returns:
        set: 异常记录的 match_key 集合
    """
    print(f"\n>>> Detecting Energy IQR outliers")
    
    outlier_keys = set()
    total_groups = 0
    groups_with_outliers = 0
    total_outliers = 0
    
    # 按 (structure, temp) 分组分析
    for (structure, temp), group in df_merged.groupby(['structure', 'temp']):
        total_groups += 1
        
        # 至少需要5个数据点才能进行IQR分析
        if len(group) < 5:
            continue
        
        energy_values = group['energy_cluster'].values
        
        # 计算 IQR
        Q1 = np.percentile(energy_values, 25)
        Q3 = np.percentile(energy_values, 75)
        IQR = Q3 - Q1
        
        # 定义异常值边界
        lower_bound = Q1 - iqr_factor * IQR
        upper_bound = Q3 + iqr_factor * IQR
        
        # 识别异常值
        outlier_mask = (energy_values < lower_bound) | (energy_values > upper_bound)
        
        if outlier_mask.any():
            groups_with_outliers += 1
            group_outliers = group[outlier_mask]
            total_outliers += len(group_outliers)
            
            # 添加到异常集合
            outlier_keys.update(group_outliers['match_key'].values)
            
            # 详细日志 (仅显示前5个异常组)
            if groups_with_outliers <= 5:
                print(f"    {structure} @ {temp}K: {len(group_outliers)} outliers")
                print(f"      IQR range: [{lower_bound:.2f}, {upper_bound:.2f}] eV")
                print(f"      Outlier E: {group_outliers['energy_cluster'].values}")
    
    print(f"\n    Summary:")
    print(f"      Total groups analyzed: {total_groups}")
    print(f"      Groups with outliers: {groups_with_outliers}")
    print(f"      Total outlier records: {total_outliers}")
    print(f"      Unique match_keys: {len(outlier_keys)}")
    
    return outlier_keys


def normalize_path(path, system_type='Cv', structure_name=''):
    """
    DEPRECATED: Use extract_path_signature instead!
    
    This function is kept for backward compatibility but should use
    the same 4-level path signature algorithm as extract_path_signature
    to ensure consistency with MSD outlier filtering.
    
    Args:
        path: file path string
        system_type: system type (ignored, kept for compatibility)
        structure_name: structure name (ignored, kept for compatibility)
    
    Returns:
        str: 4-level path signature like "run3/o2/o2pt4sn6/t1000.r24.gpu0"
    """
    # Use the same algorithm as Step 6 MSD outlier filtering!
    # Determine if this is an MSD path or energy/Lindemann path
    is_msd_path = 'msd_' in path.lower() or re.search(r'\d+K[/\\]', path) is not None
    
    return extract_path_signature(path, is_msd_path=is_msd_path)


def load_energy_data(energy_file, system_filter=None, file_type='cluster'):
    """
    Load energy data at individual run level
    
    Args:
        energy_file: path to energy CSV file
        system_filter: list of system types to filter (e.g., ['Cv', 'Pt8SnX'])
                      None means load all systems
        file_type: 'cluster' or 'support'
    """
    print(f"\n>>> Loading {file_type} energy data: {energy_file.name}")
    df = pd.read_csv(energy_file, encoding='utf-8')
    
    # Handle different column formats
    if '结构' in df.columns:
        df.columns = ['path', 'structure', 'temp', 'run_num', 'total_steps', 'sample_steps', 
                      'avg_energy', 'std', 'min', 'max', 'sample_interval', 
                      'skip_steps', 'full_path']
    else:
        # Already in English
        expected_cols = ['path', 'structure', 'temp', 'run_num', 'total_steps', 'sample_steps', 
                        'avg_energy', 'std', 'min', 'max', 'sample_interval', 
                        'skip_steps', 'full_path']
        df.columns = expected_cols[:len(df.columns)]
    
    # Detect system types
    df[['system_type', 'system_id']] = df['structure'].apply(
        lambda x: pd.Series(detect_system_type(x))
    )
    
    # Filter by system if specified
    if system_filter:
        df = df[df['system_type'].isin(system_filter)].copy()
        print(f"    Filter applied: {system_filter}")
    
    # Create matching key
    df['match_key'] = df.apply(
        lambda row: normalize_path(row['full_path'], row['system_type'], row['structure']), axis=1
    )
    df = df[df['match_key'].notna()]
    
    # Print system distribution
    system_counts = df['system_type'].value_counts()
    print(f"    Loaded: {len(df)} records")
    print(f"    System distribution:")
    for sys_type, count in system_counts.items():
        print(f"      - {sys_type}: {count} records")
    
    return df


def load_lindemann_individual_runs(system_filter=None):
    """
    Load Lindemann raw data at individual run level
    
    Args:
        system_filter: list of system types to filter
    """
    print(f"\n>>> Loading Lindemann individual run data")
    
    # Find all lindemann_master_run_*.csv files
    lindemann_files = sorted(LINDEMANN_DIR.glob('lindemann_master_run_*.csv'))
    
    if not lindemann_files:
        print(f"    Error: No Lindemann files found in {LINDEMANN_DIR}")
        return None
    
    print(f"    Found {len(lindemann_files)} files")
    
    # Read and merge
    dfs = []
    for f in lindemann_files:
        df_temp = pd.read_csv(f, encoding='utf-8')
        dfs.append(df_temp)
        print(f"      - {f.name}: {len(df_temp)} records")
    
    df = pd.concat(dfs, ignore_index=True)
    print(f"    Total concatenated: {len(df)} records")
    
    # Remove duplicates (based on directory + structure + temperature)
    # Keep the latest record (assuming later files have more recent data)
    key_cols = ['目录', '结构', '温度(K)'] if '目录' in df.columns else ['directory', 'structure', 'temp']
    df_before_dedup = len(df)
    df = df.drop_duplicates(subset=key_cols, keep='last')
    duplicates_removed = df_before_dedup - len(df)
    
    if duplicates_removed > 0:
        print(f"    [OK] Removed {duplicates_removed} duplicate records ({duplicates_removed/df_before_dedup*100:.1f}%)")
        print(f"    After deduplication: {len(df)} unique records")
    
    # Map Chinese column names
    col_mapping = {
        '目录': 'directory',
        '结构': 'structure',
        '温度(K)': 'temp',
        'Lindemann指数': 'delta'
    }
    
    # Handle both Chinese and English column names
    if '目录' in df.columns:
        df.rename(columns=col_mapping, inplace=True)
    
    # Detect system types
    df[['system_type', 'system_id']] = df['structure'].apply(
        lambda x: pd.Series(detect_system_type(x))
    )
    
    # Filter by system if specified
    if system_filter:
        df = df[df['system_type'].isin(system_filter)].copy()
        print(f"    Filter applied: {system_filter}")
    
    # Create matching key
    df['match_key'] = df.apply(
        lambda row: normalize_path(row['directory'], row['system_type'], row['structure']), axis=1
    )
    df = df[df['match_key'].notna()]
    
    # Print system distribution
    system_counts = df['system_type'].value_counts()
    print(f"    Loaded: {len(df)} records")
    print(f"    System distribution:")
    for sys_type, count in system_counts.items():
        print(f"      - {sys_type}: {count} records")
    
    return df


def classify_single_run(delta):
    """Simple classification with hard thresholds"""
    if delta < LINDEMANN_THRESHOLDS['solid']:
        return 'solid'
    elif delta < LINDEMANN_THRESHOLDS['melting']:
        return 'premelting'
    else:
        return 'liquid'


def merge_energy_lindemann(df_energy, df_lindemann, outlier_signatures=None, 
                          apply_iqr_filter=False, iqr_factor=3.0):
    """
    Merge energy and Lindemann data with optional filtering
    
    Args:
        df_energy: Energy DataFrame
        df_lindemann: Lindemann DataFrame  
        outlier_signatures: Set of path signatures to exclude (from Step 1 MSD outliers)
        apply_iqr_filter: Whether to apply IQR filtering (both Lindemann and Energy)
        iqr_factor: IQR multiplier for outlier detection (default 3.0)
    
    Returns:
        tuple: (df_merged_filtered, df_merged_original)
            - df_merged_filtered: Final filtered data
            - df_merged_original: Original merged data (before any filtering)
    """
    print(f"\n{'='*80}")
    print("Merging energy and Lindemann data")
    print("="*80)
    print(f"    Energy records: {len(df_energy)}")
    print(f"    Lindemann records: {len(df_lindemann)}")
    
    # First, create unfiltered merged data for comparison
    df_e_orig = df_energy[['match_key', 'structure', 'system_type', 'system_id', 
                            'temp', 'avg_energy', 'std']].copy()
    df_e_orig.rename(columns={'std': 'energy_std'}, inplace=True)
    df_l_orig = df_lindemann[['match_key', 'delta']].copy()
    
    df_merged_original = pd.merge(df_e_orig, df_l_orig, on='match_key', how='inner')
    df_merged_original['phase'] = df_merged_original['delta'].apply(classify_single_run)
    
    def extract_run(key):
        match = re.search(r'_r(\d+)$', key)
        return int(match.group(1)) if match else None
    
    df_merged_original['run_id'] = df_merged_original['match_key'].apply(extract_run)
    
    print(f"    Original merged: {len(df_merged_original)} records")
    
    # Now apply Step 1 filtering using path signatures (same as Step 6)
    if outlier_signatures and len(outlier_signatures) > 0:
        print(f"\n    [Filter Method 1] Step 1 MSD outliers (path signature matching)")
        print(f"    Outlier signatures to exclude: {len(outlier_signatures)}")
        
        # Filter energy data
        energy_before = len(df_energy)
        df_energy['path_signature'] = df_energy['full_path'].apply(
            lambda x: extract_path_signature(x, is_msd_path=False) if pd.notna(x) else None
        )
        df_energy['is_outlier'] = df_energy['path_signature'].isin(outlier_signatures)
        energy_filtered = df_energy['is_outlier'].sum()
        df_energy = df_energy[~df_energy['is_outlier']].copy()
        df_energy.drop(columns=['path_signature', 'is_outlier'], inplace=True)
        
        # Filter Lindemann data
        lindemann_before = len(df_lindemann)
        if 'directory' in df_lindemann.columns:
            df_lindemann['path_signature'] = df_lindemann['directory'].apply(
                lambda x: extract_path_signature(x, is_msd_path=False) if pd.notna(x) else None
            )
            df_lindemann['is_outlier'] = df_lindemann['path_signature'].isin(outlier_signatures)
            lindemann_filtered = df_lindemann['is_outlier'].sum()
            df_lindemann = df_lindemann[~df_lindemann['is_outlier']].copy()
            df_lindemann.drop(columns=['path_signature', 'is_outlier'], inplace=True)
        else:
            lindemann_filtered = 0
        
        print(f"    Energy filtered: {energy_filtered} records ({energy_filtered/energy_before*100:.1f}%)")
        print(f"    Lindemann filtered: {lindemann_filtered} records ({lindemann_filtered/lindemann_before*100:.1f}%)")
        print(f"    Remaining - Energy: {len(df_energy)}, Lindemann: {len(df_lindemann)}")
    
    # Select needed columns (keep avg_energy name for consistency)
    df_e = df_energy[['match_key', 'structure', 'system_type', 'system_id', 
                      'temp', 'avg_energy', 'std']].copy()
    df_e.rename(columns={'std': 'energy_std'}, inplace=True)
    
    df_l = df_lindemann[['match_key', 'delta']].copy()
    
    # Inner join
    df_merged = pd.merge(df_e, df_l, on='match_key', how='inner')
    
    # Classify
    df_merged['phase'] = df_merged['delta'].apply(classify_single_run)
    
    # Extract run number
    def extract_run(key):
        match = re.search(r'_r(\d+)$', key)
        return int(match.group(1)) if match else None
    
    df_merged['run_id'] = df_merged['match_key'].apply(extract_run)
    
    print(f"    Success: {len(df_merged)} records matched ({len(df_merged)/len(df_energy)*100:.1f}%)")
    
    # Apply IQR filtering if requested (both Lindemann and Energy)
    if apply_iqr_filter:
        print(f"\n    [Filter Method 2] IQR outlier detection (Lindemann + Energy, factor={iqr_factor})")
        merged_before = len(df_merged)
        
        # Temporarily rename for IQR functions
        df_merged['energy_cluster'] = df_merged['avg_energy']
        
        # Detect Lindemann outliers
        lindemann_outlier_keys = detect_lindemann_outliers_iqr(df_merged)
        
        # Detect Energy outliers
        energy_outlier_keys = detect_energy_outliers_iqr(df_merged, iqr_factor)
        
        # Drop temporary column
        df_merged.drop(columns=['energy_cluster'], inplace=True)
        
        # Combine both sets of outliers
        all_iqr_outliers = lindemann_outlier_keys | energy_outlier_keys
        
        print(f"\n    Combined IQR outliers:")
        print(f"      Lindemann outliers: {len(lindemann_outlier_keys)} unique match_keys")
        print(f"      Energy outliers: {len(energy_outlier_keys)} unique match_keys")
        print(f"      Total unique outliers: {len(all_iqr_outliers)} match_keys")
        print(f"      Overlap: {len(lindemann_outlier_keys & energy_outlier_keys)} match_keys")
        
        # Remove outliers
        df_merged = df_merged[~df_merged['match_key'].isin(all_iqr_outliers)].copy()
        iqr_filtered = merged_before - len(df_merged)
        
        print(f"\n    IQR filtered: {iqr_filtered} records ({iqr_filtered/merged_before*100:.1f}%)")
        print(f"    Remaining: {len(df_merged)} records")
    
    # System distribution
    print(f"\n    Final system-wise distribution:")
    for sys_type in sorted(df_merged['system_type'].unique()):
        count = (df_merged['system_type'] == sys_type).sum()
        print(f"      {sys_type}: {count} records")
    
    # Overall phase distribution
    phase_counts = df_merged['phase'].value_counts()
    print(f"\n    Overall phase distribution:")
    for phase, count in sorted(phase_counts.items()):
        print(f"      {phase}: {count} points ({count/len(df_merged)*100:.1f}%)")
    
    # Return both filtered and original for reporting
    return df_merged, df_merged_original


def generate_filtering_report(df_original, df_filtered, method1_applied, method2_applied, iqr_factor=3.0):
    """
    Generate detailed filtering statistics report
    
    Args:
        df_original: Original merged data before IQR filtering
        df_filtered: Final filtered data
        method1_applied: Whether Method 1 (MSD outliers) was applied
        method2_applied: Whether Method 2 (IQR) was applied
        iqr_factor: IQR factor used
    """
    print(f"\n{'='*80}")
    print("Generating Filtering Statistics Report")
    print("="*80)
    
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    report_file = RESULTS_DIR / 'filtering_statistics_report.md'
    csv_file = RESULTS_DIR / 'filtering_statistics_by_structure_temp.csv'
    
    # Overall statistics
    n_original = len(df_original)
    n_filtered = len(df_filtered)
    n_removed = n_original - n_filtered
    pct_removed = (n_removed / n_original * 100) if n_original > 0 else 0
    
    # Filtering methods used
    methods_used = []
    if method1_applied:
        methods_used.append("Method 1 (Step 1 MSD outliers)")
    if method2_applied:
        methods_used.append(f"Method 2 (IQR, factor={iqr_factor})")
    if not methods_used:
        methods_used.append("None (No filtering)")
    
    # Statistics by system_type
    system_stats = []
    for sys_type in sorted(df_original['system_type'].unique()):
        orig_count = (df_original['system_type'] == sys_type).sum()
        filt_count = (df_filtered['system_type'] == sys_type).sum()
        removed = orig_count - filt_count
        pct = (removed / orig_count * 100) if orig_count > 0 else 0
        
        system_stats.append({
            'system_type': sys_type,
            'original': orig_count,
            'filtered': filt_count,
            'removed': removed,
            'percent_removed': pct
        })
    
    # Statistics by (structure, temperature)
    detailed_stats = []
    
    # Get all unique (structure, temp) combinations from original data
    for (structure, temp) in df_original.groupby(['structure', 'temp']).size().index:
        orig_group = df_original[(df_original['structure'] == structure) & (df_original['temp'] == temp)]
        filt_group = df_filtered[(df_filtered['structure'] == structure) & (df_filtered['temp'] == temp)]
        
        orig_count = len(orig_group)
        filt_count = len(filt_group)
        removed = orig_count - filt_count
        pct = (removed / orig_count * 100) if orig_count > 0 else 0
        
        # Get system_type from original data
        sys_type = orig_group['system_type'].iloc[0] if len(orig_group) > 0 else 'Unknown'
        
        detailed_stats.append({
            'system_type': sys_type,
            'structure': structure,
            'temperature': temp,
            'original_points': orig_count,
            'filtered_points': filt_count,
            'removed_points': removed,
            'percent_removed': pct
        })
    
    df_detailed = pd.DataFrame(detailed_stats)
    df_detailed = df_detailed.sort_values(['system_type', 'structure', 'temperature'])
    
    # Add a note row at the top explaining the data
    # Create a comment file alongside the CSV
    csv_readme = csv_file.parent / (csv_file.stem + '_README.txt')
    with open(csv_readme, 'w', encoding='utf-8') as f:
        f.write("=" * 80 + "\n")
        f.write("Filtering Statistics CSV - Data Explanation\n")
        f.write("=" * 80 + "\n\n")
        f.write("DATA SOURCES:\n\n")
        f.write("  Energy:      ~3262 unique records (LAMMPS MD simulations)\n")
        f.write("  Lindemann:   4012 raw records from 3 CSV files\n")
        f.write("               -> 3262 unique after deduplication (by directory+structure+temp)\n")
        f.write("  Merged:      ~8,030 records (inner join on match_key)\n")
        f.write("               Note: Average ~2.5 runs per (structure, temperature) combination\n\n")
        f.write("COLUMN DEFINITIONS:\n\n")
        f.write("  original_points:  Number of data points in the ORIGINAL MERGED dataset\n")
        f.write("                    (Energy + Lindemann data, before any filtering)\n\n")
        f.write("  filtered_points:  Number of data points AFTER applying filtering methods\n\n")
        f.write("  removed_points:   Number of data points removed by filtering\n\n")
        f.write("  percent_removed:  Percentage of data removed (removed/original * 100)\n\n")
        f.write("IMPORTANT NOTES:\n\n")
        f.write("1. 'Original' refers to the merged Energy-Lindemann dataset, NOT raw LAMMPS data\n\n")
        f.write("2. Step 1 MSD filtering is based on Pt and Sn ELEMENTAL MSD:\n")
        f.write("   - msd_Pt.xvg:  Pt atoms' mean square displacement\n")
        f.write("   - msd_Sn.xvg:  Sn atoms' mean square displacement\n")
        f.write("   - NOT msd_Pt-Sn.xvg (Pt-Sn relative distance)\n\n")
        f.write("3. When EITHER Pt or Sn MSD is flagged as anomalous for a simulation run,\n")
        f.write("   that ENTIRE run is excluded from both Energy and Lindemann analyses.\n\n")
        f.write("4. This ensures consistent data quality across all analysis types.\n\n")
        f.write("=" * 80 + "\n")
    
    # Save CSV
    df_detailed.to_csv(csv_file, index=False, encoding='utf-8-sig')
    print(f"    [OK] Detailed CSV saved: {csv_file.name}")
    print(f"    [OK] CSV explanation saved: {csv_readme.name}")
    
    # Write markdown report
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write("# Data Filtering Statistics Report\n\n")
        f.write(f"**Generated**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write("="*80 + "\n\n")
        
        # Important note about data source
        f.write("## Data Source Explanation\n\n")
        f.write("**Important**: The 'Original' data refers to the **merged energy-Lindemann dataset** ")
        f.write("before applying the filtering methods listed below. This merged dataset combines:\n\n")
        f.write("1. **Energy data**: Average cluster energy from LAMMPS MD simulations\n")
        f.write("   - Raw source: ~3262 records (energy_master_*.csv)\n")
        f.write("2. **Lindemann index data**: Calculated from MSD (Mean Square Displacement) analysis\n")
        f.write("   - Raw source: 4012 records from 3 lindemann_master_run_*.csv files\n")
        f.write("   - After deduplication: ~3262 unique (directory, structure, temperature) combinations\n")
        f.write("   - Note: Multiple files contain overlapping data; duplicates are removed before merging\n\n")
        f.write("**Note about MSD types**: Step 1 filtering is based on **Pt and Sn elemental MSD** ")
        f.write("(`msd_Pt.xvg` and `msd_Sn.xvg`), NOT the Pt-Sn distance MSD. When either Pt or Sn ")
        f.write("MSD shows anomalous behavior for a simulation run, that run is excluded from ")
        f.write("**both** energy and Lindemann analyses.\n\n")
        f.write("="*80 + "\n\n")
        
        # Filtering methods
        f.write("## 1. Filtering Methods Applied\n\n")
        for i, method in enumerate(methods_used, 1):
            f.write(f"{i}. {method}\n")
        f.write("\n")
        
        # Overall summary
        f.write("## 2. Overall Statistics\n\n")
        f.write(f"- **Original merged data points**: {n_original:,}\n")
        f.write(f"- **Filtered data points**: {n_filtered:,}\n")
        f.write(f"- **Removed data points**: {n_removed:,} ({pct_removed:.2f}%)\n")
        f.write(f"- **Retention rate**: {100-pct_removed:.2f}%\n\n")
        
        # By system type
        f.write("## 3. Statistics by System Type\n\n")
        f.write("| System Type | Original | Filtered | Removed | % Removed |\n")
        f.write("|-------------|----------|----------|---------|----------|\n")
        for stat in system_stats:
            f.write(f"| {stat['system_type']} | {stat['original']:,} | {stat['filtered']:,} | "
                   f"{stat['removed']:,} | {stat['percent_removed']:.2f}% |\n")
        f.write("\n")
        
        # By structure (summary)
        f.write("## 4. Statistics by Structure (Summary)\n\n")
        structure_summary = df_detailed.groupby('structure').agg({
            'original_points': 'sum',
            'filtered_points': 'sum',
            'removed_points': 'sum'
        }).reset_index()
        structure_summary['percent_removed'] = (
            structure_summary['removed_points'] / structure_summary['original_points'] * 100
        )
        structure_summary = structure_summary.sort_values('removed_points', ascending=False)
        
        f.write("| Structure | Original | Filtered | Removed | % Removed |\n")
        f.write("|-----------|----------|----------|---------|----------|\n")
        for _, row in structure_summary.head(20).iterrows():
            f.write(f"| {row['structure']} | {int(row['original_points']):,} | "
                   f"{int(row['filtered_points']):,} | {int(row['removed_points']):,} | "
                   f"{row['percent_removed']:.2f}% |\n")
        
        if len(structure_summary) > 20:
            f.write(f"\n*Showing top 20 structures by removed points. See CSV for complete data.*\n")
        f.write("\n")
        
        # Temperature-wise analysis
        f.write("## 5. Statistics by Temperature (All Structures)\n\n")
        temp_summary = df_detailed.groupby('temperature').agg({
            'original_points': 'sum',
            'filtered_points': 'sum',
            'removed_points': 'sum'
        }).reset_index()
        temp_summary['percent_removed'] = (
            temp_summary['removed_points'] / temp_summary['original_points'] * 100
        )
        temp_summary = temp_summary.sort_values('temperature')
        
        f.write("| Temperature (K) | Original | Filtered | Removed | % Removed |\n")
        f.write("|-----------------|----------|----------|---------|----------|\n")
        for _, row in temp_summary.iterrows():
            f.write(f"| {int(row['temperature'])} | {int(row['original_points']):,} | "
                   f"{int(row['filtered_points']):,} | {int(row['removed_points']):,} | "
                   f"{row['percent_removed']:.2f}% |\n")
        f.write("\n")
        
        # Detailed breakdown note
        f.write("## 6. Detailed Breakdown\n\n")
        f.write("For complete structure-by-structure, temperature-by-temperature breakdown, ")
        f.write(f"see: `{csv_file.name}`\n\n")
        
        # High outlier structures
        high_outlier_structures = structure_summary[structure_summary['percent_removed'] > 50]
        if len(high_outlier_structures) > 0:
            f.write("## 7. High Outlier Structures (>50% removed)\n\n")
            f.write("| Structure | Original | Removed | % Removed |\n")
            f.write("|-----------|----------|---------|----------|\n")
            for _, row in high_outlier_structures.iterrows():
                f.write(f"| {row['structure']} | {int(row['original_points']):,} | "
                       f"{int(row['removed_points']):,} | {row['percent_removed']:.2f}% |\n")
            f.write("\n")
    
    print(f"    [OK] Markdown report saved: {report_file.name}")
    print(f"    Total removed: {n_removed}/{n_original} ({pct_removed:.2f}%)")


def calculate_support_cv():
    """Calculate support heat capacity"""
    print(f"\n>>> Calculating support heat capacity")
    
    if not SUPPORT_ENERGY_FILE.exists():
        print(f"    Warning: Support file not found, using default 38.2151 meV/K")
        return 38.2151
    
    df_sup = pd.read_csv(SUPPORT_ENERGY_FILE, encoding='utf-8')
    
    # Handle column names
    if '结构' in df_sup.columns:
        df_sup.columns = ['path', 'structure', 'temp', 'run_num', 'total_steps', 'sample_steps', 
                          'avg_energy', 'std', 'min', 'max', 'sample_interval', 
                          'skip_steps', 'full_path']
    
    # Filter support data
    df_sup = df_sup[df_sup['structure'].str.contains('sup-1|sup-2', na=False)]
    
    if len(df_sup) == 0:
        print(f"    Warning: No support data, using default 38.2151 meV/K")
        return 38.2151
    
    # Average by temperature
    df_sup_avg = df_sup.groupby('temp')['avg_energy'].mean().reset_index()
    df_sup_avg = df_sup_avg.sort_values('temp')
    
    # Linear fit
    T = df_sup_avg['temp'].values
    E = df_sup_avg['avg_energy'].values
    
    if len(T) < 2:
        print(f"    Warning: Insufficient support data, using default 38.2151 meV/K")
        return 38.2151
    
    slope, intercept, r_value, p_value, std_err = linregress(T, E)
    
    Cv_support = slope * 1000  # eV/K -> meV/K
    R2 = r_value ** 2
    
    print(f"    Support Cv: {Cv_support:.4f} meV/K")
    print(f"    Fit R2: {R2:.6f}")
    
    return Cv_support


def perform_lindemann_clustering(df_structure, structure_name, n_clusters=3):
    """
    Perform K-means clustering on Lindemann index to auto-detect phase boundaries
    
    Args:
        df_structure: DataFrame for single structure
        structure_name: Name of the structure
        n_clusters: Number of clusters (default 3 for solid/premelting/liquid)
    
    Returns:
        dict: Clustering results with thresholds and labels
    """
    print(f"\n{'='*80}")
    print(f"Lindemann Index Clustering Analysis: {structure_name}")
    print("="*80)
    
    if len(df_structure) < 10:
        print(f"  Warning: Insufficient data ({len(df_structure)} points), need at least 10")
        return None
    
    # Prepare data: (temperature, lindemann_index)
    X = df_structure[['temp', 'delta']].values
    
    # Standardize features
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # K-means clustering
    kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
    df_structure = df_structure.copy()
    df_structure['cluster'] = kmeans.fit_predict(X_scaled)
    
    # Analyze clusters
    print(f"\n>>> Cluster Statistics:")
    cluster_stats = []
    for i in range(n_clusters):
        cluster_data = df_structure[df_structure['cluster'] == i]
        delta_mean = cluster_data['delta'].mean()
        delta_std = cluster_data['delta'].std()
        temp_mean = cluster_data['temp'].mean()
        n_points = len(cluster_data)
        
        cluster_stats.append({
            'cluster_id': i,
            'n_points': n_points,
            'delta_mean': delta_mean,
            'delta_std': delta_std,
            'temp_mean': temp_mean
        })
        
        print(f"  Cluster {i}: n={n_points}, δ={delta_mean:.4f}±{delta_std:.4f}, T_avg={temp_mean:.1f}K")
    
    # Sort clusters by mean delta (solid < premelting < liquid)
    cluster_stats = sorted(cluster_stats, key=lambda x: x['delta_mean'])
    
    # Assign phase labels
    phase_labels = ['solid', 'premelting', 'liquid']
    cluster_to_phase = {stat['cluster_id']: phase_labels[i] for i, stat in enumerate(cluster_stats)}
    
    df_structure['phase_clustered'] = df_structure['cluster'].map(cluster_to_phase)
    
    # Calculate thresholds (boundaries between clusters)
    thresholds = []
    if n_clusters >= 2:
        for i in range(n_clusters - 1):
            lower_cluster = cluster_stats[i]
            upper_cluster = cluster_stats[i + 1]
            # Threshold = midpoint between cluster means
            threshold = (lower_cluster['delta_mean'] + upper_cluster['delta_mean']) / 2
            thresholds.append(threshold)
            print(f"\n  Threshold {phase_labels[i]}/{phase_labels[i+1]}: δ = {threshold:.4f}")
    
    # Visualization
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    cluster_dir = RESULTS_DIR / 'clustering_analysis'
    cluster_dir.mkdir(parents=True, exist_ok=True)
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    fig.suptitle(f'Lindemann Index Clustering Analysis - {structure_name}', 
                 fontsize=14, fontweight='bold')
    
    colors = {'solid': '#3498db', 'premelting': '#e67e22', 'liquid': '#e74c3c'}
    
    # Plot 1: Temperature vs Lindemann with clusters
    ax1 = axes[0]
    for phase in ['solid', 'premelting', 'liquid']:
        df_phase = df_structure[df_structure['phase_clustered'] == phase]
        if len(df_phase) > 0:
            ax1.scatter(df_phase['temp'], df_phase['delta'], 
                       c=colors[phase], alpha=0.6, s=80, 
                       label=f'{phase} (n={len(df_phase)})',
                       edgecolors='black', linewidths=0.5)
    
    # Add threshold lines
    for i, thresh in enumerate(thresholds):
        ax1.axhline(y=thresh, color='red', linestyle='--', linewidth=2, 
                   label=f'Threshold {i+1}: δ={thresh:.4f}', alpha=0.7)
    
    ax1.set_xlabel('Temperature (K)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Lindemann Index δ', fontsize=12, fontweight='bold')
    ax1.set_title('(a) Clustered Data with Auto-Detected Thresholds', 
                  fontsize=12, fontweight='bold')
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Comparison with fixed thresholds
    ax2 = axes[1]
    
    # Original phase classification (fixed thresholds)
    for phase in ['solid', 'premelting', 'liquid']:
        df_phase = df_structure[df_structure['phase'] == phase]
        if len(df_phase) > 0:
            ax2.scatter(df_phase['temp'], df_phase['delta'], 
                       c=colors[phase], alpha=0.6, s=80,
                       label=f'{phase} (n={len(df_phase)})',
                       edgecolors='black', linewidths=0.5)
    
    # Fixed thresholds
    ax2.axhline(y=0.1, color='gray', linestyle='--', linewidth=2, 
               label='Fixed: δ=0.1', alpha=0.7)
    ax2.axhline(y=0.15, color='gray', linestyle='--', linewidth=2, 
               label='Fixed: δ=0.15', alpha=0.7)
    
    ax2.set_xlabel('Temperature (K)', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Lindemann Index δ', fontsize=12, fontweight='bold')
    ax2.set_title('(b) Fixed Thresholds (0.1, 0.15)', 
                  fontsize=12, fontweight='bold')
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    output_file = cluster_dir / f'{structure_name}_clustering_analysis.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n  [OK] Clustering plot saved: {output_file.name}")
    plt.close()
    
    # Save results
    
    
    
    results = {
        'structure': structure_name,
        'n_clusters': n_clusters,
        'thresholds': thresholds,
        'cluster_to_phase': cluster_to_phase,
        'cluster_stats': cluster_stats,
        'df_clustered': df_structure
    }
    
    # Save comparison CSV
    comparison_df = df_structure[['temp', 'delta', 'phase', 'phase_clustered', 'cluster']].copy()
    comparison_csv = cluster_dir / f'{structure_name}_clustering_comparison.csv'
    comparison_df.to_csv(comparison_csv, index=False, encoding='utf-8-sig')
    print(f"  [OK] Comparison CSV saved: {comparison_csv.name}")
    
    return results


def fit_regional_heat_capacity(df_data, Cv_support, structure_name=None):
    """
    Fit three-region heat capacity using individual runs
    
    Args:
        df_data: dataframe to analyze (can be filtered by structure)
        Cv_support: support heat capacity
        structure_name: specific structure name for display (None for all)
    """
    if structure_name:
        print(f"\n{'='*80}")
        print(f"Three-region heat capacity calculation: {structure_name}")
        print("="*80)
    else:
        print(f"\n{'='*80}")
        print("Three-region heat capacity calculation (all data)")
        print("="*80)
    
    regions = ['solid', 'premelting', 'liquid']
    results = {}
    
    for region_name in regions:
        print(f"\n>>> {region_name.capitalize()} region")
        
        df_region = df_data[df_data['phase'] == region_name].copy()
        
        if len(df_region) < 5:
            print(f"    Warning: Insufficient data (n={len(df_region)} < 5), skipping")
            continue
        
        # Extract temperature and energy
        T = df_region['temp'].values
        E = df_region['avg_energy'].values
        
        # Check if all temperatures are the same
        if len(np.unique(T)) < 2:
            print(f"    Warning: Only {len(np.unique(T))} unique temperature(s), cannot fit, skipping")
            continue
        
        # Linear fit
        slope, intercept, r_value, p_value, std_err = linregress(T, E)
        
        # Heat capacity calculation
        Cv_total = slope * 1000  # eV/K -> meV/K
        Cv_cluster = Cv_total - Cv_support
        R2 = r_value ** 2
        
        # Temperature range
        T_min, T_max = T.min(), T.max()
        
        # Save results
        results[region_name] = {
            'n_points': len(df_region),
            'T_range': (T_min, T_max),
            'slope': slope,
            'slope_err': std_err,
            'intercept': intercept,
            'Cv_total': Cv_total,
            'Cv_cluster': Cv_cluster,
            'R2': R2,
            'p_value': p_value,
            'data': df_region[['temp', 'avg_energy', 'delta']].copy()
        }
        
        print(f"    Data points: {len(df_region)}")
        print(f"    Temp range: {T_min:.0f}-{T_max:.0f} K")
        print(f"    Linear fit: E = {slope:.6f} * T + {intercept:.3f}")
        print(f"    Cv_total = {Cv_total:.4f} +/- {std_err*1000:.4f} meV/K")
        print(f"    Cv_cluster = {Cv_cluster:.4f} meV/K")
        
        # Grade marker
        if R2 > 0.95:
            mark = "Excellent"
        elif R2 > 0.90:
            mark = "Good"
        else:
            mark = "Fair"
        
        print(f"    R2 = {R2:.6f} [{mark}]")
        print(f"    p-value = {p_value:.2e}")
    
    return results


def plot_individual_structure_analysis(df_structure, results, Cv_support, structure_name, output_dir):
    """
    Generate individual structure analysis plot (7.3 style)
    
    Args:
        df_structure: DataFrame for this structure
        results: Fitting results dictionary
        Cv_support: Support heat capacity
        structure_name: Name of the structure (e.g., 'Pt8Sn0', 'Cv-1')
        output_dir: Output directory for plots
    """
    
    # Create figure with 2x2 subplots (exactly like 7.3)
    fig, axes = plt.subplots(2, 2, figsize=(18, 14))
    fig.suptitle(f'Step 7.4: Individual Run Analysis - {structure_name}', 
                 fontsize=16, fontweight='bold', y=0.995)
    
    # Color scheme (exactly like 7.3)
    colors = {
        'solid': '#3498db',      # Blue
        'premelting': '#e67e22', # Orange
        'liquid': '#e74c3c'      # Red
    }
    
    phase_labels = {
        'solid': '固态 Solid',
        'premelting': '预熔化 Premelting',
        'liquid': '液态 Liquid'
    }
    
    # Calculate relative energy (subtract minimum)
    E_min = df_structure['avg_energy'].min()
    df_structure = df_structure.copy()
    df_structure['relative_energy'] = df_structure['avg_energy'] - E_min
    
    # ===== Plot 1 (Top Left): Scatter plot with fit lines (a) =====
    ax1 = axes[0, 0]
    
    for phase in ['solid', 'premelting', 'liquid']:
        df_phase = df_structure[df_structure['phase'] == phase]
        
        if len(df_phase) > 0:
            # Scatter plot with relative energy
            ax1.scatter(df_phase['temp'], df_phase['relative_energy'], 
                       c=colors[phase], alpha=0.5, s=50, edgecolors='black', linewidths=0.5,
                       label=f'{phase_labels[phase]} (n={len(df_phase)})',
                       zorder=3)
            
            # Fit line (convert to relative)
            if phase in results:
                res = results[phase]
                T_min, T_max = res['T_range']
                T_fit = np.linspace(T_min, T_max, 100)
                E_fit = res['slope'] * T_fit + res['intercept'] - E_min
                
                ax1.plot(T_fit, E_fit, color=colors[phase], linewidth=3, 
                        linestyle='--', alpha=0.8,
                        label=f'{phase_labels[phase]} 拟合 (R²={res["R2"]:.4f})',
                        zorder=2)
    
    ax1.set_xlabel('温度 Temperature (K)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('相对能量 Relative Energy (eV)', fontsize=12, fontweight='bold')
    ax1.set_title('(a) 单次模拟能量分布与线性拟合\nIndividual Run Energy vs Temperature', 
                  fontsize=13, fontweight='bold', pad=10)
    ax1.legend(fontsize=9, loc='upper left', framealpha=0.95, ncol=2)
    ax1.grid(True, alpha=0.3, linestyle=':', linewidth=0.8)
    
    # ===== Plot 2 (Top Right): Heat capacity bar chart (b) =====
    ax2 = axes[0, 1]
    
    regions = []
    cv_values = []
    cv_errors = []
    r2_values = []
    point_counts = []
    
    for phase in ['solid', 'premelting', 'liquid']:
        if phase in results:
            regions.append(phase_labels[phase])
            cv_values.append(results[phase]['Cv_cluster'])
            cv_errors.append(results[phase]['slope_err'] * 1000)
            r2_values.append(results[phase]['R2'])
            point_counts.append(results[phase]['n_points'])
    
    x = np.arange(len(regions))
    width = 0.6
    
    bars = ax2.bar(x, cv_values, width, 
                   color=[colors[p] for p in ['solid', 'premelting', 'liquid'] if p in results], 
                   alpha=0.8, edgecolor='black', linewidth=1.5,
                   yerr=cv_errors, capsize=8, error_kw={'linewidth': 2})
    
    # Annotate bars
    for i, (bar, cv, cv_err, r2, n) in enumerate(zip(bars, cv_values, cv_errors, r2_values, point_counts)):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height + cv_err + 0.2,
                f'{cv:.3f}±{cv_err:.3f}\nmeV/K\nR²={r2:.4f}\n(n={n})',
                ha='center', va='bottom', fontsize=10, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
    
    ax2.set_xlabel('相态区域 Phase Region', fontsize=12, fontweight='bold')
    ax2.set_ylabel('团簇热容 Cluster Heat Capacity Cv (meV/K)', fontsize=12, fontweight='bold')
    ax2.set_title('(b) 三段分区热容结果\nThree-Region Heat Capacity Results', 
                  fontsize=13, fontweight='bold', pad=10)
    ax2.set_xticks(x)
    ax2.set_xticklabels(regions, fontsize=10)
    ax2.grid(True, alpha=0.3, linestyle=':', axis='y', linewidth=0.8)
    if len(cv_values) > 0:
        ax2.set_ylim(0, max(cv_values) * 1.6)
    
    # ===== Plot 3 (Bottom Left): Temperature-phase distribution (c) =====
    ax3 = axes[1, 0]
    
    temp_phase = df_structure.groupby(['temp', 'phase']).size().unstack(fill_value=0)
    temp_sorted = sorted(df_structure['temp'].unique())
    
    # Stacked bar chart
    solid_counts = [temp_phase.loc[t, 'solid'] if t in temp_phase.index and 'solid' in temp_phase.columns else 0 for t in temp_sorted]
    pre_counts = [temp_phase.loc[t, 'premelting'] if t in temp_phase.index and 'premelting' in temp_phase.columns else 0 for t in temp_sorted]
    liquid_counts = [temp_phase.loc[t, 'liquid'] if t in temp_phase.index and 'liquid' in temp_phase.columns else 0 for t in temp_sorted]
    
    x_pos = np.arange(len(temp_sorted))
    
    ax3.bar(x_pos, solid_counts, label=phase_labels['solid'], 
            color=colors['solid'], alpha=0.8, edgecolor='black', linewidth=0.5)
    ax3.bar(x_pos, pre_counts, bottom=solid_counts, label=phase_labels['premelting'], 
            color=colors['premelting'], alpha=0.8, edgecolor='black', linewidth=0.5)
    ax3.bar(x_pos, liquid_counts, bottom=np.array(solid_counts)+np.array(pre_counts), 
            label=phase_labels['liquid'], color=colors['liquid'], alpha=0.8, edgecolor='black', linewidth=0.5)
    
    ax3.set_xlabel('温度 Temperature (K)', fontsize=12, fontweight='bold')
    ax3.set_ylabel('模拟次数 Number of Runs', fontsize=12, fontweight='bold')
    ax3.set_title('(c) 各温度点的相态分布\nPhase Distribution at Each Temperature', 
                  fontsize=13, fontweight='bold', pad=10)
    
    # Adjust x-axis labels based on number of temperatures
    if len(temp_sorted) > 15:
        step = max(1, len(temp_sorted) // 10)
        ax3.set_xticks(x_pos[::step])
        ax3.set_xticklabels([f'{int(t)}' for t in temp_sorted[::step]], rotation=45)
    else:
        ax3.set_xticks(x_pos)
        ax3.set_xticklabels([f'{int(t)}' for t in temp_sorted], rotation=45)
    
    ax3.legend(fontsize=10, loc='upper left', framealpha=0.95)
    ax3.grid(True, alpha=0.3, linestyle=':', axis='y', linewidth=0.8)
    
    # ===== Plot 4 (Bottom Right): Lindemann index distribution (d) =====
    ax4 = axes[1, 1]
    
    for phase in ['solid', 'premelting', 'liquid']:
        df_phase = df_structure[df_structure['phase'] == phase]
        if len(df_phase) > 0:
            ax4.scatter(df_phase['temp'], df_phase['delta'], 
                       c=colors[phase], alpha=0.6, s=60, 
                       edgecolors='black', linewidths=0.8,
                       label=f'{phase_labels[phase]} (n={len(df_phase)})',
                       zorder=3)
    
    # Threshold lines
    ax4.axhline(y=0.1, color='gray', linestyle='--', linewidth=2.5, 
                label='固态/预熔化阈值 δ=0.1', alpha=0.7, zorder=1)
    ax4.axhline(y=0.15, color='red', linestyle='--', linewidth=2.5, 
                label='预熔化/液态阈值 δ=0.15', alpha=0.7, zorder=1)
    
    # Shade regions
    if len(df_structure) > 0:
        delta_max = max(df_structure['delta']) * 1.1
        ax4.axhspan(0, 0.1, alpha=0.1, color=colors['solid'], zorder=0)
        ax4.axhspan(0.1, 0.15, alpha=0.1, color=colors['premelting'], zorder=0)
        ax4.axhspan(0.15, delta_max, alpha=0.1, color=colors['liquid'], zorder=0)
        ax4.set_ylim(0, delta_max)
    
    ax4.set_xlabel('温度 Temperature (K)', fontsize=12, fontweight='bold')
    ax4.set_ylabel('Lindemann 指数 δ', fontsize=12, fontweight='bold')
    ax4.set_title('(d) 单次模拟 Lindemann 指数分布\nIndividual Run Lindemann Index Distribution', 
                  fontsize=13, fontweight='bold', pad=10)
    ax4.legend(fontsize=9, loc='upper left', framealpha=0.95, ncol=2)
    ax4.grid(True, alpha=0.3, linestyle=':', linewidth=0.8)
    
    plt.tight_layout(rect=[0, 0, 1, 0.99])
    
    # Save figure
    output_file = output_dir / f'{structure_name}_individual_runs_analysis.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"      > Plot saved: {structure_name}_individual_runs_analysis.png")
    
    plt.close()


def generate_system_comparison_report(df_merged, Cv_support):
    """Generate comparison report for all systems - BY INDIVIDUAL STRUCTURE"""
    print(f"\n{'='*80}")
    print("Generating individual structure comparison report")
    print("="*80)
    
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    report_file = RESULTS_DIR / 'step7_4_multi_system_comparison.md'
    
    # Create subdirectory for individual structure plots
    individual_plots_dir = RESULTS_DIR / 'individual_structure_plots'
    individual_plots_dir.mkdir(parents=True, exist_ok=True)
    
    # Analyze each INDIVIDUAL STRUCTURE (system_id), not system_type
    structures = sorted(df_merged['system_id'].unique())
    system_results = {}
    
    print(f"\n  Total structures to analyze: {len(structures)}")
    
    for system_id in structures:
        print(f"\n  Analyzing {system_id}...")
        # Filter by system_id instead of system_type
        df_structure = df_merged[df_merged['system_id'] == system_id].copy()
        
        # Only analyze if enough data points
        if len(df_structure) < 10:
            print(f"    Warning: Only {len(df_structure)} points, skipping")
            continue
        
        results = fit_regional_heat_capacity(df_structure, Cv_support, system_id)
        system_results[system_id] = results
        
        # Generate individual visualization for this structure
        if results:
            plot_individual_structure_analysis(df_structure, results, Cv_support, system_id, individual_plots_dir)
    
    # Write report
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write("# Step 7.4: 多体系热容分析对比报告\n\n")
        f.write("Multi-System Heat Capacity Analysis Comparison Report\n\n")
        f.write("=" * 80 + "\n\n")
        
        f.write(f"**报告生成时间 Report Date**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        # Overview
        f.write("## 1. 分析概览 Analysis Overview\n\n")
        f.write(f"- **总数据点 Total Points**: {len(df_merged)}\n")
        f.write(f"- **结构数量 Number of Structures**: {len(structures)}\n")
        f.write(f"- **已分析结构 Analyzed Structures**: {len(system_results)}\n")
        f.write(f"- **Support 热容 Support Cv**: {Cv_support:.4f} meV/K\n\n")
        
        f.write("### 1.1 结构列表 Structure List\n\n")
        
        # Group by system_type for display
        system_type_groups = df_merged.groupby('system_type')['system_id'].apply(lambda x: sorted(x.unique())).to_dict()
        
        for sys_type, ids in sorted(system_type_groups.items()):
            count = sum(df_merged['system_id'] == sid for sid in ids).sum()
            f.write(f"- **{sys_type}**: {len(ids)} 个结构 ({', '.join(ids[:5])}{'...' if len(ids) > 5 else ''})\n")
        
        # System-wise results
        f.write("\n## 2. 各结构热容结果 Structure-wise Heat Capacity Results\n\n")
        
        for system_id in sorted(system_results.keys()):
            results = system_results[system_id]
            
            if not results:
                f.write(f"### {system_id}\n\n")
                f.write("⚠️ 数据不足，无法完成分析 Insufficient data for analysis\n\n")
                continue
            
            f.write(f"### {system_id}\n\n")
            
            # Summary table
            f.write("| 区域 Region | 温度范围 Temp Range | 数据点 Points | Cv_cluster (meV/K) | R² |\n")
            f.write("|-------------|---------------------|--------------|-------------------|----||\n")
            
            phase_map = {'solid': '固态 Solid', 'premelting': '预熔化 Premelting', 'liquid': '液态 Liquid'}
            
            for phase in ['solid', 'premelting', 'liquid']:
                if phase in results:
                    res = results[phase]
                    T_range = f"{res['T_range'][0]:.0f}-{res['T_range'][1]:.0f} K"
                    cv_str = f"{res['Cv_cluster']:.4f} ± {res['slope_err']*1000:.4f}"
                    f.write(f"| {phase_map[phase]} | {T_range} | {res['n_points']} | {cv_str} | {res['R2']:.6f} |\n")
            
            f.write("\n")
        
        # Comparison tables - group by system_type
        f.write("## 3. 体系对比 System Comparison\n\n")
        
        for sys_type in sorted(system_type_groups.keys()):
            structure_ids = system_type_groups[sys_type]
            
            f.write(f"### {sys_type} 系列\n\n")
            
            # Solid comparison
            f.write("#### 固态热容 Solid Heat Capacity\n\n")
            f.write("| 结构 Structure | Cv_cluster (meV/K) | R² | 数据点 Points |\n")
            f.write("|----------------|-------------------|----|---------------|\n")
            
            for system_id in structure_ids:
                if system_id in system_results and 'solid' in system_results[system_id]:
                    res = system_results[system_id]['solid']
                    cv_str = f"{res['Cv_cluster']:.4f} ± {res['slope_err']*1000:.4f}"
                    f.write(f"| {system_id} | {cv_str} | {res['R2']:.6f} | {res['n_points']} |\n")
            
            f.write("\n")
            
            # Liquid comparison
            f.write("#### 液态热容 Liquid Heat Capacity\n\n")
            f.write("| 结构 Structure | Cv_cluster (meV/K) | R² | 数据点 Points |\n")
            f.write("|----------------|-------------------|----|---------------|\n")
            
            for system_id in structure_ids:
                if system_id in system_results and 'liquid' in system_results[system_id]:
                    res = system_results[system_id]['liquid']
                    cv_str = f"{res['Cv_cluster']:.4f} ± {res['slope_err']*1000:.4f}"
                    f.write(f"| {system_id} | {cv_str} | {res['R2']:.6f} | {res['n_points']} |\n")
            
            f.write("\n")
        
        # Conclusions - find best structures overall
        f.write("\n## 4. 结论 Conclusions\n\n")
        
        # Find best R2 for each region across ALL structures
        for region in ['solid', 'premelting', 'liquid']:
            best_r2 = 0
            best_structure = None
            for system_id in system_results:
                results = system_results[system_id]
                if region in results and results[region]['R2'] > best_r2:
                    best_r2 = results[region]['R2']
                    best_structure = system_id
            
            if best_structure:
                f.write(f"- **{region.capitalize()} 区域最佳拟合 Best fit**: {best_structure} (R² = {best_r2:.6f})\n")
        
        f.write(f"\n---\n")
        f.write(f"**脚本版本 Script Version**: step7_4_multi_system_heat_capacity.py v1.0\n")
    
    print(f"    [OK] Report saved: {report_file}")
    
    # Save merged data
    # IMPORTANT: Replace 'structure' with 'system_id' for consistency
    # This ensures Cv-1/Cv-2/.../Cv-5 all show as 'Cv' in saved CSV
    df_merged_export = df_merged.copy()
    df_merged_export['structure'] = df_merged_export['system_id']
    
    csv_file = RESULTS_DIR / 'step7_4_all_systems_data.csv'
    df_merged_export.to_csv(csv_file, index=False, encoding='utf-8-sig')
    print(f"    [OK] Data saved: {csv_file}")
    
    return system_results


def plot_multi_system_comparison(df_merged, system_results, Cv_support):
    """Generate multi-system comparison visualization - BY INDIVIDUAL STRUCTURE"""
    print(f"\n>>> Generating multi-structure comparison plots")
    
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    
    # Get list of structures (system_id)
    structures = sorted(df_merged['system_id'].unique())
    
    # Group by system_type for coloring
    system_type_map = df_merged[['system_id', 'system_type']].drop_duplicates().set_index('system_id')['system_type'].to_dict()
    
    # Create color map by system_type
    system_types = sorted(df_merged['system_type'].unique())
    cmap = plt.cm.get_cmap('tab10')
    type_colors = {sys_type: cmap(i) for i, sys_type in enumerate(system_types)}
    
    # Map structure to color via its type
    structure_colors = {sid: type_colors[system_type_map[sid]] for sid in structures}
    
    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(18, 14))
    fig.suptitle('Step 7.4: Multi-Structure Heat Capacity Comparison (Individual Analysis)', 
                 fontsize=16, fontweight='bold', y=0.995)
    
    # Plot 1: All structures scatter (color by system_type)
    ax1 = axes[0, 0]
    
    # Plot by system_type for legend
    for sys_type in system_types:
        df_type = df_merged[df_merged['system_type'] == sys_type]
        ax1.scatter(df_type['temp'], df_type['avg_energy'], 
                   c=[type_colors[sys_type]], alpha=0.5, s=30,
                   label=f'{sys_type}', edgecolors='black', linewidths=0.3)
    
    ax1.set_xlabel('温度 Temperature (K)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('能量 Energy (eV)', fontsize=12, fontweight='bold')
    ax1.set_title('(a) 所有结构能量分布 (按体系类型着色)\nAll Structures Energy Distribution (Colored by System Type)', 
                  fontsize=13, fontweight='bold', pad=10)
    ax1.legend(fontsize=9, loc='upper left', ncol=2)
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Solid Cv comparison (top structures with R²>0.9)
    ax2 = axes[0, 1]
    solid_data = []
    
    for system_id in structures:
        if system_id in system_results and 'solid' in system_results[system_id]:
            res = system_results[system_id]['solid']
            if res['R2'] > 0.9:  # Only show high-quality results
                solid_data.append({
                    'id': system_id,
                    'cv': res['Cv_cluster'],
                    'err': res['slope_err'] * 1000,
                    'r2': res['R2'],
                    'color': structure_colors[system_id]
                })
    
    if solid_data:
        # Sort by R² descending
        solid_data.sort(key=lambda x: x['r2'], reverse=True)
        # Take top 15
        solid_data = solid_data[:15]
        
        x = np.arange(len(solid_data))
        bars = ax2.bar(x, [d['cv'] for d in solid_data], 
                      color=[d['color'] for d in solid_data],
                      alpha=0.8, edgecolor='black', linewidth=1.0,
                      yerr=[d['err'] for d in solid_data], capsize=4)
        
        ax2.set_xlabel('结构 Structure', fontsize=12, fontweight='bold')
        ax2.set_ylabel('固态热容 Solid Cv (meV/K)', fontsize=12, fontweight='bold')
        ax2.set_title(f'(b) 固态热容对比 (R²>0.9, Top 15)\nSolid Heat Capacity Comparison', 
                      fontsize=13, fontweight='bold', pad=10)
        ax2.set_xticks(x)
        ax2.set_xticklabels([d['id'] for d in solid_data], rotation=45, ha='right', fontsize=8)
        ax2.grid(True, alpha=0.3, axis='y')
    else:
        ax2.text(0.5, 0.5, 'No high-quality solid data\n(R² > 0.9)', 
                ha='center', va='center', fontsize=14, transform=ax2.transAxes)
    
    # Plot 3: Liquid Cv comparison (top structures with R²>0.9)
    ax3 = axes[1, 0]
    liquid_data = []
    
    for system_id in structures:
        if system_id in system_results and 'liquid' in system_results[system_id]:
            res = system_results[system_id]['liquid']
            if res['R2'] > 0.9:
                liquid_data.append({
                    'id': system_id,
                    'cv': res['Cv_cluster'],
                    'err': res['slope_err'] * 1000,
                    'r2': res['R2'],
                    'color': structure_colors[system_id]
                })
    
    if liquid_data:
        liquid_data.sort(key=lambda x: x['r2'], reverse=True)
        liquid_data = liquid_data[:15]
        
        x = np.arange(len(liquid_data))
        bars = ax3.bar(x, [d['cv'] for d in liquid_data], 
                      color=[d['color'] for d in liquid_data],
                      alpha=0.8, edgecolor='black', linewidth=1.0,
                      yerr=[d['err'] for d in liquid_data], capsize=4)
        
        ax3.set_xlabel('结构 Structure', fontsize=12, fontweight='bold')
        ax3.set_ylabel('液态热容 Liquid Cv (meV/K)', fontsize=12, fontweight='bold')
        ax3.set_title(f'(c) 液态热容对比 (R²>0.9, Top 15)\nLiquid Heat Capacity Comparison', 
                      fontsize=13, fontweight='bold', pad=10)
        ax3.set_xticks(x)
        ax3.set_xticklabels([d['id'] for d in liquid_data], rotation=45, ha='right', fontsize=8)
        ax3.grid(True, alpha=0.3, axis='y')
    else:
        ax3.text(0.5, 0.5, 'No high-quality liquid data\n(R² > 0.9)', 
                ha='center', va='center', fontsize=14, transform=ax3.transAxes)
    
    # Plot 4: R² distribution by system_type
    ax4 = axes[1, 1]
    
    # Collect R² values by system_type
    r2_by_type = {sys_type: {'solid': [], 'liquid': []} for sys_type in system_types}
    
    for system_id in structures:
        sys_type = system_type_map[system_id]
        if system_id in system_results:
            if 'solid' in system_results[system_id]:
                r2_by_type[sys_type]['solid'].append(system_results[system_id]['solid']['R2'])
            if 'liquid' in system_results[system_id]:
                r2_by_type[sys_type]['liquid'].append(system_results[system_id]['liquid']['R2'])
    
    # Box plot
    solid_r2_data = []
    liquid_r2_data = []
    labels = []
    
    for sys_type in system_types:
        if r2_by_type[sys_type]['solid']:
            solid_r2_data.append(r2_by_type[sys_type]['solid'])
            liquid_r2_data.append(r2_by_type[sys_type]['liquid'])
            labels.append(sys_type)
    
    if solid_r2_data:
        positions_solid = np.arange(len(labels)) * 2
        positions_liquid = positions_solid + 0.8
        
        bp1 = ax4.boxplot(solid_r2_data, positions=positions_solid, widths=0.6,
                          patch_artist=True, labels=labels, showfliers=False)
        bp2 = ax4.boxplot(liquid_r2_data, positions=positions_liquid, widths=0.6,
                          patch_artist=True, labels=[''] * len(labels), showfliers=False)
        
        for patch in bp1['boxes']:
            patch.set_facecolor('skyblue')
        for patch in bp2['boxes']:
            patch.set_facecolor('salmon')
        
        ax4.axhline(y=0.9, color='red', linestyle='--', linewidth=2, alpha=0.7, label='R²=0.9 threshold')
        ax4.set_xlabel('体系类型 System Type', fontsize=12, fontweight='bold')
        ax4.set_ylabel('R² 值', fontsize=12, fontweight='bold')
        ax4.set_title('(d) R² 分布 (固态 vs 液态)\nR² Distribution (Solid vs Liquid)', 
                      fontsize=13, fontweight='bold', pad=10)
        ax4.set_xticks(positions_solid + 0.4)
        ax4.set_xticklabels(labels, rotation=45, ha='right')
        ax4.legend(['R²=0.9', '固态 Solid', '液态 Liquid'], fontsize=9)
        ax4.grid(True, alpha=0.3, axis='y')
        ax4.set_ylim(-0.1, 1.1)
    
    plt.tight_layout(rect=[0, 0, 1, 0.99])
    
    # Save
    output_file = RESULTS_DIR / 'step7_4_multi_system_comparison.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"    [OK] Figure saved: {output_file}")
    
    plt.close()


def main():
    """Main function with command-line argument support"""
    
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='Step 7.4: Multi-System Heat Capacity Analysis with Optional Filtering',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Default: No filtering (use all data)
  python step7_4_multi_system_heat_capacity.py

  # Apply Step 1 MSD outlier filtering only
  python step7_4_multi_system_heat_capacity.py --msd-filter

  # Apply IQR filtering only (Lindemann + Energy)
  python step7_4_multi_system_heat_capacity.py --iqr-filter

  # Apply both filtering methods
  python step7_4_multi_system_heat_capacity.py --msd-filter --iqr-filter

  # Adjust IQR threshold (default 3.0)
  python step7_4_multi_system_heat_capacity.py --iqr-filter --iqr-factor 2.5
        """
    )
    
    parser.add_argument(
        '--msd-filter',
        action='store_true',
        help='Enable Step 1 MSD outlier filtering (path signature matching)'
    )
    
    parser.add_argument(
        '--iqr-filter',
        action='store_true',
        help='Enable IQR-based outlier detection (both Lindemann and Energy)'
    )
    
    parser.add_argument(
        '--iqr-factor',
        type=float,
        default=3.0,
        help='IQR multiplier for outlier detection (default: 3.0, stricter than 1.5)'
    )
    
    parser.add_argument(
        '--cluster-analysis',
        type=str,
        default=None,
        help='Structure name for clustering analysis (e.g., Pt6Sn8). Uses K-means to auto-detect phase boundaries.'
    )
    
    args = parser.parse_args()
    
    print(f"\n{'='*80}")
    print("Step 7.4: Multi-Structure Individual Run Heat Capacity Analysis")
    print("="*80)
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    
    # Display filtering configuration
    print(f"[Filtering Configuration]")
    if args.msd_filter:
        print(f"  [ON]  Filter Method 1: Step 1 MSD outliers (path signature matching)")
    else:
        print(f"  [OFF] Filter Method 1: Disabled")
    
    if args.iqr_filter:
        print(f"  [ON]  Filter Method 2: IQR outliers (Lindemann + Energy, factor={args.iqr_factor})")
    else:
        print(f"  [OFF] Filter Method 2: Disabled")
    
    if not args.msd_filter and not args.iqr_filter:
        print(f"  [!!!] WARNING: No filtering applied - using all raw data!")
    
    print()
    
    # User can specify which systems to analyze
    # Set to None to analyze all systems
    # Examples:
    #   system_filter = ['Cv']  # Only Cv series
    #   system_filter = ['Cv', 'Pt8SnX', 'Pt6SnX']  # Multiple series
    #   system_filter = None  # All systems
    # 
    # NOTE: This filters by system_type (Cv, Pt8SnX, etc.)
    #       But analysis is done PER STRUCTURE (Cv-1, Pt8Sn3, etc.)
    
    system_filter = None  # Change this to filter specific system types
    
    print(f"System type filter: {system_filter if system_filter else 'All system types'}")
    print(f"Analysis level: Individual structures (system_id)")
    print(f"Note: Each structure (e.g., Cv-1, Pt8Sn3) analyzed separately\n")
    
    # 1. Load data
    df_energy = load_energy_data(CLUSTER_ENERGY_FILE, system_filter, 'cluster')
    df_lindemann = load_lindemann_individual_runs(system_filter)
    
    if df_energy is None or df_lindemann is None:
        print("\nError: Data loading failed")
        return
    
    if len(df_energy) == 0 or len(df_lindemann) == 0:
        print("\nError: No data after filtering")
        return
    
    # 1.5. Load Step 1 outliers for data quality filtering (if enabled)
    outliers = set()
    if args.msd_filter:
        print("\n" + "="*80)
        print("Loading Step 1 Outlier Filtering Results")
        print("="*80)
        outliers = load_outliers()
    
    # 2. Merge data with optional filtering
    print("\n" + "="*80)
    print("Merging Energy and Lindemann Data")
    print("="*80)
    df_merged, df_merged_original = merge_energy_lindemann(
        df_energy, 
        df_lindemann, 
        outlier_signatures=outliers if args.msd_filter else None,
        apply_iqr_filter=args.iqr_filter,
        iqr_factor=args.iqr_factor
    )
    
    if len(df_merged) == 0:
        print("\nError: Data matching failed")
        return
    
    # 2.5. Generate filtering statistics report
    if args.msd_filter or args.iqr_filter:
        generate_filtering_report(
            df_merged_original,  # Original unfiltered data
            df_merged,           # Final filtered data
            method1_applied=args.msd_filter,
            method2_applied=args.iqr_filter,
            iqr_factor=args.iqr_factor
        )
    
    # 3. Calculate support Cv
    Cv_support = calculate_support_cv()
    
    # 3.5. Optional: Clustering analysis for specific structure
    if args.cluster_analysis:
        print(f"\n{'='*80}")
        print(f"Performing Clustering Analysis for: {args.cluster_analysis}")
        print("="*80)
        
        # Filter data for the specified structure
        df_cluster_target = df_merged[df_merged['structure'] == args.cluster_analysis].copy()
        
        if len(df_cluster_target) == 0:
            print(f"  [ERROR] Structure '{args.cluster_analysis}' not found in data!")
            print(f"  Available structures: {sorted(df_merged['structure'].unique())[:10]}...")
        else:
            print(f"  Found {len(df_cluster_target)} data points for {args.cluster_analysis}")
            clustering_results = perform_lindemann_clustering(df_cluster_target, args.cluster_analysis, n_clusters=3)
            
            if clustering_results:
                # Optionally: Re-analyze heat capacity using clustered phases
                print(f"\n  Analyzing heat capacity with clustered phase boundaries...")
                df_cluster_target['phase_original'] = df_cluster_target['phase']
                df_cluster_target['phase'] = clustering_results['df_clustered']['phase_clustered']
                
                cluster_results = fit_regional_heat_capacity(df_cluster_target, Cv_support, args.cluster_analysis)
                
                if cluster_results:
                    print(f"\n  >>> Comparison: Fixed vs Clustered Thresholds")
                    print(f"  Fixed thresholds: δ=0.10, 0.15")
                    print(f"  Clustered thresholds: {', '.join([f'δ={t:.4f}' for t in clustering_results['thresholds']])}")
                    
                    # Compare Cv values
                    for phase in ['solid', 'premelting', 'liquid']:
                        if phase in cluster_results:
                            print(f"  {phase}: Cv_cluster = {cluster_results[phase]['Cv_cluster']:.4f} meV/K "
                                  f"(R2={cluster_results[phase]['R2']:.4f}, n={cluster_results[phase]['n_points']})")
    
    # 4. Multi-system analysis
    system_results = generate_system_comparison_report(df_merged, Cv_support)
    
    # 5. Visualization
    plot_multi_system_comparison(df_merged, system_results, Cv_support)
    
    # 6. Summary
    print(f"\n{'='*80}")
    print("Step 7.4 Multi-System Analysis Complete")
    print("="*80)
    
    print(f"\n[Analysis Summary]")
    print(f"  Total data points: {len(df_merged)}")
    print(f"  Number of system types: {len(df_merged['system_type'].unique())}")
    print(f"  Number of structures: {len(df_merged['system_id'].unique())}")
    print(f"  Structures analyzed: {len(system_results)}")
    
    print(f"\n[Structure Distribution]")
    for sys_type in sorted(df_merged['system_type'].unique()):
        structures = df_merged[df_merged['system_type'] == sys_type]['system_id'].unique()
        print(f"  {sys_type}: {len(structures)} structures")
    
    print(f"\nEnd time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*80)


if __name__ == '__main__':
    main()
