#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
Step6: 能量分析 (改进版) - 多体系支持 + 4级路径签名算法
================================================================================

作者: GitHub Copilot
日期: 2025-10-18
版本: v2.1

❓❓❓  关键问题：载体热容未知  ❓❓❓
================================================================================
本脚本计算的热容基于LAMMPS的etotal输出，该值包含：
  - 团簇部分：~60个原子 (Pt + Sn + O，实际数量因体系而异)
  - 载体部分：240个原子 (Al₂O₃)
  - 总计：~300个原子

热容计算公式：
  Cv_total = dE_total/dT
           = Cv_cluster + Cv_support
           = Cv_cluster + C (常数，❓未知)

❓ 我们不知道载体热容 Cv_support 的确切值！
  - 载体是固定的 240 个 Al₂O₃ 原子
  - Cv_support 应该是常数（不随团簇大小变化）
  - 但需要单独模拟才能确定其数值

输出数据的含义：
  - Cv_total_meV_K：团簇 + 载体的总热容（包含未知常数）
  - Cv_per_atom_meV_K：Cv_total / N_cluster（表观值，含载体）

获得准确团簇热容的方法：
  1. ✅ [推荐] 单独模拟纯 Al₂O₃ 体系（240 atoms）
     测量 Cv_support，然后扣除
  2. ⚠️  [替代] 使用拟合估计值（~18-21 meV/K）
     python subtract_support_v2.py
     注意：这只是估计值，需要实验验证
  
详细说明请参阅：CRITICAL_SUPPORT_HEAT_CAPACITY_UNKNOWN.md
================================================================================

功能概述
========
本脚本用于分析LAMMPS模拟的能量数据,支持多体系分析、MSD异常筛选、
热容计算和熔化温度识别。

核心功能
========
1. **智能数据加载与筛选**
   - 加载LAMMPS能量数据 (energy_master_*.csv)
   - 支持正则表达式体系筛选 (O1, O2, Pt8SnX等)
   - 可选MSD异常筛选 (与step1交叉验证)
   - 4级路径签名算法,精确匹配批次信息 (run3, run2等)

2. **能量计算与分析**
   - 自动识别化学式中的Pt、Sn、O原子数
   - 计算每原子能量 (eV/atom) - 注意：能量总和包含载体！
   - 对每个体系计算相对能量 (E - E_min)
   - 按体系系列分类 (O1, O2, Pt8SnX, Pt6SnX等)

3. **热容计算** ❓ 重要：包含未知的载体贡献
   - 计算总热容: Cv_total = dE_total/dT
   - 计算每原子热容: Cv_per_atom = d(E/N)/dT
   - ❓ E_total = E_cluster + E_support，其中 E_support 的热容未知
   - Cv_support 是常数（240个Al₂O₃），但数值待测量
   - 热容跃变检测，识别熔化温度 Tm

4. **数据可视化**
   - 能量-温度曲线 (按Sn含量分组)
   - 热容对比图 (按系列分组)
   - 热容热力图 (结构 vs 温度)
   - 氧系列综合对比
   - 能量波动分析

5. **4级路径签名算法** (核心改进)
   - 解决批次重复问题 (run3, run2, run4, run5)
   - MSD路径: batch/parent/composition/run
     例: run3/o2/o2pt4sn6/t1000.r24.gpu0
   - 能量路径: batch/parent/composition/run
     例: run3/o2/o2pt7sn7/t200.r0.gpu0
   - 自动检测批次标识符,无批次时使用3级签名
   - 精确1:1匹配,避免误判

命令行参数
==========
--no-filter : 不使用MSD异常筛选,分析所有能量数据
              (适用于某些体系数据全被筛选时,如Pt3Sn5 @ 1100K)

使用示例
========
# 默认模式: 使用MSD异常筛选
python step6_energy_analysis_v2.py

# 不筛选模式: 使用所有数据
python step6_energy_analysis_v2.py --no-filter

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
1. energy_master_*.csv - 能量数据主文件
   列: 路径, 结构, 温度, 模拟序号, 总步数, 采样步数, 
       平均能量, 标准差, 最小值, 最大值, 采样间隔, 跳过步数, 完整路径

2. large_D_outliers.csv - MSD异常记录 (可选)
   列: composition, filepath, element, temperature, D_value, ...

输出文件
========
results/energy_analysis_v2/
├── energy_per_system.csv              # 每体系能量数据 (包含载体)
├── heat_capacity_per_system.csv       # 每体系热容数据 ⚠️ 含载体贡献！
├── melting_temperatures.csv           # 熔化温度汇总
├── energy_diffusion_merged.csv        # 能量-扩散系数关联
├── energy_analysis_v2_report.txt      # 分析报告
├── Pt8SnX_Energy_vs_T_Sn*.png        # 能量-温度曲线
├── HeatCapacity_comparison_*.png      # 热容对比图 ⚠️ 含载体！
├── HeatCapacity_heatmap_*.png         # 热容热力图 ⚠️ 含载体！
├── HeatCapacity_Oxygen_Comprehensive.png  # 氧系列综合图
└── EnergyFluctuation_vs_Temperature.png   # 能量波动图

⚠️  注意：heat_capacity_per_system.csv 中的热容值包含载体贡献！
    - Cv_meV_per_K = (dE_total/dT) / N_cluster (表观值)
    - 获得纯团簇热容，请运行：python subtract_support_v2.py
    - 扣除载体后数据位于：results/cluster_only_analysis/

关键改进 (v2.1)
===============
1. **4级路径签名算法**
   - 修改前: 890签名 → 900记录 (批次重复,10个run被重复计数)
   - 修改后: 892签名 → 892记录 (完美1:1匹配)
   - 最终结果: 3262 - 892 = 2370 有效记录 (精确过滤)

2. **通用路径提取函数**
   - extract_path_signature(filepath, is_msd_path)
   - 支持MSD路径(有温度目录)和能量路径(无温度目录)
   - 自动检测批次标识符(run3, run2, run4, run5)
   - 向上搜索最多3级目录

3. **命令行参数支持**
   - --no-filter: 跳过MSD异常筛选
   - 解决Pt3Sn5 @ 1100K等体系全被筛选的问题

技术细节
========
路径签名算法逻辑:
1. 从文件路径提取run信息 (T1000.r24.gpu0)
2. 定位关键目录索引 (温度目录或run所在位置)
3. 提取composition和parent目录
4. 向上搜索批次标识符 (最多3级)
5. 构建签名:
   - 有批次: batch/parent/composition/run (4级)
   - 无批次: parent/composition/run (3级)

批次标识符: ['run3', 'run2', 'run4', 'run5']

注意事项
========
1. 路径签名算法要求MSD和能量数据使用相同的目录结构
2. 化学式自动从路径提取,支持多种命名格式
3. Cv系列需要特殊处理 (从g-xxxx-化学式提取)
4. 热容计算使用滑动窗口,至少需要3个温度点

依赖库
======
- pandas: 数据处理
- numpy: 数值计算
- matplotlib: 绘图
- scipy: 插值和统计
- argparse: 命令行参数解析

相关脚本
========
- step7_lindemann_analysis.py: 林德曼指数分析 (使用相同的路径签名算法)
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
- 热容跃变检测
- 系列分类功能

v1.0:
- 基础能量分析
- 热容计算
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from datetime import datetime
from scipy import stats
from scipy.interpolate import UnivariateSpline
import warnings
import re
import argparse
warnings.filterwarnings('ignore')


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
        
        能量路径:
        >>> extract_path_signature(
        ...     "/home/data/run3/o-sorted/o2/O2Pt7Sn7/T200.r0.gpu0",
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
        # 能量路径: 找run所在位置
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


# ===== 全局配置 =====
BASE_DIR = Path(__file__).parent
ENERGY_MASTER = BASE_DIR / 'data' / 'lammps_energy' / 'energy_master_20251016_121110.csv'
OUTLIERS_FILE = BASE_DIR / 'results' / 'large_D_outliers.csv'
D_VALUES_FILE = BASE_DIR / 'results' / 'ensemble_D_analysis' / 'ensemble_D_values.csv'

OUTPUT_DIR = BASE_DIR / 'results' / 'energy_analysis_v2'
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

# ===== 体系过滤配置 (与step3保持一致) =====
# 使用正则表达式过滤要分析的体系
SYSTEM_FILTER = {
    'include_patterns': [
        # ===== 当前分析: 含氧体系 (Pt-Sn-O三元系统) =====
        # 考虑Sn-O成键特性,按O原子数分类
        # 支持O1开头(O1Sn4Pt3)或O1结尾(Pt2Sn2O1)的命名
        r'[Oo]1',          # O1系列: O1xxx, Pt2Sn2O1, Pt3Sn2O1等
        r'[Oo]2',          # O2系列: O2Pt4Sn6, O2Pt7Sn7, Pt3Sn3O2等
        r'[Oo]3',          # O3系列: O3Pt5Sn7, O3Sn4Pt2, Sn7Pt4O3等
        r'[Oo]4',          # O4系列: O4Pt3Sn6, Sn10Pt7O4等
        
        # ===== Pt(8-x)SnX体系 (总原子数=8,变Pt/Sn比) =====
        r'^pt7sn1',        # Pt7Sn1: 7+1=8
        r'^pt6sn2',        # Pt6Sn2: 6+2=8
        r'^pt5sn3',        # Pt5Sn3: 5+3=8
        r'^pt4sn4',        # Pt4Sn4: 4+4=8
        r'^pt3sn5',        # Pt3Sn5: 3+5=8
        r'^pt2sn6',        # Pt2Sn6: 2+6=8
        r'^pt8$',          # Pt8: 纯Pt (用于对比)
        
        # ===== Pt8SnX体系 (固定8个Pt原子,变Sn含量) =====
        r'^pt8sn',         # Pt8SnX体系: pt8sn0, pt8sn1, ..., pt8sn10
        
        # ===== Pt6SnX体系 (固定6个Pt原子,变Sn含量) =====
        r'^pt6',           # Pt6SnX体系: pt6, pt6sn1, pt6sn2, ..., pt6sn9
        
        # ===== Cv系列 (特殊命名规则) =====
        r'^Cv-',           # Cv系列: Cv-1, Cv-2, Cv-3, Cv-4, Cv-5
    ],
    'exclude_patterns': [
        # 含氧分析时不排除
        # r'[Oo]\d+',        # 排除所有含氧系统
    ]
}

# 颜色方案 (11种Sn含量)
SN_COLORS = plt.cm.viridis(np.linspace(0, 1, 11))
SN_COLOR_MAP = {i: SN_COLORS[i] for i in range(11)}

plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.dpi'] = 150


def extract_chemical_formula_from_path(full_path):
    """
    从完整路径中提取化学式
    
    Examples:
    - /home/.../g-1535-Sn8Pt6O4/Cv-1/T1000.r16.gpu0 -> Sn8Pt6O4
    - /home/.../GPU-Pt8/Pt8/pt8sn5-1-best/T800.r14.gpu0 -> pt8sn5-1-best
    - /home/.../pt6/pt6sn3/T700.r5.gpu0 -> pt6sn3
    
    策略: 
    1. 对于Cv系列，从上两级目录提取g-xxxx-化学式
    2. 其他情况，路径倒数第二层通常是化学式目录
    """
    if not full_path:
        return None
    
    # 分割路径
    parts = full_path.replace('\\', '/').split('/')
    
    # 去除空字符串
    parts = [p for p in parts if p]
    
    if len(parts) < 2:
        return None
    
    # 特殊处理Cv系列: .../g-1535-Sn8Pt6O4/Cv-1/...
    # Cv系列的特征：倒数第二层是Cv-X
    if len(parts) >= 3 and parts[-2].startswith('Cv-'):
        parent_dir = parts[-3]  # g-1535-Sn8Pt6O4
        if parent_dir.startswith('g-'):
            # 提取化学式: g-1535-Sn8Pt6O4 -> Sn8Pt6O4
            match = re.search(r'g-\d+-(.+)', parent_dir, re.IGNORECASE)
            if match:
                return match.group(1)
    
    # 倒数第二层: .../化学式目录/Cv-1/...
    parent_dir = parts[-2]
    
    # 检查是否是g-开头的目录(如g-1535-Sn8Pt6O4)
    if parent_dir.startswith('g-'):
        # 提取化学式: g-1535-Sn8Pt6O4 -> Sn8Pt6O4
        match = re.search(r'g-\d+-(.+)', parent_dir, re.IGNORECASE)
        if match:
            return match.group(1)
    
    # 否则返回倒数第二层作为结构名
    return parent_dir


def extract_sn_content(structure_name):
    """
    从结构名提取Sn含量
    
    Examples:
    - pt8sn5-2-best -> 5
    - pt6sn3 -> 3
    - O2Pt4Sn6 -> 6
    - Sn8Pt6O4 -> 8
    - pt6 -> 0 (纯Pt)
    - Cv-1 -> None (无法提取)
    """
    # 匹配 ptXsnY 格式 (pt8sn5, pt6sn3等)
    match = re.search(r'pt\d+sn(\d+)', structure_name, re.IGNORECASE)
    if match:
        return int(match.group(1))
    
    # 匹配 SnXPtY 或 OXSnYPtZ 格式 (Sn8Pt6O4, O2Sn8Pt7等)
    match = re.search(r'sn(\d+)', structure_name, re.IGNORECASE)
    if match:
        return int(match.group(1))
    
    # 匹配纯Pt系列 (pt6, pt7等, Sn含量为0)
    match = re.search(r'^pt\d+$', structure_name, re.IGNORECASE)
    if match:
        return 0
    
    return None


def extract_pt_sn_o_atoms(structure_name):
    """
    从结构名提取Pt、Sn、O原子数
    
    Returns:
    --------
    (n_pt, n_sn, n_o) : tuple or (None, None, None)
    
    Examples:
    - pt8sn5 -> (8, 5, 0)
    - pt6sn2 -> (6, 2, 0)
    - pt6 -> (6, 0, 0)
    - Sn8Pt6O4 -> (6, 8, 4)
    - O2Pt4Sn6 -> (4, 6, 2)
    - Pt5Sn3O1 -> (5, 3, 1)
    """
    structure_lower = structure_name.lower()
    
    # 匹配含氧格式: O?Pt?Sn?O? (各种组合)
    # 例如: O2Pt4Sn6, Sn8Pt6O4, Pt5Sn3O1
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
    
    # 如果没有提取到任何原子数,返回None
    if n_pt == 0 and n_sn == 0 and n_o == 0:
        return (None, None, None)
    
    return (n_pt, n_sn, n_o)


def extract_pt_sn_atoms(structure_name):
    """
    从结构名提取Pt和Sn原子数 (向后兼容,不包含O)
    
    Returns:
    --------
    (n_pt, n_sn) : tuple or (None, None)
    """
    n_pt, n_sn, n_o = extract_pt_sn_o_atoms(structure_name)
    return (n_pt, n_sn)


def classify_system_series(structure_name):
    """
    分类体系所属系列
    
    特殊分类:
    - Pt8SnX: 固定8个Pt原子,改变Sn含量 (总原子数8~18)
    - Pt6SnX: 固定6个Pt原子,改变Sn含量 (总原子数6~16)
    - Pt(8-x)SnX: 固定总原子数=8 (如Pt6Sn2, Pt5Sn3, Pt4Sn4等)
    - O1-O4: 含氧系列,按O原子数分类
    
    Returns:
    --------
    series : str
        系列名称
    """
    structure_lower = structure_name.lower()
    
    # 优先检查含氧系列 (因为含氧体系名称中也包含Pt和Sn)
    # 修改: 支持O在开头或结尾的命名方式
    # O1可能出现在: O1Sn4Pt3, Pt2Sn2O1, Sn1Pt2O1等位置
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
    
    # 提取Pt和Sn原子数 (用于PtSn体系分类)
    n_pt, n_sn = extract_pt_sn_atoms(structure_name)
    
    if n_pt is not None and n_sn is not None:
        total_atoms = n_pt + n_sn
        
        # Pt8SnX体系: 固定8个Pt
        if n_pt == 8:
            return 'Pt8SnX'
        
        # Pt6SnX体系: 固定6个Pt
        if n_pt == 6:
            return 'Pt6SnX'
        
        # Pt(8-x)SnX体系: 总原子数=8
        if total_atoms == 8:
            return 'Pt(8-x)SnX'
        
        # 其他PtSn体系
        if n_pt == 7:
            return 'Pt7SnX'
        if n_pt == 5:
            return 'Pt5SnX'
        if n_pt == 4:
            return 'Pt4SnX'
        if n_pt == 3:
            return 'Pt3SnX'
    
    # 纯Pt系列
    if re.search(r'^pt\d+$', structure_lower):
        return 'Pt-pure'
    
    return 'Other'


def filter_structures_by_pattern(structures, include_patterns=None, exclude_patterns=None):
    """
    根据正则表达式过滤结构列表 (与step3保持一致)
    
    Parameters:
    -----------
    structures : list
        结构名称列表
    include_patterns : list of str
        包含模式列表 (正则表达式)
        如果为空,则包含所有结构
    exclude_patterns : list of str
        排除模式列表 (正则表达式)
    
    Returns:
    --------
    filtered_structures : list
        过滤后的结构列表
    """
    if include_patterns is None:
        include_patterns = []
    if exclude_patterns is None:
        exclude_patterns = []
    
    filtered = []
    
    for structure in structures:
        # 检查是否匹配include模式
        if include_patterns:
            included = any(re.search(pattern, structure, re.IGNORECASE) 
                          for pattern in include_patterns)
            if not included:
                continue
        
        # 检查是否匹配exclude模式
        if exclude_patterns:
            excluded = any(re.search(pattern, structure, re.IGNORECASE) 
                          for pattern in exclude_patterns)
            if excluded:
                continue
        
        # 通过过滤
        filtered.append(structure)
    
    return filtered


def build_path_filter_set(outliers_df):
    r"""
    从outliers构建路径过滤集合 (使用通用路径签名算法)
    
    MSD路径示例:
    1. ...\more\run3\o2\O2Pt4Sn6\1000K\T1000.r24.gpu0_msd_Pt.xvg
       → run3/o2/o2pt4sn6/t1000.r24.gpu0
       
    2. ...\o68\o2-3\O2Pt4Sn6\1000K\T1000.r8.gpu0_msd_Sn.xvg
       → o2-3/o2pt4sn6/t1000.r8.gpu0
       
    3. ...\GPU-Pt8\Pt8\pt8sn0-2-best\700K\T700.r4.gpu0_msd_Pt.xvg
       → pt8/pt8sn0-2-best/t700.r4.gpu0
    
    ✅ 使用通用extract_path_signature函数
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


def should_filter_energy_record(row, filter_path_signature_set):
    """
    判断能量记录是否应该被筛选掉 (使用通用路径签名算法)
    
    ✅ 使用通用extract_path_signature函数
    
    能量记录路径示例:
    1. /home/.../more/run3/o-sorted-dell-o68/o2/O2Pt7Sn7/T200.r0.gpu0
       → run3/o2/o2pt7sn7/t200.r0.gpu0
       
    2. /home/.../o68/o-sorted-dell-o68/o2-3/O2Pt4Sn6/T500.r3.gpu0
       → o2-3/o2pt4sn6/t500.r3.gpu0
       
    3. /home/.../GPU-Pt8/Pt8/pt8sn0-2-best/T700.r4.gpu0
       → pt8/pt8sn0-2-best/t700.r4.gpu0
    """
    full_path = row.get('完整路径', '')
    if not full_path:
        return False
    
    # 使用通用路径签名函数 (能量路径无温度目录)
    path_signature = extract_path_signature(full_path, is_msd_path=False)
    if not path_signature:
        return False
    
    return path_signature in filter_path_signature_set


def load_energy_data(no_filter=False):
    """
    加载能量数据并使用正则表达式过滤体系 (与step3保持一致)
    
    Args:
        no_filter: bool, 如果为True,跳过MSD异常筛选,使用所有数据
    """
    print(f"\n[*] Loading energy data...")
    
    # 1. 加载能量数据
    df_energy = pd.read_csv(ENERGY_MASTER, encoding='utf-8')
    df_energy.columns = ['路径', '结构', '温度', '模拟序号', '总步数', '采样步数', 
                        '平均能量', '标准差', '最小值', '最大值', '采样间隔', 
                        '跳过步数', '完整路径']
    
    print(f"  [OK] Loaded {len(df_energy)} energy records")
    
    # 2. 应用正则表达式过滤
    include_patterns = SYSTEM_FILTER.get('include_patterns', [])
    exclude_patterns = SYSTEM_FILTER.get('exclude_patterns', [])
    
    if include_patterns or exclude_patterns:
        print(f"\n[*] Applying system filters...")
        if include_patterns:
            print(f"    Include patterns: {include_patterns}")
        if exclude_patterns:
            print(f"    Exclude patterns: {exclude_patterns}")
        
        all_structures = df_energy['结构'].unique().tolist()
        filtered_structures = filter_structures_by_pattern(
            all_structures, include_patterns, exclude_patterns
        )
        
        df_energy = df_energy[df_energy['结构'].isin(filtered_structures)].copy()
        print(f"    Filtered: {len(filtered_structures)}/{len(all_structures)} structures selected")
        print(f"    Records: {len(df_energy)}")
    else:
        print(f"  [INFO] No system filter applied, using all structures")
    
    # 3. 提取真实化学式 (从路径上一级目录)
    df_energy['化学式'] = df_energy['完整路径'].apply(extract_chemical_formula_from_path)
    
    # 如果化学式为空,使用结构名
    df_energy['化学式'] = df_energy['化学式'].fillna(df_energy['结构'])
    
    # 3.1 提取Sn含量 (从化学式提取)
    df_energy['Sn含量'] = df_energy['化学式'].apply(extract_sn_content)
    
    # 3.2 分类体系系列
    df_energy['系列'] = df_energy['化学式'].apply(classify_system_series)
    
    # 3.3 过滤掉无法提取Sn含量的结构
    n_total = len(df_energy)
    df_energy = df_energy[df_energy['Sn含量'].notna()].copy()
    n_filtered = len(df_energy)
    if n_total > n_filtered:
        print(f"  [INFO] Filtered out {n_total - n_filtered} records with unknown Sn content")
    
    print(f"\n[*] System series classification:")
    series_counts = df_energy['系列'].value_counts().sort_index()
    for series, count in series_counts.items():
        print(f"    - {series}: {count} records")
    
    # 4. 加载outliers并构建过滤集 (可选)
    if no_filter:
        print("\n[*] [!] Skipping MSD outlier filtering (--no-filter mode)")
        print(f"  [INFO] Using all {len(df_energy)} records without filtering")
    elif OUTLIERS_FILE.exists():
        print("\n[*] Loading outliers for filtering...")
        df_outliers = pd.read_csv(OUTLIERS_FILE)
        print(f"  [OK] Loaded {len(df_outliers)} outlier runs")
        
        filter_path_signature_set = build_path_filter_set(df_outliers)
        
        # 5. 应用过滤 (精确路径匹配)
        df_energy['is_outlier'] = df_energy.apply(
            lambda row: should_filter_energy_record(row, filter_path_signature_set), 
            axis=1
        )
        
        n_outliers = df_energy['is_outlier'].sum()
        print(f"  [FILTER] Identified {n_outliers} outlier records in energy data")
        print(f"    - Using precise path signature matching (避免误判)")
        
        # 移除outliers
        df_energy = df_energy[~df_energy['is_outlier']].copy()
        print(f"  [OK] After filtering: {len(df_energy)} valid records")
    else:
        print(f"  [WARNING] Outliers file not found: {OUTLIERS_FILE}")
        print(f"  [WARNING] Proceeding without outlier filtering")
    
    # 6. 统计信息
    print(f"\n[*] Data summary:")
    print(f"  Structures: {df_energy['结构'].nunique()}")
    print(f"  Sn content range: {df_energy['Sn含量'].min():.0f} - {df_energy['Sn含量'].max():.0f}")
    print(f"  Temperature range: {df_energy['温度'].min():.0f} - {df_energy['温度'].max():.0f} K")
    print(f"  Run IDs: {df_energy['模拟序号'].nunique()}")
    
    return df_energy


def calculate_per_atom_energy(df):
    """计算每原子能量 (智能识别Pt、Sn、O原子数)"""
    print("\n[*] Calculating per-atom energy...")
    
    # 提取Pt、Sn、O原子数
    df[['Pt原子数', 'Sn原子数', 'O原子数']] = df['化学式'].apply(
        lambda x: pd.Series(extract_pt_sn_o_atoms(x))
    )
    
    # 计算总原子数 (包含O)
    df['原子数'] = df['Pt原子数'] + df['Sn原子数'] + df['O原子数']
    
    # 计算每原子能量
    df['每原子能量'] = df['平均能量'] / df['原子数']
    df['每原子标准差'] = df['标准差'] / df['原子数']
    
    # 计算Pt+Sn比例 (不含O)
    df['PtSn原子数'] = df['Pt原子数'] + df['Sn原子数']
    df['Pt比例'] = df['Pt原子数'] / df['PtSn原子数']
    df['Sn比例'] = df['Sn原子数'] / df['PtSn原子数']
    
    print(f"  [OK] Total atoms: {df['原子数'].min():.0f} - {df['原子数'].max():.0f}")
    print(f"  [OK] Pt atoms: {df['Pt原子数'].min():.0f} - {df['Pt原子数'].max():.0f}")
    print(f"  [OK] Sn atoms: {df['Sn原子数'].min():.0f} - {df['Sn原子数'].max():.0f}")
    if df['O原子数'].max() > 0:
        print(f"  [OK] O atoms: {df['O原子数'].min():.0f} - {df['O原子数'].max():.0f}")
    
    return df


def calculate_relative_energy_per_system(df):
    """对每个体系计算相对能量 (E - E_min)"""
    print("\n[*] Calculating relative energy for each system...")
    
    results = []
    
    for structure in df['结构'].unique():
        df_sys = df[df['结构'] == structure].copy()
        
        # 找到最低能量(通常在最低温度)
        E_min = df_sys['每原子能量'].min()
        
        # 计算相对能量
        df_sys['相对能量'] = df_sys['每原子能量'] - E_min
        
        results.append(df_sys)
    
    df_result = pd.concat(results, ignore_index=True)
    print(f"  [OK] Calculated relative energy for {len(df['结构'].unique())} systems")
    
    return df_result


def calculate_heat_capacity_per_system(df):
    """
    对每个体系计算热容 Cv = dE/dT
    
    ⚠️  重要说明：
    ==============
    本函数计算整个体系的总热容（含团簇+载体）：
    
    计算方式：
      1. E_total = LAMMPS的etotal (团簇~60原子 + 载体240原子)
      2. Cv_total = dE_total/dT  (单位：eV/K 或 meV/K)
    
    输出列说明：
      - Cv_total_meV_K: 总热容 (meV/K) - 推荐使用
      - Cv_total_eV_K: 总热容 (eV/K)
      - Cv_per_atom_meV_K: 每团簇原子热容 (表观值，含载体贡献)
      - N_cluster: 团簇原子数
    
    如何获得纯团簇热容：
      使用 subtract_support_v2.py 扣除载体贡献
    
    使用数值微分,单位说明:
    - 输入: eV/atom, K
    - 输出: meV/K (总热容), eV/K (总热容)
    - 主要单位: meV/K (纳米团簇研究常用单位)
    - 转换: 1 eV = 1000 meV
    """
    print("\n[*] Calculating heat capacity for each system...")
    
    cv_results = []
    
    for structure in sorted(df['结构'].unique()):
        df_sys = df[df['结构'] == structure].copy()
        
        # 对每个温度求平均能量(多个run_id)
        agg_dict = {
            '平均能量': 'mean',      # ⚠️ 使用总能量计算总热容
            '每原子能量': 'mean',    # 保留，用于计算每原子热容
            'Sn含量': 'first'
        }
        
        # 添加系列、化学式、原子数
        if '系列' in df_sys.columns:
            agg_dict['系列'] = 'first'
        if '化学式' in df_sys.columns:
            agg_dict['化学式'] = 'first'
        if 'Pt原子数' in df_sys.columns:
            agg_dict['Pt原子数'] = 'first'
        if 'Sn原子数' in df_sys.columns:
            agg_dict['Sn原子数'] = 'first'
        if 'O原子数' in df_sys.columns:
            agg_dict['O原子数'] = 'first'
        
        df_avg = df_sys.groupby('温度').agg(agg_dict).reset_index().sort_values('温度')
        
        temperatures = df_avg['温度'].values
        energies_total = df_avg['平均能量'].values      # 总能量 (eV)
        energies_per_atom = df_avg['每原子能量'].values # 每原子能量 (eV/atom)
        sn_content = df_avg['Sn含量'].iloc[0]
        
        # 获取系列和化学式
        series = df_avg['系列'].iloc[0] if '系列' in df_avg.columns else 'Unknown'
        chemical_formula = df_avg['化学式'].iloc[0] if '化学式' in df_avg.columns else structure
        
        # 获取Pt、Sn、O原子数
        pt_atoms = df_avg['Pt原子数'].iloc[0] if 'Pt原子数' in df_avg.columns else None
        sn_atoms = df_avg['Sn原子数'].iloc[0] if 'Sn原子数' in df_avg.columns else None
        o_atoms = df_avg['O原子数'].iloc[0] if 'O原子数' in df_avg.columns else 0
        
        # 计算团簇总原子数
        n_cluster = 0
        if pt_atoms is not None and sn_atoms is not None:
            n_cluster = pt_atoms + sn_atoms + (o_atoms if o_atoms is not None else 0)
        
        # 需要至少3个温度点才能计算热容
        if len(temperatures) < 3:
            print(f"  [SKIP] {structure}: only {len(temperatures)} temperature points")
            continue
        
        # 使用中心差分计算dE/dT
        cv_total_values = []      # 总热容 (eV/K)
        cv_per_atom_values = []   # 每原子热容 (eV/K/atom)
        cv_temps = []
        
        for i in range(1, len(temperatures) - 1):
            dT_left = temperatures[i] - temperatures[i-1]
            dT_right = temperatures[i+1] - temperatures[i]
            
            # 总能量的变化率 → 总热容
            dE_total_left = energies_total[i] - energies_total[i-1]
            dE_total_right = energies_total[i+1] - energies_total[i]
            
            # 每原子能量的变化率 → 每原子热容
            dE_per_atom_left = energies_per_atom[i] - energies_per_atom[i-1]
            dE_per_atom_right = energies_per_atom[i+1] - energies_per_atom[i]
            
            # 中心差分
            if dT_left > 0 and dT_right > 0:
                cv_total = (dE_total_left/dT_left + dE_total_right/dT_right) / 2
                cv_per_atom = (dE_per_atom_left/dT_left + dE_per_atom_right/dT_right) / 2
                
                cv_total_values.append(cv_total)
                cv_per_atom_values.append(cv_per_atom)
                cv_temps.append(temperatures[i])
        
        # 端点使用前向/后向差分
        if len(temperatures) >= 2:
            # 第一个点 (前向差分)
            cv_total_first = (energies_total[1] - energies_total[0]) / (temperatures[1] - temperatures[0])
            cv_per_atom_first = (energies_per_atom[1] - energies_per_atom[0]) / (temperatures[1] - temperatures[0])
            
            cv_total_values.insert(0, cv_total_first)
            cv_per_atom_values.insert(0, cv_per_atom_first)
            cv_temps.insert(0, temperatures[0])
            
            # 最后一个点 (后向差分)
            cv_total_last = (energies_total[-1] - energies_total[-2]) / (temperatures[-1] - temperatures[-2])
            cv_per_atom_last = (energies_per_atom[-1] - energies_per_atom[-2]) / (temperatures[-1] - temperatures[-2])
            
            cv_total_values.append(cv_total_last)
            cv_per_atom_values.append(cv_per_atom_last)
            cv_temps.append(temperatures[-1])
        
        
        # 保存结果
        for temp, cv_total, cv_per_atom in zip(cv_temps, cv_total_values, cv_per_atom_values):
            # cv_total 是总热容 dE_total/dT，单位 eV/K
            # cv_per_atom 是每原子热容 d(E/N)/dT，单位 eV/K/atom
            # ⚠️ 注意：E_total 包含团簇+载体，约300个原子
            
            result_dict = {
                '结构': structure,
                '化学式': chemical_formula,
                '系列': series,
                'Sn含量': sn_content,
                '温度': temp,
                # ✅ 整个体系的总热容（含团簇+载体，约300个原子）
                'Cv_total_meV_K': cv_total * 1000,           # 总热容 (meV/K)
                'Cv_total_eV_K': cv_total,                   # 总热容 (eV/K)
                # ✅ 每团簇原子的热容（表观值，含载体贡献）
                'Cv_per_atom_meV_K': cv_per_atom * 1000,     # 每原子 (meV/K/atom)
                'Cv_per_atom_eV_K': cv_per_atom,             # 每原子 (eV/K/atom)
                # 为了兼容性，保留旧列名
                'Cv_meV_per_K': cv_per_atom * 1000,          # [已废弃] 每原子热容
                'Cv_eV_per_K': cv_per_atom,                  # [已废弃] 每原子热容
            }
            
            # 添加原子数信息
            if n_cluster > 0:
                result_dict['N_cluster'] = n_cluster         # 团簇原子数
            
            # 添加Pt、Sn、O原子数 (如果有)
            if pt_atoms is not None:
                result_dict['Pt原子数'] = pt_atoms
            if sn_atoms is not None:
                result_dict['Sn原子数'] = sn_atoms
            if o_atoms is not None and o_atoms > 0:
                result_dict['O原子数'] = o_atoms
            
            cv_results.append(result_dict)
    
    df_cv = pd.DataFrame(cv_results)
    print(f"  [OK] Calculated {len(df_cv)} heat capacity values")
    
    return df_cv


def detect_melting_temperature(df_cv, threshold_factor=1.5):
    """
    检测熔化温度 Tm (通过热容跃变)
    
    判断标准:
    1. 热容突然增大 (Cv > threshold_factor * baseline_Cv)
    2. 出现在中高温区域 (T > 500K)
    3. 热容跃变后保持在高值
    
    Parameters:
    -----------
    df_cv : DataFrame
        热容数据,包含列: 结构, 温度, Cv_meV_per_K
    threshold_factor : float
        热容跃变阈值倍数 (默认1.5倍)
    
    Returns:
    --------
    df_tm : DataFrame
        熔化温度结果,包含列: 结构, Sn含量, Tm, 跃变前Cv, 跃变后Cv
    """
    print(f"\n[*] Detecting melting temperature (Cv jump analysis)...")
    print(f"    - Threshold factor: {threshold_factor}x baseline Cv")
    print(f"    - Temperature range: > 500K")
    
    tm_results = []
    
    for structure in sorted(df_cv['结构'].unique()):
        df_sys = df_cv[df_cv['结构'] == structure].copy().sort_values('温度')
        
        if len(df_sys) < 4:
            continue
        
        temperatures = df_sys['温度'].values
        cv_values = df_sys['Cv_meV_per_K'].values
        sn_content = df_sys['Sn含量'].iloc[0]
        
        # 只在T>500K范围内寻找熔化点
        high_temp_mask = temperatures > 500
        if high_temp_mask.sum() < 3:
            continue
        
        temps_high = temperatures[high_temp_mask]
        cv_high = cv_values[high_temp_mask]
        
        # 计算低温区基线Cv (T<500K的平均值)
        low_temp_mask = temperatures <= 500
        if low_temp_mask.sum() > 0:
            baseline_cv = cv_values[low_temp_mask].mean()
        else:
            # 如果没有低温数据,用前3个点的平均
            baseline_cv = cv_values[:3].mean()
        
        # 寻找热容跃变点
        tm_detected = None
        cv_before = None
        cv_after = None
        
        for i in range(1, len(cv_high)):
            if cv_high[i] > threshold_factor * baseline_cv:
                # 找到跃变点
                tm_detected = temps_high[i]
                cv_before = cv_high[i-1]
                cv_after = cv_high[i]
                break
        
        if tm_detected:
            tm_results.append({
                '结构': structure,
                'Sn含量': sn_content,
                'Tm(K)': tm_detected,
                'Cv_before(J/mol·K)': cv_before,
                'Cv_after(J/mol·K)': cv_after,
                'Cv_jump_ratio': cv_after / cv_before if cv_before > 0 else np.nan,
                'baseline_Cv(J/mol·K)': baseline_cv
            })
    
    df_tm = pd.DataFrame(tm_results)
    
    if len(df_tm) > 0:
        print(f"  [OK] Detected melting temperature for {len(df_tm)} systems:")
        for _, row in df_tm.iterrows():
            print(f"    - {row['结构']:20s} Tm = {row['Tm(K)']:4.0f}K "
                  f"(Cv: {row['Cv_before(J/mol·K)']:.2f} -> {row['Cv_after(J/mol·K)']:.2f} J/mol·K, "
                  f"ratio={row['Cv_jump_ratio']:.2f}x)")
    else:
        print(f"  [WARNING] No melting temperature detected")
    
    return df_tm


def plot_energy_vs_temperature_by_system(df, output_dir):
    """按Sn含量分组,绘制每个体系的能量-温度曲线"""
    print("\n[*] Plotting energy vs temperature by system...")
    
    # 判断系列类型和总原子数
    series_name = df['系列'].iloc[0] if '系列' in df.columns else 'Unknown'
    total_atoms = df['原子数'].iloc[0] if '原子数' in df.columns else None
    
    # 构建文件名前缀
    if series_name == 'Pt8SnX':
        file_prefix = 'Pt8SnX'
        title_prefix = 'Pt8Sn'
    elif series_name == 'Pt(8-x)SnX':
        file_prefix = 'Pt(8-x)SnX_TotalAtoms8'
        title_prefix = 'Pt(8-x)Sn'
    else:
        file_prefix = series_name
        title_prefix = series_name
    
    sn_contents = sorted(df['Sn含量'].unique())
    
    for sn in sn_contents:
        df_sn = df[df['Sn含量'] == sn]
        structures = sorted(df_sn['结构'].unique())
        
        if len(structures) == 0:
            continue
        
        # 对于Pt(8-x)SnX体系,提取Pt原子数用于标题
        if series_name == 'Pt(8-x)SnX' and 'Pt原子数' in df_sn.columns:
            pt_count = int(df_sn['Pt原子数'].iloc[0])
            title = f'Pt{pt_count}Sn{int(sn)} (总原子数=8) 能量 vs 温度'
        else:
            title = f'{title_prefix}{int(sn)} 能量 vs 温度 (各体系独立)'
        
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        fig.suptitle(title, fontsize=14, fontweight='bold')
        
        # 左图: 绝对能量
        ax1 = axes[0]
        for struct in structures:
            df_struct = df_sn[df_sn['结构'] == struct].sort_values('温度')
            ax1.plot(df_struct['温度'], df_struct['每原子能量'], 
                    marker='o', label=struct, alpha=0.7)
        
        ax1.set_xlabel('温度 (K)', fontsize=12)
        ax1.set_ylabel('每原子能量 (eV/atom)', fontsize=12)
        ax1.set_title('绝对能量', fontsize=11)
        ax1.legend(fontsize=8, loc='best')
        ax1.grid(True, alpha=0.3)
        
        # 右图: 相对能量
        ax2 = axes[1]
        for struct in structures:
            df_struct = df_sn[df_sn['结构'] == struct].sort_values('温度')
            ax2.plot(df_struct['温度'], df_struct['相对能量'], 
                    marker='s', label=struct, alpha=0.7)
        
        ax2.set_xlabel('温度 (K)', fontsize=12)
        ax2.set_ylabel('相对能量 (eV/atom)', fontsize=12)
        ax2.set_title('相对能量 (E - E_min)', fontsize=11)
        ax2.legend(fontsize=8, loc='best')
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        filename = output_dir / f'{file_prefix}_Energy_vs_T_Sn{int(sn)}.png'
        plt.savefig(filename, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  [OK] Saved: {filename.name}")


def plot_heat_capacity_comparison(df_cv, output_dir, df_tm=None):
    """
    对比不同Sn含量的热容 (按系列分组)
    
    Parameters:
    -----------
    df_cv : DataFrame
        热容数据 (包含'系列'列)
    output_dir : Path
        输出目录
    df_tm : DataFrame, optional
        熔化温度数据,如果提供则在图上标记Tm
    """
    print("\n[*] Plotting heat capacity comparison by series...")
    
    # 获取所有系列
    if '系列' not in df_cv.columns:
        print("  [WARNING] '系列' column not found, skipping series-based plots")
        return
    
    series_list = sorted(df_cv['系列'].unique())
    print(f"    Found series: {series_list}")
    
    # 为每个系列绘制独立的图
    for series in series_list:
        df_series = df_cv[df_cv['系列'] == series].copy()
        
        if len(df_series) == 0:
            continue
        
        # 对于Pt(8-x)SnX体系,按Pt含量分组;其他体系按Sn含量分组
        if series == 'Pt(8-x)SnX':
            # 提取Pt原子数作为分组依据
            if 'Pt原子数' in df_series.columns:
                group_key = 'Pt原子数'
                group_values = sorted(df_series['Pt原子数'].dropna().unique())
                x_label = 'Pt含量'
                title_suffix = 'Pt含量'
            else:
                print(f"  [SKIP] {series}: 'Pt原子数' column not found")
                continue
        elif series in ['O1', 'O2', 'O3', 'O4']:
            # 含氧系列: 按结构分组，展示每个独立的Pt:Sn:O配比
            # 按照总原子数和Pt比例排序
            group_key = '结构'
            
            # 提取结构信息并排序
            struct_info = []
            for struct in df_series['结构'].unique():
                row = df_series[df_series['结构'] == struct].iloc[0]
                pt = row['Pt原子数']
                sn = row['Sn原子数']
                o = row['O原子数']
                total = pt + sn + o
                pt_ratio = pt / (pt + sn) if (pt + sn) > 0 else 0
                struct_info.append((struct, pt, sn, o, total, pt_ratio))
            
            # 排序: 按总原子数, 然后按Pt比例
            struct_info_sorted = sorted(struct_info, key=lambda x: (x[4], x[5]))
            group_values = [x[0] for x in struct_info_sorted]
            
            x_label = '组分配比'
            title_suffix = '不同组分'
        else:
            # 其他体系按Sn含量分组
            group_key = 'Sn含量'
            group_values = sorted(df_series['Sn含量'].unique())
            x_label = 'Sn含量'
            title_suffix = 'Sn含量'
        
        if len(group_values) < 2:
            print(f"  [SKIP] {series}: only {len(group_values)} {x_label}, need at least 2")
            continue
        
        # ========== 创建两张图：总热容 和 每原子热容 ==========
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle(f'热容对比 - {series}系列', fontsize=16, fontweight='bold')
        
        # 上排: 总热容 (Cv_total_meV_K)
        # 下排: 每原子热容 (Cv_per_atom_meV_K)
        
        # 左上图: 总热容 vs 温度
        ax1 = axes[0, 0]
        
        # 为该系列创建颜色映射
        colors = plt.cm.viridis(np.linspace(0, 1, len(group_values)))
        color_map = {val: colors[i] for i, val in enumerate(group_values)}
        
        # 含氧系列: 准备同分异构体标记
        formula_count = {}
        formula_index = {}
        if series in ['O1', 'O2', 'O3', 'O4']:
            # 统计每个化学式的出现次数
            for group_val in group_values:
                row = df_series[df_series[group_key] == group_val].iloc[0]
                pt = int(row['Pt原子数'])
                sn = int(row['Sn原子数'])
                formula = f'Pt{pt}Sn{sn}'
                formula_count[formula] = formula_count.get(formula, 0) + 1
        
        for group_val in group_values:
            df_group = df_series[df_series[group_key] == group_val]
            
            # 对每个分组值,可能有多个体系
            for structure in df_group['结构'].unique():
                df_struct = df_group[df_group['结构'] == structure].sort_values('温度')
                
                # 构建标签 (显示Pt和Sn原子数)
                if 'Pt原子数' in df_struct.columns and 'Sn原子数' in df_struct.columns:
                    n_pt = int(df_struct['Pt原子数'].iloc[0])
                    n_sn = int(df_struct['Sn原子数'].iloc[0])
                    
                    # 含氧系列: 添加同分异构体标记
                    if series in ['O1', 'O2', 'O3', 'O4']:
                        n_o = int(df_struct['O原子数'].iloc[0])
                        formula = f'Pt{n_pt}Sn{n_sn}'
                        
                        # 如果是同分异构体，添加星号
                        if formula_count.get(formula, 0) > 1:
                            idx = formula_index.get(formula, 0)
                            formula_index[formula] = idx + 1
                            markers = ['', '*', '**', '***', '****', '*****']  # 扩展标记列表
                            # 确保idx不超出范围
                            marker = markers[min(idx, len(markers)-1)]
                            label = f'Pt{n_pt}Sn{n_sn}O{n_o}{marker}'
                        else:
                            label = f'Pt{n_pt}Sn{n_sn}O{n_o}'
                    else:
                        label = f'Pt{n_pt}Sn{n_sn}'
                else:
                    label = structure
                
                # 绘制总热容 (上排左图)
                ax1.plot(df_struct['温度'], df_struct['Cv_total_meV_K'], 
                        marker='o', color=color_map[group_val], 
                        alpha=0.6, linewidth=1.5, label=label)
                
                # 标记熔化温度 (如果有)
                if df_tm is not None and len(df_tm) > 0:
                    tm_row = df_tm[df_tm['结构'] == structure]
                    if len(tm_row) > 0:
                        tm = tm_row['Tm(K)'].values[0]
                        ax1.axvline(tm, color=color_map[group_val], 
                                   linestyle='--', alpha=0.3, linewidth=1.5)
                        ax1.text(tm, ax1.get_ylim()[1]*0.95, f'Tm', 
                                fontsize=8, color=color_map[group_val], 
                                rotation=90, va='top', ha='right')
        
        ax1.set_xlabel('温度 (K)', fontsize=12)
        ax1.set_ylabel('总热容 Cv_total [meV/K]', fontsize=12, fontweight='bold', color='darkred')
        ax1.set_title(f'{series}系列: 总热容 vs 温度 (含团簇+载体)', fontsize=11, fontweight='bold')
        ax1.legend(fontsize=8, loc='best', ncol=2)
        ax1.grid(True, alpha=0.3)
        
        # 左下图: 每原子热容 vs 温度
        ax3 = axes[1, 0]
        
        # 重新绘制每原子热容
        for group_val in group_values:
            df_group = df_series[df_series[group_key] == group_val]
            
            for structure in df_group['结构'].unique():
                df_struct = df_group[df_group['结构'] == structure].sort_values('温度')
                
                # 构建标签 (与上图相同)
                if 'Pt原子数' in df_struct.columns and 'Sn原子数' in df_struct.columns:
                    n_pt = int(df_struct['Pt原子数'].iloc[0])
                    n_sn = int(df_struct['Sn原子数'].iloc[0])
                    
                    if series in ['O1', 'O2', 'O3', 'O4']:
                        n_o = int(df_struct['O原子数'].iloc[0])
                        formula = f'Pt{n_pt}Sn{n_sn}'
                        if formula_count.get(formula, 0) > 1:
                            idx = formula_index.get(formula, 0)
                            markers = ['', '*', '**', '***', '****', '*****']
                            marker = markers[min(idx, len(markers)-1)]
                            label = f'Pt{n_pt}Sn{n_sn}O{n_o}{marker}'
                        else:
                            label = f'Pt{n_pt}Sn{n_sn}O{n_o}'
                    else:
                        label = f'Pt{n_pt}Sn{n_sn}'
                else:
                    label = structure
                
                # 绘制每原子热容
                ax3.plot(df_struct['温度'], df_struct['Cv_per_atom_meV_K'], 
                        marker='s', color=color_map[group_val], 
                        alpha=0.6, linewidth=1.5, label=label)
        
        ax3.set_xlabel('温度 (K)', fontsize=12)
        ax3.set_ylabel('每原子热容 Cv/N [meV/K/atom]', fontsize=12, fontweight='bold', color='darkblue')
        ax3.set_title(f'{series}系列: 每原子热容 vs 温度 (表观值，含载体)', fontsize=11, fontweight='bold')
        ax3.legend(fontsize=8, loc='best', ncol=2)
        ax3.grid(True, alpha=0.3)
        
        # 右上图: 平均总热容 vs 分组键 (不同温度)
        ax2 = axes[0, 1]
        
        # 选择几个代表性温度
        selected_temps = [300, 500, 700, 900, 1100]
        temp_colors = plt.cm.coolwarm(np.linspace(0, 1, len(selected_temps)))
        
        # 右上图: 总热容
        for i, temp in enumerate(selected_temps):
            avg_cv_total = []
            group_list = []
            x_positions = []
            
            for idx, group_val in enumerate(group_values):
                df_temp = df_series[(df_series[group_key] == group_val) & 
                               (np.abs(df_series['温度'] - temp) <= 50)]
                if len(df_temp) > 0:
                    avg_cv_total.append(df_temp['Cv_total_meV_K'].mean())
                    group_list.append(group_val)
                    
                    if series in ['O1', 'O2', 'O3', 'O4']:
                        x_positions.append(idx)
                    else:
                        x_positions.append(group_val)
            
            if len(group_list) > 0:
                ax2.plot(x_positions, avg_cv_total, marker='o', color=temp_colors[i],
                        linewidth=2, markersize=8, label=f'{temp}K')
        
        # 含氧系列: 设置x轴刻度标签
        if series in ['O1', 'O2', 'O3', 'O4']:
            ax2.set_xticks(range(len(group_values)))
            # 构建标签: PtXSnY，同分异构体用*标记
            x_tick_labels = []
            formula_count = {}  # 统计每个化学式的出现次数
            
            # 第一遍: 统计化学式
            for struct in group_values:
                row = df_series[df_series['结构'] == struct].iloc[0]
                pt = int(row['Pt原子数'])
                sn = int(row['Sn原子数'])
                formula = f'Pt{pt}Sn{sn}'
                formula_count[formula] = formula_count.get(formula, 0) + 1
            
            # 第二遍: 生成标签
            formula_index = {}  # 记录每个化学式的当前索引
            for struct in group_values:
                row = df_series[df_series['结构'] == struct].iloc[0]
                pt = int(row['Pt原子数'])
                sn = int(row['Sn原子数'])
                formula = f'Pt{pt}Sn{sn}'
                
                # 如果是同分异构体，添加标记
                if formula_count[formula] > 1:
                    idx = formula_index.get(formula, 0)
                    formula_index[formula] = idx + 1
                    # 使用罗马数字或字母标记
                    markers = ['', '*', '**', '***', '****', '*****', '******']
                    # 确保idx不超出范围
                    marker = markers[min(idx, len(markers)-1)]
                    x_tick_labels.append(f'{formula}{marker}')
                else:
                    x_tick_labels.append(formula)
            
            ax2.set_xticklabels(x_tick_labels, rotation=45, ha='right', fontsize=9)
        
        ax2.set_xlabel(x_label, fontsize=12)
        ax2.set_ylabel('平均总热容 [meV/K]', fontsize=12, fontweight='bold', color='darkred')
        ax2.set_title(f'{series}系列: 平均总热容 vs {title_suffix}', fontsize=11, fontweight='bold')
        ax2.legend(fontsize=9, loc='best')
        ax2.grid(True, alpha=0.3)
        
        # 右下图: 每原子热容平均值
        ax4 = axes[1, 1]
        
        for i, temp in enumerate(selected_temps):
            avg_cv_per_atom = []
            group_list = []
            x_positions = []
            
            for idx, group_val in enumerate(group_values):
                df_temp = df_series[(df_series[group_key] == group_val) & 
                               (np.abs(df_series['温度'] - temp) <= 50)]
                if len(df_temp) > 0:
                    avg_cv_per_atom.append(df_temp['Cv_per_atom_meV_K'].mean())
                    group_list.append(group_val)
                    
                    if series in ['O1', 'O2', 'O3', 'O4']:
                        x_positions.append(idx)
                    else:
                        x_positions.append(group_val)
            
            if len(group_list) > 0:
                ax4.plot(x_positions, avg_cv_per_atom, marker='s', color=temp_colors[i],
                        linewidth=2, markersize=8, label=f'{temp}K')
        
        # 含氧系列: 设置x轴刻度标签 (右下图)
        if series in ['O1', 'O2', 'O3', 'O4']:
            ax4.set_xticks(range(len(group_values)))
            ax4.set_xticklabels(x_tick_labels, rotation=45, ha='right', fontsize=9)
        
        ax4.set_xlabel(x_label, fontsize=12)
        ax4.set_ylabel('平均每原子热容 [meV/K/atom]', fontsize=12, fontweight='bold', color='darkblue')
        ax4.set_title(f'{series}系列: 平均每原子热容 vs {title_suffix}', fontsize=11, fontweight='bold')
        ax4.legend(fontsize=9, loc='best')
        ax4.grid(True, alpha=0.3)
        
        plt.tight_layout()
        filename = output_dir / f'HeatCapacity_comparison_{series}.png'
        plt.savefig(filename, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  [OK] Saved: {filename.name}")


def plot_heat_capacity_heatmap(df_cv, output_dir):
    """
    绘制热容热熔图(Heatmap): 展示热容随温度和组分的二维分布
    
    Parameters:
    -----------
    df_cv : DataFrame
        热容数据,包含列: 系列, 结构, 温度, Cv_meV_per_K, Pt原子数, Sn原子数
    output_dir : Path
        输出目录
    """
    print("\n[*] Plotting heat capacity heatmap...")
    
    if len(df_cv) == 0:
        print("  [WARNING] No heat capacity data available for heatmap")
        return
    
    # 按系列分组绘制
    series_list = df_cv['系列'].unique()
    
    for series in series_list:
        df_series = df_cv[df_cv['系列'] == series].copy()
        
        if len(df_series) == 0:
            continue
        
        print(f"  [*] Generating heatmap for series: {series}")
        
        # 确定分组键 (与plot_heat_capacity_comparison保持一致)
        if series == 'Pt(8-x)SnX':
            group_key = 'Pt原子数'
            x_label = 'Pt含量'
            title_suffix = 'Pt含量'
        elif series == 'Pt8SnX':
            group_key = 'Sn原子数'
            x_label = 'Sn含量'
            title_suffix = 'Sn含量'
        elif series in ['O1', 'O2', 'O3', 'O4']:
            # 含氧系列: 按结构分组，展示完整的Pt:Sn:O配比
            # 因为同一Sn含量可能有不同的Pt含量，需要区分
            group_key = '结构'
            x_label = '组分(Pt-Sn-O)'
            title_suffix = '组分配比'
        else:
            # 默认按Sn含量
            if 'Sn原子数' in df_series.columns:
                group_key = 'Sn原子数'
                x_label = 'Sn含量'
                title_suffix = 'Sn含量'
            else:
                print(f"  [WARNING] Cannot determine group key for series {series}")
                continue
        
        # 准备数据透视表: 行=温度, 列=组分
        # 对每个(温度, 组分)组合,计算平均Cv (两种：总热容和每原子热容)
        pivot_data_total = df_series.groupby(['温度', group_key])['Cv_total_meV_K'].mean().reset_index()
        heatmap_data_total = pivot_data_total.pivot(index='温度', columns=group_key, values='Cv_total_meV_K')
        
        pivot_data_per_atom = df_series.groupby(['温度', group_key])['Cv_per_atom_meV_K'].mean().reset_index()
        heatmap_data_per_atom = pivot_data_per_atom.pivot(index='温度', columns=group_key, values='Cv_per_atom_meV_K')
        
        # 含氧系列: 按Pt:Sn比例排序列
        if series in ['O1', 'O2', 'O3', 'O4']:
            # 获取每个结构的Pt、Sn、O数量
            struct_info = []
            for struct in heatmap_data_total.columns:
                row = df_series[df_series['结构'] == struct].iloc[0]
                pt = int(row['Pt原子数'])
                sn = int(row['Sn原子数'])
                o = int(row['O原子数'])
                total = pt + sn + o
                pt_ratio = pt / (pt + sn) if (pt + sn) > 0 else 0  # Pt在金属中的比例
                struct_info.append((struct, pt, sn, o, total, pt_ratio))
            
            # 排序: 优先按总原子数，其次按Pt比例
            struct_info_sorted = sorted(struct_info, key=lambda x: (x[4], x[5]))
            sorted_columns = [s[0] for s in struct_info_sorted]
            
            # 重新排列列顺序
            heatmap_data_total = heatmap_data_total[sorted_columns]
            heatmap_data_per_atom = heatmap_data_per_atom[sorted_columns]
        
        # 检查数据
        if heatmap_data_total.empty:
            print(f"  [WARNING] No valid data for heatmap in series {series}")
            continue
        
        # 创建图形：两个子图，上方为总热容，下方为每原子热容
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 14))
        fig.suptitle(f'{series}系列 热容热力图对比', fontsize=16, fontweight='bold')
        
        # 绘制热力图
        import matplotlib.colors as mcolors
        
        # 自定义colormap: 蓝色(低Cv) -> 绿色 -> 黄色 -> 红色(高Cv)
        cmap = plt.cm.RdYlBu_r  # 反转: 红=高温/高Cv, 蓝=低温/低Cv
        
        # === 上图: 总热容 ===
        im1 = ax1.imshow(heatmap_data_total.values, aspect='auto', cmap=cmap, 
                        interpolation='bilinear', origin='lower')
        
        # 设置坐标轴刻度和标签
        ax1.set_xticks(np.arange(len(heatmap_data_total.columns)))
        ax1.set_yticks(np.arange(len(heatmap_data_total.index)))
        
        # 根据group_key构建标签
        if group_key == 'Pt原子数':
            # Pt(8-x)SnX: 显示PtXSnY
            x_labels = []
            for pt_num in heatmap_data_total.columns:
                row = df_series[df_series['Pt原子数'] == pt_num].iloc[0]
                sn_num = int(row['Sn原子数'])
                x_labels.append(f'Pt{int(pt_num)}Sn{sn_num}')
        elif group_key == '结构':
            # 含氧系列: 按结构显示组分 PtXSnY (同分异构体用*标记)
            x_labels = []
            formula_count = {}
            
            # 第一遍: 统计化学式
            for struct_name in heatmap_data_total.columns:
                row = df_series[df_series['结构'] == struct_name].iloc[0]
                pt_num = int(row['Pt原子数'])
                sn_num = int(row['Sn原子数'])
                formula = f'Pt{pt_num}Sn{sn_num}'
                formula_count[formula] = formula_count.get(formula, 0) + 1
            
            # 第二遍: 生成标签
            formula_index = {}
            for struct_name in heatmap_data_total.columns:
                row = df_series[df_series['结构'] == struct_name].iloc[0]
                pt_num = int(row['Pt原子数'])
                sn_num = int(row['Sn原子数'])
                formula = f'Pt{pt_num}Sn{sn_num}'
                
                if formula_count[formula] > 1:
                    idx = formula_index.get(formula, 0)
                    formula_index[formula] = idx + 1
                    markers = ['', '*', '**', '***', '****', '*****', '******']
                    marker = markers[min(idx, len(markers)-1)]
                    x_labels.append(f'{formula}{marker}')
                else:
                    x_labels.append(formula)
        else:
            # 其他系列: 直接显示Sn含量
            x_labels = [f'Sn{int(val)}' for val in heatmap_data_total.columns]
        
        # === 上图: 总热容热力图 ===
        ax1.set_xticklabels(x_labels, fontsize=10)
        ax1.set_yticklabels([f'{int(temp)}K' for temp in heatmap_data_total.index], fontsize=10)
        plt.setp(ax1.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
        
        # 标注数值
        for i in range(len(heatmap_data_total.index)):
            for j in range(len(heatmap_data_total.columns)):
                value = heatmap_data_total.values[i, j]
                if not np.isnan(value):
                    text_color = 'white' if value < heatmap_data_total.values.mean() else 'black'
                    ax1.text(j, i, f'{value:.1f}',
                           ha="center", va="center", color=text_color,
                           fontsize=8, fontweight='bold')
        
        ax1.set_xlabel(x_label, fontsize=12, fontweight='bold')
        ax1.set_ylabel('温度', fontsize=12, fontweight='bold')
        ax1.set_title(f'总热容 Cv_total [meV/K] (含团簇+载体)', fontsize=13, fontweight='bold', color='darkred')
        
        cbar1 = plt.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
        cbar1.set_label('Cv_total [meV/K]', fontsize=11, fontweight='bold')
        
        # === 下图: 每原子热容热力图 ===
        im2 = ax2.imshow(heatmap_data_per_atom.values, aspect='auto', cmap=cmap, 
                        interpolation='bilinear', origin='lower')
        
        ax2.set_xticks(np.arange(len(heatmap_data_per_atom.columns)))
        ax2.set_yticks(np.arange(len(heatmap_data_per_atom.index)))
        ax2.set_xticklabels(x_labels, fontsize=10)
        ax2.set_yticklabels([f'{int(temp)}K' for temp in heatmap_data_per_atom.index], fontsize=10)
        plt.setp(ax2.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
        
        # 标注数值
        for i in range(len(heatmap_data_per_atom.index)):
            for j in range(len(heatmap_data_per_atom.columns)):
                value = heatmap_data_per_atom.values[i, j]
                if not np.isnan(value):
                    text_color = 'white' if value < heatmap_data_per_atom.values.mean() else 'black'
                    ax2.text(j, i, f'{value:.2f}',
                           ha="center", va="center", color=text_color,
                           fontsize=8, fontweight='bold')
        
        ax2.set_xlabel(x_label, fontsize=12, fontweight='bold')
        ax2.set_ylabel('温度', fontsize=12, fontweight='bold')
        ax2.set_title(f'每原子热容 Cv/N [meV/K/atom] (表观值，含载体)', fontsize=13, fontweight='bold', color='darkblue')
        
        cbar2 = plt.colorbar(im2, ax=ax2, fraction=0.046, pad=0.04)
        cbar2.set_label('Cv/N [meV/K/atom]', fontsize=11, fontweight='bold')
        
        plt.tight_layout()
        
        # 保存文件
        filename = output_dir / f'HeatCapacity_heatmap_{series}.png'
        plt.savefig(filename, dpi=200, bbox_inches='tight')
        plt.close()
        print(f"  [OK] Saved: {filename.name}")


def plot_heat_capacity_heatmap_cluster_only(df_cv, output_dir):
    """
    绘制团簇热容热熔图(扣除载体后): 展示纯团簇热容随温度和组分的二维分布
    
    Parameters:
    -----------
    df_cv : DataFrame
        热容数据,包含列: 系列, 结构, 温度, Cv_total_meV_K, N_cluster
    output_dir : Path
        输出目录
    """
    print("\n[*] Plotting cluster heat capacity heatmap (support subtracted)...")
    
    # 载体热容 (更新为CP2K直接计算值)
    SUPPORT_CV_TOTAL = 38.17  # meV/K (from CP2K calculation, R²=0.999998)
    
    # 计算团簇热容
    df = df_cv.copy()
    df['Cv_cluster_meV_K'] = df['Cv_total_meV_K'] - SUPPORT_CV_TOTAL
    df['Cv_cluster_per_atom_meV_K'] = df['Cv_cluster_meV_K'] / df['N_cluster']
    
    # 过滤负值
    df_valid = df[df['Cv_cluster_meV_K'] > 0].copy()
    
    if len(df_valid) == 0:
        print("  [WARNING] No valid data for heatmap after support subtraction")
        return
    
    print(f"  [INFO] Valid records: {len(df_valid)}/{len(df)} for heatmap")
    
    # 经典极限
    kB = 0.0862  # meV/K
    CLASSICAL_CV_PER_ATOM = 3 * kB
    
    # 按系列分组绘制
    series_list = df_valid['系列'].unique()
    
    for series in series_list:
        df_series = df_valid[df_valid['系列'] == series].copy()
        
        if len(df_series) == 0:
            continue
        
        print(f"  [*] Generating cluster heatmap for series: {series}")
        
        # 确定分组键 (与原函数保持一致)
        if series == 'Pt(8-x)SnX':
            group_key = 'Pt原子数'
            x_label = 'Pt含量'
        elif series == 'Pt8SnX':
            group_key = 'Sn原子数'
            x_label = 'Sn含量'
        elif series in ['O1', 'O2', 'O3', 'O4']:
            group_key = '结构'
            x_label = '组分(Pt-Sn-O)'
        else:
            if 'Sn原子数' in df_series.columns:
                group_key = 'Sn原子数'
                x_label = 'Sn含量'
            else:
                continue
        
        # 准备数据透视表
        pivot_data_total = df_series.groupby(['温度', group_key])['Cv_cluster_meV_K'].mean().reset_index()
        heatmap_data_total = pivot_data_total.pivot(index='温度', columns=group_key, values='Cv_cluster_meV_K')
        
        pivot_data_per_atom = df_series.groupby(['温度', group_key])['Cv_cluster_per_atom_meV_K'].mean().reset_index()
        heatmap_data_per_atom = pivot_data_per_atom.pivot(index='温度', columns=group_key, values='Cv_cluster_per_atom_meV_K')
        
        # 含氧系列: 按Pt:Sn比例排序列
        if series in ['O1', 'O2', 'O3', 'O4']:
            struct_info = []
            for struct in heatmap_data_total.columns:
                row = df_series[df_series['结构'] == struct].iloc[0]
                pt = int(row['Pt原子数'])
                sn = int(row['Sn原子数'])
                o = int(row['O原子数'])
                total = pt + sn + o
                pt_ratio = pt / (pt + sn) if (pt + sn) > 0 else 0
                struct_info.append((struct, pt, sn, o, total, pt_ratio))
            
            struct_info_sorted = sorted(struct_info, key=lambda x: (x[4], x[5]))
            sorted_columns = [s[0] for s in struct_info_sorted]
            
            heatmap_data_total = heatmap_data_total[sorted_columns]
            heatmap_data_per_atom = heatmap_data_per_atom[sorted_columns]
        
        # 创建图表 (1×2)
        fig, axes = plt.subplots(1, 2, figsize=(16, 8))
        fig.suptitle(f'团簇热容热力图（已扣除载体 {SUPPORT_CV_TOTAL} meV/K）- {series}系列', 
                     fontsize=16, fontweight='bold')
        
        ax1, ax2 = axes
        
        # 准备X轴标签
        if series in ['O1', 'O2', 'O3', 'O4']:
            x_labels = []
            for struct in heatmap_data_total.columns:
                row = df_series[df_series['结构'] == struct].iloc[0]
                pt = int(row['Pt原子数'])
                sn = int(row['Sn原子数'])
                o = int(row['O原子数'])
                x_labels.append(f'Pt{pt}Sn{sn}O{o}')
        else:
            x_labels = [str(c) for c in heatmap_data_total.columns]
        
        # 颜色映射
        cmap = 'YlOrRd'
        
        # === 左图: 团簇总热容热力图 ===
        im1 = ax1.imshow(heatmap_data_total.values, aspect='auto', cmap=cmap,
                        interpolation='bilinear', origin='lower')
        
        ax1.set_xticks(np.arange(len(heatmap_data_total.columns)))
        ax1.set_yticks(np.arange(len(heatmap_data_total.index)))
        ax1.set_xticklabels(x_labels, fontsize=10)
        ax1.set_yticklabels([f'{int(temp)}K' for temp in heatmap_data_total.index], fontsize=10)
        plt.setp(ax1.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
        
        # 标注数值
        for i in range(len(heatmap_data_total.index)):
            for j in range(len(heatmap_data_total.columns)):
                value = heatmap_data_total.values[i, j]
                if not np.isnan(value):
                    text_color = 'white' if value < heatmap_data_total.values.mean() else 'black'
                    ax1.text(j, i, f'{value:.1f}',
                           ha="center", va="center", color=text_color,
                           fontsize=8, fontweight='bold')
        
        ax1.set_xlabel(x_label, fontsize=12, fontweight='bold')
        ax1.set_ylabel('温度', fontsize=12, fontweight='bold')
        ax1.set_title(f'团簇总热容 Cv_cluster [meV/K]', fontsize=13, fontweight='bold', color='darkred')
        
        cbar1 = plt.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
        cbar1.set_label('Cv_cluster [meV/K]', fontsize=11, fontweight='bold')
        
        # === 右图: 团簇每原子热容热力图 ===
        im2 = ax2.imshow(heatmap_data_per_atom.values, aspect='auto', cmap='RdYlGn_r',
                        interpolation='bilinear', origin='lower',
                        vmin=0, vmax=CLASSICAL_CV_PER_ATOM*2)
        
        ax2.set_xticks(np.arange(len(heatmap_data_per_atom.columns)))
        ax2.set_yticks(np.arange(len(heatmap_data_per_atom.index)))
        ax2.set_xticklabels(x_labels, fontsize=10)
        ax2.set_yticklabels([f'{int(temp)}K' for temp in heatmap_data_per_atom.index], fontsize=10)
        plt.setp(ax2.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
        
        # 标注数值和经典极限比较
        for i in range(len(heatmap_data_per_atom.index)):
            for j in range(len(heatmap_data_per_atom.columns)):
                value = heatmap_data_per_atom.values[i, j]
                if not np.isnan(value):
                    # 根据与经典极限的比较选择文字颜色
                    ratio = value / CLASSICAL_CV_PER_ATOM
                    text_color = 'white' if ratio > 1.5 else 'black'
                    ax2.text(j, i, f'{value:.2f}',
                           ha="center", va="center", color=text_color,
                           fontsize=8, fontweight='bold')
        
        ax2.set_xlabel(x_label, fontsize=12, fontweight='bold')
        ax2.set_ylabel('温度', fontsize=12, fontweight='bold')
        ax2.set_title(f'团簇每原子热容 Cv/N [meV/K/atom]\n经典极限={CLASSICAL_CV_PER_ATOM:.3f}', 
                     fontsize=13, fontweight='bold', color='darkblue')
        
        cbar2 = plt.colorbar(im2, ax=ax2, fraction=0.046, pad=0.04)
        cbar2.set_label('Cv/N [meV/K/atom]', fontsize=11, fontweight='bold')
        
        plt.tight_layout()
        
        # 保存文件
        filename = output_dir / f'ClusterHeatCapacity_heatmap_{series}.png'
        plt.savefig(filename, dpi=200, bbox_inches='tight')
        plt.close()
        print(f"  [OK] Saved: {filename.name}")


def plot_oxygen_series_comprehensive(df_cv, output_dir):
    """
    绘制O1-O4系列的综合对比图
    展示不同O含量对Cv的影响，考虑Pt:Sn:O比值
    
    Parameters:
    -----------
    df_cv : DataFrame
        热容数据
    output_dir : Path
        输出目录
    """
    print("\n[*] Plotting comprehensive oxygen series comparison...")
    
    # 筛选含氧系列
    oxygen_series = df_cv[df_cv['系列'].isin(['O1', 'O2', 'O3', 'O4'])].copy()
    
    if len(oxygen_series) == 0:
        print("  [WARNING] No oxygen series data available")
        return
    
    # 创建2x2子图
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('含氧体系综合对比 (O1-O4系列)\n考虑Pt-Sn-O三元配比', 
                 fontsize=16, fontweight='bold', y=0.995)
    
    # 图1: Cv vs 总原子数 (不同O含量)
    ax1 = axes[0, 0]
    o_colors = {'O1': '#1f77b4', 'O2': '#ff7f0e', 'O3': '#2ca02c', 'O4': '#d62728'}
    
    for o_series in ['O1', 'O2', 'O3', 'O4']:
        df_o = oxygen_series[oxygen_series['系列'] == o_series]
        if len(df_o) == 0:
            continue
        
        # 计算总原子数
        df_o['总原子数'] = df_o['Pt原子数'] + df_o['Sn原子数'] + df_o['O原子数']
        
        # 按总原子数分组，计算平均Cv (300K附近)
        df_300k = df_o[np.abs(df_o['温度'] - 300) <= 50]
        avg_data = df_300k.groupby('总原子数')['Cv_meV_per_K'].mean()
        
        ax1.plot(avg_data.index, avg_data.values, 'o-', 
                color=o_colors[o_series], linewidth=2, markersize=8, 
                label=f'{o_series} ({int(o_series[1])}个O)')
    
    ax1.set_xlabel('总原子数 (Pt+Sn+O)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Cv [meV/K per atom] @ 300K', fontsize=12, fontweight='bold')
    ax1.set_title('热容 vs 总原子数', fontsize=13, fontweight='bold')
    ax1.legend(fontsize=10, loc='best')
    ax1.grid(True, alpha=0.3)
    
    # 图2: Cv vs Pt比例 (Pt/(Pt+Sn))
    ax2 = axes[0, 1]
    
    for o_series in ['O1', 'O2', 'O3', 'O4']:
        df_o = oxygen_series[oxygen_series['系列'] == o_series]
        if len(df_o) == 0:
            continue
        
        # 计算Pt比例
        df_o['Pt比例'] = df_o['Pt原子数'] / (df_o['Pt原子数'] + df_o['Sn原子数'])
        
        # 300K附近的数据
        df_300k = df_o[np.abs(df_o['温度'] - 300) <= 50]
        
        # 按结构分组
        for struct in df_300k['结构'].unique():
            df_struct = df_300k[df_300k['结构'] == struct]
            pt_ratio = df_struct['Pt比例'].iloc[0]
            cv = df_struct['Cv_meV_per_K'].mean()
            ax2.scatter(pt_ratio, cv, color=o_colors[o_series], s=100, alpha=0.7)
    
    ax2.set_xlabel('Pt比例 (Pt/(Pt+Sn))', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Cv [meV/K per atom] @ 300K', fontsize=12, fontweight='bold')
    ax2.set_title('热容 vs Pt/Sn比例', fontsize=13, fontweight='bold')
    ax2.legend([plt.Line2D([0], [0], marker='o', color='w', 
                           markerfacecolor=o_colors[s], markersize=10) 
               for s in ['O1', 'O2', 'O3', 'O4']], 
              ['O1', 'O2', 'O3', 'O4'], fontsize=10, loc='best')
    ax2.grid(True, alpha=0.3)
    
    # 图3: Cv vs 温度 (不同O含量的平均趋势)
    ax3 = axes[1, 0]
    
    for o_series in ['O1', 'O2', 'O3', 'O4']:
        df_o = oxygen_series[oxygen_series['系列'] == o_series]
        if len(df_o) == 0:
            continue
        
        # 按温度计算平均Cv
        avg_cv = df_o.groupby('温度')['Cv_meV_per_K'].mean()
        std_cv = df_o.groupby('温度')['Cv_meV_per_K'].std()
        
        ax3.plot(avg_cv.index, avg_cv.values, 'o-', 
                color=o_colors[o_series], linewidth=2, markersize=6, 
                label=f'{o_series}')
        ax3.fill_between(avg_cv.index, 
                        avg_cv.values - std_cv.values, 
                        avg_cv.values + std_cv.values, 
                        color=o_colors[o_series], alpha=0.2)
    
    ax3.set_xlabel('温度 (K)', fontsize=12, fontweight='bold')
    ax3.set_ylabel('平均Cv [meV/K per atom]', fontsize=12, fontweight='bold')
    ax3.set_title('热容 vs 温度 (各系列平均)', fontsize=13, fontweight='bold')
    ax3.legend(fontsize=10, loc='best')
    ax3.grid(True, alpha=0.3)
    
    # 图4: Cv热熔图 (O含量 vs 总原子数)
    ax4 = axes[1, 1]
    
    # 准备数据: 每个(O含量, 总原子数)组合的平均Cv
    oxygen_series['O含量'] = oxygen_series['系列'].str[1].astype(int)
    oxygen_series['总原子数'] = (oxygen_series['Pt原子数'] + 
                               oxygen_series['Sn原子数'] + 
                               oxygen_series['O原子数'])
    
    # 在300K附近计算
    df_300k = oxygen_series[np.abs(oxygen_series['温度'] - 300) <= 50]
    pivot_comprehensive = df_300k.groupby(['O含量', '总原子数'])['Cv_meV_per_K'].mean().reset_index()
    heatmap_comprehensive = pivot_comprehensive.pivot(index='O含量', columns='总原子数', values='Cv_meV_per_K')
    
    # 绘制热熔图
    im4 = ax4.imshow(heatmap_comprehensive.values, aspect='auto', cmap=plt.cm.RdYlBu_r,
                     interpolation='bilinear', origin='lower')
    
    ax4.set_xticks(np.arange(len(heatmap_comprehensive.columns)))
    ax4.set_yticks(np.arange(len(heatmap_comprehensive.index)))
    ax4.set_xticklabels([f'{int(x)}' for x in heatmap_comprehensive.columns], fontsize=10)
    ax4.set_yticklabels([f'O{int(y)}' for y in heatmap_comprehensive.index], fontsize=10)
    
    # 标注数值
    for i in range(len(heatmap_comprehensive.index)):
        for j in range(len(heatmap_comprehensive.columns)):
            value = heatmap_comprehensive.values[i, j]
            if not np.isnan(value):
                text_color = 'white' if value < heatmap_comprehensive.values[~np.isnan(heatmap_comprehensive.values)].mean() else 'black'
                ax4.text(j, i, f'{value:.2f}', ha="center", va="center", 
                        color=text_color, fontsize=9, fontweight='bold')
    
    ax4.set_xlabel('总原子数', fontsize=12, fontweight='bold')
    ax4.set_ylabel('O原子数', fontsize=12, fontweight='bold')
    ax4.set_title('Cv热熔图 @ 300K', fontsize=13, fontweight='bold')
    
    # 添加颜色条
    cbar4 = plt.colorbar(im4, ax=ax4, fraction=0.046, pad=0.04)
    cbar4.set_label('Cv [J/(mol·K)]', fontsize=11, fontweight='bold')
    
    plt.tight_layout()
    
    # 保存
    filename = output_dir / 'HeatCapacity_Oxygen_Comprehensive.png'
    plt.savefig(filename, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  [OK] Saved: {filename.name}")


def plot_energy_fluctuation(df, output_dir):
    """绘制能量波动(标准差) vs 温度"""
    print("\n[*] Plotting energy fluctuation...")
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    sn_contents = sorted(df['Sn含量'].unique())
    
    for sn in sn_contents:
        df_sn = df[df['Sn含量'] == sn]
        
        # 对每个Sn含量,按温度分组求平均标准差
        fluctuation_data = []
        for temp in sorted(df_sn['温度'].unique()):
            df_temp = df_sn[df_sn['温度'] == temp]
            avg_std = df_temp['每原子标准差'].mean()
            fluctuation_data.append((temp, avg_std))
        
        if len(fluctuation_data) > 0:
            temps, stds = zip(*fluctuation_data)
            ax.plot(temps, stds, marker='o', color=SN_COLOR_MAP[sn],
                   linewidth=2, markersize=6, label=f'Sn{sn}')
    
    ax.set_xlabel('温度 (K)', fontsize=12)
    ax.set_ylabel('能量标准差 (eV/atom)', fontsize=12)
    ax.set_title('能量波动 vs 温度 (各Sn含量)', fontsize=13, fontweight='bold')
    ax.legend(fontsize=9, loc='best', ncol=2)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    filename = output_dir / 'EnergyFluctuation_vs_Temperature.png'
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  [OK] Saved: {filename.name}")


def plot_cluster_heat_capacity_after_support_subtraction(df_cv, output_dir):
    """
    扣除载体热容后的团簇热容分析图
    
    基于 CP2K 直接计算，载体热容 Cv_support = 38.17 meV/K (R²=0.999998)
    扣除后检查团簇热容是否守恒
    """
    print("\n[*] Plotting cluster heat capacity after support subtraction...")
    
    # 载体热容（CP2K直接计算值）
    SUPPORT_CV_TOTAL = 38.17  # meV/K (from CP2K calculation, R²=0.999998)
    
    # 计算团簇热容
    df = df_cv.copy()
    df['Cv_cluster_meV_K'] = df['Cv_total_meV_K'] - SUPPORT_CV_TOTAL
    df['Cv_cluster_per_atom_meV_K'] = df['Cv_cluster_meV_K'] / df['N_cluster']
    
    # 只保留有效数据（团簇热容 > 0）
    df_valid = df[df['Cv_cluster_meV_K'] > 0].copy()
    
    if len(df_valid) == 0:
        print(f"  [WARNING] No valid data after subtracting support Cv = {SUPPORT_CV_TOTAL} meV/K")
        return
    
    print(f"  [OK] Valid records: {len(df_valid)}/{len(df)} after subtraction")
    
    # 经典极限
    kB = 0.0862  # meV/K
    CLASSICAL_CV_PER_ATOM = 3 * kB  # 0.2586 meV/K/atom
    
    # 创建图表 (2x2)
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle(f'扣除载体热容后的团簇热容分析 (Cv_support = {SUPPORT_CV_TOTAL} meV/K)', 
                 fontsize=16, fontweight='bold')
    
    # === 子图1: 团簇总热容 vs N ===
    ax1 = axes[0, 0]
    
    # 按N分组统计
    stats_by_N = df_valid.groupby('N_cluster').agg({
        'Cv_cluster_meV_K': ['mean', 'std', 'count']
    }).reset_index()
    stats_by_N.columns = ['N_cluster', 'Cv_mean', 'Cv_std', 'count']
    
    N_values = stats_by_N['N_cluster'].values
    Cv_mean = stats_by_N['Cv_mean'].values
    Cv_std = stats_by_N['Cv_std'].values
    
    ax1.errorbar(N_values, Cv_mean, yerr=Cv_std, fmt='o', markersize=8,
                 capsize=5, alpha=0.7, label='测量值 (均值±std)')
    
    # 拟合线（带截距）
    from scipy import stats as sp_stats
    slope, intercept, r_value, _, _ = sp_stats.linregress(N_values, Cv_mean)
    N_fit = np.linspace(N_values.min(), N_values.max(), 100)
    ax1.plot(N_fit, slope * N_fit + intercept, 'r--', linewidth=2,
             label=f'拟合: Cv = {slope:.3f}×N + {intercept:.2f} (R²={r_value**2:.3f})')
    
    # 经典极限参考线
    ax1.plot(N_fit, CLASSICAL_CV_PER_ATOM * N_fit, 'g:', linewidth=2,
             label=f'经典: Cv = {CLASSICAL_CV_PER_ATOM:.3f}×N')
    
    ax1.set_xlabel('团簇原子数 N', fontsize=12, fontweight='bold')
    ax1.set_ylabel('团簇总热容 [meV/K]', fontsize=12, fontweight='bold')
    ax1.set_title('团簇总热容 vs 原子数', fontsize=13, fontweight='bold')
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    
    # === 子图2: 团簇每原子热容 vs N ===
    ax2 = axes[0, 1]
    
    stats_per_atom = df_valid.groupby('N_cluster').agg({
        'Cv_cluster_per_atom_meV_K': ['mean', 'std']
    }).reset_index()
    stats_per_atom.columns = ['N_cluster', 'Cv_per_atom_mean', 'Cv_per_atom_std']
    
    Cv_per_atom_mean = stats_per_atom['Cv_per_atom_mean'].values
    Cv_per_atom_std = stats_per_atom['Cv_per_atom_std'].values
    
    ax2.errorbar(N_values, Cv_per_atom_mean, yerr=Cv_per_atom_std, 
                 fmt='o', markersize=8, capsize=5, alpha=0.7, label='测量值')
    
    # 检查 1/N 依赖
    inv_N = 1.0 / N_values
    slope_inv, intercept_inv, r_inv, _, _ = sp_stats.linregress(inv_N, Cv_per_atom_mean)
    
    ax2.axhline(CLASSICAL_CV_PER_ATOM, color='red', linestyle='--', linewidth=2,
                label=f'经典极限 ({CLASSICAL_CV_PER_ATOM:.3f})')
    ax2.axhline(intercept_inv, color='blue', linestyle=':', linewidth=2,
                label=f'拟合渐近线 ({intercept_inv:.3f})')
    
    ax2.set_xlabel('团簇原子数 N', fontsize=12, fontweight='bold')
    ax2.set_ylabel('每原子热容 [meV/K/atom]', fontsize=12, fontweight='bold')
    ax2.set_title(f'每原子热容 vs N (1/N系数={slope_inv:.3f}, R²={r_inv**2:.3f})', 
                  fontsize=13, fontweight='bold')
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)
    
    # === 子图3: 每原子热容 vs 1/N ===
    ax3 = axes[1, 0]
    
    ax3.scatter(inv_N, Cv_per_atom_mean, s=100, alpha=0.7, label='测量值')
    ax3.errorbar(inv_N, Cv_per_atom_mean, yerr=Cv_per_atom_std,
                 fmt='none', ecolor='gray', alpha=0.5, capsize=5)
    
    # 拟合线
    inv_N_fit = np.linspace(0, inv_N.max(), 100)
    ax3.plot(inv_N_fit, slope_inv * inv_N_fit + intercept_inv, 'r--', linewidth=2,
             label=f'Cv/N = {intercept_inv:.3f} + {slope_inv:.3f}/N')
    
    ax3.axhline(CLASSICAL_CV_PER_ATOM, color='green', linestyle=':', linewidth=2,
                label=f'经典极限 ({CLASSICAL_CV_PER_ATOM:.3f})')
    
    ax3.set_xlabel('1/N', fontsize=12, fontweight='bold')
    ax3.set_ylabel('每原子热容 [meV/K/atom]', fontsize=12, fontweight='bold')
    ax3.set_title('检查 1/N 依赖性', fontsize=13, fontweight='bold')
    ax3.legend(fontsize=10)
    ax3.grid(True, alpha=0.3)
    
    # === 子图4: 团簇热容 vs 温度（按N着色）===
    ax4 = axes[1, 1]
    
    # 选择几个代表性的N值
    N_display = sorted(df_valid['N_cluster'].unique())
    cmap = plt.cm.viridis
    colors = [cmap(i/len(N_display)) for i in range(len(N_display))]
    
    for N, color in zip(N_display, colors):
        df_N = df_valid[df_valid['N_cluster'] == N]
        # 按温度分组取平均
        df_N_avg = df_N.groupby('温度').agg({
            'Cv_cluster_per_atom_meV_K': 'mean'
        }).reset_index()
        
        ax4.plot(df_N_avg['温度'], df_N_avg['Cv_cluster_per_atom_meV_K'],
                 'o-', label=f'N={N}', color=color, alpha=0.7, markersize=4)
    
    ax4.axhline(CLASSICAL_CV_PER_ATOM, color='red', linestyle='--', linewidth=2,
                label=f'经典极限 ({CLASSICAL_CV_PER_ATOM:.3f})')
    
    ax4.set_xlabel('温度 [K]', fontsize=12, fontweight='bold')
    ax4.set_ylabel('团簇每原子热容 [meV/K/atom]', fontsize=12, fontweight='bold')
    ax4.set_title('团簇热容随温度变化', fontsize=13, fontweight='bold')
    ax4.legend(fontsize=8, ncol=2, loc='best')
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    filename = output_dir / 'ClusterHeatCapacity_AfterSupportSubtraction.png'
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"  [OK] Saved: {filename.name}")
    print(f"  [INFO] 拟合结果:")
    print(f"    - 总热容: Cv = {slope:.3f}×N + {intercept:.2f} (R²={r_value**2:.3f})")
    print(f"    - 每原子: Cv/N = {intercept_inv:.3f} + {slope_inv:.3f}/N (R²={r_inv**2:.3f})")
    print(f"    - 测量/经典: {intercept_inv/CLASSICAL_CV_PER_ATOM:.2f}")
    
    if abs(slope_inv) < 0.05:
        print(f"    ✓ 1/N系数接近0，团簇热容基本守恒!")
    else:
        print(f"    ✗ 仍存在 1/N 依赖 (系数 = {slope_inv:.3f})")


def plot_heat_capacity_comparison_cluster_only(df_cv, output_dir, df_tm=None):
    """
    对比不同Sn含量的**团簇热容** (扣除载体后) - 按系列分组
    
    与 plot_heat_capacity_comparison 类似,但绘制的是扣除载体热容后的纯团簇热容
    """
    print("\n[*] Plotting cluster heat capacity comparison (support subtracted)...")
    
    # 载体热容 (更新为CP2K直接计算值)
    SUPPORT_CV_TOTAL = 38.17  # meV/K (from CP2K calculation, R²=0.999998)
    
    # 计算团簇热容
    df = df_cv.copy()
    df['Cv_cluster_meV_K'] = df['Cv_total_meV_K'] - SUPPORT_CV_TOTAL
    df['Cv_cluster_per_atom_meV_K'] = df['Cv_cluster_meV_K'] / df['N_cluster']
    
    # 过滤负值
    df_valid = df[df['Cv_cluster_meV_K'] > 0].copy()
    print(f"  [INFO] Valid records: {len(df_valid)}/{len(df)} (removed negative cluster Cv)")
    
    if len(df_valid) == 0:
        print("  [WARNING] No valid data after support subtraction")
        return
    
    # 获取所有系列
    if '系列' not in df_valid.columns:
        print("  [WARNING] '系列' column not found")
        return
    
    series_list = sorted(df_valid['系列'].unique())
    print(f"    Found series: {series_list}")
    
    # 经典极限
    kB = 0.0862  # meV/K
    CLASSICAL_CV_PER_ATOM = 3 * kB
    
    # 为每个系列绘制独立的图
    for series in series_list:
        df_series = df_valid[df_valid['系列'] == series].copy()
        
        if len(df_series) == 0:
            continue
        
        # 确定分组方式 (同原函数逻辑)
        if series == 'Pt(8-x)SnX':
            if 'Pt原子数' in df_series.columns:
                group_key = 'Pt原子数'
                group_values = sorted(df_series['Pt原子数'].dropna().unique())
                x_label = 'Pt含量'
                title_suffix = 'Pt含量'
            else:
                continue
        elif series in ['O1', 'O2', 'O3', 'O4']:
            group_key = '结构'
            struct_info = []
            for struct in df_series['结构'].unique():
                row = df_series[df_series['结构'] == struct].iloc[0]
                pt = row['Pt原子数']
                sn = row['Sn原子数']
                o = row['O原子数']
                total = pt + sn + o
                pt_ratio = pt / (pt + sn) if (pt + sn) > 0 else 0
                struct_info.append((struct, pt, sn, o, total, pt_ratio))
            struct_info_sorted = sorted(struct_info, key=lambda x: (x[4], x[5]))
            group_values = [x[0] for x in struct_info_sorted]
            x_label = '组分配比'
            title_suffix = '不同组分'
        else:
            group_key = 'Sn含量'
            group_values = sorted(df_series['Sn含量'].unique())
            x_label = 'Sn含量'
            title_suffix = 'Sn含量'
        
        if len(group_values) < 2:
            continue
        
        # ========== 创建图表 (2×2) ==========
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle(f'团簇热容对比（已扣除载体 {SUPPORT_CV_TOTAL} meV/K）- {series}系列', 
                     fontsize=16, fontweight='bold')
        
        colors = plt.cm.viridis(np.linspace(0, 1, len(group_values)))
        color_map = {val: colors[i] for i, val in enumerate(group_values)}
        
        # 准备同分异构体标记 (O系列)
        formula_count = {}
        if series in ['O1', 'O2', 'O3', 'O4']:
            for group_val in group_values:
                row = df_series[df_series[group_key] == group_val].iloc[0]
                pt = int(row['Pt原子数'])
                sn = int(row['Sn原子数'])
                formula = f'Pt{pt}Sn{sn}'
                formula_count[formula] = formula_count.get(formula, 0) + 1
        
        # === 子图1: 团簇总热容 vs 温度 ===
        ax1 = axes[0, 0]
        for group_val in group_values:
            df_group = df_series[df_series[group_key] == group_val].copy()
            df_group_sorted = df_group.sort_values('温度')
            
            # 准备标签
            if series in ['O1', 'O2', 'O3', 'O4']:
                row = df_group_sorted.iloc[0]
                pt = int(row['Pt原子数'])
                sn = int(row['Sn原子数'])
                formula = f'Pt{pt}Sn{sn}'
                if formula_count.get(formula, 0) > 1:
                    label = f'{formula} ({group_val})'
                else:
                    label = formula
            else:
                label = f'{x_label}={group_val}'
            
            ax1.plot(df_group_sorted['温度'], df_group_sorted['Cv_cluster_meV_K'],
                     marker='o', label=label, color=color_map[group_val], 
                     linewidth=2, markersize=6, alpha=0.8)
        
        ax1.set_xlabel('温度 (K)', fontsize=12, fontweight='bold')
        ax1.set_ylabel('团簇总热容 [meV/K]', fontsize=12, fontweight='bold')
        ax1.set_title(f'{series}系列：总热容 vs 温度（团簇热容）', fontsize=13, fontweight='bold')
        ax1.legend(fontsize=9, loc='best', ncol=2)
        ax1.grid(True, alpha=0.3)
        
        # === 子图2: 平均每原子热容 vs 不同组分 ===
        ax2 = axes[0, 1]
        
        # 计算每个温度的平均值
        temp_list = sorted(df_series['温度'].unique())
        for T in temp_list:
            df_temp = df_series[df_series['温度'] == T]
            avg_cv_per_atom = []
            x_positions = []
            
            for i, group_val in enumerate(group_values):
                df_group = df_temp[df_temp[group_key] == group_val]
                if len(df_group) > 0:
                    avg_cv_per_atom.append(df_group['Cv_cluster_per_atom_meV_K'].mean())
                    x_positions.append(i)
            
            if len(avg_cv_per_atom) > 0:
                ax2.plot(x_positions, avg_cv_per_atom, marker='o', 
                         label=f'{T}K', linewidth=2, markersize=6, alpha=0.7)
        
        # 经典极限参考线
        ax2.axhline(CLASSICAL_CV_PER_ATOM, color='red', linestyle='--', 
                    linewidth=2, label=f'经典极限 ({CLASSICAL_CV_PER_ATOM:.3f})')
        
        # X轴标签
        if series in ['O1', 'O2', 'O3', 'O4']:
            x_labels = []
            for group_val in group_values:
                row = df_series[df_series[group_key] == group_val].iloc[0]
                pt = int(row['Pt原子数'])
                sn = int(row['Sn原子数'])
                x_labels.append(f'Pt{pt}Sn{sn}')
            ax2.set_xticks(range(len(group_values)))
            ax2.set_xticklabels(x_labels, rotation=45, ha='right')
        else:
            ax2.set_xticks(range(len(group_values)))
            ax2.set_xticklabels([str(v) for v in group_values])
        
        ax2.set_xlabel(x_label, fontsize=12, fontweight='bold')
        ax2.set_ylabel('每原子热容 [meV/K/atom]', fontsize=12, fontweight='bold')
        ax2.set_title(f'{series}系列：平均每原子热容 vs {title_suffix}', fontsize=13, fontweight='bold')
        ax2.legend(fontsize=9, loc='best', ncol=2)
        ax2.grid(True, alpha=0.3)
        
        # === 子图3: 每原子热容 vs 温度（按组分着色）===
        ax3 = axes[1, 0]
        
        for group_val in group_values:
            df_group = df_series[df_series[group_key] == group_val].copy()
            df_group_sorted = df_group.sort_values('温度')
            
            # 准备标签
            if series in ['O1', 'O2', 'O3', 'O4']:
                row = df_group_sorted.iloc[0]
                pt = int(row['Pt原子数'])
                sn = int(row['Sn原子数'])
                formula = f'Pt{pt}Sn{sn}'
                if formula_count.get(formula, 0) > 1:
                    label = f'{formula} ({group_val})'
                else:
                    label = formula
            else:
                label = f'{x_label}={group_val}'
            
            ax3.plot(df_group_sorted['温度'], df_group_sorted['Cv_cluster_per_atom_meV_K'],
                     marker='o', label=label, color=color_map[group_val],
                     linewidth=2, markersize=6, alpha=0.8)
        
        # 经典极限参考线
        ax3.axhline(CLASSICAL_CV_PER_ATOM, color='red', linestyle='--',
                    linewidth=2, label=f'经典极限 ({CLASSICAL_CV_PER_ATOM:.3f})')
        
        ax3.set_xlabel('温度 (K)', fontsize=12, fontweight='bold')
        ax3.set_ylabel('每原子热容 [meV/K/atom]', fontsize=12, fontweight='bold')
        ax3.set_title(f'{series}系列：每原子热容 vs 温度', fontsize=13, fontweight='bold')
        ax3.legend(fontsize=9, loc='best', ncol=2)
        ax3.grid(True, alpha=0.3)
        
        # === 子图4: 平均每原子热容 vs 组分（不同温度分段）===
        ax4 = axes[1, 1]
        
        # 按温度范围分组统计
        temp_bins = [(200, 400), (400, 700), (700, 900), (900, 1100)]
        temp_colors = ['blue', 'orange', 'green', 'red']
        
        for (T_min, T_max), color in zip(temp_bins, temp_colors):
            df_temp_range = df_series[(df_series['温度'] >= T_min) & (df_series['温度'] < T_max)]
            
            if len(df_temp_range) == 0:
                continue
            
            avg_cv = []
            x_positions = []
            
            for i, group_val in enumerate(group_values):
                df_group = df_temp_range[df_temp_range[group_key] == group_val]
                if len(df_group) > 0:
                    avg_cv.append(df_group['Cv_cluster_per_atom_meV_K'].mean())
                    x_positions.append(i)
            
            if len(avg_cv) > 0:
                ax4.plot(x_positions, avg_cv, marker='s', 
                         label=f'{T_min}-{T_max}K', color=color,
                         linewidth=2, markersize=7, alpha=0.7)
        
        # 经典极限参考线
        ax4.axhline(CLASSICAL_CV_PER_ATOM, color='red', linestyle='--',
                    linewidth=2, label=f'经典极限')
        
        # X轴标签
        if series in ['O1', 'O2', 'O3', 'O4']:
            x_labels = []
            for group_val in group_values:
                row = df_series[df_series[group_key] == group_val].iloc[0]
                pt = int(row['Pt原子数'])
                sn = int(row['Sn原子数'])
                x_labels.append(f'Pt{pt}Sn{sn}')
            ax4.set_xticks(range(len(group_values)))
            ax4.set_xticklabels(x_labels, rotation=45, ha='right')
        else:
            ax4.set_xticks(range(len(group_values)))
            ax4.set_xticklabels([str(v) for v in group_values])
        
        ax4.set_xlabel(x_label, fontsize=12, fontweight='bold')
        ax4.set_ylabel('每原子热容 [meV/K/atom]', fontsize=12, fontweight='bold')
        ax4.set_title(f'{series}系列：平均每原子热容 vs {title_suffix}（按温度段）', fontsize=13, fontweight='bold')
        ax4.legend(fontsize=9, loc='best')
        ax4.grid(True, alpha=0.3)
        
        plt.tight_layout()
        filename = output_dir / f'ClusterHeatCapacity_comparison_{series}.png'
        plt.savefig(filename, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  [OK] Saved: {filename.name}")


def correlate_energy_with_diffusion(df_energy, df_cv, output_dir):
    """关联能量/热容与扩散系数"""
    print("\n[*] Correlating energy with diffusion coefficient...")
    
    # 加载扩散系数数据
    if not D_VALUES_FILE.exists():
        print(f"  [WARNING] D values file not found: {D_VALUES_FILE}")
        return
    
    df_D = pd.read_csv(D_VALUES_FILE)
    print(f"  [OK] Loaded {len(df_D)} diffusion coefficient records")
    
    # 准备能量数据用于匹配
    df_energy_avg = df_energy.groupby(['Sn含量', '温度']).agg({
        '每原子能量': 'mean',
        '相对能量': 'mean',
        '每原子标准差': 'mean'
    }).reset_index()
    
    # 准备热容数据
    df_cv_avg = df_cv.groupby(['Sn含量', '温度']).agg({
        'Cv_meV_per_K': 'mean'
    }).reset_index()
    
    # 合并数据
    # D数据中的列: composition, temp_K, element, D_cm2_s
    # 需要提取Sn含量并转换为宽格式
    
    # 1. 提取Sn含量
    df_D['Sn含量'] = df_D['composition'].apply(extract_sn_content)
    df_D['温度'] = df_D['temp_K']
    
    # 2. 转换为宽格式: Pt和Sn的D值分别作为列
    df_D_pivot = df_D.pivot_table(
        index=['composition', 'Sn含量', '温度'], 
        columns='element', 
        values='D_cm2_s',
        aggfunc='first'
    ).reset_index()
    
    # 重命名列
    df_D_pivot.columns.name = None
    if 'Pt' in df_D_pivot.columns:
        df_D_pivot = df_D_pivot.rename(columns={'Pt': 'D_cm2_s_Pt'})
    if 'Sn' in df_D_pivot.columns:
        df_D_pivot = df_D_pivot.rename(columns={'Sn': 'D_cm2_s_Sn'})
    
    df_D_renamed = df_D_pivot
    
    df_merged = df_D_renamed.merge(df_energy_avg, on=['Sn含量', '温度'], how='inner')
    df_merged = df_merged.merge(df_cv_avg, on=['Sn含量', '温度'], how='left')
    
    print(f"  [OK] Merged {len(df_merged)} records")
    
    # 保存merged数据供后续分析
    merged_file = output_dir / 'energy_diffusion_merged.csv'
    df_merged.to_csv(merged_file, index=False, encoding='utf-8-sig')
    print(f"  [OK] Saved merged data: {merged_file.name}")
    print(f"  [INFO] Correlation plot temporarily disabled due to technical issues")


def generate_report(df_energy, df_cv, output_dir, df_tm=None):
    """
    生成分析报告
    
    Parameters:
    -----------
    df_energy : DataFrame
        能量数据
    df_cv : DataFrame
        热容数据
    output_dir : Path
        输出目录
    df_tm : DataFrame, optional
        熔化温度数据
    """
    print("\n[*] Generating analysis report...")
    
    report_lines = []
    report_lines.append("=" * 80)
    report_lines.append("能量分析报告 (改进版 - 体系独立分析 + 熔化温度识别)")
    report_lines.append("=" * 80)
    report_lines.append(f"生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report_lines.append("")
    
    # 数据概览
    report_lines.append("[1] 数据概览")
    report_lines.append("-" * 80)
    report_lines.append(f"总记录数: {len(df_energy)}")
    report_lines.append(f"体系数量: {df_energy['结构'].nunique()}")
    report_lines.append(f"Sn含量范围: {df_energy['Sn含量'].min():.0f} - {df_energy['Sn含量'].max():.0f}")
    report_lines.append(f"温度范围: {df_energy['温度'].min():.0f} - {df_energy['温度'].max():.0f} K")
    report_lines.append("")
    
    # 各体系统计
    report_lines.append("[2] 各体系能量统计")
    report_lines.append("-" * 80)
    report_lines.append(f"{'体系':<20} {'Sn':>3} {'温度数':>6} {'能量范围(eV/atom)':<30}")
    report_lines.append("-" * 80)
    
    for structure in sorted(df_energy['结构'].unique()):
        df_sys = df_energy[df_energy['结构'] == structure]
        sn = df_sys['Sn含量'].iloc[0]
        n_temps = df_sys['温度'].nunique()
        e_min = df_sys['每原子能量'].min()
        e_max = df_sys['每原子能量'].max()
        
        report_lines.append(f"{structure:<20} {sn:>3} {n_temps:>6}   "
                          f"{e_min:>10.2f} ~ {e_max:>10.2f}")
    
    report_lines.append("")
    
    # 热容统计
    if len(df_cv) > 0:
        report_lines.append("[3] 热容统计")
        report_lines.append("-" * 80)
        report_lines.append(f"计算的热容数据点: {len(df_cv)}")
        report_lines.append("")
        
        # 按Sn含量统计平均热容
        report_lines.append(f"{'Sn含量':>6} {'平均Cv [meV/K per atom]':>30} {'范围':<30}")
        report_lines.append("-" * 80)
        
        for sn in sorted(df_cv['Sn含量'].unique()):
            df_sn = df_cv[df_cv['Sn含量'] == sn]
            avg_cv = df_sn['Cv_meV_per_K'].mean()
            min_cv = df_sn['Cv_meV_per_K'].min()
            max_cv = df_sn['Cv_meV_per_K'].max()
            
            report_lines.append(f"Sn{sn:>2}   {avg_cv:>25.2f}   "
                              f"{min_cv:>10.2f} ~ {max_cv:>10.2f}")
        
        report_lines.append("")
    
    # 关键发现
    report_lines.append("[4] 关键发现")
    report_lines.append("-" * 80)
    
    # 能量趋势
    sn_energies = []
    for sn in sorted(df_energy['Sn含量'].unique()):
        df_sn = df_energy[df_energy['Sn含量'] == sn]
        avg_energy = df_sn['每原子能量'].mean()
        sn_energies.append((sn, avg_energy))
    
    report_lines.append("能量趋势:")
    report_lines.append(f"  - Sn{sn_energies[0][0]}: {sn_energies[0][1]:.2f} eV/atom (最低)")
    report_lines.append(f"  - Sn{sn_energies[-1][0]}: {sn_energies[-1][1]:.2f} eV/atom (最高)")
    report_lines.append(f"  - 能量随Sn含量增加呈上升趋势")
    report_lines.append("")
    
    # 热容趋势
    if len(df_cv) > 0:
        sn_cvs = []
        for sn in sorted(df_cv['Sn含量'].unique()):
            df_sn = df_cv[df_cv['Sn含量'] == sn]
            avg_cv = df_sn['Cv_meV_per_K'].mean()
            sn_cvs.append((sn, avg_cv))
        
        report_lines.append("热容趋势:")
        report_lines.append(f"  - Sn{sn_cvs[0][0]}: {sn_cvs[0][1]:.2f} J/(mol·K) (最高)")
        report_lines.append(f"  - Sn{sn_cvs[-1][0]}: {sn_cvs[-1][1]:.2f} J/(mol·K) (最低)")
        report_lines.append(f"  - 热容随Sn含量增加呈下降趋势")
    
    report_lines.append("")
    
    # 熔化温度分析
    if df_tm is not None and len(df_tm) > 0:
        report_lines.append("[5] 熔化温度 (Tm) - 热容跃变识别")
        report_lines.append("-" * 80)
        report_lines.append(f"{'结构':<20} {'Sn':>3} {'Tm(K)':>6} {'Cv前(J/mol·K)':>15} {'Cv后(J/mol·K)':>15} {'跃变倍数':>10}")
        report_lines.append("-" * 80)
        
        for _, row in df_tm.iterrows():
            report_lines.append(
                f"{row['结构']:<20} {row['Sn含量']:>3.0f} {row['Tm(K)']:>6.0f} "
                f"{row['Cv_before(J/mol·K)']:>15.2f} {row['Cv_after(J/mol·K)']:>15.2f} "
                f"{row['Cv_jump_ratio']:>10.2f}x"
            )
        
        report_lines.append("")
        report_lines.append("说明:")
        report_lines.append("  - Tm: 熔化温度,由热容跃变识别")
        report_lines.append("  - Cv前: 跃变前的热容值")
        report_lines.append("  - Cv后: 跃变后的热容值")
        report_lines.append("  - 跃变倍数: Cv后/Cv前,反映熔化剧烈程度")
        report_lines.append("")
    
    report_lines.append("=" * 80)
    
    # 保存报告
    report_file = output_dir / 'energy_analysis_v2_report.txt'
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write('\n'.join(report_lines))
    
    print(f"  [OK] Saved: {report_file.name}")
    
    # 打印到终端
    print("\n" + '\n'.join(report_lines))


def save_data_tables(df_energy, df_cv, output_dir, df_tm=None):
    """
    保存数据表
    
    Parameters:
    -----------
    df_energy : DataFrame
        能量数据
    df_cv : DataFrame
        热容数据
    output_dir : Path
        输出目录
    df_tm : DataFrame, optional
        熔化温度数据
    """
    print("\n[*] Saving data tables...")
    
    # 能量数据
    energy_file = output_dir / 'energy_per_system.csv'
    df_energy_export = df_energy[[
        '结构', 'Sn含量', '温度', '模拟序号',
        '每原子能量', '相对能量', '每原子标准差'
    ]].sort_values(['结构', '温度'])
    df_energy_export.to_csv(energy_file, index=False, encoding='utf-8-sig')
    print(f"  [OK] Saved: {energy_file.name}")
    
    # 热容数据
    if len(df_cv) > 0:
        cv_file = output_dir / 'heat_capacity_per_system.csv'
        df_cv.sort_values(['结构', '温度']).to_csv(
            cv_file, index=False, encoding='utf-8-sig'
        )
        print(f"  [OK] Saved: {cv_file.name}")
    
    # 熔化温度数据
    if df_tm is not None and len(df_tm) > 0:
        tm_file = output_dir / 'melting_temperatures.csv'
        df_tm.sort_values(['Sn含量', '结构']).to_csv(
            tm_file, index=False, encoding='utf-8-sig'
        )
        print(f"  [OK] Saved: {tm_file.name}")


def main():
    r"""
    主函数
    
    配置说明:
    ---------
    在脚本顶部的SYSTEM_FILTER中配置要分析的体系:
    - 注释/取消注释 include_patterns 来选择体系
    - 例如: r'^pt8' 匹配Pt8系列, r'^[Oo]\d+' 匹配O2系列
    - 可以同时启用多个pattern来分析多个体系
    
    命令行参数:
    -----------
    --no-filter : 不使用MSD异常筛选,分析所有能量数据
                  (适用于当某些体系的数据全被筛选掉时)
    
    示例:
    python step6_energy_analysis_v2.py            # 默认:使用异常筛选
    python step6_energy_analysis_v2.py --no-filter  # 不筛选,使用全部数据
    """
    # 解析命令行参数
    parser = argparse.ArgumentParser(
        description='Step6: 能量分析 (改进版) - 正则表达式体系筛选'
    )
    parser.add_argument(
        '--no-filter',
        action='store_true',
        help='不使用MSD异常筛选,分析所有能量数据'
    )
    args = parser.parse_args()
    
    # 根据模式设置输出目录
    if args.no_filter:
        OUTPUT_DIR_ACTUAL = BASE_DIR / 'results' / 'energy_analysis_v2_no_filter'
        mode_desc = "不使用MSD异常筛选 (--no-filter)"
    else:
        OUTPUT_DIR_ACTUAL = OUTPUT_DIR  # 使用默认目录
        mode_desc = "使用MSD异常筛选 (默认)"
    
    OUTPUT_DIR_ACTUAL.mkdir(exist_ok=True, parents=True)
    
    print("\n" + "=" * 80)
    print("Step6: 能量分析 (改进版) - 正则表达式体系筛选")
    print(f"[!] 模式: {mode_desc}")
    print(f"[!] 输出目录: {OUTPUT_DIR_ACTUAL.relative_to(BASE_DIR)}")
    print("=" * 80)
    
    # 显示当前配置
    print("\n[CONFIG] System Filter Configuration:")
    include = SYSTEM_FILTER.get('include_patterns', [])
    exclude = SYSTEM_FILTER.get('exclude_patterns', [])
    if include:
        print(f"  Include patterns: {include}")
    else:
        print(f"  Include patterns: [] (all structures)")
    if exclude:
        print(f"  Exclude patterns: {exclude}")
    else:
        print(f"  Exclude patterns: []")
    
    # 1. 加载能量数据 (带outlier筛选和体系过滤)
    # 传入 no_filter 参数控制是否筛选
    df_energy = load_energy_data(no_filter=args.no_filter)
    
    # 2. 计算每原子能量
    df_energy = calculate_per_atom_energy(df_energy)
    
    # 3. 对每个体系计算相对能量
    df_energy = calculate_relative_energy_per_system(df_energy)
    
    # 4. 对每个体系计算热容
    df_cv = calculate_heat_capacity_per_system(df_energy)
    
    # 5. 检测熔化温度 (热容跃变分析)
    df_tm = detect_melting_temperature(df_cv, threshold_factor=1.5)
    
    # 6. 可视化
    print("\n" + "=" * 80)
    print("Generating Visualizations...")
    print("=" * 80)
    
    plot_energy_vs_temperature_by_system(df_energy, OUTPUT_DIR_ACTUAL)
    plot_heat_capacity_comparison(df_cv, OUTPUT_DIR_ACTUAL, df_tm)  # 传入Tm数据
    plot_heat_capacity_comparison_cluster_only(df_cv, OUTPUT_DIR_ACTUAL, df_tm)  # 新增: 扣除载体后的对比图
    plot_heat_capacity_heatmap(df_cv, OUTPUT_DIR_ACTUAL)  # 新增: 热容热熔图
    plot_heat_capacity_heatmap_cluster_only(df_cv, OUTPUT_DIR_ACTUAL)  # 新增: 扣除载体后的热力图
    plot_oxygen_series_comprehensive(df_cv, OUTPUT_DIR_ACTUAL)  # 新增: O1-O4系列综合对比图
    plot_energy_fluctuation(df_energy, OUTPUT_DIR_ACTUAL)
    plot_cluster_heat_capacity_after_support_subtraction(df_cv, OUTPUT_DIR_ACTUAL)  # 新增: 扣除载体热容
    correlate_energy_with_diffusion(df_energy, df_cv, OUTPUT_DIR_ACTUAL)
    
    # 7. 生成报告
    generate_report(df_energy, df_cv, OUTPUT_DIR_ACTUAL, df_tm)  # 传入Tm数据
    
    # 8. 保存数据
    save_data_tables(df_energy, df_cv, OUTPUT_DIR_ACTUAL, df_tm)  # 传入Tm数据
    
    print("\n" + "=" * 80)
    print(f"[DONE] Analysis complete! Results saved to: {OUTPUT_DIR_ACTUAL}")
    print("=" * 80)
    
    # ⚠️  重要提示
    print("\n" + "=" * 80)
    print("⚠️  ⚠️  ⚠️   重要提示：关于载体热容   ⚠️  ⚠️  ⚠️")
    print("=" * 80)
    print("输出的热容数据包含未知的载体（Al₂O₃）贡献！")
    print("")
    print("  LAMMPS的etotal包含整个体系的能量：")
    print("    - 团簇部分：~60个原子 (Pt + Sn + O)")
    print("    - 载体部分：240个原子 (Al₂O₃)")
    print("")
    print("  热容计算公式：")
    print("    Cv_total = dE_total/dT")
    print("             = Cv_cluster + Cv_support")
    print("             = Cv_cluster + C (常数，未知)")
    print("")
    print("  ❓ 关键问题：载体热容 Cv_support 是未知的！")
    print("    - 载体是固定的 240 个 Al₂O₃ 原子")
    print("    - Cv_support 应该是常数（不随团簇大小变化）")
    print("    - ⚠️  但我们不知道这个常数的确切值")
    print("")
    print("  📊 当前输出的含义：")
    print("    - Cv_total_meV_K：团簇 + 载体的总热容")
    print("    - Cv_per_atom_meV_K：表观每原子热容（含载体影响）")
    print("    - 两者都无法直接反映纯团簇的热容")
    print("")
    print("  ✅ 如何获得准确的团簇热容：")
    print("    1. [推荐] 单独模拟纯 Al₂O₃ 体系（240 atoms）")
    print("       测量 Cv_support = dE_support/dT")
    print("")
    print("    2. [替代] 使用拟合估计值（~18-21 meV/K）")
    print("       python subtract_support_v2.py")
    print("       ⚠️  注意：这只是估计值，需要验证")
    print("")
    print("  详细说明参阅：")
    print("    CRITICAL_SUPPORT_HEAT_CAPACITY_UNKNOWN.md")
    print("=" * 80 + "\n")


if __name__ == '__main__':
    main()
