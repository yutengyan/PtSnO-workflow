#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
Step 7.2: Per-Atom 动力学分析 - 林德曼指数 + 扩散系数
================================================================================

作者: GitHub Copilot
日期: 2025-12-03
版本: v1.0

功能概述
========
结合每个原子的林德曼指数(δ)和扩散系数(D)，分析：
1. 哪些原子在动？哪个元素先动？
2. 高振动原子是否也是高扩散原子？
3. 表面vs核心的动力学差异
4. 温度依赖性分析

核心分析
========
- δ vs D 散点图 (按元素着色)
- 元素级统计: mean(δ_Pt), mean(D_Pt), std(δ_Pt)
- 高扩散原子的元素分布
- 低温vs高温的原子活跃程度变化
- δ-D相关性分析

输入数据
========
1. Per-atom Lindemann: data/lindemann/per-atoms/per_atom_master_run_*.csv
   列: 目录, 结构, 温度(K), atom_id, element, lindemann_index, ...
   
2. Per-atom D: data/gmx_msd/per-atom/collected_gmx_per_atom_msd/per_atom_diffusion_coefficients_*.csv
   列: 结构, 温度(K), atom_id, 元素, D(1e-5 cm²/s), ...

输出
====
results/per_atom_dynamics/
├── merged_per_atom_data.csv           # 合并后的数据
├── element_statistics.csv             # 元素级统计
├── correlation_analysis.csv           # δ-D相关性
├── delta_vs_D_scatter_*.png           # δ vs D 散点图
├── element_comparison_*.png           # 元素对比图
├── temperature_evolution_*.png        # 温度演化图
└── analysis_report.txt                # 分析报告

使用示例
========
# 基础分析
python step7_2_per_atom_dynamics_analysis.py

# 指定结构
python step7_2_per_atom_dynamics_analysis.py --structure pt8sn6

# 指定温度范围
python step7_2_per_atom_dynamics_analysis.py --temp-range 300 900

# 交互模式
python step7_2_per_atom_dynamics_analysis.py --interactive
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns
from pathlib import Path
import re
from scipy import stats
import argparse
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# 配置部分
# ============================================================================

# 字体设置 - Times New Roman for publication
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman', 'DejaVu Serif']
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['font.size'] = 12

# 基础路径
BASE_DIR = Path(__file__).parent
DATA_DIR = BASE_DIR / 'data'
OUTPUT_DIR = BASE_DIR / 'results' / 'per_atom_dynamics'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# 数据文件路径
# 默认读取 sup86-50k（50K步长，19温度点）；旧数据: DATA_DIR / 'lindemann' / 'per-atoms'
LINDEMANN_DIR = DATA_DIR / 'lindemann' / 'per-atoms' / 'sup86-50k'
MSD_DIR = DATA_DIR / 'gmx_msd' / 'per-atom' / 'collected_gmx_per_atom_msd'

# 元素颜色方案
ELEMENT_COLORS = {
    'Pt': '#1f77b4',   # 蓝色
    'Sn': '#ff7f0e',   # 橙色
    'O': '#2ca02c',    # 绿色
    'Al': '#d62728',   # 红色
}

# 元素标记
ELEMENT_MARKERS = {
    'Pt': 'o',
    'Sn': 's',
    'O': '^',
    'Al': 'd',
}

# 阈值
HIGH_LINDEMANN_THRESHOLD = 0.10  # δ > 0.10 为高振动
HIGH_D_THRESHOLD = 0.01          # D > 0.01 (1e-5 cm²/s) 为高扩散

# GMX ensemble MSD 数据路径
GMX_MSD_DIR = DATA_DIR / 'gmx_msd' / 'unwrap' / 'gmx_msd_results_20251118_152614'

# Structure-level Lindemann 数据路径 (与Step 7.8.5一致)
STRUCTURE_LEVEL_LINDEMANN_FILE = BASE_DIR / 'results' / 'step7_8_2_alloy_series' / 'structure_level_lindemann.csv'

# 体系分类颜色 - 与 step6_3_1 一致
SERIES_COLORS = {
    'air': '#1f77b4',        # 蓝色 - 气相合金
    'sum8': '#ff7f0e',       # 橙色 - Pt+Sn=8
    'pt8snx': '#2ca02c',     # 绿色 - Pt=8
    'pt6snx': '#d62728',     # 红色 - Pt=6
    'o1': '#9467bd',         # 紫色 - O=1
    'o2': '#8c564b',         # 棕色 - O=2
    'o3': '#e377c2',         # 粉色 - O=3
    'o4': '#7f7f7f',         # 灰色 - O=4
    'other': '#bcbd22',      # 黄绿 - 其他
}

# 特征温度
FEATURE_TEMPS = [300, 900]  # K


# ============================================================================
# 体系分类函数 - 与 step6_3_1_cv_scatter_plot.py 保持一致
# ============================================================================

def classify_system(structure_name, path=None):
    """
    分类体系类型 - 使用与 step6_3_1 一致的分类逻辑
    
    分组:
    - air: 气相合金 (路径含 air/airmd)
    - sum8: Pt+Sn=8 的无氧体系 (如 Pt6Sn2, Pt4Sn4)
    - pt8snx: Pt=8 的无氧体系 (如 Pt8Sn2, Pt8Sn4)
    - pt6snx: Pt=6 的无氧体系 (如 Pt6Sn2, Pt6Sn4, Pt6Sn8)
    - o1, o2, o3, o4: 按氧原子数分类
    
    Returns:
        dict with keys: series, has_oxide, n_pt, n_sn, n_o, composition
    """
    if pd.isna(structure_name):
        structure_name = ''
    
    name_lower = str(structure_name).lower()
    path_str = str(path).lower() if path else ''
    
    # 提取原子数
    pt_match = re.search(r'pt(\d+)', name_lower)
    sn_match = re.search(r'sn(\d+)', name_lower)
    o_match = re.search(r'o(\d+)', name_lower)
    
    n_pt = int(pt_match.group(1)) if pt_match else 0
    n_sn = int(sn_match.group(1)) if sn_match else 0
    n_o = int(o_match.group(1)) if o_match else 0
    
    # 判断是否含氧
    has_oxide = n_o > 0
    
    # 判断是否为气相 (air)
    is_air = 'air' in name_lower or 'air' in path_str or 'airmd' in path_str
    
    # 确定系列
    series = 'other'
    
    if has_oxide:
        # 按氧原子数分类
        if n_o == 1:
            series = 'o1'
        elif n_o == 2:
            series = 'o2'
        elif n_o == 3:
            series = 'o3'
        elif n_o >= 4:
            series = 'o4'
    elif is_air:
        series = 'air'
    else:
        # 无氧体系，按 Pt+Sn 总数或 Pt 数量分类
        total = n_pt + n_sn
        if total == 8:
            series = 'sum8'
        elif n_pt == 8:
            series = 'pt8snx'
        elif n_pt == 6:
            series = 'pt6snx'
    
    return {
        'series': series,
        'has_oxide': has_oxide,
        'n_pt': n_pt,
        'n_sn': n_sn,
        'n_o': n_o,
        'composition': f'Pt{n_pt}Sn{n_sn}' + (f'O{n_o}' if n_o > 0 else '')
    }


def add_system_classification(df):
    """为数据添加体系分类"""
    print("\n[*] 添加体系分类...")
    
    # 获取path列（如果存在）
    path_col = 'path' if 'path' in df.columns else None
    
    classifications = df.apply(
        lambda row: classify_system(row['structure'], 
                                   row[path_col] if path_col else None), 
        axis=1
    )
    
    df['series'] = [c['series'] for c in classifications]
    df['has_oxide'] = [c['has_oxide'] for c in classifications]
    df['n_pt'] = [c['n_pt'] for c in classifications]
    df['n_sn'] = [c['n_sn'] for c in classifications]
    df['n_o'] = [c['n_o'] for c in classifications]
    df['composition'] = [c['composition'] for c in classifications]
    
    # 打印分类统计
    print(f"  体系分类统计:")
    series_order = ['air', 'sum8', 'pt8snx', 'pt6snx', 'o1', 'o2', 'o3', 'o4', 'other']
    for series in series_order:
        if series in df['series'].values:
            count = len(df[df['series'] == series])
            pct = count / len(df) * 100
            n_struct = df[df['series'] == series]['structure'].nunique()
            print(f"    {series}: {count} ({pct:.1f}%), {n_struct} structures")
    
    print(f"  含氧体系: {df['has_oxide'].sum()} ({df['has_oxide'].mean()*100:.1f}%)")
    
    return df


# ============================================================================
# 数据加载函数
# ============================================================================

def load_lindemann_data(data_dir):
    """加载 per-atom Lindemann 数据"""
    print("\n[1/4] 加载 Per-Atom Lindemann 数据...")
    
    files = sorted(Path(data_dir).glob('per_atom_master_*.csv'))
    # 同时搜索子目录（兼容旧目录结构和 sup86-50k 等子目录）
    if not files:
        files = sorted(Path(data_dir).glob('*/per_atom_master_*.csv'))
    if not files:
        print(f"  [WARNING] 未找到 Lindemann 数据文件: {data_dir}")
        return None
    
    dfs = []
    for f in files:
        try:
            df = pd.read_csv(f, encoding='utf-8')
        except:
            df = pd.read_csv(f, encoding='gbk')
        dfs.append(df)
        print(f"  - {f.name}: {len(df)} records")
    
    df_all = pd.concat(dfs, ignore_index=True)
    
    # 标准化列名
    col_map = {
        '目录': 'path',
        '结构': 'structure', 
        '温度(K)': 'temp',
        'atom_id': 'atom_id',
        'element': 'element',
        'lindemann_index': 'delta'
    }
    
    # 查找实际列名并重命名
    for old, new in col_map.items():
        for col in df_all.columns:
            if old in col or old == col:
                df_all = df_all.rename(columns={col: new})
                break
    
    # 从路径中提取结构名（如果结构列为空）
    if 'structure' in df_all.columns:
        if df_all['structure'].isna().all() or (df_all['structure'] == '').all():
            df_all['structure'] = df_all['path'].apply(extract_structure_from_path)
    
    print(f"  [OK] 合并: {len(df_all)} atom records, {df_all['structure'].nunique()} structures")
    
    return df_all


def load_msd_data(data_dir):
    """加载 per-atom MSD/D 数据"""
    print("\n[2/4] 加载 Per-Atom MSD/D 数据...")
    
    files = sorted(Path(data_dir).glob('per_atom_diffusion_coefficients_*.csv'))
    if not files:
        print(f"  [WARNING] 未找到 MSD 数据文件: {data_dir}")
        return None
    
    dfs = []
    for f in files:
        try:
            df = pd.read_csv(f, encoding='utf-8')
        except:
            df = pd.read_csv(f, encoding='gbk')
        dfs.append(df)
        print(f"  - {f.name}: {len(df)} records")
    
    df_all = pd.concat(dfs, ignore_index=True)
    
    # 标准化列名
    col_map = {
        '结构': 'structure',
        '温度(K)': 'temp',
        'atom_id': 'atom_id',
        '元素': 'element',
        'D(1e-5 cm': 'D',  # D值列可能有不完整名称
        '完整目录路径': 'full_path',  # 用于精确匹配run
    }
    
    for old, new in col_map.items():
        for col in df_all.columns:
            if old in col:
                df_all = df_all.rename(columns={col: new})
                break
    
    # 确保有D列
    if 'D' not in df_all.columns:
        for col in df_all.columns:
            if 'D(' in col or col.startswith('D'):
                df_all = df_all.rename(columns={col: 'D'})
                break
    
    # 确保有full_path列
    if 'full_path' not in df_all.columns:
        # 尝试从其他列获取
        if '完整目录路径' in df_all.columns:
            df_all['full_path'] = df_all['完整目录路径']
        elif 'XVG 文件路径' in df_all.columns:
            # 从XVG路径提取目录
            df_all['full_path'] = df_all['XVG 文件路径'].apply(lambda x: str(x).rsplit('/', 1)[0] if pd.notna(x) else None)
    
    print(f"  [OK] 合并: {len(df_all)} atom records, {df_all['structure'].nunique()} structures")
    
    return df_all


def extract_structure_from_path(path):
    """从路径中提取结构名"""
    if pd.isna(path):
        return 'Unknown'
    parts = str(path).replace('\\', '/').split('/')
    for part in reversed(parts):
        part_lower = part.lower()
        if 'pt' in part_lower and ('sn' in part_lower or part_lower.startswith('pt')):
            # 移除后缀如 .r1.gpu0
            name = re.sub(r'\.[rR]\d+\.gpu\d+.*$', '', part)
            return name
    return 'Unknown'


def merge_data(df_lindemann, df_msd, structure_filter=None, temp_range=None):
    """合并 Lindemann 和 MSD 数据"""
    print("\n[3/4] 合并数据...")
    
    if df_lindemann is None or df_msd is None:
        print("  [ERROR] 数据缺失，无法合并")
        return None
    
    # 标准化结构名
    df_lindemann['structure_std'] = df_lindemann['structure'].str.lower().str.replace('-', '').str.replace('_', '')
    df_msd['structure_std'] = df_msd['structure'].str.lower().str.replace('-', '').str.replace('_', '')
    
    # 确保类型一致
    df_lindemann['temp'] = pd.to_numeric(df_lindemann['temp'], errors='coerce')
    df_msd['temp'] = pd.to_numeric(df_msd['temp'], errors='coerce')
    df_lindemann['atom_id'] = pd.to_numeric(df_lindemann['atom_id'], errors='coerce')
    df_msd['atom_id'] = pd.to_numeric(df_msd['atom_id'], errors='coerce')
    
    # 重要: 按完整路径+temp+atom_id合并,避免多个run交叉匹配
    # 检查是否有path列(lindemann)和full_path列(msd)
    if 'path' in df_lindemann.columns and 'full_path' in df_msd.columns:
        print("  [INFO] 使用完整路径匹配 (避免多run交叉)")
        df_merged = pd.merge(
            df_lindemann[['structure', 'structure_std', 'path', 'temp', 'atom_id', 'element', 'delta']],
            df_msd[['structure_std', 'full_path', 'temp', 'atom_id', 'D']],
            left_on=['path', 'temp', 'atom_id'],
            right_on=['full_path', 'temp', 'atom_id'],
            how='inner'
        )
    else:
        # 降级方案: 按 structure_std + temp + atom_id 合并 (可能产生重复!)
        print("  [WARNING] 无完整路径信息,使用structure+temp+atom_id合并 (可能重复)")
        df_merged = pd.merge(
            df_lindemann[['structure', 'structure_std', 'temp', 'atom_id', 'element', 'delta']],
            df_msd[['structure_std', 'temp', 'atom_id', 'D']],
            on=['structure_std', 'temp', 'atom_id'],
            how='inner'
        )
    
    print(f"  [INFO] 合并结果: {len(df_merged)} matched records")
    
    if len(df_merged) == 0:
        print("  [WARNING] 无匹配数据，尝试宽松合并...")
        # 尝试只按 temp + atom_id 合并
        df_merged = pd.merge(
            df_lindemann[['structure', 'temp', 'atom_id', 'element', 'delta']],
            df_msd[['temp', 'atom_id', 'D']],
            on=['temp', 'atom_id'],
            how='inner'
        )
        print(f"  [INFO] 宽松合并结果: {len(df_merged)} records")
    
    # 应用筛选
    if structure_filter and len(df_merged) > 0:
        df_merged = df_merged[df_merged['structure'].str.contains(structure_filter, case=False, na=False)]
        print(f"  [FILTER] 结构筛选 '{structure_filter}': {len(df_merged)} records")
    
    if temp_range and len(df_merged) > 0:
        df_merged = df_merged[(df_merged['temp'] >= temp_range[0]) & (df_merged['temp'] <= temp_range[1])]
        print(f"  [FILTER] 温度范围 {temp_range}: {len(df_merged)} records")
    
    # 添加分类标签
    if len(df_merged) > 0:
        df_merged['high_delta'] = df_merged['delta'] > HIGH_LINDEMANN_THRESHOLD
        df_merged['high_D'] = df_merged['D'] > HIGH_D_THRESHOLD
        df_merged['mobility_class'] = df_merged.apply(classify_mobility, axis=1)
        
        # 添加体系分类
        df_merged = add_system_classification(df_merged)
    
    return df_merged


def classify_mobility(row):
    """分类原子活跃程度"""
    if row['high_delta'] and row['high_D']:
        return 'Active'      # 高振动 + 高扩散
    elif row['high_delta'] and not row['high_D']:
        return 'Vibrating'   # 高振动 + 低扩散 (局域振动)
    elif not row['high_delta'] and row['high_D']:
        return 'Diffusing'   # 低振动 + 高扩散 (跳跃式)
    else:
        return 'Stable'      # 低振动 + 低扩散


# ============================================================================
# 统计分析函数
# ============================================================================

def calculate_ensemble_D(df_group):
    """
    计算系综平均扩散系数
    
    ⚠️ 已弃用此方法! 
    现在使用 calculate_ensemble_D_from_gmx_msd() 从GMX ensemble MSD拟合
    
    旧方法: 对每个元素,先对所有run的原子D值求平均(得到每个run的平均D),
         然后对所有run的平均D再求平均,得到系综平均D
         
    这比直接对所有原子的D取平均更合理,因为考虑了run间的统计独立性
    
    参数:
        df_group: 包含structure, temp, element, D列的DataFrame
        
    返回:
        D_ensemble: 系综平均扩散系数
        D_std_runs: run间的标准差
        n_runs: run数量
    """
    if 'full_path' not in df_group.columns or df_group['full_path'].isna().all():
        # 如果没有run信息,退化为简单平均
        return df_group['D'].mean(), df_group['D'].std(), 1
    
    # 按run分组,计算每个run的平均D
    run_D_means = df_group.groupby('full_path')['D'].mean()
    
    # 计算系综平均
    D_ensemble = run_D_means.mean()
    D_std_runs = run_D_means.std() if len(run_D_means) > 1 else 0
    n_runs = len(run_D_means)
    
    return D_ensemble, D_std_runs, n_runs


def read_gmx_msd_xvg(filepath):
    """读取GMX MSD .xvg文件"""
    time_data = []
    msd_data = []
    
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#') or line.startswith('@'):
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        t = float(parts[0])
                        msd_nm2 = float(parts[1])
                        msd_a2 = msd_nm2 * 100  # nm² -> Ų
                        time_data.append(t)
                        msd_data.append(msd_a2)
                    except ValueError:
                        continue
    except Exception as e:
        return None, None
    
    if len(time_data) == 0:
        return None, None
    
    return np.array(time_data), np.array(msd_data)


def calculate_ensemble_D_from_gmx_msd(structure, temp, element, 
                                       fit_range=(20, 140), 
                                       window_size=120, step=5):
    """
    从GMX ensemble MSD计算扩散系数 (与step7_5一致的方法)
    
    方法:
    1. 加载该structure/temp/element的所有runs的MSD
    2. 计算平均MSD曲线
    3. 使用滑动窗口拟合 (window=120ps, step=5ps)
    4. 返回平均D值和CV
    
    参数:
        structure: 如 'pt8sn6'
        temp: 温度(K), 如 900
        element: 元素, 如 'Pt', 'Sn'
        fit_range: 拟合起始范围 (ps)
        window_size: 滑动窗口宽度 (ps)
        step: 滑动步长 (ps)
        
    返回:
        D_mean: 平均扩散系数 (×10⁻⁵ cm²/s)
        D_std: 标准差 (从滑动窗口)
        CV: 变异系数 (%)
        r2_mean: 平均R²
        n_runs: run数量
    """
    from scipy import stats
    
    # 构建温度字符串
    temp_str = f'{temp}K'
    
    # 扫描文件
    msd_list = []
    for xvg_file in GMX_MSD_DIR.rglob("*_msd_*.xvg"):
        try:
            parts = xvg_file.parts
            filename = xvg_file.stem
            
            # 检查元素
            if '_msd_' not in filename:
                continue
            file_element = filename.split('_msd_')[-1]
            if file_element != element:
                continue
            
            # 检查温度和结构
            temperature = None
            composition = None
            for i in range(len(parts)-1, 0, -1):
                if parts[i].endswith('K'):
                    temperature = parts[i]
                    composition = parts[i-1]
                    break
            
            if temperature != temp_str:
                continue
            if not composition.lower().startswith(structure.lower()):
                continue
            
            # 读取数据
            time, msd = read_gmx_msd_xvg(xvg_file)
            if time is not None:
                msd_list.append((time, msd))
        except:
            continue
    
    if not msd_list:
        return None, None, None, None, 0
    
    # 计算平均MSD
    min_len = min(len(msd) for _, msd in msd_list)
    msd_aligned = np.array([msd[:min_len] for _, msd in msd_list])
    time_common = msd_list[0][0][:min_len]
    msd_mean = np.mean(msd_aligned, axis=0)
    
    # 滑动窗口拟合 (与step7_5一致)
    actual_max = time_common[-1]
    slide_start = fit_range[0]
    slide_end = actual_max - window_size
    
    if slide_end < slide_start + 2 * step:
        # 回退到单窗口
        mask = (time_common >= fit_range[0]) & (time_common <= min(fit_range[1], actual_max))
        if mask.sum() < 10:
            return None, None, None, None, len(msd_list)
        
        time_fit = time_common[mask]
        msd_fit = msd_mean[mask]
        
        slope, intercept, r_value, p_value, std_err = stats.linregress(time_fit, msd_fit)
        D_A2_per_ps = slope / 6.0
        D_cm2_per_s = D_A2_per_ps * 1e-4
        
        return D_cm2_per_s * 1e5, 0.0, 0.0, r_value**2, len(msd_list)
    
    # 滑动窗口拟合
    D_values = []
    r2_values = []
    
    for t_start in np.arange(slide_start, slide_end + step, step):
        t_end = t_start + window_size
        mask = (time_common >= t_start) & (time_common <= t_end)
        
        if mask.sum() < 10:
            continue
        
        time_fit = time_common[mask]
        msd_fit = msd_mean[mask]
        
        slope, intercept, r_value, p_value, std_err = stats.linregress(time_fit, msd_fit)
        D_A2_per_ps = slope / 6.0
        D_cm2_per_s = D_A2_per_ps * 1e-4
        
        D_values.append(D_cm2_per_s * 1e5)  # 转换为10^-5 cm²/s
        r2_values.append(r_value ** 2)
    
    if not D_values:
        return None, None, None, None, len(msd_list)
    
    D_values = np.array(D_values)
    D_mean = np.mean(D_values)
    D_std = np.std(D_values)
    CV = D_std / D_mean * 100 if D_mean != 0 else 0
    r2_mean = np.mean(r2_values)
    
    return D_mean, D_std, CV, r2_mean, len(msd_list)


def calculate_ensemble_D_with_filtering(structure, temp, element, 
                                         fit_range=(20, 140), 
                                         window_size=120, step=5,
                                         filter_method='r2',
                                         r2_threshold=0.90,
                                         outlier_sigma=2.0):
    """
    从GMX ensemble MSD计算扩散系数 - 带质量筛选
    
    筛选方法:
    1. 'r2': 筛除R²低于阈值的run
    2. 'outlier': 筛除D值偏离均值超过N个标准差的run
    3. 'both': 同时使用R²和离群值筛选
    4. None: 不筛选(默认)
    
    参数:
        structure: 结构名称
        temp: 温度(K)
        element: 元素
        fit_range: 拟合范围 (ps)
        window_size: 滑动窗口大小 (ps)
        step: 滑动步长 (ps)
        filter_method: 筛选方法 ('r2', 'outlier', 'both', None)
        r2_threshold: R²阈值 (默认0.90)
        outlier_sigma: 离群值判断标准差倍数 (默认2.0)
    
    返回:
        D_mean, D_std, CV, r2_mean, n_runs, n_filtered, filter_info
    """
    from scipy import stats
    
    temp_str = f'{temp}K'
    
    # 1. 收集所有MSD文件
    msd_list = []
    msd_files = []  # 保存文件路径
    
    for xvg_file in GMX_MSD_DIR.rglob("*_msd_*.xvg"):
        try:
            parts = xvg_file.parts
            filename = xvg_file.stem
            
            if '_msd_' not in filename:
                continue
            file_element = filename.split('_msd_')[-1]
            if file_element != element:
                continue
            
            temperature = None
            composition = None
            for i in range(len(parts)-1, 0, -1):
                if parts[i].endswith('K'):
                    temperature = parts[i]
                    composition = parts[i-1]
                    break
            
            if temperature != temp_str:
                continue
            if not composition.lower().startswith(structure.lower()):
                continue
            
            time, msd = read_gmx_msd_xvg(xvg_file)
            if time is not None:
                msd_list.append((time, msd))
                msd_files.append(str(xvg_file))
        except:
            continue
    
    if not msd_list:
        return None, None, None, None, 0, 0, {}
    
    n_total = len(msd_list)
    
    # 2. 如果需要筛选,先对每个run拟合计算D和R²
    if filter_method is not None:
        run_D_values = []
        run_r2_values = []
        
        for time, msd in msd_list:
            # 简单拟合整个fit_range
            mask = (time >= fit_range[0]) & (time <= fit_range[1])
            if mask.sum() < 10:
                run_D_values.append(np.nan)
                run_r2_values.append(0.0)
                continue
            
            time_fit = time[mask]
            msd_fit = msd[mask]
            
            slope, intercept, r_value, p_value, std_err = stats.linregress(time_fit, msd_fit)
            D_A2_per_ps = slope / 6.0
            D_cm2_per_s = D_A2_per_ps * 1e-4
            
            run_D_values.append(D_cm2_per_s * 1e5)
            run_r2_values.append(r_value ** 2)
        
        run_D_values = np.array(run_D_values)
        run_r2_values = np.array(run_r2_values)
        
        # 3. 应用筛选
        keep_mask = np.ones(n_total, dtype=bool)
        filter_reasons = []
        
        if filter_method in ['r2', 'both']:
            # 筛除低R²
            r2_mask = run_r2_values >= r2_threshold
            keep_mask &= r2_mask
            n_r2_filtered = (~r2_mask).sum()
            if n_r2_filtered > 0:
                filter_reasons.append(f"R²<{r2_threshold}: {n_r2_filtered} runs")
        
        if filter_method in ['outlier', 'both']:
            # 筛除离群值 (使用有效D值)
            valid_D = run_D_values[~np.isnan(run_D_values)]
            if len(valid_D) > 2:
                D_mean_temp = np.mean(valid_D)
                D_std_temp = np.std(valid_D)
                
                outlier_mask = np.ones(n_total, dtype=bool)
                for i in range(n_total):
                    if not np.isnan(run_D_values[i]):
                        z_score = abs(run_D_values[i] - D_mean_temp) / D_std_temp
                        if z_score > outlier_sigma:
                            outlier_mask[i] = False
                
                keep_mask &= outlier_mask
                n_outlier_filtered = (~outlier_mask).sum()
                if n_outlier_filtered > 0:
                    filter_reasons.append(f"D outlier (>{outlier_sigma}σ): {n_outlier_filtered} runs")
        
        # 4. 应用筛选
        n_filtered = (~keep_mask).sum()
        
        if n_filtered > 0:
            msd_list = [msd_list[i] for i in range(n_total) if keep_mask[i]]
            msd_files = [msd_files[i] for i in range(n_total) if keep_mask[i]]
            
            filter_info = {
                'n_total': n_total,
                'n_kept': len(msd_list),
                'n_filtered': n_filtered,
                'filter_reasons': filter_reasons,
                'kept_files': msd_files,
                'run_D_before': run_D_values,
                'run_r2_before': run_r2_values
            }
        else:
            filter_info = {
                'n_total': n_total,
                'n_kept': n_total,
                'n_filtered': 0,
                'filter_reasons': ['No runs filtered']
            }
    else:
        n_filtered = 0
        filter_info = {}
    
    if not msd_list:
        return None, None, None, None, n_total, n_filtered, filter_info
    
    # 5. 计算筛选后的ensemble D (使用滑动窗口)
    min_len = min(len(msd) for _, msd in msd_list)
    msd_aligned = np.array([msd[:min_len] for _, msd in msd_list])
    time_common = msd_list[0][0][:min_len]
    msd_mean = np.mean(msd_aligned, axis=0)
    
    # 滑动窗口拟合
    actual_max = time_common[-1]
    slide_start = fit_range[0]
    slide_end = actual_max - window_size
    
    if slide_end < slide_start + 2 * step:
        # 回退到单窗口
        mask = (time_common >= fit_range[0]) & (time_common <= min(fit_range[1], actual_max))
        if mask.sum() < 10:
            return None, None, None, None, len(msd_list), n_filtered, filter_info
        
        time_fit = time_common[mask]
        msd_fit = msd_mean[mask]
        
        slope, intercept, r_value, p_value, std_err = stats.linregress(time_fit, msd_fit)
        D_A2_per_ps = slope / 6.0
        D_cm2_per_s = D_A2_per_ps * 1e-4
        
        return D_cm2_per_s * 1e5, 0.0, 0.0, r_value**2, len(msd_list), n_filtered, filter_info
    
    # 滑动窗口拟合
    D_values = []
    r2_values = []
    
    for t_start in np.arange(slide_start, slide_end + step, step):
        t_end = t_start + window_size
        mask = (time_common >= t_start) & (time_common <= t_end)
        
        if mask.sum() < 10:
            continue
        
        time_fit = time_common[mask]
        msd_fit = msd_mean[mask]
        
        slope, intercept, r_value, p_value, std_err = stats.linregress(time_fit, msd_fit)
        D_A2_per_ps = slope / 6.0
        D_cm2_per_s = D_A2_per_ps * 1e-4
        
        D_values.append(D_cm2_per_s * 1e5)
        r2_values.append(r_value ** 2)
    
    if not D_values:
        return None, None, None, None, len(msd_list), n_filtered, filter_info
    
    D_values = np.array(D_values)
    D_mean = np.mean(D_values)
    D_std = np.std(D_values)
    CV = D_std / D_mean * 100 if D_mean != 0 else 0
    r2_mean = np.mean(r2_values)
    
    return D_mean, D_std, CV, r2_mean, len(msd_list), n_filtered, filter_info


def calculate_element_statistics(df, use_ensemble_D=True, use_gmx_msd=True):
    """
    按元素计算统计量 (包含系列信息)
    
    参数:
        df: 合并后的数据
        use_ensemble_D: 是否使用系综平均计算D (推荐True)
        use_gmx_msd: 是否从GMX ensemble MSD计算D (推荐True, 与step7_5一致)
    """
    print(f"\n[4/4] 计算元素统计 (系综平均D: {use_ensemble_D}, 使用GMX MSD: {use_gmx_msd})...")
    
    stats_list = []
    
    # 获取系列信息
    has_series = 'series' in df.columns
    
    for (structure, temp), group in df.groupby(['structure', 'temp']):
        series = group['series'].iloc[0] if has_series else 'unknown'
        has_oxide = group['has_oxide'].iloc[0] if 'has_oxide' in group.columns else False
        
        for element in group['element'].unique():
            elem_data = group[group['element'] == element]
            
            # 计算扩散系数
            if use_gmx_msd:
                # 新方法: 从GMX ensemble MSD拟合 (与step7_5一致)
                D_mean, D_std, CV, r2_mean, n_runs = calculate_ensemble_D_from_gmx_msd(
                    structure, temp, element
                )
                
                if D_mean is None:
                    print(f"  [WARNING] 无法从GMX MSD计算 {structure} {temp}K {element}, 使用per-atom平均")
                    # 回退到旧方法
                    if use_ensemble_D:
                        D_ensemble, D_std_runs, n_runs = calculate_ensemble_D(elem_data)
                        D_mean = D_ensemble
                        D_std = D_std_runs
                        CV = None
                        r2_mean = None
                    else:
                        D_mean = elem_data['D'].mean()
                        D_std = elem_data['D'].std()
                        n_runs = elem_data['full_path'].nunique() if 'full_path' in elem_data.columns else 1
                        CV = None
                        r2_mean = None
            elif use_ensemble_D:
                # 旧方法1: 系综平均 (per-atom平均)
                D_ensemble, D_std_runs, n_runs = calculate_ensemble_D(elem_data)
                D_mean = D_ensemble
                D_std = D_std_runs
                CV = None
                r2_mean = None
            else:
                # 旧方法2: 直接对所有原子取平均
                D_mean = elem_data['D'].mean()
                D_std = elem_data['D'].std()
                n_runs = elem_data['full_path'].nunique() if 'full_path' in elem_data.columns else 1
                CV = None
                r2_mean = None
            
            stats_dict = {
                'structure': structure,
                'temp': temp,
                'series': series,
                'has_oxide': has_oxide,
                'element': element,
                'n_atoms': len(elem_data),
                'n_runs': n_runs,  # run数量
                # Lindemann 统计
                'delta_mean': elem_data['delta'].mean(),
                'delta_std': elem_data['delta'].std(),
                'delta_min': elem_data['delta'].min(),
                'delta_max': elem_data['delta'].max(),
                'delta_median': elem_data['delta'].median(),
                # D 统计
                'D_mean': D_mean,
                'D_std': D_std,
                'D_min': elem_data['D'].min(),
                'D_max': elem_data['D'].max(),
                'D_median': elem_data['D'].median(),
                # 高活性比例
                'high_delta_ratio': (elem_data['delta'] > HIGH_LINDEMANN_THRESHOLD).mean(),
                'high_D_ratio': (elem_data['D'] > HIGH_D_THRESHOLD).mean(),
                'active_ratio': ((elem_data['delta'] > HIGH_LINDEMANN_THRESHOLD) & 
                                (elem_data['D'] > HIGH_D_THRESHOLD)).mean(),
            }
            
            # 添加GMX MSD特有的指标
            if use_gmx_msd and CV is not None:
                stats_dict['CV'] = CV
                stats_dict['r2_mean'] = r2_mean
            
            stats_list.append(stats_dict)
    
    df_stats = pd.DataFrame(stats_list)
    print(f"  [OK] 生成 {len(df_stats)} 条元素统计")
    
    return df_stats


def load_structure_level_lindemann():
    """
    读取structure-level的Lindemann指数 (与Step 7.8.5一致)
    
    这个方法更符合系综平均的物理意义:
    1. 每次run内: 先对元素内所有原子求平均 (保持该run的特征)
    2. 多次run间: 再对多次run求平均 (真正的系综平均)
    
    Returns:
        pd.DataFrame with columns: structure, temp, delta_Pt, delta_Sn, delta_Pt_std, delta_Sn_std
    """
    if not STRUCTURE_LEVEL_LINDEMANN_FILE.exists():
        print(f"  [WARNING] Structure-level文件不存在: {STRUCTURE_LEVEL_LINDEMANN_FILE}")
        return None
    
    df = pd.read_csv(STRUCTURE_LEVEL_LINDEMANN_FILE, encoding='utf-8-sig')
    
    # 只保留需要的列
    df_lindemann = df[['structure', 'temperature', 'delta_Pt', 'delta_Sn', 
                       'delta_Pt_std', 'delta_Sn_std']].copy()
    df_lindemann = df_lindemann.rename(columns={'temperature': 'temp'})
    
    print(f"  [OK] 读取structure-level Lindemann: {len(df_lindemann)} 条记录")
    return df_lindemann


def calculate_series_statistics(df, use_ensemble_D=True):
    """
    按系列(气相/负载/含氧)计算统计量
    
    参数:
        df: 合并后的数据
        use_ensemble_D: 是否使用系综平均计算D (推荐True)
    """
    print(f"\n[*] 计算分系列统计 (系综平均D: {use_ensemble_D})...")
    
    if 'series' not in df.columns:
        print("  [WARNING] 无系列分类信息")
        return None
    
    stats_list = []
    
    for series in df['series'].unique():
        series_data = df[df['series'] == series]
        
        for element in series_data['element'].unique():
            elem_data = series_data[series_data['element'] == element]
            
            # 按温度计算
            for temp in sorted(elem_data['temp'].unique()):
                temp_data = elem_data[elem_data['temp'] == temp]
                
                # 计算系综平均D
                if use_ensemble_D:
                    D_ensemble, D_std_runs, n_runs = calculate_ensemble_D(temp_data)
                    D_mean = D_ensemble
                    D_std = D_std_runs
                else:
                    D_mean = temp_data['D'].mean()
                    D_std = temp_data['D'].std()
                    n_runs = temp_data['full_path'].nunique() if 'full_path' in temp_data.columns else 1
                
                stats_list.append({
                    'series': series,
                    'element': element,
                    'temp': temp,
                    'n_atoms': len(temp_data),
                    'n_structures': temp_data['structure'].nunique(),
                    'n_runs': n_runs,  # 新增
                    # δ 统计
                    'delta_mean': temp_data['delta'].mean(),
                    'delta_std': temp_data['delta'].std(),
                    'delta_median': temp_data['delta'].median(),
                    # D 统计 (系综平均)
                    'D_mean': D_mean,
                    'D_std': D_std,
                    'D_median': temp_data['D'].median(),
                    # 活性比例
                    'high_delta_ratio': (temp_data['delta'] > HIGH_LINDEMANN_THRESHOLD).mean(),
                    'high_D_ratio': (temp_data['D'] > HIGH_D_THRESHOLD).mean(),
                    'active_ratio': ((temp_data['delta'] > HIGH_LINDEMANN_THRESHOLD) & 
                                    (temp_data['D'] > HIGH_D_THRESHOLD)).mean(),
                })
    
    df_stats = pd.DataFrame(stats_list)
    print(f"  [OK] 生成 {len(df_stats)} 条分系列统计")
    
    return df_stats


def calculate_correlation(df):
    """计算 δ-D 相关性"""
    print("\n[*] 计算 δ-D 相关性...")
    
    corr_list = []
    
    for (structure, temp), group in df.groupby(['structure', 'temp']):
        if len(group) < 5:
            continue
        
        # 整体相关性
        r_all, p_all = stats.pearsonr(group['delta'], group['D'])
        
        corr_list.append({
            'structure': structure,
            'temp': temp,
            'element': 'All',
            'n_atoms': len(group),
            'pearson_r': r_all,
            'p_value': p_all,
        })
        
        # 按元素分
        for element in group['element'].unique():
            elem_data = group[group['element'] == element]
            if len(elem_data) < 3:
                continue
            
            try:
                r, p = stats.pearsonr(elem_data['delta'], elem_data['D'])
                corr_list.append({
                    'structure': structure,
                    'temp': temp,
                    'element': element,
                    'n_atoms': len(elem_data),
                    'pearson_r': r,
                    'p_value': p,
                })
            except:
                pass
    
    df_corr = pd.DataFrame(corr_list)
    print(f"  [OK] 计算 {len(df_corr)} 条相关性")
    
    return df_corr


# ============================================================================
# 可视化函数
# ============================================================================

def plot_series_comparison(df_series_stats, output_dir):
    """绘制分系列对比图 - 核心分析图"""
    print("\n[*] 绘制分系列对比图 (气相 vs 负载 vs 含氧)...")
    
    if df_series_stats is None or len(df_series_stats) == 0:
        print("  [WARNING] 无分系列统计数据")
        return
    
    # 创建大图: 3行 x 3列
    fig, axes = plt.subplots(3, 3, figsize=(18, 15))
    
    series_list = df_series_stats['series'].unique()
    series_colors = {s: SERIES_COLORS.get(s, 'gray') for s in series_list}
    
    # ===== 第一行: Pt 元素对比 =====
    for col_idx, metric in enumerate(['delta_mean', 'D_mean', 'high_delta_ratio']):
        ax = axes[0, col_idx]
        
        for series in series_list:
            series_data = df_series_stats[(df_series_stats['series'] == series) & 
                                          (df_series_stats['element'] == 'Pt')]
            if len(series_data) == 0:
                continue
            series_data = series_data.sort_values('temp')
            
            color = series_colors[series]
            
            if metric == 'high_delta_ratio':
                y_data = series_data[metric] * 100
                ylabel = 'High δ Ratio (%)'
            else:
                y_data = series_data[metric]
                ylabel = 'Mean δ' if metric == 'delta_mean' else 'Mean D (10⁻⁵ cm²/s)'
            
            ax.plot(series_data['temp'], y_data, 'o-', color=color, 
                   linewidth=2, markersize=6, label=series)
        
        ax.set_xlabel('Temperature (K)', fontsize=11)
        ax.set_ylabel(ylabel, fontsize=11)
        ax.set_title(f'Pt: {ylabel}', fontsize=12, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        if metric == 'D_mean':
            ax.set_yscale('log')
    
    # ===== 第二行: Sn 元素对比 =====
    for col_idx, metric in enumerate(['delta_mean', 'D_mean', 'high_delta_ratio']):
        ax = axes[1, col_idx]
        
        for series in series_list:
            series_data = df_series_stats[(df_series_stats['series'] == series) & 
                                          (df_series_stats['element'] == 'Sn')]
            if len(series_data) == 0:
                continue
            series_data = series_data.sort_values('temp')
            
            color = series_colors[series]
            
            if metric == 'high_delta_ratio':
                y_data = series_data[metric] * 100
                ylabel = 'High δ Ratio (%)'
            else:
                y_data = series_data[metric]
                ylabel = 'Mean δ' if metric == 'delta_mean' else 'Mean D (10⁻⁵ cm²/s)'
            
            ax.plot(series_data['temp'], y_data, 's-', color=color, 
                   linewidth=2, markersize=6, label=series)
        
        ax.set_xlabel('Temperature (K)', fontsize=11)
        ax.set_ylabel(ylabel, fontsize=11)
        ax.set_title(f'Sn: {ylabel}', fontsize=12, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        if metric == 'D_mean':
            ax.set_yscale('log')
    
    # ===== 第三行: O 元素 + Pt-Sn差异 + 活跃比例 =====
    # O 元素
    ax = axes[2, 0]
    o_data = df_series_stats[df_series_stats['element'] == 'O']
    if len(o_data) > 0:
        for series in o_data['series'].unique():
            series_data = o_data[o_data['series'] == series].sort_values('temp')
            color = series_colors.get(series, 'green')
            ax.plot(series_data['temp'], series_data['delta_mean'], '^-', 
                   color=color, linewidth=2, markersize=6, label=f'{series}')
        ax.set_xlabel('Temperature (K)', fontsize=11)
        ax.set_ylabel('Mean δ', fontsize=11)
        ax.set_title('O: Mean Lindemann Index', fontsize=12, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
    else:
        ax.text(0.5, 0.5, 'No O data', ha='center', va='center', fontsize=14)
        ax.set_title('O: Mean Lindemann Index', fontsize=12, fontweight='bold')
    
    # Pt vs Sn 差异
    ax = axes[2, 1]
    for series in series_list:
        pt_data = df_series_stats[(df_series_stats['series'] == series) & 
                                  (df_series_stats['element'] == 'Pt')].set_index('temp')
        sn_data = df_series_stats[(df_series_stats['series'] == series) & 
                                  (df_series_stats['element'] == 'Sn')].set_index('temp')
        
        common_temps = pt_data.index.intersection(sn_data.index)
        if len(common_temps) == 0:
            continue
        
        diff = pt_data.loc[common_temps, 'delta_mean'] - sn_data.loc[common_temps, 'delta_mean']
        color = series_colors[series]
        ax.plot(common_temps, diff, 'o-', color=color, linewidth=2, markersize=6, label=series)
    
    ax.axhline(0, color='black', linestyle='--', alpha=0.5)
    ax.set_xlabel('Temperature (K)', fontsize=11)
    ax.set_ylabel('δ(Pt) - δ(Sn)', fontsize=11)
    ax.set_title('Pt-Sn Lindemann Difference', fontsize=12, fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 活跃比例对比
    ax = axes[2, 2]
    for series in series_list:
        series_data = df_series_stats[df_series_stats['series'] == series]
        temp_active = series_data.groupby('temp')['active_ratio'].mean() * 100
        color = series_colors[series]
        ax.plot(temp_active.index, temp_active.values, 'D-', 
               color=color, linewidth=2, markersize=6, label=series)
    
    ax.set_xlabel('Temperature (K)', fontsize=11)
    ax.set_ylabel('Active Atom Ratio (%)', fontsize=11)
    ax.set_title('Overall Active Ratio by Series', fontsize=12, fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.suptitle('Per-Atom Dynamics: Series Comparison\n(Air vs Supported vs Oxide)', 
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    
    output_path = output_dir / 'series_comparison.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"  [SAVED] {output_path}")
    
    return fig


def plot_oxide_vs_nonoxide(df, output_dir):
    """绘制含氧 vs 无氧体系对比图"""
    print("\n[*] 绘制含氧 vs 无氧对比图...")
    
    if 'has_oxide' not in df.columns:
        print("  [WARNING] 无含氧分类信息")
        return
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # (1) δ 分布箱线图
    ax1 = axes[0, 0]
    data_to_plot = []
    labels = []
    for has_ox, label in [(False, 'Non-oxide'), (True, 'Oxide')]:
        subset = df[df['has_oxide'] == has_ox]
        for elem in ['Pt', 'Sn']:
            elem_data = subset[subset['element'] == elem]['delta']
            if len(elem_data) > 0:
                data_to_plot.append(elem_data)
                labels.append(f'{label}\n{elem}')
    
    if data_to_plot:
        bp = ax1.boxplot(data_to_plot, labels=labels, patch_artist=True)
        colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728'][:len(data_to_plot)]
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.6)
    
    ax1.axhline(HIGH_LINDEMANN_THRESHOLD, color='red', linestyle='--', alpha=0.5)
    ax1.set_ylabel('Lindemann Index (δ)', fontsize=12)
    ax1.set_title('(a) δ Distribution: Oxide vs Non-oxide', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    
    # (2) D 分布箱线图
    ax2 = axes[0, 1]
    data_to_plot = []
    labels = []
    for has_ox, label in [(False, 'Non-oxide'), (True, 'Oxide')]:
        subset = df[df['has_oxide'] == has_ox]
        for elem in ['Pt', 'Sn']:
            elem_data = subset[subset['element'] == elem]['D']
            if len(elem_data) > 0:
                data_to_plot.append(elem_data)
                labels.append(f'{label}\n{elem}')
    
    if data_to_plot:
        bp = ax2.boxplot(data_to_plot, labels=labels, patch_artist=True)
        colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728'][:len(data_to_plot)]
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.6)
    
    ax2.axhline(HIGH_D_THRESHOLD, color='blue', linestyle='--', alpha=0.5)
    ax2.set_ylabel('D (10⁻⁵ cm²/s)', fontsize=12)
    ax2.set_title('(b) D Distribution: Oxide vs Non-oxide', fontsize=12, fontweight='bold')
    ax2.set_yscale('log')
    ax2.grid(True, alpha=0.3)
    
    # (3) 活跃比例 vs 温度
    ax3 = axes[1, 0]
    for has_ox, label, color in [(False, 'Non-oxide', '#1f77b4'), (True, 'Oxide', '#2ca02c')]:
        subset = df[df['has_oxide'] == has_ox]
        active_by_temp = subset.groupby('temp').apply(
            lambda x: ((x['delta'] > HIGH_LINDEMANN_THRESHOLD) & (x['D'] > HIGH_D_THRESHOLD)).mean() * 100
        )
        ax3.plot(active_by_temp.index, active_by_temp.values, 'o-', 
                color=color, linewidth=2, markersize=8, label=label)
    
    ax3.set_xlabel('Temperature (K)', fontsize=12)
    ax3.set_ylabel('Active Atom Ratio (%)', fontsize=12)
    ax3.set_title('(c) Active Ratio vs Temperature', fontsize=12, fontweight='bold')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # (4) O 原子动力学
    ax4 = axes[1, 1]
    o_data = df[(df['has_oxide'] == True) & (df['element'] == 'O')]
    if len(o_data) > 0:
        o_by_temp = o_data.groupby('temp').agg({'delta': 'mean', 'high_delta': 'mean'})
        ax4_twin = ax4.twinx()
        l1, = ax4.plot(o_by_temp.index, o_by_temp['delta'], 'o-', color='green', 
                       linewidth=2, markersize=8, label='O: Mean δ')
        l2, = ax4_twin.plot(o_by_temp.index, o_by_temp['high_delta'] * 100, 's--', 
                           color='orange', linewidth=2, markersize=8, label='O: High δ %')
        ax4.set_xlabel('Temperature (K)', fontsize=12)
        ax4.set_ylabel('Mean δ (O atoms)', fontsize=12, color='green')
        ax4_twin.set_ylabel('High δ Ratio (%)', fontsize=12, color='orange')
        ax4.set_title('(d) Oxygen Atom Dynamics', fontsize=12, fontweight='bold')
        ax4.legend(handles=[l1, l2], loc='upper left')
        ax4.grid(True, alpha=0.3)
    else:
        ax4.text(0.5, 0.5, 'No O data available', ha='center', va='center', fontsize=14)
        ax4.set_title('(d) Oxygen Atom Dynamics', fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    output_path = output_dir / 'oxide_vs_nonoxide.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"  [SAVED] {output_path}")
    return fig


def plot_delta_vs_D_scatter(df, output_dir, structure=None, temp=None):
    """绘制 δ vs D 散点图"""
    print("\n[*] 绘制 δ vs D 散点图...")
    
    # 筛选数据
    df_plot = df.copy()
    if structure:
        df_plot = df_plot[df_plot['structure'].str.contains(structure, case=False)]
    if temp:
        df_plot = df_plot[df_plot['temp'] == temp]
    
    if len(df_plot) == 0:
        print("  [WARNING] 无数据可绘制")
        return
    
    # 创建图形
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # ===== 左图: 按元素着色 =====
    ax1 = axes[0]
    
    for element in df_plot['element'].unique():
        elem_data = df_plot[df_plot['element'] == element]
        color = ELEMENT_COLORS.get(element, 'gray')
        marker = ELEMENT_MARKERS.get(element, 'o')
        
        ax1.scatter(elem_data['delta'], elem_data['D'],
                   c=color, marker=marker, s=50, alpha=0.6,
                   label=f'{element} (n={len(elem_data)})')
    
    # 添加阈值线
    ax1.axvline(HIGH_LINDEMANN_THRESHOLD, color='red', linestyle='--', alpha=0.5, label=f'δ={HIGH_LINDEMANN_THRESHOLD}')
    ax1.axhline(HIGH_D_THRESHOLD, color='blue', linestyle='--', alpha=0.5, label=f'D={HIGH_D_THRESHOLD}')
    
    ax1.set_xlabel('Lindemann Index (δ)', fontsize=14)
    ax1.set_ylabel('Diffusion Coefficient D (10⁻⁵ cm²/s)', fontsize=14)
    ax1.set_title('δ vs D by Element', fontsize=14, fontweight='bold')
    ax1.legend(loc='upper left')
    ax1.set_yscale('log')
    ax1.grid(True, alpha=0.3)
    
    # ===== 右图: 按温度着色 =====
    ax2 = axes[1]
    
    temps = sorted(df_plot['temp'].unique())
    cmap = plt.cm.coolwarm
    norm = plt.Normalize(min(temps), max(temps))
    
    for temp_val in temps:
        temp_data = df_plot[df_plot['temp'] == temp_val]
        color = cmap(norm(temp_val))
        ax2.scatter(temp_data['delta'], temp_data['D'],
                   c=[color], s=50, alpha=0.6, label=f'{int(temp_val)}K')
    
    ax2.axvline(HIGH_LINDEMANN_THRESHOLD, color='red', linestyle='--', alpha=0.5)
    ax2.axhline(HIGH_D_THRESHOLD, color='blue', linestyle='--', alpha=0.5)
    
    ax2.set_xlabel('Lindemann Index (δ)', fontsize=14)
    ax2.set_ylabel('Diffusion Coefficient D (10⁻⁵ cm²/s)', fontsize=14)
    ax2.set_title('δ vs D by Temperature', fontsize=14, fontweight='bold')
    ax2.legend(loc='upper left', ncol=2, fontsize=9)
    ax2.set_yscale('log')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # 保存
    suffix = ''
    if structure:
        suffix += f'_{structure}'
    if temp:
        suffix += f'_{int(temp)}K'
    
    output_path = output_dir / f'delta_vs_D_scatter{suffix}.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"  [SAVED] {output_path}")
    
    return fig


def plot_element_comparison(df_stats, output_dir, structure=None):
    """绘制元素对比图"""
    print("\n[*] 绘制元素对比图...")
    
    df_plot = df_stats.copy()
    if structure:
        df_plot = df_plot[df_plot['structure'].str.contains(structure, case=False)]
    
    if len(df_plot) == 0:
        print("  [WARNING] 无数据可绘制")
        return
    
    # 读取structure-level的Lindemann指数 (与Step 7.8.5一致)
    df_lindemann = load_structure_level_lindemann()
    if df_lindemann is not None and structure:
        # 筛选相同结构
        df_lindemann = df_lindemann[df_lindemann['structure'].str.contains(structure, case=False)]
        print(f"  [INFO] 使用structure-level Lindemann数据 (系综平均方法)")
    
    # 改为3x2布局,添加e、f两个比值图
    fig, axes = plt.subplots(3, 2, figsize=(14, 18))
    
    # ===== (1) δ mean vs Temperature =====
    ax1 = axes[0, 0]
    
    # 如果有structure_level数据,优先使用
    if df_lindemann is not None and len(df_lindemann) > 0:
        for element in ['Pt', 'Sn']:
            if element == 'Pt':
                temps = df_lindemann['temp'].values
                delta_mean = df_lindemann['delta_Pt'].values
                delta_std = df_lindemann['delta_Pt_std'].values
            else:
                temps = df_lindemann['temp'].values
                delta_mean = df_lindemann['delta_Sn'].values
                delta_std = df_lindemann['delta_Sn_std'].values
            
            color = ELEMENT_COLORS.get(element, 'gray')
            marker = ELEMENT_MARKERS.get(element, 'o')
            
            ax1.plot(temps, delta_mean,
                    marker=marker, color=color, linewidth=2, markersize=8,
                    label=f'{element} (structure-level)')
            ax1.fill_between(temps,
                            delta_mean - delta_std,
                            delta_mean + delta_std,
                            color=color, alpha=0.2)
    else:
        # 降级使用per-atom平均
        for element in df_plot['element'].unique():
            elem_data = df_plot[df_plot['element'] == element].sort_values('temp')
            color = ELEMENT_COLORS.get(element, 'gray')
            marker = ELEMENT_MARKERS.get(element, 'o')
            
            ax1.plot(elem_data['temp'], elem_data['delta_mean'],
                    marker=marker, color=color, linewidth=2, markersize=8,
                    label=f'{element} (per-atom avg)')
            ax1.fill_between(elem_data['temp'],
                            elem_data['delta_mean'] - elem_data['delta_std'],
                            elem_data['delta_mean'] + elem_data['delta_std'],
                            color=color, alpha=0.2)
    
    ax1.axhline(HIGH_LINDEMANN_THRESHOLD, color='red', linestyle='--', alpha=0.5)
    ax1.set_xlabel('Temperature (K)', fontsize=12)
    ax1.set_ylabel('Mean Lindemann Index (δ)', fontsize=12)
    ax1.set_title('(a) δ vs Temperature', fontsize=13, fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # ===== (2) D mean vs Temperature (改进: 突出Pt vs Sn差异) =====
    ax2 = axes[0, 1]
    for element in df_plot['element'].unique():
        elem_data = df_plot[df_plot['element'] == element].sort_values('temp')
        color = ELEMENT_COLORS.get(element, 'gray')
        marker = ELEMENT_MARKERS.get(element, 'o')
        
        # 绘制主线
        ax2.plot(elem_data['temp'], elem_data['D_mean'],
                marker=marker, color=color, linewidth=2.5, markersize=10,
                label=element, zorder=3)
        
        # 添加误差带(使用D_std)
        if 'D_std' in elem_data.columns:
            ax2.fill_between(elem_data['temp'],
                            elem_data['D_mean'] - elem_data['D_std'],
                            elem_data['D_mean'] + elem_data['D_std'],
                            color=color, alpha=0.15, zorder=1)
    
    # 不使用对数坐标,使用线性坐标以突出差异
    ax2.set_xlabel('Temperature (K)', fontsize=13, fontweight='bold')
    ax2.set_ylabel('Mean D (×10⁻⁵ cm²/s)', fontsize=13, fontweight='bold')
    ax2.set_title('(b) Diffusion Coefficient vs Temperature\n(GMX Ensemble MSD Method)', 
                  fontsize=14, fontweight='bold')
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3, linestyle=':', linewidth=0.8)
    
    # 使用对数坐标Y轴 (模仿d图)
    ax2.set_yscale('log')
    # 设置Y轴范围: 自动范围,留余量
    y_min = df_plot['D_mean'].min() * 0.5
    y_max = df_plot['D_mean'].max() * 2.0
    ax2.set_ylim(y_min, y_max)
    
    # 添加文本说明
    ax2.text(0.02, 0.98, 'Note: D calculation uses same method\n(no difference for structure-level vs per-atom)',
            transform=ax2.transAxes, fontsize=9, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))
    
    # ===== (3) High activity ratio vs Temperature =====
    ax3 = axes[1, 0]
    for element in df_plot['element'].unique():
        elem_data = df_plot[df_plot['element'] == element].sort_values('temp')
        color = ELEMENT_COLORS.get(element, 'gray')
        marker = ELEMENT_MARKERS.get(element, 'o')
        
        ax3.plot(elem_data['temp'], elem_data['high_delta_ratio'] * 100,
                marker=marker, color=color, linewidth=2, markersize=8,
                label=f'{element} (δ>{HIGH_LINDEMANN_THRESHOLD})')
    
    ax3.set_xlabel('Temperature (K)', fontsize=12)
    ax3.set_ylabel('High Vibration Ratio (%)', fontsize=12)
    ax3.set_title('(c) High δ Ratio vs Temperature', fontsize=13, fontweight='bold')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    ax3.set_ylim(0, 105)
    
    # ===== (4) δ-D correlation (改进: Sn方块+Pt圆点,颜色深浅表温度) =====
    ax4 = axes[1, 1]
    
    # 获取温度范围,用于颜色映射
    temps = sorted(df_plot['temp'].unique())
    temp_min = min(temps)
    temp_max = max(temps)
    
    # 使用颜色映射 (冷色到暖色: 低温到高温)
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable
    norm = Normalize(vmin=temp_min, vmax=temp_max)
    cmap = plt.cm.RdYlBu_r  # 蓝色(冷)到红色(热)
    
    # 按元素绘制,使用不同marker
    for element in sorted(df_plot['element'].unique()):
        elem_data = df_plot[df_plot['element'] == element].sort_values('temp')
        
        # Sn用方块,Pt用圆点
        if element == 'Sn':
            marker = 's'  # 方块
            marker_size = 150
        elif element == 'Pt':
            marker = 'o'  # 圆点
            marker_size = 150
        else:
            marker = '^'  # 三角形(其他元素)
            marker_size = 120
        
        # 为每个温度点绘制不同颜色
        for _, row in elem_data.iterrows():
            temp = row['temp']
            color = cmap(norm(temp))
            
            ax4.scatter(row['delta_mean'], row['D_mean'],
                       c=[color], marker=marker, s=marker_size, 
                       alpha=0.8, edgecolors='black', linewidths=1.5,
                       zorder=3)
    
    # 添加图例 (元素类型)
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='gray',
               markersize=10, label='Pt', markeredgecolor='black', markeredgewidth=1.5),
        Line2D([0], [0], marker='s', color='w', markerfacecolor='gray',
               markersize=10, label='Sn', markeredgecolor='black', markeredgewidth=1.5),
    ]
    ax4.legend(handles=legend_elements, loc='upper left', fontsize=11)
    
    # 添加颜色条 (温度)
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax4)
    cbar.set_label('Temperature (K)', fontsize=11)
    
    ax4.set_xlabel('Mean Lindemann Index (δ)', fontsize=12)
    ax4.set_ylabel('Mean D (10⁻⁵ cm²/s)', fontsize=12)
    ax4.set_title('(d) δ-D Correlation (color=T, shape=element)', fontsize=13, fontweight='bold')
    ax4.set_yscale('log')
    ax4.grid(True, alpha=0.3, linestyle=':', linewidth=0.8)
    
    # ===== (e) Diffusion Coefficient Ratio (D_Sn/D_Pt) vs Temperature =====
    ax5 = axes[2, 0]
    
    if 'Pt' in df_plot['element'].values and 'Sn' in df_plot['element'].values:
        temps = sorted(df_plot['temp'].unique())
        D_ratios = []
        D_ratio_errors = []
        
        for temp in temps:
            pt_data = df_plot[(df_plot['element']=='Pt') & (df_plot['temp']==temp)]
            sn_data = df_plot[(df_plot['element']=='Sn') & (df_plot['temp']==temp)]
            
            if len(pt_data) > 0 and len(sn_data) > 0:
                pt_D = pt_data['D_mean'].values[0]
                sn_D = sn_data['D_mean'].values[0]
                ratio = sn_D / pt_D
                D_ratios.append(ratio)
                
                # 误差传播: ratio_error = ratio * sqrt((sn_std/sn_D)^2 + (pt_std/pt_D)^2)
                if 'D_std' in pt_data.columns and 'D_std' in sn_data.columns:
                    pt_std = pt_data['D_std'].values[0]
                    sn_std = sn_data['D_std'].values[0]
                    ratio_error = ratio * np.sqrt((sn_std/sn_D)**2 + (pt_std/pt_D)**2)
                    D_ratio_errors.append(ratio_error)
                else:
                    D_ratio_errors.append(0)
        
        # 绘制比值曲线
        ax5.errorbar(temps, D_ratios, yerr=D_ratio_errors,
                    marker='o', color='purple', linewidth=2.5, markersize=10,
                    capsize=5, capthick=2, label='D_Sn/D_Pt')
        
        # 添加参考线 y=1
        ax5.axhline(y=1.0, color='gray', linestyle='--', linewidth=1.5, 
                   label='Equal diffusion (ratio=1)', alpha=0.7)
        
        # 标注平均比值
        mean_ratio = np.mean(D_ratios)
        ax5.text(0.05, 0.95, f'Mean ratio = {mean_ratio:.3f}',
                transform=ax5.transAxes, fontsize=12, fontweight='bold',
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    ax5.set_xlabel('Temperature (K)', fontsize=13, fontweight='bold')
    ax5.set_ylabel('D_Sn / D_Pt', fontsize=13, fontweight='bold')
    ax5.set_title('(e) Diffusion Coefficient Ratio vs Temperature', fontsize=14, fontweight='bold')
    ax5.legend(fontsize=11)
    ax5.grid(True, alpha=0.3, linestyle=':', linewidth=0.8)
    ax5.set_ylim(0.8, max(D_ratios)*1.15 if D_ratios else 1.5)
    
    # ===== (f) Lindemann Index Ratio (δ_Pt/δ_Sn) vs Temperature =====
    ax6 = axes[2, 1]
    
    # 优先使用structure_level数据计算比值
    if df_lindemann is not None and len(df_lindemann) > 0:
        temps = sorted(df_lindemann['temp'].unique())
        delta_ratios = []
        delta_ratio_errors = []
        
        for temp in temps:
            temp_data = df_lindemann[df_lindemann['temp'] == temp]
            if len(temp_data) > 0:
                pt_delta = temp_data['delta_Pt'].values[0]
                sn_delta = temp_data['delta_Sn'].values[0]
                ratio = pt_delta / sn_delta  # Pt/Sn (Pt活性更高)
                delta_ratios.append(ratio)
                
                # 误差传播
                pt_std = temp_data['delta_Pt_std'].values[0]
                sn_std = temp_data['delta_Sn_std'].values[0]
                ratio_error = ratio * np.sqrt((pt_std/pt_delta)**2 + (sn_std/sn_delta)**2)
                delta_ratio_errors.append(ratio_error)
        
        # 绘制比值曲线
        ax6.errorbar(temps, delta_ratios, yerr=delta_ratio_errors,
                    marker='s', color='teal', linewidth=2.5, markersize=10,
                    capsize=5, capthick=2, label='δ_Pt/δ_Sn (structure-level)')
        
        # 添加参考线 y=1
        ax6.axhline(y=1.0, color='gray', linestyle='--', linewidth=1.5,
                   label='Equal mobility (ratio=1)', alpha=0.7)
        
        # 标注平均比值
        mean_ratio = np.mean(delta_ratios)
        ax6.text(0.05, 0.95, f'Mean ratio = {mean_ratio:.3f}',
                transform=ax6.transAxes, fontsize=12, fontweight='bold',
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
    
    elif 'Pt' in df_plot['element'].values and 'Sn' in df_plot['element'].values:
        # 降级使用per-atom平均数据
        temps = sorted(df_plot['temp'].unique())
        delta_ratios = []
        delta_ratio_errors = []
        
        for temp in temps:
            pt_data = df_plot[(df_plot['element']=='Pt') & (df_plot['temp']==temp)]
            sn_data = df_plot[(df_plot['element']=='Sn') & (df_plot['temp']==temp)]
            
            if len(pt_data) > 0 and len(sn_data) > 0:
                pt_delta = pt_data['delta_mean'].values[0]
                sn_delta = sn_data['delta_mean'].values[0]
                ratio = pt_delta / sn_delta  # 改为 Pt/Sn (Pt活性更高)
                delta_ratios.append(ratio)
                
                # 误差传播
                if 'delta_std' in pt_data.columns and 'delta_std' in sn_data.columns:
                    pt_std = pt_data['delta_std'].values[0]
                    sn_std = sn_data['delta_std'].values[0]
                    ratio_error = ratio * np.sqrt((pt_std/pt_delta)**2 + (sn_std/sn_delta)**2)
                    delta_ratio_errors.append(ratio_error)
                else:
                    delta_ratio_errors.append(0)
        
        # 绘制比值曲线
        ax6.errorbar(temps, delta_ratios, yerr=delta_ratio_errors,
                    marker='s', color='orange', linewidth=2.5, markersize=10,
                    capsize=5, capthick=2, label='δ_Pt/δ_Sn (per-atom avg)')
        
        # 添加参考线 y=1
        ax6.axhline(y=1.0, color='gray', linestyle='--', linewidth=1.5,
                   label='Equal mobility (ratio=1)', alpha=0.7)
        
        # 标注平均比值
        mean_ratio = np.mean(delta_ratios)
        ax6.text(0.05, 0.95, f'Mean ratio = {mean_ratio:.3f}',
                transform=ax6.transAxes, fontsize=12, fontweight='bold',
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
    
    ax6.set_xlabel('Temperature (K)', fontsize=13, fontweight='bold')
    ax6.set_ylabel('δ_Pt / δ_Sn', fontsize=13, fontweight='bold')
    ax6.set_title('(f) Lindemann Index Ratio vs Temperature', fontsize=14, fontweight='bold')
    ax6.legend(fontsize=11)
    ax6.grid(True, alpha=0.3, linestyle=':', linewidth=0.8)
    # Pt活性更高,所以比值>1,范围调整为1.0-1.6
    if delta_ratios:
        ax6.set_ylim(1.0, max(delta_ratios)*1.1)
    else:
        ax6.set_ylim(1.0, 1.6)
    
    plt.tight_layout()
    
    suffix = f'_{structure}' if structure else ''
    output_path = output_dir / f'element_comparison{suffix}.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"  [SAVED] {output_path}")
    
    return fig


def plot_mobility_distribution(df, output_dir, temp=None):
    """绘制活跃度分布图"""
    print("\n[*] 绘制活跃度分布图...")
    
    df_plot = df.copy()
    if temp:
        df_plot = df_plot[df_plot['temp'] == temp]
    
    if len(df_plot) == 0:
        print("  [WARNING] 无数据可绘制")
        return
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # ===== 左图: 活跃度分类饼图 =====
    ax1 = axes[0]
    
    mobility_counts = df_plot['mobility_class'].value_counts()
    colors_map = {
        'Active': '#d62728',     # 红色 - 活跃
        'Vibrating': '#ff7f0e',  # 橙色 - 局域振动
        'Diffusing': '#1f77b4',  # 蓝色 - 跳跃扩散
        'Stable': '#2ca02c',     # 绿色 - 稳定
    }
    colors = [colors_map.get(c, 'gray') for c in mobility_counts.index]
    
    wedges, texts, autotexts = ax1.pie(mobility_counts, labels=mobility_counts.index,
                                        autopct='%1.1f%%', colors=colors,
                                        explode=[0.02]*len(mobility_counts))
    ax1.set_title(f'Atom Mobility Classification\n(T={int(temp)}K)' if temp else 'Atom Mobility Classification',
                  fontsize=13, fontweight='bold')
    
    # ===== 右图: 按元素的活跃度分布 =====
    ax2 = axes[1]
    
    # 计算每个元素的活跃度分布
    mobility_by_element = df_plot.groupby(['element', 'mobility_class']).size().unstack(fill_value=0)
    mobility_by_element_pct = mobility_by_element.div(mobility_by_element.sum(axis=1), axis=0) * 100
    
    # 按顺序绘制
    class_order = ['Stable', 'Vibrating', 'Diffusing', 'Active']
    available_classes = [c for c in class_order if c in mobility_by_element_pct.columns]
    
    mobility_by_element_pct[available_classes].plot(kind='bar', stacked=True, ax=ax2,
                                                     color=[colors_map.get(c, 'gray') for c in available_classes])
    
    ax2.set_xlabel('Element', fontsize=12)
    ax2.set_ylabel('Percentage (%)', fontsize=12)
    ax2.set_title('Mobility Distribution by Element', fontsize=13, fontweight='bold')
    ax2.legend(title='Mobility Class', bbox_to_anchor=(1.02, 1))
    ax2.set_xticklabels(ax2.get_xticklabels(), rotation=0)
    ax2.set_ylim(0, 100)
    
    plt.tight_layout()
    
    suffix = f'_{int(temp)}K' if temp else ''
    output_path = output_dir / f'mobility_distribution{suffix}.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"  [SAVED] {output_path}")
    
    return fig


def plot_temperature_heatmap(df_stats, output_dir):
    """绘制温度-元素热力图"""
    print("\n[*] 绘制温度-元素热力图...")
    
    if len(df_stats) == 0:
        print("  [WARNING] 无数据可绘制")
        return
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # ===== 左图: δ_mean 热力图 =====
    ax1 = axes[0]
    
    pivot_delta = df_stats.pivot_table(index='element', columns='temp', values='delta_mean', aggfunc='mean')
    sns.heatmap(pivot_delta, annot=True, fmt='.3f', cmap='YlOrRd', ax=ax1,
                cbar_kws={'label': 'Mean δ'})
    ax1.set_title('(a) Mean Lindemann Index (δ)', fontsize=13, fontweight='bold')
    ax1.set_xlabel('Temperature (K)', fontsize=12)
    ax1.set_ylabel('Element', fontsize=12)
    
    # ===== 右图: D_mean 热力图 =====
    ax2 = axes[1]
    
    pivot_D = df_stats.pivot_table(index='element', columns='temp', values='D_mean', aggfunc='mean')
    # 使用对数刻度
    pivot_D_log = np.log10(pivot_D + 1e-10)
    sns.heatmap(pivot_D_log, annot=True, fmt='.2f', cmap='YlOrRd', ax=ax2,
                cbar_kws={'label': 'log₁₀(D)'})
    ax2.set_title('(b) Mean Diffusion Coefficient (log₁₀D)', fontsize=13, fontweight='bold')
    ax2.set_xlabel('Temperature (K)', fontsize=12)
    ax2.set_ylabel('Element', fontsize=12)
    
    plt.tight_layout()
    
    output_path = output_dir / 'temperature_element_heatmap.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"  [SAVED] {output_path}")
    
    return fig


# ============================================================================
# 特征温度对比 (300K vs 900K)
# ============================================================================

def plot_feature_temperature_comparison(df, output_dir):
    """绘制特征温度对比图 (300K vs 900K) - 按分组分析"""
    print("\n[*] 绘制特征温度对比图 (300K vs 900K)...")
    
    if 'series' not in df.columns:
        print("  [WARNING] 无分系列信息")
        return None
    
    # 筛选特征温度数据
    df_feature = df[df['temp'].isin(FEATURE_TEMPS)]
    if len(df_feature) == 0:
        print(f"  [WARNING] 无特征温度 {FEATURE_TEMPS} 的数据")
        # 尝试找到最接近的温度
        available_temps = sorted(df['temp'].unique())
        print(f"  可用温度: {available_temps}")
        return None
    
    print(f"  特征温度数据: {len(df_feature)} records")
    
    # 定义分组
    # 主要分组: air, sum8, pt8snx, pt6snx
    # 含氧分组: o1, o2, o3, o4
    main_series = ['air', 'sum8', 'pt8snx', 'pt6snx']
    oxide_series = ['o1', 'o2', 'o3', 'o4']
    
    available_main = [s for s in main_series if s in df_feature['series'].unique()]
    available_oxide = [s for s in oxide_series if s in df_feature['series'].unique()]
    
    print(f"  主分组: {available_main}")
    print(f"  含氧分组: {available_oxide}")
    
    # ===== 图1: 主分组 300K vs 900K 对比 =====
    if len(available_main) > 0:
        fig1, axes1 = plt.subplots(2, 3, figsize=(18, 12))
        
        # 第一行: 300K
        # 第二行: 900K
        
        for row_idx, temp in enumerate(FEATURE_TEMPS):
            df_temp = df_feature[df_feature['temp'] == temp]
            
            # (1) δ 分布
            ax = axes1[row_idx, 0]
            data_to_plot = []
            labels = []
            colors = []
            for series in available_main:
                for elem in ['Pt', 'Sn']:
                    subset = df_temp[(df_temp['series'] == series) & (df_temp['element'] == elem)]
                    if len(subset) > 0:
                        data_to_plot.append(subset['delta'].values)
                        labels.append(f'{series}\n{elem}')
                        colors.append(SERIES_COLORS.get(series, 'gray'))
            
            if data_to_plot:
                bp = ax.boxplot(data_to_plot, labels=labels, patch_artist=True)
                for i, (patch, c) in enumerate(zip(bp['boxes'], colors)):
                    patch.set_facecolor(c)
                    patch.set_alpha(0.6)
            
            ax.axhline(HIGH_LINDEMANN_THRESHOLD, color='red', linestyle='--', alpha=0.5)
            ax.set_ylabel('Lindemann Index (δ)', fontsize=11)
            ax.set_title(f'{int(temp)}K: δ Distribution', fontsize=12, fontweight='bold')
            ax.grid(True, alpha=0.3)
            ax.tick_params(axis='x', rotation=45)
            
            # (2) D 分布
            ax = axes1[row_idx, 1]
            data_to_plot = []
            labels = []
            colors = []
            for series in available_main:
                for elem in ['Pt', 'Sn']:
                    subset = df_temp[(df_temp['series'] == series) & (df_temp['element'] == elem)]
                    if len(subset) > 0:
                        data_to_plot.append(subset['D'].values)
                        labels.append(f'{series}\n{elem}')
                        colors.append(SERIES_COLORS.get(series, 'gray'))
            
            if data_to_plot:
                bp = ax.boxplot(data_to_plot, labels=labels, patch_artist=True)
                for i, (patch, c) in enumerate(zip(bp['boxes'], colors)):
                    patch.set_facecolor(c)
                    patch.set_alpha(0.6)
            
            ax.axhline(HIGH_D_THRESHOLD, color='blue', linestyle='--', alpha=0.5)
            ax.set_ylabel('D (10⁻⁵ cm²/s)', fontsize=11)
            ax.set_title(f'{int(temp)}K: D Distribution', fontsize=12, fontweight='bold')
            ax.set_yscale('log')
            ax.grid(True, alpha=0.3)
            ax.tick_params(axis='x', rotation=45)
            
            # (3) 活跃比例柱状图
            ax = axes1[row_idx, 2]
            series_labels = []
            pt_active = []
            sn_active = []
            
            for series in available_main:
                pt_data = df_temp[(df_temp['series'] == series) & (df_temp['element'] == 'Pt')]
                sn_data = df_temp[(df_temp['series'] == series) & (df_temp['element'] == 'Sn')]
                
                if len(pt_data) > 0 or len(sn_data) > 0:
                    series_labels.append(series)
                    if len(pt_data) > 0:
                        pt_active.append((pt_data['delta'] > HIGH_LINDEMANN_THRESHOLD).mean() * 100)
                    else:
                        pt_active.append(0)
                    if len(sn_data) > 0:
                        sn_active.append((sn_data['delta'] > HIGH_LINDEMANN_THRESHOLD).mean() * 100)
                    else:
                        sn_active.append(0)
            
            x = np.arange(len(series_labels))
            width = 0.35
            ax.bar(x - width/2, pt_active, width, label='Pt', color=ELEMENT_COLORS['Pt'], alpha=0.8)
            ax.bar(x + width/2, sn_active, width, label='Sn', color=ELEMENT_COLORS['Sn'], alpha=0.8)
            
            ax.set_xticks(x)
            ax.set_xticklabels(series_labels)
            ax.set_ylabel('High δ Ratio (%)', fontsize=11)
            ax.set_title(f'{int(temp)}K: High Vibration Ratio', fontsize=12, fontweight='bold')
            ax.legend()
            ax.grid(True, alpha=0.3, axis='y')
            ax.set_ylim(0, 105)
        
        plt.suptitle('Feature Temperature Comparison: Main Series (Air, Sum8, Pt8SnX, Pt6SnX)\n'
                     'Upper: 300K (Low T), Lower: 900K (High T)', 
                     fontsize=14, fontweight='bold', y=1.02)
        plt.tight_layout()
        
        output_path = output_dir / 'feature_temp_main_series.png'
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"  [SAVED] {output_path}")
    
    # ===== 图2: 含氧分组 300K vs 900K 对比 =====
    if len(available_oxide) > 0:
        fig2, axes2 = plt.subplots(2, 3, figsize=(18, 12))
        
        for row_idx, temp in enumerate(FEATURE_TEMPS):
            df_temp = df_feature[df_feature['temp'] == temp]
            
            # (1) δ 分布
            ax = axes2[row_idx, 0]
            data_to_plot = []
            labels = []
            colors = []
            for series in available_oxide:
                for elem in ['Pt', 'Sn', 'O']:
                    subset = df_temp[(df_temp['series'] == series) & (df_temp['element'] == elem)]
                    if len(subset) > 0:
                        data_to_plot.append(subset['delta'].values)
                        labels.append(f'{series}\n{elem}')
                        colors.append(SERIES_COLORS.get(series, 'gray'))
            
            if data_to_plot:
                bp = ax.boxplot(data_to_plot, labels=labels, patch_artist=True)
                for i, (patch, c) in enumerate(zip(bp['boxes'], colors)):
                    patch.set_facecolor(c)
                    patch.set_alpha(0.6)
            
            ax.axhline(HIGH_LINDEMANN_THRESHOLD, color='red', linestyle='--', alpha=0.5)
            ax.set_ylabel('Lindemann Index (δ)', fontsize=11)
            ax.set_title(f'{int(temp)}K: δ Distribution (Oxide)', fontsize=12, fontweight='bold')
            ax.grid(True, alpha=0.3)
            ax.tick_params(axis='x', rotation=45)
            
            # (2) D 分布
            ax = axes2[row_idx, 1]
            data_to_plot = []
            labels = []
            colors = []
            for series in available_oxide:
                for elem in ['Pt', 'Sn', 'O']:
                    subset = df_temp[(df_temp['series'] == series) & (df_temp['element'] == elem)]
                    if len(subset) > 0:
                        data_to_plot.append(subset['D'].values)
                        labels.append(f'{series}\n{elem}')
                        colors.append(SERIES_COLORS.get(series, 'gray'))
            
            if data_to_plot:
                bp = ax.boxplot(data_to_plot, labels=labels, patch_artist=True)
                for i, (patch, c) in enumerate(zip(bp['boxes'], colors)):
                    patch.set_facecolor(c)
                    patch.set_alpha(0.6)
            
            ax.axhline(HIGH_D_THRESHOLD, color='blue', linestyle='--', alpha=0.5)
            ax.set_ylabel('D (10⁻⁵ cm²/s)', fontsize=11)
            ax.set_title(f'{int(temp)}K: D Distribution (Oxide)', fontsize=12, fontweight='bold')
            ax.set_yscale('log')
            ax.grid(True, alpha=0.3)
            ax.tick_params(axis='x', rotation=45)
            
            # (3) 元素对比柱状图
            ax = axes2[row_idx, 2]
            series_labels = []
            pt_delta = []
            sn_delta = []
            o_delta = []
            
            for series in available_oxide:
                pt_data = df_temp[(df_temp['series'] == series) & (df_temp['element'] == 'Pt')]
                sn_data = df_temp[(df_temp['series'] == series) & (df_temp['element'] == 'Sn')]
                o_data = df_temp[(df_temp['series'] == series) & (df_temp['element'] == 'O')]
                
                if len(pt_data) > 0 or len(sn_data) > 0 or len(o_data) > 0:
                    series_labels.append(series)
                    pt_delta.append(pt_data['delta'].mean() if len(pt_data) > 0 else 0)
                    sn_delta.append(sn_data['delta'].mean() if len(sn_data) > 0 else 0)
                    o_delta.append(o_data['delta'].mean() if len(o_data) > 0 else 0)
            
            x = np.arange(len(series_labels))
            width = 0.25
            ax.bar(x - width, pt_delta, width, label='Pt', color=ELEMENT_COLORS['Pt'], alpha=0.8)
            ax.bar(x, sn_delta, width, label='Sn', color=ELEMENT_COLORS['Sn'], alpha=0.8)
            ax.bar(x + width, o_delta, width, label='O', color=ELEMENT_COLORS['O'], alpha=0.8)
            
            ax.set_xticks(x)
            ax.set_xticklabels(series_labels)
            ax.set_ylabel('Mean δ', fontsize=11)
            ax.set_title(f'{int(temp)}K: Mean δ by Element', fontsize=12, fontweight='bold')
            ax.legend()
            ax.grid(True, alpha=0.3, axis='y')
            ax.axhline(HIGH_LINDEMANN_THRESHOLD, color='red', linestyle='--', alpha=0.5)
        
        plt.suptitle('Feature Temperature Comparison: Oxide Series (O1, O2, O3, O4)\n'
                     'Upper: 300K (Low T), Lower: 900K (High T)', 
                     fontsize=14, fontweight='bold', y=1.02)
        plt.tight_layout()
        
        output_path = output_dir / 'feature_temp_oxide_series.png'
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"  [SAVED] {output_path}")
    
    # ===== 图3: 汇总对比 - 雷达图/热力图 =====
    fig3, axes3 = plt.subplots(1, 2, figsize=(16, 7))
    
    # (1) 300K vs 900K 热力图
    all_series = available_main + available_oxide
    if len(all_series) > 0:
        ax = axes3[0]
        
        # 构建数据: 行=series, 列=[Pt-300K, Sn-300K, Pt-900K, Sn-900K]
        cols = ['Pt@300K', 'Sn@300K', 'Pt@900K', 'Sn@900K']
        data_matrix = []
        series_labels = []
        
        for series in all_series:
            row = []
            for temp in FEATURE_TEMPS:
                for elem in ['Pt', 'Sn']:
                    subset = df_feature[(df_feature['series'] == series) & 
                                        (df_feature['temp'] == temp) & 
                                        (df_feature['element'] == elem)]
                    if len(subset) > 0:
                        row.append(subset['delta'].mean())
                    else:
                        row.append(np.nan)
            if not all(np.isnan(row)):
                data_matrix.append(row)
                series_labels.append(series)
        
        if data_matrix:
            df_heatmap = pd.DataFrame(data_matrix, index=series_labels, columns=cols)
            sns.heatmap(df_heatmap, annot=True, fmt='.3f', cmap='RdYlGn_r', ax=ax,
                       vmin=0, vmax=0.15, linewidths=0.5)
            ax.set_title('Mean δ: 300K vs 900K Comparison\n(Green=Stable, Red=Active)', 
                        fontsize=12, fontweight='bold')
            ax.set_ylabel('Series', fontsize=11)
    
    # (2) 活跃度变化量柱状图
    ax = axes3[1]
    
    series_labels = []
    delta_change_pt = []  # 900K - 300K
    delta_change_sn = []
    
    for series in all_series:
        pt_300 = df_feature[(df_feature['series'] == series) & 
                            (df_feature['temp'] == 300) & 
                            (df_feature['element'] == 'Pt')]['delta'].mean()
        pt_900 = df_feature[(df_feature['series'] == series) & 
                            (df_feature['temp'] == 900) & 
                            (df_feature['element'] == 'Pt')]['delta'].mean()
        sn_300 = df_feature[(df_feature['series'] == series) & 
                            (df_feature['temp'] == 300) & 
                            (df_feature['element'] == 'Sn')]['delta'].mean()
        sn_900 = df_feature[(df_feature['series'] == series) & 
                            (df_feature['temp'] == 900) & 
                            (df_feature['element'] == 'Sn')]['delta'].mean()
        
        if not (np.isnan(pt_300) or np.isnan(pt_900) or np.isnan(sn_300) or np.isnan(sn_900)):
            series_labels.append(series)
            delta_change_pt.append(pt_900 - pt_300)
            delta_change_sn.append(sn_900 - sn_300)
    
    if series_labels:
        x = np.arange(len(series_labels))
        width = 0.35
        ax.bar(x - width/2, delta_change_pt, width, label='Pt', color=ELEMENT_COLORS['Pt'], alpha=0.8)
        ax.bar(x + width/2, delta_change_sn, width, label='Sn', color=ELEMENT_COLORS['Sn'], alpha=0.8)
        
        ax.set_xticks(x)
        ax.set_xticklabels(series_labels, rotation=45, ha='right')
        ax.set_ylabel('Δδ (900K - 300K)', fontsize=11)
        ax.set_title('Temperature-Induced Mobility Change\n(Higher = More Temperature-Sensitive)', 
                    fontsize=12, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')
        ax.axhline(0, color='black', linestyle='-', alpha=0.3)
    
    plt.tight_layout()
    
    output_path = output_dir / 'feature_temp_summary.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"  [SAVED] {output_path}")
    
    return fig3


# ============================================================================
# 报告生成
# ============================================================================

def generate_report(df, df_stats, df_corr, output_dir, df_series_stats=None):
    """生成分析报告 (包含分系列分析)"""
    print("\n[*] 生成分析报告...")
    
    report_lines = []
    report_lines.append("=" * 80)
    report_lines.append("Step 7.2: Per-Atom 动力学分析报告")
    report_lines.append("=" * 80)
    report_lines.append("")
    
    # 1. 数据概况
    report_lines.append("1. 数据概况")
    report_lines.append("-" * 40)
    report_lines.append(f"   总原子记录数: {len(df)}")
    report_lines.append(f"   结构数: {df['structure'].nunique()}")
    report_lines.append(f"   温度范围: {df['temp'].min():.0f}K - {df['temp'].max():.0f}K")
    report_lines.append(f"   元素: {', '.join(df['element'].unique())}")
    
    # 分系列统计
    if 'series' in df.columns:
        report_lines.append("\n   体系分类:")
        for series in df['series'].unique():
            count = len(df[df['series'] == series])
            pct = count / len(df) * 100
            n_struct = df[df['series'] == series]['structure'].nunique()
            report_lines.append(f"      {series}: {count} records ({pct:.1f}%), {n_struct} structures")
    
    report_lines.append("")
    
    # 2. 分系列元素统计 (核心！)
    report_lines.append("2. 分系列统计 (气相 vs 负载 vs 含氧)")
    report_lines.append("-" * 40)
    
    if 'series' in df.columns:
        for series in df['series'].unique():
            series_data = df[df['series'] == series]
            report_lines.append(f"\n   【{series.upper()}】")
            
            for element in ['Pt', 'Sn', 'O']:
                elem_data = series_data[series_data['element'] == element]
                if len(elem_data) == 0:
                    continue
                
                report_lines.append(f"      {element}: n={len(elem_data)}")
                report_lines.append(f"         δ = {elem_data['delta'].mean():.4f} ± {elem_data['delta'].std():.4f}")
                report_lines.append(f"         D = {elem_data['D'].mean():.4f} ± {elem_data['D'].std():.4f}")
                report_lines.append(f"         高振动比例: {(elem_data['delta'] > HIGH_LINDEMANN_THRESHOLD).mean()*100:.1f}%")
                report_lines.append(f"         高扩散比例: {(elem_data['D'] > HIGH_D_THRESHOLD).mean()*100:.1f}%")
    
    report_lines.append("")
    
    # 3. 全局元素统计
    report_lines.append("3. 全局元素统计 (所有体系合并)")
    report_lines.append("-" * 40)
    
    for element in df['element'].unique():
        elem_data = df[df['element'] == element]
        report_lines.append(f"\n   {element}:")
        report_lines.append(f"      原子数: {len(elem_data)}")
        report_lines.append(f"      δ: {elem_data['delta'].mean():.4f} ± {elem_data['delta'].std():.4f}")
        report_lines.append(f"      D: {elem_data['D'].mean():.4f} ± {elem_data['D'].std():.4f} (1e-5 cm²/s)")
        report_lines.append(f"      高振动比例 (δ>{HIGH_LINDEMANN_THRESHOLD}): {(elem_data['delta'] > HIGH_LINDEMANN_THRESHOLD).mean()*100:.1f}%")
        report_lines.append(f"      高扩散比例 (D>{HIGH_D_THRESHOLD}): {(elem_data['D'] > HIGH_D_THRESHOLD).mean()*100:.1f}%")
    
    report_lines.append("")
    
    # 4. 原子活跃度分类
    report_lines.append("4. 原子活跃度分类")
    report_lines.append("-" * 40)
    
    mobility_counts = df['mobility_class'].value_counts()
    for cls, count in mobility_counts.items():
        report_lines.append(f"   {cls}: {count} ({count/len(df)*100:.1f}%)")
    
    report_lines.append("")
    report_lines.append("   分类定义:")
    report_lines.append(f"   - Active: δ>{HIGH_LINDEMANN_THRESHOLD} AND D>{HIGH_D_THRESHOLD} (高振动+高扩散)")
    report_lines.append(f"   - Vibrating: δ>{HIGH_LINDEMANN_THRESHOLD} AND D<={HIGH_D_THRESHOLD} (局域振动)")
    report_lines.append(f"   - Diffusing: δ<={HIGH_LINDEMANN_THRESHOLD} AND D>{HIGH_D_THRESHOLD} (跳跃扩散)")
    report_lines.append(f"   - Stable: δ<={HIGH_LINDEMANN_THRESHOLD} AND D<={HIGH_D_THRESHOLD} (稳定)")
    report_lines.append("")
    
    # 5. 关键发现
    report_lines.append("5. 关键发现")
    report_lines.append("-" * 40)
    
    # 分系列发现
    if 'series' in df.columns:
        report_lines.append("\n   【分系列对比】")
        for series in df['series'].unique():
            series_data = df[df['series'] == series]
            
            # 找出该系列中哪个元素先活跃
            for element in ['Pt', 'Sn']:
                elem_data = series_data[series_data['element'] == element]
                if len(elem_data) == 0:
                    continue
                
                # 按温度统计高振动比例
                high_by_temp = elem_data.groupby('temp')['high_delta'].mean()
                first_high = high_by_temp[high_by_temp > 0.5].index
                if len(first_high) > 0:
                    report_lines.append(f"   {series}: {element} 首次>50%高振动 @ {first_high[0]:.0f}K")
    
    # 找出哪个元素先"动起来"
    if len(df_stats) > 0:
        # 按温度找第一个高活性元素
        high_temps = df_stats[df_stats['high_delta_ratio'] > 0.5].sort_values('temp')
        if len(high_temps) > 0:
            first_active = high_temps.iloc[0]
            report_lines.append(f"   首先出现大量高振动(>50%)的元素: {first_active['element']} at {first_active['temp']:.0f}K")
    
    # δ-D 相关性
    if len(df_corr) > 0:
        avg_corr = df_corr[df_corr['element'] == 'All']['pearson_r'].mean()
        report_lines.append(f"   delta-D avg correlation: r = {avg_corr:.3f}")
        if avg_corr > 0.5:
            report_lines.append("   -> High vibration atoms are usually high diffusion (consistent dynamics)")
        elif avg_corr < 0:
            report_lines.append("   -> Some atoms show anomalous behavior (vibration-diffusion decoupled)")
    
    report_lines.append("")
    report_lines.append("=" * 80)
    
    # 保存报告
    report_path = output_dir / 'analysis_report.txt'
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(report_lines))
    
    print(f"  [SAVED] {report_path}")
    
    # 打印简化版本
    print("\n" + "=" * 60)
    print("Report saved to:", report_path)
    print("=" * 60)


# ============================================================================
# 主函数
# ============================================================================

def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='Step 7.2: Per-Atom 动力学分析')
    parser.add_argument('--structure', type=str, default=None,
                       help='筛选特定结构 (支持部分匹配)')
    parser.add_argument('--temp-range', type=float, nargs=2, default=None,
                       help='温度范围 (min max)')
    parser.add_argument('--interactive', action='store_true',
                       help='交互模式 (显示图形)')
    parser.add_argument('--lindemann-dir', type=str, default=None,
                       help='Lindemann数据目录')
    parser.add_argument('--msd-dir', type=str, default=None,
                       help='MSD数据目录')
    parser.add_argument('--no-ensemble-D', action='store_true',
                       help='不使用系综平均计算D (使用旧的直接平均方法)')
    
    args = parser.parse_args()
    
    # 确定是否使用系综平均D
    use_ensemble_D = not args.no_ensemble_D
    
    print("=" * 80)
    print("Step 7.2: Per-Atom 动力学分析 - 林德曼指数 + 扩散系数")
    print("=" * 80)
    print(f"\n计算方法:")
    print(f"  扩散系数D: {'系综平均 (推荐)' if use_ensemble_D else '直接平均 (旧方法)'}")
    if use_ensemble_D:
        print(f"    - 先对每个run的原子D求平均")
        print(f"    - 再对所有run的平均D求平均")
        print(f"    - 更合理,考虑了run间的统计独立性")
    
    # 确定数据目录
    lindemann_dir = Path(args.lindemann_dir) if args.lindemann_dir else LINDEMANN_DIR
    msd_dir = Path(args.msd_dir) if args.msd_dir else MSD_DIR
    
    print(f"\n数据目录:")
    print(f"  Lindemann: {lindemann_dir}")
    print(f"  MSD: {msd_dir}")
    
    # 加载数据
    df_lindemann = load_lindemann_data(lindemann_dir)
    df_msd = load_msd_data(msd_dir)
    
    # 合并数据
    df_merged = merge_data(df_lindemann, df_msd, 
                          structure_filter=args.structure,
                          temp_range=args.temp_range)
    
    if df_merged is None or len(df_merged) == 0:
        print("\n[ERROR] 无法合并数据或数据为空")
        print("\n可能原因:")
        print("  1. per-atom Lindemann 和 MSD 数据的结构/温度/atom_id不匹配")
        print("  2. 数据文件不存在或格式错误")
        print("\n建议:")
        print("  - 检查两个数据源的列名和内容")
        print("  - 使用 --lindemann-dir 和 --msd-dir 指定正确路径")
        return
    
    # 保存合并数据
    merged_path = OUTPUT_DIR / 'merged_per_atom_data.csv'
    df_merged.to_csv(merged_path, index=False, encoding='utf-8-sig')
    print(f"\n[SAVED] {merged_path}")
    
    # ========== 输出详细数值 ==========
    print("\n" + "=" * 80)
    print("Per-Atom 数据详情")
    print("=" * 80)
    
    # 按温度和元素汇总
    for temp in sorted(df_merged['temp'].unique()):
        temp_data = df_merged[df_merged['temp'] == temp]
        print(f"\n【温度: {temp}K】")
        
        for element in sorted(temp_data['element'].unique()):
            elem_data = temp_data[temp_data['element'] == element]
            print(f"\n  {element} 元素 ({len(elem_data)} 个原子×run记录):")
            print(f"    Lindemann Index (δ):")
            print(f"      平均值: {elem_data['delta'].mean():.6f}")
            print(f"      标准差: {elem_data['delta'].std():.6f}")
            print(f"      中位数: {elem_data['delta'].median():.6f}")
            print(f"      范围:   [{elem_data['delta'].min():.6f}, {elem_data['delta'].max():.6f}]")
            
            print(f"    Diffusion Coefficient (D, 10^-5 cm^2/s):")
            print(f"      平均值: {elem_data['D'].mean():.6f}")
            print(f"      标准差: {elem_data['D'].std():.6f}")
            print(f"      中位数: {elem_data['D'].median():.6f}")
            print(f"      范围:   [{elem_data['D'].min():.6f}, {elem_data['D'].max():.6f}]")
            
            # 相关性
            if len(elem_data) >= 3:
                r, p = stats.pearsonr(elem_data['delta'], elem_data['D'])
                significance = "显著" if p < 0.05 else "不显著"
                print(f"    δ-D 相关性: r = {r:.4f}, p = {p:.4f} ({significance})")
        
        # 元素间对比
        pt_data = temp_data[temp_data['element'] == 'Pt']
        sn_data = temp_data[temp_data['element'] == 'Sn']
        if len(pt_data) > 0 and len(sn_data) > 0:
            print(f"\n  【Pt vs Sn 对比】")
            print(f"    D_Sn/D_Pt  = {sn_data['D'].mean() / pt_data['D'].mean():.4f}")
            print(f"    δ_Sn/δ_Pt  = {sn_data['delta'].mean() / pt_data['delta'].mean():.4f}")
    
    print("\n" + "=" * 80)
    
    # 计算统计 (使用GMX ensemble MSD, 与step7_5一致)
    df_stats = calculate_element_statistics(df_merged, use_ensemble_D=use_ensemble_D, use_gmx_msd=True)
    stats_path = OUTPUT_DIR / 'element_statistics.csv'
    df_stats.to_csv(stats_path, index=False, encoding='utf-8-sig')
    print(f"\n[SAVED] {stats_path}")
    
    # ========== 输出拟合质量和变异系数 (CV & R²) ==========
    if 'CV' in df_stats.columns and 'r2_mean' in df_stats.columns:
        print("\n" + "=" * 80)
        print("Element Statistics (GMX Ensemble MSD Method)")
        print("=" * 80)
        
        for temp in sorted(df_stats['temp'].unique()):
            temp_stats = df_stats[df_stats['temp'] == temp]
            print(f"\n[{int(temp)}K]")
            
            for element in sorted(temp_stats['element'].unique()):
                elem_stats = temp_stats[temp_stats['element'] == element]
                if len(elem_stats) > 0:
                    row = elem_stats.iloc[0]
                    print(f"  {element}:")
                    print(f"    D_mean = {row['D_mean']:.6f} +/- {row['D_std']:.6f} (x10^-5 cm^2/s)")
                    print(f"    CV = {row['CV']:.2f}% (coefficient of variation)")
                    print(f"    R^2 = {row['r2_mean']:.4f} (fit quality)")
                    print(f"    n_runs = {int(row['n_runs'])}")
            
            # 对比
            pt_stats = temp_stats[temp_stats['element'] == 'Pt']
            sn_stats = temp_stats[temp_stats['element'] == 'Sn']
            if len(pt_stats) > 0 and len(sn_stats) > 0:
                pt = pt_stats.iloc[0]
                sn = sn_stats.iloc[0]
                print(f"\n  [Quality Assessment]")
                
                # CV 评价
                if pt['CV'] < 5 and sn['CV'] < 5:
                    cv_grade = "Excellent (CV<5%)"
                elif pt['CV'] < 10 and sn['CV'] < 10:
                    cv_grade = "Good (CV<10%)"
                elif pt['CV'] < 20 and sn['CV'] < 20:
                    cv_grade = "Fair (CV<20%)"
                else:
                    cv_grade = "WARNING: Poor (CV>=20%)"
                
                # R² 评价
                if pt['r2_mean'] > 0.98 and sn['r2_mean'] > 0.98:
                    r2_grade = "Excellent (R^2>0.98)"
                elif pt['r2_mean'] > 0.95 and sn['r2_mean'] > 0.95:
                    r2_grade = "Good (R^2>0.95)"
                elif pt['r2_mean'] > 0.90 and sn['r2_mean'] > 0.90:
                    r2_grade = "Fair (R^2>0.90)"
                else:
                    r2_grade = "WARNING: Poor (R^2<=0.90)"
                
                print(f"    CV: {cv_grade}")
                print(f"    R^2: {r2_grade}")
        
        print("\n" + "=" * 80)
    
    # 计算分系列统计
    df_series_stats = calculate_series_statistics(df_merged, use_ensemble_D=use_ensemble_D)
    if df_series_stats is not None:
        series_stats_path = OUTPUT_DIR / 'series_statistics.csv'
        df_series_stats.to_csv(series_stats_path, index=False, encoding='utf-8-sig')
        print(f"[SAVED] {series_stats_path}")
    
    # 计算相关性
    df_corr = calculate_correlation(df_merged)
    corr_path = OUTPUT_DIR / 'correlation_analysis.csv'
    df_corr.to_csv(corr_path, index=False, encoding='utf-8-sig')
    print(f"[SAVED] {corr_path}")
    
    # ========== 绘图 ==========
    
    # 1. 分系列对比图 (核心图！)
    plot_series_comparison(df_series_stats, OUTPUT_DIR)
    
    # 2. 含氧 vs 无氧对比
    plot_oxide_vs_nonoxide(df_merged, OUTPUT_DIR)
    
    # 3. δ vs D 散点图
    plot_delta_vs_D_scatter(df_merged, OUTPUT_DIR, structure=args.structure)
    
    # 4. 元素对比图
    plot_element_comparison(df_stats, OUTPUT_DIR, structure=args.structure)
    
    # 5. 按温度绘制活跃度分布
    temps = sorted(df_merged['temp'].unique())
    if len(temps) > 0:
        # 选择几个代表性温度
        selected_temps = [temps[0], temps[len(temps)//2], temps[-1]]
        for t in selected_temps:
            plot_mobility_distribution(df_merged, OUTPUT_DIR, temp=t)
    
    # 6. 温度-元素热力图
    plot_temperature_heatmap(df_stats, OUTPUT_DIR)
    
    # 7. 特征温度对比 (300K vs 900K) - 核心分析图!
    plot_feature_temperature_comparison(df_merged, OUTPUT_DIR)
    
    # 生成报告 (包含分系列信息)
    generate_report(df_merged, df_stats, df_corr, OUTPUT_DIR, df_series_stats)
    
    print("\n" + "=" * 80)
    print("分析完成！输出目录:", OUTPUT_DIR)
    print("=" * 80)
    
    if args.interactive:
        plt.show()


if __name__ == '__main__':
    main()
