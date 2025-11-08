#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
统一多温度分析脚本 - v625数据 (2025-10-26)

支持三大系列分析,通过参数指定:
1. Pt8Snx: Pt8+0-10Sn系列 (原子数8-18)
2. PtxSn8-x: 总原子数=8的系列 (Pt8→Pt3Sn5)
3. Pt6Snx: Pt6+0-9Sn系列 (原子数6-15,缺pt6sn10)

用法:
    python step7-5-unified_multi_temp_v625_analysis.py --series Pt8Snx
    python step7-5-unified_multi_temp_v625_analysis.py --series PtxSn8-x
    python step7-5-unified_multi_temp_v625_analysis.py --series Pt6Snx
    python step7-5-unified_multi_temp_v625_analysis.py --all  # 运行所有系列

功能:
- 自动检测运行文件夹(使用V625DataLocator)
- 多运行平均(4-8次)
- 完整可视化(综合图3×10 + 热图2×3 + Q6对比2×3)
- 键类型统计(Pt-Pt, Pt-Sn, Sn-Sn)
- 温度效应和组分效应分析

Author: AI Assistant
Date: 2025-10-26
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager
from scipy.stats import pearsonr
import os
import sys
import warnings
import argparse
import re
from pathlib import Path
from v625_data_locator import V625DataLocator
import seaborn as sns

# 设置控制台输出编码 - 必须在最早导入后立即设置
if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8')
    # 设置默认编码
    import locale
    locale.setlocale(locale.LC_ALL, 'zh_CN.UTF-8')

warnings.filterwarnings('ignore')

# 配置中文字体
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'Arial']
plt.rcParams['axes.unicode_minus'] = False

def extract_path_signature(filepath, is_msd_path=True):
    """
    从文件路径提取4级路径签名 (与Step 7.4完全一致)
    
    Args:
        filepath: 完整文件路径
        is_msd_path: True=MSD路径(有温度目录), False=能量路径(无温度目录)
    
    Returns:
        path_signature: 4级路径签名,如 "run3/o2/o2pt4sn6/t1000.r24.gpu0"
                       或 "parent/composition/t1000.r24.gpu0" (无批次时)
    
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
    if not filepath or pd.isna(filepath):
        return None
    
    filepath = str(filepath).replace('\\', '/')
    
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

def load_msd_outliers():
    """
    加载Step 1的MSD异常值列表
    返回路径签名集合
    """
    outlier_file = Path(__file__).parent.parent / 'results' / 'large_D_outliers.csv'
    
    if not outlier_file.exists():
        print(f"\n[警告] MSD异常值文件不存在: {outlier_file}")
        print("  将不应用MSD过滤")
        return set()
    
    print(f"\n[MSD过滤] 加载异常值文件: {outlier_file}")
    
    try:
        df_outliers = pd.read_csv(outlier_file)
        
        if 'filepath' not in df_outliers.columns:
            print(f"  [错误] 文件缺少filepath列")
            return set()
        
        # 提取路径签名
        outlier_signatures = set()
        for filepath in df_outliers['filepath'].dropna():
            sig = extract_path_signature(filepath, is_msd_path=True)
            if sig:
                outlier_signatures.add(sig)
        
        print(f"  加载完成: {len(df_outliers)}条异常记录")
        print(f"  唯一路径签名: {len(outlier_signatures)}个")
        
        # 显示异常原因分布
        if 'reason' in df_outliers.columns:
            print(f"\n  异常原因分布:")
            reason_counts = df_outliers['reason'].value_counts()
            for reason, count in reason_counts.items():
                print(f"    {reason}: {count}条 ({count/len(df_outliers)*100:.1f}%)")
        
        return outlier_signatures
        
    except Exception as e:
        print(f"  [错误] 加载失败: {e}")
        return set()

class UnifiedMultiTempAnalyzer:
    """统一的多温度分析器"""
    
    # 系列配置
    SERIES_CONFIGS = {
        'Pt8Snx': {
            'name': 'Pt8Snx',
            'display_name': 'Pt8+Sn系列',
            'run_folder': 'Pt8',
            'systems': [
                ('pt8sn0-2-best', 0, 8, 0, 'Pt8'),
                ('pt8sn1-2-best', 1, 8, 1, 'Pt8Sn1'),
                ('pt8sn2-1-best', 2, 8, 2, 'Pt8Sn2'),
                ('pt8sn3-1-best', 3, 8, 3, 'Pt8Sn3'),
                ('pt8sn4-1-best', 4, 8, 4, 'Pt8Sn4'),
                ('pt8sn5-1-best', 5, 8, 5, 'Pt8Sn5'),
                ('pt8sn6-1-best', 6, 8, 6, 'Pt8Sn6'),
                ('pt8sn7-1-best', 7, 8, 7, 'Pt8Sn7'),
                ('pt8sn8-1-best', 8, 8, 8, 'Pt8Sn8'),
                ('pt8sn9-1-best', 9, 8, 9, 'Pt8Sn9'),
                ('pt8sn10-2-best', 10, 8, 10, 'Pt8Sn10')
            ],
            'output_subdir': 'step7.5.unified',  # 统一输出目录
            'exclude_from_heatmap': []  # 默认不屏蔽,使用--filter参数自定义
        },
        'PtxSn8-x': {
            'name': 'PtxSn8-x',
            'display_name': 'PtxSn8-x系列(总原子数=8)',
            'run_folders': {
                'Pt8': ['Pt8'],  # Pt8系列
                'PtxSn8-x': ['PtxSn8-x']  # PtxSn8-x系列
            },
            'systems': [
                ('pt8sn0-2-best', 0, 8, 0, 'Pt8', 'Pt8'),
                ('pt7sn1-1', 1, 7, 1, 'Pt7Sn1', 'PtxSn8-x'),
                ('pt6sn2', 2, 6, 2, 'Pt6Sn2', 'PtxSn8-x'),  # 修正文件夹名
                ('pt5sn3-1-best', 3, 5, 3, 'Pt5Sn3', 'PtxSn8-x'),  # 修正文件夹名
                ('pt4sn4-2', 4, 4, 4, 'Pt4Sn4', 'PtxSn8-x'),  # 修正文件夹名
                ('pt3sn5', 5, 3, 5, 'Pt3Sn5', 'PtxSn8-x')  # 修正文件夹名
            ],
            'output_subdir': 'step7.5.unified'  # 统一输出目录
        },
        'Pt6Snx': {
            'name': 'Pt6Snx',
            'display_name': 'Pt6+Sn系列',
            'run_folder': 'Pt6',
            'systems': [
                ('pt6', 0, 6, 0, 'Pt6'),
                ('pt6sn1', 1, 6, 1, 'Pt6Sn1'),
                ('pt6sn2', 2, 6, 2, 'Pt6Sn2'),
                ('pt6sn3', 3, 6, 3, 'Pt6Sn3'),
                ('pt6sn4', 4, 6, 4, 'Pt6Sn4'),
                ('pt6sn5-2', 5, 6, 5, 'Pt6Sn5'),
                ('pt6sn6-2', 6, 6, 6, 'Pt6Sn6'),
                ('pt6sn7', 7, 6, 7, 'Pt6Sn7'),
                ('pt6sn8', 8, 6, 8, 'Pt6Sn8'),
                ('pt6sn9-2', 9, 6, 9, 'Pt6Sn9')
            ],
            'output_subdir': 'step7.5.unified',  # 统一输出目录
            'exclude_from_heatmap': []  # 不屏蔽
        }
    }
    
    def __init__(self, base_path, output_base_dir, series_name, enable_msd_filter=False):
        """
        初始化分析器
        
        Args:
            base_path: v626数据根目录
            output_base_dir: 输出基础目录
            series_name: 系列名称 ('Pt8Snx', 'PtxSn8-x', 'Pt6Snx')
            enable_msd_filter: 是否启用MSD异常值筛选（默认False）
        """
        self.base_path = Path(base_path)
        self.output_base_dir = Path(output_base_dir)
        self.series_name = series_name
        self.enable_msd_filter = enable_msd_filter
        
        # 获取系列配置
        if series_name not in self.SERIES_CONFIGS:
            raise ValueError(f"未知系列: {series_name}, 可选: {list(self.SERIES_CONFIGS.keys())}")
        
        self.config = self.SERIES_CONFIGS[series_name]
        
        # 创建输出目录
        self.output_dir = self.output_base_dir / self.config['output_subdir']
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # 初始化数据定位器
        self.locator = V625DataLocator(base_path)
        
        # 加载系列配置
        self.systems = self.config['systems']
        self.temperatures = ['200K', '300K', '400K', '500K', '600K', 
                            '700K', '800K', '900K', '1000K', '1100K']
        self.results = {}
        
        # MSD过滤统计
        self.msd_filter_stats = {
            'total_runs': 0,
            'filtered_runs': 0,
            'filtered_data_points': 0
        }
        
        # 加载MSD异常路径签名（如果启用）
        if enable_msd_filter:
            self.msd_outliers = self.load_msd_outliers()
        else:
            self.msd_outliers = set()
        
        print(f"\n{'='*100}")
        print(f"初始化分析器: {self.config['display_name']}")
        print(f"输出目录: {self.output_dir}")
        if self.enable_msd_filter:
            print(f"✅ MSD过滤: 已启用 ({len(self.msd_outliers)}个异常路径签名)")
        else:
            print(f"⚠️ MSD过滤: 未启用，使用全部数据")
        print(f"{'='*100}")
    
    def load_msd_outliers(self):
        """
        加载Step 1的MSD异常路径签名
        
        Returns:
            set: 异常路径签名集合
        """
        # 使用绝对路径（相对于scripts目录的上级）
        script_dir = Path(__file__).parent
        outliers_file = script_dir.parent / 'results' / 'large_D_outliers.csv'
        
        if not outliers_file.exists():
            print(f"\n⚠️ 未找到MSD异常文件: {outliers_file}")
            print(f"   将使用全部数据，不进行筛选")
            return set()
        
        try:
            df_outliers = pd.read_csv(outliers_file, encoding='utf-8')
            print(f"\n>>> 加载MSD异常数据:")
            print(f"    异常记录数: {len(df_outliers)}")
            
            # 统计异常原因
            if 'reason' in df_outliers.columns:
                for reason, count in df_outliers['reason'].value_counts().items():
                    pct = count / len(df_outliers) * 100
                    print(f"      - {reason}: {count} ({pct:.1f}%)")
            
            # 提取路径签名
            filter_signatures = set()
            for _, row in df_outliers.iterrows():
                filepath = row.get('filepath', '')
                if filepath:
                    sig = extract_path_signature(filepath, is_msd_path=True)
                    if sig:
                        filter_signatures.add(sig)
            
            print(f"    唯一路径签名: {len(filter_signatures)}")
            return filter_signatures
            
        except Exception as e:
            print(f"\n⚠️ 加载MSD异常文件出错: {e}")
            print(f"   将使用全部数据，不进行筛选")
            return set()
    
    def load_coordination_time_series(self, run_path, sys_name, temp):
        """
        加载coordination时间序列数据
        支持v625格式(300K)和v626格式(T300.r3.gpu0)
        v626会自动加载该温度的所有重复运行并返回列表
        
        Returns:
            list[DataFrame] 或 None: 返回所有运行的DataFrame列表
        """
        sys_path = run_path / sys_name
        if not sys_path.exists():
            return None
        
        # 尝试v625格式: 300K (单次运行)
        csv_path_v625 = sys_path / temp / 'coordination_time_series.csv'
        if csv_path_v625.exists():
            try:
                df = pd.read_csv(csv_path_v625)
                return [df]  # 返回列表格式以统一处理
            except Exception as e:
                return None
        
        # 尝试v626格式: T300.r*.gpu* (多次运行)
        temp_value = temp.replace('K', '')
        temp_pattern = f"T{temp_value}.*"
        
        # 查找匹配的所有目录
        matching_dirs = sorted(sys_path.glob(temp_pattern))
        
        if matching_dirs:
            # 加载所有运行的数据
            all_dfs = []
            for temp_dir in matching_dirs:
                csv_path = temp_dir / 'coordination_time_series.csv'
                if csv_path.exists():
                    try:
                        df = pd.read_csv(csv_path)
                        all_dfs.append(df)
                    except Exception as e:
                        pass
            
            if all_dfs:
                return all_dfs  # 返回所有运行的列表
        
        return None
    
    def load_q6_time_series(self, run_path, sys_name, temp):
        """
        加载Q6时间序列数据
        支持v625格式(300K)和v626格式(T300.r3.gpu0)
        v626会自动加载该温度的所有重复运行并返回列表
        
        Returns:
            list[DataFrame] 或 None: 返回所有运行的DataFrame列表
        """
        sys_path = run_path / sys_name
        if not sys_path.exists():
            return None
        
        # 尝试v625格式: 300K (单次运行)
        csv_path_v625 = sys_path / temp / 'cluster_global_q6_time_series.csv'
        if csv_path_v625.exists():
            try:
                df = pd.read_csv(csv_path_v625)
                if 'cluster_metal_q6_global' not in df.columns:
                    return None
                return [df]  # 返回列表格式以统一处理
            except Exception as e:
                return None
        
        # 尝试v626格式: T300.r*.gpu* (多次运行)
        temp_value = temp.replace('K', '')
        temp_pattern = f"T{temp_value}.*"
        
        # 查找匹配的所有目录
        matching_dirs = sorted(sys_path.glob(temp_pattern))
        
        if matching_dirs:
            # 加载所有运行的数据
            all_dfs = []
            for temp_dir in matching_dirs:
                csv_path = temp_dir / 'cluster_global_q6_time_series.csv'
                if csv_path.exists():
                    try:
                        df = pd.read_csv(csv_path)
                        if 'cluster_metal_q6_global' not in df.columns:
                            continue
                        all_dfs.append(df)
                    except Exception as e:
                        pass
            
            if all_dfs:
                return all_dfs  # 返回所有运行的列表
        
        return None
    
    def load_geometry_time_series(self, run_path, sys_name, temp):
        """
        加载几何数据(回转半径、质心距离等)
        支持v625格式(300K)和v626格式(T300.r3.gpu0)
        v626会自动加载该温度的所有重复运行并返回列表
        
        Returns:
            list[DataFrame] 或 None: 返回所有运行的DataFrame列表
        """
        sys_path = run_path / sys_name
        if not sys_path.exists():
            return None
        
        # 尝试v625格式: 300K (单次运行)
        csv_path_v625 = sys_path / temp / 'cluster_geometry_time_series.csv'
        if csv_path_v625.exists():
            try:
                df = pd.read_csv(csv_path_v625)
                return [df]  # 返回列表格式以统一处理
            except Exception as e:
                return None
        
        # 尝试v626格式: T300.r*.gpu* (多次运行)
        temp_value = temp.replace('K', '')
        temp_pattern = f"T{temp_value}.*"
        
        # 查找匹配的所有目录
        matching_dirs = sorted(sys_path.glob(temp_pattern))
        
        if matching_dirs:
            # 加载所有运行的数据
            all_dfs = []
            for temp_dir in matching_dirs:
                csv_path = temp_dir / 'cluster_geometry_time_series.csv'
                if csv_path.exists():
                    try:
                        df = pd.read_csv(csv_path)
                        all_dfs.append(df)
                    except Exception as e:
                        pass
            
            if all_dfs:
                return all_dfs  # 返回所有运行的列表
        
        return None
    
    def load_element_comparison(self, run_path, sys_name, temp):
        """
        加载元素对比数据(Q4等)
        支持v625格式(300K)和v626格式(T300.r3.gpu0)
        v626会自动加载该温度的所有重复运行并返回列表
        
        Returns:
            list[DataFrame] 或 None: 返回所有运行的DataFrame列表
        """
        sys_path = run_path / sys_name
        if not sys_path.exists():
            return None
        
        # 尝试v625格式: 300K (单次运行)
        csv_path_v625 = sys_path / temp / 'element_comparison.csv'
        if csv_path_v625.exists():
            try:
                df = pd.read_csv(csv_path_v625)
                return [df]  # 返回列表格式以统一处理
            except Exception as e:
                return None
        
        # 尝试v626格式: T300.r*.gpu* (多次运行)
        temp_value = temp.replace('K', '')
        temp_pattern = f"T{temp_value}.*"
        
        # 查找匹配的所有目录
        matching_dirs = sorted(sys_path.glob(temp_pattern))
        
        if matching_dirs:
            # 加载所有运行的数据
            all_dfs = []
            for temp_dir in matching_dirs:
                csv_path = temp_dir / 'element_comparison.csv'
                if csv_path.exists():
                    try:
                        df = pd.read_csv(csv_path)
                        all_dfs.append(df)
                    except Exception as e:
                        pass
            
            if all_dfs:
                return all_dfs  # 返回所有运行的列表
        
        return None
    
    def _get_run_paths_for_system(self, system_info):
        """获取特定体系的运行路径"""
        if self.series_name == 'PtxSn8-x':
            # PtxSn8-x系列需要区分Pt8和PtxSn8-x
            folder_key = system_info[5]  # 'Pt8' 或 'PtxSn8-x'
            return self.locator.find_all_runs(folder_key)
        else:
            # Pt8Snx和Pt6Snx系列
            return self.locator.find_all_runs(self.config['run_folder'])
    
    def load_system_data(self, system_info):
        """
        加载单个体系在特定温度下的数据
        
        Args:
            system_info: (sys_name, sn_num, n_pt, n_sn, display_name, [folder_key])
        
        Returns:
            dict: {temp: data_dict}
        """
        sys_name = system_info[0]
        sn_num = system_info[1]
        n_pt = system_info[2]
        n_sn = system_info[3]
        display_name = system_info[4]
        
        # 获取该体系的运行路径
        run_paths = self._get_run_paths_for_system(system_info)
        
        print(f"\n处理 {display_name} ({sys_name}, {n_pt}Pt+{n_sn}Sn={n_pt+n_sn}原子)...")
        print(f"  使用{len(run_paths)}个运行文件夹")
        
        temp_results = {}
        
        for temp in self.temperatures:
            run_data_list = []
            
            # 读取所有运行的数据
            for run_path in run_paths:
                try:
                    # 加载coordination数据 - 返回列表
                    cn_dfs_list = self.load_coordination_time_series(run_path, sys_name, temp)
                    if cn_dfs_list is None:
                        continue
                    
                    # 加载Q6数据 - 返回列表
                    q6_dfs_list = self.load_q6_time_series(run_path, sys_name, temp)
                    if q6_dfs_list is None:
                        continue
                    
                    # 加载几何数据 - 返回列表
                    geo_dfs_list = self.load_geometry_time_series(run_path, sys_name, temp)
                    
                    # 加载元素对比数据(Q4等) - 返回列表
                    elem_dfs_list = self.load_element_comparison(run_path, sys_name, temp)
                    
                    # 获取温度目录列表用于路径签名检查
                    sys_path = run_path / sys_name
                    temp_value = temp.replace('K', '')
                    temp_dirs = sorted(sys_path.glob(f"T{temp_value}.*"))
                    
                    # 处理每个重复运行
                    num_runs = len(cn_dfs_list)
                    for idx in range(num_runs):
                        # MSD筛选: 检查当前运行的路径签名
                        if self.enable_msd_filter and idx < len(temp_dirs):
                            current_dir = temp_dirs[idx]
                            path_signature = extract_path_signature(str(current_dir), is_msd_path=False)
                            
                            if path_signature and path_signature in self.msd_outliers:
                                self.msd_filter_stats['filtered_runs'] += 1
                                continue  # 跳过该次运行
                        
                        df_cn = cn_dfs_list[idx]
                        df_q6 = q6_dfs_list[idx]
                        
                        # 计算统计量
                        pt_cn_total = df_cn['Pt_cn_total'].mean()
                        pt_pt_bonds = df_cn['Pt_cn_Pt_Pt'].mean() if 'Pt_cn_Pt_Pt' in df_cn.columns else 0
                        pt_sn_bonds = df_cn['Pt_cn_Pt_Sn'].mean() if 'Pt_cn_Pt_Sn' in df_cn.columns and n_sn > 0 else 0
                        
                        sn_cn_total = df_cn['Sn_cn_total'].mean() if 'Sn_cn_total' in df_cn.columns and n_sn > 0 else 0
                        sn_sn_bonds = df_cn['Sn_cn_Sn_Sn'].mean() if 'Sn_cn_Sn_Sn' in df_cn.columns and n_sn > 0 else 0
                        sn_pt_bonds = df_cn['Sn_cn_Sn_Pt'].mean() if 'Sn_cn_Sn_Pt' in df_cn.columns and n_sn > 0 else 0
                        
                        # 归一化的键密度
                        pt_pt_bonds_per_pt = pt_pt_bonds / n_pt if n_pt > 0 else 0
                        pt_sn_bonds_per_pt = pt_sn_bonds / n_pt if n_pt > 0 else 0
                        sn_sn_bonds_per_sn = sn_sn_bonds / n_sn if n_sn > 0 else 0
                        sn_pt_bonds_per_sn = sn_pt_bonds / n_sn if n_sn > 0 else 0
                        
                        run_data = {
                            'q6': df_q6['cluster_metal_q6_global'].mean(),
                            'q6_std': df_q6['cluster_metal_q6_global'].std(),
                            'pt_q6': df_q6['Pt_q6'].mean() if 'Pt_q6' in df_q6.columns else 0,
                            'sn_q6': df_q6['Sn_q6'].mean() if 'Sn_q6' in df_q6.columns and n_sn > 0 else 0,
                            'pt_cn_total': pt_cn_total,
                            'pt_pt_bonds': pt_pt_bonds,
                            'pt_sn_bonds': pt_sn_bonds,
                            'pt_pt_bonds_per_pt': pt_pt_bonds_per_pt,
                            'pt_sn_bonds_per_pt': pt_sn_bonds_per_pt,
                            'sn_cn_total': sn_cn_total,
                            'sn_sn_bonds': sn_sn_bonds,
                            'sn_pt_bonds': sn_pt_bonds,
                            'sn_sn_bonds_per_sn': sn_sn_bonds_per_sn,
                            'sn_pt_bonds_per_sn': sn_pt_bonds_per_sn,
                            'rg': 0, 'rg_std': 0,
                            'pt_dist': 0, 'sn_dist': 0, 'd_sn_pt': 0,
                            'pt_q4': 0, 'sn_q4': 0, 'cluster_q4': 0,
                            'n_pt': n_pt, 'n_sn': n_sn,
                            'total_atoms': n_pt + n_sn
                        }
                        
                        # 读取几何数据
                        if geo_dfs_list and idx < len(geo_dfs_list):
                            df_geo = geo_dfs_list[idx]
                            try:
                                run_data['rg'] = df_geo['gyration_radius'].mean()
                                run_data['rg_std'] = df_geo['gyration_radius'].std()
                                run_data['pt_dist'] = df_geo['pt_avg_dist_to_center'].mean()
                                if 'sn_avg_dist_to_center' in df_geo.columns and n_sn > 0:
                                    run_data['sn_dist'] = df_geo['sn_avg_dist_to_center'].mean()
                                    run_data['d_sn_pt'] = run_data['sn_dist'] - run_data['pt_dist']
                            except:
                                pass
                        
                        # 读取Q4数据
                        if elem_dfs_list and idx < len(elem_dfs_list):
                            df_elem = elem_dfs_list[idx]
                            try:
                                # element_comparison.csv格式: Element,CN_total,wGCN,Sn_wGCN,Q6,Q4
                                pt_row = df_elem[df_elem['Element'] == 'Pt']
                                if not pt_row.empty and 'Q4' in df_elem.columns:
                                    run_data['pt_q4'] = pt_row['Q4'].values[0]
                                
                                sn_row = df_elem[df_elem['Element'] == 'Sn']
                                if not sn_row.empty and 'Q4' in df_elem.columns and n_sn > 0:
                                    run_data['sn_q4'] = sn_row['Q4'].values[0]
                                
                                # 团簇整体Q4 (加权平均)
                                if n_sn > 0:
                                    run_data['cluster_q4'] = (run_data['pt_q4'] * n_pt + run_data['sn_q4'] * n_sn) / (n_pt + n_sn)
                                else:
                                    run_data['cluster_q4'] = run_data['pt_q4']
                            except:
                                pass
                        
                        run_data_list.append(run_data)
                    
                except Exception as e:
                    continue
            
            # 如果有成功的运行,计算平均值
            if run_data_list:
                avg_data = {}
                for key in run_data_list[0].keys():
                    values = [rd[key] for rd in run_data_list]
                    if isinstance(values[0], (int, float, np.number)):
                        avg_data[key] = np.mean(values)
                    else:
                        avg_data[key] = values[0]
                
                avg_data['n_runs'] = len(run_data_list)
                temp_results[temp] = avg_data
                
                # 输出信息
                info_parts = [
                    f"Q6={avg_data['q6']:.3f}",
                    f"Pt-Pt={avg_data['pt_pt_bonds_per_pt']:.3f}",
                    f"Pt-Sn={avg_data['pt_sn_bonds_per_pt']:.3f}"
                ]
                
                # 如果有几何数据，添加Rg和Δd
                if avg_data['rg'] > 0:
                    info_parts.append(f"Rg={avg_data['rg']:.3f}")
                    if n_sn > 0 and avg_data['d_sn_pt'] != 0:
                        info_parts.append(f"Δd={avg_data['d_sn_pt']:.3f}")
                
                # 如果有Q4数据，添加Q4
                if avg_data['cluster_q4'] > 0:
                    info_parts.append(f"Q4={avg_data['cluster_q4']:.3f}")
                
                info_parts.append(f"(平均{len(run_data_list)}次运行)")
                print(f"  [OK] {temp}: {', '.join(info_parts)}")
        
        return temp_results
    
    def collect_all_data(self):
        """收集所有体系所有温度的数据"""
        print("\n" + "="*100)
        print(f"收集{self.config['display_name']}数据 - 多运行平均")
        print("="*100)
        
        total_expected = len(self.systems) * len(self.temperatures)
        total_success = 0
        total_runs_used = 0
        run_distribution = {}
        
        for system_info in self.systems:
            sn_num = system_info[1]
            temp_results = self.load_system_data(system_info)
            
            if temp_results:
                self.results[sn_num] = temp_results
                for temp, data in temp_results.items():
                    total_success += 1
                    n_runs = data.get('n_runs', 1)
                    total_runs_used += n_runs
                    
                    if n_runs not in run_distribution:
                        run_distribution[n_runs] = 0
                    run_distribution[n_runs] += 1
        
        print("\n" + "="*100)
        print("数据收集完成!")
        print("="*100)
        print(f"\n数据统计:")
        print(f"  预期数据点: {total_expected} ({len(self.systems)}体系 × {len(self.temperatures)}温度)")
        print(f"  成功读取: {total_success} ({total_success/total_expected*100:.1f}%)")
        print(f"  失败: {total_expected - total_success}")
        print(f"  总运行次数: {total_runs_used}")
        if total_success > 0:
            print(f"  平均每点运行数: {total_runs_used/total_success:.2f}")
        
        if run_distribution:
            print(f"\n运行数分布:")
            for n_runs in sorted(run_distribution.keys()):
                count = run_distribution[n_runs]
                pct = count / total_success * 100
                print(f"    {n_runs}次运行: {count}个数据点 ({pct:.1f}%)")
        
        # 显示MSD筛选统计
        if self.enable_msd_filter and self.msd_filter_stats['filtered_runs'] > 0:
            print(f"\n[MSD筛选统计]:")
            total_attempted = total_runs_used + self.msd_filter_stats['filtered_runs']
            filter_rate = self.msd_filter_stats['filtered_runs'] / total_attempted * 100 if total_attempted > 0 else 0
            print(f"  尝试读取运行: {total_attempted}")
            print(f"  保留运行: {total_runs_used} ({100-filter_rate:.1f}%)")
            print(f"  过滤运行: {self.msd_filter_stats['filtered_runs']} ({filter_rate:.1f}%)")
        
        print("="*100)
    
    def save_data_table(self):
        """保存数据表"""
        rows = []
        for sn_num in sorted(self.results.keys()):
            for temp in self.temperatures:
                if temp in self.results[sn_num]:
                    data = self.results[sn_num][temp]
                    row = {
                        'series': self.series_name,
                        'sn_num': sn_num,
                        'temperature': temp,
                        'temp_value': int(temp.replace('K', '')),
                        'n_pt': data['n_pt'],
                        'n_sn': data['n_sn'],
                        'total_atoms': data['total_atoms'],
                        'q6': data['q6'],
                        'q6_std': data['q6_std'],
                        'pt_q6': data.get('pt_q6', 0),
                        'sn_q6': data.get('sn_q6', 0),
                        'pt_cn_total': data['pt_cn_total'],
                        'pt_pt_bonds_per_pt': data['pt_pt_bonds_per_pt'],
                        'pt_sn_bonds_per_pt': data['pt_sn_bonds_per_pt'],
                        'sn_cn_total': data.get('sn_cn_total', 0),
                        'sn_sn_bonds_per_sn': data.get('sn_sn_bonds_per_sn', 0),
                        'sn_pt_bonds_per_sn': data.get('sn_pt_bonds_per_sn', 0),
                        'rg': data.get('rg', 0),
                        'd_sn_pt': data.get('d_sn_pt', 0),
                        'n_runs': data.get('n_runs', 1)
                    }
                    rows.append(row)
        
        df = pd.DataFrame(rows)
        csv_file = self.output_dir / f"{self.series_name.lower()}_multi_temp_data.csv"
        df.to_csv(csv_file, index=False, encoding='utf-8-sig')
        print(f"\n💾 数据表已保存: {csv_file}")
    
    def plot_comprehensive_analysis(self):
        """绘制综合分析图(6行×10列),支持过滤指定Sn含量"""
        fig = plt.figure(figsize=(24, 18))
        gs = fig.add_gridspec(6, 10, hspace=0.35, wspace=0.3)
        
        # 获取过滤配置
        exclude_sn = self.config.get('exclude_from_heatmap', [])
        
        # 为每个温度绘制6个子图
        for temp_idx, temp in enumerate(self.temperatures):
            sn_nums = []
            q6_vals = []
            pt_pt_vals = []
            pt_sn_vals = []
            rg_vals = []
            d_sn_pt_vals = []
            q4_vals = []
            
            for sn_num in sorted(self.results.keys()):
                # 应用过滤
                if sn_num in exclude_sn:
                    continue
                    
                if temp in self.results[sn_num]:
                    data = self.results[sn_num][temp]
                    sn_nums.append(sn_num)
                    q6_vals.append(data['q6'])
                    pt_pt_vals.append(data['pt_pt_bonds_per_pt'])
                    pt_sn_vals.append(data['pt_sn_bonds_per_pt'])
                    rg_vals.append(data.get('rg', 0))
                    d_sn_pt_vals.append(data.get('d_sn_pt', 0))
                    q4_vals.append(data.get('cluster_q4', 0))
            
            if not sn_nums:
                continue
            
            # Row 1: Q6 vs Sn含量
            ax1 = fig.add_subplot(gs[0, temp_idx])
            ax1.plot(sn_nums, q6_vals, 'o-', linewidth=2, markersize=6, color='blue')
            ax1.set_xlabel('Sn含量', fontsize=8)
            ax1.set_ylabel('Q6', fontsize=8)
            ax1.set_title(f'{temp}', fontsize=9, fontweight='bold')
            ax1.grid(True, alpha=0.3)
            ax1.tick_params(labelsize=7)
            
            # Row 2: Pt-Pt键 vs Sn含量
            ax2 = fig.add_subplot(gs[1, temp_idx])
            ax2.plot(sn_nums, pt_pt_vals, 's-', linewidth=2, markersize=6, color='green')
            ax2.set_xlabel('Sn含量', fontsize=8)
            ax2.set_ylabel('Pt-Pt键/Pt', fontsize=8)
            ax2.grid(True, alpha=0.3)
            ax2.tick_params(labelsize=7)
            
            # Row 3: Pt-Sn键 vs Sn含量
            ax3 = fig.add_subplot(gs[2, temp_idx])
            ax3.plot(sn_nums, pt_sn_vals, '^-', linewidth=2, markersize=6, color='red')
            ax3.set_xlabel('Sn含量', fontsize=8)
            ax3.set_ylabel('Pt-Sn键/Pt', fontsize=8)
            ax3.grid(True, alpha=0.3)
            ax3.tick_params(labelsize=7)
            
            # Row 4: 回转半径Rg vs Sn含量
            ax4 = fig.add_subplot(gs[3, temp_idx])
            if any(v > 0 for v in rg_vals):
                ax4.plot(sn_nums, rg_vals, 'd-', linewidth=2, markersize=6, color='purple')
                ax4.set_xlabel('Sn含量', fontsize=8)
                ax4.set_ylabel('Rg (Å)', fontsize=8)
                ax4.grid(True, alpha=0.3)
                ax4.tick_params(labelsize=7)
            else:
                ax4.text(0.5, 0.5, 'No Rg data', ha='center', va='center', transform=ax4.transAxes)
                ax4.axis('off')
            
            # Row 5: 距离差Δd vs Sn含量 (核壳指标)
            ax5 = fig.add_subplot(gs[4, temp_idx])
            if any(v != 0 for v in d_sn_pt_vals):
                ax5.plot(sn_nums, d_sn_pt_vals, 'v-', linewidth=2, markersize=6, color='orange')
                ax5.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
                ax5.set_xlabel('Sn含量', fontsize=8)
                ax5.set_ylabel('Δd (Sn-Pt) (Å)', fontsize=8)
                ax5.grid(True, alpha=0.3)
                ax5.tick_params(labelsize=7)
            else:
                ax5.text(0.5, 0.5, 'No Δd data', ha='center', va='center', transform=ax5.transAxes)
                ax5.axis('off')
            
            # Row 6: Q4四次对称性 vs Sn含量
            ax6 = fig.add_subplot(gs[5, temp_idx])
            if any(v > 0 for v in q4_vals):
                ax6.plot(sn_nums, q4_vals, 'h-', linewidth=2, markersize=6, color='brown')
                ax6.set_xlabel('Sn含量', fontsize=8)
                ax6.set_ylabel('Q4', fontsize=8)
                ax6.grid(True, alpha=0.3)
                ax6.tick_params(labelsize=7)
            else:
                ax6.text(0.5, 0.5, 'No Q4 data', ha='center', va='center', transform=ax6.transAxes)
                ax6.axis('off')
        
        plt.suptitle(f'{self.config["display_name"]} - 多温度综合分析 (v626完整版)', 
                     fontsize=16, fontweight='bold', y=0.995)
        
        output_file = self.output_dir / f'{self.series_name.lower()}_comprehensive_analysis.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"📊 综合分析图已保存: {output_file}")
    
    def plot_heatmaps(self):
        """绘制热图,支持屏蔽指定Sn含量"""
        fig, axes = plt.subplots(3, 3, figsize=(22, 16))
        
        # 准备数据 - 过滤掉需要屏蔽的Sn含量
        exclude_sn = self.config.get('exclude_from_heatmap', [])
        sn_nums = [sn for sn in sorted(self.results.keys()) if sn not in exclude_sn]
        temps = self.temperatures
        
        if exclude_sn:
            print(f"\n[热图] 屏蔽Sn含量: {exclude_sn}")
        
        def create_matrix(field):
            matrix = []
            for temp in temps:
                row = []
                for sn_num in sn_nums:
                    if temp in self.results[sn_num]:
                        value = self.results[sn_num][temp].get(field, 0)
                        # 对于新增字段，如果为0可能表示没有数据
                        if field in ['rg', 'd_sn_pt', 'cluster_q4'] and value == 0:
                            row.append(np.nan)
                        else:
                            row.append(value)
                    else:
                        row.append(np.nan)
                matrix.append(row)
            return np.array(matrix)
        
        # 9个热图
        fields = [
            ('q6', 'Q6 六次对称性'),
            ('pt_pt_bonds_per_pt', 'Pt-Pt键/Pt'),
            ('pt_sn_bonds_per_pt', 'Pt-Sn键/Pt'),
            ('rg', '回转半径 Rg (Å)'),
            ('d_sn_pt', '距离差 Δd(Sn-Pt) (Å)'),
            ('cluster_q4', 'Q4 四次对称性'),
            ('sn_sn_bonds_per_sn', 'Sn-Sn键/Sn'),
            ('pt_dist', 'Pt到质心距离 (Å)'),
            ('sn_dist', 'Sn到质心距离 (Å)')
        ]
        
        for idx, (field, title) in enumerate(fields):
            ax = axes[idx // 3, idx % 3]
            matrix = create_matrix(field)
            
            # 检查是否有有效数据
            if np.all(np.isnan(matrix)):
                ax.text(0.5, 0.5, f'No {field} data', 
                       ha='center', va='center', fontsize=12, transform=ax.transAxes)
                ax.set_title(title, fontsize=11, fontweight='bold')
                ax.axis('off')
                continue
            
            # 特殊处理Δd的colormap (中心为0)
            if field == 'd_sn_pt':
                vmax = np.nanmax(np.abs(matrix))
                im = ax.imshow(matrix, aspect='auto', cmap='RdBu_r', 
                              interpolation='nearest', vmin=-vmax, vmax=vmax)
            else:
                im = ax.imshow(matrix, aspect='auto', cmap='RdYlBu_r', interpolation='nearest')
            
            ax.set_xticks(range(len(sn_nums)))
            ax.set_xticklabels([f'Sn{s}' for s in sn_nums], fontsize=8)
            ax.set_yticks(range(len(temps)))
            ax.set_yticklabels(temps, fontsize=8)
            ax.set_xlabel('Sn含量', fontsize=9)
            ax.set_ylabel('温度', fontsize=9)
            ax.set_title(title, fontsize=10, fontweight='bold')
            
            # 添加数值标注
            for i in range(len(temps)):
                for j in range(len(sn_nums)):
                    if not np.isnan(matrix[i, j]):
                        text = ax.text(j, i, f'{matrix[i, j]:.2f}',
                                     ha="center", va="center", color="black", fontsize=6)
            
            plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        
        plt.suptitle(f'{self.config["display_name"]} - 热图分析 (v626完整版)', 
                    fontsize=14, fontweight='bold')
        plt.tight_layout()
        
        output_file = self.output_dir / f'{self.series_name.lower()}_heatmap.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"📊 热图已保存: {output_file}")
    
    def run_analysis(self):
        """运行完整分析流程"""
        print(f"\n{'='*100}")
        print(f"开始{self.config['display_name']}完整分析")
        print(f"{'='*100}")
        
        # 1. 收集数据
        self.collect_all_data()
        
        # 2. 保存数据表
        self.save_data_table()
        
        # 3. 绘制图表
        self.plot_comprehensive_analysis()
        self.plot_heatmaps()
        
        print(f"\n{'='*100}")
        print(f"✅ {self.config['display_name']}分析完成!")
        print(f"{'='*100}")


def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='统一多温度分析脚本 - v625数据')
    parser.add_argument('--series', type=str, choices=['Pt8Snx', 'PtxSn8-x', 'Pt6Snx'],
                       help='指定要分析的系列')
    parser.add_argument('--all', action='store_true',
                       help='运行所有系列的分析')
    parser.add_argument('--filter', type=str, default='',
                       help='自定义过滤Sn含量(逗号分隔),例如: --filter 0,1,2 屏蔽Sn0/Sn1/Sn2')
    parser.add_argument('--enable-msd-filter', action='store_true',
                       help='启用Step 1的MSD异常值筛选(移除数据质量差的运行)')
    
    args = parser.parse_args()
    
    # 路径配置 (使用v626数据)
    base_path = Path(__file__).parent / 'data' / 'coordination' / 'coordination_time_series_results_sample_20251106_214943'
    output_base_dir = Path(__file__).parent / 'results'
    
    # 确定要运行的系列
    if args.all:
        series_list = ['Pt8Snx', 'PtxSn8-x', 'Pt6Snx']
    elif args.series:
        series_list = [args.series]
    else:
        print("错误: 请指定 --series 或 --all")
        print("用法示例:")
        print("  python step7-5-unified_multi_temp_v626_analysis.py --series Pt8Snx")
        print("  python step7-5-unified_multi_temp_v626_analysis.py --all")
        print("  python step7-5-unified_multi_temp_v626_analysis.py --series Pt8Snx --enable-msd-filter")
        return
    
    # 运行分析
    for series in series_list:
        try:
            analyzer = UnifiedMultiTempAnalyzer(
                base_path=base_path,
                output_base_dir=output_base_dir,
                series_name=series,
                enable_msd_filter=args.enable_msd_filter
            )
            
            # 应用自定义过滤
            if args.filter:
                filter_sn = [int(x.strip()) for x in args.filter.split(',')]
                analyzer.config['exclude_from_heatmap'] = filter_sn
                print(f"\n[自定义过滤] 屏蔽Sn含量: {filter_sn}")
            
            analyzer.run_analysis()
                
        except Exception as e:
            print(f"\n❌ {series}分析失败: {e}")
            import traceback
            traceback.print_exc()
            continue


if __name__ == '__main__':
    main()
