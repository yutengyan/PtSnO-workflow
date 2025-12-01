#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 7.6 新服务器数据 - 温度对比分析 (完整整合版)

整合原有3个脚本的功能，同时支持 Q6 和 OP2 分析：

  7.6.1: 温度并排对比 (sidebyside)
         → 对比两个温度的 Q6/OP2 时间演化
         → 所有体系在同一张图上并排展示
         
  7.6.2: 单体系温度对比 (individual)
         → 每个体系单独一张图
         → 4行(cluster_q6, pt_q6, sn_q6, op2) × 2列(低温, 高温)
         → 统计盒子: 均值、标准差、变异系数(CV)
         
  7.6.3: 统计对比 (stats)
         → 柱状图对比不同体系的 Q6/OP2 均值
         → 变异系数(CV)对比 → 识别温度稳定性
         → 温度敏感性散点图
         
  额外: 统计热图 (heatmap) + Q6 vs OP2 相关性

数据结构:
dp-md/4090-ustc/more/Pt6/pt6sn8/T500.r3.gpu0/
dp-md/4090-ustc/o68/g-1535-Sn8Pt6O4/Cv-1/T500.r6.gpu0/

用法:
    # 运行所有分析
    python step7-6-new_server_temp_comparison.py
    
    # 模式1: 并排对比两个温度
    python step7-6-new_server_temp_comparison.py --mode sidebyside --temps 300,900

    # 模式2: 单个系统详细分析
    python step7-6-new_server_temp_comparison.py --mode individual --temps 300,900

    # 模式3: 统计对比 (柱状图+散点图)
    python step7-6-new_server_temp_comparison.py --mode stats --temps 300,900
    
    # 模式4: 统计热图
    python step7-6-new_server_temp_comparison.py --mode heatmap

Author: AI Assistant
Date: 2025-11-30
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from scipy.signal import savgol_filter
import os
import sys
import warnings
import argparse
import re
from pathlib import Path
import seaborn as sns

# 设置控制台输出编码
if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8')

warnings.filterwarnings('ignore')

# 配置中文字体
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'Arial']
plt.rcParams['axes.unicode_minus'] = False

# 中文字体属性
try:
    chinese_font = FontProperties(family='Microsoft YaHei', size=9)
except:
    chinese_font = FontProperties(size=9)


class NewServerTempComparisonAnalyzer:
    """新服务器数据温度对比分析器"""
    
    def __init__(self, data_root, output_dir=None, unified_ylim=True):
        """
        初始化分析器
        
        Args:
            data_root: 解压后的数据根目录
            output_dir: 输出目录
            unified_ylim: 是否统一Y轴范围（同一行的左右两列使用相同Y轴）
        """
        self.data_root = Path(data_root)
        if output_dir:
            self.output_dir = Path(output_dir)
        else:
            self.output_dir = Path(__file__).parent / 'results' / 'step7.6_new_server'
        
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # 存储扫描到的数据结构
        self.systems = {}
        
        # Y轴对齐选项
        self.unified_ylim = unified_ylim
        
        print(f"\n{'='*80}")
        print(f"Step 7.6 新服务器数据温度对比分析")
        print(f"数据目录: {self.data_root}")
        print(f"输出目录: {self.output_dir}")
        print(f"{'='*80}")
    
    def scan_data_structure(self):
        """扫描数据目录结构"""
        print("\n>>> 扫描数据结构...")
        
        csv_files = list(self.data_root.rglob('coordination_time_series.csv'))
        print(f"    找到 {len(csv_files)} 个数据点")
        
        for csv_path in csv_files:
            run_dir = csv_path.parent
            system_dir = run_dir.parent
            
            run_name = run_dir.name
            match = re.match(r'T(\d+)\.r(\d+)\.gpu\d+', run_name)
            if not match:
                continue
            
            temp = int(match.group(1))
            run_num = int(match.group(2))
            
            system_name = system_dir.name.lower()
            path_parts = csv_path.relative_to(self.data_root).parts
            
            if len(path_parts) >= 4:
                system_key = '/'.join(path_parts[:-2])
            else:
                system_key = system_name
            
            if system_key not in self.systems:
                self.systems[system_key] = {
                    'name': system_name,
                    'path': system_dir,
                    'temps': {},
                    'display_name': self._get_display_name(system_name, path_parts)
                }
            
            if temp not in self.systems[system_key]['temps']:
                self.systems[system_key]['temps'][temp] = []
            
            self.systems[system_key]['temps'][temp].append({
                'run_num': run_num,
                'run_dir': run_dir,
                'run_name': run_name
            })
        
        print(f"    发现 {len(self.systems)} 个系统")
        
        # 创建 ensemble 平均系统
        self._create_ensemble_systems()
        
        return self.systems
    
    def _create_ensemble_systems(self):
        """
        创建 ensemble 平均系统:
        - Pt6Sn8_ensemble: Pt6Sn8 + Pt6-2 + Pt6-3 的平均 (3个初始结构，纯金属团簇)
        - Cv_ensemble: Cv-1 ~ Cv-5 的平均 (5个初始结构，含氧化物)
        - Air68_ensemble: air/68 + air-2/68 的平均 (Pt6Sn8 气相)
        - Air86_ensemble: air/86 + air-2/86 的平均 (Pt8Sn6 气相)
        
        关键区分:
        - Pt6Sn8 系列路径: dp-md/4090-ustc/more/Pt6*/pt6sn8/
        - Cv 系列路径: dp-md/4090-ustc/o68/g-1535-Sn8Pt6O4/Cv-*/
        - Air68 系列路径: air*/68/
        - Air86 系列路径: air*/86/
        """
        # === 1. 识别 Pt6Sn8 系列 (纯金属，路径中有 /more/Pt6) ===
        pt6sn8_keys = [k for k in self.systems.keys() 
                       if '/more/Pt6' in k and 'Cv-' not in k]
        
        if len(pt6sn8_keys) >= 2:
            self._add_ensemble('Pt6Sn8', pt6sn8_keys, 'Pt6Sn8 系综平均')
        
        # === 2. 识别 Cv 系列 (Pt6Sn8O4，路径中有 /Cv-) ===
        cv_keys = [k for k in self.systems.keys() if '/Cv-' in k]
        if len(cv_keys) >= 2:
            self._add_ensemble('Cv', cv_keys, 'Pt6Sn8O4 系综平均')
        
        # === 3. 识别 Air-68 系列 (Pt6Sn8 气相，路径中有 /68) ===
        # 路径格式: air/68 或 air-2/68
        air68_keys = [k for k in self.systems.keys() 
                      if '/68' in k and ('air' in k.lower() or 'Air' in k)]
        if len(air68_keys) >= 2:
            self._add_ensemble('Air68', air68_keys, 'Air-68(Pt6Sn8) 系综平均')
        
        # === 4. 识别 Air-86 系列 (Pt8Sn6 气相，路径中有 /86) ===
        # 路径格式: air/86 或 air-2/86
        air86_keys = [k for k in self.systems.keys() 
                      if '/86' in k and ('air' in k.lower() or 'Air' in k)]
        if len(air86_keys) >= 2:
            self._add_ensemble('Air86', air86_keys, 'Air-86(Pt8Sn6) 系综平均')
    
    def _add_ensemble(self, name, sys_keys, display_prefix):
        """
        通用的系综创建辅助函数
        
        Args:
            name: 系综名称（如 'Pt6Sn8', 'Cv', 'Air68', 'Air86'）
            sys_keys: 属于该系综的系统键列表
            display_prefix: 显示名称前缀
        """
        # 合并所有变体
        ensemble_temps = {}
        for sys_key in sys_keys:
            for temp, runs in self.systems[sys_key]['temps'].items():
                if temp not in ensemble_temps:
                    ensemble_temps[temp] = []
                for run in runs:
                    run_with_source = run.copy()
                    run_with_source['source_system'] = sys_key
                    ensemble_temps[temp].append(run_with_source)
        
        # 只保留所有系统都有的温度
        common_temps = {}
        for temp, runs in ensemble_temps.items():
            if len(runs) >= len(sys_keys):  # 至少每个系统都有数据
                common_temps[temp] = runs
        
        if common_temps:
            ensemble_key = f'{name}_ensemble'
            self.systems[ensemble_key] = {
                'name': f'{name.lower()}_ensemble',
                'path': None,
                'temps': common_temps,
                'display_name': f'{display_prefix} (n={len(sys_keys)})',
                'is_ensemble': True,
                'source_systems': sys_keys
            }
            print(f"    ✅ 创建 {name} 系综: {len(sys_keys)} 个初始结构")
    
    def _get_display_name(self, system_name, path_parts):
        """生成显示名称"""
        # Cv 系列
        if system_name.startswith('cv-'):
            return f"Cv{system_name[-1]}"
        
        # Pt6Sn8 负载系列
        if 'pt6sn8' in system_name:
            for part in path_parts:
                if part.startswith('Pt6-'):
                    return f"Pt6Sn8({part})"
                if part == 'Pt6':
                    return "Pt6Sn8"
            return "Pt6Sn8"
        
        # Air 气相系列
        # 路径格式: air/68 或 air-2/86
        path_str = '/'.join(path_parts)
        if 'air' in path_str.lower():
            # 找到 air 或 air-2 部分
            air_part = None
            cluster_type = None
            for i, part in enumerate(path_parts):
                if part.lower().startswith('air'):
                    air_part = part
                if part in ['68', '86']:
                    cluster_type = part
            
            if cluster_type == '68':
                if air_part and air_part != 'air':
                    return f"Air-68({air_part})"  # Air-68(air-2)
                return "Air-68"
            elif cluster_type == '86':
                if air_part and air_part != 'air':
                    return f"Air-86({air_part})"  # Air-86(air-2)
                return "Air-86"
        
        return system_name.upper()
    
    def load_time_series(self, run_dir, data_type='q6'):
        """加载时间序列数据"""
        if data_type == 'q6':
            csv_file = run_dir / 'cluster_global_q6_time_series.csv'
        elif data_type == 'op2':
            csv_file = run_dir / 'cluster_op2_time_series.csv'
        elif data_type == 'coordination':
            csv_file = run_dir / 'coordination_time_series.csv'
        elif data_type == 'geometry':
            csv_file = run_dir / 'cluster_geometry_time_series.csv'
        else:
            return None
        
        if csv_file.exists():
            try:
                return pd.read_csv(csv_file)
            except:
                return None
        return None
    
    def calculate_statistics(self, values):
        """计算统计指标"""
        mean = np.mean(values)
        std = np.std(values)
        cv = std / mean if mean != 0 else 0
        
        return {
            'mean': mean,
            'std': std,
            'cv': cv,
            'min': np.min(values),
            'max': np.max(values),
        }
    
    def plot_side_by_side_comparison(self, temps):
        """
        模式1: 并排对比两个温度的Q6/OP2演化（所有系统在同一张图）
        """
        print(f"\n>>> 绘制并排温度对比: {temps[0]}K vs {temps[1]}K")
        
        # 创建大图: 4行 × 2列
        # 行1: cluster_metal_q6, 行2: pt_q6, 行3: sn_q6, 行4: op2_all_metal
        fig, axes = plt.subplots(4, 2, figsize=(16, 16))
        
        fields = [
            ('q6', 'cluster_metal_q6_global', 'Q6 整体金属团簇'),
            ('q6', 'pt_q6_global', 'Q6 Pt原子'),
            ('q6', 'sn_q6_global', 'Q6 Sn原子'),
            ('op2', 'op2_all_metal', 'OP2 取向参数'),
        ]
        
        colors = plt.cm.tab10(np.linspace(0, 1, len(self.systems)))
        
        for row_idx, (data_type, field, title) in enumerate(fields):
            for col_idx, temp in enumerate(temps):
                ax = axes[row_idx, col_idx]
                
                valid_count = 0
                for sys_idx, (sys_key, sys_info) in enumerate(self.systems.items()):
                    if temp not in sys_info['temps']:
                        continue
                    
                    # 加载该温度所有运行的数据
                    all_dfs = []
                    for run_info in sys_info['temps'][temp]:
                        df = self.load_time_series(run_info['run_dir'], data_type)
                        if df is not None and field in df.columns:
                            all_dfs.append(df)
                    
                    if not all_dfs:
                        continue
                    
                    valid_count += 1
                    
                    # 计算平均
                    if len(all_dfs) > 1:
                        avg_df = pd.concat(all_dfs).groupby('frame').mean().reset_index()
                    else:
                        avg_df = all_dfs[0]
                    
                    time_ps = avg_df['time_ps'].values
                    values = avg_df[field].values
                    
                    # 绘制
                    ax.plot(time_ps, values, linewidth=1.5, 
                           color=colors[sys_idx], 
                           label=sys_info['display_name'], alpha=0.8)
                
                # 设置标题和标签
                if row_idx == 0:
                    ax.set_title(f'{temp}K', fontsize=14, fontweight='bold')
                
                if col_idx == 0:
                    ax.set_ylabel(title, fontsize=11, fontweight='bold')
                
                if row_idx == 3:
                    ax.set_xlabel('时间 (ps)', fontsize=11)
                
                ax.grid(True, alpha=0.3)
                
                # 图例(只在右上角第一行显示)
                if row_idx == 0 and col_idx == 1:
                    ax.legend(loc='upper right', fontsize=8, ncol=2)
                
                ax.text(0.02, 0.98, f'n={valid_count}系统', 
                       transform=ax.transAxes, fontsize=9, 
                       verticalalignment='top',
                       bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        plt.suptitle(f'温度对比: {temps[0]}K vs {temps[1]}K (Q6 & OP2 时间演化)', 
                    fontsize=16, fontweight='bold', y=0.995)
        plt.tight_layout()
        
        output_file = self.output_dir / f'sidebyside_{temps[0]}K_vs_{temps[1]}K.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"📊 并排对比图已保存: {output_file}")
    
    def plot_individual_system_comparison(self, temps, system_key=None):
        """
        模式2: 单个系统的详细温度对比（每个系统一张图）
        """
        print(f"\n>>> 绘制单系统温度对比: {temps[0]}K vs {temps[1]}K")
        
        systems_to_plot = [system_key] if system_key else list(self.systems.keys())
        
        for sys_key in systems_to_plot:
            if sys_key not in self.systems:
                continue
            
            sys_info = self.systems[sys_key]
            
            # 检查两个温度是否都有数据
            if temps[0] not in sys_info['temps'] or temps[1] not in sys_info['temps']:
                print(f"  ⚠️ {sys_info['display_name']} 缺少温度数据，跳过")
                continue
            
            self._plot_single_system(sys_info, temps)
    
    def _plot_single_system(self, sys_info, temps):
        """绘制单个系统的详细对比图（带Y轴对齐）"""
        is_ensemble = sys_info.get('is_ensemble', False)
        
        fig, axes = plt.subplots(4, 2, figsize=(14, 14))
        
        # 标题区分 ensemble 和单系统
        if is_ensemble:
            title = f"{sys_info['display_name']}\n温度对比: {temps[0]}K vs {temps[1]}K"
        else:
            title = f"{sys_info['display_name']} 温度对比: {temps[0]}K vs {temps[1]}K"
        
        fig.suptitle(title, fontsize=16, fontweight='bold', y=0.995)
        
        fields = [
            ('q6', 'cluster_metal_q6_global', 'Q6 整体金属'),
            ('q6', 'pt_q6_global', 'Q6 Pt原子'),
            ('q6', 'sn_q6_global', 'Q6 Sn原子'),
            ('op2', 'op2_all_metal', 'OP2 取向参数'),
        ]
        
        # 为 ensemble 的不同来源系统分配颜色
        if is_ensemble:
            source_systems = sys_info.get('source_systems', [])
            source_colors = plt.cm.tab10(np.linspace(0, 1, len(source_systems)))
            source_color_map = {sys: source_colors[i] for i, sys in enumerate(source_systems)}
        
        # ========== 第一遍: 收集所有数据以计算统一Y轴范围 ==========
        y_limits = {}
        all_data = {}  # 存储加载的数据，避免重复读取
        
        for row_idx, (data_type, field, title_label) in enumerate(fields):
            all_values = []
            all_data[row_idx] = {}
            
            for col_idx, temp in enumerate(temps):
                all_data[row_idx][col_idx] = {'dfs': [], 'avg_df': None, 'avg_values': None, 'source_info': []}
                
                if temp not in sys_info['temps']:
                    continue
                
                # 加载所有运行的数据
                all_dfs = []
                source_info = []  # 记录每个df的来源
                for run_info in sys_info['temps'][temp]:
                    df = self.load_time_series(run_info['run_dir'], data_type)
                    if df is not None and field in df.columns:
                        all_dfs.append(df)
                        all_values.extend(df[field].values)
                        # 记录来源系统（用于ensemble着色）
                        source_info.append(run_info.get('source_system', 'unknown'))
                
                if all_dfs:
                    all_data[row_idx][col_idx]['dfs'] = all_dfs
                    all_data[row_idx][col_idx]['source_info'] = source_info
                    
                    # 计算平均
                    if len(all_dfs) > 1:
                        avg_df = pd.concat(all_dfs).groupby('frame').mean().reset_index()
                    else:
                        avg_df = all_dfs[0]
                    
                    all_data[row_idx][col_idx]['avg_df'] = avg_df
                    all_data[row_idx][col_idx]['avg_values'] = avg_df[field].values
            
            # 计算该行的Y轴范围（两个温度统一）
            if all_values and self.unified_ylim:
                y_min = np.min(all_values)
                y_max = np.max(all_values)
                y_range = y_max - y_min
                # 添加5%的边距
                y_limits[row_idx] = (y_min - 0.05*y_range, y_max + 0.05*y_range)
        
        # ========== 第二遍: 绘图 ==========
        for row_idx, (data_type, field, title_label) in enumerate(fields):
            for col_idx, temp in enumerate(temps):
                ax = axes[row_idx, col_idx]
                
                if temp not in sys_info['temps']:
                    ax.text(0.5, 0.5, 'No Data', ha='center', va='center', 
                           transform=ax.transAxes, fontsize=12)
                    continue
                
                data = all_data[row_idx][col_idx]
                all_dfs = data['dfs']
                avg_df = data['avg_df']
                avg_values = data['avg_values']
                source_info = data['source_info']
                
                if not all_dfs or avg_df is None:
                    ax.text(0.5, 0.5, f'No {field} data', ha='center', va='center', 
                           transform=ax.transAxes, fontsize=12)
                    continue
                
                # 绘制原始数据（半透明）
                if is_ensemble:
                    # ensemble: 不同来源用不同颜色
                    for df, src in zip(all_dfs, source_info):
                        color = source_color_map.get(src, 'gray')
                        ax.plot(df['time_ps'], df[field], color=color, alpha=0.3, linewidth=0.8)
                else:
                    # 单系统: 用灰色
                    for df in all_dfs:
                        ax.plot(df['time_ps'], df[field], color='gray', alpha=0.2, linewidth=0.5)
                
                time_ps = avg_df['time_ps'].values
                
                # 绘制平均曲线
                main_color = 'blue' if col_idx == 0 else 'red'
                ax.plot(time_ps, avg_values, color=main_color, linewidth=2.5, 
                       label=f'平均 (n={len(all_dfs)})', alpha=0.8)
                
                # 平滑曲线
                if len(avg_values) > 51:
                    smoothed = savgol_filter(avg_values, window_length=51, polyorder=3)
                    ax.plot(time_ps, smoothed, color=main_color, linewidth=1.5, 
                           linestyle='--', label='平滑', alpha=0.6)
                
                # 统计信息
                stats = self.calculate_statistics(avg_values)
                stats_text = (
                    f"均值: {stats['mean']:.4f}\n"
                    f"标准差: {stats['std']:.4f}\n"
                    f"CV: {stats['cv']:.3f}\n"
                    f"范围: [{stats['min']:.3f}, {stats['max']:.3f}]"
                )
                
                # 根据CV选择背景色
                if stats['cv'] < 0.1:
                    box_color = '#e8f5e9'  # 浅绿 - 很稳定
                elif stats['cv'] < 0.2:
                    box_color = '#fff9c4'  # 浅黄 - 较稳定
                else:
                    box_color = '#ffebee'  # 浅红 - 不稳定
                
                ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
                       verticalalignment='top', fontsize=9,
                       bbox=dict(boxstyle='round', facecolor=box_color, alpha=0.8))
                
                # 设置标签
                if row_idx == 0:
                    ax.set_title(f'{temp}K', fontsize=12, fontweight='bold')
                if col_idx == 0:
                    ax.set_ylabel(title_label, fontsize=10, fontweight='bold')
                if row_idx == 3:
                    ax.set_xlabel('时间 (ps)', fontsize=10)
                
                ax.grid(True, alpha=0.3, linestyle='--')
                ax.legend(loc='upper right', fontsize=8)
                
                # 应用统一Y轴范围
                if self.unified_ylim and row_idx in y_limits:
                    ax.set_ylim(y_limits[row_idx])
        
        # 为 ensemble 图添加颜色图例
        if is_ensemble and len(source_systems) > 0:
            # 在图的底部添加来源系统图例
            legend_elements = []
            for sys_key in source_systems:
                color = source_color_map[sys_key]
                display = self.systems[sys_key]['display_name'] if sys_key in self.systems else sys_key
                legend_elements.append(plt.Line2D([0], [0], color=color, linewidth=2, label=display))
            
            fig.legend(handles=legend_elements, loc='lower center', ncol=len(source_systems),
                      fontsize=9, framealpha=0.9, bbox_to_anchor=(0.5, 0.01))
        
        plt.tight_layout()
        if is_ensemble:
            plt.subplots_adjust(bottom=0.06)  # 为底部图例留空间
        
        # 保存
        safe_name = sys_info['display_name'].replace('/', '_').replace(' ', '_').replace('(', '_').replace(')', '_')
        output_file = self.output_dir / f'individual_{safe_name}_{temps[0]}K_vs_{temps[1]}K.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"📊 {sys_info['display_name']} 对比图已保存: {output_file.name}")
    
    def plot_statistics_heatmap(self):
        """
        模式3: Q6和OP2统计热图（所有系统 × 所有温度）
        """
        print("\n>>> 绘制统计热图...")
        
        # 收集所有温度
        all_temps = set()
        for sys_info in self.systems.values():
            all_temps.update(sys_info['temps'].keys())
        temps = sorted(all_temps)
        
        systems = list(self.systems.keys())
        
        # 创建图表
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        
        metrics = [
            ('q6', 'cluster_metal_q6_global', 'Q6 整体金属'),
            ('op2', 'op2_all_metal', 'OP2 取向参数'),
            ('geometry', 'gyration_radius', '回转半径 Rg'),
            ('q6', 'pt_q6_global', 'Q6 Pt原子'),
            ('q6', 'sn_q6_global', 'Q6 Sn原子'),
            ('op2', 'op2_pt', 'OP2 Pt原子'),
        ]
        
        for idx, (data_type, field, title) in enumerate(metrics):
            ax = axes[idx // 3, idx % 3]
            
            # 创建矩阵
            matrix = np.full((len(temps), len(systems)), np.nan)
            
            for j, sys_key in enumerate(systems):
                sys_info = self.systems[sys_key]
                for i, temp in enumerate(temps):
                    if temp not in sys_info['temps']:
                        continue
                    
                    # 加载并计算平均值
                    all_values = []
                    for run_info in sys_info['temps'][temp]:
                        df = self.load_time_series(run_info['run_dir'], data_type)
                        if df is not None and field in df.columns:
                            all_values.extend(df[field].values)
                    
                    if all_values:
                        matrix[i, j] = np.mean(all_values)
            
            # 检查是否有数据
            if np.all(np.isnan(matrix)):
                ax.text(0.5, 0.5, f'No {field} data', ha='center', va='center',
                       transform=ax.transAxes, fontsize=12)
                ax.set_title(title, fontsize=11, fontweight='bold')
                ax.axis('off')
                continue
            
            # 绘制热图
            im = ax.imshow(matrix, aspect='auto', cmap='RdYlBu_r', interpolation='nearest')
            
            # 设置标签
            ax.set_yticks(range(len(temps)))
            ax.set_yticklabels([f'{t}K' for t in temps], fontsize=8)
            ax.set_xticks(range(len(systems)))
            x_labels = [self.systems[s]['display_name'][:12] for s in systems]
            ax.set_xticklabels(x_labels, fontsize=7, rotation=45, ha='right')
            ax.set_ylabel('温度')
            ax.set_title(title, fontsize=10, fontweight='bold')
            
            # 添加数值
            for i in range(len(temps)):
                for j in range(len(systems)):
                    if not np.isnan(matrix[i, j]):
                        ax.text(j, i, f'{matrix[i, j]:.2f}',
                               ha='center', va='center', color='black', fontsize=6)
            
            plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        
        plt.suptitle('多系统多温度统计热图', fontsize=14, fontweight='bold')
        plt.tight_layout()
        
        output_file = self.output_dir / 'statistics_heatmap.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"📊 统计热图已保存: {output_file}")
    
    def plot_q6_vs_op2_scatter(self, temps=None):
        """
        额外分析: Q6 vs OP2 散点图（探索相关性）
        """
        print("\n>>> 绘制Q6 vs OP2相关性分析...")
        
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        
        # 收集所有数据点
        all_q6 = []
        all_op2 = []
        all_temps_list = []
        all_systems_list = []
        
        for sys_key, sys_info in self.systems.items():
            for temp, runs in sys_info['temps'].items():
                for run_info in runs:
                    df_q6 = self.load_time_series(run_info['run_dir'], 'q6')
                    df_op2 = self.load_time_series(run_info['run_dir'], 'op2')
                    
                    if df_q6 is not None and df_op2 is not None:
                        if 'cluster_metal_q6_global' in df_q6.columns and 'op2_all_metal' in df_op2.columns:
                            q6_mean = df_q6['cluster_metal_q6_global'].mean()
                            op2_mean = df_op2['op2_all_metal'].mean()
                            all_q6.append(q6_mean)
                            all_op2.append(op2_mean)
                            all_temps_list.append(temp)
                            all_systems_list.append(sys_info['display_name'])
        
        if not all_q6:
            print("  ⚠️ 无足够数据绘制Q6 vs OP2图")
            plt.close()
            return
        
        # 左图: 按温度着色
        ax1 = axes[0]
        scatter1 = ax1.scatter(all_q6, all_op2, c=all_temps_list, cmap='coolwarm', 
                               s=50, alpha=0.7, edgecolors='black', linewidths=0.5)
        plt.colorbar(scatter1, ax=ax1, label='温度 (K)')
        
        # 计算相关系数
        from scipy.stats import pearsonr
        r, p = pearsonr(all_q6, all_op2)
        
        ax1.set_xlabel('Q6 (六次对称性)', fontsize=11)
        ax1.set_ylabel('OP2 (取向参数)', fontsize=11)
        ax1.set_title(f'Q6 vs OP2 相关性 (r={r:.3f}, p={p:.2e})', fontsize=12, fontweight='bold')
        ax1.grid(True, alpha=0.3)
        
        # 右图: 按系统类型分组
        ax2 = axes[1]
        unique_systems = list(set(all_systems_list))
        colors = plt.cm.tab10(np.linspace(0, 1, len(unique_systems)))
        
        for i, sys_name in enumerate(unique_systems):
            mask = [s == sys_name for s in all_systems_list]
            q6_subset = [all_q6[j] for j in range(len(mask)) if mask[j]]
            op2_subset = [all_op2[j] for j in range(len(mask)) if mask[j]]
            ax2.scatter(q6_subset, op2_subset, c=[colors[i]], label=sys_name, 
                       s=50, alpha=0.7, edgecolors='black', linewidths=0.5)
        
        ax2.set_xlabel('Q6 (六次对称性)', fontsize=11)
        ax2.set_ylabel('OP2 (取向参数)', fontsize=11)
        ax2.set_title('按系统分类', fontsize=12, fontweight='bold')
        ax2.legend(fontsize=8, loc='best')
        ax2.grid(True, alpha=0.3)
        
        plt.suptitle('Q6 与 OP2 相关性分析', fontsize=14, fontweight='bold')
        plt.tight_layout()
        
        output_file = self.output_dir / 'q6_vs_op2_correlation.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"📊 Q6 vs OP2相关性图已保存: {output_file}")
    
    def plot_stats_comparison(self, temps):
        """
        模式3 (7.6.3): Q6和OP2统计对比
        
        包含6个子图:
        1. Q6均值对比 (柱状图)
        2. Q6 CV对比 (柱状图)
        3. OP2均值对比 (柱状图)
        4. OP2 CV对比 (柱状图)
        5. 温度敏感性散点图 (ΔQ6 vs ΔCV)
        6. 稳定性评分散点图 (Q6 vs CV)
        """
        print(f"\n>>> 绘制统计对比图: {temps[0]}K vs {temps[1]}K")
        
        # 收集所有系统的统计数据
        stats_data = []
        
        for sys_key, sys_info in self.systems.items():
            if temps[0] not in sys_info['temps'] or temps[1] not in sys_info['temps']:
                continue
            
            row = {
                'system': sys_info['display_name'],
                'sys_key': sys_key,
            }
            
            # 计算两个温度下的Q6和OP2统计
            for temp in temps:
                # 加载Q6数据
                q6_values = []
                op2_values = []
                
                for run_info in sys_info['temps'][temp]:
                    df_q6 = self.load_time_series(run_info['run_dir'], 'q6')
                    df_op2 = self.load_time_series(run_info['run_dir'], 'op2')
                    
                    if df_q6 is not None and 'cluster_metal_q6_global' in df_q6.columns:
                        q6_values.extend(df_q6['cluster_metal_q6_global'].values)
                    
                    if df_op2 is not None and 'op2_all_metal' in df_op2.columns:
                        op2_values.extend(df_op2['op2_all_metal'].values)
                
                if q6_values:
                    row[f'q6_mean_{temp}'] = np.mean(q6_values)
                    row[f'q6_std_{temp}'] = np.std(q6_values)
                    row[f'q6_cv_{temp}'] = np.std(q6_values) / np.mean(q6_values) if np.mean(q6_values) != 0 else 0
                else:
                    row[f'q6_mean_{temp}'] = np.nan
                    row[f'q6_std_{temp}'] = np.nan
                    row[f'q6_cv_{temp}'] = np.nan
                
                if op2_values:
                    row[f'op2_mean_{temp}'] = np.mean(op2_values)
                    row[f'op2_std_{temp}'] = np.std(op2_values)
                    row[f'op2_cv_{temp}'] = np.std(op2_values) / np.mean(op2_values) if np.mean(op2_values) != 0 else 0
                else:
                    row[f'op2_mean_{temp}'] = np.nan
                    row[f'op2_std_{temp}'] = np.nan
                    row[f'op2_cv_{temp}'] = np.nan
            
            stats_data.append(row)
        
        if not stats_data:
            print("  ⚠️ 无足够数据进行统计对比")
            return
        
        df = pd.DataFrame(stats_data)
        
        # 计算变化量
        df['delta_q6_mean'] = df[f'q6_mean_{temps[0]}'] - df[f'q6_mean_{temps[1]}']
        df['delta_q6_cv'] = df[f'q6_cv_{temps[1]}'] - df[f'q6_cv_{temps[0]}']
        df['delta_op2_mean'] = df[f'op2_mean_{temps[0]}'] - df[f'op2_mean_{temps[1]}']
        df['delta_op2_cv'] = df[f'op2_cv_{temps[1]}'] - df[f'op2_cv_{temps[0]}']
        
        # 创建图表
        fig = plt.figure(figsize=(18, 16))
        gs = fig.add_gridspec(4, 2, hspace=0.35, wspace=0.3)
        
        colors_t1 = 'steelblue'
        colors_t2 = 'orangered'
        
        x = np.arange(len(df))
        width = 0.35
        
        # === 子图1: Q6均值对比 ===
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.bar(x - width/2, df[f'q6_mean_{temps[0]}'], width, 
               label=f'{temps[0]}K', color=colors_t1, alpha=0.8)
        ax1.bar(x + width/2, df[f'q6_mean_{temps[1]}'], width, 
               label=f'{temps[1]}K', color=colors_t2, alpha=0.8)
        ax1.set_xlabel('系统', fontsize=10)
        ax1.set_ylabel('Q6均值', fontsize=10, fontweight='bold')
        ax1.set_title('Q6均值对比', fontsize=12, fontweight='bold')
        ax1.set_xticks(x)
        ax1.set_xticklabels(df['system'], fontsize=8, rotation=45, ha='right')
        ax1.legend()
        ax1.grid(True, alpha=0.3, axis='y')
        
        # === 子图2: Q6 CV对比 ===
        ax2 = fig.add_subplot(gs[0, 1])
        ax2.bar(x - width/2, df[f'q6_cv_{temps[0]}'], width, 
               label=f'{temps[0]}K', color=colors_t1, alpha=0.8)
        ax2.bar(x + width/2, df[f'q6_cv_{temps[1]}'], width, 
               label=f'{temps[1]}K', color=colors_t2, alpha=0.8)
        ax2.axhline(y=0.1, color='green', linestyle='--', linewidth=1, alpha=0.5, label='CV=0.1 (稳定)')
        ax2.axhline(y=0.2, color='orange', linestyle='--', linewidth=1, alpha=0.5, label='CV=0.2 (临界)')
        ax2.set_xlabel('系统', fontsize=10)
        ax2.set_ylabel('Q6 CV', fontsize=10, fontweight='bold')
        ax2.set_title('Q6变异系数(CV)对比', fontsize=12, fontweight='bold')
        ax2.set_xticks(x)
        ax2.set_xticklabels(df['system'], fontsize=8, rotation=45, ha='right')
        ax2.legend(fontsize=8)
        ax2.grid(True, alpha=0.3, axis='y')
        
        # === 子图3: OP2均值对比 ===
        ax3 = fig.add_subplot(gs[1, 0])
        ax3.bar(x - width/2, df[f'op2_mean_{temps[0]}'], width, 
               label=f'{temps[0]}K', color=colors_t1, alpha=0.8)
        ax3.bar(x + width/2, df[f'op2_mean_{temps[1]}'], width, 
               label=f'{temps[1]}K', color=colors_t2, alpha=0.8)
        ax3.set_xlabel('系统', fontsize=10)
        ax3.set_ylabel('OP2均值', fontsize=10, fontweight='bold')
        ax3.set_title('OP2均值对比', fontsize=12, fontweight='bold')
        ax3.set_xticks(x)
        ax3.set_xticklabels(df['system'], fontsize=8, rotation=45, ha='right')
        ax3.legend()
        ax3.grid(True, alpha=0.3, axis='y')
        
        # === 子图4: OP2 CV对比 ===
        ax4 = fig.add_subplot(gs[1, 1])
        ax4.bar(x - width/2, df[f'op2_cv_{temps[0]}'], width, 
               label=f'{temps[0]}K', color=colors_t1, alpha=0.8)
        ax4.bar(x + width/2, df[f'op2_cv_{temps[1]}'], width, 
               label=f'{temps[1]}K', color=colors_t2, alpha=0.8)
        ax4.axhline(y=0.1, color='green', linestyle='--', linewidth=1, alpha=0.5, label='CV=0.1 (稳定)')
        ax4.axhline(y=0.2, color='orange', linestyle='--', linewidth=1, alpha=0.5, label='CV=0.2 (临界)')
        ax4.set_xlabel('系统', fontsize=10)
        ax4.set_ylabel('OP2 CV', fontsize=10, fontweight='bold')
        ax4.set_title('OP2变异系数(CV)对比', fontsize=12, fontweight='bold')
        ax4.set_xticks(x)
        ax4.set_xticklabels(df['system'], fontsize=8, rotation=45, ha='right')
        ax4.legend(fontsize=8)
        ax4.grid(True, alpha=0.3, axis='y')
        
        # === 子图5: Q6温度敏感性 (ΔQ6 vs ΔCV) ===
        ax5 = fig.add_subplot(gs[2, 0])
        valid_mask = ~(df['delta_q6_mean'].isna() | df['delta_q6_cv'].isna())
        if valid_mask.any():
            scatter = ax5.scatter(df.loc[valid_mask, 'delta_q6_mean'], 
                                 df.loc[valid_mask, 'delta_q6_cv'],
                                 c=range(valid_mask.sum()), cmap='viridis',
                                 s=200, alpha=0.7, edgecolors='black', linewidth=1.5)
            for i, (idx, row) in enumerate(df[valid_mask].iterrows()):
                ax5.annotate(row['system'][:8], 
                            (row['delta_q6_mean'], row['delta_q6_cv']),
                            fontsize=8, ha='center', va='center')
        ax5.axhline(y=0, color='gray', linestyle='-', linewidth=1, alpha=0.5)
        ax5.axvline(x=0, color='gray', linestyle='-', linewidth=1, alpha=0.5)
        ax5.set_xlabel(f'ΔQ6 (Q6@{temps[0]}K - Q6@{temps[1]}K)', fontsize=10)
        ax5.set_ylabel(f'ΔCV (CV@{temps[1]}K - CV@{temps[0]}K)', fontsize=10)
        ax5.set_title('Q6温度敏感性分析', fontsize=12, fontweight='bold')
        ax5.grid(True, alpha=0.3)
        ax5.text(0.02, 0.98, '左下: 温度稳定\n右上: 温度敏感', 
                transform=ax5.transAxes, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5), fontsize=9)
        
        # === 子图6: OP2温度敏感性 (ΔOP2 vs ΔCV) ===
        ax6 = fig.add_subplot(gs[2, 1])
        valid_mask = ~(df['delta_op2_mean'].isna() | df['delta_op2_cv'].isna())
        if valid_mask.any():
            scatter = ax6.scatter(df.loc[valid_mask, 'delta_op2_mean'], 
                                 df.loc[valid_mask, 'delta_op2_cv'],
                                 c=range(valid_mask.sum()), cmap='viridis',
                                 s=200, alpha=0.7, edgecolors='black', linewidth=1.5)
            for i, (idx, row) in enumerate(df[valid_mask].iterrows()):
                ax6.annotate(row['system'][:8], 
                            (row['delta_op2_mean'], row['delta_op2_cv']),
                            fontsize=8, ha='center', va='center')
        ax6.axhline(y=0, color='gray', linestyle='-', linewidth=1, alpha=0.5)
        ax6.axvline(x=0, color='gray', linestyle='-', linewidth=1, alpha=0.5)
        ax6.set_xlabel(f'ΔOP2 (OP2@{temps[0]}K - OP2@{temps[1]}K)', fontsize=10)
        ax6.set_ylabel(f'ΔCV (CV@{temps[1]}K - CV@{temps[0]}K)', fontsize=10)
        ax6.set_title('OP2温度敏感性分析', fontsize=12, fontweight='bold')
        ax6.grid(True, alpha=0.3)
        ax6.text(0.02, 0.98, '左下: 温度稳定\n右上: 温度敏感', 
                transform=ax6.transAxes, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5), fontsize=9)
        
        # === 子图7: Q6稳定性评分 ===
        ax7 = fig.add_subplot(gs[3, 0])
        valid_mask = ~(df[f'q6_mean_{temps[0]}'].isna() | df[f'q6_cv_{temps[0]}'].isna())
        if valid_mask.any():
            scatter1 = ax7.scatter(df.loc[valid_mask, f'q6_mean_{temps[0]}'], 
                                  df.loc[valid_mask, f'q6_cv_{temps[0]}'],
                                  c='steelblue', s=150, alpha=0.7, marker='o',
                                  edgecolors='black', linewidth=1, label=f'{temps[0]}K')
        valid_mask = ~(df[f'q6_mean_{temps[1]}'].isna() | df[f'q6_cv_{temps[1]}'].isna())
        if valid_mask.any():
            scatter2 = ax7.scatter(df.loc[valid_mask, f'q6_mean_{temps[1]}'], 
                                  df.loc[valid_mask, f'q6_cv_{temps[1]}'],
                                  c='orangered', s=150, alpha=0.7, marker='s',
                                  edgecolors='black', linewidth=1, label=f'{temps[1]}K')
        ax7.axhline(y=0.1, color='green', linestyle='--', linewidth=1, alpha=0.5)
        ax7.axhline(y=0.2, color='orange', linestyle='--', linewidth=1, alpha=0.5)
        ax7.set_xlabel('Q6均值', fontsize=10, fontweight='bold')
        ax7.set_ylabel('CV', fontsize=10, fontweight='bold')
        ax7.set_title('Q6稳定性评分 (理想: 高Q6 + 低CV)', fontsize=12, fontweight='bold')
        ax7.legend()
        ax7.grid(True, alpha=0.3)
        
        # === 子图8: OP2稳定性评分 ===
        ax8 = fig.add_subplot(gs[3, 1])
        valid_mask = ~(df[f'op2_mean_{temps[0]}'].isna() | df[f'op2_cv_{temps[0]}'].isna())
        if valid_mask.any():
            scatter1 = ax8.scatter(df.loc[valid_mask, f'op2_mean_{temps[0]}'], 
                                  df.loc[valid_mask, f'op2_cv_{temps[0]}'],
                                  c='steelblue', s=150, alpha=0.7, marker='o',
                                  edgecolors='black', linewidth=1, label=f'{temps[0]}K')
        valid_mask = ~(df[f'op2_mean_{temps[1]}'].isna() | df[f'op2_cv_{temps[1]}'].isna())
        if valid_mask.any():
            scatter2 = ax8.scatter(df.loc[valid_mask, f'op2_mean_{temps[1]}'], 
                                  df.loc[valid_mask, f'op2_cv_{temps[1]}'],
                                  c='orangered', s=150, alpha=0.7, marker='s',
                                  edgecolors='black', linewidth=1, label=f'{temps[1]}K')
        ax8.axhline(y=0.1, color='green', linestyle='--', linewidth=1, alpha=0.5)
        ax8.axhline(y=0.2, color='orange', linestyle='--', linewidth=1, alpha=0.5)
        ax8.set_xlabel('OP2均值', fontsize=10, fontweight='bold')
        ax8.set_ylabel('CV', fontsize=10, fontweight='bold')
        ax8.set_title('OP2稳定性评分', fontsize=12, fontweight='bold')
        ax8.legend()
        ax8.grid(True, alpha=0.3)
        
        plt.suptitle(f'Q6 & OP2 统计对比: {temps[0]}K vs {temps[1]}K', 
                    fontsize=16, fontweight='bold')
        
        # 保存图表
        output_file = self.output_dir / f'stats_comparison_{temps[0]}K_vs_{temps[1]}K.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"📊 统计对比图已保存: {output_file}")
        
        # 保存CSV
        csv_file = self.output_dir / f'stats_comparison_{temps[0]}K_vs_{temps[1]}K.csv'
        df.to_csv(csv_file, index=False, float_format='%.6f', encoding='utf-8-sig')
        print(f"📄 统计数据已保存: {csv_file}")
        
        # 打印统计摘要
        print(f"\n>>> 统计摘要:")
        for _, row in df.iterrows():
            print(f"  {row['system']}:")
            print(f"    Q6: {temps[0]}K={row[f'q6_mean_{temps[0]}']:.4f}(CV={row[f'q6_cv_{temps[0]}']:.3f}), "
                  f"{temps[1]}K={row[f'q6_mean_{temps[1]}']:.4f}(CV={row[f'q6_cv_{temps[1]}']:.3f})")
            if not np.isnan(row[f'op2_mean_{temps[0]}']):
                print(f"    OP2: {temps[0]}K={row[f'op2_mean_{temps[0]}']:.4f}(CV={row[f'op2_cv_{temps[0]}']:.3f}), "
                      f"{temps[1]}K={row[f'op2_mean_{temps[1]}']:.4f}(CV={row[f'op2_cv_{temps[1]}']:.3f})")
    
    def run_analysis(self, mode='all', temps=None):
        """运行分析"""
        # 扫描数据
        self.scan_data_structure()
        
        if not self.systems:
            print("❌ 未找到任何数据!")
            return
        
        # 默认温度对比
        if temps is None:
            all_temps = set()
            for sys_info in self.systems.values():
                all_temps.update(sys_info['temps'].keys())
            temps_sorted = sorted(all_temps)
            if len(temps_sorted) >= 2:
                temps = [temps_sorted[0], temps_sorted[-1]]  # 最低温和最高温
            else:
                temps = [temps_sorted[0], temps_sorted[0]]
        
        print(f"\n>>> 温度对比: {temps[0]}K vs {temps[1]}K")
        
        # 根据模式运行
        if mode == 'sidebyside' or mode == 'all':
            self.plot_side_by_side_comparison(temps)
        
        if mode == 'individual' or mode == 'all':
            self.plot_individual_system_comparison(temps)
        
        if mode == 'stats' or mode == 'all':
            self.plot_stats_comparison(temps)
        
        if mode == 'heatmap' or mode == 'all':
            self.plot_statistics_heatmap()
        
        if mode == 'all':
            self.plot_q6_vs_op2_scatter()
        
        print("\n" + "="*80)
        print("✅ Step 7.6 分析完成!")
        print(f"结果保存在: {self.output_dir}")
        print("="*80)


def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description='Step 7.6 新服务器数据温度对比分析',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例用法:

1. 运行所有分析模式:
   python step7-6-new_server_temp_comparison.py

2. 并排对比两个温度 (7.6.1):
   python step7-6-new_server_temp_comparison.py --mode sidebyside --temps 300,900

3. 单个系统详细分析 (7.6.2):
   python step7-6-new_server_temp_comparison.py --mode individual --temps 200,1100

4. 统计对比柱状图和散点图 (7.6.3):
   python step7-6-new_server_temp_comparison.py --mode stats --temps 300,900

5. 统计热图:
   python step7-6-new_server_temp_comparison.py --mode heatmap

6. 禁用Y轴对齐（每列独立调整Y轴范围）:
   python step7-6-new_server_temp_comparison.py --mode individual --temps 300,900 --no-unified-ylim

Y轴对齐说明:
- 默认启用: 同一行的左右两列使用相同Y轴范围，便于直观对比温度差异
- 使用 --no-unified-ylim 禁用: 每列独立调整Y轴，可显示更多细节

变异系数(CV)说明:
- CV < 0.1  → 非常稳定 (绿色背景)
- 0.1 < CV < 0.2 → 较稳定 (黄色背景)
- CV > 0.2  → 不稳定 (红色背景)
        """
    )
    
    parser.add_argument('--data', type=str, 
                       default=r'C:\Users\11207\OneDrive\02_Code\work1-PtSnO\workflow\data\coordination\seletion\coordination_time_series_results_sample_20251130_193923',
                       help='数据根目录')
    parser.add_argument('--output', type=str, default=None,
                       help='输出目录')
    parser.add_argument('--mode', type=str, default='all',
                       choices=['all', 'sidebyside', 'individual', 'stats', 'heatmap'],
                       help='分析模式: all=全部, sidebyside=并排对比, individual=单系统, stats=统计对比, heatmap=热图')
    parser.add_argument('--temps', type=str, default=None,
                       help='温度对比(逗号分隔), 如: 300,900')
    parser.add_argument('--unified-ylim', action='store_true', default=True,
                       help='统一Y轴范围(默认启用，同一行左右两列使用相同Y轴)')
    parser.add_argument('--no-unified-ylim', action='store_true',
                       help='禁用统一Y轴范围(每列独立调整Y轴)')
    
    args = parser.parse_args()
    
    # 解析温度
    temps = None
    if args.temps:
        temps = [int(t.strip()) for t in args.temps.split(',')]
        if len(temps) != 2:
            print("❌ 错误: 必须指定恰好2个温度")
            return
    
    # 确定Y轴对齐选项
    unified_ylim = True  # 默认启用
    if args.no_unified_ylim:
        unified_ylim = False
    
    # 创建分析器并运行
    analyzer = NewServerTempComparisonAnalyzer(args.data, args.output, unified_ylim=unified_ylim)
    analyzer.run_analysis(mode=args.mode, temps=temps)


if __name__ == '__main__':
    main()
