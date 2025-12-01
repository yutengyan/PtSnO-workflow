#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 7.1.2: Air-68 vs Air-86 系综平均对比分析

对比两种气相团簇在相同温度下的Q6/OP2时间演化：
- Air-68 (Pt6Sn8): 6个Pt + 8个Sn = 14原子
- Air-86 (Pt8Sn6): 8个Pt + 6个Sn = 14原子

功能：
1. 统一Y轴范围便于对比
2. 系综平均（多个初始结构平均）
3. 计算并显示统计量

Author: AI Assistant
Date: 2025-12-01
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
from pathlib import Path

# 设置控制台输出编码
if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8')

warnings.filterwarnings('ignore')

# 配置中文字体
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'Arial']
plt.rcParams['axes.unicode_minus'] = False

try:
    chinese_font = FontProperties(family='Microsoft YaHei', size=9)
except:
    chinese_font = FontProperties(size=9)


class AirEnsembleComparison:
    """Air系列系综平均对比分析器"""
    
    def __init__(self, data_root, output_dir=None):
        self.data_root = Path(data_root)
        if output_dir:
            self.output_dir = Path(output_dir)
        else:
            self.output_dir = Path(__file__).parent / 'results' / 'step7.1.2_air_comparison'
        
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # 存储数据
        self.air68_data = {}  # {temp: [df1, df2, ...]}
        self.air86_data = {}
        
    def scan_data(self):
        """扫描并分类数据"""
        print(">>> 扫描数据结构...")
        
        for root, dirs, files in os.walk(self.data_root):
            root_path = Path(root)
            
            # 查找Q6时间序列文件
            q6_file = root_path / 'cluster_global_q6_time_series.csv'
            op2_file = root_path / 'op2_time_series.csv'
            
            if q6_file.exists():
                # 解析路径确定系统类型和温度
                path_str = str(root_path)
                
                # 判断是68还是86
                is_68 = '/68' in path_str or '\\68' in path_str
                is_86 = '/86' in path_str or '\\86' in path_str
                
                if not (is_68 or is_86):
                    continue
                
                # 提取温度
                temp = None
                for part in root_path.parts:
                    if part.startswith('T') and part[1:].split('.')[0].isdigit():
                        temp = int(part[1:].split('.')[0])
                        break
                
                if temp is None:
                    continue
                
                # 加载数据
                try:
                    df_q6 = pd.read_csv(q6_file)
                    df_op2 = pd.read_csv(op2_file) if op2_file.exists() else None
                    
                    data_entry = {
                        'q6': df_q6,
                        'op2': df_op2,
                        'path': str(root_path)
                    }
                    
                    if is_68:
                        if temp not in self.air68_data:
                            self.air68_data[temp] = []
                        self.air68_data[temp].append(data_entry)
                    else:
                        if temp not in self.air86_data:
                            self.air86_data[temp] = []
                        self.air86_data[temp].append(data_entry)
                        
                except Exception as e:
                    print(f"  ⚠️ 加载失败: {root_path}: {e}")
        
        print(f"    Air-68 数据: {len(self.air68_data)} 个温度")
        for temp, data_list in sorted(self.air68_data.items()):
            print(f"      {temp}K: {len(data_list)} 个结构")
        
        print(f"    Air-86 数据: {len(self.air86_data)} 个温度")
        for temp, data_list in sorted(self.air86_data.items()):
            print(f"      {temp}K: {len(data_list)} 个结构")
    
    def compute_ensemble_average(self, data_list, field, source='q6'):
        """计算系综平均"""
        all_values = []
        
        for data_entry in data_list:
            df = data_entry[source]
            if df is not None and field in df.columns:
                all_values.append(df[field].values)
        
        if not all_values:
            return None, None, None
        
        # 找到最短长度
        min_len = min(len(v) for v in all_values)
        all_values = [v[:min_len] for v in all_values]
        
        # 计算平均和标准差
        values_array = np.array(all_values)
        mean = np.mean(values_array, axis=0)
        std = np.std(values_array, axis=0)
        
        # 平滑处理
        if len(mean) > 21:
            mean_smooth = savgol_filter(mean, min(21, len(mean)//2*2+1), 3)
        else:
            mean_smooth = mean
        
        return mean_smooth, std, len(all_values)
    
    def plot_comparison(self, temp, fields_config):
        """
        绘制Air-68 vs Air-86对比图
        
        Args:
            temp: 温度 (K)
            fields_config: 字段配置列表 [(field_name, source, title, ylabel), ...]
        """
        if temp not in self.air68_data or temp not in self.air86_data:
            print(f"  ⚠️ 温度 {temp}K 数据不完整")
            return
        
        n_fields = len(fields_config)
        fig, axes = plt.subplots(n_fields, 2, figsize=(14, 4*n_fields))
        
        if n_fields == 1:
            axes = axes.reshape(1, -1)
        
        colors = {'air68': '#E74C3C', 'air86': '#3498DB'}  # 红色和蓝色
        
        for row, (field, source, title, ylabel) in enumerate(fields_config):
            ax_68 = axes[row, 0]
            ax_86 = axes[row, 1]
            
            # 计算系综平均
            mean_68, std_68, n_68 = self.compute_ensemble_average(
                self.air68_data[temp], field, source)
            mean_86, std_86, n_86 = self.compute_ensemble_average(
                self.air86_data[temp], field, source)
            
            # 确定统一Y轴范围
            y_values = []
            if mean_68 is not None and std_68 is not None:
                y_values.extend(mean_68 - std_68)
                y_values.extend(mean_68 + std_68)
            if mean_86 is not None and std_86 is not None:
                y_values.extend(mean_86 - std_86)
                y_values.extend(mean_86 + std_86)
            
            if y_values:
                y_min, y_max = np.min(y_values), np.max(y_values)
                y_margin = (y_max - y_min) * 0.15
                y_min -= y_margin
                y_max += y_margin
            else:
                y_min, y_max = 0, 1
            
            # 绘制Air-68
            if mean_68 is not None:
                frames_68 = np.arange(len(mean_68))
                ax_68.fill_between(frames_68, mean_68 - std_68, mean_68 + std_68, 
                                  alpha=0.3, color=colors['air68'])
                ax_68.plot(frames_68, mean_68, color=colors['air68'], linewidth=1.5,
                          label=f'Air-68 系综平均 (n={n_68})')
                
                # 统计信息
                stats_68 = {
                    'mean': np.mean(mean_68),
                    'std': np.std(mean_68),
                    'cv': np.std(mean_68) / np.mean(mean_68) * 100 if np.mean(mean_68) != 0 else 0
                }
                stats_text_68 = f"均值: {stats_68['mean']:.4f}\n标准差: {stats_68['std']:.4f}\nCV: {stats_68['cv']:.1f}%"
                ax_68.text(0.98, 0.98, stats_text_68, transform=ax_68.transAxes,
                          verticalalignment='top', horizontalalignment='right',
                          bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
                          fontsize=9, fontproperties=chinese_font)
            
            ax_68.set_xlim(0, len(mean_68) if mean_68 is not None else 100)
            ax_68.set_ylim(y_min, y_max)
            ax_68.set_xlabel('Frame', fontsize=10)
            ax_68.set_ylabel(ylabel, fontsize=10, fontweight='bold')
            ax_68.set_title(f'Air-68 (Pt₆Sn₈) - {title}', fontsize=11, fontweight='bold')
            ax_68.legend(loc='upper left', fontsize=9)
            ax_68.grid(True, alpha=0.3)
            
            # 绘制Air-86
            if mean_86 is not None:
                frames_86 = np.arange(len(mean_86))
                ax_86.fill_between(frames_86, mean_86 - std_86, mean_86 + std_86, 
                                  alpha=0.3, color=colors['air86'])
                ax_86.plot(frames_86, mean_86, color=colors['air86'], linewidth=1.5,
                          label=f'Air-86 系综平均 (n={n_86})')
                
                # 统计信息
                stats_86 = {
                    'mean': np.mean(mean_86),
                    'std': np.std(mean_86),
                    'cv': np.std(mean_86) / np.mean(mean_86) * 100 if np.mean(mean_86) != 0 else 0
                }
                stats_text_86 = f"均值: {stats_86['mean']:.4f}\n标准差: {stats_86['std']:.4f}\nCV: {stats_86['cv']:.1f}%"
                ax_86.text(0.98, 0.98, stats_text_86, transform=ax_86.transAxes,
                          verticalalignment='top', horizontalalignment='right',
                          bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8),
                          fontsize=9, fontproperties=chinese_font)
            
            ax_86.set_xlim(0, len(mean_86) if mean_86 is not None else 100)
            ax_86.set_ylim(y_min, y_max)
            ax_86.set_xlabel('Frame', fontsize=10)
            ax_86.set_ylabel(ylabel, fontsize=10, fontweight='bold')
            ax_86.set_title(f'Air-86 (Pt₈Sn₆) - {title}', fontsize=11, fontweight='bold')
            ax_86.legend(loc='upper left', fontsize=9)
            ax_86.grid(True, alpha=0.3)
        
        plt.suptitle(f'Air-68 vs Air-86 系综平均对比 @ {temp}K', 
                    fontsize=14, fontweight='bold')
        plt.tight_layout()
        
        output_file = self.output_dir / f'air68_vs_air86_ensemble_{temp}K.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"📊 对比图已保存: {output_file.name}")
        
        return output_file
    
    def plot_overlay_comparison(self, temp, fields_config):
        """
        绘制Air-68和Air-86叠加对比图（同一坐标系）
        """
        if temp not in self.air68_data or temp not in self.air86_data:
            print(f"  ⚠️ 温度 {temp}K 数据不完整")
            return
        
        n_fields = len(fields_config)
        fig, axes = plt.subplots(n_fields, 1, figsize=(12, 4*n_fields))
        
        if n_fields == 1:
            axes = [axes]
        
        colors = {'air68': '#E74C3C', 'air86': '#3498DB'}
        
        stats_data = []
        
        for row, (field, source, title, ylabel) in enumerate(fields_config):
            ax = axes[row]
            
            # 计算系综平均
            mean_68, std_68, n_68 = self.compute_ensemble_average(
                self.air68_data[temp], field, source)
            mean_86, std_86, n_86 = self.compute_ensemble_average(
                self.air86_data[temp], field, source)
            
            # 绘制Air-68
            if mean_68 is not None:
                frames_68 = np.arange(len(mean_68))
                ax.fill_between(frames_68, mean_68 - std_68, mean_68 + std_68, 
                               alpha=0.2, color=colors['air68'])
                ax.plot(frames_68, mean_68, color=colors['air68'], linewidth=2,
                       label=f'Air-68 (Pt₆Sn₈, n={n_68})')
                
                stats_68 = {
                    'system': 'Air-68',
                    'field': field,
                    'mean': np.mean(mean_68),
                    'std': np.std(mean_68),
                    'cv': np.std(mean_68) / np.mean(mean_68) * 100
                }
                stats_data.append(stats_68)
            
            # 绘制Air-86
            if mean_86 is not None:
                frames_86 = np.arange(len(mean_86))
                ax.fill_between(frames_86, mean_86 - std_86, mean_86 + std_86, 
                               alpha=0.2, color=colors['air86'])
                ax.plot(frames_86, mean_86, color=colors['air86'], linewidth=2,
                       label=f'Air-86 (Pt₈Sn₆, n={n_86})')
                
                stats_86 = {
                    'system': 'Air-86',
                    'field': field,
                    'mean': np.mean(mean_86),
                    'std': np.std(mean_86),
                    'cv': np.std(mean_86) / np.mean(mean_86) * 100
                }
                stats_data.append(stats_86)
            
            # 添加统计信息框
            if mean_68 is not None and mean_86 is not None:
                stats_text = (
                    f"Air-68: μ={stats_68['mean']:.4f}, CV={stats_68['cv']:.1f}%\n"
                    f"Air-86: μ={stats_86['mean']:.4f}, CV={stats_86['cv']:.1f}%\n"
                    f"Δ(68-86): {stats_68['mean'] - stats_86['mean']:+.4f}"
                )
                ax.text(0.98, 0.98, stats_text, transform=ax.transAxes,
                       verticalalignment='top', horizontalalignment='right',
                       bbox=dict(boxstyle='round', facecolor='white', alpha=0.9),
                       fontsize=10, family='monospace')
            
            ax.set_xlabel('Frame', fontsize=11)
            ax.set_ylabel(ylabel, fontsize=11, fontweight='bold')
            ax.set_title(f'{title} @ {temp}K', fontsize=12, fontweight='bold')
            ax.legend(loc='upper left', fontsize=10)
            ax.grid(True, alpha=0.3)
        
        plt.suptitle(f'Air-68 (Pt₆Sn₈) vs Air-86 (Pt₈Sn₆) 系综平均叠加对比 @ {temp}K', 
                    fontsize=14, fontweight='bold')
        plt.tight_layout()
        
        output_file = self.output_dir / f'air68_vs_air86_overlay_{temp}K.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"📊 叠加对比图已保存: {output_file.name}")
        
        # 保存统计数据
        if stats_data:
            stats_df = pd.DataFrame(stats_data)
            csv_file = self.output_dir / f'air68_vs_air86_stats_{temp}K.csv'
            stats_df.to_csv(csv_file, index=False, float_format='%.6f')
            print(f"📄 统计数据已保存: {csv_file.name}")
        
        return output_file
    
    def run_analysis(self, temps=None):
        """运行完整分析"""
        self.scan_data()
        
        if temps is None:
            # 找出共同温度
            common_temps = set(self.air68_data.keys()) & set(self.air86_data.keys())
            temps = sorted(common_temps)
        
        if not temps:
            print("❌ 没有找到共同的温度数据")
            return
        
        print(f"\n>>> 分析温度: {temps}")
        
        # 定义要分析的字段
        fields_config = [
            ('cluster_metal_q6_global', 'q6', 'Q6 (全局)', 'Q6'),
            ('pt_q6_global', 'q6', 'Pt-Q6', 'Pt Q6'),
            ('sn_q6_global', 'q6', 'Sn-Q6', 'Sn Q6'),
            ('op2_all_metal', 'op2', 'OP2 (全部金属)', 'OP2'),
        ]
        
        for temp in temps:
            print(f"\n>>> 处理 {temp}K...")
            # 并排对比
            self.plot_comparison(temp, fields_config)
            # 叠加对比
            self.plot_overlay_comparison(temp, fields_config)
        
        print(f"\n{'='*60}")
        print(f"✅ 分析完成! 结果保存在: {self.output_dir}")
        print(f"{'='*60}")


def main():
    parser = argparse.ArgumentParser(
        description='Air-68 vs Air-86 系综平均对比分析',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('--data', type=str, 
                       default=str(Path(__file__).parent / 'data' / 'coordination' / 'air' / 
                                  'coordination_time_series_results_air-sample_20251130_211818'),
                       help='数据目录')
    parser.add_argument('--output', type=str, default=None,
                       help='输出目录')
    parser.add_argument('--temps', type=str, default='300',
                       help='温度列表，逗号分隔 (默认: 300)')
    
    args = parser.parse_args()
    
    temps = [int(t.strip()) for t in args.temps.split(',')]
    
    print(f"\n{'='*60}")
    print(f"Step 7.1.2: Air-68 vs Air-86 系综平均对比分析")
    print(f"数据目录: {args.data}")
    print(f"温度: {temps}")
    print(f"{'='*60}")
    
    analyzer = AirEnsembleComparison(args.data, args.output)
    analyzer.run_analysis(temps)


if __name__ == '__main__':
    main()
