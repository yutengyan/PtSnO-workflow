#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
新服务器数据分析脚本 - 分析coordination_time_series_results数据

支持新的服务器输出格式:
dp-md/4090-ustc/more/Pt6/pt6sn8/T500.r3.gpu0/
dp-md/4090-ustc/o68/g-1535-Sn8Pt6O4/Cv-1/T500.r6.gpu0/

数据包含:
- coordination_time_series.csv
- cluster_global_q6_time_series.csv
- cluster_op2_time_series.csv
- cluster_geometry_time_series.csv
- element_comparison.csv
- analysis_info.txt

Author: AI Assistant
Date: 2025-11-30
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager
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


class NewServerDataAnalyzer:
    """新服务器数据分析器"""
    
    def __init__(self, data_root, output_dir=None):
        """
        初始化分析器
        
        Args:
            data_root: 解压后的数据根目录
            output_dir: 输出目录
        """
        self.data_root = Path(data_root)
        if output_dir:
            self.output_dir = Path(output_dir)
        else:
            self.output_dir = Path(__file__).parent / 'results' / 'step7.5.new_server'
        
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # 存储发现的系统
        self.systems = {}
        self.results = {}
        
        print(f"\n{'='*80}")
        print(f"新服务器数据分析器")
        print(f"数据目录: {self.data_root}")
        print(f"输出目录: {self.output_dir}")
        print(f"{'='*80}")
    
    def scan_data_structure(self):
        """扫描数据目录结构，发现所有可用的系统和温度"""
        print("\n>>> 扫描数据结构...")
        
        # 递归查找所有包含 coordination_time_series.csv 的目录
        csv_files = list(self.data_root.rglob('coordination_time_series.csv'))
        
        print(f"    找到 {len(csv_files)} 个数据点")
        
        for csv_path in csv_files:
            run_dir = csv_path.parent  # T500.r3.gpu0
            system_dir = run_dir.parent  # pt6sn8 或 Cv-1
            
            # 提取温度和运行信息
            run_name = run_dir.name  # T500.r3.gpu0
            match = re.match(r'T(\d+)\.r(\d+)\.gpu\d+', run_name)
            if not match:
                continue
            
            temp = int(match.group(1))
            run_num = int(match.group(2))
            
            # 提取系统名称
            system_name = system_dir.name.lower()
            
            # 获取完整路径信息用于分组
            path_parts = csv_path.relative_to(self.data_root).parts
            
            # 创建唯一的系统键
            # 例如: "more/Pt6/pt6sn8" 或 "o68/g-1535-Sn8Pt6O4/Cv-1"
            if len(path_parts) >= 4:
                system_key = '/'.join(path_parts[:-2])  # 去掉温度目录和文件名
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
        
        # 打印发现的系统
        print(f"\n    发现 {len(self.systems)} 个系统:")
        for sys_key, sys_info in self.systems.items():
            temps = sorted(sys_info['temps'].keys())
            total_runs = sum(len(runs) for runs in sys_info['temps'].values())
            print(f"      {sys_info['display_name']}: {len(temps)}个温度, {total_runs}次运行")
            print(f"        温度范围: {min(temps)}K - {max(temps)}K")
    
    def _get_display_name(self, system_name, path_parts):
        """生成显示名称"""
        # 检测是否是Cv系列
        if system_name.startswith('cv-'):
            return f"Cv系列-{system_name}"
        
        # 检测是否是pt6sn8系列
        if 'pt6sn8' in system_name:
            # 查看路径确定是哪个运行
            for part in path_parts:
                if part.startswith('Pt6-'):
                    return f"Pt6Sn8 ({part})"
                if part == 'Pt6':
                    return "Pt6Sn8 (Pt6)"
            return f"Pt6Sn8"
        
        return system_name.upper()
    
    def load_run_data(self, run_dir):
        """加载单次运行的所有数据"""
        data = {}
        
        # 加载配位数时间序列
        cn_file = run_dir / 'coordination_time_series.csv'
        if cn_file.exists():
            data['coordination'] = pd.read_csv(cn_file)
        
        # 加载Q6时间序列
        q6_file = run_dir / 'cluster_global_q6_time_series.csv'
        if q6_file.exists():
            data['q6'] = pd.read_csv(q6_file)
        
        # 加载OP2时间序列
        op2_file = run_dir / 'cluster_op2_time_series.csv'
        if op2_file.exists():
            data['op2'] = pd.read_csv(op2_file)
        
        # 加载几何时间序列
        geo_file = run_dir / 'cluster_geometry_time_series.csv'
        if geo_file.exists():
            data['geometry'] = pd.read_csv(geo_file)
        
        # 加载元素对比数据
        elem_file = run_dir / 'element_comparison.csv'
        if elem_file.exists():
            data['element'] = pd.read_csv(elem_file)
        
        return data
    
    def analyze_system(self, sys_key, sys_info):
        """分析单个系统的所有数据"""
        print(f"\n>>> 分析系统: {sys_info['display_name']}")
        
        temp_results = {}
        
        for temp in sorted(sys_info['temps'].keys()):
            runs = sys_info['temps'][temp]
            run_data_list = []
            
            for run_info in runs:
                data = self.load_run_data(run_info['run_dir'])
                if data:
                    run_data_list.append(data)
            
            if not run_data_list:
                continue
            
            # 计算该温度的统计量
            stats = self._calculate_statistics(run_data_list)
            stats['n_runs'] = len(run_data_list)
            temp_results[temp] = stats
            
            # 打印信息
            info = f"  [OK] {temp}K: Q6={stats['q6']:.3f}"
            if stats.get('op2_all_metal', 0) > 0:
                info += f", OP2={stats['op2_all_metal']:.3f}"
            if stats.get('rg', 0) > 0:
                info += f", Rg={stats['rg']:.3f}Å"
            info += f" ({len(run_data_list)}次运行)"
            print(info)
        
        return temp_results
    
    def _calculate_statistics(self, run_data_list):
        """计算多次运行的统计量"""
        stats = {
            'q6': [], 'q6_pt': [], 'q6_sn': [],
            'op2_all_metal': [], 'op2_pt': [], 'op2_sn': [],
            'rg': [], 'pt_dist': [], 'sn_dist': [],
            'pt_cn_total': [], 'pt_pt_bonds': [], 'pt_sn_bonds': [],
            'sn_cn_total': [], 'sn_sn_bonds': [], 'sn_pt_bonds': [],
        }
        
        for data in run_data_list:
            # Q6数据
            if 'q6' in data:
                df = data['q6']
                if 'cluster_metal_q6_global' in df.columns:
                    stats['q6'].append(df['cluster_metal_q6_global'].mean())
                if 'pt_q6_global' in df.columns:
                    stats['q6_pt'].append(df['pt_q6_global'].mean())
                if 'sn_q6_global' in df.columns:
                    stats['q6_sn'].append(df['sn_q6_global'].mean())
            
            # OP2数据
            if 'op2' in data:
                df = data['op2']
                if 'op2_all_metal' in df.columns:
                    stats['op2_all_metal'].append(df['op2_all_metal'].mean())
                if 'op2_pt' in df.columns:
                    stats['op2_pt'].append(df['op2_pt'].mean())
                if 'op2_sn' in df.columns:
                    stats['op2_sn'].append(df['op2_sn'].mean())
            
            # 几何数据
            if 'geometry' in data:
                df = data['geometry']
                if 'gyration_radius' in df.columns:
                    stats['rg'].append(df['gyration_radius'].mean())
                if 'pt_avg_dist_to_center' in df.columns:
                    stats['pt_dist'].append(df['pt_avg_dist_to_center'].mean())
                if 'sn_avg_dist_to_center' in df.columns:
                    stats['sn_dist'].append(df['sn_avg_dist_to_center'].mean())
            
            # 配位数数据
            if 'coordination' in data:
                df = data['coordination']
                if 'Pt_cn_total' in df.columns:
                    stats['pt_cn_total'].append(df['Pt_cn_total'].mean())
                if 'Pt_cn_Pt_Pt' in df.columns:
                    stats['pt_pt_bonds'].append(df['Pt_cn_Pt_Pt'].mean())
                if 'Pt_cn_Pt_Sn' in df.columns:
                    stats['pt_sn_bonds'].append(df['Pt_cn_Pt_Sn'].mean())
                if 'Sn_cn_total' in df.columns:
                    stats['sn_cn_total'].append(df['Sn_cn_total'].mean())
                if 'Sn_cn_Sn_Sn' in df.columns:
                    stats['sn_sn_bonds'].append(df['Sn_cn_Sn_Sn'].mean())
                if 'Sn_cn_Sn_Pt' in df.columns:
                    stats['sn_pt_bonds'].append(df['Sn_cn_Sn_Pt'].mean())
        
        # 计算平均值
        result = {}
        for key, values in stats.items():
            if values:
                result[key] = np.mean(values)
                result[f'{key}_std'] = np.std(values)
            else:
                result[key] = 0
                result[f'{key}_std'] = 0
        
        return result
    
    def run_analysis(self):
        """运行完整分析"""
        # 1. 扫描数据结构
        self.scan_data_structure()
        
        if not self.systems:
            print("\n❌ 未找到任何数据!")
            return
        
        # 2. 分析每个系统
        print("\n" + "="*80)
        print("开始数据分析")
        print("="*80)
        
        for sys_key, sys_info in self.systems.items():
            temp_results = self.analyze_system(sys_key, sys_info)
            if temp_results:
                self.results[sys_key] = {
                    'info': sys_info,
                    'data': temp_results
                }
        
        # 3. 保存数据
        self.save_data_tables()
        
        # 4. 绘制图表
        self.plot_temperature_trends()
        self.plot_comparison_heatmaps()
        self.plot_op2_analysis()
        
        print("\n" + "="*80)
        print("✅ 分析完成!")
        print(f"结果保存在: {self.output_dir}")
        print("="*80)
    
    def save_data_tables(self):
        """保存数据表"""
        rows = []
        for sys_key, sys_data in self.results.items():
            sys_info = sys_data['info']
            for temp, stats in sys_data['data'].items():
                row = {
                    'system': sys_info['display_name'],
                    'system_key': sys_key,
                    'temperature': temp,
                    'n_runs': stats['n_runs'],
                    'q6': stats['q6'],
                    'q6_std': stats.get('q6_std', 0),
                    'q6_pt': stats.get('q6_pt', 0),
                    'q6_sn': stats.get('q6_sn', 0),
                    'op2_all_metal': stats.get('op2_all_metal', 0),
                    'op2_pt': stats.get('op2_pt', 0),
                    'op2_sn': stats.get('op2_sn', 0),
                    'rg': stats.get('rg', 0),
                    'pt_dist': stats.get('pt_dist', 0),
                    'sn_dist': stats.get('sn_dist', 0),
                    'pt_cn_total': stats.get('pt_cn_total', 0),
                    'pt_pt_bonds': stats.get('pt_pt_bonds', 0),
                    'pt_sn_bonds': stats.get('pt_sn_bonds', 0),
                    'sn_cn_total': stats.get('sn_cn_total', 0),
                    'sn_sn_bonds': stats.get('sn_sn_bonds', 0),
                    'sn_pt_bonds': stats.get('sn_pt_bonds', 0),
                }
                rows.append(row)
        
        df = pd.DataFrame(rows)
        csv_file = self.output_dir / 'all_systems_data.csv'
        df.to_csv(csv_file, index=False, encoding='utf-8-sig')
        print(f"\n💾 数据表已保存: {csv_file}")
        
        return df
    
    def plot_temperature_trends(self):
        """绘制温度趋势图"""
        n_systems = len(self.results)
        if n_systems == 0:
            return
        
        # 计算子图布局
        n_cols = min(3, n_systems)
        n_rows = (n_systems + n_cols - 1) // n_cols
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(6*n_cols, 5*n_rows))
        if n_systems == 1:
            axes = np.array([axes])
        axes = axes.flatten()
        
        for idx, (sys_key, sys_data) in enumerate(self.results.items()):
            ax = axes[idx]
            sys_info = sys_data['info']
            data = sys_data['data']
            
            temps = sorted(data.keys())
            q6_vals = [data[t]['q6'] for t in temps]
            op2_vals = [data[t].get('op2_all_metal', 0) for t in temps]
            
            ax.plot(temps, q6_vals, 'o-', label='Q6', linewidth=2, markersize=6)
            if any(v > 0 for v in op2_vals):
                ax2 = ax.twinx()
                ax2.plot(temps, op2_vals, 's--', label='OP2', color='red', linewidth=2, markersize=6)
                ax2.set_ylabel('OP2', color='red')
                ax2.tick_params(axis='y', labelcolor='red')
            
            ax.set_xlabel('温度 (K)')
            ax.set_ylabel('Q6', color='blue')
            ax.tick_params(axis='y', labelcolor='blue')
            ax.set_title(sys_info['display_name'])
            ax.grid(True, alpha=0.3)
        
        # 隐藏多余的子图
        for idx in range(len(self.results), len(axes)):
            axes[idx].axis('off')
        
        plt.suptitle('温度对Q6和OP2的影响', fontsize=14, fontweight='bold')
        plt.tight_layout()
        
        output_file = self.output_dir / 'temperature_trends.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"📊 温度趋势图已保存: {output_file}")
    
    def plot_comparison_heatmaps(self):
        """绘制对比热图"""
        if not self.results:
            return
        
        # 收集所有温度
        all_temps = set()
        for sys_data in self.results.values():
            all_temps.update(sys_data['data'].keys())
        temps = sorted(all_temps)
        
        systems = list(self.results.keys())
        
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        
        fields = [
            ('q6', 'Q6 六次对称性'),
            ('op2_all_metal', 'OP2 取向参数'),
            ('rg', '回转半径 Rg (Å)'),
            ('pt_pt_bonds', 'Pt-Pt键数'),
            ('pt_sn_bonds', 'Pt-Sn键数'),
            ('sn_sn_bonds', 'Sn-Sn键数'),
        ]
        
        for idx, (field, title) in enumerate(fields):
            ax = axes[idx // 3, idx % 3]
            
            # 创建矩阵
            matrix = []
            for temp in temps:
                row = []
                for sys_key in systems:
                    if temp in self.results[sys_key]['data']:
                        value = self.results[sys_key]['data'][temp].get(field, 0)
                        row.append(value if value else np.nan)
                    else:
                        row.append(np.nan)
                matrix.append(row)
            
            matrix = np.array(matrix)
            
            if np.all(np.isnan(matrix)) or np.nanmax(matrix) == 0:
                ax.text(0.5, 0.5, f'No {field} data', 
                       ha='center', va='center', fontsize=12, transform=ax.transAxes)
                ax.set_title(title, fontsize=11, fontweight='bold')
                ax.axis('off')
                continue
            
            im = ax.imshow(matrix, aspect='auto', cmap='RdYlBu_r', interpolation='nearest')
            
            # 设置标签
            y_labels = [f'{t}K' for t in temps]
            x_labels = [self.results[s]['info']['display_name'][:15] for s in systems]
            
            ax.set_yticks(range(len(temps)))
            ax.set_yticklabels(y_labels, fontsize=8)
            ax.set_xticks(range(len(systems)))
            ax.set_xticklabels(x_labels, fontsize=7, rotation=45, ha='right')
            ax.set_ylabel('温度')
            ax.set_title(title, fontsize=10, fontweight='bold')
            
            # 添加数值
            for i in range(len(temps)):
                for j in range(len(systems)):
                    if not np.isnan(matrix[i, j]):
                        ax.text(j, i, f'{matrix[i, j]:.2f}',
                               ha="center", va="center", color="black", fontsize=6)
            
            plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        
        plt.suptitle('多系统对比热图', fontsize=14, fontweight='bold')
        plt.tight_layout()
        
        output_file = self.output_dir / 'comparison_heatmaps.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"📊 对比热图已保存: {output_file}")
    
    def plot_op2_analysis(self):
        """绘制OP2专项分析图"""
        # 检查是否有OP2数据
        has_op2 = False
        for sys_data in self.results.values():
            for temp_data in sys_data['data'].values():
                if temp_data.get('op2_all_metal', 0) > 0:
                    has_op2 = True
                    break
            if has_op2:
                break
        
        if not has_op2:
            print("  (跳过OP2分析 - 无OP2数据)")
            return
        
        fig, axes = plt.subplots(2, 2, figsize=(14, 12))
        
        # 1. OP2 vs 温度 (所有系统)
        ax1 = axes[0, 0]
        for sys_key, sys_data in self.results.items():
            temps = sorted(sys_data['data'].keys())
            op2_vals = [sys_data['data'][t].get('op2_all_metal', 0) for t in temps]
            if any(v > 0 for v in op2_vals):
                ax1.plot(temps, op2_vals, 'o-', label=sys_data['info']['display_name'][:20])
        ax1.set_xlabel('温度 (K)')
        ax1.set_ylabel('OP2 (All Metal)')
        ax1.set_title('OP2 vs 温度')
        ax1.legend(fontsize=8)
        ax1.grid(True, alpha=0.3)
        
        # 2. OP2_Pt vs OP2_Sn
        ax2 = axes[0, 1]
        all_op2_pt = []
        all_op2_sn = []
        all_temps = []
        for sys_data in self.results.values():
            for temp, data in sys_data['data'].items():
                op2_pt = data.get('op2_pt', 0)
                op2_sn = data.get('op2_sn', 0)
                if op2_pt > 0 and op2_sn > 0:
                    all_op2_pt.append(op2_pt)
                    all_op2_sn.append(op2_sn)
                    all_temps.append(temp)
        
        if all_op2_pt and all_op2_sn:
            scatter = ax2.scatter(all_op2_pt, all_op2_sn, c=all_temps, cmap='coolwarm', s=50)
            plt.colorbar(scatter, ax=ax2, label='温度 (K)')
            ax2.plot([0, max(all_op2_pt)], [0, max(all_op2_pt)], 'k--', alpha=0.5, label='y=x')
            ax2.set_xlabel('OP2 (Pt)')
            ax2.set_ylabel('OP2 (Sn)')
            ax2.set_title('Pt vs Sn 取向参数对比')
            ax2.legend()
            ax2.grid(True, alpha=0.3)
        else:
            ax2.text(0.5, 0.5, 'No Pt/Sn OP2 data', ha='center', va='center', transform=ax2.transAxes)
        
        # 3. Q6 vs OP2 相关性
        ax3 = axes[1, 0]
        all_q6 = []
        all_op2 = []
        all_temps_q6 = []
        for sys_data in self.results.values():
            for temp, data in sys_data['data'].items():
                q6 = data.get('q6', 0)
                op2 = data.get('op2_all_metal', 0)
                if q6 > 0 and op2 > 0:
                    all_q6.append(q6)
                    all_op2.append(op2)
                    all_temps_q6.append(temp)
        
        if all_q6 and all_op2:
            scatter = ax3.scatter(all_q6, all_op2, c=all_temps_q6, cmap='coolwarm', s=50)
            plt.colorbar(scatter, ax=ax3, label='温度 (K)')
            # 计算相关系数
            from scipy.stats import pearsonr
            r, p = pearsonr(all_q6, all_op2)
            ax3.set_xlabel('Q6 六次对称性')
            ax3.set_ylabel('OP2 取向参数')
            ax3.set_title(f'Q6 vs OP2 相关性 (r={r:.3f}, p={p:.3e})')
            ax3.grid(True, alpha=0.3)
        else:
            ax3.text(0.5, 0.5, 'No Q6/OP2 data', ha='center', va='center', transform=ax3.transAxes)
        
        # 4. OP2物理意义说明
        ax4 = axes[1, 1]
        ax4.axis('off')
        explanation = """
OP2 (Orientation Parameter 2) 物理意义:

采用Steinhardt型序参量公式计算:
OP_l = √(4π/(2l+1) × Σ|q_lm|²)

其中 q_lm = (1/N_bond) × Σ w(r) × Y_lm(n̂)
权重函数: w(r) = exp(-0.5×((r-rc)/rc)²)

物理解释:
• OP2 高值 (>0.4): 2D平面/分散结构
  - 邻居分布在同一平面
  - 结构较为扁平

• OP2 低值 (<0.3): 3D紧凑/对称结构
  - 邻居分布在3D空间
  - 接近正四面体或正二十面体

温度效应:
• 升温 → OP2可能增加 (结构更分散)
• 降温 → OP2可能降低 (结构更紧凑)

元素差异:
• Pt通常倾向于形成紧凑结构
• Sn可能导致结构变化
        """
        ax4.text(0.1, 0.95, explanation, transform=ax4.transAxes, 
                fontsize=10, verticalalignment='top', fontfamily='monospace')
        
        plt.suptitle('OP2 取向参数详细分析', fontsize=14, fontweight='bold')
        plt.tight_layout()
        
        output_file = self.output_dir / 'op2_analysis.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"📊 OP2分析图已保存: {output_file}")


def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='分析新服务器数据')
    parser.add_argument('--data', type=str, 
                       default=r'C:\Users\11207\OneDrive\02_Code\work1-PtSnO\workflow\data\coordination\seletion\coordination_time_series_results_sample_20251130_193923',
                       help='数据根目录')
    parser.add_argument('--output', type=str, default=None,
                       help='输出目录')
    
    args = parser.parse_args()
    
    analyzer = NewServerDataAnalyzer(args.data, args.output)
    analyzer.run_analysis()


if __name__ == '__main__':
    main()
