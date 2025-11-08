#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step7-6-3: Q6统计对比可视化

对比不同体系在300K vs 900K的Q6均值和变异系数(CV)
生成柱状图和散点图，帮助识别温度稳定性
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

sys.path.append(str(Path(__file__).parent))
from v625_data_locator import V625DataLocator

plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False


class Q6StatsComparison:
    """Q6统计对比分析器"""
    
    def __init__(self, base_path, output_dir):
        self.locator = V625DataLocator(base_path)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def load_q6_time_series(self, run_path, sys_name, temp):
        """加载Q6时间序列数据（支持v625和v626）"""
        sys_path = run_path / sys_name
        if not sys_path.exists():
            return None
        
        # v625格式
        csv_path_v625 = sys_path / temp / 'cluster_global_q6_time_series.csv'
        if csv_path_v625.exists():
            try:
                df = pd.read_csv(csv_path_v625)
                return [df]
            except:
                return None
        
        # v626格式
        temp_value = temp.replace('K', '')
        temp_pattern = f"T{temp_value}.*"
        matching_dirs = sorted(sys_path.glob(temp_pattern))
        
        if matching_dirs:
            all_dfs = []
            for temp_dir in matching_dirs:
                csv_path = temp_dir / 'cluster_global_q6_time_series.csv'
                if csv_path.exists():
                    try:
                        df = pd.read_csv(csv_path)
                        all_dfs.append(df)
                    except:
                        pass
            if all_dfs:
                return all_dfs
        
        return None
    
    def calculate_stats(self, run_paths, system_name, temp, field='cluster_metal_q6_global'):
        """计算Q6统计量"""
        all_dfs = []
        for run_path in run_paths:
            dfs_list = self.load_q6_time_series(run_path, system_name, temp)
            if dfs_list is not None:
                all_dfs.extend(dfs_list)
        
        if not all_dfs:
            return None
        
        # 计算平均
        if len(all_dfs) > 1:
            avg_df = pd.concat(all_dfs).groupby('frame').mean().reset_index()
        else:
            avg_df = all_dfs[0]
        
        values = avg_df[field].values
        mean = np.mean(values)
        std = np.std(values)
        cv = std / mean if mean != 0 else 0
        
        return {
            'mean': mean,
            'std': std,
            'cv': cv,
            'n_runs': len(all_dfs)
        }
    
    def plot_comparison(self, series_name, run_paths, systems, temps):
        """
        生成Q6均值和CV的对比图
        
        包含4个子图:
        1. Q6均值对比 (柱状图)
        2. CV对比 (柱状图)
        3. 温度敏感性散点图 (ΔQ6 vs ΔCV)
        4. 稳定性评分散点图 (Q6 vs CV)
        """
        # 提取Sn编号
        sn_nums = []
        for sys_name in systems:
            import re
            match = re.search(r'sn(\d+)', sys_name)
            sn_num = int(match.group(1)) if match else 0
            sn_nums.append(sn_num)
        
        # 收集数据
        data = {
            'sn': sn_nums,
            'system': systems,
            f'mean_{temps[0]}': [],
            f'mean_{temps[1]}': [],
            f'cv_{temps[0]}': [],
            f'cv_{temps[1]}': [],
        }
        
        print(f"\n收集统计数据...")
        for sys_name in systems:
            # 300K统计
            stats_t1 = self.calculate_stats(run_paths, sys_name, temps[0])
            # 900K统计
            stats_t2 = self.calculate_stats(run_paths, sys_name, temps[1])
            
            if stats_t1 and stats_t2:
                data[f'mean_{temps[0]}'].append(stats_t1['mean'])
                data[f'mean_{temps[1]}'].append(stats_t2['mean'])
                data[f'cv_{temps[0]}'].append(stats_t1['cv'])
                data[f'cv_{temps[1]}'].append(stats_t2['cv'])
                print(f"  {sys_name}: Q6({temps[0]})={stats_t1['mean']:.4f}, CV={stats_t1['cv']:.3f} | Q6({temps[1]})={stats_t2['mean']:.4f}, CV={stats_t2['cv']:.3f}")
            else:
                data[f'mean_{temps[0]}'].append(np.nan)
                data[f'mean_{temps[1]}'].append(np.nan)
                data[f'cv_{temps[0]}'].append(np.nan)
                data[f'cv_{temps[1]}'].append(np.nan)
                print(f"  {sys_name}: 数据缺失")
        
        df = pd.DataFrame(data)
        
        # 计算温度差异
        df['delta_mean'] = df[f'mean_{temps[0]}'] - df[f'mean_{temps[1]}']
        df['delta_cv'] = df[f'cv_{temps[1]}'] - df[f'cv_{temps[0]}']
        
        # 创建图表
        fig = plt.figure(figsize=(16, 12))
        gs = fig.add_gridspec(3, 2, hspace=0.3, wspace=0.3)
        
        colors_t1 = 'steelblue'
        colors_t2 = 'orangered'
        
        # === 子图1: Q6均值对比 ===
        ax1 = fig.add_subplot(gs[0, 0])
        x = np.arange(len(df))
        width = 0.35
        
        ax1.bar(x - width/2, df[f'mean_{temps[0]}'], width, 
               label=f'{temps[0]}', color=colors_t1, alpha=0.8)
        ax1.bar(x + width/2, df[f'mean_{temps[1]}'], width, 
               label=f'{temps[1]}', color=colors_t2, alpha=0.8)
        
        ax1.set_xlabel('Sn含量', fontsize=12, fontweight='bold')
        ax1.set_ylabel('Q6均值', fontsize=12, fontweight='bold')
        ax1.set_title('Q6均值对比', fontsize=14, fontweight='bold')
        ax1.set_xticks(x)
        ax1.set_xticklabels([f'Sn{s}' for s in df['sn']])
        ax1.legend()
        ax1.grid(True, alpha=0.3, axis='y')
        
        # === 子图2: CV对比 ===
        ax2 = fig.add_subplot(gs[0, 1])
        
        ax2.bar(x - width/2, df[f'cv_{temps[0]}'], width, 
               label=f'{temps[0]}', color=colors_t1, alpha=0.8)
        ax2.bar(x + width/2, df[f'cv_{temps[1]}'], width, 
               label=f'{temps[1]}', color=colors_t2, alpha=0.8)
        
        # CV参考线
        ax2.axhline(y=0.1, color='green', linestyle='--', linewidth=1, alpha=0.5, label='CV=0.1 (稳定)')
        ax2.axhline(y=0.2, color='orange', linestyle='--', linewidth=1, alpha=0.5, label='CV=0.2 (临界)')
        
        ax2.set_xlabel('Sn含量', fontsize=12, fontweight='bold')
        ax2.set_ylabel('变异系数 (CV)', fontsize=12, fontweight='bold')
        ax2.set_title('变异系数(CV)对比', fontsize=14, fontweight='bold')
        ax2.set_xticks(x)
        ax2.set_xticklabels([f'Sn{s}' for s in df['sn']])
        ax2.legend()
        ax2.grid(True, alpha=0.3, axis='y')
        
        # === 子图3: 温度敏感性 (ΔQ6 vs ΔCV) ===
        ax3 = fig.add_subplot(gs[1, :])
        
        scatter = ax3.scatter(df['delta_mean'], df['delta_cv'], 
                             c=df['sn'], cmap='viridis', 
                             s=200, alpha=0.7, edgecolors='black', linewidth=1.5)
        
        # 添加Sn标签
        for i, row in df.iterrows():
            ax3.annotate(f"Sn{row['sn']}", 
                        (row['delta_mean'], row['delta_cv']),
                        fontsize=10, fontweight='bold',
                        ha='center', va='center')
        
        # 参考线
        ax3.axhline(y=0, color='gray', linestyle='-', linewidth=1, alpha=0.5)
        ax3.axvline(x=0, color='gray', linestyle='-', linewidth=1, alpha=0.5)
        
        ax3.set_xlabel(f'ΔQ6 (Q6@{temps[0]} - Q6@{temps[1]})', fontsize=12, fontweight='bold')
        ax3.set_ylabel(f'ΔCV (CV@{temps[1]} - CV@{temps[0]})', fontsize=12, fontweight='bold')
        ax3.set_title('温度敏感性分析', fontsize=14, fontweight='bold')
        ax3.grid(True, alpha=0.3)
        
        # 添加说明文本
        text = (
            "理想区域 (左下): ΔQ6小且ΔCV小 → 温度稳定\n"
            "敏感区域 (右上): ΔQ6大且ΔCV大 → 温度敏感"
        )
        ax3.text(0.02, 0.98, text, transform=ax3.transAxes,
                verticalalignment='top', bbox=dict(boxstyle='round', 
                facecolor='wheat', alpha=0.5), fontsize=10)
        
        cbar = plt.colorbar(scatter, ax=ax3)
        cbar.set_label('Sn含量', fontsize=10)
        
        # === 子图4: 稳定性评分散点图 ===
        ax4_1 = fig.add_subplot(gs[2, 0])
        ax4_2 = fig.add_subplot(gs[2, 1])
        
        # 300K: Q6 vs CV
        scatter1 = ax4_1.scatter(df[f'mean_{temps[0]}'], df[f'cv_{temps[0]}'], 
                                c=df['sn'], cmap='viridis', 
                                s=200, alpha=0.7, edgecolors='black', linewidth=1.5)
        
        for i, row in df.iterrows():
            ax4_1.annotate(f"Sn{row['sn']}", 
                          (row[f'mean_{temps[0]}'], row[f'cv_{temps[0]}']),
                          fontsize=10, fontweight='bold',
                          ha='center', va='center')
        
        ax4_1.axhline(y=0.1, color='green', linestyle='--', linewidth=1, alpha=0.5)
        ax4_1.axhline(y=0.2, color='orange', linestyle='--', linewidth=1, alpha=0.5)
        ax4_1.set_xlabel('Q6均值', fontsize=12, fontweight='bold')
        ax4_1.set_ylabel('CV', fontsize=12, fontweight='bold')
        ax4_1.set_title(f'{temps[0]}稳定性评分', fontsize=14, fontweight='bold')
        ax4_1.grid(True, alpha=0.3)
        ax4_1.text(0.02, 0.98, '理想: 高Q6+低CV', transform=ax4_1.transAxes,
                  verticalalignment='top', bbox=dict(boxstyle='round', 
                  facecolor='lightgreen', alpha=0.5), fontsize=9)
        
        # 900K: Q6 vs CV
        scatter2 = ax4_2.scatter(df[f'mean_{temps[1]}'], df[f'cv_{temps[1]}'], 
                                c=df['sn'], cmap='viridis', 
                                s=200, alpha=0.7, edgecolors='black', linewidth=1.5)
        
        for i, row in df.iterrows():
            ax4_2.annotate(f"Sn{row['sn']}", 
                          (row[f'mean_{temps[1]}'], row[f'cv_{temps[1]}']),
                          fontsize=10, fontweight='bold',
                          ha='center', va='center')
        
        ax4_2.axhline(y=0.1, color='green', linestyle='--', linewidth=1, alpha=0.5)
        ax4_2.axhline(y=0.2, color='orange', linestyle='--', linewidth=1, alpha=0.5)
        ax4_2.set_xlabel('Q6均值', fontsize=12, fontweight='bold')
        ax4_2.set_ylabel('CV', fontsize=12, fontweight='bold')
        ax4_2.set_title(f'{temps[1]}稳定性评分', fontsize=14, fontweight='bold')
        ax4_2.grid(True, alpha=0.3)
        ax4_2.text(0.02, 0.98, '理想: 高Q6+低CV', transform=ax4_2.transAxes,
                  verticalalignment='top', bbox=dict(boxstyle='round', 
                  facecolor='lightcoral', alpha=0.5), fontsize=9)
        
        plt.suptitle(f'{series_name}系列 Q6统计对比: {temps[0]} vs {temps[1]}',
                    fontsize=16, fontweight='bold')
        
        # 保存图表
        output_file = self.output_dir / f'q6_stats_comparison_{temps[0]}_vs_{temps[1]}.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"\n📊 统计对比图已保存: {output_file.name}")
        
        # 保存CSV
        csv_file = self.output_dir / f'q6_stats_comparison_{temps[0]}_vs_{temps[1]}.csv'
        df.to_csv(csv_file, index=False, float_format='%.6f')
        print(f"📄 统计数据已保存: {csv_file.name}")
        
        return df


def main():
    parser = argparse.ArgumentParser(
        description='Q6均值和CV的可视化对比分析',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例用法:

python step7-6-3_q6_stats_comparison.py --series Pt8 --temps "300K,900K" --systems "pt8sn1-2-best,pt8sn2-1-best,pt8sn3-1-best,pt8sn4-1-best,pt8sn5-1-best,pt8sn6-1-best,pt8sn7-1-best,pt8sn8-1-best,pt8sn9-1-best,pt8sn10-2-best"
        """
    )
    
    parser.add_argument('--series', type=str, required=True,
                       help='系列名称(如Pt8)')
    parser.add_argument('--temps', type=str, required=True,
                       help='温度对比,逗号分隔(如"300K,900K")')
    parser.add_argument('--systems', type=str, required=True,
                       help='体系列表,逗号分隔')
    
    args = parser.parse_args()
    
    # 自动修正为 workflow 目录下统一保存
    base_path = Path(__file__).parent / 'data' / 'coordination' / 'coordination_time_series_results_sample_20251106_214943'
    output_dir = Path(__file__).parent / 'results' / 'step7.6_q6_stats'
    
    temps = [t.strip() for t in args.temps.split(',')]
    systems = [s.strip() for s in args.systems.split(',')]
    
    if len(temps) != 2:
        print("❌ 错误: 必须指定恰好2个温度")
        sys.exit(1)
    
    print(f"\n{'='*80}")
    print(f"Q6统计对比分析")
    print(f"{'='*80}")
    print(f"系列: {args.series}")
    print(f"温度对比: {temps[0]} vs {temps[1]}")
    print(f"体系数量: {len(systems)}")
    print(f"{'='*80}")
    
    analyzer = Q6StatsComparison(base_path, output_dir)
    
    run_paths = analyzer.locator.find_all_runs(args.series)
    if not run_paths:
        print(f"❌ 错误: 未找到{args.series}的运行数据")
        sys.exit(1)
    
    print(f"找到{len(run_paths)}个运行文件夹")
    
    df = analyzer.plot_comparison(args.series, run_paths, systems, temps)
    
    print(f"\n{'='*80}")
    print(f"✅ 分析完成!")
    print(f"{'='*80}")
    print(f"\n💡 图表说明:")
    print(f"  1. Q6均值对比: 高温下Q6通常降低")
    print(f"  2. CV对比: 高温下CV通常增大(波动增强)")
    print(f"  3. 温度敏感性: ΔQ6和ΔCV都小的体系温度稳定")
    print(f"  4. 稳定性评分: 理想体系应高Q6+低CV")


if __name__ == '__main__':
    main()
