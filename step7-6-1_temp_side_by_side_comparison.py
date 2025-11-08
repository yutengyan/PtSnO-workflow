#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
温度对比分析 - 并排展示多个温度的Q6演化

专门用于对比分析,例如300K vs 900K

用法:
    python step7-6-1_temp_side_by_side_comparison.py \
        --series Pt8Snx \
        --temps 300K,900K \
        --systems pt8sn1-2-best,pt8sn2-1-best,pt8sn3-1-best,pt8sn4-1-best,pt8sn5-1-best,pt8sn6-1-best,pt8sn7-1-best,pt8sn8-1-best,pt8sn9-1-best,pt8sn10-2-best

Author: AI Assistant
Date: 2025-10-27
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager
from scipy.signal import savgol_filter
import argparse
from pathlib import Path
import warnings
from v625_data_locator import V625DataLocator

warnings.filterwarnings('ignore')

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False


class TempSideBySideAnalyzer:
    """温度并排对比分析器"""
    
    def __init__(self, base_path, output_dir):
        self.base_path = Path(base_path)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.locator = V625DataLocator(base_path)
    
    def load_q6_time_series(self, run_path, sys_name, temp):
        """
        加载Q6时间序列数据
        支持v625格式(300K)和v626格式(T300.r3.gpu0)
        v626会自动加载该温度的所有重复运行并返回列表
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
    
    def plot_side_by_side_comparison(self, series_name, run_paths, systems, temps):
        """
        并排对比图: 左边第一个温度,右边第二个温度
        
        3行(cluster_metal_q6, pt_q6, sn_q6) × 2列(temp1, temp2)
        """
        print(f"\n{'='*80}")
        print(f"温度并排对比: {temps[0]} vs {temps[1]}")
        print(f"体系数量: {len(systems)}")
        print(f"{'='*80}")
        
        fields = ['cluster_metal_q6_global', 'pt_q6_global', 'sn_q6_global']
        titles = ['整体金属团簇Q6 (主分析)', 'Pt原子Q6 (辅助)', 'Sn原子Q6 (辅助)']
        
        fig, axes = plt.subplots(3, 2, figsize=(16, 12))
        
        colors = plt.cm.tab10(np.linspace(0, 1, len(systems)))
        
        for row_idx, (field, title) in enumerate(zip(fields, titles)):
            for col_idx, temp in enumerate(temps):
                ax = axes[row_idx, col_idx]
                
                valid_systems = []
                
                for sys_idx, sys_name in enumerate(systems):
                    # 加载数据
                    all_dfs = []
                    for run_path in run_paths:
                        dfs_list = self.load_q6_time_series(run_path, sys_name, temp)
                        if dfs_list is not None:
                            # load_q6_time_series返回列表，需要展开
                            for df in dfs_list:
                                if field in df.columns:
                                    all_dfs.append(df)
                    
                    if not all_dfs:
                        continue
                    
                    valid_systems.append(sys_name)
                    
                    # 计算平均
                    if len(all_dfs) > 1:
                        avg_df = pd.concat(all_dfs).groupby('frame').mean().reset_index()
                    else:
                        avg_df = all_dfs[0]
                    
                    time_ps = avg_df['time_ps'].values
                    avg_values = avg_df[field].values
                    
                    # 提取Sn编号用于标签
                    sn_num = sys_name.split('sn')[1].split('-')[0]
                    
                    # 绘制
                    ax.plot(time_ps, avg_values, linewidth=1.5, 
                           color=colors[sys_idx], 
                           label=f'Sn{sn_num}', alpha=0.8)
                
                # 设置标题和标签
                if row_idx == 0:
                    ax.set_title(f'{temp}', fontsize=14, fontweight='bold')
                
                if col_idx == 0:
                    ax.set_ylabel(title, fontsize=11, fontweight='bold')
                
                if row_idx == 2:
                    ax.set_xlabel('时间 (ps)', fontsize=11)
                
                ax.grid(True, alpha=0.3)
                
                # 图例(只在右上角显示一次)
                if row_idx == 0 and col_idx == 1:
                    ax.legend(loc='upper right', fontsize=9, ncol=2,
                             title='Sn含量', title_fontsize=10)
                
                # 显示有效体系数
                ax.text(0.02, 0.98, f'n={len(valid_systems)}体系', 
                       transform=ax.transAxes, fontsize=9, 
                       verticalalignment='top',
                       bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        plt.suptitle(f'Pt8Snx系列温度对比: {temps[0]} vs {temps[1]} (Q6时间演化)', 
                    fontsize=16, fontweight='bold', y=0.995)
        plt.tight_layout()
        
        # 保存
        temp_str = f'{temps[0]}_vs_{temps[1]}'
        output_file = self.output_dir / f'q6_time_comparison_{temp_str}_Pt8Snx_all.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"📊 温度对比图已保存: {output_file}")
        print(f"✅ 包含{len(systems)}个体系的完整对比")


def main():
    parser = argparse.ArgumentParser(description='温度并排对比分析')
    parser.add_argument('--series', type=str, required=True,
                       choices=['Pt8Snx', 'PtxSn8-x', 'Pt6Snx'],
                       help='指定系列')
    parser.add_argument('--temps', type=str, required=True,
                       help='两个温度,用逗号分隔,例如: 300K,900K')
    parser.add_argument('--systems', type=str, required=True,
                       help='体系列表,用逗号分隔')
    
    args = parser.parse_args()
    
    # 解析参数
    temps = [t.strip() for t in args.temps.split(',')]
    systems = [s.strip() for s in args.systems.split(',')]
    
    if len(temps) != 2:
        print("错误: 必须指定恰好2个温度")
        return
    
    # 路径配置 (使用v626数据)
    base_path = Path(__file__).parent / 'data' / 'coordination' / 'coordination_time_series_results_sample_20251106_214943'
    output_dir = Path(__file__).parent / 'results' / 'step7.6_q6_time'
    
    # 初始化分析器
    analyzer = TempSideBySideAnalyzer(base_path, output_dir)
    
    # 获取运行路径
    if args.series == 'PtxSn8-x':
        run_paths = analyzer.locator.find_all_runs('Pt8')  # 假设主要是Pt8
    else:
        series_folder = 'Pt8' if args.series == 'Pt8Snx' else 'Pt6'
        run_paths = analyzer.locator.find_all_runs(series_folder)
    
    print(f"\n找到{len(run_paths)}个运行文件夹")
    
    # 生成对比图
    analyzer.plot_side_by_side_comparison(args.series, run_paths, systems, temps)
    
    print(f"\n{'='*80}")
    print("✅ 分析完成!")
    print(f"{'='*80}")


if __name__ == '__main__':
    main()
