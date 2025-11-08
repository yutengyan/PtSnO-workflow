#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step7-6-2: 单个体系的温度对比分析

为每个体系单独生成300K vs 900K对比图
每张图包含:
- 3行(cluster_metal_q6, pt_q6, sn_q6) × 2列(300K, 900K)
- 统计信息盒子(均值、标准差、变异系数CV)
- 多运行平均+原始数据展示
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.signal import savgol_filter

# 添加父目录到路径
sys.path.append(str(Path(__file__).parent))
from v625_data_locator import V625DataLocator

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['font.size'] = 10
# 确保图例也使用正确的字体
from matplotlib.font_manager import FontProperties
try:
    chinese_font = FontProperties(family='Microsoft YaHei', size=9)
except:
    try:
        chinese_font = FontProperties(family='SimHei', size=9)
    except:
        chinese_font = FontProperties(size=9)


class IndividualSystemTempAnalyzer:
    """单个体系温度对比分析器"""
    
    def __init__(self, base_path, output_dir, unified_ylim=False):
        self.locator = V625DataLocator(base_path)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.unified_ylim = unified_ylim
        
        # 设置中文字体
        try:
            self.chinese_font = FontProperties(family='Microsoft YaHei', size=9)
        except:
            try:
                self.chinese_font = FontProperties(family='SimHei', size=9)
            except:
                self.chinese_font = FontProperties(size=9)
    
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
                return [df]  # 返回列表格式以统一处理
            except Exception as e:
                print(f"⚠️ 加载失败 {csv_path_v625}: {e}")
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
                        all_dfs.append(df)
                    except Exception as e:
                        print(f"⚠️ 加载失败 {csv_path}: {e}")
            
            if all_dfs:
                return all_dfs  # 返回所有运行的列表
        
        return None
    
    def calculate_statistics(self, values):
        """
        计算统计指标
        
        变异系数(Coefficient of Variation, CV):
        - 定义: CV = 标准差 / 均值
        - 含义: 衡量数据相对离散程度
        - 特点: 无量纲,便于比较不同量级的数据
        - 解读:
          * CV < 0.1: 波动很小,非常稳定
          * 0.1 < CV < 0.2: 波动适中,较稳定
          * CV > 0.2: 波动较大,不稳定
        """
        mean = np.mean(values)
        std = np.std(values)
        cv = std / mean if mean != 0 else 0
        
        # 线性趋势拟合
        x = np.arange(len(values))
        slope, intercept = np.polyfit(x, values, 1)
        
        return {
            'mean': mean,
            'std': std,
            'cv': cv,  # 变异系数
            'min': np.min(values),
            'max': np.max(values),
            'slope': slope,
        }
    
    def plot_single_system_comparison(self, series_name, run_paths, system_name, temps):
        """
        为单个体系绘制温度对比图
        
        Args:
            series_name: 系列名称(如Pt8Snx)
            run_paths: 运行路径列表
            system_name: 体系名称(如pt8sn2-1-best)
            temps: 温度列表[temp1, temp2]
        """
        # 提取Sn含量
        import re
        match = re.search(r'sn(\d+)', system_name)
        sn_num = int(match.group(1)) if match else 0
        
        # 创建图表
        fig, axes = plt.subplots(3, 2, figsize=(14, 10))
        fig.suptitle(f'{series_name} {system_name} 温度对比: {temps[0]} vs {temps[1]}',
                     fontsize=16, fontweight='bold', y=0.995)
        
        # Q6字段
        fields = ['cluster_metal_q6_global', 'pt_q6_global', 'sn_q6_global']
        field_titles = {
            'cluster_metal_q6_global': 'Cluster Metal Q6 (主要指标)',
            'pt_q6_global': 'Pt Q6 (辅助)',
            'sn_q6_global': 'Sn Q6 (辅助)'
        }
        
        # 如果需要统一Y轴，先收集所有数据的范围
        y_limits = {}
        if self.unified_ylim:
            for field in fields:
                all_values = []
                for temp in temps:
                    for run_path in run_paths:
                        dfs_list = self.load_q6_time_series(run_path, system_name, temp)
                        if dfs_list is not None:
                            for df in dfs_list:
                                if field in df.columns:
                                    all_values.extend(df[field].values)
                
                if all_values:
                    y_min = np.min(all_values)
                    y_max = np.max(all_values)
                    y_range = y_max - y_min
                    # 添加5%的边距
                    y_limits[field] = (y_min - 0.05*y_range, y_max + 0.05*y_range)
        
        # 为每个Q6字段和温度绘图
        for row_idx, field in enumerate(fields):
            for col_idx, temp in enumerate(temps):
                ax = axes[row_idx, col_idx]
                
                # 加载所有运行的数据
                all_dfs = []
                for run_path in run_paths:
                    dfs_list = self.load_q6_time_series(run_path, system_name, temp)
                    if dfs_list is not None:
                        # load_q6_time_series返回列表，需要展开
                        all_dfs.extend(dfs_list)
                
                if not all_dfs:
                    ax.text(0.5, 0.5, 'No Data', 
                           ha='center', va='center', transform=ax.transAxes)
                    continue
                
                # 绘制原始数据(半透明)
                for df in all_dfs:
                    time_ps = df['time_ps'].values
                    values = df[field].values
                    ax.plot(time_ps, values, color='gray', alpha=0.2, linewidth=0.5)
                
                # 计算平均值
                if len(all_dfs) > 1:
                    # 多运行平均
                    avg_df = pd.concat(all_dfs).groupby('frame').mean().reset_index()
                else:
                    avg_df = all_dfs[0]
                
                time_ps = avg_df['time_ps'].values
                avg_values = avg_df[field].values
                
                # 绘制平均曲线(粗线)
                color = 'blue' if col_idx == 0 else 'red'
                ax.plot(time_ps, avg_values, color=color, linewidth=2.5, 
                       label=f'平均 (n={len(all_dfs)})', alpha=0.8)
                
                # 平滑曲线
                if len(avg_values) > 51:
                    smoothed = savgol_filter(avg_values, window_length=51, polyorder=3)
                    ax.plot(time_ps, smoothed, color=color, linewidth=1.5, 
                           linestyle='--', label='平滑', alpha=0.6)
                
                # 计算统计信息
                stats = self.calculate_statistics(avg_values)
                
                # 添加统计信息盒子
                stats_text = (
                    f"均值: {stats['mean']:.4f}\n"
                    f"标准差: {stats['std']:.4f}\n"
                    f"CV: {stats['cv']:.3f}\n"
                    f"范围: [{stats['min']:.4f}, {stats['max']:.4f}]"
                )
                
                # 根据CV值选择背景色
                if stats['cv'] < 0.1:
                    box_color = '#e8f5e9'  # 浅绿色 - 很稳定
                elif stats['cv'] < 0.2:
                    box_color = '#fff9c4'  # 浅黄色 - 较稳定
                else:
                    box_color = '#ffebee'  # 浅红色 - 不稳定
                
                ax.text(0.02, 0.98, stats_text,
                       transform=ax.transAxes,
                       verticalalignment='top',
                       bbox=dict(boxstyle='round', facecolor=box_color, alpha=0.8),
                       fontproperties=self.chinese_font)
                
                # 设置标题
                if row_idx == 0:
                    ax.set_title(f'{temp}', fontsize=12, fontweight='bold')
                
                # 设置Y轴标签
                if col_idx == 0:
                    ax.set_ylabel(field_titles[field], fontsize=10, fontweight='bold')
                
                # 设置X轴标签
                if row_idx == 2:
                    ax.set_xlabel('Time (ps)', fontsize=10)
                
                # 网格和图例
                ax.grid(True, alpha=0.3, linestyle='--')
                ax.legend(loc='upper right', prop=self.chinese_font, framealpha=0.9)
                
                # 如果启用统一Y轴，应用统一的范围
                if self.unified_ylim and field in y_limits:
                    ax.set_ylim(y_limits[field])
        
        plt.tight_layout()
        
        # 保存图表
        output_file = self.output_dir / f'q6_comparison_{system_name}_{temps[0]}_vs_{temps[1]}.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"📊 {system_name} 温度对比图已保存: {output_file.name}")
        
        return output_file


def main():
    parser = argparse.ArgumentParser(
        description='单个体系的温度对比分析(每个体系一张图)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例用法:

1. 为Pt8Snx所有体系生成300K vs 900K对比图:
   python step7-6-2_individual_system_temp_comparison.py --series Pt8Snx --temps "300K,900K" --systems "pt8sn1-2-best,pt8sn2-1-best,pt8sn3-1-best,pt8sn4-1-best,pt8sn5-1-best,pt8sn6-1-best,pt8sn7-1-best,pt8sn8-1-best,pt8sn9-1-best,pt8sn10-2-best"

2. 使用统一Y轴范围(便于对比左右两列):
   python step7-6-2_individual_system_temp_comparison.py --series Pt8Snx --temps "300K,900K" --systems "pt8sn2-1-best,pt8sn4-1-best" --unified-ylim

3. 只对比几个关键体系:
   python step7-6-2_individual_system_temp_comparison.py --series Pt8Snx --temps "300K,900K" --systems "pt8sn2-1-best,pt8sn4-1-best,pt8sn6-1-best"

4. 其他温度对比:
   python step7-6-2_individual_system_temp_comparison.py --series Pt8Snx --temps "200K,600K" --systems "pt8sn2-1-best"

变异系数(CV)解读:
- CV < 0.1  (绿色背景): 波动很小,非常稳定
- 0.1 < CV < 0.2 (黄色背景): 波动适中,较稳定  
- CV > 0.2  (红色背景): 波动较大,不稳定

Y轴设置:
- 默认: 每个子图自动调整Y轴范围(更清晰显示各自细节)
- --unified-ylim: 同一行的两列使用相同Y轴(便于直接对比温度差异)
        """
    )
    
    parser.add_argument('--series', type=str, required=True,
                       help='系列名称(如Pt8Snx)')
    parser.add_argument('--temps', type=str, required=True,
                       help='温度对比,逗号分隔(如"300K,900K")')
    parser.add_argument('--systems', type=str, required=True,
                       help='体系列表,逗号分隔(如"pt8sn2-1-best,pt8sn4-1-best")')
    parser.add_argument('--base-path', type=str,
                       default=r'D:\OneDrive\py\Cv\lin\MSD_Analysis_Collection\v3_simplified_workflow\files\q6_cn\v626\coordination_time_series_results_sample_20251026_200908',
                       help='数据根目录(默认v626)')
    parser.add_argument('--output-dir', type=str,
                       default=r'D:\OneDrive\py\Cv\lin\MSD_Analysis_Collection\v3_simplified_workflow\results\step7.6_q6_time',
                       help='输出目录')
    parser.add_argument('--unified-ylim', action='store_true',
                       help='统一Y轴范围(每行的左右两列使用相同Y轴)')
    
    args = parser.parse_args()
    
    # 解析参数
    temps = [t.strip() for t in args.temps.split(',')]
    systems = [s.strip() for s in args.systems.split(',')]
    
    if len(temps) != 2:
        print("❌ 错误: 必须指定恰好2个温度")
        sys.exit(1)
    
    print(f"\n{'='*80}")
    print(f"单个体系温度对比分析")
    print(f"{'='*80}")
    print(f"系列: {args.series}")
    print(f"温度对比: {temps[0]} vs {temps[1]}")
    print(f"体系数量: {len(systems)}")
    print(f"统一Y轴: {'是' if args.unified_ylim else '否'}")
    print(f"{'='*80}\n")
    
    # 创建分析器
    # 自动修正为 workflow 目录下统一保存
    base_path = Path(__file__).parent / 'data' / 'coordination' / 'coordination_time_series_results_sample_20251106_214943'
    output_dir = Path(__file__).parent / 'results' / 'step7.6_individual_system'
    analyzer = IndividualSystemTempAnalyzer(base_path, output_dir, unified_ylim=args.unified_ylim)
    
    # 获取运行路径
    run_paths = analyzer.locator.find_all_runs(args.series)
    if not run_paths:
        print(f"❌ 错误: 未找到{args.series}的运行数据")
        sys.exit(1)
    
    print(f"找到{len(run_paths)}个运行文件夹\n")
    
    # 为每个体系生成对比图
    success_count = 0
    for system_name in systems:
        try:
            analyzer.plot_single_system_comparison(
                args.series, run_paths, system_name, temps
            )
            success_count += 1
        except Exception as e:
            print(f"❌ {system_name} 分析失败: {e}")
    
    print(f"\n{'='*80}")
    print(f"✅ 完成! 成功生成 {success_count}/{len(systems)} 个体系的对比图")
    print(f"{'='*80}")
    print(f"\n📁 输出目录: {args.output_dir}")
    print(f"📊 文件命名: q6_comparison_<体系名>_{temps[0]}_vs_{temps[1]}.png")
    
    print("\n💡 变异系数(CV)说明:")
    print("  - CV = 标准差 / 均值 (无量纲)")
    print("  - CV < 0.1  → 非常稳定 (绿色背景)")
    print("  - 0.1 < CV < 0.2 → 较稳定 (黄色背景)")
    print("  - CV > 0.2  → 不稳定 (红色背景)")
    print("  - CV越小说明Q6波动越小,结构越稳定")


if __name__ == '__main__':
    main()
