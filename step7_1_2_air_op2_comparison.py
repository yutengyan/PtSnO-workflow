#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 7.1.2: Air-68 vs Air-86 OP2 系综平均对比图

绘制一张简洁的对比图：
- 上: Air-86 (Pt8Sn6) OP2
- 下: Air-68 (Pt6Sn8) OP2
- 统一Y轴便于对比

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

# 配置字体
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False


def scan_air_data(data_root):
    """扫描并分类Air数据"""
    data_root = Path(data_root)
    air68_data = {}
    air86_data = {}
    
    for root, dirs, files in os.walk(data_root):
        root_path = Path(root)
        # 尝试两种可能的文件名
        op2_file = root_path / 'cluster_op2_time_series.csv'
        if not op2_file.exists():
            op2_file = root_path / 'op2_time_series.csv'
        
        if op2_file.exists():
            path_str = str(root_path)
            
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
            
            try:
                df = pd.read_csv(op2_file)
                if is_68:
                    if temp not in air68_data:
                        air68_data[temp] = []
                    air68_data[temp].append(df)
                else:
                    if temp not in air86_data:
                        air86_data[temp] = []
                    air86_data[temp].append(df)
            except Exception as e:
                pass
    
    return air68_data, air86_data


def compute_ensemble_average(data_list, field='op2_all_metal'):
    """计算系综平均"""
    all_values = []
    
    for df in data_list:
        if df is not None and field in df.columns:
            all_values.append(df[field].values)
    
    if not all_values:
        return None, None, 0
    
    min_len = min(len(v) for v in all_values)
    all_values = [v[:min_len] for v in all_values]
    
    values_array = np.array(all_values)
    mean = np.mean(values_array, axis=0)
    std = np.std(values_array, axis=0)
    
    # 平滑
    if len(mean) > 21:
        mean = savgol_filter(mean, 21, 3)
    
    return mean, std, len(all_values)


def plot_op2_comparison(air68_data, air86_data, temp, output_dir):
    """
    绘制OP2对比图
    上: Air-86 (Pt8Sn6)
    下: Air-68 (Pt6Sn8)
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if temp not in air68_data or temp not in air86_data:
        print(f"  ⚠️ 温度 {temp}K 数据不完整")
        return
    
    # 计算系综平均
    mean_68, std_68, n_68 = compute_ensemble_average(air68_data[temp])
    mean_86, std_86, n_86 = compute_ensemble_average(air86_data[temp])
    
    if mean_68 is None or mean_86 is None:
        print(f"  ⚠️ OP2数据缺失")
        return
    
    # 统一Y轴范围
    y_values = []
    y_values.extend(mean_68 - std_68)
    y_values.extend(mean_68 + std_68)
    y_values.extend(mean_86 - std_86)
    y_values.extend(mean_86 + std_86)
    y_min, y_max = np.min(y_values), np.max(y_values)
    y_margin = (y_max - y_min) * 0.1
    y_min -= y_margin
    y_max += y_margin
    
    # 时间轴：175 ps
    total_time_ps = 175.0
    n_frames_86 = len(mean_86)
    n_frames_68 = len(mean_68)
    time_86 = np.linspace(0, total_time_ps, n_frames_86)
    time_68 = np.linspace(0, total_time_ps, n_frames_68)
    
    # 创建图表 - 科研绘图尺寸 (11.5cm x 9cm)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11.5/2.54, 9/2.54), sharex=True)
    
    # 黑灰配色
    line_color = '#333333'  # 深灰/黑色
    fill_color = '#888888'  # 中灰色
    
    # === 上图: Air-86 ===
    ax1.fill_between(time_86, mean_86 - std_86, mean_86 + std_86, 
                     alpha=0.3, color=fill_color)
    ax1.plot(time_86, mean_86, color=line_color, linewidth=1.5)
    
    ax1.set_ylim(y_min, y_max)
    ax1.set_ylabel('OP2', fontsize=10, fontweight='bold')
    ax1.set_title(r'Pt$_8$Sn$_6$ @ ' + f'{temp} K', fontsize=10, fontweight='bold')
    ax1.tick_params(axis='both', which='major', labelsize=9)
    
    # === 下图: Air-68 ===
    ax2.fill_between(time_68, mean_68 - std_68, mean_68 + std_68, 
                     alpha=0.3, color=fill_color)
    ax2.plot(time_68, mean_68, color=line_color, linewidth=1.5)
    
    ax2.set_ylim(y_min, y_max)
    ax2.set_xlabel('Time (ps)', fontsize=10, fontweight='bold')
    ax2.set_ylabel('OP2', fontsize=10, fontweight='bold')
    ax2.set_title(r'Pt$_6$Sn$_8$ @ ' + f'{temp} K', fontsize=10, fontweight='bold')
    ax2.tick_params(axis='both', which='major', labelsize=9)
    
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.3)
    
    output_file = output_dir / f'air_op2_comparison_{temp}K.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"📊 OP2对比图已保存: {output_file}")
    
    # 打印对比
    mean_86_val = np.mean(mean_86)
    mean_68_val = np.mean(mean_68)
    print(f"\n>>> OP2 统计对比 @ {temp}K:")
    print(f"    Pt8Sn6: OP2 = {mean_86_val:.4f}")
    print(f"    Pt6Sn8: OP2 = {mean_68_val:.4f}")
    print(f"    Δ = {mean_86_val - mean_68_val:+.4f}")
    
    return output_file


def main():
    parser = argparse.ArgumentParser(description='Air-68 vs Air-86 OP2对比图')
    
    parser.add_argument('--data', type=str, 
                       default=str(Path(__file__).parent / 'data' / 'coordination' / 'air' / 
                                  'coordination_time_series_results_air-sample_20251130_211818'),
                       help='数据目录')
    parser.add_argument('--output', type=str, 
                       default=str(Path(__file__).parent / 'results' / 'step7.1.2_air_comparison'),
                       help='输出目录')
    parser.add_argument('--temp', type=int, default=300, help='温度 (默认: 300)')
    
    args = parser.parse_args()
    
    print(f"\n{'='*60}")
    print(f"Step 7.1.2: Air OP2 对比图")
    print(f"温度: {args.temp}K")
    print(f"{'='*60}")
    
    print("\n>>> 扫描数据...")
    air68_data, air86_data = scan_air_data(args.data)
    
    print(f"    Air-68: {len(air68_data.get(args.temp, []))} 个结构")
    print(f"    Air-86: {len(air86_data.get(args.temp, []))} 个结构")
    
    plot_op2_comparison(air68_data, air86_data, args.temp, args.output)
    
    print(f"\n{'='*60}")
    print(f"✅ 完成!")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
