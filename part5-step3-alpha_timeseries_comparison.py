#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Part 5 Step 3: Alpha时间序列分析与对比图

功能:
1. 批量处理XYZ轨迹文件，计算alpha时间序列
2. 生成类似OP2的对比图，对比不同温度的alpha演化
3. 支持系综平均和标准差计算

Author: AI Assistant
Date: 2025-12-28
"""

import argparse
import csv
import math
import sys
from pathlib import Path
from typing import List, Tuple, Iterator
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
import warnings

warnings.filterwarnings('ignore')

# 设置高质量论文图样式 - 与 step6_1_3_lindemann_only 保持一致
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
plt.rcParams['mathtext.fontset'] = 'dejavusans'
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['xtick.major.width'] = 1.5
plt.rcParams['ytick.major.width'] = 1.5
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# 字体大小常量 - 适配大图尺寸 (10x8 英寸)，与 step6_1_3 保持一致
FONT_TICK = 28
FONT_LABEL = 34
FONT_LEGEND = 26
FONT_ANNOT = 26  # inset 图内标注文字大小

# ==================== 工具函数 ====================

def convert_subscript(text: str) -> str:
    """
    将 _{数字} 格式转换为 matplotlib LaTeX 下标
    
    例如: "Pt_{8}Sn_{6}" -> "Pt$_8$Sn$_6$"
    
    这种格式可以被 matplotlib 正确渲染为下标
    """
    import re
    # 将 _{数字} 转换为 $_数字$
    result = re.sub(r'_\{(\d+)\}', r'$_{\1}$', text)
    
    return result

# ==================== CSV文件读取功能 ====================

def read_alpha_csv(csv_path: Path) -> Tuple[np.ndarray, np.ndarray]:
    """
    读取alpha_timeseries_*.csv文件
    
    Returns:
        frames: 帧编号数组
        alphas: alpha值数组
    """
    df = pd.read_csv(csv_path)
    
    # 支持两种列名格式
    if 'frame' in df.columns and 'alpha' in df.columns:
        frames = df['frame'].values
        alphas = df['alpha'].values
    elif 'Frame' in df.columns and 'Alpha' in df.columns:
        frames = df['Frame'].values
        alphas = df['Alpha'].values
    else:
        raise ValueError(f"CSV文件格式不正确: {csv_path}")
    
    return frames, alphas

def scan_csv_files(data_dir: Path, temperatures: List[int]) -> dict:
    """
    扫描目录中的alpha_timeseries_*.csv文件
    
    支持文件名格式:
    - alpha_timeseries_300K_r5.csv
    - alpha_timeseries_300K_run5.csv  
    - alpha_timeseries_300K_best2_r5.csv (带系列标识)
    
    Args:
        data_dir: CSV文件所在目录
        temperatures: 要查找的温度列表
    
    Returns:
        {temperature: [csv_paths]}
    """
    temp_data = {t: [] for t in temperatures}
    
    # 查找所有CSV文件
    for csv_file in data_dir.glob("alpha_timeseries_*.csv"):
        filename = csv_file.stem  # 去除.csv后缀
        
        for temp in temperatures:
            # 匹配 300K, 900K 等 (更宽松的匹配模式)
            if f"_{temp}K" in filename or f"_{temp}k" in filename:
                temp_data[temp].append(csv_file)
                print(f"    匹配: {csv_file.name} -> {temp}K")
                break
    
    return temp_data
plt.rcParams['axes.unicode_minus'] = False

# ==================== Alpha计算核心函数 ====================

COV_RAD_ANG = {
    "Pt": 1.36,
    "Sn": 1.39,
    "Ga": 1.22,
    "Pd": 1.39,
    "Ni": 1.24,
}

def rational_switch(r: float, r0: float, nn: int, nd: int, d0: float = 0.0) -> float:
    """PLUMED RATIONAL switching function"""
    if r0 <= 0.0:
        raise ValueError("r0 must be > 0")
    x = (r - d0) / r0
    if abs(x - 1.0) < 1e-10:
        return float(nn) / float(nd)
    try:
        xn = x ** nn
        xd = x ** nd
    except OverflowError:
        return 0.0
    denom = 1.0 - xd
    if denom == 0.0:
        return float(nn) / float(nd)
    return (1.0 - xn) / denom

def dist(a: Tuple[float,float,float], b: Tuple[float,float,float]) -> float:
    dx = a[0]-b[0]; dy = a[1]-b[1]; dz = a[2]-b[2]
    return math.sqrt(dx*dx + dy*dy + dz*dz)

def read_xyz_frames(path: str) -> Iterator[Tuple[List[str], List[Tuple[float,float,float]]]]:
    """读取XYZ轨迹的每一帧"""
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        while True:
            line = f.readline()
            if not line:
                return
            line = line.strip()
            if not line:
                continue
            try:
                n = int(line)
            except ValueError:
                continue
            _comment = f.readline()
            symbols: List[str] = []
            coords: List[Tuple[float,float,float]] = []
            for _ in range(n):
                parts = f.readline().split()
                if len(parts) < 4:
                    continue
                sym = parts[0]
                x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                symbols.append(sym)
                coords.append((x,y,z))
            yield symbols, coords

def cn_xy(
    X: List[Tuple[float,float,float]],
    Y: List[Tuple[float,float,float]],
    r0: float,
    nn: int,
    nd: int,
    same_set: bool = False
) -> float:
    """计算配位数CN(X-Y)"""
    nx = len(X)
    if nx == 0:
        return float("nan")
    total = 0.0
    if not same_set:
        for i in range(nx):
            ai = X[i]
            for j in range(len(Y)):
                total += rational_switch(dist(ai, Y[j]), r0, nn, nd)
    else:
        for i in range(nx):
            ai = X[i]
            for j in range(nx):
                if j == i:
                    continue
                total += rational_switch(dist(ai, Y[j]), r0, nn, nd)
    return total / float(nx)

def compute_alpha_timeseries(xyz_file, metal='Pt', second='Sn', nn=12, nd=24, 
                            r0_ms=None, r0_mm=None, auto_r0=True):
    """
    计算单个XYZ文件的alpha时间序列
    
    Returns:
        list of alpha values (one per frame)
    """
    if auto_r0:
        r0_ms = 1.3 * (COV_RAD_ANG[metal] + COV_RAD_ANG[second])
        r0_mm = 1.3 * (COV_RAD_ANG[metal] + COV_RAD_ANG[metal])
    
    if r0_ms is None or r0_mm is None:
        raise ValueError("Must specify r0_ms and r0_mm or use auto_r0=True")
    
    alphas = []
    
    for symbols, coords in read_xyz_frames(str(xyz_file)):
        pt = [coords[i] for i,s in enumerate(symbols) if s == metal]
        sn = [coords[i] for i,s in enumerate(symbols) if s == second]
        
        if len(pt) == 0 or len(sn) == 0:
            alphas.append(float('nan'))
            continue
        
        cn_ms = cn_xy(pt, sn, r0_ms, nn, nd, same_set=False)
        cn_mm = cn_xy(pt, pt, r0_mm, nn, nd, same_set=True)
        
        denom = (cn_ms + cn_mm)
        alpha = float("nan") if (not math.isfinite(denom) or denom == 0.0) else (1.0 - 2.0*cn_ms/denom)
        
        alphas.append(alpha)
    
    return alphas

# ==================== 数据扫描和处理 ====================

def scan_xyz_data(data_root, temp_list=[300, 900]):
    """
    扫描数据目录，查找不同温度的XYZ轨迹文件
    
    预期目录结构:
    data_root/
    ├── 300K/
    │   ├── traj_1.xyz
    │   ├── traj_2.xyz
    │   └── ...
    └── 900K/
        ├── traj_1.xyz
        └── ...
    
    Returns:
        dict: {temp: [list of xyz files]}
    """
    data_root = Path(data_root)
    temp_data = {temp: [] for temp in temp_list}
    
    for temp in temp_list:
        # 查找包含温度标识的目录
        for pattern in [f'{temp}K', f'T{temp}', f'{temp}']:
            temp_dirs = list(data_root.rglob(pattern))
            for temp_dir in temp_dirs:
                if temp_dir.is_dir():
                    # 查找XYZ文件
                    xyz_files = list(temp_dir.glob('*.xyz'))
                    xyz_files.extend(temp_dir.glob('**/*.xyz'))
                    
                    for xyz_file in xyz_files:
                        if xyz_file.is_file() and xyz_file.stat().st_size > 0:
                            temp_data[temp].append(xyz_file)
    
    return temp_data

def compute_ensemble_average(alpha_list_of_arrays, smooth=True):
    """
    计算alpha时间序列的系综平均
    
    Parameters:
        alpha_list_of_arrays: list of numpy arrays (each array is one trajectory)
        smooth: whether to apply Savitzky-Golay filter
    
    Returns:
        mean, std, n_trajectories
    """
    if not alpha_list_of_arrays:
        return None, None, 0
    
    # 找到最短的序列长度
    min_len = min(len(arr) for arr in alpha_list_of_arrays)
    
    # 截断所有序列到相同长度
    truncated = [arr[:min_len] for arr in alpha_list_of_arrays]
    
    # 转换为numpy数组
    values_array = np.array(truncated)
    
    # 计算均值和标准差
    mean = np.nanmean(values_array, axis=0)
    std = np.nanstd(values_array, axis=0)
    
    # 平滑处理
    if smooth and len(mean) > 21:
        mean = savgol_filter(mean, 21, 3)
    
    return mean, std, len(alpha_list_of_arrays)

# ==================== 绘图函数 ====================

def plot_alpha_comparison(temp1_data, temp2_data, temp1, temp2, output_dir, 
                         total_time_ps=175.0, system_name='Supported Pt₈Sn₆',
                         color_mode='gray', hide_title=False, figsize=(10, 8),
                         hide_temp_label=False):
    """
    绘制两个温度的alpha对比图
    类似OP2对比图的风格,上下堆叠布局
    
    Parameters:
        temp1_data: list of alpha arrays for temperature 1
        temp2_data: list of alpha arrays for temperature 2
        temp1, temp2: temperature values (K)
        output_dir: output directory
        total_time_ps: total simulation time in picoseconds
        system_name: system name for title (使用下标格式 Pt₈Sn₆, 不加粗)
        color_mode: 'gray' (黑灰色) 或 'color' (彩色: 300K蓝, 900K红)
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 计算系综平均
    mean1, std1, n1 = compute_ensemble_average(temp1_data)
    mean2, std2, n2 = compute_ensemble_average(temp2_data)
    
    if mean1 is None or mean2 is None:
        print(f"  ⚠️ Alpha数据缺失")
        return
    
    print(f"  📊 {temp1}K: {n1} 条轨迹, {len(mean1)} 帧")
    print(f"  📊 {temp2}K: {n2} 条轨迹, {len(mean2)} 帧")
    
    # 统一Y轴范围
    y_values = []
    y_values.extend(mean1 - std1)
    y_values.extend(mean1 + std1)
    y_values.extend(mean2 - std2)
    y_values.extend(mean2 + std2)
    y_min, y_max = np.nanmin(y_values), np.nanmax(y_values)
    y_margin = (y_max - y_min) * 0.1
    y_min -= y_margin
    y_max += y_margin
    
    # 时间轴
    n_frames1 = len(mean1)
    n_frames2 = len(mean2)
    time1 = np.linspace(0, total_time_ps, n_frames1)
    time2 = np.linspace(0, total_time_ps, n_frames2)
    
    # 创建图表 - 大图尺寸 (10x8 英寸)，无间隙紧贴布局
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, sharex=True,
                                   gridspec_kw={'hspace': 0})
    
    # 配色方案
    if color_mode == 'color':
        # 彩色模式: 300K蓝色, 900K红色
        color1 = '#1f77b4'  # 蓝色
        color2 = '#d62728'  # 红色
        line_color1 = color1
        line_color2 = color2
        fill_color1 = color1
        fill_color2 = color2
    else:
        # 黑灰模式 (默认, 与OP2保持一致)
        line_color1 = '#333333'  # 深灰/黑色
        line_color2 = '#333333'
        fill_color1 = '#888888'  # 中灰色
        fill_color2 = '#888888'
    
    # === 上图: temp1 (300K) ===
    ax1.fill_between(time1, mean1 - std1, mean1 + std1, 
                     alpha=0.3, color=fill_color1)
    ax1.plot(time1, mean1, color=line_color1, linewidth=4)
    
    ax1.set_ylim(y_min, y_max)
    ax1.set_ylabel('Alpha', fontsize=FONT_LABEL)
    ax1.tick_params(axis='y', labelsize=FONT_TICK, width=1.5, length=6)
    ax1.tick_params(axis='x', which='both', bottom=False, labelbottom=False)  # 隐藏上图X轴刻度

    # inset 标注：左上角系统名，右上角温度
    if not hide_title:
        ax1.text(0.03, 0.95, system_name, transform=ax1.transAxes,
                 fontsize=FONT_ANNOT, va='top', ha='left')
    if not hide_temp_label:
        ax1.text(0.97, 0.95, f'{temp1} K', transform=ax1.transAxes,
                 fontsize=FONT_ANNOT, va='top', ha='right')
    
    # === 下图: temp2 (900K) ===
    ax2.fill_between(time2, mean2 - std2, mean2 + std2, 
                     alpha=0.3, color=fill_color2)
    ax2.plot(time2, mean2, color=line_color2, linewidth=4)
    
    ax2.set_ylim(y_min, y_max)
    ax2.set_xlabel('Time (ps)', fontsize=FONT_LABEL)
    ax2.set_ylabel('Alpha', fontsize=FONT_LABEL)
    ax2.tick_params(axis='both', labelsize=FONT_TICK, width=1.5, length=6)

    # inset 标注：右上角温度
    if not hide_temp_label:
        ax2.text(0.97, 0.95, f'{temp2} K', transform=ax2.transAxes,
                 fontsize=FONT_ANNOT, va='top', ha='right')
    
    plt.tight_layout()
    plt.subplots_adjust(hspace=0)  # tight_layout后强制归零，消除上下子图间隙
    
    output_file = output_dir / f'alpha_comparison_{temp1}K_vs_{temp2}K.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"\n📊 Alpha对比图已保存: {output_file}")
    
    # 打印统计信息
    mean1_val = np.nanmean(mean1)
    mean2_val = np.nanmean(mean2)
    print(f"\n>>> Alpha 统计对比:")
    print(f"    {temp1}K: Alpha = {mean1_val:.4f} ± {np.nanmean(std1):.4f}")
    print(f"    {temp2}K: Alpha = {mean2_val:.4f} ± {np.nanmean(std2):.4f}")
    print(f"    Δ(Alpha) = {mean2_val - mean1_val:+.4f}")
    
    return output_file

# ==================== 主函数 ====================

def main():
    parser = argparse.ArgumentParser(
        description='Alpha时间序列分析与温度对比',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例用法:
  # 直接指定两个XYZ文件
  python part5-step3-alpha_timeseries_comparison.py --xyz1 path/to/300K.xyz --xyz2 path/to/900K.xyz
  
  # 指定温度标签
  python part5-step3-alpha_timeseries_comparison.py --xyz1 file1.xyz --xyz2 file2.xyz --temp1 300 --temp2 1000
  
  # 自动扫描目录（旧模式）
  python part5-step3-alpha_timeseries_comparison.py --data data/xyz_trajs --scan-mode
        """
    )
    
    parser.add_argument('--xyz1', type=str, default=None,
                       help='温度1的XYZ文件路径')
    parser.add_argument('--xyz2', type=str, default=None,
                       help='温度2的XYZ文件路径')
    
    parser.add_argument('--data', type=str, 
                       default=str(Path(__file__).parent / 'data' / 'xyz_trajs'),
                       help='数据根目录（scan-mode/csv-mode时使用）')
    parser.add_argument('--scan-mode', action='store_true',
                       help='扫描模式：自动在data目录中查找XYZ文件')
    parser.add_argument('--csv-mode', action='store_true',
                       help='CSV模式：读取已计算的alpha_timeseries_*.csv文件')
    
    parser.add_argument('--output', type=str,
                       default=str(Path(__file__).parent / 'results' / 'alpha_comparison'),
                       help='输出目录 (默认: results/alpha_comparison)')
    
    parser.add_argument('--temp1', type=int, default=300, help='温度1标签 (K, 默认: 300)')
    parser.add_argument('--temp2', type=int, default=900, help='温度2标签 (K, 默认: 900)')
    
    parser.add_argument('--metal', type=str, default='Pt', help='金属元素 (默认: Pt)')
    parser.add_argument('--second', type=str, default='Sn', help='第二元素 (默认: Sn)')
    
    parser.add_argument('--nn', type=int, default=12, help='NN指数 (默认: 12)')
    parser.add_argument('--nd', type=int, default=None, help='ND指数 (默认: 2*NN)')
    
    parser.add_argument('--total-time', type=float, default=175.0, 
                       help='总模拟时间(ps) (默认: 175.0)')
    parser.add_argument('--system-name', type=str, default='Supported Pt_{8}Sn_{6}', 
                       help='体系名称，支持 _{数字} 格式自动转为下标 (默认: Supported Pt_{8}Sn_{6})')
    
    parser.add_argument('--color-mode', type=str, default='gray', choices=['gray', 'color'],
                       help='配色模式: gray=黑灰色, color=彩色(300K蓝/900K红) (默认: gray)')
    
    parser.add_argument('--figsize', type=str, default='10x8',
                       help='图像尺寸，格式 WxH（英寸），默认: 10x8')
    parser.add_argument('--hide-title', action='store_true',
                       help='隐藏图内左上角系统名标注（温度标注保留）')
    parser.add_argument('--hide-temp-label', action='store_true',
                       help='隐藏图内右上角温度标注（300 K / 900 K）')
    
    parser.add_argument('--no-auto-r0', action='store_true',
                       help='禁用自动R0计算（需手动指定--r0-ms和--r0-mm）')
    parser.add_argument('--r0-ms', type=float, default=None, help='Metal-Second的R0 (Å)')
    parser.add_argument('--r0-mm', type=float, default=None, help='Metal-Metal的R0 (Å)')
    
    args = parser.parse_args()
    
    # 转换下标格式
    args.system_name = convert_subscript(args.system_name)
    
    # 解析 figsize
    try:
        fw, fh = [float(x) for x in args.figsize.lower().split('x')]
    except Exception:
        print(f"⚠️ --figsize 格式错误，应为 WxH（如 10x8），使用默认值 10x8")
        fw, fh = 10.0, 8.0
    
    nd = args.nd if args.nd is not None else 2 * args.nn
    auto_r0 = not args.no_auto_r0
    
    print(f"\n{'='*70}")
    print(f"Part 5 Step 3: Alpha时间序列对比分析")
    print(f"{'='*70}")
    print(f"体系名称: {args.system_name}")
    print(f"配色模式: {args.color_mode}")
    print(f"金属: {args.metal}, 第二元素: {args.second}")
    print(f"温度对比: {args.temp1}K vs {args.temp2}K")
    print(f"NN={args.nn}, ND={nd}")
    if auto_r0:
        r0_ms = 1.3 * (COV_RAD_ANG[args.metal] + COV_RAD_ANG[args.second])
        r0_mm = 1.3 * (COV_RAD_ANG[args.metal] + COV_RAD_ANG[args.metal])
        print(f"R0 (auto): Metal-Second={r0_ms:.4f} Å, Metal-Metal={r0_mm:.4f} Å")
    else:
        print(f"R0 (manual): Metal-Second={args.r0_ms:.4f} Å, Metal-Metal={args.r0_mm:.4f} Å")
    print(f"{'='*70}\n")
    
    # 判断使用哪种模式
    if args.xyz1 and args.xyz2:
        # 直接文件模式
        print(">>> 模式: 直接文件输入")
        print(f"  温度1 ({args.temp1}K): {args.xyz1}")
        print(f"  温度2 ({args.temp2}K): {args.xyz2}")
        
        xyz1_path = Path(args.xyz1)
        xyz2_path = Path(args.xyz2)
        
        if not xyz1_path.exists():
            print(f"\n❌ 错误: 文件不存在 - {args.xyz1}")
            return
        if not xyz2_path.exists():
            print(f"\n❌ 错误: 文件不存在 - {args.xyz2}")
            return
        
        temp_data = {
            args.temp1: [xyz1_path],
            args.temp2: [xyz2_path]
        }
        
    elif args.scan_mode:
        # 扫描目录模式 (XYZ文件)
        print(">>> 模式: 目录扫描 (XYZ)")
        print(f"  扫描目录: {args.data}")
        
        temp_data = scan_xyz_data(args.data, [args.temp1, args.temp2])
        
        n_temp1 = len(temp_data[args.temp1])
        n_temp2 = len(temp_data[args.temp2])
        
        print(f"  找到 {args.temp1}K: {n_temp1} 个XYZ文件")
        print(f"  找到 {args.temp2}K: {n_temp2} 个XYZ文件")
        
        if n_temp1 == 0 or n_temp2 == 0:
            print("\n❌ 错误: 某个温度没有找到XYZ文件")
            print(f"请检查数据目录: {args.data}")
            return
    
    elif args.csv_mode:
        # CSV模式 (读取预计算的alpha时间序列)
        print(">>> 模式: CSV读取 (预计算的Alpha时间序列)")
        print(f"  扫描目录: {args.data}")
        
        data_dir = Path(args.data)
        if not data_dir.exists():
            print(f"\n❌ 错误: 目录不存在 - {args.data}")
            return
        
        temp_data = scan_csv_files(data_dir, [args.temp1, args.temp2])
        
        n_temp1 = len(temp_data[args.temp1])
        n_temp2 = len(temp_data[args.temp2])
        
        print(f"  找到 {args.temp1}K: {n_temp1} 个CSV文件")
        print(f"  找到 {args.temp2}K: {n_temp2} 个CSV文件")
        
        if n_temp1 == 0 or n_temp2 == 0:
            print("\n❌ 错误: 某个温度没有找到CSV文件")
            print(f"请检查数据目录: {args.data}")
            print(f"期望文件名格式: alpha_timeseries_300K_*.csv")
            return
    
    else:
        print("\n❌ 错误: 请指定输入模式")
        print("  选项1: 使用 --xyz1 和 --xyz2 直接指定两个XYZ文件")
        print("  选项2: 使用 --scan-mode --data <目录> 扫描XYZ文件")
        print("  选项3: 使用 --csv-mode --data <目录> 读取CSV文件")
        return
    
    # 获取文件列表
    n_temp1 = len(temp_data[args.temp1])
    n_temp2 = len(temp_data[args.temp2])
    
    # 2. 计算/读取alpha时间序列
    if args.csv_mode:
        # CSV模式: 直接读取
        print(f"\n>>> 步骤2: 读取Alpha时间序列...")
        
        alpha_temp1 = []
        print(f"  读取 {args.temp1}K ({n_temp1} 个CSV文件)...")
        for i, csv_file in enumerate(temp_data[args.temp1], 1):
            print(f"    [{i}/{n_temp1}] {csv_file.name}...", end=' ')
            try:
                frames, alphas = read_alpha_csv(csv_file)
                alpha_temp1.append(alphas)
                print(f"✓ ({len(alphas)} 帧)")
            except Exception as e:
                print(f"✗ 错误: {e}")
        
        alpha_temp2 = []
        print(f"  读取 {args.temp2}K ({n_temp2} 个CSV文件)...")
        for i, csv_file in enumerate(temp_data[args.temp2], 1):
            print(f"    [{i}/{n_temp2}] {csv_file.name}...", end=' ')
            try:
                frames, alphas = read_alpha_csv(csv_file)
                alpha_temp2.append(alphas)
                print(f"✓ ({len(alphas)} 帧)")
            except Exception as e:
                print(f"✗ 错误: {e}")
    
    else:
        # XYZ模式: 计算alpha
        print(f"\n>>> 步骤2: 计算Alpha时间序列...")
        
        alpha_temp1 = []
        print(f"  处理 {args.temp1}K ({n_temp1} 个文件)...")
        for i, xyz_file in enumerate(temp_data[args.temp1], 1):
            print(f"    [{i}/{n_temp1}] {xyz_file.name}...", end=' ')
            try:
                alphas = compute_alpha_timeseries(
                    xyz_file, 
                    metal=args.metal, 
                    second=args.second,
                    nn=args.nn, 
                    nd=nd,
                    r0_ms=args.r0_ms,
                    r0_mm=args.r0_mm,
                    auto_r0=auto_r0
                )
                alpha_temp1.append(np.array(alphas))
                print(f"✓ ({len(alphas)} 帧)")
            except Exception as e:
                print(f"✗ 错误: {e}")
        
        alpha_temp2 = []
        print(f"  处理 {args.temp2}K ({n_temp2} 个文件)...")
        for i, xyz_file in enumerate(temp_data[args.temp2], 1):
            print(f"    [{i}/{n_temp2}] {xyz_file.name}...", end=' ')
            try:
                alphas = compute_alpha_timeseries(
                    xyz_file,
                    metal=args.metal,
                    second=args.second,
                    nn=args.nn,
                    nd=nd,
                    r0_ms=args.r0_ms,
                    r0_mm=args.r0_mm,
                    auto_r0=auto_r0
                )
                alpha_temp2.append(np.array(alphas))
                print(f"✓ ({len(alphas)} 帧)")
            except Exception as e:
                print(f"✗ 错误: {e}")
    
    # 3. 绘制对比图
    print(f"\n>>> 步骤3: 绘制Alpha对比图...")
    plot_alpha_comparison(
        alpha_temp1, alpha_temp2,
        args.temp1, args.temp2,
        args.output,
        total_time_ps=args.total_time,
        system_name=args.system_name,
        color_mode=args.color_mode,
        hide_title=args.hide_title,
        figsize=(fw, fh),
        hide_temp_label=args.hide_temp_label
    )
    
    print(f"\n{'='*70}")
    print(f"✅ 分析完成!")
    print(f"输出目录: {args.output}")
    print(f"{'='*70}\n")

if __name__ == '__main__':
    main()
