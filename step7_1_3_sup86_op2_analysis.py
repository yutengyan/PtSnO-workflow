#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 7.1.3: Supported Pt8Sn6 OP2 分析

数据来源:
  data/coordination/sup/coordination_time_series_results_sample_20260310_191856/
    4090-ustc/more/run3/Pt8/pt8sn6-1-best/   → 30 runs (sample模式, 每温度3个run)
    dp-md/GPU-Pt8/Pt8{,-2,-3}/pt8sn6-1-best/  → 27 runs (dp-md原始)
    dp-md/more/Pt8{,-2,-3,-4}/pt8sn6-1-best/  → 40 runs (more补充)

功能:
  1. --mode scan   : OP2 均值 vs 温度 散点/折线图 (所有来源系综平均)
  2. --mode ts     : 指定温度的 OP2 时间序列 (可叠加多个来源)
  3. --mode compare: 两个温度的时间序列上下对比 (与 step7_1_2 一致的双子图布局)

Author: AI Assistant
Date: 2026-03-10
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
import os
import sys
import warnings
import argparse
from pathlib import Path
import matplotlib.ticker as ticker

# 控制台编码
if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8')

warnings.filterwarnings('ignore')

# 样式 - 与 step7_1_2 / part5-step3 完全一致
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

FONT_TICK = 28
FONT_LABEL = 34
FONT_ANNOT = 26

# 数据根目录
DEFAULT_DATA = str(
    Path(__file__).parent / 'data' / 'coordination' / 'sup' /
    'coordination_time_series_results_sample_20260310_191856'
)
DEFAULT_OUTPUT = str(Path(__file__).parent / 'results' / 'step7.1.3_sup86_op2')


# ──────────────────────────────────────────────
# 数据扫描
# ──────────────────────────────────────────────

def scan_sup86_data(data_root, sources=None):
    """
    扫描所有 supported Pt8Sn6 的 cluster_op2_time_series.csv

    返回: dict  { temp_K: [df, df, ...] }
    sources: None → 全部; 或 ['4090', 'dp-md-gpu', 'dp-md-more'] 的子集
    """
    data_root = Path(data_root)
    temp_data = {}  # {temp: [df, ...]}

    for root, dirs, files in os.walk(data_root):
        root_path = Path(root)

        op2_file = root_path / 'cluster_op2_time_series.csv'
        if not op2_file.exists():
            continue

        path_str = str(root_path).replace('\\', '/')

        # 判断来源
        if '4090-ustc/more/run3' in path_str or '4090-ustc\\more\\run3' in path_str:
            src_tag = '4090'
        elif 'dp-md' in path_str and 'GPU-Pt8' in path_str:
            src_tag = 'dp-md-gpu'
        elif 'dp-md' in path_str and '/more/' in path_str:
            src_tag = 'dp-md-more'
        elif 'dp-md' in path_str:
            src_tag = 'dp-md-gpu'
        else:
            src_tag = 'other'

        if sources is not None and src_tag not in sources:
            continue

        # 提取温度
        temp = None
        for part in root_path.parts:
            if part.startswith('T') and '.' in part:
                try:
                    temp = int(part[1:].split('.')[0])
                    break
                except ValueError:
                    pass
            elif part.startswith('T') and part[1:].isdigit():
                temp = int(part[1:])
                break

        if temp is None:
            continue

        try:
            df = pd.read_csv(op2_file)
            if temp not in temp_data:
                temp_data[temp] = []
            temp_data[temp].append(df)
        except Exception:
            pass

    return temp_data


def compute_ensemble(data_list, field='op2_all_metal', run_index=None, smooth=True):
    """系综平均或单次运行"""
    values = []
    for df in data_list:
        if df is not None and field in df.columns:
            values.append(df[field].values)

    if not values:
        return None, None, 0

    min_len = min(len(v) for v in values)
    values = [v[:min_len] for v in values]

    if run_index is not None:
        if run_index >= len(values):
            print(f"  ⚠️ run_index={run_index} 超出范围 (共{len(values)}次运行)")
            return None, None, 0
        v = values[run_index]
        if smooth and len(v) > 21:
            v = savgol_filter(v, 21, 3)
        return v, None, 1

    arr = np.array(values)
    mean = np.mean(arr, axis=0)
    std = np.std(arr, axis=0)
    if smooth and len(mean) > 21:
        mean = savgol_filter(mean, 21, 3)
    return mean, std, len(values)


# ──────────────────────────────────────────────
# 模式 1: OP2 均值 vs 温度
# ──────────────────────────────────────────────

def plot_op2_vs_temp(temp_data, field, output_dir, hide_title=False, figsize=(10, 6)):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    temps = sorted(temp_data.keys())
    means = []
    stds = []
    ns = []

    for t in temps:
        m, s, n = compute_ensemble(temp_data[t], field=field)
        if m is not None:
            means.append(np.mean(m))
            stds.append(np.std(m) if s is None else np.mean(s))
            ns.append(n)
        else:
            means.append(np.nan)
            stds.append(np.nan)
            ns.append(0)

    temps = np.array(temps)
    means = np.array(means)
    stds = np.array(stds)

    fig, ax = plt.subplots(figsize=figsize)
    line_color = '#333333'
    fill_color = '#888888'

    ax.fill_between(temps, means - stds, means + stds,
                    alpha=0.25, color=fill_color)
    ax.plot(temps, means, color=line_color, linewidth=3,
            marker='o', markersize=8, markerfacecolor='white',
            markeredgecolor=line_color, markeredgewidth=2)

    ax.set_xlabel('Temperature (K)', fontsize=FONT_LABEL)
    ax.set_ylabel('OP$_2$', fontsize=FONT_LABEL)
    ax.tick_params(axis='both', labelsize=FONT_TICK, width=1.5, length=6)

    if not hide_title:
        ax.text(0.03, 0.95, r'Pt$_8$Sn$_6$/SnO$_2$',
                transform=ax.transAxes, fontsize=FONT_ANNOT,
                va='top', ha='left')

    plt.tight_layout()
    out_file = output_dir / f'sup86_op2_vs_temp_{field}.png'
    plt.savefig(out_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"📊 OP2 vs T 图已保存: {out_file}")

    print(f"\n>>> Supported Pt8Sn6 OP2 ({field}) 温度依赖:")
    for t, m, s, n in zip(temps, means, stds, ns):
        print(f"    {t:5d} K  OP2 = {m:.4f} ± {s:.4f}  (n={n})")

    return out_file


# ──────────────────────────────────────────────
# 模式 2: 单温度时间序列
# ──────────────────────────────────────────────

def plot_ts_single(temp_data, temp, field, output_dir, run_index=None,
                   hide_title=False, hide_temp_label=False, figsize=(10, 4)):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if temp not in temp_data:
        print(f"  ⚠️ 温度 {temp}K 无数据")
        return

    mean, std, n = compute_ensemble(temp_data[temp], field=field, run_index=run_index)
    if mean is None:
        return

    total_time_ps = 175.0
    time = np.linspace(0, total_time_ps, len(mean))

    fig, ax = plt.subplots(figsize=figsize)
    line_color = '#333333'
    fill_color = '#888888'

    if std is not None:
        ax.fill_between(time, mean - std, mean + std,
                        alpha=0.3, color=fill_color)
    ax.plot(time, mean, color=line_color, linewidth=4)

    ax.set_xlabel('Time (ps)', fontsize=FONT_LABEL)
    ax.set_ylabel('OP$_2$', fontsize=FONT_LABEL)
    ax.tick_params(axis='both', labelsize=FONT_TICK, width=1.5, length=6)

    if not hide_title:
        ax.text(0.03, 0.95, r'Pt$_8$Sn$_6$/SnO$_2$',
                transform=ax.transAxes, fontsize=FONT_ANNOT, va='top', ha='left')
    if not hide_temp_label:
        ax.text(0.97, 0.95, f'{temp} K',
                transform=ax.transAxes, fontsize=FONT_ANNOT, va='top', ha='right')

    plt.tight_layout()
    run_str = f'_run{run_index}' if run_index is not None else ''
    out_file = output_dir / f'sup86_op2_ts_{temp}K{run_str}.png'
    plt.savefig(out_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"📊 时间序列图已保存: {out_file}")
    mode_str = f'(Run #{run_index})' if run_index is not None else f'(Ensemble n={n})'
    print(f"    {temp}K {mode_str}: mean OP2 = {np.mean(mean):.4f}")
    return out_file


# ──────────────────────────────────────────────
# 模式 3: 两温度双子图对比 (与 step7_1_2 一致)
# ──────────────────────────────────────────────

def plot_ts_compare(temp_data, temps, field, output_dir, run_index=None,
                    hide_title=False, hide_temp_label=False, figsize=(10, 8),
                    y_ticks=None):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    total_time_ps = 175.0
    data_dict = {}
    for t in temps:
        if t not in temp_data:
            print(f"  ⚠️ 温度 {t}K 无数据，跳过")
            continue
        mean, std, n = compute_ensemble(temp_data[t], field=field, run_index=run_index)
        if mean is not None:
            time = np.linspace(0, total_time_ps, len(mean))
            data_dict[t] = {'mean': mean, 'std': std, 'time': time, 'n': n}

    if not data_dict:
        print("  ⚠️ 无有效数据")
        return

    # 统一 Y 轴范围
    y_vals = []
    for d in data_dict.values():
        y_vals.extend(d['mean'])
        if d['std'] is not None:
            y_vals.extend(d['mean'] - d['std'])
            y_vals.extend(d['mean'] + d['std'])
    y_min, y_max = np.min(y_vals), np.max(y_vals)
    margin = (y_max - y_min) * 0.1
    y_min -= margin
    y_max += margin

    n_panels = len(data_dict)
    fig, axes = plt.subplots(n_panels, 1, figsize=figsize, sharex=True,
                             gridspec_kw={'hspace': 0})
    if n_panels == 1:
        axes = [axes]

    line_color = '#333333'
    fill_color = '#888888'

    for idx, (t, ax) in enumerate(zip(sorted(data_dict.keys()), axes)):
        d = data_dict[t]
        if d['std'] is not None:
            ax.fill_between(d['time'], d['mean'] - d['std'], d['mean'] + d['std'],
                            alpha=0.3, color=fill_color)
        ax.plot(d['time'], d['mean'], color=line_color, linewidth=4)

        ax.set_ylim(y_min, y_max)
        ax.set_ylabel('OP$_2$', fontsize=FONT_LABEL)
        ax.tick_params(axis='y', labelsize=FONT_TICK, width=1.5, length=6)
        # Y 轴刻度：手动指定或自动选取
        if y_ticks is not None:
            ax.set_yticks(y_ticks)
        else:
            ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=4, prune='both'))

        if not hide_title:
            ax.text(0.03, 0.95, r'Pt$_8$Sn$_6$/SnO$_2$',
                    transform=ax.transAxes, fontsize=FONT_ANNOT, va='top', ha='left')
        if not hide_temp_label:
            ax.text(0.97, 0.95, f'{t} K',
                    transform=ax.transAxes, fontsize=FONT_ANNOT, va='top', ha='right')

        if idx < n_panels - 1:
            ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
        else:
            ax.set_xlabel('Time (ps)', fontsize=FONT_LABEL)
            ax.tick_params(axis='x', labelsize=FONT_TICK, width=1.5, length=6)

    plt.tight_layout()
    plt.subplots_adjust(hspace=0)

    temps_str = '_'.join(map(str, sorted(data_dict.keys())))
    run_str = f'_run{run_index}' if run_index is not None else ''
    out_file = output_dir / f'sup86_op2_compare_{temps_str}K{run_str}.png'
    plt.savefig(out_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"📊 温度对比图已保存: {out_file}")

    print(f"\n>>> Supported Pt8Sn6 OP2 ({field}) 温度对比:")
    for t in sorted(data_dict.keys()):
        d = data_dict[t]
        mode_str = f'(Run #{run_index})' if run_index is not None else f'(n={d["n"]})'
        print(f"    {t}K {mode_str}: mean OP2 = {np.mean(d['mean']):.4f}")

    return out_file


# ──────────────────────────────────────────────
# main
# ──────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Supported Pt8Sn6 OP2 分析',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例用法:
  # 1. OP2 均值 vs 温度 扫描图
  python step7_1_3_sup86_op2_analysis.py --mode scan

  # 2. 指定温度时间序列 (系综平均)
  python step7_1_3_sup86_op2_analysis.py --mode ts --temp 900

  # 3. 指定温度时间序列 (单次运行)
  python step7_1_3_sup86_op2_analysis.py --mode ts --temp 900 --run 0

  # 4. 两温度双子图对比 (系综平均)
  python step7_1_3_sup86_op2_analysis.py --mode compare --temps 300 900

  # 5. 两温度双子图对比 (单次运行, 隐藏标注)
  python step7_1_3_sup86_op2_analysis.py --mode compare --temps 300 900 --run 0 --hide-title --hide-temp-label

  # 6. 只使用特定来源
  python step7_1_3_sup86_op2_analysis.py --mode scan --sources dp-md-more
        """
    )
    parser.add_argument('--data', type=str, default=DEFAULT_DATA, help='数据根目录')
    parser.add_argument('--output', type=str, default=DEFAULT_OUTPUT, help='输出目录')
    parser.add_argument('--mode', type=str, default='scan',
                        choices=['scan', 'ts', 'compare'],
                        help='scan: OP2 vs 温度 | ts: 单温度时间序列 | compare: 双温度双子图')
    parser.add_argument('--field', type=str, default='op2_all_metal',
                        choices=['op2_all_metal', 'op2_pt', 'op2_sn', 'op2_global'],
                        help='OP2 字段 (默认: op2_all_metal)')
    parser.add_argument('--temp', type=int, default=900,
                        help='单温度 (ts模式, 默认: 900)')
    parser.add_argument('--temps', type=int, nargs='+', default=[300, 900],
                        help='温度列表 (compare模式, 默认: 300 900)')
    parser.add_argument('--run', type=int, default=None,
                        help='指定运行索引 (0,1,2...), 不指定则系综平均')
    parser.add_argument('--sources', type=str, nargs='+',
                        choices=['4090', 'dp-md-gpu', 'dp-md-more'],
                        default=None,
                        help='数据来源过滤 (默认: 全部)')
    parser.add_argument('--hide-title', action='store_true',
                        help='隐藏图内结构名标注')
    parser.add_argument('--hide-temp-label', action='store_true',
                        help='隐藏图内温度标注')
    parser.add_argument('--y-ticks', type=str, default=None,
                        help='手动指定 Y 轴刻度，逗号分隔，如 0.33,0.35,0.37')
    parser.add_argument('--figsize', type=str, default=None,
                        help='图像尺寸 WxH 英寸 (scan默认10x6, ts默认10x4, compare默认10x8)')

    args = parser.parse_args()

    # 默认 figsize
    default_fs = {'scan': (10, 6), 'ts': (10, 4), 'compare': (10, 8)}
    if args.figsize:
        try:
            fw, fh = [float(x) for x in args.figsize.lower().split('x')]
            figsize = (fw, fh)
        except Exception:
            print("⚠️ --figsize 格式错误，使用默认值")
            figsize = default_fs[args.mode]
    else:
        figsize = default_fs[args.mode]

    # 解析 y_ticks
    y_ticks = None
    if args.y_ticks:
        try:
            y_ticks = [float(v) for v in args.y_ticks.split(',')]
        except Exception:
            print("⚠️ --y-ticks 格式错误，使用自动刻度")

    print(f"\n{'='*60}")
    print(f"Step 7.1.3: Supported Pt8Sn6 OP2 分析")
    print(f"{'='*60}")
    print(f"数据目录: {args.data}")
    print(f"模式: {args.mode}, 字段: {args.field}")
    if args.sources:
        print(f"来源过滤: {args.sources}")

    print("\n>>> 扫描数据...")
    temp_data = scan_sup86_data(args.data, sources=args.sources)
    print(f"    发现温度: {sorted(temp_data.keys())}")
    for t in sorted(temp_data.keys()):
        print(f"    {t}K: {len(temp_data[t])} 个 run")

    if args.mode == 'scan':
        plot_op2_vs_temp(temp_data, args.field, args.output,
                         hide_title=args.hide_title, figsize=figsize)

    elif args.mode == 'ts':
        plot_ts_single(temp_data, args.temp, args.field, args.output,
                       run_index=args.run,
                       hide_title=args.hide_title,
                       hide_temp_label=args.hide_temp_label,
                       figsize=figsize)

    elif args.mode == 'compare':
        plot_ts_compare(temp_data, args.temps, args.field, args.output,
                        run_index=args.run,
                        hide_title=args.hide_title,
                        hide_temp_label=args.hide_temp_label,
                        figsize=figsize,
                        y_ticks=y_ticks)

    print(f"\n{'='*60}")
    print(f"✅ 完成!")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
