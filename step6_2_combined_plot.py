#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 6.2: Lindemann + 能量/热容 上下组合图

将 step6_1_3 (Lindemann 指数图) 与 step6_1_1 (能量-热容图) 合并为一张
上下两个子图共享 X 轴的组合图，方便直接对比熔化行为。

  上子图 (ax_lind)  — Lindemann 指数随温度变化
  下子图 (ax_e/ax_cv) — 能量随温度变化 + 右轴 Cv 曲线

用法:
  # 基本用法（复现两个脚本各自的默认输出）
  python step6_2_combined_plot.py

  # 完整参数示例（对应用户原始两条命令）
  python step6_2_combined_plot.py \\
      --partitions-sup86 200-600,650-1100 --x-ticks 200,400,600,800,1000 \\
      --lind-y-ticks 0.1,0.2,0.3 --show-error-bars \\
      --structure Pt8sn6 --partitions 200-550,600-1100 \\
      --exclude "300K:5,9" "400K:0,1,3,7,9" "600K:6" "700K:2,8" "800K:0" "900K:8" "1000K:5" "1100K:2" \\
      --exclude-sort-by energy --cv-y-ticks 0,2,4 --cv-ticks 3,4,5,6 \\
      --peak-method partition --use-lindemann-threshold 0.1 --clustering-method lindemann-threshold \\
      --override-energy "400K:-95384.8107±0.1287" "750K:-95370.0535±0.1804" \\
      --use-lines --hide-override-marker \\
      --figsize 10x12 --height-ratio 1,1 --font-scale 0.75

作者: AI Assistant
日期: 2026-03-11
"""

import os
import sys
import glob
import argparse
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import linregress
from scipy.interpolate import CubicSpline
from pathlib import Path

# ─────────────────────────── 全局绘图样式 ───────────────────────────
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
plt.rcParams['mathtext.fontset'] = 'dejavusans'
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['xtick.major.width'] = 1.5
plt.rcParams['ytick.major.width'] = 1.5

# 字体大小常量（可被 --font-scale 缩放）
FONT_TICK  = 28
FONT_LABEL = 34
FONT_LEGEND = 26

# 分区配色（与两个原始脚本保持一致）
PARTITION_COLORS = {
    1: '#0173B2',   # 蓝色 — 固相
    2: '#DE4343',   # 红色 — 液相
}

# 载体默认热容 (meV/K)
CV_SUPPORT = 38.2151


# ══════════════════════════════════════════════════════════════════════
#  公共工具函数
# ══════════════════════════════════════════════════════════════════════

def parse_ticks(tick_str, dtype=float):
    """解析逗号分隔的刻度字符串，例如 '0,0.1,0.2' → [0.0, 0.1, 0.2]"""
    if not tick_str:
        return None
    try:
        return [dtype(x) for x in tick_str.split(',')]
    except Exception as e:
        print(f"  ⚠ 无法解析刻度 '{tick_str}': {e}")
        return None


def parse_partitions(partition_str):
    """解析分区字符串，例如 '200-550,600-1100' → [(200, 550), (600, 1100)]"""
    if not partition_str:
        return None
    try:
        return [tuple(map(int, p.split('-'))) for p in partition_str.split(',')]
    except Exception as e:
        print(f"  ⚠ 无法解析分区 '{partition_str}': {e}")
        return None


def parse_exclude_points(exclude_args):
    """解析排除点参数，例如 ['300K:5,9', '400K:0'] → {300: [5,9], 400: [0]}"""
    exclude_dict = {}
    if not exclude_args:
        return exclude_dict
    for arg in exclude_args:
        try:
            temp_str, indices_str = arg.split(':')
            temp = int(temp_str.replace('K', '').replace('k', ''))
            indices = [int(i) for i in indices_str.split(',')]
            exclude_dict[temp] = indices
        except Exception as e:
            print(f"  ⚠ 无法解析排除点 '{arg}': {e}")
    return exclude_dict


def parse_override_energy(override_args):
    """
    解析能量覆盖参数，支持 '400K:-95384.715' 或 '400K:-95384.715±0.20' 或 '+/-' 格式
    返回 {temp: (mean_eV, std_eV)}
    """
    override_dict = {}
    if not override_args:
        return override_dict
    for arg in override_args:
        try:
            colon_idx = arg.index(':')
            temp_str = arg[:colon_idx]
            energy_str = arg[colon_idx + 1:]
            temp = int(temp_str.replace('K', '').replace('k', ''))
            std_val = 0.0
            for sep in ['±', '+/-']:
                if sep in energy_str:
                    parts = energy_str.split(sep, 1)
                    energy_str = parts[0]
                    std_val = float(parts[1])
                    break
            override_dict[temp] = (float(energy_str), std_val)
        except Exception as e:
            print(f"  ⚠ 无法解析能量覆盖 '{arg}': {e}")
    return override_dict


def load_support_energy_data():
    """加载载体能量数据，返回 (slope, intercept, R2) 或 None"""
    support_csv = 'data/lammps_energy/sup/energy_master_20251021_151520.csv'
    if not os.path.exists(support_csv):
        return None
    try:
        df = pd.read_csv(support_csv)
        if 'temp' in df.columns and 'avg_energy' in df.columns:
            T = df['temp'].values
            E = df['avg_energy'].values
            slope, intercept, r, _, _ = linregress(T, E)
            return slope, intercept, r**2
    except Exception as e:
        print(f"  ⚠ 读取载体能量失败: {e}")
    return None


def find_clustering_csv(structure, base_dir='results/step6_1_clustering', method='auto'):
    """
    根据结构名查找聚类 CSV 文件。
    method: 'auto'（优先 lindemann-threshold），'lindemann-threshold'，'kmeans'
    """
    # 规范化结构名（文件名大小写不敏感搜索）
    patterns_priority = []
    if method in ('auto', 'lindemann-threshold'):
        patterns_priority.append(f'{structure}_lindemann-threshold_n2_clustered_data.csv')
    if method in ('auto', 'kmeans'):
        patterns_priority.append(f'{structure}_kmeans_n2_clustered_data.csv')

    for pattern in patterns_priority:
        p = Path(base_dir) / pattern
        if p.exists():
            return str(p)
        # 大小写不敏感搜索
        candidates = list(Path(base_dir).glob('*.csv'))
        for c in candidates:
            if c.name.lower() == pattern.lower():
                return str(c)
    return None


def filter_data_by_exclusion(df, exclude_dict, sort_by='delta'):
    """按排除规则过滤 DataFrame"""
    if not exclude_dict:
        return df
    print(f"\n  [排除点过滤] 排序依据: {sort_by}")
    mask = np.ones(len(df), dtype=bool)
    for temp, indices in exclude_dict.items():
        temp_mask = df['temp'] == temp
        temp_indices = np.where(temp_mask)[0]
        if len(temp_indices) == 0:
            print(f"    警告: 温度 {temp}K 没有数据点")
            continue
        temp_df = df[temp_mask].copy()
        temp_df['original_idx'] = temp_indices
        if sort_by == 'energy':
            temp_df_sorted = temp_df.sort_values('avg_energy', ascending=False)
        else:
            temp_df_sorted = temp_df.sort_values('delta', ascending=False)
        for idx in indices:
            if idx < len(temp_df_sorted):
                original_idx = int(temp_df_sorted.iloc[idx]['original_idx'])
                mask[original_idx] = False
                val = (temp_df_sorted.iloc[idx]['avg_energy'] if sort_by == 'energy'
                       else temp_df_sorted.iloc[idx]['delta'])
                label = 'energy' if sort_by == 'energy' else 'delta'
                print(f"    排除: {temp}K 第{idx}个点 ({label}={val:.4f})")
    filtered = df[mask].copy()
    print(f"  过滤结果: {len(df)}条 → {len(filtered)}条 (排除{len(df)-len(filtered)}条)")
    return filtered


# ══════════════════════════════════════════════════════════════════════
#  上子图：Lindemann 指数
# ══════════════════════════════════════════════════════════════════════

def draw_lindemann_ax(ax, df, font_tick, font_label, font_legend,
                      y_ticks=None, x_ticks=None, custom_partitions=None,
                      show_error_bars=True, y_lim=None,
                      hide_x_label=True, partition_labels=None):
    """
    在给定 ax 上绘制 Lindemann 指数图。
    （从 step6_1_3.plot_lindemann_single 提取，去掉 figure 创建）
    """
    if custom_partitions is not None:
        df = df.copy()
        for i, (T_min, T_max) in enumerate(custom_partitions):
            mask = (df['temp'] >= T_min) & (df['temp'] <= T_max)
            df.loc[mask, 'phase_clustered'] = f'partition{i+1}'

    temps = sorted(df['temp'].unique())

    # 散点（不显示误差棒时才绘制原始散点）
    if not show_error_bars:
        for phase in sorted(df['phase_clustered'].unique()):
            df_phase = df[df['phase_clustered'] == phase]
            phase_num = int(phase.replace('partition', '')) if isinstance(phase, str) else int(phase)
            color = PARTITION_COLORS.get(phase_num, '#999999')
            label = (partition_labels[phase_num - 1]
                     if partition_labels and phase_num - 1 < len(partition_labels)
                     else f'Partition {phase_num}')
            ax.scatter(df_phase['temp'], df_phase['delta'],
                       c=color, s=80, alpha=0.7, edgecolors='black', linewidth=0.5,
                       label=label, zorder=3)

    # 误差棒
    if show_error_bars:
        for temp in temps:
            df_temp = df[df['temp'] == temp]
            mean_d = df_temp['delta'].mean()
            std_d  = df_temp['delta'].std()
            phase = df_temp['phase_clustered'].iloc[0]
            phase_num = int(phase.replace('partition', '')) if isinstance(phase, str) else int(phase)
            color = PARTITION_COLORS.get(phase_num, '#999999')
            ax.errorbar(temp, mean_d, yerr=std_d,
                        fmt='none', ecolor=color, elinewidth=2.5,
                        capsize=6, capthick=2.5, alpha=0.5, zorder=1)

    # 分区内均值连线（误差棒模式或 show_average_line 时自动绘制）
    all_phase_data = {}
    for phase in sorted(df['phase_clustered'].unique()):
        df_phase = df[df['phase_clustered'] == phase]
        phase_num = int(phase.replace('partition', '')) if isinstance(phase, str) else int(phase)
        t_means = [(t, df_phase[df_phase['temp'] == t]['delta'].mean())
                   for t in sorted(df_phase['temp'].unique())]
        t_vals = [x[0] for x in t_means]
        m_vals = [x[1] for x in t_means]
        all_phase_data[phase_num] = (t_vals, m_vals)
        color = PARTITION_COLORS.get(phase_num, '#999999')
        label = (partition_labels[phase_num - 1]
                 if partition_labels and phase_num - 1 < len(partition_labels)
                 else f'Partition {phase_num}') if show_error_bars else None
        ax.plot(t_vals, m_vals,
                color=color, linestyle='-', linewidth=2.5, alpha=0.8, zorder=3,
                marker='o', markersize=14.14,
                markerfacecolor=color, markeredgecolor=color, markeredgewidth=0,
                label=label)

    # 跨分区渐变连接线
    if 1 in all_phase_data and 2 in all_phase_data:
        t1, m1 = all_phase_data[1]
        t2, m2 = all_phase_data[2]
        x_s, x_e = t1[-1], t2[0]
        y_s, y_e = m1[-1], m2[0]
        n = 100
        xi = np.linspace(x_s, x_e, n)
        yi = np.linspace(y_s, y_e, n)
        pts = np.array([xi, yi]).T.reshape(-1, 1, 2)
        segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
        cmap = LinearSegmentedColormap.from_list('br', [PARTITION_COLORS[1], PARTITION_COLORS[2]])
        lc = LineCollection(segs, cmap=cmap, linewidth=2.5, alpha=0.8, zorder=3)
        lc.set_array(np.linspace(0, 1, len(segs)))
        ax.add_collection(lc)

    # 轴标签
    if hide_x_label:
        ax.set_xlabel('', fontsize=font_label)
        ax.tick_params(axis='x', which='both', bottom=False, top=False,
                       labelbottom=False, labelsize=font_tick)
    else:
        ax.set_xlabel('Temperature (K)', fontsize=font_label)
        ax.tick_params(axis='x', labelsize=font_tick, width=1.5, length=6)
    ax.set_ylabel('Lindemann Index', fontsize=font_label)
    ax.tick_params(axis='y', labelsize=font_tick, width=1.5, length=6)

    # X 轴范围和刻度
    ax.set_xlim(150, 1150)
    if x_ticks is not None:
        ax.set_xticks(x_ticks)
    else:
        ax.set_xticks(temps)

    # Y 轴范围和刻度
    if y_lim is not None:
        ax.set_ylim(y_lim)
    if y_ticks is not None:
        ax.set_yticks(y_ticks)

    # 图例
    ax.legend(loc='lower right', fontsize=font_legend, frameon=False)


# ══════════════════════════════════════════════════════════════════════
#  下子图：能量 + 热容（从 step6_1_1.plot_partition_cv 提取核心逻辑）
# ══════════════════════════════════════════════════════════════════════

def _determine_partitions_by_lindemann(df, threshold=0.1):
    """按 Lindemann 阈值自动分区，返回 (custom_partitions, temp_to_partition)"""
    temp_delta = df.groupby('temp')['delta'].mean().reset_index()
    temp_delta = temp_delta.sort_values('temp')
    temp_delta['partition'] = temp_delta['delta'].apply(lambda x: 1 if x < threshold else 2)

    p1_temps = temp_delta[temp_delta['partition'] == 1]['temp'].values
    p2_temps = temp_delta[temp_delta['partition'] == 2]['temp'].values

    if len(p1_temps) == 0 or len(p2_temps) == 0:
        return None, None

    custom_partitions = [(p1_temps.min(), p1_temps.max()),
                         (p2_temps.min(), p2_temps.max())]
    temp_to_partition = {t: (1 if t in p1_temps else 2) for t in temp_delta['temp']}
    return custom_partitions, temp_to_partition


def _calc_transition_cv(df, temp_to_partition, phase_fits, temps_unique, E_rel,
                        T_left, T_right, ph_left, ph_right,
                        peak_method, is_air, slope_sup, int_sup):
    """计算过渡区热容（对应 step6_1_1.calculate_transition_cv）"""
    result = {}

    idx_l = np.where(temps_unique == T_left)[0]
    idx_r = np.where(temps_unique == T_right)[0]
    if len(idx_l) and len(idx_r):
        result['Cv_data'] = (E_rel[idx_r[0]] - E_rel[idx_l[0]]) / (T_right - T_left) * 1000
    else:
        result['Cv_data'] = (phase_fits[ph_left]['Cv'] + phase_fits[ph_right]['Cv']) / 2

    df_l = df[(df['temp'] == T_left)  & (df['phase_clustered'] == ph_left)]
    df_r = df[(df['temp'] == T_right) & (df['phase_clustered'] == ph_right)]
    if len(df_l) and len(df_r):
        if is_air:
            El = df_l['avg_energy'].mean()
            Er = df_r['avg_energy'].mean()
        else:
            El = df_l['avg_energy'].mean() - (slope_sup * T_left  + int_sup)
            Er = df_r['avg_energy'].mean() - (slope_sup * T_right + int_sup)
        result['Cv_partition'] = (Er - El) / (T_right - T_left) * 1000
    else:
        result['Cv_partition'] = result['Cv_data']

    fl = phase_fits[ph_left]
    fr = phase_fits[ph_right]
    result['Cv_fit'] = ((fr['slope'] * T_right + fr['intercept']) -
                        (fl['slope'] * T_left  + fl['intercept'])) / (T_right - T_left) * 1000
    return result


def draw_energy_cv_ax(ax_e, ax_cv, df, font_tick, font_label,
                      custom_partitions=None,
                      peak_method='partition',
                      use_lines=False,
                      override_energy=None,
                      hide_override_marker=False,
                      y_ticks=None, cv_ticks=None,
                      x_ticks=None,
                      lindemann_threshold=None,
                      clustering_method='auto'):
    """
    在给定的 ax_e（能量左轴）和 ax_cv（热容右轴，twin）上绘制能量-热容图。
    （从 step6_1_1.plot_partition_cv 提取，去掉 figure/保存逻辑）

    返回 phase_fits dict（供外部使用）。
    """
    is_air = df.attrs.get('is_air', False)

    # ── 载体能量基线 ──
    if is_air:
        slope_sup, int_sup = 0.0, 0.0
    else:
        sup = load_support_energy_data()
        if sup is not None:
            slope_sup, int_sup, _ = sup
        else:
            slope_sup = CV_SUPPORT / 1000
            T_min = df['temp'].min()
            int_sup = df[df['temp'] == T_min]['avg_energy'].mean() * 0.9 - slope_sup * T_min
            print("  [警告] 使用默认Cv_support估算载体能量")

    # ── 按温度分组计算团簇能量均值/std ──
    temps_unique, E_mean, E_std = [], [], []
    for temp, grp in df.groupby('temp'):
        E_cl = grp['avg_energy'].values - (0 if is_air else slope_sup * temp + int_sup)
        temps_unique.append(temp)
        E_mean.append(E_cl.mean())
        E_std.append(E_cl.std() if len(E_cl) > 1 else 0.0)

    temps_unique = np.array(temps_unique)
    E_mean  = np.array(E_mean)
    E_std   = np.array(E_std)

    # ── 能量覆盖 ──
    override_mask = np.zeros(len(temps_unique), dtype=bool)
    if override_energy:
        for i, temp in enumerate(temps_unique):
            t_int = int(round(temp))
            if t_int in override_energy:
                raw_e, raw_std = override_energy[t_int]
                cluster_e   = raw_e if is_air else raw_e - (slope_sup * temp + int_sup)
                E_mean[i]   = cluster_e
                E_std[i]    = raw_std
                override_mask[i] = True

    E_ref  = E_mean.min()
    E_rel  = E_mean - E_ref

    # ── Lindemann 阈值自动分区（可选） ──
    if lindemann_threshold is not None and custom_partitions is None:
        cp, _tp = _determine_partitions_by_lindemann(df, lindemann_threshold)
        if cp is not None:
            custom_partitions = cp

    # ── 重新分类：根据 custom_partitions 更新 phase_clustered ──
    # 注意：这里需要区分两种情况：
    #  1. 按温度范围（custom_partitions）→ 用于确定能量拟合区间
    #  2. 按 Lindemann 阈值逐点分类    → 用于 Cv_partition 子集过滤
    # 对应 step6_1_1 中 "根据δ阈值重新分配数据点" 的逻辑
    if custom_partitions is not None:
        df = df.copy()
        if lindemann_threshold is not None:
            # 按每点 delta 值逐点分类（还原 step6_1_1 的行为）
            print(f"\n  [重新分类] 根据δ阈值={lindemann_threshold}逐点重新分配数据点...")
            df.loc[df['delta'] <  lindemann_threshold, 'phase_clustered'] = 'partition1'
            df.loc[df['delta'] >= lindemann_threshold, 'phase_clustered'] = 'partition2'
        else:
            # 按温度范围分类
            print(f"\n  [重新分类] 根据自定义分区按温度范围重新分配数据点...")
            for i, (T_min, T_max) in enumerate(custom_partitions):
                mask = (df['temp'] >= T_min) & (df['temp'] <= T_max)
                df.loc[mask, 'phase_clustered'] = f'partition{i+1}'

    # ── 每个温度的分区归属 ──
    temp_to_partition = {}
    if custom_partitions is not None:
        for i, (T_min, T_max) in enumerate(custom_partitions):
            pname = f'partition{i+1}'
            for temp in temps_unique:
                if T_min <= temp <= T_max:
                    temp_to_partition[temp] = pname
        unassigned = [t for t in temps_unique if t not in temp_to_partition]
        for temp in unassigned:
            min_dist, nearest = float('inf'), None
            for i, (T_min, T_max) in enumerate(custom_partitions):
                d = min(abs(temp - T_min), abs(temp - T_max))
                if d < min_dist:
                    min_dist, nearest = d, f'partition{i+1}'
            temp_to_partition[temp] = nearest
    else:
        for temp in temps_unique:
            grp = df[df['temp'] == temp]
            temp_to_partition[temp] = grp['phase_clustered'].value_counts().idxmax()

    # ── 分区拟合 ──
    phases_sorted = sorted(set(temp_to_partition.values()),
                           key=lambda x: int(x.replace('partition', '')))
    phase_fits = {}
    for phase in phases_sorted:
        ph_temps = sorted(t for t, p in temp_to_partition.items() if p == phase)
        if len(ph_temps) < 2:
            continue
        mask_ph = np.isin(temps_unique, ph_temps)
        T_ph  = temps_unique[mask_ph]
        E_ph  = E_rel[mask_ph]
        std_ph = E_std[mask_ph]
        sl, ic, rv, _, se = linregress(T_ph, E_ph)
        phase_fits[phase] = {
            'slope': sl, 'intercept': ic,
            'R2': rv**2, 'Cv': sl * 1000, 'Cv_err': se * 1000,
            'n_temps': len(T_ph),
            'T_range': (T_ph.min(), T_ph.max()),
            'T_data': T_ph, 'E_data': E_ph, 'E_std': std_ph
        }
        print(f"  {phase}: Cv={sl*1000:.4f}+/-{se*1000:.4f} meV/K, "
              f"R2={rv**2:.4f}, n={len(T_ph)}, T={T_ph.min():.0f}-{T_ph.max():.0f}K")

    # ── 绘制能量数据点 ──
    ax_e.errorbar(temps_unique, E_rel, yerr=E_std,
                  fmt='o', markersize=10, color='black',
                  ecolor='gray', elinewidth=2, capsize=4, capthick=2,
                  zorder=5, label='Data')
    if override_mask.any() and not hide_override_marker:
        ax_e.scatter(temps_unique[override_mask], E_rel[override_mask],
                     marker='*', s=240, color='red', zorder=7,
                     label='Overridden energy')

    # ── 绘制拟合线 / 连接线 ──
    for phase in phases_sorted:
        fit = phase_fits[phase]
        if use_lines:
            ph_mask = ((temps_unique >= fit['T_range'][0]) &
                       (temps_unique <= fit['T_range'][1]))
            ax_e.plot(temps_unique[ph_mask], E_rel[ph_mask],
                      '-', color='black', linewidth=2.5, zorder=4)
        else:
            T_fit = np.linspace(fit['T_range'][0], fit['T_range'][1], 50)
            ax_e.plot(T_fit, fit['slope'] * T_fit + fit['intercept'],
                      '-', color='black', linewidth=2.5, zorder=4)

    # 跨分区连接线
    for i in range(len(phases_sorted) - 1):
        fl = phase_fits[phases_sorted[i]]
        fr = phase_fits[phases_sorted[i+1]]
        T_le = fl['T_range'][1];  idx_l = np.where(temps_unique == T_le)[0]
        T_rs = fr['T_range'][0];  idx_r = np.where(temps_unique == T_rs)[0]
        E_le = E_rel[idx_l[0]] if len(idx_l) else fl['slope']*T_le + fl['intercept']
        E_rs = E_rel[idx_r[0]] if len(idx_r) else fr['slope']*T_rs + fr['intercept']
        ax_e.plot([T_le, T_rs], [E_le, E_rs], '-', color='black', linewidth=2.5, zorder=4)

    # 轴标签和刻度
    ax_e.set_xlabel('Temperature (K)', fontsize=font_label)
    ax_e.set_ylabel('Total Energy (eV)', fontsize=font_label)
    ax_e.tick_params(axis='both', labelsize=font_tick)

    if x_ticks is not None:
        ax_e.set_xticks(x_ticks)
        ax_e.set_xlim(150, 1150)

    if y_ticks is not None:
        ax_e.set_yticks(y_ticks)

    # ── 热容曲线 ──
    Cv_list = [phase_fits[p]['Cv'] for p in phases_sorted]
    boundaries = []

    if len(phases_sorted) >= 2:
        for i in range(len(phases_sorted) - 1):
            ph_cur_temps = [t for t, p in temp_to_partition.items() if p == phases_sorted[i]]
            ph_nxt_temps = [t for t, p in temp_to_partition.items() if p == phases_sorted[i+1]]
            T_cl = max(ph_cur_temps)
            T_nf = min(ph_nxt_temps)
            T_bd = (T_cl + T_nf) / 2

            cv_res = _calc_transition_cv(df, temp_to_partition, phase_fits,
                                         temps_unique, E_rel,
                                         T_cl, T_nf,
                                         phases_sorted[i], phases_sorted[i+1],
                                         peak_method, is_air, slope_sup, int_sup)
            Cv_trans = cv_res[f'Cv_{peak_method}'] if f'Cv_{peak_method}' in cv_res else cv_res['Cv_data']
            has_peak = Cv_trans > max(Cv_list[i], Cv_list[i+1])
            boundaries.append({'T_boundary': T_bd, 'T_left': T_cl, 'T_right': T_nf,
                                'Cv_trans': Cv_trans, 'has_peak': has_peak,
                                'idx_left': i, 'idx_right': i+1})
            print(f"  ★ 边界{i+1}: Cv_peak={Cv_trans:.2f} meV/K ({'有峰' if has_peak else '无峰'})")

        # 绘制热容曲线
        T_plot, Cv_plot = [temps_unique.min()], [Cv_list[0]]
        for i, bd in enumerate(boundaries):
            T_plot.append(bd['T_left'])
            Cv_plot.append(Cv_list[i])
            if bd['has_peak']:
                cs = CubicSpline([bd['T_left'], bd['T_boundary'], bd['T_right']],
                                 [Cv_list[i], bd['Cv_trans'], Cv_list[i+1]],
                                 bc_type='clamped')
                T_sm = np.linspace(bd['T_left'], bd['T_right'], 50)
                T_plot.extend(T_sm[1:])
                Cv_plot.extend(cs(T_sm)[1:])
            else:
                T_plot.append(bd['T_right'])
                Cv_plot.append(Cv_list[i+1])
        T_plot.append(temps_unique.max())
        Cv_plot.append(Cv_list[-1])
        ax_cv.plot(T_plot, Cv_plot, 'r-', linewidth=2.5, zorder=3)
    else:
        Cv_single = Cv_list[0]
        ax_cv.plot([temps_unique.min(), temps_unique.max()],
                   [Cv_single, Cv_single], 'r-', linewidth=2.5, zorder=3)

    ax_cv.set_ylabel(r'$C_v$ (meV/K)', fontsize=font_label, color='red')
    ax_cv.tick_params(axis='y', labelcolor='red', labelsize=font_tick, color='red')
    ax_cv.spines['right'].set_color('red')

    # Cv 轴范围
    cv_vals = Cv_list.copy()
    for bd in boundaries:
        if bd.get('has_peak', False):
            cv_vals.append(bd['Cv_trans'])
    if cv_vals:
        ax_cv.set_ylim(min(cv_vals) * 0.85, max(cv_vals) * 1.1)

    if cv_ticks is not None:
        ax_cv.set_yticks(cv_ticks)

    return phase_fits


# ══════════════════════════════════════════════════════════════════════
#  主流程
# ══════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description='Step 6.2: Lindemann + 能量/热容 上下组合图',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # 快速上手（复现 step6_1_3 + step6_1_1 的典型输出）
  python step6_2_combined_plot.py \\
      --partitions-sup86 200-600,650-1100 --x-ticks 200,400,600,800,1000 \\
      --lind-y-ticks 0.1,0.2,0.3 --show-error-bars \\
      --structure Pt8sn6 --partitions 200-550,600-1100 \\
      --exclude "300K:5,9" "400K:0,1,3,7,9" "600K:6" "700K:2,8" \\
             "800K:0" "900K:8" "1000K:5" "1100K:2" \\
      --exclude-sort-by energy --cv-y-ticks 0,2,4 --cv-ticks 3,4,5,6 \\
      --peak-method partition --use-lindemann-threshold 0.1 \\
      --clustering-method lindemann-threshold \\
      --override-energy "400K:-95384.8107±0.1287" "750K:-95370.0535±0.1804" \\
      --use-lines --hide-override-marker \\
      --figsize 10x12 --height-ratio 1,1 --font-scale 0.75
""")

    # ── 通用参数 ──
    parser.add_argument('--figsize', type=str, default='10x12',
                        help='整体图片尺寸，格式 WxH (默认: 10x12)')
    parser.add_argument('--dpi', type=int, default=300,
                        help='分辨率 (默认: 300)')
    parser.add_argument('--height-ratio', type=str, default='1,1',
                        metavar='R1,R2',
                        help='上下子图高度比，例如 1,1 (等高) 或 1,1.5 (默认: 1,1)')
    parser.add_argument('--font-scale', type=float, default=1.0,
                        metavar='SCALE',
                        help='字体整体缩放比例 (默认: 1.0)，例如 0.75 缩小到 75%%')
    parser.add_argument('--x-ticks', type=str, metavar='X1,X2,...',
                        help='X轴刻度（两个子图共享），例如 200,400,600,800,1000')
    parser.add_argument('--output-dir', type=str, default='results/step6_2_combined',
                        help='输出目录 (默认: results/step6_2_combined)')
    parser.add_argument('--format', type=str, default='png', choices=['png', 'pdf', 'svg'],
                        help='输出格式 (默认: png)')

    # ── Lindemann 子图参数（对应 step6_1_3） ──
    lind_grp = parser.add_argument_group('Lindemann 子图参数（上子图）')
    lind_grp.add_argument('--lind-y-ticks', type=str, metavar='Y1,Y2,...',
                           help='Lindemann Y轴刻度，例如 0.1,0.2,0.3')
    lind_grp.add_argument('--show-error-bars', action='store_true',
                           help='显示 Lindemann 误差棒（标准差）')
    lind_grp.add_argument('--partitions-sup86', type=str, metavar='T1-T2,T3-T4',
                           help='sup86 的 Lindemann 分区，例如 200-600,650-1100')
    lind_grp.add_argument('--exclude-sup86', nargs='+', metavar='TEMP:INDICES',
                           help='sup86 Lindemann 排除点，格式: "300K:0,1"')
    lind_grp.add_argument('--partition-labels', type=str, metavar='NAME1,NAME2',
                           help='Lindemann 图例名，逗号分隔，例如 "Solid,Liquid"')
    lind_grp.add_argument('--lind-align-y', action='store_true',
                           help='对齐 Lindemann Y轴到数据范围（自动设置 y_lim）')

    # ── 能量/热容子图参数（对应 step6_1_1） ──
    cv_grp = parser.add_argument_group('能量/热容子图参数（下子图）')
    cv_grp.add_argument('--structure', type=str, default='Pt8sn6',
                         help='结构名称 (默认: Pt8sn6)')
    cv_grp.add_argument('--clustering-method', type=str, default='auto',
                         choices=['auto', 'lindemann-threshold', 'kmeans'],
                         help='聚类文件优先级 (默认: auto)')
    cv_grp.add_argument('--partitions', type=str, metavar='T1-T2,T3-T4',
                         help='能量/Cv 手动分区，例如 200-550,600-1100')
    cv_grp.add_argument('--use-lindemann-threshold', type=float, default=None,
                         metavar='THRESHOLD',
                         help='用 Lindemann 阈值自动分区（覆盖 --partitions）')
    cv_grp.add_argument('--exclude', nargs='+', metavar='TEMP:INDICES',
                         help='能量/Cv 排除点，格式: "300K:5,9"')
    cv_grp.add_argument('--exclude-sort-by', type=str, default='delta',
                         choices=['delta', 'energy'],
                         help='排除点排序依据 (默认: delta)')
    cv_grp.add_argument('--peak-method', type=str, default='partition',
                         choices=['data', 'partition', 'fit'],
                         help='热容峰计算方法 (默认: partition)')
    cv_grp.add_argument('--use-lines', action='store_true',
                         help='用直线连接各温度均值点（不绘制拟合回归线）')
    cv_grp.add_argument('--override-energy', nargs='+', metavar='TEMPK:ENERGY_EV',
                         help='覆盖指定温度的平均能量，支持 ±std，例如 "400K:-95384.8107±0.1287"')
    cv_grp.add_argument('--hide-override-marker', action='store_true',
                         help='隐藏 override-energy 红色★标记')
    cv_grp.add_argument('--cv-y-ticks', type=str, metavar='Y1,Y2,...',
                         help='能量 Y轴刻度，例如 0,2,4')
    cv_grp.add_argument('--cv-ticks', type=str, metavar='C1,C2,...',
                         help='Cv 右轴刻度，例如 3,4,5,6')

    args = parser.parse_args()

    print("=" * 70)
    print("Step 6.2: Lindemann + 能量/热容 上下组合图")
    print("=" * 70)

    # ── 解析通用参数 ──
    try:
        fw, fh = map(float, args.figsize.lower().split('x'))
        figsize = (fw, fh)
    except Exception:
        print(f"⚠ 无法解析 figsize '{args.figsize}'，使用默认 10x12")
        figsize = (10, 12)

    try:
        hr = [float(r) for r in args.height_ratio.split(',')]
        if len(hr) != 2:
            raise ValueError
    except Exception:
        print(f"⚠ 无法解析 height-ratio '{args.height_ratio}'，使用默认 1,1")
        hr = [1, 1]

    font_tick   = int(round(FONT_TICK  * args.font_scale))
    font_label  = int(round(FONT_LABEL * args.font_scale))
    font_legend = int(round(FONT_LEGEND * args.font_scale))

    x_ticks    = parse_ticks(args.x_ticks)
    lind_y_ticks = parse_ticks(args.lind_y_ticks)
    cv_y_ticks = parse_ticks(args.cv_y_ticks)
    cv_ticks   = parse_ticks(args.cv_ticks)
    custom_partitions_sup86 = parse_partitions(args.partitions_sup86)
    custom_partitions_cv    = parse_partitions(args.partitions)
    exclude_sup86  = parse_exclude_points(args.exclude_sup86)
    exclude_cv     = parse_exclude_points(args.exclude)
    override_energy = parse_override_energy(args.override_energy)
    partition_labels = ([s.strip() for s in args.partition_labels.split(',')]
                        if args.partition_labels else None)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # ── 加载 Lindemann 数据（sup86） ──
    print(f"\n[Lindemann 子图] 加载 sup86 数据...")
    lind_csv = find_clustering_csv(args.structure,
                                   method=args.clustering_method)
    if lind_csv is None:
        print(f"✗ 找不到 {args.structure} 的聚类 CSV，请检查 results/step6_1_clustering/")
        return 1

    df_lind = pd.read_csv(lind_csv)
    print(f"  ✓ 已加载: {lind_csv}  ({len(df_lind)} 条)")
    if exclude_sup86:
        df_lind = filter_data_by_exclusion(df_lind, exclude_sup86, sort_by='delta')

    # ── 加载 能量/Cv 数据 ──
    print(f"\n[能量/Cv 子图] 加载 {args.structure} 数据...")
    cv_csv = find_clustering_csv(args.structure,
                                 method=args.clustering_method)
    if cv_csv is None:
        print(f"✗ 找不到 {args.structure} 的聚类 CSV")
        return 1

    df_cv_raw = pd.read_csv(cv_csv)
    print(f"  ✓ 已加载: {cv_csv}  ({len(df_cv_raw)} 条)")
    is_air = args.structure.lower().startswith('air')
    df_cv_raw.attrs['is_air'] = is_air
    if exclude_cv:
        df_cv_raw = filter_data_by_exclusion(df_cv_raw, exclude_cv,
                                              sort_by=args.exclude_sort_by)
    df_cv_raw.attrs['is_air'] = is_air  # attrs 在 filter 后可能丢失

    # ── 构建图形 ──
    print(f"\n>>> 构建组合图  figsize={figsize}  height_ratio={hr}  font_scale={args.font_scale}")
    fig = plt.figure(figsize=figsize, dpi=args.dpi)
    gs  = GridSpec(2, 1, figure=fig,
                   height_ratios=hr,
                   hspace=0.0)          # 子图间无间距（共享 X 轴）

    ax_lind = fig.add_subplot(gs[0])
    ax_e    = fig.add_subplot(gs[1], sharex=ax_lind)  # 共享 X 轴
    ax_cv   = ax_e.twinx()

    # ── 绘制上子图：Lindemann ──
    print("\n── 上子图: Lindemann 指数 ──")
    y_lim_lind = None
    if args.lind_align_y:
        d_min = df_lind['delta'].min()
        d_max = df_lind['delta'].max()
        margin = (d_max - d_min) * 0.05
        y_lim_lind = (d_min - margin, d_max + margin)

    draw_lindemann_ax(ax_lind, df_lind,
                      font_tick, font_label, font_legend,
                      y_ticks=lind_y_ticks,
                      x_ticks=x_ticks,
                      custom_partitions=custom_partitions_sup86,
                      show_error_bars=args.show_error_bars,
                      y_lim=y_lim_lind,
                      hide_x_label=True,          # 上图隐藏 X 轴，由下图显示
                      partition_labels=partition_labels)

    # 上图的 X 轴刻度线也隐藏（sharex 已共享轴，但刻度标签由下图显示）
    plt.setp(ax_lind.get_xticklabels(), visible=False)
    ax_lind.tick_params(axis='x', which='both', bottom=True, labelbottom=False)

    # ── 绘制下子图：能量 + Cv ──
    print("\n── 下子图: 能量 + Cv ──")
    draw_energy_cv_ax(ax_e, ax_cv, df_cv_raw,
                      font_tick, font_label,
                      custom_partitions=custom_partitions_cv,
                      peak_method=args.peak_method,
                      use_lines=args.use_lines,
                      override_energy=override_energy,
                      hide_override_marker=args.hide_override_marker,
                      y_ticks=cv_y_ticks,
                      cv_ticks=cv_ticks,
                      x_ticks=x_ticks,
                      lindemann_threshold=args.use_lindemann_threshold,
                      clustering_method=args.clustering_method)

    # ── 保存 ──
    plt.tight_layout()
    out_name = f'{args.structure}_combined_lind_cv.{args.format}'
    out_path = output_dir / out_name
    fig.savefig(out_path, dpi=args.dpi, bbox_inches='tight', facecolor='white')
    plt.close(fig)

    print(f"\n✅ 组合图已保存: {out_path}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
