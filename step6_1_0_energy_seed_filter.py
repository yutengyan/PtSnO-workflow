#!/usr/bin/env python3
"""
step6_1_0_energy_seed_filter.py
===============================
评估不同随机种子的LAMMPS能量数据质量，与两段分区（P1固相/P2液相）的
相对能量基线比较，筛选合理 run 并给出 --override-energy 建议值。

用法：
    python step6_1_0_energy_seed_filter.py \\
        --energy-table results/step6_1_1_partition_cv/Pt8Sn6_energy_table.csv \\
        --random-csv data/lammps_energy_analysis/.../energy_master_20260311_143941.csv \\
        --p1-range 200-550 --p2-range 600-1100 \\
        --threshold 0.15

也可直接硬编码运行（无参数时使用默认路径）。

输出：
    - 控制台：每个温度的诊断表、多阈值筛选结果、建议命令
    - 可与 step6_1_1_partition_cv_plot.py --override-energy 联动
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import linregress


# ============================================================
# 默认路径（硬编码，无参数时使用）
# ============================================================
DEFAULT_ETABLE = os.path.join(
    os.path.dirname(__file__),
    "results", "step6_1_1_partition_cv", "Pt8Sn6_energy_table.csv"
)
DEFAULT_RANDOM_CSV = os.path.join(
    os.path.dirname(__file__),
    "data", "lammps_energy_analysis", "lammps_energy_analysis",
    "lammps_energy_analysis-random", "energy_master_20260311_143941.csv"
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="评估不同随机种子的能量数据质量（vs 分区基线）",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("--energy-table", type=str, default=DEFAULT_ETABLE,
                        help="step6_1_1 导出的能量表 CSV\n"
                             "须含列: 温度(K), 绝对能量_eV, Ecluster_eV, Ecluster_rel_eV")
    parser.add_argument("--random-csv", type=str, default=DEFAULT_RANDOM_CSV,
                        help="新随机种子的 energy_master CSV\n"
                             "须含列: 温度(K), 平均总能(eV)")
    parser.add_argument("--p1-range", type=str, default="200-550",
                        help="P1(固相)温度范围, 如 200-550")
    parser.add_argument("--p2-range", type=str, default="600-1100",
                        help="P2(液相)温度范围, 如 600-1100")
    parser.add_argument("--threshold", type=float, default=0.15,
                        help="筛选合理 run 的偏差阈值 (eV), 默认 0.15")
    parser.add_argument("--thresholds", type=str, default="0.30,0.25,0.20,0.15,0.12,0.10,0.08,0.05",
                        help="多阈值扫描列表, 逗号分隔")
    return parser.parse_args()


def load_energy_table(path):
    """加载 step6_1_1 导出的能量表"""
    df = pd.read_csv(path)
    # 自动适配中文列名（Windows GBK 编码可能乱码）
    col_map = {}
    for c in df.columns:
        cl = c.strip().lower()
        if "温度" in c or cl == "t":
            col_map[c] = "T"
        elif "绝对" in c or "abs" in cl:
            col_map[c] = "E_abs"
        elif "Ecluster_eV" == c.strip():
            col_map[c] = "Ecluster"
        elif "rel" in cl:
            col_map[c] = "E_rel"
    df = df.rename(columns=col_map)
    return df


def load_random_runs(path):
    """加载新随机种子的 energy_master CSV，返回 {T: [E_abs_list]}"""
    df = pd.read_csv(path)
    # 自动适配列名
    col_T = [c for c in df.columns if "温度" in c or c.lower() == "t"][0]
    col_E = [c for c in df.columns if "平均总能" in c or "avg_energy" in c.lower()][0]

    result = {}
    for T, grp in df.groupby(col_T):
        T_int = int(round(T))
        result[T_int] = grp[col_E].values.tolist()
    return result


def derive_support_baseline(etable):
    """从能量表反推载体能量基线: E_abs = Ecluster + (s_sup*T + i_sup)"""
    if "Ecluster" not in etable.columns:
        return None, None
    E_sup = etable["E_abs"] - etable["Ecluster"]
    s_sup, i_sup, _, _, _ = linregress(etable["T"], E_sup)
    return s_sup, i_sup


def abs2rel(E_abs, T, s_sup, i_sup, E_cluster_ref):
    """系统绝对能量 -> 相对团簇能量（与图 Y 轴一致）"""
    E_cluster = E_abs - (s_sup * T + i_sup)
    return E_cluster - E_cluster_ref


def main():
    args = parse_args()

    # ── 加载能量表 ──────────────────────────────────────────────────
    print("=" * 70)
    print("  step6_1_0 : 随机种子能量数据质量筛选")
    print("=" * 70)

    if not os.path.exists(args.energy_table):
        print(f"  [错误] 能量表不存在: {args.energy_table}")
        sys.exit(1)
    etable = load_energy_table(args.energy_table)
    print(f"  能量表: {args.energy_table}  ({len(etable)} 个温度点)")

    # ── 载体基线 & 相对能量参考 ──────────────────────────────────────
    s_sup, i_sup = derive_support_baseline(etable)
    if s_sup is None:
        print("  [错误] 能量表缺少 Ecluster_eV 列，无法反推载体基线")
        sys.exit(1)

    if "Ecluster" in etable.columns:
        E_cluster_ref = etable["Ecluster"].min()
    else:
        E_cluster_ref = etable["E_rel"].min()
    print(f"  载体基线: slope={s_sup:.6f} eV/K  intercept={i_sup:.4f}")
    print(f"  Ecluster 参考值: {E_cluster_ref:.4f} eV")

    # 验证转换精度
    check = etable.iloc[0]
    converted = abs2rel(check["E_abs"], check["T"], s_sup, i_sup, E_cluster_ref)
    error = abs(converted - check["E_rel"])
    print(f"  转换验证 (T={check['T']:.0f}K): 转换={converted:.4f}  表={check['E_rel']:.4f}  误差={error:.6f}")
    if error > 0.01:
        print(f"  [警告] 转换误差较大 ({error:.4f})，请检查载体基线")

    # ── 两段基线拟合（相对能量坐标）────────────────────────────────
    p1_lo, p1_hi = map(float, args.p1_range.split("-"))
    p2_lo, p2_hi = map(float, args.p2_range.split("-"))

    p1 = etable[(etable["T"] >= p1_lo) & (etable["T"] <= p1_hi)]
    p2 = etable[(etable["T"] >= p2_lo) & (etable["T"] <= p2_hi)]

    s1, i1, r1, _, _ = linregress(p1["T"], p1["E_rel"])
    s2, i2, r2, _, _ = linregress(p2["T"], p2["E_rel"])

    print(f"\n  P1 ({p1_lo:.0f}-{p1_hi:.0f}K): slope={s1:.6f} eV/K  R2={r1**2:.6f}")
    print(f"  P2 ({p2_lo:.0f}-{p2_hi:.0f}K): slope={s2:.6f} eV/K  R2={r2**2:.6f}")

    boundary = (p1_hi + p2_lo) / 2

    def baseline_rel(T):
        return s1 * T + i1 if T <= boundary else s2 * T + i2

    # ── 全温度段偏差总览 ─────────────────────────────────────────
    print(f"\n  === 全温度段 vs 两段基线（相对能量）===")
    print(f"  {'T':>5}  {'E_rel':>8}  {'基线':>8}  {'偏差':>8}  分区")
    print(f"  {'-'*48}")
    for _, row in etable.iterrows():
        T = row["T"]
        E = row["E_rel"]
        pred = baseline_rel(T)
        dev = E - pred
        flag = "  <--!!!" if abs(dev) > 0.15 else ""
        pt = "P1" if T <= boundary else "P2"
        print(f"  {T:5.0f}  {E:8.4f}  {pred:8.4f}  {dev:+8.4f}{flag}  {pt}")

    # ── 加载新随机数据 ───────────────────────────────────────────
    if not os.path.exists(args.random_csv):
        print(f"\n  [错误] 随机种子数据不存在: {args.random_csv}")
        sys.exit(1)
    random_runs = load_random_runs(args.random_csv)
    print(f"\n  随机种子数据: {args.random_csv}")
    print(f"  温度: {sorted(random_runs.keys())}  共 {sum(len(v) for v in random_runs.values())} 条")

    # ── 逐温度分析 ──────────────────────────────────────────────
    thresholds = [float(x) for x in args.thresholds.split(",")]
    main_threshold = args.threshold
    suggestions = []

    for T_key in sorted(random_runs.keys()):
        runs_abs = np.array(random_runs[T_key])
        runs_rel = np.array([abs2rel(e, T_key, s_sup, i_sup, E_cluster_ref) for e in runs_abs])
        pred = baseline_rel(T_key)
        pt = "P1" if T_key <= boundary else "P2"

        # 旧数据
        old_row = etable[etable["T"] == T_key]
        if len(old_row) == 0:
            print(f"\n  [警告] T={T_key}K 不在能量表中，跳过")
            continue
        old_E_rel = float(old_row["E_rel"].iloc[0])
        old_E_abs = float(old_row["E_abs"].iloc[0])
        old_dev = old_E_rel - pred

        print(f"\n  {'='*60}")
        print(f"  T={T_key}K [{pt}]  基线={pred:.4f} eV (相对)")
        print(f"  旧数据: E_rel={old_E_rel:.4f}  偏差={old_dev:+.4f} eV")
        print(f"  新数据: {len(runs_abs)} runs  "
              f"E_rel均值={runs_rel.mean():.4f}  偏差={runs_rel.mean()-pred:+.4f}  "
              f"std={runs_rel.std():.4f}")

        # 逐 run 列表
        print(f"  {'run':>5}  {'abs(eV)':>14}  {'E_rel':>8}  {'vs基线':>8}  状态")
        for i in np.argsort(runs_rel):
            e_a = runs_abs[i]
            e_r = runs_rel[i]
            d = e_r - pred
            st = "OK" if abs(d) <= main_threshold else ("偏高" if d > 0 else "偏低")
            bar = "#" * min(int(abs(d) / 0.05), 22)
            print(f"  run{i:>2}  {e_a:14.4f}  {e_r:8.4f}  {d:+8.4f}  {st:4s}  {bar}")

        # 多阈值扫描
        print(f"\n  --- 阈值扫描 ---")
        print(f"  {'阈值':>8}  {'通过':>4}  {'均值abs':>14}  {'E_rel':>8}  {'vs基线':>8}  {'std':>8}")
        for th in thresholds:
            mask = np.abs(runs_rel - pred) <= th
            n = mask.sum()
            if n > 0:
                gm_abs = runs_abs[mask].mean()
                gm_rel = runs_rel[mask].mean()
                d = gm_rel - pred
                std_abs = runs_abs[mask].std(ddof=1) if n > 1 else 0
                marker = "  <-- *" if abs(th - main_threshold) < 0.001 else ""
                print(f"  {th:8.2f}  {n:>4}  {gm_abs:14.4f}  {gm_rel:8.4f}  {d:+8.4f}  {std_abs:8.4f}{marker}")
            else:
                print(f"  {th:8.2f}  {n:>4}  ---")

        # 主阈值筛选结果
        mask_main = np.abs(runs_rel - pred) <= main_threshold
        n_main = mask_main.sum()
        if n_main > 0:
            gm = runs_abs[mask_main].mean()
            gs = runs_abs[mask_main].std(ddof=1) if n_main > 1 else 0.0
            gm_rel = runs_rel[mask_main].mean()
            gdev = gm_rel - pred
            std_str = f"\u00b1{gs:.4f}" if gs > 1e-6 else ""
            suggestions.append({
                "T": T_key, "mean_abs": gm, "std_abs": gs,
                "mean_rel": gm_rel, "dev": gdev, "n": n_main,
                "old_abs": old_E_abs, "old_dev": old_dev,
            })
        else:
            suggestions.append({
                "T": T_key, "mean_abs": None, "std_abs": None,
                "mean_rel": None, "dev": None, "n": 0,
                "old_abs": old_E_abs, "old_dev": old_dev,
            })

    # ── 汇总建议 ────────────────────────────────────────────────
    print(f"\n{'='*70}")
    print(f"  汇总（主阈值 = {main_threshold} eV）")
    print(f"{'='*70}")

    override_parts = []
    for s in suggestions:
        T = s["T"]
        n = s["n"]
        if n >= 5:
            quality = "OK"
        elif n >= 3:
            quality = "CAUTION"
        elif n >= 1:
            quality = "FEW"
        else:
            quality = "NONE"

        print(f"\n  T={T}K  [{n}/10 合理]  {quality}")
        if s["mean_abs"] is not None:
            gs = s["std_abs"]
            std_part = f"\u00b1{gs:.4f}" if gs > 1e-6 else ""
            print(f"    新数据均值: {s['mean_abs']:.4f} eV  std={gs:.4f}")
            print(f"    E_rel={s['mean_rel']:.4f}  vs基线={s['dev']:+.4f} eV")
            print(f"    旧数据偏差: {s['old_dev']:+.4f} eV")
            override_parts.append(f'"{T}K:{s["mean_abs"]:.4f}{std_part}"')
        else:
            print(f"    无合理run，保留旧值: {s['old_abs']:.4f} eV (偏差={s['old_dev']:+.4f})")
            override_parts.append(f'"{T}K:{s["old_abs"]:.4f}"')

    print(f"\n  建议命令:")
    print(f"    --override-energy {' '.join(override_parts)}")
    print()


if __name__ == "__main__":
    main()
