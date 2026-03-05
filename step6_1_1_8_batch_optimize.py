#!/usr/bin/env python3
"""
批量优化排除点分析工具
对多个结构进行自动化的排除点优化分析

分区策略:
- T1(固相边界): Lindemann 阈值 (T1_lindemann)
- T2(液相边界): O 活化温度 (T2_B 或 T2_Bprime)
- 数据来源: T1_T2_summary.csv  (自动读取，无需手动维护硬编码字典)
"""

import os
import sys
import re
import subprocess
import pandas as pd
from pathlib import Path


# ============================================================================
# 从 T1_T2_summary.csv 自动加载分区配置
# ============================================================================

def load_partitions_from_csv(csv_file='T1_T2_summary.csv', T2_col='T2_B'):
    """
    从 T1_T2_summary.csv 读取每个结构的分区边界，自动生成三分区配置。

    分区规则:
      partition1 (固态):   200K → round(T1, -2)
      partition2 (熔化中): round(T1, -2) → T2
      partition3 (液态):   T2 → 1800K

    注意:
      - T1_lindemann 会被对齐到最近的 100K 整数（向下取整）用作边界
      - T2 直接使用 CSV 中的整数值
      - 若 T1_boundary >= T2，则退化为 2 分区（无中间相）

    Parameters
    ----------
    csv_file : str
        T1_T2_summary.csv 的路径
    T2_col : str
        使用哪列作为 T2：'T2_B' 或 'T2_Bprime'

    Returns
    -------
    dict: {structure_name: {'partitions': [(T_lo, T_hi), ...], 'T1': float, 'T2': float}}
    """
    if not Path(csv_file).exists():
        raise FileNotFoundError(f"找不到分区配置文件: {csv_file}")

    df = pd.read_csv(csv_file)
    # 列名映射（兼容带括号的列名）
    col_map = {}
    for col in df.columns:
        if 'T1_lindemann' in col:
            col_map['T1'] = col
        if 'T2_B(' in col and 'prime' not in col:
            col_map['T2_B'] = col
        if 'T2_Bprime' in col or "T3_onset" in col:
            col_map['T2_Bprime'] = col

    T1_col_name = col_map.get('T1', 'T1_lindemann(K)')
    T2_col_name = col_map.get(T2_col, col_map.get('T2_B'))

    result = {}
    for _, row in df.iterrows():
        structure = str(row['case']).strip()
        T1_raw = float(row[T1_col_name])
        T2_val = float(row[T2_col_name])

        # T1 向下对齐到 100K 整数（例如 1492 → 1400）
        T1_boundary = int(T1_raw // 100) * 100

        # 避免 T1_boundary >= T2 的情况（退化为 2 分区）
        if T1_boundary >= T2_val:
            partitions = [
                (200, T2_val),
                (T2_val, 1800),
            ]
        else:
            partitions = [
                (200, T1_boundary),
                (T1_boundary, T2_val),
                (T2_val, 1800),
            ]

        result[structure] = {
            'partitions': partitions,
            'T1_raw': T1_raw,
            'T1_boundary': T1_boundary,
            'T2': T2_val,
            'n_partitions': len(partitions),
            # 数据文件可能使用不带 "g-xxx-" 前缀的短名，自动提取备用
            'data_name': re.sub(r'^g-\d+-', '', structure),
        }

    return result


# 默认 CSV 文件和 T2 列
DEFAULT_CSV = 'T1_T2_summary.csv'
DEFAULT_T2_COL = 'T2_Bprime'


def get_partition_str(partitions):
    """将 [(200,1400),(1400,1700),(1700,1800)] 转成 '200-1400,1400-1700,1700-1800'"""
    return ','.join(f"{int(lo)}-{int(hi)}" for lo, hi in partitions)


def extract_plot_cmd_from_report(structure, data_name=None):
    """
    从 *_EXCLUDE_RECOMMENDATIONS.md 中提取 --partitions 和 --exclude 参数。
    先用 structure 找，找不到再用 data_name（短名）。

    返回 (report_name, partitions_str, exclude_args_list) 或 None。
    """
    for name in dict.fromkeys([structure, data_name]):
        if not name:
            continue
        report_file = f"{name}_EXCLUDE_RECOMMENDATIONS.md"
        if os.path.exists(report_file):
            break
    else:
        return None

    try:
        content = Path(report_file).read_text(encoding='utf-8')
    except Exception:
        return None

    # 在 powershell 代码块里找完整命令（跨行，用反引号续行）
    # 找 ```powershell ... ``` 块
    ps_block = re.search(r'```powershell\s*(.*?)```', content, re.DOTALL)
    if not ps_block:
        return None

    cmd_text = ps_block.group(1)
    # 去掉续行符，拼成一行
    cmd_text = cmd_text.replace('`\n', ' ').replace('\\\n', ' ')
    cmd_text = ' '.join(cmd_text.split())  # 压缩多余空白

    # 提取 --partitions
    m_part = re.search(r'--partitions\s+(\S+)', cmd_text)
    partitions_str = m_part.group(1) if m_part else None

    # 提取所有 --exclude 参数（可能有多个 "xxxK:0,1,2" 片段）
    # 找 --exclude 后面的所有带引号的参数
    m_excl = re.search(r'--exclude\s+((?:"[^"]*"\s*)+)', cmd_text)
    if m_excl:
        exclude_args = re.findall(r'"([^"]*)"', m_excl.group(1))
    else:
        exclude_args = []

    return name, partitions_str, exclude_args


def run_plot_for_structure(structure, figsize="10x8", extra_args=None, data_name=None):
    """
    读取报告里的筛选参数，调用 step6_1_1_partition_cv_plot.py 绘图。

    Parameters
    ----------
    structure  : str  CSV 中的全名（用于日志显示）
    data_name  : str  数据文件/报告的短名（去掉 g-xxx- 前缀）
    figsize    : str  传给 --figsize，默认 '10x8'
    extra_args : list  附加参数
    """
    parsed = extract_plot_cmd_from_report(structure, data_name)
    if parsed is None:
        print(f"  [SKIP] {structure}: 未找到报告或解析失败，跳过绘图")
        return False

    report_name, partitions_str, exclude_args = parsed
    # 绘图脚本用报告里的实际名字（可能是短名）
    plot_structure = report_name

    cmd = ["python", "step6_1_1_partition_cv_plot.py",
           "--structure", plot_structure,       # 用实际文件名（短名）
           "--figsize", figsize,
           "--peak-method", "partition"]

    if partitions_str:
        cmd += ["--partitions", partitions_str]

    for ex in exclude_args:
        cmd += ["--exclude", ex]

    if exclude_args:
        cmd += ["--exclude-sort-by", "energy"]

    if extra_args:
        cmd += extra_args

    print(f"  绘图命令: {' '.join(cmd)}")

    try:
        import locale
        encoding = locale.getpreferredencoding() if sys.platform == 'win32' else 'utf-8'
        env = os.environ.copy()
        env.setdefault('PYTHONIOENCODING', 'utf-8')
        result = subprocess.run(cmd, capture_output=True, text=True,
                                encoding=encoding, errors='replace', env=env)
        if result.returncode == 0:
            # 推算图片路径（与 step6_1_1_partition_cv_plot.py 保持一致）
            if partitions_str:
                range_str = '_'.join(f"{seg}K" for seg in partitions_str.split(','))
                img_path = f"results/step6_1_1_partition_cv/{plot_structure}_partition_{range_str}_cv.png"
            else:
                img_path = f"results/step6_1_1_partition_cv/{plot_structure}_partition_cv.png"
            img_exists = os.path.exists(img_path)
            status = "[OK]" if img_exists else "[?]"
            print(f"  {status} 绘图成功: {img_path}")
            return img_path if img_exists else True
        else:
            print(f"  [FAIL] 绘图失败: {structure}")
            if result.stderr:
                print(f"  错误: {result.stderr[:400].strip()}")
            return False
    except Exception as e:
        print(f"  [ERROR] 绘图异常: {e}")
        return False


def check_data_file(structure, data_name=None):
    """
    检查数据文件是否存在。
    先用 structure 名查找，找不到再用 data_name（去掉 g-xxx- 前缀的短名）。
    返回实际找到的文件名（str），找不到返回 None。
    """
    base = f"results/step6_1_clustering/{{}}_lindemann-threshold_n2_clustered_data.csv"
    for name in dict.fromkeys([structure, data_name]):  # 去重，保持顺序
        if name and os.path.exists(base.format(name)):
            return name
    return None


def run_optimize_for_structure(structure, partitions, mode="suggest", threshold=60,
                               platform="windows", data_name=None):
    """
    对单个结构运行优化分析，自动附加 --partitions 参数。

    Parameters
    ----------
    structure  : str  CSV 中的结构名（用于报告/日志文件名）
    partitions : list[(float,float)]
    data_name  : str  数据文件实际使用的短名（去掉 g-xxx- 前缀），None 时与 structure 相同
    """
    print("\n" + "="*80)
    print(f"处理: {structure}")
    print("="*80)

    # 检查数据文件（先尝试 structure，再尝试 data_name）
    found_name = check_data_file(structure, data_name)
    if not found_name:
        tried = list(dict.fromkeys([structure, data_name]))
        print(f"  [SKIP] 数据文件不存在，跳过: {tried}")
        return False

    if found_name != structure:
        print(f"  [注意] 数据文件使用短名: {found_name}  (CSV名: {structure})")

    partitions_str = get_partition_str(partitions)
    print(f"  分区: {partitions_str}")

    # 构建命令：--structure 用实际文件名，报告/日志仍用 structure 全名
    cmd = [
        "python", "step6_1_1_8_auto_optimize_exclude.py",
        "--structure", found_name,
        "--mode", mode,
        "--threshold", str(threshold),
        "--platform", platform,
        "--partitions", partitions_str,
    ]

    print(f"  命令: {' '.join(cmd)}")

    # 运行命令并记录日志
    try:
        import locale
        encoding = locale.getpreferredencoding() if sys.platform == 'win32' else 'utf-8'
        env = os.environ.copy()
        env.setdefault('PYTHONIOENCODING', 'utf-8')
        result = subprocess.run(cmd, capture_output=True, text=True,
                                encoding=encoding, errors='replace', env=env)

        logs_dir = Path("logs")
        logs_dir.mkdir(parents=True, exist_ok=True)
        try:
            (logs_dir / f"{structure}.out.txt").write_text(
                result.stdout or "", encoding='utf-8', errors='replace')
            (logs_dir / f"{structure}.err.txt").write_text(
                result.stderr or "", encoding='utf-8', errors='replace')
        except Exception as we:
            print(f"  警告: 无法写入日志文件 ({we})")

        if result.returncode == 0:
            print(f"  [OK] 成功")
            # 打印子进程全部输出（对 GBK 终端做安全编码）
            if result.stdout:
                safe_out = result.stdout.encode(
                    sys.stdout.encoding or 'utf-8', errors='replace'
                ).decode(sys.stdout.encoding or 'utf-8', errors='replace')
                print(safe_out)
            return True
        else:
            print(f"  [FAIL] 失败 (查看 logs/{structure}.err.txt)")
            if result.stderr:
                print(f"  错误(摘要): {result.stderr[:400].strip()}")
            return False
    except Exception as e:
        print(f"  [ERROR] 异常: {e}")
        try:
            Path("logs").mkdir(parents=True, exist_ok=True)
            (Path("logs") / f"{structure}.err.txt").write_text(str(e), encoding='utf-8')
        except Exception:
            pass
        return False


def batch_process(partition_map, structures=None, mode="suggest", threshold=60,
                  platform="windows", auto_plot=True, figsize="10x8"):
    """
    批量处理多个结构。

    Parameters
    ----------
    partition_map : dict  由 load_partitions_from_csv() 返回的分区配置
    structures    : list  要处理的结构名列表；None 表示处理 partition_map 中所有结构
    auto_plot     : bool  筛选完成后自动调用绘图脚本（默认 True）
    figsize       : str   绘图尺寸，传给 --figsize（默认 '10x8'）
    """
    if structures is None:
        structures = list(partition_map.keys())

    # 过滤掉没有数据文件的结构
    available = [s for s in structures if s in partition_map]
    missing   = [s for s in structures if s not in partition_map]
    if missing:
        print(f"  [警告] 以下结构在 CSV 中找不到分区配置，跳过: {missing}")

    print("\n" + "="*80)
    print("批量优化排除点分析 (分区由 T1_T2_summary.csv 自动读取)")
    print("="*80)
    print(f"结构数量: {len(available)}")
    print(f"模式: {mode}  阈值: {threshold} meV  平台: {platform}")
    print("="*80)

    results = []
    for i, structure in enumerate(available, 1):
        print(f"\n[{i}/{len(available)}] {structure}")
        cfg = partition_map[structure]
        success = run_optimize_for_structure(
            structure, cfg['partitions'], mode, threshold, platform,
            data_name=cfg.get('data_name'),
        )
        # 报告可能用短名生成（g-xxx-前缀被去掉），两个都尝试
        data_name = cfg.get('data_name', structure)
        report_file = (f"{data_name}_EXCLUDE_RECOMMENDATIONS.md"
                       if os.path.exists(f"{data_name}_EXCLUDE_RECOMMENDATIONS.md")
                       else f"{structure}_EXCLUDE_RECOMMENDATIONS.md")
        results.append({
            "structure": structure,
            "data_name": data_name,
            "partitions": get_partition_str(cfg['partitions']),
            "success": success,
            "report": report_file if os.path.exists(report_file) else "",
            "out_log": f"logs/{structure}.out.txt",
            "err_log": f"logs/{structure}.err.txt",
        })

    # 汇总
    print("\n" + "="*80)
    print("批量处理汇总")
    print("="*80)
    success_count = sum(1 for r in results if r['success'])
    fail_count = len(results) - success_count
    print(f"总数: {len(results)}  成功: {success_count}  失败: {fail_count}")
    if fail_count > 0:
        print("\n失败的结构:")
        for r in results:
            if not r['success']:
                print(f"  - {r['structure']}")

    # 保存汇总 CSV
    try:
        import csv
        with open('batch_run_results.csv', 'w', newline='', encoding='utf-8-sig') as cf:
            writer = csv.DictWriter(cf, fieldnames=['structure', 'data_name', 'partitions', 'success', 'report', 'out_log', 'err_log'],
                                    extrasaction='ignore')
            writer.writeheader()
            writer.writerows(results)
        print("\n[OK] 批量运行结果已保存: batch_run_results.csv")
    except Exception as e:
        print(f"  无法写入汇总: {e}")

    print("\n[OK] 批量处理完成!")

    # ── 自动绘图 ──────────────────────────────────────────────────────
    if auto_plot:
        print("\n" + "="*80)
        print("批量绘图 (使用筛选报告中的参数)")
        print("="*80)
        plot_ok, plot_fail = 0, 0
        plot_images = []
        for r in results:
            structure = r['structure']
            if not r['success']:
                print(f"  [SKIP] {structure}: 筛选失败，跳过绘图")
                continue
            print(f"\n  → {structure}")
            img = run_plot_for_structure(structure, figsize=figsize,
                                         data_name=r.get('data_name'))
            if img:
                plot_ok += 1
                if isinstance(img, str):
                    plot_images.append((structure, img))
            else:
                plot_fail += 1

        print(f"\n{'='*80}")
        print(f"绘图汇总: 共 {plot_ok + plot_fail} 个结构  [成功 {plot_ok}]  [失败 {plot_fail}]")
        if plot_images:
            print(f"\n生成的图片 ({len(plot_images)} 个):")
            for name, path in plot_images:
                print(f"  {name:<25}  {path}")
        print("="*80)
    else:
        print("\n[提示] 绘图已跳过（使用 --no-plot 禁用，去掉可自动绘图）")


def generate_summary_table(partition_map):
    """生成所有结构的优化效果汇总表"""
    print("\n" + "="*80)
    print("生成优化效果汇总表")
    print("="*80)

    summary_data = []
    for structure, cfg in partition_map.items():
        # 先用全名找，找不到再用短名（data_name，去掉 g-xxx- 前缀）
        data_name = cfg.get('data_name', structure)
        for name in dict.fromkeys([structure, data_name]):
            report_file = f"{name}_EXCLUDE_RECOMMENDATIONS.md"
            if os.path.exists(report_file):
                break
        else:
            print(f"  跳过 {structure}: 报告不存在 (尝试过: {structure}, {data_name})")
            continue
        
        # 读取实际数据文件，获取温度范围
        actual_t_range = "N/A"
        runs_per_t_str = "N/A"
        data_sources   = "N/A"
        try:
            found_name = check_data_file(structure, data_name)
            if found_name:
                data_csv = (f"results/step6_1_clustering/"
                            f"{found_name}_lindemann-threshold_n2_clustered_data.csv")
                df_data = pd.read_csv(data_csv)
                # 温度列可能叫 temp / T / temperature
                t_col = next(
                    (c for c in df_data.columns
                     if c.strip().lower() in ('t', 'temperature', 'temp')),
                    None
                )
                if t_col:
                    t_min = int(df_data[t_col].min())
                    t_max = int(df_data[t_col].max())
                    n_pts = len(df_data)
                    actual_t_range = f"{t_min}-{t_max}K ({n_pts}pts)"

                    # 每个温度点模拟次数（min-max 或单值）
                    runs_per_t = df_data.groupby(t_col).size()
                    r_min, r_max = int(runs_per_t.min()), int(runs_per_t.max())
                    runs_per_t_str = str(r_min) if r_min == r_max else f"{r_min}-{r_max}"

                    # 上游数据来源（data_source 列）
                    if 'data_source' in df_data.columns:
                        src_list = sorted(df_data['data_source'].dropna().unique())
                        data_sources = ','.join(str(s) for s in src_list)
        except Exception:
            pass  # 找不到数据文件时保持 N/A

        # 读取报告，提取关键指标
        try:
            with open(report_file, 'r', encoding='utf-8') as f:
                content = f.read()

            # 从报告内容中用正则提取 Cv 数值（可选，找不到就记 N/A）
            cv_vals = re.findall(r'Cv_net\s*[=≈:]\s*([\d.]+)', content)
            cv_summary = ', '.join(cv_vals) if cv_vals else "N/A"

            summary_data.append({
                "结构": structure,
                "分区数": cfg['n_partitions'],
                "分区范围": get_partition_str(cfg['partitions']),
                "T1_lindemann": cfg['T1_raw'],
                "T2_B": cfg['T2'],
                "数据温度范围": actual_t_range,
                "runs/T点": runs_per_t_str,
                "数据来源": data_sources,
                "Cv_net(meV/K)": cv_summary,
                "报告": "[OK]",
            })
        except Exception as e:
            print(f"  错误 {structure}: {e}")

    # 生成 CSV
    if summary_data:
        df = pd.DataFrame(summary_data)
        output_file = "batch_optimization_summary.csv"
        df.to_csv(output_file, index=False, encoding='utf-8-sig')
        print(f"\n[OK] 汇总表已保存: {output_file}")
        print(df.to_string(index=False))
    else:
        print("  无数据生成汇总表")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="批量优化排除点分析工具（分区自动从 T1_T2_summary.csv 读取）"
    )
    parser.add_argument('--structures', nargs='+',
                        help='指定要处理的结构列表（空格分隔）；默认处理 CSV 中全部结构')
    parser.add_argument('--mode', default='suggest',
                        choices=['suggest', 'test', 'compare', 'all'],
                        help='分析模式 (默认: suggest)')
    parser.add_argument('--threshold', type=int, default=60,
                        help='残差阈值(meV) (默认: 60)')
    parser.add_argument('--platform', default='windows',
                        choices=['windows', 'linux', 'mac'],
                        help='目标平台 (默认: windows)')
    parser.add_argument('--partition-type',
                        choices=['2', '3', 'all'], default='all',
                        help='按分区数筛选 (2=2分区, 3=3分区, all=全部)')
    parser.add_argument('--csv-file', default=DEFAULT_CSV,
                        help=f'T1/T2 汇总 CSV 文件路径 (默认: {DEFAULT_CSV})')
    parser.add_argument('--t2-col', default=DEFAULT_T2_COL,
                        choices=['T2_B', 'T2_Bprime'],
                        help=f'使用哪一列作为 T2 (默认: {DEFAULT_T2_COL})')
    parser.add_argument('--summary-only', action='store_true',
                        help='仅生成汇总表，不运行分析')
    parser.add_argument('--show-config', action='store_true',
                        help='显示从 CSV 读取的分区配置信息后退出')
    parser.add_argument('--no-plot', action='store_true',
                        help='筛选完成后不自动绘图（默认筛选后自动绘图）')
    parser.add_argument('--figsize', default='10x8',
                        help='绘图尺寸，传给 --figsize（默认: 10x8）')

    args = parser.parse_args()

    # ── 从 CSV 加载分区配置 ──────────────────────────────────────────
    partition_map = load_partitions_from_csv(args.csv_file, args.t2_col)
    structures_2 = [s for s, c in partition_map.items() if c['n_partitions'] == 2]
    structures_3 = [s for s, c in partition_map.items() if c['n_partitions'] == 3]

    print(f"\n[配置] CSV: {args.csv_file}  T2列: {args.t2_col}")
    print(f"  2分区结构: {len(structures_2)}个  3分区结构: {len(structures_3)}个")

    # ── --show-config：打印分区配置后退出 ────────────────────────────
    if args.show_config:
        print("\n" + "="*90)
        print(f"分区配置 (来源: {args.csv_file}, T2列: {args.t2_col})")
        print("="*90)
        print(f"{'结构':<18} {'分区数':>5}  {'T1_raw':>8}  {'T1_boundary':>11}  {'T2':>6}  {'分区范围'}")
        print("-"*90)
        for s, c in partition_map.items():
            print(f"{s:<18} {c['n_partitions']:>5}  {c['T1_raw']:>8.0f}  "
                  f"{c['T1_boundary']:>11}  {c['T2']:>6}  "
                  f"{get_partition_str(c['partitions'])}")
        print("="*90)
        sys.exit(0)

    # ── --summary-only：只生成汇总表 ─────────────────────────────────
    if args.summary_only:
        generate_summary_table(partition_map)
        sys.exit(0)

    # ── 确定目标结构列表 ──────────────────────────────────────────────
    if args.structures:
        structures = args.structures
    elif args.partition_type == '2':
        structures = structures_2
    elif args.partition_type == '3':
        structures = structures_3
    else:
        structures = list(partition_map.keys())

    # ── 运行批量处理 ──────────────────────────────────────────────────
    batch_process(partition_map, structures, args.mode, args.threshold, args.platform,
                  auto_plot=not args.no_plot, figsize=args.figsize)

    # ── 生成汇总表 ────────────────────────────────────────────────────
    if args.mode in ('all', 'suggest'):
        generate_summary_table(partition_map)
