#!/usr/bin/env python3
"""
批量优化排除点分析工具
对多个结构进行自动化的排除点优化分析

分区策略:
- T1(固相边界): 可选 Lindemann阈值(δ=0.1) 或 kmeans聚类
- T2(液相边界): 基于kmeans聚类分析确定
- 数据来源: results/step6_1_clustering/t1_t2_kmeans_comparison_19structures.csv
"""

import os
import sys
import subprocess
import pandas as pd
from pathlib import Path

# ============================================================================
# 分区配置 - 支持两种边界类型: kmeans (默认) 和 lindemann
# ============================================================================

# 边界类型配置
# T1_kmeans, T2_kmeans: 基于kmeans聚类的边界温度
# T1_lindemann: 基于Lindemann指数δ=0.1的边界温度
# 注意: T2只有kmeans版本，因为Lindemann主要用于确定熔点（T1）

PARTITION_CONFIG_FULL = {
    # ========== 2分区结构 (6个): 使用kmeans时有T2，使用Lindemann时无T2 ==========
    "Sn1Pt2O1": {
        "T1_kmeans": 750, "T2_kmeans": 1450, "T1_lindemann": 1800,
        "type_kmeans": "3-partition", "type_lindemann": "2-partition"
    },
    "Pt2Sn2O1": {
        "T1_kmeans": 750, "T2_kmeans": 1350, "T1_lindemann": 1500,
        "type_kmeans": "3-partition", "type_lindemann": "2-partition"
    },
    "Sn6Pt5O2": {
        "T1_kmeans": 750, "T2_kmeans": 1350, "T1_lindemann": 600,
        "type_kmeans": "3-partition", "type_lindemann": "2-partition"
    },
    "Pt6Sn5O2": {
        "T1_kmeans": 550, "T2_kmeans": 750, "T1_lindemann": 800,
        "type_kmeans": "3-partition", "type_lindemann": "2-partition"
    },
    "Pt6Sn6O3": {
        "T1_kmeans": 450, "T2_kmeans": 650, "T1_lindemann": 700,
        "type_kmeans": "3-partition", "type_lindemann": "2-partition"
    },
    "Sn7Pt6O4": {
        "T1_kmeans": 400, "T2_kmeans": 700, "T1_lindemann": 800,
        "type_kmeans": "3-partition", "type_lindemann": "2-partition"
    },
    
    # ========== 3分区结构 (13个): 两种边界类型都有T2 ==========
    "Pt3Sn2O1": {
        "T1_kmeans": 750, "T2_kmeans": 1200, "T1_lindemann": 1100, "T2_lindemann": 1200,
        "type_kmeans": "3-partition", "type_lindemann": "3-partition"
    },
    "Sn3O2Pt2": {
        "T1_kmeans": 750, "T2_kmeans": 1400, "T1_lindemann": 1300, "T2_lindemann": 1500,
        "type_kmeans": "3-partition", "type_lindemann": "3-partition"
    },
    "O3Sn4Pt2": {
        "T1_kmeans": 700, "T2_kmeans": 1250, "T1_lindemann": 900, "T2_lindemann": 1200,
        "type_kmeans": "3-partition", "type_lindemann": "3-partition"
    },
    "Pt3Sn3O2": {
        "T1_kmeans": 750, "T2_kmeans": 1300, "T1_lindemann": 1100, "T2_lindemann": 1200,
        "type_kmeans": "3-partition", "type_lindemann": "3-partition"
    },
    "Sn3Pt4O1": {
        "T1_kmeans": 700, "T2_kmeans": 1250, "T1_lindemann": 600, "T2_lindemann": 900,
        "type_kmeans": "3-partition", "type_lindemann": "3-partition"
    },
    "Pt5Sn3O1": {
        "T1_kmeans": 700, "T2_kmeans": 1200, "T1_lindemann": 600, "T2_lindemann": 800,
        "type_kmeans": "3-partition", "type_lindemann": "3-partition"
    },
    "Pt5Sn4O1": {
        "T1_kmeans": 750, "T2_kmeans": 1250, "T1_lindemann": 800, "T2_lindemann": 900,
        "type_kmeans": "3-partition", "type_lindemann": "3-partition"
    },
    "O2Pt4Sn6": {
        "T1_kmeans": 750, "T2_kmeans": 1300, "T1_lindemann": 600, "T2_lindemann": 900,
        "type_kmeans": "3-partition", "type_lindemann": "3-partition"
    },
    "Sn7Pt4O3": {
        "T1_kmeans": 800, "T2_kmeans": 1300, "T1_lindemann": 700, "T2_lindemann": 900,
        "type_kmeans": "3-partition", "type_lindemann": "3-partition"
    },
    "O3Pt5Sn7": {
        "T1_kmeans": 700, "T2_kmeans": 1200, "T1_lindemann": 800, "T2_lindemann": 900,
        "type_kmeans": "3-partition", "type_lindemann": "3-partition"
    },
    "Pt7Sn5O1": {
        "T1_kmeans": 650, "T2_kmeans": 1150, "T1_lindemann": 600, "T2_lindemann": 800,
        "type_kmeans": "3-partition", "type_lindemann": "3-partition"
    },
    "Pt7Sn6O1": {
        "T1_kmeans": 600, "T2_kmeans": 1150, "T1_lindemann": 600, "T2_lindemann": 900,
        "type_kmeans": "3-partition", "type_lindemann": "3-partition"
    },
    "O2Pt7Sn7": {
        "T1_kmeans": 750, "T2_kmeans": 1275, "T1_lindemann": 600, "T2_lindemann": 750,
        "type_kmeans": "3-partition", "type_lindemann": "3-partition"
    },
}

# 默认使用 kmeans 边界
DEFAULT_BOUNDARY_TYPE = "kmeans"


def get_partition_config(boundary_type="kmeans"):
    """
    根据边界类型生成分区配置
    
    Args:
        boundary_type: "kmeans" (默认) 或 "lindemann"
    
    Returns:
        PARTITION_CONFIG 字典
        
    分区策略 (3分区结构):
        - T1分区 (固相): 200K → T1边界
        - T_mid分区 (过渡区): T1边界 → T2边界
        - T2分区 (液相): T2边界 → 1800K
    
    分区策略 (2分区结构, 无T2):
        - T1分区 (固相): 200K → T1边界
        - T_mid分区: T1边界 → 1800K
    """
    config = {}
    
    for structure, full_config in PARTITION_CONFIG_FULL.items():
        T1_key = f"T1_{boundary_type}"
        T2_key = f"T2_{boundary_type}" if boundary_type == "lindemann" else "T2_kmeans"
        type_key = f"type_{boundary_type}"
        
        T1_val = full_config.get(T1_key)
        T2_val = full_config.get(T2_key)
        partition_type = full_config.get(type_key, "3-partition")
        
        # T1范围 (固相): 200K → T1边界
        # 使用T1边界本身作为范围终点（包含边界点）
        T1_range = (200, T1_val) if T1_val else None
        
        # T_mid范围 (过渡区): T1边界 → T2边界
        if T1_val and T2_val:
            T_mid_range = (T1_val, T2_val)
        elif T1_val:
            # 无T2时，过渡区延伸到1800K
            T_mid_range = (T1_val, 1800)
        else:
            T_mid_range = None
        
        # T2范围 (液相): T2边界 → 1800K
        if T2_val:
            T2_range = (T2_val, 1800)
        else:
            T2_range = None
        
        config[structure] = {
            "T1": T1_range,           # 固相: 200K → T1
            "T_mid": T_mid_range,     # 过渡区: T1 → T2
            "T2": T2_range,           # 液相: T2 → 1800K
            "type": partition_type,
            "T1_boundary": T1_val,
            "T2_boundary": T2_val
        }
    
    return config


# 生成默认配置（兼容旧代码）
PARTITION_CONFIG = get_partition_config(DEFAULT_BOUNDARY_TYPE)

# 结构列表
STRUCTURES_2_PARTITION = [s for s, c in PARTITION_CONFIG.items() if c["type"] == "2-partition"]
STRUCTURES_3_PARTITION = [s for s, c in PARTITION_CONFIG.items() if c["type"] == "3-partition"]
ALL_STRUCTURES = list(PARTITION_CONFIG.keys())


def check_data_file(structure):
    """检查数据文件是否存在"""
    data_file = f"results/step6_1_clustering/{structure}_lindemann-threshold_n2_clustered_data.csv"
    return os.path.exists(data_file)


def run_optimize_for_structure(structure, mode="all", threshold=60, platform="windows"):
    """对单个结构运行优化分析"""
    print("\n" + "="*80)
    print(f"处理: {structure}")
    print("="*80)
    
    # 检查数据文件
    if not check_data_file(structure):
        print(f"  [SKIP] 数据文件不存在,跳过: {structure}")
        return False
    
    # 构建命令
    cmd = [
        "python", "step6_1_1_8_auto_optimize_exclude.py",
        "--structure", structure,
        "--mode", mode,
        "--threshold", str(threshold),
        "--platform", platform
    ]
    
    print(f"  命令: {' '.join(cmd)}")
    
    # 运行命令并记录日志(始终继续，不抛出)
    try:
        # Windows中文环境下使用系统默认编码
        import sys
        import locale
        encoding = locale.getpreferredencoding() if sys.platform == 'win32' else 'utf-8'

        # Ensure child Python uses UTF-8 for stdout/stderr to avoid encoding errors on Windows
        env = os.environ.copy()
        env.setdefault('PYTHONIOENCODING', 'utf-8')
        result = subprocess.run(cmd, capture_output=True, text=True, encoding=encoding, errors='replace', env=env)

        # 写入日志文件
        logs_dir = Path("logs")
        logs_dir.mkdir(parents=True, exist_ok=True)
        out_file = logs_dir / f"{structure}.out.txt"
        err_file = logs_dir / f"{structure}.err.txt"
        try:
            with open(out_file, 'w', encoding='utf-8') as fo:
                fo.write(result.stdout or "")
            with open(err_file, 'w', encoding='utf-8') as fe:
                fe.write(result.stderr or "")
        except Exception:
            # 最坏情况：仍然不要中断批量流程
            print("  警告: 无法写入日志文件")

        if result.returncode == 0:
            print(f"  [OK] 成功")
            return True
        else:
            print(f"  [FAIL] 失败 (查看 logs/{structure}.err.txt 获取详细信息)")
            # 打印一小段 stderr 便于实时查看
            if result.stderr:
                print(f"  错误(摘要): {result.stderr[:400].strip()}")
            return False
    except Exception as e:
        print(f"  [ERROR] 异常: {str(e)}")
        # 记录异常到日志
        try:
            logs_dir = Path("logs")
            logs_dir.mkdir(parents=True, exist_ok=True)
            with open(logs_dir / f"{structure}.err.txt", 'w', encoding='utf-8') as fe:
                fe.write(str(e))
        except Exception:
            pass
        return False


def generate_partition_commands(structure, exclude_dict, config):
    """生成分区绘图命令"""
    partition_info = PARTITION_CONFIG.get(structure)
    if not partition_info:
        print(f"  警告: {structure} 没有分区配置")
        return []
    
    commands = []
    
    # T1分区命令
    if partition_info["T1"]:
        T1_min, T1_max = partition_info["T1"]
        exclude_args = " ".join([f'"{k}K:{",".join(map(str, v))}"' for k, v in exclude_dict.items()
                                 if T1_min <= k <= T1_max])
        
        cmd = f"""python step6_1_1_partition_cv_plot.py `
    --structure {structure} `
    --partitions {T1_min}-{T1_max}"""
        
        if exclude_args:
            cmd += f""" `
    --exclude {exclude_args} `
    --exclude-sort-by energy"""
        
        cmd += """ `
    --y-ticks 0,2,4 --cv-ticks 3,4,5,6 --figsize 10x8 --peak-method partition"""
        
        commands.append(("T1", cmd))
    
    # T2分区命令
    if partition_info["T2"]:
        T2_min, T2_max = partition_info["T2"]
        exclude_args = " ".join([f'"{k}K:{",".join(map(str, v))}"' for k, v in exclude_dict.items()
                                 if T2_min <= k <= T2_max])
        
        cmd = f"""python step6_1_1_partition_cv_plot.py `
    --structure {structure} `
    --partitions {T2_min}-{T2_max}"""
        
        if exclude_args:
            cmd += f""" `
    --exclude {exclude_args} `
    --exclude-sort-by energy"""
        
        cmd += """ `
    --y-ticks 0,2,4 --cv-ticks 3,4,5,6 --figsize 10x8 --peak-method partition"""
        
        commands.append(("T2", cmd))
    
    return commands


def batch_process(structures, mode="all", threshold=60, platform="windows"):
    """批量处理多个结构"""
    print("\n" + "="*80)
    print(f"批量优化排除点分析")
    print("="*80)
    print(f"结构数量: {len(structures)}")
    print(f"模式: {mode}")
    print(f"阈值: {threshold} meV")
    print(f"平台: {platform}")
    print("="*80)
    
    results = []

    for i, structure in enumerate(structures, 1):
        print(f"\n[{i}/{len(structures)}] {structure}")
        success = run_optimize_for_structure(structure, mode, threshold, platform)
        # 记录结果并关联日志/报告路径
        report_file = f"{structure}_EXCLUDE_RECOMMENDATIONS.md"
        out_log = f"logs/{structure}.out.txt"
        err_log = f"logs/{structure}.err.txt"
        results.append({
            "structure": structure,
            "success": success,
            "report": report_file if os.path.exists(report_file) else "",
            "out_log": out_log if os.path.exists(out_log) else "",
            "err_log": err_log if os.path.exists(err_log) else "",
        })
    
    # 汇总
    print("\n" + "="*80)
    print("批量处理汇总")
    print("="*80)
    success_count = sum(1 for r in results if r['success'])
    fail_count = len(results) - success_count
    
    print(f"总数: {len(results)}")
    print(f"成功: {success_count}")
    print(f"失败: {fail_count}")
    
    if fail_count > 0:
        print("\n失败的结构:")
        for r in results:
            if not r['success']:
                print(f"  - {r['structure']} (查看 {r['err_log'] or 'logs/<structure>.err.txt'})")

    # 保存批量运行日志汇总CSV
    try:
        import csv
        summary_file = 'batch_run_results.csv'
        with open(summary_file, 'w', newline='', encoding='utf-8-sig') as cf:
            writer = csv.DictWriter(cf, fieldnames=['structure', 'success', 'report', 'out_log', 'err_log'])
            writer.writeheader()
            for r in results:
                writer.writerow(r)
        print(f"\n[OK] 批量运行结果已保存: {summary_file}")
    except Exception as e:
        print(f"  无法写入批量运行汇总: {e}")
    
    print("\n[OK] 批量处理完成!")
    print(f"报告目录: 当前工作目录下生成 *_EXCLUDE_RECOMMENDATIONS.md")


def generate_summary_table():
    """生成所有结构的优化效果汇总表"""
    print("\n" + "="*80)
    print("生成优化效果汇总表")
    print("="*80)
    
    summary_data = []
    
    for structure in list(PARTITION_CONFIG.keys()):
        report_file = f"{structure}_EXCLUDE_RECOMMENDATIONS.md"
        
        if not os.path.exists(report_file):
            print(f"  跳过 {structure}: 报告不存在")
            continue
        
        # 读取报告提取关键指标
        try:
            with open(report_file, 'r', encoding='utf-8') as f:
                content = f.read()
                
                # 简单提取(实际应该用正则表达式)
                partition_type = PARTITION_CONFIG[structure]["type"]
                T1_range = PARTITION_CONFIG[structure]["T1"]
                T2_range = PARTITION_CONFIG[structure]["T2"]
                
                summary_data.append({
                    "结构": structure,
                    "分区类型": partition_type,
                    "T1范围": f"{T1_range[0]}-{T1_range[1]}K" if T1_range else "-",
                    "T2范围": f"{T2_range[0]}-{T2_range[1]}K" if T2_range else "-",
                    "报告": "[OK]"
                })
        except Exception as e:
            print(f"  错误 {structure}: {e}")
    
    # 生成CSV
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
    
    parser = argparse.ArgumentParser(description="批量优化排除点分析工具")
    parser.add_argument('--structures', nargs='+', 
                        help='指定要处理的结构列表(空格分隔)')
    parser.add_argument('--mode', default='all', 
                        choices=['suggest', 'test', 'compare', 'all'],
                        help='分析模式 (默认: all)')
    parser.add_argument('--threshold', type=int, default=60,
                        help='残差阈值(meV) (默认: 60)')
    parser.add_argument('--platform', default='windows', 
                        choices=['windows', 'linux', 'mac'],
                        help='目标平台 (默认: windows)')
    parser.add_argument('--partition-type', 
                        choices=['2', '3', 'all'],
                        help='按分区类型筛选 (2=2分区, 3=3分区, all=全部)')
    parser.add_argument('--boundary-type', default='kmeans',
                        choices=['kmeans', 'lindemann'],
                        help='边界类型: kmeans(默认,基于聚类) 或 lindemann(基于δ=0.1)')
    parser.add_argument('--summary-only', action='store_true',
                        help='仅生成汇总表,不运行分析')
    parser.add_argument('--show-config', action='store_true',
                        help='显示当前分区配置信息')
    
    args = parser.parse_args()
    
    # 根据边界类型生成配置
    partition_config = get_partition_config(args.boundary_type)
    structures_2_partition = [s for s, c in partition_config.items() if c["type"] == "2-partition"]
    structures_3_partition = [s for s, c in partition_config.items() if c["type"] == "3-partition"]
    all_structures = list(partition_config.keys())
    
    print(f"\n[配置] 边界类型: {args.boundary_type.upper()}")
    print(f"  - 2分区结构: {len(structures_2_partition)}个")
    print(f"  - 3分区结构: {len(structures_3_partition)}个")
    
    # 显示配置信息
    if args.show_config:
        print("\n" + "="*110)
        print(f"分区配置详情 (当前边界类型: {args.boundary_type})")
        print("="*110)
        print(f"{'结构':<15} {'类型':<12} {'T1(固相)':<15} {'T_mid(过渡区)':<18} {'T2(液相)':<15}")
        print("-"*110)
        for s, c in partition_config.items():
            t1_range = f"{c['T1'][0]}-{c['T1'][1]}K" if c['T1'] else "N/A"
            t_mid_range = f"{c['T_mid'][0]}-{c['T_mid'][1]}K" if c['T_mid'] else "N/A"
            t2_range = f"{c['T2'][0]}-{c['T2'][1]}K" if c['T2'] else "N/A"
            print(f"{s:<15} {c['type']:<12} {t1_range:<15} {t_mid_range:<18} {t2_range:<15}")
        
        # 显示两种边界类型的对比
        print("\n" + "="*110)
        print("两种边界类型对比 (都基于Lindemann指数，但计算方法不同)")
        print("="*110)
        print(f"{'结构':<15} | {'T1_kmeans':>9} | {'T1_lindemann':>12} | {'T2_kmeans':>9} | {'T2_lindemann':>12} | {'T1差异':>6}")
        print("-"*110)
        for structure, config in PARTITION_CONFIG_FULL.items():
            t1_k = config.get('T1_kmeans', 0)
            t1_l = config.get('T1_lindemann', 0)
            t2_k = config.get('T2_kmeans', 0)
            t2_l = config.get('T2_lindemann', 0)
            diff = abs(t1_k - t1_l)
            t2_l_str = f"{t2_l}K" if t2_l else "N/A"
            print(f"{structure:<15} | {t1_k:>8}K | {t1_l:>11}K | {t2_k:>8}K | {t2_l_str:>12} | {diff:>5}K")
        
        print("\n" + "="*110)
        print("说明:")
        print("  - kmeans: 聚类算法自动识别相变温度 (可能是非100K整数倍，如700K、650K)")
        print("  - lindemann: Lindemann指数首次超过δ=0.1阈值的温度")
        print("  - T1差异: 两种方法计算的T1边界差值")
        print("  - T2_lindemann为N/A: 该结构在lindemann模式下为2分区(无过渡区)")
        print("="*110)
        sys.exit(0)
    
    # 确定要处理的结构列表
    if args.summary_only:
        generate_summary_table()
        sys.exit(0)
    
    if args.structures:
        structures = args.structures
    elif args.partition_type == '2':
        structures = structures_2_partition
    elif args.partition_type == '3':
        structures = structures_3_partition
    else:
        structures = structures_2_partition + structures_3_partition
    
    # 运行批量处理
    batch_process(structures, args.mode, args.threshold, args.platform)
    
    # 生成汇总表
    if args.mode == 'all':
        generate_summary_table()
