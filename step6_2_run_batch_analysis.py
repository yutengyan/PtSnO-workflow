#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Author: Yutengyan 1120798743@qq.com
Date: 2025-10-28 11:59:33
LastEditors: yutengyan 1120798743@qq.com
LastEditTime: 2025-11-29 10:31:34
FilePath: step6_2_run_batch_analysis.py
Description: 批量运行所有Pt-Sn结构的聚类分析
'''
"""
批量运行所有Pt-Sn结构的聚类分析
调用 step6_1_clustering_analysis.py 对每个结构进行分析

使用方法:
  默认数据:  python step6_2_run_batch_analysis.py
  50K数据:   python step6_2_run_batch_analysis.py --data-source 50K
"""
import subprocess
import sys
import argparse
from datetime import datetime

# 默认结构列表 (100K温度间隔数据)
structures_default = [
    # Air系列 (气相纳米团簇) - 不需要扣除载体热容
    "Air68", "Air86",
    
    # # Cv系列 (Pt6Sn8O4, 即 g-1535-Sn8Pt6O4)
    "Cv",
    # PtxSnyOz系列
    "O2pt4sn6", "O2pt7sn7", "O2sn8pt7", "O3pt5sn7", "O3sn4pt2",
    "O3sn5pt3", "O4pt3sn6", "Pt2sn2o1", "Pt3sn2o1", "Pt3sn3o2",
    "Pt3sn5", "Pt4sn4", "Pt5sn3", "Pt5sn3o1", "Pt5sn4o1",
    "Pt6sn1", "Pt6sn2", "Pt6sn3", "Pt6sn4", "Pt6sn5",
    "Pt6sn6", "Pt6sn7", "Pt6sn8", "Pt6sn9", "Pt7sn1",
    "Pt7sn5o1", "Pt7sn6o1", "Pt7sn9o4", "Pt8sn0", "Pt8sn1",
    "Pt8sn10", "Pt8sn2", "Pt8sn3", "Pt8sn4", "Pt8sn5",
    "Pt8sn6", "Pt8sn7", "Pt8sn8", "Pt8sn9", "Sn10pt7o4",
    "Sn1pt2o1", "Sn3o2pt2", "Sn3pt4o1", "Sn4pt3o1", "Sn5o2pt4",
    "Sn5o4pt2", "Sn6pt5o2", "Sn7pt4o3"
]

# 50K数据结构列表 (50K温度间隔数据)
structures_50K = [
    "Pt6Sn8",   # 160点, 16温度 × 10 runs
    "Pt8Sn6",   # 83点, 9温度 × ~10 runs  
    "O2Pt7Sn7", # 120点, 18+温度 × ~6-7 runs
]

def main():
    # 解析命令行参数
    parser = argparse.ArgumentParser(description='批量运行聚类分析')
    parser.add_argument('--data-source', type=str, default='default',
                        choices=['default', '50K'],
                        help='数据源: default (100K间隔) 或 50K (50K间隔)')
    parser.add_argument('--n-partitions', type=int, default=2,
                        help='分区数 (默认: 2)')
    parser.add_argument('--use-energy', action='store_true', default=True,
                        help='使用能量特征 (默认开启)')
    parser.add_argument('--use-d-value', action='store_true', default=False,
                        help='使用扩散系数D值特征')
    args = parser.parse_args()
    
    # 根据数据源选择结构列表
    if args.data_source == '50K':
        structures = structures_50K
        feature_desc = "Temperature + Lindemann-δ + Energy"
    else:
        structures = structures_default
        feature_desc = "Temperature + Lindemann-δ + Energy + Diffusion-D" if args.use_d_value else "Temperature + Lindemann-δ + Energy"
    
    total = len(structures)
    success = []
    failed = []
    
    print("=" * 80)
    print(f"开始批量分析 {total} 个结构")
    print(f"数据源: {args.data_source}")
    print(f"特征空间: {feature_desc}")
    print(f"分区数: {args.n_partitions}")
    print(f"开始时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 80)
    print()
    
    for idx, structure in enumerate(structures, 1):
        print(f"[{idx}/{total}] 正在分析: {structure}")
        print("-" * 80)
        
        try:
            cmd = [
                sys.executable,
                "step6_1_clustering_analysis.py",
                "--structure", structure,
                "--n-partitions", str(args.n_partitions),
            ]
            
            if args.use_energy:
                cmd.append("--use-energy")
            if args.use_d_value and args.data_source == 'default':
                cmd.append("--use-d-value")  # 50K数据暂不支持D值
            
            result = subprocess.run(cmd, capture_output=False, text=True)
            
            if result.returncode == 0:
                success.append(structure)
                print(f"[✓] {structure} 分析完成\n")
            else:
                failed.append(structure)
                print(f"[✗] {structure} 分析失败 (退出码: {result.returncode})\n")
                
        except Exception as e:
            failed.append(structure)
            print(f"[✗] {structure} 分析出错: {e}\n")
    
    print("=" * 80)
    print("分析完成统计")
    print("=" * 80)
    print(f"总计: {total} 个结构")
    print(f"成功: {len(success)} 个")
    print(f"失败: {len(failed)} 个")
    print(f"结束时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    if failed:
        print("\n失败的结构列表:")
        for f in failed:
            print(f"  - {f}")
    
    print("\n所有结果保存在: results/step6_1_clustering/")
    print("=" * 80)

if __name__ == "__main__":
    main()
