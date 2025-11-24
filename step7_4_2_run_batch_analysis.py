#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Author: Yutengyan 1120798743@qq.com
Date: 2025-10-28 11:59:33
LastEditors: Yutengyan 1120798743@qq.com
LastEditTime: 2025-11-21 
FilePath: step7_4_2_run_batch_analysis.py
Description: 批量运行所有Pt-Sn结构的聚类分析
'''
"""
批量运行所有Pt-Sn结构的聚类分析
"""
import subprocess
import sys
from datetime import datetime

structures = [
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

def main():
    total = len(structures)
    success = []
    failed = []
    
    print("=" * 80)
    print(f"开始批量分析 {total} 个结构（使用4D特征）")
    print(f"特征空间: Temperature + Lindemann-δ + Energy + Diffusion-D")
    print(f"开始时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 80)
    print()
    
    for idx, structure in enumerate(structures, 1):
        print(f"[{idx}/{total}] 正在分析: {structure}")
        print("-" * 80)
        
        try:
            cmd = [
                sys.executable,
                "step7_4_2_clustering_analysis.py",
                "--structure", structure,
                "--n-partitions", "2",
                "--use-energy",      # 添加能量特征
                "--use-d-value"      # 添加扩散系数D值特征
            ]
            
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
    
    print("\n所有结果保存在: results/step7_4_2_clustering/")
    print("=" * 80)

if __name__ == "__main__":
    main()
