#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
快速验证workflow数据路径配置
检查所有数据文件是否可访问
"""

from pathlib import Path
import sys

def check_path(path, description):
    """检查路径是否存在"""
    if path.exists():
        if path.is_file():
            size_mb = path.stat().st_size / (1024 * 1024)
            print(f"  ✅ {description}: {path.name} ({size_mb:.2f} MB)")
        else:
            count = len(list(path.rglob('*'))) if path.is_dir() else 0
            print(f"  ✅ {description}: {path.name}/ ({count} 项)")
        return True
    else:
        print(f"  ❌ {description}: {path} (不存在!)")
        return False

def main():
    print("=" * 70)
    print("workflow 数据路径验证")
    print("=" * 70)
    print()
    
    BASE_DIR = Path(__file__).parent
    all_ok = True
    
    # 1. 检查GMX MSD数据
    print("📁 Step 1-5: GMX MSD数据")
    gmx_dirs = [
        BASE_DIR / 'data' / 'gmx_msd' / 'collected_gmx_msd',
        BASE_DIR / 'data' / 'gmx_msd' / 'gmx_msd_results_20251015_184626_collected'
    ]
    for gmx_dir in gmx_dirs:
        all_ok &= check_path(gmx_dir, "GMX目录")
    print()
    
    # 2. 检查LAMMPS能量数据
    print("📁 Step 6: LAMMPS能量数据")
    energy_files = [
        BASE_DIR / 'data' / 'lammps_energy' / 'energy_master_20251016_121110.csv',
        BASE_DIR / 'data' / 'lammps_energy' / 'energy_master_20251021_134929.csv',
        BASE_DIR / 'data' / 'lammps_energy' / 'sup' / 'energy_master_20251021_151520.csv',
    ]
    for energy_file in energy_files:
        all_ok &= check_path(energy_file, "能量文件")
    print()
    
    # 3. 检查Lindemann数据
    print("📁 Step 7: Lindemann指数数据")
    lindemann_dir = BASE_DIR / 'data' / 'lindemann'
    all_ok &= check_path(lindemann_dir, "Lindemann目录")
    if lindemann_dir.exists():
        for lind_file in lindemann_dir.glob('lindemann_*.csv'):
            all_ok &= check_path(lind_file, "  Lindemann文件")
    print()
    
    # 4. 检查配位数数据
    print("📁 Step 7.5: 配位数/Q6数据")
    coord_dir = BASE_DIR / 'data' / 'coordination' / 'coordination_time_series_results_sample_20251026_200908'
    all_ok &= check_path(coord_dir, "配位数目录")
    print()
    
    # 5. 检查results目录
    print("📁 输出目录")
    results_dir = BASE_DIR / 'results'
    if not results_dir.exists():
        print(f"  ℹ️  results/目录不存在，将在首次运行时创建")
    else:
        all_ok &= check_path(results_dir, "结果目录")
    print()
    
    # 总结
    print("=" * 70)
    if all_ok:
        print("✅ 所有数据路径验证通过！")
        print()
        print("🎯 下一步测试:")
        print("  1. 测试Step 1: python step1_detect_outliers.py")
        print("  2. 测试Step 6: python step6_energy_analysis_v2.py")
        print("  3. 测试Step 7: python step7_lindemann_analysis.py")
    else:
        print("❌ 部分路径验证失败，请检查数据迁移")
        sys.exit(1)
    print("=" * 70)

if __name__ == "__main__":
    main()
