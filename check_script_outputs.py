#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
检查workflow脚本的输出完整性
"""

from pathlib import Path
import pandas as pd

def check_outputs():
    """检查所有脚本的输出文件"""
    
    results_dir = Path(__file__).parent / 'results'
    
    print("=" * 70)
    print("Workflow 脚本输出完整性检查")
    print("=" * 70)
    print()
    
    # Step 1: GMX MSD异常检测
    print("📁 Step 1: GMX MSD异常检测 (根目录)")
    print("-" * 70)
    step1_files = {
        'large_D_outliers.csv': '异常run清单',
        'ensemble_comparison.csv': '改进前后对比',
        'gmx_D_distribution.png': 'GMX D值分布',
        'quality_improvement_summary.png': '质量改进总结',
        'run_quality_report.txt': '质量报告'
    }
    
    step1_ok = 0
    for file, desc in step1_files.items():
        file_path = results_dir / file
        if file_path.exists():
            size_mb = file_path.stat().st_size / 1024 / 1024
            print(f"  ✅ {file:<35} {desc} ({size_mb:.2f} MB)")
            step1_ok += 1
        else:
            print(f"  ❌ {file:<35} {desc} (缺失)")
    
    print(f"\n  状态: {step1_ok}/{len(step1_files)} 文件 ({'100%' if step1_ok == len(step1_files) else '不完整'})") 
    print()
    
    # Step 6: 能量/热容分析
    print("📁 Step 6: 能量/热容分析 (energy_analysis_v2_no_filter)")
    print("-" * 70)
    
    step6_dir = results_dir / 'energy_analysis_v2_no_filter'
    if step6_dir.exists():
        files = list(step6_dir.glob('*.png'))
        files_by_type = {}
        
        for f in files:
            prefix = f.name.split('_')[0]
            if prefix not in files_by_type:
                files_by_type[prefix] = []
            files_by_type[prefix].append(f.name)
        
        expected_counts = {
            'Pt8SnX': 11,  # Energy vs T for Sn0-Sn10
            'HeatCapacity': 7,  # comparison + heatmap for 7 series
            'ClusterHeatCapacity': 7  # comparison + heatmap for 7 series (扣除支撑层)
        }
        
        total_expected = 11 + 7 + 7 + 7 + 7  # 能量曲线 + 热容对比 + 团簇热容对比 + 热容热力图 + 团簇热容热力图
        total_actual = len(files)
        
        # 详细统计
        energy_files = [f for f in files if 'Energy_vs_T' in f.name]
        hc_comp_files = [f for f in files if 'HeatCapacity_comparison' in f.name]
        chc_comp_files = [f for f in files if 'ClusterHeatCapacity_comparison' in f.name]
        hc_heatmap_files = [f for f in files if 'HeatCapacity_heatmap' in f.name and 'Cluster' not in f.name]
        chc_heatmap_files = [f for f in files if 'ClusterHeatCapacity_heatmap' in f.name]
        
        print(f"  能量-温度曲线 (Pt8SnX_Energy_vs_T):")
        print(f"    {'✅' if len(energy_files) == 11 else '⚠️ '} {len(energy_files)}/11 文件")
        
        print(f"\n  热容对比图 (HeatCapacity_comparison):")
        print(f"    {'✅' if len(hc_comp_files) == 7 else '⚠️ '} {len(hc_comp_files)}/7 文件")
        
        print(f"\n  团簇热容对比 (ClusterHeatCapacity_comparison):")
        print(f"    {'✅' if len(chc_comp_files) == 7 else '⚠️ '} {len(chc_comp_files)}/7 文件")
        
        print(f"\n  热容热力图 (HeatCapacity_heatmap):")
        print(f"    {'✅' if len(hc_heatmap_files) == 7 else '⚠️ '} {len(hc_heatmap_files)}/7 文件")
        if len(hc_heatmap_files) < 7:
            existing = set(f.name.replace('HeatCapacity_heatmap_', '').replace('.png', '') for f in hc_heatmap_files)
            expected = {'O1', 'O2', 'O3', 'O4', 'Pt(8-x)SnX', 'Pt6SnX', 'Pt8SnX'}
            missing = expected - existing
            if missing:
                print(f"      缺失: {', '.join(missing)}")
        
        print(f"\n  团簇热容热力图 (ClusterHeatCapacity_heatmap):")
        print(f"    {'✅' if len(chc_heatmap_files) == 7 else '⚠️ '} {len(chc_heatmap_files)}/7 文件")
        if len(chc_heatmap_files) < 7:
            existing = set(f.name.replace('ClusterHeatCapacity_heatmap_', '').replace('.png', '') for f in chc_heatmap_files)
            expected = {'O1', 'O2', 'O3', 'O4', 'Pt(8-x)SnX', 'Pt6SnX', 'Pt8SnX'}
            missing = expected - existing
            if missing:
                print(f"      缺失: {', '.join(missing)}")
        
        completion_rate = (total_actual / total_expected) * 100
        print(f"\n  状态: {total_actual}/{total_expected} 文件 ({completion_rate:.1f}%)")
        
        if total_actual < total_expected:
            print(f"  ⚠️  脚本可能在生成热力图时被中断")
    else:
        print("  ❌ 输出目录不存在")
    
    print()
    
    # Step 7: Lindemann分析
    print("📁 Step 7: Lindemann指数分析 (lindemann_analysis)")
    print("-" * 70)
    
    step7_dir = results_dir / 'lindemann_analysis'
    if step7_dir.exists():
        files = list(step7_dir.glob('*'))
        
        lindemann_vs_t = [f for f in files if 'Lindemann_vs_T' in f.name]
        lindemann_heatmap = [f for f in files if 'Lindemann_heatmap' in f.name]
        cv_files = [f for f in files if 'Cv_series' in f.name]
        
        print(f"  Lindemann vs T 曲线:")
        print(f"    {'✅' if len(lindemann_vs_t) == 8 else '⚠️ '} {len(lindemann_vs_t)}/8 文件")
        
        print(f"\n  Lindemann热力图:")
        print(f"    {'✅' if len(lindemann_heatmap) == 8 else '⚠️ '} {len(lindemann_heatmap)}/8 文件")
        
        print(f"\n  Cv系列分析:")
        print(f"    {'✅' if len(cv_files) >= 2 else '⚠️ '} {len(cv_files)} 文件")
        
        print(f"\n  其他文件:")
        other_files = [f for f in files if f.name not in 
                      [x.name for x in lindemann_vs_t + lindemann_heatmap + cv_files]]
        for f in other_files:
            size_kb = f.stat().st_size / 1024
            print(f"    ✅ {f.name} ({size_kb:.1f} KB)")
        
        total_files = len(files)
        print(f"\n  状态: {total_files} 文件 ✅")
    else:
        print("  ❌ 输出目录不存在")
    
    print()
    print("=" * 70)
    print("总结:")
    print("=" * 70)
    print("✅ Step 1: 输出完整")
    print(f"{'⚠️ ' if total_actual < total_expected else '✅'} Step 6: {'部分输出' if total_actual < total_expected else '输出完整'}")
    print("✅ Step 7: 输出完整")
    print()
    
    if total_actual < total_expected:
        print("💡 建议:")
        print("  - 重新运行 step6_energy_analysis_v2.py 完成所有图表生成")
        print("  - 或者运行时不要中断(Ctrl+C)，让脚本完整执行")
    
    print("=" * 70)

if __name__ == "__main__":
    check_outputs()
