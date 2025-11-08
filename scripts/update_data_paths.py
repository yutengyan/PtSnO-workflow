#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
批量更新workflow脚本中的数据路径
将原始路径替换为workflow/data/下的相对路径
"""

from pathlib import Path
import re

# 脚本路径映射
PATH_REPLACEMENTS = {
    # Step 1-5: GMX MSD路径
    "step1_detect_outliers.py": [
        {
            "old": "GMX_DATA_DIRS = [\n    'd:/OneDrive/py/Cv/lin/MSD_Analysis_Collection/test-unwrap-new/file/collected_gmx_msd',\n    'd:/OneDrive/py/Cv/lin/MSD_Analysis_Collection/test-unwrap-new/file/gmx_msd_results_20251015_184626_collected'\n]",
            "new": "GMX_DATA_DIRS = [\n    BASE_DIR / 'data' / 'gmx_msd' / 'collected_gmx_msd',\n    BASE_DIR / 'data' / 'gmx_msd' / 'gmx_msd_results_20251015_184626_collected'\n]"
        },
        {
            "old": "BASE_DIR = Path(__file__).parent.parent",
            "new": "BASE_DIR = Path(__file__).parent"
        }
    ],
    
    "step2_ensemble_analysis.py": [
        {
            "old": "GMX_DATA_DIRS = [\n    'd:/OneDrive/py/Cv/lin/MSD_Analysis_Collection/test-unwrap-new/file/collected_gmx_msd',\n    'd:/OneDrive/py/Cv/lin/MSD_Analysis_Collection/test-unwrap-new/file/gmx_msd_results_20251015_184626_collected'\n]",
            "new": "GMX_DATA_DIRS = [\n    BASE_DIR / 'data' / 'gmx_msd' / 'collected_gmx_msd',\n    BASE_DIR / 'data' / 'gmx_msd' / 'gmx_msd_results_20251015_184626_collected'\n]"
        },
        {
            "old": "BASE_DIR = Path(__file__).parent.parent",
            "new": "BASE_DIR = Path(__file__).parent"
        }
    ],
    
    # Step 6: LAMMPS能量路径
    "step6_energy_analysis_v2.py": [
        {
            "old": "BASE_DIR = Path(__file__).parent.parent",
            "new": "BASE_DIR = Path(__file__).parent"
        },
        {
            "old": "ENERGY_MASTER = BASE_DIR / 'files' / 'lammps_energy_analysis' / 'energy_master_20251016_121110.csv'",
            "new": "ENERGY_MASTER = BASE_DIR / 'data' / 'lammps_energy' / 'energy_master_20251016_121110.csv'"
        }
    ],
    
    "step5.9calculate_support_heat_capacity.py": [
        {
            "old": "BASE_DIR = Path(__file__).parent.parent",
            "new": "BASE_DIR = Path(__file__).parent"
        },
        {
            "old": "ENERGY_FILE = BASE_DIR / 'files' / 'lammps_energy_analysis' / 'sup' / 'energy_master_20251021_151520.csv'",
            "new": "ENERGY_FILE = BASE_DIR / 'data' / 'lammps_energy' / 'sup' / 'energy_master_20251021_151520.csv'"
        }
    ],
    
    # Step 7: Lindemann路径
    "step7_lindemann_analysis.py": [
        {
            "old": "BASE_DIR = Path(r'd:\\OneDrive\\py\\Cv\\lin\\MSD_Analysis_Collection\\v3_simplified_workflow')",
            "new": "BASE_DIR = Path(__file__).parent"
        },
        {
            "old": "DATA_DIR = BASE_DIR / 'files' / 'takeit'",
            "new": "DATA_DIR = BASE_DIR / 'data' / 'lindemann'"
        }
    ],
    
    "step7_4_multi_system_heat_capacity.py": [
        {
            "old": "BASE_DIR = Path(__file__).parent.parent",
            "new": "BASE_DIR = Path(__file__).parent"
        },
        {
            "old": "ENERGY_MASTER = BASE_DIR / 'files' / 'lammps_energy_analysis' / 'energy_master_20251016_121110.csv'",
            "new": "ENERGY_MASTER = BASE_DIR / 'data' / 'lammps_energy' / 'energy_master_20251016_121110.csv'"
        },
        {
            "old": "DATA_DIR = BASE_DIR / 'files' / 'takeit'",
            "new": "DATA_DIR = BASE_DIR / 'data' / 'lindemann'"
        }
    ],
    
    # Step 7.5: 配位数路径
    "step7-5-unified_multi_temp_v626_analysis.py": [
        {
            "old": 'base_path = r"D:\\OneDrive\\py\\Cv\\lin\\MSD_Analysis_Collection\\v3_simplified_workflow\\files\\q6_cn\\v626\\coordination_time_series_results_sample_20251026_200908"',
            "new": "base_path = Path(__file__).parent / 'data' / 'coordination' / 'coordination_time_series_results_sample_20251026_200908'"
        }
    ],
}

def update_script(script_name, replacements):
    """更新单个脚本的路径"""
    script_path = Path(__file__).parent / script_name
    
    if not script_path.exists():
        print(f"❌ 未找到: {script_name}")
        return False
    
    try:
        # 读取脚本内容
        content = script_path.read_text(encoding='utf-8')
        original_content = content
        
        # 应用所有替换
        for replacement in replacements:
            old_text = replacement["old"]
            new_text = replacement["new"]
            
            if old_text in content:
                content = content.replace(old_text, new_text)
                print(f"  ✅ 替换: {old_text[:50]}... → {new_text[:50]}...")
            else:
                print(f"  ⚠️  未找到: {old_text[:50]}...")
        
        # 如果有修改，写回文件
        if content != original_content:
            script_path.write_text(content, encoding='utf-8')
            print(f"✅ 更新成功: {script_name}\n")
            return True
        else:
            print(f"⚠️  无需修改: {script_name}\n")
            return False
            
    except Exception as e:
        print(f"❌ 错误: {script_name} - {e}\n")
        return False

def main():
    """批量更新所有脚本"""
    print("=" * 70)
    print("开始批量更新workflow脚本的数据路径")
    print("=" * 70)
    print()
    
    success_count = 0
    total_count = len(PATH_REPLACEMENTS)
    
    for script_name, replacements in PATH_REPLACEMENTS.items():
        print(f"📝 处理: {script_name}")
        if update_script(script_name, replacements):
            success_count += 1
    
    print("=" * 70)
    print(f"完成! 成功更新 {success_count}/{total_count} 个脚本")
    print("=" * 70)
    print()
    print("📋 更新的脚本:")
    print("  Step 1-5: step1_detect_outliers.py, step2_ensemble_analysis.py")
    print("  Step 6:   step6_energy_analysis_v2.py, step5.9calculate_support_heat_capacity.py")
    print("  Step 7:   step7_lindemann_analysis.py, step7_4_multi_system_heat_capacity.py")
    print("  Step 7.5: step7-5-unified_multi_temp_v626_analysis.py")
    print()
    print("🎯 下一步:")
    print("  1. 检查step3, step4, step5的输出路径")
    print("  2. 检查step6.2, step6_3的输入路径")
    print("  3. 检查step7_4_2, step7-6系列的输入路径")
    print("  4. 运行测试: python step1_detect_outliers.py")

if __name__ == "__main__":
    main()
