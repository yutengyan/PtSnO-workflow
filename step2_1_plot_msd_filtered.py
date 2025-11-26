#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
绘制特定体系的MSD曲线 - 仅300K和900K (可选过滤)
====================================================
创建时间: 2025-10-16
最后更新: 2025-11-25

支持的数据集:
    1. 'pt8sn6' - Pt8Sn6 负载型数据
       路径: data/gmx_msd/unwrap/gmx_msd_results_20251118_152614
       输出: results/msd_curves_pt8sn6_loaded/
    
    2. 'air_86' - 气象数据 86 系统
       路径: data/gmx_msd/unwrap/air/gmx_msd_results_20251124_170114
       输出: results/msd_curves_air_86/

功能:
1. 支持多个数据集切换 (通过 DATASET 配置)
2. 仅绘制300K和900K两个温度
3. 可选择是否过滤D值错误的runs (通过ENABLE_FILTERING配置)
4. 统一Y轴坐标范围
5. 简洁的图表显示

配置方法:
---------
1. 切换数据集:
   修改 DATASET 变量:
   DATASET = 'pt8sn6'    # 使用 Pt8Sn6 负载型数据
   DATASET = 'air_86'    # 使用气象数据 86 系统

2. 启用/禁用过滤:
   ENABLE_FILTERING = True   # 过滤异常runs (默认)
   ENABLE_FILTERING = False  # 绘制所有runs

使用示例:
---------
# 绘制 Pt8Sn6 数据
DATASET = 'pt8sn6'
ENABLE_FILTERING = True
python step3_1_plot_msd_filtered.py

# 绘制气象数据 86 系统
DATASET = 'air_86'
ENABLE_FILTERING = True
python step3_1_plot_msd_filtered.py
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from datetime import datetime
import re  # 添加正则表达式支持
import warnings
from collections import defaultdict
warnings.filterwarnings('ignore')

# ===== 数据集选择配置 =====
# 可选值: 'pt8sn6' 或 'air_86'
DATASET = 'pt8sn6'  # 🔧 修改此处切换数据集

# 数据集配置字典
DATASET_CONFIGS = {
    'pt8sn6': {
        'name': 'Pt8Sn6负载型',
        'data_dir': 'data/gmx_msd/unwrap/gmx_msd_results_20251118_152614',
        'output_dir': 'results/msd_curves_pt8sn6_loaded',
        'system_pattern': r'^pt8sn6',
        'target_temps': ['300K', '900K'],
        'description': 'Pt8Sn6 负载型数据 (unwrap per-atom MSD)'
    },
    'air_86': {
        'name': '气象数据86',
        'data_dir': 'data/gmx_msd/unwrap/air/gmx_msd_results_20251124_170114',
        'output_dir': 'results/msd_curves_air_86',
        'system_pattern': r'^86$',
        'target_temps': ['300K', '900K'],
        'description': '气象数据 86 系统 (atmospheric conditions)'
    }
}

# ===== 全局配置 =====
BASE_DIR = Path(__file__).parent  # workflow目录
OUTLIERS_CSV = BASE_DIR / 'results' / 'large_D_outliers.csv'

# 根据选择的数据集设置配置
if DATASET not in DATASET_CONFIGS:
    raise ValueError(f"未知的数据集: {DATASET}. 可选值: {list(DATASET_CONFIGS.keys())}")

current_config = DATASET_CONFIGS[DATASET]
GMX_DATA_DIRS = [BASE_DIR / current_config['data_dir']]
OUTPUT_DIR = BASE_DIR / current_config['output_dir']
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

# ===== 专用配置 =====
TARGET_SYSTEM_PATTERN = current_config['system_pattern']  # 目标体系匹配模式
TARGET_TEMPS = current_config['target_temps']  # 目标温度
ENABLE_FILTERING = True  # 是否启用D值过滤 (True=过滤异常runs, False=绘制所有runs)

COLORS = {
    'Pt': '#E74C3C',
    'Sn': '#3498DB',
    'PtSn': '#2ECC71'
}

plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False


def load_large_D_outliers():
    """
    加载大D值异常run清单
    
    Returns:
    --------
    outliers : set
        异常文件路径集合
    """
    if not ENABLE_FILTERING:
        print(f"\n[!] 过滤功能已禁用 - 将绘制所有runs")
        return set()
    
    try:
        df_outliers = pd.read_csv(OUTLIERS_CSV)
        outlier_files = set(df_outliers['filepath'].values)
        print(f"\n[√] 已加载异常run清单: {len(outlier_files)} 个文件")
        print(f"    这些runs将被自动过滤掉")
        return outlier_files
    except FileNotFoundError:
        print(f"\n[!] 未找到异常文件清单: {OUTLIERS_CSV}")
        print(f"    将绘制所有runs")
        return set()


def build_file_index_for_target(outlier_files=None):
    """
    构建目标体系的文件索引 (带联动过滤)
    
    Parameters:
    -----------
    outlier_files : set, optional
        异常文件路径集合
    
    Returns:
    --------
    file_index : dict
        {(composition, temperature, element): [file_path1, file_path2, ...]}
    stats : dict
        筛选统计
    compositions_found : list
        找到的所有匹配composition列表
    
    联动过滤逻辑:
        如果同一个模拟的任何元素被标记为异常，则该模拟的所有元素都被过滤
    """
    if outlier_files is None:
        outlier_files = set()
    
    print(f"\n[*] 正在为匹配 '{TARGET_SYSTEM_PATTERN}' 的体系构建文件索引...")
    print(f"    目标温度: {', '.join(TARGET_TEMPS)}")
    
    # ========== 第一遍扫描: 收集所有文件并识别异常模拟 ==========
    all_files = []  # [(xvg_file, composition, temperature, element, sim_id), ...]
    bad_simulations = set()  # 有异常元素的模拟ID集合
    
    for gmx_dir in GMX_DATA_DIRS:
        if not gmx_dir.exists():
            continue
            
        print(f"  扫描目录: {gmx_dir.name}...")
        
        for xvg_file in gmx_dir.rglob("*_msd_*.xvg"):
            try:
                parts = xvg_file.parts
                filename = xvg_file.stem
                if '_msd_' not in filename:
                    continue
                element = filename.split('_msd_')[-1]
                
                # 提取 sim_id (如 T300.r1.gpu0)
                sim_id = filename.split('_msd_')[0]  # T300.r1.gpu0
                
                temperature = None
                composition = None
                for i in range(len(parts)-1, 0, -1):
                    if parts[i].endswith('K'):
                        temperature = parts[i]
                        composition = parts[i-1]
                        break
                
                if not temperature or not composition:
                    continue
                
                if not re.match(TARGET_SYSTEM_PATTERN, composition, re.IGNORECASE):
                    continue
                if temperature not in TARGET_TEMPS:
                    continue
                
                # 完整的模拟标识 (composition + temperature + sim_id)
                full_sim_id = f"{composition}_{temperature}_{sim_id}"
                
                all_files.append((xvg_file, composition, temperature, element, full_sim_id))
                
                # 检查是否是异常run - 如果是，标记整个模拟
                if str(xvg_file) in outlier_files:
                    bad_simulations.add(full_sim_id)
                    
            except Exception as e:
                continue
    
    print(f"  [INFO] 识别到 {len(bad_simulations)} 个异常模拟 (联动过滤)")
    
    # ========== 第二遍: 按联动规则过滤 ==========
    file_index = defaultdict(list)
    compositions_found = set()
    stats = {
        'total_found': 0,
        'kept': 0,
        'filtered': 0,
        'by_temp': {},
        'by_element': {}
    }
    
    for temp in TARGET_TEMPS:
        stats['by_temp'][temp] = {'total': 0, 'kept': 0, 'filtered': 0}
    for elem in ['Pt', 'Sn', 'PtSn']:
        stats['by_element'][elem] = {'total': 0, 'kept': 0, 'filtered': 0}
    
    for xvg_file, composition, temperature, element, full_sim_id in all_files:
        compositions_found.add(composition)
        stats['total_found'] += 1
        stats['by_temp'][temperature]['total'] += 1
        if element in stats['by_element']:
            stats['by_element'][element]['total'] += 1
        
        # 联动过滤: 如果该模拟有任何异常元素，则全部过滤
        if full_sim_id in bad_simulations:
            stats['filtered'] += 1
            stats['by_temp'][temperature]['filtered'] += 1
            if element in stats['by_element']:
                stats['by_element'][element]['filtered'] += 1
        else:
            key = (composition, temperature, element)
            file_index[key].append(xvg_file)
            stats['kept'] += 1
            stats['by_temp'][temperature]['kept'] += 1
            if element in stats['by_element']:
                stats['by_element'][element]['kept'] += 1
    
    compositions_found = sorted(compositions_found)
    
    print(f"\n  [√] 索引构建完成 (联动过滤已启用):")
    print(f"      找到的composition: {', '.join(compositions_found) if compositions_found else '无'}")
    print(f"      总文件数: {stats['total_found']}")
    print(f"      保留: {stats['kept']} 个有效文件")
    print(f"      过滤: {stats['filtered']} 个异常文件 (联动)")
    
    for temp in TARGET_TEMPS:
        temp_stats = stats['by_temp'][temp]
        print(f"      {temp}: {temp_stats['kept']}/{temp_stats['total']} 个有效文件 "
              f"(过滤 {temp_stats['filtered']} 个)")
    
    # 按元素统计 (验证联动一致性)
    print(f"      按元素统计:")
    for elem in ['Pt', 'Sn', 'PtSn']:
        elem_stats = stats['by_element'].get(elem, {'kept': 0, 'total': 0})
        print(f"        {elem}: {elem_stats['kept']}/{elem_stats['total']} 个有效文件")
    
    return file_index, stats, compositions_found


def read_gmx_msd_xvg(filepath):
    """读取GMX MSD .xvg文件"""
    time_data = []
    msd_data = []
    
    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#') or line.startswith('@'):
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        t = float(parts[0])
                        msd_nm2 = float(parts[1])
                        msd_a2 = msd_nm2 * 100  # nm^2 -> A^2
                        time_data.append(t)
                        msd_data.append(msd_a2)
                    except ValueError:
                        continue
    except:
        return None, None
    
    if len(time_data) == 0:
        return None, None
    
    return np.array(time_data), np.array(msd_data)


def plot_pt8sn6_300k_900k(file_index, stats, compositions_found):
    """
    绘制Pt8Sn6的300K和900K MSD曲线
    
    Parameters:
    -----------
    file_index : dict
        文件索引
    stats : dict
        统计信息
    compositions_found : list
        找到的composition列表
    """
    comp_str = ', '.join(compositions_found) if len(compositions_found) > 1 else compositions_found[0]
    print(f"\n{'='*80}")
    print(f"绘制: {comp_str} - 300K & 900K (过滤版)")
    print(f"{'='*80}")
    
    # 创建2列布局 (300K | 900K)
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # 第1步: 加载所有数据并找全局最大MSD值
    print("\n[*] 加载MSD数据...")
    
    msd_cache = {}  # {(temp, element): [(time, msd), ...]}
    global_max_msd = 0
    
    for temp in TARGET_TEMPS:
        for element in ['Pt', 'Sn', 'PtSn']:
            msd_list = []
            
            # 遍历所有找到的composition
            for comp in compositions_found:
                key = (comp, temp, element)
                files = file_index.get(key, [])
                
                for filepath in files:
                    time, msd = read_gmx_msd_xvg(filepath)
                    if time is not None:
                        msd_list.append((time, msd))
                        
                        # 更新全局最大值
                        max_val = np.max(msd)
                        if max_val > global_max_msd:
                            global_max_msd = max_val
            
            if msd_list:
                msd_cache[(temp, element)] = msd_list
    
    total_curves = sum(len(v) for v in msd_cache.values())
    print(f"  [√] 已加载 {total_curves} 条有效MSD曲线")
    print(f"  [√] 全局最大MSD值: {global_max_msd:.2f} Å²")
    
    # 设置统一的Y轴上限 (留10%余量)
    unified_ylim = global_max_msd * 1.1
    
    # 第2步: 绘制每个温度
    print("\n[*] 绘制中...")
    
    for idx, temp in enumerate(TARGET_TEMPS):
        ax = axes[idx]
        has_data = False
        
        temp_stats = stats['by_temp'][temp]
        
        # 绘制三种元素
        for element in ['Pt', 'Sn', 'PtSn']:
            msd_list = msd_cache.get((temp, element), [])
            
            if not msd_list:
                continue
            
            has_data = True
            color = COLORS.get(element, '#95a5a6')
            
            # 对齐时间轴并绘制所有runs
            min_len = min(len(msd) for _, msd in msd_list)
            for time, msd in msd_list:
                ax.plot(time[:min_len], msd[:min_len], 
                       color=color, alpha=0.3, linewidth=1)
            
            # 计算平均曲线
            msd_aligned = np.array([msd[:min_len] for _, msd in msd_list])
            time_common = msd_list[0][0][:min_len]
            msd_mean = np.mean(msd_aligned, axis=0)
            
            # 绘制平均曲线 (粗线)
            ax.plot(time_common, msd_mean, 
                   color=color, linewidth=3, alpha=0.9,
                   label=f'{element} (n={len(msd_list)})')
        
        if has_data:
            # 统一Y轴范围
            ax.set_ylim(0, unified_ylim)
            
            ax.set_xlabel('Time (ps)', fontsize=12, fontweight='bold')
            ax.set_ylabel(r'MSD ($\AA^2$)', fontsize=12, fontweight='bold')
            
            # 标题：根据是否过滤显示不同信息
            if ENABLE_FILTERING and temp_stats['filtered'] > 0:
                title_text = (f'{comp_str} @ {temp}\n'
                             f'有效runs: {temp_stats["kept"]}')
            else:
                title_text = f'{comp_str} @ {temp}\nRuns: {temp_stats["kept"]}'
            
            ax.set_title(title_text, fontsize=13, fontweight='bold')
            ax.legend(fontsize=10, loc='upper left', framealpha=0.9)
            ax.grid(True, alpha=0.3, linestyle=':', linewidth=1)
        else:
            ax.text(0.5, 0.5, f'无数据\n{temp}',
                   ha='center', va='center',
                   fontsize=14, fontweight='bold',
                   transform=ax.transAxes)
            ax.axis('off')
    
    # 总标题
    if ENABLE_FILTERING and stats['filtered'] > 0:
        title_suffix = f'(过滤版 - 已过滤{stats["filtered"]}个异常runs)'
    else:
        title_suffix = '(所有runs)'
    
    fig.suptitle(f'{comp_str.upper()} MSD曲线对比 - 300K vs 900K {title_suffix}',
                 fontsize=15, fontweight='bold', y=0.98)
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    # 保存图片
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    filename_base = compositions_found[0].replace('-', '_') if compositions_found else 'pt8sn6'
    filter_tag = 'filtered' if ENABLE_FILTERING else 'all'
    output_file = OUTPUT_DIR / f'{filename_base}_300K_900K_{filter_tag}_{timestamp}.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n[√] 图片已保存: {output_file}")
    
    plt.show()
    plt.close()


def main():
    """主函数"""
    print("\n" + "="*80)
    print(f"MSD曲线绘制 - {current_config['name']} - 300K & 900K")
    print("="*80)
    print(f"\n数据集: {DATASET}")
    print(f"  描述: {current_config['description']}")
    print(f"  数据路径: {current_config['data_dir']}")
    print(f"\n配置:")
    print(f"  目标体系模式: {TARGET_SYSTEM_PATTERN}")
    print(f"  目标温度: {', '.join(TARGET_TEMPS)}")
    print(f"  D值过滤: {'启用' if ENABLE_FILTERING else '禁用'}")
    print(f"  输出目录: {OUTPUT_DIR}")
    
    # 1. 加载异常run清单
    outlier_files = load_large_D_outliers()
    
    # 2. 构建文件索引 (仅目标体系和温度)
    file_index, stats, compositions_found = build_file_index_for_target(outlier_files)
    
    if stats['kept'] == 0:
        print(f"\n[X] 错误: 没有找到有效的MSD数据文件!")
        print(f"    请检查:")
        print(f"    1. 数据目录是否正确")
        print(f"    2. 匹配 '{TARGET_SYSTEM_PATTERN}' 的体系是否存在")
        print(f"    3. 300K和900K温度点是否存在")
        return
    
    # 3. 绘制
    plot_pt8sn6_300k_900k(file_index, stats, compositions_found)
    
    print("\n" + "="*80)
    print("完成!")
    print("="*80)


if __name__ == '__main__':
    main()
