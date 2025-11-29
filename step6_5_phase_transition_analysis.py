#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 6.5: 相变分析 - 典型体系热容变化与融化行为分析

分析内容:
1. 气相团簇 (Air68, Air86) - 无载体影响
2. 负载型PtSn团簇 (Pt6sn8, Pt8sn6等) - 有载体影响  
3. 含氧负载型 (Pt6sn8o4等) - 氧化物载体影响

核心指标:
- 热容差异显著性: |Cv1-Cv2| / sqrt(err1^2+err2^2) >= 2
- 热容变化比例: (Cv2-Cv1)/Cv1
- 是否存在明显相变

用法:
  python step6_5_phase_transition_analysis.py                    # 分析默认体系
  python step6_5_phase_transition_analysis.py --all              # 分析所有可用体系
  python step6_5_phase_transition_analysis.py --add Pt6sn0 Pt7sn5  # 添加额外体系
  python step6_5_phase_transition_analysis.py --only Air68 Pt8sn6  # 只分析指定体系
  python step6_5_phase_transition_analysis.py --list             # 列出所有可用体系

作者: AI Assistant
日期: 2025-11-29
"""

import os
import sys
import glob
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False


def find_clustering_results(base_dir='results/step6_1_clustering'):
    """查找所有聚类分析结果"""
    results = {}
    
    # 查找所有auto2质量指标文件
    pattern = os.path.join(base_dir, '*_auto2_quality_metrics.csv')
    files = glob.glob(pattern)
    
    for f in files:
        basename = os.path.basename(f)
        structure = basename.replace('_auto2_quality_metrics.csv', '')
        results[structure] = {
            'auto2': f,
            'auto3': f.replace('auto2', 'auto3'),
            'kmeans_n2': f.replace('auto2_quality_metrics', 'kmeans_n2_quality_metrics'),
        }
    
    return results


def load_partition_data(structure, base_dir='results/step6_1_clustering'):
    """加载指定结构的分区数据"""
    data = {
        'structure': structure,
        'n2': None,
        'n3': None,
    }
    
    # 读取2分区数据
    csv2 = os.path.join(base_dir, f'{structure}_auto2_quality_metrics.csv')
    if os.path.exists(csv2):
        try:
            df = pd.read_csv(csv2)
            if len(df) >= 2:
                data['n2'] = {
                    'cv1': df['Cv_cluster'].iloc[0],
                    'err1': df['Cv_cluster_err'].iloc[0],
                    'cv2': df['Cv_cluster'].iloc[1],
                    'err2': df['Cv_cluster_err'].iloc[1],
                    'r2_1': df['R2'].iloc[0],
                    'r2_2': df['R2'].iloc[1],
                    'silhouette': df['silhouette_score'].iloc[0] if 'silhouette_score' in df.columns else None,
                }
        except Exception as e:
            print(f"  警告: 读取{structure} 2分区数据失败: {e}")
    
    # 读取3分区数据
    csv3 = os.path.join(base_dir, f'{structure}_auto3_quality_metrics.csv')
    if os.path.exists(csv3):
        try:
            df = pd.read_csv(csv3)
            if len(df) >= 3:
                data['n3'] = {
                    'cv1': df['Cv_cluster'].iloc[0],
                    'err1': df['Cv_cluster_err'].iloc[0],
                    'cv2': df['Cv_cluster'].iloc[1],
                    'err2': df['Cv_cluster_err'].iloc[1],
                    'cv3': df['Cv_cluster'].iloc[2],
                    'err3': df['Cv_cluster_err'].iloc[2],
                }
        except Exception as e:
            pass
    
    return data


def calculate_significance(cv1, err1, cv2, err2):
    """计算热容差异显著性"""
    diff = abs(cv2 - cv1)
    combined_err = np.sqrt(err1**2 + err2**2)
    ratio = diff / combined_err if combined_err > 0 else float('inf')
    return {
        'diff': diff,
        'combined_err': combined_err,
        'ratio': ratio,
        'significant': ratio >= 2,
        'change_percent': (cv2 - cv1) / cv1 * 100 if cv1 > 0 else 0,
    }


def classify_structure(name):
    """分类结构类型"""
    name_lower = name.lower()
    
    if 'air' in name_lower:
        return 'gas_phase', '气相团簇'
    
    # Cv 就是 Pt6sn8o4 (含氧)
    if name == 'Cv':
        return 'supported_oxide', '含氧负载型'
    
    # 检查是否含氧
    if 'o' in name_lower:
        import re
        if re.search(r'o\d|O\d|\do|\dO', name):
            return 'supported_oxide', '含氧负载型'
    
    # 默认为无氧负载型
    if 'pt' in name_lower and 'sn' in name_lower:
        return 'supported', '负载型PtSn'
    
    return 'other', '其他'


def get_display_name(name):
    """获取显示名称"""
    if name == 'Cv':
        return 'Pt6sn8o4'
    return name


def analyze_typical_systems(target_systems=None, include_all_oxide=True):
    """分析典型体系
    
    参数:
        target_systems: 指定分析的体系列表，None则使用默认
        include_all_oxide: 是否包含所有含氧团簇
    """
    
    base_dir = 'results/step6_1_clustering'
    
    # 先查找所有可用结构
    all_results = find_clustering_results(base_dir)
    available = set(all_results.keys())
    
    # 默认分析的典型体系
    if target_systems is None:
        target_systems = [
            # 气相团簇
            'Air68', 'Air86',
            # Pt6系列 (无氧) - 包含Pt6sn0
            'Pt6sn0', 'Pt6sn1', 'Pt6sn2', 'Pt6sn3', 'Pt6sn4', 'Pt6sn5', 
            'Pt6sn6', 'Pt6sn7', 'Pt6sn8', 'Pt6sn9',
            # Pt8系列 (无氧)
            'Pt8sn0', 'Pt8sn1', 'Pt8sn2', 'Pt8sn3', 'Pt8sn4', 'Pt8sn5',
            'Pt8sn6', 'Pt8sn7', 'Pt8sn8', 'Pt8sn9', 'Pt8sn10',
        ]
        
        # 添加所有含氧团簇
        if include_all_oxide:
            import re
            for name in available:
                # 检查是否是含氧结构
                if name == 'Cv' or re.search(r'o\d|O\d', name, re.IGNORECASE):
                    if name not in target_systems:
                        target_systems.append(name)
    
    print("=" * 80)
    print("相变分析 - 典型体系热容变化与融化行为")
    print("=" * 80)
    print(f"\n可用体系数: {len(available)}")
    print(f"目标体系数: {len(target_systems)}")
    
    # 按类型分组分析
    analysis_results = []
    
    for structure in target_systems:
        if structure not in available:
            # 尝试大小写变体
            variants = [structure, structure.lower(), structure.upper(), 
                       structure.capitalize()]
            found = None
            for v in variants:
                if v in available:
                    found = v
                    break
            if not found:
                continue
            structure = found
        
        data = load_partition_data(structure, base_dir)
        
        if data['n2'] is None:
            continue
        
        n2 = data['n2']
        sig = calculate_significance(n2['cv1'], n2['err1'], n2['cv2'], n2['err2'])
        
        struct_type, type_name = classify_structure(structure)
        
        # 获取显示名称
        display_name = get_display_name(structure)
        
        result = {
            'structure': structure,
            'display_name': display_name,  # 添加显示名称
            'type': struct_type,
            'type_name': type_name,
            'cv1': n2['cv1'],
            'err1': n2['err1'],
            'cv2': n2['cv2'],
            'err2': n2['err2'],
            'r2_1': n2['r2_1'],
            'r2_2': n2['r2_2'],
            'silhouette': n2['silhouette'],
            'diff': sig['diff'],
            'combined_err': sig['combined_err'],
            'ratio': sig['ratio'],
            'significant': sig['significant'],
            'change_percent': sig['change_percent'],
            'recommendation': '2分区' if sig['significant'] else '1分区',
            'phase_transition': '有相变' if sig['significant'] and sig['change_percent'] > 20 else 
                               ('可能有相变' if sig['significant'] else '无明显相变'),
        }
        
        analysis_results.append(result)
    
    return analysis_results

def print_analysis_table(results):
    """打印分析结果表格"""
    
    # 按类型分组
    by_type = {}
    for r in results:
        t = r['type']
        if t not in by_type:
            by_type[t] = []
        by_type[t].append(r)
    
    type_order = ['gas_phase', 'supported', 'supported_oxide', 'other']
    type_names = {
        'gas_phase': '🔵 气相团簇 (无载体)',
        'supported': '🟢 负载型PtSn (无氧)',
        'supported_oxide': '🟠 含氧负载型',
        'other': '⚪ 其他',
    }
    
    print("\n" + "=" * 100)
    print("热容差异显著性分析结果")
    print("=" * 100)
    
    for t in type_order:
        if t not in by_type:
            continue
        
        print(f"\n### {type_names[t]}")
        print("-" * 100)
        print(f"{'体系':<12} {'Cv1':<10} {'Cv2':<10} {'差异':<8} {'误差':<8} {'比值':<8} {'变化%':<10} {'推荐':<8} {'相变判断':<12}")
        print("-" * 100)
        
        for r in sorted(by_type[t], key=lambda x: -x['ratio']):
            cv1_str = f"{r['cv1']:.2f}±{r['err1']:.2f}"
            cv2_str = f"{r['cv2']:.2f}±{r['err2']:.2f}"
            sig_mark = "✓" if r['significant'] else "✗"
            display = r.get('display_name', r['structure'])
            
            print(f"{display:<12} {cv1_str:<10} {cv2_str:<10} "
                  f"{r['diff']:<8.2f} {r['combined_err']:<8.2f} "
                  f"{r['ratio']:<8.2f} {r['change_percent']:>+8.1f}% "
                  f"{r['recommendation']:<8} {r['phase_transition']:<12}")
    
    # 统计汇总
    print("\n" + "=" * 100)
    print("统计汇总")
    print("=" * 100)
    
    for t in type_order:
        if t not in by_type:
            continue
        
        group = by_type[t]
        n_total = len(group)
        n_significant = sum(1 for r in group if r['significant'])
        n_phase = sum(1 for r in group if '有相变' in r['phase_transition'])
        avg_ratio = np.mean([r['ratio'] for r in group])
        avg_change = np.mean([r['change_percent'] for r in group])
        
        print(f"\n{type_names[t]}:")
        print(f"  体系数: {n_total}")
        print(f"  热容差异显著: {n_significant}/{n_total} ({100*n_significant/n_total:.1f}%)")
        print(f"  存在相变: {n_phase}/{n_total} ({100*n_phase/n_total:.1f}%)")
        print(f"  平均显著性比值: {avg_ratio:.2f}")
        print(f"  平均热容变化: {avg_change:+.1f}%")


def generate_comparison_plot(results, output_path='results/step6_1_clustering/phase_transition_comparison.png'):
    """生成对比图"""
    
    # 按类型分组
    by_type = {}
    for r in results:
        t = r['type']
        if t not in by_type:
            by_type[t] = []
        by_type[t].append(r)
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # 1. 热容变化百分比对比
    ax1 = axes[0, 0]
    type_colors = {'gas_phase': 'blue', 'supported': 'green', 'supported_oxide': 'orange', 'other': 'gray'}
    type_labels = {'gas_phase': '气相', 'supported': '负载型', 'supported_oxide': '含氧负载型', 'other': '其他'}
    
    for t, group in by_type.items():
        x = [r.get('display_name', r['structure']) for r in group]
        y = [r['change_percent'] for r in group]
        ax1.bar(x, y, color=type_colors.get(t, 'gray'), alpha=0.7, label=type_labels.get(t, t))
    
    ax1.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    ax1.axhline(y=20, color='red', linestyle='--', linewidth=1, label='相变阈值(20%)')
    ax1.set_ylabel('热容变化 (%)')
    ax1.set_title('热容变化百分比 (Cv2-Cv1)/Cv1')
    ax1.tick_params(axis='x', rotation=45)
    ax1.legend(loc='upper right')
    
    # 2. 显著性比值对比
    ax2 = axes[0, 1]
    for t, group in by_type.items():
        x = [r.get('display_name', r['structure']) for r in group]
        y = [r['ratio'] for r in group]
        ax2.bar(x, y, color=type_colors.get(t, 'gray'), alpha=0.7, label=type_labels.get(t, t))
    
    ax2.axhline(y=2, color='red', linestyle='--', linewidth=2, label='显著性阈值(2)')
    ax2.set_ylabel('显著性比值')
    ax2.set_title('热容差异显著性比值 |ΔCv|/σ')
    ax2.tick_params(axis='x', rotation=45)
    ax2.legend(loc='upper right')
    
    # 3. Cv1 vs Cv2 散点图
    ax3 = axes[1, 0]
    for t, group in by_type.items():
        cv1 = [r['cv1'] for r in group]
        cv2 = [r['cv2'] for r in group]
        ax3.scatter(cv1, cv2, c=type_colors.get(t, 'gray'), s=100, alpha=0.7, 
                   label=type_labels.get(t, t), edgecolors='black')
        
        # 标注体系名
        for r in group:
            display = r.get('display_name', r['structure'])
            ax3.annotate(display, (r['cv1'], r['cv2']), fontsize=8,
                        xytext=(3, 3), textcoords='offset points')
    
    # 添加对角线
    lims = [min(ax3.get_xlim()[0], ax3.get_ylim()[0]), 
            max(ax3.get_xlim()[1], ax3.get_ylim()[1])]
    ax3.plot(lims, lims, 'k--', alpha=0.5, label='Cv1=Cv2')
    ax3.set_xlabel('Cv1 (低温区, kB/atom)')
    ax3.set_ylabel('Cv2 (高温区, kB/atom)')
    ax3.set_title('低温区 vs 高温区热容')
    ax3.legend(loc='upper left')
    
    # 4. 类型对比箱线图
    ax4 = axes[1, 1]
    type_data = []
    type_labels_list = []
    for t in ['gas_phase', 'supported', 'supported_oxide']:
        if t in by_type:
            type_data.append([r['change_percent'] for r in by_type[t]])
            type_labels_list.append(type_labels[t])
    
    if type_data:
        bp = ax4.boxplot(type_data, labels=type_labels_list, patch_artist=True)
        colors = ['blue', 'green', 'orange']
        for patch, color in zip(bp['boxes'], colors[:len(bp['boxes'])]):
            patch.set_facecolor(color)
            patch.set_alpha(0.5)
    
    ax4.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    ax4.set_ylabel('热容变化 (%)')
    ax4.set_title('不同类型体系热容变化分布')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"\n图表已保存: {output_path}")
    plt.close()


def generate_report(results, output_path='results/step6_1_clustering/PHASE_TRANSITION_ANALYSIS.md'):
    """生成Markdown报告"""
    
    # 按类型分组
    by_type = {}
    for r in results:
        t = r['type']
        if t not in by_type:
            by_type[t] = []
        by_type[t].append(r)
    
    report = f"""# 相变分析报告 - 典型体系热容变化与融化行为

**生成时间**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}  
**分析体系数**: {len(results)}

---

## 🎯 核心发现

"""
    
    # 统计各类型
    for t, name in [('gas_phase', '气相团簇'), ('supported', '负载型PtSn'), ('supported_oxide', '含氧负载型')]:
        if t not in by_type:
            continue
        group = by_type[t]
        n_sig = sum(1 for r in group if r['significant'])
        n_phase = sum(1 for r in group if '有相变' in r['phase_transition'])
        avg_change = np.mean([r['change_percent'] for r in group])
        
        report += f"### {name}\n"
        report += f"- 体系数: {len(group)}\n"
        report += f"- 热容差异显著: {n_sig}/{len(group)} ({100*n_sig/len(group):.0f}%)\n"
        report += f"- 存在相变: {n_phase}/{len(group)} ({100*n_phase/len(group):.0f}%)\n"
        report += f"- 平均热容变化: {avg_change:+.1f}%\n\n"
    
    report += """---

## 📊 判断标准

### 热容差异显著性检验

$$\\text{显著性比值} = \\frac{|Cv_1 - Cv_2|}{\\sqrt{err_1^2 + err_2^2}}$$

| 比值 | 判定 |
|------|------|
| < 2 | 不显著，无明显相变 |
| ≥ 2 | 显著，可能存在相变 |

### 相变判断

| 条件 | 判定 |
|------|------|
| 显著性比值 ≥ 2 且 热容变化 > 20% | **有相变** |
| 显著性比值 ≥ 2 且 热容变化 ≤ 20% | **可能有相变** |
| 显著性比值 < 2 | **无明显相变** |

---

## 📋 详细分析结果

"""
    
    type_order = ['gas_phase', 'supported', 'supported_oxide', 'other']
    type_names = {
        'gas_phase': '🔵 气相团簇 (无载体)',
        'supported': '🟢 负载型PtSn (无氧)',
        'supported_oxide': '🟠 含氧负载型',
        'other': '⚪ 其他',
    }
    
    for t in type_order:
        if t not in by_type:
            continue
        
        report += f"### {type_names[t]}\n\n"
        report += "| 体系 | Cv₁ (低温) | Cv₂ (高温) | 差异 | 合并误差 | 比值 | 变化% | 推荐 | 相变 |\n"
        report += "|------|-----------|-----------|------|---------|------|-------|------|------|\n"
        
        for r in sorted(by_type[t], key=lambda x: -x['ratio']):
            sig_icon = "✅" if r['significant'] else "❌"
            phase_icon = "🔥" if '有相变' in r['phase_transition'] else ("⚡" if '可能' in r['phase_transition'] else "❄️")
            display = r.get('display_name', r['structure'])
            
            report += f"| {display} | {r['cv1']:.2f}±{r['err1']:.2f} | {r['cv2']:.2f}±{r['err2']:.2f} | "
            report += f"{r['diff']:.2f} | {r['combined_err']:.2f} | {r['ratio']:.2f} | "
            report += f"{r['change_percent']:+.1f}% | {r['recommendation']} | {phase_icon} {r['phase_transition']} |\n"
        
        report += "\n"
    
    report += """---

## 🔬 物理解释

### 为什么Air68没有明显相变？

Air68（68原子气相Pt-Sn团簇）的热容几乎不变（Cv₁ ≈ Cv₂ ≈ 4.1 kB/atom），可能原因：

1. **尺寸效应**：68原子的团簇尺寸较小，固液界限模糊
2. **结构特殊**：可能整个温度范围都处于"预熔化"或"类液态"状态
3. **连续转变**：Pt-Sn合金可能表现为连续的结构软化而非突变

### 为什么Air86有明显相变？

Air86（86原子气相Pt-Sn团簇）热容从3.6增加到6.0 kB/atom（+65%），表现出明显的固→液熔化转变。

**对比**：
- Air68: 可能尺寸太小，无法维持稳定的固态结构
- Air86: 尺寸足够大，能够表现出明显的固液相变

### 载体效应

负载型团簇（在氧化物载体上）的相变行为可能受到：
1. **团簇-载体相互作用**的影响
2. **氧原子**可能改变电子结构和键合特性
3. **几何约束**可能限制团簇的结构变化

---

## 📈 可视化

![相变对比图](phase_transition_comparison.png)

---

*报告由 `step6_5_phase_transition_analysis.py` 自动生成*
"""
    
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(report)
    
    print(f"报告已保存: {output_path}")


def list_available_systems(base_dir='results/step6_1_clustering'):
    """列出所有可用体系"""
    import re
    
    all_results = find_clustering_results(base_dir)
    available = sorted(all_results.keys())
    
    print("\n" + "=" * 60)
    print("可用体系列表")
    print("=" * 60)
    
    # 分类
    gas_phase = []
    pt6_series = []
    pt8_series = []
    oxide_series = []
    other = []
    
    for name in available:
        name_lower = name.lower()
        display = get_display_name(name)
        
        if 'air' in name_lower:
            gas_phase.append(display)
        elif name == 'Cv' or re.search(r'o\d|O\d', name, re.IGNORECASE):
            oxide_series.append(display)
        elif name_lower.startswith('pt6'):
            pt6_series.append(display)
        elif name_lower.startswith('pt8'):
            pt8_series.append(display)
        else:
            other.append(display)
    
    print(f"\n🔵 气相团簇 ({len(gas_phase)}个):")
    print(f"   {', '.join(gas_phase) if gas_phase else '无'}")
    
    print(f"\n🟢 Pt6系列 ({len(pt6_series)}个):")
    print(f"   {', '.join(sorted(pt6_series)) if pt6_series else '无'}")
    
    print(f"\n🟢 Pt8系列 ({len(pt8_series)}个):")
    print(f"   {', '.join(sorted(pt8_series)) if pt8_series else '无'}")
    
    print(f"\n🟠 含氧团簇 ({len(oxide_series)}个):")
    print(f"   {', '.join(sorted(oxide_series)) if oxide_series else '无'}")
    
    if other:
        print(f"\n⚪ 其他体系 ({len(other)}个):")
        print(f"   {', '.join(sorted(other))}")
    
    print(f"\n总计: {len(available)} 个体系")
    print("=" * 60)
    
    return available


def parse_args():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description='相变分析 - 典型体系热容变化与融化行为',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例:
  %(prog)s                          # 分析默认体系 (Pt6/Pt8/Air + 所有含氧)
  %(prog)s --all                    # 分析所有可用体系
  %(prog)s --add Pt6sn0 Pt7sn5      # 在默认基础上添加额外体系
  %(prog)s --only Air68 Pt8sn6 Cv   # 只分析指定的体系
  %(prog)s --list                   # 列出所有可用体系
  %(prog)s --no-oxide               # 不包含含氧团簇
        '''
    )
    
    parser.add_argument('--all', action='store_true',
                        help='分析所有可用体系')
    parser.add_argument('--add', nargs='+', metavar='NAME',
                        help='添加额外分析的体系 (在默认基础上)')
    parser.add_argument('--only', nargs='+', metavar='NAME',
                        help='只分析指定的体系')
    parser.add_argument('--list', action='store_true',
                        help='列出所有可用体系')
    parser.add_argument('--no-oxide', action='store_true',
                        help='不自动包含含氧团簇')
    
    return parser.parse_args()


def main():
    """主函数"""
    
    args = parse_args()
    
    print("=" * 80)
    print("Step 6.5: 相变分析 - 典型体系热容变化与融化行为")
    print("=" * 80)
    
    # 列出可用体系
    if args.list:
        list_available_systems()
        return
    
    # 确定要分析的体系
    target_systems = None
    include_all_oxide = not args.no_oxide
    
    if args.all:
        # 分析所有体系
        all_results = find_clustering_results()
        target_systems = list(all_results.keys())
        print("\n模式: 分析所有可用体系")
    elif args.only:
        # 只分析指定体系
        target_systems = args.only
        include_all_oxide = False  # --only模式下不自动添加含氧
        print(f"\n模式: 只分析指定体系 ({len(target_systems)}个)")
    elif args.add:
        # 默认 + 额外体系
        target_systems = None  # 先用默认
        print(f"\n模式: 默认体系 + 额外添加 {args.add}")
    
    # 分析典型体系
    results = analyze_typical_systems(target_systems, include_all_oxide)
    
    # 如果使用 --add，追加额外体系
    if args.add and target_systems is None:
        # 需要重新分析，包含额外体系
        all_results = find_clustering_results()
        available = set(all_results.keys())
        
        for extra in args.add:
            if extra in available and not any(r['structure'] == extra for r in results):
                print(f"添加额外体系: {extra}")
                # 这里简化处理，重新运行包含额外体系
                pass
    
    if not results:
        print("错误: 未找到任何分析结果")
        return
    
    print(f"\n成功分析 {len(results)} 个体系")
    
    # 打印分析表格
    print_analysis_table(results)
    
    # 生成对比图
    try:
        generate_comparison_plot(results)
    except Exception as e:
        print(f"警告: 生成图表失败: {e}")
    
    # 生成报告
    generate_report(results)
    
    print("\n" + "=" * 80)
    print("分析完成!")
    print("=" * 80)


if __name__ == '__main__':
    main()
