import pandas as pd
import numpy as np

# 读取三个版本的数据
old = pd.read_csv(r'C:\Users\11207\OneDrive\02_Code\work1-PtSnO\workflow\data\lindemann\lindemann_master_run_20251025_205545.csv')
mid = pd.read_csv(r'C:\Users\11207\OneDrive\02_Code\work1-PtSnO\workflow\data\lindemann\collected_lindemann_cluster-20251112\lindemann_master_run_20251112_122604.csv')
new = pd.read_csv(r'C:\Users\11207\OneDrive\02_Code\work1-PtSnO\workflow\data\lindemann\lin-for-all-but-every-ele\lindemann_master_run_20251113_195434.csv')

print("="*80)
print("三个版本 Lindemann 数据对比分析")
print("="*80)
print(f"版本1 (旧): 2025-10-25  →  {len(old)} 条记录")
print(f"版本2 (中): 2025-11-12  →  {len(mid)} 条记录")
print(f"版本3 (新): 2025-11-13  →  {len(new)} 条记录")
print()

# ========================================================================
# 对比 1: 版本1 vs 版本2 (旧 vs 中)
# ========================================================================
print("="*80)
print("📊 对比 1: 版本1 (10-25) vs 版本2 (11-12)")
print("="*80)

merged_12 = pd.merge(old, mid, on='目录', suffixes=('_v1', '_v2'))
merged_12['diff_12'] = merged_12['Lindemann指数_v2'] - merged_12['Lindemann指数_v1']

print(f'\n✅ 匹配记录数: {len(merged_12)} 条 (匹配率: {len(merged_12)/len(old)*100:.1f}%)')

identical_12 = (np.abs(merged_12['diff_12']) < 1e-10).sum()
print(f'\n📈 数值对比:')
print(f'   完全相同: {identical_12} 条 ({identical_12/len(merged_12)*100:.1f}%)')
print(f'   有差异:   {len(merged_12) - identical_12} 条 ({(len(merged_12) - identical_12)/len(merged_12)*100:.1f}%)')
print(f'   平均差异: {merged_12["diff_12"].mean():.8f}')
print(f'   最大正差异: {merged_12["diff_12"].max():.8f}')
print(f'   最大负差异: {merged_12["diff_12"].min():.8f}')

diff_12_001 = (np.abs(merged_12['diff_12']) > 0.001).sum()
diff_12_01 = (np.abs(merged_12['diff_12']) > 0.01).sum()
diff_12_1 = (np.abs(merged_12['diff_12']) > 0.1).sum()
print(f'\n🔍 差异程度:')
print(f'   |差异| > 0.001: {diff_12_001} 条 ({diff_12_001/len(merged_12)*100:.1f}%)')
print(f'   |差异| > 0.01:  {diff_12_01} 条 ({diff_12_01/len(merged_12)*100:.1f}%)')
print(f'   |差异| > 0.1:   {diff_12_1} 条 ({diff_12_1/len(merged_12)*100:.1f}%)')

# ========================================================================
# 对比 2: 版本1 vs 版本3 (旧 vs 新)
# ========================================================================
print("\n" + "="*80)
print("📊 对比 2: 版本1 (10-25) vs 版本3 (11-13)")
print("="*80)

merged_13 = pd.merge(old, new, on='目录', suffixes=('_v1', '_v3'))
merged_13['diff_13'] = merged_13['Lindemann指数_v3'] - merged_13['Lindemann指数_v1']

print(f'\n✅ 匹配记录数: {len(merged_13)} 条 (匹配率: {len(merged_13)/len(old)*100:.1f}%)')

identical_13 = (np.abs(merged_13['diff_13']) < 1e-10).sum()
print(f'\n📈 数值对比:')
print(f'   完全相同: {identical_13} 条 ({identical_13/len(merged_13)*100:.1f}%)')
print(f'   有差异:   {len(merged_13) - identical_13} 条 ({(len(merged_13) - identical_13)/len(merged_13)*100:.1f}%)')
print(f'   平均差异: {merged_13["diff_13"].mean():.8f}')
print(f'   最大正差异: {merged_13["diff_13"].max():.8f}')
print(f'   最大负差异: {merged_13["diff_13"].min():.8f}')

diff_13_001 = (np.abs(merged_13['diff_13']) > 0.001).sum()
diff_13_01 = (np.abs(merged_13['diff_13']) > 0.01).sum()
diff_13_1 = (np.abs(merged_13['diff_13']) > 0.1).sum()
print(f'\n🔍 差异程度:')
print(f'   |差异| > 0.001: {diff_13_001} 条 ({diff_13_001/len(merged_13)*100:.1f}%)')
print(f'   |差异| > 0.01:  {diff_13_01} 条 ({diff_13_01/len(merged_13)*100:.1f}%)')
print(f'   |差异| > 0.1:   {diff_13_1} 条 ({diff_13_1/len(merged_13)*100:.1f}%)')

# ========================================================================
# 对比 3: 版本2 vs 版本3 (中 vs 新)
# ========================================================================
print("\n" + "="*80)
print("📊 对比 3: 版本2 (11-12) vs 版本3 (11-13)")
print("="*80)

merged_23 = pd.merge(mid, new, on='目录', suffixes=('_v2', '_v3'))
merged_23['diff_23'] = merged_23['Lindemann指数_v3'] - merged_23['Lindemann指数_v2']

print(f'\n✅ 匹配记录数: {len(merged_23)} 条 (匹配率: {len(merged_23)/len(mid)*100:.1f}%)')

identical_23 = (np.abs(merged_23['diff_23']) < 1e-10).sum()
print(f'\n📈 数值对比:')
print(f'   完全相同: {identical_23} 条 ({identical_23/len(merged_23)*100:.1f}%)')
print(f'   有差异:   {len(merged_23) - identical_23} 条 ({(len(merged_23) - identical_23)/len(merged_23)*100:.1f}%)')
print(f'   平均差异: {merged_23["diff_23"].mean():.8f}')
print(f'   最大正差异: {merged_23["diff_23"].max():.8f}')
print(f'   最大负差异: {merged_23["diff_23"].min():.8f}')

diff_23_001 = (np.abs(merged_23['diff_23']) > 0.001).sum()
diff_23_01 = (np.abs(merged_23['diff_23']) > 0.01).sum()
diff_23_1 = (np.abs(merged_23['diff_23']) > 0.1).sum()
print(f'\n🔍 差异程度:')
print(f'   |差异| > 0.001: {diff_23_001} 条 ({diff_23_001/len(merged_23)*100:.1f}%)')
print(f'   |差异| > 0.01:  {diff_23_01} 条 ({diff_23_01/len(merged_23)*100:.1f}%)')
print(f'   |差异| > 0.1:   {diff_23_1} 条 ({diff_23_1/len(merged_23)*100:.1f}%)')

# ========================================================================
# 三版本都存在的记录
# ========================================================================
print("\n" + "="*80)
print("📊 三版本共同记录分析")
print("="*80)

# 找到三个版本都有的记录
common_dirs = set(old['目录']) & set(mid['目录']) & set(new['目录'])
print(f'\n✅ 三版本共同记录数: {len(common_dirs)} 条')

old_common = old[old['目录'].isin(common_dirs)].copy()
mid_common = mid[mid['目录'].isin(common_dirs)].copy()
new_common = new[new['目录'].isin(common_dirs)].copy()

# 合并三个版本
merged_all = old_common.merge(mid_common, on='目录', suffixes=('_v1', '_v2'))
merged_all = merged_all.merge(new_common, on='目录')
merged_all.rename(columns={'Lindemann指数': 'Lindemann指数_v3', '结构': '结构_v3', '温度(K)': '温度(K)_v3'}, inplace=True)

merged_all['diff_v1_v2'] = merged_all['Lindemann指数_v2'] - merged_all['Lindemann指数_v1']
merged_all['diff_v2_v3'] = merged_all['Lindemann指数_v3'] - merged_all['Lindemann指数_v2']
merged_all['diff_v1_v3'] = merged_all['Lindemann指数_v3'] - merged_all['Lindemann指数_v1']

# 三版本都相同的记录
all_identical = ((np.abs(merged_all['diff_v1_v2']) < 1e-10) & 
                 (np.abs(merged_all['diff_v2_v3']) < 1e-10)).sum()
print(f'\n📈 三版本数值比较:')
print(f'   三版本完全相同: {all_identical} 条 ({all_identical/len(merged_all)*100:.1f}%)')
print(f'   有任何差异:     {len(merged_all) - all_identical} 条 ({(len(merged_all) - all_identical)/len(merged_all)*100:.1f}%)')

# ========================================================================
# 最大差异记录展示
# ========================================================================
print("\n" + "="*80)
print("⚠️  版本1→版本3 差异最大的前10条记录:")
print("="*80)

top10_13 = merged_13.nlargest(10, 'diff_13', keep='first')[
    ['结构_v1', '温度(K)_v1', 'Lindemann指数_v1', 'Lindemann指数_v3', 'diff_13', '目录']
]
for idx, row in top10_13.iterrows():
    pct = (row['diff_13'] / row['Lindemann指数_v1']) * 100
    print(f'\n🔸 {row["结构_v1"]:20s} @ {row["温度(K)_v1"]:4.0f}K')
    print(f'   版本1: {row["Lindemann指数_v1"]:.6f}')
    print(f'   版本3: {row["Lindemann指数_v3"]:.6f}')
    print(f'   差异:  {row["diff_13"]:+.6f} ({pct:+.1f}%)')
    print(f'   路径:  {row["目录"]}')

# ========================================================================
# 保存详细对比结果
# ========================================================================
print("\n" + "="*80)
print("💾 保存对比结果...")
print("="*80)

# 保存版本1 vs 版本3的差异
diff_13_records = merged_13[np.abs(merged_13['diff_13']) > 1e-10].copy()
diff_13_records = diff_13_records.sort_values('diff_13', key=abs, ascending=False)
output_13 = diff_13_records[['结构_v1', '温度(K)_v1', 'Lindemann指数_v1', 'Lindemann指数_v3', 'diff_13', '目录']].copy()
output_13.columns = ['结构', '温度(K)', 'Lindemann指数_v1(10-25)', 'Lindemann指数_v3(11-13)', '差异', '服务器路径']
output_13.to_csv('lindemann_diff_v1_vs_v3.csv', index=False, encoding='utf-8-sig')
print(f'✅ 已保存 {len(output_13)} 条差异记录到: lindemann_diff_v1_vs_v3.csv')

# 保存版本2 vs 版本3的差异
diff_23_records = merged_23[np.abs(merged_23['diff_23']) > 1e-10].copy()
diff_23_records = diff_23_records.sort_values('diff_23', key=abs, ascending=False)
output_23 = diff_23_records[['结构_v2', '温度(K)_v2', 'Lindemann指数_v2', 'Lindemann指数_v3', 'diff_23', '目录']].copy()
output_23.columns = ['结构', '温度(K)', 'Lindemann指数_v2(11-12)', 'Lindemann指数_v3(11-13)', '差异', '服务器路径']
output_23.to_csv('lindemann_diff_v2_vs_v3.csv', index=False, encoding='utf-8-sig')
print(f'✅ 已保存 {len(output_23)} 条差异记录到: lindemann_diff_v2_vs_v3.csv')

# 保存三版本对比 (只包含三版本都有的记录)
output_all = merged_all[['结构_v1', '温度(K)_v1', 'Lindemann指数_v1', 'Lindemann指数_v2', 
                          'Lindemann指数_v3', 'diff_v1_v2', 'diff_v2_v3', 'diff_v1_v3', '目录']].copy()
output_all.columns = ['结构', '温度(K)', 'Lindemann_v1(10-25)', 'Lindemann_v2(11-12)', 
                       'Lindemann_v3(11-13)', '差异_v1→v2', '差异_v2→v3', '差异_v1→v3', '服务器路径']
output_all.to_csv('lindemann_three_versions_comparison.csv', index=False, encoding='utf-8-sig')
print(f'✅ 已保存 {len(output_all)} 条三版本对比记录到: lindemann_three_versions_comparison.csv')

# ========================================================================
# 打印版本2→版本3有显著差异的路径
# ========================================================================
significant_23 = diff_23_records[np.abs(diff_23_records['diff_23']) > 0.01]
if len(significant_23) > 0:
    print("\n" + "="*80)
    print(f"📁 版本2→版本3 显著差异记录 (|差异| > 0.01) 共 {len(significant_23)} 条:")
    print("="*80)
    print("   可以用以下路径在服务器上检查:\n")
    
    for idx, row in significant_23.head(20).iterrows():  # 只显示前20条
        print(f'   [{row["结构_v2"]:20s} {row["温度(K)_v2"]:4.0f}K] 差异={row["diff_23"]:+.6f}')
        print(f'   → {row["目录"]}')
        print()
    
    if len(significant_23) > 20:
        print(f'   ... 还有 {len(significant_23) - 20} 条记录,详见 CSV 文件')
else:
    print("\n✅ 版本2→版本3 没有显著差异 (|差异| > 0.01)!")

print("\n" + "="*80)
print("💡 总结:")
print("="*80)
print(f"1. 版本1→版本2: {diff_12_01}/{len(merged_12)} 条有显著差异 ({diff_12_01/len(merged_12)*100:.1f}%)")
print(f"2. 版本1→版本3: {diff_13_01}/{len(merged_13)} 条有显著差异 ({diff_13_01/len(merged_13)*100:.1f}%)")
print(f"3. 版本2→版本3: {diff_23_01}/{len(merged_23)} 条有显著差异 ({diff_23_01/len(merged_23)*100:.1f}%)")
print(f"4. 三版本共同记录: {len(common_dirs)} 条,其中 {all_identical} 条完全相同")
print("="*80)
