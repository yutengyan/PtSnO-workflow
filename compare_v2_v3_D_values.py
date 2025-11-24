#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
版本2 vs 版本3 D值详细对比
"""

import pandas as pd
import numpy as np
from pathlib import Path

BASE_DIR = Path(__file__).parent

# 读取数据
v2_file = BASE_DIR / 'data' / 'gmx_msd' / 'unwrap' / 'collected_gmx_per_atom_msd' / 'diffusion_coefficients_gmx_gmx_msd_20251118_152614.csv'
v3_file = BASE_DIR / 'data' / 'gmx_msd' / 'unwrap' / 'collected_gmx_per_atom_msd' / 'diffusion_coefficients_gmx_gmx_msd_20251118_151543.csv'

print("=" * 100)
print("版本2 vs 版本3 D值详细对比".center(100))
print("=" * 100)

df_v2 = pd.read_csv(v2_file)
df_v3 = pd.read_csv(v3_file)

print(f"\n📁 数据加载:")
print(f"   版本2: {len(df_v2)} 条记录")
print(f"   版本3: {len(df_v3)} 条记录")

# 创建完整键
df_v2['full_key'] = (df_v2['完整目录路径'].astype(str) + '||' + 
                     df_v2['结构'].astype(str) + '||' + 
                     df_v2['温度(K)'].astype(str) + '||' + 
                     df_v2['原子组'].astype(str))

df_v3['full_key'] = (df_v3['完整目录路径'].astype(str) + '||' + 
                     df_v3['结构'].astype(str) + '||' + 
                     df_v3['温度(K)'].astype(str) + '||' + 
                     df_v3['原子组'].astype(str))

print("\n" + "=" * 100)
print("📊 样例数据对比 (前10条)")
print("=" * 100)

print("\n【版本2】前10条:")
print(df_v2[['结构', '温度(K)', '原子组', 'D(1e-5 cm²/s)', 'D_err(1e-5 cm²/s)']].head(10).to_string(index=True))

print("\n【版本3】前10条:")
print(df_v3[['结构', '温度(K)', '原子组', 'D(1e-5 cm²/s)', 'D_err(1e-5 cm²/s)']].head(10).to_string(index=True))

# 合并对比
print("\n" + "=" * 100)
print("🔍 逐条D值对比")
print("=" * 100)

merged = df_v2.merge(df_v3, on='full_key', suffixes=('_v2', '_v3'))
print(f"\n成功匹配: {len(merged)} 条记录")

# 计算D值差异
merged['D_diff'] = merged['D(1e-5 cm²/s)_v3'] - merged['D(1e-5 cm²/s)_v2']
merged['D_diff_abs'] = merged['D_diff'].abs()
merged['D_diff_pct'] = (merged['D_diff'] / merged['D(1e-5 cm²/s)_v2'].replace(0, np.nan) * 100)

# 计算D_err差异
merged['D_err_diff'] = merged['D_err(1e-5 cm²/s)_v3'] - merged['D_err(1e-5 cm²/s)_v2']

# 统计
print(f"\n📈 D值差异统计:")
print(f"   最小差异: {merged['D_diff'].min():.10f}")
print(f"   最大差异: {merged['D_diff'].max():.10f}")
print(f"   平均差异: {merged['D_diff'].mean():.10f}")
print(f"   差异标准差: {merged['D_diff'].std():.10f}")
print(f"   差异绝对值平均: {merged['D_diff_abs'].mean():.10f}")

# 不同阈值下的差异统计
thresholds = [1e-10, 1e-8, 1e-6, 1e-4, 1e-2]
print(f"\n📊 不同阈值下的差异分布:")
for thresh in thresholds:
    count = (merged['D_diff_abs'] >= thresh).sum()
    pct = count / len(merged) * 100
    print(f"   |差异| ≥ {thresh:.0e}: {count:5d} 条 ({pct:5.2f}%)")

# 完全相同的记录
identical = (merged['D_diff_abs'] < 1e-10).sum()
print(f"\n✅ 完全相同 (差异<1e-10): {identical} 条 ({identical/len(merged)*100:.2f}%)")

# D_err差异统计
print(f"\n📈 D_err差异统计:")
print(f"   最小差异: {merged['D_err_diff'].min():.10f}")
print(f"   最大差异: {merged['D_err_diff'].max():.10f}")
print(f"   平均差异: {merged['D_err_diff'].mean():.10f}")

# 差异最大的记录
if merged['D_diff_abs'].max() > 1e-10:
    print("\n" + "=" * 100)
    print("⚠️  差异最大的前20条记录")
    print("=" * 100)
    
    top_diff = merged.nlargest(20, 'D_diff_abs')
    for i, (idx, row) in enumerate(top_diff.iterrows(), 1):
        print(f"\n[{i}] {row['结构_v2']} @ {row['温度(K)_v2']}K - {row['原子组_v2']}")
        print(f"    路径: {row['完整目录路径_v2']}")
        print(f"    D_v2  = {row['D(1e-5 cm²/s)_v2']:12.8f}")
        print(f"    D_v3  = {row['D(1e-5 cm²/s)_v3']:12.8f}")
        print(f"    差异   = {row['D_diff']:+12.8f}  ({row['D_diff_pct']:.4f}%)")
        print(f"    D_err_v2 = {row['D_err(1e-5 cm²/s)_v2']:12.8f}")
        print(f"    D_err_v3 = {row['D_err(1e-5 cm²/s)_v3']:12.8f}")

# 随机抽样10条对比
print("\n" + "=" * 100)
print("🎲 随机抽样10条对比")
print("=" * 100)

sample = merged.sample(min(10, len(merged)))
for i, (idx, row) in enumerate(sample.iterrows(), 1):
    print(f"\n[{i}] {row['结构_v2']} @ {row['温度(K)_v2']}K - {row['原子组_v2']}")
    print(f"    D_v2 = {row['D(1e-5 cm²/s)_v2']:12.8f}  ±{row['D_err(1e-5 cm²/s)_v2']:12.8f}")
    print(f"    D_v3 = {row['D(1e-5 cm²/s)_v3']:12.8f}  ±{row['D_err(1e-5 cm²/s)_v3']:12.8f}")
    print(f"    差异  = {row['D_diff']:+12.8f}")

# 保存详细对比
output_file = BASE_DIR / 'v2_v3_D_value_comparison.csv'
merged_output = merged[['结构_v2', '温度(K)_v2', '原子组_v2', '完整目录路径_v2',
                        'D(1e-5 cm²/s)_v2', 'D(1e-5 cm²/s)_v3', 'D_diff', 'D_diff_abs',
                        'D_err(1e-5 cm²/s)_v2', 'D_err(1e-5 cm²/s)_v3', 'D_err_diff']].copy()
merged_output.columns = ['结构', '温度(K)', '原子组', '路径', 
                         'D_v2', 'D_v3', 'D差异', 'D差异绝对值',
                         'D_err_v2', 'D_err_v3', 'D_err差异']
merged_output.to_csv(output_file, index=False)

print("\n" + "=" * 100)
print("💾 保存结果")
print("=" * 100)
print(f"详细对比已保存至: {output_file}")

print("\n" + "=" * 100)
print("📝 总结")
print("=" * 100)

print(f"""
版本2 vs 版本3 D值对比结果:
----------------------------
• 总记录数: {len(merged)}
• 完全相同 (差异<1e-10): {identical} ({identical/len(merged)*100:.2f}%)
• 有差异 (差异≥1e-10): {len(merged)-identical} ({(len(merged)-identical)/len(merged)*100:.2f}%)
• 最大差异: {merged['D_diff_abs'].max():.10f}
• 平均差异: {merged['D_diff'].mean():.10f}

结论: {'两版本D值完全一致!' if identical == len(merged) else f'存在 {len(merged)-identical} 条记录有差异'}
""")

print("=" * 100)
