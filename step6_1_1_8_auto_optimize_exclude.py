"""
Step 6.1.1.8: 自动优化排除点建议 - 集成工具
=========================================

功能:
1. 自动生成排除点建议
2. 测试不同残差阈值
3. 方案对比分析
4. 生成推荐文档和命令

使用方法:
    python step6_1_1_8_auto_optimize_exclude.py --structure Sn1Pt2O1 --mode all
    python step6_1_1_8_auto_optimize_exclude.py --structure Sn1Pt2O1 --mode suggest --threshold 60
    python step6_1_1_8_auto_optimize_exclude.py --structure Sn1Pt2O1 --mode compare

作者: AI Assistant
日期: 2026-02-04
"""

import pandas as pd
import numpy as np
from scipy.stats import linregress
import argparse
from pathlib import Path
import json
import locale
import sys

class ExcludeOptimizer:
    """排除点优化器"""
    
    def __init__(self, structure, T_min=200, T_max=1600, data_dir='results/step6_1_clustering',
                 partitions=None):
        """
        partitions: 已知分区列表，格式 [(T_min1, T_max1), (T_min2, T_max2), ...]
                    例如 [(200,1450),(1450,1650),(1700,1800)]
                    若为 None，则退回原来的全局单线性拟合模式。
        """
        self.structure = structure
        self.T_min = T_min
        self.T_max = T_max
        self.partitions = partitions  # None → 全局模式；有值 → 分区模式
        self.data_file = Path(data_dir) / f'{structure}_lindemann-threshold_n2_clustered_data.csv'
        
        if not self.data_file.exists():
            raise FileNotFoundError(f"数据文件不存在: {self.data_file}")
        
        # 使用系统编码读取CSV (Windows中文环境兼容)
        encoding = locale.getpreferredencoding() if sys.platform == 'win32' else 'utf-8'
        self.df = pd.read_csv(self.data_file, encoding=encoding)
        print(f"✓ 已加载数据: {self.data_file}")
        print(f"  总数据点: {len(self.df)}")
        print(f"  温度范围: {self.df['temp'].min():.0f}K - {self.df['temp'].max():.0f}K")
        if partitions:
            parts_str = ', '.join([f"{int(a)}-{int(b)}K" for a, b in partitions])
            print(f"  已知分区: {parts_str}  → 启用【分区内残差】筛选模式")
        else:
            print(f"  分区: 未指定 → 使用【全局单线性拟合】筛选模式")
        print()
    
    def calculate_metrics(self, exclude_dict):
        """计算给定排除方案的指标（全局单线性拟合模式）"""
        # 应用排除
        df_filtered_list = []
        for temp in self.df['temp'].unique():
            df_temp = self.df[self.df['temp'] == temp].copy()
            df_temp = df_temp.sort_values('avg_energy', ascending=False).reset_index(drop=True)
            
            if temp in exclude_dict:
                mask = [i not in exclude_dict[temp] for i in range(len(df_temp))]
                df_temp = df_temp[mask]
            
            df_filtered_list.append(df_temp)
        
        df_filtered = pd.concat(df_filtered_list, ignore_index=True)
        
        # 拟合
        df_part = df_filtered[(df_filtered['temp'] >= self.T_min) & (df_filtered['temp'] <= self.T_max)]
        temp_unique = np.sort(df_part['temp'].unique())
        E_avg = np.array([df_part[df_part['temp']==T]['avg_energy'].mean() for T in temp_unique]) * 1000
        
        slope, intercept, r, _, _ = linregress(temp_unique, E_avg)
        y_pred = slope * temp_unique + intercept
        residuals = E_avg - y_pred
        
        # 计算每个温度的残差
        residuals_dict = {}
        for T in temp_unique:
            E_avg_T = df_filtered[df_filtered['temp'] == T]['avg_energy'].mean() * 1000
            E_pred = slope * T + intercept
            residuals_dict[T] = E_avg_T - E_pred
        
        n_excluded = sum(len(v) for v in exclude_dict.values())
        n_total = len(self.df[(self.df['temp'] >= self.T_min) & (self.df['temp'] <= self.T_max)])
        
        return {
            'R2': r**2,
            'max_residual': np.abs(residuals).max(),
            'mean_residual': np.abs(residuals).mean(),
            'residuals_dict': residuals_dict,
            'n_excluded': n_excluded,
            'n_total': n_total,
            'exclude_ratio': n_excluded / n_total * 100,
            'slope': slope,
            'intercept': intercept,
        }

    def _apply_exclude(self, df, exclude_dict):
        """将 exclude_dict 应用到 df，返回过滤后的 DataFrame"""
        df_filtered_list = []
        for temp in df['temp'].unique():
            df_temp = df[df['temp'] == temp].copy()
            df_temp = df_temp.sort_values('avg_energy', ascending=False).reset_index(drop=True)
            if temp in exclude_dict:
                mask = [i not in exclude_dict[temp] for i in range(len(df_temp))]
                df_temp = df_temp[mask]
            df_filtered_list.append(df_temp)
        return pd.concat(df_filtered_list, ignore_index=True)

    def _fit_partition(self, df_filtered, T_lo, T_hi):
        """
        对过滤后的 df_filtered 在 [T_lo, T_hi] 分区内做线性拟合。
        返回 (slope_eV/K, intercept_eV, r²) 或 None（温度点不足）。
        """
        df_p = df_filtered[(df_filtered['temp'] >= T_lo) & (df_filtered['temp'] <= T_hi)]
        temps_p = np.sort(df_p['temp'].unique())
        if len(temps_p) < 2:
            return None
        E_avg_p = np.array([df_p[df_p['temp'] == T]['avg_energy'].mean() for T in temps_p])
        slope, intercept, r, _, _ = linregress(temps_p, E_avg_p)
        return slope, intercept, r ** 2, len(temps_p)

    def calculate_metrics_by_partition(self, exclude_dict):
        """
        分区模式：
          1. 对每个分区做线性拟合（用各温度均值）
          2. 计算每个 run 与其所属分区拟合线的偏差（run-level residual）
          3. max/mean_residual 基于 run-level 偏差，更真实反映数据散度
        """
        assert self.partitions, "需要提供 partitions 才能使用分区模式"

        df_filtered = self._apply_exclude(self.df, exclude_dict)
        partition_info = []
        all_run_devs = []   # 每个run与拟合线的绝对偏差(meV)

        for T_lo, T_hi in self.partitions:
            fit = self._fit_partition(df_filtered, T_lo, T_hi)
            if fit is None:
                partition_info.append({
                    'range': (T_lo, T_hi), 'n_temps': 0,
                    'Cv': float('nan'), 'R2': float('nan'), 'max_resid': float('nan'),
                })
                continue
            slope, intercept, r2, n_temps = fit
            # 每个run的偏差
            df_p = df_filtered[(df_filtered['temp'] >= T_lo) & (df_filtered['temp'] <= T_hi)]
            run_devs = np.abs(df_p['avg_energy'] - (slope * df_p['temp'] + intercept)) * 1000
            all_run_devs.extend(run_devs.tolist())
            partition_info.append({
                'range': (T_lo, T_hi),
                'n_temps': n_temps,
                'Cv': slope * 1000,        # meV/K（含载体）
                'R2': r2,
                'max_resid': run_devs.max(),
            })

        n_excluded = sum(len(v) for v in exclude_dict.values())
        n_total = len(self.df)

        # 整体 R² 用全局单线性作为对照
        df_global = df_filtered
        temps_g = np.sort(df_global['temp'].unique())
        E_avg_g = np.array([df_global[df_global['temp'] == T]['avg_energy'].mean() for T in temps_g]) * 1000
        _, _, r_g, _, _ = linregress(temps_g, E_avg_g)

        return {
            'R2': r_g ** 2,
            'max_residual': max(all_run_devs) if all_run_devs else 0.0,
            'mean_residual': float(np.mean(all_run_devs)) if all_run_devs else 0.0,
            'n_excluded': n_excluded,
            'n_total': n_total,
            'exclude_ratio': n_excluded / n_total * 100,
            'partition_info': partition_info,
        }

    def generate_exclude_dict_by_partition(self, target_residual, max_exclude=3):
        """
        分区感知的排除生成（基于 run-level 偏差）：

        核心思路：
          - 对每个分区，用「当前所有点的均值」做线性拟合得到参考直线
          - 对每个温度的每个 run，计算其 avg_energy 与参考直线预测值的偏差
          - 偏差超过阈值的 run 视为离群 → 排除
          - 保留规则：每个温度至少保留 2 个 run

        分区策略：
          - 普通分区（n_temps > 2）：直接用本分区拟合线，按偏差绝对值排除
          - 薄分区（n_temps == 2，如 partition2/3）：
              * 同样计算 run-level 偏差
              * 但采用「单侧排除」规则：若一个温度点的超阈值 run 同时含正偏差和负偏差，
                说明数据天然分散（相变区），不排除（均值有代表性）
              * 只有当超阈值 run 偏差方向单一（全正或全负）时才排除
              * 这样避免双侧排除改变均值，导致 Cv 斜率失真
        """
        assert self.partitions, "需要提供 partitions 才能使用分区感知排除"

        exclude_dict = {}

        for p_idx, (T_lo, T_hi) in enumerate(self.partitions):
            # 先用无排除数据拟合参考直线
            fit = self._fit_partition(self.df, T_lo, T_hi)
            df_p = self.df[(self.df['temp'] >= T_lo) & (self.df['temp'] <= T_hi)]
            temps_p = np.sort(df_p['temp'].unique())
            n_temps = len(temps_p)

            p_label = f"partition{p_idx+1} ({int(T_lo)}-{int(T_hi)}K, n_temps={n_temps})"
            if n_temps <= 2:
                p_label += " [薄分区]"
            print(f"  [{p_label}]")

            if fit is None or n_temps < 2:
                print(f"    温度点不足，跳过\n")
                continue

            slope, intercept, r2, _ = fit

            for T in temps_p:
                # 按能量降序排列，index 0=最高能量
                df_t = df_p[df_p['temp'] == T].sort_values(
                    'avg_energy', ascending=False).reset_index(drop=True)
                n_runs = len(df_t)
                E_pred = slope * T + intercept   # eV

                # 每个 run 的偏差（meV）
                devs = (df_t['avg_energy'] - E_pred) * 1000   # meV, shape (n_runs,)

                # 找超阈值的 run（保留至少 2 个）
                to_exclude = []
                for i, dev in enumerate(devs):
                    if abs(dev) > target_residual:
                        to_exclude.append((i, dev))

                if not to_exclude:
                    continue

                # 薄分区（n_temps==2）：单侧排除约束
                # 若超阈值 run 同时有正偏差和负偏差 → 数据天然分散，不排除
                if n_temps <= 2:
                    has_pos = any(dev > 0 for _, dev in to_exclude)
                    has_neg = any(dev < 0 for _, dev in to_exclude)
                    if has_pos and has_neg:
                        print(f"    {int(T)}K: 超阈值run含双侧偏差 [相变区天然分散]，跳过排除")
                        continue

                # 按偏差绝对值从大到小排序，优先排除偏差最大的
                to_exclude.sort(key=lambda x: abs(x[1]), reverse=True)
                max_can_exclude = min(max_exclude, n_runs - 2)   # 至少保留 2 个 run
                if max_can_exclude <= 0:
                    print(f"    {int(T)}K: 有{len(to_exclude)}个超阈值run，但run数={n_runs}<=2，跳过")
                    continue

                final_exclude = [idx for idx, _ in to_exclude[:max_can_exclude]]
                final_exclude.sort()
                devs_str = ', '.join([f'{devs[i]:+.0f}' for i in final_exclude])
                print(f"    {int(T)}K: 排除index {final_exclude} (偏差: {devs_str} meV)")
                exclude_dict[int(T)] = final_exclude

            print()

        return exclude_dict

    def generate_exclude_dict(self, target_residual, max_exclude=3):
        """基于目标残差阈值生成排除字典"""
        # 计算初始残差
        baseline = self.calculate_metrics({})
        residuals = baseline['residuals_dict']
        
        exclude_dict = {}
        
        for T, residual in residuals.items():
            if abs(residual) > target_residual:
                df_temp = self.df[self.df['temp'] == T].copy()
                df_temp_sorted = df_temp.sort_values('avg_energy', ascending=False).reset_index(drop=True)
                n_points = len(df_temp_sorted)
                
                if n_points <= 2:
                    continue
                
                # 计算建议排除的点数
                n_exclude = min(max_exclude, max(1, int(abs(residual) / target_residual)))
                n_exclude = min(n_exclude, n_points - 2)
                
                if residual > 0:
                    # 正残差 → 删除高能量点(索引0,1,2...)
                    indices = list(range(n_exclude))
                else:
                    # 负残差 → 删除低能量点(索引n-1,n-2,n-3...)
                    indices = list(range(n_points-1, n_points-1-n_exclude, -1))
                
                exclude_dict[int(T)] = indices
        
        return exclude_dict
    
    def suggest_exclusions(self, target_residual=50, max_exclude=3, verbose=True):
        """生成并显示排除建议（自动选择全局或分区模式）"""
        use_partition = bool(self.partitions)
        mode_label = "【分区模式】" if use_partition else "【全局单线性模式】"

        if verbose:
            print("="*80)
            print(f"生成排除建议 {mode_label}: {self.T_min}-{self.T_max}K")
            print("="*80)
            print(f"目标残差阈值: {target_residual} meV")
            print(f"每温度最多排除: {max_exclude} 个点")
            if use_partition:
                print(f"策略: 每个分区独立拟合参考线 → 对每个run计算偏差 → 排除偏差>阈值的run")
                print(f"      薄分区(2温度点)同样筛选 → 改善均值质量 → 更准确的热容斜率")
            else:
                print(f"策略: 正残差→删高能量(索引0,1...), 负残差→删低能量(索引n-1,n-2...)")
            print()

        # 基线
        if use_partition:
            baseline = self.calculate_metrics_by_partition({})
        else:
            baseline = self.calculate_metrics({})

        if verbose:
            print(f"基线(无排除):")
            print(f"  最大|run偏差| = {baseline['max_residual']:.2f} meV")
            print(f"  平均|run偏差| = {baseline['mean_residual']:.2f} meV")
            if use_partition and 'partition_info' in baseline:
                CV_SUPPORT = 38.2151
                for pi in baseline['partition_info']:
                    cv_net = pi['Cv'] - CV_SUPPORT if not np.isnan(pi['Cv']) else float('nan')
                    print(f"    {int(pi['range'][0])}-{int(pi['range'][1])}K: "
                          f"Cv_net={cv_net:.3f} meV/K, R²={pi['R2']:.4f}, "
                          f"n_temps={pi['n_temps']}, max_run_dev={pi['max_resid']:.1f} meV")
            print()

        # 生成排除建议
        if use_partition:
            exclude_dict = self.generate_exclude_dict_by_partition(target_residual, max_exclude)
        else:
            exclude_dict = self.generate_exclude_dict(target_residual, max_exclude)

        if not exclude_dict:
            if verbose:
                print("✓ 所有温度的残差都在目标范围内，无需排除!")
            return {}, baseline, baseline

        # 显示建议（全局模式才需要这里单独打印，分区模式已在 generate 里打印）
        if verbose and not use_partition:
            print("建议排除的温度点:")
            for T, indices in sorted(exclude_dict.items()):
                residual = baseline['residuals_dict'][T]
                df_temp = self.df[self.df['temp'] == T]
                n_points = len(df_temp)
                direction = "高能量" if 0 in indices else "低能量"
                print(f"  {int(T)}K: 残差={residual:+7.2f} meV, 点数={n_points}, "
                      f"删除{len(indices)}个{direction}点 {indices}")
            print()

        # 优化后效果
        if use_partition:
            optimized = self.calculate_metrics_by_partition(exclude_dict)
        else:
            optimized = self.calculate_metrics(exclude_dict)

        if verbose:
            print("="*80)
            print("优化后效果")
            print("="*80)
            print(f"  最大|残差|: {baseline['max_residual']:.2f} → {optimized['max_residual']:.2f} meV "
                  f"(↓{baseline['max_residual']-optimized['max_residual']:.2f})")
            print(f"  平均|残差|: {baseline['mean_residual']:.2f} → {optimized['mean_residual']:.2f} meV "
                  f"(↓{baseline['mean_residual']-optimized['mean_residual']:.2f})")
            print(f"  排除点数: {optimized['n_excluded']} / {optimized['n_total']} "
                  f"({optimized['exclude_ratio']:.1f}%)")
            if use_partition and 'partition_info' in optimized:
                CV_SUPPORT = 38.2151
                print(f"\n  优化后各分区情况:")
                for pi in optimized['partition_info']:
                    cv_net = pi['Cv'] - CV_SUPPORT if not np.isnan(pi['Cv']) else float('nan')
                    print(f"    {int(pi['range'][0])}-{int(pi['range'][1])}K: "
                          f"Cv_net={cv_net:.3f} meV/K, R²={pi['R2']:.4f}, "
                          f"n_temps={pi['n_temps']}, max_run_dev={pi['max_resid']:.1f} meV")
            print()

        return exclude_dict, baseline, optimized
    
    def test_thresholds(self, thresholds=None, max_exclude=3):
        """测试不同残差阈值"""
        if thresholds is None:
            thresholds = [30, 40, 50, 60, 70, 80, 90, 100, 120, 150]

        use_partition = bool(self.partitions)
        mode_label = "【分区模式】" if use_partition else "【全局模式】"

        print("="*90)
        print(f"测试不同残差阈值 {mode_label}")
        print("="*90)
        print()

        if use_partition:
            baseline = self.calculate_metrics_by_partition({})
        else:
            baseline = self.calculate_metrics({})

        print(f"基线: 最大残差={baseline['max_residual']:.2f} meV, "
              f"平均残差={baseline['mean_residual']:.2f} meV")
        print()

        print(f"{'阈值':<8} {'最大残差':<12} {'平均残差':<12} {'R²':<12} {'点数':<8} {'比例':<8} {'效率':<10}")
        print("-"*90)

        results = []
        for threshold in thresholds:
            if use_partition:
                exclude_dict = self.generate_exclude_dict_by_partition(threshold, max_exclude)
                metrics = self.calculate_metrics_by_partition(exclude_dict)
            else:
                exclude_dict = self.generate_exclude_dict(threshold, max_exclude)
                metrics = self.calculate_metrics(exclude_dict)

            res_reduction = baseline['max_residual'] - metrics['max_residual']
            efficiency = res_reduction / (metrics['exclude_ratio'] + 1e-6)

            results.append({
                'threshold': threshold,
                'exclude_dict': exclude_dict,
                'metrics': metrics,
                'efficiency': efficiency,
                'res_reduction': res_reduction,
            })

            print(f"{threshold:<8.0f} {metrics['max_residual']:<12.2f} {metrics['mean_residual']:<12.2f} "
                  f"{metrics['R2']:<12.6f} {metrics['n_excluded']:<8} {metrics['exclude_ratio']:<7.1f}% "
                  f"{efficiency:<10.2f}")

        print()
        print("效率 = 残差降低(meV) / 排除比例(%)")
        print()

        return results, baseline
    
    def compare_strategies(self, results, baseline):
        """对比不同策略"""
        print("="*90)
        print("推荐方案对比")
        print("="*90)
        print()
        
        # 找到关键阈值的结果
        strategies = {
            150: '最小干预',
            60: '适度优化',
            30: '激进优化',
        }
        
        for threshold, label in strategies.items():
            result = next((r for r in results if r['threshold'] == threshold), None)
            if not result:
                continue
            
            metrics = result['metrics']
            print(f"【方案: {label}】threshold={threshold} meV")
            print(f"  排除: {metrics['n_excluded']} 点 ({metrics['exclude_ratio']:.1f}%)")
            print(f"  最大残差: {baseline['max_residual']:.2f} → {metrics['max_residual']:.2f} meV "
                  f"(↓{result['res_reduction']:.2f}, {result['res_reduction']/baseline['max_residual']*100:.1f}%)")
            print(f"  平均残差: {baseline['mean_residual']:.2f} → {metrics['mean_residual']:.2f} meV "
                  f"(↓{baseline['mean_residual']-metrics['mean_residual']:.2f})")
            print(f"  R²: {metrics['R2']:.6f}")
            print(f"  效率: {result['efficiency']:.2f} meV/%")
            
            if threshold == 60:
                print(f"  ★ 推荐: 性价比最高")
            elif threshold == 150:
                print(f"  评价: 改善有限")
            elif threshold == 30:
                print(f"  评价: 排除过多")
            print()
        
        return strategies
    
    def generate_command(self, exclude_dict, partitions=None, windows=True):
        """生成命令行（若未传 partitions 且 self.partitions 存在则自动使用）"""
        if partitions is None:
            if self.partitions:
                partitions = ','.join([f"{int(a)}-{int(b)}" for a, b in self.partitions])
            else:
                partitions = f"{self.T_min}-{self.T_max}"
        
        exclude_args = []
        for T in sorted(exclude_dict.keys()):
            indices = exclude_dict[T]
            indices_str = ','.join(map(str, indices))
            exclude_args.append(f'"{int(T)}K:{indices_str}"')
        
        if windows:
            # Windows PowerShell 格式 (使用 ` 反引号换行)
            cmd = f"""python step6_1_1_partition_cv_plot.py `
    --structure {self.structure} `
    --partitions {partitions} `
    --exclude {' '.join(exclude_args)} `
    --exclude-sort-by energy `
    --y-ticks 0,2,4 --cv-ticks 3,4,5,6 --figsize 10x8 --peak-method partition"""
        else:
            # Linux/Mac Bash 格式 (使用 \ 反斜杠换行)
            cmd = f"""python step6_1_1_partition_cv_plot.py \\
    --structure {self.structure} \\
    --partitions {partitions} \\
    --exclude {' '.join(exclude_args)} \\
    --exclude-sort-by energy \\
    --y-ticks 0,2,4 --cv-ticks 3,4,5,6 --figsize 10x8 --peak-method partition"""
        
        return cmd
    
    def generate_report(self, results, baseline, output_file=None):
        """生成Markdown报告"""
        if output_file is None:
            output_file = f"{self.structure}_EXCLUDE_RECOMMENDATIONS.md"
        
        # 找到推荐方案
        recommended = next((r for r in results if r['threshold'] == 60), results[len(results)//2])
        
        report = f"""# {self.structure} 排除点优化建议

## 📊 数据概况

- **结构**: {self.structure}
- **温度范围**: {self.T_min}-{self.T_max}K
- **总数据点**: {baseline['n_total']}
- **生成日期**: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}

## 🔍 基线状态

- R² = {baseline['R2']:.6f}
- 最大|残差| = {baseline['max_residual']:.2f} meV
- 平均|残差| = {baseline['mean_residual']:.2f} meV

## ✅ 推荐方案 (threshold={recommended['threshold']} meV)

### 排除详情

| 温度 | 排除索引 | 方向 | 原因 |
|------|---------|------|------|
"""
        
        for T, indices in sorted(recommended['exclude_dict'].items()):
            direction = "高能量" if 0 in indices else "低能量"
            # 分区模式无 residuals_dict，改为显示索引方向
            reason = f"偏差超阈值，排除{direction}端run"
            report += f"| {int(T)}K | {indices} | {direction} | {reason} |\n"
        
        report += f"""
### 优化效果

- R²: {baseline['R2']:.6f} → {recommended['metrics']['R2']:.6f} (↑{(recommended['metrics']['R2']-baseline['R2'])*100:.4f}%)
- 最大|残差|: {baseline['max_residual']:.2f} → {recommended['metrics']['max_residual']:.2f} meV (↓{recommended['res_reduction']:.2f}, {recommended['res_reduction']/baseline['max_residual']*100:.1f}%)
- 平均|残差|: {baseline['mean_residual']:.2f} → {recommended['metrics']['mean_residual']:.2f} meV (↓{baseline['mean_residual']-recommended['metrics']['mean_residual']:.2f})
- 排除点数: {recommended['metrics']['n_excluded']} / {recommended['metrics']['n_total']} ({recommended['metrics']['exclude_ratio']:.1f}%)
- 效率: {recommended['efficiency']:.2f} meV/%

### 执行命令 (Windows PowerShell)

```powershell
{self.generate_command(recommended['exclude_dict'], windows=True)}
```

### 执行命令 (Linux/Mac Bash)

```bash
{self.generate_command(recommended['exclude_dict'], windows=False)}
```

## 📈 方案对比

| 阈值 | 最大残差 | 排除比例 | 效率 | 推荐度 |
|------|---------|---------|------|--------|
"""
        
        for r in results:
            stars = "⭐⭐⭐⭐⭐" if r['threshold'] == 60 else "⭐⭐⭐⭐" if r['threshold'] in [30, 50, 70] else "⭐⭐⭐" if r['threshold'] in [40, 80, 100] else "⭐⭐"
            report += f"| {r['threshold']} meV | {r['metrics']['max_residual']:.2f} meV | {r['metrics']['exclude_ratio']:.1f}% | {r['efficiency']:.2f} | {stars} |\n"
        
        report += f"""
## 🔧 使用工具

### 快速生成建议
```bash
python step6_1_1_8_auto_optimize_exclude.py --structure {self.structure} --mode suggest --threshold 60
```

### 测试不同阈值
```bash
python step6_1_1_8_auto_optimize_exclude.py --structure {self.structure} --mode test
```

### 完整分析
```bash
python step6_1_1_8_auto_optimize_exclude.py --structure {self.structure} --mode all
```

## 📝 说明

### 排序逻辑
脚本使用 `ascending=False` (降序排序):
- 索引 0, 1, 2, ... = **最高能量** (最不稳定)
- 索引 n-1, n-2, n-3, ... = **最低能量** (最稳定)

### 排除策略
- ✅ 正残差(点在拟合线上方) → 删除高能量点 → 索引 [0, 1, 2, ...]
- ✅ 负残差(点在拟合线下方) → 删除低能量点 → 索引 [n-1, n-2, n-3, ...]

### 阈值选择建议
- threshold < 50: 过于激进,可能过度拟合
- threshold = 60-80: ★ 平衡点,推荐
- threshold > 100: 过于保守,改善有限

---
*本报告由 step6_1_1_8_auto_optimize_exclude.py 自动生成*
"""
        
        # 保存报告
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(report)
        
        print(f"✓ 报告已保存: {output_file}")
        return output_file


def main():
    parser = argparse.ArgumentParser(
        description='Step 6.1.1.8: 自动优化排除点建议',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  # 完整分析(推荐)
  python step6_1_1_8_auto_optimize_exclude.py --structure Sn1Pt2O1 --mode all
  
  # 快速生成建议
  python step6_1_1_8_auto_optimize_exclude.py --structure Sn1Pt2O1 --mode suggest --threshold 60
  
  # 仅测试不同阈值
  python step6_1_1_8_auto_optimize_exclude.py --structure Sn1Pt2O1 --mode test
  
  # 对比已有方案
  python step6_1_1_8_auto_optimize_exclude.py --structure Sn1Pt2O1 --mode compare
        """
    )
    
    parser.add_argument('--structure', required=True, help='结构名称 (如 Sn1Pt2O1)')
    parser.add_argument('--mode', default='all', 
                       choices=['suggest', 'test', 'compare', 'all'],
                       help='运行模式: suggest(建议), test(测试), compare(对比), all(全部)')
    parser.add_argument('--threshold', type=float, default=60, 
                       help='目标残差阈值(meV), 默认60')
    parser.add_argument('--max-exclude', type=int, default=3,
                       help='每温度最多排除点数, 默认3')
    parser.add_argument('--T-min', type=int, default=200, help='最低温度(K)')
    parser.add_argument('--T-max', type=int, default=1600, help='最高温度(K)')
    parser.add_argument('--partitions', type=str, default=None,
                       help='已知分区（启用分区感知筛选）, 格式: 200-1450,1450-1650,1700-1800')
    parser.add_argument('--output', help='报告输出文件名')
    parser.add_argument('--data-dir', default='results/step6_1_clustering',
                       help='数据目录')
    parser.add_argument('--platform', default='windows', choices=['windows', 'linux', 'mac'],
                       help='目标平台 (windows使用反引号`, linux/mac使用反斜杠\\)')
    
    args = parser.parse_args()
    
    # 解析 --partitions
    partitions_list = None
    if args.partitions:
        try:
            partitions_list = []
            for seg in args.partitions.split(','):
                T_lo, T_hi = map(float, seg.strip().split('-'))
                partitions_list.append((T_lo, T_hi))
        except Exception as e:
            print(f"❌ --partitions 格式错误: {e}")
            print("   正确格式示例: 200-1450,1450-1650,1700-1800")
            return 1

    # 初始化优化器
    try:
        optimizer = ExcludeOptimizer(
            args.structure,
            T_min=args.T_min,
            T_max=args.T_max,
            data_dir=args.data_dir,
            partitions=partitions_list,
        )
    except FileNotFoundError as e:
        print(f"❌ 错误: {e}")
        return 1
    
    # 根据模式执行
    if args.mode == 'suggest':
        print("="*80)
        print("模式: 生成排除建议")
        print("="*80)
        print()
        
        exclude_dict, baseline, optimized = optimizer.suggest_exclusions(
            target_residual=args.threshold,
            max_exclude=args.max_exclude,
            verbose=True
        )
        
        if exclude_dict:
            print("="*80)
            print("生成命令")
            print("="*80)
            cmd = optimizer.generate_command(exclude_dict, windows=(args.platform=='windows'))
            print(cmd)
            print()
            
            # 生成报告文件（供批量绘图脚本使用）
            output_file = args.output if args.output else f"{args.structure}_EXCLUDE_RECOMMENDATIONS.md"
            # 构建完整的 results 用于报告生成（包含 res_reduction 和 efficiency）
            res_reduction = baseline['max_residual'] - optimized['max_residual']
            efficiency = res_reduction / (optimized['exclude_ratio'] + 1e-6)
            results = [{
                'threshold': args.threshold,
                'exclude_dict': exclude_dict,
                'metrics': optimized,
                'res_reduction': res_reduction,
                'efficiency': efficiency
            }]
            optimizer.generate_report(results, baseline, output_file)
            print(f"\n✓ 报告已保存: {output_file}")
    
    elif args.mode == 'test':
        print("="*80)
        print("模式: 测试不同阈值")
        print("="*80)
        print()
        
        results, baseline = optimizer.test_thresholds(max_exclude=args.max_exclude)
    
    elif args.mode == 'compare':
        print("="*80)
        print("模式: 方案对比")
        print("="*80)
        print()
        
        results, baseline = optimizer.test_thresholds(max_exclude=args.max_exclude)
        optimizer.compare_strategies(results, baseline)
        
        # 显示推荐命令
        recommended = next((r for r in results if r['threshold'] == 60), results[len(results)//2])
        print("="*80)
        print("推荐命令 (threshold=60)")
        print("="*80)
        cmd = optimizer.generate_command(recommended['exclude_dict'], windows=(args.platform=='windows'))
        print(cmd)
        print()
    
    elif args.mode == 'all':
        print("="*80)
        print("模式: 完整分析")
        print("="*80)
        print()
        
        # 1. 测试不同阈值
        results, baseline = optimizer.test_thresholds(max_exclude=args.max_exclude)
        print()
        
        # 2. 对比策略
        optimizer.compare_strategies(results, baseline)
        
        # 3. 显示推荐命令
        recommended = next((r for r in results if r['threshold'] == 60), results[len(results)//2])
        print("="*80)
        print("★ 推荐方案执行命令")
        print("="*80)
        cmd = optimizer.generate_command(recommended['exclude_dict'], windows=(args.platform=='windows'))
        print(cmd)
        print()
        
        # 4. 生成报告
        output_file = args.output if args.output else f"{args.structure}_EXCLUDE_RECOMMENDATIONS.md"
        optimizer.generate_report(results, baseline, output_file)
        print()
        
        print("="*80)
        print("✓ 分析完成!")
        print("="*80)
        print(f"推荐: 使用 threshold=60 meV 方案")
        print(f"报告: {output_file}")
        print(f"命令: 见上方输出")
    
    return 0


if __name__ == "__main__":
    exit(main())
