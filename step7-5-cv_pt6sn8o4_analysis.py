#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pt6Sn8O4体系热容模拟分析 - Cv-1到Cv-5

专门用于分析Pt6Sn8O4氧化物体系的5次热容模拟
温度范围: 200K-1100K (间隔50K, 19个温度点)
重复次数: 5次 (Cv-1, Cv-2, Cv-3, Cv-4, Cv-5)

用法:
    python step7-5-cv_pt6sn8o4_analysis.py
    python step7-5-cv_pt6sn8o4_analysis.py --enable-msd-filter

Author: AI Assistant
Date: 2025-10-27
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager
from scipy.stats import pearsonr
import os
import sys
import warnings
import argparse
import re
from pathlib import Path
import seaborn as sns

# 设置控制台输出编码 - 必须在最早导入后立即设置
if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8')
    # 设置默认编码
    import locale
    try:
        locale.setlocale(locale.LC_ALL, 'zh_CN.UTF-8')
    except:
        pass

warnings.filterwarnings('ignore')

# 配置中文字体
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'Arial']
plt.rcParams['axes.unicode_minus'] = False

def extract_path_signature(filepath, is_msd_path=True):
    """
    从文件路径提取4级路径签名
    
    Args:
        filepath: 完整文件路径
        is_msd_path: True=MSD路径(有温度目录), False=能量路径(无温度目录)
    
    Returns:
        path_signature: 4级路径签名
    """
    if not filepath or pd.isna(filepath):
        return None
    
    filepath = str(filepath).replace('\\', '/')
    
    # 1. 提取run信息 (T1000.r24.gpu0)
    run_match = re.search(r'(T\d+\.r\d+\.gpu\d+)', filepath, re.IGNORECASE)
    if not run_match:
        return None
    run_info = run_match.group(1).lower()
    
    # 2. 分割路径
    parts = re.split(r'[\\/]', filepath)
    
    # 3. 找到关键目录的索引
    if is_msd_path:
        # MSD路径: 找温度目录 (1000K)
        key_idx = None
        for i, part in enumerate(parts):
            if re.match(r'\d+K$', part, re.IGNORECASE):
                key_idx = i
                break
    else:
        # 能量/Lindemann路径: 找run所在位置
        key_idx = None
        for i, part in enumerate(parts):
            if re.search(r'T\d+\.r\d+\.gpu\d+', part, re.IGNORECASE):
                key_idx = i
                break
    
    if key_idx is None or key_idx < 2:
        # 无法提取足够的层级,返回简化签名
        return run_info
    
    # 4. 提取目录层级
    composition_dir = parts[key_idx - 1].lower()  # Cv-1
    parent_dir = parts[key_idx - 2].lower()       # g-1535-sn8pt6o4
    
    # 5. 检查批次标识符
    batch_keywords = ['run3', 'run2', 'run4', 'run5']
    path_signature = f"{parent_dir}/{composition_dir}/{run_info}"
    
    # 向上搜索批次标识符
    if key_idx >= 3:
        for check_idx in range(key_idx - 3, max(-1, key_idx - 6), -1):
            if check_idx < 0 or check_idx >= len(parts):
                break
            check_dir = parts[check_idx].lower()
            if check_dir in batch_keywords:
                path_signature = f"{check_dir}/{parent_dir}/{composition_dir}/{run_info}"
                break
    
    return path_signature


class Pt6Sn8O4Analyzer:
    """Pt6Sn8O4氧化物体系分析器"""
    
    def __init__(self, base_path, output_dir, enable_msd_filter=False):
        """
        初始化分析器
        
        Args:
            base_path: v626数据根目录
            output_dir: 输出目录
            enable_msd_filter: 是否启用MSD异常值筛选
        """
        self.base_path = Path(base_path)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.enable_msd_filter = enable_msd_filter
        
        # 体系路径
        self.system_path = self.base_path / 'dp-md' / '4090-ustc' / 'o68' / 'g-1535-Sn8Pt6O4'
        
        # 温度列表 (200K-1100K, 间隔50K)
        self.temperatures = [f'{t}K' for t in range(200, 1150, 50)]
        
        # Cv模拟次数
        self.cv_runs = ['Cv-1', 'Cv-2', 'Cv-3', 'Cv-4', 'Cv-5']
        
        # 原子数
        self.n_pt = 6
        self.n_sn = 8
        self.n_o = 4
        self.total_atoms = self.n_pt + self.n_sn + self.n_o  # 18原子
        
        # 结果存储
        self.results = {}
        
        # MSD过滤统计
        self.msd_filter_stats = {
            'total_runs': 0,
            'filtered_runs': 0
        }
        
        # 加载MSD异常路径签名（如果启用）
        if enable_msd_filter:
            self.msd_outliers = self.load_msd_outliers()
        else:
            self.msd_outliers = set()
        
        print(f"\n{'='*100}")
        print(f"初始化Pt6Sn8O4分析器")
        print(f"体系路径: {self.system_path}")
        print(f"输出目录: {self.output_dir}")
        print(f"温度点数: {len(self.temperatures)} (200K-1100K, 间隔50K)")
        print(f"重复次数: {len(self.cv_runs)} (Cv-1到Cv-5)")
        print(f"原子组成: Pt{self.n_pt}Sn{self.n_sn}O{self.n_o} ({self.total_atoms}原子)")
        if self.enable_msd_filter:
            print(f"✅ MSD过滤: 已启用 ({len(self.msd_outliers)}个异常路径签名)")
        else:
            print(f"⚠️ MSD过滤: 未启用，使用全部数据")
        print(f"{'='*100}")
    
    def load_msd_outliers(self):
        """加载MSD异常路径签名"""
        script_dir = Path(__file__).parent
        outliers_file = script_dir.parent / 'results' / 'large_D_outliers.csv'
        
        if not outliers_file.exists():
            print(f"\n⚠️ 未找到MSD异常文件: {outliers_file}")
            print(f"   将使用全部数据，不进行筛选")
            return set()
        
        try:
            df_outliers = pd.read_csv(outliers_file, encoding='utf-8')
            print(f"\n>>> 加载MSD异常数据:")
            print(f"    异常记录数: {len(df_outliers)}")
            
            # 统计异常原因
            if 'reason' in df_outliers.columns:
                for reason, count in df_outliers['reason'].value_counts().items():
                    pct = count / len(df_outliers) * 100
                    print(f"      - {reason}: {count} ({pct:.1f}%)")
            
            # 提取路径签名
            filter_signatures = set()
            for _, row in df_outliers.iterrows():
                filepath = row.get('filepath', '')
                if filepath:
                    sig = extract_path_signature(filepath, is_msd_path=True)
                    if sig:
                        filter_signatures.add(sig)
            
            print(f"    唯一路径签名: {len(filter_signatures)}")
            return filter_signatures
            
        except Exception as e:
            print(f"\n⚠️ 加载MSD异常文件出错: {e}")
            print(f"   将使用全部数据，不进行筛选")
            return set()
    
    def load_data_for_temp(self, cv_run, temp):
        """
        加载指定Cv运行和温度的数据
        
        Args:
            cv_run: Cv运行名称 (e.g., 'Cv-1')
            temp: 温度 (e.g., '300K')
        
        Returns:
            dict: 数据字典，如果失败返回None
        """
        # 构建路径
        temp_value = temp.replace('K', '')
        cv_path = self.system_path / cv_run
        
        # 查找匹配的温度目录
        temp_dirs = sorted(cv_path.glob(f"T{temp_value}.*"))
        
        if not temp_dirs:
            return None
        
        temp_dir = temp_dirs[0]  # 每个温度只有一个目录
        
        # MSD筛选检查
        if self.enable_msd_filter:
            path_signature = extract_path_signature(str(temp_dir), is_msd_path=False)
            if path_signature and path_signature in self.msd_outliers:
                self.msd_filter_stats['filtered_runs'] += 1
                return None
            self.msd_filter_stats['total_runs'] += 1
        
        # 读取数据文件
        try:
            # 1. coordination_time_series.csv
            cn_file = temp_dir / 'coordination_time_series.csv'
            if not cn_file.exists():
                return None
            df_cn = pd.read_csv(cn_file)
            
            # 2. cluster_global_q6_time_series.csv
            q6_file = temp_dir / 'cluster_global_q6_time_series.csv'
            if not q6_file.exists():
                return None
            df_q6 = pd.read_csv(q6_file)
            
            # 3. cluster_geometry_time_series.csv
            geo_file = temp_dir / 'cluster_geometry_time_series.csv'
            df_geo = pd.read_csv(geo_file) if geo_file.exists() else None
            
            # 4. element_comparison.csv
            elem_file = temp_dir / 'element_comparison.csv'
            df_elem = pd.read_csv(elem_file) if elem_file.exists() else None
            
            # 提取统计量
            data = {
                'temp': temp,
                'cv_run': cv_run,
                'q6': df_q6['cluster_metal_q6_global'].mean(),
                'q6_std': df_q6['cluster_metal_q6_global'].std(),
            }
            
            # Pt配位数
            data['pt_cn_total'] = df_cn['Pt_cn_total'].mean()
            data['pt_pt_bonds'] = df_cn['Pt_cn_Pt_Pt'].mean() if 'Pt_cn_Pt_Pt' in df_cn.columns else 0
            data['pt_sn_bonds'] = df_cn['Pt_cn_Pt_Sn'].mean() if 'Pt_cn_Pt_Sn' in df_cn.columns else 0
            data['pt_o_bonds'] = df_cn['Pt_cn_Pt_O'].mean() if 'Pt_cn_Pt_O' in df_cn.columns else 0
            
            # Sn配位数
            data['sn_cn_total'] = df_cn['Sn_cn_total'].mean() if 'Sn_cn_total' in df_cn.columns else 0
            data['sn_sn_bonds'] = df_cn['Sn_cn_Sn_Sn'].mean() if 'Sn_cn_Sn_Sn' in df_cn.columns else 0
            data['sn_pt_bonds'] = df_cn['Sn_cn_Sn_Pt'].mean() if 'Sn_cn_Sn_Pt' in df_cn.columns else 0
            data['sn_o_bonds'] = df_cn['Sn_cn_Sn_O'].mean() if 'Sn_cn_Sn_O' in df_cn.columns else 0
            
            # O配位数
            data['o_cn_total'] = df_cn['O_cn_total'].mean() if 'O_cn_total' in df_cn.columns else 0
            data['o_pt_bonds'] = df_cn['O_cn_O_Pt'].mean() if 'O_cn_O_Pt' in df_cn.columns else 0
            data['o_sn_bonds'] = df_cn['O_cn_O_Sn'].mean() if 'O_cn_O_Sn' in df_cn.columns else 0
            
            # 归一化键密度
            data['pt_pt_bonds_per_pt'] = data['pt_pt_bonds'] / self.n_pt
            data['pt_sn_bonds_per_pt'] = data['pt_sn_bonds'] / self.n_pt
            data['pt_o_bonds_per_pt'] = data['pt_o_bonds'] / self.n_pt
            
            data['sn_sn_bonds_per_sn'] = data['sn_sn_bonds'] / self.n_sn
            data['sn_pt_bonds_per_sn'] = data['sn_pt_bonds'] / self.n_sn
            data['sn_o_bonds_per_sn'] = data['sn_o_bonds'] / self.n_sn
            
            data['o_pt_bonds_per_o'] = data['o_pt_bonds'] / self.n_o
            data['o_sn_bonds_per_o'] = data['o_sn_bonds'] / self.n_o
            
            # 几何数据
            if df_geo is not None:
                data['rg'] = df_geo['gyration_radius'].mean()
                data['rg_std'] = df_geo['gyration_radius'].std()
                data['pt_dist'] = df_geo['pt_avg_dist_to_center'].mean() if 'pt_avg_dist_to_center' in df_geo.columns else 0
                data['sn_dist'] = df_geo['sn_avg_dist_to_center'].mean() if 'sn_avg_dist_to_center' in df_geo.columns else 0
                
                # O到质心距离
                if 'o_avg_dist_to_center' in df_geo.columns:
                    data['o_dist'] = df_geo['o_avg_dist_to_center'].mean()
                else:
                    data['o_dist'] = 0
            else:
                data['rg'] = 0
                data['rg_std'] = 0
                data['pt_dist'] = 0
                data['sn_dist'] = 0
                data['o_dist'] = 0
            
            # Q4数据
            if df_elem is not None and 'Q4' in df_elem.columns:
                pt_row = df_elem[df_elem['Element'] == 'Pt']
                sn_row = df_elem[df_elem['Element'] == 'Sn']
                o_row = df_elem[df_elem['Element'] == 'O']
                
                data['pt_q4'] = pt_row['Q4'].values[0] if not pt_row.empty else 0
                data['sn_q4'] = sn_row['Q4'].values[0] if not sn_row.empty else 0
                data['o_q4'] = o_row['Q4'].values[0] if not o_row.empty else 0
                
                # Q6数据
                data['pt_q6'] = pt_row['Q6'].values[0] if not pt_row.empty else 0
                data['sn_q6'] = sn_row['Q6'].values[0] if not sn_row.empty else 0
                data['o_q6'] = o_row['Q6'].values[0] if not o_row.empty else 0
            else:
                data['pt_q4'] = data['sn_q4'] = data['o_q4'] = 0
                data['pt_q6'] = data['sn_q6'] = data['o_q6'] = 0
            
            return data
            
        except Exception as e:
            print(f"  [错误] {cv_run} @ {temp}: {e}")
            return None
    
    def collect_all_data(self):
        """收集所有Cv运行的所有温度数据"""
        print(f"\n{'='*100}")
        print(f"开始收集Pt6Sn8O4数据")
        print(f"{'='*100}")
        
        all_data = []
        
        for temp in self.temperatures:
            print(f"\n处理温度: {temp}")
            temp_data = []
            
            for cv_run in self.cv_runs:
                data = self.load_data_for_temp(cv_run, temp)
                if data is not None:
                    temp_data.append(data)
                    all_data.append(data)
                    print(f"  ✓ {cv_run}: Q6={data['q6']:.3f}, Rg={data['rg']:.3f}Å")
                else:
                    print(f"  ✗ {cv_run}: 数据缺失或被过滤")
            
            # 计算该温度的平均值
            if temp_data:
                avg_data = {}
                for key in temp_data[0].keys():
                    if key in ['temp', 'cv_run']:
                        avg_data[key] = temp_data[0][key]
                    else:
                        values = [d[key] for d in temp_data]
                        if isinstance(values[0], (int, float)):
                            avg_data[key] = np.mean(values)
                        else:
                            avg_data[key] = values[0]
                
                avg_data['n_runs'] = len(temp_data)
                self.results[temp] = avg_data
                
                print(f"  → 平均 ({len(temp_data)}次): Q6={avg_data['q6']:.3f}, Rg={avg_data['rg']:.3f}Å")
        
        # 创建DataFrame
        self.df_all = pd.DataFrame(all_data)
        
        print(f"\n{'='*100}")
        print(f"数据收集完成!")
        print(f"{'='*100}")
        print(f"\n数据统计:")
        print(f"  预期数据点: {len(self.temperatures) * len(self.cv_runs)} ({len(self.temperatures)}温度 × {len(self.cv_runs)}次Cv)")
        print(f"  成功读取: {len(all_data)}")
        print(f"  失败: {len(self.temperatures) * len(self.cv_runs) - len(all_data)}")
        
        if self.enable_msd_filter and self.msd_filter_stats['filtered_runs'] > 0:
            total = self.msd_filter_stats['total_runs'] + self.msd_filter_stats['filtered_runs']
            filter_rate = self.msd_filter_stats['filtered_runs'] / total * 100
            print(f"\n[MSD筛选统计]:")
            print(f"  尝试读取: {total}")
            print(f"  保留: {self.msd_filter_stats['total_runs']} ({100-filter_rate:.1f}%)")
            print(f"  过滤: {self.msd_filter_stats['filtered_runs']} ({filter_rate:.1f}%)")
        
        print(f"{'='*100}")
    
    def save_data_table(self):
        """保存数据表"""
        # 保存原始数据
        csv_file = self.output_dir / 'pt6sn8o4_cv_all_data.csv'
        self.df_all.to_csv(csv_file, index=False, encoding='utf-8-sig')
        print(f"\n💾 原始数据已保存: {csv_file}")
        
        # 保存平均数据
        avg_rows = []
        for temp in sorted(self.results.keys(), key=lambda x: int(x.replace('K', ''))):
            data = self.results[temp]
            avg_rows.append({
                'temperature': temp,
                'temp_value': int(temp.replace('K', '')),
                'n_runs': data.get('n_runs', 0),
                'q6': data['q6'],
                'q6_std': data['q6_std'],
                'pt_q6': data.get('pt_q6', 0),
                'sn_q6': data.get('sn_q6', 0),
                'o_q6': data.get('o_q6', 0),
                'pt_pt_bonds_per_pt': data['pt_pt_bonds_per_pt'],
                'pt_sn_bonds_per_pt': data['pt_sn_bonds_per_pt'],
                'pt_o_bonds_per_pt': data['pt_o_bonds_per_pt'],
                'sn_pt_bonds_per_sn': data['sn_pt_bonds_per_sn'],
                'sn_o_bonds_per_sn': data['sn_o_bonds_per_sn'],
                'o_pt_bonds_per_o': data['o_pt_bonds_per_o'],
                'o_sn_bonds_per_o': data['o_sn_bonds_per_o'],
                'rg': data['rg'],
                'pt_dist': data['pt_dist'],
                'sn_dist': data['sn_dist'],
                'o_dist': data['o_dist']
            })
        
        df_avg = pd.DataFrame(avg_rows)
        avg_csv = self.output_dir / 'pt6sn8o4_cv_averaged_data.csv'
        df_avg.to_csv(avg_csv, index=False, encoding='utf-8-sig')
        print(f"💾 平均数据已保存: {avg_csv}")
    
    def plot_temperature_trends(self):
        """绘制温度趋势图"""
        fig, axes = plt.subplots(3, 3, figsize=(18, 15))
        fig.suptitle('Pt6Sn8O4 温度效应分析 (Cv-1到Cv-5平均)', fontsize=16, fontweight='bold')
        
        # 提取温度值
        temps = sorted([int(t.replace('K', '')) for t in self.results.keys()])
        
        # 9个物理量
        properties = [
            ('q6', 'Q6 六次对称性', 'blue'),
            ('rg', '回转半径 Rg (Å)', 'green'),
            ('pt_pt_bonds_per_pt', 'Pt-Pt键/Pt', 'red'),
            ('pt_sn_bonds_per_pt', 'Pt-Sn键/Pt', 'orange'),
            ('pt_o_bonds_per_pt', 'Pt-O键/Pt', 'purple'),
            ('sn_pt_bonds_per_sn', 'Sn-Pt键/Sn', 'brown'),
            ('sn_o_bonds_per_sn', 'Sn-O键/Sn', 'pink'),
            ('o_pt_bonds_per_o', 'O-Pt键/O', 'cyan'),
            ('o_sn_bonds_per_o', 'O-Sn键/O', 'magenta')
        ]
        
        for idx, (prop, label, color) in enumerate(properties):
            ax = axes[idx // 3, idx % 3]
            
            values = [self.results[f'{t}K'][prop] for t in temps]
            
            ax.plot(temps, values, 'o-', color=color, linewidth=2, markersize=6)
            ax.set_xlabel('温度 (K)', fontsize=10)
            ax.set_ylabel(label, fontsize=10)
            ax.set_title(label, fontsize=11, fontweight='bold')
            ax.grid(True, alpha=0.3)
            ax.tick_params(labelsize=9)
        
        plt.tight_layout()
        
        output_file = self.output_dir / 'pt6sn8o4_temperature_trends.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"📊 温度趋势图已保存: {output_file}")
    
    def plot_cv_comparison(self):
        """绘制不同Cv运行的对比"""
        # 选择几个关键温度点
        key_temps = ['300K', '600K', '900K']
        
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))
        fig.suptitle('Pt6Sn8O4 不同Cv运行对比 (关键温度)', fontsize=14, fontweight='bold')
        
        for idx, temp in enumerate(key_temps):
            ax = axes[idx]
            
            # 提取该温度的所有Cv数据
            temp_df = self.df_all[self.df_all['temp'] == temp]
            
            if len(temp_df) == 0:
                continue
            
            cv_labels = temp_df['cv_run'].values
            q6_values = temp_df['q6'].values
            
            colors = plt.cm.Set3(range(len(cv_labels)))
            bars = ax.bar(cv_labels, q6_values, color=colors, edgecolor='black', linewidth=1.5)
            
            # 添加数值标签
            for bar, val in zip(bars, q6_values):
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2., height,
                       f'{val:.3f}', ha='center', va='bottom', fontsize=9)
            
            ax.set_xlabel('Cv运行', fontsize=10)
            ax.set_ylabel('Q6', fontsize=10)
            ax.set_title(f'{temp}', fontsize=11, fontweight='bold')
            ax.grid(True, alpha=0.3, axis='y')
            ax.tick_params(labelsize=9)
            
            # 添加平均线
            avg_q6 = temp_df['q6'].mean()
            ax.axhline(y=avg_q6, color='red', linestyle='--', linewidth=2, label=f'平均: {avg_q6:.3f}')
            ax.legend(fontsize=9)
        
        plt.tight_layout()
        
        output_file = self.output_dir / 'pt6sn8o4_cv_comparison.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"📊 Cv对比图已保存: {output_file}")
    
    def run_analysis(self):
        """运行完整分析流程"""
        print(f"\n{'='*100}")
        print(f"开始Pt6Sn8O4完整分析")
        print(f"{'='*100}")
        
        # 1. 收集数据
        self.collect_all_data()
        
        # 2. 保存数据表
        self.save_data_table()
        
        # 3. 绘制图表
        self.plot_temperature_trends()
        self.plot_cv_comparison()
        
        print(f"\n{'='*100}")
        print(f"✅ Pt6Sn8O4分析完成!")
        print(f"{'='*100}")


def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='Pt6Sn8O4氧化物体系热容模拟分析')
    parser.add_argument('--enable-msd-filter', action='store_true',
                       help='启用Step 1的MSD异常值筛选')
    
    args = parser.parse_args()
    
    # 路径配置
    base_path = r"D:\OneDrive\py\Cv\lin\MSD_Analysis_Collection\v3_simplified_workflow\files\q6_cn\v626\coordination_time_series_results_sample_20251026_200908"
    output_dir = r"D:\OneDrive\py\Cv\lin\MSD_Analysis_Collection\v3_simplified_workflow\results\pt6sn8o4_cv"
    
    try:
        analyzer = Pt6Sn8O4Analyzer(
            base_path=base_path,
            output_dir=output_dir,
            enable_msd_filter=args.enable_msd_filter
        )
        
        analyzer.run_analysis()
            
    except Exception as e:
        print(f"\n❌ 分析失败: {e}")
        import traceback
        traceback.print_exc()


if __name__ == '__main__':
    main()
