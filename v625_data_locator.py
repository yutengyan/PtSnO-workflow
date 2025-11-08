"""
v625数据定位器 - 智能查找和加载分散的运行文件夹
支持新v625数据结构(2025-10-26)

新结构特点:
- 运行文件夹分散在多个子目录中
- Pt8: 8个运行 (4090-ustc/more/run3/, dp-md/4090-ustc/GPU-Pt8/, dp-md/4090-ustc/more/)
- Pt6: 3个运行 (dp-md/4090-ustc/more/)
- PtxSn8-x: 4个运行 (4090-ustc/more/run3/, dp-md/4090-ustc/more/)
"""
import os
import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Dict, Optional, Tuple
import warnings
warnings.filterwarnings('ignore')

class V625DataLocator:
    """v625数据定位和加载工具"""
    
    def __init__(self, base_path: str):
        """
        初始化数据定位器
        
        Args:
            base_path: v625数据根目录
                例如: .../coordination_time_series_results_sample_20251026_170451
        """
        self.base_path = Path(base_path)
        if not self.base_path.exists():
            raise ValueError(f"数据根目录不存在: {base_path}")
        
        # 缓存已找到的运行路径
        self._run_paths_cache = {}
        
        print(f"[V625DataLocator] 初始化完成: {self.base_path.name}")
    
    def find_all_runs(self, series_name: str) -> List[Path]:
        """
        递归查找指定系列的所有运行文件夹
        
        Args:
            series_name: 系列名称 ('Pt8', 'Pt6', 'PtxSn8-x')
        
        Returns:
            List[Path]: 运行文件夹的绝对路径列表
        
        Examples:
            >>> locator.find_all_runs('Pt8')
            [Path('.../dp-md/4090-ustc/more/Pt8'),
             Path('.../dp-md/4090-ustc/more/Pt8-2'),
             ...]
        """
        # 检查缓存
        if series_name in self._run_paths_cache:
            return self._run_paths_cache[series_name]
        
        run_paths = []
        
        # 递归搜索所有目录
        for root, dirs, files in os.walk(self.base_path):
            for dir_name in dirs:
                # 匹配模式: Pt8, Pt8-2, Pt8-3等 或 PtxSn8-x, PtxSn8-x-2等
                if dir_name == series_name or dir_name.startswith(f"{series_name}-"):
                    full_path = Path(root) / dir_name
                    run_paths.append(full_path)
        
        # 排序: 主文件夹在前,然后按数字排序
        def sort_key(path):
            name = path.name
            if name == series_name:
                return (0, 0)
            elif '-' in name:
                try:
                    num = int(name.split('-')[-1])
                    return (1, num)
                except:
                    return (2, name)
            return (3, name)
        
        run_paths = sorted(run_paths, key=sort_key)
        
        # 缓存结果
        self._run_paths_cache[series_name] = run_paths
        
        print(f"[{series_name}] 找到{len(run_paths)}个运行文件夹:")
        for i, path in enumerate(run_paths, 1):
            rel_path = path.relative_to(self.base_path)
            print(f"  {i}. {rel_path}")
        
        return run_paths
    
    def load_coordination_data(self, run_path: Path, system_name: str, 
                               temperature: str, load_all_runs: bool = True) -> Optional[pd.DataFrame]:
        """
        加载coordination数据
        支持v625格式(200K)和v626格式(T200.r*.gpu*)
        
        Args:
            run_path: 运行文件夹路径
            system_name: 体系名称 (如'pt8sn7-1-best')
            temperature: 温度 (如'200K')
            load_all_runs: True=加载并平均所有重复运行, False=只取第一个
        
        Returns:
            DataFrame或None: 成功返回平均数据,失败返回None
        """
        sys_path = run_path / system_name
        
        # v625格式: 200K/coordination_time_series.csv
        csv_path_v625 = sys_path / temperature / 'coordination_time_series.csv'
        if csv_path_v625.exists():
            try:
                return pd.read_csv(csv_path_v625)
            except Exception as e:
                print(f"    ⚠️ 读取失败: {csv_path_v625.name} - {e}")
                return None
        
        # v626格式: T200.r*.gpu*/coordination_time_series.csv
        temp_value = temperature.replace('K', '')
        temp_pattern = f"T{temp_value}.*"
        matching_dirs = sorted(sys_path.glob(temp_pattern))
        
        if matching_dirs:
            if not load_all_runs:
                # 只取第一个
                csv_path = matching_dirs[0] / 'coordination_time_series.csv'
                if csv_path.exists():
                    try:
                        return pd.read_csv(csv_path)
                    except Exception as e:
                        print(f"    ⚠️ 读取失败: {csv_path.name} - {e}")
                        return None
            else:
                # 加载所有重复运行并平均
                all_dfs = []
                for temp_dir in matching_dirs:
                    csv_path = temp_dir / 'coordination_time_series.csv'
                    if csv_path.exists():
                        try:
                            df = pd.read_csv(csv_path)
                            all_dfs.append(df)
                        except Exception as e:
                            continue
                
                if all_dfs:
                    if len(all_dfs) == 1:
                        return all_dfs[0]
                    else:
                        # 按frame分组求平均
                        combined = pd.concat(all_dfs)
                        avg_df = combined.groupby('frame').mean().reset_index()
                        return avg_df
        
        return None
    
    def load_q6_data(self, run_path: Path, system_name: str, 
                     temperature: str, load_all_runs: bool = True) -> Optional[pd.DataFrame]:
        """
        加载Q6数据
        支持v625格式(200K)和v626格式(T200.r*.gpu*)
        
        Args:
            run_path: 运行文件夹路径
            system_name: 体系名称
            temperature: 温度
            load_all_runs: True=加载并平均所有重复运行, False=只取第一个
        
        Returns:
            DataFrame或None: 成功返回平均数据,失败返回None
        """
        sys_path = run_path / system_name
        
        # v625格式: 200K/cluster_global_q6_time_series.csv
        csv_path_v625 = sys_path / temperature / 'cluster_global_q6_time_series.csv'
        if csv_path_v625.exists():
            try:
                return pd.read_csv(csv_path_v625)
            except Exception as e:
                print(f"    ⚠️ 读取Q6失败: {csv_path_v625.name} - {e}")
                return None
        
        # v626格式: T200.r*.gpu*/cluster_global_q6_time_series.csv
        temp_value = temperature.replace('K', '')
        temp_pattern = f"T{temp_value}.*"
        matching_dirs = sorted(sys_path.glob(temp_pattern))
        
        if matching_dirs:
            if not load_all_runs:
                # 只取第一个
                csv_path = matching_dirs[0] / 'cluster_global_q6_time_series.csv'
                if csv_path.exists():
                    try:
                        return pd.read_csv(csv_path)
                    except Exception as e:
                        print(f"    ⚠️ 读取Q6失败: {csv_path.name} - {e}")
                        return None
            else:
                # 加载所有重复运行并平均
                all_dfs = []
                for temp_dir in matching_dirs:
                    csv_path = temp_dir / 'cluster_global_q6_time_series.csv'
                    if csv_path.exists():
                        try:
                            df = pd.read_csv(csv_path)
                            all_dfs.append(df)
                        except Exception as e:
                            continue
                
                if all_dfs:
                    if len(all_dfs) == 1:
                        return all_dfs[0]
                    else:
                        # 按frame分组求平均
                        combined = pd.concat(all_dfs)
                        avg_df = combined.groupby('frame').mean().reset_index()
                        return avg_df
        
        return None
    
    def load_geometry_data(self, run_path: Path, system_name: str, 
                          temperature: str, load_all_runs: bool = True) -> Optional[pd.DataFrame]:
        """
        加载几何数据(回转半径、质心距离等)
        支持v625格式(200K)和v626格式(T200.r*.gpu*)
        
        Args:
            run_path: 运行文件夹路径
            system_name: 体系名称
            temperature: 温度
            load_all_runs: True=加载并平均所有重复运行, False=只取第一个
        
        Returns:
            DataFrame或None: 包含sn_avg_dist_to_center, pt_avg_dist_to_center, gyration_radius
        """
        sys_path = run_path / system_name
        
        # v625格式: 200K/cluster_geometry_time_series.csv
        csv_path_v625 = sys_path / temperature / 'cluster_geometry_time_series.csv'
        if csv_path_v625.exists():
            try:
                return pd.read_csv(csv_path_v625)
            except Exception as e:
                print(f"    ⚠️ 读取几何数据失败: {csv_path_v625.name} - {e}")
                return None
        
        # v626格式: T200.r*.gpu*/cluster_geometry_time_series.csv
        temp_value = temperature.replace('K', '')
        temp_pattern = f"T{temp_value}.*"
        matching_dirs = sorted(sys_path.glob(temp_pattern))
        
        if matching_dirs:
            if not load_all_runs:
                # 只取第一个
                csv_path = matching_dirs[0] / 'cluster_geometry_time_series.csv'
                if csv_path.exists():
                    try:
                        return pd.read_csv(csv_path)
                    except Exception as e:
                        print(f"    ⚠️ 读取几何数据失败: {csv_path.name} - {e}")
                        return None
            else:
                # 加载所有重复运行并平均
                all_dfs = []
                for temp_dir in matching_dirs:
                    csv_path = temp_dir / 'cluster_geometry_time_series.csv'
                    if csv_path.exists():
                        try:
                            df = pd.read_csv(csv_path)
                            all_dfs.append(df)
                        except Exception as e:
                            continue
                
                if all_dfs:
                    if len(all_dfs) == 1:
                        return all_dfs[0]
                    else:
                        # 按frame分组求平均
                        combined = pd.concat(all_dfs)
                        avg_df = combined.groupby('frame').mean().reset_index()
                        return avg_df
        
        return None
    
    def load_multi_run_average(self, series_name: str, system_name: str, 
                               temperature: str, data_type: str = 'coordination'
                               ) -> Optional[Dict]:
        """
        加载多个运行并求平均
        
        Args:
            series_name: 系列名称 ('Pt8', 'Pt6', 'PtxSn8-x')
            system_name: 体系名称 (如'pt8sn7-1-best')
            temperature: 温度 (如'200K')
            data_type: 数据类型 ('coordination' 或 'q6')
        
        Returns:
            dict: {
                'mean': 平均值字典,
                'std': 标准差字典,
                'n_runs': 成功运行数,
                'failed_runs': 失败运行列表
            }
            或None: 所有运行都失败
        """
        run_paths = self.find_all_runs(series_name)
        
        if not run_paths:
            return None
        
        # 选择加载函数
        load_func = self.load_coordination_data if data_type == 'coordination' else self.load_q6_data
        
        # 加载所有运行的数据
        all_data = []
        failed_runs = []
        
        for run_path in run_paths:
            df = load_func(run_path, system_name, temperature)
            if df is not None:
                # 时间平均
                time_avg = df.mean(numeric_only=True)
                all_data.append(time_avg)
            else:
                failed_runs.append(run_path.name)
        
        if not all_data:
            return None
        
        # 集合平均(多次运行平均)
        all_df = pd.DataFrame(all_data)
        mean_values = all_df.mean()
        std_values = all_df.std()
        
        return {
            'mean': mean_values.to_dict(),
            'std': std_values.to_dict(),
            'n_runs': len(all_data),
            'failed_runs': failed_runs
        }
    
    def get_field_name(self, base_name: str, check_df: Optional[pd.DataFrame] = None) -> str:
        """
        获取正确的字段名(处理大小写不一致)
        
        Args:
            base_name: 基础字段名 (如'pt_cn_total' 或 'Pt_cn_total')
            check_df: 可选的DataFrame用于验证字段是否存在
        
        Returns:
            str: 正确的字段名
        """
        # v625统一使用大写开头
        if base_name.startswith('pt_'):
            correct_name = 'Pt_' + base_name[3:]
        elif base_name.startswith('sn_'):
            correct_name = 'Sn_' + base_name[3:]
        else:
            correct_name = base_name
        
        # 如果提供了DataFrame,验证字段是否存在
        if check_df is not None and correct_name not in check_df.columns:
            # 尝试小写版本
            if base_name in check_df.columns:
                return base_name
        
        return correct_name
    
    def extract_values(self, data_dict: Dict, field_mappings: Dict[str, str]) -> Dict:
        """
        从数据字典中提取指定字段的值
        
        Args:
            data_dict: load_multi_run_average返回的字典
            field_mappings: 字段映射 {'pt_cn_total': 'Pt_cn_total', ...}
        
        Returns:
            dict: 提取的值
        """
        if data_dict is None:
            return None
        
        result = {}
        mean_data = data_dict['mean']
        std_data = data_dict['std']
        
        for output_key, field_name in field_mappings.items():
            result[output_key] = mean_data.get(field_name, np.nan)
            result[f"{output_key}_std"] = std_data.get(field_name, np.nan)
        
        result['n_runs'] = data_dict['n_runs']
        
        return result


def test_locator():
    """测试数据定位器"""
    base_path = r"D:\OneDrive\py\Cv\lin\MSD_Analysis_Collection\v3_simplified_workflow\files\q6_cn\v625\coordination_time_series_results_sample_20251026_170451"
    
    locator = V625DataLocator(base_path)
    
    print("\n" + "="*60)
    print("测试1: 查找Pt8系列运行")
    print("="*60)
    pt8_runs = locator.find_all_runs('Pt8')
    print(f"✓ 找到{len(pt8_runs)}个Pt8运行\n")
    
    print("="*60)
    print("测试2: 查找Pt6系列运行")
    print("="*60)
    pt6_runs = locator.find_all_runs('Pt6')
    print(f"✓ 找到{len(pt6_runs)}个Pt6运行\n")
    
    print("="*60)
    print("测试3: 查找PtxSn8-x系列运行")
    print("="*60)
    ptxsn8x_runs = locator.find_all_runs('PtxSn8-x')
    print(f"✓ 找到{len(ptxsn8x_runs)}个PtxSn8-x运行\n")
    
    print("="*60)
    print("测试4: 加载Pt8的pt8sn7-1-best @ 200K数据(多运行平均)")
    print("="*60)
    result = locator.load_multi_run_average('Pt8', 'pt8sn7-1-best', '200K')
    if result:
        print(f"✓ 成功加载并平均")
        print(f"  - 成功运行数: {result['n_runs']}")
        print(f"  - 失败运行: {result['failed_runs'] if result['failed_runs'] else '无'}")
        print(f"  - Pt配位数均值: {result['mean'].get('Pt_cn_total', 'N/A'):.3f}")
        print(f"  - Pt配位数标准差: {result['std'].get('Pt_cn_total', 'N/A'):.3f}")
    else:
        print("✗ 加载失败")
    
    print("\n" + "="*60)
    print("测试完成!")
    print("="*60)


if __name__ == '__main__':
    test_locator()
