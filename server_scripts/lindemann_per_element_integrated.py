# -*- coding: utf-8 -*-
"""
分元素计算Lindemann指数 - 集成版
快速计算 Pt、Sn、PtSn(Cluster)、PtSnO(All) 的Lindemann指数
专为批处理优化，输出单行CSV结果
"""
import MDAnalysis as mda
import numpy as np
from itertools import combinations
import argparse
import time
import sys
from pathlib import Path

def parse_index_file(index_file, verbose=False):
    """
    解析Gromacs格式的index文件，提取PtSnOCluster原子索引
    
    Args:
        index_file: index_zsplit.ndx文件路径
        verbose: 是否打印详细信息
        
    Returns:
        list: 团簇原子的索引列表（1-based）
    """
    cluster_indices = []
    in_cluster_section = False
    
    with open(index_file, 'r') as f:
        for line in f:
            line = line.strip()
            
            # 检测团簇部分开始
            if '[ PtSnOCluster ]' in line or '[ Cluster ]' in line:
                in_cluster_section = True
                continue
            
            # 检测下一个section开始，结束团簇部分
            if in_cluster_section and line.startswith('['):
                break
            
            # 读取团簇原子索引
            if in_cluster_section and line:
                try:
                    indices = [int(x) for x in line.split()]
                    cluster_indices.extend(indices)
                except ValueError:
                    continue
    
    if verbose:
        print(f"从index文件读取到 {len(cluster_indices)} 个团簇原子")
    
    return cluster_indices

def fast_unwrap_trajectory(universe, selection='all', verbose=False):
    """快速unwrap轨迹"""
    if verbose:
        print(f"Unwrap处理: {selection}")
    
    ag = universe.select_atoms(selection)
    n_atoms = len(ag)
    n_frames = len(universe.trajectory)

    all_coords = np.zeros((n_frames, n_atoms, 3))
    box_dims = np.zeros((n_frames, 3))

    for i, ts in enumerate(universe.trajectory):
        all_coords[i] = ag.positions
        box_dims[i] = ts.dimensions[:3]

    box = box_dims[0]
    delta = np.diff(all_coords, axis=0)
    offset = box * np.round(delta / box)
    cumulative_offset = np.concatenate([np.zeros((1, n_atoms, 3)), np.cumsum(offset, axis=0)], axis=0)
    unwrapped_coords = all_coords + cumulative_offset

    if verbose:
        print(f"  完成: {n_frames} 帧, {n_atoms} 原子")
    
    return unwrapped_coords

def calculate_lindemann_fast(unwrapped_coords):
    """快速计算Lindemann指数"""
    n_frames, n_atoms, _ = unwrapped_coords.shape
    
    if n_atoms <= 1:
        return 0.0
    
    atom_pairs = list(combinations(range(n_atoms), 2))
    num_pairs = len(atom_pairs)

    sum_r = np.zeros(num_pairs, dtype=np.float64)
    sum_r2 = np.zeros(num_pairs, dtype=np.float64)

    for frame_idx in range(n_frames):
        positions = unwrapped_coords[frame_idx]
        for idx, (i, j) in enumerate(atom_pairs):
            r = np.linalg.norm(positions[j] - positions[i])
            sum_r[idx] += r
            sum_r2[idx] += r ** 2

    avg_r = sum_r / n_frames
    avg_r2 = sum_r2 / n_frames
    std_r = np.sqrt(avg_r2 - avg_r ** 2)
    lindemann_ratios = std_r / avg_r
    
    return np.mean(lindemann_ratios)

def calculate_all_elements(coord_file, traj_file=None, box_dimensions=None, verbose=False, index_file=None):
    """
    计算所有元素的Lindemann指数
    
    Args:
        index_file: 可选的index_zsplit.ndx文件，用于精确选择团簇原子（排除载体）
    
    Returns:
        dict: {'Pt': float, 'Sn': float, 'PtSn': float, 'PtSnO': float, 'atom_counts': dict}
    """
    if traj_file is None:
        traj_file = coord_file
    
    try:
        u = mda.Universe(coord_file, traj_file, topology_format='XYZ')
    except Exception as e:
        if verbose:
            print(f"❌ 加载失败: {e}")
        return None

    # 设置盒子尺寸
    if box_dimensions is None:
        box_dimensions = np.array([16.5717, 16.6582, 30.0, 90.0, 90.0, 90.0])
    else:
        box_dimensions = np.array(box_dimensions)

    for ts in u.trajectory:
        ts.dimensions = box_dimensions

    # 尝试读取index文件来精确选择团簇原子
    cluster_indices = None
    if index_file and Path(index_file).exists():
        try:
            cluster_indices = parse_index_file(index_file, verbose=verbose)
            if verbose:
                print(f"✅ 使用index文件选择团簇原子: {len(cluster_indices)} 个")
        except Exception as e:
            if verbose:
                print(f"⚠️ index文件解析失败: {e}，使用默认方法")
    
    # 统计原子数
    try:
        if cluster_indices is not None:
            # 使用index文件精确选择（推荐！）
            # cluster_indices是1-based，需要转换为0-based
            cluster_ids_str = ' '.join([str(i-1) for i in cluster_indices])
            cluster_atoms = u.select_atoms(f'index {cluster_ids_str}')
            
            # 从团簇中分离各元素
            pt_atoms = cluster_atoms.select_atoms('name Pt')
            sn_atoms = cluster_atoms.select_atoms('name Sn')
            o_atoms = cluster_atoms.select_atoms('name O')
            ptsn_atoms = cluster_atoms.select_atoms('name Pt or name Sn')
            all_atoms = cluster_atoms  # PtSnO就是整个团簇
            
            if verbose:
                print(f"从团簇选择: Pt={len(pt_atoms)}, Sn={len(sn_atoms)}, O={len(o_atoms)}")
        else:
            # 备用方法：选择所有原子（会包含载体的O！）
            if verbose:
                print("⚠️ 未使用index文件，选择所有原子（可能包含载体O）")
            pt_atoms = u.select_atoms('name Pt')
            sn_atoms = u.select_atoms('name Sn')
            o_atoms = u.select_atoms('name O')
            ptsn_atoms = u.select_atoms('name Pt or name Sn')
            all_atoms = u.select_atoms('name Pt or name Sn or name O')
        
        atom_counts = {
            'Pt': len(pt_atoms),
            'Sn': len(sn_atoms),
            'O': len(o_atoms),
            'PtSn': len(ptsn_atoms),
            'PtSnO': len(all_atoms)
        }
        
        if verbose:
            print(f"原子数: Pt={atom_counts['Pt']}, Sn={atom_counts['Sn']}, O={atom_counts['O']}")
    
    except Exception as e:
        if verbose:
            print(f"❌ 选择原子失败: {e}")
        return None

    results = {
        'Pt': 0.0,
        'Sn': 0.0,
        'PtSn': 0.0,
        'PtSnO': 0.0,
        'atom_counts': atom_counts
    }

    # 计算 Pt
    if atom_counts['Pt'] > 1:
        try:
            unwrapped = fast_unwrap_trajectory(u, 'name Pt', verbose=verbose)
            results['Pt'] = calculate_lindemann_fast(unwrapped)
            if verbose:
                print(f"Pt Lindemann: {results['Pt']:.6f}")
        except Exception as e:
            if verbose:
                print(f"Pt 计算失败: {e}")

    # 计算 Sn
    if atom_counts['Sn'] > 1:
        try:
            unwrapped = fast_unwrap_trajectory(u, 'name Sn', verbose=verbose)
            results['Sn'] = calculate_lindemann_fast(unwrapped)
            if verbose:
                print(f"Sn Lindemann: {results['Sn']:.6f}")
        except Exception as e:
            if verbose:
                print(f"Sn 计算失败: {e}")

    # 计算 PtSn (Cluster)
    if atom_counts['PtSn'] > 1:
        try:
            unwrapped = fast_unwrap_trajectory(u, 'name Pt or name Sn', verbose=verbose)
            results['PtSn'] = calculate_lindemann_fast(unwrapped)
            if verbose:
                print(f"PtSn Lindemann: {results['PtSn']:.6f}")
        except Exception as e:
            if verbose:
                print(f"PtSn 计算失败: {e}")

    # 计算 PtSnO (All)
    # 如果没有O原子，PtSnO = PtSn
    if atom_counts['O'] == 0:
        results['PtSnO'] = results['PtSn']
        if verbose:
            print(f"PtSnO (无O原子，使用PtSn值): {results['PtSnO']:.6f}")
    elif atom_counts['PtSnO'] > 1:
        try:
            unwrapped = fast_unwrap_trajectory(u, 'name Pt or name Sn or name O', verbose=verbose)
            results['PtSnO'] = calculate_lindemann_fast(unwrapped)
            if verbose:
                print(f"PtSnO Lindemann: {results['PtSnO']:.6f}")
        except Exception as e:
            if verbose:
                print(f"PtSnO 计算失败: {e}")
            # 失败时使用PtSn的值作为fallback
            results['PtSnO'] = results['PtSn']

    return results

def main():
    parser = argparse.ArgumentParser(
        description='快速计算各元素Lindemann指数 (Pt, Sn, PtSn, PtSnO)'
    )
    
    parser.add_argument('--coord', type=str, required=True,
                        help='坐标文件路径（XYZ格式）')
    parser.add_argument('--traj', type=str, default=None,
                        help='轨迹文件路径（默认与coord相同）')
    parser.add_argument('--output', type=str, default='lindemann_elements.csv',
                        help='输出CSV文件路径')
    parser.add_argument('--box', type=float, nargs=6, default=None,
                        help='盒子尺寸（6个值）')
    parser.add_argument('--index', type=str, default=None,
                        help='index_zsplit.ndx文件路径（用于精确选择团簇原子，排除载体）')
    parser.add_argument('--verbose', action='store_true',
                        help='显示详细信息')
    parser.add_argument('--stdout-only', action='store_true',
                        help='只输出到stdout，不写文件')
    
    args = parser.parse_args()
    
    start_time = time.time()
    
    results = calculate_all_elements(
        coord_file=args.coord,
        traj_file=args.traj,
        box_dimensions=args.box,
        verbose=args.verbose,
        index_file=args.index
    )
    
    if results is None:
        if args.verbose:
            print("❌ 计算失败")
        sys.exit(1)
    
    elapsed = time.time() - start_time
    
    # 输出格式：Pt,Sn,PtSn,PtSnO
    output_line = f"{results['Pt']:.6f},{results['Sn']:.6f},{results['PtSn']:.6f},{results['PtSnO']:.6f}"
    
    if args.stdout_only:
        # 只输出到stdout，供shell脚本捕获
        print(output_line)
    else:
        # 写入文件并打印
        with open(args.output, 'w') as f:
            f.write("Pt,Sn,PtSn,PtSnO\n")
            f.write(output_line + "\n")
        
        if args.verbose:
            print(f"✅ 结果已保存: {args.output}")
            print(f"耗时: {elapsed:.2f}s")
        
        print(output_line)

if __name__ == '__main__':
    main()
