#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 6.3.1: Heat Capacity Scatter Plot - Cv1 vs Cv2

Creates a scatter plot comparing low-T and high-T heat capacity.
Labels use short format: (x,y) or (x,y,z) for Pt_x Sn_y O_z

Author: AI Assistant
Date: 2025-12-01

================================================================================
Usage (use vscode-1 environment):
================================================================================

# ★★★ Recommended commands ★★★

# With labels (recommended for publication)
& C:/Users/11207/.conda/envs/vscode-1/python.exe step6_3_1_cv_scatter_plot.py --no-errorbars --fontscale 1.5 --markerscale 2

# Scatter only (no labels, cleaner)
& C:/Users/11207/.conda/envs/vscode-1/python.exe step6_3_1_cv_scatter_plot.py --no-labels --no-errorbars --markerscale 2 --fontscale 1.5

# ★★★ All options ★★★

# Plot all data (default)
& C:/Users/11207/.conda/envs/vscode-1/python.exe step6_3_1_cv_scatter_plot.py

# Plot without error bars (cleaner, less overlap)
& C:/Users/11207/.conda/envs/vscode-1/python.exe step6_3_1_cv_scatter_plot.py --no-errorbars

# Plot without labels (scatter only)
& C:/Users/11207/.conda/envs/vscode-1/python.exe step6_3_1_cv_scatter_plot.py --no-labels

# Scale fonts and markers:
# --fontscale X   : scale all fonts (default 1.0, try 1.5 for larger)
# --markerscale X : scale marker sizes (default 1.0, try 2 for larger)

# Select which types to plot:
# --no-air        : exclude gas-phase
# --no-supported  : exclude supported Pt-Sn (no oxygen)
# --no-oxide      : exclude supported Pt-Sn-O

# Only plot gas-phase
& C:/Users/11207/.conda/envs/vscode-1/python.exe step6_3_1_cv_scatter_plot.py --no-supported --no-oxide

# Only plot supported Pt-Sn (no oxygen)
& C:/Users/11207/.conda/envs/vscode-1/python.exe step6_3_1_cv_scatter_plot.py --no-air --no-oxide

# Only plot supported Pt-Sn-O
& C:/Users/11207/.conda/envs/vscode-1/python.exe step6_3_1_cv_scatter_plot.py --no-air --no-supported

# Exclude specific compositions using --exclude "x,y,z;x,y,z":
# Exclude (3,5,3) and (3,4,1)
& C:/Users/11207/.conda/envs/vscode-1/python.exe step6_3_1_cv_scatter_plot.py --exclude "3,5,3;3,4,1"

# Exclude (6,8) - affects both air(6,8) and supported (6,8)
& C:/Users/11207/.conda/envs/vscode-1/python.exe step6_3_1_cv_scatter_plot.py --exclude "6,8"

# Combine filters: only oxide, no error bars, exclude (3,5,3)
& C:/Users/11207/.conda/envs/vscode-1/python.exe step6_3_1_cv_scatter_plot.py --no-air --no-supported --no-errorbars --exclude "3,5,3"

# Custom output path
& C:/Users/11207/.conda/envs/vscode-1/python.exe step6_3_1_cv_scatter_plot.py -o "results/my_plot.png"

================================================================================
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
from adjustText import adjust_text

# Font settings for English journals
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['font.size'] = 10


def parse_composition(name):
    """
    Parse structure name to extract Pt, Sn, O numbers
    
    Handles various formats:
    - Pt6sn8, Pt6sn8o4
    - O2pt4sn6, O2sn8pt7
    - Sn6pt5o2, Sn3o2pt2
    - Air68, Air86
    - Cv
    
    Returns: (pt_num, sn_num, o_num) or None
    """
    import re
    name_lower = name.lower()
    
    # Special case: Cv = Pt6Sn8O4
    if name == 'Cv':
        return (6, 8, 4)
    
    # Air series: Air68 = Pt6Sn8
    if 'air' in name_lower:
        match = re.search(r'air(\d+)(\d+)', name_lower)
        if match:
            return (int(match.group(1)), int(match.group(2)), 0)
        return None
    
    # Extract all components using flexible regex
    pt_match = re.search(r'pt(\d+)', name_lower)
    sn_match = re.search(r'sn(\d+)', name_lower)
    o_match = re.search(r'o(\d+)', name_lower)
    
    pt_num = int(pt_match.group(1)) if pt_match else 0
    sn_num = int(sn_match.group(1)) if sn_match else 0
    o_num = int(o_match.group(1)) if o_match else 0
    
    if pt_num > 0 or sn_num > 0:
        return (pt_num, sn_num, o_num)
    
    return None


def format_label_short(name):
    """
    Format structure name in short format for plot labels
    
    Rules:
    - Gas-phase: air(x,y) e.g., air(6,8)
    - Supported Pt-Sn: (x,y) e.g., (6,8)
    - Supported Pt-Sn-O: (x,y,z) e.g., (6,8,4)
    """
    comp = parse_composition(name)
    if comp is None:
        return name
    
    pt, sn, o = comp
    
    # Air series
    if 'air' in name.lower():
        return f'air({pt},{sn})'
    
    # With oxygen
    if o > 0:
        return f'({pt},{sn},{o})'
    else:
        return f'({pt},{sn})'


def classify_structure(name):
    """
    Classify structure type
    
    Returns: (type_key, display_name)
    """
    name_lower = name.lower()
    
    if 'air' in name_lower:
        return 'air', 'Gas-phase'
    
    # Use composition parser
    comp = parse_composition(name)
    if comp:
        pt, sn, o = comp
        if o > 0:
            return 'oxide', 'Supported Pt-Sn-O'
        else:
            return 'supported', 'Supported Pt-Sn'
    
    return 'other', 'Other'


def load_partition_data(structure, base_dir='results/step6_1_clustering'):
    """Load Cv1 and Cv2 from quality metrics file"""
    filepath = os.path.join(base_dir, f'{structure}_auto2_quality_metrics.csv')
    
    if not os.path.exists(filepath):
        return None
    
    try:
        df = pd.read_csv(filepath)
        
        if len(df) < 2:
            return None
        
        return {
            'cv1': df['Cv_cluster'].iloc[0],
            'err1': df['Cv_cluster_err'].iloc[0],
            'cv2': df['Cv_cluster'].iloc[1],
            'err2': df['Cv_cluster_err'].iloc[1],
        }
    except Exception as e:
        print(f"  Warning: Failed to read {structure}: {e}")
        return None


def find_all_structures(base_dir='results/step6_1_clustering'):
    """Find all available structures"""
    pattern = os.path.join(base_dir, '*_auto2_quality_metrics.csv')
    files = glob.glob(pattern)
    
    structures = []
    for f in files:
        basename = os.path.basename(f)
        structure = basename.replace('_auto2_quality_metrics.csv', '')
        structures.append(structure)
    
    return sorted(structures)


def should_exclude(name, exclude_list):
    """
    Check if a structure should be excluded based on composition
    
    Args:
        name: structure name
        exclude_list: list of tuples like [(3,5,3), (3,4,1)] to exclude
    
    Returns: True if should be excluded
    """
    if not exclude_list:
        return False
    
    comp = parse_composition(name)
    if comp is None:
        return False
    
    pt, sn, o = comp
    
    for exc in exclude_list:
        if len(exc) == 2:
            # (x, y) format - match Pt and Sn
            if pt == exc[0] and sn == exc[1]:
                return True
        elif len(exc) == 3:
            # (x, y, z) format - match Pt, Sn, O
            if pt == exc[0] and sn == exc[1] and o == exc[2]:
                return True
    
    return False


def create_cv_scatter_plot(output_path='results/step6_1_clustering/cv1_vs_cv2_scatter.png',
                           show_air=True, show_supported=True, show_oxide=True,
                           exclude_list=None, show_errorbars=True, show_labels=True,
                           fontscale=1.0, markerscale=1.0):
    """
    Create Cv1 vs Cv2 scatter plot with non-overlapping labels
    
    Args:
        output_path: output file path
        show_air: whether to show gas-phase data
        show_supported: whether to show supported Pt-Sn (no oxygen) data
        show_oxide: whether to show supported Pt-Sn-O data
        exclude_list: list of compositions to exclude, e.g., [(3,5,3), (3,4,1)]
        show_errorbars: whether to show error bars (default True)
        show_labels: whether to show data point labels (default True)
        fontscale: scale factor for all fonts (default 1.0)
        markerscale: scale factor for marker sizes (default 1.0)
    """
    
    # Base font sizes (will be scaled)
    FONT_LABEL = 8 * fontscale      # data point labels
    FONT_AXIS = 12 * fontscale      # axis labels
    FONT_TITLE = 13 * fontscale     # title
    FONT_LEGEND = 9 * fontscale     # legend
    FONT_TICK = 10 * fontscale      # tick labels
    STROKE_WIDTH = 2 * fontscale    # text stroke width
    
    # Base marker sizes (will be scaled)
    MARKER_AIR = 50 * markerscale
    MARKER_SUPPORTED = 40 * markerscale
    MARKER_OXIDE = 45 * markerscale
    
    base_dir = 'results/step6_1_clustering'
    structures = find_all_structures(base_dir)
    
    print(f"Found {len(structures)} systems")
    
    # Print filter settings
    print(f"\nFilter settings:")
    print(f"  Show Gas-phase: {show_air}")
    print(f"  Show Supported Pt-Sn: {show_supported}")
    print(f"  Show Supported Pt-Sn-O: {show_oxide}")
    print(f"  Show Error bars: {show_errorbars}")
    print(f"  Show Labels: {show_labels}")
    print(f"  Font scale: {fontscale}")
    print(f"  Marker scale: {markerscale}")
    if exclude_list:
        print(f"  Excluded compositions: {exclude_list}")
    
    # Load and classify data
    data_by_type = {
        'air': [],
        'supported': [],
        'oxide': [],
    }
    
    excluded_count = 0
    
    for struct in structures:
        # Check exclusion list
        if should_exclude(struct, exclude_list):
            comp = parse_composition(struct)
            print(f"    Excluding: {struct} -> {comp}")
            excluded_count += 1
            continue
        
        d = load_partition_data(struct, base_dir)
        if d is None:
            continue
        
        type_key, _ = classify_structure(struct)
        short_name = format_label_short(struct)
        
        data_by_type[type_key].append({
            'name': struct,
            'display': short_name,
            'cv1': d['cv1'],
            'cv2': d['cv2'],
            'err1': d['err1'],
            'err2': d['err2'],
        })
    
    if excluded_count > 0:
        print(f"  Total excluded: {excluded_count}")
    
    # Print statistics
    print(f"\nClassification statistics (after filtering):")
    print(f"  Gas-phase: {len(data_by_type['air'])}")
    print(f"  Supported Pt-Sn: {len(data_by_type['supported'])}")
    print(f"  Supported Pt-Sn-O: {len(data_by_type['oxide'])}")
    
    # Style settings (using scaled marker sizes)
    style = {
        'air': {'color': '#1E90FF', 'marker': 'o', 'label': 'Gas: air(x,y)', 's': MARKER_AIR},
        'supported': {'color': '#2E8B57', 'marker': 's', 'label': 'Supp.: (x,y)', 's': MARKER_SUPPORTED},
        'oxide': {'color': '#FF6347', 'marker': '^', 'label': 'Supp.-O: (x,y,z)', 's': MARKER_OXIDE},
    }
    
    # Determine which types to plot
    types_to_plot = []
    if show_air:
        types_to_plot.append('air')
    if show_supported:
        types_to_plot.append('supported')
    if show_oxide:
        types_to_plot.append('oxide')
    
    # Create figure - larger for better spacing
    fig, ax = plt.subplots(figsize=(24/2.54, 22/2.54))
    
    # Increase tick label font size (scaled)
    ax.tick_params(axis='both', labelsize=FONT_TICK)
    
    all_points = []
    
    # Plot data points
    for type_key in types_to_plot:
        data = data_by_type[type_key]
        if not data:
            continue
        
        s = style[type_key]
        
        cv1_vals = [d['cv1'] for d in data]
        cv2_vals = [d['cv2'] for d in data]
        err1_vals = [d['err1'] for d in data]
        err2_vals = [d['err2'] for d in data]
        
        # Error bars (optional)
        if show_errorbars:
            ax.errorbar(cv1_vals, cv2_vals, xerr=err1_vals, yerr=err2_vals,
                        fmt='none', ecolor=s['color'], alpha=0.3, capsize=2, zorder=5)
        
        # Scatter points
        ax.scatter(cv1_vals, cv2_vals, c=s['color'], marker=s['marker'],
                   s=s['s'], alpha=0.9, label=s['label'], edgecolors='white',
                   linewidths=0.5, zorder=15)
        
        # Collect label data
        for d in data:
            all_points.append((d['cv1'], d['cv2'], d['display'], s['color']))
    
    # Add diagonal line (Cv1 = Cv2)
    lims = [
        min(ax.get_xlim()[0], ax.get_ylim()[0]),
        max(ax.get_xlim()[1], ax.get_ylim()[1])
    ]
    ax.plot(lims, lims, 'k--', alpha=0.4, linewidth=1, label=r'$C_{v,1} = C_{v,2}$', zorder=1)
    
    # Set axis limits with padding
    padding = (lims[1] - lims[0]) * 0.05
    ax.set_xlim(lims[0] - padding, lims[1] + padding)
    ax.set_ylim(lims[0] - padding, lims[1] + padding)
    
    # Add labels with adjustText (optional)
    if show_labels:
        texts = []
        for x, y, label, color in all_points:
            # Initial offset away from point
            txt = ax.text(x, y, label, fontsize=FONT_LABEL, color='black', ha='center', va='center',
                          path_effects=[path_effects.withStroke(linewidth=STROKE_WIDTH, foreground='white')],
                          zorder=20)
            texts.append(txt)
        
        # Get scatter point positions for adjustText to avoid
        x_points = [p[0] for p in all_points]
        y_points = [p[1] for p in all_points]
        
        # Adjust text positions with strong parameters
        adjust_text(texts, 
                    x=x_points, y=y_points,
                    ax=ax,
                    arrowprops=dict(arrowstyle='-', color='gray', alpha=0.4, lw=0.3),
                    expand_points=(2.5, 2.5),
                    expand_text=(1.8, 1.8),
                    force_text=(1.0, 1.0),
                    force_points=(0.8, 0.8),
                    lim=500,
                    only_move={'points': 'y', 'texts': 'xy'})
    
    # Labels and title (using scaled font sizes)
    ax.set_xlabel(r'$C_{v,1}$ (Low-T region) [meV/K]', fontsize=FONT_AXIS, fontweight='bold')
    ax.set_ylabel(r'$C_{v,2}$ (High-T region) [meV/K]', fontsize=FONT_AXIS, fontweight='bold')
    ax.set_title('Heat Capacity: Low-T vs High-T Region', fontsize=FONT_TITLE, fontweight='bold')
    
    # Legend - compact (only show types that are plotted)
    ax.legend(loc='lower right', fontsize=FONT_LEGEND, framealpha=0.9, 
              handlelength=1, handletextpad=0.3, borderpad=0.4)
    
    ax.set_aspect('equal', adjustable='box')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"\n✅ Figure saved: {output_path}")
    
    return output_path


def parse_exclude_arg(exclude_str):
    """
    Parse exclude argument string to list of tuples
    
    Examples:
        "3,5,3" -> [(3,5,3)]
        "3,5,3;3,4,1" -> [(3,5,3), (3,4,1)]
        "6,8" -> [(6,8)]
    """
    if not exclude_str:
        return None
    
    exclude_list = []
    for item in exclude_str.split(';'):
        parts = item.strip().split(',')
        try:
            nums = tuple(int(p.strip()) for p in parts)
            if len(nums) in [2, 3]:
                exclude_list.append(nums)
        except ValueError:
            print(f"  Warning: Could not parse '{item}'")
    
    return exclude_list if exclude_list else None


def main():
    """Main function with command-line arguments"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Step 6.3.1: Heat Capacity Scatter Plot - Cv1 vs Cv2',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Plot all data
  python step6_3_1_cv_scatter_plot.py
  
  # Plot without error bars (cleaner)
  python step6_3_1_cv_scatter_plot.py --no-errorbars
  
  # Only plot gas-phase and supported Pt-Sn (no oxide)
  python step6_3_1_cv_scatter_plot.py --no-oxide
  
  # Only plot supported Pt-Sn-O
  python step6_3_1_cv_scatter_plot.py --no-air --no-supported
  
  # Exclude specific compositions: (3,5,3) and (3,4,1)
  python step6_3_1_cv_scatter_plot.py --exclude "3,5,3;3,4,1"
  
  # Exclude (6,8) - affects both air(6,8) and supported (6,8)
  python step6_3_1_cv_scatter_plot.py --exclude "6,8"
  
  # Combine filters
  python step6_3_1_cv_scatter_plot.py --no-air --no-errorbars --exclude "3,5,3;3,4,1"
        """)
    
    parser.add_argument('--no-air', action='store_true',
                        help='Do not plot gas-phase data')
    parser.add_argument('--no-supported', action='store_true',
                        help='Do not plot supported Pt-Sn (no oxygen) data')
    parser.add_argument('--no-oxide', action='store_true',
                        help='Do not plot supported Pt-Sn-O data')
    parser.add_argument('--no-errorbars', action='store_true',
                        help='Do not show error bars (cleaner plot)')
    parser.add_argument('--no-labels', action='store_true',
                        help='Do not show data point labels (scatter only)')
    parser.add_argument('--exclude', type=str, default=None,
                        help='Compositions to exclude, e.g., "3,5,3;3,4,1" or "6,8"')
    parser.add_argument('--fontscale', type=float, default=1.0,
                        help='Scale factor for all fonts (default 1.0, try 1.5 for larger)')
    parser.add_argument('--markerscale', type=float, default=1.0,
                        help='Scale factor for marker sizes (default 1.0, try 1.5 for larger)')
    parser.add_argument('-o', '--output', type=str, 
                        default='results/step6_1_clustering/cv1_vs_cv2_scatter.png',
                        help='Output file path')
    
    args = parser.parse_args()
    
    print("=" * 60)
    print("Step 6.3.1: Heat Capacity Scatter Plot - Cv1 vs Cv2")
    print("=" * 60)
    
    # Parse exclude list
    exclude_list = parse_exclude_arg(args.exclude)
    
    output_path = create_cv_scatter_plot(
        output_path=args.output,
        show_air=not args.no_air,
        show_supported=not args.no_supported,
        show_oxide=not args.no_oxide,
        exclude_list=exclude_list,
        show_errorbars=not args.no_errorbars,
        show_labels=not args.no_labels,
        fontscale=args.fontscale,
        markerscale=args.markerscale
    )
    
    print("\n" + "=" * 60)
    print("Done!")
    print("=" * 60)


if __name__ == '__main__':
    main()
