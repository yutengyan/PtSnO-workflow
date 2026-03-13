#!/usr/bin/env python3
"""
分析粘附能与分区边界温度(T1/T2)及氧迁移起始温度(T_onset_O)的相关性

三种粘附能类型:
1. Pt-Sn 和 O-Al2O3 的粘附能 (金属团簇与氧化铝载体)
2. Pt-Sn 和 SnO+Al2O3 的粘附能 (金属团簇与锡氧化物修饰的载体)  
3. SnO 和 Al2O3 的粘附能 (锡氧化物与氧化铝载体)

分析目标:
- 粘附能强度是否影响相变温度(T1, T2)?
- 粘附能越强，是否熔点(T1)越高?
- 粘附能与氧迁移起始温度(T_onset_O)是否相关?
  → T_onset_O 定义: AIMD轨迹中首次出现O迁移事件的温度
  → 物理意义: 粘附能越强(界面结合越紧密) → O越难迁移 → T_onset_O越高?
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib
import argparse

matplotlib.rcParams['font.family'] = ['Arial', 'SimHei']
matplotlib.rcParams['axes.unicode_minus'] = False

# ============================================================================
# 命令行参数 (模仿 plot_eadh_tm_correlation.py)
# ============================================================================
parser = argparse.ArgumentParser(
    description='分析粘附能与分区边界温度的相关性 (含偏相关/离群点检测)',
    formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('--threshold', type=float, default=2.0,
                    help='离群点检测阈值 (标准化残差距离σ, 默认: 2.0)\n'
                         '  1.5 → 更严格, 排除更多点\n'
                         '  2.0 → 默认\n'
                         '  2.5 → 更宽松, 排除更少点\n'
                         '  999 → 实际不排除任何点')
parser.add_argument('--no-labels', action='store_true',
                    help='不显示数据点标签')
parser.add_argument('--no-r', action='store_true',
                    help='图例中不显示 Pearson r 值')
parser.add_argument('--no-r2', action='store_true',
                    help='图例中不显示 R² 值')
parser.add_argument('--exclude', nargs='*', default=None,
                    help='手动指定额外排除的体系 (空格分隔)\n'
                         '  例: --exclude Sn1Pt2O1 O3Sn4Pt2 Pt3Sn3O2 Sn7Pt6O4')
parser.add_argument('--consistent-outliers', action='store_true',
                    help='A/B使用一致的离群点集合 (取并集)\n'
                         '  默认: A/B各自独立检测离群点\n'
                         '  开启后: 先分别检测, 再取并集统一排除')
parser.add_argument('--plot-clean', action='store_true',
                    help='绘制去离群点后的出版级独立图 (每组3张: a/b/c)\n'
                         '  a = 简单相关, b = 尺寸着色, c = 偏相关残差\n'
                         '  格式: 10x8in, Arial 28/34, 透明, 无标题/中文/图注框')
parser.add_argument('--clean-xticks', nargs='*', default=None,
                    help='手动指定 --plot-clean 各面板的X轴刻度值\n'
                         '  格式: PANEL:v1,v2,v3,...  (PANEL = Aa/Ab/Ac/Ba/Bb/Bc)\n'
                         '  例: --clean-xticks Aa:-0.6,-0.4,-0.2,0 Ac:-0.1,0,0.1')
parser.add_argument('--clean-yticks', nargs='*', default=None,
                    help='手动指定 --plot-clean 各面板的Y轴刻度值\n'
                         '  格式: PANEL:v1,v2,v3,...  (PANEL = Aa/Ab/Ac/Ba/Bb/Bc)\n'
                         '  例: --clean-yticks Aa:500,600,700,800 Ba:300,400,500')
parser.add_argument('--clean-panel', nargs='*', default=None,
                    help='只输出指定面板 (默认全部)\n'
                         '  例: --clean-panel Bc     → 只输出 hypothesisB (c)\n'
                         '       --clean-panel Aa Bf  → 只输出 A(a) 和 B(f) Bc2 size-effect')
parser.add_argument('--interactive', nargs='*', default=None,
                    help='交互模式: 可拖动标签调整位置\n'
                         '  不带参数 → 所有面板全部交互\n'
                         '  带参数   → 只对指定面板交互, 例:\n'
                         '    --interactive Ag-adh          → 只交互 hypothesisA (g_adh)\n'
                         '    --interactive Aa Bb Bg-size   → 多面板\n'
                         '  面板名规则: 假设字母(A/B/C) + 面板字母(a-g) + 可选后缀(-adh/-size/-f)\n'
                         '  拖动后控制台打印偏移量, 复制到 --clean-offsets 参数\n'
                         '  建议配合 --clean-panel 只打开需要的面板')
parser.add_argument('--clean-offsets', nargs='*', default=None,
                    help='手动指定标签偏移 (来自 --interactive 输出)\n'
                         '  格式: PANEL@Structure:dx,dy\n'
                         '  例: --clean-offsets Bc@Pt8Sn3O2:0.5,-20 Bc@Pt6Sn7O1:-0.3,15')
parser.add_argument('--no-stars', action='store_true',
                    help='图例中不显示显著性星号 (***/**/*), 只显示 r = +0.xx')
parser.add_argument('--heatmap-vmax', type=float, default=1.0,
                    help='summary heatmap 色标上限 (默认: 1.0, 即截断 |β|>1)\n'
                         '  1.0 → 截断到1, 色阶均匀, 视觉干净\n'
                         '  1.5 → 容纳 |β|>1, 但整体颜色偏浅')
parser.add_argument('--heatmap-r2', action='store_true', default=False,
                    help='在 summary heatmap 最前面补充两列单变量 R²\n'
                         '  R²(adh only) 和 R²(size only) (默认: 不显示)\n'
                         '  加上后共7列, 可直观对比粘附能与尺寸的独立解释力')
parser.add_argument('--exclude-c', nargs='*', default=None,
                    help='手动排除预期C (Pt8Snx) 中的特定结构 (空格分隔)\n'
                         '  例: --exclude-c Pt8Sn0\n'
                         '       --exclude-c Pt8Sn0 Pt4Sn4')
parser.add_argument('--offsets-file', default='label_offsets.json',
                    help='标签偏移量持久化文件 (默认: label_offsets.json)\n'
                         '  --interactive 模式下拖动后自动保存到该文件\n'
                         '  下次运行时自动读取, 无需再粘贴 --clean-offsets\n'
                         '  --clean-offsets 的值可覆盖文件中同名键')
parser.add_argument('--export-temperatures', action='store_true',
                    help='导出所有结构的温度汇总表到 T1_T2_summary.csv\n'
                         '  T1 = T1_lindemann (Lindemann δ=0.1 阈值温度)\n'
                         '  T2_B = T_onset_O (假说B用, 阈值 2.5/ps, 统一新阈值)\n'
                         '  T2_Bprime = T3_onset_O=T_onset_O_perO (假说B\'用, 阈值 2.5/nO /ps)')
args = parser.parse_args()

# --------------------------------------------------------------------------
# 辅助: 判断某个面板是否需要交互
#   --interactive          (args.interactive == [])   → 全部面板交互
#   --interactive Ag-adh   (args.interactive == ['Ag-adh']) → 只指定面板
#   未传 --interactive      (args.interactive is None) → 不交互
#
# panel_key 规则 (与 --interactive 参数保持一致):
#   单字母面板:  'Aa' 'Bb' 'Cc' 'Ad' 'Be' 'Bf'
#   g 面板后缀:  'Ag-adh' 'Ag-size' 'Bg-adh' 'Bg-size' 'Cg-adh' 'Cg-size'
# --------------------------------------------------------------------------
def _is_interactive(panel_key: str) -> bool:
    """返回 True 表示该面板应进入交互模式。"""
    if args.interactive is None:
        return False          # 未传 --interactive
    if len(args.interactive) == 0:
        return True           # --interactive 不带参数 → 全部
    return panel_key in args.interactive



# 分区边界温度 (来自kmeans聚类和Lindemann阈值δ=0.1) + T_onset_O (氧迁移起始温度)
PARTITION_DATA = {
    # T1_kmeans: kmeans聚类确定的固相边界
    # T1_lindemann: Lindemann指数首次超过δ=0.1的温度
    # T2_kmeans: kmeans聚类确定的液相边界
    # T_onset_O / T_onset_O_perO:
    #   M1 active_frac >= 80% (at least 4/5 independent runs show O migration, freq>0)
    #   由 m1_onset_inspector.py --af 0.80 计算，2026-03-09
    #   R²_B基线=0.248, R²_B(k=4)=0.856 (排除 O2Pt4Sn6, Pt7Sn5O1, Sn1Pt2O1, Sn3O2Pt2)
    "Sn1Pt2O1": {"T1_kmeans": 750, "T2_kmeans": 1450, "T1_lindemann": 1719.47, "T_onset_O": 1000, "T_onset_O_perO": 1000},
    "Pt2Sn2O1": {"T1_kmeans": 750, "T2_kmeans": 1350, "T1_lindemann": 1492.00, "T_onset_O": 1500, "T_onset_O_perO": 1500},
    "Pt3Sn2O1": {"T1_kmeans": 750, "T2_kmeans": 1200, "T1_lindemann": 1063.90, "T_onset_O": 1500, "T_onset_O_perO": 1500},
    "Sn3O2Pt2": {"T1_kmeans": 750, "T2_kmeans": 1400, "T1_lindemann": 1230.55, "T_onset_O": 1400, "T_onset_O_perO": 1400},
    "O3Sn4Pt2": {"T1_kmeans": 700, "T2_kmeans": 1250, "T1_lindemann":  861.84, "T_onset_O": 1100, "T_onset_O_perO": 1100},
    "Pt3Sn3O2": {"T1_kmeans": 750, "T2_kmeans": 1300, "T1_lindemann":  849.52, "T_onset_O": 1200, "T_onset_O_perO": 1200},
    "Sn3Pt4O1": {"T1_kmeans": 700, "T2_kmeans": 1250, "T1_lindemann":  617.20, "T_onset_O": 1400, "T_onset_O_perO": 1400},
    "Pt5Sn3O1": {"T1_kmeans": 700, "T2_kmeans": 1200, "T1_lindemann":  585.42, "T_onset_O": 1200, "T_onset_O_perO": 1200},
    "Pt5Sn4O1": {"T1_kmeans": 750, "T2_kmeans": 1250, "T1_lindemann":  727.74, "T_onset_O": 1400, "T_onset_O_perO": 1400},
    "O2Pt4Sn6": {"T1_kmeans": 750, "T2_kmeans": 1300, "T1_lindemann":  490.69, "T_onset_O": 1300, "T_onset_O_perO": 1300},
    "Sn6Pt5O2": {"T1_kmeans": 750, "T2_kmeans": 1350, "T1_lindemann":  537.64, "T_onset_O": 1100, "T_onset_O_perO": 1100},
    "Sn7Pt4O3": {"T1_kmeans": 800, "T2_kmeans": 1300, "T1_lindemann":  623.71, "T_onset_O": 1000, "T_onset_O_perO": 1000},
    "O3Pt5Sn7": {"T1_kmeans": 700, "T2_kmeans": 1200, "T1_lindemann":  619.82, "T_onset_O": 1200, "T_onset_O_perO": 1200},
    "Pt7Sn5O1": {"T1_kmeans": 650, "T2_kmeans": 1150, "T1_lindemann":  542.90, "T_onset_O": 1200, "T_onset_O_perO": 1200},
    "Pt7Sn6O1": {"T1_kmeans": 600, "T2_kmeans": 1150, "T1_lindemann":  558.02, "T_onset_O": 1300, "T_onset_O_perO": 1300},
    "Pt6Sn5O2": {"T1_kmeans": 550, "T2_kmeans": 750, "T1_lindemann":  739.98, "T_onset_O": 1100, "T_onset_O_perO": 1100},   # g-948-Pt6Sn5O2
    "Pt6Sn6O3": {"T1_kmeans": 450, "T2_kmeans": 650, "T1_lindemann":  649.13, "T_onset_O": 1100, "T_onset_O_perO": 1100},   # g-948-Pt6Sn6O3
    "Sn7Pt6O4": {"T1_kmeans": 400, "T2_kmeans": 700, "T1_lindemann":  710.76, "T_onset_O": 1100, "T_onset_O_perO": 1100},   # g-1051-Sn7Pt6O4
    "O2Pt7Sn7": {"T1_kmeans": 750, "T2_kmeans": 1275, "T1_lindemann":  562.89, "T_onset_O": 1250, "T_onset_O_perO": 1250},
}

# 粘附能数据 (单位: eV)
# Eadh_first: 团簇/SnO层放置到载体上但未弛豫的单点能 (初始构型, 仅作参考)
# Eadh_last:  完全弛豫后的真实吸附能 (用于 per_atom 计算)

# 类型1: Pt-Sn 和 O-Al2O3 的粘附能
ADHESION_TYPE1 = {
    "Sn1Pt2O1": {"Eadh_first": -7.626978, "Eadh_last": -6.313275, "nPt": 2, "nSn": 1, "nO": 1, "nMetal": 3},
    "Pt2Sn2O1": {"Eadh_first": -7.057111, "Eadh_last": -6.896038, "nPt": 2, "nSn": 2, "nO": 1, "nMetal": 4},
    "Pt3Sn2O1": {"Eadh_first": -8.317548, "Eadh_last": -7.698159, "nPt": 3, "nSn": 2, "nO": 1, "nMetal": 5},
    "Sn3Pt4O1": {"Eadh_first": -7.761824, "Eadh_last": -6.994111, "nPt": 4, "nSn": 3, "nO": 1, "nMetal": 7},
    "Pt5Sn3O1": {"Eadh_first": -7.441544, "Eadh_last": -6.612727, "nPt": 5, "nSn": 3, "nO": 1, "nMetal": 8},
    "Pt5Sn4O1": {"Eadh_first": -7.80499, "Eadh_last": -6.177435, "nPt": 5, "nSn": 4, "nO": 1, "nMetal": 9},
    "Pt7Sn5O1": {"Eadh_first": -7.67398, "Eadh_last": -5.716766, "nPt": 7, "nSn": 5, "nO": 1, "nMetal": 12},
    "Pt7Sn6O1": {"Eadh_first": -6.820869, "Eadh_last": -5.675522, "nPt": 7, "nSn": 6, "nO": 1, "nMetal": 13},
    "Sn3O2Pt2": {"Eadh_first": -12.622551, "Eadh_last": -11.635394, "nPt": 2, "nSn": 3, "nO": 2, "nMetal": 5},
    "Pt3Sn3O2": {"Eadh_first": -12.376116, "Eadh_last": -11.090447, "nPt": 3, "nSn": 3, "nO": 2, "nMetal": 6},
    "O2Pt4Sn6": {"Eadh_first": -11.53096, "Eadh_last": -11.018506, "nPt": 4, "nSn": 6, "nO": 2, "nMetal": 10},
    "Sn6Pt5O2": {"Eadh_first": -12.358586, "Eadh_last": -10.206829, "nPt": 5, "nSn": 6, "nO": 2, "nMetal": 11},
    "O3Sn4Pt2": {"Eadh_first": -19.457509, "Eadh_last": -14.163041, "nPt": 2, "nSn": 4, "nO": 3, "nMetal": 6},
    "Sn7Pt4O3": {"Eadh_first": -17.902354, "Eadh_last": -15.362398, "nPt": 4, "nSn": 7, "nO": 3, "nMetal": 11},
    "O3Pt5Sn7": {"Eadh_first": -18.935104, "Eadh_last": -14.9003, "nPt": 5, "nSn": 7, "nO": 3, "nMetal": 12},
    "Pt6Sn5O2": {"Eadh_first": -12.218972, "Eadh_last": -10.353433, "nPt": 6, "nSn": 5, "nO": 2, "nMetal": 11},
    "Pt6Sn6O3": {"Eadh_first": -17.810794, "Eadh_last": -13.848144, "nPt": 6, "nSn": 6, "nO": 3, "nMetal": 12},
    "Sn7Pt6O4": {"Eadh_first": -24.781255, "Eadh_last": -19.804273, "nPt": 6, "nSn": 7, "nO": 4, "nMetal": 13},
    "O2Pt7Sn7": {"Eadh_first": -12.903869, "Eadh_last": -10.041334, "nPt": 7, "nSn": 7, "nO": 2, "nMetal": 14},
}

# 类型2: Pt-Sn 和 SnO+Al2O3 的粘附能
ADHESION_TYPE2 = {
    "Sn1Pt2O1": {"Eadh_first": -6.933193, "Eadh_last": -5.4368, "nPt": 2, "nSn": 1, "nO": 1, "nMetal": 3},
    "Pt2Sn2O1": {"Eadh_first": -7.941399, "Eadh_last": -6.972672, "nPt": 2, "nSn": 2, "nO": 1, "nMetal": 4},
    "Pt3Sn2O1": {"Eadh_first": -8.329207, "Eadh_last": -7.385584, "nPt": 3, "nSn": 2, "nO": 1, "nMetal": 5},
    "Sn3O2Pt2": {"Eadh_first": -7.84662, "Eadh_last": -6.384237, "nPt": 2, "nSn": 3, "nO": 2, "nMetal": 5},
    "O3Sn4Pt2": {"Eadh_first": -7.754083, "Eadh_last": -6.982616, "nPt": 2, "nSn": 4, "nO": 3, "nMetal": 6},
    "Pt3Sn3O2": {"Eadh_first": -9.328695, "Eadh_last": -7.183718, "nPt": 3, "nSn": 3, "nO": 2, "nMetal": 6},
    "Sn3Pt4O1": {"Eadh_first": -8.046561, "Eadh_last": -7.235841, "nPt": 4, "nSn": 3, "nO": 1, "nMetal": 7},
    "Pt5Sn3O1": {"Eadh_first": -9.744485, "Eadh_last": -6.679104, "nPt": 5, "nSn": 3, "nO": 1, "nMetal": 8},
    "Pt5Sn4O1": {"Eadh_first": -7.701398, "Eadh_last": -6.30215, "nPt": 5, "nSn": 4, "nO": 1, "nMetal": 9},
    "O2Pt4Sn6": {"Eadh_first": -4.469224, "Eadh_last": -4.581913, "nPt": 4, "nSn": 6, "nO": 2, "nMetal": 10},
    "Sn6Pt5O2": {"Eadh_first": -8.228137, "Eadh_last": -6.203692, "nPt": 5, "nSn": 6, "nO": 2, "nMetal": 11},
    "Sn7Pt4O3": {"Eadh_first": -5.693145, "Eadh_last": -4.868677, "nPt": 4, "nSn": 7, "nO": 3, "nMetal": 11},
    "O3Pt5Sn7": {"Eadh_first": -7.438638, "Eadh_last": -6.085726, "nPt": 5, "nSn": 7, "nO": 3, "nMetal": 12},
    "Pt7Sn5O1": {"Eadh_first": -9.118131, "Eadh_last": -6.220163, "nPt": 7, "nSn": 5, "nO": 1, "nMetal": 12},
    "Pt7Sn6O1": {"Eadh_first": -3.684547, "Eadh_last": -1.977057, "nPt": 7, "nSn": 6, "nO": 1, "nMetal": 13},
    "Pt6Sn5O2": {"Eadh_first": -9.887237, "Eadh_last": -6.96025, "nPt": 6, "nSn": 5, "nO": 2, "nMetal": 11},
    "Pt6Sn6O3": {"Eadh_first": -8.198437, "Eadh_last": -6.793535, "nPt": 6, "nSn": 6, "nO": 3, "nMetal": 12},
    "Sn7Pt6O4": {"Eadh_first": -8.28393, "Eadh_last": -4.277785, "nPt": 6, "nSn": 7, "nO": 4, "nMetal": 13},
    "O2Pt7Sn7": {"Eadh_first": -8.043795, "Eadh_last": -6.059779, "nPt": 7, "nSn": 7, "nO": 2, "nMetal": 14},
}

# 类型3: SnO 和 Al2O3 的粘附能
ADHESION_TYPE3 = {
    "Sn1Pt2O1": {"Eadh_first": -10.968696, "Eadh_last": -11.128002, "nPt": 2, "nSn": 1, "nO": 1, "nMetal": 3},
    "Pt2Sn2O1": {"Eadh_first": -14.082388, "Eadh_last": -13.519749, "nPt": 2, "nSn": 2, "nO": 1, "nMetal": 4},
    "Pt3Sn2O1": {"Eadh_first": -18.248562, "Eadh_last": -18.075406, "nPt": 3, "nSn": 2, "nO": 1, "nMetal": 5},
    "Sn3O2Pt2": {"Eadh_first": -14.910821, "Eadh_last": -14.282319, "nPt": 2, "nSn": 3, "nO": 2, "nMetal": 5},
    "O3Sn4Pt2": {"Eadh_first": -15.837873, "Eadh_last": -15.169723, "nPt": 2, "nSn": 4, "nO": 3, "nMetal": 6},
    "Pt3Sn3O2": {"Eadh_first": -18.622631, "Eadh_last": -18.094711, "nPt": 3, "nSn": 3, "nO": 2, "nMetal": 6},
    "Sn3Pt4O1": {"Eadh_first": -27.560365, "Eadh_last": -27.486836, "nPt": 4, "nSn": 3, "nO": 1, "nMetal": 7},
    "Pt5Sn3O1": {"Eadh_first": -32.202601, "Eadh_last": -31.956066, "nPt": 5, "nSn": 3, "nO": 1, "nMetal": 8},
    "Pt5Sn4O1": {"Eadh_first": -36.731596, "Eadh_last": -36.98543, "nPt": 5, "nSn": 4, "nO": 1, "nMetal": 9},
    "O2Pt4Sn6": {"Eadh_first": -36.855768, "Eadh_last": -36.852494, "nPt": 4, "nSn": 6, "nO": 2, "nMetal": 10},
    "Sn6Pt5O2": {"Eadh_first": -42.86392, "Eadh_last": -42.552825, "nPt": 5, "nSn": 6, "nO": 2, "nMetal": 11},
    "Sn7Pt4O3": {"Eadh_first": -38.530871, "Eadh_last": -37.858729, "nPt": 4, "nSn": 7, "nO": 3, "nMetal": 11},
    "Pt6Sn5O2": {"Eadh_first": -42.197058, "Eadh_last": -42.444754, "nPt": 6, "nSn": 5, "nO": 2, "nMetal": 11},
    "O3Pt5Sn7": {"Eadh_first": -44.680551, "Eadh_last": -43.143246, "nPt": 5, "nSn": 7, "nO": 3, "nMetal": 12},
    "Pt7Sn5O1": {"Eadh_first": -51.484049, "Eadh_last": -51.692154, "nPt": 7, "nSn": 5, "nO": 1, "nMetal": 12},
    "Pt6Sn6O3": {"Eadh_first": -44.164425, "Eadh_last": -41.404068, "nPt": 6, "nSn": 6, "nO": 3, "nMetal": 12},
    "Pt7Sn6O1": {"Eadh_first": -58.765302, "Eadh_last": -59.469702, "nPt": 7, "nSn": 6, "nO": 1, "nMetal": 13},
    "Sn7Pt6O4": {"Eadh_first": -50.695858, "Eadh_last": -48.893116, "nPt": 6, "nSn": 7, "nO": 4, "nMetal": 13},
    "O2Pt7Sn7": {"Eadh_first": -57.88468, "Eadh_last": -57.527967, "nPt": 7, "nSn": 7, "nO": 2, "nMetal": 14},
}

# 类型4: 整个 PtSnO 簇 与 Al2O3 的吸附能 (E4)
# 物理含义: 金属-氧化物整体界面粘附强度 (含O效应), 归一化至 nO → 反映每个O原子对界面结合的贡献
# 数据来源: 用户 2026-03 新增 (与 Type2 同批次 AIMD, 不同分割面)
ADHESION_TYPE4 = {
    "Sn1Pt2O1": {"Eadh_last": -3.268905, "nPt": 2, "nSn": 1, "nO": 1, "nMetal": 3},
    "Pt2Sn2O1": {"Eadh_last": -2.684449, "nPt": 2, "nSn": 2, "nO": 1, "nMetal": 4},
    "Pt3Sn2O1": {"Eadh_last": -3.320033, "nPt": 3, "nSn": 2, "nO": 1, "nMetal": 5},
    "Sn3O2Pt2": {"Eadh_last": -3.174117, "nPt": 2, "nSn": 3, "nO": 2, "nMetal": 5},
    "O3Sn4Pt2": {"Eadh_last": -4.039789, "nPt": 2, "nSn": 4, "nO": 3, "nMetal": 6},
    "Pt3Sn3O2": {"Eadh_last": -3.682498, "nPt": 3, "nSn": 3, "nO": 2, "nMetal": 6},
    "Sn3Pt4O1": {"Eadh_last": -2.400166, "nPt": 4, "nSn": 3, "nO": 1, "nMetal": 7},
    "Pt5Sn3O1": {"Eadh_last": -1.968795, "nPt": 5, "nSn": 3, "nO": 1, "nMetal": 8},
    "Pt5Sn4O1": {"Eadh_last": -2.966004, "nPt": 5, "nSn": 4, "nO": 1, "nMetal": 9},
    "O2Pt4Sn6": {"Eadh_last": -2.643805, "nPt": 4, "nSn": 6, "nO": 2, "nMetal": 10},
    "Sn6Pt5O2": {"Eadh_last": -2.964427, "nPt": 5, "nSn": 6, "nO": 2, "nMetal": 11},
    "Sn7Pt4O3": {"Eadh_last": -2.233587, "nPt": 4, "nSn": 7, "nO": 3, "nMetal": 11},
    "Pt6Sn5O2": {"Eadh_last": -3.426205, "nPt": 6, "nSn": 5, "nO": 2, "nMetal": 11},
    "O3Pt5Sn7": {"Eadh_last": -3.692219, "nPt": 5, "nSn": 7, "nO": 3, "nMetal": 12},
    "Pt7Sn5O1": {"Eadh_last": -2.740339, "nPt": 7, "nSn": 5, "nO": 1, "nMetal": 12},
    "Pt6Sn6O3": {"Eadh_last": -3.214217, "nPt": 6, "nSn": 6, "nO": 3, "nMetal": 12},
    "Pt7Sn6O1": {"Eadh_last": -2.442727, "nPt": 7, "nSn": 6, "nO": 1, "nMetal": 13},
    "Sn7Pt6O4": {"Eadh_last": -6.249845, "nPt": 6, "nSn": 7, "nO": 4, "nMetal": 13},
    "O2Pt7Sn7": {"Eadh_last": -3.965053, "nPt": 7, "nSn": 7, "nO": 2, "nMetal": 14},
}

# SnO 层的原子组成 (来自 MS 结构文件中实际的 SnO 团簇)
# nSn_sno: SnO层中的Sn原子数; nO_sno: SnO层中的O原子数
# n_SnO = nSn_sno + nO_sno (Type3 每原子粘附能的分母)
# n_PtSn = nPt + nSn - nSn_sno (Type2 每原子粘附能的分母, 即 Pt-Sn 团簇原子数)
SNO_COMPOSITION = {
    "Sn1Pt2O1": {"nSn_sno": 1, "nO_sno": 1},   # n_SnO=2,  n_PtSn=2
    "Pt2Sn2O1": {"nSn_sno": 2, "nO_sno": 1},   # n_SnO=3,  n_PtSn=2
    "Pt3Sn2O1": {"nSn_sno": 2, "nO_sno": 1},   # n_SnO=3,  n_PtSn=3
    "Sn3O2Pt2": {"nSn_sno": 3, "nO_sno": 2},   # n_SnO=5,  n_PtSn=2
    "O3Sn4Pt2": {"nSn_sno": 4, "nO_sno": 3},   # n_SnO=7,  n_PtSn=2
    "Pt3Sn3O2": {"nSn_sno": 3, "nO_sno": 2},   # n_SnO=5,  n_PtSn=3
    "Sn3Pt4O1": {"nSn_sno": 2, "nO_sno": 1},   # n_SnO=3,  n_PtSn=5
    "Pt5Sn3O1": {"nSn_sno": 2, "nO_sno": 1},   # n_SnO=3,  n_PtSn=6
    "Pt5Sn4O1": {"nSn_sno": 2, "nO_sno": 1},   # n_SnO=3,  n_PtSn=7
    "O2Pt4Sn6": {"nSn_sno": 3, "nO_sno": 2},   # n_SnO=5,  n_PtSn=7
    "Sn6Pt5O2": {"nSn_sno": 3, "nO_sno": 2},   # n_SnO=5,  n_PtSn=8
    "Sn7Pt4O3": {"nSn_sno": 4, "nO_sno": 3},   # n_SnO=7,  n_PtSn=7
    "O3Pt5Sn7": {"nSn_sno": 4, "nO_sno": 3},   # n_SnO=7,  n_PtSn=8
    "Pt7Sn5O1": {"nSn_sno": 2, "nO_sno": 1},   # n_SnO=3,  n_PtSn=10
    "Pt7Sn6O1": {"nSn_sno": 1, "nO_sno": 1},   # n_SnO=2,  n_PtSn=12
    "Pt6Sn5O2": {"nSn_sno": 3, "nO_sno": 2},   # n_SnO=5,  n_PtSn=8
    "Pt6Sn6O3": {"nSn_sno": 4, "nO_sno": 3},   # n_SnO=7,  n_PtSn=8
    "Sn7Pt6O4": {"nSn_sno": 4, "nO_sno": 4},   # n_SnO=8,  n_PtSn=9
    "O2Pt7Sn7": {"nSn_sno": 3, "nO_sno": 2},   # n_SnO=5,  n_PtSn=11
}

# ============================================================================
# 预期C 数据: Pt8Snx 系列 (来自 plot_eadh_tm_correlation.py)
#   Eadh = Type2 每原子粘附能 (eV/atom), Tm = T1_lindemann (K)
#   两个子系列: Pt8Snx (固定 nPt=8, nSn=0~10) 和 PtxSn8-x (nPt+nSn=8)
# ============================================================================
PT8SNX_DATA = {
    # series, nPt, nSn, Eadh/atom, Tm
    'Pt3Sn5':    {'series': 'PtxSn8-x', 'nPt': 3, 'nSn': 5, 'Eadh': -0.1161645,    'Tm': 535.8},
    'Pt4Sn4':    {'series': 'PtxSn8-x', 'nPt': 4, 'nSn': 4, 'Eadh':  0.000542875,  'Tm': 505.3},
    'Pt5Sn3':    {'series': 'PtxSn8-x', 'nPt': 5, 'nSn': 3, 'Eadh': -0.319300375,  'Tm': 665.1},
    'Pt6Sn2':    {'series': 'PtxSn8-x', 'nPt': 6, 'nSn': 2, 'Eadh': -0.2375385,    'Tm': 602.6},
    'Pt7Sn1':    {'series': 'PtxSn8-x', 'nPt': 7, 'nSn': 1, 'Eadh': -0.293225125,  'Tm': 641.3},
    'Pt8Sn0':    {'series': 'Pt8Snx',   'nPt': 8, 'nSn': 0, 'Eadh': -0.661488,     'Tm': 800.0},
    'Pt8Sn1':    {'series': 'Pt8Snx',   'nPt': 8, 'nSn': 1, 'Eadh': -0.315793444,  'Tm': 604.0},
    'Pt8Sn2':    {'series': 'Pt8Snx',   'nPt': 8, 'nSn': 2, 'Eadh': -0.2898329,    'Tm': 705.0},
    'Pt8Sn3':    {'series': 'Pt8Snx',   'nPt': 8, 'nSn': 3, 'Eadh': -0.201466091,  'Tm': 574.0},
    'Pt8Sn4':    {'series': 'Pt8Snx',   'nPt': 8, 'nSn': 4, 'Eadh': -0.136585333,  'Tm': 566.3},
    'Pt8Sn5':    {'series': 'Pt8Snx',   'nPt': 8, 'nSn': 5, 'Eadh': -0.114243077,  'Tm': 575.5},
    'Pt8Sn6':    {'series': 'Pt8Snx',   'nPt': 8, 'nSn': 6, 'Eadh': -0.135212286,  'Tm': 560.2367},
    'Pt8Sn7':    {'series': 'Pt8Snx',   'nPt': 8, 'nSn': 7, 'Eadh': -0.152302333,  'Tm': 577.1},
    'Pt8Sn8':    {'series': 'Pt8Snx',   'nPt': 8, 'nSn': 8, 'Eadh': -0.08049725,   'Tm': 477.0},
    'Pt8Sn9':    {'series': 'Pt8Snx',   'nPt': 8, 'nSn': 9, 'Eadh': -0.040358529,  'Tm': 504.5},
    'Pt8Sn10':   {'series': 'Pt8Snx',   'nPt': 8, 'nSn': 10,'Eadh':  0.000571111,  'Tm': 472.2},
}


def build_pt8snx_dataframe():
    """构建 Pt8Snx 系列的 DataFrame (用于预期C偏相关分析)"""
    rows = []
    for name, d in PT8SNX_DATA.items():
        rows.append({
            'Structure': name,
            'series': d['series'],
            'nPt': d['nPt'],
            'nSn': d['nSn'],
            'nMetal': d['nPt'] + d['nSn'],
            'Eadh_per_atom': d['Eadh'],   # = Type2/atom
            'T1_lindemann': d['Tm'],
        })
    return pd.DataFrame(rows)


def build_dataframe():
    """构建包含所有数据的DataFrame"""
    data = []
    
    for structure in PARTITION_DATA.keys():
        row = {
            "Structure": structure,
            "T1_kmeans": PARTITION_DATA[structure]["T1_kmeans"],
            "T2_kmeans": PARTITION_DATA[structure]["T2_kmeans"],
            "T1_lindemann": PARTITION_DATA[structure].get("T1_lindemann", 0),
            "T_onset_O": PARTITION_DATA[structure].get("T_onset_O", np.nan),
            "T_range": PARTITION_DATA[structure]["T2_kmeans"] - PARTITION_DATA[structure]["T1_kmeans"],
        }
        
        # 添加粘附能数据
        if structure in ADHESION_TYPE1:
            row["Type1_Eadh_first"] = ADHESION_TYPE1[structure]["Eadh_first"]
            row["Type1_Eadh_last"] = ADHESION_TYPE1[structure]["Eadh_last"]
            row["Type1_Eadh_avg"] = (row["Type1_Eadh_first"] + row["Type1_Eadh_last"]) / 2
            row["nPt"] = ADHESION_TYPE1[structure]["nPt"]
            row["nSn"] = ADHESION_TYPE1[structure]["nSn"]
            row["nO"] = ADHESION_TYPE1[structure]["nO"]
            row["nO_total"] = ADHESION_TYPE1[structure]["nO"]
            row["nMetal"] = ADHESION_TYPE1[structure]["nMetal"]
        
        if structure in ADHESION_TYPE2:
            row["Type2_Eadh_first"] = ADHESION_TYPE2[structure]["Eadh_first"]
            row["Type2_Eadh_last"] = ADHESION_TYPE2[structure]["Eadh_last"]
            row["Type2_Eadh_avg"] = (row["Type2_Eadh_first"] + row["Type2_Eadh_last"]) / 2
        
        if structure in ADHESION_TYPE3:
            row["Type3_Eadh_first"] = ADHESION_TYPE3[structure]["Eadh_first"]
            row["Type3_Eadh_last"] = ADHESION_TYPE3[structure]["Eadh_last"]
            row["Type3_Eadh_avg"] = (row["Type3_Eadh_first"] + row["Type3_Eadh_last"]) / 2
        
        # 计算每原子粘附能 (均使用 Eadh_last — 弛豫后真实吸附能)
        # Type1: Pt-Sn / O-Al2O3   → Eadh_last / nMetal       (nMetal = nPt + nSn)
        # Type2: Pt-Sn / SnO+Al2O3 → Eadh_last / n_PtSn       (n_PtSn = nPt + nSn - nSn_sno, 即纯Pt-Sn团簇原子数)
        # Type3: SnO / Al2O3        → Eadh_last / n_SnO        (n_SnO  = nSn_sno + nO_sno, 即SnO层总原子数)
        if "nMetal" in row and row["nMetal"] > 0:
            if "Type1_Eadh_last" in row:
                row["Type1_per_atom"] = row["Type1_Eadh_last"] / row["nMetal"]
            
            if structure in SNO_COMPOSITION:
                sno = SNO_COMPOSITION[structure]
                n_SnO = sno["nSn_sno"] + sno["nO_sno"]
                n_PtSn = row["nPt"] + row["nSn"] - sno["nSn_sno"]  # 不在SnO层的Sn + 所有Pt
                row["nSn_sno"] = sno["nSn_sno"]
                row["nO_sno"] = sno["nO_sno"]
                row["nSn_SnO"] = sno["nSn_sno"]
                row["nO_SnO"] = sno["nO_sno"]
                row["n_SnO"] = n_SnO
                row["n_PtSn"] = n_PtSn
                row["nAtoms_SnO"] = n_SnO
                row["nAtoms_PtSn"] = n_PtSn
                
                if "Type2_Eadh_last" in row and n_PtSn > 0:
                    row["Type2_per_atom"] = row["Type2_Eadh_last"] / n_PtSn
                if "Type3_Eadh_last" in row and n_SnO > 0:
                    row["Type3_per_atom"] = row["Type3_Eadh_last"] / n_SnO

        # Type4: 整个 PtSnO 簇 / Al2O3 吸附能, 归一化至 nO
        # Type4_per_nO = Eadh_last / nO  (每O原子的整体界面粘附能)
        if structure in ADHESION_TYPE4:
            t4 = ADHESION_TYPE4[structure]
            row["Type4_Eadh_last"] = t4["Eadh_last"]
            nO_t4 = t4["nO"]
            if nO_t4 > 0:
                row["Type4_per_nO"] = t4["Eadh_last"] / nO_t4

        # T3_onset_O: 每个 O 原子归一化的迁移起始温度
        # 来源: PARTITION_DATA["T_onset_O_perO"] — 由 process_melting_summary.py 以阈值 2.5/nO /ps 重新计算
        # 物理含义: 体系总迁移频率等效于每个O原子独立达到2.5/ps时的最低温度
        # 与 T_onset_O/nO 的区别: 本字段是对原始AIMD数据重新判定，而非简单数学除法
        t_perO = PARTITION_DATA.get(structure, {}).get("T_onset_O_perO", np.nan)
        if not np.isnan(float(t_perO)) if t_perO is not None else True:
            row["T3_onset_O"] = float(t_perO) if t_perO is not None else np.nan
        else:
            row["T3_onset_O"] = np.nan

        data.append(row)
    
    return pd.DataFrame(data)


def analyze_correlation(df):
    """分析粘附能与T1/T2的相关性"""
    
    print("\n")
    print("=" * 80)
    print("   粘附能与分区边界温度(T1/T2)的相关性分析")
    print("=" * 80)
    
    # 定义要分析的列组合
    adhesion_types = [
        ("Type1_Eadh_avg", "Type1_per_atom", "Pt-Sn / O-Al2O3", "金属团簇与氧化铝载体"),
        ("Type2_Eadh_avg", "Type2_per_atom", "Pt-Sn / SnO+Al2O3", "金属团簇与SnO修饰载体"),
        ("Type3_Eadh_avg", "Type3_per_atom", "SnO / Al2O3", "SnO层与氧化铝载体"),
    ]
    
    temp_cols = [
        ("T1_kmeans", "T1_kmeans (聚类法)"),
        ("T1_lindemann", "T1_lindemann (δ=0.1)"),
        ("T2_kmeans", "T2_kmeans (液相边界)"),
        ("T_onset_O", "T_onset_O (O迁移起始)"),
        ("T3_onset_O", "T3_onset_O (O迁移/nO归一化)"),
    ]
    
    results = []
    
    # 分块输出，每种粘附能类型单独一个表格
    for total_col, per_atom_col, type_name, description in adhesion_types:
        print(f"\n{'─' * 80}")
        print(f"  {type_name}")
        print(f"  {description}")
        print(f"{'─' * 80}")
        print(f"  {'温度参数':<25} │ {'类型':<8} │ {'r':>7} │ {'R²':>7} │ {'p值':>12} │ {'显著性':>6}")
        print(f"  {'─' * 70}")
        
        for temp_col, temp_name in temp_cols:
            for col, col_type in [(total_col, "总能"), (per_atom_col, "每原子")]:
                valid_df = df[[col, temp_col]].dropna()
                
                if len(valid_df) >= 3:
                    r, p = stats.pearsonr(valid_df[col], valid_df[temp_col])
                    r2 = r ** 2
                    
                    # 显著性标记
                    if p < 0.001:
                        sig = "***"
                        highlight = " <=="
                    elif p < 0.01:
                        sig = "**"
                        highlight = " <="
                    elif p < 0.05:
                        sig = "*"
                        highlight = ""
                    else:
                        sig = "ns"
                        highlight = ""
                    
                    print(f"  {temp_name:<25} │ {col_type:<8} │ {r:>+7.3f} │ {r2:>7.3f} │ {p:>12.2e} │ {sig:>4}{highlight}")
                    
                    results.append({
                        "粘附能类型": type_name,
                        "数据类型": col_type,
                        "温度参数": temp_name,
                        "r": r,
                        "R2": r2,
                        "p": p,
                        "显著性": sig
                    })
    
    # 打印关键发现摘要
    print("\n")
    print("=" * 80)
    print("   关键发现摘要")
    print("=" * 80)
    
    # 找出最显著的相关性
    results_df = pd.DataFrame(results)
    significant = results_df[results_df['p'] < 0.05].sort_values('p')
    
    if len(significant) > 0:
        print("\n  显著相关性 (p < 0.05):")
        print(f"  {'─' * 75}")
        for _, row in significant.head(6).iterrows():
            direction = "↑" if row['r'] > 0 else "↓"
            meaning = "粘附能越强→温度越高" if row['r'] < 0 else "粘附能越强→温度越低"
            符合预期 = "✓ 符合预期" if row['r'] < 0 else "✗ 异常"
            print(f"  {direction} {row['粘附能类型']:<20} ({row['数据类型']}) vs {row['温度参数']}")
            print(f"      r = {row['r']:+.3f}, R² = {row['R2']:.3f}, p = {row['p']:.2e} {row['显著性']}")
            print(f"      解释: {meaning} ({符合预期})")
            print()
    
    print(f"  {'─' * 75}")
    print("  注: * p<0.05, ** p<0.01, *** p<0.001, ns=不显著")
    print("  注: 负相关(r<0)表示粘附能越强(越负)，熔点越高 → 符合物理预期")
    print()
    
    return results_df


def analyze_by_oxygen_content(df):
    """按氧含量分组分析"""
    
    print("\n" + "="*100)
    print("按氧含量(nO)分组的相关性分析")
    print("="*100)
    
    for nO in sorted(df['nO'].unique()):
        sub_df = df[df['nO'] == nO]
        if len(sub_df) >= 3:
            print(f"\n--- nO = {nO} (n={len(sub_df)}) ---")
            
            # Type1 vs T1
            valid = sub_df[['Type1_Eadh_avg', 'T1_kmeans']].dropna()
            if len(valid) >= 3:
                r, p = stats.pearsonr(valid['Type1_Eadh_avg'], valid['T1_kmeans'])
                sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
                print(f"  Type1(Pt-Sn/O-Al2O3) vs T1: r={r:+.4f}, p={p:.4e} {sig}")
            
            # Type3 vs T1
            valid = sub_df[['Type3_Eadh_avg', 'T1_kmeans']].dropna()
            if len(valid) >= 3:
                r, p = stats.pearsonr(valid['Type3_Eadh_avg'], valid['T1_kmeans'])
                sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
                print(f"  Type3(SnO/Al2O3) vs T1: r={r:+.4f}, p={p:.4e} {sig}")


def plot_correlations(df, output_dir="results/adhesion_analysis"):
    """绘制相关性图"""
    import os
    os.makedirs(output_dir, exist_ok=True)
    
    fig, axes = plt.subplots(3, 3, figsize=(15, 15))
    
    adhesion_cols = [
        ("Type1_Eadh_avg", "Pt-Sn / O-Al2O3", "tab:blue"),
        ("Type2_Eadh_avg", "Pt-Sn / SnO+Al2O3", "tab:green"),
        ("Type3_Eadh_avg", "SnO / Al2O3", "tab:red"),
    ]
    
    # 第一行: vs T1
    for i, (col, name, color) in enumerate(adhesion_cols):
        ax = axes[0, i]
        valid = df[[col, 'T1_kmeans', 'nO']].dropna()
        
        # 按nO着色
        scatter = ax.scatter(valid[col], valid['T1_kmeans'], 
                            c=valid['nO'], cmap='viridis', 
                            s=80, alpha=0.7, edgecolors='black')
        
        # 添加拟合线
        if len(valid) >= 3:
            z = np.polyfit(valid[col], valid['T1_kmeans'], 1)
            p = np.poly1d(z)
            x_line = np.linspace(valid[col].min(), valid[col].max(), 100)
            ax.plot(x_line, p(x_line), '--', color=color, alpha=0.8, lw=2)
            
            r, pval = stats.pearsonr(valid[col], valid['T1_kmeans'])
            ax.set_title(f'{name}\nr={r:.3f}, p={pval:.3e}')
        
        ax.set_xlabel(f'Eadh (eV)')
        ax.set_ylabel('T1 (K)')
        ax.grid(True, alpha=0.3)
    
    # 第二行: vs T2
    for i, (col, name, color) in enumerate(adhesion_cols):
        ax = axes[1, i]
        valid = df[[col, 'T2_kmeans', 'nO']].dropna()
        
        scatter = ax.scatter(valid[col], valid['T2_kmeans'], 
                            c=valid['nO'], cmap='viridis', 
                            s=80, alpha=0.7, edgecolors='black')
        
        if len(valid) >= 3:
            z = np.polyfit(valid[col], valid['T2_kmeans'], 1)
            p = np.poly1d(z)
            x_line = np.linspace(valid[col].min(), valid[col].max(), 100)
            ax.plot(x_line, p(x_line), '--', color=color, alpha=0.8, lw=2)
            
            r, pval = stats.pearsonr(valid[col], valid['T2_kmeans'])
            ax.set_title(f'{name}\nr={r:.3f}, p={pval:.3e}')
        
        ax.set_xlabel(f'Eadh (eV)')
        ax.set_ylabel('T2 (K)')
        ax.grid(True, alpha=0.3)
    
    # 第三行: vs T_onset_O
    for i, (col, name, color) in enumerate(adhesion_cols):
        ax = axes[2, i]
        valid = df[[col, 'T_onset_O', 'nO']].dropna()
        
        scatter = ax.scatter(valid[col], valid['T_onset_O'], 
                            c=valid['nO'], cmap='viridis', 
                            s=80, alpha=0.7, edgecolors='black')
        
        if len(valid) >= 3:
            z = np.polyfit(valid[col], valid['T_onset_O'], 1)
            p = np.poly1d(z)
            x_line = np.linspace(valid[col].min(), valid[col].max(), 100)
            ax.plot(x_line, p(x_line), '--', color=color, alpha=0.8, lw=2)
            
            r, pval = stats.pearsonr(valid[col], valid['T_onset_O'])
            ax.set_title(f'{name}\nr={r:.3f}, p={pval:.3e}')
        
        ax.set_xlabel(f'Eadh (eV)')
        ax.set_ylabel('T_onset_O (K)')
        ax.grid(True, alpha=0.3)
    
    # 添加colorbar
    cbar = fig.colorbar(scatter, ax=axes, label='nO (氧原子数)', shrink=0.6)
    
    plt.suptitle('粘附能 vs 分区边界温度 (T1/T2/T_onset_O)\n颜色表示氧原子数nO', fontsize=14)
    plt.tight_layout(rect=[0, 0, 0.92, 0.95])
    
    output_path = f"{output_dir}/adhesion_vs_partition_correlation.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"\n[OK] 相关性图已保存: {output_path}")
    
    plt.close()


def plot_per_atom_adhesion(df, output_dir="results/adhesion_analysis"):
    """绘制每原子粘附能的相关性"""
    import os
    os.makedirs(output_dir, exist_ok=True)
    
    # 使用 build_dataframe 中已正确计算的每原子粘附能
    # Type1_per_atom = Eadh_avg / nMetal (全部金属原子)
    # Type2_per_atom = Eadh_avg / n_PtSn (Pt-Sn团簇原子数)
    # Type3_per_atom = Eadh_avg / n_SnO  (SnO层原子数)
    df['Type1_Eadh_per_atom'] = df['Type1_per_atom']
    df['Type2_Eadh_per_atom'] = df['Type2_per_atom']
    df['Type3_Eadh_per_atom'] = df['Type3_per_atom']
    # Type1_total: 总吸附能（非每原子），即 Eadh_last 原始值
    if 'Type1_Eadh_last' in df.columns:
        df['Type1_total'] = df['Type1_Eadh_last']
    # Type1_per_nO: Type1总吸附能 / nO（每O原子的Type1粘附能贡献）
    if 'Type1_Eadh_last' in df.columns and 'nO' in df.columns:
        df['Type1_per_nO'] = df['Type1_Eadh_last'] / df['nO']
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    
    per_atom_cols = [
        ("Type1_Eadh_per_atom", "Type1: Pt-Sn/O-Al2O3", "tab:blue"),
        ("Type2_Eadh_per_atom", "Type2: Pt-Sn/SnO+Al2O3", "tab:green"),
        ("Type3_Eadh_per_atom", "Type3: SnO/Al2O3", "tab:red"),
    ]
    
    # 第一行: vs T1_lindemann
    for i, (col, name, color) in enumerate(per_atom_cols):
        ax = axes[0, i]
        valid = df[[col, 'T1_lindemann', 'nO']].dropna()
        scatter = ax.scatter(valid[col], valid['T1_lindemann'],
                            c=valid['nO'], cmap='viridis', s=80, alpha=0.7, edgecolors='black')
        if len(valid) >= 3:
            z = np.polyfit(valid[col], valid['T1_lindemann'], 1)
            p_fit = np.poly1d(z)
            x_line = np.linspace(valid[col].min(), valid[col].max(), 100)
            ax.plot(x_line, p_fit(x_line), '--', color=color, alpha=0.8, lw=2)
            r, pval = stats.pearsonr(valid[col], valid['T1_lindemann'])
            ax.set_title(f'{name}\nr={r:.3f}, p={pval:.3e}')
        ax.set_xlabel('Eadh per metal atom (eV/atom)')
        ax.set_ylabel('T1_lindemann (K)')
        ax.grid(True, alpha=0.3)
    
    # 第二行: vs T_onset_O
    for i, (col, name, color) in enumerate(per_atom_cols):
        ax = axes[1, i]
        valid = df[[col, 'T_onset_O', 'nO']].dropna()
        scatter = ax.scatter(valid[col], valid['T_onset_O'],
                            c=valid['nO'], cmap='viridis', s=80, alpha=0.7, edgecolors='black')
        if len(valid) >= 3:
            z = np.polyfit(valid[col], valid['T_onset_O'], 1)
            p_fit = np.poly1d(z)
            x_line = np.linspace(valid[col].min(), valid[col].max(), 100)
            ax.plot(x_line, p_fit(x_line), '--', color=color, alpha=0.8, lw=2)
            r, pval = stats.pearsonr(valid[col], valid['T_onset_O'])
            ax.set_title(f'{name}\nr={r:.3f}, p={pval:.3e}')
        ax.set_xlabel('Eadh per metal atom (eV/atom)')
        ax.set_ylabel('T_onset_O (K)')
        ax.grid(True, alpha=0.3)
    
    plt.colorbar(scatter, ax=axes, label='nO', shrink=0.8)
    plt.suptitle('每金属原子粘附能 vs T1_lindemann / T_onset_O', fontsize=12)
    plt.tight_layout(rect=[0, 0, 0.92, 0.95])
    
    output_path = f"{output_dir}/adhesion_per_atom_vs_T1.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"[OK] 每原子粘附能图已保存: {output_path}")
    
    plt.close()


def print_summary_table(df):
    """打印数据汇总表"""
    
    print("\n")
    print("=" * 100)
    print("   数据汇总表")
    print("=" * 100)
    
    # 表格1: 基本信息和温度边界
    print("\n  [1] 结构信息与分区边界温度")
    print(f"  {'─' * 105}")
    print(f"  {'结构':<12} │ {'nO':>3} │ {'nMetal':>6} │ {'T1_kmeans':>9} │ {'T1_lind':>8} │ {'T2_kmeans':>9} │ {'T_onset_O':>9} │ {'T3_onset_O':>11} │ {'ΔT1':>6}")
    print(f"  {'─' * 105}")
    
    for _, row in df.sort_values('T1_lindemann', ascending=False).iterrows():
        t1_km = row['T1_kmeans']
        t1_lind = row.get('T1_lindemann', 0)
        t2_km = row['T2_kmeans']
        t_onset = row.get('T_onset_O', np.nan)
        t3_onset = row.get('T3_onset_O', np.nan)
        n_O = row.get('nO', np.nan)
        delta_t1 = t1_lind - t1_km
        t_onset_str = f"{t_onset:>8.0f}K" if pd.notna(t_onset) else f"{'---':>9}"
        t3_onset_str = f"{t3_onset:>9.1f}K" if pd.notna(t3_onset) else f"{'---':>11}"
        n_O_str = f"{n_O:>3.0f}" if pd.notna(n_O) else f"{'?':>3}"
        print(f"  {row['Structure']:<12} │ {n_O_str} │ {row['nMetal']:>6.0f} │ {t1_km:>8.0f}K │ {t1_lind:>7.0f}K │ {t2_km:>8.0f}K │ {t_onset_str} │ {t3_onset_str} │ {delta_t1:>+5.0f}K")
    
    # 表格2: 每原子粘附能 (使用正确的分母)
    print(f"\n  [2] 每原子粘附能 (eV/atom, 均使用 Eadh_last — 弛豫后真实吸附能)")
    print(f"       Type1/at = Eadh_last / nMetal  (nMetal = nPt+nSn, 总金属原子数)")
    print(f"       Type2/at = Eadh_last / n_PtSn  (n_PtSn = nPt+nSn-nSn_sno, 纯Pt-Sn团簇原子数)")
    print(f"       Type3/at = Eadh_last / n_SnO   (n_SnO  = nSn_sno+nO_sno, SnO层总原子数)")
    print(f"  {'─' * 130}")
    print(f"  {'结构':<12} │ {'nO':>3} │ {'nMetal':>6} │ {'n_PtSn':>6} │ {'n_SnO':>5} │ {'T1_lind':>7} │ {'T_onset_O':>9} │ {'T3_onset_O':>11} │ {'Type1/at':>9} │ {'Type2/at':>9} │ {'Type3/at':>9} │ 备注")
    print(f"  {'─' * 130}")
    
    for _, row in df.sort_values('T1_lindemann', ascending=False).iterrows():
        t1_lind = row.get('T1_lindemann', 0)
        t_onset = row.get('T_onset_O', np.nan)
        t3_onset = row.get('T3_onset_O', np.nan)
        t1_per = row.get('Type1_per_atom', 0)
        t2_per = row.get('Type2_per_atom', 0)
        t3_per = row.get('Type3_per_atom', 0)
        n_PtSn = row.get('n_PtSn', row.get('nMetal', 0))
        n_SnO = row.get('n_SnO', 0)
        n_O = row.get('nO', np.nan)
        
        # 根据T1_lindemann添加备注
        if t1_lind >= 1500:
            note = "高熔点"
        elif t1_lind >= 1000:
            note = "中高熔点"
        elif t1_lind >= 700:
            note = "中等熔点"
        else:
            note = "低熔点"
        
        t_onset_str = f"{t_onset:>8.0f}K" if pd.notna(t_onset) else f"{'---':>9}"
        t3_onset_str = f"{t3_onset:>9.1f}K" if pd.notna(t3_onset) else f"{'---':>11}"
        n_O_str = f"{n_O:>3.0f}" if pd.notna(n_O) else f"{'?':>3}"
        print(f"  {row['Structure']:<12} │ {n_O_str} │ {row['nMetal']:>6.0f} │ {n_PtSn:>6.0f} │ {n_SnO:>5.0f} │ {t1_lind:>6.0f}K │ {t_onset_str} │ {t3_onset_str} │ {t1_per:>9.3f} │ {t2_per:>9.3f} │ {t3_per:>9.3f} │ {note}")


# ============================================================================
# 偏相关/残差分析 + 离群点检测
# ============================================================================

def _sig_label(p):
    """p值显著性标记"""
    if p < 0.001: return '***'
    if p < 0.01:  return '**'
    if p < 0.05:  return '*'
    return 'ns'


def _calc_residuals(v, adh_col, temp_col, size_col):
    """
    计算残差: 从粘附能和温度中分别线性回归掉尺寸变量的影响。
    
    原理:
      1. 拟合 adh = a*size + b → 残差 = 粘附能中"非尺寸"的部分
      2. 拟合 temp = c*size + d → 残差 = 温度中"非尺寸"的部分
      3. 残差之间的相关性 = 控制尺寸后的净粘附能效应 (等价于偏相关)
    """
    v = v.copy()
    sl_a, int_a, _, _, _ = stats.linregress(v[size_col], v[adh_col])
    sl_t, int_t, _, _, _ = stats.linregress(v[size_col], v[temp_col])
    v['adh_pred'] = sl_a * v[size_col] + int_a
    v['T_pred'] = sl_t * v[size_col] + int_t
    v['adh_res'] = v[adh_col] - v['adh_pred']
    v['T_res'] = v[temp_col] - v['T_pred']
    return v, (sl_a, int_a, sl_t, int_t)


def _detect_outliers(v, threshold=2.0, manual_exclude=None):
    """
    基于标准化残差距离检测离群点。
    
    方法: 将粘附能残差和温度残差分别标准化为z分数,
          计算每个点到原点的欧氏距离, 距离 > threshold 视为离群。
    
    Parameters:
        v:               DataFrame (需含 adh_res, T_res, Structure)
        threshold:       标准化残差距离阈值 (σ), 越小越严格
        manual_exclude:  手动额外排除的体系名列表
    
    Returns: (outlier_names_list, v_with_dist, auto_outliers, manual_outliers)
    """
    v = v.copy()
    adh_std = v['adh_res'].std()
    T_std = v['T_res'].std()
    if adh_std == 0 or T_std == 0:
        v['dist_z'] = 0
        return [], v, [], []
    v['adh_res_z'] = (v['adh_res'] - v['adh_res'].mean()) / adh_std
    v['T_res_z'] = (v['T_res'] - v['T_res'].mean()) / T_std
    v['dist_z'] = np.sqrt(v['adh_res_z']**2 + v['T_res_z']**2)
    auto_outliers = v[v['dist_z'] > threshold]['Structure'].tolist()
    
    # 合并手动排除
    manual_added = []
    if manual_exclude:
        for s in manual_exclude:
            if s in v['Structure'].values and s not in auto_outliers:
                manual_added.append(s)
    
    all_outliers = list(dict.fromkeys(auto_outliers + manual_added))  # 保序去重
    return all_outliers, v, auto_outliers, manual_added


def analyze_partial_correlation(df, output_dir="results/adhesion_analysis"):
    """
    偏相关分析: 控制尺寸(nMetal)后, 分析粘附能与温度的净相关性。
    
    包含:
    1. 两组关键关系的偏相关 (残差分析)
    2. 离群点检测 (标准化残差距离 > 2σ)
    3. 去掉离群点后重新计算
    
    两组关键关系:
      预期A: Type2/at → T1_lindemann (金属粘附 → 金属熔化)
      预期B: Type3/at → T_onset_O   (SnO粘附 → O迁移)
    """
    import os
    os.makedirs(output_dir, exist_ok=True)
    
    print("\n")
    print("=" * 90)
    print("   偏相关分析: 控制尺寸后的净粘附能—温度关系")
    print("=" * 90)
    print("   原理: 简单÷n不能完全消除尺寸效应 (因为 Eadh ∝ n^(2/3), 不是 ∝ n)")
    print("         偏相关 = 先从粘附能和温度中分别去除尺寸的线性影响, 再看残差相关")
    print("   控制变量: A→n_PtSn, B→n_SnO, C→nMetal (各自与粘附能分母口径一致)")
    print()
    
    # 定义分析组
    # 预期C 使用独立的 Pt8Snx 数据集 (通过 df_override 传入)
    df_c = build_pt8snx_dataframe()
    # --exclude-c: 手动排除 Pt8Snx 系列中的特定结构
    if args.exclude_c:
        before = len(df_c)
        df_c = df_c[~df_c['Structure'].isin(args.exclude_c)].copy()
        excluded = [s for s in args.exclude_c if s in PT8SNX_DATA]
        print(f"  ⚙ --exclude-c: 从预期C排除 {excluded}  ({before} → {len(df_c)} 个点)")
    analyses = [
        {
            'label': '预期A',
            'adh_col': 'Type2_per_atom', 'temp_col': 'T1_lindemann',
            'size_col': 'n_PtSn',
            'adh_name': 'Type2/at (÷n_PtSn)',
            'temp_name': 'T1_lindemann',
            'physical': '金属粘附决定金属熔化温度',
            'color': 'steelblue', 'color_clean': 'darkorange',
        },
        {
            'label': '预期B',
            'adh_col': 'Type3_per_atom', 'temp_col': 'T_onset_O',
            'size_col': 'n_SnO',
            'adh_name': 'Type3/at (÷n_SnO)',
            'temp_name': 'T_onset_O',
            'physical': 'SnO界面粘附决定O迁移温度',
            'color': 'seagreen', 'color_clean': 'teal',
        },
        {
            'label': "预期B'",
            'adh_col': 'Type3_per_atom', 'temp_col': 'T3_onset_O',
            'size_col': 'n_SnO',
            'adh_name': 'Type3/at (÷n_SnO)',
            'temp_name': 'T3_onset_O (÷nO)',
            'physical': "SnO粘附→O迁移(T3=T_onset/nO归一化，消除nO尺寸效应)",
            'color': 'mediumpurple', 'color_clean': 'rebeccapurple',
        },
        {
            'label': "预期B''",
            'adh_col': 'Type1_total', 'temp_col': 'T3_onset_O',
            'size_col': 'nMetal',
            'adh_name': 'Type1_total (总吸附能, 非每原子)',
            'temp_name': 'T3_onset_O (÷nO)',
            'physical': "总粘附能(非归一化)→O迁移温度T3: 粘附能总量主导氧活化",
            'color': 'darkorange', 'color_clean': 'saddlebrown',
        },
        {
            'label': '预期B3',
            'adh_col': 'Type1_per_nO', 'temp_col': 'T3_onset_O',
            'size_col': 'nO',
            'adh_name': 'Type1/nO (÷nO)',
            'temp_name': 'T3_onset_O (÷nO)',
            'physical': "每O的Type1粘附能→T3: 界面粘附强度(归一化至O数)决定O迁移难易",
            'color': 'crimson', 'color_clean': 'darkred',
        },
        {
            'label': '预期D',
            'adh_col': 'Type4_per_nO', 'temp_col': 'T3_onset_O',
            'size_col': 'nO',
            'adh_name': 'Type4/nO (÷nO)',
            'temp_name': 'T3_onset_O (÷nO)',
            'physical': "整体PtSnO/Al2O3粘附能(每O归一化)→T3: 含O的整体界面结合强度决定O迁移",
            'color': 'chocolate', 'color_clean': 'saddlebrown',
            'axis_adh_label': r'$E_{adh}^2$ (eV/atom)',
            'axis_adh_res_label': r'$E_{adh}^2$ residual (eV/atom)',
            'axis_temp_label': r'$T_2$ (K)',
            'axis_temp_res_label': r'$T_2$ residual (K)',
        },
        {
            'label': '预期C',
            'adh_col': 'Eadh_per_atom', 'temp_col': 'T1_lindemann',
            'size_col': 'nMetal',
            'adh_name': 'Eadh/atom (Pt₈Snₓ)',
            'temp_name': 'T_m',
            'physical': 'Pt8Snx系列: 粘附能→熔点 (无O, 纯尺寸/组分效应)',
            'color': '#1f77b4', 'color_clean': '#1f77b4',
            'df_override': df_c,      # 使用独立数据集
            'consistent_group': False,  # 不参与 A∪B 并集
        },
    ]
    
    results_partial = []
    
    # ====================================================================
    # --consistent-outliers: 预扫描, 收集A/B各自离群点后取并集
    #   注意: 预期C 有独立数据集, 不参与 A∪B 并集
    # ====================================================================
    unified_outliers = None  # None = 独立模式; list = 一致模式
    if args.consistent_outliers:
        _pre_outliers = {}  # label -> outlier list
        for a in analyses:
            if a.get('consistent_group') is False:
                continue  # C 不参与并集
            cols_pre = ['Structure', a['adh_col'], a['temp_col'], a['size_col']]
            v_pre = df[cols_pre].dropna().copy()
            if len(v_pre) < 5:
                continue
            v_pre, _ = _calc_residuals(v_pre, a['adh_col'], a['temp_col'], a['size_col'])
            out_pre, _, auto_pre, _ = _detect_outliers(
                v_pre, threshold=args.threshold, manual_exclude=args.exclude)
            _pre_outliers[a['label']] = out_pre
        
        # 取并集 (A ∪ B)
        all_sets = [set(v) for v in _pre_outliers.values()]
        unified_outliers = list(dict.fromkeys(
            s for st in all_sets for s in st))  # 保序去重
        
        print(f"  ⚙ --consistent-outliers 模式: A/B 使用统一排除集 (C独立)")
        for lbl, outs in _pre_outliers.items():
            print(f"    {lbl} 单独检测: {', '.join(outs) if outs else '无'}")
        print(f"    并集 ({len(unified_outliers)}个): {', '.join(unified_outliers) if unified_outliers else '无'}")
        print()
    
    for a in analyses:
        # 选择数据源: df_override 或默认 df
        src_df = a.get('df_override', df)
        use_consistent = (unified_outliers is not None
                          and a.get('consistent_group') is not False)

        cols = ['Structure', a['adh_col'], a['temp_col'], a['size_col']]
        # 添加可能需要的额外尺寸列
        for col_extra in ['n_PtSn', 'n_SnO', 'nPt', 'nSn', 'nO']:
            if col_extra in src_df.columns and col_extra not in cols:
                cols.append(col_extra)
        
        v = src_df[cols].dropna().copy()
        if len(v) < 5:
            print(f"  {a['label']}: 数据不足 (n={len(v)}), 跳过")
            continue
        
        # 简单相关 (全样本)
        r_simple, p_simple = stats.pearsonr(v[a['adh_col']], v[a['temp_col']])
        
        # ---- 确定排除列表 ----
        if use_consistent:
            outliers = [s for s in unified_outliers if s in v['Structure'].values]
            auto_out = list(outliers)
            manual_out = []
            outlier_mode = 'consistent'
        else:
            # 独立模式: 先在全样本上计算残差做自动检测
            v_pre, _ = _calc_residuals(v, a['adh_col'], a['temp_col'], a['size_col'])
            outliers, v_pre, auto_out, manual_out = _detect_outliers(
                v_pre, threshold=args.threshold, manual_exclude=args.exclude)
            outlier_mode = 'independent'
        
        # ---- 先排除, 再做所有残差分析 ----
        v_clean = v[~v['Structure'].isin(outliers)].copy()
        
        if len(v_clean) < 5:
            print(f"  {a['label']}: 去离群后数据不足 (n={len(v_clean)}), 跳过")
            continue
        
        # 偏相关 (全样本, 仅用于对比报告)
        v_full_res, slopes_full = _calc_residuals(v, a['adh_col'], a['temp_col'], a['size_col'])
        r_partial, p_partial = stats.pearsonr(v_full_res['adh_res'], v_full_res['T_res'])
        
        # 偏相关 (clean 样本 — 这是最终使用的结果)
        v_clean, slopes_clean = _calc_residuals(v_clean, a['adh_col'], a['temp_col'], a['size_col'])
        r_clean, p_clean = stats.pearsonr(v_clean['adh_res'], v_clean['T_res'])
        
        # 在 clean 数据上计算 dist_z (用于打印表排序)
        _out_dummy, v_clean, _, _ = _detect_outliers(
            v_clean, threshold=999, manual_exclude=None)
        
        # ---- 辅助: 获取结构的 (nPt,nSn,nO) 标签 ----
        def _comp_label(struct_name):
            """返回 '(nPt,nSn,nO)' 或 '(nPt,nSn)' 字符串"""
            # 先在 v_clean 中找, 再在 v 中找
            for df_src in [v_clean, v]:
                r = df_src[df_src['Structure'] == struct_name]
                if len(r) > 0 and 'nPt' in r.columns:
                    r = r.iloc[0]
                    if 'nO' in r.index and pd.notna(r.get('nO', np.nan)):
                        return f"({int(r['nPt'])},{int(r['nSn'])},{int(r['nO'])})"
                    else:
                        return f"({int(r['nPt'])},{int(r['nSn'])})"
            return ""
        
        # 打印结果
        mode_tag = " [A∪B一致]" if outlier_mode == 'consistent' else ""
        print(f"  ┌─── {a['label']}: {a['adh_name']} vs {a['temp_name']}{mode_tag} ───┐")
        print(f"  │ 物理预期: {a['physical']}")
        print(f"  │ 离群阈值: {args.threshold}σ{mode_tag}")
        print(f"  │")
        print(f"  │ 简单相关(全={len(v)}点, 已含手动排除):       r = {r_simple:+.3f}, p = {p_simple:.3e} {_sig_label(p_simple)}")
        print(f"  │ 偏相关(全={len(v)}点, 控制{a['size_col']}): r = {r_partial:+.3f}, p = {p_partial:.3e} {_sig_label(p_partial)}")
        
        # 离群点详细信息
        if outliers:
            if outlier_mode == 'consistent':
                out_info = [f"{s} {_comp_label(s)}" for s in outliers]
                print(f"  │ 统一排除 ({len(outliers)}个, A∪B并集): {', '.join(out_info)}")
            else:
                auto_info = [f"{s} {_comp_label(s)}" for s in auto_out]
                print(f"  │ 自动排除 ({len(auto_out)}个, >{args.threshold}σ): {', '.join(auto_info) if auto_info else '无'}")
                if manual_out:
                    manual_info = [f"{s} {_comp_label(s)}" for s in manual_out]
                    print(f"  │ 手动排除 ({len(manual_out)}个, --exclude): {', '.join(manual_info)}")
                print(f"  │ 排除合计 ({len(outliers)}个):  {', '.join(outliers)}")
        else:
            print(f"  │ 离群点: 无 (阈值 {args.threshold}σ)")
        
        print(f"  │ clean={len(v_clean)}点 (在全基础上再去自动离群) 偏相关:  r = {r_clean:+.3f}, p = {p_clean:.3e} {_sig_label(p_clean)}")
        print(f"  │")
        # ---- 单变量 R²: adh-only 和 size-only ----
        from sklearn.linear_model import LinearRegression as _LR
        _x_adh  = v_clean[a['adh_col']].values
        _x_size = v_clean[a['size_col']].values
        _y_T    = v_clean[a['temp_col']].values
        _R2_adh_only  = _LR().fit(_x_adh.reshape(-1,1),  _y_T).score(_x_adh.reshape(-1,1),  _y_T)
        _R2_size_only = _LR().fit(_x_size.reshape(-1,1), _y_T).score(_x_size.reshape(-1,1), _y_T)
        # 单变量 Pearson r（有符号）及 p 值
        _r_adh_univ,  _p_adh_univ  = stats.pearsonr(_x_adh,  _y_T)
        _r_size_univ, _p_size_univ = stats.pearsonr(_x_size, _y_T)
        print(f"  │ 单变量 R²(adh only)  = {_R2_adh_only:.3f}  (r = {_r_adh_univ:+.3f} {_sig_label(_p_adh_univ)})")
        print(f"  │ 单变量 R²(size only) = {_R2_size_only:.3f}  (r = {_r_size_univ:+.3f} {_sig_label(_p_size_univ)})")
        print(f"  │")
        # ---- 多元回归: T = a*adh + b*size + c ----
        _X_mr = np.column_stack([_x_adh, _x_size])
        _reg_mr = _LR().fit(_X_mr, _y_T)
        _R2_mr = _reg_mr.score(_X_mr, _y_T)
        _coef_adh, _coef_size = _reg_mr.coef_
        _intercept = _reg_mr.intercept_
        # t-values & partial r for both predictors
        _n_mr = len(_y_T); _k_mr = 2
        _resid_mr = _y_T - _reg_mr.predict(_X_mr)
        _MSE_mr = np.sum(_resid_mr**2) / (_n_mr - _k_mr - 1)
        _X_c = np.column_stack([np.ones(_n_mr), _X_mr])
        _XtX_inv = np.linalg.inv(_X_c.T @ _X_c)
        _se_mr = np.sqrt(_MSE_mr * np.diag(_XtX_inv))[1:]
        _t_mr = _reg_mr.coef_ / _se_mr
        _df_mr = _n_mr - _k_mr - 1
        _rp_adh  = _t_mr[0] / np.sqrt(_t_mr[0]**2 + _df_mr)
        _rp_size = _t_mr[1] / np.sqrt(_t_mr[1]**2 + _df_mr)
        from scipy.stats import t as _t_dist
        _p_adh  = 2 * (1 - _t_dist.cdf(abs(_t_mr[0]), df=_df_mr))
        _p_size = 2 * (1 - _t_dist.cdf(abs(_t_mr[1]), df=_df_mr))
        # 标准化 β
        from sklearn.preprocessing import StandardScaler as _SS
        _X_std = _SS().fit_transform(_X_mr)
        _y_std = (_y_T - _y_T.mean()) / _y_T.std()
        _reg_s = _LR(fit_intercept=False).fit(_X_std, _y_std)
        _beta_adh, _beta_size = _reg_s.coef_
        print(f"  │ 多元回归(clean {len(v_clean)}点):")
        print(f"  │   {a['temp_name']} = {_coef_adh:+.1f}·{a['adh_name']} {_coef_size:+.1f}·{a['size_col']} {_intercept:+.1f}")
        print(f"  │   R² = {_R2_mr:.3f}")
        print(f"  │   r_partial(adh)  = {_rp_adh:+.3f}, p = {_p_adh:.3e} {_sig_label(_p_adh)},  β = {_beta_adh:+.3f}")
        print(f"  │   r_partial(size) = {_rp_size:+.3f}, p = {_p_size:.3e} {_sig_label(_p_size)},  β = {_beta_size:+.3f}")
        print(f"  └{'─' * 70}┘")
        
        # 逐点数据 — 基于 clean 子集的残差 (与 Ac/Ad 图完全一致)
        print(f"\n  {'结构':>12} │ (Pt,Sn,O) │ {a['size_col']:>6} │ {a['adh_col']:>12} │ adh残差 │ {a['temp_col']:>12} │  T残差 │ 距离z")
        print(f"  {'─' * 100}")
        for _, row in v_clean.sort_values('dist_z', ascending=False).iterrows():
            comp = _comp_label(row['Structure'])
            print(f"  {row['Structure']:>12} │ {comp:>9} │ {row[a['size_col']]:>6.0f} │ {row[a['adh_col']]:>12.2f} │ {row['adh_res']:>+7.2f} │ {row[a['temp_col']]:>12.0f} │ {row['T_res']:>+6.0f} │ {row['dist_z']:>5.2f}")
        
        # 单独列出被排除的结构
        if outliers:
            print(f"\n  已排除 ({len(outliers)}个):")
            for s in outliers:
                comp = _comp_label(s)
                r_excl = v[v['Structure'] == s]
                if len(r_excl) > 0:
                    r_excl = r_excl.iloc[0]
                    print(f"  {'✗ ' + s:>14} │ {comp:>9} │ {r_excl[a['size_col']]:>6.0f} │ {r_excl[a['adh_col']]:>12.2f} │    --- │ {r_excl[a['temp_col']]:>12.0f} │   --- │  ---")
        print()
        
        # 保存结果用于返回
        results_partial.append({
            'label': a['label'], 'adh': a['adh_name'], 'temp': a['temp_name'],
            'r_simple': r_simple, 'p_simple': p_simple,
            'r_partial': r_partial, 'p_partial': p_partial,
            'outliers': outliers, 'n_outliers': len(outliers),
            'r_clean': r_clean, 'p_clean': p_clean,
            'n_total': len(v), 'n_clean': len(v_clean),
            'v_all': v_full_res, 'v_clean': v_clean, 'analysis': a,
            # 多元回归信息
            'R2_mr': _R2_mr,
            'R2_adh_only': _R2_adh_only, 'R2_size_only': _R2_size_only,
            'r_adh_univ': _r_adh_univ, 'p_adh_univ': _p_adh_univ,
            'r_size_univ': _r_size_univ, 'p_size_univ': _p_size_univ,
            'coef_adh': _coef_adh, 'coef_size': _coef_size, 'intercept': _intercept,
            'r_partial_adh': _rp_adh, 'p_partial_adh': _p_adh, 'beta_adh': _beta_adh,
            'r_partial_size': _rp_size, 'p_partial_size': _p_size, 'beta_size': _beta_size,
        })
    
    # 总结对比表
    print("  " + "=" * 120)
    print("  偏相关分析总结")
    print("  " + "=" * 120)
    print("  口径说明: 全 = 当前可用全集(已包含手动排除如 --exclude/--exclude-c); clean = 全基础上再做自动离群清洗")
    print(f"  {'关系':<28} │ {'简单r(全,手动后)':>14} │ {'偏r(全,手动后)':>14} │ {'偏r(clean,再去离群)':>18} │ {'R²':>5} │ {'r_adh':>7} │ {'r_size':>7} │ {'β_adh':>6} │ {'β_size':>6} │ 结论")
    print(f"  {'─' * 120}")
    for res in results_partial:
        s1 = f"{res['r_simple']:+.3f}{_sig_label(res['p_simple']):>3}"
        s2 = f"{res['r_partial']:+.3f}{_sig_label(res['p_partial']):>3}"
        s3 = f"{res['r_clean']:+.3f}{_sig_label(res['p_clean']):>3}"
        s_R2 = f"{res['R2_mr']:.2f}"
        s_ra = f"{res['r_partial_adh']:+.3f}"
        s_rs = f"{res['r_partial_size']:+.3f}"
        s_ba = f"{res['beta_adh']:+.2f}"
        s_bs = f"{res['beta_size']:+.2f}"
        # 判定
        if _sig_label(res['p_simple']) != 'ns' and _sig_label(res['p_partial']) == 'ns':
            conclusion = "尺寸驱动的虚假显著"
        elif _sig_label(res['p_simple']) == 'ns' and _sig_label(res['p_partial']) != 'ns':
            conclusion = "被尺寸遮盖的真实关系"
        elif _sig_label(res['p_clean']) != 'ns' and _sig_label(res['p_partial']) == 'ns':
            conclusion = "去离群后显著"
        else:
            conclusion = "一致"
        print(f"  {res['label']+': '+res['adh']+' → '+res['temp']:<28} │ {s1:>10} │ {s2:>10} │ {s3:>10} │ {s_R2:>5} │ {s_ra:>7} │ {s_rs:>7} │ {s_ba:>6} │ {s_bs:>6} │ {conclusion}")
    print()
    
    return results_partial


def _confidence_band(x, y, x_plot, confidence=0.95):
    """计算线性回归的置信区间带"""
    from scipy.stats import t as t_dist
    slope, intercept, _, _, _ = stats.linregress(x, y)
    y_pred = slope * x + intercept
    residuals = y - y_pred
    n = len(x)
    df = n - 2
    se = np.sqrt(np.sum(residuals**2) / df)
    t_val = t_dist.ppf((1 + confidence) / 2, df)
    x_mean = np.mean(x)
    x_var = np.sum((x - x_mean)**2)
    y_plot = slope * x_plot + intercept
    ci = t_val * se * np.sqrt(1/n + (x_plot - x_mean)**2 / x_var)
    return y_plot, ci


def _style_ax(ax):
    """统一坐标轴样式: 与 plot_eadh_tm_correlation.py 一致"""
    ax.tick_params(axis='both', direction='out', length=6, width=2, labelsize=18)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(2)
    ax.patch.set_alpha(0)


def plot_size_deconvolution(df, results_partial, output_dir="results/adhesion_analysis"):
    """
    绘制 deconvolution 拆分图 — 对每个假设生成两张 3-panel 图：
      *_adh_deconvolution.png  — 以粘附能为 X 轴, 控制 size
      *_size_deconvolution.png — 以尺寸变量为 X 轴, 控制 adhesion

    每张图 3 列:
      (a) 简单散点 + 置信区间
      (b) 着色散点 (颜色=被控制的另一变量)
      (c) 偏相关残差散点 + 置信区间

    共 3 假设 × 2 视角 = 6 张图.
    """
    import os
    from scipy.stats import linregress as _lr
    os.makedirs(output_dir, exist_ok=True)

    for res in results_partial:
        a   = res['analysis']
        hk  = a['label'].replace('预期', 'hypothesis')   # e.g. 'hypothesisA'

        for perspective in ('adh', 'size'):
            # ---- 每次循环重新拷贝, 避免残留 ----
            v_all   = res['v_all'].copy()
            v_clean = res['v_clean'].copy()
            outliers = res['outliers']

            # ---- 选择 X 轴 / 颜色 / 残差方向 ----
            if perspective == 'adh':
                # X=粘附能, 颜色=size, 残差去 size
                x_col     = a['adh_col']
                color_col = a['size_col']
                ctrl_col  = a['size_col']        # 被控制 (被去除) 的变量
            else:
                # X=size, 颜色=粘附能, 残差去 adh
                x_col     = a['size_col']
                color_col = a['adh_col']
                ctrl_col  = a['adh_col']

            # ---- 英文标签 ----
            def _adh_lbl(col):
                if 'Type2' in col:
                    return r'$E_{\mathrm{adh}}^{\mathrm{Type2}}$/atom'
                elif 'Type3' in col:
                    return r'$E_{\mathrm{adh}}^{\mathrm{Type3}}$/atom'
                elif 'Eadh' in col:
                    return r'$E_{\mathrm{adh}}$/atom'
                return col

            def _size_lbl(col):
                if col == 'n_SnO':
                    return r'$n_{\mathrm{SnO}}$'
                elif col == 'n_PtSn':
                    return r'$n_{\mathrm{PtSn}}$'
                elif col == 'nMetal':
                    return r'$n_{\mathrm{Metal}}$'
                return col

            def _col_label(col, unit=True):
                """为任意列生成美观标签"""
                if 'Type' in col or 'Eadh' in col:
                    lbl = _adh_lbl(col)
                    return f'{lbl} (eV/atom)' if unit else lbl
                else:
                    return _size_lbl(col)

            x_label     = _col_label(x_col)
            color_label = _col_label(color_col)
            temp_label  = (r'$T_{\mathrm{m}}$' if 'lindemann' in a['temp_col']
                           else r'$T_{\mathrm{onset,O}}$')

            # ---- 离群掩码 ----
            is_outlier = v_all['Structure'].isin(outliers)
            normal_all = v_all[~is_outlier]
            out_all    = v_all[is_outlier]

            # ---- 标签生成 ----
            def _point_label(row):
                if 'nPt' in row.index:
                    if 'nO' in row.index and pd.notna(row.get('nO', np.nan)):
                        return f"({int(row['nPt'])},{int(row['nSn'])},{int(row['nO'])})"
                    else:
                        return f"({int(row['nPt'])},{int(row['nSn'])})"
                return row['Structure']

            show_labels = not args.no_labels

            # ---- 简单相关 ----
            r_all_s, p_all_s = (stats.pearsonr(v_all[x_col].astype(float),
                                               v_all[a['temp_col']].astype(float))
                                if len(v_all) >= 5 else (np.nan, np.nan))
            r_cln_s, p_cln_s = (stats.pearsonr(normal_all[x_col].astype(float),
                                               normal_all[a['temp_col']].astype(float))
                                if len(normal_all) >= 5 else (np.nan, np.nan))

            # ---- 残差计算: 从 x_col 和 temp 中去掉 ctrl_col 的线性影响 ----
            # v_all
            sl_x, it_x, *_ = _lr(v_all[ctrl_col].astype(float), v_all[x_col].astype(float))
            sl_t, it_t, *_ = _lr(v_all[ctrl_col].astype(float), v_all[a['temp_col']].astype(float))
            v_all['_xres'] = v_all[x_col] - (sl_x * v_all[ctrl_col] + it_x)
            v_all['_yres'] = v_all[a['temp_col']] - (sl_t * v_all[ctrl_col] + it_t)
            # v_clean
            sl_x2, it_x2, *_ = _lr(v_clean[ctrl_col].astype(float), v_clean[x_col].astype(float))
            sl_t2, it_t2, *_ = _lr(v_clean[ctrl_col].astype(float), v_clean[a['temp_col']].astype(float))
            v_clean['_xres'] = v_clean[x_col] - (sl_x2 * v_clean[ctrl_col] + it_x2)
            v_clean['_yres'] = v_clean[a['temp_col']] - (sl_t2 * v_clean[ctrl_col] + it_t2)

            # 偏相关 r
            r_cln_c, p_cln_c = stats.pearsonr(v_clean['_xres'], v_clean['_yres'])
            r_all_c, p_all_c = stats.pearsonr(v_all['_xres'],   v_all['_yres'])

            # 重建 outlier masks (v_all 已被修改)
            is_outlier = v_all['Structure'].isin(outliers)
            normal_all = v_all[~is_outlier]
            out_all    = v_all[is_outlier]

            # ================================================================
            fig, axes = plt.subplots(1, 3, figsize=(30, 9))
            fig.patch.set_alpha(0)

            # ================================================================
            # (a) 简单散点 + 置信区间
            # ================================================================
            ax = axes[0]
            if len(normal_all) >= 5:
                xa = normal_all[x_col].values.astype(float)
                ya = normal_all[a['temp_col']].values.astype(float)
                xp = np.linspace(xa.min(), xa.max(), 100)
                yf, ci = _confidence_band(xa, ya, xp)
                ax.fill_between(xp, yf - ci, yf + ci, alpha=0.18, color=a['color'])
                ax.plot(xp, yf, '-', color=a['color'], lw=4)

            # 全数据拟合线
            z = np.polyfit(v_all[x_col].astype(float), v_all[a['temp_col']].astype(float), 1)
            pf = np.poly1d(z)
            xl = np.linspace(v_all[x_col].min(), v_all[x_col].max(), 100)
            ax.plot(xl, pf(xl), ':', color='gray', lw=2, alpha=0.5)

            ax.scatter(normal_all[x_col], normal_all[a['temp_col']],
                       s=200, c=a['color'], edgecolors='black', linewidths=2,
                       alpha=0.85, zorder=5)
            if len(out_all) > 0:
                ax.scatter(out_all[x_col], out_all[a['temp_col']],
                           s=250, c='red', marker='X', edgecolors='darkred',
                           linewidths=2, zorder=6, alpha=0.9)

            if show_labels:
                for _, row in v_all.iterrows():
                    is_out = row['Structure'] in outliers
                    ax.annotate(_point_label(row),
                               (row[x_col], row[a['temp_col']]),
                               fontsize=11, ha='center', va='bottom',
                               xytext=(0, 10), textcoords='offset points',
                               fontweight='bold' if is_out else 'normal',
                               color='red' if is_out else 'black', alpha=0.8,
                               arrowprops=dict(arrowstyle='-', color='gray',
                                               lw=0.8, alpha=0.4))

            ax.text(0.05, 0.95,
                    f'All (n={len(v_all)}):   r = {r_all_s:+.3f} ({_sig_label(p_all_s)})\n'
                    f'Clean (n={len(normal_all)}): r = {r_cln_s:+.3f} ({_sig_label(p_cln_s)})',
                    transform=ax.transAxes, fontsize=16, va='top', ha='left',
                    bbox=dict(boxstyle='round,pad=0.4', fc='white', ec='gray', alpha=0.9))

            ax.set_xlabel(x_label, fontsize=22)
            ax.set_ylabel(f'{temp_label} (K)', fontsize=22)
            ax.set_title('(a) Simple correlation', fontsize=22, fontweight='bold', loc='left')
            _style_ax(ax)

            # ================================================================
            # (b) 着色散点  — 颜色 = 被控制的变量
            # ================================================================
            ax = axes[1]
            # 选 cmap: 粘附能→ RdYlBu, 尺寸 → RdYlGn_r
            cmap_b = 'RdYlBu' if ('Type' in color_col or 'Eadh' in color_col) else 'RdYlGn_r'

            scatter = ax.scatter(normal_all[x_col], normal_all[a['temp_col']], s=200,
                                 c=normal_all[color_col], cmap=cmap_b, edgecolors='black',
                                 linewidths=2, zorder=5, alpha=0.85,
                                 vmin=v_all[color_col].min(), vmax=v_all[color_col].max())
            if len(out_all) > 0:
                ax.scatter(out_all[x_col], out_all[a['temp_col']], s=250,
                           c='red', marker='X', edgecolors='darkred', linewidths=2,
                           zorder=6, alpha=0.9)
            cbar = plt.colorbar(scatter, ax=ax, shrink=0.85, pad=0.02)
            cbar.set_label(color_label, fontsize=20)
            cbar.ax.tick_params(labelsize=16)

            if show_labels:
                for _, row in v_all.iterrows():
                    is_out = row['Structure'] in outliers
                    ax.annotate(_point_label(row),
                               (row[x_col], row[a['temp_col']]),
                               fontsize=11, ha='center', va='bottom',
                               xytext=(0, 10), textcoords='offset points',
                               fontweight='bold' if is_out else 'normal',
                               color='red' if is_out else 'black', alpha=0.8,
                               arrowprops=dict(arrowstyle='-', color='gray',
                                               lw=0.8, alpha=0.4))

            # 混淆相关文本
            ctrl_label_short = _col_label(ctrl_col, unit=False)
            r_ctrl_x, _ = stats.pearsonr(normal_all[ctrl_col].astype(float),
                                         normal_all[x_col].astype(float))
            r_ctrl_T, _ = stats.pearsonr(normal_all[ctrl_col].astype(float),
                                         normal_all[a['temp_col']].astype(float))
            confound_name = 'Size' if perspective == 'adh' else 'Adhesion'
            ax.text(0.05, 0.95,
                    f'{confound_name} confounding:\n'
                    f'{ctrl_label_short} → {_col_label(x_col, unit=False)}: r = {r_ctrl_x:+.2f}\n'
                    f'{ctrl_label_short} → {temp_label}: r = {r_ctrl_T:+.2f}',
                    transform=ax.transAxes, fontsize=14, va='top', ha='left',
                    bbox=dict(boxstyle='round,pad=0.4', fc='lightyellow', ec='orange', alpha=0.9))

            ax.set_xlabel(x_label, fontsize=22)
            ax.set_ylabel(f'{temp_label} (K)', fontsize=22)
            title_b = f'(b) Colored by {_col_label(color_col, unit=False)}'
            ax.set_title(title_b, fontsize=22, fontweight='bold', loc='left')
            _style_ax(ax)

            # ================================================================
            # (c) 偏相关残差散点 + 置信区间
            # ================================================================
            ax = axes[2]
            normal = v_all[~is_outlier]
            out_pts = v_all[is_outlier]

            if len(v_clean) >= 5:
                x_arr = v_clean['_xres'].values
                y_arr = v_clean['_yres'].values
                xp = np.linspace(x_arr.min(), x_arr.max(), 100)
                yf, ci = _confidence_band(x_arr, y_arr, xp)
                ax.fill_between(xp, yf - ci, yf + ci, alpha=0.18, color=a['color_clean'])
                ax.plot(xp, yf, '-', color=a['color_clean'], lw=4,
                        label=f'Clean (n={len(v_clean)}): r = {r_cln_c:+.3f} {_sig_label(p_cln_c)}')

            z_c = np.polyfit(v_all['_xres'], v_all['_yres'], 1)
            pf_c = np.poly1d(z_c)
            xl_c = np.linspace(v_all['_xres'].min(), v_all['_xres'].max(), 100)
            ax.plot(xl_c, pf_c(xl_c), ':', color='gray', lw=2, alpha=0.5,
                    label=f'All (n={len(v_all)}): r = {r_all_c:+.3f} {_sig_label(p_all_c)}')

            ax.scatter(normal['_xres'], normal['_yres'], s=200, c=a['color_clean'],
                       edgecolors='black', linewidths=2, zorder=5, alpha=0.85)
            if len(out_pts) > 0:
                ax.scatter(out_pts['_xres'], out_pts['_yres'], s=250, c='red',
                           marker='X', edgecolors='darkred', linewidths=2, zorder=6, alpha=0.9,
                           label=f'Outliers (n={len(out_pts)})')

            if show_labels:
                for _, row in v_all.iterrows():
                    is_out = row['Structure'] in outliers
                    ax.annotate(_point_label(row),
                               (row['_xres'], row['_yres']),
                               fontsize=11, ha='center', va='bottom',
                               xytext=(0, 10), textcoords='offset points',
                               fontweight='bold' if is_out else 'normal',
                               color='red' if is_out else 'black', alpha=0.8,
                               arrowprops=dict(arrowstyle='-', color='gray',
                                               lw=0.8, alpha=0.4))

            ax.axhline(y=0, color='gray', ls=':', alpha=0.4, lw=1)
            ax.axvline(x=0, color='gray', ls=':', alpha=0.4, lw=1)

            ctrl_what = _col_label(ctrl_col, unit=False)
            c_xlabel = f'{_col_label(x_col, unit=False)} residual'
            if perspective == 'adh':
                c_xlabel += ' (eV/atom)'
            ax.set_xlabel(c_xlabel, fontsize=22)
            ax.set_ylabel(f'{temp_label} residual (K)', fontsize=22)
            corrected_by = 'size' if perspective == 'adh' else 'adhesion'
            ax.set_title(f'(c) Partial correlation ({corrected_by}-corrected)',
                         fontsize=22, fontweight='bold', loc='left')
            ax.legend(loc='lower right', fontsize=14, frameon=True, fancybox=True,
                      framealpha=0.9, edgecolor='gray')
            _style_ax(ax)

            # ---- 保存 ----
            plt.tight_layout(w_pad=4)
            fname = f'{hk}_{perspective}_deconvolution.png'
            fig.savefig(f"{output_dir}/{fname}", dpi=300, bbox_inches='tight',
                        transparent=True)
            print(f"  [OK] {a['label']} {perspective} 拆分图: {output_dir}/{fname}")
            plt.close()


# ============================================================================
# --plot-clean: 出版级独立图 (去离群后, 每组3张)
# ============================================================================

def _nice_ticks(vmin, vmax, target_n=5):
    """Generate 4-7 nice (symmetric, integer-ish) tick values.
    
    优先选取 round step (1,2,5 × 10^k), 绝不使用 linspace 以避免
    产生 -6.2, 0.125 之类的丑陋小数.
    """
    import math
    data_range = vmax - vmin
    if data_range == 0:
        return np.array([vmin - 1, vmin, vmin + 1])
    raw_step = data_range / target_n
    mag = 10 ** math.floor(math.log10(raw_step))
    candidates = [1, 2, 5, 10]

    # --- 第一轮: 严格 4-7 个 ticks ---
    best_step = None
    best_n = 999
    for c in candidates:
        step = c * mag
        lo = math.floor(vmin / step) * step
        hi = math.ceil(vmax / step) * step
        n = round((hi - lo) / step) + 1
        if 4 <= n <= 7 and abs(n - target_n) < abs(best_n - target_n):
            best_step = step
            best_n = n

    # --- 第二轮: 放宽到 3-8 个 ticks ---
    if best_step is None:
        best_n = 999
        for c in candidates:
            step = c * mag
            lo = math.floor(vmin / step) * step
            hi = math.ceil(vmax / step) * step
            n = round((hi - lo) / step) + 1
            if 3 <= n <= 8 and abs(n - target_n) < abs(best_n - target_n):
                best_step = step
                best_n = n

    # --- 第三轮: 再放宽, 接受任意 >= 3 ---
    if best_step is None:
        best_n = 999
        for c in candidates:
            step = c * mag
            lo = math.floor(vmin / step) * step
            hi = math.ceil(vmax / step) * step
            n = round((hi - lo) / step) + 1
            if n >= 3 and abs(n - target_n) < abs(best_n - target_n):
                best_step = step
                best_n = n

    if best_step is None:
        # 终极 fallback: 用 raw_step (几乎不会到这里)
        best_step = raw_step
    lo = math.floor(vmin / best_step) * best_step
    hi = math.ceil(vmax / best_step) * best_step
    all_ticks = np.arange(lo, hi + best_step * 0.5, best_step)
    # 只保留数据范围附近的刻度 (用较宽的 margin 避免误删)
    margin = best_step * 0.6
    ticks = all_ticks[(all_ticks >= vmin - margin) & (all_ticks <= vmax + margin)]
    # 如果过滤太狠, 回退到全部 ticks (不用 linspace!)
    if len(ticks) < 3:
        ticks = all_ticks
    return ticks


def _fmt_tick(val):
    """Format tick value to 2 significant figures, prefer integer if close."""
    if val == 0:
        return '0'
    import math
    magnitude = math.floor(math.log10(abs(val)))
    if magnitude >= 1:  # >= 10, use integer
        rounded = round(val)
        return f'{rounded:g}'
    elif magnitude == 0:  # 1-9.x
        rounded = round(val, 1)
        return f'{rounded:g}'
    else:
        rounded = round(val, -magnitude + 1)
        return f'{rounded:g}'


def _style_clean_ax(ax, xlabel, ylabel, x_data, y_data,
                    user_xticks=None, user_yticks=None):
    """Apply strict publication style to a single 10x8 panel.
    
    与 plot_eadh_tm_correlation.py 完全一致:
      - 4 spines visible, lw=2
      - 仅左/下刻度线 (右/上无刻度线)
      - 刻度朝外, length=6, width=2
      - 可手动指定 xticks / yticks, 否则自动生成
    """

    # Font sizes — labels 34pt, NOT bold
    ax.set_xlabel(xlabel, fontsize=34, fontweight='normal', fontfamily='Arial')
    ax.set_ylabel(ylabel, fontsize=34, fontweight='normal', fontfamily='Arial')

    # Spines: all 4 visible, lw=2
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(2)

    # Ticks: outward, only left & bottom (no right/top ticks)
    ax.tick_params(axis='both', direction='out', length=6, width=2,
                   labelsize=28, top=False, right=False)

    # --- X ticks ---
    if user_xticks is not None:
        xt = np.array(user_xticks)
    else:
        xt = _nice_ticks(x_data.min(), x_data.max())
    ax.set_xticks(xt)
    ax.set_xticklabels([_fmt_tick(v) for v in xt], fontfamily='Arial')

    # --- Y ticks ---
    if user_yticks is not None:
        yt = np.array(user_yticks)
    else:
        yt = _nice_ticks(y_data.min(), y_data.max())
    ax.set_yticks(yt)
    ax.set_yticklabels([_fmt_tick(v) for v in yt], fontfamily='Arial')

    # Axis limits with padding
    x_pad = (x_data.max() - x_data.min()) * 0.12
    y_pad = (y_data.max() - y_data.min()) * 0.12
    ax.set_xlim(x_data.min() - x_pad, x_data.max() + x_pad)
    ax.set_ylim(y_data.min() - y_pad, y_data.max() + y_pad)

    # No grid, no title
    ax.grid(False)
    ax.set_title('')
    ax.patch.set_alpha(0)


def plot_clean_panels(df, results_partial, output_dir="results/adhesion_analysis"):
    """
    --plot-clean: 出版级独立图 (去离群点后的clean数据).

    每组关系4张独立图 (10x8 in):
      (a) 简单散点 + 拟合线 + 95% CI
      (b) 尺寸着色散点 (colorbar = nMetal)
      (c) 偏相关残差散点 + 拟合线 + 95% CI
      (d) 多元回归 T_predicted vs T_actual + y=x + β标注
    
    风格与 plot_eadh_tm_correlation.py 完全一致.
    """
    import os
    os.makedirs(output_dir, exist_ok=True)

    plt.rcParams.update({
        'font.family': 'Arial',
        'font.size': 28,
        'mathtext.default': 'regular',
    })

    # ---- 解析用户指定刻度 ----
    # 格式: --clean-xticks Aa:-0.6,-0.4,-0.2,0  Bc:100,200,300
    # 映射到 tick_map[('A','a','x')] = [-0.6, -0.4, -0.2, 0.0]
    def _parse_tick_specs(spec_list):
        """Parse 'Aa:v1,v2,...' specs into dict keyed by ('hypothesisX', panel_letter).
        支持短格式 'Cg' 和扩展格式 'Cg-adh'/'Cg-size'（g 面板两子图共用同一刻度）。
        """
        result = {}
        if spec_list is None:
            return result
        for spec in spec_list:
            if ':' not in spec:
                print(f"  [WARN] 无效刻度格式 (缺少':'): {spec}")
                continue
            key, vals = spec.split(':', 1)
            # 兼容扩展格式: Cg-adh / Cg-size → 取前两个字符
            key_base = key.split('-')[0]  # 去掉 -adh / -size 后缀
            if len(key_base) != 2 or key_base[0] not in 'ABC' or key_base[1] not in 'abcdefg':
                print(f"  [WARN] 无效面板标识 '{key}' (应为 Aa-Ag/Ba-Bg/.../Cg-adh)")
                continue
            try:
                tick_vals = [float(v) for v in vals.split(',')]
                # 转换为长格式 key，与绘图代码中 xt_map.get((hk, panel_letter)) 匹配
                hyp_long = f'hypothesis{key_base[0]}'   # 'C' → 'hypothesisC'
                panel_letter = key_base[1]               # 'g'
                result[(hyp_long, panel_letter)] = tick_vals
            except ValueError:
                print(f"  [WARN] 无法解析刻度值: {vals}")
        return result

    xt_map = _parse_tick_specs(args.clean_xticks)
    yt_map = _parse_tick_specs(args.clean_yticks)

    # hypothesis index: first = 'A', second = 'B', third = 'C'
    hyp_keys = ['A', 'B', 'C']

    # ---- 标签标注辅助函数 (模仿 plot_eadh_tm_correlation.py) ----
    # 解析用户偏移 --clean-offsets  格式: PANEL@Structure:dx,dy
    offset_map = {}  # (panel_key, structure) -> (dx, dy)

    # 1) 先从 JSON 文件加载历史偏移
    import json as _json, os as _os
    _offsets_path = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)),
                                  args.offsets_file)
    if _os.path.exists(_offsets_path):
        try:
            with open(_offsets_path, 'r', encoding='utf-8') as _f:
                _saved = _json.load(_f)
            for _k, _v in _saved.items():
                _panel, _struct = _k.split('@', 1)
                offset_map[(_panel, _struct)] = tuple(_v)
            print(f'  [INFO] 已从 {args.offsets_file} 加载 {len(_saved)} 条标签偏移')
        except Exception as _e:
            print(f'  [WARN] 读取 {args.offsets_file} 失败: {_e}')

    # 2) --clean-offsets 命令行参数覆盖文件中的同名键
    if args.clean_offsets:
        for spec in args.clean_offsets:
            try:
                left, coords = spec.split(':', 1)
                panel_part, struct = left.split('@', 1)
                dx, dy = [float(v) for v in coords.split(',')]
                offset_map[(panel_part, struct)] = (dx, dy)
            except Exception:
                print(f"  [WARN] 无法解析偏移: {spec}  (格式: PANEL@Structure:dx,dy)")

    def _auto_layout(points, x_range, y_range, n_iter=80):
        """Repulsion-based label placement to avoid overlaps.
        
        points: list of (x, y) data coordinates
        Returns: list of (dx, dy) offsets in data coordinates
        """
        n = len(points)
        if n == 0:
            return []

        # Approximate label bbox in data units
        lbl_w = x_range * 0.12   # estimated label width
        lbl_h = y_range * 0.055  # estimated label height

        # Initial offsets: spread around the point at 8 compass directions
        angles = [45, -45, 135, -135, 90, -90, 0, 180]
        r0_x = x_range * 0.06
        r0_y = y_range * 0.06
        import math
        offsets = []
        for i in range(n):
            ang = math.radians(angles[i % len(angles)])
            offsets.append([r0_x * math.cos(ang), r0_y * math.sin(ang)])

        # Iterative repulsion
        for iteration in range(n_iter):
            for i in range(n):
                fx, fy = 0.0, 0.0
                # label center position
                lx_i = points[i][0] + offsets[i][0]
                ly_i = points[i][1] + offsets[i][1]

                for j in range(n):
                    if i == j:
                        continue
                    lx_j = points[j][0] + offsets[j][0]
                    ly_j = points[j][1] + offsets[j][1]

                    # Overlap in bbox?
                    dx_ij = lx_i - lx_j
                    dy_ij = ly_i - ly_j
                    overlap_x = lbl_w - abs(dx_ij)
                    overlap_y = lbl_h - abs(dy_ij)

                    if overlap_x > 0 and overlap_y > 0:
                        # Push apart
                        push = 0.3
                        if dx_ij == 0 and dy_ij == 0:
                            dx_ij, dy_ij = 0.01 * x_range, 0.01 * y_range
                        fx += push * (1 if dx_ij >= 0 else -1) * overlap_x
                        fy += push * (1 if dy_ij >= 0 else -1) * overlap_y

                # Also repel from data points (avoid covering other dots)
                for j in range(n):
                    dx_p = lx_i - points[j][0]
                    dy_p = ly_i - points[j][1]
                    if abs(dx_p) < lbl_w * 0.5 and abs(dy_p) < lbl_h * 0.5:
                        push = 0.15
                        if dx_p == 0 and dy_p == 0:
                            dx_p, dy_p = 0.005 * x_range, 0.005 * y_range
                        fx += push * (1 if dx_p >= 0 else -1) * lbl_w * 0.3
                        fy += push * (1 if dy_p >= 0 else -1) * lbl_h * 0.3

                # Spring: don't stray too far from own point
                dist_x = offsets[i][0]
                dist_y = offsets[i][1]
                spring = 0.02
                fx -= spring * dist_x
                fy -= spring * dist_y

                offsets[i][0] += fx
                offsets[i][1] += fy

                # Clamp: stay within reasonable range
                max_dx = x_range * 0.25
                max_dy = y_range * 0.25
                offsets[i][0] = max(-max_dx, min(max_dx, offsets[i][0]))
                offsets[i][1] = max(-max_dy, min(max_dy, offsets[i][1]))

        return [(o[0], o[1]) for o in offsets]

    def _annotate_clean(ax, v, xcol, ycol, panel_key=None, x_override=None):
        """为每个数据点添加 (nPt,nSn,nO) 标签, 自动避重叠.
        
        Parameters:
            x_override: 如提供, 用此数组代替 v[xcol] 作为 x 坐标 (用于 Ad 面板).
        
        Returns: (annotations_list, structure_names_list) for interactive mode.
        """
        ann_list = []
        struct_names = []
        if args.no_labels:
            return ann_list, struct_names

        # x values: from override array or column
        if x_override is not None:
            x_vals = np.asarray(x_override)
        else:
            x_vals = v[xcol].values
        y_vals = v[ycol].values

        x_range = x_vals.max() - x_vals.min()
        y_range = y_vals.max() - y_vals.min()
        if x_range == 0:
            x_range = 1.0
        if y_range == 0:
            y_range = 1.0

        # Collect points and labels
        rows_data = []
        for i_row, (_, row) in enumerate(v.iterrows()):
            if 'nPt' in row.index:
                nPt = int(row['nPt']); nSn = int(row['nSn'])
                if 'nO' in row.index and pd.notna(row.get('nO', np.nan)):
                    lbl = f"({nPt},{nSn},{int(row['nO'])})"
                else:
                    lbl = f"({nPt},{nSn})"
            else:
                lbl = row['Structure']
            struct = row['Structure'] if 'Structure' in row.index else lbl
            rows_data.append((x_vals[i_row], y_vals[i_row], lbl, struct))

        points = [(r[0], r[1]) for r in rows_data]

        # Auto-layout offsets
        auto_offsets = _auto_layout(points, x_range, y_range)

        for i, (x, y, lbl, struct) in enumerate(rows_data):
            # Check user override first
            user_key = (panel_key, struct) if panel_key else None
            if user_key and user_key in offset_map:
                dx, dy = offset_map[user_key]
            else:
                dx, dy = auto_offsets[i] if i < len(auto_offsets) else (x_range * 0.04, y_range * 0.04)

            va = 'bottom' if dy > 0 else 'top'
            # 交互模式下关闭裁剪以便拖动查看超界标签; 非交互模式裁掉以保持图幅固定
            _clip = (args.interactive is None)
            ann = ax.annotate(lbl,
                              xy=(x, y),
                              xytext=(x + dx, y + dy),
                              fontsize=20, fontfamily='Arial',
                              ha='left', va=va, color='black',
                              clip_on=_clip,
                              bbox=dict(boxstyle='round,pad=0.5',
                                        facecolor='lightyellow' if args.interactive is not None else 'white',
                                        edgecolor='gray' if args.interactive is not None else 'none',
                                        alpha=0.7 if args.interactive is not None else 0),
                              arrowprops=dict(arrowstyle='-',
                                              connectionstyle='arc3,rad=0',
                                              color='gray', linewidth=1, alpha=0.5),
                              picker=True)
            ann_list.append(ann)
            struct_names.append(struct)

        return ann_list, struct_names

    # ---- panel filter ----
    panel_filter = None
    draw_panel_D = False          # 大写 D: 跨假设总结热图
    if args.clean_panel:
        panel_filter = set()
        for p in args.clean_panel:
            if p == 'D':
                draw_panel_D = True
            else:
                # 支持两种格式:
                #   短格式(原有): Aa / Bc / Cg  (单字母假说 + 单字母面板)
                #   扩展格式: B3g-adh / B'g-adh / B''g-adh (多字符假说 + g-adh/g-size)
                # 提取规则: 末尾 g-adh / g-size 或 单字母 a-g
                import re as _re
                _m = _re.match(r'^([A-C][^a-g]*)([a-g](?:-adh|-size)?)$', p)
                if _m:
                    panel_filter.add(p)
                elif len(p) == 2 and p[0] in 'ABC' and p[1] in 'abcdefg':
                    panel_filter.add(p)
                else:
                    print(f"  [WARN] 无效面板标识 '{p}' (应为 Aa-Ag/Ba-Bg/Ca-Cg/B3g-adh/D)")

    def _should_draw(hk, panel_letter):
        """判断是否需要绘制指定面板。
        hk: 'hypothesisA' / 'hypothesisB3' 等
        panel_letter: 'a'-'g'
        panel_filter 中的条目格式: 'Ag' / 'B3g-adh' 等（去掉 'hypothesis' 前缀）
        """
        if panel_filter is None:
            return True
        # hk_short: 去掉 'hypothesis' 前缀，如 'A' / "B'" / 'B3'
        hk_short = hk.replace('hypothesis', '')
        # 精确匹配: 短格式 'Ag', 或 g 面板扩展格式 'B3g-adh' / 'B3g-size' / 'B3g'
        if f'{hk_short}{panel_letter}' in panel_filter:
            return True
        # g 面板还接受 'B3g-adh' / 'B3g-size'
        if panel_letter == 'g':
            if f'{hk_short}g-adh' in panel_filter or f'{hk_short}g-size' in panel_filter:
                return True
        return False

    # ---- DraggableAnnotation for --interactive ----
    class DraggableAnnotation:
        def __init__(self, annotation, structure_name, panel_key):
            self.annotation = annotation
            self.structure_name = structure_name
            self.panel_key = panel_key
            self.dragging = False

        def connect(self):
            canvas = self.annotation.figure.canvas
            canvas.mpl_connect('button_press_event', self.on_press)
            canvas.mpl_connect('button_release_event', self.on_release)
            canvas.mpl_connect('motion_notify_event', self.on_motion)

        def _event_to_data(self, event):
            """将任意位置的鼠标事件转换为本 axes 的数据坐标.
            
            即使鼠标不在 axes 内 (event.inaxes is None), 也能通过
            display→data 变换获取正确坐标, 从而支持拖动到边框之外.
            """
            ax = self.annotation.axes
            try:
                return ax.transData.inverted().transform((event.x, event.y))
            except Exception:
                return None, None

        def on_press(self, event):
            contains, _ = self.annotation.contains(event)
            if not contains:
                return
            self.dragging = True
            xd, yd = self._event_to_data(event)
            self.x0 = xd
            self.y0 = yd
            print(f"  📍 选中: {self.structure_name}")

        def on_motion(self, event):
            if not self.dragging:
                return
            xd, yd = self._event_to_data(event)
            if xd is None or yd is None:
                return
            dx = xd - self.x0
            dy = yd - self.y0
            cur = self.annotation.get_position()
            self.annotation.set_position((cur[0] + dx, cur[1] + dy))
            self.x0 = xd
            self.y0 = yd
            self.annotation.figure.canvas.draw_idle()

        def on_release(self, event):
            if not self.dragging:
                return
            self.dragging = False
            xy = self.annotation.xy
            xytext = self.annotation.get_position()
            dx = xytext[0] - xy[0]
            dy = xytext[1] - xy[1]
            key = f'{self.panel_key}@{self.structure_name}'
            print(f"  ✅ {key}:{dx:.4f},{dy:.2f}")
            # 实时写回 JSON 文件
            offset_map[(self.panel_key, self.structure_name)] = (dx, dy)
            try:
                _data = {f'{p}@{s}': list(v) for (p, s), v in offset_map.items()}
                with open(_offsets_path, 'w', encoding='utf-8') as _f:
                    _json.dump(_data, _f, ensure_ascii=False, indent=2)
            except Exception as _e:
                print(f"  [WARN] 保存偏移失败: {_e}")

    # ---- 面板绘制总数统计 ----
    n_drawn = 0

    for idx, res in enumerate(results_partial):
        a = res['analysis']
        v_clean = res['v_clean'].copy()
        v_all = res['v_all']
        # hk 与 plot_size_deconvolution 保持一致：由 label 派生，而非固定列表索引
        hk = a['label'].replace('预期', 'hypothesis')   # e.g. 'hypothesisA', "hypothesisB'"

        # English axis labels (LaTeX)
        # 横坐标: A→E1adh, B→E2adh, C→Eadh
        # 纵坐标: A→T1, B→T2, C→Tm
        if 'Type2' in a['adh_col']:
            adh_label = r'$E_{adh}^1$ (eV/atom)'
            adh_res_label = r'$E_{adh}^1$ residual (eV/atom)'
        elif 'Type3' in a['adh_col']:
            adh_label = r'$E_{adh}^2$ (eV/atom)'
            adh_res_label = r'$E_{adh}^2$ residual (eV/atom)'
        elif 'Type4' in a['adh_col']:
            adh_label = r'$E_{adh}^4$ (eV/O)'
            adh_res_label = r'$E_{adh}^4$ residual (eV/O)'
        else:
            # 预期C: 通用 Eadh/atom
            adh_label = r'$E_{adh}$ (eV/atom)'
            adh_res_label = r'$E_{adh}$ residual (eV/atom)'

        if 'lindemann' in a['temp_col']:
            if 'Type2' in a['adh_col']:
                # A组: Type2 → T1
                temp_label = r'$T_1$ (K)'
                temp_res_label = r'$T_1$ residual (K)'
            else:
                # C组: Eadh → Tm
                temp_label = r'$T_m$ (K)'
                temp_res_label = r'$T_m$ residual (K)'
        elif 'T3' in a['temp_col']:
            # D组: Type4 → T3 (T_onset_O, 整体界面)
            temp_label = r'$T_3$ (K)'
            temp_res_label = r'$T_3$ residual (K)'
        else:
            # B组: Type3 → T2 (T_onset_O)
            temp_label = r'$T_2$ (K)'
            temp_res_label = r'$T_2$ residual (K)'

        # 用 analyses 配置里的自定义标签覆盖自动推断（优先级最高）
        adh_label      = a.get('axis_adh_label',      adh_label)
        adh_res_label  = a.get('axis_adh_res_label',  adh_res_label)
        temp_label     = a.get('axis_temp_label',     temp_label)
        temp_res_label = a.get('axis_temp_res_label', temp_res_label)

        htag = a['label'].replace('预期', 'hypothesis')
        # ============================================================
        # (a) Simple correlation — clean data only
        # ============================================================
        if _should_draw(hk, 'a'):
            fig, ax = plt.subplots(figsize=(10, 8))
            fig.patch.set_alpha(0)

            x_a = v_clean[a['adh_col']].values
            y_a = v_clean[a['temp_col']].values

            x_plot = np.linspace(x_a.min(), x_a.max(), 200)
            y_fit, ci = _confidence_band(x_a, y_a, x_plot)
            ax.fill_between(x_plot, y_fit - ci, y_fit + ci,
                            alpha=0.18, color=a['color'])
            r_val, p_val = stats.pearsonr(x_a, y_a)
            if args.no_stars:
                fit_label = f'r = {r_val:+.2f}'
            else:
                fit_label = f'r = {r_val:+.2f} ({_sig_label(p_val)})'
            ax.plot(x_plot, y_fit, '-', color=a['color'], lw=4,
                    label=fit_label)

            ax.scatter(x_a, y_a, s=200, c=a['color'], edgecolors='black',
                       linewidths=2, alpha=0.85, zorder=5)

            anns, names = _annotate_clean(ax, v_clean, a['adh_col'],
                                          a['temp_col'], panel_key=f'{hk}a')

            ax.legend(loc='best', fontsize=26, frameon=False,
                      prop={'family': 'Arial'})
            _style_clean_ax(ax, adh_label, temp_label, x_a, y_a,
                            user_xticks=xt_map.get((hk, 'a')),
                            user_yticks=yt_map.get((hk, 'a')))

            plt.tight_layout(pad=0.5)
            fname_a = f'{htag}_clean_a.png'

            if _is_interactive(f'{hk}a'):
                # Save initial version first
                fig.savefig(f'{output_dir}/{fname_a}', dpi=300,
                            bbox_inches='tight', transparent=True)
                _drags = []                       # keep refs → prevent GC
                for ann, sn in zip(anns, names):
                    da = DraggableAnnotation(ann, sn, f'{hk}a')
                    da.connect()
                    _drags.append(da)
                print(f'\n  🎯 交互模式: {hk}(a) — 拖动标签, 关闭窗口后自动保存')
                print(f'     复制 ✅ 行到 --clean-offsets 参数')
                plt.show()
                # Save again with adjusted positions (override bbox for clean output)
                for ann in anns:
                    ann.get_bbox_patch().set_alpha(0)
                    ann.get_bbox_patch().set_edgecolor('none')
                fig.savefig(f'{output_dir}/{fname_a}', dpi=300,
                            bbox_inches='tight', transparent=True)
            else:
                fig.savefig(f'{output_dir}/{fname_a}', dpi=300,
                            bbox_inches='tight', transparent=True)
            n_drawn += 1

            plt.close()
            print(f'  [OK] {htag} (a) simple clean: {output_dir}/{fname_a}')

        # ============================================================
        # (b) Colored by nMetal — clean data only
        # ============================================================
        if _should_draw(hk, 'b'):
            fig, ax = plt.subplots(figsize=(10, 8))
            fig.patch.set_alpha(0)

            x_b = v_clean[a['adh_col']].values
            y_b = v_clean[a['temp_col']].values
            c_b = v_clean[a['size_col']].values

            scatter = ax.scatter(x_b, y_b, s=200, c=c_b, cmap='RdYlGn_r',
                                 edgecolors='black', linewidths=2, zorder=5,
                                 alpha=0.85,
                                 vmin=v_all[a['size_col']].min(),
                                 vmax=v_all[a['size_col']].max())

            cbar = plt.colorbar(scatter, ax=ax, shrink=0.85, pad=0.02)
            cbar.set_label(r'$n_{\mathrm{Metal}}$', fontsize=28,
                           fontfamily='Arial')
            cbar.ax.tick_params(labelsize=22, direction='out', length=4, width=1.5)
            cbar.outline.set_linewidth(1.5)

            anns, names = _annotate_clean(ax, v_clean, a['adh_col'],
                                          a['temp_col'], panel_key=f'{hk}b')

            _style_clean_ax(ax, adh_label, temp_label, x_b, y_b,
                            user_xticks=xt_map.get((hk, 'b')),
                            user_yticks=yt_map.get((hk, 'b')))

            plt.tight_layout(pad=0.5)
            fname_b = f'{htag}_clean_b.png'

            if _is_interactive(f'{hk}b'):
                fig.savefig(f'{output_dir}/{fname_b}', dpi=300,
                            bbox_inches='tight', transparent=True)
                _drags = []                       # keep refs → prevent GC
                for ann, sn in zip(anns, names):
                    da = DraggableAnnotation(ann, sn, f'{hk}b')
                    da.connect()
                    _drags.append(da)
                print(f'\n  🎯 交互模式: {hk}(b) — 拖动标签, 关闭窗口后自动保存')
                print(f'     复制 ✅ 行到 --clean-offsets 参数')
                plt.show()
                for ann in anns:
                    ann.get_bbox_patch().set_alpha(0)
                    ann.get_bbox_patch().set_edgecolor('none')
                fig.savefig(f'{output_dir}/{fname_b}', dpi=300,
                            bbox_inches='tight', transparent=True)
            else:
                fig.savefig(f'{output_dir}/{fname_b}', dpi=300,
                            bbox_inches='tight', transparent=True)
            n_drawn += 1

            plt.close()
            print(f'  [OK] {htag} (b) size-colored clean: {output_dir}/{fname_b}')

        # ============================================================
        # (c) Partial correlation residuals — clean data only
        # ============================================================
        if _should_draw(hk, 'c'):
            fig, ax = plt.subplots(figsize=(10, 8))
            fig.patch.set_alpha(0)

            x_c = v_clean['adh_res'].values
            y_c = v_clean['T_res'].values

            x_plot = np.linspace(x_c.min(), x_c.max(), 200)
            y_fit, ci = _confidence_band(x_c, y_c, x_plot)
            ax.fill_between(x_plot, y_fit - ci, y_fit + ci,
                            alpha=0.18, color=a['color_clean'])
            r_c, p_c = stats.pearsonr(x_c, y_c)
            if args.no_stars:
                fit_label_c = f'r = {r_c:+.2f}'
            else:
                fit_label_c = f'r = {r_c:+.2f} ({_sig_label(p_c)})'
            ax.plot(x_plot, y_fit, '-', color=a['color_clean'], lw=4,
                    label=fit_label_c)

            ax.scatter(x_c, y_c, s=200, c=a['color_clean'], edgecolors='black',
                       linewidths=2, alpha=0.85, zorder=5)

            anns, names = _annotate_clean(ax, v_clean, 'adh_res', 'T_res',
                                          panel_key=f'{hk}c')

            ax.legend(loc='best', fontsize=26, frameon=False,
                      prop={'family': 'Arial'})
            _style_clean_ax(ax, adh_res_label, temp_res_label, x_c, y_c,
                            user_xticks=xt_map.get((hk, 'c')),
                            user_yticks=yt_map.get((hk, 'c')))

            plt.tight_layout(pad=0.5)
            fname_c = f'{htag}_clean_c.png'

            if _is_interactive(f'{hk}c'):
                fig.savefig(f'{output_dir}/{fname_c}', dpi=300,
                            bbox_inches='tight', transparent=True)
                _drags = []                       # keep refs → prevent GC
                for ann, sn in zip(anns, names):
                    da = DraggableAnnotation(ann, sn, f'{hk}c')
                    da.connect()
                    _drags.append(da)
                print(f'\n  🎯 交互模式: {hk}(c) — 拖动标签, 关闭窗口后自动保存')
                print(f'     复制 ✅ 行到 --clean-offsets 参数')
                plt.show()
                for ann in anns:
                    ann.get_bbox_patch().set_alpha(0)
                    ann.get_bbox_patch().set_edgecolor('none')
                fig.savefig(f'{output_dir}/{fname_c}', dpi=300,
                            bbox_inches='tight', transparent=True)
            else:
                fig.savefig(f'{output_dir}/{fname_c}', dpi=300,
                            bbox_inches='tight', transparent=True)
            n_drawn += 1

            plt.close()
            print(f'  [OK] {htag} (c) partial residual clean: {output_dir}/{fname_c}')

        # ============================================================
        # (d) Multivariate regression: predicted vs actual T
        # ============================================================
        if _should_draw(hk, 'd'):
            from sklearn.linear_model import LinearRegression
            from sklearn.preprocessing import StandardScaler

            fig, ax = plt.subplots(figsize=(10, 8))
            fig.patch.set_alpha(0)

            x_adh = v_clean[a['adh_col']].values
            x_size = v_clean[a['size_col']].values
            y_obs = v_clean[a['temp_col']].values

            # --- OLS multivariate regression ---
            X_raw = np.column_stack([x_adh, x_size])
            reg = LinearRegression().fit(X_raw, y_obs)
            y_pred = reg.predict(X_raw)
            R2 = reg.score(X_raw, y_obs)

            # --- Standardised β ---
            scaler = StandardScaler()
            X_std = scaler.fit_transform(X_raw)
            y_std = (y_obs - y_obs.mean()) / y_obs.std()
            reg_std = LinearRegression(fit_intercept=False).fit(X_std, y_std)
            beta_adh = reg_std.coef_[0]
            beta_size = reg_std.coef_[1]

            # --- p-values via t-test for each coefficient ---
            n_pts = len(y_obs)
            k_vars = 2
            residuals = y_obs - y_pred
            MSE = np.sum(residuals ** 2) / (n_pts - k_vars - 1)
            # 必须用包含截距列的设计矩阵计算 SE
            X_with_const = np.column_stack([np.ones(n_pts), X_raw])
            XtX_inv_full = np.linalg.inv(X_with_const.T @ X_with_const)
            # coef 顺序: [intercept, adh, size]
            se_full = np.sqrt(MSE * np.diag(XtX_inv_full))
            se_b = se_full[1:]  # 只取 adh, size 的 SE
            t_vals = reg.coef_ / se_b
            from scipy.stats import t as t_dist
            p_vals = 2 * (1 - t_dist.cdf(np.abs(t_vals), df=n_pts - k_vars - 1))

            # 从 t 统计量推算偏相关系数 (与 Ac 残差法数学等价)
            df_res = n_pts - k_vars - 1
            r_partial_adh  = t_vals[0] / np.sqrt(t_vals[0]**2 + df_res)
            r_partial_size = t_vals[1] / np.sqrt(t_vals[1]**2 + df_res)

            # --- Diagonal reference y = x ---
            all_vals = np.concatenate([y_obs, y_pred])
            lo, hi = all_vals.min(), all_vals.max()
            pad = (hi - lo) * 0.08
            diag = np.linspace(lo - pad, hi + pad, 50)
            ax.plot(diag, diag, '--', color='gray', lw=2, alpha=0.6,
                    zorder=1, label='y = x')

            # --- Scatter: predicted vs actual ---
            ax.scatter(y_pred, y_obs, s=200, c=a['color'], edgecolors='black',
                       linewidths=2, alpha=0.85, zorder=5)

            # --- Annotate data points ---
            anns, names = _annotate_clean(ax, v_clean, None, a['temp_col'],
                                          panel_key=f'{hk}d',
                                          x_override=y_pred)

            # --- Annotation box: β coefficients + R² ---
            # Determine nice names for the two predictors
            if 'Type2' in a['adh_col']:
                adh_sym = r'$E_{adh}^1$'
            elif 'Type3' in a['adh_col']:
                adh_sym = r'$E_{adh}^2$'
            else:
                adh_sym = r'$E_{adh}$'

            sig_adh = _sig_label(p_vals[0])
            sig_size = _sig_label(p_vals[1])

            if args.no_stars:
                info_text = (
                    f'$R^2$ = {R2:.2f}\n'
                    f'$r_{{adh}}$ = {r_partial_adh:+.2f}\n'
                    f'$r_{{size}}$ = {r_partial_size:+.2f}\n'
                    f'$\\beta_{{adh}}$ = {beta_adh:+.2f}\n'
                    f'$\\beta_{{size}}$ = {beta_size:+.2f}'
                )
            else:
                info_text = (
                    f'$R^2$ = {R2:.2f}\n'
                    f'$r_{{adh}}$ = {r_partial_adh:+.2f} ({sig_adh})\n'
                    f'$r_{{size}}$ = {r_partial_size:+.2f} ({sig_size})\n'
                    f'$\\beta_{{adh}}$ = {beta_adh:+.2f} ({sig_adh})\n'
                    f'$\\beta_{{size}}$ = {beta_size:+.2f} ({sig_size})'
                )
            ax.text(0.05, 0.95, info_text, transform=ax.transAxes,
                    fontsize=24, fontfamily='Arial', verticalalignment='top',
                    bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                              alpha=0.8, edgecolor='gray', linewidth=1))

            # --- Legend ---
            ax.legend(loc='lower right', fontsize=24, frameon=False,
                      prop={'family': 'Arial'})

            # --- Axis labels ---
            pred_label = temp_label.replace('(K)', 'predicted (K)')
            _style_clean_ax(ax, pred_label, temp_label, y_pred, y_obs,
                            user_xticks=xt_map.get((hk, 'd')),
                            user_yticks=yt_map.get((hk, 'd')))

            plt.tight_layout(pad=0.5)
            fname_d = f'{htag}_clean_d.png'

            if _is_interactive(f'{hk}d'):
                fig.savefig(f'{output_dir}/{fname_d}', dpi=300,
                            bbox_inches='tight', transparent=True)
                _drags = []
                for ann, sn in zip(anns, names):
                    da = DraggableAnnotation(ann, sn, f'{hk}d')
                    da.connect()
                    _drags.append(da)
                print(f'\n  🎯 交互模式: {hk}(d) — 拖动标签, 关闭窗口后自动保存')
                print(f'     复制 ✅ 行到 --clean-offsets 参数')
                plt.show()
                for ann in anns:
                    ann.get_bbox_patch().set_alpha(0)
                    ann.get_bbox_patch().set_edgecolor('none')
                fig.savefig(f'{output_dir}/{fname_d}', dpi=300,
                            bbox_inches='tight', transparent=True)
            else:
                fig.savefig(f'{output_dir}/{fname_d}', dpi=300,
                            bbox_inches='tight', transparent=True)
            n_drawn += 1

            plt.close()
            print(f'  [OK] {htag} (d) multivariate regression: {output_dir}/{fname_d}')

        # ============================================================
        # (e) Contour heatmap: predicted T on adh × size grid
        # ============================================================
        if _should_draw(hk, 'e'):
            from sklearn.linear_model import LinearRegression as _LR_e

            fig, ax = plt.subplots(figsize=(10, 8))
            fig.patch.set_alpha(0)

            x_adh_e = v_clean[a['adh_col']].values
            x_size_e = v_clean[a['size_col']].values
            y_T_e = v_clean[a['temp_col']].values

            # 多元回归
            X_e = np.column_stack([x_adh_e, x_size_e])
            reg_e = _LR_e().fit(X_e, y_T_e)
            R2_e = reg_e.score(X_e, y_T_e)

            # 创建网格
            adh_pad = (x_adh_e.max() - x_adh_e.min()) * 0.10
            size_pad = (x_size_e.max() - x_size_e.min()) * 0.15
            adh_range = np.linspace(x_adh_e.min() - adh_pad,
                                     x_adh_e.max() + adh_pad, 120)
            size_range = np.linspace(x_size_e.min() - size_pad,
                                      x_size_e.max() + size_pad, 120)
            adh_grid, size_grid = np.meshgrid(adh_range, size_range)
            T_grid = (reg_e.intercept_
                      + reg_e.coef_[0] * adh_grid
                      + reg_e.coef_[1] * size_grid)

            # 等高线填充
            contourf = ax.contourf(adh_grid, size_grid, T_grid,
                                   levels=15, cmap='RdYlBu_r', alpha=0.85)
            # 等高线
            contour = ax.contour(adh_grid, size_grid, T_grid,
                                  levels=15, colors='black',
                                  linewidths=0.6, alpha=0.35)
            ax.clabel(contour, inline=True, fontsize=12, fmt='%d')

            # 实测散点
            ax.scatter(x_adh_e, x_size_e,
                       c=y_T_e, cmap='RdYlBu_r',
                       vmin=T_grid.min(), vmax=T_grid.max(),
                       s=200, edgecolors='black', linewidths=2,
                       zorder=10, alpha=1.0)

            # 标注结构名
            anns, names = _annotate_clean(ax, v_clean,
                                          a['adh_col'], a['size_col'],
                                          panel_key=f'{hk}e')

            # 色标
            cbar = plt.colorbar(contourf, ax=ax, pad=0.02)
            cbar.set_label(temp_label, fontsize=28, fontfamily='Arial')
            cbar.ax.tick_params(labelsize=22)
            cbar.outline.set_linewidth(1.5)

            # 注释框: 回归方程 + R²
            eq_txt = (
                f'$R^2$ = {R2_e:.2f}\n'
                f'{a["temp_name"]} = '
                f'{reg_e.coef_[0]:+.1f}·{a["adh_col"]} '
                f'{reg_e.coef_[1]:+.1f}·{a["size_col"]} '
                f'{reg_e.intercept_:+.0f}'
            )
            ax.text(0.03, 0.97, eq_txt, transform=ax.transAxes,
                    fontsize=18, fontfamily='Arial',
                    verticalalignment='top',
                    bbox=dict(boxstyle='round,pad=0.4',
                              facecolor='white', alpha=0.85,
                              edgecolor='gray', linewidth=1))

            # 坐标轴 — Y轴标签随 size_col 变化
            size_label_map = {
                'nMetal': '$n_{Metal}$',
                'n_PtSn': '$n_{PtSn}$',
                'n_SnO':  '$n_{SnO}$',
            }
            size_axis_label = size_label_map.get(a['size_col'], a['size_col'])
            _style_clean_ax(ax, adh_label, size_axis_label,
                            x_adh_e, x_size_e,
                            user_xticks=xt_map.get((hk, 'e')),
                            user_yticks=yt_map.get((hk, 'e')))

            plt.tight_layout(pad=0.5)
            fname_e = f'{htag}_clean_e.png'

            if _is_interactive(f'{hk}e'):
                fig.savefig(f'{output_dir}/{fname_e}', dpi=300,
                            bbox_inches='tight', transparent=True)
                _drags = []
                for ann, sn in zip(anns, names):
                    da = DraggableAnnotation(ann, sn, f'{hk}e')
                    da.connect()
                    _drags.append(da)
                print(f'\n  🎯 交互模式: {hk}(e) — 拖动标签')
                plt.show()
                for ann in anns:
                    ann.get_bbox_patch().set_alpha(0)
                    ann.get_bbox_patch().set_edgecolor('none')
                fig.savefig(f'{output_dir}/{fname_e}', dpi=300,
                            bbox_inches='tight', transparent=True)
            else:
                fig.savefig(f'{output_dir}/{fname_e}', dpi=300,
                            bbox_inches='tight', transparent=True)
            n_drawn += 1

            plt.close()
            print(f'  [OK] {htag} (e) contour heatmap: {output_dir}/{fname_e}')

        # ============================================================
        # (f) Bc2: n_SnO (x) vs T_onset_O (y) — size-effect direct plot
        #     Only drawn for hypothesis B (hk == 'B')
        # ============================================================
        if hk == 'B' and _should_draw(hk, 'f'):
            fig, ax = plt.subplots(figsize=(10, 8))
            fig.patch.set_alpha(0)

            # B 的 size_col = 'n_SnO', temp_col = 'T_onset_O'
            x_f = v_clean[a['size_col']].values.astype(float)
            y_f = v_clean[a['temp_col']].values.astype(float)

            # 简单线性拟合 + 95% CI
            x_plot_f = np.linspace(x_f.min(), x_f.max(), 200)
            y_fit_f, ci_f = _confidence_band(x_f, y_f, x_plot_f)
            ax.fill_between(x_plot_f, y_fit_f - ci_f, y_fit_f + ci_f,
                            alpha=0.18, color=a['color_clean'])
            r_f, p_f = stats.pearsonr(x_f, y_f)
            if args.no_stars:
                fit_label_f = f'r = {r_f:+.2f}'
            else:
                fit_label_f = f'r = {r_f:+.2f} ({_sig_label(p_f)})'
            ax.plot(x_plot_f, y_fit_f, '-', color=a['color_clean'], lw=4,
                    label=fit_label_f)

            ax.scatter(x_f, y_f, s=200, c=a['color_clean'], edgecolors='black',
                       linewidths=2, alpha=0.85, zorder=5)

            anns, names = _annotate_clean(ax, v_clean, a['size_col'],
                                          a['temp_col'], panel_key=f'{hk}f')

            ax.legend(loc='best', fontsize=26, frameon=False,
                      prop={'family': 'Arial'})
            _style_clean_ax(ax, r'$n_{SnO}$', r'$T_2$ (K)', x_f, y_f,
                            user_xticks=xt_map.get((hk, 'f')),
                            user_yticks=yt_map.get((hk, 'f')))

            plt.tight_layout(pad=0.5)
            fname_f = f'{htag}_clean_f_size_effect.png'

            if _is_interactive(f'{hk}f'):
                fig.savefig(f'{output_dir}/{fname_f}', dpi=300,
                            bbox_inches='tight', transparent=True)
                _drags = []
                for ann, sn in zip(anns, names):
                    da = DraggableAnnotation(ann, sn, f'{hk}f')
                    da.connect()
                    _drags.append(da)
                print(f'\n  🎯 交互模式: {hk}(f) Bc2 — 拖动标签')
                plt.show()
                for ann in anns:
                    ann.get_bbox_patch().set_alpha(0)
                    ann.get_bbox_patch().set_edgecolor('none')
                fig.savefig(f'{output_dir}/{fname_f}', dpi=300,
                            bbox_inches='tight', transparent=True)
            else:
                fig.savefig(f'{output_dir}/{fname_f}', dpi=300,
                            bbox_inches='tight', transparent=True)
            n_drawn += 1

            plt.close()
            print(f'  [OK] {htag} (f) Bc2 size-effect: {output_dir}/{fname_f}')

        # ============================================================
        # (g) 单变量 R² 散点图 — 两张独立图
        #     g_adh:  adh → T   文件名: hypothesisX_clean_g_adh.png
        #     g_size: size → T  文件名: hypothesisX_clean_g_size.png
        #   R² = r²（单变量OLS），显示在图例中，无文本框
        # ============================================================
        if _should_draw(hk, 'g'):
            x_adh_g  = v_clean[a['adh_col']].values.astype(float)
            x_size_g = v_clean[a['size_col']].values.astype(float)
            y_g      = v_clean[a['temp_col']].values.astype(float)

            size_label_map_g = {
                'nMetal': '$n_{Metal}$',
                'n_PtSn': '$n_{PtSn}$',
                'n_SnO':  '$n_{SnO}$',
            }
            size_axis_label_g = size_label_map_g.get(a['size_col'], a['size_col'])

            # ---- 辅助: 画单变量散点图并保存 ----
            def _draw_univar(x_vals, y_vals, x_label, panel_suffix, xcol, ycol):
                fig_u, ax_u = plt.subplots(figsize=(10, 8))
                fig_u.patch.set_alpha(0)
                x_plot_u = np.linspace(x_vals.min(), x_vals.max(), 200)
                y_fit_u, ci_u = _confidence_band(x_vals, y_vals, x_plot_u)
                r_u, p_u = stats.pearsonr(x_vals, y_vals)
                r2_u = r_u ** 2  # 单变量 R² = r²，符号唯一
                ax_u.fill_between(x_plot_u, y_fit_u - ci_u, y_fit_u + ci_u,
                                  alpha=0.18, color=a['color_clean'])
                _parts = []
                if not args.no_r:
                    _parts.append(f'r = {r_u:+.2f}')
                if not args.no_r2:
                    _parts.append(f'$R^2$ = {r2_u:.2f}')
                fit_lbl_u = ',  '.join(_parts) if _parts else '_nolegend_'
                _fit_line, = ax_u.plot(x_plot_u, y_fit_u, '-', color=a['color_clean'],
                                       lw=4, label=fit_lbl_u)
                ax_u.scatter(x_vals, y_vals, s=200, c=a['color_clean'],
                             edgecolors='black', linewidths=2, alpha=0.85, zorder=5)
                anns_u, names_u = _annotate_clean(ax_u, v_clean, xcol, ycol,
                                                  panel_key=f'{hk}g_{panel_suffix}')
                # 散点说明图例项: 根据数据集是否含 O 自动选择化学式
                import matplotlib.lines as _mlines
                _has_O = ('nO' in v_clean.columns and v_clean['nO'].notna().any()
                          and (v_clean['nO'] > 0).any())
                _dot_label = (r'$\mathrm{Pt}_x\mathrm{Sn}_y\mathrm{O}_z$' if _has_O
                              else r'$\mathrm{Pt}_x\mathrm{Sn}_y$')
                _dot_handle = _mlines.Line2D(
                    [], [], linestyle='none',
                    marker='o', markersize=14,
                    markerfacecolor=a['color_clean'],
                    markeredgecolor='black', markeredgewidth=1.5,
                    label=_dot_label)
                _leg_handles = [_dot_handle]
                if _parts:   # 有统计量: 追加拟合线句柄
                    _leg_handles.append(_fit_line)
                ax_u.legend(handles=_leg_handles,
                            loc='best', fontsize=22, frameon=False,
                            prop={'family': 'Arial'})
                _style_clean_ax(ax_u, x_label, temp_label, x_vals, y_vals,
                                user_xticks=xt_map.get((hk, 'g')),
                                user_yticks=yt_map.get((hk, 'g')))
                # tight_layout 根据轴标签/刻度自动分配边距，与图例条数无关
                plt.tight_layout(pad=0.5)
                fname_u = f'{htag}_clean_g_{panel_suffix}.png'
                # 同时支持短格式 (Ag-adh) 和长格式 (hypothesisAg-adh / "hypothesisB'g-size")
                _hk_short = hk.replace('hypothesis', '')   # 'hypothesisA'→'A', "hypothesisB'"→"B'"
                _ikey_long  = f'{hk}g-{panel_suffix}'      # "hypothesisAg-adh"
                _ikey_short = f'{_hk_short}g-{panel_suffix}'  # "Ag-adh"
                if _is_interactive(_ikey_long) or _is_interactive(_ikey_short):
                    fig_u.savefig(f'{output_dir}/{fname_u}', dpi=300,
                                  bbox_inches='tight', transparent=True)
                    _drags = []
                    for ann, sn in zip(anns_u, names_u):
                        da = DraggableAnnotation(ann, sn, f'{hk}g_{panel_suffix}')
                        da.connect()
                        _drags.append(da)
                    print(f'\n  🎯 交互模式: {hk}(g_{panel_suffix}) — 拖动标签')
                    plt.show()
                    for ann in anns_u:
                        ann.get_bbox_patch().set_alpha(0)
                        ann.get_bbox_patch().set_edgecolor('none')
                    fig_u.savefig(f'{output_dir}/{fname_u}', dpi=300,
                                  bbox_inches='tight', transparent=True)
                else:
                    fig_u.savefig(f'{output_dir}/{fname_u}', dpi=300,
                                  bbox_inches='tight', transparent=True)
                plt.close()
                print(f'  [OK] {htag} (g_{panel_suffix}) univariate R²: {output_dir}/{fname_u}')
                return fname_u

            _draw_univar(x_adh_g,  y_g, adh_label,         'adh',  a['adh_col'],  a['temp_col'])
            _draw_univar(x_size_g, y_g, size_axis_label_g, 'size', a['size_col'], a['temp_col'])
            n_drawn += 2


    if (panel_filter is None or draw_panel_D) and len(results_partial) >= 2:
        import matplotlib.colors as mcolors

        # ---- 构建数据矩阵 ----
        # 行: 假设 A/B/C   列: 指标 (去掉前三列简单/偏相关, 只保留回归指标)
        # --heatmap-r2: 在最前面插入两列单变量 R²
        if args.heatmap_r2:
            col_keys = ['R2_adh_only', 'R2_size_only',
                        'R2_mr', 'r_partial_adh', 'r_partial_size',
                        'beta_adh', 'beta_size']
            col_labels = [
                '$R^2_{adh}$', '$R^2_{size}$',
                '$R^2$', '$|r_{adh}|$', '$|r_{size}|$',
                r'$|\beta_{adh}|$', r'$|\beta_{size}|$',
            ]
            p_keys = [None, None,
                      None, 'p_partial_adh', 'p_partial_size',
                      'p_partial_adh', 'p_partial_size']
        else:
            col_keys = ['R2_mr', 'r_partial_adh', 'r_partial_size',
                        'beta_adh', 'beta_size']
            col_labels = [
                '$R^2$', '$|r_{adh}|$', '$|r_{size}|$',
                r'$|\beta_{adh}|$', r'$|\beta_{size}|$',
            ]
            # 对应 p 值的键 (仅用于内部, 不显示星号)
            p_keys = [None, 'p_partial_adh', 'p_partial_size',
                      'p_partial_adh', 'p_partial_size']

        n_cols = len(col_keys)

        # ================================================================
        # 两套行配置: 图1=全部5行, 图2=3行(C/A/B3)
        # ================================================================
        ROW_ORDER_FULL = {'预期C': 0, '预期A': 1, "预期B'": 2, "预期B''": 3, '预期B3': 4, '预期D': 5}
        ROW_ORDER_3    = {'预期C': 0, '预期A': 1, '预期D': 2}

        ROW_LABEL_MAP = {
            '预期C':   r'$T_m$ (PtSn–AlO)',
            '预期A':   r'$T_1$ (PtSn–SnO)',
            "预期B'":  r"$T_2$ (SnO–AlO, $E^3_{adh}$)",
            "预期B''": r"$T_2$ (SnO–AlO, $E^1_{adh}$ total)",
            '预期B3':  r'$T_2$ (SnO–AlO)',
            '预期D':   r'$T_2$ (PtSnO–AlO)',
        }

        heatmap_configs = [
            # (文件名后缀, ROW_ORDER字典)
            ('partial_correlation_summary_heatmap.png',   ROW_ORDER_FULL),
            ('partial_correlation_summary_heatmap2.png',  ROW_ORDER_3),
        ]

        for fname_D, row_order in heatmap_configs:
            row_labels = []
            results_sorted = sorted(
                [r for r in results_partial if r.get('label', '') in row_order],
                key=lambda r: row_order.get(r.get('label', ''), 99))
            n_rows = len(results_sorted)
            mat = np.full((n_rows, n_cols), np.nan)
            sig_mat = [[''] * n_cols for _ in range(n_rows)]

            for i, res in enumerate(results_sorted):
                lbl = res.get('label', '')
                row_labels.append(ROW_LABEL_MAP.get(lbl, lbl))
                for j, ck in enumerate(col_keys):
                    val = res.get(ck, np.nan)
                    _is_r2_col = ck.startswith('R2_')
                    if not _is_r2_col and not np.isnan(val):
                        val = abs(val)
                    mat[i, j] = val
                    pk = p_keys[j]
                    if pk and pk in res:
                        sig_mat[i][j] = _sig_label(res[pk])

            # ---- 选择颜色映射 ----
            # 全部列统一用 sequential (YlOrRd)，范围 [0, VMAX]
            # 默认 VMAX=1.0 截断 |β|>1；可用 --heatmap-vmax 1.5 扩展
            cmap_seq = plt.cm.YlOrRd
            VMAX = args.heatmap_vmax

            # ---- 自适应图宽/高: 5列→10, 7列→13; 每行1.4英寸 ----
            fig_width  = 10 + 3 * int(args.heatmap_r2)  # 10 or 13
            fig_height = max(2.5, n_rows * 1.4)          # 3行→4.2, 5行→7.0
            val_fs     = 18 if not args.heatmap_r2 else 15  # 数值字号
            col_lbl_fs = 16 if not args.heatmap_r2 else 13  # 列标题字号

            # ---- 绘图 ----
            fig, ax = plt.subplots(figsize=(fig_width, fig_height))
            fig.patch.set_alpha(0)

            for i in range(n_rows):
                for j in range(n_cols):
                    val = mat[i, j]
                    if np.isnan(val):
                        color = 'white'
                    else:
                        color = cmap_seq(min(val / VMAX, 1.0))

                    rect = plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                         facecolor=color, edgecolor='white',
                                         linewidth=2)
                    ax.add_patch(rect)

                    # 文字颜色: 深色背景用白字 (阈值按 VMAX 比例换算)
                    text_color = 'white' if val > VMAX * 0.55 else 'black'

                    # 数值文本 (全部正值，不显示符号)
                    txt = f'{val:.2f}'

                    ax.text(j, i, txt, ha='center', va='center',
                            fontsize=val_fs, fontfamily='Arial', fontweight='bold',
                            color=text_color)

            # --heatmap-r2: 在第2列与第3列之间画分隔竖线 (单变量 | 多元)
            if args.heatmap_r2:
                ax.axvline(x=1.5, color='#555555', linewidth=2.0, linestyle='--', alpha=0.7)

            # 坐标轴设置
            ax.set_xlim(-0.5, n_cols - 0.5)
            ax.set_ylim(-0.5, n_rows - 0.5)
            ax.set_xticks(range(n_cols))
            ax.set_xticklabels(col_labels, fontsize=col_lbl_fs, fontfamily='Arial',
                               ha='center')
            ax.xaxis.set_ticks_position('top')
            ax.xaxis.set_label_position('top')
            ax.set_yticks(range(n_rows))
            ax.set_yticklabels(row_labels, fontsize=20, fontfamily='Arial')
            ax.invert_yaxis()

            # 边框
            for spine in ax.spines.values():
                spine.set_linewidth(2)
                spine.set_color('black')
            ax.tick_params(axis='both', length=0)

            # ---- 单色标 (统一 sequential) ----
            from mpl_toolkits.axes_grid1 import make_axes_locatable
            divider = make_axes_locatable(ax)

            cax_right = divider.append_axes("right", size="3%", pad=0.3)
            norm_seq = mcolors.Normalize(vmin=0, vmax=VMAX)
            sm_seq = plt.cm.ScalarMappable(cmap=cmap_seq, norm=norm_seq)
            sm_seq.set_array([])
            cbar_seq = plt.colorbar(sm_seq, cax=cax_right, orientation='vertical')
            cbar_seq.set_label('strength', fontsize=16, fontfamily='Arial')
            cbar_seq.ax.tick_params(labelsize=14)

            plt.tight_layout(pad=0.5)
            fig.savefig(f'{output_dir}/{fname_D}', dpi=300,
                        bbox_inches='tight', transparent=True)
            n_drawn += 1
            plt.close()
            print(f'  [OK] (D) summary heatmap: {output_dir}/{fname_D}')

            # ================================================================
            # (D2) 单变量 r + R² heatmap
            #      列: R²_adh | r_adh(univ) | R²_size | r_size(univ)
            #      R² 列: sequential (YlOrRd, 0→1)
            #      行: 同 (D), 由 row_order 决定
            # ================================================================
            col_keys_r  = ['R2_adh_only', 'r_adh_univ',
                           'R2_size_only', 'r_size_univ']
            col_labels_r = [r'$R^2_{adh}$',  r'$|r_{adh}|$',
                            r'$R^2_{size}$', r'$|r_{size}|$']
            p_keys_r    = [None, None, None, None]   # 不显示星号
            # 标记哪些列是 R²（sequential 色），哪些是 r（分歧色）
            is_r2_col_r = [True, False, True, False]

            n_cols_r = len(col_keys_r)
            mat_r    = np.full((n_rows, n_cols_r), np.nan)
            sig_mat_r = [[''] * n_cols_r for _ in range(n_rows)]

            for i, res in enumerate(results_sorted):
                for j, ck in enumerate(col_keys_r):
                    val = res.get(ck, np.nan)
                    mat_r[i, j] = val          # R² 本身≥0；r 保留正负号
                    pk = p_keys_r[j]
                    if pk and pk in res:
                        sig_mat_r[i][j] = _sig_label(res[pk])

            # 两套色标
            cmap_seq_r = plt.cm.YlOrRd          # R² 列: 0→1 sequential
            VLIM_R = 1.0

            # 图宽: 4列 → 9，R² 和 |r| 统一用 YlOrRd sequential
            fig_r, ax_r2 = plt.subplots(figsize=(9, fig_height))
            fig_r.patch.set_alpha(0)

            for i in range(n_rows):
                for j in range(n_cols_r):
                    val = mat_r[i, j]
                    # r 列取绝对值参与着色
                    disp_val = val if is_r2_col_r[j] else abs(val)
                    if np.isnan(disp_val):
                        color = 'white'
                    else:
                        color = cmap_seq_r(min(disp_val / VLIM_R, 1.0))

                    rect = plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                         facecolor=color, edgecolor='white', linewidth=2)
                    ax_r2.add_patch(rect)

                    text_color = 'white' if disp_val > 0.55 else 'black'
                    # R² 和 |r| 均显示绝对值，无符号，2位小数
                    txt = f'{disp_val:.2f}'
                    ax_r2.text(j, i, txt, ha='center', va='center',
                               fontsize=17, fontfamily='Arial', fontweight='bold',
                               color=text_color)

            # 虚线分隔: adh组(R²+|r|) | size组(R²+|r|)
            ax_r2.axvline(x=1.5, color='#444444', linewidth=2.0, linestyle='--', alpha=0.8)

            ax_r2.set_xlim(-0.5, n_cols_r - 0.5)
            ax_r2.set_ylim(-0.5, n_rows - 0.5)
            ax_r2.set_xticks(range(n_cols_r))
            ax_r2.set_xticklabels(col_labels_r, fontsize=15, fontfamily='Arial', ha='center')
            ax_r2.xaxis.set_ticks_position('top')
            ax_r2.xaxis.set_label_position('top')
            ax_r2.set_yticks(range(n_rows))
            ax_r2.set_yticklabels(row_labels, fontsize=20, fontfamily='Arial')
            ax_r2.invert_yaxis()
            for spine in ax_r2.spines.values():
                spine.set_linewidth(2); spine.set_color('black')
            ax_r2.tick_params(axis='both', length=0)

            # 单色标 (统一 sequential)
            from mpl_toolkits.axes_grid1 import make_axes_locatable as _mad
            div_r = _mad(ax_r2)
            cax_r = div_r.append_axes("right", size="3%", pad=0.3)
            norm_r = mcolors.Normalize(vmin=0, vmax=VLIM_R)
            sm_r   = plt.cm.ScalarMappable(cmap=cmap_seq_r, norm=norm_r)
            sm_r.set_array([])
            cbar_r = plt.colorbar(sm_r, cax=cax_r, orientation='vertical')
            cbar_r.set_label('strength', fontsize=14, fontfamily='Arial')
            cbar_r.ax.tick_params(labelsize=12)

            plt.tight_layout(pad=0.5)
            fname_D2 = fname_D.replace('.png', '_r.png')
            fig_r.savefig(f'{output_dir}/{fname_D2}', dpi=300,
                          bbox_inches='tight', transparent=True)
            n_drawn += 1
            plt.close()
            print(f'  [OK] (D2) univariate-r heatmap: {output_dir}/{fname_D2}')

    if args.interactive is not None:
        print('\n' + '=' * 60)
        print(f'  ✅ 交互完成! 偏移量已自动保存到 {args.offsets_file}')
        print(f'     下次运行时将自动加载, 无需 --clean-offsets')
        print(f'     如需手动覆盖某个点: --clean-offsets PANEL@Structure:dx,dy')
        print('=' * 60)

    print(f'\n  [OK] --plot-clean: 共生成 {n_drawn} 张独立图 (PNG)')


def main():
    """主函数"""
    
    # 构建数据框
    df = build_dataframe()
    
    # 打印数据汇总
    print_summary_table(df)
    
    # 相关性分析
    results_df = analyze_correlation(df)
    
    # 按氧含量分组分析
    # analyze_by_oxygen_content(df)  # 暂时注释，简化输出
    
    # 绘图
    plot_correlations(df)
    plot_per_atom_adhesion(df)
    
    # 偏相关分析 + 离群点检测
    results_partial = analyze_partial_correlation(df)
    
    # 尺寸拆分绘图
    plot_size_deconvolution(df, results_partial)
    
    # --plot-clean: 出版级独立图 (去离群后)
    if args.plot_clean:
        plot_clean_panels(df, results_partial)
    
    # 保存结果
    output_dir = "results/adhesion_analysis"
    import os
    os.makedirs(output_dir, exist_ok=True)
    
    df.to_csv(f"{output_dir}/adhesion_partition_data.csv", index=False)
    results_df.to_csv(f"{output_dir}/correlation_results.csv", index=False)

    # --export-temperatures: 导出 T1/T2/T2' 温度汇总表
    if args.export_temperatures:
        DISPLAY_NAME = {'Pt6Sn5O2': 'g-948-Pt6Sn5O2',
                        'Pt6Sn6O3': 'g-948-Pt6Sn6O3',
                        'Sn7Pt6O4': 'g-1051-Sn7Pt6O4'}
        t_rows = []
        for case, pd_data in PARTITION_DATA.items():
            a1 = ADHESION_TYPE1.get(case, {})
            T2_B  = pd_data.get('T_onset_O', None)
            T2_Bp = pd_data.get('T_onset_O_perO', None)
            t_rows.append({
                'case':                    DISPLAY_NAME.get(case, case),
                'nPt':                     a1.get('nPt', ''),
                'nSn':                     a1.get('nSn', ''),
                'nO':                      a1.get('nO', ''),
                'nMetal':                  a1.get('nMetal', ''),
                'T1_lindemann(K)':         pd_data.get('T1_lindemann', None),
                'T2_B(T_onset_O,K)':       T2_B,
                "T2_Bprime(T3_onset_O,K)": T2_Bp,
                'same?': ('Y' if T2_B == T2_Bp else f'diff({T2_Bp - T2_B:+d})')
                         if (T2_B is not None and T2_Bp is not None) else 'N/A',
            })
        df_temps = pd.DataFrame(t_rows)
        out_path = 'T1_T2_summary.csv'
        df_temps.to_csv(out_path, index=False, encoding='utf-8-sig')
        print(f"\n  [OK] 温度汇总表已保存: {out_path}")
        print(df_temps.to_string(index=False))
    
    # 保存偏相关结果
    partial_rows = []
    for res in results_partial:
        partial_rows.append({
            'label': res['label'], 'adhesion': res['adh'], 'temperature': res['temp'],
            'r_simple': res['r_simple'], 'p_simple': res['p_simple'],
            'r_partial': res['r_partial'], 'p_partial': res['p_partial'],
            'outliers': '; '.join(res['outliers']), 'n_outliers': res['n_outliers'],
            'r_clean': res['r_clean'], 'p_clean': res['p_clean'],
            'n_total': res['n_total'], 'n_clean': res['n_clean'],
        })
    pd.DataFrame(partial_rows).to_csv(f"{output_dir}/partial_correlation_results.csv", index=False)
    
    print(f"\n  [OK] 数据已保存到: {output_dir}/")
    
    # 结论
    print("\n")
    print("=" * 90)
    print("   物理图像: 偏相关分析揭示的粘附能—温度关系")
    print("=" * 90)

    # ---------- 从 results_partial 提取最新数据 ----------
    res_map = {r['label']: r for r in results_partial}
    _A = res_map.get('预期A', {})
    _B = res_map.get('预期B', {})
    _C = res_map.get('预期C', {})

    def _fmt(v, width=6):
        return f"{v:+.2f}" if isinstance(v, (int, float)) else str(v)

    print(f"""
  ══════════════════════════════════════════════════════════════════════════════
   I. 三条假说定义
  ══════════════════════════════════════════════════════════════════════════════
  ┌──────────────────────────────────────────────────────────────────────────┐
  │ 假说C (基准):  Eadh/atom (Pt₈Snₓ) → T_m                               │
  │   无SnO界面层的纯金属团簇, 粘附能与熔点的直接关系                        │
  │   控制变量: nMetal                                                       │
  │                                                                          │
  │ 假说A (氧化物延伸): Type2/at (÷n_PtSn) → T1_lindemann                  │
  │   含SnO的体系中, Pt-Sn团簇对SnO修饰载体的粘附 → 金属熔化温度            │
  │   控制变量: n_PtSn (Pt-Sn团簇原子数, 与Type2分母一致)                   │
  │                                                                          │
  │ 假说B (O迁移): Type3/at (÷n_SnO) → T_onset_O                           │
  │   SnO层对Al₂O₃载体的粘附 → O从SnO层向团簇迁移的起始温度                 │
  │   控制变量: n_SnO (SnO层原子数, 与Type3分母一致)                         │
  └──────────────────────────────────────────────────────────────────────────┘

  ══════════════════════════════════════════════════════════════════════════════
   II. 偏相关核心结果 (控制尺寸后)
  ══════════════════════════════════════════════════════════════════════════════
  ┌──────────────────────────────────────────────────────────────────────────┐
  │ C:  r_partial(adh) = {_fmt(_C.get('r_partial_adh',0))},  β_adh = {_fmt(_C.get('beta_adh',0))},  R² = {_C.get('R2_mr',0):.2f}  │
  │     β_size = {_fmt(_C.get('beta_size',0))} (ns)                                             │
  │     ✓ 粘附能是熔点的唯一决定因素, 尺寸几乎无独立贡献                     │
  │                                                                          │
  │ A:  r_partial(adh) = {_fmt(_A.get('r_partial_adh',0))},  β_adh = {_fmt(_A.get('beta_adh',0))},  R² = {_A.get('R2_mr',0):.2f}  │
  │     β_size = {_fmt(_A.get('beta_size',0))} (ns)                                             │
  │     ✓ 含SnO体系中粘附能仍主导金属熔化, C→A的物理图像成立                 │
  │                                                                          │
  │ B:  r_partial(adh) = {_fmt(_B.get('r_partial_adh',0))},  β_adh = {_fmt(_B.get('beta_adh',0))},  R² = {_B.get('R2_mr',0):.2f}  │
  │     β_size = {_fmt(_B.get('beta_size',0))} (***)                                            │
  │     ✗ 粘附能无独立贡献; T_onset_O 由 n_SnO 尺寸主导                      │
  └──────────────────────────────────────────────────────────────────────────┘

  ══════════════════════════════════════════════════════════════════════════════
   III. 物理图像
  ══════════════════════════════════════════════════════════════════════════════

  1. 粘附能→金属稳定性 (C+A: 已验证)
  ┌──────────────────────────────────────────────────────────────────────────┐
  │ C 建立基准: 对于无氧化物的 Pt₈Snₓ 系列:                                 │
  │   每原子粘附能越强 → 团簇与载体结合越紧密 → 熔点越高                     │
  │   偏相关控制nMetal后 r=-0.95***, 表明这是粘附能的本征效应                 │
  │                                                                          │
  │ A 延伸验证: 引入SnO界面层后, 粘附能→T1_lindemann 的关系依然成立          │
  │   r_partial=-0.85***, β_adh=-1.18 (超标准化系数: Pt-Sn团簇对               │
  │   SnO修饰载体的粘附强度是金属熔化温度的核心决定因素)                      │
  │                                                                          │
  │ → 结论: 粘附能→金属熔化 是跨体系的普适物理规律                           │
  └──────────────────────────────────────────────────────────────────────────┘

  2. O迁移温度的尺寸效应 (B: 新发现)
  ┌──────────────────────────────────────────────────────────────────────────┐
  │ 简单相关 Type3/at vs T_onset_O: r=-0.19 (ns), 本身就不显著              │
  │ 偏相关控制 n_SnO 后: r_partial=+0.20 (ns), 粘附能仍无独立贡献           │
  │                                                                          │
  │ 但 β_size(n_SnO) = -0.92***, R²=0.77:                                   │
  │   T_onset_O 几乎完全由 SnO层原子数 决定                                  │
  │                                                                          │
  │ 物理解释:                                                                │
  │   n_SnO 越大 → SnO层中可供迁移的O原子越多                                │
  │              → 存在更多低势垒迁移路径 (统计效应)                          │
  │              → 某个O原子以较低温度即可完成迁移                             │
  │              → T_onset_O 越低                                             │
  │                                                                          │
  │   每原子SnO粘附能(Type3/at)量化的是 平均 SnO-载体键强度,                │
  │   但O迁移是 局部 事件: 只需一个O原子克服势垒即可触发.                    │
  │   因此 "平均粘附强度" 不如 "可迁移O数目" 重要.                           │
  │                                                                          │
  │ → 结论: T_onset_O 的决定因素是 SnO层规模(n_SnO), 而非每原子粘附强度      │
  │         这是一个尺寸驱动的统计学效应, 而非热力学效应                      │
  └──────────────────────────────────────────────────────────────────────────┘

  3. 统一图景
  ┌──────────────────────────────────────────────────────────────────────────┐
  │                                                                          │
  │  金属熔化 (T1_lindemann):  热力学过程, 需要集体失序                      │
  │    → 粘附能(每原子)越强, 越难集体失序 → 粘附能主导 (C, A 验证)           │
  │                                                                          │
  │  O迁移起始 (T_onset_O):   局部动力学事件, 只需单个O跨越势垒              │
  │    → SnO层越大, 出现低势垒路径的概率越高 → 尺寸主导 (B 发现)             │
  │                                                                          │
  │  两种温度参数反映截然不同的物理机制:                                      │
  │    T1 = f(粘附强度)    ← 集体热力学                                      │
  │    T_onset_O = f(n_SnO) ← 局部统计力学                                   │
  │                                                                          │
  └──────────────────────────────────────────────────────────────────────────┘

  注: T_onset_O 数据来自 oxygen_migration_master CSV (AIMD轨迹分析)
      18/19 个canonical结构有T_onset_O值 (Pt6Sn5O2缺失)
    """)


if __name__ == "__main__":
    main()

