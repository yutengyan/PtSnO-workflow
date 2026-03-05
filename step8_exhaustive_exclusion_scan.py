"""
step8_exhaustive_exclusion_scan.py
====================================
穷举排除扫描: 同时优化两个假说的 R²
  目标 A: E2/n_PtSn (Type2_per_atom)  vs  T1_lindemann
  目标 B: E4/nO     (整体吸附能/O数)   vs  T3_perO (T_onset_O_perO)

用法:
    python step8_exhaustive_exclusion_scan.py            # 扫描 k=1~4, 每组显示 Top8
    python step8_exhaustive_exclusion_scan.py --topn 12  # 每组显示 Top12
    python step8_exhaustive_exclusion_scan.py --kmax 3   # 只扫到 k=3
"""
import argparse
import numpy as np
from itertools import combinations
from scipy import stats

# ══════════════════════════════════════════════════════════════════
#  数据
# ══════════════════════════════════════════════════════════════════
# T1_lindemann: Lindemann指数熔化温度 (K)
# T3_perO:      每O原子归一化O迁移起始温度 T_onset_O_perO (K)
TEMP = {
    'Sn1Pt2O1':  {'T1_lind': 1719.47, 'T3_perO': 1200},
    'Pt2Sn2O1':  {'T1_lind': 1492.00, 'T3_perO': 1700},
    'Pt3Sn2O1':  {'T1_lind': 1063.90, 'T3_perO': 1600},
    'Sn3O2Pt2':  {'T1_lind': 1230.55, 'T3_perO': 1500},
    'O3Sn4Pt2':  {'T1_lind':  861.84, 'T3_perO': 1300},
    'Pt3Sn3O2':  {'T1_lind':  849.52, 'T3_perO': 1600},
    'Sn3Pt4O1':  {'T1_lind':  617.20, 'T3_perO': 1700},
    'Pt5Sn3O1':  {'T1_lind':  585.42, 'T3_perO': 1600},
    'Pt5Sn4O1':  {'T1_lind':  727.74, 'T3_perO': 1600},
    'O2Pt4Sn6':  {'T1_lind':  490.69, 'T3_perO': 1300},
    'Sn6Pt5O2':  {'T1_lind':  537.64, 'T3_perO': 1300},
    'Sn7Pt4O3':  {'T1_lind':  623.71, 'T3_perO': 1100},
    'Pt6Sn5O2':  {'T1_lind':  739.98, 'T3_perO': 1300},
    'O3Pt5Sn7':  {'T1_lind':  619.82, 'T3_perO': 1100},
    'Pt7Sn5O1':  {'T1_lind':  542.90, 'T3_perO': 1600},
    'Pt6Sn6O3':  {'T1_lind':  649.13, 'T3_perO': 1200},
    'Pt7Sn6O1':  {'T1_lind':  558.02, 'T3_perO': 1600},
    'Sn7Pt6O4':  {'T1_lind':  710.76, 'T3_perO': 1200},
    'O2Pt7Sn7':  {'T1_lind':  562.89, 'T3_perO': 1300},
}

# E2 = Type2_per_atom = ADHESION_TYPE2.Eadh_last / n_PtSn
# n_PtSn = nPt + nSn - nSn_sno  (来自 SNO_COMPOSITION)
E2_RAW = {
    #  case            Eadh_last(Type2)  nO  n_PtSn
    'Sn1Pt2O1':  {'E': -5.436800,  'nO': 1, 'nPtSn': 2},
    'Pt2Sn2O1':  {'E': -6.972672,  'nO': 1, 'nPtSn': 2},
    'Pt3Sn2O1':  {'E': -7.385584,  'nO': 1, 'nPtSn': 3},
    'Sn3O2Pt2':  {'E': -6.384237,  'nO': 2, 'nPtSn': 2},
    'O3Sn4Pt2':  {'E': -6.982616,  'nO': 3, 'nPtSn': 2},
    'Pt3Sn3O2':  {'E': -7.183718,  'nO': 2, 'nPtSn': 3},
    'Sn3Pt4O1':  {'E': -7.235841,  'nO': 1, 'nPtSn': 5},
    'Pt5Sn3O1':  {'E': -6.679104,  'nO': 1, 'nPtSn': 6},
    'Pt5Sn4O1':  {'E': -6.302150,  'nO': 1, 'nPtSn': 7},
    'O2Pt4Sn6':  {'E': -4.581913,  'nO': 2, 'nPtSn': 7},
    'Sn6Pt5O2':  {'E': -6.203692,  'nO': 2, 'nPtSn': 8},
    'Sn7Pt4O3':  {'E': -4.868677,  'nO': 3, 'nPtSn': 7},
    'Pt6Sn5O2':  {'E': -6.960250,  'nO': 2, 'nPtSn': 8},
    'O3Pt5Sn7':  {'E': -6.085726,  'nO': 3, 'nPtSn': 8},
    'Pt7Sn5O1':  {'E': -6.220163,  'nO': 1, 'nPtSn': 10},
    'Pt6Sn6O3':  {'E': -6.793535,  'nO': 3, 'nPtSn': 8},
    'Pt7Sn6O1':  {'E': -1.977057,  'nO': 1, 'nPtSn': 12},
    'Sn7Pt6O4':  {'E': -4.277785,  'nO': 4, 'nPtSn': 9},
    'O2Pt7Sn7':  {'E': -6.059779,  'nO': 2, 'nPtSn': 11},
}

# E4 = 整个 PtSnO 簇 / Al2O3 的吸附能 (用户新提供)
E4_RAW = {
    'Sn1Pt2O1': {'E': -3.268905, 'nO': 1},
    'Pt2Sn2O1': {'E': -2.684449, 'nO': 1},
    'Pt3Sn2O1': {'E': -3.320033, 'nO': 1},
    'Sn3O2Pt2': {'E': -3.174117, 'nO': 2},
    'O3Sn4Pt2': {'E': -4.039789, 'nO': 3},
    'Pt3Sn3O2': {'E': -3.682498, 'nO': 2},
    'Sn3Pt4O1': {'E': -2.400166, 'nO': 1},
    'Pt5Sn3O1': {'E': -1.968795, 'nO': 1},
    'Pt5Sn4O1': {'E': -2.966004, 'nO': 1},
    'O2Pt4Sn6': {'E': -2.643805, 'nO': 2},
    'Sn6Pt5O2': {'E': -2.964427, 'nO': 2},
    'Sn7Pt4O3': {'E': -2.233587, 'nO': 3},
    'Pt6Sn5O2': {'E': -3.426205, 'nO': 2},
    'O3Pt5Sn7': {'E': -3.692219, 'nO': 3},
    'Pt7Sn5O1': {'E': -2.740339, 'nO': 1},
    'Pt6Sn6O3': {'E': -3.214217, 'nO': 3},
    'Pt7Sn6O1': {'E': -2.442727, 'nO': 1},
    'Sn7Pt6O4': {'E': -6.249845, 'nO': 4},
    'O2Pt7Sn7': {'E': -3.965053, 'nO': 2},
}

# ══════════════════════════════════════════════════════════════════
#  构建数组
# ══════════════════════════════════════════════════════════════════
ALL_CASES = sorted(set(TEMP) & set(E2_RAW) & set(E4_RAW))


def build_arrays(exclude=None):
    """返回 (cases, xA, yA, xB, yB)"""
    if exclude is None:
        exclude = set()
    cases, xA, yA, xB, yB = [], [], [], [], []
    for c in ALL_CASES:
        if c in exclude:
            continue
        e2 = E2_RAW[c]; e4 = E4_RAW[c]; t = TEMP[c]
        cases.append(c)
        xA.append(e2['E'] / e2['nPtSn'])   # Type2_per_atom
        yA.append(t['T1_lind'])             # T1_lindemann
        xB.append(e4['E'] / e4['nO'])       # E4/nO
        yB.append(t['T3_perO'])             # T3_perO
    return (np.array(cases),
            np.array(xA), np.array(yA),
            np.array(xB), np.array(yB))


def r2_pair(excl):
    """计算排除 excl 后的 (R²_A, R²_B)，若样本不足返回 (None, None)"""
    _, xA, yA, xB, yB = build_arrays(exclude=excl)
    if len(xA) < 4:
        return None, None
    rA, _ = stats.pearsonr(xA, yA)
    rB, _ = stats.pearsonr(xB, yB)
    return rA**2, rB**2


def sig(p):
    return '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'


# ══════════════════════════════════════════════════════════════════
#  主程序
# ══════════════════════════════════════════════════════════════════
def main():
    parser = argparse.ArgumentParser(description='穷举排除扫描: 同时优化 R²_A 和 R²_B')
    parser.add_argument('--topn', type=int, default=8,  help='每个 k 显示前 N 名 (默认 8)')
    parser.add_argument('--kmax', type=int, default=4,  help='最大排除点数 (默认 4)')
    parser.add_argument('--sort', choices=['joint','A','B'], default='joint',
                        help='排序依据: joint=R²_A+R²_B, A=仅R²_A, B=仅R²_B (默认 joint)')
    args = parser.parse_args()

    # 基线
    _, xA0, yA0, xB0, yB0 = build_arrays()
    r2A0 = stats.pearsonr(xA0, yA0)[0] ** 2
    r2B0 = stats.pearsonr(xB0, yB0)[0] ** 2
    rA0, pA0 = stats.pearsonr(xA0, yA0)
    rB0, pB0 = stats.pearsonr(xB0, yB0)

    print()
    print('=' * 84)
    print('  穷举排除扫描  |  目标A: E2/n_PtSn vs T1_lind   目标B: E4/nO vs T3_perO')
    print('=' * 84)
    print(f'  基线 (n={len(xA0)}):')
    print(f'    A: r={rA0:+.3f}{sig(pA0):3}  R²={r2A0:.3f}')
    print(f'    B: r={rB0:+.3f}{sig(pB0):3}  R²={r2B0:.3f}')
    print(f'    综合 = {r2A0 + r2B0:.3f}')
    print('=' * 84)

    sort_key = {
        'joint': lambda x: -(x[0] + x[1]),
        'A':     lambda x: -x[0],
        'B':     lambda x: -x[1],
    }[args.sort]

    ALL = list(ALL_CASES)
    for k in range(1, args.kmax + 1):
        results = []
        for combo in combinations(ALL, k):
            excl = set(combo)
            rA2, rB2 = r2_pair(excl)
            if rA2 is None:
                continue
            results.append((rA2, rB2, excl))
        results.sort(key=sort_key)
        total = len(results)
        show  = min(args.topn, total)

        print(f'\n── k={k}  排除 {k} 个点  共 {total} 组  → Top {show}  (排序: {args.sort}) ──')
        hdr = f"  {'排除集合':<44} {'R²_A':>6} {'R²_B':>6} {'综合':>6}  {'ΔR²_A':>7} {'ΔR²_B':>7}"
        print(hdr)
        print('  ' + '-' * (len(hdr) - 2))
        for rA2, rB2, excl in results[:show]:
            excl_str = '{' + ', '.join(sorted(excl)) + '}'
            dA = rA2 - r2A0
            dB = rB2 - r2B0
            print(f"  {excl_str:<44} {rA2:6.3f} {rB2:6.3f} {rA2+rB2:6.3f}"
                  f"  {dA:+7.3f} {dB:+7.3f}")

    print()
    print('=' * 84)
    print('  说明:')
    print('    R²_A  = E2/n_PtSn (Type2_per_atom) vs T1_lindemann 的决定系数')
    print('    R²_B  = E4/nO                       vs T3_perO      的决定系数')
    print('    综合  = R²_A + R²_B  (等权重联合优化指标)')
    print('    ΔR²   = 相对基线的提升量')
    print('=' * 84)


if __name__ == '__main__':
    main()
