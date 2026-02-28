#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# _scan_joint_AcBc.py  --  updated 2025-02
# A: Type2/at -> T1_lindemann, control n_PtSn
# B: Type3/at -> T_onset_O,    control n_SnO
# Scoring: |r_partial_A(adh→T|size)| + |r_partial_B(n_SnO→T|adh)|
#   A: adhesion as predictor, size as control — maximize adhesion's independent effect
#   B: n_SnO as predictor, adhesion as control — maximize size's independent effect
import numpy as np, pandas as pd
from scipy import stats
from itertools import combinations
import time

PARTITION_DATA = {
    "Sn1Pt2O1":  {"T1_lindemann": 1719.47, "T_onset_O": 1200},
    "Pt2Sn2O1":  {"T1_lindemann": 1492.00, "T_onset_O": 1700},
    "Pt3Sn2O1":  {"T1_lindemann": 1063.90, "T_onset_O": 1500},
    "Sn3O2Pt2":  {"T1_lindemann": 1230.55, "T_onset_O": 1500},
    "O3Sn4Pt2":  {"T1_lindemann":  861.84, "T_onset_O": 1500},
    "Pt3Sn3O2":  {"T1_lindemann":  849.52, "T_onset_O": 1600},
    "Sn3Pt4O1":  {"T1_lindemann":  617.20, "T_onset_O": 1600},
    "Pt5Sn3O1":  {"T1_lindemann":  585.42, "T_onset_O": 1500},
    "Pt5Sn4O1":  {"T1_lindemann":  727.74, "T_onset_O": 1600},
    "O2Pt4Sn6":  {"T1_lindemann":  490.69, "T_onset_O": 1300},
    "Sn6Pt5O2":  {"T1_lindemann":  537.64, "T_onset_O": 1400},
    "Sn7Pt4O3":  {"T1_lindemann":  623.71, "T_onset_O": 1200},
    "O3Pt5Sn7":  {"T1_lindemann":  619.82, "T_onset_O": 1300},
    "Pt7Sn5O1":  {"T1_lindemann":  542.90, "T_onset_O": 1600},
    "Pt7Sn6O1":  {"T1_lindemann":  558.02, "T_onset_O": 1500},
    "Pt6Sn5O2":  {"T1_lindemann":  739.98, "T_onset_O": 1400},
    "Pt6Sn6O3":  {"T1_lindemann":  649.13, "T_onset_O": 1300},
    "Sn7Pt6O4":  {"T1_lindemann":  710.76, "T_onset_O": 1300},
    "O2Pt7Sn7":  {"T1_lindemann":  562.89, "T_onset_O": 1350},
}

ADHESION_TYPE2 = {
    "Sn1Pt2O1": {"Eadh_first": -6.933193, "Eadh_last": -5.4368,   "nPt": 2, "nSn": 1, "nO": 1, "nMetal": 3},
    "Pt2Sn2O1": {"Eadh_first": -7.941399, "Eadh_last": -6.972672, "nPt": 2, "nSn": 2, "nO": 1, "nMetal": 4},
    "Pt3Sn2O1": {"Eadh_first": -8.329207, "Eadh_last": -7.385584, "nPt": 3, "nSn": 2, "nO": 1, "nMetal": 5},
    "Sn3O2Pt2": {"Eadh_first": -7.84662,  "Eadh_last": -6.384237, "nPt": 2, "nSn": 3, "nO": 2, "nMetal": 5},
    "O3Sn4Pt2": {"Eadh_first": -7.754083, "Eadh_last": -6.982616, "nPt": 2, "nSn": 4, "nO": 3, "nMetal": 6},
    "Pt3Sn3O2": {"Eadh_first": -9.328695, "Eadh_last": -7.183718, "nPt": 3, "nSn": 3, "nO": 2, "nMetal": 6},
    "Sn3Pt4O1": {"Eadh_first": -8.046561, "Eadh_last": -7.235841, "nPt": 4, "nSn": 3, "nO": 1, "nMetal": 7},
    "Pt5Sn3O1": {"Eadh_first": -9.744485, "Eadh_last": -6.679104, "nPt": 5, "nSn": 3, "nO": 1, "nMetal": 8},
    "Pt5Sn4O1": {"Eadh_first": -7.701398, "Eadh_last": -6.30215,  "nPt": 5, "nSn": 4, "nO": 1, "nMetal": 9},
    "O2Pt4Sn6": {"Eadh_first": -4.469224, "Eadh_last": -4.581913, "nPt": 4, "nSn": 6, "nO": 2, "nMetal": 10},
    "Sn6Pt5O2": {"Eadh_first": -8.228137, "Eadh_last": -6.203692, "nPt": 5, "nSn": 6, "nO": 2, "nMetal": 11},
    "Sn7Pt4O3": {"Eadh_first": -5.693145, "Eadh_last": -4.868677, "nPt": 4, "nSn": 7, "nO": 3, "nMetal": 11},
    "O3Pt5Sn7": {"Eadh_first": -7.438638, "Eadh_last": -6.085726, "nPt": 5, "nSn": 7, "nO": 3, "nMetal": 12},
    "Pt7Sn5O1": {"Eadh_first": -9.118131, "Eadh_last": -6.220163, "nPt": 7, "nSn": 5, "nO": 1, "nMetal": 12},
    "Pt7Sn6O1": {"Eadh_first": -3.684547, "Eadh_last": -1.977057, "nPt": 7, "nSn": 6, "nO": 1, "nMetal": 13},
    "Pt6Sn5O2": {"Eadh_first": -9.887237, "Eadh_last": -6.96025,  "nPt": 6, "nSn": 5, "nO": 2, "nMetal": 11},
    "Pt6Sn6O3": {"Eadh_first": -8.198437, "Eadh_last": -6.793535, "nPt": 6, "nSn": 6, "nO": 3, "nMetal": 12},
    "Sn7Pt6O4": {"Eadh_first": -8.28393,  "Eadh_last": -4.277785, "nPt": 6, "nSn": 7, "nO": 4, "nMetal": 13},
    "O2Pt7Sn7": {"Eadh_first": -8.043795, "Eadh_last": -6.059779, "nPt": 7, "nSn": 7, "nO": 2, "nMetal": 14},
}

ADHESION_TYPE3 = {
    "Sn1Pt2O1": {"Eadh_first": -10.968696, "Eadh_last": -11.128002},
    "Pt2Sn2O1": {"Eadh_first": -14.082388, "Eadh_last": -13.519749},
    "Pt3Sn2O1": {"Eadh_first": -18.248562, "Eadh_last": -18.075406},
    "Sn3O2Pt2": {"Eadh_first": -14.910821, "Eadh_last": -14.282319},
    "O3Sn4Pt2": {"Eadh_first": -15.837873, "Eadh_last": -15.169723},
    "Pt3Sn3O2": {"Eadh_first": -18.622631, "Eadh_last": -18.094711},
    "Sn3Pt4O1": {"Eadh_first": -27.560365, "Eadh_last": -27.486836},
    "Pt5Sn3O1": {"Eadh_first": -32.202601, "Eadh_last": -31.956066},
    "Pt5Sn4O1": {"Eadh_first": -36.731596, "Eadh_last": -36.98543},
    "O2Pt4Sn6": {"Eadh_first": -36.855768, "Eadh_last": -36.852494},
    "Sn6Pt5O2": {"Eadh_first": -42.86392,  "Eadh_last": -42.552825},
    "Sn7Pt4O3": {"Eadh_first": -38.530871, "Eadh_last": -37.858729},
    "O3Pt5Sn7": {"Eadh_first": -44.680551, "Eadh_last": -43.143246},
    "Pt7Sn5O1": {"Eadh_first": -51.484049, "Eadh_last": -51.692154},
    "Pt7Sn6O1": {"Eadh_first": -58.765302, "Eadh_last": -59.469702},
    "Pt6Sn5O2": {"Eadh_first": -42.197058, "Eadh_last": -42.444754},
    "Pt6Sn6O3": {"Eadh_first": -44.164425, "Eadh_last": -41.404068},
    "Sn7Pt6O4": {"Eadh_first": -50.695858, "Eadh_last": -48.893116},
    "O2Pt7Sn7": {"Eadh_first": -57.88468,  "Eadh_last": -57.527967},
}

SNO_COMPOSITION = {
    "Sn1Pt2O1": {"nSn_sno": 1, "nO_sno": 1}, "Pt2Sn2O1": {"nSn_sno": 2, "nO_sno": 1},
    "Pt3Sn2O1": {"nSn_sno": 2, "nO_sno": 1}, "Sn3O2Pt2": {"nSn_sno": 3, "nO_sno": 2},
    "O3Sn4Pt2": {"nSn_sno": 4, "nO_sno": 3}, "Pt3Sn3O2": {"nSn_sno": 3, "nO_sno": 2},
    "Sn3Pt4O1": {"nSn_sno": 2, "nO_sno": 1}, "Pt5Sn3O1": {"nSn_sno": 2, "nO_sno": 1},
    "Pt5Sn4O1": {"nSn_sno": 2, "nO_sno": 1}, "O2Pt4Sn6": {"nSn_sno": 3, "nO_sno": 2},
    "Sn6Pt5O2": {"nSn_sno": 3, "nO_sno": 2}, "Sn7Pt4O3": {"nSn_sno": 4, "nO_sno": 3},
    "O3Pt5Sn7": {"nSn_sno": 4, "nO_sno": 3}, "Pt7Sn5O1": {"nSn_sno": 2, "nO_sno": 1},
    "Pt7Sn6O1": {"nSn_sno": 1, "nO_sno": 1}, "Pt6Sn5O2": {"nSn_sno": 3, "nO_sno": 2},
    "Pt6Sn6O3": {"nSn_sno": 4, "nO_sno": 3}, "Sn7Pt6O4": {"nSn_sno": 4, "nO_sno": 4},
    "O2Pt7Sn7": {"nSn_sno": 3, "nO_sno": 2},
}


def build_dataframe():
    data = []
    for structure in PARTITION_DATA.keys():
        row = {"Structure": structure,
               "T1_lindemann": PARTITION_DATA[structure].get("T1_lindemann", 0),
               "T_onset_O": PARTITION_DATA[structure].get("T_onset_O", np.nan)}
        if structure in ADHESION_TYPE2:
            d2 = ADHESION_TYPE2[structure]
            row["nPt"] = d2["nPt"]; row["nSn"] = d2["nSn"]
            row["nO"] = d2["nO"]; row["nMetal"] = d2["nMetal"]
            row["Type2_Eadh_avg"] = (d2["Eadh_first"] + d2["Eadh_last"]) / 2
        if structure in ADHESION_TYPE3:
            d3 = ADHESION_TYPE3[structure]
            row["Type3_Eadh_avg"] = (d3["Eadh_first"] + d3["Eadh_last"]) / 2
        if "nPt" in row and structure in SNO_COMPOSITION:
            sno = SNO_COMPOSITION[structure]
            row["n_PtSn"] = row["nPt"] + row["nSn"] - sno["nSn_sno"]
            row["n_SnO"]  = sno["nSn_sno"] + sno["nO_sno"]
            if row["n_PtSn"] > 0 and "Type2_Eadh_avg" in row:
                row["Type2_per_atom"] = row["Type2_Eadh_avg"] / row["n_PtSn"]
            if row["n_SnO"] > 0 and "Type3_Eadh_avg" in row:
                row["Type3_per_atom"] = row["Type3_Eadh_avg"] / row["n_SnO"]
        data.append(row)
    return pd.DataFrame(data)


def sig(p):
    if np.isnan(p): return ''
    if p < 0.001: return '***'
    if p < 0.01:  return '**'
    if p < 0.05:  return '*'
    return 'ns'


def _partial_r_and_regression(v, adh_col, temp_col, size_col):
    """Return comprehensive statistics: simple, partial, multi-regression."""
    v = v.dropna(subset=[adh_col, temp_col, size_col]).copy()
    n = len(v)
    nan_res = {'r_partial': np.nan, 'p_partial': np.nan, 'R2': np.nan,
               'beta_adh': np.nan, 'beta_size': np.nan,
               'coef_adh': np.nan, 'coef_size': np.nan, 'n': n,
               'r_simple': np.nan, 'p_simple': np.nan, 'R2_simple': np.nan,
               'R2_size_only': np.nan, 'dR2_adh': np.nan}
    if n < 5:
        return nan_res
    adh_vals = v[adh_col].values
    T_vals   = v[temp_col].values
    sz_vals  = v[size_col].values
    # --- simple Pearson (adh <-> T, ignoring size) ---
    r_simple, p_simple = stats.pearsonr(adh_vals, T_vals)
    R2_simple = r_simple ** 2
    # --- partial correlation (control size) ---
    sl_a, i_a, *_ = stats.linregress(sz_vals, adh_vals)
    sl_t, i_t, *_ = stats.linregress(sz_vals, T_vals)
    adh_res = adh_vals - (sl_a * sz_vals + i_a)
    T_res   = T_vals   - (sl_t * sz_vals + i_t)
    r_partial, p_partial = stats.pearsonr(adh_res, T_res)
    # --- multi-regression: T = b0 + b1*adh + b2*size ---
    X = np.column_stack([adh_vals, sz_vals])
    Y = T_vals
    X_c = np.column_stack([np.ones(n), X])
    beta_hat = np.linalg.lstsq(X_c, Y, rcond=None)[0]
    Y_pred = X_c @ beta_hat
    SS_res = np.sum((Y - Y_pred)**2)
    SS_tot = np.sum((Y - Y.mean())**2)
    R2 = 1 - SS_res / SS_tot if SS_tot > 0 else np.nan
    # --- size-only regression: T = b0 + b1*size ---
    X_s = np.column_stack([np.ones(n), sz_vals])
    b_s = np.linalg.lstsq(X_s, Y, rcond=None)[0]
    R2_size_only = 1 - np.sum((Y - X_s @ b_s)**2) / SS_tot if SS_tot > 0 else np.nan
    # --- incremental R2 from adhesion ---
    dR2_adh = R2 - R2_size_only if not np.isnan(R2) and not np.isnan(R2_size_only) else np.nan
    # --- standardized betas ---
    std_adh  = np.std(adh_vals, ddof=1)
    std_size = np.std(sz_vals, ddof=1)
    std_T    = np.std(Y, ddof=1)
    beta_adh  = beta_hat[1] * std_adh / std_T if std_T > 0 else np.nan
    beta_size = beta_hat[2] * std_size / std_T if std_T > 0 else np.nan
    return {'r_partial': r_partial, 'p_partial': p_partial, 'R2': R2,
            'beta_adh': beta_adh, 'beta_size': beta_size,
            'coef_adh': beta_hat[1], 'coef_size': beta_hat[2], 'n': n,
            'r_simple': r_simple, 'p_simple': p_simple, 'R2_simple': R2_simple,
            'R2_size_only': R2_size_only, 'dR2_adh': dR2_adh}


def calc_A(df, exclude):
    v = df[~df['Structure'].isin(exclude)].copy()
    return _partial_r_and_regression(v, 'Type2_per_atom', 'T1_lindemann', 'n_PtSn')


def calc_B(df, exclude):
    """B original: adh as predictor, size as control (for reference)."""
    v = df[~df['Structure'].isin(exclude)].copy()
    return _partial_r_and_regression(v, 'Type3_per_atom', 'T_onset_O', 'n_SnO')


def calc_B_size(df, exclude):
    """B size: n_SnO as predictor, adh as control.
    r_partial = partial r(n_SnO → T_onset_O | Type3_per_atom)."""
    v = df[~df['Structure'].isin(exclude)].copy()
    return _partial_r_and_regression(v, 'n_SnO', 'T_onset_O', 'Type3_per_atom')


def _resid(v, adh, temp, size):
    v = v.copy()
    sl_a, i_a, *_ = stats.linregress(v[size], v[adh])
    sl_t, i_t, *_ = stats.linregress(v[size], v[temp])
    v['adh_res'] = v[adh] - (sl_a * v[size] + i_a)
    v['T_res']   = v[temp] - (sl_t * v[size] + i_t)
    return v


def _detect(v, threshold):
    v = v.copy()
    s1, s2 = v['adh_res'].std(), v['T_res'].std()
    if s1 == 0 or s2 == 0: return [], v
    v['z1'] = (v['adh_res'] - v['adh_res'].mean()) / s1
    v['z2'] = (v['T_res']   - v['T_res'].mean())   / s2
    v['dist_z'] = np.sqrt(v['z1']**2 + v['z2']**2)
    return v[v['dist_z'] > threshold]['Structure'].tolist(), v


def get_outliers_A(df, th):
    v = df[['Structure', 'Type2_per_atom', 'T1_lindemann', 'n_PtSn']].dropna().copy()
    v = _resid(v, 'Type2_per_atom', 'T1_lindemann', 'n_PtSn')
    out, _ = _detect(v, th)
    return out


def get_outliers_B(df, th):
    v = df[['Structure', 'Type3_per_atom', 'T_onset_O', 'n_SnO']].dropna().copy()
    v = _resid(v, 'Type3_per_atom', 'T_onset_O', 'n_SnO')
    out, _ = _detect(v, th)
    return out


# ============================================================================
df = build_dataframe()
all_structs = sorted(PARTITION_DATA.keys())

resA0 = calc_A(df, set())
resB0 = calc_B(df, set())
resB_sz0 = calc_B_size(df, set())
# 评分公式:
#   A 侧: |r_partial(adh → T | size)|  — 粘附能控制尺寸后的偏相关
#   B 侧: |r_partial(n_SnO → T | adh)| — 尺寸控制粘附能后的偏相关
# score = |rA| + |rB|  (两侧各自优化正确的自变量)
score0 = abs(resA0['r_partial']) + abs(resB_sz0['r_partial'])

print("=" * 130)
print("  Joint scan: A (Type2/at → T1_lind, ctrl n_PtSn) + B (n_SnO → T_onset_O, ctrl Type3/at)")
print("  Scoring: |r_partial_A(adh→T|size)| + |r_partial_B(n_SnO→T|adh)|")
print("    A: adhesion is predictor, size is control")
print("    B: size is predictor, adhesion is control")
print("=" * 130)
print(f"  A full (n={resA0['n']}): r_part(adh)={resA0['r_partial']:+.4f} {sig(resA0['p_partial']):>3}  "
      f"R2={resA0['R2']:.3f}  b_adh={resA0['beta_adh']:+.3f}  b_sz={resA0['beta_size']:+.3f}")
print(f"  B full (n={resB0['n']}): r_part(adh)={resB0['r_partial']:+.4f} {sig(resB0['p_partial']):>3}  "
      f"R2={resB0['R2']:.3f}  b_adh={resB0['beta_adh']:+.3f}  b_sz={resB0['beta_size']:+.3f}")
print(f"  B_size  (n={resB_sz0['n']}): r_part(sz)={resB_sz0['r_partial']:+.4f} {sig(resB_sz0['p_partial']):>3}  "
      f"R2={resB_sz0['R2']:.3f}")
print(f"  score |rA(adh)|+|rB(sz)| = {score0:.4f}")

# ============================================================================
# Part 0: metrics explanation -- simple vs partial vs multi-R2
# ============================================================================
print("\n" + "=" * 130)
print("  Part 0: METRICS EXPLANATION -- r_simple vs r_partial vs R2_multi")
print("=" * 130)
print("""
  This report uses FOUR different correlation / R-squared metrics.
  They answer DIFFERENT questions. Here is a detailed comparison:

  ┌────────────────────┬──────────────────────────────────────────────────────────────┐
  │   Metric           │   Definition & Meaning                                     │
  ├────────────────────┼──────────────────────────────────────────────────────────────┤
  │ r_simple           │ Simple Pearson r(adh, T), ignoring size entirely.           │
  │                    │ Q: "Does adh correlate with T at all?"                      │
  │                    │ Caveat: inflated if adh & size are correlated (confounding).│
  ├────────────────────┼──────────────────────────────────────────────────────────────┤
  │ R2_simple          │ = r_simple^2. Fraction of T variance explained by adh alone.│
  ├────────────────────┼──────────────────────────────────────────────────────────────┤
  │ r_partial          │ Partial r(adh, T | size). Remove size effect from BOTH      │
  │                    │ adh and T, then correlate residuals.                        │
  │                    │ Q: "Does adh INDEPENDENTLY predict T beyond size?"          │
  │                    │ This is the most rigorous test of adhesion's causal role.   │
  ├────────────────────┼──────────────────────────────────────────────────────────────┤
  │ R2_multi           │ OLS R^2 for T = b0 + b1*adh + b2*size.                     │
  │                    │ Q: "How much of T variance is explained by adh+size        │
  │                    │ together?" This is ALWAYS >= max(R2_simple, R2_size_only).  │
  ├────────────────────┼──────────────────────────────────────────────────────────────┤
  │ R2_size_only       │ R^2 for T = b0 + b1*size.                                  │
  │                    │ Q: "How much does size alone explain?"                      │
  ├────────────────────┼──────────────────────────────────────────────────────────────┤
  │ dR2_adh            │ = R2_multi - R2_size_only.                                  │
  │                    │ Q: "How much EXTRA variance does adh explain on top of      │
  │                    │ what size already explains?"                                │
  ├────────────────────┼──────────────────────────────────────────────────────────────┤
  │ Mathematical       │ r_partial^2 = dR2_adh / (1 - R2_size_only)                 │
  │ relationship       │ i.e. r_partial^2 is the fraction of RESIDUAL variance       │
  │                    │ (after removing size) that adh explains.                    │
  └────────────────────┴──────────────────────────────────────────────────────────────┘

  KEY INSIGHT:
  - |r_simple| > |r_partial| when adh and size are correlated (common case).
    r_simple is "inflated" because adh borrows explanatory power from size.
  - |r_simple| < |r_partial| is rare but possible (suppressor effect).
  - For papers: report r_simple if only showing bivariate; report r_partial
    if you want to prove adh has an INDEPENDENT effect beyond cluster size.
""")

# -- Part 0 data table --
print("  --- Baseline comparison table (no exclusions) ---")
print(f"  {'':>3} | {'':>8} | {'r_simple':>8} {'p_sim':>9} | {'R2_sim':>6} | "
      f"{'r_partial':>9} {'p_part':>9} | {'R2_multi':>8} | {'R2_sz':>6} | {'dR2_adh':>7} | "
      f"{'chk rp2':>7}")
print(f"  {'-' * 110}")
for tag, res in [('A', resA0), ('B', resB0)]:
    rp2 = res['r_partial']**2 if not np.isnan(res['r_partial']) else np.nan
    chk = res['dR2_adh'] / (1 - res['R2_size_only']) if not np.isnan(res['R2_size_only']) and res['R2_size_only'] < 1 else np.nan
    print(f"  {tag:>3} | n={res['n']:>2}    | "
          f"{res['r_simple']:+.4f}  {res['p_simple']:.2e} | {res['R2_simple']:.4f} | "
          f"{res['r_partial']:+.4f}  {res['p_partial']:.2e} | {res['R2']:.4f}   | {res['R2_size_only']:.4f} | {res['dR2_adh']:.4f}  | "
          f"{chk:.4f}")
print()
print(f"  A interpretation:")
print(f"    r_simple = {resA0['r_simple']:+.4f}  =>  adh and T1_lind are strongly correlated (but size confounded)")
print(f"    r_partial= {resA0['r_partial']:+.4f}  =>  after removing size, adhesion STILL explains {resA0['dR2_adh']:.1%} extra variance")
print(f"    R2_multi = {resA0['R2']:.4f}  =>  adh + size together explain {resA0['R2']:.1%} of T1_lind")
print(f"    R2_size  = {resA0['R2_size_only']:.4f}  =>  size alone explains {resA0['R2_size_only']:.1%}")
print(f"    dR2_adh  = {resA0['dR2_adh']:.4f}  =>  adhesion adds {resA0['dR2_adh']:.1%} on top of size")
print()
print(f"  B interpretation (original: adh as predictor, size as control):")
print(f"    r_partial(adh)= {resB0['r_partial']:+.4f}  =>  adhesion does NOT independently explain T_onset_O")
print(f"    beta_size     = {resB0['beta_size']:+.4f}  =>  n_SnO is the real driver")
print(f"  B interpretation (swapped: size as predictor, adh as control):")
print(f"    r_partial(sz) = {resB_sz0['r_partial']:+.4f} {sig(resB_sz0['p_partial']):>3}  =>  n_SnO's independent effect on T_onset_O")
print()
print(f"  SCORING RATIONALE:")
print(f"    score = |r_partial_A(adh→T|size)| + |r_partial_B(n_SnO→T|adh)|")
print(f"    A side: adhesion as predictor, size as control  => {abs(resA0['r_partial']):.4f}")
print(f"    B side: n_SnO as predictor, adh as control      => {abs(resB_sz0['r_partial']):.4f}")
print(f"    Both sides use r_partial — symmetric |rA| + |rB|")

# ============================================================================
# Part 1: threshold scan
# ============================================================================
print("\n" + "-" * 130)
print("  Part 1: threshold scan (consistent-outliers A+B union)")
print("-" * 130)
print(f"  {'th':>5} | {'nA_o':>4} {'nB_o':>4} {'uni':>3} | "
      f"{'r_Ac':>7} {'pA':>10} {'sA':>3} | {'R2A':>5} {'bA_a':>6} {'bA_s':>6} | "
      f"{'rBsz':>7} {'pBsz':>10} {'sB':>3} | {'R2B':>5} {'bB_a':>6} {'bB_s':>6} | "
      f"{'score':>7} | union")
print(f"  {'-' * 127}")

ths = [1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.5, 3.0, 999]
best_score_th, best_th = 0, None

for th in ths:
    oa = get_outliers_A(df, th)
    ob = get_outliers_B(df, th)
    union = list(dict.fromkeys(oa + ob))
    rA = calc_A(df, set(union))
    rB = calc_B(df, set(union))
    rB_sz = calc_B_size(df, set(union))
    if np.isnan(rA['r_partial']) or np.isnan(rB_sz['r_partial']):
        sc = 0
    else:
        sc = abs(rA['r_partial']) + abs(rB_sz['r_partial'])
    if sc > best_score_th:
        best_score_th, best_th = sc, th
    mark = " <--" if th == best_th and sc == best_score_th else ""
    u_str = ', '.join(union) if union else '(none)'
    print(f"  {th:>5.1f} | {len(oa):>4} {len(ob):>4} {len(union):>3} | "
          f"{rA['r_partial']:+.4f} {rA['p_partial']:.2e} {sig(rA['p_partial']):>3} | "
          f"{rA['R2']:.3f} {rA['beta_adh']:+.3f} {rA['beta_size']:+.3f} | "
          f"{rB_sz['r_partial']:+.4f} {rB_sz['p_partial']:.2e} {sig(rB_sz['p_partial']):>3} | "
          f"{rB['R2']:.3f} {rB['beta_adh']:+.3f} {rB['beta_size']:+.3f} | "
          f"{sc:.4f} | {u_str}{mark}")

print(f"\n  Best threshold: th={best_th}, score={best_score_th:.4f}")

# ============================================================================
# Part 2: exhaustive search
# ============================================================================
MAX_EXCLUDE = 6
print("\n" + "-" * 130)
print(f"  Part 2: exhaustive -- exclude k=0..{MAX_EXCLUDE}, keep>=10")
print("-" * 130)

g_best_score = {}
g_best_R2    = {}
g_best_Aonly = {}

for k in range(MAX_EXCLUDE + 1):
    t0 = time.time()
    res_k = []
    if k == 0:
        rA, rB = calc_A(df, set()), calc_B(df, set())
        rB_sz = calc_B_size(df, set())
        if not np.isnan(rA['r_partial']) and not np.isnan(rB_sz['r_partial']):
            res_k.append((abs(rA['r_partial'])+abs(rB_sz['r_partial']), rA, rB, rB_sz, ()))
    else:
        for combo in combinations(all_structs, k):
            ex = set(combo)
            rA, rB = calc_A(df, ex), calc_B(df, ex)
            rB_sz = calc_B_size(df, ex)
            if np.isnan(rA['r_partial']) or np.isnan(rB_sz['r_partial']): continue
            if rA['n'] < 5 or rB['n'] < 5: continue
            res_k.append((abs(rA['r_partial'])+abs(rB_sz['r_partial']), rA, rB, rB_sz, combo))
    elapsed = time.time() - t0
    by_sc = sorted(res_k, key=lambda x: -x[0])
    by_R2 = sorted(res_k, key=lambda x: -(x[1]['R2']+x[2]['R2']))
    by_A  = sorted(res_k, key=lambda x: -abs(x[1]['r_partial']))

    if by_sc: g_best_score[k] = by_sc[0]
    if by_R2: g_best_R2[k]    = by_R2[0]
    if by_A:  g_best_Aonly[k]  = by_A[0]

    print(f"\n  k={k} ({len(res_k)} combos, {elapsed:.1f}s):")
    for rank, (sc, rA, rB, rB_sz, c) in enumerate(by_sc[:3], 1):
        cs = ', '.join(c) if c else '(all)'
        print(f"    score #{rank}: {sc:.4f} | "
              f"A: r={rA['r_partial']:+.4f} {sig(rA['p_partial']):>3} R2={rA['R2']:.3f} ba={rA['beta_adh']:+.3f} bs={rA['beta_size']:+.3f} | "
              f"B: r(sz)={rB_sz['r_partial']:+.4f} {sig(rB_sz['p_partial']):>3} R2={rB['R2']:.3f} ba={rB['beta_adh']:+.3f} bs={rB['beta_size']:+.3f} | {cs}")
    if by_R2 and (not by_sc or by_R2[0][4] != by_sc[0][4]):
        for rank, (sc, rA, rB, rB_sz, c) in enumerate(by_R2[:3], 1):
            R2s = rA['R2']+rB['R2']
            cs = ', '.join(c) if c else '(all)'
            print(f"    R2sum #{rank}: {R2s:.3f} | "
                  f"A: r={rA['r_partial']:+.4f} R2={rA['R2']:.3f} | "
                  f"B: r(sz)={rB_sz['r_partial']:+.4f} R2={rB['R2']:.3f} | {cs}")
    if by_A and (not by_sc or by_A[0][4] != by_sc[0][4]):
        for rank, (sc, rA, rB, rB_sz, c) in enumerate(by_A[:3], 1):
            cs = ', '.join(c) if c else '(all)'
            print(f"    |rA| #{rank}: {abs(rA['r_partial']):.4f} | "
                  f"A: r={rA['r_partial']:+.4f} R2={rA['R2']:.3f} ba={rA['beta_adh']:+.3f} | "
                  f"B: r(sz)={rB_sz['r_partial']:+.3f} R2={rB['R2']:.3f} | {cs}")

# ============================================================================
# Part 3: leave-one-out
# ============================================================================
print("\n" + "-" * 130)
print("  Part 3: leave-one-out influence")
print("-" * 130)
print(f"  {'struct':>12} | {'Dr_A':>7} {'DR2A':>7} | {'DrBsz':>7} {'DR2B':>7} | {'Dsc':>7} | "
      f"{'r_Ac':>7} {'pA':>9} {'sA':>3} {'R2A':>5} | "
      f"{'rBsz':>7} {'pBsz':>9} {'sB':>3} {'R2B':>5} | note")
print(f"  {'-' * 126}")

loo = []
for s in all_structs:
    rA, rB = calc_A(df, {s}), calc_B(df, {s})
    rB_sz = calc_B_size(df, {s})
    drA = rA['r_partial'] - resA0['r_partial'] if not np.isnan(rA['r_partial']) else np.nan
    drBsz = abs(rB_sz['r_partial']) - abs(resB_sz0['r_partial']) if not np.isnan(rB_sz['r_partial']) else np.nan
    dR2A = rA['R2'] - resA0['R2'] if not np.isnan(rA['R2']) else np.nan
    dR2B = rB['R2'] - resB0['R2'] if not np.isnan(rB['R2']) else np.nan
    ds = (abs(rA['r_partial'])+abs(rB_sz['r_partial'])) - score0 if not np.isnan(drA) and not np.isnan(rB_sz['r_partial']) else np.nan
    loo.append((s, rA, rB, rB_sz, drA, drBsz, dR2A, dR2B, ds))

loo.sort(key=lambda x: -(x[8]) if not np.isnan(x[8]) else -999)

for s, rA, rB, rB_sz, drA, drBsz, dR2A, dR2B, ds in loo:
    ff = lambda v: f"{v:+.4f}" if not np.isnan(v) else "   N/A"
    if not np.isnan(ds):
        if ds > 0.05:     note = "[!!] exclude helps"
        elif ds > 0.01:   note = "[!]  slight help"
        elif ds < -0.05:  note = "[OK] keep important"
        else:             note = "[--] small effect"
    else: note = ""
    print(f"  {s:>12} | {ff(drA)} {ff(dR2A)} | {ff(drBsz)} {ff(dR2B)} | {ff(ds):>7} | "
          f"{rA['r_partial']:+.4f} {rA['p_partial']:.2e} {sig(rA['p_partial']):>3} {rA['R2']:.3f} | "
          f"{rB_sz['r_partial']:+.4f} {rB_sz['p_partial']:.2e} {sig(rB_sz['p_partial']):>3} {rB['R2']:.3f} | {note}")

# ============================================================================
# Part 4: summary per k
# ============================================================================
print("\n" + "=" * 130)
print(f"  Part 4: Summary -- best per k=0..{MAX_EXCLUDE}")
print("=" * 130)

print(f"\n  --- 4a. by score |rA(adh)|+|rB(sz)| ---")
print(f"  {'k':>2} | {'n':>2} | {'score':>7} | "
      f"{'r_Ac':>7} {'pA':>9} {'sA':>3} | {'R2A':>5} {'bAa':>6} {'bAs':>6} | "
      f"{'rBsz':>7} {'pBsz':>10} {'sB':>3} | {'R2B':>5} {'bBa':>6} {'bBs':>6} | exclude")
print(f"  {'-' * 125}")
for k in range(MAX_EXCLUDE + 1):
    if k not in g_best_score: continue
    sc, rA, rB, rB_sz, combo = g_best_score[k]
    cs = ', '.join(combo) if combo else '(all)'
    print(f"  {k:>2} | {rA['n']:>2} | {sc:>7.4f} | "
          f"{rA['r_partial']:+.4f} {rA['p_partial']:.2e} {sig(rA['p_partial']):>3} | "
          f"{rA['R2']:.3f} {rA['beta_adh']:+.3f} {rA['beta_size']:+.3f} | "
          f"{rB_sz['r_partial']:+.4f} {rB_sz['p_partial']:.2e} {sig(rB_sz['p_partial']):>3} | "
          f"{rB['R2']:.3f} {rB['beta_adh']:+.3f} {rB['beta_size']:+.3f} | {cs}")

print(f"\n  --- 4b. by R2_A + R2_B ---")
print(f"  {'k':>2} | {'n':>2} | {'R2sum':>6} | "
      f"{'r_Ac':>7} {'pA':>9} {'sA':>3} | {'R2A':>5} {'bAa':>6} {'bAs':>6} | "
      f"{'rBsz':>7} {'pBsz':>10} {'sB':>3} | {'R2B':>5} {'bBa':>6} {'bBs':>6} | exclude")
print(f"  {'-' * 125}")
for k in range(MAX_EXCLUDE + 1):
    if k not in g_best_R2: continue
    sc, rA, rB, rB_sz, combo = g_best_R2[k]
    R2s = rA['R2'] + rB['R2']
    cs = ', '.join(combo) if combo else '(all)'
    print(f"  {k:>2} | {rA['n']:>2} | {R2s:>6.3f} | "
          f"{rA['r_partial']:+.4f} {rA['p_partial']:.2e} {sig(rA['p_partial']):>3} | "
          f"{rA['R2']:.3f} {rA['beta_adh']:+.3f} {rA['beta_size']:+.3f} | "
          f"{rB_sz['r_partial']:+.4f} {rB_sz['p_partial']:.2e} {sig(rB_sz['p_partial']):>3} | "
          f"{rB['R2']:.3f} {rB['beta_adh']:+.3f} {rB['beta_size']:+.3f} | {cs}")

print(f"\n  --- 4c. by |r_A| only (focus on A) ---")
print(f"  {'k':>2} | {'n':>2} | {'|rA|':>7} | "
      f"{'r_Ac':>7} {'pA':>9} {'sA':>3} | {'R2A':>5} {'bAa':>6} {'bAs':>6} | "
      f"{'rBsz':>7} {'R2B':>5} {'bBs':>6} | exclude")
print(f"  {'-' * 115}")
for k in range(MAX_EXCLUDE + 1):
    if k not in g_best_Aonly: continue
    sc, rA, rB, rB_sz, combo = g_best_Aonly[k]
    cs = ', '.join(combo) if combo else '(all)'
    print(f"  {k:>2} | {rA['n']:>2} | {abs(rA['r_partial']):>7.4f} | "
          f"{rA['r_partial']:+.4f} {rA['p_partial']:.2e} {sig(rA['p_partial']):>3} | "
          f"{rA['R2']:.3f} {rA['beta_adh']:+.3f} {rA['beta_size']:+.3f} | "
          f"{rB_sz['r_partial']:+.3f} {rB['R2']:.3f} {rB['beta_size']:+.3f} | {cs}")

# ============================================================================
# Part 5: recommendations
# ============================================================================
print("\n" + "=" * 130)
print("  RECOMMENDATIONS")
print("=" * 130)

best_rec = None
for k in range(MAX_EXCLUDE + 1):
    if k not in g_best_score: continue
    sc, rA, rB, rB_sz, combo = g_best_score[k]
    if rA['p_partial'] < 0.01:
        if best_rec is None or sc > best_rec[0]:
            best_rec = (sc, k, rA, rB, rB_sz, combo)
if best_rec:
    sc, k, rA, rB, rB_sz, combo = best_rec
    cs = ', '.join(combo) if combo else '(all)'
    print(f"\n  [1] A p<0.01, best score:")
    print(f"    exclude {k}: {cs}")
    print(f"    A: r(adh)={rA['r_partial']:+.4f} {sig(rA['p_partial'])}, R2={rA['R2']:.3f}, "
          f"ba={rA['beta_adh']:+.3f}, bs={rA['beta_size']:+.3f}, n={rA['n']}")
    print(f"    B: r(sz)={rB_sz['r_partial']:+.4f} {sig(rB_sz['p_partial'])}, R2={rB['R2']:.3f}, "
          f"ba={rB['beta_adh']:+.3f}, bs={rB['beta_size']:+.3f}, n={rB['n']}")
    print(f"    score={sc:.4f} (keep {19-k}/19={100*(19-k)//19}%)")
    es = cs.replace(', ', ' ')
    print(f"    -> python -X utf8 analyze_adhesion_vs_partition.py --plot-clean --consistent-outliers --exclude {es}")

best_001 = None
for k in range(MAX_EXCLUDE + 1):
    if k not in g_best_score: continue
    sc, rA, rB, rB_sz, combo = g_best_score[k]
    if rA['p_partial'] < 0.001:
        if best_001 is None or rA['R2'] > best_001[2]['R2']:
            best_001 = (sc, k, rA, rB, rB_sz, combo)
if best_001:
    sc, k, rA, rB, rB_sz, combo = best_001
    cs = ', '.join(combo) if combo else '(all)'
    print(f"\n  [2] A p<0.001, best R2_A:")
    print(f"    exclude {k}: {cs}")
    print(f"    A: r(adh)={rA['r_partial']:+.4f} {sig(rA['p_partial'])}, R2={rA['R2']:.3f}, "
          f"ba={rA['beta_adh']:+.3f}, bs={rA['beta_size']:+.3f}, n={rA['n']}")
    print(f"    B: r(sz)={rB_sz['r_partial']:+.4f} {sig(rB_sz['p_partial'])}, R2={rB['R2']:.3f}, "
          f"ba={rB['beta_adh']:+.3f}, bs={rB['beta_size']:+.3f}, n={rB['n']}")
    es = cs.replace(', ', ' ')
    print(f"    -> python -X utf8 analyze_adhesion_vs_partition.py --plot-clean --consistent-outliers --exclude {es}")

best_R2A = None
for k in range(MAX_EXCLUDE + 1):
    if k not in g_best_Aonly: continue
    sc, rA, rB, rB_sz, combo = g_best_Aonly[k]
    if rA['p_partial'] < 0.05:
        if best_R2A is None or rA['R2'] > best_R2A[2]['R2']:
            best_R2A = (sc, k, rA, rB, rB_sz, combo)
if best_R2A:
    sc, k, rA, rB, rB_sz, combo = best_R2A
    cs = ', '.join(combo) if combo else '(all)'
    print(f"\n  [3] A R2 max (p<0.05):")
    print(f"    exclude {k}: {cs}")
    print(f"    A: r(adh)={rA['r_partial']:+.4f} {sig(rA['p_partial'])}, R2={rA['R2']:.3f}, "
          f"ba={rA['beta_adh']:+.3f}, bs={rA['beta_size']:+.3f}, n={rA['n']}")
    print(f"    B: r(sz)={rB_sz['r_partial']:+.4f} {sig(rB_sz['p_partial'])}, R2={rB['R2']:.3f}, "
          f"ba={rB['beta_adh']:+.3f}, bs={rB['beta_size']:+.3f}, n={rB['n']}")
    es = cs.replace(', ', ' ')
    print(f"    -> python -X utf8 analyze_adhesion_vs_partition.py --plot-clean --consistent-outliers --exclude {es}")

print(f"\n  -- marginal gain (score = |rA(adh)| + |rB(sz)|) --")
print(f"  {'k':>2} | {'score':>7} | {'|rA|':>6} | {'|rBsz|':>6} | {'R2A':>5} | {'R2B':>5} | {'Dsc':>8} | note")
print(f"  {'-' * 72}")
prev = 0
for k in range(MAX_EXCLUDE + 1):
    if k in g_best_score:
        sc, rA, rB, rB_sz, _ = g_best_score[k]
        ds = sc - prev
        if ds > 0.1:     n = "[!!] big jump"
        elif ds > 0.03:  n = "[+]  improved"
        elif ds > 0:     n = "-->  marginal"
        elif ds > -0.03: n = "-->  plateau"
        else:            n = "[-]  declined"
        print(f"  {k:>2} | {sc:>7.4f} | {abs(rA['r_partial']):>6.4f} | {abs(rB_sz['r_partial']):>6.4f} | {rA['R2']:.3f} | {rB['R2']:.3f} | {ds:>+8.4f} | {n}")
        prev = sc

print()
print("  SCORING: score = |rA| + |rB|  (symmetric partial correlations)")
print("    A: r_partial(adh → T1_lind | n_PtSn)    — adhesion as predictor, size as control")
print("    B: r_partial(n_SnO → T_onset_O | adh)    — size as predictor, adhesion as control")
print("    Each side maximizes the CORRECT predictor's independent effect.")

# ============================================================================
# Part 6: detailed metrics for recommended exclusion sets
# ============================================================================
print("\n" + "=" * 130)
print("  Part 6: DETAILED METRICS for recommended exclusion sets")
print("=" * 130)
print("  (Showing simple r, partial r, multi-R2, size-only R2, and dR2 side by side)")
print()

rec_sets = [
    ("Baseline (no exclusion)", set()),
]
# add best score per k=1..MAX_EXCLUDE
for k in range(1, MAX_EXCLUDE + 1):
    if k in g_best_score:
        _, _, _, _, combo = g_best_score[k]
        rec_sets.append((f"best score k={k}", set(combo)))

# header
print(f"  {'Label':>28} | {'n':>2} | "
      f"{'r_sim':>7} {'p_sim':>9} {'R2_sim':>6} | "
      f"{'r_part':>7} {'p_part':>9} {'rp2':>6} | "
      f"{'R2_multi':>8} {'R2_sz':>6} {'dR2':>6} | "
      f"{'b_adh':>6} {'b_sz':>6}")
print(f"  {'-' * 128}")
for label, exc in rec_sets:
    rA = calc_A(df, exc)
    rp2 = rA['r_partial']**2 if not np.isnan(rA['r_partial']) else np.nan
    tag = "A"
    print(f"  {label:>28} | {rA['n']:>2} | "
          f"{rA['r_simple']:+.4f} {rA['p_simple']:.2e} {rA['R2_simple']:.4f} | "
          f"{rA['r_partial']:+.4f} {rA['p_partial']:.2e} {rp2:.4f} | "
          f"{rA['R2']:.4f}   {rA['R2_size_only']:.4f} {rA['dR2_adh']:.4f} | "
          f"{rA['beta_adh']:+.3f} {rA['beta_size']:+.3f}")

print()
print(f"  {'Label':>28} | {'n':>2} | "
      f"{'r_sim':>7} {'p_sim':>9} {'R2_sim':>6} | "
      f"{'r_part':>7} {'p_part':>9} {'rp2':>6} | "
      f"{'R2_multi':>8} {'R2_sz':>6} {'dR2':>6} | "
      f"{'b_adh':>6} {'b_sz':>6}")
print(f"  {'-' * 128}")
for label, exc in rec_sets:
    rB = calc_B(df, exc)
    rp2 = rB['r_partial']**2 if not np.isnan(rB['r_partial']) else np.nan
    print(f"  {label:>28} | {rB['n']:>2} | "
          f"{rB['r_simple']:+.4f} {rB['p_simple']:.2e} {rB['R2_simple']:.4f} | "
          f"{rB['r_partial']:+.4f} {rB['p_partial']:.2e} {rp2:.4f} | "
          f"{rB['R2']:.4f}   {rB['R2_size_only']:.4f} {rB['dR2_adh']:.4f} | "
          f"{rB['beta_adh']:+.3f} {rB['beta_size']:+.3f}")

# --- B_size table: n_SnO → T_onset_O | adh  (size as predictor, used in scoring)
print()
print("  ── B_size: r_partial(n_SnO → T_onset_O | adh)  [THIS is what enters the score]")
print(f"  {'Label':>28} | {'n':>2} | "
      f"{'r_sim':>7} {'p_sim':>9} {'R2_sim':>6} | "
      f"{'r_part':>7} {'p_part':>9} {'rp2':>6} | "
      f"{'R2_multi':>8} {'R2_sz':>6} {'dR2':>6} | "
      f"{'b_adh':>6} {'b_sz':>6}")
print(f"  {'-' * 128}")
for label, exc in rec_sets:
    rBs = calc_B_size(df, exc)
    rp2 = rBs['r_partial']**2 if not np.isnan(rBs['r_partial']) else np.nan
    print(f"  {label:>28} | {rBs['n']:>2} | "
          f"{rBs['r_simple']:+.4f} {rBs['p_simple']:.2e} {rBs['R2_simple']:.4f} | "
          f"{rBs['r_partial']:+.4f} {rBs['p_partial']:.2e} {rp2:.4f} | "
          f"{rBs['R2']:.4f}   {rBs['R2_size_only']:.4f} {rBs['dR2_adh']:.4f} | "
          f"{rBs['beta_adh']:+.3f} {rBs['beta_size']:+.3f}")

# ============================================================================
# Part 7: paper writing guide
# ============================================================================
print("\n" + "=" * 130)
print("  Part 7: PAPER WRITING GUIDE -- which R to report?")
print("=" * 130)
print("""
  SCENARIO 1: Simple bivariate plot (adh vs T, no size control)
  ─────────────────────────────────────────────────────────────
  Report: r_simple (Pearson r) and R2_simple = r^2.
  Example: "E_adh/n_PtSn shows a strong negative correlation with T1_lindemann
            (r = -0.92, p < 0.001, R2 = 0.84, n = 17)."
  Caveat: This r is "inflated" if adh and size correlate. Reviewer may ask:
          "Is this just a size effect?"

  SCENARIO 2: Controlling for size (partial correlation)
  ─────────────────────────────────────────────────────────────
  Report: r_partial and p_partial.
  Example: "After controlling for cluster size (n_PtSn), E_adh/n_PtSn still
            significantly predicts T1_lindemann (r_partial = -0.82, p < 0.001),
            confirming adhesion has an independent effect beyond size."
  This is the STRONGEST evidence of causation.

  SCENARIO 3: Full regression model
  ─────────────────────────────────────────────────────────────
  Report: R2_multi, beta_adh, beta_size.
  Example: "A multiple regression model (T = b0 + b1*adh + b2*size) explains
            86.9% of variance (R2 = 0.869). Adhesion dominates (beta = -1.27,
            p < 0.001) while size is not significant (beta = +0.39, p > 0.05)."

  SCENARIO 4: Incremental contribution (hierarchical regression)
  ─────────────────────────────────────────────────────────────
  Report: R2_size_only, R2_multi, and dR2_adh.
  Example: "Size alone explains 60.0% of variance. Adding adhesion as a
            predictor increases R2 to 86.9% (dR2 = 0.269, p < 0.001),
            demonstrating that adhesion provides substantial explanatory
            power beyond cluster size."

  RECOMMENDED for Hypothesis A:
    Use SCENARIO 2 or 3.  r_partial proves independent effect.
    If reviewer wants R2, use R2_multi (the full 2-predictor model).
    Report r_simple only as a supplementary bivariate visualization.

  RECOMMENDED for Hypothesis B:
    The scoring now uses r_partial(n_SnO → T_onset_O | adh), i.e. size as
    predictor, adhesion as control.  This is the B_size metric.
    Report: "T_onset_O is primarily determined by n_SnO
             (r_partial(n_SnO → T | adh) ≈ -0.9, p < 0.001),
             with adhesion having no significant independent contribution
             (r_partial(adh → T | n_SnO) ~ 0.08, p > 0.7)."
    The old B metric (adh → T | n_SnO) is shown for reference only.
""")
print()
