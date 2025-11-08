# Data Filtering Statistics Report

**Generated**: 2025-11-08 11:50:04

================================================================================

## Data Source Explanation

**Important**: The 'Original' data refers to the **merged energy-Lindemann dataset** before applying the filtering methods listed below. This merged dataset combines:

1. **Energy data**: Average cluster energy from LAMMPS MD simulations
   - Raw source: ~3262 records (energy_master_*.csv)
2. **Lindemann index data**: Calculated from MSD (Mean Square Displacement) analysis
   - Raw source: 4012 records from 3 lindemann_master_run_*.csv files
   - After deduplication: ~3262 unique (directory, structure, temperature) combinations
   - Note: Multiple files contain overlapping data; duplicates are removed before merging

**Note about MSD types**: Step 1 filtering is based on **Pt and Sn elemental MSD** (`msd_Pt.xvg` and `msd_Sn.xvg`), NOT the Pt-Sn distance MSD. When either Pt or Sn MSD shows anomalous behavior for a simulation run, that run is excluded from **both** energy and Lindemann analyses.

================================================================================

## 1. Filtering Methods Applied

1. Method 1 (Step 1 MSD outliers)

## 2. Overall Statistics

- **Original merged data points**: 3,262
- **Filtered data points**: 2,692
- **Removed data points**: 570 (17.47%)
- **Retention rate**: 82.53%

## 3. Statistics by System Type

| System Type | Original | Filtered | Removed | % Removed |
|-------------|----------|----------|---------|----------|
| Cv | 95 | 95 | 0 | 0.00% |
| Pt6 | 30 | 30 | 0 | 0.00% |
| Pt6SnX | 330 | 327 | 3 | 0.91% |
| Pt8SnX | 1,067 | 880 | 187 | 17.53% |
| PtxSny | 240 | 222 | 18 | 7.50% |
| PtxSnyOz | 1,500 | 1,138 | 362 | 24.13% |

## 4. Statistics by Structure (Summary)

| Structure | Original | Filtered | Removed | % Removed |
|-----------|----------|----------|---------|----------|
| pt8sn1-2-best | 97 | 65 | 32 | 32.99% |
| pt8sn2-1-best | 97 | 72 | 25 | 25.77% |
| pt8sn10-2-best | 97 | 73 | 24 | 24.74% |
| g-1-O1Sn4Pt3 | 60 | 38 | 22 | 36.67% |
| Sn3Pt4O1 | 60 | 39 | 21 | 35.00% |
| Pt7Sn6O1 | 60 | 40 | 20 | 33.33% |
| Sn4Pt3O1 | 60 | 41 | 19 | 31.67% |
| O2Pt4Sn6 | 60 | 41 | 19 | 31.67% |
| Pt3Sn2O1 | 60 | 42 | 18 | 30.00% |
| pt8sn3-1-best | 97 | 80 | 17 | 17.53% |
| O2Sn8Pt7 | 60 | 43 | 17 | 28.33% |
| pt8sn4-1-best | 97 | 81 | 16 | 16.49% |
| Sn1Pt2O1 | 60 | 44 | 16 | 26.67% |
| Sn7Pt4O3 | 60 | 44 | 16 | 26.67% |
| Pt3Sn3O2 | 60 | 44 | 16 | 26.67% |
| Pt7Sn5O1 | 60 | 45 | 15 | 25.00% |
| Sn6Pt5O2 | 60 | 45 | 15 | 25.00% |
| pt8sn5-1-best | 97 | 82 | 15 | 15.46% |
| pt8sn8-1-best | 97 | 82 | 15 | 15.46% |
| Sn5O2Pt4 | 60 | 46 | 14 | 23.33% |

*Showing top 20 structures by removed points. See CSV for complete data.*

## 5. Statistics by Temperature (All Structures)

| Temperature (K) | Original | Filtered | Removed | % Removed |
|-----------------|----------|----------|---------|----------|
| 200 | 292 | 269 | 23 | 7.88% |
| 250 | 5 | 5 | 0 | 0.00% |
| 300 | 325 | 294 | 31 | 9.54% |
| 350 | 5 | 5 | 0 | 0.00% |
| 400 | 325 | 291 | 34 | 10.46% |
| 450 | 5 | 5 | 0 | 0.00% |
| 500 | 325 | 277 | 48 | 14.77% |
| 550 | 5 | 5 | 0 | 0.00% |
| 600 | 325 | 279 | 46 | 14.15% |
| 650 | 5 | 5 | 0 | 0.00% |
| 700 | 325 | 280 | 45 | 13.85% |
| 750 | 5 | 5 | 0 | 0.00% |
| 800 | 325 | 275 | 50 | 15.38% |
| 850 | 5 | 5 | 0 | 0.00% |
| 900 | 325 | 250 | 75 | 23.08% |
| 950 | 5 | 5 | 0 | 0.00% |
| 1000 | 325 | 235 | 90 | 27.69% |
| 1050 | 5 | 5 | 0 | 0.00% |
| 1100 | 325 | 197 | 128 | 39.38% |

## 6. Detailed Breakdown

For complete structure-by-structure, temperature-by-temperature breakdown, see: `filtering_statistics_by_structure_temp.csv`

