================================================================================
Filtering Statistics CSV - Data Explanation
================================================================================

DATA SOURCES:

  Energy:      ~3262 unique records (LAMMPS MD simulations)
  Lindemann:   4012 raw records from 3 CSV files
               -> 3262 unique after deduplication (by directory+structure+temp)
  Merged:      ~8,030 records (inner join on match_key)
               Note: Average ~2.5 runs per (structure, temperature) combination

COLUMN DEFINITIONS:

  original_points:  Number of data points in the ORIGINAL MERGED dataset
                    (Energy + Lindemann data, before any filtering)

  filtered_points:  Number of data points AFTER applying filtering methods

  removed_points:   Number of data points removed by filtering

  percent_removed:  Percentage of data removed (removed/original * 100)

IMPORTANT NOTES:

1. 'Original' refers to the merged Energy-Lindemann dataset, NOT raw LAMMPS data

2. Step 1 MSD filtering is based on Pt and Sn ELEMENTAL MSD:
   - msd_Pt.xvg:  Pt atoms' mean square displacement
   - msd_Sn.xvg:  Sn atoms' mean square displacement
   - NOT msd_Pt-Sn.xvg (Pt-Sn relative distance)

3. When EITHER Pt or Sn MSD is flagged as anomalous for a simulation run,
   that ENTIRE run is excluded from both Energy and Lindemann analyses.

4. This ensures consistent data quality across all analysis types.

================================================================================
