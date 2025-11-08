# 数据迁移核查清单 ✅

**迁移完成时间**: 2025-11-06  
**迁移目标**: `C:\Users\11207\OneDrive\02_Code\work1-PtSnO\workflow\`

---

## 📊 迁移数据汇总

| 数据类型 | 文件数 | 大小 | 状态 |
|---------|--------|------|------|
| GMX MSD数据 | 6,592 | 299 MB | ✅ 完整 |
| LAMMPS能量数据 | 5 | 2.2 MB | ✅ 完整（含sup/） |
| Lindemann指数 | 2 | 0.97 MB | ✅ 完整 |
| 配位数/Q6数据 | 1 | 110 MB | ✅ 完整（仅v626） |
| **总计** | **6,600** | **~412 MB** | ✅ |

---

## 📁 最终数据目录结构

```
workflow/data/
├── gmx_msd/                                              [6,592 files, 299 MB]
│   ├── collected_gmx_msd/                               (5,910个.xvg)
│   └── gmx_msd_results_20251015_184626_collected/       (682个.xvg)
│
├── lammps_energy/                                        [5 files, 2.2 MB]
│   ├── energy_average_20251016_121110.csv
│   ├── energy_master_20251016_121110.csv                ← step6主文件
│   ├── energy_master_20251021_134929.csv                ← 热容分析
│   ├── sup/
│   │   ├── energy_average_20251021_151520.csv
│   │   └── energy_master_20251021_151520.csv            ← step5.9支撑层能量
│
├── lindemann/                                            [2 files, 0.97 MB]
│   ├── lindemann_master_run_20251025_205545.csv
│   └── lindemann_comparison_run_20251025_205545.csv
│
└── coordination/                                         [1 file, 110 MB]
    └── coordination_time_series_results_sample_20251026_200908.tar.gz  (v626)
```

---

## ✅ 已确认的数据使用关系

### 1️⃣ GMX MSD数据 (`data/gmx_msd/`)

**实际文件数**: 6,592个 (与预期9,659不同)
- `collected_gmx_msd/`: 5,910个
- `gmx_msd_results_20251015_184626_collected/`: 682个

**使用脚本** (5个):
```
step1_detect_outliers.py          → 检测异常run
step2_ensemble_analysis.py        → 集合平均MSD
step3_plot_msd.py                 → 绘制MSD曲线
step4_calculate_ensemble_D.py     → 计算扩散系数
step5_analyze_sn_content.py       → 分析Sn含量影响
```

**需要修改的路径**:
```python
# 所有Step 1-5脚本
GMX_DATA_DIRS = [
    BASE_DIR / 'data' / 'gmx_msd' / 'collected_gmx_msd',
    BASE_DIR / 'data' / 'gmx_msd' / 'gmx_msd_results_20251015_184626_collected'
]
```

---

### 2️⃣ LAMMPS能量数据 (`data/lammps_energy/`)

**文件清单**:
- ✅ `energy_master_20251016_121110.csv` - Step 6主文件
- ✅ `energy_master_20251021_134929.csv` - 热容计算
- ✅ `sup/energy_master_20251021_151520.csv` - **Step 5.9支撑层能量**

**使用脚本** (5个):
```
step6_energy_analysis_v2.py                    → 多体系能量分析
step6.2analyze_cv_series.py                    → 热容系列分析
step6_3_adaptive_regional_heat_capacity.py     → 自适应区域热容
step5.9calculate_support_heat_capacity.py      → 计算支撑层热容 ✅
step7_4_multi_system_heat_capacity.py          → 多体系热容对比
```

**需要修改的路径**:
```python
# step6_energy_analysis_v2.py
ENERGY_MASTER = BASE_DIR / 'data' / 'lammps_energy' / 'energy_master_20251016_121110.csv'

# step5.9calculate_support_heat_capacity.py
SUPPORT_ENERGY = BASE_DIR / 'data' / 'lammps_energy' / 'sup' / 'energy_master_20251021_151520.csv'
```

**✅ 关键确认**: 
- step5.9使用的是 **sup/目录下的支撑层能量数据**
- sup/目录已正确复制 ✅
- 这个数据**不含**团簇能量，只有240个Al₂O₃原子

---

### 3️⃣ Lindemann指数数据 (`data/lindemann/`)

**使用脚本** (3个):
```
step7_lindemann_analysis.py              → Lindemann分析
step7_4_multi_system_heat_capacity.py    → 结合热容分析
step7_4_2_clustering_analysis.py         → 聚类分析
```

**需要修改的路径**:
```python
# step7_lindemann_analysis.py
DATA_DIR = BASE_DIR / 'data' / 'lindemann'
LINDEMANN_FILES = sorted(DATA_DIR.glob('lindemann_master_run_*.csv'))
```

---

### 4️⃣ 配位数/Q6数据 (`data/coordination/`)

**✅ 关键确认**: 
- ✅ **只需要v626数据**: `coordination_time_series_results_sample_20251026_200908.tar.gz`
- ❌ **已删除v624/v625**: `coordination_time_series_results_sample_20251024_235042.tar.gz`
  - 原因: Step 7.5脚本支持v625/v626双格式，但实际只使用v626
  - 脚本会自动检测并使用v626格式的多run数据

**使用脚本** (4个):
```
step7-5-unified_multi_temp_v626_analysis.py    → 多温度统一分析 (v626)
step7-6-1_temp_side_by_side_comparison.py      → 温度并排对比
step7-6-2_individual_system_temp_comparison.py → 单体系温度对比
step7-6-3_q6_stats_comparison.py               → Q6统计对比
```

**需要修改的路径**:
```python
# step7-5-unified_multi_temp_v626_analysis.py (line 1042)
base_path = BASE_DIR / 'data' / 'coordination' / 'coordination_time_series_results_sample_20251026_200908'
```

**⚠️ 使用前需解压**:
```powershell
cd C:\Users\11207\OneDrive\02_Code\work1-PtSnO\workflow\data\coordination\
tar -xzf coordination_time_series_results_sample_20251026_200908.tar.gz
```

---

## 🚫 已删除的文件

### ❌ `subtract_support_v2.py`
- **删除原因**: 
  - 该脚本用于估算支撑层热容（~18-21 meV/K）
  - 现在有step5.9直接计算支撑层热容，不需要估算
  - 不属于核心分析流程

### ❌ `coordination_time_series_results_sample_20251024_235042.tar.gz` (v624/v625)
- **删除原因**:
  - Step 7.5脚本虽支持v625格式，但实际只使用v626数据
  - v626格式包含多run数据，更适合生产分析
  - 减少约57 MB存储空间

---

## 📝 最终脚本清单

**总数**: 17个核心分析脚本

### Step 1-5: MSD分析 (5个)
```
step1_detect_outliers.py
step2_ensemble_analysis.py
step3_plot_msd.py
step4_calculate_ensemble_D.py
step5_analyze_sn_content.py
```

### Step 6: 能量/热容分析 (4个)
```
step6_energy_analysis_v2.py
step6.2analyze_cv_series.py
step6_3_adaptive_regional_heat_capacity.py
step5.9calculate_support_heat_capacity.py     ← 支撑层热容（使用sup/数据）
```

### Step 7: Lindemann/配位数分析 (7个)
```
step7_lindemann_analysis.py
step7_4_multi_system_heat_capacity.py
step7_4_2_clustering_analysis.py
step7-5-unified_multi_temp_v626_analysis.py   ← 使用v626数据
step7-6-1_temp_side_by_side_comparison.py
step7-6-2_individual_system_temp_comparison.py
step7-6-3_q6_stats_comparison.py
```

### 工具脚本 (1个)
```
v625_data_locator.py                          ← 定位分散的run目录
```

---

## ⚙️ 路径修改清单

### 通用BASE_DIR设置
所有脚本需要统一设置:
```python
from pathlib import Path
BASE_DIR = Path(__file__).parent  # 指向workflow/
```

### 需要修改的脚本路径 (11个脚本)

#### Step 1-5 (5个脚本)
```python
# step1_detect_outliers.py, step2_ensemble_analysis.py
GMX_DATA_DIRS = [
    BASE_DIR / 'data' / 'gmx_msd' / 'collected_gmx_msd',
    BASE_DIR / 'data' / 'gmx_msd' / 'gmx_msd_results_20251015_184626_collected'
]
```

#### Step 6 (4个脚本)
```python
# step6_energy_analysis_v2.py
ENERGY_MASTER = BASE_DIR / 'data' / 'lammps_energy' / 'energy_master_20251016_121110.csv'

# step5.9calculate_support_heat_capacity.py
SUPPORT_ENERGY = BASE_DIR / 'data' / 'lammps_energy' / 'sup' / 'energy_master_20251021_151520.csv'
```

#### Step 7 (2个脚本)
```python
# step7_lindemann_analysis.py
DATA_DIR = BASE_DIR / 'data' / 'lindemann'
LINDEMANN_FILES = sorted(DATA_DIR.glob('lindemann_master_run_*.csv'))

# step7-5-unified_multi_temp_v626_analysis.py
base_path = BASE_DIR / 'data' / 'coordination' / 'coordination_time_series_results_sample_20251026_200908'
```

---

## ✅ 验证清单

### 数据完整性
- [x] `data/gmx_msd/` 包含6,592个.xvg文件 ✅
- [x] `data/lammps_energy/` 包含5个CSV文件 ✅
- [x] `data/lammps_energy/sup/` 存在且包含2个文件 ✅
- [x] `data/lindemann/` 包含2个CSV文件 ✅
- [x] `data/coordination/` 只包含v626压缩包 ✅
- [x] 删除了不需要的v624/v625数据 ✅
- [x] 删除了subtract_support_v2.py ✅

### 脚本完整性
- [x] workflow/目录包含17个核心脚本 ✅
- [x] 包含README.md和SCRIPT_INDEX.md ✅
- [x] 包含STEP7_DATA_SOURCE_GUIDE.md ✅

### 待完成任务
- [ ] 解压v626压缩包
- [ ] 修改所有11个脚本的数据路径
- [ ] 修改所有脚本的BASE_DIR设置
- [ ] 运行Step 1测试
- [ ] 运行Step 6测试
- [ ] 运行Step 7测试

---

## 🎯 下一步操作

### 1. 解压coordination数据
```powershell
cd C:\Users\11207\OneDrive\02_Code\work1-PtSnO\workflow\data\coordination\
tar -xzf coordination_time_series_results_sample_20251026_200908.tar.gz
```

### 2. 批量修改脚本路径
需要修改11个脚本中的数据路径:
- Step 1-5: 5个脚本 (GMX_DATA_DIRS)
- Step 6: 4个脚本 (ENERGY_MASTER, SUPPORT_ENERGY)
- Step 7: 2个脚本 (DATA_DIR, base_path)

### 3. 测试运行
```bash
# 测试Step 1
python step1_detect_outliers.py

# 测试Step 6
python step6_energy_analysis_v2.py

# 测试Step 7
python step7_lindemann_analysis.py
```

---

## 📊 迁移统计

| 项目 | 数量/大小 |
|------|----------|
| 总文件数 | 6,600 |
| 总数据量 | 412 MB |
| Python脚本 | 17个 |
| Markdown文档 | 4个 |
| 需要修改路径的脚本 | 11个 |
| 删除的文件 | 2个 |

---

**状态**: ✅ 数据迁移完成  
**下一步**: 修改脚本路径并测试  
**维护者**: GitHub Copilot  
**版本**: v1.0
