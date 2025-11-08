# 数据迁移文档 (Data Migration Guide)

## 📊 数据迁移概览

本文档记录了从原始分析目录到 `workflow/` 文件夹的数据迁移情况。

### 迁移状态

✅ **已完成**: 关键数据已复制到 `workflow/data/` 目录  
⏳ **待处理**: 脚本路径需要修改以使用本地数据  
📍 **位置**: `C:\Users\11207\OneDrive\02_Code\work1-PtSnO\workflow\data\`

---

## 📁 数据目录结构

```
workflow/data/
├── gmx_msd/                    # GMX MSD数据 (Step 1-5使用)
│   ├── collected_gmx_msd/      # 5,910个.xvg文件
│   └── gmx_msd_results_20251015_184626_collected/  # 3,749个.xvg文件
│
├── lammps_energy/              # LAMMPS能量数据 (Step 6, 7.4使用)
│   ├── energy_master_20251016_121110.csv
│   ├── energy_master_20251021_134929.csv
│   ├── sup/
│   │   └── energy_master_20251021_151520.csv
│   └── (其他能量CSV文件)
│
├── lindemann/                  # Lindemann指数数据 (Step 7.4使用)
│   ├── lindemann_master_run_20251025_205545.csv
│   └── lindemann_comparison_run_20251025_205545.csv
│
└── coordination/               # 配位数/Q6数据 (Step 7.5, 7.6使用)
    ├── coordination_time_series_results_sample_20251026_200908.tar.gz  (v626)
    └── coordination_time_series_results_sample_20251024_235042.tar.gz  (v624/v625)
```

---

## 🔍 数据源与使用关系

### 1️⃣ GMX MSD数据 (`data/gmx_msd/`)

**使用脚本**: Step 1-5  
**数据格式**: `.xvg` (GROMACS输出)  
**文件数量**: 9,659个 (5,910 + 3,749)  
**数据大小**: ~299 MB

#### 原始位置
```
源目录1: d:/OneDrive/py/Cv/lin/MSD_Analysis_Collection/test-unwrap-new/file/collected_gmx_msd
源目录2: d:/OneDrive/py/Cv/lin/MSD_Analysis_Collection/test-unwrap-new/file/gmx_msd_results_20251015_184626_collected
```

#### 使用的脚本
| 脚本 | 读取路径变量 | 作用 |
|------|-------------|------|
| `step1_detect_outliers.py` | `GMX_DATA_DIRS` (line 73-76) | 检测异常run，提取GMX D值 |
| `step2_ensemble_analysis.py` | `GMX_DATA_DIRS` | 多run集合平均，计算MSD |
| `step3_plot_msd.py` | 使用step2输出 | 绘制MSD曲线 |
| `step4_calculate_ensemble_D.py` | 使用step2输出 | 计算扩散系数 |
| `step5_analyze_sn_content.py` | 使用step4输出 | 分析Sn含量影响 |

#### 数据内容
- **文件命名格式**: `msd_<element>_<composition>_<temperature>_<run>.xvg`
  - 例: `msd_Pt_pt8_t300_r0.xvg`
- **数据列**: 时间(ps), MSD(nm²)
- **注释行**: 包含GMX计算的扩散系数D值
  ```
  # D[Pt] = 0.2585 (+/- 0.7552) (1e-5 cm^2/s)
  ```

---

### 2️⃣ LAMMPS能量数据 (`data/lammps_energy/`)

**使用脚本**: Step 6系列, Step 7.4  
**数据格式**: `.csv`  
**文件数量**: 5个主文件  
**数据大小**: ~2.2 MB

#### 原始位置
```
源目录: C:\Users\11207\OneDrive\02_Code\work1-PtSnO\Cv\lin\MSD_Analysis_Collection\v3_simplified_workflow\files\lammps_energy_analysis\
附加源: C:\Users\11207\OneDrive\02_Code\work1-PtSnO\Cv\lin\MSD_Analysis_Collection\v3_simplified_workflow\files\heat_capacit\
```

#### 使用的脚本
| 脚本 | 读取路径变量 | 作用 |
|------|-------------|------|
| `step6_energy_analysis_v2.py` | `ENERGY_MASTER` (line 309) | 计算热容，识别熔化温度 |
| `step6.2analyze_cv_series.py` | 使用step6输出 | 分析热容系列趋势 |
| `step6_3_adaptive_regional_heat_capacity.py` | 使用step6输出 | 自适应区域热容分析 |
| `step5.9calculate_support_heat_capacity.py` | `HEAT_CAPACITY_FILE` | 计算支撑层热容 |
| `subtract_support_v2.py` | `HEAT_CAPACITY_FILE` | 扣除支撑层热容影响 |
| `step7_4_multi_system_heat_capacity.py` | `ENERGY_MASTER` | 多体系热容对比 |

#### 关键文件
- **主数据文件**: `energy_master_20251016_121110.csv`
  - 列: 结构, 温度(K), 平均能量(eV), 能量标准差, 目录
  - ⚠️ **重要**: 能量包含团簇(~60原子) + 支撑层(240个Al₂O₃原子)
- **热容文件**: `energy_master_20251021_134929.csv`
  - 用于热容计算和熔化温度识别
- **支撑层数据**: `sup/energy_master_20251021_151520.csv`
  - 用于扣除支撑层热容贡献

#### ❓ 关键问题
```
热容包含未知的支撑层贡献:
  Cv_total = Cv_cluster + Cv_support
           = Cv_cluster + C (常数，❓未知)

支撑层热容 Cv_support:
  - 240个Al₂O₃原子
  - 估计值: ~18-21 meV/K (来自拟合)
  - ⚠️ 需要单独模拟验证
```

---

### 3️⃣ Lindemann指数数据 (`data/lindemann/`)

**使用脚本**: Step 7.4  
**数据格式**: `.csv`  
**文件数量**: 2个  
**数据大小**: ~0.97 MB

#### 原始位置
```
源目录: C:\Users\11207\OneDrive\02_Code\work1-PtSnO\Cv\lin\MSD_Analysis_Collection\v3_simplified_workflow\files\takeit\
```

#### 使用的脚本
| 脚本 | 读取路径变量 | 作用 |
|------|-------------|------|
| `step7_lindemann_analysis.py` | `LINDEMANN_FILES` (line 298) | 分析Lindemann指数，识别相变 |
| `step7_4_multi_system_heat_capacity.py` | `LINDEMANN_FILES` | 结合热容分析相变 |
| `step7_4_2_clustering_analysis.py` | 使用step7输出 | 聚类分析相态分布 |

#### 关键文件
- **主数据文件**: `lindemann_master_run_20251025_205545.csv`
  - 列: 结构, 温度(K), Lindemann指数, 目录
  - 覆盖范围: r0-r29 (30个run)
- **对比文件**: `lindemann_comparison_run_20251025_205545.csv`
  - 用于交叉验证

#### Lindemann指数物理意义
```
δ = Lindemann Index = <√(⟨r²⟩ - ⟨r⟩²)> / ⟨r⟩

判据:
  δ < 0.1  → 固态 (原子振动幅度小)
  δ ≥ 0.1  → 液态 (原子自由移动)
```

---

### 4️⃣ 配位数/Q6数据 (`data/coordination/`)

**使用脚本**: Step 7.5, 7.6  
**数据格式**: `.tar.gz` (压缩包) + 内部CSV文件  
**文件数量**: 2个压缩包 (解压后数千个CSV)  
**数据大小**: ~167 MB (压缩)

#### 原始位置
```
v626格式: C:\Users\11207\OneDrive\02_Code\work1-PtSnO\Cv\lin\MSD_Analysis_Collection\v3_simplified_workflow\files\q6_cn\v626\
v625格式: C:\Users\11207\OneDrive\02_Code\work1-PtSnO\Cv\lin\MSD_Analysis_Collection\v3_simplified_workflow\files\q6_cn\v624\
```

#### 使用的脚本
| 脚本 | 读取路径变量 | 作用 |
|------|-------------|------|
| `step7-5-unified_multi_temp_v626_analysis.py` | `base_path` (line 1042) | 多温度配位数/Q6分析 |
| `step7-6-1_temp_side_by_side_comparison.py` | 使用step7-5输出 | 不同温度并排对比 |
| `step7-6-2_individual_system_temp_comparison.py` | 使用step7-5输出 | 单体系温度对比 |
| `step7-6-3_q6_stats_comparison.py` | 使用step7-5输出 | Q6统计量对比 |
| `v625_data_locator.py` | N/A | 工具: 定位分散的run目录 |

#### 压缩包内容
- **v626格式**: `coordination_time_series_results_sample_20251026_200908.tar.gz`
  ```
  PtX-Y/
  └── composition/
      └── T<temp>.r<run>.gpu<gpu>/
          ├── coordination_time_series.csv
          ├── cluster_global_q6_time_series.csv
          ├── cluster_geometry_time_series.csv
          └── element_comparison.csv
  ```

- **v625格式**: `coordination_time_series_results_sample_20251024_235042.tar.gz`
  ```
  PtX-Y/
  └── composition/
      └── <temp>K/
          ├── coordination_time_series.csv
          ├── cluster_global_q6_time_series.csv
          └── ...
  ```

#### 格式差异 (v625 vs v626)
| 特性 | v625 | v626 |
|------|------|------|
| Run数量 | 单run | 多run (r0, r1, r2...) |
| 温度目录 | `300K/`, `400K/` | `T300.r0.gpu0/` |
| 数据量 | 较小 | 较大 (多run) |
| 使用场景 | 快速原型 | 生产分析 |

#### 主要数据列
- **coordination_time_series.csv**: 时间, 配位数(CN)
- **cluster_global_q6_time_series.csv**: 时间, Q6参数
- **element_comparison.csv**: 元素对比数据

---

## ⚠️ 待修改路径清单

### Step 1-5 (GMX MSD)
```python
# step1_detect_outliers.py (line 73-76)
旧路径:
GMX_DATA_DIRS = [
    'd:/OneDrive/py/Cv/lin/MSD_Analysis_Collection/test-unwrap-new/file/collected_gmx_msd',
    'd:/OneDrive/py/Cv/lin/MSD_Analysis_Collection/test-unwrap-new/file/gmx_msd_results_20251015_184626_collected'
]

新路径:
GMX_DATA_DIRS = [
    BASE_DIR / 'data' / 'gmx_msd' / 'collected_gmx_msd',
    BASE_DIR / 'data' / 'gmx_msd' / 'gmx_msd_results_20251015_184626_collected'
]
```

### Step 6 (LAMMPS能量)
```python
# step6_energy_analysis_v2.py (line 309)
旧路径:
ENERGY_MASTER = BASE_DIR / 'files' / 'lammps_energy_analysis' / 'energy_master_20251016_121110.csv'

新路径:
ENERGY_MASTER = BASE_DIR / 'data' / 'lammps_energy' / 'energy_master_20251016_121110.csv'
```

### Step 7.4 (Lindemann)
```python
# step7_lindemann_analysis.py (line 298)
旧路径:
DATA_DIR = BASE_DIR / 'files' / 'takeit'
LINDEMANN_FILES = sorted(DATA_DIR.glob('lindemann_master_run_*.csv'))

新路径:
DATA_DIR = BASE_DIR / 'data' / 'lindemann'
LINDEMANN_FILES = sorted(DATA_DIR.glob('lindemann_master_run_*.csv'))
```

### Step 7.5-7.6 (配位数/Q6)
```python
# step7-5-unified_multi_temp_v626_analysis.py (line 1042)
旧路径:
base_path = r"D:\OneDrive\py\Cv\lin\MSD_Analysis_Collection\v3_simplified_workflow\files\q6_cn\v626\coordination_time_series_results_sample_20251026_200908"

新路径:
base_path = BASE_DIR / 'data' / 'coordination' / 'coordination_time_series_results_sample_20251026_200908'

⚠️ 注意: 需要先解压 .tar.gz 文件
```

---

## 📦 压缩包解压说明

### 需要解压的文件
1. `coordination_time_series_results_sample_20251026_200908.tar.gz` (v626)
2. `coordination_time_series_results_sample_20251024_235042.tar.gz` (v625)

### 解压命令 (PowerShell)
```powershell
# 进入coordination目录
cd C:\Users\11207\OneDrive\02_Code\work1-PtSnO\workflow\data\coordination\

# 解压v626数据 (需要tar工具或7-Zip)
tar -xzf coordination_time_series_results_sample_20251026_200908.tar.gz

# 解压v625数据
tar -xzf coordination_time_series_results_sample_20251024_235042.tar.gz
```

### 或使用Python解压
```python
import tarfile
from pathlib import Path

coord_dir = Path('C:/Users/11207/OneDrive/02_Code/work1-PtSnO/workflow/data/coordination')

# 解压v626
with tarfile.open(coord_dir / 'coordination_time_series_results_sample_20251026_200908.tar.gz', 'r:gz') as tar:
    tar.extractall(coord_dir)

# 解压v625
with tarfile.open(coord_dir / 'coordination_time_series_results_sample_20251024_235042.tar.gz', 'r:gz') as tar:
    tar.extractall(coord_dir)
```

---

## 🔧 BASE_DIR统一设置

所有脚本应使用统一的BASE_DIR设置:

```python
from pathlib import Path

# 自动检测脚本位置
BASE_DIR = Path(__file__).parent.parent  # 指向workflow/

# 或使用绝对路径
BASE_DIR = Path(r'C:\Users\11207\OneDrive\02_Code\work1-PtSnO\workflow')
```

---

## ✅ 验证清单

迁移完成后，请验证:

- [ ] `data/gmx_msd/` 包含9,659个.xvg文件
- [ ] `data/lammps_energy/` 包含5个CSV文件
- [ ] `data/lindemann/` 包含2个CSV文件
- [ ] `data/coordination/` 解压后包含完整目录结构
- [ ] 所有脚本的BASE_DIR已修改
- [ ] 所有数据路径已更新为 `data/` 子目录
- [ ] Step 1-5运行正常 (测试step1)
- [ ] Step 6运行正常 (测试step6)
- [ ] Step 7.4运行正常 (测试step7_lindemann)
- [ ] Step 7.5/7.6运行正常 (测试step7-5)

---

## 📊 数据统计汇总

| 数据类型 | 文件数量 | 数据大小 | 用于步骤 |
|---------|---------|---------|---------|
| GMX MSD | 9,659 | 299 MB | Step 1-5 |
| LAMMPS能量 | 5 | 2.2 MB | Step 6, 7.4 |
| Lindemann | 2 | 0.97 MB | Step 7.4 |
| 配位数/Q6 | 2压缩包 | 167 MB | Step 7.5, 7.6 |
| **总计** | **9,668** | **~469 MB** | **20个脚本** |

---

## 📝 注意事项

1. **数据完整性**: 已复制的数据与原始数据目录内容一致
2. **路径兼容性**: 使用 `Path` 对象确保Windows/Linux兼容
3. **相对路径**: 所有路径基于 `BASE_DIR` 相对设置
4. **压缩包**: coordination数据需要手动解压才能使用
5. **支撑层热容**: Step 6的热容结果包含未知的支撑层贡献，需要后续校正

---

## 🔗 相关文档

- **STEP7_DATA_SOURCE_GUIDE.md**: Step 7数据源详细说明
- **README.md**: 工作流程总览和使用指南
- **SCRIPT_INDEX.md**: 脚本快速索引

---

**更新时间**: 2025-11-06  
**维护者**: GitHub Copilot  
**版本**: v1.0
