# 项目重构指南 (Project Refactoring Guide)

📅 **创建日期**: 2025年11月8日  
🎯 **目标**: 改善代码组织、统一命名规范、提升可维护性

---

## 📋 目录

1. [当前问题分析](#当前问题分析)
2. [建议的目标结构](#建议的目标结构)
3. [命名规范](#命名规范)
4. [重构步骤计划](#重构步骤计划)
5. [代码重构清单](#代码重构清单)
6. [实施优先级](#实施优先级)
7. [风险评估与备份](#风险评估与备份)

---

## 🔍 当前问题分析

### 问题1: 命名不统一
- ❌ **混用下划线和连字符**: `step7_4` vs `step7-5` vs `step7-6-1`
- ❌ **版本号混乱**: `v2`, `v625`, `v626` 散落在不同脚本中
- ❌ **临时命名未清理**: `step5.9`, `step6.2`, `step6_3` 等过渡版本

### 问题2: 目录结构扁平
- ❌ **21个脚本全在根目录**: 难以快速定位相关功能
- ❌ **缺少代码复用**: 每个脚本独立实现相似功能
- ❌ **配置分散**: 路径、参数硬编码在各个脚本中

### 问题3: 文档管理
- ✅ **已有文档**: README.md, WORKFLOW_OVERVIEW.md, SCRIPT_INDEX.md
- ⚠️ **文档同步**: 代码更新后文档需要手动同步
- ⚠️ **文档冗余**: 三个文档有部分重叠内容

---

## 🏗️ 建议的目标结构

```
workflow/
│
├── README.md                          # 项目总览(保留)
├── WORKFLOW_GUIDE.md                  # 工作流程指南(合并+简化)
├── CHANGELOG.md                       # 变更日志(新增)
│
├── config/                            # 📁 配置文件夹(新增)
│   ├── __init__.py
│   ├── paths.py                       # 统一路径配置
│   ├── analysis_params.py             # 分析参数配置
│   └── plotting_styles.py             # 绘图风格配置
│
├── scripts/                           # 📁 分析脚本(重组)
│   ├── __init__.py
│   │
│   ├── step01_diffusion/              # 步骤1: 扩散分析
│   │   ├── __init__.py
│   │   ├── detect_outliers.py         # 原 step1_detect_outliers.py
│   │   └── README.md
│   │
│   ├── step02_ensemble/               # 步骤2: 系综分析
│   │   ├── __init__.py
│   │   ├── ensemble_analysis.py       # 原 step2_ensemble_analysis.py
│   │   └── README.md
│   │
│   ├── step03_msd/                    # 步骤3: MSD分析
│   │   ├── __init__.py
│   │   ├── plot_msd.py                # 原 step3_plot_msd.py
│   │   └── README.md
│   │
│   ├── step04_diffusion_coeff/        # 步骤4: 扩散系数
│   │   ├── __init__.py
│   │   ├── calculate_ensemble_D.py    # 原 step4_calculate_ensemble_D.py
│   │   └── README.md
│   │
│   ├── step05_composition/            # 步骤5: 成分分析
│   │   ├── __init__.py
│   │   ├── analyze_sn_content.py      # 原 step5_analyze_sn_content.py
│   │   └── README.md
│   │
│   ├── step06_energy/                 # 步骤6: 能量分析
│   │   ├── __init__.py
│   │   ├── energy_analysis.py         # 原 step6_energy_analysis_v2.py
│   │   ├── analyze_cv_series.py       # 原 step6.2analyze_cv_series.py
│   │   ├── adaptive_regional_cv.py    # 原 step6_3_adaptive_regional_heat_capacity.py
│   │   ├── calculate_support_cv.py    # 原 step5.9calculate_support_heat_capacity.py
│   │   └── README.md
│   │
│   └── step07_structure/              # 步骤7: 结构分析
│       ├── __init__.py
│       ├── lindemann_analysis.py      # 原 step7_lindemann_analysis.py
│       ├── clustering_analysis.py     # 原 step7_4_2_clustering_analysis.py
│       ├── multi_system_cv.py         # 原 step7_4_multi_system_heat_capacity.py
│       ├── unified_multi_temp.py      # 原 step7-5-unified_multi_temp_v626_analysis.py
│       ├── temp_comparison_parallel.py    # 原 step7-6-1_temp_side_by_side_comparison.py
│       ├── temp_comparison_single.py      # 原 step7-6-2_individual_system_temp_comparison.py
│       ├── q6_stats_comparison.py         # 原 step7-6-3_q6_stats_comparison.py
│       └── README.md
│
├── utils/                             # 📁 工具函数(新增)
│   ├── __init__.py
│   ├── data_loader.py                 # 数据读取通用函数
│   ├── path_helper.py                 # 路径处理通用函数
│   ├── plotting.py                    # 绘图通用函数
│   ├── statistics.py                  # 统计分析通用函数
│   └── v625_data_locator.py           # 原有定位工具(保留)
│
├── tools/                             # 📁 辅助工具(新增)
│   ├── __init__.py
│   ├── run_script.py                  # 统一脚本启动器
│   ├── check_dependencies.py          # 依赖检查
│   ├── update_data_paths.py           # 原 update_data_paths.py
│   ├── verify_data_paths.py           # 原 verify_data_paths.py
│   └── check_script_outputs.py        # 原 check_script_outputs.py
│
├── data/                              # 📁 数据文件夹(保留结构)
│   ├── coordination/
│   ├── gmx_msd/
│   ├── lammps_energy/
│   └── lindemann/
│
├── results/                           # 📁 结果文件夹(保留结构)
│   ├── step01_outliers/
│   ├── step02_ensemble/
│   ├── step03_msd_curves/
│   ├── step04_diffusion_coeff/
│   ├── step05_composition/
│   ├── step06_energy/
│   └── step07_structure/
│
└── tests/                             # 📁 单元测试(新增)
    ├── __init__.py
    ├── test_data_loader.py
    ├── test_path_helper.py
    └── test_statistics.py
```

---

## 📝 命名规范

### 1. 脚本命名规则

**格式**: `<动词>_<对象>.py`

| 原文件名 | 新文件名 | 命名逻辑 |
|---------|---------|---------|
| `step1_detect_outliers.py` | `detect_outliers.py` | 动词+对象 |
| `step2_ensemble_analysis.py` | `ensemble_analysis.py` | 对象+动作 |
| `step6_energy_analysis_v2.py` | `energy_analysis.py` | 去除版本号 |
| `step7-5-unified_multi_temp_v626_analysis.py` | `unified_multi_temp.py` | 简化名称 |
| `step7-6-1_temp_side_by_side_comparison.py` | `temp_comparison_parallel.py` | 描述性命名 |

**核心原则**:
- ✅ 使用下划线 `_` (不用连字符 `-`)
- ✅ 去除步骤编号(用文件夹区分)
- ✅ 去除版本号(用Git管理)
- ✅ 使用完整单词(避免缩写如 `sn`, `cv`)
- ✅ 长度控制在3-4个词以内

### 2. 文件夹命名规则

**格式**: `step<NN>_<功能描述>`

```
step01_diffusion        # 两位数编号
step02_ensemble         # 全小写
step06_energy           # 下划线分隔
step07_structure        # 功能性描述
```

### 3. 变量命名统一

| 概念 | 旧命名(混乱) | 新命名(统一) |
|------|------------|------------|
| 基础路径 | `base_path`, `base_dir`, `data_root` | `base_path` |
| 输出路径 | `output_dir`, `output_base`, `result_dir` | `output_dir` |
| 数据框 | `df`, `data`, `dataset` | `df_<描述>` |
| 温度 | `temp`, `T`, `temperature` | `temperature` |
| 系统名 | `sys`, `system`, `sys_name` | `system_name` |

---

## 🚀 重构步骤计划

### 阶段1: 准备工作(30分钟)
1. ✅ **创建Git分支**
   ```powershell
   git checkout -b refactor/project-reorganization
   ```

2. ✅ **备份当前状态**
   ```powershell
   # 压缩整个workflow文件夹
   Compress-Archive -Path "C:\Users\11207\OneDrive\02_Code\work1-PtSnO\workflow" `
                    -DestinationPath "C:\Users\11207\OneDrive\02_Code\work1-PtSnO\workflow_backup_20251108.zip"
   ```

3. ✅ **创建新文件夹结构**
   ```powershell
   # 创建所有新文件夹
   $folders = @(
       "config",
       "scripts\step01_diffusion",
       "scripts\step02_ensemble",
       "scripts\step03_msd",
       "scripts\step04_diffusion_coeff",
       "scripts\step05_composition",
       "scripts\step06_energy",
       "scripts\step07_structure",
       "utils",
       "tools",
       "tests"
   )
   foreach ($folder in $folders) {
       New-Item -ItemType Directory -Path $folder -Force
   }
   ```

### 阶段2: 代码迁移(1-2小时)

#### 步骤2.1: 创建配置模块
**文件**: `config/paths.py`

```python
"""统一路径配置"""
from pathlib import Path

# 项目根目录
PROJECT_ROOT = Path(__file__).parent.parent

# 数据文件夹
DATA_DIR = PROJECT_ROOT / "data"
DATA_COORDINATION = DATA_DIR / "coordination"
DATA_GMX_MSD = DATA_DIR / "gmx_msd"
DATA_LAMMPS_ENERGY = DATA_DIR / "lammps_energy"
DATA_LINDEMANN = DATA_DIR / "lindemann"

# 结果文件夹
RESULTS_DIR = PROJECT_ROOT / "results"

# 当前数据版本(集中管理)
COORDINATION_RUN = "coordination_time_series_results_sample_20251106_214943"
ENERGY_MASTER = "energy_master_20251016_121110.csv"
LINDEMANN_RUN = "lindemann_master_run_20251025_205545.csv"

def get_coordination_data_path():
    """获取最新配位数据路径"""
    return DATA_COORDINATION / COORDINATION_RUN

def get_energy_data_path():
    """获取最新能量数据路径"""
    return DATA_LAMMPS_ENERGY / ENERGY_MASTER

def get_lindemann_data_path():
    """获取最新Lindemann数据路径"""
    return DATA_LINDEMANN / LINDEMANN_RUN
```

**文件**: `config/analysis_params.py`

```python
"""统一分析参数配置"""

# 温度范围
TEMPERATURE_RANGE = {
    'low': [300, 400, 500, 600],
    'medium': [700, 800, 900, 1000],
    'high': [1100, 1200, 1300, 1400]
}

# 系统系列
SYSTEM_SERIES = {
    'Pt8': ['pt8sn1-2-best', 'pt8sn2-1-best', 'pt8sn3-1-best', 
            'pt8sn4-1-best', 'pt8sn5-1-best', 'pt8sn6-1-best',
            'pt8sn7-1-best', 'pt8sn8-1-best', 'pt8sn9-1-best', 'pt8sn10-2-best'],
    'Pt6': ['pt6sn2-1-best', 'pt6sn4-1-best', 'pt6sn6-1-best', 'pt6sn8-1-best'],
    'PtSn8': ['pt1sn8-1-best', 'pt2sn8-2-best', 'pt3sn8-1-best', 
              'pt4sn8-1-best', 'pt6sn8-1-best']
}

# 分析窗口(ps)
ANALYSIS_WINDOWS = {
    'equilibration': (0, 1000),
    'production': (1000, 10000)
}

# 统计阈值
THRESHOLDS = {
    'outlier_std': 3.0,           # 异常值标准差
    'min_samples': 5,             # 最小样本数
    'cv_threshold': 0.15,         # 热容变异系数阈值
    'lindemann_melting': 0.10     # Lindemann熔化阈值
}
```

#### 步骤2.2: 创建工具模块
**文件**: `utils/data_loader.py`

```python
"""数据读取通用函数"""
import pandas as pd
from pathlib import Path
from typing import Optional, List

def load_csv_safe(file_path: Path, **kwargs) -> Optional[pd.DataFrame]:
    """安全读取CSV文件(自动处理编码)"""
    encodings = ['utf-8', 'gbk', 'gb2312', 'latin1']
    
    for encoding in encodings:
        try:
            df = pd.read_csv(file_path, encoding=encoding, **kwargs)
            print(f"✅ 成功读取: {file_path.name} (编码: {encoding})")
            return df
        except UnicodeDecodeError:
            continue
        except Exception as e:
            print(f"❌ 读取失败: {file_path.name} - {e}")
            return None
    
    print(f"❌ 所有编码尝试失败: {file_path.name}")
    return None

def load_coordination_data(base_path: Path, system: str, temp: str) -> Optional[pd.DataFrame]:
    """读取配位数时间序列数据"""
    file_path = base_path / system / f"{system}-{temp}-coord_time_series.csv"
    return load_csv_safe(file_path)

def filter_by_time_range(df: pd.DataFrame, time_col: str, 
                         start: float, end: float) -> pd.DataFrame:
    """按时间范围过滤数据"""
    return df[(df[time_col] >= start) & (df[time_col] <= end)].copy()
```

**文件**: `utils/plotting.py`

```python
"""绘图通用函数"""
import matplotlib.pyplot as plt
from pathlib import Path
from typing import List, Tuple

# 统一绘图风格
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['figure.dpi'] = 300

def setup_figure(nrows: int = 1, ncols: int = 1, 
                 figsize: Tuple[int, int] = (10, 6)) -> Tuple:
    """创建标准图形布局"""
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, dpi=300)
    return fig, axes

def save_figure(fig, output_path: Path, formats: List[str] = ['png']):
    """保存图形(多种格式)"""
    for fmt in formats:
        save_path = output_path.with_suffix(f'.{fmt}')
        fig.savefig(save_path, bbox_inches='tight', dpi=300)
        print(f"✅ 已保存: {save_path}")
```

#### 步骤2.3: 移动并更新脚本

**示例迁移**: `step7-6-3_q6_stats_comparison.py` → `scripts/step07_structure/q6_stats_comparison.py`

```python
# 旧导入方式
import sys
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

# 新导入方式(使用配置和工具模块)
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))  # 添加项目根目录到路径

from config.paths import get_coordination_data_path, RESULTS_DIR
from config.analysis_params import SYSTEM_SERIES, TEMPERATURE_RANGE
from utils.data_loader import load_coordination_data, load_csv_safe
from utils.plotting import setup_figure, save_figure

# 旧路径定义
base_path = Path(__file__).parent / "data" / "coordination" / "coordination_time_series_results_sample_20251106_214943"
output_dir = Path(__file__).parent / "results" / "step7.6_q6_stats"

# 新路径定义(使用配置)
base_path = get_coordination_data_path()
output_dir = RESULTS_DIR / "step07_structure" / "q6_stats"
```

### 阶段3: 测试验证(30分钟)

#### 测试清单
```powershell
# 测试所有步骤脚本是否能正常运行
cd scripts/step07_structure
python q6_stats_comparison.py --series Pt8 --temps "300K,900K" --systems "pt8sn1-2-best,pt8sn2-1-best"

# 检查输出文件
Test-Path "..\..\results\step07_structure\q6_stats\*"

# 验证所有导入
cd ..\..\tools
python check_dependencies.py
```

### 阶段4: 文档更新(30分钟)

1. **更新 README.md** (简化为项目总览)
2. **合并文档** → **WORKFLOW_GUIDE.md** (保留核心内容)
3. **创建 CHANGELOG.md** (记录重构历史)
4. **为每个脚本文件夹创建 README.md** (说明该步骤功能)

---

## ✅ 代码重构清单

### 提取公共函数示例

#### 重复代码1: 数据读取
**出现位置**: step7-5, step7-6-1, step7-6-2, step7-6-3

```python
# 旧代码(每个脚本重复)
def load_data(file_path):
    encodings = ['utf-8', 'gbk', 'latin1']
    for enc in encodings:
        try:
            return pd.read_csv(file_path, encoding=enc)
        except:
            continue
    return None

# 新代码(统一到 utils/data_loader.py)
from utils.data_loader import load_csv_safe
df = load_csv_safe(file_path)
```

#### 重复代码2: 路径处理
**出现位置**: 所有脚本

```python
# 旧代码(硬编码)
base_path = Path(__file__).parent / "data" / "coordination" / "coordination_time_series_results_sample_20251106_214943"

# 新代码(配置化)
from config.paths import get_coordination_data_path
base_path = get_coordination_data_path()
```

#### 重复代码3: 图形保存
**出现位置**: 所有绘图脚本

```python
# 旧代码(每个脚本重复)
output_file = output_dir / f"{prefix}.png"
fig.savefig(output_file, bbox_inches='tight', dpi=300)
print(f"保存图形: {output_file}")

# 新代码(统一到 utils/plotting.py)
from utils.plotting import save_figure
save_figure(fig, output_dir / prefix, formats=['png', 'pdf'])
```

---

## 🎯 实施优先级

### 🔴 高优先级(必须做)
1. ✅ **创建备份** (5分钟) - 防止数据丢失
2. ✅ **创建config/paths.py** (15分钟) - 统一路径管理
3. ✅ **迁移step07脚本** (30分钟) - 最常用的脚本
4. ✅ **测试step07功能** (15分钟) - 确保可用性

### 🟡 中优先级(建议做)
5. ⚠️ **创建utils模块** (1小时) - 提取公共函数
6. ⚠️ **迁移step01-06脚本** (1小时) - 完整迁移
7. ⚠️ **更新文档** (30分钟) - 同步文档
8. ⚠️ **创建run_script.py** (30分钟) - 统一启动器

### 🟢 低优先级(可选)
9. 💡 **编写单元测试** (2小时) - 提升质量
10. 💡 **配置CI/CD** (1小时) - 自动化测试
11. 💡 **添加类型提示** (1小时) - 代码规范

---

## ⚠️ 风险评估与备份

### 风险点
1. **路径错误**: 迁移后导入路径可能失效
   - **缓解**: 在每个__init__.py中添加路径调整代码
   
2. **依赖丢失**: 脚本间隐式依赖可能断裂
   - **缓解**: 先测试单个脚本,再测试完整流程
   
3. **数据覆盖**: 新旧结果文件夹混淆
   - **缓解**: 重命名旧results为results_old_20251108

### 备份策略
```powershell
# 三重备份
# 1. 压缩备份
Compress-Archive -Path "workflow" -DestinationPath "workflow_backup_$(Get-Date -Format 'yyyyMMdd_HHmmss').zip"

# 2. Git提交
git add .
git commit -m "备份: 重构前状态"

# 3. 复制到其他位置
Copy-Item -Path "workflow" -Destination "D:\Backups\workflow_$(Get-Date -Format 'yyyyMMdd')" -Recurse
```

### 回滚计划
```powershell
# 如果重构失败,快速回滚
git checkout main
git branch -D refactor/project-reorganization
```

---

## 📌 快速启动命令

### 创建脚本启动器: `tools/run_script.py`

```python
#!/usr/bin/env python3
"""统一脚本启动器"""
import sys
import subprocess
from pathlib import Path

# 脚本别名映射
SCRIPT_ALIASES = {
    's01': 'scripts/step01_diffusion/detect_outliers.py',
    's02': 'scripts/step02_ensemble/ensemble_analysis.py',
    's03': 'scripts/step03_msd/plot_msd.py',
    's04': 'scripts/step04_diffusion_coeff/calculate_ensemble_D.py',
    's05': 'scripts/step05_composition/analyze_sn_content.py',
    's06': 'scripts/step06_energy/energy_analysis.py',
    's07': 'scripts/step07_structure/lindemann_analysis.py',
    's07.5': 'scripts/step07_structure/unified_multi_temp.py',
    's07.6.1': 'scripts/step07_structure/temp_comparison_parallel.py',
    's07.6.2': 'scripts/step07_structure/temp_comparison_single.py',
    's07.6.3': 'scripts/step07_structure/q6_stats_comparison.py',
}

def main():
    if len(sys.argv) < 2:
        print("用法: python run_script.py <别名> [参数...]")
        print("\n可用别名:")
        for alias, path in SCRIPT_ALIASES.items():
            print(f"  {alias:8} -> {path}")
        return
    
    alias = sys.argv[1]
    if alias not in SCRIPT_ALIASES:
        print(f"❌ 未知别名: {alias}")
        return
    
    script_path = Path(__file__).parent.parent / SCRIPT_ALIASES[alias]
    args = sys.argv[2:]
    
    cmd = ['python', str(script_path)] + args
    print(f"🚀 运行: {' '.join(cmd)}")
    subprocess.run(cmd)

if __name__ == '__main__':
    main()
```

**使用示例**:
```powershell
# 旧方式(繁琐)
python step7-6-3_q6_stats_comparison.py --series Pt8 --temps "300K,900K" --systems "..."

# 新方式(简洁)
python tools/run_script.py s07.6.3 --series Pt8 --temps "300K,900K" --systems "..."
```

---

## 📖 参考资料

### 推荐阅读
1. **Python项目结构**: [The Hitchhiker's Guide to Python](https://docs.python-guide.org/writing/structure/)
2. **代码重构**: [Refactoring: Improving the Design of Existing Code](https://refactoring.com/)
3. **命名规范**: [PEP 8 -- Style Guide for Python Code](https://peps.python.org/pep-0008/)

### 工具推荐
1. **代码格式化**: `black` (自动格式化)
2. **代码检查**: `pylint` / `flake8`
3. **类型检查**: `mypy`
4. **文档生成**: `Sphinx`

---

## ✨ 预期收益

### 可维护性提升
- ✅ **查找脚本时间**: 从30秒 → 5秒(文件夹直接定位)
- ✅ **添加新功能**: 从复制粘贴 → 调用utils函数
- ✅ **修改配置**: 从全局搜索 → 修改config/paths.py

### 代码质量提升
- ✅ **代码复用率**: 从10% → 60%
- ✅ **命名一致性**: 从混乱 → 统一规范
- ✅ **文档覆盖率**: 从50% → 90%

### 协作效率提升
- ✅ **新人上手时间**: 从2小时 → 30分钟(文档清晰)
- ✅ **Bug修复时间**: 从30分钟 → 10分钟(结构清晰)
- ✅ **功能扩展时间**: 从1小时 → 20分钟(工具完善)

---

## 🎉 结语

这份指南为您的项目重构提供了完整的路线图。建议分阶段实施:

1. **第一次**: 只做高优先级任务(2小时)
2. **第二次**: 完成中优先级任务(3小时)
3. **第三次**: 根据需要添加低优先级功能(弹性)

**记住**: 重构不是一次性工作,而是持续改进的过程。每次修改都记得:
- 📝 更新文档
- ✅ 运行测试
- 💾 提交Git
- 🏷️ 打版本标签

祝重构顺利! 🚀

---

**维护者**: AI Assistant  
**最后更新**: 2025-11-08  
**版本**: v1.0
