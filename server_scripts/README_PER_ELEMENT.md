# 分元素Lindemann指数计算 - 使用指南

## 📋 概述

本方案在原有批处理脚本基础上，增加了**分元素Lindemann指数计算**功能。

### 新增功能

在 `lindemann_comparison_run_*.csv` 中增加以下列：

| 原有列 | 新增列 |
|--------|--------|
| Cluster_Unwrapped (PtSn团簇) | **Pt** - Pt元素单独的Lindemann指数 |
| Cluster_Wrapped | **Sn** - Sn元素单独的Lindemann指数 |
| 差异、差异% | **PtSn** - PtSn团簇整体（验证用） |
| | **PtSnO** - 包含O的整体系统 |

### 输出CSV格式

```csv
目录,结构,温度(K),Cluster_Unwrapped,Cluster_Wrapped,差异,差异%,Pt,Sn,PtSn,PtSnO,时间戳
```

## 📦 文件列表

生成的文件：
1. **lindemann_per_element_integrated.py** - 分元素计算核心脚本
2. **run_lindemann_batch_cluster_v8_per_element.sh** - v8批处理脚本（集成分元素）
3. **test_per_element_single.sh** - 单目录测试脚本

## 🚀 使用步骤

### 步骤1：上传脚本到服务器

```bash
# 在本地执行，上传到服务器
scp -P 2002 lindemann_per_element_integrated.py \
    jychen@211.86.151.148:/home/scms/jychen/tools/cp2k/md/msd/nnmD/dp-md/20251009/lin-new/

scp -P 2002 run_lindemann_batch_cluster_v8_per_element.sh \
    jychen@211.86.151.148:/home/scms/jychen/tools/cp2k/md/msd/nnmD/dp-md/20251009/lin-new/

scp -P 2002 test_per_element_single.sh \
    jychen@211.86.151.148:/home/scms/jychen/tools/cp2k/md/msd/nnmD/dp-md/20251009/lin-new/
```

### 步骤2：服务器上设置权限

```bash
ssh -p 2002 jychen@211.86.151.148
cd /home/scms/jychen/tools/cp2k/md/msd/nnmD/dp-md/20251009/lin-new

chmod +x run_lindemann_batch_cluster_v8_per_element.sh
chmod +x test_per_element_single.sh
```

### 步骤3：单目录测试

```bash
# 测试单个目录，验证功能
bash test_per_element_single.sh
```

测试成功后会看到：
```
✅ 计算成功！

结果文件内容:
Pt,Sn,PtSn,PtSnO
0.045123,0.052341,0.048567,0.051234
```

### 步骤4：批量运行

#### 方案A：直接运行（前台）

```bash
bash run_lindemann_batch_cluster_v8_per_element.sh
```

#### 方案B：后台运行（推荐）

修改 `start_full_batch_on_node30.sh` 中的脚本名称：

```bash
# 原来：
nohup bash run_lindemann_batch_cluster_v7_fixed.sh > run_v7_full_$(date +%Y%m%d_%H%M%S).log 2>&1 &

# 改为：
nohup bash run_lindemann_batch_cluster_v8_per_element.sh > run_v8_full_$(date +%Y%m%d_%H%M%S).log 2>&1 &
```

然后启动：
```bash
bash start_full_batch_on_node30.sh
```

#### 方案C：小批量测试（推荐先做）

创建一个测试用的目录列表：
```bash
# 从原列表中提取前10个目录
head -10 dir_list2.txt > dir_list_test10.txt

# 修改脚本中的 DIR_LIST_FILE
# DIR_LIST_FILE="dir_list_test10.txt"

bash run_lindemann_batch_cluster_v8_per_element.sh
```

## 📊 结果查看

### 查看实时进度

```bash
tail -f run_v8_full_*.log
```

### 查看comparison结果（含分元素）

```bash
# 查看最新的comparison文件
ls -lt collected_lindemann_cluster/lindemann_comparison_*.csv | head -1

# 查看前10行
head -11 collected_lindemann_cluster/lindemann_comparison_run_*.csv
```

### 统计分析

```bash
# 计算各元素平均Lindemann指数
tail -n +2 collected_lindemann_cluster/lindemann_comparison_run_*.csv | awk -F',' '
BEGIN {
    print "元素统计:"
}
{
    cluster+=$4; pt+=$8; sn+=$9; ptsn+=$10; ptsno+=$11; n++
}
END {
    print "样本数: " n
    print "Cluster(原方法): " cluster/n
    print "Pt:             " pt/n
    print "Sn:             " sn/n
    print "PtSn:           " ptsn/n
    print "PtSnO:          " ptsno/n
}'
```

## 🔍 性能优化

分元素计算采用了以下优化：
1. **并行unwrap** - 每个元素独立unwrap，充分利用numpy向量化
2. **快速计算** - 简化输出，只保留核心数值
3. **stdout输出** - 避免文件I/O开销，直接传递给shell

预计额外耗时：**每个目录增加 10-20秒**（取决于原子数）

## ⚙️ 配置选项

在 `run_lindemann_batch_cluster_v8_per_element.sh` 中：

```bash
# 是否启用分元素计算（可设置为false关闭）
ENABLE_PER_ELEMENT=true

# 分元素计算脚本路径
PER_ELEMENT_SCRIPT="/path/to/lindemann_per_element_integrated.py"
```

## 🐛 故障排查

### 问题1：分元素脚本未找到

**症状：** 日志中显示 "分元素计算失败，使用默认值0"

**解决：**
```bash
# 检查脚本是否存在
ls -l /home/scms/jychen/tools/cp2k/md/msd/nnmD/dp-md/20251009/lin-new/lindemann_per_element_integrated.py

# 检查权限
chmod +x lindemann_per_element_integrated.py
```

### 问题2：MDAnalysis错误

**症状：** ImportError: No module named 'MDAnalysis'

**解决：**
```bash
# 确认conda环境
conda activate mda_env

# 检查MDAnalysis
python3 -c "import MDAnalysis; print(MDAnalysis.__version__)"
```

### 问题3：分元素结果全为0

**可能原因：**
- 元素选择表达式错误
- XYZ文件中元素名称不匹配（如 'PT' vs 'Pt'）

**调试：**
```bash
# 查看XYZ文件前几行
head -20 sampling-simply.xyz

# 手动测试
python3 lindemann_per_element_integrated.py \
    --coord sampling-simply.xyz \
    --verbose
```

## 📈 与下游分析集成

生成的CSV文件可以直接用于后续分析：

```python
import pandas as pd

# 读取结果
df = pd.read_csv('collected_lindemann_cluster/lindemann_comparison_run_*.csv')

# 提取分元素数据
df_elements = df[['结构', '温度(K)', 'Pt', 'Sn', 'PtSn', 'PtSnO']]

# 按温度分组统计
df_elements.groupby('温度(K)').mean()

# 绘图对比
import matplotlib.pyplot as plt
df_elements.groupby('温度(K)')[['Pt', 'Sn', 'PtSn', 'PtSnO']].mean().plot()
plt.ylabel('Lindemann Index')
plt.title('Lindemann Index by Element vs Temperature')
plt.savefig('lindemann_by_element.png')
```

## 📝 技术细节

### 计算方法差异

| 方法 | 原子选择 | 说明 |
|------|----------|------|
| **Cluster_Unwrapped** | `name Pt or name Sn` | 原方法，团簇整体 |
| **Pt** | `name Pt` | 只计算Pt原子之间的距离涨落 |
| **Sn** | `name Sn` | 只计算Sn原子之间的距离涨落 |
| **PtSn** | `name Pt or name Sn` | 与Cluster相同（验证用） |
| **PtSnO** | `name Pt or name Sn or name O` | 包含O原子的全系统 |

### Unwrap说明

每个元素独立进行unwrap处理，确保边界条件正确处理。

## ✅ 检查清单

部署前确认：
- [ ] 已上传 `lindemann_per_element_integrated.py`
- [ ] 已上传 `run_lindemann_batch_cluster_v8_per_element.sh`
- [ ] 已设置可执行权限
- [ ] 已在单个目录测试成功
- [ ] 已在10个目录小批量测试
- [ ] 已修改 `start_full_batch_on_node30.sh`（如需要）
- [ ] 已确认磁盘空间充足

## 🎯 预期结果

成功运行后，comparison CSV 会包含完整的分元素信息：

```csv
目录,结构,温度(K),Cluster_Unwrapped,Cluster_Wrapped,差异,差异%,Pt,Sn,PtSn,PtSnO,时间戳
/path/to/dir,Pt6Sn8O4,300,0.048567,0.051234,-0.002667,-5.21,0.045123,0.052341,0.048567,0.051234,2025-11-11 19:30:15
```

可用于：
- ✅ 对比不同元素的流动性差异
- ✅ 分析温度对各元素的影响
- ✅ 研究Pt-Sn相互作用
- ✅ 评估O原子对整体动力学的贡献
