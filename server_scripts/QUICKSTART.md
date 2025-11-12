# 🚀 快速启动指南 - 分元素Lindemann指数计算

## 一键部署和运行

### 步骤1️⃣：部署到服务器（本地PowerShell执行）

```powershell
# 进入脚本目录
cd c:\Users\11207\OneDrive\02_Code\work1-PtSnO\workflow\server_scripts

# 一键部署（需要bash环境，如Git Bash）
bash deploy_per_element_to_server.sh
```

或者手动上传：

```powershell
# 上传核心脚本
scp -P 2002 lindemann_per_element_integrated.py jychen@211.86.151.148:/home/scms/jychen/tools/cp2k/md/msd/nnmD/dp-md/20251009/lin-new/

# 上传批处理脚本
scp -P 2002 run_lindemann_batch_cluster_v8_per_element.sh jychen@211.86.151.148:/home/scms/jychen/tools/cp2k/md/msd/nnmD/dp-md/20251009/lin-new/

# 上传测试脚本
scp -P 2002 test_per_element_single.sh jychen@211.86.151.148:/home/scms/jychen/tools/cp2k/md/msd/nnmD/dp-md/20251009/lin-new/
```

### 步骤2️⃣：连接服务器测试

```bash
# SSH连接
ssh -p 2002 jychen@211.86.151.148

# 进入工作目录
cd /home/scms/jychen/tools/cp2k/md/msd/nnmD/dp-md/20251009/lin-new

# 设置权限
chmod +x run_lindemann_batch_cluster_v8_per_element.sh
chmod +x test_per_element_single.sh

# 激活conda环境
conda activate mda_env

# 测试单个目录
bash test_per_element_single.sh
```

**预期看到：**
```
✅ 计算成功！

结果文件内容:
Pt,Sn,PtSn,PtSnO
0.045123,0.052341,0.048567,0.051234
```

### 步骤3️⃣：小批量测试（10个目录）

```bash
# 创建测试列表（前10个目录）
head -10 dir_list2.txt > dir_list_test10.txt

# 编辑v8脚本，改为使用测试列表
nano run_lindemann_batch_cluster_v8_per_element.sh
# 修改第3行：DIR_LIST_FILE="dir_list_test10.txt"

# 运行小批量测试
bash run_lindemann_batch_cluster_v8_per_element.sh
```

**预期：** 约5-10分钟完成10个目录

**检查结果：**
```bash
# 查看生成的comparison文件
ls -lh collected_lindemann_cluster/lindemann_comparison_run_*.csv

# 查看前几行（包含表头和数据）
head -5 collected_lindemann_cluster/lindemann_comparison_run_*.csv
```

应该看到：
```csv
目录,结构,温度(K),Cluster_Unwrapped,Cluster_Wrapped,差异,差异%,Pt,Sn,PtSn,PtSnO,时间戳
/path/to/dir1,Pt6Sn8O4,300,0.048567,0.051234,-0.002667,-5.21,0.045123,0.052341,0.048567,0.051234,2025-11-11 19:30:15
```

✅ 看到 **Pt, Sn, PtSn, PtSnO** 列就成功了！

### 步骤4️⃣：完整批量运行（3273个目录）

```bash
# 方案A：修改启动脚本（推荐）
nano start_full_batch_on_node30.sh

# 找到这行：
nohup bash run_lindemann_batch_cluster_v7_fixed.sh > run_v7_full_$(date +%Y%m%d_%H%M%S).log 2>&1 &

# 改为：
nohup bash run_lindemann_batch_cluster_v8_per_element.sh > run_v8_full_$(date +%Y%m%d_%H%M%S).log 2>&1 &

# 保存后启动
bash start_full_batch_on_node30.sh
```

```bash
# 方案B：直接在v8脚本中改回完整列表
nano run_lindemann_batch_cluster_v8_per_element.sh
# 确保第3行：DIR_LIST_FILE="dir_list2.txt"

# 后台运行
nohup bash run_lindemann_batch_cluster_v8_per_element.sh > run_v8_full_$(date +%Y%m%d_%H%M%S).log 2>&1 &

# 记录进程ID
echo $! > .batch_v8.pid
```

### 步骤5️⃣：监控进度

```bash
# 实时查看日志
tail -f run_v8_full_*.log

# 查看进度（每10秒更新）
# 日志中会显示：
# 进度: 150/3273 (4.6%) | ✅ 148 | ⏭️ 0 | ❌ 2

# 检查进程是否还在运行
ps aux | grep run_lindemann_batch

# 查看已完成数量
wc -l collected_lindemann_cluster/lindemann_comparison_run_*.csv
```

### 步骤6️⃣：查看结果

```bash
# 完成后查看统计
tail -50 run_v8_full_*.log

# 查看分元素平均值
tail -n +2 collected_lindemann_cluster/lindemann_comparison_run_*.csv | awk -F',' '
{
    cluster+=$4; pt+=$8; sn+=$9; ptsn+=$10; ptsno+=$11; n++
}
END {
    printf "样本数: %d\n", n
    printf "Cluster: %.6f\n", cluster/n
    printf "Pt:      %.6f\n", pt/n
    printf "Sn:      %.6f\n", sn/n
    printf "PtSn:    %.6f\n", ptsn/n
    printf "PtSnO:   %.6f\n", ptsno/n
}'
```

## 📊 数据下载到本地分析

### 下载CSV文件

```powershell
# 在本地PowerShell执行
scp -P 2002 jychen@211.86.151.148:/home/scms/jychen/tools/cp2k/md/msd/nnmD/dp-md/20251009/lin-new/collected_lindemann_cluster/lindemann_comparison_run_*.csv c:\Users\11207\OneDrive\02_Code\work1-PtSnO\workflow\results\
```

### Python分析

```python
import pandas as pd
import matplotlib.pyplot as plt

# 读取数据
df = pd.read_csv('results/lindemann_comparison_run_20251111_190000.csv')

# 基本统计
print("分元素Lindemann指数统计:")
print(df[['Pt', 'Sn', 'PtSn', 'PtSnO']].describe())

# 按温度分组
temp_analysis = df.groupby('温度(K)')[['Pt', 'Sn', 'PtSn', 'PtSnO']].mean()
print("\n各温度下的平均Lindemann指数:")
print(temp_analysis)

# 绘图
temp_analysis.plot(kind='line', marker='o', figsize=(12, 6))
plt.ylabel('Lindemann Index')
plt.xlabel('Temperature (K)')
plt.title('Lindemann Index by Element vs Temperature')
plt.legend(['Pt only', 'Sn only', 'PtSn cluster', 'PtSnO all'])
plt.grid(True, alpha=0.3)
plt.savefig('lindemann_elements_vs_temp.png', dpi=300, bbox_inches='tight')
plt.show()

# 相关性分析
print("\n元素间相关性:")
print(df[['Pt', 'Sn', 'PtSn', 'PtSnO']].corr())
```

## 🎯 关键输出文件

运行完成后会生成：

```
collected_lindemann_cluster/
├── lindemann_master_run_20251111_190000.csv      # 主表（简化）
├── lindemann_comparison_run_20251111_190000.csv  # 对比表（含分元素）⭐
├── convergence_master_run_20251111_190000.csv    # 收敛性分析
├── summary_20251111_190000.log                    # 详细日志
└── error_20251111_190000.log                      # 错误日志
```

**最重要的是 `lindemann_comparison_run_*.csv`** - 包含所有分元素数据！

## ⏱️ 预期时间

- **单个目录测试：** ~1分钟
- **10个目录测试：** ~5-10分钟
- **完整3273个目录：** ~12-15小时（并行16任务）

## ✅ 成功标志

1. ✅ 测试脚本输出4个Lindemann值
2. ✅ comparison CSV包含 Pt, Sn, PtSn, PtSnO 列
3. ✅ 数值合理（通常 0.01-0.15 范围）
4. ✅ Pt ≈ Sn（数量级相近）
5. ✅ PtSn ≈ Cluster_Unwrapped（验证正确）

## ❌ 常见问题

### 问题1：ImportError: No module named 'MDAnalysis'

```bash
# 确认conda环境
conda activate mda_env
python3 -c "import MDAnalysis; print('OK')"
```

### 问题2：分元素结果全为0

```bash
# 检查XYZ文件
head -20 sampling-simply.xyz

# 确认元素名称大小写（Pt vs PT）
grep -E "^(Pt|Sn|O) " sampling-simply.xyz | head -5
```

### 问题3：脚本找不到

```bash
# 确认文件存在
ls -l /home/scms/jychen/tools/cp2k/md/msd/nnmD/dp-md/20251009/lin-new/lindemann_per_element_integrated.py

# 确认权限
chmod +x *.sh
```

## 📞 需要帮助？

1. 查看详细文档：`README_PER_ELEMENT.md`
2. 查看实现总结：`IMPLEMENTATION_SUMMARY.md`
3. 检查日志文件：`error_*.log` 和 `summary_*.log`

---

**准备好了就开始吧！** 🚀

第一步：`bash deploy_per_element_to_server.sh`
