# 🤖 代码 AI 处理最佳实践指南

**基于实际项目经验总结**  
**日期**: 2025-11-08  
**场景**: Python 多脚本项目迁移与路径重构

---

## 📋 核心原则

### 1. **文件编码问题 - 最高优先级 ⚠️**

#### ❌ 常见问题
```python
# 复制脚本时可能出现的编码混乱:
# ===== 鍏ㄥ眬閰嶇疆 =====  # ❌ 乱码 (应该是"全局配置")
BASE_DIR = Path(__file__).parent  # workflow鐩綍  # ❌ 乱码 (应该是"workflow目录")
```

#### ✅ 正确做法
1. **复制前检查编码**:
   ```powershell
   # 使用 Python 检测编码
   python -c "import chardet; print(chardet.detect(open('script.py','rb').read()))"
   ```

2. **保持 UTF-8 编码**:
   ```python
   # 文件头部必须包含
   # -*- coding: utf-8 -*-
   ```

3. **避免编码转换**:
   - PowerShell `Get-Content` 默认使用 GBK，会破坏 UTF-8 中文
   - 使用 `-Encoding UTF8` 参数或直接用 Python 处理

#### 🔧 修复乱码脚本
```powershell
# PowerShell 批量修复
$content = Get-Content "script.py" -Raw -Encoding UTF8
$fixed = $content -replace '鍏ㄥ眬閰嶇疆', '全局配置' `
                  -replace 'workflow鐩綍', 'workflow目录'
[System.IO.File]::WriteAllText("script.py", $fixed, [System.Text.UTF8Encoding]::new($false))
```

---

## 📂 路径配置问题

### 2. **相对路径优于绝对路径**

#### ❌ 错误示例
```python
# 使用硬编码绝对路径 - 不可移植
DATA_DIR = Path(r'd:\OneDrive\py\Cv\lin\data')
BASE_DIR = Path(__file__).parent.parent  # 假设错误的目录层级
```

#### ✅ 正确做法
```python
# 使用相对路径 - 可移植
BASE_DIR = Path(__file__).parent  # 脚本所在目录
DATA_DIR = BASE_DIR / 'data'      # 数据子目录
RESULTS_DIR = BASE_DIR / 'results'  # 输出子目录
```

### 3. **目录结构变化需同步更新**

#### 场景: 从旧结构迁移
```
旧结构:                      新结构:
v3_simplified_workflow/      workflow/
├── scripts/                 ├── step*.py (直接在根目录)
│   └── step*.py            ├── data/
└── files/                   └── results/
    └── data/
```

#### 提示词模板
```
当我将脚本从子目录移到根目录时:
1. 检查所有 Path(__file__).parent 的层级是否正确
2. 将 parent.parent 改为 parent (如果脚本层级减少)
3. 更新所有硬编码的数据路径为相对路径
4. 验证 BASE_DIR 定义在使用前 (避免 NameError)
```

---

## 🔍 AI 交互提示词

### 4. **路径诊断提示词**

```
请帮我分析以下问题:
1. 检查所有脚本中的 BASE_DIR 定义是否指向正确目录
2. 列出所有使用绝对路径的地方
3. 验证数据路径是否存在
4. 比较旧结构和新结构的差异

项目结构:
[粘贴目录树]

错误信息:
[粘贴完整错误]
```

### 5. **批量修复提示词**

```
请帮我批量修复以下问题:

问题类型: [路径配置错误/编码问题/导入错误]

受影响的文件:
- script1.py: [具体问题]
- script2.py: [具体问题]

修复要求:
1. 保持代码逻辑不变
2. 只修改配置部分
3. 生成修复前后对比
4. 提供验证命令

示例:
旧: BASE_DIR = Path(__file__).parent.parent
新: BASE_DIR = Path(__file__).parent
```

### 6. **输出验证提示词**

```
请创建一个脚本来验证:
1. 所有脚本是否生成了预期的输出文件
2. 输出文件的完整性 (大小、数量)
3. 生成检查报告并标记缺失/异常的输出

预期输出清单:
[列出所有应生成的文件]
```

---

## 🛠️ 实战技巧

### 7. **分阶段验证**

```
Step 1: 路径验证
✅ 数据文件存在性检查
✅ BASE_DIR 指向正确
✅ 所有导入路径正确

Step 2: 小范围测试
✅ 运行 1-2 个脚本验证
✅ 检查输出是否符合预期

Step 3: 完整测试
✅ 运行所有脚本
✅ 检查输出完整性
```

### 8. **错误信息完整性**

#### ❌ 不完整的问题描述
```
"脚本报错了，怎么办？"
```

#### ✅ 完整的问题描述
```
脚本: step3_plot_msd.py
错误: FileNotFoundError: 'c:\...\results\ensemble_analysis_results.csv'
当前目录: C:\Users\...\workflow
已检查:
- 文件不存在于该路径
- BASE_DIR = Path(__file__).parent
期望: 读取 workflow/results/ 下的文件
```

---

## 📝 AI 工具选择建议

### 9. **何时使用不同工具**

| 任务类型 | 推荐工具 | 原因 |
|---------|---------|------|
| 路径检查 | `file_search` + `grep_search` | 快速定位所有路径配置 |
| 批量修改 | `replace_string_in_file` | 精确替换，避免破坏代码 |
| 编码问题 | Python 脚本 | PowerShell 可能破坏 UTF-8 |
| 输出验证 | 自定义脚本 | 自动化检查更可靠 |
| 测试运行 | `run_in_terminal` | 实时查看输出和错误 |

### 10. **PowerShell vs Python 选择**

```
使用 PowerShell:
✅ 文件操作 (复制、移动、删除)
✅ 目录遍历和统计
✅ 简单的文本替换 (ASCII)

使用 Python:
✅ UTF-8 中文处理
✅ 复杂的文本解析
✅ 编码检测和转换
✅ 数据验证和分析
```

---

## 🚨 常见陷阱

### 11. **需要特别注意的问题**

1. **BASE_DIR 定义顺序**
   ```python
   # ❌ 错误: 使用前未定义
   DATA_DIR = BASE_DIR / 'data'  # NameError!
   BASE_DIR = Path(__file__).parent
   
   # ✅ 正确: 先定义后使用
   BASE_DIR = Path(__file__).parent
   DATA_DIR = BASE_DIR / 'data'
   ```

2. **PowerShell 目录切换**
   ```powershell
   # ❌ cd 不会持久化到下一条命令
   cd workflow
   python script.py  # 仍在原目录运行!
   
   # ✅ 使用绝对路径或组合命令
   cd workflow; python script.py
   # 或
   python C:\full\path\to\script.py
   ```

3. **输出目录自动创建**
   ```python
   # ✅ 确保输出目录存在
   OUTPUT_DIR = BASE_DIR / 'results'
   OUTPUT_DIR.mkdir(exist_ok=True, parents=True)
   ```

---

## 🎯 提示词模板库

### 模板 1: 路径重构
```
我需要将以下脚本从 [旧路径] 迁移到 [新路径]:

脚本列表:
- script1.py
- script2.py

当前问题:
- BASE_DIR 指向错误
- 数据路径使用绝对路径

请:
1. 分析所有路径配置
2. 生成修复方案
3. 批量修改文件
4. 提供验证命令
```

### 模板 2: 编码修复
```
以下脚本出现中文乱码:

症状: "鍏ㄥ眬閰嶇疆" 应该是 "全局配置"

请:
1. 检测文件实际编码
2. 识别所有乱码模式
3. 批量修复为正确中文
4. 确保保存为 UTF-8 BOM
```

### 模板 3: 输出检查
```
请创建验证脚本检查:

预期输出:
- Step 1: 5个文件 (CSV + PNG)
- Step 2: 2个CSV文件
- Step 3: 10+ PNG图表

当前状态: [未知/部分完成]

生成报告显示:
- 完整性百分比
- 缺失的文件列表
- 异常大小的文件
```

---

## 📖 参考资源

### 有用的验证命令

```powershell
# 1. 检查文件编码
python -c "import chardet; print(chardet.detect(open('file.py','rb').read()))"

# 2. 统计输出文件
Get-ChildItem results -Recurse -File | Measure-Object | Select-Object Count

# 3. 查找路径配置
Get-Content script.py | Select-String "BASE_DIR|Path.*parent"

# 4. 验证数据文件存在
Test-Path "workflow\data\gmx_msd\*" -PathType Leaf

# 5. 检查脚本语法
python -m py_compile script.py
```

---

## ✅ 最佳实践清单

在请求 AI 处理代码时，确保提供:

- [ ] 完整的错误信息 (包括堆栈跟踪)
- [ ] 当前目录结构 (`tree /F` 或 `ls -R`)
- [ ] 相关文件的关键代码片段 (不只是出错行)
- [ ] 已尝试的解决方案
- [ ] 预期行为 vs 实际行为
- [ ] 环境信息 (Python版本、OS、编码设置)

在批量修改前，记得:

- [ ] 备份原始文件
- [ ] 小范围测试修复方案
- [ ] 逐步验证每个修改
- [ ] 保留修复记录文档

---

## 🎓 经验总结

从本次项目学到的关键教训:

1. **编码问题比路径问题更隐蔽** - 看起来正常但实际已损坏
2. **PowerShell 不是处理 UTF-8 的好工具** - 优先使用 Python
3. **相对路径 + 自包含结构 = 可移植性** - 避免硬编码
4. **分步验证比一次性修复更安全** - 降低风险
5. **自动化检查脚本很重要** - 节省人工验证时间

---

**记住**: AI 是强大的助手，但最终验证和决策需要人类完成！

**建议工作流**:
```
AI 分析 → 人工审查 → 小范围测试 → AI 批量修复 → 全面验证 → 文档记录
```
