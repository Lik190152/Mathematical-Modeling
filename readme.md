# 赛题A分析与建模工程

这个工程把赛题 A 的三问串成一条可复现流水线：

`原始 xlsx -> 数据清洗与配准 -> Q1 温度预测 -> Q2 CO 预测 -> Q3 负压优化 -> 报告表格/图表`

## 目录结构

- `src/sinter_co/`: 数据处理、特征工程、模型训练、优化核心逻辑
- `scripts/`: 六个按顺序执行的入口脚本
- `tests/`: 最小回归测试
- `artifacts/`: 中间数据、模型和指标输出
- `reports/`: 论文可直接使用的图表和表格

## 依赖

推荐使用带 Anaconda 科学计算栈的 Python 环境。当前实现依赖见 [requirements.txt](/d:/数学建模/赛题A/requirements.txt:1)。

如果 `to_parquet()` 缺少 `pyarrow`/`fastparquet`，脚本会自动回退为 `.csv` 输出，不会中断流程。

## 运行顺序

在仓库根目录执行：

```powershell
$env:Path='D:\Anaconda3;D:\Anaconda3\Library\bin;D:\Anaconda3\Library\mingw-w64\bin;D:\Anaconda3\Library\usr\bin;D:\Anaconda3\Scripts;' + $env:Path
D:\Anaconda3\python.exe scripts\01_profile_data.py
D:\Anaconda3\python.exe scripts\02_build_features.py
D:\Anaconda3\python.exe scripts\03_train_q1.py
D:\Anaconda3\python.exe scripts\04_train_q2.py
D:\Anaconda3\python.exe scripts\05_optimize_q3.py
D:\Anaconda3\python.exe scripts\06_generate_report_tables.py
```

所有脚本都支持：

- `--input`: 原始 Excel 路径，默认 `附件1.原始数据.xlsx`
- `--artifacts-dir`: 中间产物目录，默认 `artifacts`
- `--seed`: 随机种子，默认 `42`

`05_optimize_q3.py` 额外支持：

- `--row`: 选择优化的样本行，默认 `-1` 表示最后一行

## 主要输出

- `artifacts/processed/`: 基础数据与清洗后数据
- `artifacts/features/`: Q1/Q2 特征和 lag 摘要
- `artifacts/models/`: Q1 风箱模型和 Q2 CO 模型
- `artifacts/metrics/`: 训练指标和 Q3 优化结果
- `reports/figures/`: 论文插图
- `reports/tables/`: 数据字典、lag 表、指标表、优化表

## 测试

```powershell
$env:Path='D:\Anaconda3;D:\Anaconda3\Library\bin;D:\Anaconda3\Library\mingw-w64\bin;D:\Anaconda3\Library\usr\bin;D:\Anaconda3\Scripts;' + $env:Path
D:\Anaconda3\python.exe -m pytest tests\test_pipeline.py -q
```
