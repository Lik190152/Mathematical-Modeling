# 冠脉狭窄分析项目

这是一个面向赛题 C 的可复现项目，用于从三维冠脉二值 `mask` 数据中完成以下工作：

- 建立血管结构模型
- 识别并量化狭窄候选病变
- 设计局部与全局观察策略
- 生成论文所需的表格、图片和正文材料

## 项目目录

```text
.
|-- case1/ , case2/                  输入的 NIfTI 冠脉 mask
|-- coronary_analysis/               核心分析代码
|-- tools/                           可执行入口脚本
|-- tests/                           单元测试与集成测试
|-- analysis_results/                生成的结果表和图片
|-- manuscript/                      论文正文源稿与图表清单
|-- requirements.txt                 Python 依赖
|-- run.bat                          Windows 一键启动脚本
`-- run_analysis.bat                 启动脚本别名
```

## 依赖安装

先安装项目依赖：

```bash
pip install -r requirements.txt
```

`requirements.txt` 中固定了当前项目已验证可运行的版本。

## 一键运行

### 方式 1：Windows 直接双击

直接双击：

```text
run_analysis.bat
```

或双击：

```text
run.bat
```

### 方式 2：命令行运行

```cmd
cmd /c run_analysis.bat --no-pause
```

常用参数：

```cmd
cmd /c run_analysis.bat --dry-run
cmd /c run_analysis.bat --no-pause
```

这几个用法的区别如下：

| 用法 | 会不会真正开始分析 | 执行结束后是否暂停 | 适合什么场景 |
| --- | --- | --- | --- |
| `run_analysis.bat` 或 `run.bat`（双击/不带参数） | 会 | 会 | 适合直接双击运行，便于查看执行结果和报错信息 |
| `cmd /c run_analysis.bat --dry-run` | 不会 | 不涉及 | 只打印将要执行的 Python 命令，不真正执行分析，用来检查脚本是否找到了正确的解释器和分析入口 |
| `cmd /c run_analysis.bat --no-pause` | 会 | 不会 | 执行完整分析，执行结束后不暂停，适合命令行、脚本调用或自动化运行 |

如果你只是想确认脚本会调用哪个 Python 和哪个入口文件，先用 `--dry-run`。  
如果你要真正跑完整分析并且不想在结束后手动按键关闭窗口，用 `--no-pause`。  
如果你是双击运行，希望在执行结束后保留窗口查看日志，则直接双击 `run_analysis.bat` 或 `run.bat` 即可。

## 启动脚本的 Python 查找逻辑

为了避免把解释器路径写死在脚本里，`run.bat` 采用以下优先级查找 Python：

1. 优先使用项目根目录下的 `.venv\Scripts\python.exe`
2. 如果没有 `.venv`，则使用系统 `PATH` 中的 `python`
3. 如果没有 `python`，则尝试使用 `py -3`

因此，不再依赖固定路径如 `D:\Anaconda3\python.exe`。

如果你使用自己的虚拟环境，推荐在项目根目录创建 `.venv`，这样双击脚本即可直接运行。

## 直接运行分析入口

如果不想走批处理脚本，也可以直接运行：

```bash
python tools/run_coronary_analysis.py
```

该脚本会自动扫描工作区下的 `case*/**/mask.nii.gz`，执行完整分析流程，并将输出写入 `analysis_results`。

## 运行后会生成什么

执行 `tools/run_coronary_analysis.py`、`run.bat` 或 `run_analysis.bat` 后，会生成：

### 表格

- `analysis_results/tables/basic_structure_stats.csv`
- `analysis_results/tables/lesion_candidates.csv`
- `analysis_results/tables/local_best_views.csv`
- `analysis_results/tables/global_view_plans.csv`
- `analysis_results/tables/sensitivity_analysis.csv`

### 图片

- `analysis_results/figures/figure1_workflow.png`
- `analysis_results/figures/figure2_mip_overview.png`
- `analysis_results/figures/figure3_centerline_overview.png`
- `analysis_results/figures/figure4_radius_curves.png`
- `analysis_results/figures/figure5_local_view_comparison.png`
- `analysis_results/figures/figure6_global_heatmap.png`

## 复现当前结果的步骤

1. 将题目给定的 NIfTI 数据放到 `case1/` 和 `case2/` 下
2. 安装依赖：

```bash
pip install -r requirements.txt
```

3. 运行完整分析：

```bash
python tools/run_coronary_analysis.py
```

或者：

```cmd
cmd /c run_analysis.bat --no-pause
```

4. 在 `analysis_results` 中查看生成的表格和图
5. 在 `manuscript` 中查看论文材料

## 论文材料位置

论文相关文件位于：

- `manuscript/赛题C_正文源稿.md`
- `manuscript/赛题C_图表与附录清单.md`

这些文件说明了如何将生成的图表与结果组织到最终论文中。

## 运行测试

如果你想验证代码是否正常，可执行：

```bash
python -m unittest tests.test_coronary_analysis -v
```

## 核心入口说明

- `tools/run_coronary_analysis.py`
  - 完整分析入口
- `tools/generate_nifti_mips.py`
  - 生成冠脉三视图 MIP 图
- `coronary_analysis/pipeline.py`
  - 分析主流程：预处理、中心线提取、病变识别、视角优化

## 补充说明

- 本项目处理的是三维二值冠脉 `mask`，不是原始 DSA 灰度造影序列。
- 当前模型评价侧重于几何一致性、可解释性、参数敏感性和视角覆盖效率，因为题目没有提供临床真值病变标注。
