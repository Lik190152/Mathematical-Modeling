# 赛题C 图表与附录清单

## 一、图件清单与图题建议

### 图1
- 文件：`analysis_results/figures/figure1_workflow.png`
- 建议图题：**图1 冠脉狭窄识别与观察策略优化总体流程**
- 正文落位：第 6 节“求解流程”首次引用
- 说明要点：
  - 从 `mask` 预处理出发；
  - 中间经过中心线与半径建模；
  - 再进行狭窄评分；
  - 最后输出局部与全局视角优化结果。

### 图2
- 文件：`analysis_results/figures/figure2_mip_overview.png`
- 建议图题：**图2 四支冠脉的三视图最大投影结果**
- 正文落位：第 7.1 节“血管结构建模结果分析”
- 说明要点：
  - 展示 `case1 left / case1 right / case2 left / case2 right` 的 `XY、XZ、YZ` 三视图；
  - 用于说明四支冠脉在空间形态、弯曲度和分支复杂性上的差异。

### 图3
- 文件：`analysis_results/figures/figure3_centerline_overview.png`
- 建议图题：**图3 冠脉中心线与原始血管结构叠加结果**
- 正文落位：第 7.1 节
- 说明要点：
  - 白色为血管轮廓，红色为中心线；
  - 用于说明中心线提取结果与原始血管主体吻合良好；
  - 证明血管拓扑模型具有可解释性。

### 图4
- 文件：`analysis_results/figures/figure4_radius_curves.png`
- 建议图题：**图4 典型病变分支的半径—狭窄率联合曲线**
- 正文落位：第 7.2 节“狭窄识别与量化结果分析”
- 说明要点：
  - 深蓝实线为局部半径；
  - 红色虚线为狭窄率；
  - 粉色阴影为最终识别出的病变区段；
  - 无主要病变的冠脉可在图中保留“no major lesion”说明。

### 图5
- 文件：`analysis_results/figures/figure5_local_view_comparison.png`
- 建议图题：**图5 主要病变的局部最佳观察视角对比**
- 正文落位：第 7.3 节“局部最佳观察策略结果分析”
- 说明要点：
  - 蓝色线段为分支中心线投影；
  - 红色线段为对应病变投影；
  - 每列对应一个局部最优视角；
  - 图题或正文中需解释 `az` 为方位角、`el` 为俯仰角。

### 图6
- 文件：`analysis_results/figures/figure6_global_heatmap.png`
- 建议图题：**图6 病例级平均局部视角评分热图**
- 正文落位：第 7.4 节“全局观察策略结果分析”
- 说明要点：
  - 横轴为方位角索引（每步 `5°`）；
  - 纵轴为俯仰角索引（每步 `5°`）；
  - 颜色越亮表示平均局部观察效果越好；
  - 用于说明优良视角呈连续角度带分布，而非单个孤立点。

## 二、表格清单与表题建议

### 表1
- 来源文件：`analysis_results/tables/basic_structure_stats.csv`
- 建议表题：**表1 四支冠脉的基础结构统计结果**
- 正文落位：第 7.1 节
- 建议保留列：
  - 病例、侧别；
  - 原始体积、重采样后体积、体积变化率；
  - 中心线节点数、边数、分支数；
  - 主要病变数。

### 表2
- 来源文件：`analysis_results/tables/lesion_candidates.csv`
- 建议表题：**表2 主要病变候选的识别与量化结果**
- 正文落位：第 7.2 节
- 建议保留列：
  - 病例、侧别、病变编号；
  - 分支编号、分支层级；
  - 中心位置、病变长度；
  - 最小半径、参考半径、狭窄率、严重度。

### 表3
- 来源文件：`analysis_results/tables/local_best_views.csv`
- 建议表题：**表3 主要病变的局部最佳观察视角结果**
- 正文落位：第 7.3 节
- 建议保留列：
  - 病变编号、排名；
  - 方位角、俯仰角；
  - 综合评分、缩短评分、重叠评分。

### 表4
- 来源文件：`analysis_results/tables/global_view_plans.csv`
- 建议表题：**表4 病例级全局观察方案**
- 正文落位：第 7.4 节
- 建议正文中优先展示 `side = case` 的记录；
- 建议保留列：
  - 病例、方案类型（`2-view` 或 `3-view`）、排名；
  - 方位角、俯仰角；
  - 平均评分；
  - 覆盖病变集合。

### 表5
- 来源文件：`analysis_results/tables/sensitivity_analysis.csv`
- 建议表题：**表5 狭窄率阈值敏感性分析结果**
- 正文落位：第 8 节“稳定性与敏感性分析”
- 建议保留列：
  - 病例、侧别、阈值；
  - 候选病变数；
  - 形成病变的分支数；
  - 最大狭窄率。

### 表6
- 来源：由 `basic_structure_stats.csv`、`lesion_candidates.csv`、`local_best_views.csv` 和 `global_view_plans.csv` 汇总整理
- 建议表题：**表6 病例级综合分析结果汇总**
- 正文落位：第 7.5 节“不同病例与不同冠脉的对比分析”
- 建议保留列：
  - 病例、侧别；
  - 分支数；
  - 是否检出主要病变；
  - 最强病变狭窄率；
  - 最佳局部视角评分；
  - 推荐全局 `2-view` 方案。

## 三、附录建议结构

### 附录A 参数设置

建议在附录中集中列出核心参数，避免在正文中反复拆散：

- 最大连通域保留：启用
- 重采样目标间距：`0.4 mm`
- 半径平滑窗口：`5` 点
- 近端/远端参考窗口：`5 mm`
- 狭窄率阈值：`0.30`
- 病变最小长度：`2 mm`
- 病变最大长度：`20 mm`
- 病变合并间隔：`3 mm`
- 方位角搜索范围：`0° ~ 355°`，步长 `5°`
- 俯仰角搜索范围：`-60° ~ 60°`，步长 `5°`
- 局部视角评分：`0.6 × 缩短评分 + 0.4 × 重叠评分`
- 全局覆盖阈值：局部评分达到最优评分的 `70%`

### 附录B 结果文件映射

建议在附录中列出论文与结果文件的对应关系，便于复核：

- `basic_structure_stats.csv` → 表1
- `lesion_candidates.csv` → 表2
- `local_best_views.csv` → 表3
- `global_view_plans.csv` → 表4
- `sensitivity_analysis.csv` → 表5
- 汇总整理后的病例级对比表 → 表6
- `figure1_workflow.png` → 图1
- `figure2_mip_overview.png` → 图2
- `figure3_centerline_overview.png` → 图3
- `figure4_radius_curves.png` → 图4
- `figure5_local_view_comparison.png` → 图5
- `figure6_global_heatmap.png` → 图6

### 附录C 代码与实现说明

如需说明可复现性，可在附录中简要列出：

- 主分析脚本：`tools/run_coronary_analysis.py`
- 核心模块：
  - `coronary_analysis/geometry.py`
  - `coronary_analysis/centerline.py`
  - `coronary_analysis/lesion.py`
  - `coronary_analysis/view_selection.py`
  - `coronary_analysis/pipeline.py`
- 单元与集成测试：`tests/test_coronary_analysis.py`

## 四、Word 迁移建议

1. 先将 [赛题C_正文源稿.md](</D:/aaaaaaaaaaaaa/manuscript/赛题C_正文源稿.md>) 作为内容底稿导入或复制到 Word。
2. 按本清单中的图号、表号顺序插入图表，保持编号与正文引用一致。
3. 图题建议放在图下方，表题建议放在表上方，格式以竞赛模板为准。
4. 在最终排版阶段，不再修改模型口径、数值结果和图表顺序，只处理样式、字体、页边距、标题编号和公式编号。
