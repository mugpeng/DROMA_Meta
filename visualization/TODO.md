# DROMA.Meta Visualization — TODO

> Status: **Phase 1-4 implemented** — Phase 5 (docs/tests) pending
> Updated: 2026-04-16
> Decisions: subdirectory (future merge into DROMA.Meta) | Tableau 10 colors | getwd() paths | all phases

---

## Phase 1: 包基础设施 ✅

- [x] `visualization/DESCRIPTION` — 声明 Package, Imports (ggplot2, ggrepel, patchwork, ComplexHeatmap, circlize, grid, ggsci, stats), Depends (DROMA.Meta)
- [x] `visualization/NAMESPACE` — 手写，20 个 export
- [x] `getVisWorkflowDefaults()` 返回 paths，与父包 `getMetaWorkflowDefaults()` 风格一致
- [x] workflow 脚本统一 pattern: `library(DROMA.Meta)` → `getVisWorkflowDefaults()` → `dir.create()`
- [x] `sanitizeName()` 通过 `DROMA.Meta::sanitizeName()` 或 `library(DROMA.Meta)` 引入
- [x] `ggrepel` 所有调用均有 `requireNamespace` guard

## Phase 2: 代码质量 ✅

- [x] `plotPipelineFunnel()` 修复: 4 阶段 (All genes → Intersection → TCGA-AD → Clinical)，不再重复 key
- [x] `collectAllBatchResults()` 改用 `basename(dirname(f))` 提取 drug/tumor_type
- [x] 每个 plot 函数有 `required_cols` check + `stop(call.=FALSE)`
- [x] 空数据用 `warning()` + 返回空 data.table（收集函数），plot 函数用 `stop()`（调用方负责 tryCatch）
- [x] roxygen: 所有公开函数完整 `@title`, `@description`, `@param`, `@return`, `@export`
- [x] 增加 `collectAllAdStats()`, `collectAllSigFeatures()` 收集函数

## Phase 3: 发表级图表 ✅

- [x] `themeMetaPaper()` — base_size=11, sans, black axis text, no minor grid, grey92 major grid, 10pt margins
- [x] 配色: Tableau 10 colorblind-safe (Up=#D62728, Down=#1F77B4, NS=#C7C7C7)
- [x] PDF 默认 7.2×5.6 inches (Nature double-column)
- [x] Volcano: 背景着色象限 + direction legend + 方向标注 + ggrepel label
- [x] Pipeline funnel: 阶段色条 + 流失百分比标注
- [x] Heatmap: diverging RdBu 色阶 + direction top annotation + cell text option + `saveComplexHeatmapPdf()`
- [x] Concordance: LM 95% CI + zero reference lines + correlation annotation + density2d (n>50) + gene labels (n<=50)
- [x] Upset: left annotation (set size bar) + stage colors + empty set guard
- [x] Summary: direction bar (带百分比) + dumbbell attrition chart

## Phase 4: 新增图表 ✅

- [x] `FuncVisTcgaAD.R` — AD p-value 直方图 + AD statistic density overlay
- [x] `FuncVisForest.R` — stage-wise forest plot (cell → pdcpdx → clinical per gene)
- [x] `19-Figure1_Panel.R` — patchwork 4-panel (funnel + volcano + concordance + direction)

## Phase 5: 文档与测试 (pending)

- [ ] 每个 workflow 脚本顶部注释块补充 input/output/depends 说明
- [ ] 空 data 边界测试: 确认所有 plot 函数对 `data.table(name=character(0))` 报错优雅
- [ ] 合并入 DROMA.Meta 主包时的 NAMESPACE 合并计划

---

## 文件清单

```
visualization/
  DESCRIPTION                     # 包元数据
  NAMESPACE                       # 20 个 export
  TODO.md                         # 本文件
  R/
    FuncVisHelper.R               # 主题、颜色、路径、数据收集、PDF导出
    FuncVisPipeline.R             # Pipeline funnel + pair count bar
    FuncVisVolcano.R              # Volcano plot (单张 + 批量)
    FuncVisHeatmap.R              # ComplexHeatmap 效应量热图
    FuncVisConcordance.R          # Preclinical vs Clinical 散点
    FuncVisUpset.R                # Upset 多层交集
    FuncVisSummary.R              # 方向分布 + dumbbell attrition
    FuncVisTcgaAD.R               # TCGA AD p-value 直方图 + density
    FuncVisForest.R               # Stage-wise forest plot
  workflow/
    10-Collect_Results.R          # 汇总 → 7 张 collected CSV
    11-Pipeline_Overview.R        # funnel + pair count PDF
    12-Volcano_Plots.R            # 批量 volcano PDF
    13-Heatmap_Biomarkers.R       # heatmap PDF
    14-Concordance_Scatter.R      # concordance PDF
    15-Upset_Intersection.R       # 批量 upset PDF
    16-Summary_Plots.R            # direction + dumbbell PDF
    17-TCGA_AD_Plots.R            # AD histogram + density PDF
    18-Forest_Plots.R             # forest PDF
    19-Figure1_Panel.R            # 4-panel composite PDF
```

## 运行方式

```r
# 在 DROMA_Meta 项目根目录 (与 Meta_Example/ 同级)
setwd("/Users/peng/Desktop/Project/DROMA/Meta_project/DROMA_Meta")

# Step 1: 收集数据 (必须先运行)
source("visualization/workflow/10-Collect_Results.R")

# Step 2-9: 生成图表 (可并行，顺序无关)
source("visualization/workflow/11-Pipeline_Overview.R")
source("visualization/workflow/12-Volcano_Plots.R")
source("visualization/workflow/13-Heatmap_Biomarkers.R")
source("visualization/workflow/14-Concordance_Scatter.R")
source("visualization/workflow/15-Upset_Intersection.R")
source("visualization/workflow/16-Summary_Plots.R")
source("visualization/workflow/17-TCGA_AD_Plots.R")
source("visualization/workflow/18-Forest_Plots.R")
source("visualization/workflow/19-Figure1_Panel.R")
```

输出目录: `Meta_Example/Output/visualization/*.pdf`
