# scRDSreport

`scRDSreport` 把单细胞 RDS 整理成一个可浏览、可追溯、可下载的 Quarto 报告。入口函数只有 `running()`；用户提供输入 RDS 和输出文件夹即可：

```r
library(scRDSreport)

running(
  input = "path/to/object.rds",
  output = "path/to/result"
)
```

报告固定写入 `output/report.html`。原始 RDS、分析对象、表达矩阵、细胞元数据、输入对象中已有的注释、结果表和图形都放在同一个 `output` 文件夹下。

请为 `output` 使用一个新的或空的专用结果目录。首次运行会写入所有权标记；以后只有带该标记（或可验证旧版 manifest）的 scRDSreport 目录才能用 `overwrite = TRUE` 更新。程序会拒绝 HOME、文件系统根目录、非空的普通项目目录，以及把输入 RDS 放在同一输出目录的 `downloads/`、`analysis/` 等受管理路径中的用法，避免误删输入或无关文件。

默认行为会先检查 RDS：raw/partial 对象使用 `profile = "full"` 补齐分析；已经含有降维和 cluster/注释的对象在 `analyze = "auto"` 下自动采用 `report_only`，只整理并完整导出现有结果。若确实希望对已分析对象追加高级模块，请显式传入 `report_config(profile = "full")`。默认注释模式始终是 `preserve`，不会擅自改写细胞注释、样本设计、轨迹起点或 CNV 参考。

## 安装

需要 R >= 4.1 和 [Quarto CLI >= 1.3](https://quarto.org/docs/get-started/)。以下是生成报告、运行 Seurat/SCP 核心流程所需的安装层：

```r
install.packages(c("BiocManager", "remotes", "quarto"))

BiocManager::install("BiocParallel", ask = FALSE, update = FALSE)
remotes::install_github("zhanghao-njmu/SCP", dependencies = NA)
remotes::install_github(
  "Bio-fanglab/scRDSreport",
  dependencies = NA,
  upgrade = "never"
)
```

`full` 中的高级章节使用可选依赖。希望尽可能运行全部模块时，再安装下面这一层：

```r
install.packages(c(
  "circlize", "data.table", "dplyr", "future", "ggrepel",
  "htmlwidgets", "igraph", "plotly", "scales", "stringr", "tidyr"
))

BiocManager::install(c(
  "AnnotationDbi", "celldex", "clusterProfiler", "ComplexHeatmap",
  "edgeR", "enrichplot", "GenomicFeatures", "GSVA", "infercnv", "rtracklayer",
  "monocle3", "org.Hs.eg.db", "org.Mm.eg.db",
  "SingleCellExperiment", "SingleR", "SummarizedExperiment",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "TxDb.Mmusculus.UCSC.mm10.knownGene"
), ask = FALSE, update = FALSE)

remotes::install_github("jinworks/CellChat", dependencies = TRUE)
remotes::install_github("satijalab/seurat-wrappers", dependencies = TRUE)
```

`KEGG.db` 仅用于仍依赖它的旧版 Bioconductor 兼容路径，可在当前 Bioconductor 仓库仍提供时另行安装。开发和本地检查再安装 `devtools` 与 `testthat`，普通用户不需要这两个包。

不同 R/Bioconductor 版本可提供的高级包可能不同。缺少某个可选包只会使对应模块被跳过，不影响核心报告和数据导出；实际使用的包版本、警告和跳过原因会写入 manifest。

## 三种运行档位

```r
# 完整档：请求 12 个模块，默认值
cfg_full <- report_config(profile = "full")
running("object.rds", "result_full", config = cfg_full)

# 核心档：QC、降维、聚类、注释/组成和下载
cfg_core <- report_config(profile = "core")
running("object.rds", "result_core", config = cfg_core)

# 只报告：不补跑 SCP 或高级分析，只整理输入对象中已有内容
cfg_report <- report_config(profile = "report_only")
running("object.rds", "result_report", config = cfg_report)
```

完整计划固定保留以下 12 个模块 ID：

```text
qc, reduction, cluster, celltype, differential, enrichment,
pseudotime, communication, cell_cycle, tf, cnv, downloads
```

可以用 `modules` 只选择一部分；`downloads` 始终保留：

```r
cfg <- report_config(
  profile = "full",
  modules = c("qc", "reduction", "cluster", "celltype")
)
running("object.rds", "result_selected", config = cfg)
```

每个模块都会记录是否请求、是否具备前提、最终状态、原因、分析引擎与版本、参数、随机种子、耗时、警告、错误和产物 ID。完整说明见 [完整分析与资源配置](docs/full-analysis.md)。

## 注释：默认只保留 RDS 原始内容

默认的 `annotation_mode = "preserve"` 会寻找并原样导出 RDS 元数据中的已有注释，不创建新的 cell type，也不覆盖已有列。已知列名不标准时可显式指定：

```r
cfg <- report_config(annotation_mode = "preserve")

running(
  "object.rds",
  "result",
  annotation_col = "celltype_manual",
  config = cfg
)
```

只有用户明确选择 `annotation_mode = "auto"` 或 `"manual"` 时，程序才会尝试产生新的注释列；对应参考数据、物种映射或 marker 不完整时，注释模块会说明原因并跳过，不会用包内自定义标签填充。

## 样本分组和重复

程序会从样本名提出保守的 `group`/`replicate` 建议。自动结果中的 `confidence` 和 `needs_review` 用于提醒人工复核，`grouping_rule` 是短规则代码，不是结果正文。正式组间比较建议始终传入真实实验设计：

```r
sample_map <- data.frame(
  sample_id = c("Ctrl_1", "Ctrl_2", "Drug_1", "Drug_2"),
  group = c("Control", "Control", "Drug", "Drug"),
  replicate = c("1", "2", "1", "2")
)

running(
  "object.rds",
  "result",
  sample_map = sample_map
)
```

每个比较组至少有两个独立生物学样本时，`differential = "auto"` 才能进行重复感知的 pseudobulk 推断。没有生物学重复时，细胞不能冒充重复：报告最多给出描述性表达量和效应大小，不制造 P 值。显式选择 cell-level Wilcoxon 的结果也只能作为探索性结果解释。

## 物种和资源

`species = "auto"` 会根据 feature ID/基因名尝试识别物种；证据含糊时保持 `unknown`，不会回退到小鼠资源。当前内置的完整资源映射只有：

- `human`：人，包含人源 OrgDb、KEGG、CellChat、hg38 TxDb 等资源入口；
- `mouse`：小鼠，包含鼠源 OrgDb、KEGG、CellChat、mm10 TxDb 等资源入口。

```r
species_resources("human")
species_resources("mouse")
species_resources("unknown")  # supported = FALSE，不套用人/鼠资源

running("object.rds", "result_human", species = "human")
running("object.rds", "result_mouse", species = "mouse")
```

其他物种可以读取和导出，但物种特异模块需要用户通过 `resource_overrides` 提供与其基因 ID 和基因组版本匹配的 OrgDb、KEGG 代码、CellChat 数据库、TxDb/GTF、基因集或 marker。提供少量字段不等于该物种的所有模块都可运行；每个模块仍会单独检查前提。

## 轨迹起点和 CNV 参考

轨迹起点和正常 CNV 参考都属于生物学假设，程序不会根据图形自动猜测：

```r
cfg <- report_config(
  profile = "full",
  trajectory_root = "Y_1",                 # Monocle principal node
  cnv_reference = c("T_cell", "B_cell")  # 注释列中的正常参考组
)

running("object.rds", "result", config = cfg)
```

未提供 `trajectory_root` 时，轨迹模块最多输出无方向的轨迹几何和候选起点，不生成定向 pseudotime 或动态基因结论；未提供 `cnv_reference` 时，`cnv` 章节显示 `cnv_reference_required`。程序不会悄悄选择起点或正常参考组。

## 超大 raw 对象

默认保留全部分析细胞和 feature。对超大 raw 10x 对象，建议限制计算和图形规模，同时保留原始输入下载：

```r
cfg <- report_config(
  profile = "full",
  limits = list(
    analysis_max_cells = 20000,
    analysis_max_features = 15000,
    plot_max_cells = 50000,
    marker_max_cells_per_ident = 1500,
    workers = 4,
    embed_max_mb = 50
  )
)

running(
  "large_raw.rds",
  "large_result",
  config = cfg,
  filter_raw_barcodes = "auto",
  filter_low_expression_features = "auto",
  scp_args = list(nHVF = 2000, linear_reduction_dims = 30),
  embed_downloads = "auto"
)
```

有限的 `analysis_max_cells`/`analysis_max_features` 会建立确定性分析子集；原始 RDS 和原始 counts 仍单独导出，manifest 会记录抽样前后规模。报告不会把子集分析标成“全数据分析”。大型 RDS 不建议使用 `embed_downloads = "always"`，因为 Base64 内置会明显增加 HTML 体积和浏览器内存占用。

## 单文件 HTML 和下载

页面样式、图形、DT 表格与脚本内置在 `report.html`，因此报告离线可打开。下载文件有三种策略：

```r
# 默认：在预算内优先把模块结果、注释和说明表内置到 HTML
running("object.rds", "result", embed_downloads = "auto", embed_max_mb = 50)

# 强制所有下载文件内置；只适合确认体积可控的小对象
running("object.rds", "result_one_file", embed_downloads = "always")

# HTML 只保留指向 output 文件夹的下载入口
running("object.rds", "result_folder", embed_downloads = "never")
```

`auto` 超出预算的 RDS/大型矩阵仍放在 output 文件夹。因此需要向别人完整交付大型项目时，推荐压缩整个 output 文件夹；如果只发送 `report.html`，只能下载其中标为“HTML 内置”的文件。

表达矩阵使用单细胞通用方向：行是 feature/基因，列是细胞；Matrix Market 文件配套的 features 和 barcodes 文件严格保持相同顺序。注释表是一行一个细胞，第一列为 cell barcode，后续列来自 RDS 中已有注释或用户明确请求的注释结果。每张报告表下方都会说明行、列、字段和统计单位。

## 输出目录

```text
output/
├── report.html
├── analysis/
│   ├── qc/
│   ├── reduction/
│   ├── cluster/
│   ├── celltype/
│   ├── differential/
│   └── ...                         # 仅有产物的高级模块
├── downloads/
│   ├── original_*.rds
│   └── analysis_object.rds
├── matrices/
│   ├── original_*_counts.mtx.gz
│   ├── original_*_counts_features.tsv.gz
│   ├── original_*_counts_barcodes.tsv.gz
│   └── analysis_*_{counts,data}.*
├── tables/
│   ├── original_cell_metadata.csv.gz
│   ├── original_cell_annotations.csv.gz   # 输入 RDS 有注释时
│   ├── analysis_cell_metadata.csv.gz
│   ├── analysis_cell_annotations.csv.gz   # 分析对象有注释时
│   ├── sample_design.csv.gz
│   └── embedding_*.csv.gz
├── figures/
└── .report/
    ├── manifest.rds
    └── manifest.json
```

没有运行前提的模块可能没有对应产物文件夹，但仍会在报告的模块状态和 manifest 中保留跳过原因。

## 开发检查

```r
devtools::test()
devtools::check()
```

问题反馈：[Bio-fanglab/scRDSreport issues](https://github.com/Bio-fanglab/scRDSreport/issues)
