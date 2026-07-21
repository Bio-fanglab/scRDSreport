# scRDSreport

`scRDSreport` 把单细胞 RDS 整理成一个可浏览、可追溯、可下载的 Quarto 报告。同时可审计多物种资源系统,分析入口只有 `running()`；用户提供输入 RDS 和输出文件夹即可：

```r
library(scRDSreport)

running(
  input = "path/to/object.rds",
  output = "path/to/result"
)
```

报告固定写入 `output/report.html`。原始 RDS、分析对象、表达矩阵、细胞元数据、输入对象中已有的注释、结果表和图形都放在同一个 `output` 文件夹下。

请为 `output` 使用一个新的或空的专用结果目录。首次运行会写入所有权标记；以后只有带该标记（或可验证旧版 manifest）的 scRDSreport 目录才能用 `overwrite = TRUE` 更新。程序会拒绝 HOME、文件系统根目录、非空的普通项目目录，以及把输入 RDS 放在同一输出目录的 `downloads/`、`analysis/` 等受管理路径中的用法，避免误删输入或无关文件。

默认行为会先检查 RDS：raw/partial 对象使用 `profile = "full"` 补齐分析；已经含有降维和 cluster/注释的对象在 `analyze = "auto"` 下自动采用 `report_only`，只整理并完整导出现有结果。若确实希望对已分析对象追加高级模块，请显式传入 `report_config(profile = "full")`。

默认 `annotation_mode = "auto_if_missing"`：RDS 已有注释时原样保留；没有注释时，只在当前物种注册了可信参考时尝试参考注释。自动结果写入独立列，不覆盖原始元数据，也不会拿另一个物种的参考填充。样本设计、轨迹起点和 CNV 正常参考仍不会被擅自猜测。

完整档默认尽量为每个章节给出可解释结果，但不会为了“有图有表”制造统计结论：没有可靠生物学重复时，差异模块按注释或 cluster 输出有数量上限的 one-vs-rest 描述性效应量和基因排序，不生成 P 值/FDR；富集模块可据此输出基因集方向与效应摘要，不声称显著性；没有轨迹起点、CellChat 推断条件或 CNV 正常参考时，分别只输出无方向几何、分组/诊断信息或 readiness 表。

## 与 FASTQ/nf-core 流程的边界

原 `run_scrnaseq.sh` 包含 FASTQ 到表达矩阵的上游流程，本包从 RDS 开始：

```text
FASTQ + samplesheet
  └─ run_scrnaseq.sh / nf-core/scrnaseq / Cell Ranger → 单细胞 RDS
       └─ scRDSreport::running() → 模块化分析结果 + report.html
```

`scRDSreport` 不执行 FASTQ 比对、定量、Cell Ranger 或 nf-core/Nextflow。FASTQ 用户应先用原脚本或其他流程生成 Seurat、SingleCellExperiment、表达矩阵或含 counts 的 list RDS，再把该 RDS 交给 `running()`。这样上游流程可以独立升级，R 包只负责从 RDS 开始的检查、补分析、报告和下载。

## 安装

需要 R >= 4.2。安装器会自动处理 R 包，但不会静默安装操作系统软件；开始前请按实际模块准备以下项目：

- 生成 HTML 需要按官方说明安装 [Quarto CLI >= 1.3](https://quarto.org/docs/get-started/)；安装名为 `quarto` 的 R 包并不等于安装了 Quarto CLI。
- `infercnv` 模块需要 JAGS 4.x。
- 当前 Monocle 3/BPCells 自动安装路径需要可由 `h5cc` 或 `pkg-config` 检出的 HDF5 开发库；安装最新版还需要 R >= 4.4.1 和 Bioconductor >= 3.21。较旧 R 上已能加载的兼容 Monocle 3 仍可继续使用。

### 方式 A：新环境一键安装（推荐）

仓库中的 bootstrap 把“安装 scRDSreport 及硬依赖 → 按 profile/物种安装分析依赖 → 实际加载并严格验收”合成一次调用。`species = "auto"` 在 full/core 档必须同时给出输入 RDS。建议先打开脚本链接检查内容，再运行：

```r
source("https://raw.githubusercontent.com/Bio-fanglab/scRDSreport/main/inst/install_scRDSreport.R")

install_scRDSreport(
  input = "path/to/object.rds",
  profile = "full",
  species = "auto"
)
```

不需要 HTML 时可设 `render = FALSE`。如果只想安装核心或报告依赖，可把 `profile` 改为 `"core"` 或 `"report_only"`。

### 方式 B：分步安装

先安装包本身。这里把 Bioconductor 仓库显式交给 `remotes`，因此 SCP 的 CRAN/Bioconductor 硬依赖会一起解析和安装：

```r
install.packages(c("BiocManager", "remotes"))

bio_repos <- BiocManager::repositories()
remotes::install_github(
  "Bio-fanglab/scRDSreport",
  dependencies = NA,
  upgrade = "never",
  repos = bio_repos
)
```

再针对本次 RDS 安装完整分析依赖并做真实加载验收。`species = "auto"` 会从输入 RDS 的稳定 feature ID 选择对应的 OrgDb/TxDb；不能可靠判断时会停止并要求显式物种，不会猜成人或小鼠：

```r
input_rds <- "path/to/object.rds"

scRDSreport::install_dependencies(
  profile = "full",
  input = input_rds,
  species = "auto"
)

scRDSreport::check_dependencies(
  profile = "full",
  input = input_rds,
  species = "auto"
)
```

已知是小鼠时不必读取 RDS 来判断物种：

```r
scRDSreport::install_dependencies(profile = "full", species = "mouse")
```

默认只装本次物种；确实需要在同一环境处理全部九个内置物种时，显式使用 `species = "all"`。这会安装九套 OrgDb，并可能占用较多磁盘空间。安装前可先查看计划，不会写入 library：

```r
plan <- scRDSreport::install_dependencies(
  profile = "full",
  species = "all",
  dry_run = TRUE
)
plan[plan$action != "none", c("component", "source", "action")]
```

`dependency_status()` 和 `check_dependencies()` 使用 namespace 实际加载结果，不会把“目录存在但因传递依赖缺失而无法加载”的包误判为可用。安装器会自动处理 CRAN、Bioconductor，以及 SCP、CellChat、Monocle 3 的 GitHub 来源。系统项缺失时，它会先完成其余 R 包，再汇总给出未就绪清单。包安装成功也不代表所有生物学前提都成立：轨迹起点、CNV 正常参考、可信注释和实验设计仍需用户提供。开发和本地检查再安装 `testthat`/`devtools`，普通用户不需要。

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

## 注释：已有结果优先，缺失时才用参考

默认的 `annotation_mode = "auto_if_missing"` 会先寻找 RDS 元数据中的已有注释。找到时原样保存和导出，不重新注释；没找到时，人和小鼠可在参考依赖可用时分别尝试注册的 `celldex` 参考，其他内置物种保持 cluster-only 并显示 `needs_input`，直到用户提供同物种参考。

四种模式的区别如下：

| 模式 | 行为 |
|---|---|
| `auto_if_missing` | 默认；已有注释就保留，没有注释才尝试物种匹配参考 |
| `preserve` | 只使用并导出 RDS 已有注释；缺失时不新增 |
| `auto` | 显式请求参考注释；结果写入新列并与原注释并列 |
| `manual` | 使用用户提供的 mapping/marker；结果写入新列 |

已知列名不标准时可显式指定并只保留它：

```r
cfg <- report_config(annotation_mode = "preserve")

running(
  "object.rds",
  "result",
  annotation_col = "celltype_manual",
  config = cfg
)
```

任何模式都不会覆盖已有元数据列，也不会使用包作者自定义的 cell-type 标签填空。参考、label 字段、feature overlap、置信标签数和输出列都会写入产物或 manifest。cluster-level SingleR 缺少可选的 `scrapper` 时，小于配置上限的对象可退回 cell-level SingleR，并在 provenance 中记录该兼容路径；不确定的 `pruned.labels` 仍保持缺失。

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

默认描述性兜底还会优先使用对象中可信的注释列，其次使用 cluster，生成 one-vs-rest marker 排序和效应量图；默认最多展示 20 个组，样本层比较最多 6 个，避免组数较多时产生平方级结果。需要改变这些限制时，可在 `module_options$differential` 中设置 `fallback_grouping`、`run_group_markers`、`max_marker_groups` 和 `max_contrasts`。

## 物种和资源

`species = "auto"` 会优先根据稳定 feature ID 前缀识别物种；只有大写 gene symbol 时无法区分人、猪、牛、狗、鸡和猕猴，因此保持 `unknown`，不会武断地套用人或小鼠资源。显式指定物种时，该选择控制资源，自动检测结果和冲突提示仍会写入 manifest。

```r
supported_species()              # 九个内置物种及当前依赖可用性
species_resources("mouse")      # 单个物种的完整资源对象
species_resources("unknown")    # supported = FALSE，不套用人/鼠资源
```

| 支持层级 | 物种 | 内置资源和边界 |
|---|---|---|
| `full` | human、mouse | OrgDb、KEGG、MSigDB、可信自动注释入口、CellChat 和 TxDb 入口；模块仍逐项检查依赖与数据前提 |
| `core` | rat、zebrafish、pig、cattle、chicken、dog、macaque | 物种匹配的 ID/QC/OrgDb/KEGG/MSigDB 元数据；没有内置可信注释、CellChat DB 或 TxDb，相关模块需要用户资源 |
| custom/unregistered | 其他物种 | 通用读取、QC、降维、聚类和下载；物种特异模块需要 `resource_overrides` |

人使用 human MSigDB（`db_species = "HS"`、Hallmark collection `H`），小鼠使用 Mouse MSigDB（`db_species = "MM"`、Hallmark collection `MH`）。其余七个非人内置物种使用 `msigdbr` 将 human MSigDB 正交投影到目标物种；报告会明确记录 `db_species = "HS"`、目标物种、数据库版本和可用的 ortholog evidence，不能把它表述成目标物种原生数据库。GO/KEGG 与 TF catalog 使用相应 OrgDb；默认 TF catalog 来自 `GOALL:GO:0003700`。非人细胞周期基因由 `babelgene` 把 Seurat 的 human S/G2M 集合映射到目标物种，并导出实际映射表。

### 指定小鼠并继续运行完整分析

对于 raw/partial RDS，以下调用会按需用 SCP 补齐基础分析，并选择小鼠资源。对于已经有降维和 cluster/注释的 RDS，显式的 `profile = "full"` 很重要；否则默认自动判断可能只生成 `report_only`。

```r
cfg_mouse <- report_config(
  profile = "full",
  annotation_mode = "auto_if_missing"
)

running(
  input = "mouse_object.rds",
  output = "mouse_result",
  species = "mouse",
  config = cfg_mouse
)
```

`species = "mouse"` 会选择 `org.Mm.eg.db`、`mmu` KEGG、小鼠 QC 基因模式、Mouse MSigDB 的 `MH` Hallmark collection、`CellChatDB.mouse`、mm10 TxDb 入口，以及经 `babelgene` 映射的小鼠 cell-cycle 基因。RDS 没有注释时，默认参考为通用小鼠 `celldex::MouseRNAseqData`；免疫数据可显式改用 `celldex::ImmGenData`：

```r
cfg_mouse_immune <- report_config(
  profile = "full",
  annotation_mode = "auto_if_missing",
  module_options = list(
    celltype = list(reference = "celldex::ImmGenData")
  )
)
```

若 RDS 已有注释，`auto_if_missing` 不运行 SingleR。若确实希望生成一个独立的参考注释列用于比较，可显式使用 `annotation_mode = "auto"`。无论物种如何，程序仍不会自动选择轨迹起点或 CNV 正常参考。

其他内置物种没有注册的可信自动注释参考。用户可以提供同物种参考：

```r
cfg_rat <- report_config(
  profile = "full",
  annotation_mode = "auto_if_missing",
  module_options = list(
    celltype = list(
      reference = rat_reference,
      reference_labels = rat_reference_labels
    )
  )
)
running("rat.rds", "rat_result", species = "rat", config = cfg_rat)
```

未注册物种或缺失资源可通过 `resource_overrides` 提供与 feature ID 和 genome assembly 匹配的 OrgDb、KEGG code、CellChat 数据库、TxDb/GTF、基因集或 marker。提供少量字段不等于所有模块都能运行；每个模块仍会独立报告前提和状态。

## 轨迹起点和 CNV 参考

轨迹起点和正常 CNV 参考都属于生物学假设，程序不会根据图形自动猜测：

```r
cfg <- report_config(
  profile = "full",
  trajectory_root = "Y_1",                 # Monocle principal node
  cnv_reference = c("T_cell", "B_cell"),  # 注释列中的正常参考组
  module_options = list(
    cnv = list(object_genome_assembly = "GRCm38/mm10")
  )
)

running("mouse_object.rds", "result", species = "mouse", config = cfg)
```

未提供 `trajectory_root` 时，轨迹模块最多输出无方向的轨迹几何和候选起点，不生成定向 pseudotime 或动态基因结论；未提供 `cnv_reference` 时，`cnv` 章节显示 `cnv_reference_required`，并可附带输入就绪情况，但绝不生成 CNV 信号。即使选择了物种和正常参考，内置 TxDb 也只有在 `object_genome_assembly` 与注册版本匹配，或用户核验后显式设置 `genome_assembly_confirmed = TRUE` 时才会用于推断；也可以直接提供与对象匹配的 `gene_order`、`gtf` 或 `txdb`。程序不会悄悄选择起点、正常参考组或基因组版本。

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

页面样式、图形、DT 表格与脚本内置在 `report.html`，因此报告离线可打开。左侧目录默认完整展开；结果表继续使用 DT 的搜索、复制、CSV 和 Excel 控件。完整文件索引的说明列采用固定高度的独立滚动区域，避免长说明把整行撑得过高。下载文件有三种策略：

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
