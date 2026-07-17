# 完整分析、统计边界与资源配置

本文说明 `scRDSreport` 0.2 的完整分析计划。最重要的原则是：每个章节先检查数据和生物学前提，再决定运行、部分运行或跳过。没有前提时，报告保留章节状态和原因，不用虚构注释、重复、P 值、轨迹起点或 CNV 参考来填满报告。

## 1. 一条命令和三个 profile

最简调用会自动检查输入状态：raw/partial 对象进入完整分析，已有降维和 cluster/注释的对象只导出现有分析：

```r
library(scRDSreport)
running("object.rds", "result")
```

对 raw/partial 对象，等价的显式完整分析写法是：

```r
cfg <- report_config(
  profile = "full",
  annotation_mode = "preserve",
  differential = "auto"
)

running(
  input = "object.rds",
  output = "result",
  species = "auto",
  config = cfg
)
```

三个 profile 的差异如下：

| profile | 请求的模块 | 适用情况 |
|---|---|---|
| `full` | 全部 12 个模块 | 希望尽可能完成标准分析，同时接受不具备前提的章节被明确跳过 |
| `core` | `qc`、`reduction`、`cluster`、`celltype`、`downloads` | 只需要基础分析、原注释/组成和数据导出 |
| `report_only` | `downloads` | 输入对象已有分析，只整理和展示，不补跑 SCP 或高级分析 |

`profile` 控制报告模块；`analyze = "auto"/"always"/"never"` 控制是否用 SCP 补全 raw 或 partial 对象。没有显式 `config` 时，`analyze = "auto"` 会让已分析 RDS 自动进入 `report_only`。显式传入 `profile = "full"` 才会在已有分析之上继续尝试高级模块。

## 2. 12 个模块

模块 ID 是固定的。模块未请求、缺少输入、缺少可选包或运行失败时，其他模块仍可继续。

| ID | 主要内容 | 必要前提和边界 |
|---|---|---|
| `qc` | 稀疏 counts 上的细胞/样本 QC、阈值、通过状态、QC 图 | 需要可读取的表达层；默认只记录指标，只有显式配置才按高级模块阈值修改分析对象 |
| `reduction` | HVG/PCA 方差、PCA、UMAP、t-SNE 等已有降维的坐标与图 | raw/partial 对象可先由 SCP 补全；已有完整降维时优先保留 |
| `cluster` | cluster 元数据、数量/比例和降维图 | 需要降维/邻居图或输入对象已有 cluster |
| `celltype` | 原注释、细胞组成、样本组成、注释降维图 | `preserve` 需要 RDS 已有注释；`auto`/`manual` 必须由用户显式启用并具备参考资源 |
| `differential` | 按 cell type/cluster 聚合的 pseudobulk、比较表、效应大小 | 需要明确的 sample、group 和分层字段；正式推断要求每组至少两个独立生物学样本 |
| `enrichment` | GO、KEGG、GSEA/GSVA 及其表图 | 需要可用差异排序和匹配物种的 OrgDb、KEGG/基因集资源 |
| `pseudotime` | 轨迹几何、候选起点，以及有显式 root 时的 pseudotime、动态基因和富集 | 需要足够细胞和可用降维；定向结论必须由用户提供 `trajectory_root` |
| `communication` | CellChat 网络、强度、配体受体/通路、发送者-接收者 | 需要可信细胞注释和匹配物种的 CellChat 数据库 |
| `cell_cycle` | S/G2M 得分、phase、组成和降维图 | 需要与物种及 feature 命名匹配的 cell-cycle 基因 |
| `tf` | TF ID 映射、平均表达、热图和摘要 | 需要匹配物种的 TF catalog/OrgDb；表达结果不是调控活性或因果推断 |
| `cnv` | inferCNV 输入、参考、信号和染色体摘要/热图 | 必须有用户指定的正常参考组和匹配 genome assembly 的 TxDb/GTF |
| `downloads` | 原始/分析 RDS、矩阵、元数据、产物、manifest 和 session info | 始终保留 |

例如只运行基础模块和差异分析：

```r
cfg <- report_config(
  profile = "full",
  modules = c(
    "qc", "reduction", "cluster", "celltype", "differential"
  )
)
```

也可以在某个 profile 上增减模块：

```r
cfg <- report_config(
  profile = "core",
  modules = list(
    include = "differential",
    exclude = "celltype"
  )
)
```

`downloads` 即使被排除也会重新加入，因为一次 `running()` 必须保留结果和可追溯信息。

## 3. 状态、原因和失败隔离

报告和 manifest 为每个模块保存下列信息：

- `requested`：profile/模块选择是否请求该模块；
- `eligible`：输入、物种资源和生物学前提是否满足；
- `status`：例如 `planned`、`completed`、`partial`、`skipped`、`needs_input`、`failed` 或 `not_requested`；
- `reason`/`reason_code`：适合程序读取的短原因代码；
- `message`：给用户阅读的说明；
- `engine` 和 `version`：实际分析引擎及版本；
- `parameters`、`seed`、开始/结束时间和耗时；
- 捕获的 `warnings`、`error` 和 `artifact_ids`。

`reason` 不是长篇方法说明。例如 `trajectory_root_required` 表示缺少轨迹起点，`cnv_reference_required` 表示缺少 CNV 参考。完整解释放在 `message` 和章节说明中，表格中的原因列保持为短代码。

以下情况是正常的可审计结果，不应视为程序偷偷漏跑：

- 未选择模块：`not_requested`；
- 没有 RDS 原始注释且使用默认 `preserve`：`celltype` 跳过或部分完成，并说明 `annotation_missing`；
- 物种无法安全识别：相应物种特异模块跳过；
- 没有轨迹起点：`pseudotime` 只输出无方向轨迹几何和候选起点，不生成定向 pseudotime；
- 没有正常参考组：`cnv` 显示 `cnv_reference_required`；
- 某个可选依赖未安装：只跳过对应模块；
- 模块内部报错：该模块记为 `failed`，其余安全模块继续。

## 4. annotation_mode

### 4.1 preserve：默认且推荐

```r
cfg <- report_config(annotation_mode = "preserve")
running(
  "object.rds",
  "result",
  annotation_col = "celltype_final",
  config = cfg
)
```

行为如下：

1. 在任何过滤或分析子集之前读取输入 RDS 的元数据；
2. 原样保存和导出已有 annotation 列；
3. 不用 cluster 名替换 cell type；
4. 不覆盖已有注释，也不增加包作者自定义标签；
5. 没找到注释时明确说明“无已有注释”。

### 4.2 auto：只有显式请求才运行

```r
cfg <- report_config(annotation_mode = "auto")
running("object.rds", "result_auto", species = "human", config = cfg)
```

自动注释需要 `SingleR`、`celldex`/用户参考对象、匹配的基因 ID 映射和足够 feature overlap。任何一项不满足时会显示跳过或需要输入。自动结果应写入独立列并与原注释并列展示，不能覆盖原列。

### 4.3 manual：用户提供 marker

```r
cfg <- report_config(
  annotation_mode = "manual",
  resource_overrides = list(
    manual_markers = list(
      T_cell = c("CD3D", "CD3E"),
      B_cell = c("MS4A1", "CD79A")
    )
  )
)
```

marker 必须来自本实验和对应物种。只设置 `manual` 而不提供 marker 时，模块显示 `manual_markers_required`，不会使用隐藏的默认 marker。

## 5. 样本设计、重复和统计限制

### 5.1 自动分组只是建议

程序从每个细胞的 sample ID 提取样本级设计。优先级是：

1. 用户提供的 `sample_map`；
2. 明确的 `rep1`、`replicate2`、`r3` 等末尾 token；
3. donor/subject/patient/mouse token；
4. 同一基名重复出现的数字后缀；
5. 无可靠模式时，每个样本保持自己的 group，replicate 留空。

自动设计中的列含义：

| 列 | 含义 |
|---|---|
| `sample_id` | 输入对象中的样本标识 |
| `group` | 建议的实验组；低置信度时可能等于 sample ID |
| `replicate` | 建议的重复编号；无法确认时为缺失 |
| `confidence` | `high`、`medium`、`low` 或用户映射的 `user` |
| `needs_review` | 是否必须人工复核 |
| `grouping_rule` | 产生该建议的短规则代码 |
| `n_cells_post_qc` | 该样本通过细胞级 QC 后的细胞数 |
| `n_cells_analysis` | 该样本实际进入分析对象（包括确定性子集）的细胞数 |

样本名中的数字可能代表时间、剂量、批次或受试者，程序不能从字符串知道真实实验设计。正式分析请传入：

```r
sample_map <- data.frame(
  sample_id = c("C01", "C02", "T01", "T02"),
  group = c("Control", "Control", "Treatment", "Treatment"),
  replicate = c("1", "2", "1", "2")
)

running("object.rds", "result", sample_map = sample_map)
```

`sample_map` 必须覆盖对象中的所有 sample，每个 `sample_id` 只能出现一次，每个 sample 只能映射到一个 group。

### 5.2 什么是可用于推断的重复

独立生物学样本是统计单位。多个细胞、同一样本的多个 10x lane 或重复测序都不能自动当作独立生物学重复。

`differential = "auto"` 的原则是：

- 每个比较组至少两个独立 sample：可以进行 replicate-aware pseudobulk（例如 edgeR quasi-likelihood）；
- 任一组少于两个独立 sample：只输出 pseudobulk CPM、log2 fold change 等描述性结果，不生成正式 P 值；
- sample/group 含义不明确：模块显示 `needs_input`，等待用户提供设计；
- 用户显式选择 `wilcox`：可以给出 cell-level 探索结果，但细胞不是生物学重复，不能把其 P 值解释成组间生物学推断。

同样地，样本细胞组成图默认是描述性图，不等于有重复的组成差异检验。

## 6. 物种资源

### 6.1 自动识别

```r
running("object.rds", "result", species = "auto")
```

程序优先使用稳定的 feature ID 前缀，再参考基因符号。符号大小写无法区分小鼠、大鼠和若干其他物种时，结果保持 `unknown`。`unknown` 不会隐式套用小鼠数据库、marker、染色体或 genome assembly。

### 6.2 当前内置范围

```r
human_resources <- species_resources("human")
mouse_resources <- species_resources("mouse")
unknown_resources <- species_resources("unknown")
```

| 资源字段 | human | mouse | unknown/其他未注册物种 |
|---|---|---|---|
| OrgDb | `org.Hs.eg.db` | `org.Mm.eg.db` | 无 |
| KEGG code | `hsa` | `mmu` | 无 |
| CellChat DB | human | mouse | 无 |
| genome/TxDb | GRCh38/hg38 | GRCm38/mm10 | 无 |
| cell-cycle/TF/gene-set strategy | 人源策略 | 鼠源策略 | `user_supplied` |
| 自动注释入口 | human celldex reference | mouse celldex reference | 无 |

内置的是“资源映射”，不代表相关数据库包已安装。模块运行前仍会检查 package、基因 ID overlap 和 genome assembly。

### 6.3 其他物种

其他物种仍可执行不依赖物种数据库的读取、QC、降维、聚类和下载。高级模块需要显式覆盖资源，例如：

```r
cfg <- report_config(
  resource_overrides = list(
    species = "your_species",
    orgdb = "your.OrgDb.package",
    feature_keytype = "ENSEMBL",
    kegg_code = "your_kegg_code",
    cellchat_db = NULL,
    txdb = NULL,
    gtf = "/absolute/path/to/matching_annotation.gtf",
    gene_sets = your_gene_sets,
    manual_markers = your_marker_list
  )
)

running("object.rds", "result", species = "your_species", config = cfg)
```

资源必须与 feature ID 和 genome assembly 一致。给出 OrgDb 并不会自动补齐 CellChat、轨迹、CNV 和基因集资源；每个模块会独立报告是否可运行。不要把本文的人/鼠示例直接用于其他物种。

## 7. trajectory_root

轨迹起点决定 pseudotime 的方向，必须有实验或生物学依据。最简单的配置是 Monocle principal node：

```r
cfg <- report_config(
  trajectory_root = "Y_1"
)
```

也可以传入命名 list 描述 root cells、marker 或元数据组，例如：

```r
cfg <- report_config(
  trajectory_root = list(
    type = "metadata_group",
    column = "celltype_final",
    value = "Stem_like"
  )
)
```

起点不能在对象中解析、细胞太少、降维/图结构不可用或 `monocle3` 未安装时，`pseudotime` 会显示具体原因。程序不会选择“看起来最早”的 cluster 来代替用户决定。

## 8. cnv_reference

inferCNV 的 reference 必须是真实的正常细胞群名称，并且能在选定注释列中找到：

```r
cfg <- report_config(
  cnv_reference = c("T_cell", "B_cell", "Endothelial")
)
running("tumor.rds", "tumor_result", species = "human", config = cfg)
```

还需要与对象 feature ID 和 genome assembly 匹配的 TxDb 或 GTF。未提供 reference 时状态为 `cnv_reference_required`；reference 不在注释中、基因坐标无法匹配或 `infercnv` 未安装时，章节也会明确跳过/失败。CNV 信号是表达推断结果，不应直接等同于 DNA 层面的拷贝数验证。

## 9. 资源限制和超大 raw 对象

`report_config(limits = ...)` 支持：

| 字段 | 默认值 | 含义 |
|---|---:|---|
| `analysis_max_cells` | `Inf` | 送入 SCP/高级分析的最大细胞数 |
| `analysis_max_features` | `Inf` | 送入分析的最大 feature 数 |
| `plot_max_cells` | `100000` | 单张散点图最多绘制的细胞数 |
| `marker_max_cells_per_ident` | `2000` | 每个 identity 用于 marker 计算的最大细胞数 |
| `min_cells_per_group` | `20` | 模块级分组的最低细胞数 |
| `workers` | `1` | 允许的 worker 数；并行是否生效由具体引擎决定 |
| `embed_max_mb` | `50` | `embed_downloads = "auto"` 的总内置预算 |

大对象示例：

```r
cfg <- report_config(
  profile = "full",
  limits = list(
    analysis_max_cells = 20000,
    analysis_max_features = 15000,
    plot_max_cells = 50000,
    marker_max_cells_per_ident = 1000,
    workers = 4,
    embed_max_mb = 50
  )
)

running(
  input = "combined_raw_matrix.seurat.rds",
  output = "result",
  config = cfg,
  filter_raw_barcodes = "auto",
  filter_low_expression_features = "auto",
  integration_method = "none",
  scp_args = list(nHVF = 2000, linear_reduction_dims = 30),
  embed_downloads = "auto"
)
```

安全边界：

- `filter_raw_barcodes = "auto"` 对需要补分析的 raw/partial 对象会在 SCP 之前使用给定 QC 阈值，保证降维、marker 和最终对象使用相同细胞集合；已有分析对象不会被自动过滤；
- 已有分析的对象默认不因 feature 过滤而被改变；
- 分析细胞/feature 上限只限制计算对象，不删除原始输入下载；
- 图形抽样只影响显示，不改变下载表和分析矩阵；
- manifest 记录过滤、分析子集和 feature 子集的前后规模；
- 强制把大型 RDS/矩阵 Base64 内置进 HTML 可能耗尽内存，应优先打包 output 文件夹；
- 内存紧张时可先设置 `render = FALSE` 完成分析/导出，再在资源充足的环境渲染。

默认 `integration_method = "none"`。sample 不一定等于 batch，程序不会仅凭样本名自动进行整合；只有用户明确指定方法时才交给 SCP。

## 10. 表、矩阵和解释单位

报告遵循下列方向，避免把细胞、样本和 feature 混在一起：

- 表达矩阵：行是 feature/基因，列是 cell barcode；
- features 文件：顺序对应矩阵行；
- barcodes 文件：顺序对应矩阵列；
- 细胞元数据/注释表：一行一个细胞，cell barcode 是键，其他列是元数据字段；
- sample design：一行一个生物学样本；
- pseudobulk counts：行是 feature，列是生物学样本；
- 差异结果：一行一个 feature-contrast 记录；
- 降维坐标：一行一个细胞，列为坐标和可用的注释/cluster 字段。

每个表格产物都应有：

1. 表的用途；
2. 一行表示什么；
3. 一列表示什么或列单位；
4. 关键字段的数据字典；
5. 结果是完整表、显示预览还是分析子集；
6. 对应下载文件和校验/来源信息。

DT 表用于浏览、搜索和小表下载。超大矩阵不直接塞入浏览器，报告显示预览并提供压缩 Matrix Market/TSV 文件；预览行数不代表下载文件行数。

## 11. standalone HTML 和 output 目录

`report.html` 内置页面 CSS、JavaScript、图形和 DT 组件，可离线打开。可下载文件是否也内置由 `embed_downloads` 决定：

| 值 | 行为 | 推荐用途 |
|---|---|---|
| `auto` | 在 `embed_max_mb` 总预算内优先内置注释和结果表；大文件保留在 output | 默认，兼顾可转发和体积 |
| `always` | 尝试把所有下载文件内置 | 只用于确认很小的对象；HTML 会膨胀且占浏览器内存 |
| `never` | 下载项指向 output 中的文件 | 服务器或整个文件夹一起交付 |

标准目录为：

```text
output/
├── report.html
├── analysis/<module>/        # 模块表、矩阵和图形
├── downloads/                # 原始和分析 RDS
├── matrices/                 # 原始/分析表达矩阵及行列名
├── tables/                   # 元数据、注释、设计、降维和其他通用表
├── figures/                  # 通用图形
└── .report/
    ├── manifest.rds
    └── manifest.json
```

`output` 是完整交付单元。只发送 `report.html` 时，报告仍能打开，但只能下载已内置的文件；大型项目应压缩并发送整个 output 文件夹。输入原始 RDS、原始 counts 与分析对象/分析矩阵分别导出，避免用户误把过滤或分析子集当成原始数据。

## 12. 依赖分层

完整安装命令见项目 [README](../README.md#安装)。模块和主要可选包的对应关系如下：

| 层/模块 | 主要包 |
|---|---|
| 核心读取、SCP、报告和下载 | `Seurat`、`SeuratObject`、`SCP`、`Matrix`、`BiocParallel`、`DT`、`ggplot2`、`ggsci`、`quarto`、`knitr`、`htmltools`、`jsonlite`、`digest` |
| 自动注释 | `SingleR`、`celldex`、`AnnotationDbi`、对应 `org.*.eg.db` |
| 重复感知差异 | `edgeR` |
| 富集/基因集 | `clusterProfiler`、`enrichplot`、`GSVA`、`ComplexHeatmap`、对应 OrgDb |
| 轨迹 | `monocle3`、`igraph`，必要时 `SeuratWrappers` |
| 通讯 | `CellChat`、`circlize`、`ComplexHeatmap` |
| TF/热图 | `AnnotationDbi`、对应 OrgDb、`ComplexHeatmap` |
| CNV | `infercnv`、`GenomicFeatures`、匹配的 TxDb/GTF |
| 可选交互图和数据整理 | `plotly`、`htmlwidgets`、`data.table`、`dplyr`、`tidyr`、`stringr`、`ggrepel`、`scales`、`future` |

可选依赖缺失不会触发静默替代。例如没有 `edgeR` 时不会伪造 edgeR P 值，没有 `CellChat` 数据库时不会生成通讯网络，没有 `infercnv` 或参考组时不会生成 CNV 热图。报告应把这些情况作为状态和原因呈现。

## 13. 发布前最小核对

对新数据运行后至少检查：

1. `report.html` 能在断网环境打开；
2. manifest 中 12 个模块都存在，未运行模块有原因；
3. sample design 与真实实验设计一致，`needs_review` 已处理；
4. species、OrgDb 和 genome assembly 匹配；
5. 原注释未被覆盖，自动/手动注释有单独来源列；
6. 表达矩阵是 feature × cell，features/barcodes 顺序一致；
7. 差异分析的统计单位是生物学样本；无重复结果没有正式 P 值；
8. pseudotime 起点和 CNV reference 有生物学依据；
9. 分析子集、图形抽样和完整下载的边界写入报告；
10. 大文件是 HTML 内置还是 output 文件，交付方式与报告一致。
