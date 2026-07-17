# scRDSreport

`scRDSreport` 将单细胞 RDS 文件整理为一个可浏览、可追溯、可下载的 Quarto HTML 报告。入口函数只有 `running()`，用户只需提供输入 RDS 和输出文件夹。

```r
scRDSreport::running(
  input = "path/to/object.rds",
  output = "path/to/result"
)
```

报告固定写入 `output/report.html`，矩阵、RDS 和表格都保存在同一个 `output` 文件夹下。默认会在 50 MB 总预算内优先把注释与结果表内置进 HTML；因此单独发送 HTML 仍可打开，并可下载标为“HTML 内置”的文件。

## 主要行为

- 支持 Seurat、SingleCellExperiment、普通表达矩阵，以及含 `counts` 的 list。
- 自动合并 Seurat v5 的 split assay layers，避免只读取第一个样本层。
- 同时存在降维结果和 cluster/已有注释时，默认保留原分析，不重复计算。
- raw/partial 对象通过 `SCP::Standard_SCP()` 完成标准化、HVG、PCA、邻居图、聚类和 UMAP。
- cluster marker 优先使用 `SCP::RunDEtest()`；SCP 旧接口与 Seurat 5 不兼容时使用 `Seurat::FindAllMarkers()` 并在 manifest 中记录。
- 只识别并原样保留输入 RDS 中已有的细胞类型注释。包不会自行猜测、覆盖或创建 cell type。
- `species = "auto"` 根据稳定的 feature ID 前缀识别人、小鼠、大鼠、斑马鱼、猪、牛、鸡、猕猴、犬等物种；也可填写任意物种名称。
- 原始输入 RDS、原始 counts、分析对象 RDS、分析 counts/data、细胞元数据、已有注释、降维坐标和 marker 表均可下载；原始 RDS 元数据与分析对象元数据分开导出，不会因分析子集丢失原始细胞记录。
- 表达矩阵方向固定为：行是 feature/基因，列是细胞；features 和 barcodes 文件严格保持相同顺序。
- HTML 使用 Lumen 全宽样式、左侧目录、DT 下载按钮和 ggsci 科研配色；每张表后都有行、列与字段解释。

## 安装

先安装 Bioconductor 基础依赖和 SCP，再安装本包：

```r
install.packages(c("BiocManager", "remotes", "quarto"))

BiocManager::install(c(
  "BiocParallel",
  "SingleCellExperiment",
  "SummarizedExperiment"
))

remotes::install_github("zhanghao-njmu/SCP")
remotes::install_github("Bio-fanglab/scRDSreport")
```

还需要安装 [Quarto CLI](https://quarto.org/docs/get-started/)。SCP 的 PAGA、scVelo 等 Python 功能不是默认报告的必需步骤。

## 使用方法

### 自动处理

```r
library(scRDSreport)

running(
  input = "object.rds",
  output = "sc_report"
)
```

### 明确样本设计

自动分组是保守建议，正式组间比较推荐传入 `sample_map`：

```r
sample_map <- data.frame(
  sample_id = c("S1", "S2", "S3", "S4"),
  group = c("Control", "Control", "Drug", "Drug"),
  replicate = c("1", "2", "1", "2")
)

running("object.rds", "sc_report", sample_map = sample_map)
```

`grouping_rule` 只是程序说明“为什么这样建议分组”的短规则代码。`needs_review = TRUE` 表示需要人工复核，不代表分析错误。

### 指定已有注释列和物种

```r
running(
  "object.rds",
  "sc_report",
  annotation_col = "celltype_manual",
  species = "human"
)
```

不提供 `annotation_col` 时，程序会搜索常见的 `celltype`、`annotation`、`predicted.id` 等元数据列。未找到时报告会明确写“无已有注释”，不会生成自定义注释。

### 只导出已有结果

```r
running("object.rds", "sc_report", analyze = "never")
```

### 大对象资源控制

默认使用全部 QC 后细胞和 feature。资源受限时可显式建立分析子集；完整输入 RDS 和原始 counts 仍会导出：

```r
running(
  "object.rds",
  "sc_report",
  analysis_max_cells = 20000,
  analysis_max_features = 15000,
  scp_args = list(nHVF = 2000, linear_reduction_dims = 30)
)
```

### 单个 HTML 内置下载

默认 `embed_downloads = "auto"`：在 `embed_max_mb = 50` 的总预算内，优先内置注释、元数据和结果表，超出预算的 RDS/大型矩阵仍保留在 output 文件夹。若必须只发送一个 HTML，并确认磁盘和浏览器内存足够，可强制内置全部文件：

```r
running(
  "object.rds",
  "sc_report",
  embed_downloads = "always"
)
```

大型单细胞 RDS 和矩阵会使 HTML 体积增加约三分之一，并在浏览器下载时占用额外内存，因此常规共享仍建议打包整个 output 文件夹。

### 其他常用参数

```r
# 不运行 marker
running("object.rds", "sc_report", run_markers = FALSE)

# 自定义 raw barcode QC
running("object.rds", "sc_report", min_features = 300, min_counts = 500)

# 显式指定 SCP integration；默认不自动整合
running("object.rds", "sc_report", integration_method = "Harmony")

# 仅导出文件，暂不渲染 HTML
running("object.rds", "sc_report", render = FALSE)
```

## 输出结构

```text
output/
├── report.html
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
│   ├── original_cell_annotations.csv.gz   # 仅当输入 RDS 已有注释
│   ├── analysis_cell_metadata.csv.gz
│   ├── analysis_cell_annotations.csv.gz   # 分析对象仍含已有注释时
│   ├── sample_design.csv.gz
│   ├── matrix_preview_*.csv.gz
│   └── embedding_*.csv.gz
└── .report/
    ├── manifest.rds
    └── manifest.json
```

`report.html` 内嵌页面资源、图形和表格预览，并按 `embed_downloads` 策略内置完整下载文件。下载表中的 `storage` 列会明确标出“HTML 内置”或“output 文件”。

## 开发检查

```r
devtools::test()
devtools::check()
```

问题反馈：[Bio-fanglab/scRDSreport issues](https://github.com/Bio-fanglab/scRDSreport/issues)
