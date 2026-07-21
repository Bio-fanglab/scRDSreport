# scRDSreport

[![R-CMD-check](https://github.com/Biobunengsi/scRDSreport/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Biobunengsi/scRDSreport/actions/workflows/R-CMD-check.yaml)
[![Example report](https://img.shields.io/badge/online-example-1177A3)](https://biobunengsi.github.io/scRDSreport/report.html)
[![License: MIT](https://img.shields.io/badge/license-MIT-1F6F5C.svg)](LICENSE)

**一个函数，把单细胞 RDS 整理成可浏览、可追溯、可下载的 Quarto HTML 报告。**

`scRDSreport` 会先检查输入对象：已有降维、聚类或细胞注释时优先保留并导出；只有原始或部分处理数据时，可调用 [SCP](https://github.com/zhanghao-njmu/SCP) 补齐基础分析。所有结果写入用户指定的输出文件夹，报告入口固定为 `report.html`。

**[查看 10x Genomics 官方示例数据生成的在线报告](https://biobunengsi.github.io/scRDSreport/report.html)** · [快速开始](#快速开始) · [安装](#安装) · [完整分析说明](docs/full-analysis.md)

```r
library(scRDSreport)

running(
  input = "path/to/object.rds",
  output = "path/to/result"
)
```

## 主要功能

- 保留并导出输入 RDS 中已有的降维、聚类、注释和分析结果，不覆盖原始元数据列。
- 对 raw/partial 对象按需使用 SCP 完成 QC、标准化、降维、聚类和 marker 分析。
- 以独立模块组织细胞注释、差异表达、富集、轨迹、细胞通讯、细胞周期、转录因子和 CNV；缺少生物学前提时明确标记 `skipped` 或 `needs_input`。
- 导出原始 RDS、分析对象、Matrix Market 表达矩阵、细胞元数据、注释、结果表和图形，同时记录参数、软件版本、校验值与模块状态。
- 生成可独立打开的 Quarto HTML；DT 表格支持搜索、复制以及 CSV/Excel 下载。
- 内置人、小鼠、大鼠、斑马鱼、猪、牛、鸡、狗和猕猴九种常见物种的分层资源配置，并允许用户提供自定义物种资源。

## 在线示例

在线示例由以下三个 10x Genomics 官方公开的 GEM-X Universal 5′ 小鼠数据集生成：

- [5k Mouse PBMCs](https://www.10xgenomics.com/datasets/5k_Mouse_PBMCs_5p_gem-x)
- [5k Mouse Splenocytes](https://www.10xgenomics.com/datasets/5k_Mouse_Splenocytes_5p_gem-x)
- [20k Mouse PBMCs and Mouse Splenocytes Multiplex Sample](https://www.10xgenomics.com/datasets/20k_Mouse_PBMC_Splenocyte_5p_gem-x_multiplex)

三个公开数据集被整理为一个 raw Seurat RDS，再交给 `scRDSreport` 完成检查、过滤、基础分析、参考注释和报告导出。为控制公开演示的计算量，下游展示使用固定抽取的 2,500 个细胞和 8,000 个 features；原始 RDS、分析 RDS 及原始/分析表达矩阵仍可从报告内下载。

这三个数据集不是经过验证的生物学重复，因此示例中的组间差异仅作为描述性效应量和基因排序解读，不作为重复感知的统计推断。示例报告由 scRDSreport 0.3.1 生成；仓库当前代码为更新后的 0.3.2。

数据版权归 10x Genomics，原始数据按 [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) 发布；本项目对数据进行了重新整理和分析。10x Genomics 未参与或认可本 R 包。

## 安装

需要 R >= 4.2。建议在新的 R 会话中使用仓库提供的一键安装器，它会安装包本身、当前分析档位和物种需要的 R 依赖，并在结束时实际加载依赖进行验收：

```r
source(
  "https://raw.githubusercontent.com/Biobunengsi/scRDSreport/main/inst/install_scRDSreport.R"
)

install_scRDSreport(
  input = "path/to/object.rds",
  profile = "full",
  species = "auto"
)
```

已经明确物种时不必读取 RDS 来选择依赖，例如：

```r
install_scRDSreport(profile = "full", species = "mouse")
```

安装器会处理 CRAN、Bioconductor 和 GitHub 上的 R 包。Quarto CLI、JAGS 和 HDF5 开发库属于操作系统软件，安装器会检测并给出提示，但不会跨平台静默安装：

- 生成 HTML 需要 [Quarto CLI >= 1.3](https://quarto.org/docs/get-started/)。
- inferCNV 需要 JAGS 4.x。
- 当前 Monocle 3/BPCells 安装路径需要可检测的 HDF5 开发库。

<details>
<summary>分步安装</summary>

```r
install.packages(c("BiocManager", "remotes"))

remotes::install_github(
  "Biobunengsi/scRDSreport",
  dependencies = NA,
  upgrade = "never",
  repos = BiocManager::repositories()
)

scRDSreport::install_dependencies(
  profile = "full",
  input = "path/to/object.rds",
  species = "auto"
)

scRDSreport::check_dependencies(
  profile = "full",
  input = "path/to/object.rds",
  species = "auto"
)
```

</details>

## 快速开始

最简单的调用只需要输入 RDS 和一个新的或空的输出目录：

```r
library(scRDSreport)

running(
  input = "object.rds",
  output = "result"
)
```

对于需要明确运行完整模块的小鼠对象：

```r
cfg <- report_config(
  profile = "full",
  annotation_mode = "auto_if_missing",
  differential = "auto"
)

running(
  input = "mouse_object.rds",
  output = "mouse_result",
  species = "mouse",
  config = cfg
)
```

`output` 应使用专用结果目录。程序只允许覆盖带有 scRDSreport 所有权标记的既有输出目录，并拒绝 HOME、文件系统根目录和无法识别的非空项目目录，避免误删无关文件。

## 支持的 RDS 输入

- `Seurat` 对象。
- `SingleCellExperiment` 对象。
- base R matrix 或 `Matrix` 稀疏矩阵。
- 含 `counts`、`count`、`raw_counts`、`matrix`、`expr` 或 `expression` 的 list，可同时包含 `metadata`、`meta.data`、`colData` 或 `cell_metadata`。
- 恰好包含一个 Seurat 对象的 list。

直接输入矩阵时必须是“行 = feature/基因，列 = 细胞”。不确定对象状态时可先运行：

```r
inspect_rds("object.rds")
```

## 自动判断与关键参数

| 参数 | 默认值 | 作用 |
|---|---|---|
| `analyze` | `"auto"` | 自动判断是否需要 SCP；也可设为 `"always"` 或 `"never"` |
| `config` | 自动选择 | raw/partial 对象默认 `full`；已有 reduction 且有 cluster/annotation 的对象默认 `report_only` |
| `species` | `"auto"` | 自动检测或显式指定九个内置物种之一；不能可靠判断时停止而不是猜测 |
| `annotation_mode` | `"auto_if_missing"` | 已有注释就保留；缺失时才尝试同物种参考 |
| `embed_downloads` | `"auto"` | 在默认预算内把结果文件嵌入 HTML；也可使用 `"always"` 或 `"never"` |
| `sample_map` | `NULL` | 提供真实的样本、实验组和生物学重复关系 |

如果必须只发送一个 HTML 并让所有导出文件都内置，可使用 `embed_downloads = "always"`。原始 RDS 和大型表达矩阵会让 HTML 变得非常大，并显著增加渲染时间、磁盘占用和浏览器内存消耗。

## 分析档位与模块

```r
report_config(profile = "full")        # 请求全部模块
report_config(profile = "core")        # QC、降维、聚类、注释、组成和下载
report_config(profile = "report_only") # 只整理输入对象已有结果
```

完整档包含 12 个固定模块：

```text
qc, reduction, cluster, celltype, differential, enrichment,
pseudotime, communication, cell_cycle, tf, cnv, downloads
```

模块不会为了填满报告而伪造结果。例如，轨迹分析没有起点时只导出无方向几何，CNV 没有正常参考时只输出 readiness，缺少可靠生物学重复时差异模块只给描述性效应量。

## 输出目录

```text
result/
├── report.html          # 固定报告入口
├── downloads/           # 原始 RDS 与分析后 RDS
├── matrices/            # Matrix Market 表达矩阵及行列索引
├── tables/              # 细胞元数据、注释和对象内已有结果表
└── analysis/            # 各分析模块生成的表、图和中间结果
```

表达矩阵保持“行 = feature/基因，列 = 细胞”；细胞元数据和注释表保持“一行一个细胞”。报告内每张表都说明行、列、统计量、阅读方式和注意事项，并提供可展开的数据字典。

## 样本分组与统计边界

程序可以从样本名提出保守的 `group`/`replicate` 建议，但自动结果不能代替实验设计。正式组间比较应提供 `sample_map`：

```r
sample_map <- data.frame(
  sample_id = c("Ctrl_1", "Ctrl_2", "Drug_1", "Drug_2"),
  group = c("Control", "Control", "Drug", "Drug"),
  replicate = c("1", "2", "1", "2")
)

running("object.rds", "result", sample_map = sample_map)
```

每组至少两个独立生物学样本时，`differential = "auto"` 才能选择重复感知的 pseudobulk 推断。没有生物学重复时，细胞不会被当作重复，报告也不会制造 P 值或 FDR。

## 物种资源

```r
supported_species()
species_resources("mouse")
```

- human、mouse：内置完整资源入口，包括可信参考注释、CellChat 和基因组坐标资源。
- rat、zebrafish、pig、cattle、chicken、dog、macaque：提供物种匹配的核心 ID、QC、OrgDb、KEGG/MSigDB 和细胞周期映射；高级模块需要用户提供经过验证的同物种资源。
- 其他物种：仍可进行通用读取、QC、降维、聚类和下载，物种特异模块通过 `resource_overrides` 配置。

`species = "auto"` 只有在证据可靠时才选择物种，不会把不确定对象默认套用人或小鼠资源。

## 文档

- [完整分析、资源和模块配置](docs/full-analysis.md)
- [样本分组算法与复核规则](docs/sample-grouping-algorithm.md)
- R 内帮助：`?running`、`?report_config`、`?install_dependencies`

## License

scRDSreport 使用 [MIT License](LICENSE)。
