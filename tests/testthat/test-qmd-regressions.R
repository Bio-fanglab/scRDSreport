test_that("Quarto embedding plots do not override mapped group colours", {
  template <- system.file("quarto", "report.qmd", package = "scRDSreport")
  expect_true(nzchar(template))
  code <- paste(readLines(template, warn = FALSE, encoding = "UTF-8"), collapse = "\n")

  expect_match(code, "Do not pass color = NULL here", fixed = TRUE)
  expect_false(grepl(
    "geom_point(size = 0.5, alpha = 0.78, color = if (is.null(color_col))",
    code,
    fixed = TRUE
  ))
  expect_match(code, "emit_html(figure_interpretation", fixed = TRUE)
  expect_match(code, "emit_html(interpretation_panel", fixed = TRUE)
  expect_false(grepl("print(figure_interpretation", code, fixed = TRUE))
  expect_false(grepl("print(interpretation_panel", code, fixed = TRUE))
  expect_match(code, "generated_annotation_sources", fixed = TRUE)
  expect_false(grepl("grepl(\"^explicit_\"", code, fixed = TRUE))
})

test_that("report HTML loops render as markup and file notes remain readable", {
  template <- system.file("quarto", "report.qmd", package = "scRDSreport")
  stylesheet <- system.file("quarto", "report.scss", package = "scRDSreport")
  expect_true(nzchar(template))
  expect_true(nzchar(stylesheet))

  code <- paste(readLines(template, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  css <- paste(readLines(stylesheet, warn = FALSE, encoding = "UTF-8"), collapse = "\n")

  expect_match(code, "toc-expand: true", fixed = TRUE)
  expect_false(grepl("collapsed: true", code, fixed = TRUE))
  expect_match(code, "```{r annotation-and-metadata, results='asis'}", fixed = TRUE)
  expect_match(code, "```{r inventory-before, results='asis'}", fixed = TRUE)
  expect_match(code, "class = \"file-index-view\"", fixed = TRUE)
  expect_match(code, "class='file-note-scroll'", fixed = TRUE)
  expect_match(code, "htmltools::htmlEscape(note)", fixed = TRUE)
  expect_match(code, "differential_statistics_text <- function", fixed = TRUE)
  expect_match(code, "statistics_text(label, module_id, x)", fixed = TRUE)
  expect_match(code, "没有执行显著性检验，也没有生成 P 值或 FDR", fixed = TRUE)
  expect_match(code, "它不是差异检验结果，本身不包含效应量、P 值或 FDR", fixed = TRUE)
  expect_match(code, "GSVA 分数", fixed = TRUE)
  expect_match(code, "因此不生成火山图；请查看下方效应量图和差异结果表", fixed = TRUE)
  expect_match(code, "basename(as.character(files$path[candidates]))", fixed = TRUE)
  expect_match(code, "if (has_fc && has_p) return(index)", fixed = TRUE)
  expect_match(code, "descriptive_index", fixed = TRUE)
  expect_false(grepl(
    'paste(files$label[candidates], files$path[candidates])', code,
    fixed = TRUE
  ))

  statistics_start <- regexpr("statistics_text <- function", code, fixed = TRUE)[[1L]]
  statistics_end <- regexpr("preview_notice <- function", code, fixed = TRUE)[[1L]]
  statistics_block <- substr(code, statistics_start, statistics_end - 1L)
  differential_branch <- regexpr('identical(module_id, "differential")', statistics_block, fixed = TRUE)[[1L]]
  marker_branch <- regexpr('grepl("marker", value)', statistics_block, fixed = TRUE)[[1L]]
  expect_gt(differential_branch, 0L)
  expect_gt(marker_branch, differential_branch)

  expect_match(css, ".file-index-view table.dataTable td:nth-child(7)", fixed = TRUE)
  expect_match(css, "min-width: 24rem !important", fixed = TRUE)
  expect_match(css, "height: 5.2rem", fixed = TRUE)
  expect_match(css, "overflow-y: auto", fixed = TRUE)
})

test_that("fallback artifact explanations do not claim analyses that were not run", {
  template <- system.file("quarto", "report.qmd", package = "scRDSreport")
  expect_true(nzchar(template))
  code <- paste(readLines(template, warn = FALSE, encoding = "UTF-8"), collapse = "\n")

  expect_match(code, "fallback_statistics_text <- function", fixed = TRUE)
  expect_match(code, "descriptive_gene_set_effect_summary", fixed = TRUE)
  expect_match(code, "它不是富集显著性检验，不包含 P 值或 FDR", fixed = TRUE)
  expect_match(code, "没有定向伪时间、起点或先后顺序", fixed = TRUE)
  expect_match(code, "它不包含配体-受体关系、通讯概率或细胞通讯推断", fixed = TRUE)
  expect_match(code, "空结果也不是生物学上不存在通讯的证据", fixed = TRUE)
  expect_match(code, "没有运行 CNV 推断，也不包含任何 CNV 信号", fixed = TRUE)
  expect_match(code, "how_to_read_text(files$label[index], id, table)", fixed = TRUE)
  expect_false(grepl(
    'how = "先核对标识列和分组列，再查看效应方向、统计量与校正显著性。"',
    code,
    fixed = TRUE
  ))
})
