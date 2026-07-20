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
