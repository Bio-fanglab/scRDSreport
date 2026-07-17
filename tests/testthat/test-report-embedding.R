test_that("deferred download payloads are selected and injected before body close", {
  output <- scRDSreport:::.ensure_output(tempfile("scrdsreport-embed-"), FALSE)
  first <- file.path(output, "downloads", "first.csv")
  second <- file.path(output, "downloads", "second.csv")
  writeLines(c("x", "1"), first, useBytes = TRUE)
  writeLines(c("y", "2"), second, useBytes = TRUE)
  sizes <- as.numeric(file.info(c(first, second))$size)
  files <- data.frame(
    section = c("module_qc", "expression_matrix"),
    path = c("downloads/first.csv", "downloads/second.csv"),
    bytes = sizes,
    embed_priority = c(1, 100),
    stringsAsFactors = FALSE
  )
  manifest <- list(
    files = files,
    download_embedding = list(mode = "auto", max_bytes = sizes[[1L]])
  )
  manifest_path <- file.path(output, ".report", "manifest.rds")
  saveRDS(manifest, manifest_path)
  report <- file.path(output, ".report", "report.html")
  writeLines("<html><body><p>正文</p></body></html>", report, useBytes = TRUE)

  expect_equal(scRDSreport:::.embedded_download_indices(manifest, output), 1L)
  expect_identical(
    scRDSreport:::.inject_embedded_downloads(report, manifest_path, output, verbose = FALSE),
    report
  )
  html <- paste(readLines(report, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  expect_match(html, "data-scrds-file=\"scrds-embedded-1\"", fixed = TRUE)
  expect_match(html, "data-filename=\"first.csv\"", fixed = TRUE)
  expect_false(grepl("data-filename=\"second.csv\"", html, fixed = TRUE))
  payload_position <- regexpr("data-scrds-file", html, fixed = TRUE)
  expect_lt(regexpr("<body>", html, fixed = TRUE), payload_position)
  expect_lt(payload_position, regexpr("</body>", html, fixed = TRUE))
  expect_match(html, "<p>正文</p>", fixed = TRUE)
})

test_that("always mode selects every existing download regardless of budget", {
  output <- scRDSreport:::.ensure_output(tempfile("scrdsreport-embed-all-"), FALSE)
  paths <- file.path(output, "downloads", c("a.txt", "b.txt"))
  writeLines("a", paths[[1L]])
  writeLines("b", paths[[2L]])
  manifest <- list(
    files = data.frame(
      path = c("downloads/a.txt", "downloads/b.txt"),
      bytes = as.numeric(file.info(paths)$size),
      stringsAsFactors = FALSE
    ),
    download_embedding = list(mode = "always", max_bytes = 1)
  )

  expect_equal(scRDSreport:::.embedded_download_indices(manifest, output), 1:2)
})

test_that("binary injection preserves UTF-8 bytes when body close crosses a chunk", {
  output <- scRDSreport:::.ensure_output(tempfile("scrdsreport-embed-boundary-"), FALSE)
  downloadable <- file.path(output, "downloads", "result.txt")
  writeLines("result", downloadable, useBytes = TRUE)
  manifest <- list(
    files = data.frame(
      path = "downloads/result.txt",
      bytes = as.numeric(file.info(downloadable)$size),
      stringsAsFactors = FALSE
    ),
    download_embedding = list(mode = "always", max_bytes = 1)
  )
  manifest_path <- file.path(output, ".report", "manifest.rds")
  saveRDS(manifest, manifest_path)

  report <- file.path(output, ".report", "report.html")
  opening <- charToRaw(enc2utf8("<html><body><p>跨块正文</p>"))
  prefix_length <- 1024L * 1024L - 3L
  original_prefix <- c(opening, rep(as.raw(32L), prefix_length - length(opening)))
  connection <- file(report, open = "wb")
  writeBin(c(original_prefix, charToRaw("</body></html>")), connection)
  close(connection)

  scRDSreport:::.inject_embedded_downloads(report, manifest_path, output, verbose = FALSE)
  connection <- file(report, open = "rb")
  rendered <- readBin(connection, what = "raw", n = file.info(report)$size)
  close(connection)

  expect_identical(rendered[seq_along(original_prefix)], original_prefix)
  payload <- scRDSreport:::.raw_fixed_match(rendered, charToRaw("data-scrds-file"))
  closing <- scRDSreport:::.raw_fixed_match(rendered, charToRaw("</body>"))
  expect_false(is.na(payload))
  expect_false(is.na(closing))
  expect_lt(payload, closing)
  expect_false(is.na(scRDSreport:::.raw_fixed_match(rendered, charToRaw(enc2utf8("跨块正文")))))
})
