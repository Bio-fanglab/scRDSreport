test_that("overwrite refuses unowned non-empty directories without deleting files", {
  output <- tempfile("scrdsreport-unowned-")
  dir.create(file.path(output, "analysis"), recursive = TRUE)
  sentinel <- file.path(output, "analysis", "keep.txt")
  writeLines("unrelated", sentinel)

  expect_error(
    scRDSreport:::.ensure_output(output, overwrite = TRUE),
    "not recognized as an scRDSreport output"
  )
  expect_true(file.exists(sentinel))
  expect_equal(readLines(sentinel), "unrelated")
})

test_that("owned outputs clean only managed paths and preserve unrelated files", {
  output <- scRDSreport:::.ensure_output(tempfile("scrdsreport-owned-"), FALSE)
  sentinel <- file.path(output, "keep.txt")
  managed_file <- file.path(output, "analysis", "old.txt")
  writeLines("preserve", sentinel)
  writeLines("replace", managed_file)

  result <- scRDSreport:::.ensure_output(output, overwrite = TRUE)

  expect_identical(result, output)
  expect_true(file.exists(file.path(output, ".scRDSreport-output")))
  expect_true(file.exists(sentinel))
  expect_false(file.exists(managed_file))
  expect_true(dir.exists(file.path(output, "analysis")))
})

test_that("legacy manifests are accepted only for their recorded output path", {
  skip_if_not_installed("jsonlite")
  output <- tempfile("scrdsreport-legacy-")
  dir.create(file.path(output, ".report"), recursive = TRUE)
  output <- normalizePath(output, winslash = "/")
  jsonlite::write_json(
    list(
      created_at = "2026-01-01T00:00:00+0000",
      input = "/data/object.rds",
      output = output,
      status_before = "raw",
      status_after = "analyzed",
      files = list(list(path = "report.html"))
    ),
    file.path(output, ".report", "manifest.json"),
    auto_unbox = TRUE
  )
  writeLines("old", file.path(output, "report.html"))

  expect_silent(scRDSreport:::.ensure_output(output, overwrite = TRUE))
  expect_true(file.exists(file.path(output, ".scRDSreport-output")))
})

test_that("running protects an input stored under managed output paths", {
  output <- scRDSreport:::.ensure_output(tempfile("scrdsreport-input-guard-"), FALSE)
  input <- file.path(output, "downloads", "analysis_object.rds")
  saveRDS(list(value = 1), input)

  expect_error(
    running(input, output, analyze = "never", render = FALSE, overwrite = TRUE),
    "Input RDS is inside a managed scRDSreport output path"
  )
  expect_true(file.exists(input))
})

test_that("broad output targets are rejected before any cleanup", {
  expect_error(
    scRDSreport:::.ensure_output(path.expand("~"), overwrite = TRUE),
    "Refusing broad output directory"
  )

  missing <- file.path(
    path.expand("~"),
    paste0("scrdsreport-path-guard-", Sys.getpid())
  )
  expect_false(dir.exists(missing))
  expect_error(
    scRDSreport:::.ensure_output(file.path(missing, ".."), overwrite = TRUE),
    "must not contain a '..'"
  )
  expect_false(dir.exists(missing))
})
