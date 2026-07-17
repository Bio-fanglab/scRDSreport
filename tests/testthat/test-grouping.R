test_that("explicit replicate suffixes are high confidence", {
  x <- c("Ctrl_rep1", "Ctrl_rep2", "Drug_rep1", "Drug_rep2")
  result <- infer_sample_design(x)$design
  expect_equal(result$group, c("Ctrl", "Ctrl", "Drug", "Drug"))
  expect_equal(result$replicate, c("1", "2", "1", "2"))
  expect_true(all(result$confidence == "high"))
  expect_false(any(result$needs_review))
  expect_true(all(result$grouping_rule == "explicit_replicate_suffix"))
})

test_that("bare numeric suffixes are proposed but require review", {
  result <- infer_sample_design(c("WT_1", "WT_2", "KO_1", "KO_2"))$design
  expect_equal(result$group, c("WT", "WT", "KO", "KO"))
  expect_true(all(result$confidence == "medium"))
  expect_true(all(result$needs_review))
})

test_that("unreplicated samples remain independent", {
  x <- c("Control", "Treatment")
  result <- infer_sample_design(x)$design
  expect_equal(result$group, x)
  expect_true(all(is.na(result$replicate)))
  expect_true(all(result$needs_review))
})

test_that("user maps override inference", {
  x <- c("A", "B")
  map <- data.frame(sample_id = x, group = c("control", "drug"), replicate = c("1", "1"))
  result <- infer_sample_design(x, sample_map = map)$design
  expect_equal(result$group, c("control", "drug"))
  expect_equal(result$confidence, c("user", "user"))
})

test_that("subject tokens are proposed as replicate blocks", {
  x <- c("Patient1_Tumor", "Patient1_Normal", "Patient2_Tumor", "Patient2_Normal")
  result <- infer_sample_design(x)$design
  expect_equal(result$group, c("Tumor", "Normal", "Tumor", "Normal"))
  expect_true(all(result$confidence == "medium"))
  expect_true(all(result$needs_review))
})

test_that("SCP-style cluster columns are recognized", {
  meta <- data.frame(scRDSreportpcaclusters = c("0", "1"), check.names = FALSE)
  expect_equal(scRDSreport:::.cluster_columns(meta), "scRDSreportpcaclusters")
})
