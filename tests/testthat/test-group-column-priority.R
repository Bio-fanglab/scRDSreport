.group_priority_object <- function() {
  skip_if_not_installed("SeuratObject")
  counts <- Matrix::Matrix(
    matrix(
      seq_len(96L) %% 5L,
      nrow = 12L,
      dimnames = list(paste0("gene", seq_len(12L)), paste0("cell", seq_len(8L)))
    ),
    sparse = TRUE
  )
  object <- SeuratObject::CreateSeuratObject(counts = counts)
  object$group <- rep(c("rds_control", "rds_treatment"), each = 4L)
  object$condition <- rep(c("condition_control", "condition_treatment"), each = 4L)
  object$treatment <- rep(c("vehicle", "drug"), each = 4L)
  object
}

.attach_group_priority_design <- function(object, sample_map = NULL) {
  sample_ids <- rep(c("WT_1", "WT_2"), each = 4L)
  inferred <- infer_sample_design(sample_ids, sample_map = sample_map)
  scRDSreport:::.attach_design(object, inferred)
}

test_that("an explicit module group column has highest priority", {
  object <- .group_priority_object()
  sample_map <- data.frame(
    sample_id = c("WT_1", "WT_2"),
    group = c("mapped_control", "mapped_treatment"),
    stringsAsFactors = FALSE
  )
  object <- .attach_group_priority_design(object, sample_map)

  expect_equal(
    scRDSreport:::.fa_group_column(object, list(group_column = "treatment")),
    "treatment"
  )
})

test_that("a user sample map outranks original RDS group fields", {
  object <- .group_priority_object()
  sample_map <- data.frame(
    sample_id = c("WT_1", "WT_2"),
    group = c("mapped_control", "mapped_treatment"),
    stringsAsFactors = FALSE
  )
  object <- .attach_group_priority_design(object, sample_map)
  metadata <- scRDSreport:::.seurat_metadata(object)

  expect_true(all(metadata$.scRDSreport_grouping_rule == "user_map"))
  expect_true(all(metadata$.scRDSreport_design_confidence == "user"))
  expect_equal(scRDSreport:::.fa_group_column(object), ".scRDSreport_group")

  object$.scRDSreport_grouping_rule <- NULL
  expect_equal(scRDSreport:::.fa_group_column(object), ".scRDSreport_group")
})

test_that("original RDS group fields outrank automatic name inference", {
  object <- .attach_group_priority_design(.group_priority_object())
  metadata <- scRDSreport:::.seurat_metadata(object)

  expect_true(all(metadata$.scRDSreport_grouping_rule == "numeric_suffix_candidate"))
  expect_equal(scRDSreport:::.fa_group_column(object), "group")

  object$group <- NULL
  expect_equal(scRDSreport:::.fa_group_column(object), "condition")

  object$condition <- NULL
  expect_equal(scRDSreport:::.fa_group_column(object), "treatment")
})

test_that("automatic inferred groups remain the final fallback", {
  object <- .group_priority_object()
  object$group <- NULL
  object$condition <- NULL
  object$treatment <- NULL
  object <- .attach_group_priority_design(object)

  expect_equal(scRDSreport:::.fa_group_column(object), ".scRDSreport_group")
})
