test_that("cluster detection ignores all-NA result columns", {
  meta <- data.frame(
    brokenclusters = factor(c(NA, NA, NA)),
    graph_SNN_res.0.6 = factor(c("0", "1", "1"))
  )
  expect_equal(scRDSreport:::.cluster_columns(meta), "graph_SNN_res.0.6")
})
