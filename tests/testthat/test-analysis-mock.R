# tests/testthat/test-analysis-mock.R — mock-only tests for analysis.R internal functions

test_that(".plot_motif_rank_scatter generates scatter with mock data", {
  skip_if_not_installed("ggplot2")

  mock_res <- data.frame(
    MotifID = paste0("MA", sprintf("%04d", 1:20), ".1"),
    MotifName = paste("TF", LETTERS[1:20]),
    Pvalue = 10^-(1:20),
    FDR = c(10^-(1:10), rep(0.5, 10)),
    OddsRatio = runif(20, 0.8, 1.6),
    Family = c(
      rep("bZIP", 5), rep("Homeobox", 5), rep("Unknown", 5),
      rep("bHLH", 5)
    ),
    stringsAsFactors = FALSE
  )

  p <- looplook:::.plot_motif_rank_scatter(mock_res, "Test_Motif", fdr_thresh = 0.05)
  expect_s3_class(p, "ggplot")

  # empty data
  p_null <- looplook:::.plot_motif_rank_scatter(
    data.frame(
      MotifID = character(0), MotifName = character(0),
      Pvalue = numeric(0), FDR = numeric(0), OddsRatio = numeric(0),
      Family = character(0), stringsAsFactors = FALSE
    ),
    "Empty"
  )
  expect_null(p_null)


  # no Family column
  res_no_fam <- mock_res[, setdiff(colnames(mock_res), "Family")]
  p2 <- looplook:::.plot_motif_rank_scatter(res_no_fam, "NoFam")
  expect_s3_class(p2, "ggplot")
})

test_that(".plot_save_motif generates barplot with mock data", {
  skip_if_not_installed("ggplot2")

  mock_res <- data.frame(
    MotifID = paste0("MA", sprintf("%04d", 1:20), ".1"),
    MotifName = paste("TF", LETTERS[1:20]),
    Pvalue = 10^-(1:20),
    FDR = 10^-(1:20),
    OddsRatio = rep(1.2, 20),
    stringsAsFactors = FALSE
  )

  p <- looplook:::.plot_save_motif(mock_res, "Test_Bar")
  expect_s3_class(p, "ggplot")

  # NULL returns NULL
  expect_null(looplook:::.plot_save_motif(NULL, "Null"))

  # empty data
  expect_null(looplook:::.plot_save_motif(
    mock_res[integer(0), ], "Empty"
  ))
})

test_that("run_heatmap_and_connectivity handles distal-only mode with mock data", {
  skip_if_not_installed("ggpointdensity")
  skip_if_not_installed("viridis")

  genes <- paste0("Gene", 1:30)
  tpm_mat <- as.data.frame(matrix(rnorm(30 * 4), 30, 4,
    dimnames = list(genes, c("s1", "s2", "s3", "s4"))
  ))
  meta_raw <- data.frame(
    SampleID = c("s1", "s2", "s3", "s4"),
    Group = c("A", "A", "B", "B"), stringsAsFactors = FALSE
  )

  loop_stats <- data.frame(
    Gene = genes[1:15],
    Total_Loops = sample(1:10, 15, replace = TRUE),
    n_Linked_Distal = sample(0:5, 15, replace = TRUE),
    stringsAsFactors = FALSE
  )
  global_glist <- setNames(rnorm(30), genes)
  targets <- genes[1:10]

  plots <- looplook:::run_heatmap_and_connectivity(
    targets, tpm_mat, meta_raw, loop_stats, global_glist,
    heatmap_ntop = 50, cor_method = "pearson",
    current_proj_name = "Test_Distal", source_type = "targets",
    target_col = "n_Linked_Distal", skip_heatmap = TRUE
  )
  expect_type(plots, "list")
})

test_that(".annotate_motif_families handles empty input", {
  res_empty <- looplook:::.annotate_motif_families(NULL)
  expect_null(res_empty)

  res_no_rows <- looplook:::.annotate_motif_families(
    data.frame(
      MotifID = character(0), MotifName = character(0),
      Pvalue = numeric(0), stringsAsFactors = FALSE
    )
  )
  expect_equal(nrow(res_no_rows), 0)
})
