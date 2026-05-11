# tests/testthat/test-analysis-modules.R — unit tests for analysis pipeline modules

# Shared setup: load pre-computed annotation and expression data
rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
has_data <- file.exists(rdata_path) && rdata_path != ""
if (has_data) {
  temp_env <- new.env()
  load(rdata_path, envir = temp_env)
  raw_annotation <- temp_env[[ls(temp_env)[1]]]

  gene_sets <- looplook:::extract_target_gene_sets(raw_annotation,
    src = "targets",
    include_Filled = TRUE
  )
  target_genes <- head(intersect(
    gene_sets[["Target_Genes"]],
    rownames(looplook:::read_robust_general(
      system.file("extdata", "example_tpm.txt", package = "looplook"),
      header = TRUE, row_name = 1, min_cols = 1
    ))
  ), 12)

  diff_df <- looplook:::read_robust_general(
    system.file("extdata", "example_deg.txt", package = "looplook"),
    header = TRUE, row_name = 1, min_cols = 1
  )
  global_glist <- sort(setNames(
    diff_df[["log2FoldChange"]],
    rownames(diff_df)
  ), decreasing = TRUE)
}

test_that("run_lfc_violin returns ggplot with valid input", {
  skip_if_not(has_data, "Pre-computed RData not available")
  skip_if_not_installed("ggplot2")

  p <- looplook:::run_lfc_violin(target_genes, global_glist, "wilcox.test", "Test_Violin")
  expect_s3_class(p, "ggplot")
  # fewer than 3 genes returns NULL
  p_null <- looplook:::run_lfc_violin(
    head(target_genes, 2), global_glist,
    "wilcox.test", "TooFew"
  )
  expect_null(p_null)
})

test_that("run_gsea_analysis returns list with result and plot", {
  skip_if_not(has_data, "Pre-computed RData not available")
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("enrichplot")

  out <- looplook:::run_gsea_analysis(
    target_genes, global_glist, 50,
    "Test_GSEA"
  )
  expect_type(out, "list")
  expect_true("result" %in% names(out))
  expect_true("plot" %in% names(out))
  # fewer than 2 overlap should return NULL
  out_null <- looplook:::run_gsea_analysis(
    c("FAKE1", "FAKE2"), global_glist, 50,
    "Empty"
  )
  expect_null(out_null$result)
})

test_that("run_go_enrichment returns list with result and plot", {
  skip_if_not(has_data, "Pre-computed RData not available")
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Hs.eg.db")
  library(org.Hs.eg.db)

  out <- looplook:::run_go_enrichment(target_genes, "org.Hs.eg.db", global_glist,
    cnet_nSample = 10, project_name = "Test_GO"
  )
  expect_type(out, "list")
  expect_true("result" %in% names(out))
  # result should be a data.frame with expected columns
  if (!is.null(out$result) && nrow(out$result) > 0) {
    expect_s3_class(out$result, "data.frame")
  }
})

test_that("run_ppi_analysis returns ggplot for valid input", {
  skip_if_not(has_data, "Pre-computed RData not available")
  skip_if_not_installed("STRINGdb")
  skip_if_not_installed("ggraph")
  run_network_tests <- identical(
    tolower(Sys.getenv("LOOPLOOK_RUN_NETWORK_TESTS", unset = "false")),
    "true"
  )
  skip_if_not(
    run_network_tests,
    "External STRINGdb integration test disabled; set LOOPLOOK_RUN_NETWORK_TESTS=true to run."
  )
  # STRINGdb requires external network access even when explicitly enabled.
  has_net <- tryCatch(
    {
      con <- url("https://string-db.org", open = "r")
      close(con)
      TRUE
    },
    error = function(e) FALSE
  )
  skip_if_not(has_net, "STRINGdb not reachable")

  p <- looplook:::run_ppi_analysis(target_genes, global_glist, "org.Hs.eg.db",
    ppi_score = 700, ppi_ntop = 20, "Test_PPI"
  )
  # may return NULL if no interactions found — not a failure
  if (!is.null(p)) {
    expect_s3_class(p, "gg")
  }
})

test_that("run_heatmap_and_connectivity returns plot list", {
  skip_if_not(has_data, "Pre-computed RData not available")
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")

  tpm_mat <- looplook:::read_robust_general(
    system.file("extdata", "example_tpm.txt", package = "looplook"),
    header = TRUE, row_name = 1, min_cols = 1
  )
  meta_raw <- data.frame(
    SampleID = colnames(tpm_mat),
    Group = rep(c("A", "B"), length.out = ncol(tpm_mat)),
    stringsAsFactors = FALSE
  )

  plots <- looplook:::run_heatmap_and_connectivity(target_genes, tpm_mat, meta_raw,
    loop_stats_df = NULL, global_glist,
    heatmap_ntop = 50, cor_method = "pearson",
    current_proj_name = "Test_Heat", source_type = "targets"
  )
  expect_type(plots, "list")
})
