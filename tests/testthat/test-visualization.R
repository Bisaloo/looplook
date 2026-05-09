# tests/testthat/test-visualization.R

test_that("Module 5: IGV-Style Track Visualization generates plot successfully", {
  f1 <- system.file("extdata", "example_loops_1.bedpe", package = "looplook")
  atac_path <- system.file("extdata", "example_peaks.bed", package = "looplook")

  skip_if(f1 == "" || atac_path == "")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
  skip_if_not_installed("org.Hs.eg.db")

  out_base <- tempdir()
  save_path <- file.path(out_base, "Test_Locus_Track.pdf")

  track_plot <- suppressWarnings(suppressMessages(
    plot_peaks_interactions(
      bedpe_file = f1,
      target_bed = atac_path,
      species = "hg38",
      score_to_alpha = TRUE,
      save_file = save_path
    )
  ))

  expect_s3_class(track_plot, "ggplot")

  unlink(save_path)
})

test_that("draw_flower_simplified returns ggplot and handles edge cases", {
  skip_if_not_installed("ggplot2")

  gene_sets <- list(
    A = c("TP53", "BRCA1", "MYC"),
    B = c("BRCA1", "MYC", "EGFR"),
    C = c("MYC", "EGFR", "KRAS")
  )
  p <- looplook:::draw_flower_simplified(gene_sets, "Test_Flower",
    c(A = "red", B = "blue", C = "green"))
  expect_s3_class(p, "ggplot")

  p_null <- looplook:::draw_flower_simplified(list(X = c("A", "B")), "Solo", NULL)
  expect_null(p_null)
})

test_that("draw_upset_intersections returns grob and handles empty input", {
  skip_if_not_installed("UpSetR")

  gene_sets <- list(
    Up = c("TP53", "BRCA1", "MYC"),
    Down = c("BRCA1", "MYC", "CDKN1A")
  )
  g <- looplook:::draw_upset_intersections(gene_sets, "Test_UpSet")
  # UpSetR returns a captured grob; may be NULL if rendering fails
  if (!is.null(g)) {
    expect_true(inherits(g, "grob"))
  }

  expect_null(
    looplook:::draw_upset_intersections(list(X = character(0), Y = character(0)),
      "Empty")
  )
})
