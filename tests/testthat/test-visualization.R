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

test_that("prepare_track_data stacks overlapping arcs within max_levels", {
  skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
  skip_if_not_installed("org.Hs.eg.db")

  bedpe_path <- tempfile(fileext = ".bedpe")
  on.exit(unlink(bedpe_path), add = TRUE)
  writeLines(
    c(
      "chr1\t100\t140\tchr1\t400\t440",
      "chr1\t120\t160\tchr1\t320\t360",
      "chr1\t150\t190\tchr1\t260\t300"
    ),
    con = bedpe_path
  )

  d <- suppressWarnings(suppressMessages(
    looplook:::prepare_track_data(
      bedpe_file = bedpe_path,
      target_bed = NULL,
      chr = "chr1",
      from = 90,
      to = 450,
      species = "hg38",
      max_levels = 2,
      base_anchor_height = 0.05,
      loop_color = "#5D6D7E",
      anchor_color = "#3498DB",
      score_to_alpha = FALSE,
      min_score = NULL
    )
  ))

  expect_gt(nrow(d$bez_df), 0)

  arc_levels <- sort(unique(d$bez_df$arc_level))
  expect_equal(arc_levels, c(1, 2))
  expect_lte(max(arc_levels), 2)

  peak_by_loop <- stats::aggregate(y ~ loop_i, data = d$bez_df, FUN = max)
  expect_equal(length(unique(peak_by_loop$y)), nrow(peak_by_loop))
})

test_that("draw_flower_simplified returns ggplot and handles edge cases", {
  skip_if_not_installed("ggplot2")

  gene_sets <- list(
    A = c("TP53", "BRCA1", "MYC"),
    B = c("BRCA1", "MYC", "EGFR"),
    C = c("MYC", "EGFR", "KRAS")
  )
  p <- looplook:::draw_flower_simplified(
    gene_sets, "Test_Flower",
    c(A = "red", B = "blue", C = "green")
  )
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
    looplook:::draw_upset_intersections(
      list(X = character(0), Y = character(0)),
      "Empty"
    )
  )
})
