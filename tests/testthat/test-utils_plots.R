# tests/testthat/test-utils_plots.R

test_that("Utility drawing functions render successfully with mock data", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("scales")
  library(ggplot2)
  library(dplyr)

  mock_loop <- data.frame(
    loop_type = c("E-P", "P-P", "E-E", "E-P", "P-P", "E-G"),
    functional_anchor1_type = c("E", "P", "E", "E", "P", "E"),
    anchor1_source = c("Native", "Native", "Native", "Promoter-derived enhancer", "Native", "Gene_body-derived enhancer"),
    functional_anchor2_type = c("P", "P", "E", "P", "P", "G"),
    anchor2_source = c("Native", "Native", "Native", "Native", "Native", "Native")
  )

  mock_target <- data.frame(
    annotation = c("Promoter (<=1kb)", "Distal Intergenic", "Intron", "Promoter (1-2kb)"),
    loop_genes_Total = c("GeneA", NA, "GeneB;GeneC", "GeneD")
  )

  mock_exp <- data.frame(
    loop_type = c("E-P", "E-P", "P-P", "E-E"),
    loop_genes = c("GeneA", "GeneB", "GeneC", "GeneD"),
    expression_value = c(10.5, 2.1, 0.5, 30.2)
  )

  mock_colors <- c("E-P" = "red", "P-P" = "blue", "E-E" = "green", "E-G" = "purple")

  proj <- "Test_Project"

  looplook:::.with_known_upstream_noise_suppressed({
    expect_s3_class(looplook:::draw_target_loop_donut(mock_loop, proj, NULL, mock_colors), "ggplot")
    expect_s3_class(looplook:::draw_target_annotation_pie(mock_target, proj, NULL), "ggplot")
    expect_s3_class(looplook:::draw_rose_plot(mock_loop, proj, NULL, mock_colors), "ggplot")
    expect_s3_class(looplook:::draw_target_connectivity_bar(mock_target, NULL, proj, NULL), "ggplot")
    expect_s3_class(looplook:::draw_enhancer_source_distribution(mock_loop, proj, NULL), "ggplot")
    expect_s3_class(looplook:::draw_expression_violin(mock_exp, proj, NULL, "TPM", mock_colors), "ggplot")
    expect_s3_class(looplook:::draw_comparison_bar(mock_loop, mock_loop, NULL, mock_colors), "ggplot")
  })
})

test_that("draw_circular_bar_plot returns ggplot with valid input", {
  skip_if_not_installed("ggplot2")
  mock <- data.frame(
    loop_type = c("E-P", "P-P", "E-E"),
    loop_genes = c("GeneA;GeneB", "GeneC", "GeneD;GeneE"),
    stringsAsFactors = FALSE
  )
  p <- looplook:::draw_circular_bar_plot(mock, "Test", NULL, c("E-P" = "red", "P-P" = "blue", "E-E" = "green"))
  expect_s3_class(p, "ggplot")
})

test_that("draw_pie_with_outside_labels handles empty and valid data", {
  skip_if_not_installed("ggplot2")
  expect_null(looplook:::draw_pie_with_outside_labels(data.frame(), "anno", "Empty", "Set2"))

  df <- data.frame(
    annotation = c("Promoter (<=1kb)", "Distal Intergenic", "Intron (1-2kb)"),
    stringsAsFactors = FALSE
  )
  p <- looplook:::draw_pie_with_outside_labels(df, "annotation", "Test", "Set2")
  expect_s3_class(p, "ggplot")
})
