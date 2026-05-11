# tests/testthat/test-annotation.R

test_that("packaged annotation example keeps the expected output contract", {
  rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
  skip_if(rdata_path == "", "Pre-computed test data not available")

  temp_env <- new.env()
  load(rdata_path, envir = temp_env)
  res_integrated <- temp_env[[ls(temp_env)[1]]]

  expect_type(res_integrated, "list")
  expect_true(all(c("target_annotation", "loop_annotation", "anchor_loci_annotation", "anchor_annotation", "plots") %in% names(res_integrated)))
  expect_type(res_integrated$plots, "list")
  expect_gt(nrow(res_integrated$loop_annotation), 0)
  expect_gt(nrow(res_integrated$target_annotation), 0)
  expect_true("Assigned_Target_Genes_Filled" %in% colnames(res_integrated$target_annotation))
})

test_that("annotate_peaks_and_loops respects quiet and write_output flags", {
  skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
  skip_if_not_installed("org.Hs.eg.db")
  invisible(requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE))
  invisible(requireNamespace("org.Hs.eg.db", quietly = TRUE))
  tiny_bedpe <- tempfile(fileext = ".bedpe")
  writeLines("chr1\t0\t100\tchr1\t200\t300", tiny_bedpe)

  out_base <- tempfile(pattern = "looplook_anno_nowrite_")
  unlink(out_base, recursive = TRUE, force = TRUE)
  expect_false(dir.exists(out_base))

  expect_no_message({
    res_integrated <- suppressPackageStartupMessages(
      suppressWarnings(
        annotate_peaks_and_loops(
          bedpe_file = tiny_bedpe,
          txdb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
          org_db = org.Hs.eg.db::org.Hs.eg.db,
          out_dir = out_base,
          project_name = "Tiny_NoWrite_Test",
          write_output = FALSE,
          quiet = TRUE
        )
      )
    )
  })

  expect_type(res_integrated, "list")
  expect_false(dir.exists(out_base))
})
