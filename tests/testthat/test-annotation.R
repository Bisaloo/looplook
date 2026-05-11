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
  skip_if_not_installed("org.Hs.eg.db")
  sample_txdb_path <- system.file(
    "extdata", "hg19_knownGene_sample.sqlite",
    package = "GenomicFeatures"
  )
  skip_if(sample_txdb_path == "", "Sample TxDb not available")
  txdb_obj <- AnnotationDbi::loadDb(sample_txdb_path)
  tiny_bedpe <- tempfile(fileext = ".bedpe")
  writeLines("chr6\t10412000\t10412600\tchr6\t10415000\t10415600", tiny_bedpe)

  out_base <- tempfile(pattern = "looplook_anno_nowrite_")
  unlink(out_base, recursive = TRUE, force = TRUE)
  expect_false(dir.exists(out_base))

  expect_no_message({
    res_integrated <- looplook:::.with_known_upstream_noise_suppressed(
      looplook::annotate_peaks_and_loops(
        bedpe_file = tiny_bedpe,
        txdb = txdb_obj,
        org_db = "org.Hs.eg.db",
        species = "hg19",
        out_dir = out_base,
        project_name = "Tiny_NoWrite_Test",
        write_output = FALSE,
        quiet = TRUE
      )
    )
  })

  expect_type(res_integrated, "list")
  expect_false(dir.exists(out_base))
})
