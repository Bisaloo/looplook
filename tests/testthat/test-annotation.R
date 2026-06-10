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

# --- Parameter validation ---
test_that("annotate_peaks_and_loops validates new parameters", {
  skip_if_not_installed("org.Hs.eg.db")
  sample_txdb_path <- system.file("extdata", "hg19_knownGene_sample.sqlite", package = "GenomicFeatures")
  skip_if(sample_txdb_path == "", "Sample TxDb not available")
  txdb_obj <- AnnotationDbi::loadDb(sample_txdb_path)
  tiny_bedpe <- tempfile(fileext = ".bedpe")
  writeLines("chr6\t10412000\t10412600\tchr6\t10415000\t10415600", tiny_bedpe)
  base_args <- list(
    bedpe_file = tiny_bedpe, txdb = txdb_obj, org_db = "org.Hs.eg.db",
    species = "hg19", out_dir = tempdir(), write_output = FALSE, quiet = TRUE
  )

  expect_error(do.call(annotate_peaks_and_loops, c(base_args, list(anchor_gap = -2))))
  expect_error(do.call(annotate_peaks_and_loops, c(base_args, list(anchor_min_overlap = 0))))
  expect_error(do.call(annotate_peaks_and_loops, c(base_args, list(anchor_min_frac = 1.5))))
  expect_error(do.call(annotate_peaks_and_loops, c(base_args, list(anchor_min_frac = -0.1))))
  expect_error(do.call(annotate_peaks_and_loops, c(base_args, list(hub_percentile = 0))))
  expect_error(do.call(annotate_peaks_and_loops, c(base_args, list(hub_percentile = 1.5))))
  expect_error(do.call(annotate_peaks_and_loops, c(base_args, list(neighbor_hop = -1))))
  expect_error(do.call(annotate_peaks_and_loops, c(base_args, list(neighbor_hop = 1.5))))
  expect_error(do.call(annotate_peaks_and_loops, c(base_args, list(karyo_bin_size = 0))))

  # Valid values should not error
  expect_no_error(do.call(annotate_peaks_and_loops, c(base_args, list(anchor_gap = 200L))))
  expect_no_error(do.call(annotate_peaks_and_loops, c(base_args, list(anchor_min_overlap = 10L))))
  expect_no_error(do.call(annotate_peaks_and_loops, c(base_args, list(anchor_min_frac = 0.5))))
})

# --- anchor_gap proximity matching ---
test_that("anchor_gap > 0 allows proximity-based peak-anchor linking", {
  skip_if_not_installed("org.Hs.eg.db")
  sample_txdb_path <- system.file("extdata", "hg19_knownGene_sample.sqlite", package = "GenomicFeatures")
  skip_if(sample_txdb_path == "", "Sample TxDb not available")
  txdb_obj <- AnnotationDbi::loadDb(sample_txdb_path)

  target_bed <- tempfile(fileext = ".bed")
  loop_bedpe <- tempfile(fileext = ".bedpe")
  # Peak 200bp away from anchor, 0bp actual overlap
  writeLines("chr6\t10412800\t10413000", target_bed)
  writeLines("chr6\t10412000\t10412600\tchr6\t10415000\t10415600", loop_bedpe)

  # anchor_gap=200: proximity hit (peak within 200bp of anchor) SHOULD link
  res <- annotate_peaks_and_loops(
    bedpe_file = loop_bedpe, target_bed = target_bed,
    txdb = txdb_obj, org_db = "org.Hs.eg.db", species = "hg19",
    anchor_gap = 200L, anchor_min_overlap = 1L,
    out_dir = tempdir(), write_output = FALSE, quiet = TRUE
  )
  has_loop <- !is.na(res$target_annotation$Linked_Loop_IDs) &
              res$target_annotation$Linked_Loop_IDs != ""
  expect_true(any(has_loop),
    info = "Peak within anchor_gap should be linked via proximity matching")

  # Default (anchor_gap=-1L): strict overlap only, should NOT link
  res2 <- annotate_peaks_and_loops(
    bedpe_file = loop_bedpe, target_bed = target_bed,
    txdb = txdb_obj, org_db = "org.Hs.eg.db", species = "hg19",
    out_dir = tempdir(), write_output = FALSE, quiet = TRUE
  )
  has_loop2 <- !is.na(res2$target_annotation$Linked_Loop_IDs) &
               res2$target_annotation$Linked_Loop_IDs != ""
  expect_false(any(has_loop2),
    info = "Peak without physical overlap should NOT be linked with default strict mode")

  unlink(c(target_bed, loop_bedpe))
})

# --- validate_epeG_by_chromatin: evidence/missing-data logic ---
# Uses pre-computed annotation results with guaranteed P/G anchors.
test_that("validate_epeG_by_chromatin: no marks → all uncertain", {
  rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
  skip_if(rdata_path == "", "Pre-computed test data not available")
  tmp <- new.env()
  load(rdata_path, envir = tmp)
  res <- tmp[[ls(tmp)[1]]]

  val <- validate_epeG_by_chromatin(res, chromatin_beds = list(), quiet = TRUE)
  skip_if(nrow(val) == 0, "No P/G anchors in pre-computed data")
  expect_true(all(val$confidence == "uncertain"))
  expect_true(all(is.na(val$H3K4me1) & is.na(val$H3K27ac) &
                  is.na(val$ATAC) & is.na(val$H3K27me3) & is.na(val$H3K4me3)))
})

test_that("validate_epeG_by_chromatin: missing negative marks → not weak", {
  rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
  skip_if(rdata_path == "", "Pre-computed test data not available")
  tmp <- new.env()
  load(rdata_path, envir = tmp)
  res <- tmp[[ls(tmp)[1]]]

  first_anchor <- res$loop_annotation[1, ]
  h3k4me1 <- tempfile(fileext = ".bed")
  writeLines(sprintf("%s\t%d\t%d", first_anchor$chr1, first_anchor$start1,
                     first_anchor$end1), h3k4me1)

  val <- validate_epeG_by_chromatin(res, chromatin_beds = list(
    H3K4me1 = h3k4me1
  ), quiet = TRUE)
  skip_if(nrow(val) == 0, "No P/G anchors in pre-computed data")
  expect_true(all(is.na(val$H3K27me3)))
  expect_true(all(is.na(val$H3K4me3)))
  expect_false(any(val$confidence == "weak"),
    info = "Missing negative marks should not produce 'weak' classification")
  unlink(h3k4me1)
})

# --- .assign_chromatin_confidence: all confidence levels ---
test_that(".assign_chromatin_confidence: gold_standard (all 5 marks aligned)", {
  anchors <- data.frame(anchor_id = "A1", chr = "chr1", start = 1L, end = 100L,
    anchor_type = "eP", anchor_gene = "GENE1", cluster_id = "C1",
    stringsAsFactors = FALSE)
  known_marks <- c("H3K4me1", "H3K27ac", "ATAC", "H3K27me3", "H3K4me3")
  mm <- as.data.frame(matrix(NA, nrow = 1, ncol = 5, dimnames = list(NULL, known_marks)))
  mm$H3K4me1 <- TRUE; mm$H3K27ac <- TRUE; mm$ATAC <- TRUE
  mm$H3K27me3 <- FALSE; mm$H3K4me3 <- FALSE

  res <- looplook:::.assign_chromatin_confidence(anchors, mm, known_marks, known_marks)
  expect_equal(as.character(res$confidence), "gold_standard")
  expect_false(any(is.na(res[, known_marks])))
})

test_that(".assign_chromatin_confidence: high_confidence (H3K4me1 + H3K27ac)", {
  anchors <- data.frame(anchor_id = "A1", chr = "chr1", start = 1L, end = 100L,
    anchor_type = "eP", anchor_gene = "GENE1", cluster_id = "C1",
    stringsAsFactors = FALSE)
  known_marks <- c("H3K4me1", "H3K27ac", "ATAC", "H3K27me3", "H3K4me3")
  mm <- as.data.frame(matrix(NA, nrow = 1, ncol = 5, dimnames = list(NULL, known_marks)))
  mm$H3K4me1 <- TRUE; mm$H3K27ac <- TRUE
  # Only 2 marks provided → all_five = FALSE → gold_standard impossible
  # H3K4me1+ and H3K27ac+ → high_confidence

  res <- looplook:::.assign_chromatin_confidence(anchors, mm,
    c("H3K4me1", "H3K27ac"), known_marks)
  expect_equal(as.character(res$confidence), "high_confidence")
})

test_that(".assign_chromatin_confidence: supported (one positive mark)", {
  anchors <- data.frame(anchor_id = "A1", chr = "chr1", start = 1L, end = 100L,
    anchor_type = "eP", anchor_gene = "GENE1", cluster_id = "C1",
    stringsAsFactors = FALSE)
  known_marks <- c("H3K4me1", "H3K27ac", "ATAC", "H3K27me3", "H3K4me3")
  mm <- as.data.frame(matrix(NA, nrow = 1, ncol = 5, dimnames = list(NULL, known_marks)))
  mm$ATAC <- TRUE  # only ATAC positive, no H3K4me1

  res <- looplook:::.assign_chromatin_confidence(anchors, mm, c("ATAC"), known_marks)
  expect_equal(as.character(res$confidence), "supported")
})

test_that(".assign_chromatin_confidence: weak (negative marks tested, no positives)", {
  anchors <- data.frame(anchor_id = "A1", chr = "chr1", start = 1L, end = 100L,
    anchor_type = "eP", anchor_gene = "GENE1", cluster_id = "C1",
    stringsAsFactors = FALSE)
  known_marks <- c("H3K4me1", "H3K27ac", "ATAC", "H3K27me3", "H3K4me3")
  mm <- as.data.frame(matrix(NA, nrow = 1, ncol = 5, dimnames = list(NULL, known_marks)))
  mm$H3K27me3 <- FALSE; mm$H3K4me3 <- FALSE  # negative marks tested, absent

  res <- looplook:::.assign_chromatin_confidence(anchors, mm,
    c("H3K27me3", "H3K4me3"), known_marks)
  expect_equal(as.character(res$confidence), "weak")
})

test_that(".assign_chromatin_confidence: uncertain (all marks absent, NA cols OK)", {
  anchors <- data.frame(anchor_id = "A1", chr = "chr1", start = 1L, end = 100L,
    anchor_type = "eP", anchor_gene = "GENE1", cluster_id = "C1",
    stringsAsFactors = FALSE)
  known_marks <- c("H3K4me1", "H3K27ac", "ATAC", "H3K27me3", "H3K4me3")
  mm <- as.data.frame(matrix(NA, nrow = 1, ncol = 5, dimnames = list(NULL, known_marks)))

  res <- looplook:::.assign_chromatin_confidence(anchors, mm, character(0), known_marks)
  expect_equal(as.character(res$confidence), "uncertain")
  # All mark columns should be NA
  expect_true(all(is.na(res[, known_marks])))
})

test_that(".assign_chromatin_confidence: NA negative marks do not produce weak", {
  anchors <- data.frame(anchor_id = "A1", chr = "chr1", start = 1L, end = 100L,
    anchor_type = "eP", anchor_gene = "GENE1", cluster_id = "C1",
    stringsAsFactors = FALSE)
  known_marks <- c("H3K4me1", "H3K27ac", "ATAC", "H3K27me3", "H3K4me3")
  mm <- as.data.frame(matrix(NA, nrow = 1, ncol = 5, dimnames = list(NULL, known_marks)))
  mm$H3K4me1 <- FALSE; mm$H3K27ac <- FALSE; mm$ATAC <- FALSE  # all active marks tested & absent
  # negative marks are NA → should NOT trigger "weak"

  res <- looplook:::.assign_chromatin_confidence(anchors, mm,
    c("H3K4me1", "H3K27ac", "ATAC"), known_marks)
  expect_equal(as.character(res$confidence), "uncertain")
  expect_true(all(is.na(res$H3K27me3)))
  expect_true(all(is.na(res$H3K4me3)))
})

test_that(".assign_chromatin_confidence: gold_standard fails with missing mark", {
  anchors <- data.frame(anchor_id = "A1", chr = "chr1", start = 1L, end = 100L,
    anchor_type = "eP", anchor_gene = "GENE1", cluster_id = "C1",
    stringsAsFactors = FALSE)
  known_marks <- c("H3K4me1", "H3K27ac", "ATAC", "H3K27me3", "H3K4me3")
  mm <- as.data.frame(matrix(NA, nrow = 1, ncol = 5, dimnames = list(NULL, known_marks)))
  mm$H3K4me1 <- TRUE; mm$H3K27ac <- TRUE; mm$ATAC <- TRUE
  mm$H3K27me3 <- FALSE; mm$H3K4me3 <- NA  # H3K4me3 NOT tested → all_five fails

  res <- looplook:::.assign_chromatin_confidence(anchors, mm,
    c("H3K4me1", "H3K27ac", "ATAC", "H3K27me3"), known_marks)
  # Not gold_standard (only 4 marks provided)
  expect_equal(as.character(res$confidence), "high_confidence")
})

# --- .record_database_versions ---
test_that(".record_database_versions returns expected structure", {
  dbv <- looplook:::.record_database_versions("hg38")
  expect_type(dbv, "list")
  expect_true(all(c("TxDb","OrgDb","BSgenome","JASPAR","clusterProfiler",
    "txdb_pkg","orgdb_pkg") %in% names(dbv)))
})

test_that(".record_database_versions handles NULL species", {
  dbv <- looplook:::.record_database_versions(NULL)
  expect_type(dbv, "list")
  expect_true(all(is.na(dbv$TxDb) & is.na(dbv$OrgDb) & is.na(dbv$BSgenome)))
})
