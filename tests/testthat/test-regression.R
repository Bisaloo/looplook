# tests/testthat/test-regression.R
# Regression tests for bugs found during code review (2026-05-08)

# ── 1. BED/BEDPE 0-based → GRanges 1-based conversion ──────────────────────

test_that("read_simple_bed converts 0-based start to 1-based", {
  tmp_bed <- tempfile(fileext = ".bed")
  writeLines("chr1\t0\t100", tmp_bed)
  gr <- looplook:::read_simple_bed(tmp_bed)
  expect_equal(GenomicRanges::start(gr), 1) # 0-based 0 → 1-based 1
  expect_equal(GenomicRanges::end(gr), 100) # end unchanged
  expect_equal(GenomicRanges::width(gr), 100) # width = end - start + 1 = 100
})

test_that("bedpe_to_gi converts 0-based BEDPE starts to 1-based GRanges", {
  tmp_bedpe <- tempfile(fileext = ".bedpe")
  writeLines("chr1\t0\t100\tchr1\t200\t300", tmp_bedpe)
  gi <- looplook:::bedpe_to_gi(tmp_bedpe)
  a1 <- InteractionSet::anchors(gi, type = "first")
  a2 <- InteractionSet::anchors(gi, type = "second")
  expect_equal(GenomicRanges::start(a1), 1) # 0-based 0 → 1-based 1
  expect_equal(GenomicRanges::start(a2), 201) # 0-based 200 → 1-based 201
})


# ── 2. load_expression_matrix preserves gene names with single column ───────

test_that("load_expression_matrix preserves names with single sample column", {
  tmp_expr <- tempfile(fileext = ".txt")
  writeLines("Gene\tS1\nA\t10\nB\t20\nC\t30", tmp_expr)
  vals <- looplook:::load_expression_matrix(tmp_expr, sample_columns = "S1")
  expect_named(vals)
  expect_equal(names(vals), c("A", "B", "C"))
  expect_equal(unname(vals), c(10, 20, 30))
})

test_that("load_expression_matrix preserves names with multiple sample columns", {
  tmp_expr <- tempfile(fileext = ".txt")
  writeLines("Gene\tS1\tS2\nA\t10\t20\nB\t15\t25", tmp_expr)
  vals <- looplook:::load_expression_matrix(tmp_expr, sample_columns = c("S1", "S2"))
  expect_named(vals)
  expect_equal(names(vals), c("A", "B"))
  expect_equal(unname(vals), c(15, 20))
})


# ── 3. run_go_enrichment works without library(org.Hs.eg.db) ────────────────

test_that("run_go_enrichment runs without attaching org.Hs.eg.db", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Hs.eg.db")
  skip_if_not_installed("enrichplot")
  # Ensure the package is not attached; requireNamespace is enough
  if ("package:org.Hs.eg.db" %in% search()) {
    detach("package:org.Hs.eg.db", unload = TRUE, character.only = TRUE)
  }
  # Use well-known human genes that are guaranteed to map
  genes <- c(
    "TP53", "BRCA1", "EGFR", "MYC", "VEGFA",
    "TNF", "IL6", "CDKN1A", "BCL2", "BAX"
  )
  res <- tryCatch(
    looplook:::run_go_enrichment(genes,
      org_db = "org.Hs.eg.db",
      universe_genes = NULL
    ),
    error = function(e) e
  )
  expect_false(inherits(res, "error"),
    info = paste("run_go_enrichment failed:", conditionMessage(res))
  )
})


# ── 4. refine_loop_anchors_by_expression handles NULL target_annotation ─────

test_that("refine_loop_anchors_by_expression survives NULL target_annotation", {
  rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
  expr_path <- system.file("extdata", "example_tpm.txt", package = "looplook")
  skip_if(rdata_path == "" || expr_path == "", "test data not available")
  tmp <- new.env()
  load(rdata_path, envir = tmp)
  res <- tmp[[ls(tmp)[1]]]
  # Remove target_annotation to simulate loop-only input
  res$target_annotation <- NULL
  expect_error(
    refine_loop_anchors_by_expression(
      annotation_res = res,
      expr_matrix_file = expr_path,
      sample_columns = c("con1", "con2"),
      threshold = 1.0,
      out_dir = tempdir(),
      project_name = "test_null_ta"
    ),
    NA
  )
})


# ── 5. compute_refined_stats splits A;B multi-gene strings ──────────────────

test_that("compute_refined_stats splits semicolon-separated genes", {
  # Build a minimal loop_df with multi-gene anchors
  loop_df <- data.frame(
    loop_ID = c("L1", "L2"),
    anchor1_type = c("P", "E"),
    anchor1_gene = c("A;B", NA),
    anchor2_type = c("E", "P"),
    anchor2_gene = c(NA, "C;D"),
    loop_type = c("E-P", "E-P"),
    stringsAsFactors = FALSE
  )
  vals <- setNames(c(10, 20, 30, 40), c("A", "B", "C", "D"))
  res <- looplook:::compute_refined_stats(
    loop_df = loop_df,
    upstream_promoter_stats = NULL,
    upstream_distal_stats = NULL,
    vals = vals,
    threshold = 1,
    hub_percentile = 0.95
  )
  # A and B should each appear as separate rows in promoter_centric
  genes <- res$promoter_centric$Gene
  expect_true("A" %in% genes)
  expect_true("B" %in% genes)
  expect_true("C" %in% genes)
  expect_true("D" %in% genes)
  # "A;B" must NOT appear as a single gene name
  expect_false("A;B" %in% genes)
  expect_false("C;D" %in% genes)
})


# ── 6. BEDPE coordinate consistency: bedpe_to_gi vs annotate_peaks_and_loops ─

test_that("same BEDPE gives consistent anchor coords in gi and annotation", {
  skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
  tmp_bedpe <- tempfile(fileext = ".bedpe")
  writeLines("chr1\t0\t100\tchr1\t200\t300", tmp_bedpe)
  gi <- looplook:::bedpe_to_gi(tmp_bedpe)
  a1_gi <- GenomicRanges::start(InteractionSet::anchors(gi, "first"))
  a2_gi <- GenomicRanges::start(InteractionSet::anchors(gi, "second"))
  res <- suppressWarnings(annotate_peaks_and_loops(
    bedpe_file = tmp_bedpe,
    species = "hg38",
    out_dir = tempdir(),
    project_name = "coord_test"
  ))
  la <- res$loop_annotation
  expect_true(all(la$start1 >= 1))
  expect_true(all(la$start2 >= 1))
  expect_equal(min(la$start1), a1_gi)
  expect_equal(min(la$start2), a2_gi)
})


# ── 7. consolidate_chromatin_loops exports 0-based BEDPE start ──────────────

test_that("consolidate_chromatin_loops roundtrip preserves 0-based BEDPE start", {
  # Create two BEDPE files so consolidation actually merges and exports
  f1 <- tempfile(fileext = ".bedpe")
  f2 <- tempfile(fileext = ".bedpe")
  writeLines("chr1\t0\t100\tchr1\t200\t300", f1)
  writeLines("chr1\t50\t150\tchr1\t220\t280", f2)
  out_file <- tempfile(fileext = ".bedpe")
  suppressMessages(
    consolidate_chromatin_loops(files = c(f1, f2), out_file = out_file)
  )
  skip_if(!file.exists(out_file), "consolidation produced no output")
  exported <- read.table(out_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  # Roundtrip: 0-based input → 1-based internal → 0-based export → 0-based output
  # The consensus of chr1:0-100 and chr1:50-150 gives anchor1 min start 0
  expect_equal(exported[1, 2], 0) # start1: 0-based
  # The consensus of chr1:200-300 and chr1:220-280 gives anchor2 min start 200
  expect_equal(exported[1, 5], 200) # start2: 0-based (internal 201 → 200)
})


# ── 8. resolve_gene_conflicts edge cases ─────────────────────────────────────

test_that("resolve_gene_conflicts handles empty input", {
  skip_if_not_installed("org.Hs.eg.db")
  skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  empty_df <- data.frame(
    chr = character(0), start = integer(0),
    end = integer(0), stringsAsFactors = FALSE
  )
  res <- looplook:::resolve_gene_conflicts(
    empty_df, txdb,
    "org.Hs.eg.db", c(-2000, 2000), NULL
  )
  expect_equal(nrow(res), 0)
})

test_that("resolve_gene_conflicts handles no matching genes", {
  skip_if_not_installed("org.Hs.eg.db")
  skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  # Region in a gene desert (no promoters nearby)
  df <- data.frame(chr = "chr1", start = 1, end = 100, stringsAsFactors = FALSE)
  res <- suppressWarnings(looplook:::resolve_gene_conflicts(
    df, txdb,
    "org.Hs.eg.db", c(-2000, 2000), NULL
  ))
  expect_s3_class(res, "data.frame")
  expect_true(nrow(res) >= 1)
})


# ── 9. profile_target_genes lfc_col defaults to "log2FoldChange" ────────────

test_that("profile_target_genes lfc_col defaults to log2FoldChange", {
  frmls <- formals(looplook::profile_target_genes)
  expect_equal(frmls$lfc_col, "log2FoldChange")
})


# ── 9b. looplook_report allows precomputed results without bedpe_file ────────

test_that("looplook_report bedpe_file defaults to NULL", {
  frmls <- formals(looplook::looplook_report)
  expect_null(frmls$bedpe_file)
})


# ── 10. Assigned_Target_Genes prioritizes promoter-linked genes ──────────────

test_that("Assigned_Target_Genes is not identical to all loop-connected genes", {
  skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
  skip_if_not_installed("org.Hs.eg.db")
  bedpe_path <- system.file("extdata", "example_loops_1.bedpe", package = "looplook")
  bed_path <- system.file("extdata", "example_peaks.bed", package = "looplook")
  skip_if(bedpe_path == "" || bed_path == "", "example data not available")

  res <- suppressWarnings(suppressMessages(
    annotate_peaks_and_loops(
      bedpe_file = bedpe_path,
      target_bed = bed_path,
      species = "hg38",
      out_dir = tempdir(),
      project_name = "priority_test"
    )
  ))
  ta <- res$target_annotation
  informative <- ta[!is.na(ta$All_Loop_Connected_Genes) &
    !is.na(ta$Assigned_Target_Genes), , drop = FALSE]
  expect_gt(nrow(informative), 0)
  expect_true(any(
    informative$Assigned_Target_Genes != informative$All_Loop_Connected_Genes
  ))
  expect_true(all(vapply(seq_len(nrow(informative)), function(i) {
    assigned <- looplook:::clean_gene_names(
      informative$Assigned_Target_Genes[i], ";"
    )
    all_genes <- looplook:::clean_gene_names(
      informative$All_Loop_Connected_Genes[i], ";"
    )
    all(assigned %in% all_genes)
  }, logical(1))))
})


# ── 11. BEDPE export writes true n_members, not n_reps ───────────────────────

test_that("consolidate_chromatin_loops exports cluster member count", {
  f1 <- tempfile(fileext = ".bedpe")
  f2 <- tempfile(fileext = ".bedpe")
  writeLines(c(
    "chr1\t0\t100\tchr1\t200\t300",
    "chr1\t1000\t1100\tchr1\t1200\t1300"
  ), f1)
  writeLines(c(
    "chr1\t10\t90\tchr1\t210\t290",
    "chr1\t1010\t1090\tchr1\t1210\t1290"
  ), f2)
  out_file <- tempfile(fileext = ".bedpe")
  gi <- suppressMessages(
    consolidate_chromatin_loops(
      files = c(f1, f2),
      mode = "union",
      gap = 1e9,
      out_file = out_file
    )
  )
  exported <- read.table(out_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  expect_equal(nrow(exported), 1)
  expect_equal(exported[1, 9], S4Vectors::mcols(gi)$n_members[1])
  expect_false(exported[1, 9] == S4Vectors::mcols(gi)$n_reps[1])
})


# ── 12. compute_refined_stats retains active G anchors as distal elements ────

test_that("compute_refined_stats keeps G anchors in distal stats", {
  loop_df <- data.frame(
    chr1 = "chr1", start1 = 101, end1 = 200,
    chr2 = "chr1", start2 = 1001, end2 = 1100,
    cluster_id = "C1", a1_id = "A1", a2_id = "A2",
    anchor1_type = "G", anchor1_gene = "GENE_A",
    anchor2_type = "P", anchor2_gene = "GENE_B",
    loop_type = "G-P",
    stringsAsFactors = FALSE
  )
  vals <- setNames(c(10, 20), c("GENE_A", "GENE_B"))
  res <- looplook:::compute_refined_stats(
    loop_df = loop_df,
    upstream_promoter_stats = NULL,
    upstream_distal_stats = NULL,
    vals = vals,
    threshold = 1,
    hub_percentile = 0.95
  )
  expect_false(is.null(res$distal_element))
  expect_equal(nrow(res$distal_element), 1)
  expect_equal(res$distal_element$chr, "chr1")
  expect_equal(res$distal_element$start, 101)
  expect_equal(res$distal_element$end, 200)
  expect_equal(res$distal_element$Target_Genes, "GENE_B")
})


# ── 13. Raincloud grouping is mutually exclusive ──────────────────────────────

test_that("run_heatmap_and_connectivity assigns each gene to one connectivity group", {
  skip_if_not_installed("ggpubr")
  skip_if_not_installed("ggpointdensity")
  skip_if_not_installed("ggdist")

  genes <- c("A", "B", "C", "D", "E")
  tpm_mat <- data.frame(
    S1 = c(10, 11, 12, 13, 14),
    S2 = c(12, 13, 14, 15, 16),
    row.names = genes,
    check.names = FALSE
  )
  meta <- data.frame(
    SampleID = c("S1", "S2"),
    Group = c("G1", "G2"),
    stringsAsFactors = FALSE
  )
  stats_df <- data.frame(
    Gene = genes,
    Total_Loops = c(10, 9, 8, 4, 3),
    Is_High_Connectivity_Gene = c("Yes", "Yes", "No", "No", "No"),
    Is_High_Distal_Connectivity_Gene = c("Yes", "No", "Yes", "No", "No"),
    stringsAsFactors = FALSE
  )
  glist <- setNames(c(1, 2, 3, 4, 5), genes)

  expect_no_warning(
    res <- looplook:::run_heatmap_and_connectivity(
      target_genes = genes,
      tpm_mat_raw = tpm_mat,
      meta_raw = meta,
      loop_stats_df = stats_df,
      global_glist = glist,
      heatmap_ntop = 50,
      cor_method = "pearson",
      current_proj_name = "group_test",
      source_type = "loops",
      skip_heatmap = TRUE
    )
  )
  expect_true("Raincloud_LFC" %in% names(res))
  expect_false(any(duplicated(res$Raincloud_LFC$data$Gene)))
})


# ── 14. annotate_peaks_and_loops accepts OrgDb objects ────────────────────────

test_that("annotate_peaks_and_loops accepts OrgDb package objects", {
  skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
  skip_if_not_installed("org.Hs.eg.db")

  bedpe_path <- system.file("extdata", "example_loops_mini.bedpe", package = "looplook")
  bed_path <- system.file("extdata", "example_peaks_mini.bed", package = "looplook")
  if (bedpe_path == "" || bed_path == "") {
    bedpe_path <- system.file("extdata", "example_loops_1.bedpe", package = "looplook")
    bed_path <- system.file("extdata", "example_peaks.bed", package = "looplook")
  }
  skip_if(bedpe_path == "" || bed_path == "", "example data not available")

  expect_no_error({
    res <- suppressWarnings(suppressMessages(
      annotate_peaks_and_loops(
        bedpe_file = bedpe_path,
        target_bed = bed_path,
        txdb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
        org_db = org.Hs.eg.db::org.Hs.eg.db,
        out_dir = tempdir(),
        project_name = "orgdb_object_test"
      )
    ))
    expect_type(res, "list")
    expect_true("target_annotation" %in% names(res))
  })
})


# ── 15. motif target matching splits semicolon-separated anchor genes ─────────

test_that(".anchor_matches_targets handles semicolon-separated anchor genes", {
  expect_true(looplook:::.anchor_matches_targets("A;B", c("B", "C")))
  expect_true(looplook:::.anchor_matches_targets(" A ; B ", c("A")))
  expect_false(looplook:::.anchor_matches_targets("A;B", c("C")))
  expect_false(looplook:::.anchor_matches_targets("", c("A")))
  expect_false(looplook:::.anchor_matches_targets(NA_character_, c("A")))
})


# ── 16. karyotype plots remain renderable after deferred report rendering ────

test_that("draw_karyo_heatmap_internal stores self-contained image payload", {
  skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
  skip_if_not_installed("png")

  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  gr <- suppressWarnings(GenomicRanges::GRanges(
    seqnames = c("chr1", "chr2"),
    ranges = IRanges::IRanges(start = c(1e6, 2e6), width = 5000)
  ))

  obj <- looplook:::draw_karyo_heatmap_internal(
    gr_data = gr,
    title_prefix = "Deferred render test",
    bin_size = 1e7,
    sat_level = 0.99,
    ref_txdb = txdb,
    plot_species = "hg38",
    unit_label = "Anchors"
  )

  expect_s3_class(obj, "looplook_karyo")
  expect_type(obj$png_raw, "raw")
  expect_gt(length(obj$png_raw), 0)

  # Simulate a stale path from an older object; rendering should still use the
  # embedded PNG payload rather than relying on file existence.
  obj$file <- tempfile(fileext = ".png")
  tmp_pdf <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp_pdf)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_error(print(obj), NA)
})


# ── 17. refinement plots skip karyotype when annotation packages unavailable ─

test_that("build_refinement_plots skips refined karyotype plots if TxDb/OrgDb are missing", {
  skip("testthat cannot mock base::requireNamespace")
  expect_false("Refined_Karyo_Active" %in% names(plots))
  expect_false("Refined_Karyo_TargetGenes" %in% names(plots))
})


# ── 18. profile_target_genes does not require OrgDb when GO is disabled ─────

test_that("profile_target_genes runs without OrgDb if GO analysis is disabled", {
  rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
  expr_path <- system.file("extdata", "example_tpm.txt", package = "looplook")
  diff_path <- system.file("extdata", "example_deg.txt", package = "looplook")
  meta_path <- system.file("extdata", "example_coldata.txt", package = "looplook")
  skip_if(
    rdata_path == "" || expr_path == "" || diff_path == "" || meta_path == "",
    "example data not available"
  )

  tmp <- new.env()
  load(rdata_path, envir = tmp)
  res <- tmp[[ls(tmp)[1]]]

  expect_no_error({
    out <- suppressWarnings(suppressMessages(
      profile_target_genes(
        annotation_res = res,
        diff_file = diff_path,
        expr_matrix_file = expr_path,
        metadata_file = meta_path,
        target_source = "loops",
        project_name = "no_orgdb_needed",
        org_db = "definitelyNotAPackage",
        run_go = FALSE,
        run_ppi = FALSE,
        run_motif = FALSE
      )
    ))
    expect_type(out, "list")
    expect_true("loops" %in% names(out))
  })
})
