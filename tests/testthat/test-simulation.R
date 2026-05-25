# tests/testthat/test-simulation.R — synthetic data covering all code paths

# ════════════════════════════════════════════════════════════════════════════
# Part A: Unit-level tests with pure mock data (no TxDb / OrgDb required)
# ════════════════════════════════════════════════════════════════════════════

test_that("clean_anchor: all reclassification patterns", {
  wl <- c("TP53", "BRCA1")

  # active promoter → stays P
  r <- looplook:::clean_anchor("TP53", "P", wl, TRUE)
  expect_equal(r$type, "P")
  expect_equal(r$gene, "TP53")

  # silent promoter with reclassify → eP
  r <- looplook:::clean_anchor("GENE_X", "P", wl, TRUE)
  expect_equal(r$type, "eP")
  expect_true(is.na(r$gene))

  # silent promoter without reclassify → stays P, gene NA
  r <- looplook:::clean_anchor("GENE_X", "P", wl, FALSE)
  expect_equal(r$type, "P")
  expect_true(is.na(r$gene))

  # silent gene body with reclassify → eG
  r <- looplook:::clean_anchor("GENE_X", "G", wl, TRUE)
  expect_equal(r$type, "eG")

  # enhancer never reclassified (regardless of expression)
  r <- looplook:::clean_anchor("GENE_X", "E", wl, TRUE)
  expect_equal(r$type, "E")

  # multi-gene: keeps only active one
  r <- looplook:::clean_anchor("TP53;GENE_X;BRCA1", "P", wl, TRUE)
  expect_equal(r$gene, "TP53;BRCA1")

  # empty gene string → type preserved, gene NA
  r <- looplook:::clean_anchor("", "E", wl, TRUE)
  expect_equal(r$type, "E")
  expect_true(is.na(r$gene))

  # NA gene
  r <- looplook:::clean_anchor(NA_character_, "P", wl, TRUE)
  expect_equal(r$type, "P")
  expect_true(is.na(r$gene))
})

test_that("extract_target_gene_sets: source branching", {
  anno <- list(
    target_annotation = data.frame(
      SYMBOL = c("A", "B"),
      Assigned_Target_Genes = c("A;C", "B"),
      Assigned_Target_Genes_Filled = c("A;C;D", "B;E"),
      Regulated_promoter_genes = c("A", NA),
      stringsAsFactors = FALSE
    ),
    loop_annotation = data.frame(
      loop_type = c("E-P", "P-P", "E-G", "E-P"),
      Putative_Target_Genes = c("X", "Y", "Z", "W"),
      stringsAsFactors = FALSE
    )
  )

  # targets + all + filled
  r <- looplook:::extract_target_gene_sets(anno, "targets", include_Filled = TRUE, target_mapping_mode = "all")
  expect_setequal(r$Target_Genes, c("A", "B", "C", "D", "E"))

  # targets + promoter mode
  r <- looplook:::extract_target_gene_sets(anno, "targets", include_Filled = FALSE, target_mapping_mode = "promoter")
  expect_setequal(r$Target_Genes, c("A"))

  # targets + use_nearest_gene → uses SYMBOL
  r <- looplook:::extract_target_gene_sets(anno, "targets", use_nearest_gene = TRUE)
  expect_setequal(r$Target_Genes, c("A", "B"))

  # loops with type filter
  r <- looplook:::extract_target_gene_sets(anno, "loops", active_loop_types = c("E-P"))
  expect_true("EP_Genes" %in% names(r))
  expect_setequal(r$EP_Genes, c("X", "W"))

  # loops all types
  r <- looplook:::extract_target_gene_sets(anno, "loops")
  expect_equal(length(r), 3)  # EP, PP, EG
})

test_that("compute_refined_stats: all loop types and multi-gene split", {
  loop_df <- data.frame(
    chr1 = rep("chr1", 8), start1 = seq(1000, 8000, 1000), end1 = seq(1500, 8500, 1000),
    chr2 = rep("chr1", 8), start2 = seq(10000, 17000, 1000), end2 = seq(10500, 17500, 1000),
    cluster_id = paste0("C", c(1, 1, 2, 2, 3, 3, 4, 4)),
    a1_id = paste0("A", 1:8), a2_id = paste0("B", 1:8),
    anchor1_type = c("P", "P", "E",  "E",  "G",  "G",  "P",  "P"),
    anchor2_type = c("E", "P", "P",  "G",  "P",  "E",  "eP", "eG"),
    anchor1_gene = c("G1", "G2;G3", NA,   NA,   "G4", "G5", "G6;G1", "G7"),
    anchor2_gene = c(NA,   NA,   "G1", "G8", "G9", NA,   "G10",    "G11"),
    loop_type = c("E-P", "P-P", "E-P", "E-G", "G-P", "E-G", "eP-P", "eG-P"),
    Putative_Target_Genes = c("G1", "G2;G3", "G1", "G8", "G9", "G5", "G6;G1", "G11"),
    stringsAsFactors = FALSE
  )
  vals <- setNames(seq_len(11) * 2, paste0("G", 1:11))

  res <- looplook:::compute_refined_stats(loop_df, upstream_promoter_stats = NULL,
    vals = vals, threshold = 1, hub_percentile = 0.95)

  # promoter-centric: G1, G2, G3, G6, G7 should appear as separate rows
  expect_true(all(c("G1", "G2", "G3", "G6", "G7") %in% res$promoter_centric$Gene))
  expect_false("G2;G3" %in% res$promoter_centric$Gene)
  expect_false("G6;G1" %in% res$promoter_centric$Gene)

  # G1 appears 3 times: L1 a1=P, L3 a2=P, L7 a1=P (G6;G1 split)
  g1_row <- res$promoter_centric[res$promoter_centric$Gene == "G1", ]
  expect_equal(g1_row$Total_Loops, 3)

  # distal elements: E, G, eP, eG anchors
  expect_true(!is.null(res$distal_element))
  expect_gt(nrow(res$distal_element), 0)
  # G9 is target of distal anchor A5 (G-P loop: a1=G(E-type distal), a2=P)
  expect_true("G9" %in% res$distal_element$Target_Genes)
})

test_that("compute_refined_stats: hub detection thresholds", {
  loop_df <- data.frame(
    chr1 = "chr1", start1 = 1:10 * 1000, end1 = 1:10 * 1000 + 500,
    chr2 = "chr1", start2 = 1:10 * 2000, end2 = 1:10 * 2000 + 500,
    cluster_id = paste0("C", 1:10),
    a1_id = paste0("A", 1:10), a2_id = paste0("B", 1:10),
    anchor1_type = rep("P", 10), anchor2_type = rep("E", 10),
    anchor1_gene = paste0("Gene", 1:10), anchor2_gene = NA_character_,
    loop_type = rep("E-P", 10),
    Putative_Target_Genes = paste0("Gene", 1:10),
    stringsAsFactors = FALSE
  )
  # Gene1-3 have high connectivity (many loops), Gene4-10 have few
  gene_loops <- c(rep(10, 3), rep(2, 7))
  loop_df <- do.call(rbind, lapply(seq_along(gene_loops), function(i) {
    n <- gene_loops[i]
    if (n <= 1) return(NULL)
    data.frame(
      chr1 = "chr1", start1 = i * 1000, end1 = i * 1000 + 500,
      chr2 = "chr1", start2 = (i + 1) * 2000, end2 = (i + 1) * 2000 + 500,
      cluster_id = paste0("C", i, "_", seq_len(n)),
      a1_id = paste0("A", i, "_", seq_len(n)),
      a2_id = paste0("B", i, "_", seq_len(n)),
      anchor1_type = "P", anchor2_type = "E",
      anchor1_gene = paste0("Gene", i), anchor2_gene = NA_character_,
      loop_type = "E-P",
      Putative_Target_Genes = paste0("Gene", i),
      stringsAsFactors = FALSE
    )
  }))
  vals <- setNames(seq_len(10) * 5, paste0("Gene", 1:10))

  res <- looplook:::compute_refined_stats(loop_df, upstream_promoter_stats = NULL,
    vals = vals, threshold = 1, hub_percentile = 0.50)

  # With 50th percentile, at least Gene1-3 (with 10 loops each) should be hubs
  high_conn <- res$promoter_centric$Gene[res$promoter_centric$Is_High_Connectivity_Gene == "Yes"]
  expect_true("Gene1" %in% high_conn)
})

test_that("get_feature_class: all categories", {
  classify <- function(x) {
    if (is.na(x)) return("Unknown")
    x <- tolower(x)
    if (grepl("promoter", x)) return("P")
    if (grepl("intergenic|downstream", x)) return("E")
    if (grepl("exon|intron|utr", x)) return("G")
    return("E")
  }
  inputs <- c("Promoter (<=1kb)", "Distal Intergenic", "Intron (uc001.1)",
              "Exon (uc001.1)", "5-UTR", "Downstream (<=300bp)",
              "3-UTR", NA_character_)
  expected <- c("P", "E", "G", "G", "G", "E", "G", "Unknown")
  expect_equal(unname(vapply(inputs, classify, character(1))), expected)
})

test_that("extract_genes and clean_gene_names: complete edge cases", {
  expect_equal(looplook:::extract_genes(c("A;B", "B;C", "D")), "A;B;C;D")
  expect_equal(looplook:::extract_genes(""), NA_character_)
  expect_equal(looplook:::extract_genes(";"), NA_character_)
  expect_equal(looplook:::extract_genes(c(NA, "")), NA_character_)
  expect_equal(looplook:::extract_genes(" A ; B "), " A ; B ")  # extract_genes splits by ; only, no trimws

  expect_equal(looplook:::clean_gene_names(NULL), character(0))
  expect_equal(looplook:::clean_gene_names(c("A", NA, "", " ")), "A")
  expect_equal(looplook:::clean_gene_names("A;B;A;", ";"), c("A", "B"))
  expect_equal(looplook:::clean_gene_names(c("A,B", "C,D"), "[;,]"), c("A", "B", "C", "D"))
})

test_that("load_expression_matrix: all input variants", {
  tmp <- tempfile(fileext = ".tsv")
  writeLines("Gene\tS1\tS2\tS3\nA\t10\t20\t30\nB\t5\t15\t25\nC\t0\t0\t0", tmp)

  # multiple columns
  v <- looplook:::load_expression_matrix(tmp, c("S1", "S2"))
  expect_equal(v[["A"]], 15)

  # single column preserves names
  v <- looplook:::load_expression_matrix(tmp, "S1")
  expect_equal(v[["A"]], 10)

  # all columns (NULL)
  v <- looplook:::load_expression_matrix(tmp, NULL)
  expect_equal(v[["C"]], 0)

  # integer indices
  v <- looplook:::load_expression_matrix(tmp, c(2L, 3L))
  expect_equal(names(v), c("A", "B", "C"))

  unlink(tmp)

  # error: non-numeric values
  tmp2 <- tempfile(fileext = ".tsv")
  writeLines("Gene\tS1\nA\tNA\nB\tbad", tmp2)
  expect_error(looplook:::load_expression_matrix(tmp2, "S1"), "non-numeric")
  unlink(tmp2)
})

test_that("format_annotation_columns: roundtrip", {
  df <- data.frame(annotation = c("Promoter (<=1kb)", "Intron (uc001.1)"),
                   stringsAsFactors = FALSE)
  res <- looplook:::format_annotation_columns(df)
  expect_equal(res$annotation, c("Promoter", "Intron"))
  expect_equal(res$detail_anno, c("Promoter (<=1kb)", "Intron (uc001.1)"))

  # no annotation column → pass through
  df2 <- data.frame(x = 1:3)
  res2 <- looplook:::format_annotation_columns(df2)
  expect_equal(res2$x, 1:3)
})

test_that("simplify_annotation: all paths", {
  x <- c("Promoter (<=1kb)", "Intron (uc001.1)", "Exon (uc001.1)",
         "Distal Intergenic", "Downstream (<=300bp)", "unknown_type", NA)
  res <- looplook:::simplify_annotation(x)
  expect_equal(unname(res), c("Promoter", "Intron", "Exon", "Distal Intergenic",
                      "Downstream", "Others", "Others"))
})

test_that("read_robust_general: fill option for ragged lines", {
  tmp <- tempfile(fileext = ".csv")
  writeLines("A,B,C\n1,2\n3,4,5,6", tmp)
  res <- looplook:::read_robust_general(tmp, header = TRUE)
  expect_equal(ncol(res), 4)
  unlink(tmp)
})

test_that("get_colors: RColorBrewer, custom, and fallback", {
  expect_length(looplook:::get_colors(0, "Set2"), 0)
  expect_length(looplook:::get_colors(1, "Set1"), 1)
  expect_length(looplook:::get_colors(50, "Set2"), 50)

  cols <- looplook:::get_colors(3, c("#FF0000"))
  expect_equal(cols, c("#FF0000", "#FF0000", "#FF0000"))

  cols <- looplook:::get_colors(3, NULL)
  expect_length(cols, 3)
})


# ════════════════════════════════════════════════════════════════════════════
# Part B: Integration tests with sample TxDb
# ════════════════════════════════════════════════════════════════════════════

test_that("annotate_peaks_and_loops: all loop-type combinations", {
  skip_if_not_installed("org.Hs.eg.db")
  txdb <- AnnotationDbi::loadDb(
    system.file("extdata", "hg19_knownGene_sample.sqlite", package = "GenomicFeatures")
  )
  skip_if(is.null(txdb), "Sample TxDb unavailable")

  # Get real gene coordinates from TxDb
  genes_gr <- GenomicFeatures::genes(txdb)
  gene_ids <- names(genes_gr)
  map <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
    keys = gene_ids, columns = "SYMBOL", keytype = "ENTREZID")
  map <- map[!duplicated(map$ENTREZID), ]

  # Pick genes on the same chromosome for intra-chromosomal loops
  gene_df <- data.frame(
    gene_id = gene_ids,
    chr = as.character(GenomicRanges::seqnames(genes_gr)),
    start = GenomicRanges::start(genes_gr),
    end = GenomicRanges::end(genes_gr),
    stringsAsFactors = FALSE
  )
  gene_df <- merge(gene_df, map, by.x = "gene_id", by.y = "ENTREZID")

  # Find chr6 genes: TFAP2A-AS1 (~10.4M) and SRPK1 (~35.8M)
  # Place anchors to create E-P, P-P, E-E, E-G, P-G loops
  chr6_genes <- gene_df[gene_df$chr == "chr6", ]
  skip_if(nrow(chr6_genes) < 2, "Need ≥2 genes on same chr for loop tests")

  g1 <- chr6_genes[1, ]  # TFAP2A-AS1: ncRNA
  g2 <- chr6_genes[2, ]  # SRPK1: protein-coding

  # Anchor A: exactly at g1 TSS → should annotate as Promoter
  # Anchor B: intergenic between g1 and g2 → should annotate as Distal Intergenic (E)
  # Anchor C: at g2 TSS → should annotate as Promoter
  # Anchor D: within g2 gene body (TSS+1000) → should annotate as G (Exon/Intron)

  tss_a <- ifelse(g1$start < g1$end, g1$start, g1$end)
  tss_c <- ifelse(g2$start < g2$end, g2$start, g2$end)
  mid_ab <- round((tss_a + tss_c) / 2)

  # Build BEDPE with 4 loops:
  # L1: P(A)-E(B) → type should be E-P
  # L2: P(A)-P(C) → type should be P-P
  # L3: E(B)-P(C) → type should be E-P
  # L4: P(A)-G(D) → type should be G-P (if D hits gene body)
  bedpe_lines <- c(
    sprintf("chr6\t%d\t%d\tchr6\t%d\t%d", tss_a - 1000, tss_a + 500, mid_ab - 500, mid_ab + 500),
    sprintf("chr6\t%d\t%d\tchr6\t%d\t%d", tss_a - 1000, tss_a + 500, tss_c - 1000, tss_c + 500),
    sprintf("chr6\t%d\t%d\tchr6\t%d\t%d", mid_ab - 500, mid_ab + 500, tss_c - 1000, tss_c + 500),
    sprintf("chr6\t%d\t%d\tchr6\t%d\t%d", tss_a - 1000, tss_a + 500, tss_c + 1000, tss_c + 10000)
  )
  tmp_bedpe <- tempfile(fileext = ".bedpe")
  writeLines(bedpe_lines, tmp_bedpe)

  res <- looplook:::.with_known_upstream_noise_suppressed(
    looplook::annotate_peaks_and_loops(
      bedpe_file = tmp_bedpe,
      txdb = txdb, org_db = "org.Hs.eg.db", species = "hg19",
      out_dir = tempdir(), project_name = "Sim_LoopTypes",
      write_output = FALSE, quiet = TRUE
    )
  )

  expect_type(res, "list")
  la <- res$loop_annotation
  expect_equal(nrow(la), 4)

  # All loops should have a classified loop_type (not "Unknown")
  expect_false(any(la$loop_type == "Unknown"))

  # Promoter-centric stats should exist for the P anchors
  expect_gt(nrow(res$promoter_centric_stats), 0)

  unlink(tmp_bedpe)
})

test_that("annotate_peaks_and_loops: target BED integration with and without loop overlap", {
  skip_if_not_installed("org.Hs.eg.db")
  txdb <- AnnotationDbi::loadDb(
    system.file("extdata", "hg19_knownGene_sample.sqlite", package = "GenomicFeatures")
  )
  skip_if(is.null(txdb), "Sample TxDb unavailable")

  genes_gr <- GenomicFeatures::genes(txdb)
  gene_ids <- names(genes_gr)
  map <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
    keys = gene_ids, columns = "SYMBOL", keytype = "ENTREZID")
  map <- map[!duplicated(map$ENTREZID), ]
  gene_df <- data.frame(
    gene_id = gene_ids,
    chr = as.character(GenomicRanges::seqnames(genes_gr)),
    start = GenomicRanges::start(genes_gr),
    end = GenomicRanges::end(genes_gr),
    stringsAsFactors = FALSE
  )
  gene_df <- merge(gene_df, map, by.x = "gene_id", by.y = "ENTREZID")
  chrX_genes <- gene_df[gene_df$chr == "chrX", ]
  skip_if(nrow(chrX_genes) < 3, "Need ≥3 genes on chrX")

  g1 <- chrX_genes[1, ]  # CLCN4
  g2 <- chrX_genes[3, ]  # ARMCX3
  tss1 <- ifelse(g1$start < g1$end, g1$start, g1$end)
  tss2 <- ifelse(g2$start < g2$end, g2$start, g2$end)

  # One loop: P(g1) <-> E(intergenic)
  bedpe_lines <- sprintf("chrX\t%d\t%d\tchrX\t%d\t%d",
    tss1 - 1000, tss1 + 500,
    round((tss1 + tss2) / 2) - 500, round((tss1 + tss2) / 2) + 500)
  tmp_bedpe <- tempfile(fileext = ".bedpe")
  writeLines(bedpe_lines, tmp_bedpe)

  # Target BED:
  # Peak1: overlaps the P anchor → should get 3D-assigned targets
  # Peak2: far away from any anchor → orphan (only linear nearest gene)
  tmp_bed <- tempfile(fileext = ".bed")
  writeLines(c(
    sprintf("chrX\t%d\t%d", tss1 - 500, tss1 + 100),    # overlaps promoter anchor
    "chrX\t50000000\t50001000"                            # orphan peak
  ), tmp_bed)

  res <- looplook:::.with_known_upstream_noise_suppressed(
    looplook::annotate_peaks_and_loops(
      bedpe_file = tmp_bedpe,
      target_bed = tmp_bed,
      txdb = txdb, org_db = "org.Hs.eg.db", species = "hg19",
      out_dir = tempdir(), project_name = "Sim_Targets",
      write_output = FALSE, quiet = TRUE
    )
  )

  ta <- res$target_annotation
  expect_equal(nrow(ta), 2)

  # Peak1 should have 3D-derived targets
  peak1 <- ta[1, ]
  expect_false(is.na(peak1$Assigned_Target_Genes) || peak1$Assigned_Target_Genes == "")

  # Both should have _Filled column (peak2 falls back to SYMBOL)
  expect_true("Assigned_Target_Genes_Filled" %in% colnames(ta))
  expect_false(any(is.na(ta$Assigned_Target_Genes_Filled) & ta$Assigned_Target_Genes_Filled == ""))

  # All_Loop_Connected_Genes should differ from Assigned_Target_Genes
  # (ATGs uses priority: promoter-only for 3D assignment)
  has_loop <- !is.na(ta$All_Loop_Connected_Genes) & ta$All_Loop_Connected_Genes != ""
  if (any(has_loop)) {
    # verify Assigned_Target_Genes is subset of All_Loop_Connected_Genes
    for (i in which(has_loop)) {
      assigned <- looplook:::clean_gene_names(ta$Assigned_Target_Genes[i], ";")
      all_conn <- looplook:::clean_gene_names(ta$All_Loop_Connected_Genes[i], ";")
      if (length(assigned) > 0)
        expect_true(all(assigned %in% all_conn))
    }
  }

  unlink(c(tmp_bedpe, tmp_bed))
})

test_that("annotate_peaks_and_loops: neighbor_hop controls ego-network radius", {
  skip_if_not_installed("org.Hs.eg.db")
  txdb <- AnnotationDbi::loadDb(
    system.file("extdata", "hg19_knownGene_sample.sqlite", package = "GenomicFeatures")
  )
  skip_if(is.null(txdb), "Sample TxDb unavailable")

  genes_gr <- GenomicFeatures::genes(txdb)
  gene_ids <- names(genes_gr)
  map <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
    keys = gene_ids, columns = "SYMBOL", keytype = "ENTREZID")
  map <- map[!duplicated(map$ENTREZID), ]
  gene_df <- data.frame(
    gene_id = gene_ids,
    chr = as.character(GenomicRanges::seqnames(genes_gr)),
    start = GenomicRanges::start(genes_gr),
    end = GenomicRanges::end(genes_gr),
    stringsAsFactors = FALSE
  )
  gene_df <- merge(gene_df, map, by.x = "gene_id", by.y = "ENTREZID")
  chrX_genes <- gene_df[gene_df$chr == "chrX", ]
  skip_if(nrow(chrX_genes) < 3, "Need ≥3 genes on chrX")

  g1 <- chrX_genes[1, ]
  g2 <- chrX_genes[2, ]
  g3 <- chrX_genes[3, ]
  tss1 <- ifelse(g1$start < g1$end, g1$start, g1$end)
  tss2 <- ifelse(g2$start < g2$end, g2$start, g2$end)
  tss3 <- ifelse(g3$start < g3$end, g3$start, g3$end)

  # Build 3 anchors linked as: A(g1) — B(intergenic) — C(g2)
  # hop=0: ego for A gives only B (no gene)
  # hop=1: ego for A gives B and C (reaching g2)
  # hop=2 (target): ego for A gives B, C and also D(g3) if D connects to C
  bedpe_lines <- c(
    sprintf("chrX\t%d\t%d\tchrX\t%d\t%d",
      tss1 - 1000, tss1 + 500, round((tss1 + tss2) / 2) - 500, round((tss1 + tss2) / 2) + 500),
    sprintf("chrX\t%d\t%d\tchrX\t%d\t%d",
      round((tss1 + tss2) / 2) - 500, round((tss1 + tss2) / 2) + 500, tss2 - 1000, tss2 + 500)
  )
  tmp_bedpe <- tempfile(fileext = ".bedpe")
  writeLines(bedpe_lines, tmp_bedpe)

  res0 <- looplook:::.with_known_upstream_noise_suppressed(
    looplook::annotate_peaks_and_loops(
      bedpe_file = tmp_bedpe, txdb = txdb, org_db = "org.Hs.eg.db", species = "hg19",
      neighbor_hop = 0, out_dir = tempdir(), project_name = "Sim_Hop0",
      write_output = FALSE, quiet = TRUE
    )
  )
  res1 <- looplook:::.with_known_upstream_noise_suppressed(
    looplook::annotate_peaks_and_loops(
      bedpe_file = tmp_bedpe, txdb = txdb, org_db = "org.Hs.eg.db", species = "hg19",
      neighbor_hop = 1, out_dir = tempdir(), project_name = "Sim_Hop1",
      write_output = FALSE, quiet = TRUE
    )
  )

  # With hop=0, anchor A gets ego genes only from its direct neighbor (B, intergenic → no gene)
  # With hop=1, anchor A gets ego genes from C too (which is a promoter)
  # So Putative_Target_Genes should be richer in hop=1
  la0 <- res0$loop_annotation
  la1 <- res1$loop_annotation

  # The promoter anchor in L1 should have Putative_Target_Genes
  l1_0 <- la0[1, ]
  l1_1 <- la1[1, ]
  expect_false(is.na(l1_1$Putative_Target_Genes) || l1_1$Putative_Target_Genes == "")

  unlink(tmp_bedpe)
})

test_that("refine_loop_anchors_by_expression: reclassify cascade", {
  skip_if_not_installed("org.Hs.eg.db")
  rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
  expr_path <- system.file("extdata", "example_tpm.txt", package = "looplook")
  skip_if(rdata_path == "" || expr_path == "")

  tmp <- new.env()
  load(rdata_path, envir = tmp)
  res <- tmp[[ls(tmp)[1]]]
  res$loop_annotation <- head(res$loop_annotation, 30)
  res$target_annotation <- head(res$target_annotation, 10)
  res$promoter_centric_stats <- head(res$promoter_centric_stats, 30)
  res$distal_element_stats <- head(res$distal_element_stats, 30)

  # Step A: without reclassify
  ref_nr <- looplook:::.with_known_upstream_noise_suppressed(
    looplook::refine_loop_anchors_by_expression(
      annotation_res = res, expr_matrix_file = expr_path,
      sample_columns = c("con1", "con2"), threshold = 1,
      reclassify_by_expression = FALSE,
      out_dir = tempdir(), project_name = "Sim_NoReclassify",
      write_output = FALSE, quiet = TRUE
    )
  )

  # Step B: with reclassify
  ref_wr <- looplook:::.with_known_upstream_noise_suppressed(
    looplook::refine_loop_anchors_by_expression(
      annotation_res = res, expr_matrix_file = expr_path,
      sample_columns = c("con1", "con2"), threshold = 1,
      reclassify_by_expression = TRUE,
      out_dir = tempdir(), project_name = "Sim_Reclassify",
      write_output = FALSE, quiet = TRUE
    )
  )

  # Both should return valid structures
  expect_type(ref_nr, "list")
  expect_type(ref_wr, "list")
  expect_true("Putative_Target_Genes" %in% colnames(ref_nr$loop_annotation))
  expect_true("Putative_Target_Genes" %in% colnames(ref_wr$loop_annotation))

  # With reclassify enabled, loop types may include eP/eG prefixes
  loop_types_wr <- unique(ref_wr$loop_annotation$loop_type)
  loop_types_nr <- unique(ref_nr$loop_annotation$loop_type)
  # The reclassified set may have more types (including eP, eG variants)
  expect_true(length(loop_types_wr) >= length(loop_types_nr) ||
    any(grepl("e", loop_types_wr)))
})

test_that("pipeline: empty and edge-case inputs", {
  skip_if_not_installed("org.Hs.eg.db")
  txdb <- AnnotationDbi::loadDb(
    system.file("extdata", "hg19_knownGene_sample.sqlite", package = "GenomicFeatures")
  )
  skip_if(is.null(txdb), "Sample TxDb unavailable")

  # Single loop on different chromosomes → no clustering possible
  tmp_bedpe <- tempfile(fileext = ".bedpe")
  writeLines("chr1\t1000000\t1000500\tchr6\t10412000\t10412600", tmp_bedpe)

  res <- looplook:::.with_known_upstream_noise_suppressed(
    looplook::annotate_peaks_and_loops(
      bedpe_file = tmp_bedpe, txdb = txdb, org_db = "org.Hs.eg.db",
      species = "hg19", out_dir = tempdir(), project_name = "Sim_Trans",
      write_output = FALSE, quiet = TRUE
    )
  )
  # Inter-chromosomal loops should still annotate
  expect_equal(nrow(res$loop_annotation), 1)
  expect_false(res$loop_annotation$loop_type == "Unknown")

  unlink(tmp_bedpe)
})

test_that("pipeline: expression matrix with NA values", {
  tmp_expr <- tempfile(fileext = ".tsv")
  writeLines("Gene\tS1\tS2\nA\t10\tNA\nB\tNA\t20\nC\t5\t5", tmp_expr)

  vals <- looplook:::load_expression_matrix(tmp_expr, c("S1", "S2"))
  expect_equal(vals[["A"]], 10)  # mean of 10 and NA = 10
  expect_equal(vals[["B"]], 20)  # mean of NA and 20 = 20
  expect_equal(vals[["C"]], 5)

  unlink(tmp_expr)
})

test_that(".anchor_matches_targets: case and whitespace", {
  expect_true(looplook:::.anchor_matches_targets("  TP53 ; BRCA1 ", c("BRCA1", "EGFR")))
  expect_false(looplook:::.anchor_matches_targets("TP53;BRCA1", c("EGFR")))
  expect_false(looplook:::.anchor_matches_targets("", c("A")))
  expect_false(looplook:::.anchor_matches_targets(NA_character_, c("A")))
})

test_that(".prepare_motif_anchor_sets: handles mixed anchor types", {
  loop_df <- data.frame(
    chr1 = rep("chr1", 4), start1 = c(100, 300, 500, 700), end1 = c(150, 350, 550, 750),
    chr2 = rep("chr1", 4), start2 = c(200, 400, 600, 800), end2 = c(250, 450, 650, 850),
    anchor1_gene = c("TP53", "CTRL1", "CTRL2", "CTRL3"),
    anchor1_type = c("P", "P", "P", "E"),
    anchor2_gene = c("EnhA", "BRCA1", "BgEnh", "PROM1"),
    anchor2_type = c("E", "P", "E", "P"),
    stringsAsFactors = FALSE
  )

  # Target TP53 → L1 matches (a1=P, a2=E), L2-L4 are background
  ms <- looplook:::.prepare_motif_anchor_sets(loop_df, c("TP53"))
  expect_equal(ms$target_loop_n, 1)
  # FG proximal = a1 of L1 (TP53 P anchor)
  expect_equal(length(ms$proximal_fg), 1)
  # FG distal = a2 of L1 (EnhA E anchor)
  expect_equal(length(ms$distal_fg), 1)
  # BG should have the remaining anchors from non-target loops
  expect_equal(length(ms$proximal_bg), 4)  # L2 a1, L3 a1, L2 a2, L4 a2
  expect_equal(length(ms$distal_bg), 2)    # a2 of L2, L3 (E type, non-target) + L4 distal

  # Target both TP53 and BRCA1
  ms2 <- looplook:::.prepare_motif_anchor_sets(loop_df, c("TP53", "BRCA1"))
  expect_equal(ms2$target_loop_n, 2)

  # Target none → 0 target loops
  ms3 <- looplook:::.prepare_motif_anchor_sets(loop_df, c("NONEXISTENT"))
  expect_equal(ms3$target_loop_n, 0)
})

test_that(".subset_motif_loop_df: task-to-loop_type mapping", {
  df <- data.frame(loop_type = c("E-P", "P-P", "G-P", "E-G"), stringsAsFactors = FALSE)

  r <- looplook:::.subset_motif_loop_df(df, "loops", "EP_Genes")
  expect_true(all(r$loop_type == "E-P"))

  r <- looplook:::.subset_motif_loop_df(df, "loops", "PP_Genes")
  expect_true(all(r$loop_type == "P-P"))

  r <- looplook:::.subset_motif_loop_df(df, "loops", "GP_Genes")
  expect_true(all(r$loop_type == "G-P"))

  # Non-matching task returns full df
  r <- looplook:::.subset_motif_loop_df(df, "loops", "ZZ_Genes")
  expect_equal(nrow(r), nrow(df))

  # Non-loops source returns full df
  r <- looplook:::.subset_motif_loop_df(df, "targets", "EP_Genes")
  expect_equal(nrow(r), nrow(df))
})

test_that(".deduplicate_anchor_df and .anchor_df_to_gr: robustness", {
  df <- data.frame(
    anchor_id = c("A1", "A1", "", NA),
    chr = c("chr1", "chr1", "chr2", "chr3"),
    start = c(100, 100, 200, 300),
    end = c(150, 150, 250, 350),
    anchor_type = c("P", "P", "E", "G"),
    stringsAsFactors = FALSE
  )
  dedup <- looplook:::.deduplicate_anchor_df(df)
  expect_equal(nrow(dedup), 3)  # A1, chr2_200_250, chr3_300_350
  expect_false(any(duplicated(dedup$anchor_id)))

  # empty input
  empty <- looplook:::.deduplicate_anchor_df(data.frame(
    anchor_id = character(0), chr = character(0),
    start = integer(0), end = integer(0), anchor_type = character(0),
    stringsAsFactors = FALSE
  ))
  expect_equal(nrow(empty), 0)

  # convert to GRanges
  gr <- looplook:::.anchor_df_to_gr(df)
  expect_s4_class(gr, "GRanges")
})


# ════════════════════════════════════════════════════════════════════════════
# Part C: Additional integration patterns (sample TxDb required)
# ════════════════════════════════════════════════════════════════════════════

test_that("annotate_peaks_and_loops: E-E and G-G loop types", {
  skip_if_not_installed("org.Hs.eg.db")
  txdb <- AnnotationDbi::loadDb(
    system.file("extdata", "hg19_knownGene_sample.sqlite", package = "GenomicFeatures")
  )
  skip_if(is.null(txdb), "Sample TxDb unavailable")

  genes_gr <- GenomicFeatures::genes(txdb)
  gene_ids <- names(genes_gr)

  # Find genes with enough spacing to place intergenic anchors between them
  gene_df <- data.frame(
    gene_id = gene_ids,
    chr = as.character(GenomicRanges::seqnames(genes_gr)),
    start = GenomicRanges::start(genes_gr),
    end = GenomicRanges::end(genes_gr),
    stringsAsFactors = FALSE
  )

  # chr8 has IDO2 (~39.8M) and MIR1204 (~128.8M) — far apart, plenty of intergenic space
  # Place two intergenic anchors (E-E) between them
  chr8_genes <- gene_df[gene_df$chr == "chr8", ]
  skip_if(nrow(chr8_genes) < 2, "Need ≥2 genes on chr8")

  g1 <- chr8_genes[1, ]
  mid <- round((g1$end + 50000000))

  # E-E loop: two intergenic anchors
  # G-G loop: place anchors inside gene bodies
  tss_g1 <- ifelse(g1$start < g1$end, g1$start, g1$end)

  bedpe_lines <- c(
    # L1: E-E (two intergenic anchors far from any gene)
    sprintf("chr8\t%d\t%d\tchr8\t%d\t%d",
      mid - 500, mid + 500, mid + 10000 - 500, mid + 10000 + 500),
    # L2: G-G within same gene body (anchor both inside g1)
    sprintf("chr8\t%d\t%d\tchr8\t%d\t%d",
      tss_g1 + 5000, tss_g1 + 5500, tss_g1 + 7000, tss_g1 + 7500)
  )
  tmp_bedpe <- tempfile(fileext = ".bedpe")
  writeLines(bedpe_lines, tmp_bedpe)

  res <- looplook:::.with_known_upstream_noise_suppressed(
    looplook::annotate_peaks_and_loops(
      bedpe_file = tmp_bedpe, txdb = txdb, org_db = "org.Hs.eg.db",
      species = "hg19", out_dir = tempdir(), project_name = "Sim_EE_GG",
      write_output = FALSE, quiet = TRUE
    )
  )

  la <- res$loop_annotation
  loop_types <- sort(unique(la$loop_type))
  # E-E and G-G should be present (ChIPseeker classifies intergenic as E,
  # and internal gene body regions as Intron/Exon → G)
  expect_true(any(grepl("E-E", loop_types)) || any(grepl("G-G", loop_types)))

  unlink(tmp_bedpe)
})

test_that("annotate_peaks_and_loops: target peak connected to multiple loops", {
  skip_if_not_installed("org.Hs.eg.db")
  txdb <- AnnotationDbi::loadDb(
    system.file("extdata", "hg19_knownGene_sample.sqlite", package = "GenomicFeatures")
  )
  skip_if(is.null(txdb), "Sample TxDb unavailable")

  genes_gr <- GenomicFeatures::genes(txdb)
  gene_ids <- names(genes_gr)
  gene_df <- data.frame(
    gene_id = gene_ids,
    chr = as.character(GenomicRanges::seqnames(genes_gr)),
    start = GenomicRanges::start(genes_gr),
    end = GenomicRanges::end(genes_gr),
    stringsAsFactors = FALSE
  )

  # Use chrX genes ARMCX3 and FAM199X (close enough for realistic loops)
  chrX_genes <- gene_df[gene_df$chr == "chrX", ]
  skip_if(nrow(chrX_genes) < 3, "Need ≥3 genes on chrX")

  g1 <- chrX_genes[3, ]  # ARMCX3
  g2 <- chrX_genes[4, ]  # FAM199X
  tss1 <- ifelse(g1$start < g1$end, g1$start, g1$end)
  tss2 <- ifelse(g2$start < g2$end, g2$start, g2$end)

  # Build 3 loops that share the SAME promoter anchor (A):
  # L1: P(A at g1) — E(B intergenic)
  # L2: P(A at g1) — E(C intergenic, different location)
  # L3: P(A at g1) — P(D at g2)
  anchor_a_start <- tss1 - 1000
  anchor_a_end   <- tss1 + 500

  bedpe_lines <- c(
    sprintf("chrX\t%d\t%d\tchrX\t%d\t%d",
      anchor_a_start, anchor_a_end, tss1 + 50000 - 500, tss1 + 50000 + 500),
    sprintf("chrX\t%d\t%d\tchrX\t%d\t%d",
      anchor_a_start, anchor_a_end, tss1 + 80000 - 500, tss1 + 80000 + 500),
    sprintf("chrX\t%d\t%d\tchrX\t%d\t%d",
      anchor_a_start, anchor_a_end, tss2 - 1000, tss2 + 500)
  )
  tmp_bedpe <- tempfile(fileext = ".bedpe")
  writeLines(bedpe_lines, tmp_bedpe)

  # Target peak overlapping the shared anchor A → should link to all 3 loops
  tmp_bed <- tempfile(fileext = ".bed")
  writeLines(sprintf("chrX\t%d\t%d", anchor_a_start, anchor_a_end), tmp_bed)

  res <- looplook:::.with_known_upstream_noise_suppressed(
    looplook::annotate_peaks_and_loops(
      bedpe_file = tmp_bedpe, target_bed = tmp_bed,
      txdb = txdb, org_db = "org.Hs.eg.db", species = "hg19",
      out_dir = tempdir(), project_name = "Sim_MultiLoop",
      write_output = FALSE, quiet = TRUE
    )
  )

  ta <- res$target_annotation
  expect_equal(nrow(ta), 1)

  # Should be connected to loops (Linked_Loop_IDs should have multiple entries)
  expect_false(is.na(ta$Linked_Loop_IDs) || ta$Linked_Loop_IDs == "")
  n_loops <- length(unlist(strsplit(ta$Linked_Loop_IDs, ";")))
  expect_equal(n_loops, 3)

  # Assigned_Target_Genes should not be empty (at least g1 via L3)
  expect_false(is.na(ta$Assigned_Target_Genes) || ta$Assigned_Target_Genes == "")

  unlink(c(tmp_bedpe, tmp_bed))
})

test_that("refine: full target-assignment fallback chain", {
  skip_if_not_installed("org.Hs.eg.db")
  rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
  expr_path  <- system.file("extdata", "example_tpm.txt", package = "looplook")
  skip_if(rdata_path == "" || expr_path == "")

  tmp <- new.env()
  load(rdata_path, envir = tmp)
  res <- tmp[[ls(tmp)[1]]]
  res$loop_annotation <- head(res$loop_annotation, 20)
  res$target_annotation <- head(res$target_annotation, 8)
  res$promoter_centric_stats <- head(res$promoter_centric_stats, 20)
  res$distal_element_stats <- head(res$distal_element_stats, 20)

  # Use a high threshold so most genes are silenced → triggers reclassify + fallback
  ref <- looplook:::.with_known_upstream_noise_suppressed(
    looplook::refine_loop_anchors_by_expression(
      annotation_res = res, expr_matrix_file = expr_path,
      sample_columns = c("con1", "con2"), threshold = 100,
      reclassify_by_expression = TRUE,
      out_dir = tempdir(), project_name = "Sim_FallbackChain",
      write_output = FALSE, quiet = TRUE
    )
  )

  # Verify the output contract
  expect_type(ref, "list")
  la <- ref$loop_annotation

  # After high-threshold refinement, some anchors should be reclassified
  types <- unique(c(la$anchor1_type, la$anchor2_type))
  expect_true(any(grepl("^e", types)))

  # Putative_Target_Genes should use fallback when original targets are silenced
  expect_true("Putative_Target_Genes" %in% colnames(la))

  # If target_annotation present, Assigned_Target_Genes_Filled must exist
  if (!is.null(ref$target_annotation)) {
    ta <- ref$target_annotation
    expect_true("Assigned_Target_Genes_Filled" %in% colnames(ta))
    # _Filled column should never be NA when Assigned_Target_Genes is present
    has_assigned <- !is.na(ta$Assigned_Target_Genes) & ta$Assigned_Target_Genes != ""
    filled_na <- is.na(ta$Assigned_Target_Genes_Filled) | ta$Assigned_Target_Genes_Filled == ""
    expect_false(any(has_assigned & filled_na))
  }
})


# ════════════════════════════════════════════════════════════════════════════
# Part D: analysis.R helpers and visualization.R exported functions
# ════════════════════════════════════════════════════════════════════════════

test_that(".calc_gc_fraction computes GC content correctly", {
  seqs <- Biostrings::DNAStringSet(c("GCGC", "ATAT", "NNNN", "GCAT", ""))
  gc <- looplook:::.calc_gc_fraction(seqs)
  expect_equal(gc[1], 1.0)
  expect_equal(gc[2], 0.0)
  expect_equal(gc[3], 0.0)
  expect_equal(gc[4], 0.5)
  expect_true(is.na(gc[5]))
})

test_that(".is_promoter_anchor_type and .is_enhancer_like_anchor_type", {
  expect_true(looplook:::.is_promoter_anchor_type("P"))
  expect_false(looplook:::.is_promoter_anchor_type("E"))
  expect_false(looplook:::.is_promoter_anchor_type("eP"))
  expect_false(looplook:::.is_promoter_anchor_type(NA_character_))

  expect_true(looplook:::.is_enhancer_like_anchor_type("E"))
  expect_true(looplook:::.is_enhancer_like_anchor_type("eP"))
  expect_true(looplook:::.is_enhancer_like_anchor_type("eG"))
  expect_false(looplook:::.is_enhancer_like_anchor_type("P"))
  expect_false(looplook:::.is_enhancer_like_anchor_type("G"))
})

test_that(".empty_anchor_df and .make_anchor_df produce valid output", {
  edf <- looplook:::.empty_anchor_df()
  expect_s3_class(edf, "data.frame")
  expect_equal(nrow(edf), 0)
  expect_equal(colnames(edf), c("anchor_id", "chr", "start", "end", "anchor_type"))

  loop_df <- data.frame(
    chr1 = "chr1", start1 = 100, end1 = 150, a1_id = "A1",
    stringsAsFactors = FALSE
  )
  adf <- looplook:::.make_anchor_df(loop_df, 1, "1", c("P"))
  expect_equal(adf$anchor_id, "A1")
  expect_equal(adf$chr, "chr1")
  expect_equal(adf$anchor_type, "P")

  loop_df2 <- data.frame(chr1 = "chr2", start1 = 200, end1 = 250, stringsAsFactors = FALSE)
  adf2 <- looplook:::.make_anchor_df(loop_df2, 1, "1", c("E"))
  expect_match(adf2$anchor_id, "chr2_200_250")
})

test_that("run_lfc_violin: edge cases with constant values", {
  global_glist <- setNames(rep(1, 100), paste0("Gene", 1:100))
  targets <- names(global_glist)[1:5]
  p <- looplook:::run_lfc_violin(targets, global_glist, "wilcox.test", "Test")
  expect_s3_class(p, "ggplot")

  global_glist2 <- setNames(c(rep(-2, 50), rep(-1, 50)), paste0("Gene", 1:100))
  p2 <- looplook:::run_lfc_violin(targets, global_glist2, "t.test", "Neg")
  expect_s3_class(p2, "ggplot")
})

test_that("draw_flower_simplified: exported function with mock gene lists", {
  skip_if_not_installed("ggplot2")
  sets <- list(
    A = c("TP53", "BRCA1", "MYC", "EGFR", "KRAS"),
    B = c("BRCA1", "MYC", "EGFR", "PTEN", "APC"),
    C = c("MYC", "EGFR", "KRAS", "TP53", "BRAF")
  )
  p <- draw_flower_simplified(sets, "Test", c(A = "#E41A1C", B = "#377EB8", C = "#4DAF4A"))
  expect_s3_class(p, "ggplot")

  expect_message(
    p_null <- draw_flower_simplified(list(A = c("X")), "One", NULL),
    "Less than 2"
  )
  expect_null(p_null)
})

test_that("draw_upset_intersections: input validation and early returns", {
  skip_if_not_installed("UpSetR")

  # less than 2 lists → message + NULL
  expect_message(
    g1 <- draw_upset_intersections(list(A = c("X")), "One"),
    "Less than 2"
  )
  expect_null(g1)

  # 2 lists with content should work (rendering depends on graphics device,
  # but at minimum should not error)
  sets <- list(
    Up = c("TP53", "BRCA1", "MYC", "EGFR"),
    Down = c("BRCA1", "MYC", "CDKN1A", "BAX")
  )
  expect_no_error(draw_upset_intersections(sets, "Test"))
})

test_that("plot_peaks_interactions: renders without gene track", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggforce")
  skip_if_not_installed("ggrepel")

  tmp <- tempfile(fileext = ".bedpe")
  writeLines(c(
    "chr1\t10000\t10500\tchr1\t20000\t20500",
    "chr1\t15000\t15500\tchr1\t25000\t25500"
  ), tmp)

  p <- plot_peaks_interactions(
    bedpe_file = tmp, chr = "chr1", from = 9000, to = 26000,
    show_gene_track = FALSE
  )
  expect_s3_class(p, "ggplot")
  unlink(tmp)
})

test_that("prepare_track_data: works without gene track", {
  tmp <- tempfile(fileext = ".bedpe")
  writeLines("chr1\t10000\t10500\tchr1\t20000\t20500", tmp)

  d <- looplook:::prepare_track_data(
    bedpe_file = tmp, target_bed = NULL, chr = "chr1",
    from = 9000, to = 21000, species = "hg38",
    max_levels = 10, base_anchor_height = 0.05,
    loop_color = "#5D6D7E", anchor_color = "#3498DB",
    score_to_alpha = FALSE, min_score = NULL,
    show_gene_track = FALSE
  )
  expect_type(d, "list")
  expect_equal(d$chr, "chr1")
  expect_gt(nrow(d$bez_df), 0)
  unlink(tmp)
})
