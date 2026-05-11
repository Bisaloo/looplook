# tests/testthat/test-helpers.R — unit tests for utility helper functions

test_that("clean_gene_names handles normal, edge, and NA cases", {
  expect_equal(looplook:::clean_gene_names(c("A", "B", "A", "C")), c("A", "B", "C"))
  expect_equal(looplook:::clean_gene_names(c("A", "", "B", NA, " ")), c("A", "B"))
  expect_equal(looplook:::clean_gene_names(character(0)), character(0))
  expect_equal(looplook:::clean_gene_names(NULL), character(0))

  # with split
  expect_equal(looplook:::clean_gene_names("TP53;BRCA1;TP53", ";"), c("TP53", "BRCA1"))
  expect_equal(looplook:::clean_gene_names(c("A;B", "B;C", "D")), c("A;B", "B;C", "D"))
  expect_equal(looplook:::clean_gene_names(c("", NA, ";"), "[;,]"), character(0))
})

test_that("extract_genes collapses delimited gene strings correctly", {
  expect_equal(looplook:::extract_genes("TP53;BRCA1;TP53"), "TP53;BRCA1")
  expect_equal(looplook:::extract_genes(c("TP53", "BRCA1")), "TP53;BRCA1")
  expect_equal(looplook:::extract_genes(NA_character_), NA_character_)
  expect_equal(looplook:::extract_genes(""), NA_character_)
  expect_equal(looplook:::extract_genes(";"), NA_character_)
  # single gene
  expect_equal(looplook:::extract_genes("MYC"), "MYC")
})

test_that("load_expression_matrix loads and averages correctly", {
  tmp <- tempfile(fileext = ".tsv")
  write.table(data.frame(
    GeneID = c("TP53", "BRCA1", "MYC"),
    con1 = c(10, 20, 30),
    con2 = c(15, 25, 35),
    trt1 = c(100, 200, 300)
  ), tmp, sep = "\t", row.names = FALSE, quote = FALSE)

  vals <- looplook:::load_expression_matrix(tmp, c("con1", "con2"))
  expect_named(vals, c("TP53", "BRCA1", "MYC"))
  expect_equal(vals[["TP53"]], 12.5)

  # error cases
  expect_error(
    looplook:::load_expression_matrix("nonexistent.tsv", "con1"),
    "not found"
  )
  expect_error(
    looplook:::load_expression_matrix(tmp, "badcol"),
    "Requested sample columns not found"
  )
  expect_error(
    looplook:::load_expression_matrix(tmp, c("con1", "con1")),
    "`sample_columns` contains duplicates"
  )

  unlink(tmp)
})

test_that("load_expression_matrix errors on duplicated expression column names", {
  tmp <- tempfile(fileext = ".tsv")
  writeLines(
    c(
      "Gene\tcon1\tcon1\ttrt1",
      "TP53\t10\t15\t100",
      "BRCA1\t20\t25\t200"
    ),
    tmp
  )

  expect_error(
    looplook:::load_expression_matrix(tmp, c("con1", "trt1")),
    "duplicated sample column names"
  )

  unlink(tmp)
})

test_that("get_colors returns expected vector lengths and fallbacks", {
  expect_no_warning(expect_length(looplook:::get_colors(0, "Set2"), 0))
  expect_no_warning(expect_length(looplook:::get_colors(1, "Set2"), 1))
  expect_no_warning(expect_length(looplook:::get_colors(2, "Set2"), 2))
  expect_no_warning(expect_length(looplook:::get_colors(5, "Set2"), 5))
  expect_no_warning(expect_length(looplook:::get_colors(10, c("#E41A1C", "#377EB8")), 10))
  # auto-generate when NULL/empty
  cols <- looplook:::get_colors(6, NULL)
  expect_length(cols, 6)
})

test_that("read_robust_general handles empty and short files", {
  tmp <- tempfile(fileext = ".tsv")
  writeLines("A\tB\tC\n1\t2\t3\n4\t5\t6", tmp)
  res <- looplook:::read_robust_general(tmp, header = TRUE, min_cols = 2)
  expect_equal(nrow(res), 2)
  expect_equal(ncol(res), 3)

  expect_error(
    looplook:::read_robust_general("", desc = "Test"),
    "path is empty"
  )
  unlink(tmp)
})

test_that("resolve_gene_conflicts returns input unchanged when empty", {
  df <- data.frame(
    chr = character(0), start = integer(0), end = integer(0),
    stringsAsFactors = FALSE
  )
  expect_equal(nrow(looplook:::resolve_gene_conflicts(
    df, NULL, NULL, c(-2000, 2000), NULL
  )), 0)
})

test_that("clean_anchor filters and reclassifies correctly", {
  res <- looplook:::clean_anchor("TP53;BRCA1", "P", c("TP53"), TRUE)
  expect_equal(res$type, "P")
  expect_equal(res$gene, "TP53")

  res2 <- looplook:::clean_anchor("GENE_X", "P", c("TP53"), TRUE)
  expect_equal(res2$type, "eP")
  expect_true(is.na(res2$gene))

  res3 <- looplook:::clean_anchor("GENE_X", "P", c("TP53"), FALSE)
  expect_equal(res3$type, "P")
  expect_true(is.na(res3$gene))

  res4 <- looplook:::clean_anchor("", "E", c("TP53"), FALSE)
  expect_equal(res4$type, "E")
  expect_true(is.na(res4$gene))
})

test_that(".map_txdb_gene_ids infers OrgDb keytypes for TxDb-like gene IDs", {
  skip_if_not_installed("org.Hs.eg.db")
  org_db <- org.Hs.eg.db::org.Hs.eg.db

  entrez_ids <- head(AnnotationDbi::keys(org_db, keytype = "ENTREZID"), 5)
  ensembl_ids <- head(AnnotationDbi::keys(org_db, keytype = "ENSEMBL"), 5)

  entrez_map <- looplook:::.map_txdb_gene_ids(
    gene_ids = entrez_ids,
    org_db = org_db,
    columns = "SYMBOL",
    warn = FALSE
  )
  ensembl_map <- looplook:::.map_txdb_gene_ids(
    gene_ids = ensembl_ids,
    org_db = org_db,
    columns = "SYMBOL",
    warn = FALSE
  )

  expect_identical(attr(entrez_map, "keytype"), "ENTREZID")
  expect_identical(attr(ensembl_map, "keytype"), "ENSEMBL")
  expect_true(any(!is.na(entrez_map$SYMBOL)))
  expect_true(any(!is.na(ensembl_map$SYMBOL)))
})

test_that(".map_txdb_gene_ids safely falls back when no OrgDb keytype matches", {
  skip_if_not_installed("org.Hs.eg.db")

  map <- looplook:::.map_txdb_gene_ids(
    gene_ids = c("not_a_real_gene_1", "still_not_real_2"),
    org_db = org.Hs.eg.db::org.Hs.eg.db,
    columns = "SYMBOL",
    warn = FALSE
  )

  expect_true(is.na(attr(map, "keytype")))
  expect_identical(attr(map, "hit_rate"), 0)
  expect_true(all(is.na(map$SYMBOL)))
})
