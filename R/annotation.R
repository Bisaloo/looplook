#' Internal: Collapse IDs
#' @keywords internal
#' @noRd
.annotation_extract_ids <- function(id_vec) {
  paste(unique(na.omit(as.character(id_vec))), collapse = ";")
}

#' Internal: Resolve AnnotationDb Package Name
#' @keywords internal
#' @noRd
.pkg_from_annotation_db <- function(x) {
  db_path <- tryCatch(AnnotationDbi::dbfile(x), error = function(e) NULL)
  if (!is.null(db_path) && length(db_path) == 1L && nzchar(db_path)) {
    return(basename(dirname(dirname(db_path))))
  }
  pkg_attr <- attr(x, "package")
  if (!is.null(pkg_attr) && length(pkg_attr) == 1L && nzchar(pkg_attr) &&
    pkg_attr != "AnnotationDbi") {
    return(pkg_attr)
  }
  NULL
}

#' Internal: Resolve Annotation Resource
#' @keywords internal
#' @noRd
.resolve_annotation_resource <- function(arg, type, desc, species) {
  if (any(inherits(arg, c("TxDb", "OrgDb", "AnnotationDb")))) {
    return(list(obj = arg, pkg = .pkg_from_annotation_db(arg)))
  }
  if (is.character(arg) && nzchar(arg)) {
    if (!requireNamespace(arg, quietly = TRUE)) stop(desc, " '", arg, "' not installed")
    return(list(obj = utils::getFromNamespace(arg, arg), pkg = arg))
  }
  pkg <- if (type == "txdb") species_txdb_pkg(species) else species_orgdb_pkg(species)
  if (!requireNamespace(pkg, quietly = TRUE)) stop(desc, " '", pkg, "' not installed")
  list(obj = utils::getFromNamespace(pkg, pkg), pkg = pkg)
}

#' Internal: Convert ChIPseeker Annotation to Anchor Class
#' @keywords internal
#' @noRd
.annotation_feature_class <- function(anno_str) {
  if (length(anno_str) == 0 || is.na(anno_str)) {
    return("Unknown")
  }
  anno_str <- tolower(anno_str)

  if (grepl("promoter", anno_str)) {
    return("P")
  }

  if (grepl("intergenic|downstream", anno_str)) {
    return("E")
  }
  if (grepl("exon|intron|utr", anno_str)) {
    return("G")
  }

  "E"
}

#' Internal: Extract Loop Locus Genes
#' @keywords internal
#' @noRd
.loop_locus_genes <- function(t1, t2, s1, s2) {
  genes <- c()
  if (!is.na(t1) && t1 %in% c("P", "G")) genes <- c(genes, s1)
  if (!is.na(t2) && t2 %in% c("P", "G")) genes <- c(genes, s2)
  paste(unique(na.omit(genes)), collapse = ";")
}

#' Internal: Build Loop Type Code
#' @keywords internal
#' @noRd
.loop_type_code <- function(t1, t2) {
  if (length(t1) == 0 || length(t2) == 0 || is.na(t1) || is.na(t2)) {
    return("Unknown")
  }
  paste(sort(c(t1, t2)), collapse = "-")
}

#' Internal: Convert Anchor IDs to Genes
#' @keywords internal
#' @noRd
.ids_to_genes_simple <- function(ids, lookup) {
  valid <- intersect(ids, names(lookup))
  if (length(valid) == 0) {
    return(NA_character_)
  }
  genes <- lookup[valid]
  genes <- genes[!is.na(genes)]
  if (length(genes) == 0) {
    return(NA_character_)
  }
  paste(sort(unique(genes)), collapse = ";")
}

#' Internal: Convert Anchor IDs to Promoter-Priority Genes
#' @keywords internal
#' @noRd
.ids_to_genes_priority <- function(ids, lookup_sym, lookup_typ) {
  valid <- intersect(ids, names(lookup_sym))
  if (length(valid) == 0) {
    return(NA_character_)
  }
  promoter_ids <- valid[!is.na(lookup_typ[valid]) & lookup_typ[valid] == "P"]
  use_ids <- if (length(promoter_ids) > 0) promoter_ids else valid
  genes_present <- lookup_sym[use_ids]
  genes_present <- genes_present[!is.na(genes_present)]
  if (length(genes_present) == 0) {
    return(NA_character_)
  }
  paste(sort(unique(genes_present)), collapse = ";")
}

#' Internal: Fill Empty Target Gene Assignments
#' @keywords internal
#' @noRd
.fill_target_gene_fallback <- function(target, fallback) {
  dplyr::case_when(
    !is.na(target) & target != "" ~ target,
    !is.na(fallback) & fallback != "" ~ fallback,
    TRUE ~ NA_character_
  )
}

#' Internal: Integrate Optional Target BED Annotation
#' @keywords internal
#' @noRd
.annotate_target_bed <- function(
  target_bed, txdb_obj, org_db_pkg, tss_region, gene_expr_map, min_expr,
  conflict_strategy,
  gr_anchors, anchor_topo_map, loop_annotation_final, map_info, ego_list_target,
  log_message
) {
  bed_target <- read_robust_general(target_bed, min_cols = 3, desc = "Target BED")
  colnames(bed_target)[c(1, 2, 3)] <- c("chr", "start", "end")
  if (nrow(bed_target) == 0) {
    warning("Target BED contains no features; skipping target annotation.")
    return(list(
      bed_info = NULL,
      target_connected_loops = NULL,
      target_gene_links = NULL
    ))
  }

  bed_target$start <- bed_target$start + 1 # BED is 0-based; GRanges is 1-based
  gr_bed <- .with_known_upstream_noise_suppressed({
    gr_bed <- GenomicRanges::makeGRangesFromDataFrame(bed_target)
    gr_bed$input_id <- paste0("Peak_", seq_len(nrow(bed_target)))
    names(gr_bed) <- gr_bed$input_id
    gr_bed
  })
  bed_annot <- .with_known_upstream_noise_suppressed(
    ChIPseeker::annotatePeak(gr_bed, TxDb = txdb_obj, tssRegion = tss_region, annoDb = org_db_pkg, verbose = FALSE)
  )
  bed_info <- format_annotation_columns(as.data.frame(bed_annot))
  if ("GENENAME" %in% colnames(bed_info)) bed_info <- bed_info %>% dplyr::rename(Gene_description = GENENAME)
  log_message("    Refining Target annotation...")
  bed_info <- resolve_gene_conflicts(bed_info, txdb_obj, org_db_pkg, tss_region, gene_expr_map, min_expr = min_expr, conflict_strategy = conflict_strategy)
  gr_bed <- .harmonize_seqlevels(gr_bed, gr_anchors, "target BED")
  hits <- GenomicRanges::findOverlaps(gr_bed, gr_anchors)

  target_connected_loops <- NULL
  target_gene_links <- NULL
  if (length(hits) > 0) {
    target_connected_loops <- loop_annotation_final %>%
      dplyr::filter(cluster_id %in% unique(gr_anchors$cluster_id[S4Vectors::subjectHits(hits)]))
    hit_df <- data.frame(qid = S4Vectors::queryHits(hits), sid = S4Vectors::subjectHits(hits))
    hit_df$anchor_id <- gr_anchors$anchor_id[hit_df$sid]
    hit_df <- hit_df %>% dplyr::left_join(anchor_topo_map, by = "anchor_id")
    anchor_loop_agg <- dplyr::bind_rows(
      loop_annotation_final %>% dplyr::select(anchor_id = a1_id, loop_ID),
      loop_annotation_final %>% dplyr::select(anchor_id = a2_id, loop_ID)
    ) %>%
      dplyr::distinct() %>%
      dplyr::group_by(anchor_id) %>%
      dplyr::summarise(linked_loops = .annotation_extract_ids(loop_ID), .groups = "drop")
    hit_df <- hit_df %>% dplyr::left_join(anchor_loop_agg, by = "anchor_id")
    summary_df <- hit_df %>%
      dplyr::group_by(qid) %>%
      dplyr::summarise(
        All_Loop_Connected_Genes = extract_genes(tgt_genes_pg),
        Regulated_promoter_genes = extract_genes(tgt_genes_p),
        Assigned_Target_Genes = extract_genes(tgt_genes_prio),
        Linked_Loop_IDs = .annotation_extract_ids(linked_loops),
        .groups = "drop"
      ) %>%
      dplyr::mutate(join_id = paste0("Peak_", qid))
    bed_info <- dplyr::left_join(bed_info, summary_df, by = c("input_id" = "join_id")) %>%
      dplyr::select(-any_of(c("join_id", "qid")))
    target_gene_links <- .build_target_gene_links(
      hit_df = hit_df,
      bed_info = bed_info,
      loop_annotation_final = loop_annotation_final,
      map_info = map_info,
      ego_list_target = ego_list_target
    )
  } else {
    bed_info$All_Loop_Connected_Genes <- NA
    bed_info$Regulated_promoter_genes <- NA
    bed_info$Assigned_Target_Genes <- NA
    bed_info$Linked_Loop_IDs <- NA
  }

  if (is.null(target_gene_links)) {
    target_gene_links <- .build_target_gene_links(
      hit_df = data.frame(qid = integer(0), anchor_id = character(0)),
      bed_info = bed_info,
      loop_annotation_final = loop_annotation_final,
      map_info = map_info,
      ego_list_target = ego_list_target
    )
  }

  evidence_df <- .summarise_regulated_promoter_evidence(target_gene_links)
  bed_info <- dplyr::left_join(bed_info, evidence_df, by = "input_id")
  bed_info$Regulated_promoter_Evidence <- ifelse(
    is.na(bed_info$Regulated_promoter_Evidence) |
      bed_info$Regulated_promoter_Evidence == "",
    "none",
    bed_info$Regulated_promoter_Evidence
  )

  fallback_col <- .target_linear_gene_column(bed_info)
  fallback_vec <- if (!is.null(fallback_col)) bed_info[[fallback_col]] else rep(NA_character_, nrow(bed_info))
  bed_info <- bed_info %>% dplyr::mutate(
    All_Loop_Connected_Genes_Filled = .fill_target_gene_fallback(All_Loop_Connected_Genes, fallback_vec),
    Regulated_promoter_genes_Filled = .fill_target_gene_fallback(Regulated_promoter_genes, fallback_vec),
    Assigned_Target_Genes_Filled = .fill_target_gene_fallback(Assigned_Target_Genes, fallback_vec),
    Regulated_promoter_Fallback_Evidence = dplyr::case_when(
      !is.na(Regulated_promoter_genes) & Regulated_promoter_genes != "" ~ "none",
      !is.na(Regulated_promoter_genes_Filled) & Regulated_promoter_genes_Filled != "" ~
        .fallback_evidence_from_annotation(annotation),
      TRUE ~ "none"
    )
  )
  target_gene_links <- .mark_target_gene_link_membership(target_gene_links, bed_info)

  if ("Linked_Loop_IDs" %in% colnames(bed_info)) {
    target_col <- if ("Gene_description" %in% colnames(bed_info)) "Gene_description" else "SYMBOL"
    if (target_col %in% colnames(bed_info)) {
      bed_info <- bed_info %>% dplyr::relocate(Linked_Loop_IDs, .after = dplyr::all_of(target_col))
    }
  }

  list(
    bed_info = bed_info,
    target_connected_loops = target_connected_loops,
    target_gene_links = target_gene_links
  )
}

#' @title Spatial mapping of 1D genomic features to 3D chromatin interaction targets
#'
#' @description
#' A dual-purpose analytical framework designed to integrate 1D genomic features with 3D chromatin architecture.
#' \enumerate{
#'   \item \strong{Loop Annotation:} Classifies 3D spatial interactions (e.g., Enhancer-Promoter, Promoter-Promoter) using a defined structural hierarchy.
#'   \item \strong{Feature-to-Target Mapping:} Links 1D genomic features (e.g., GWAS risk SNPs, ATAC-seq peaks, ChIP-seq binding sites) to putative target genes via 3D chromatin contacts, providing a spatial complement to linear proximity-based assignments.
#' }
#'
#' @details
#' \strong{Mapping Strategy and Fallback Mechanism}
#' The method prioritizes physical 3D chromatin contacts while keeping strict
#' and fallback semantics separate. \code{Regulated_promoter_genes} reports
#' promoter genes supported by loop-anchor context, \code{Assigned_Target_Genes}
#' preserves the historical promoter-first 3D assignment, and \code{*_Filled}
#' columns add a linear nearest-gene fallback only when strict 3D assignments are
#' empty. \code{Regulated_promoter_Evidence},
#' \code{Regulated_promoter_Fallback_Evidence}, and \code{target_gene_links}
#' provide row-level provenance for these decisions.
#'
#' \strong{Hierarchical Conflict Resolution}
#' To address complex loci where a single anchor overlaps multiple promoters (e.g., dense gene clusters or bidirectional promoters), the function executes a 3-step resolution:
#' \enumerate{
#'   \item \emph{Biotype Prioritization:} Selects the highest-priority candidates by functional class: \code{Protein Coding > small-ncRNA (miRNA, snoRNA, snRNA, rRNA, scaRNA) > Antisense > lncRNA/ncRNA > Pseudogene}.
#'   \item \emph{Expression Filter:} Within the selected biotype tier, excludes transcriptionally silent genes using a user-provided expression matrix. If no gene in the tier is expressed, all candidates in that tier are retained.
#'   \item \emph{Expression Tiebreaker:} Among remaining candidates of equal biotype priority, retains all genes whose expression is within one order of magnitude of the highest-expressing candidate (i.e., expression >= 10\% of the group maximum). This co-dominant rule preserves functionally redundant candidates such as bidirectional promoter pairs.
#' }
#'
#' \strong{Network Topology Analysis}
#' \itemize{
#'   \item \strong{Ego-Network Expansion (\code{neighbor_hop}):} Implements k-hop neighborhood expansion via \code{igraph::ego()}. A value of \code{0} restricts to direct contacts, while \code{1} includes secondary contacts to capture broader regulatory cliques.
#'   \item \strong{Hub Detection:} Utilizes a node-degree quantile threshold (\code{hub_percentile}) to identify highly connected regulatory elements.
#' }
#'
#' @param bedpe_file Character. Path to a BEDPE file (at least 6 columns: chr1, start1, end1, chr2, start2, end2).
#' @param target_bed Optional path to a BED file of genomic features (e.g., ChIP-seq peaks, GWAS SNPs). When provided, these 1D regions are mapped to 3D target genes. Default: \code{NULL}.
#' @param txdb A \code{\link[GenomicFeatures]{TxDb}} object, a package name string, or \code{NULL} to auto-resolve from \code{species}. Default: \code{NULL}.
#' @param org_db An \code{OrgDb} object, a package name string, or \code{NULL} to auto-resolve from \code{species}. Default: \code{NULL}.
#' @param species Character. Genome assembly used when \code{txdb} and \code{org_db} are \code{NULL}. One of \code{"hg38"}, \code{"hg19"}, \code{"mm10"}, \code{"mm9"}. Default: \code{"hg38"}.
#' @param tss_region Numeric vector of length 2. Promoter window around the TSS in bp. Default: \code{c(-2000, 2000)}.
#' @param out_dir Character. Output directory for the Excel results file. Default: \code{"./results"}.
#' @param expr_matrix_file Optional path to a normalised expression matrix (TPM/FPKM, genes × samples). Enables expression-aware conflict resolution. Default: \code{NULL}.
#' @param sample_columns Character vector or integer indices. Columns in \code{expr_matrix_file} to average for baseline expression. Default: \code{NULL}.
#' @param min_expr Numeric. Minimum expression value for a gene to be considered
#'   active during anchor-level conflict resolution. Used only when
#'   \code{expr_matrix_file} is provided. Default: \code{0} (any detectable
#'   expression qualifies). Increase to \code{1} or higher to require stronger
#'   evidence. See \code{\link{refine_loop_anchors_by_expression}} for the
#'   separate \code{threshold} parameter that controls promoter reclassification.
#' @param conflict_strategy Character. Conflict resolution order for
#'   overlapping gene assignments. \code{"biotype_first"} (default): select
#'   the best biotype tier first, then apply expression filtering within that
#'   tier. \code{"expression_first"}: apply expression filtering across all
#'   biotypes first, then pick the best biotype among expressed candidates.
#' @param project_name Character. Prefix for output files and plot titles. Default: \code{"HiChIP"}.
#' @param color_palette Character. RColorBrewer palette name. Default: \code{"Set2"}.
#' @param karyo_bin_size Integer. Bin width in bp for karyotype heatmaps. Default: \code{1e5}.
#' @param neighbor_hop Integer. k-hop ego-network expansion order via \code{igraph::ego()}. \code{0} restricts to direct contacts. Default: \code{0}.
#' @param hub_percentile Numeric (0–1). Node-degree quantile for hub detection. Default: \code{0.95}.
#' @param write_output Logical. If \code{TRUE} (default), write the Excel workbook to \code{out_dir}. If \code{FALSE}, return results without creating directories or files.
#' @param quiet Logical. If \code{TRUE}, suppress progress messages while preserving warnings. Default: \code{FALSE}.
#'
#' @return An invisible named list:
#' \itemize{
#'   \item \code{target_annotation} — Target features (peaks) with gene assignments.
#'     Key columns include:
#'     \itemize{
#'       \item \code{All_Loop_Connected_Genes}: All genes from loop-connected anchors (P/G types).
#'       \item \code{Regulated_promoter_genes}: Promoter genes supported by loop-anchor context.
#'       \item \code{Assigned_Target_Genes}: Promoter-first 3D assignment (prioritises P > G > E).
#'       \item \code{*_Filled} variants: Linear nearest-gene fallback when strict 3D assignments are empty.
#'       \item \code{Regulated_promoter_Evidence}: Provenance of \code{Regulated_promoter_genes}
#'         (e.g., \code{local_promoter_overlap}, \code{direct_opposite_promoter}).
#'         \strong{Read with} \code{Regulated_promoter_genes}; do not cross-reference
#'         with \code{Assigned_Target_Genes} or other columns.
#'       \item \code{Regulated_promoter_Fallback_Evidence}: Provenance of
#'         \code{Regulated_promoter_genes_Filled}.
#'         \strong{Read with} \code{Regulated_promoter_genes_Filled}; indicates
#'         which \code{*_Filled} column supplied the fallback gene.
#'     }
#'   \item \code{target_gene_links} — Long-format peak-gene provenance table.
#'     Each row records one peak-gene linkage with full provenance.
#'     \strong{Read} \code{evidence}, \code{anchor_role}, and \code{gene_role}
#'     \strong{together as a group} — they jointly describe how each gene was
#'     assigned to each peak; do not interpret any one column in isolation.
#'     \itemize{
#'       \item \code{input_id}, \code{loop_ID}, \code{anchor_id}: Identifiers.
#'       \item \code{gene}: Linked gene symbol.
#'       \item \code{gene_role}: \code{"promoter"}, \code{"gene_body"}, or \code{"linear_annotation"}.
#'       \item \code{source}: \code{"loop_anchor"} (3D-derived) or \code{"linear_annotation"} (nearest gene).
#'       \item \code{evidence}: Provenance label —
#'         \code{"local_promoter_overlap"} (peak overlaps anchor promoter),
#'         \code{"direct_opposite_promoter"} (opposite anchor is promoter),
#'         \code{"gene_body_context"} (gene body linkage),
#'         \code{"expanded_promoter_loop"} (via ego-network expansion),
#'         \code{"linear_annotation"} (direct nearest gene),
#'         or \code{"linear_fallback"} (filled when 3D assignment was empty).
#'       \item \code{anchor_role}: \code{"local_anchor"}, \code{"opposite_anchor"},
#'         \code{"expanded_anchor"}, or \code{"linear_annotation"}.
#'       \item \code{used_as_fallback}: Logical. \code{TRUE} when this link was added
#'         via the \code{*_Filled} linear nearest-gene fallback mechanism.
#'       \item \code{in_regulated_promoter} through \code{in_assigned_target_filled}:
#'         Logical membership flags indicating which target annotation column(s)
#'         this gene appears in.
#'     }
#'   \item \code{loop_annotation} — Annotated 3D interactome with \code{Putative_Target_Genes}.
#'   \item \code{anchor_loci_annotation} — Non-redundant anchor-locus genomic classifications after within-cluster interval reduction.
#'   \item \code{anchor_annotation} — Backward-compatible alias of \code{anchor_loci_annotation}.
#'   \item \code{promoter_centric_stats} — Gene-level connectivity statistics.
#'   \item \code{distal_element_stats} — Distal-element connectivity statistics.
#'   \item \code{plots} — Named list of ggplot objects (donut, karyotype, rose, flower).
#'   \item \code{plot_list} — Backward-compatible alias of \code{plots}.
#' }
#' If \code{write_output = TRUE}, also writes a multi-sheet Excel workbook to \code{out_dir}.
#'
#' @export
#'
#' @examples
#' # Minimal runnable example for package checks
#' if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
#'   txdb_example <- AnnotationDbi::loadDb(
#'     system.file("extdata", "hg19_knownGene_sample.sqlite", package = "GenomicFeatures")
#'   )
#'   bedpe_path <- tempfile(fileext = ".bedpe")
#'   writeLines(
#'     "chr6\t10412000\t10412600\tchr6\t10415000\t10415600",
#'     bedpe_path
#'   )
#'
#'   res <- annotate_peaks_and_loops(
#'     bedpe_file = bedpe_path,
#'     txdb = txdb_example,
#'     org_db = "org.Hs.eg.db",
#'     species = "hg19",
#'     out_dir = tempdir(),
#'     project_name = "Quick_Example",
#'     write_output = FALSE,
#'     quiet = TRUE
#'   )
#'   head(res$loop_annotation, 1)
#' }
annotate_peaks_and_loops <- function(
  bedpe_file,
  target_bed = NULL,
  txdb = NULL,
  org_db = NULL,
  species = "hg38",
  tss_region = c(-2000, 2000),
  out_dir = "./results",
  expr_matrix_file = NULL,
  sample_columns = NULL,
  project_name = "HiChIP",
  color_palette = "Set2",
  karyo_bin_size = 1e5,
  neighbor_hop = 0,
  hub_percentile = 0.95,
  min_expr = 0,
  conflict_strategy = c("biotype_first", "expression_first"),
  write_output = TRUE,
  quiet = FALSE
) {
  species <- match.arg(species, c("hg38", "hg19", "mm10", "mm9"))
  stopifnot(length(tss_region) == 2L)
  conflict_strategy <- match.arg(conflict_strategy)
  log_message <- function(...) {
    if (!quiet) message(...)
  }

  if (write_output && !dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  tx_res <- .resolve_annotation_resource(txdb, "txdb", "TxDb", species)
  org_res <- .resolve_annotation_resource(org_db, "orgdb", "OrgDb", species)
  txdb_obj <- tx_res$obj
  org_db_pkg <- org_res$pkg
  if (is.null(org_db_pkg) || !nzchar(org_db_pkg)) {
    stop(
      "Unable to resolve OrgDb package name from supplied object. ",
      "Pass a package name string such as 'org.Hs.eg.db', or an installed OrgDb package object."
    )
  }

  gene_expr_map <- NULL
  if (!is.null(expr_matrix_file) && !is.null(sample_columns)) {
    log_message("Step 0: Loading expression data...")
    gene_expr_map <- load_expression_matrix(expr_matrix_file, sample_columns)
    log_message("    >>> Expression loaded for ", length(gene_expr_map), " genes.")
  } else if (!is.null(expr_matrix_file) && is.null(sample_columns)) {
    warning("expr_matrix_file provided but sample_columns is NULL; expression data will not be used for conflict resolution.", call. = FALSE)
  }

  log_message("Step 1: Reading BEDPE file...")
  loops <- read_robust_general(bedpe_file, min_cols = 6, desc = "BEDPE")
  colnames(loops)[seq_len(6)] <- c("chr1", "start1", "end1", "chr2", "start2", "end2")

  loops$start1 <- loops$start1 + 1
  loops$start2 <- loops$start2 + 1
  anchors <- dplyr::bind_rows(loops %>% dplyr::select(chr = chr1, start = start1, end = end1), loops %>% dplyr::select(chr = chr2, start = start2, end = end2)) %>%
    dplyr::distinct() %>%
    dplyr::mutate(anchor_id = paste0("A", seq_len(dplyr::n())))
  loops <- loops %>%
    dplyr::left_join(anchors %>% dplyr::select(chr, start, end, a1_id = anchor_id), by = c("chr1" = "chr", "start1" = "start", "end1" = "end")) %>%
    dplyr::left_join(anchors %>% dplyr::select(chr, start, end, a2_id = anchor_id), by = c("chr2" = "chr", "start2" = "start", "end2" = "end"))

  log_message("Step 2: Clustering loops...")
  valid_loops <- loops %>% dplyr::filter(!is.na(a1_id) & !is.na(a2_id))
  g <- igraph::graph_from_data_frame(valid_loops[, c("a1_id", "a2_id")], directed = FALSE)
  comp <- igraph::components(g)
  anchors$cluster_id <- NA
  comm <- intersect(anchors$anchor_id, names(comp$membership))
  if (length(comm) > 0) anchors$cluster_id[match(comm, anchors$anchor_id)] <- as.character(comp$membership[comm])
  anchors <- anchors %>% dplyr::filter(!is.na(cluster_id))
  loops <- loops %>% dplyr::left_join(anchors %>% dplyr::select(anchor_id, cluster_id), by = c("a1_id" = "anchor_id"))
  gr_anchors <- .with_known_upstream_noise_suppressed({
    gr_anchors <- GenomicRanges::makeGRangesFromDataFrame(anchors, keep.extra.columns = TRUE)
    gr_anchors$anchor_id <- anchors$anchor_id
    gr_anchors
  })
  cluster_regions <- .with_known_upstream_noise_suppressed({
    gr_list <- GenomicRanges::GRangesList(split(gr_anchors, gr_anchors$cluster_id))
    cluster_regions <- unlist(GenomicRanges::reduce(gr_list))
    cluster_regions$cluster_id <- names(cluster_regions)
    names(cluster_regions) <- paste0("peak_", seq_along(cluster_regions))
    cluster_regions
  })

  log_message("Step 3: Biological Classification & Topology...")
  if (length(gr_anchors) == 0) {
    warning("No valid loop anchors found; returning empty annotation.")
    return(list(
      anchor_loci_annotation = data.frame(), anchor_annotation = data.frame(), loop_annotation = data.frame(),
      promoter_centric_stats = data.frame(), distal_element_stats = data.frame(),
      target_annotation = NULL, target_gene_links = NULL, plots = list(), plot_list = list()
    ))
  }
  anchor_anno <- .with_known_upstream_noise_suppressed(
    ChIPseeker::annotatePeak(gr_anchors, TxDb = txdb_obj, tssRegion = tss_region, annoDb = org_db_pkg, verbose = FALSE)
  )
  anchor_anno_df <- format_annotation_columns(as.data.frame(anchor_anno))
  anchor_anno_df <- resolve_gene_conflicts(anchor_anno_df, txdb_obj, org_db_pkg, tss_region, gene_expr_map, min_expr = min_expr, conflict_strategy = conflict_strategy)
  anchor_anno_df$type_code <- vapply(anchor_anno_df$annotation, .annotation_feature_class, FUN.VALUE = character(1))
  map_info <- anchor_anno_df %>% dplyr::select(anchor_id, type_code, SYMBOL)
  loops_annotated <- loops %>%
    dplyr::left_join(map_info %>% dplyr::rename(t1 = type_code, s1 = SYMBOL), by = c("a1_id" = "anchor_id")) %>%
    dplyr::left_join(map_info %>% dplyr::rename(t2 = type_code, s2 = SYMBOL), by = c("a2_id" = "anchor_id"))

  loops_annotated$loop_type <- unlist(Map(.loop_type_code, loops_annotated$t1, loops_annotated$t2), use.names = FALSE)
  locus_genes <- unlist(Map(.loop_locus_genes, loops_annotated$t1, loops_annotated$t2, loops_annotated$s1, loops_annotated$s2), use.names = FALSE)
  loops_annotated$single_loop_genes <- locus_genes
  loops_annotated$reg_loop_genes <- locus_genes

  log_message("    Calculating Topology (Hops)...")
  map_info$SYMBOL <- trimws(map_info$SYMBOL)
  valid_pg_nodes <- map_info %>% dplyr::filter(type_code %in% c("P", "G") & !is.na(SYMBOL) & SYMBOL != "")
  lookup_pg_symbol <- valid_pg_nodes$SYMBOL
  names(lookup_pg_symbol) <- valid_pg_nodes$anchor_id
  lookup_pg_type <- valid_pg_nodes$type_code
  names(lookup_pg_type) <- valid_pg_nodes$anchor_id
  lookup_p_symbol <- map_info %>%
    dplyr::filter(type_code == "P" & !is.na(SYMBOL) & SYMBOL != "") %>%
    dplyr::pull(SYMBOL)
  names(lookup_p_symbol) <- map_info %>%
    dplyr::filter(type_code == "P" & !is.na(SYMBOL) & SYMBOL != "") %>%
    dplyr::pull(anchor_id)
  nodes_in_graph <- igraph::V(g)$name
  input_hop <- if (is.null(neighbor_hop)) 0 else neighbor_hop
  ego_list_loop <- igraph::ego(g, order = input_hop, nodes = nodes_in_graph, mode = "all")
  names(ego_list_loop) <- nodes_in_graph
  ego_list_target <- igraph::ego(g, order = input_hop + 1, nodes = nodes_in_graph, mode = "all")
  names(ego_list_target) <- nodes_in_graph
  anchor_topo_map <- data.frame(anchor_id = nodes_in_graph, topo_genes_p = vapply(ego_list_loop, function(x) .ids_to_genes_simple(names(x), lookup_p_symbol), character(1)), topo_genes_pg = vapply(ego_list_loop, function(x) .ids_to_genes_simple(names(x), lookup_pg_symbol), character(1)), tgt_genes_pg = vapply(ego_list_target, function(x) .ids_to_genes_simple(names(x), lookup_pg_symbol), character(1)), tgt_genes_p = vapply(ego_list_target, function(x) .ids_to_genes_simple(names(x), lookup_p_symbol), character(1)), tgt_genes_prio = vapply(ego_list_target, function(x) .ids_to_genes_priority(names(x), lookup_pg_symbol, lookup_pg_type), character(1)), stringsAsFactors = FALSE)
  anchor_topo_map[is.na(anchor_topo_map)] <- NA_character_

  log_message("Step 4: Constructing Loop Tables...")
  loops_annotated <- loops_annotated %>%
    dplyr::left_join(anchor_topo_map %>% dplyr::select(anchor_id, pg1 = topo_genes_pg), by = c("a1_id" = "anchor_id")) %>%
    dplyr::left_join(anchor_topo_map %>% dplyr::select(anchor_id, pg2 = topo_genes_pg), by = c("a2_id" = "anchor_id")) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(proximate_loop_gene = dplyr::case_when((!is.na(t1) & t1 == "G" & !is.na(t2) & t2 == "P") ~ extract_genes(pg2), (!is.na(t1) & t1 == "P" & !is.na(t2) & t2 == "G") ~ extract_genes(pg1), TRUE ~ extract_genes(c(pg1, pg2)))) %>%
    dplyr::ungroup()
  clust_vec <- setNames(anchors$cluster_id, anchors$anchor_id)
  loops_annotated$cluster_id <- clust_vec[loops_annotated$a1_id]
  agg_cluster_reg <- loops_annotated %>%
    dplyr::filter(!is.na(cluster_id)) %>%
    dplyr::group_by(cluster_id) %>%
    dplyr::summarise(all_cluster_loop_genes = extract_genes(reg_loop_genes), .groups = "drop")
  loop_annotation_final <- loops_annotated %>%
    dplyr::left_join(agg_cluster_reg, by = "cluster_id") %>%
    dplyr::mutate(loop_ID = paste0("L", seq_len(dplyr::n()))) %>%
    dplyr::select(loop_ID, chr1, start1, end1, chr2, start2, end2, cluster_id, loop_type, anchor1_gene = s1, anchor1_type = t1, anchor2_gene = s2, anchor2_type = t2, all_cluster_loop_genes, single_loop_genes, proximate_loop_gene, a1_id, a2_id) %>%
    dplyr::rename(All_Anchor_Genes = single_loop_genes, Putative_Target_Genes = proximate_loop_gene, Cluster_All_Genes = all_cluster_loop_genes)
  agg_cluster_locus <- loop_annotation_final %>%
    dplyr::filter(!is.na(cluster_id)) %>%
    dplyr::group_by(cluster_id) %>%
    dplyr::summarise(Cluster_Locus_Genes = extract_genes(All_Anchor_Genes), .groups = "drop")
  if (length(cluster_regions) == 0) {
    warning("No cluster regions formed from loop anchors.")
    cluster_info <- data.frame()
  } else {
    gene_annot <- .with_known_upstream_noise_suppressed(
      ChIPseeker::annotatePeak(cluster_regions, TxDb = txdb_obj, tssRegion = tss_region, annoDb = org_db_pkg, verbose = FALSE)
    )
    cluster_info <- format_annotation_columns(as.data.frame(gene_annot))
  }
  if ("GENENAME" %in% colnames(cluster_info)) cluster_info <- cluster_info %>% dplyr::rename(Gene_description = GENENAME)
  cluster_info$cluster_id <- as.character(cluster_info$cluster_id)
  cluster_info <- cluster_info %>% dplyr::left_join(agg_cluster_locus, by = "cluster_id")

  log_message("    Generating Promoter Centric Stats...")
  raw_stats_df <- dplyr::bind_rows(
    loop_annotation_final %>% dplyr::filter(anchor1_type == "P" & !is.na(anchor1_gene)) %>% dplyr::select(Gene = anchor1_gene, Neighbor_Type = anchor2_type, Loop_Type = loop_type) %>% dplyr::mutate(Gene = as.character(Gene)),
    loop_annotation_final %>% dplyr::filter(anchor2_type == "P" & !is.na(anchor2_gene)) %>% dplyr::select(Gene = anchor2_gene, Neighbor_Type = anchor1_type, Loop_Type = loop_type) %>% dplyr::mutate(Gene = as.character(Gene))
  ) %>%
    tidyr::separate_rows(Gene, sep = ";") %>%
    dplyr::mutate(Gene = trimws(Gene)) %>%
    dplyr::filter(Gene != "" & !is.na(Gene)) %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(
      Total_Loops = dplyr::n(),
      n_Linked_Promoters = sum(Neighbor_Type == "P", na.rm = TRUE),
      n_Linked_Distal = sum(Neighbor_Type %in% c("E", "G"), na.rm = TRUE),
      Dominant_Interaction = names(which.max(table(Loop_Type))),
      .groups = "drop"
    )

  if (nrow(raw_stats_df) == 0) {
    promoter_centric_df <- data.frame(
      Gene = character(), Total_Loops = integer(),
      n_Linked_Promoters = integer(), n_Linked_Distal = integer(),
      Dominant_Interaction = character(),
      Is_High_Connectivity_Gene = character(),
      Is_High_Distal_Connectivity_Gene = character(),
      stringsAsFactors = FALSE
    )
  } else {
  final_cutoff <- max(quantile(raw_stats_df$Total_Loops, hub_percentile, na.rm = TRUE), 3)
  distal_cutoff <- max(quantile(raw_stats_df$n_Linked_Distal, hub_percentile, na.rm = TRUE), 2)

  promoter_centric_df <- raw_stats_df %>%
    dplyr::mutate(
      Is_High_Connectivity_Gene = dplyr::if_else(Total_Loops >= final_cutoff, "Yes", "No"),
      Is_High_Distal_Connectivity_Gene = dplyr::if_else(n_Linked_Distal >= distal_cutoff, "Yes", "No")
    ) %>%
    dplyr::arrange(dplyr::desc(n_Linked_Distal))
  }

  log_message("    Generating Distal Element Stats...")
  distal_raw_df <- dplyr::bind_rows(
    loop_annotation_final %>% dplyr::filter(anchor1_type %in% c("E", "G")) %>% dplyr::select(Distal_Anchor_ID = a1_id, Distal_Type = anchor1_type, Neighbor_Gene = anchor2_gene, Neighbor_Type = anchor2_type, Loop_Type = loop_type),
    loop_annotation_final %>% dplyr::filter(anchor2_type %in% c("E", "G")) %>% dplyr::select(Distal_Anchor_ID = a2_id, Distal_Type = anchor2_type, Neighbor_Gene = anchor1_gene, Neighbor_Type = anchor1_type, Loop_Type = loop_type)
  ) %>%
    dplyr::group_by(Distal_Anchor_ID) %>%
    dplyr::summarise(
      Total_Loops = dplyr::n(),
      n_Linked_Promoters = sum(Neighbor_Type == "P", na.rm = TRUE),
      n_Linked_Distal = sum(Neighbor_Type %in% c("E", "G"), na.rm = TRUE),
      Dominant_Interaction = names(which.max(table(Loop_Type))),
      Target_Genes = extract_genes(Neighbor_Gene[Neighbor_Type == "P"]),
      .groups = "drop"
    )

  anchor_coords_map <- anchors %>%
    dplyr::select(anchor_id, chr, start, end, cluster_id) %>%
    dplyr::distinct()

  if (nrow(distal_raw_df) > 0) {
    final_cutoff_dist <- max(quantile(distal_raw_df$Total_Loops, hub_percentile, na.rm = TRUE), 3)
    distal_element_df <- distal_raw_df %>%
      dplyr::left_join(anchor_coords_map, by = c("Distal_Anchor_ID" = "anchor_id")) %>%
      dplyr::mutate(Is_High_Connectivity_Distal_Element = dplyr::if_else(Total_Loops >= final_cutoff_dist, "Yes", "No")) %>%
      dplyr::select(chr, start, end, cluster_id, Total_Loops, n_Linked_Promoters, n_Linked_Distal, Dominant_Interaction, Is_High_Connectivity_Distal_Element, Target_Genes) %>%
      dplyr::arrange(dplyr::desc(Total_Loops))
  } else {
    distal_element_df <- data.frame(
      chr = character(), start = integer(), end = integer(),
      cluster_id = character(), Total_Loops = integer(),
      n_Linked_Promoters = integer(), n_Linked_Distal = integer(),
      Dominant_Interaction = character(),
      Is_High_Connectivity_Distal_Element = character(),
      Target_Genes = character(),
      stringsAsFactors = FALSE
    )
  }

  bed_info <- NULL
  target_connected_loops <- NULL
  target_gene_links <- NULL
  if (!is.null(target_bed)) {
    log_message("Step 5: Integrating Target Annotations...")
    target_res <- .annotate_target_bed(
      target_bed = target_bed,
      txdb_obj = txdb_obj,
      org_db_pkg = org_db_pkg,
      tss_region = tss_region,
      gene_expr_map = gene_expr_map,
      min_expr = min_expr,
      conflict_strategy = conflict_strategy,
      gr_anchors = gr_anchors,
      anchor_topo_map = anchor_topo_map,
      loop_annotation_final = loop_annotation_final,
      map_info = map_info,
      ego_list_target = ego_list_target,
      log_message = log_message
    )
    bed_info <- target_res$bed_info
    target_connected_loops <- target_res$target_connected_loops
    target_gene_links <- target_res$target_gene_links
  }


  # Step 6: Visualization
  log_message("Step 6: Generating Visualizations (Returning plot objects)...")
  plot_df <- loop_annotation_final
  plot_df$loop_genes <- plot_df$All_Anchor_Genes

  plot_list <- if (quiet) {
    .with_messages_silenced(
      .with_known_upstream_noise_suppressed(
        build_annotation_plots(
          plot_df = plot_df,
          bed_info = bed_info,
          cluster_info = cluster_info,
          target_connected_loops = target_connected_loops,
          txdb_obj = txdb_obj,
          org_db_pkg = org_db_pkg,
          species = species,
          project_name = project_name,
          color_palette = color_palette,
          karyo_bin_size = karyo_bin_size
        )
      )
    )
  } else {
    .with_known_upstream_noise_suppressed(
      build_annotation_plots(
        plot_df = plot_df,
        bed_info = bed_info,
        cluster_info = cluster_info,
        target_connected_loops = target_connected_loops,
        txdb_obj = txdb_obj,
        org_db_pkg = org_db_pkg,
        species = species,
        project_name = project_name,
        color_palette = color_palette,
        karyo_bin_size = karyo_bin_size
      )
    )
  }

  loop_annotation_clean <- loop_annotation_final %>% dplyr::select(-any_of(c("a1_id", "a2_id")))
  if (write_output) {
    log_message("Step 7: Exporting to Excel...")
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "Loop Annotation")
    openxlsx::writeData(wb, "Loop Annotation", loop_annotation_clean)
    openxlsx::addWorksheet(wb, "Anchor Loci Annotation")
    openxlsx::writeData(wb, "Anchor Loci Annotation", cluster_info)
    openxlsx::addWorksheet(wb, "Promoter Stats")
    openxlsx::writeData(wb, "Promoter Stats", promoter_centric_df)
    if (!is.null(distal_element_df)) {
      openxlsx::addWorksheet(wb, "Distal Element Stats")
      openxlsx::writeData(wb, "Distal Element Stats", distal_element_df)
    }
    if (!is.null(bed_info)) {
      openxlsx::addWorksheet(wb, "Target Annotation")
      openxlsx::writeData(wb, "Target Annotation", bed_info)
    }
    if (!is.null(target_gene_links) && nrow(target_gene_links) > 0) {
      openxlsx::addWorksheet(wb, "Target Gene Links")
      openxlsx::writeData(wb, "Target Gene Links", target_gene_links)
    }
    tryCatch(
      openxlsx::saveWorkbook(wb, file.path(out_dir, paste0(project_name, "_Basic_Results.xlsx")), overwrite = TRUE),
      error = function(e) warning("Failed to save Excel workbook: ", conditionMessage(e), call. = FALSE)
    )
    log_message("    Excel file saved.")
  }

  log_message("Analysis Complete.")
  return(list(
    anchor_loci_annotation = cluster_info,
    anchor_annotation = cluster_info,
    loop_annotation = loop_annotation_clean,
    promoter_centric_stats = promoter_centric_df,
    distal_element_stats = distal_element_df,
    target_annotation = bed_info,
    target_gene_links = target_gene_links,
    plots = plot_list,
    plot_list = plot_list
  ))
}

.empty_target_gene_links <- function() {
  data.frame(
    input_id = character(),
    loop_ID = character(),
    anchor_id = character(),
    gene = character(),
    gene_role = character(),
    source = character(),
    evidence = character(),
    anchor_role = character(),
    used_as_fallback = logical(),
    in_regulated_promoter = logical(),
    in_assigned_target = logical(),
    in_all_loop_connected = logical(),
    in_regulated_promoter_filled = logical(),
    in_assigned_target_filled = logical(),
    stringsAsFactors = FALSE
  )
}

.target_gene_link_flags <- function(n) {
  data.frame(
    used_as_fallback = rep(FALSE, n),
    in_regulated_promoter = rep(FALSE, n),
    in_assigned_target = rep(FALSE, n),
    in_all_loop_connected = rep(FALSE, n),
    in_regulated_promoter_filled = rep(FALSE, n),
    in_assigned_target_filled = rep(FALSE, n),
    stringsAsFactors = FALSE
  )
}

.target_linear_gene_column <- function(bed_info) {
  if ("SYMBOL" %in% colnames(bed_info)) {
    return("SYMBOL")
  }
  if ("geneId" %in% colnames(bed_info)) {
    return("geneId")
  }
  NULL
}

.collapse_target_values <- function(x, default = NA_character_) {
  x <- unique(trimws(as.character(na.omit(x))))
  x <- x[nzchar(x)]
  if (length(x) == 0) {
    return(default)
  }
  paste(sort(x), collapse = ";")
}

.target_anchor_gene_map <- function(map_info) {
  if (is.null(map_info) || nrow(map_info) == 0 ||
    !all(c("anchor_id", "type_code", "SYMBOL") %in% colnames(map_info))) {
    return(data.frame(
      anchor_id = character(), gene = character(),
      gene_role = character(), stringsAsFactors = FALSE
    ))
  }

  map_info %>%
    dplyr::filter(type_code %in% c("P", "G"), !is.na(SYMBOL), SYMBOL != "") %>%
    dplyr::select(anchor_id, type_code, SYMBOL) %>%
    dplyr::mutate(SYMBOL = as.character(SYMBOL)) %>%
    tidyr::separate_rows(SYMBOL, sep = ";") %>%
    dplyr::mutate(gene = trimws(SYMBOL)) %>%
    dplyr::filter(!is.na(gene), gene != "") %>%
    dplyr::transmute(
      anchor_id = as.character(anchor_id),
      gene = gene,
      gene_role = dplyr::if_else(type_code == "P", "promoter", "gene_body")
    ) %>%
    dplyr::distinct()
}

.linear_target_gene_links <- function(bed_info) {
  empty <- .empty_target_gene_links()[, c(
    "input_id", "loop_ID", "anchor_id", "gene", "gene_role",
    "source", "evidence", "anchor_role"
  )]
  linear_col <- .target_linear_gene_column(bed_info)
  if (is.null(linear_col) || is.null(bed_info) || nrow(bed_info) == 0) {
    return(empty)
  }

  input_ids <- if ("input_id" %in% colnames(bed_info)) {
    as.character(bed_info$input_id)
  } else {
    paste0("Peak_", seq_len(nrow(bed_info)))
  }

  rows <- lapply(seq_len(nrow(bed_info)), function(i) {
    genes <- clean_gene_names(bed_info[[linear_col]][i], ";")
    if (length(genes) == 0) {
      return(NULL)
    }
    data.frame(
      input_id = input_ids[[i]],
      loop_ID = NA_character_,
      anchor_id = NA_character_,
      gene = genes,
      gene_role = "linear_annotation",
      source = "linear_annotation",
      evidence = "linear_annotation",
      anchor_role = "linear_annotation",
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(out)) {
    return(empty)
  }
  out
}

.build_target_gene_links <- function(
  hit_df, bed_info, loop_annotation_final, map_info, ego_list_target
) {
  base_cols <- c(
    "input_id", "loop_ID", "anchor_id", "gene", "gene_role",
    "source", "evidence", "anchor_role"
  )
  empty <- .empty_target_gene_links()[, base_cols]
  gene_map <- .target_anchor_gene_map(map_info)
  linear_rows <- .linear_target_gene_links(bed_info)

  if (is.null(hit_df) || nrow(hit_df) == 0 || nrow(gene_map) == 0) {
    rows <- dplyr::bind_rows(linear_rows)
    if (nrow(rows) == 0) {
      return(.empty_target_gene_links())
    }
    rows <- dplyr::distinct(rows)
    return(cbind(rows, .target_gene_link_flags(nrow(rows))))
  }

  hit_base <- hit_df %>%
    dplyr::filter(!is.na(anchor_id), anchor_id != "") %>%
    dplyr::transmute(
      input_id = paste0("Peak_", qid),
      local_anchor_id = as.character(anchor_id)
    ) %>%
    dplyr::distinct()

  edge_long <- dplyr::bind_rows(
    loop_annotation_final %>% dplyr::select(loop_ID, local_anchor_id = a1_id, opposite_anchor_id = a2_id),
    loop_annotation_final %>% dplyr::select(loop_ID, local_anchor_id = a2_id, opposite_anchor_id = a1_id)
  ) %>%
    dplyr::filter(!is.na(local_anchor_id), !is.na(opposite_anchor_id)) %>%
    dplyr::mutate(
      local_anchor_id = as.character(local_anchor_id),
      opposite_anchor_id = as.character(opposite_anchor_id)
    ) %>%
    dplyr::distinct()

  local_seed <- hit_base %>%
    dplyr::inner_join(edge_long %>% dplyr::select(local_anchor_id, loop_ID) %>% dplyr::distinct(),
      by = "local_anchor_id"
    )
  if (nrow(local_seed) == 0) {
    local_seed <- hit_base %>% dplyr::mutate(loop_ID = NA_character_)
  }

  local_rows <- local_seed %>%
    dplyr::left_join(gene_map, by = c("local_anchor_id" = "anchor_id")) %>%
    dplyr::filter(!is.na(gene), gene != "") %>%
    dplyr::transmute(
      input_id,
      loop_ID,
      anchor_id = local_anchor_id,
      gene,
      gene_role,
      source = "loop_anchor",
      evidence = dplyr::if_else(gene_role == "promoter", "local_promoter_overlap", "gene_body_context"),
      anchor_role = "local_anchor"
    )

  direct_rows <- hit_base %>%
    dplyr::inner_join(edge_long, by = "local_anchor_id") %>%
    dplyr::left_join(gene_map, by = c("opposite_anchor_id" = "anchor_id")) %>%
    dplyr::filter(!is.na(gene), gene != "") %>%
    dplyr::transmute(
      input_id,
      loop_ID,
      anchor_id = opposite_anchor_id,
      gene,
      gene_role,
      source = "loop_anchor",
      evidence = dplyr::if_else(gene_role == "promoter", "direct_opposite_promoter", "gene_body_context"),
      anchor_role = "opposite_anchor"
    )

  expanded_seed <- do.call(rbind, lapply(seq_len(nrow(hit_base)), function(i) {
    local_id <- hit_base$local_anchor_id[[i]]
    ego_ids <- if (local_id %in% names(ego_list_target)) {
      names(ego_list_target[[local_id]])
    } else {
      character(0)
    }
    direct_ids <- edge_long$opposite_anchor_id[edge_long$local_anchor_id == local_id]
    expanded_ids <- setdiff(ego_ids, c(local_id, direct_ids))
    if (length(expanded_ids) == 0) {
      return(NULL)
    }
    data.frame(
      input_id = hit_base$input_id[[i]],
      anchor_id = expanded_ids,
      stringsAsFactors = FALSE
    )
  }))
  expanded_rows <- empty
  if (!is.null(expanded_seed) && nrow(expanded_seed) > 0) {
    expanded_rows <- expanded_seed %>%
      dplyr::left_join(gene_map, by = "anchor_id") %>%
      dplyr::filter(!is.na(gene), gene != "") %>%
      dplyr::transmute(
        input_id,
        loop_ID = NA_character_,
        anchor_id,
        gene,
        gene_role,
        source = "loop_anchor",
        evidence = dplyr::if_else(gene_role == "promoter", "expanded_promoter_loop", "gene_body_context"),
        anchor_role = "expanded_anchor"
      )
  }

  rows <- dplyr::bind_rows(local_rows, direct_rows, expanded_rows, linear_rows) %>%
    dplyr::distinct()
  if (nrow(rows) == 0) {
    return(.empty_target_gene_links())
  }

  cbind(rows, .target_gene_link_flags(nrow(rows)))
}

.summarise_regulated_promoter_evidence <- function(target_gene_links) {
  if (is.null(target_gene_links) || nrow(target_gene_links) == 0) {
    return(data.frame(input_id = character(), Regulated_promoter_Evidence = character()))
  }
  target_gene_links %>%
    dplyr::filter(source == "loop_anchor", gene_role == "promoter") %>%
    dplyr::group_by(input_id) %>%
    dplyr::summarise(
      Regulated_promoter_Evidence = .collapse_target_values(evidence, default = "none"),
      .groups = "drop"
    )
}

.fallback_evidence_from_annotation <- function(annotation) {
  vapply(annotation, function(x) {
    if (is.na(x) || x == "") {
      return("linear_fallback")
    }
    x <- tolower(as.character(x))
    if (grepl("promoter", x)) {
      return("local_promoter")
    }
    if (grepl("exon|intron|utr", x)) {
      return("local_gene_body")
    }
    "linear_nearest"
  }, FUN.VALUE = character(1))
}

.contains_target_gene <- function(gene, gene_string) {
  if (is.na(gene) || gene == "" || is.na(gene_string) || gene_string == "") {
    return(FALSE)
  }
  gene %in% clean_gene_names(gene_string, ";")
}

.mark_target_gene_link_membership <- function(target_gene_links, bed_info) {
  if (is.null(target_gene_links) || nrow(target_gene_links) == 0) {
    return(.empty_target_gene_links())
  }

  target_cols <- c(
    "Regulated_promoter_genes", "Assigned_Target_Genes",
    "All_Loop_Connected_Genes", "Regulated_promoter_genes_Filled",
    "Assigned_Target_Genes_Filled"
  )
  lookup_cols <- intersect(c("input_id", target_cols), colnames(bed_info))
  link_df <- target_gene_links %>%
    dplyr::select(-dplyr::any_of(c(
      "used_as_fallback", "in_regulated_promoter", "in_assigned_target",
      "in_all_loop_connected", "in_regulated_promoter_filled",
      "in_assigned_target_filled"
    ))) %>%
    dplyr::left_join(bed_info[, lookup_cols, drop = FALSE], by = "input_id")

  for (col in target_cols) {
    if (!col %in% colnames(link_df)) {
      link_df[[col]] <- NA_character_
    }
  }

  link_df$in_regulated_promoter <- unlist(Map(.contains_target_gene, link_df$gene, link_df$Regulated_promoter_genes), use.names = FALSE)
  link_df$in_assigned_target <- unlist(Map(.contains_target_gene, link_df$gene, link_df$Assigned_Target_Genes), use.names = FALSE)
  link_df$in_all_loop_connected <- unlist(Map(.contains_target_gene, link_df$gene, link_df$All_Loop_Connected_Genes), use.names = FALSE)
  link_df$in_regulated_promoter_filled <- unlist(Map(.contains_target_gene, link_df$gene, link_df$Regulated_promoter_genes_Filled), use.names = FALSE)
  link_df$in_assigned_target_filled <- unlist(Map(.contains_target_gene, link_df$gene, link_df$Assigned_Target_Genes_Filled), use.names = FALSE)

  link_df$evidence[link_df$source == "linear_annotation"] <- "linear_annotation"
  reg_empty <- is.na(link_df$Regulated_promoter_genes) | link_df$Regulated_promoter_genes == ""
  assigned_empty <- is.na(link_df$Assigned_Target_Genes) | link_df$Assigned_Target_Genes == ""
  link_df$used_as_fallback <- link_df$source == "linear_annotation" &
    ((reg_empty & link_df$in_regulated_promoter_filled) |
      (assigned_empty & link_df$in_assigned_target_filled))
  link_df$evidence[link_df$used_as_fallback] <- "linear_fallback"

  link_df %>%
    dplyr::select(
      input_id, loop_ID, anchor_id, gene, gene_role, source, evidence,
      anchor_role, used_as_fallback, in_regulated_promoter,
      in_assigned_target, in_all_loop_connected,
      in_regulated_promoter_filled, in_assigned_target_filled
    ) %>%
    dplyr::distinct()
}

#' Internal: Filter Target Gene Links After Expression-Aware Refinement
#'
#' Re-marks membership flags against expression-filtered target columns,
#' appends \code{Mean_Expression} and \code{Passes_Expression_Filter},
#' and retains only rows where at least one membership flag is still
#' \code{TRUE}. Evidence labels such as \code{local_promoter_overlap}
#' and \code{linear_fallback} are preserved.
#'
#' @param target_gene_links Data frame from
#'   \code{\link{.build_target_gene_links}}.
#' @param bed_info Target annotation data frame after expression filtering.
#' @param vals Named numeric vector of per-gene mean expression.
#' @param threshold Numeric. Expression threshold.
#' @return A data frame with columns: \code{input_id}, \code{loop_ID},
#'   \code{anchor_id}, \code{gene}, \code{Mean_Expression},
#'   \code{Passes_Expression_Filter}, \code{gene_role}, \code{source},
#'   \code{evidence}, \code{anchor_role}, \code{used_as_fallback},
#'   \code{in_regulated_promoter}, \code{in_assigned_target},
#'   \code{in_all_loop_connected}, \code{in_regulated_promoter_filled},
#'   \code{in_assigned_target_filled}.
#' @keywords internal
#' @noRd
.filter_refined_target_gene_links <- function(target_gene_links, bed_info, vals, threshold) {
  refined_cols <- c(
    "input_id", "loop_ID", "anchor_id", "gene", "Mean_Expression",
    "Passes_Expression_Filter", "gene_role", "source", "evidence",
    "anchor_role", "used_as_fallback", "in_regulated_promoter",
    "in_assigned_target", "in_all_loop_connected",
    "in_regulated_promoter_filled", "in_assigned_target_filled"
  )
  empty_refined <- function() {
    out <- .empty_target_gene_links()
    out$Mean_Expression <- numeric()
    out$Passes_Expression_Filter <- logical()
    out[, refined_cols, drop = FALSE]
  }

  link_df <- .mark_target_gene_link_membership(target_gene_links, bed_info)
  if (nrow(link_df) == 0) {
    return(empty_refined())
  }

  expr <- unname(vals[link_df$gene])
  link_df$Mean_Expression <- as.numeric(expr)
  link_df$Passes_Expression_Filter <- !is.na(link_df$Mean_Expression) &
    link_df$Mean_Expression >= threshold
  link_df$retained_after_refinement <- link_df$in_regulated_promoter |
    link_df$in_assigned_target |
    link_df$in_all_loop_connected |
    link_df$in_regulated_promoter_filled |
    link_df$in_assigned_target_filled

  link_df %>%
    dplyr::filter(retained_after_refinement) %>%
    dplyr::select(dplyr::all_of(refined_cols)) %>%
    dplyr::distinct()
}

#' Internal: Build Karyotype Gene GRanges
#'
#' Shared helper that loads the full gene catalog, maps gene IDs to SYMBOLs,
#' and subsets to a given gene list. Used by karyotype heatmap sections in
#' both \code{\link{build_annotation_plots}} and related refinement plots.
#'
#' @param gene_symbols Character vector of gene symbols to retain.
#' @param txdb_obj A \code{TxDb} object for gene coordinates.
#' @param org_db_pkg Character. OrgDb package name for symbol mapping.
#' @param context Character. Diagnostic label for mapping messages.
#' @return A \code{GRanges} object subset to matching genes.
#' @keywords internal
#' @noRd
.build_karyo_gene_gr <- function(gene_symbols, txdb_obj, org_db_pkg, context) {
  all_genes_gr <- .with_known_upstream_noise_suppressed(
    GenomicFeatures::genes(txdb_obj)
  )
  map <- .map_txdb_gene_ids(
    gene_ids = .extract_txdb_gene_ids(all_genes_gr),
    org_db = org_db_pkg,
    columns = "SYMBOL",
    context = context,
    warn = FALSE
  )
  S4Vectors::mcols(all_genes_gr)$SYMBOL <- map$SYMBOL[match(
    .extract_txdb_gene_ids(all_genes_gr), map$gene_id
  )]
  all_genes_gr[S4Vectors::mcols(all_genes_gr)$SYMBOL %in% gene_symbols]
}

#' Internal: Build Loop Type Donut Plot
#'
#' @param plot_df Loop annotation data frame.
#' @param custom_colors Named color vector keyed by loop_type.
#' @param project_name Character. Project prefix for the plot title.
#' @return A \code{ggplot} object.
#' @keywords internal
#' @noRd
.build_donut_plot <- function(plot_df, custom_colors, project_name) {
  donut_data <- plot_df %>%
    dplyr::group_by(loop_type) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(
      prop = n / sum(n),
      legend_label = paste0(loop_type, " (n=", n, ", ", round(prop * 100, 1), "%)")
    ) %>%
    dplyr::arrange(dplyr::desc(n))

  donut_data$loop_type <- factor(donut_data$loop_type, levels = donut_data$loop_type)

  ggplot2::ggplot(donut_data, ggplot2::aes(x = 2, y = n, fill = loop_type)) +
    ggplot2::geom_bar(stat = "identity", color = "white") +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::xlim(0.5, 2.9) +
    ggplot2::geom_text(ggplot2::aes(x = 2.65, label = loop_type),
      position = ggplot2::position_stack(vjust = 0.5), size = 2.5
    ) +
    ggplot2::scale_fill_manual(
      values = rev(custom_colors),
      labels = setNames(donut_data$legend_label, donut_data$loop_type)
    ) +
    ggplot2::theme_void() +
    ggplot2::labs(title = paste0(project_name, ": Loop Type Distribution")) +
    ggplot2::theme(
      legend.position = "right",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )
}

#' Internal: Build Karyotype Loop-Genes Heatmap
#'
#' @param plot_df Loop annotation data frame with \code{All_Anchor_Genes} column.
#' @param txdb_obj A \code{TxDb} object.
#' @param org_db_pkg Character. OrgDb package name.
#' @param species Character. Genome assembly.
#' @param karyo_bin_size Integer. Bin size for karyotype heatmaps.
#' @param red_palette Character vector. Red color palette.
#' @return A karyotype grob, or \code{NULL}.
#' @keywords internal
#' @noRd
.build_karyo_loop_genes_plot <- function(
  plot_df, txdb_obj, org_db_pkg, species, karyo_bin_size, red_palette
) {
  if (!"All_Anchor_Genes" %in% colnames(plot_df)) {
    return(NULL)
  }
  genes_loop <- clean_gene_names(plot_df$All_Anchor_Genes, ";")
  if (length(genes_loop) == 0) {
    return(NULL)
  }
  target_genes_gr <- .build_karyo_gene_gr(
    genes_loop, txdb_obj, org_db_pkg,
    "build_annotation_plots loop-gene karyotype"
  )
  draw_karyo_heatmap_internal(
    target_genes_gr, "Loop Genes Distribution", karyo_bin_size, 0.99,
    txdb_obj, species, "Genes",
    custom_colors = red_palette
  )
}

#' Internal: Build Karyotype Anchor-Load Heatmap
#'
#' @param plot_df Loop annotation data frame with coordinate columns.
#' @param txdb_obj A \code{TxDb} object.
#' @param species Character. Genome assembly.
#' @param karyo_bin_size Integer. Bin size for karyotype heatmaps.
#' @param blue_palette Character vector. Blue color palette.
#' @return A karyotype grob, or \code{NULL}.
#' @keywords internal
#' @noRd
.build_karyo_anchors_plot <- function(
  plot_df, txdb_obj, species, karyo_bin_size, blue_palette
) {
  all_anchors <- dplyr::bind_rows(
    plot_df %>% dplyr::select(chr = chr1, start = start1, end = end1),
    plot_df %>% dplyr::select(chr = chr2, start = start2, end = end2)
  ) %>% dplyr::distinct()
  if (nrow(all_anchors) == 0) {
    return(NULL)
  }
  draw_karyo_heatmap_internal(
    .with_known_upstream_noise_suppressed(
      GenomicRanges::makeGRangesFromDataFrame(all_anchors)
    ),
    "Loop Anchor Load", karyo_bin_size, 0.99, txdb_obj, species, "Anchors",
    custom_colors = blue_palette
  )
}

#' Internal: Build Simplified Flower Plot for Annotation
#'
#' @param plot_df Loop annotation data frame.
#' @param project_name Character. Project prefix for the plot title.
#' @param custom_colors Named color vector keyed by loop_type.
#' @return A \code{ggplot} object, or \code{NULL} if fewer than 2 gene sets.
#' @keywords internal
#' @noRd
.build_flower_plot <- function(plot_df, project_name, custom_colors) {
  temp_df_flower <- plot_df %>%
    dplyr::filter(!is.na(All_Anchor_Genes) & All_Anchor_Genes != "") %>%
    tidyr::separate_rows(All_Anchor_Genes, sep = ";") %>%
    dplyr::mutate(All_Anchor_Genes = trimws(All_Anchor_Genes)) %>%
    dplyr::filter(All_Anchor_Genes != "")
  gene_sets <- split(temp_df_flower$All_Anchor_Genes, temp_df_flower$loop_type)
  gene_sets <- lapply(gene_sets, unique)
  if (length(gene_sets) <= 1) {
    return(NULL)
  }
  draw_flower_simplified(gene_sets, project_name, custom_colors)
}

#' Internal: Build Target-Connected Rose Plot
#'
#' @param target_connected_loops Data frame of target-connected loops.
#' @param custom_colors Named color vector keyed by loop_type.
#' @param project_name Character. Project prefix for the plot title.
#' @return A \code{ggplot} object, or \code{NULL}.
#' @keywords internal
#' @noRd
.build_target_rose_plot <- function(
  target_connected_loops, custom_colors, project_name
) {
  if (is.null(target_connected_loops) || nrow(target_connected_loops) == 0) {
    return(NULL)
  }
  rose_data <- target_connected_loops %>%
    dplyr::group_by(loop_type) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(
      prop = n / sum(n),
      legend_label = paste0(loop_type, " (n=", n, ", ", round(prop * 100, 1), "%)")
    ) %>%
    dplyr::arrange(dplyr::desc(n))

  rose_data$loop_type <- factor(rose_data$loop_type, levels = rose_data$loop_type)

  ggplot2::ggplot(rose_data, ggplot2::aes(x = loop_type, y = n, fill = loop_type)) +
    ggplot2::geom_bar(stat = "identity", width = 1, color = "white") +
    ggplot2::coord_polar(theta = "x") +
    ggplot2::scale_fill_manual(
      values = custom_colors,
      labels = setNames(rose_data$legend_label, rose_data$loop_type)
    ) +
    ggplot2::theme_void() +
    ggplot2::labs(title = paste0(project_name, ": Target-Linked Loop Types")) +
    ggplot2::theme(
      legend.position = "right",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )
}

#' Internal: Build Target Genomic Distribution Pie Charts
#'
#' @param bed_info Target annotation data frame.
#' @param color_palette Character. RColorBrewer palette name.
#' @param project_name Character. Project prefix for plot titles.
#' @return A named list of \code{ggplot} objects (may be empty).
#' @keywords internal
#' @noRd
.build_target_genomic_pies <- function(bed_info, color_palette, project_name) {
  plot_list <- list()
  if (is.null(bed_info)) {
    return(plot_list)
  }
  plot_list$Target_Genomic_Distribution <- draw_pie_with_outside_labels(
    bed_info, "annotation",
    paste0(project_name, ": All Targets Genomic Distribution"), color_palette
  )
  if ("Linked_Loop_IDs" %in% colnames(bed_info)) {
    linked_rows <- bed_info[!is.na(bed_info$Linked_Loop_IDs) &
      bed_info$Linked_Loop_IDs != "", ]
    if (nrow(linked_rows) > 0) {
      plot_list$Target_Loop_Genomic_Distribution <- draw_pie_with_outside_labels(
        linked_rows, "annotation",
        paste0(project_name, ": Loop-Connected Targets Distribution"),
        color_palette
      )
    }
  }
  plot_list
}

#' Internal: Build Karyotype Target-Genes Heatmap
#'
#' @param bed_info Target annotation data frame (or NULL).
#' @param txdb_obj A \code{TxDb} object.
#' @param org_db_pkg Character. OrgDb package name.
#' @param species Character. Genome assembly.
#' @param karyo_bin_size Integer. Bin size for karyotype heatmaps.
#' @param purple_palette Character vector. Purple color palette.
#' @return A karyotype grob, or \code{NULL}.
#' @keywords internal
#' @noRd
.build_karyo_target_genes_plot <- function(
  bed_info, txdb_obj, org_db_pkg, species,
  karyo_bin_size, purple_palette
) {
  if (is.null(bed_info) ||
    !"Assigned_Target_Genes_Filled" %in% colnames(bed_info)) {
    return(NULL)
  }
  genes_target <- clean_gene_names(
    bed_info$Assigned_Target_Genes_Filled, ";"
  )
  if (length(genes_target) == 0) {
    return(NULL)
  }
  target_genes_gr <- .build_karyo_gene_gr(
    genes_target, txdb_obj, org_db_pkg,
    "build_annotation_plots target-gene karyotype"
  )
  draw_karyo_heatmap_internal(
    target_genes_gr, "Target Genes (Assigned+Local)", karyo_bin_size,
    0.99, txdb_obj, species, "Genes",
    custom_colors = purple_palette
  )
}

#' Internal: Build Annotation Visualization Suite
#'
#' Generates the complete diagnostic plot collection (donut, rose, karyotype
#' heatmaps, flower plot, pie charts) for the annotation results.
#'
#' @param plot_df Loop annotation data frame with loop_type and All_Anchor_Genes.
#' @param bed_info Target annotation data frame (optional).
#' @param cluster_info Anchor annotation data frame.
#' @param target_connected_loops Data frame of target-connected loops (optional).
#' @param txdb_obj TxDb object for gene coordinate lookup.
#' @param org_db_pkg Character. Organism database package name.
#' @param species Character. Genome assembly.
#' @param project_name Character. Project prefix for plot titles.
#' @param color_palette Character. RColorBrewer palette name.
#' @param karyo_bin_size Integer. Bin size for karyotype heatmaps.
#' @return A named list of ggplot / grob objects.
#' @keywords internal
#' @noRd
build_annotation_plots <- function(
  plot_df, bed_info, cluster_info,
  target_connected_loops, txdb_obj, org_db_pkg, species, project_name,
  color_palette, karyo_bin_size
) {
  loop_types_sorted <- sort(unique(plot_df$loop_type))
  custom_colors <- get_colors(length(loop_types_sorted), color_palette)
  names(custom_colors) <- loop_types_sorted

  red_palette <- c("#FFFFFF", "#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C", "#FC4E2A", "#E31A1C", "#BD0026", "#800026", "#000000")
  blue_palette <- c("#FFFFFF", "#E1F5FE", "#B3E5FC", "#4FC3F7", "#039BE5", "#0277BD", "#01579B", "#000000")
  purple_palette <- c("#FFFFFF", "#F3E5F5", "#E1BEE7", "#BA68C8", "#9C27B0", "#7B1FA2", "#4A148C", "#000000")

  plot_list <- list()

  plot_list$Basic_Donut <- .build_donut_plot(plot_df, custom_colors, project_name)

  plot_list$Basic_Circular <- draw_circular_bar_plot(plot_df, project_name,
    color_vec = custom_colors
  )

  plot_list$Karyo_LoopGenes <- .build_karyo_loop_genes_plot(
    plot_df, txdb_obj, org_db_pkg, species, karyo_bin_size, red_palette
  )

  plot_list$Karyo_Anchors <- .build_karyo_anchors_plot(
    plot_df, txdb_obj, species, karyo_bin_size, blue_palette
  )

  plot_list$Basic_Flower <- .build_flower_plot(
    plot_df, project_name, custom_colors
  )

  if (nrow(cluster_info) > 0) {
    plot_list$Anchor_Genomic_Distribution <- draw_pie_with_outside_labels(
      cluster_info, "annotation",
      paste0(project_name, ": Anchor Loci Genomic Distribution"), color_palette
    )
  }

  if (!is.null(bed_info)) {
    plot_list$Karyo_TargetGenes <- .build_karyo_target_genes_plot(
      bed_info, txdb_obj, org_db_pkg, species, karyo_bin_size, purple_palette
    )

    plot_list$Target_Rose <- .build_target_rose_plot(
      target_connected_loops, custom_colors, project_name
    )

    target_pies <- .build_target_genomic_pies(bed_info, color_palette, project_name)
    for (nm in names(target_pies)) plot_list[[nm]] <- target_pies[[nm]]
  }

  return(plot_list)
}

# --- Internal helpers for refine_loop_anchors_by_expression ---

#' Reclassify anchor types based on expression whitelist.
#' @keywords internal
#' @noRd
.refine_reclassify_anchors <- function(loop_df, whitelist, reclassify_by_expression) {
  a1_res <- Map(
    function(g, t) clean_anchor(g, t, allow = whitelist, down = reclassify_by_expression),
    loop_df$anchor1_gene, loop_df$anchor1_type
  )
  a2_res <- Map(
    function(g, t) clean_anchor(g, t, allow = whitelist, down = reclassify_by_expression),
    loop_df$anchor2_gene, loop_df$anchor2_type
  )
  loop_df$anchor1_type <- vapply(a1_res, function(x) x$type, character(1))
  loop_df$anchor1_gene <- vapply(a1_res, function(x) x$gene, character(1))
  loop_df$anchor2_type <- vapply(a2_res, function(x) x$type, character(1))
  loop_df$anchor2_gene <- vapply(a2_res, function(x) x$gene, character(1))
  loop_df$loop_type <- unlist(Map(
    function(t1, t2) paste(sort(c(t1, t2)), collapse = "-"),
    loop_df$anchor1_type, loop_df$anchor2_type
  ), use.names = FALSE)
  loop_df
}

#' Compute active/fallback target genes and refinement status flags.
#' @keywords internal
#' @noRd
.refine_compute_targets <- function(loop_df, original_ptg, whitelist, orig_anchor1_type = NULL, orig_anchor2_type = NULL) {
  filter_genes_wl <- function(x) {
    if (is.na(x) || x == "") {
      return(NA_character_)
    }
    gs <- trimws(unlist(strsplit(as.character(x), ";")))
    gs <- unique(gs[gs %in% whitelist])
    if (length(gs) == 0) {
      return(NA_character_)
    }
    paste(sort(gs), collapse = ";")
  }
  filtered_ptg <- vapply(original_ptg, filter_genes_wl, character(1))

  is_enh_like <- function(t) t %in% c("E", "eP", "eG")
  is_promoter <- function(t) t == "P"
  is_gene_body <- function(t) t == "G"

  loop_df <- loop_df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      .fallback_ptg = dplyr::case_when(
        (is_promoter(anchor1_type) & is_enh_like(anchor2_type)) ~ extract_genes(anchor1_gene),
        (is_enh_like(anchor1_type) & is_promoter(anchor2_type)) ~ extract_genes(anchor2_gene),
        (is_promoter(anchor1_type) & is_promoter(anchor2_type)) ~ extract_genes(c(anchor1_gene, anchor2_gene)),
        (is_promoter(anchor1_type) & is_gene_body(anchor2_type)) ~ extract_genes(anchor1_gene),
        (is_gene_body(anchor1_type) & is_promoter(anchor2_type)) ~ extract_genes(anchor2_gene),
        (is_gene_body(anchor1_type) & is_enh_like(anchor2_type)) ~ extract_genes(anchor1_gene),
        (is_enh_like(anchor1_type) & is_gene_body(anchor2_type)) ~ extract_genes(anchor2_gene),
        TRUE ~ extract_genes(c(anchor1_gene, anchor2_gene))
      )
    ) %>%
    dplyr::ungroup()

  loop_df$Active_Target_Genes <- filtered_ptg
  loop_df$Putative_Target_Genes <- filtered_ptg
  empty_idx <- is.na(loop_df$Putative_Target_Genes) | loop_df$Putative_Target_Genes == ""
  loop_df$Putative_Target_Genes[empty_idx] <- loop_df$.fallback_ptg[empty_idx]
  loop_df <- loop_df %>% dplyr::select(-.fallback_ptg)

  has_active_target <- !is.na(loop_df$Active_Target_Genes) & loop_df$Active_Target_Genes != ""
  if (!is.null(orig_anchor1_type) && !is.null(orig_anchor2_type)) {
    reclassified_a1 <- orig_anchor1_type != loop_df$anchor1_type & loop_df$anchor1_type %in% c("eP", "eG")
    reclassified_a2 <- orig_anchor2_type != loop_df$anchor2_type & loop_df$anchor2_type %in% c("eP", "eG")
    has_reclassified <- reclassified_a1 | reclassified_a2
  } else {
    has_reclassified <- loop_df$anchor1_type %in% c("eP", "eG") | loop_df$anchor2_type %in% c("eP", "eG")
  }
  had_original <- !is.na(original_ptg) & original_ptg != ""

  loop_df$Has_Active_Target <- has_active_target
  loop_df$Refinement_Action <- dplyr::case_when(
    has_active_target ~ "retained_active_target",
    !has_active_target & has_reclassified & had_original ~ "reclassified_silent_anchor",
    !has_active_target & had_original ~ "expression_filtered_no_active_target",
    TRUE ~ "structural_only_no_active_target"
  )
  loop_df$Retained_In_Functional_Network <- has_active_target
  loop_df
}

#' Refine target BED annotations by expression filtering.
#' @keywords internal
#' @noRd
.refine_target_annotations <- function(bed_info, loop_df, whitelist, target_gene_links, vals, threshold) {
  cols_to_clean <- grep("Strict|Physical|Loop_Genes|promoter|Filled|Target_Genes|Assigned", colnames(bed_info), value = TRUE)
  cols_to_clean <- cols_to_clean[!grepl("Evidence|Linked_Loop_IDs|SANKEY_RAW_GENES", cols_to_clean)]
  raw_tgt_col <- "Assigned_Target_Genes_Filled"
  if (!raw_tgt_col %in% colnames(bed_info)) {
    raw_tgt_col <- grep("Filled", cols_to_clean, value = TRUE)[1]
  }
  if (!is.na(raw_tgt_col) && raw_tgt_col %in% colnames(bed_info)) bed_info$SANKEY_RAW_GENES <- bed_info[[raw_tgt_col]]
  for (col in cols_to_clean) {
    if (col %in% colnames(bed_info)) {
      bed_info[[col]] <- vapply(as.character(bed_info[[col]]), function(x) {
        if (is.na(x) || x == "") {
          return(NA_character_)
        }
        gs <- unlist(strsplit(x, ";"))
        gs_active <- gs[trimws(gs) %in% whitelist]
        if (length(gs_active) == 0) {
          return(NA_character_)
        }
        paste(unique(sort(trimws(gs_active))), collapse = ";")
      }, FUN.VALUE = character(1))
    }
  }

  any_sym_in_wl <- function(s) {
    if (is.na(s) || s == "") {
      return(FALSE)
    }
    any(trimws(unlist(strsplit(as.character(s), ";"))) %in% whitelist)
  }
  filter_sym_expressed <- function(s) {
    if (is.na(s) || s == "") {
      return(NA_character_)
    }
    gs <- trimws(unlist(strsplit(as.character(s), ";")))
    gs <- gs[gs %in% whitelist]
    if (length(gs) == 0) {
      return(NA_character_)
    }
    paste(sort(gs), collapse = ";")
  }
  linear_col <- .target_linear_gene_column(bed_info)
  if (!is.null(linear_col)) {
    has_sym <- vapply(bed_info[[linear_col]], any_sym_in_wl, logical(1))
    sym_fill <- vapply(bed_info[[linear_col]], filter_sym_expressed, character(1))
  } else {
    has_sym <- rep(FALSE, nrow(bed_info))
    sym_fill <- rep(NA_character_, nrow(bed_info))
  }

  if ("Regulated_promoter_genes" %in% colnames(bed_info) &&
    "Regulated_promoter_genes_Filled" %in% colnames(bed_info)) {
    has_reg <- !is.na(bed_info$Regulated_promoter_genes) & bed_info$Regulated_promoter_genes != ""
    bed_info$Regulated_promoter_genes_Filled <- dplyr::case_when(
      has_reg ~ bed_info$Regulated_promoter_genes,
      has_sym ~ sym_fill,
      TRUE ~ NA_character_
    )
    if ("Regulated_promoter_Fallback_Evidence" %in% colnames(bed_info)) {
      bed_info$Regulated_promoter_Fallback_Evidence <- dplyr::case_when(
        has_reg ~ "none",
        has_sym ~ bed_info$Regulated_promoter_Fallback_Evidence,
        TRUE ~ "none"
      )
    }
  }

  if ("Linked_Loop_IDs" %in% colnames(bed_info) &&
    "loop_ID" %in% colnames(loop_df) &&
    "Active_Target_Genes" %in% colnames(loop_df)) {
    loop_tgt <- loop_df %>%
      dplyr::filter(!is.na(Active_Target_Genes) & Active_Target_Genes != "") %>%
      dplyr::select(loop_ID, Active_Target_Genes) %>%
      dplyr::distinct()
    get_loop_tgt <- function(linked) {
      if (is.na(linked) || linked == "") {
        return(NA_character_)
      }
      ids <- trimws(unlist(strsplit(as.character(linked), ";")))
      tgt <- loop_tgt$Active_Target_Genes[match(ids, loop_tgt$loop_ID)]
      genes <- unique(trimws(unlist(strsplit(na.omit(tgt), ";"))))
      genes <- genes[genes != ""]
      if (length(genes) == 0) {
        return(NA_character_)
      }
      paste(sort(genes), collapse = ";")
    }
    loop_tgt_vec <- vapply(bed_info$Linked_Loop_IDs, get_loop_tgt, FUN.VALUE = character(1))
    if ("Assigned_Target_Genes" %in% colnames(bed_info)) {
      empty <- is.na(bed_info$Assigned_Target_Genes) | bed_info$Assigned_Target_Genes == ""
      fill_ok <- !is.na(loop_tgt_vec) & loop_tgt_vec != ""
      bed_info$Assigned_Target_Genes[empty & fill_ok] <- loop_tgt_vec[empty & fill_ok]
      if ("Assigned_Target_Genes_Filled" %in% colnames(bed_info)) {
        has_tgt <- !is.na(bed_info$Assigned_Target_Genes) & bed_info$Assigned_Target_Genes != ""
        bed_info$Assigned_Target_Genes_Filled <- dplyr::case_when(
          has_tgt ~ bed_info$Assigned_Target_Genes,
          has_sym ~ sym_fill,
          TRUE ~ NA_character_
        )
      }
    }
  }

  if (!is.null(target_gene_links)) {
    target_gene_links <- .filter_refined_target_gene_links(
      target_gene_links, bed_info, vals, threshold
    )
  }

  list(bed_info = bed_info, target_gene_links = target_gene_links)
}

#' Export refined results to Excel workbook.
#' @keywords internal
#' @noRd
.refine_export_workbook <- function(loop_df, clust_info, promoter_centric_df, distal_element_df, bed_info, target_gene_links, out_dir, project_name) {
  wb <- openxlsx::createWorkbook()
  loop_export <- loop_df %>%
    dplyr::select(-any_of(c("a1_id", "a2_id", "loop_genes", "single_loop_genes", "proximate_loop_gene")))

  openxlsx::addWorksheet(wb, "Filtered Loop Annotation")
  openxlsx::writeData(wb, "Filtered Loop Annotation", loop_export)

  functional_loops <- loop_export %>% dplyr::filter(Retained_In_Functional_Network == TRUE)
  openxlsx::addWorksheet(wb, "Functional Loop Annotation")
  openxlsx::writeData(wb, "Functional Loop Annotation", functional_loops)

  openxlsx::addWorksheet(wb, "Filtered Anchor Loci Annotation")
  openxlsx::writeData(wb, "Filtered Anchor Loci Annotation", clust_info)

  if (!is.null(promoter_centric_df)) {
    openxlsx::addWorksheet(wb, "Filtered Promoter Stats")
    openxlsx::writeData(wb, "Filtered Promoter Stats", promoter_centric_df)
  }
  if (!is.null(distal_element_df)) {
    openxlsx::addWorksheet(wb, "Filtered Distal Element Stats")
    openxlsx::writeData(wb, "Filtered Distal Element Stats", distal_element_df)
  }
  if (!is.null(bed_info)) {
    bed_export <- bed_info %>% dplyr::select(-any_of("SANKEY_RAW_GENES"))
    openxlsx::addWorksheet(wb, "Filtered Target Annotation")
    openxlsx::writeData(wb, "Filtered Target Annotation", bed_export)
  }
  if (!is.null(target_gene_links) && nrow(target_gene_links) > 0) {
    openxlsx::addWorksheet(wb, "Filtered Target Gene Links")
    openxlsx::writeData(wb, "Filtered Target Gene Links", target_gene_links)
  }

  tryCatch(
    openxlsx::saveWorkbook(wb, file.path(out_dir, paste0(project_name, "_Refined_Results.xlsx")), overwrite = TRUE),
    error = function(e) warning("Failed to save refined Excel workbook: ", conditionMessage(e), call. = FALSE)
  )
}

#' @title Expression-Aware refinement of loop anchors and target linkages
#'
#' @description
#' Integrates quantitative RNA-seq data (e.g., TPM/FPKM) with 3D structural data
#' to reclassify regulatory elements and annotate functional status, deriving a
#' functionally interpretable regulatory network from physical chromatin contacts.
#' All structural loops are preserved; refinement status columns indicate which
#' loops belong to the high-confidence active subset.
#'
#' @details
#' \strong{Algorithmic Framework:}
#' \itemize{
#'   \item \strong{Target Filtration:} Parses merged gene assignments (e.g., \code{"GeneA;GeneB"}), evaluates individual genes against a defined expression threshold, and retains only transcriptionally active targets.
#'   \item \strong{Biological Reclassification:} Reclassifies physically annotated promoters (\code{P}) and gene bodies (\code{G}) lacking active transcription as enhancer-like regulatory elements (\code{eP}, \code{eG}). This adjusts the regulatory syntax to reflect functional states (e.g., reannotating a silent \code{P-P} loop to an \code{eP-P} interaction).
#'   \item \strong{Expression-Aware Connectivity Statistics:} Recomputes promoter-centric and distal-element connectivity after expression-aware anchor refinement, while preserving all structural loops in the refined loop annotation. This separates the complete physical contact map from the high-confidence active subset.
#'   \item \strong{External Target Refinement:} Filters auxiliary target mapping columns (e.g., \code{Assigned_Target_Genes_Filled}) based on expression criteria, ensuring that mapped 1D genomic features are exclusively linked to active genes.
#'   \item \strong{Target Provenance Preservation:} Recomputes \code{*_Filled}
#'   membership flags in \code{target_gene_links} after expression filtering,
#'   retains only links still used by the refined target columns, and appends
#'   \code{Mean_Expression} plus \code{Passes_Expression_Filter}. Evidence labels
#'   such as \code{local_promoter_overlap}, \code{direct_opposite_promoter}, and
#'   \code{linear_fallback} are preserved.
#' }
#'
#' \strong{Design Philosophy:}
#' This function does not discard structural loops based on expression state.
#' Hi-C, HiChIP, and PLAC-seq capture 3D chromatin contacts; RNA-seq captures
#' current transcriptional state. A silent promoter may reflect cell-state,
#' time-point, or technical factors rather than absence of physical contact.
#' All structural loops are retained with refinement status columns, and a
#' high-confidence functional subset is provided via
#' \code{Retained_In_Functional_Network} and the \emph{Functional Loop Annotation}
#' Excel sheet.
#'
#' \strong{Interpretation of eP/eG labels:}
#' The \code{eP} and \code{eG} labels capture expression-aware enhancer-like
#' regulatory states, enabling \code{looplook} to distinguish transcriptionally
#' silent reference promoters or gene bodies from putative regulatory anchors in
#' 3D chromatin space. Orthogonal chromatin evidence, including ATAC-seq
#' accessibility, H3K27ac enrichment, or H3K27me3 depletion, can further
#' strengthen biological interpretation when available. Users holding matched
#' ATAC-seq or ChIP-seq data may overlay eP/eG loci with these tracks to
#' confirm residual regulatory activity at transcriptionally silent promoters.
#'
#' @param annotation_res List. The raw foundational output object returned by \code{\link{annotate_peaks_and_loops}}.
#' @param expr_matrix_file Path to a normalised expression matrix (TPM/FPKM, genes × samples). Required for refinement. Default: \code{NULL}.
#' @param sample_columns Character vector or integer indices. Columns in \code{expr_matrix_file} to average. Default: \code{NULL}.
#' @param threshold Numeric. Minimum expression (e.g. TPM >= 1) for a gene to be considered active. Default: \code{1}.
#' @param unit_type Character. Expression unit for plot labels (e.g., \code{"TPM"}). Default: \code{"TPM"}.
#' @param species Character. Genome assembly. One of \code{"hg38"}, \code{"hg19"}, \code{"mm10"}, \code{"mm9"}. Default: \code{"hg38"}.
#' @param out_dir Character. Output directory for the Excel results file. Default: \code{"./results/filtered"}.
#' @param project_name Character. Prefix for output files (automatically appends \code{"_Filtered"}). Default: \code{"HiChIP"}.
#' @param color_palette Character. RColorBrewer palette name for loop-type colour assignments. Default: \code{"Paired"}.
#' @param karyo_bin_size Integer. Bin width in bp for karyotype heatmaps. Default: \code{1e5}.
#' @param reclassify_by_expression Logical. If \code{TRUE} (default), silent promoters (P) and gene bodies (G) are reclassified as eP/eG.
#' @param hub_percentile Numeric (0–1). Node-degree quantile for hub detection. Default: \code{0.95}.
#' @param write_output Logical. If \code{TRUE} (default), write the refined Excel workbook to \code{out_dir}. If \code{FALSE}, return results without creating directories or files.
#' @param quiet Logical. If \code{TRUE}, suppress progress messages while preserving warnings. Default: \code{FALSE}.
#'
#' @return An invisible named list:
#' \itemize{
#'   \item \code{loop_annotation} — Full refined 3D network with updated \code{loop_type}
#'     (e.g., eP-P) and two target gene columns:
#'     \itemize{
#'       \item \code{Active_Target_Genes}: Expression-filtered active-only targets (no fallback).
#'       \item \code{Putative_Target_Genes}: Display column; may include linear nearest-gene
#'         fallback when \code{Active_Target_Genes} is empty.
#'     }
#'     Refinement status columns:
#'     \code{Has_Active_Target}, \code{Retained_In_Functional_Network}, and
#'     \code{Refinement_Action} (\code{"retained_active_target"},
#'     \code{"reclassified_silent_anchor"}, \code{"expression_filtered_no_active_target"},
#'     or \code{"structural_only_no_active_target"}).
#'     All structural loops are preserved; filter on \code{Retained_In_Functional_Network}
#'     for the high-confidence active subset.
#'   \item \code{anchor_loci_annotation} — Filtered non-redundant anchor-locus annotations with expressed targets.
#'   \item \code{anchor_annotation} — Backward-compatible alias of \code{anchor_loci_annotation}.
#'   \item \code{promoter_centric_stats} — Gene-level connectivity statistics.
#'   \item \code{distal_element_stats} — Distal-element connectivity statistics.
#'   \item \code{target_annotation} — Target features (peaks) with gene assignments.
#'     Key columns include:
#'     \itemize{
#'       \item \code{All_Loop_Connected_Genes}: All genes from loop-connected anchors (P/G types).
#'       \item \code{Regulated_promoter_genes}: Promoter genes supported by loop-anchor context.
#'       \item \code{Assigned_Target_Genes}: Promoter-first 3D assignment (prioritises P > G > E).
#'       \item \code{*_Filled} variants: Linear nearest-gene fallback when strict 3D assignments are empty.
#'       \item \code{Regulated_promoter_Evidence}: Provenance of \code{Regulated_promoter_genes}
#'         (e.g., \code{local_promoter_overlap}, \code{direct_opposite_promoter}).
#'         \strong{Read with} \code{Regulated_promoter_genes}; do not cross-reference
#'         with \code{Assigned_Target_Genes} or other columns.
#'       \item \code{Regulated_promoter_Fallback_Evidence}: Provenance of
#'         \code{Regulated_promoter_genes_Filled}.
#'         \strong{Read with} \code{Regulated_promoter_genes_Filled}; indicates
#'         which \code{*_Filled} column supplied the fallback gene.
#'     }
#'   \item \code{target_gene_links} — Long-format peak-gene provenance table.
#'     Each row records one peak-gene linkage with full provenance.
#'     \strong{Read} \code{evidence}, \code{anchor_role}, and \code{gene_role}
#'     \strong{together as a group} — they jointly describe how each gene was
#'     assigned to each peak; do not interpret any one column in isolation.
#'     \itemize{
#'       \item \code{input_id}, \code{loop_ID}, \code{anchor_id}: Identifiers.
#'       \item \code{gene}: Linked gene symbol.
#'       \item \code{gene_role}: \code{"promoter"}, \code{"gene_body"}, or \code{"linear_annotation"}.
#'       \item \code{source}: \code{"loop_anchor"} (3D-derived) or \code{"linear_annotation"} (nearest gene).
#'       \item \code{evidence}: Provenance label —
#'         \code{"local_promoter_overlap"} (peak overlaps anchor promoter),
#'         \code{"direct_opposite_promoter"} (opposite anchor is promoter),
#'         \code{"gene_body_context"} (gene body linkage),
#'         \code{"expanded_promoter_loop"} (via ego-network expansion),
#'         \code{"linear_annotation"} (direct nearest gene),
#'         or \code{"linear_fallback"} (filled when 3D assignment was empty).
#'       \item \code{anchor_role}: \code{"local_anchor"}, \code{"opposite_anchor"},
#'         \code{"expanded_anchor"}, or \code{"linear_annotation"}.
#'       \item \code{used_as_fallback}: Logical. \code{TRUE} when this link was added
#'         via the \code{*_Filled} linear nearest-gene fallback mechanism.
#'       \item \code{in_regulated_promoter} through \code{in_assigned_target_filled}:
#'         Logical membership flags indicating which target annotation column(s)
#'         this gene appears in. A gene may appear in multiple columns simultaneously.
#'       \item (Refine only) \code{Mean_Expression}: Per-gene mean expression value.
#'       \item (Refine only) \code{Passes_Expression_Filter}: Logical. \code{TRUE} if
#'         \code{Mean_Expression >= threshold}.
#'     }
#'   \item \code{plots} — Named list of ggplot objects (dumbbell, rose, karyotype).
#'   \item \code{plot_list} — Backward-compatible alias of \code{plots}.
#' }
#' If \code{write_output = TRUE}, also writes \code{_Refined_Results.xlsx} to \code{out_dir}.
#' The workbook contains a \emph{Functional Loop Annotation} sheet with only
#' loops where \code{Retained_In_Functional_Network == TRUE}.
#'
#' @importFrom dplyr %>% filter group_by summarise ungroup mutate select rename left_join full_join arrange desc case_when rowwise coalesce any_of distinct pull
#' @importFrom ggplot2 ggplot aes geom_bar geom_segment geom_point geom_text scale_color_manual scale_fill_manual theme_minimal theme_void labs coord_polar xlim
#' @importFrom tidyr pivot_longer separate_rows
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps
#' @importFrom openxlsx createWorkbook addWorksheet writeData saveWorkbook
#' @export
#'
#' @examples
#' # 1. Get paths to the required example files in the package
#' rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
#' expr_path <- system.file("extdata", "example_tpm.txt", package = "looplook")
#'
#' # Safely load the pre-computed annotation result from RData
#' temp_env <- new.env()
#' load(rdata_path, envir = temp_env)
#' raw_annotation <- temp_env[[ls(temp_env)[1]]]
#' raw_annotation$loop_annotation <- head(raw_annotation$loop_annotation, 6)
#' raw_annotation$target_annotation <- head(raw_annotation$target_annotation, 3)
#' raw_annotation$promoter_centric_stats <- head(raw_annotation$promoter_centric_stats, 6)
#' raw_annotation$distal_element_stats <- head(raw_annotation$distal_element_stats, 6)
#'
#' res_reclassified <- refine_loop_anchors_by_expression(
#'   annotation_res = raw_annotation,
#'   expr_matrix_file = expr_path,
#'   sample_columns = "con1",
#'   threshold = 1.0,
#'   unit_type = "TPM",
#'   species = "hg38",
#'   out_dir = tempdir(),
#'   project_name = "Example_Reclassified",
#'   reclassify_by_expression = TRUE,
#'   write_output = FALSE,
#'   quiet = TRUE
#' )
#' print(table(res_reclassified$loop_annotation$loop_type))
refine_loop_anchors_by_expression <- function(
  annotation_res,
  expr_matrix_file = NULL,
  sample_columns = NULL,
  threshold = 1,
  unit_type = "TPM",
  species = "hg38",
  out_dir = "./results/filtered",
  project_name = "HiChIP",
  color_palette = "Paired",
  karyo_bin_size = 1e5,
  reclassify_by_expression = TRUE,
  hub_percentile = 0.95,
  write_output = TRUE,
  quiet = FALSE
) {
  species <- match.arg(species, c("hg38", "hg19", "mm10", "mm9"))
  log_message <- function(...) {
    if (!quiet) message(...)
  }

  # --- 0. Setup ---
  if (!grepl("_Filtered$", project_name)) project_name <- paste0(project_name, "_Filtered")
  if (write_output && !dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  log_message(">>> [Refinement] Project Name: ", project_name)

  # --- 1. Load Data ---
  log_message(">>> [Step 1] Loading Data & Expression Matrix...")
  if (is.null(annotation_res$loop_annotation)) stop("'loop_annotation' missing.")

  original_loop_df <- annotation_res$loop_annotation
  loop_df <- annotation_res$loop_annotation
  clust_info <- annotation_res$anchor_loci_annotation
  if (is.null(clust_info)) {
    clust_info <- annotation_res$anchor_annotation
  }
  bed_info <- annotation_res$target_annotation
  target_gene_links <- annotation_res$target_gene_links

  # Validate required columns in loop_annotation
  required_cols <- c(
    "chr1", "start1", "end1", "chr2", "start2", "end2",
    "anchor1_gene", "anchor1_type", "anchor2_gene", "anchor2_type",
    "Putative_Target_Genes"
  )
  missing_cols <- setdiff(required_cols, colnames(loop_df))
  if (length(missing_cols) > 0) {
    stop(
      "annotation_res$loop_annotation is missing required columns: ",
      paste(missing_cols, collapse = ", "), ". ",
      "Ensure the input was generated by annotate_peaks_and_loops()."
    )
  }
  if (!"cluster_id" %in% colnames(loop_df)) {
    warning("'cluster_id' not found in loop_annotation; cluster-level stats and target-donut will be skipped.", call. = FALSE)
    loop_df$cluster_id <- NA_character_
    has_cluster_id <- FALSE
  } else {
    has_cluster_id <- TRUE
  }

  # Check and reconstruct IDs if missing
  if (!"a1_id" %in% colnames(loop_df) || !"a2_id" %in% colnames(loop_df)) {
    log_message("    Reconstructing anchor IDs from coordinates (omitted from upstream loop_annotation output).")
    loop_df <- loop_df %>%
      dplyr::mutate(
        a1_id = paste(chr1, start1, end1, sep = "_"),
        a2_id = paste(chr2, start2, end2, sep = "_")
      )
  }

  upstream_promoter_stats <- annotation_res$promoter_centric_stats

  if (is.null(expr_matrix_file) || is.null(sample_columns)) {
    stop("Expression matrix file and sample columns are required for refinement.")
  }
  vals <- load_expression_matrix(expr_matrix_file, sample_columns)
  whitelist <- names(vals)[vals >= threshold & !is.na(vals) & names(vals) != ""]
  log_message(sprintf("    >>> Active Genes (>= %s %s): %d", threshold, unit_type, length(whitelist)))

  anno_genes <- unique(c(
    trimws(unlist(strsplit(na.omit(loop_df$anchor1_gene), ";"))),
    trimws(unlist(strsplit(na.omit(loop_df$anchor2_gene), ";")))
  ))
  anno_genes <- anno_genes[nzchar(anno_genes)]
  if (length(anno_genes) > 0) {
    overlap_rate <- length(intersect(whitelist, anno_genes)) / length(anno_genes)
    if (overlap_rate < 0.1) {
      warning(
        sprintf(
          "Only %.1f%% of annotation gene symbols match the expression matrix row names. ",
          overlap_rate * 100
        ),
        "Check that expression matrix row names use the same gene identifier convention (e.g., SYMBOL)."
      )
    }
  }

  # --- 2. Update Anchors & Loops ---
  log_message(">>> [Step 2] Updating Anchors & Loops...")

  orig_anchor1_type <- loop_df$anchor1_type
  orig_anchor2_type <- loop_df$anchor2_type
  loop_df <- .refine_reclassify_anchors(loop_df, whitelist, reclassify_by_expression)
  original_ptg <- loop_df$Putative_Target_Genes

  loop_df <- .refine_compute_targets(loop_df, original_ptg, whitelist, orig_anchor1_type, orig_anchor2_type)

  log_message(sprintf(
    "    Retained: %d / %d loops",
    sum(loop_df$Has_Active_Target), nrow(loop_df)
  ))

  # --- 3. Stats Update ---
  log_message(">>> [Step 3] Updating Stats...")

  if (has_cluster_id) {
    agg_cluster <- loop_df %>%
      dplyr::filter(!is.na(cluster_id)) %>%
      dplyr::group_by(cluster_id) %>%
      dplyr::summarise(
        Cluster_All_Genes = extract_genes(Putative_Target_Genes),
        Cluster_Active_Target_Genes = extract_genes(Active_Target_Genes),
        .groups = "drop"
      )
    loop_df <- loop_df %>%
      dplyr::select(-any_of(c("Cluster_All_Genes", "Cluster_Active_Target_Genes"))) %>%
      dplyr::left_join(agg_cluster, by = "cluster_id")
    if (!is.null(clust_info) && "cluster_id" %in% colnames(clust_info)) {
      clust_info <- clust_info %>%
        dplyr::select(-any_of(c("Cluster_All_Genes", "Cluster_Active_Target_Genes"))) %>%
        dplyr::left_join(agg_cluster, by = "cluster_id")
    }
  }

  stats_res <- compute_refined_stats(
    loop_df = loop_df,
    upstream_promoter_stats = upstream_promoter_stats,
    vals = vals,
    threshold = threshold,
    hub_percentile = hub_percentile
  )
  promoter_centric_df <- stats_res$promoter_centric
  distal_element_df <- stats_res$distal_element

  log_message(">>> [Step 4] Refining Target Annotations...")
  if (!is.null(bed_info)) {
    tgt_refined <- .refine_target_annotations(
      bed_info, loop_df, whitelist, target_gene_links, vals, threshold
    )
    bed_info <- tgt_refined$bed_info
    target_gene_links <- tgt_refined$target_gene_links
  }

  # --- 5. Visualization ---
  log_message(">>> [Step 5] Generating Visualizations (Returning plot objects)...")
  plot_list <- if (quiet) {
    .with_messages_silenced(
      .with_known_upstream_noise_suppressed(
        build_refinement_plots(
          original_loop_df = original_loop_df,
          loop_df = loop_df,
          bed_info = bed_info,
          whitelist = whitelist,
          project_name = project_name,
          karyo_bin_size = karyo_bin_size,
          species = species,
          color_palette = color_palette
        )
      )
    )
  } else {
    .with_known_upstream_noise_suppressed(
      build_refinement_plots(
        original_loop_df = original_loop_df,
        loop_df = loop_df,
        bed_info = bed_info,
        whitelist = whitelist,
        project_name = project_name,
        karyo_bin_size = karyo_bin_size,
        species = species,
        color_palette = color_palette
      )
    )
  }

  # --- 6. Export ---
  if (write_output) {
    log_message(">>> [Step 6] Exporting Refined Results...")
    .refine_export_workbook(
      loop_df, clust_info, promoter_centric_df, distal_element_df,
      bed_info, target_gene_links, out_dir, project_name
    )
    log_message("    Excel saved.")
  }

  log_message("Refinement Complete.")
  return(list(
    loop_annotation = loop_df,
    anchor_loci_annotation = clust_info,
    anchor_annotation = clust_info,
    promoter_centric_stats = promoter_centric_df,
    distal_element_stats = distal_element_df,
    target_annotation = bed_info,
    target_gene_links = target_gene_links,
    plots = plot_list,
    plot_list = plot_list
  ))
}

#' Internal: Build Refinement Dumbbell Comparison Plot
#'
#' Compares loop-type counts before vs. after expression-aware filtering.
#'
#' @param original_loop_df Loop annotation before refinement.
#' @param loop_df Loop annotation after refinement.
#' @param project_name Character. Project prefix for the plot title.
#' @return A \code{ggplot} object.
#' @keywords internal
#' @noRd
.build_dumbbell_plot <- function(original_loop_df, loop_df, project_name) {
  df_orig <- original_loop_df %>%
    dplyr::group_by(loop_type) %>%
    dplyr::summarise(Original = dplyr::n(), .groups = "drop")
  df_filt <- loop_df %>%
    dplyr::group_by(loop_type) %>%
    dplyr::summarise(Refined = dplyr::n(), .groups = "drop")
  df_dumbbell <- dplyr::full_join(df_orig, df_filt, by = "loop_type") %>%
    dplyr::mutate(
      Original = ifelse(is.na(Original), 0, Original),
      Refined = ifelse(is.na(Refined), 0, Refined),
      is_e_type = grepl("e", loop_type)
    ) %>%
    dplyr::arrange(is_e_type, dplyr::desc(Original))
  df_dumbbell$loop_type <- factor(df_dumbbell$loop_type,
    levels = rev(df_dumbbell$loop_type)
  )
  df_long <- df_dumbbell %>%
    tidyr::pivot_longer(
      cols = c("Original", "Refined"),
      names_to = "Source", values_to = "Count"
    )

  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = df_dumbbell,
      ggplot2::aes(y = loop_type, yend = loop_type, x = Original, xend = Refined),
      color = "#b2b2b2", linewidth = 0.8
    ) +
    ggplot2::geom_point(
      data = df_long,
      ggplot2::aes(x = Count, y = loop_type, color = Source), size = 3
    ) +
    ggplot2::scale_color_manual(values = c("Original" = "#999999", "Refined" = "#E69F00")) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = paste0(project_name, ": Refinement Effect (Dumbbell)"),
      x = "Number of Loops", y = "Loop Type"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "top"
    )
}

#' Internal: Build Refinement Target-Loop Donut Plot
#'
#' Donut chart of loop-type distribution among target-connected refined loops.
#'
#' @param bed_info Target annotation data frame (or NULL).
#' @param loop_df Refined loop annotation.
#' @param custom_colors Named color vector keyed by loop_type.
#' @param project_name Character. Project prefix for the plot title.
#' @return A \code{ggplot} object, or \code{NULL} if no target-connected loops.
#' @keywords internal
#' @noRd
.build_refinement_donut <- function(bed_info, loop_df, custom_colors, project_name) {
  if (is.null(bed_info)) {
    return(NULL)
  }
  gr_bed <- GenomicRanges::makeGRangesFromDataFrame(bed_info,
    keep.extra.columns = TRUE
  )
  if (!"cluster_id" %in% colnames(loop_df) || all(is.na(loop_df$cluster_id))) {
    return(NULL)
  }
  active_anc <- dplyr::bind_rows(
    loop_df %>% dplyr::select(chr = chr1, start = start1, end = end1, cluster_id),
    loop_df %>% dplyr::select(chr = chr2, start = start2, end = end2, cluster_id)
  ) %>%
    dplyr::filter(!is.na(cluster_id)) %>%
    dplyr::distinct()
  if (nrow(active_anc) == 0) {
    return(NULL)
  }
  gr_anc <- GenomicRanges::makeGRangesFromDataFrame(active_anc,
    keep.extra.columns = TRUE
  )
  hits <- GenomicRanges::findOverlaps(gr_bed, gr_anc)
  if (length(hits) == 0) {
    return(NULL)
  }
  hit_ids <- unique(gr_anc$cluster_id[S4Vectors::subjectHits(hits)])
  tgt_loops <- loop_df %>% dplyr::filter(cluster_id %in% hit_ids)
  if (nrow(tgt_loops) == 0) {
    return(NULL)
  }
  donut_data <- tgt_loops %>%
    dplyr::group_by(loop_type) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(
      fraction = count / sum(count),
      legend_label = paste0(
        loop_type, " (n=", count, ", ",
        round(fraction * 100, 1), "%)"
      ),
      plot_label = loop_type, is_lower_e = grepl("^e", loop_type)
    ) %>%
    dplyr::arrange(is_lower_e, dplyr::desc(count)) %>%
    dplyr::mutate(loop_type = factor(loop_type, levels = loop_type))
  ggplot2::ggplot(
    donut_data,
    ggplot2::aes(x = 2, y = count, fill = loop_type)
  ) +
    ggplot2::geom_bar(stat = "identity", width = 1, color = "white") +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::xlim(0.5, 2.9) +
    ggplot2::geom_text(ggplot2::aes(x = 2.8, label = plot_label),
      position = ggplot2::position_stack(vjust = 0.5), size = 3
    ) +
    ggplot2::scale_fill_manual(
      values = custom_colors,
      labels = setNames(donut_data$legend_label, donut_data$loop_type)
    ) +
    ggplot2::theme_void() +
    ggplot2::labs(title = paste0(project_name, ": Refined Structural Loops at Target Regions")) +
    ggplot2::theme(
      legend.position = "right",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )
}

#' Internal: Build Multi-Omics Sankey Tracking Plot
#'
#' Three-tier Sankey diagram tracing genomic features through loop
#' connection status to expression activity.
#'
#' @param bed_info Target annotation data frame.
#' @param whitelist Character vector of active gene symbols.
#' @param project_name Character. Project prefix for the plot title.
#' @return A \code{networkD3} sankey widget, or \code{NULL} if requirements
#'   are not met.
#' @keywords internal
#' @noRd
.build_sankey_plot <- function(bed_info, whitelist, project_name) {
  if (is.null(bed_info)) {
    return(NULL)
  }
  if (!"SANKEY_RAW_GENES" %in% colnames(bed_info)) {
    return(NULL)
  }
  if (!requireNamespace("networkD3", quietly = TRUE)) {
    return(NULL)
  }

  get_label_mapping <- function(vec, total_targets) {
    tab <- table(vec)
    label_map <- list()
    for (i in seq_along(tab)) {
      name <- names(tab)[i]
      ct <- as.integer(tab[i])
      pct <- round(ct / total_targets * 100, 1)
      label_map[[name]] <- sprintf("%s (n=%d, %.1f%%)", name, ct, pct)
    }
    label_map
  }
  check_status_strict <- function(s, wl) {
    if (is.na(s) || s == "") {
      return("Inactive")
    }
    gs <- trimws(unlist(strsplit(as.character(s), ";")))
    gs <- gs[gs != ""]
    if (length(gs) == 0) {
      return("Inactive")
    }
    if (any(gs %in% wl)) {
      return("Active")
    }
    return("Inactive")
  }

  bed_info$L1_Raw <- vapply(bed_info$annotation, function(a) {
    if (grepl("Intergenic", a, ignore.case = TRUE)) {
      return("Distal Intergenic")
    }
    if (grepl("Promoter", a, ignore.case = TRUE)) {
      return("Promoter")
    }
    if (grepl("Exon", a, ignore.case = TRUE)) {
      return("Exon")
    }
    if (grepl("Intron", a, ignore.case = TRUE)) {
      return("Intron")
    }
    return("Others")
  }, FUN.VALUE = character(1))

  bed_info$L2_Raw <- ifelse(
    !is.na(bed_info$Linked_Loop_IDs) & bed_info$Linked_Loop_IDs != "",
    "Connected", "Unconnected"
  )

  bed_info$L3_Raw <- vapply(bed_info$SANKEY_RAW_GENES, check_status_strict,
    wl = whitelist, FUN.VALUE = character(1)
  )

  sankey_df <- bed_info %>%
    dplyr::filter(
      !is.na(L3_Raw) & L3_Raw != "" &
        !is.na(L1_Raw) & L1_Raw != "" &
        !is.na(L2_Raw) & L2_Raw != ""
    )
  if (nrow(sankey_df) == 0) {
    return(NULL)
  }

  total_targets <- nrow(sankey_df)
  label_map_l1 <- get_label_mapping(sankey_df$L1_Raw, total_targets)
  label_map_l2 <- get_label_mapping(sankey_df$L2_Raw, total_targets)
  label_map_l3 <- get_label_mapping(sankey_df$L3_Raw, total_targets)

  sankey_df$L1_Label <- label_map_l1[sankey_df$L1_Raw]
  sankey_df$L2_Label <- label_map_l2[sankey_df$L2_Raw]
  sankey_df$L3_Label <- label_map_l3[sankey_df$L3_Raw]

  sankey_df <- sankey_df %>%
    dplyr::filter(
      !is.na(.data$L1_Label) & .data$L1_Label != "" &
        !is.na(.data$L2_Label) & .data$L2_Label != "" &
        !is.na(.data$L3_Label) & .data$L3_Label != ""
    )
  if (nrow(sankey_df) == 0) {
    return(NULL)
  }

  nodes <- unique(c(
    unlist(sankey_df$L1_Label, use.names = FALSE),
    unlist(sankey_df$L2_Label, use.names = FALSE),
    unlist(sankey_df$L3_Label, use.names = FALSE)
  ))
  nodes <- data.frame(name = nodes, stringsAsFactors = FALSE)

  get_idx <- function(label) match(label, nodes$name) - 1
  links <- data.frame(
    source = get_idx(sankey_df$L1_Label),
    target = get_idx(sankey_df$L2_Label),
    value = 1,
    stringsAsFactors = FALSE
  )
  links2 <- data.frame(
    source = get_idx(sankey_df$L2_Label),
    target = get_idx(sankey_df$L3_Label),
    value = 1,
    stringsAsFactors = FALSE
  )
  links <- rbind(links, links2)
  links <- links %>%
    dplyr::group_by(.data$source, .data$target) %>%
    dplyr::summarise(value = dplyr::n(), .groups = "drop")

  sankey_colors <- get_colors(nrow(nodes), "Paired")
  color_scale <- paste0('d3.scaleOrdinal().range(["',
    paste(sankey_colors, collapse = '","'), '"])')
  sn <- networkD3::sankeyNetwork(
    Links = links, Nodes = nodes,
    Source = "source", Target = "target",
    Value = "value", NodeID = "name",
    units = "TWh", fontSize = 12, nodeWidth = 30,
    colourScale = networkD3::JS(color_scale),
    iterations = 0
  )

  sn$sizingPolicy$defaultWidth <- "100%"
  sn$sizingPolicy$defaultHeight <- "450px"
  sn <- htmlwidgets::onRender(sn, sprintf('
	function(el, x) {
	  var svg = d3.select(el).select("svg");
	  function createValidID(name) {
	    if (!name) return "unknown";
	    return name.replace(/[^a-zA-Z0-9-]/g, "_");
	  }
	  svg.selectAll(".link").each(function(d) {
	    var gradientID = "gradient-" + createValidID(d.source.name) +
	      "-" + createValidID(d.target.name);
	    if (svg.select("#" + gradientID).empty()) {
	      var gradient = svg.append("defs")
	        .append("linearGradient")
	        .attr("id", gradientID)
	        .attr("gradientUnits", "userSpaceOnUse")
	        .attr("x1", d.source.x + d.source.dx / 2)
	        .attr("y1", d.source.y + d.source.dy / 2)
	        .attr("x2", d.target.x + d.target.dx / 2)
	        .attr("y2", d.target.y + d.target.dy / 2);
	      var sourceColor = d3.select(el).selectAll(".node")
	        .filter(function(node) { return node.name === d.source.name; })
	        .select("rect").style("fill");
	      var targetColor = d3.select(el).selectAll(".node")
	        .filter(function(node) { return node.name === d.target.name; })
	        .select("rect").style("fill");
	      gradient.append("stop").attr("offset", "0%%")
	        .attr("stop-color", sourceColor);
	      gradient.append("stop").attr("offset", "100%%")
	        .attr("stop-color", targetColor);
	    }
	    d3.select(this).style("stroke", "url(#" + gradientID + ")")
	      .style("stroke-opacity", 0.6)
	      .style("stroke-width", function(d) { return Math.max(2, d.width); });
	  });
	  svg.selectAll(".node rect")
	    .style("stroke", "#333333")
	    .style("stroke-width", "1px");
	  svg.selectAll("text")
	    .style("font-size", "12px")
	    .style("font-weight", "bold");
	}
	'))
  sn
}

#' Internal: Build Refinement Karyotype Heatmaps
#'
#' Generates \code{Refined_Karyo_Active} and \code{Refined_Karyo_TargetGenes}
#' karyotype heatmaps from refined loop and target annotation data.
#'
#' @param loop_df Refined loop annotation.
#' @param bed_info Target annotation data frame (or NULL).
#' @param species Character. Genome assembly.
#' @param karyo_bin_size Integer. Bin size for karyotype heatmaps.
#' @param red_palette Character vector. Red color palette.
#' @param purple_palette Character vector. Purple color palette.
#' @return A named list of karyotype grob objects (may be empty).
#' @keywords internal
#' @noRd
.build_refinement_karyotypes <- function(
  loop_df, bed_info, species,
  karyo_bin_size, red_palette, purple_palette
) {
  plot_list <- list()
  txdb_pkg <- tryCatch(species_txdb_pkg(species), error = function(e) NULL)
  org_db <- tryCatch(species_orgdb_pkg(species), error = function(e) NULL)
  if (is.null(txdb_pkg) || is.null(org_db) ||
    !requireNamespace(txdb_pkg, quietly = TRUE) ||
    !requireNamespace(org_db, quietly = TRUE)) {
    return(plot_list)
  }
  txdb_obj <- utils::getFromNamespace(txdb_pkg, txdb_pkg)

  if ("Active_Target_Genes" %in% colnames(loop_df)) {
    g_active <- clean_gene_names(loop_df$Active_Target_Genes, ";")
    if (length(g_active) > 0) {
      all_genes_gr <- .with_known_upstream_noise_suppressed(
        GenomicFeatures::genes(txdb_obj)
      )
      map <- .map_txdb_gene_ids(
        gene_ids = .extract_txdb_gene_ids(all_genes_gr),
        org_db = org_db,
        columns = "SYMBOL",
        context = "build_refinement_plots active-gene karyotype",
        warn = FALSE
      )
      S4Vectors::mcols(all_genes_gr)$SYMBOL <- map$SYMBOL[match(
        .extract_txdb_gene_ids(all_genes_gr), map$gene_id
      )]
      target_genes_gr <- all_genes_gr[
        S4Vectors::mcols(all_genes_gr)$SYMBOL %in% g_active
      ]
      plot_list$Refined_Karyo_Active <- draw_karyo_heatmap_internal(
        target_genes_gr,
        "Refined Active Genes", karyo_bin_size, 0.99, txdb_obj, species,
        "Genes",
        custom_colors = red_palette
      )
    }
  }

  if (!is.null(bed_info) &&
    "Assigned_Target_Genes_Filled" %in% colnames(bed_info)) {
    g_target <- clean_gene_names(bed_info$Assigned_Target_Genes_Filled, ";")
    if (length(g_target) > 0) {
      all_genes_gr <- .with_known_upstream_noise_suppressed(
        GenomicFeatures::genes(txdb_obj)
      )
      map <- .map_txdb_gene_ids(
        gene_ids = .extract_txdb_gene_ids(all_genes_gr),
        org_db = org_db,
        columns = "SYMBOL",
        context = "build_refinement_plots target-gene karyotype",
        warn = FALSE
      )
      S4Vectors::mcols(all_genes_gr)$SYMBOL <- map$SYMBOL[match(
        .extract_txdb_gene_ids(all_genes_gr), map$gene_id
      )]
      target_genes_gr <- all_genes_gr[
        S4Vectors::mcols(all_genes_gr)$SYMBOL %in% g_target
      ]
      plot_list$Refined_Karyo_TargetGenes <- draw_karyo_heatmap_internal(
        target_genes_gr,
        "Refined Target Genes", karyo_bin_size, 0.99, txdb_obj, species,
        "Genes",
        custom_colors = purple_palette
      )
    }
  }
  plot_list
}

#' Internal: Build Refinement Rose Plot
#'
#' Polar bar chart of refined loop-type distribution.
#'
#' @param loop_df Refined loop annotation.
#' @param custom_colors Named color vector keyed by loop_type.
#' @param project_name Character. Project prefix for the plot title.
#' @return A \code{ggplot} object.
#' @keywords internal
#' @noRd
.build_rose_plot <- function(loop_df, custom_colors, project_name) {
  rose_data <- loop_df %>%
    dplyr::group_by(loop_type) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(
      fraction = count / sum(count),
      legend_label = paste0(loop_type, " (n=", count, ", ",
        round(fraction * 100, 1), "%)"),
      is_lower_e = grepl("^e", loop_type)
    ) %>%
    dplyr::arrange(dplyr::desc(count))

  plot_order <- rose_data$loop_type
  rose_data$loop_type <- factor(rose_data$loop_type, levels = plot_order)
  legend_order <- rose_data %>%
    dplyr::arrange(is_lower_e, dplyr::desc(count)) %>%
    dplyr::pull(loop_type)

  ggplot2::ggplot(rose_data, ggplot2::aes(x = loop_type, y = count, fill = loop_type)) +
    ggplot2::geom_bar(stat = "identity", width = 1, color = "white") +
    ggplot2::coord_polar(theta = "x") +
    ggplot2::scale_fill_manual(
      values = custom_colors,
      labels = setNames(rose_data$legend_label, as.character(rose_data$loop_type)),
      breaks = legend_order,
      name = "Loop Type"
    ) +
    ggplot2::theme_void() +
    ggplot2::labs(title = paste0(project_name, ": Structural Loop Types After Reclassification"), subtitle = "Full refined network (all loops). For active subset, filter on Retained_In_Functional_Network.") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "right",
      legend.text = ggplot2::element_text(size = 10)
    )
}
#' Internal: Build Refinement Visualization Suite
#'
#' Generates diagnostic plots for expression-aware refinement: dumbbell
#' comparison, donut, Sankey, karyotype heatmaps, and rose plot.
#'
#' @param original_loop_df Loop annotation before refinement.
#' @param loop_df Loop annotation after refinement.
#' @param bed_info Target annotation data frame (optional).
#' @param whitelist Character vector of active gene symbols.
#' @param project_name Character. Project prefix for plot titles.
#' @param karyo_bin_size Integer. Bin size for karyotype heatmaps.
#' @param species Character. Genome assembly.
#' @return A named list of ggplot / htmlwidget / grob objects.
#' @keywords internal
#' @noRd
build_refinement_plots <- function(
  original_loop_df, loop_df, bed_info,
  whitelist, project_name, karyo_bin_size, species,
  color_palette = "Paired"
) {
  red_palette <- c("#FFFFFF", "#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C", "#FC4E2A", "#E31A1C", "#BD0026", "#800026", "#000000")
  purple_palette <- c("#FFFFFF", "#F3E5F5", "#E1BEE7", "#BA68C8", "#9C27B0", "#7B1FA2", "#4A148C", "#000000")

  # Assign colours by descending loop-type frequency (not alphabetical)
  type_counts <- loop_df %>%
    dplyr::group_by(loop_type) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(n))
  all_types <- type_counts$loop_type
  custom_colors <- get_colors(length(all_types), color_palette)
  names(custom_colors) <- all_types

  plot_list <- list()

  plot_list$Comparison_Dumbbell <- .build_dumbbell_plot(
    original_loop_df, loop_df, project_name
  )

  plot_list$Target_Loop_Donut <- .build_refinement_donut(
    bed_info, loop_df, custom_colors, project_name
  )

  plot_list$Target_Sankey <- .build_sankey_plot(
    bed_info, whitelist, project_name
  )

  karyo_plots <- .build_refinement_karyotypes(
    loop_df, bed_info, species, karyo_bin_size,
    red_palette, purple_palette
  )
  for (nm in names(karyo_plots)) plot_list[[nm]] <- karyo_plots[[nm]]

  plot_list$Rose <- .build_rose_plot(loop_df, custom_colors, project_name)

  return(plot_list)
}

#' Internal: Build Refinement Karyotype Heatmaps
#'
#' Generates \code{Refined_Karyo_Active} and \code{Refined_Karyo_TargetGenes}
#' karyotype heatmaps from refined loop and target annotation data.
#'
#' @param loop_df Refined loop annotation.
#' @param bed_info Target annotation data frame (or NULL).
#' @param species Character. Genome assembly.
#' @param karyo_bin_size Integer. Bin size for karyotype heatmaps.
#' @param red_palette Character vector. Red color palette.
#' @param purple_palette Character vector. Purple color palette.
#' @return A named list of karyotype grob objects (may be empty).
#' @keywords internal
#' @noRd
.build_refinement_karyotypes <- function(
  loop_df, bed_info, species,
  karyo_bin_size, red_palette, purple_palette
) {
  plot_list <- list()
  txdb_pkg <- tryCatch(species_txdb_pkg(species), error = function(e) NULL)
  org_db <- tryCatch(species_orgdb_pkg(species), error = function(e) NULL)
  if (is.null(txdb_pkg) || is.null(org_db) ||
    !requireNamespace(txdb_pkg, quietly = TRUE) ||
    !requireNamespace(org_db, quietly = TRUE)) {
    return(plot_list)
  }
  txdb_obj <- utils::getFromNamespace(txdb_pkg, txdb_pkg)

  if ("Active_Target_Genes" %in% colnames(loop_df)) {
    g_active <- clean_gene_names(loop_df$Active_Target_Genes, ";")
    if (length(g_active) > 0) {
      all_genes_gr <- .with_known_upstream_noise_suppressed(
        GenomicFeatures::genes(txdb_obj)
      )
      map <- .map_txdb_gene_ids(
        gene_ids = .extract_txdb_gene_ids(all_genes_gr),
        org_db = org_db,
        columns = "SYMBOL",
        context = "build_refinement_plots active-gene karyotype",
        warn = FALSE
      )
      S4Vectors::mcols(all_genes_gr)$SYMBOL <- map$SYMBOL[match(
        .extract_txdb_gene_ids(all_genes_gr), map$gene_id
      )]
      target_genes_gr <- all_genes_gr[
        S4Vectors::mcols(all_genes_gr)$SYMBOL %in% g_active
      ]
      plot_list$Refined_Karyo_Active <- draw_karyo_heatmap_internal(
        target_genes_gr,
        "Refined Active Genes", karyo_bin_size, 0.99, txdb_obj, species,
        "Genes",
        custom_colors = red_palette
      )
    }
  }

  if (!is.null(bed_info) &&
    "Assigned_Target_Genes_Filled" %in% colnames(bed_info)) {
    g_target <- clean_gene_names(bed_info$Assigned_Target_Genes_Filled, ";")
    if (length(g_target) > 0) {
      all_genes_gr <- .with_known_upstream_noise_suppressed(
        GenomicFeatures::genes(txdb_obj)
      )
      map <- .map_txdb_gene_ids(
        gene_ids = .extract_txdb_gene_ids(all_genes_gr),
        org_db = org_db,
        columns = "SYMBOL",
        context = "build_refinement_plots target-gene karyotype",
        warn = FALSE
      )
      S4Vectors::mcols(all_genes_gr)$SYMBOL <- map$SYMBOL[match(
        .extract_txdb_gene_ids(all_genes_gr), map$gene_id
      )]
      target_genes_gr <- all_genes_gr[
        S4Vectors::mcols(all_genes_gr)$SYMBOL %in% g_target
      ]
      plot_list$Refined_Karyo_TargetGenes <- draw_karyo_heatmap_internal(
        target_genes_gr,
        "Refined Target Genes", karyo_bin_size, 0.99, txdb_obj, species,
        "Genes",
        custom_colors = purple_palette
      )
    }
  }
  plot_list
}

#' Internal: Build Refinement Rose Plot
#'
#' Polar bar chart of refined loop-type distribution.
#'
#' @param loop_df Refined loop annotation.
#' @param custom_colors Named color vector keyed by loop_type.
#' @param project_name Character. Project prefix for the plot title.
#' @return A \code{ggplot} object.
#' @keywords internal
#' @noRd
.build_rose_plot <- function(loop_df, custom_colors, project_name) {
  rose_data <- loop_df %>%
    dplyr::group_by(loop_type) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(
      fraction = count / sum(count),
      legend_label = paste0(loop_type, " (n=", count, ", ",
        round(fraction * 100, 1), "%)"),
      is_lower_e = grepl("^e", loop_type)
    ) %>%
    dplyr::arrange(dplyr::desc(count))

  plot_order <- rose_data$loop_type
  rose_data$loop_type <- factor(rose_data$loop_type, levels = plot_order)
  legend_order <- rose_data %>%
    dplyr::arrange(is_lower_e, dplyr::desc(count)) %>%
    dplyr::pull(loop_type)

  ggplot2::ggplot(rose_data, ggplot2::aes(x = loop_type, y = count, fill = loop_type)) +
    ggplot2::geom_bar(stat = "identity", width = 1, color = "white") +
    ggplot2::coord_polar(theta = "x") +
    ggplot2::scale_fill_manual(
      values = custom_colors,
      labels = setNames(rose_data$legend_label, as.character(rose_data$loop_type)),
      breaks = legend_order,
      name = "Loop Type"
    ) +
    ggplot2::theme_void() +
    ggplot2::labs(title = paste0(project_name, ": Structural Loop Types After Reclassification"), subtitle = "Full refined network (all loops). For active subset, filter on Retained_In_Functional_Network.") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "right",
      legend.text = ggplot2::element_text(size = 10)
    )
}
#' Internal: Build Refinement Visualization Suite
#'
#' Generates diagnostic plots for expression-aware refinement: dumbbell
#' comparison, donut, Sankey, karyotype heatmaps, and rose plot.
#'
#' @param original_loop_df Loop annotation before refinement.
#' @param loop_df Loop annotation after refinement.
#' @param bed_info Target annotation data frame (optional).
#' @param whitelist Character vector of active gene symbols.
#' @param project_name Character. Project prefix for plot titles.
#' @param karyo_bin_size Integer. Bin size for karyotype heatmaps.
#' @param species Character. Genome assembly.
#' @return A named list of ggplot / htmlwidget / grob objects.
#' @keywords internal
#' @noRd
build_refinement_plots <- function(
  original_loop_df, loop_df, bed_info,
  whitelist, project_name, karyo_bin_size, species,
  color_palette = "Paired"
) {
  red_palette <- c("#FFFFFF", "#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C", "#FC4E2A", "#E31A1C", "#BD0026", "#800026", "#000000")
  purple_palette <- c("#FFFFFF", "#F3E5F5", "#E1BEE7", "#BA68C8", "#9C27B0", "#7B1FA2", "#4A148C", "#000000")

  # Assign colours by descending loop-type frequency (not alphabetical)
  type_counts <- loop_df %>%
    dplyr::group_by(loop_type) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(n))
  all_types <- type_counts$loop_type
  custom_colors <- get_colors(length(all_types), color_palette)
  names(custom_colors) <- all_types

  plot_list <- list()

  plot_list$Comparison_Dumbbell <- .build_dumbbell_plot(
    original_loop_df, loop_df, project_name
  )

  plot_list$Target_Loop_Donut <- .build_refinement_donut(
    bed_info, loop_df, custom_colors, project_name
  )

  plot_list$Target_Sankey <- .build_sankey_plot(
    bed_info, whitelist, project_name
  )

  karyo_plots <- .build_refinement_karyotypes(
    loop_df, bed_info, species, karyo_bin_size,
    red_palette, purple_palette
  )
  for (nm in names(karyo_plots)) plot_list[[nm]] <- karyo_plots[[nm]]

  plot_list$Rose <- .build_rose_plot(loop_df, custom_colors, project_name)

  return(plot_list)
}

#' @title Standardize and Clean ChIPseeker Annotations
#'
#' @description
#' An internal helper function that parses the verbose annotation strings generated by
#' \code{ChIPseeker::annotatePeak}. It extracts the broad genomic feature category
#' while preserving the exact spatial details.
#'
#' @details
#' \code{ChIPseeker} often outputs annotations with highly specific distance or transcript
#' information, such as \code{"Promoter (<=1kb)"} or \code{"Intron (uc001.1/exon 1)"}.
#' This function creates a clean, categorical \code{annotation} column (e.g., \code{"Promoter"}, \code{"Intron"})
#' which is strictly required for robust downstream regular expression matching and Pie/Donut chart visualizations,
#' while safely moving the original verbose string to a new \code{detail_anno} column.
#'
#' @param df A data frame representation of a \code{csAnno} object (generated by \code{as.data.frame(annotatePeak(...))}).
#'
#' @return A modified data frame where:
#' \itemize{
#'   \item \code{annotation} contains the broad feature class.
#'   \item \code{detail_anno} contains the original verbose string.
#' }
#'
#' @keywords internal
#' @noRd
format_annotation_columns <- function(df) {
  if ("annotation" %in% colnames(df)) {
    df <- df %>%
      dplyr::rename(detail_anno = annotation) %>%
      dplyr::mutate(annotation = gsub(" \\(.*", "", detail_anno)) %>%
      dplyr::relocate(annotation, .before = detail_anno)
  }
  return(df)
}
