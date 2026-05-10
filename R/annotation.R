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
#' The method prioritizes physical 3D chromatin contacts. If a genomic feature overlaps an anchor looping to a distal gene, the distal gene is assigned as the target. In the absence of spatial loop evidence, the function implements a linear proximity-based fallback to assign the nearest active local gene, ensuring continuous annotation coverage.
#'
#' \strong{Hierarchical Conflict Resolution}
#' To address complex loci where a single anchor overlaps multiple promoters (e.g., dense gene clusters or bidirectional promoters), the function executes a 3-step resolution:
#' \enumerate{
#'   \item \emph{Expression Filter:} Excludes transcriptionally silent genes using a user-provided expression matrix.
#'   \item \emph{Biotype Prioritization:} Ranks remaining candidates by functional class: \code{Protein Coding > Antisense > lncRNA > Pseudogene}.
#'   \item \emph{Expression Tiebreaker:} Resolves remaining ambiguities by designating the gene with the highest transcriptional abundance as the primary target.
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
#' @param project_name Character. Prefix for output files and plot titles. Default: \code{"HiChIP"}.
#' @param color_palette Character. RColorBrewer palette name. Default: \code{"Set2"}.
#' @param karyo_bin_size Integer. Bin width in bp for karyotype heatmaps. Default: \code{1e5}.
#' @param neighbor_hop Integer. k-hop ego-network expansion order via \code{igraph::ego()}. \code{0} restricts to direct contacts. Default: \code{0}.
#' @param hub_percentile Numeric (0–1). Node-degree quantile for hub detection. Default: \code{0.95}.
#'
#' @return An invisible named list:
#' \itemize{
#'   \item \code{target_annotation} — Peak-to-gene assignments with \code{Assigned_Target_Genes_Filled} (3D-prioritised, falling back to nearest gene).
#'   \item \code{loop_annotation} — Annotated 3D interactome with \code{Putative_Target_Genes}.
#'   \item \code{anchor_annotation} — Anchor-level genomic classifications.
#'   \item \code{promoter_centric_stats} — Gene-level connectivity statistics.
#'   \item \code{distal_element_stats} — Distal-element connectivity statistics.
#'   \item \code{plot_list} — Named list of ggplot objects (donut, karyotype, rose, flower).
#' }
#' Also writes a multi-sheet Excel workbook to \code{out_dir}.
#'
#' @export
#'
#' @examples
#' # Mini example files (smaller subset for fast package checks)
#' bedpe_path <- system.file("extdata", "example_loops_mini.bedpe", package = "looplook")
#' bed_path <- system.file("extdata", "example_peaks_mini.bed", package = "looplook")
#'
#' res <- annotate_peaks_and_loops(
#'   bedpe_file = bedpe_path,
#'   target_bed = bed_path,
#'   species = "hg38",
#'   tss_region = c(-2000, 2000),
#'   out_dir = tempdir(),
#'   neighbor_hop = 0,
#'   project_name = "Quick_Example"
#' )
#' head(res$target_annotation)
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
  hub_percentile = 0.95
) {
  # Ensure output directory exists
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  extract_ids <- function(id_vec) {
    paste(unique(na.omit(as.character(id_vec))), collapse = ";")
  }

  .resolve_db <- function(arg, species_map, desc) {
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
    if (any(inherits(arg, c("TxDb", "OrgDb", "AnnotationDb")))) {
      return(list(obj = arg, pkg = .pkg_from_annotation_db(arg)))
    }
    if (is.character(arg) && nzchar(arg)) {
      if (!requireNamespace(arg, quietly = TRUE)) stop(desc, " '", arg, "' not installed")
      return(list(obj = utils::getFromNamespace(arg, arg), pkg = arg))
    }
    pkg <- species_map[[species]]
    if (is.null(pkg)) stop("Species not supported: ", species)
    if (!requireNamespace(pkg, quietly = TRUE)) stop(desc, " '", pkg, "' not installed")
    list(obj = utils::getFromNamespace(pkg, pkg), pkg = pkg)
  }
  tx_species <- list(hg38 = "TxDb.Hsapiens.UCSC.hg38.knownGene",
    hg19 = "TxDb.Hsapiens.UCSC.hg19.knownGene",
    mm10 = "TxDb.Mmusculus.UCSC.mm10.knownGene",
    mm9 = "TxDb.Mmusculus.UCSC.mm9.knownGene")
  org_species <- list(hg38 = "org.Hs.eg.db", hg19 = "org.Hs.eg.db",
    mm10 = "org.Mm.eg.db", mm9 = "org.Mm.eg.db")

  tx_res <- .resolve_db(txdb, tx_species, "TxDb")
  org_res <- .resolve_db(org_db, org_species, "OrgDb")
  txdb_obj <- tx_res$obj
  org_db_pkg <- org_res$pkg
  if (is.null(org_db_pkg) || !nzchar(org_db_pkg)) {
    stop("Unable to resolve OrgDb package name from supplied object. ",
      "Pass a package name string such as 'org.Hs.eg.db', or an installed OrgDb package object.")
  }

  gene_expr_map <- NULL
  if (!is.null(expr_matrix_file) && !is.null(sample_columns)) {
    message("Step 0: Loading expression data...")
    gene_expr_map <- load_expression_matrix(expr_matrix_file, sample_columns)
    message("    >>> Expression loaded for ", length(gene_expr_map), " genes.")
  }

  get_feature_class <- function(anno_str) {
    if (is.na(anno_str)) {
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

    return("E")
  }

  message("Step 1: Reading BEDPE file...")
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

  message("Step 2: Clustering loops...")
  valid_loops <- loops %>% dplyr::filter(!is.na(a1_id) & !is.na(a2_id))
  g <- igraph::graph_from_data_frame(valid_loops[, c("a1_id", "a2_id")], directed = FALSE)
  comp <- igraph::components(g)
  anchors$cluster_id <- NA
  comm <- intersect(anchors$anchor_id, names(comp$membership))
  if (length(comm) > 0) anchors$cluster_id[match(comm, anchors$anchor_id)] <- as.character(comp$membership[comm])
  anchors <- anchors %>% dplyr::filter(!is.na(cluster_id))
  loops <- loops %>% dplyr::left_join(anchors %>% dplyr::select(anchor_id, cluster_id), by = c("a1_id" = "anchor_id"))
  gr_anchors <- GenomicRanges::makeGRangesFromDataFrame(anchors, keep.extra.columns = TRUE)
  gr_anchors$anchor_id <- anchors$anchor_id
  gr_list <- GenomicRanges::GRangesList(split(gr_anchors, gr_anchors$cluster_id))
  cluster_regions <- unlist(GenomicRanges::reduce(gr_list))
  cluster_regions$cluster_id <- names(cluster_regions)
  names(cluster_regions) <- paste0("peak_", seq_along(cluster_regions))

  message("Step 3: Biological Classification & Topology...")
  if (length(gr_anchors) == 0) {
    warning("No valid loop anchors found; returning empty annotation.")
    return(list(anchor_annotation = data.frame(), loop_annotation = data.frame(),
      promoter_centric_stats = data.frame(), distal_element_stats = data.frame(),
      target_annotation = NULL, plots = list()))
  }
  anchor_anno <- ChIPseeker::annotatePeak(gr_anchors, TxDb = txdb_obj, tssRegion = tss_region, annoDb = org_db_pkg, verbose = FALSE)
  anchor_anno_df <- format_annotation_columns(as.data.frame(anchor_anno))
  anchor_anno_df <- resolve_gene_conflicts(anchor_anno_df, txdb_obj, org_db_pkg, tss_region, gene_expr_map)
  anchor_anno_df$type_code <- vapply(anchor_anno_df$annotation, get_feature_class, FUN.VALUE = character(1))
  map_info <- anchor_anno_df %>% dplyr::select(anchor_id, type_code, SYMBOL)
  loops_annotated <- loops %>%
    dplyr::left_join(map_info %>% dplyr::rename(t1 = type_code, s1 = SYMBOL), by = c("a1_id" = "anchor_id")) %>%
    dplyr::left_join(map_info %>% dplyr::rename(t2 = type_code, s2 = SYMBOL), by = c("a2_id" = "anchor_id"))

  get_locus_genes <- function(t1, t2, s1, s2) {
    genes <- c()
    if (!is.na(t1) && t1 %in% c("P", "G")) genes <- c(genes, s1)
    if (!is.na(t2) && t2 %in% c("P", "G")) genes <- c(genes, s2)
    paste(unique(na.omit(genes)), collapse = ";")
  }
  get_type_code <- function(t1, t2) {
    if (is.na(t1) || is.na(t2)) {
      return("Unknown")
    }
    paste(sort(c(t1, t2)), collapse = "-")
  }
  loops_annotated$loop_type <- mapply(get_type_code, loops_annotated$t1, loops_annotated$t2)
  locus_genes <- mapply(get_locus_genes, loops_annotated$t1, loops_annotated$t2, loops_annotated$s1, loops_annotated$s2)
  loops_annotated$single_loop_genes <- locus_genes
  loops_annotated$reg_loop_genes <- locus_genes

  message("    Calculating Topology (Hops)...")
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
  ids_to_genes_simple <- function(ids, lookup) {
    valid <- intersect(ids, names(lookup))
    if (length(valid) == 0) {
      return(NA_character_)
    }
    paste(sort(unique(lookup[valid])), collapse = ";")
  }

  ids_to_genes_priority <- function(ids, lookup_sym, lookup_typ) {
    valid <- intersect(ids, names(lookup_sym))
    if (length(valid) == 0) {
      return(NA_character_)
    }
    promoter_ids <- valid[!is.na(lookup_typ[valid]) & lookup_typ[valid] == "P"]
    use_ids <- if (length(promoter_ids) > 0) promoter_ids else valid
    genes_present <- lookup_sym[use_ids]
    paste(sort(unique(genes_present)), collapse = ";")
  }

  input_hop <- if (is.null(neighbor_hop)) 0 else neighbor_hop
  ego_list_loop <- igraph::ego(g, order = input_hop, nodes = nodes_in_graph, mode = "all")
  names(ego_list_loop) <- nodes_in_graph
  ego_list_target <- igraph::ego(g, order = input_hop + 1, nodes = nodes_in_graph, mode = "all")
  names(ego_list_target) <- nodes_in_graph
  anchor_topo_map <- data.frame(anchor_id = nodes_in_graph, topo_genes_p = vapply(ego_list_loop, function(x) ids_to_genes_simple(names(x), lookup_p_symbol), character(1)), topo_genes_pg = vapply(ego_list_loop, function(x) ids_to_genes_simple(names(x), lookup_pg_symbol), character(1)), tgt_genes_pg = vapply(ego_list_target, function(x) ids_to_genes_simple(names(x), lookup_pg_symbol), character(1)), tgt_genes_p = vapply(ego_list_target, function(x) ids_to_genes_simple(names(x), lookup_p_symbol), character(1)), tgt_genes_prio = vapply(ego_list_target, function(x) ids_to_genes_priority(names(x), lookup_pg_symbol, lookup_pg_type), character(1)), stringsAsFactors = FALSE)
  anchor_topo_map[is.na(anchor_topo_map)] <- NA_character_

  message("Step 4: Constructing Loop Tables...")
  loops_annotated <- loops_annotated %>%
    dplyr::left_join(anchor_topo_map %>% dplyr::select(anchor_id, pg1 = topo_genes_pg, p1 = topo_genes_p), by = c("a1_id" = "anchor_id")) %>%
    dplyr::left_join(anchor_topo_map %>% dplyr::select(anchor_id, pg2 = topo_genes_pg, p2 = topo_genes_p), by = c("a2_id" = "anchor_id")) %>%
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
    gene_annot <- ChIPseeker::annotatePeak(cluster_regions, TxDb = txdb_obj, tssRegion = tss_region, annoDb = org_db_pkg, verbose = FALSE)
    cluster_info <- format_annotation_columns(as.data.frame(gene_annot))
  }
  if ("GENENAME" %in% colnames(cluster_info)) cluster_info <- cluster_info %>% dplyr::rename(Gene_description = GENENAME)
  cluster_info$cluster_id <- as.character(cluster_info$cluster_id)
  cluster_info <- cluster_info %>% dplyr::left_join(agg_cluster_locus, by = "cluster_id")

  message("    Generating Promoter Centric Stats...")
  raw_stats_df <- dplyr::bind_rows(
    loop_annotation_final %>% dplyr::filter(anchor1_type == "P" & !is.na(anchor1_gene)) %>% dplyr::select(Gene = anchor1_gene, Neighbor_Type = anchor2_type, Loop_Type = loop_type),
    loop_annotation_final %>% dplyr::filter(anchor2_type == "P" & !is.na(anchor2_gene)) %>% dplyr::select(Gene = anchor2_gene, Neighbor_Type = anchor1_type, Loop_Type = loop_type)
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

  final_cutoff <- max(quantile(raw_stats_df$Total_Loops, hub_percentile, na.rm = TRUE), 3)
  distal_cutoff <- max(quantile(raw_stats_df$n_Linked_Distal, hub_percentile, na.rm = TRUE), 2)

  promoter_centric_df <- raw_stats_df %>%
    dplyr::mutate(
      Is_High_Connectivity_Gene = dplyr::if_else(Total_Loops >= final_cutoff, "Yes", "No"),
      Is_High_Distal_Connectivity_Gene = dplyr::if_else(n_Linked_Distal >= distal_cutoff, "Yes", "No")
    ) %>%
    dplyr::arrange(dplyr::desc(n_Linked_Distal))

  message("    Generating Distal Element Stats...")
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
      dplyr::arrange(dplyr::desc(n_Linked_Promoters))
  } else {
    distal_element_df <- NULL
  }

  bed_info <- NULL
  target_connected_loops <- NULL
  if (!is.null(target_bed)) {
    message("Step 5: Integrating Target Annotations...")
    bed_target <- read_robust_general(target_bed, min_cols = 3, desc = "Target BED")
    colnames(bed_target)[c(1, 2, 3)] <- c("chr", "start", "end")
    if (nrow(bed_target) == 0) {
      warning("Target BED contains no features; skipping target annotation.")
      bed_info <- NULL
    } else {
      bed_target$start <- bed_target$start + 1  # BED is 0-based; GRanges is 1-based
      gr_bed <- GenomicRanges::makeGRangesFromDataFrame(bed_target)
      gr_bed$input_id <- paste0("Peak_", seq_len(nrow(bed_target)))
      names(gr_bed) <- gr_bed$input_id
      bed_annot <- ChIPseeker::annotatePeak(gr_bed, TxDb = txdb_obj, tssRegion = tss_region, annoDb = org_db_pkg, verbose = FALSE)
      bed_info <- format_annotation_columns(as.data.frame(bed_annot))
      if ("GENENAME" %in% colnames(bed_info)) bed_info <- bed_info %>% dplyr::rename(Gene_description = GENENAME)
      message("    Refining Target annotation...")
      bed_info <- resolve_gene_conflicts(bed_info, txdb_obj, org_db_pkg, tss_region, gene_expr_map)
      hits <- GenomicRanges::findOverlaps(gr_bed, gr_anchors)
      if (length(hits) > 0) {
        target_connected_loops <- loop_annotation_final %>% dplyr::filter(cluster_id %in% unique(gr_anchors$cluster_id[S4Vectors::subjectHits(hits)]))
        hit_df <- data.frame(qid = S4Vectors::queryHits(hits), sid = S4Vectors::subjectHits(hits))
        hit_df$anchor_id <- gr_anchors$anchor_id[hit_df$sid]
        hit_df <- hit_df %>% dplyr::left_join(anchor_topo_map, by = "anchor_id")
        anchor_loop_agg <- dplyr::bind_rows(loop_annotation_final %>% dplyr::select(anchor_id = a1_id, loop_ID), loop_annotation_final %>% dplyr::select(anchor_id = a2_id, loop_ID)) %>%
          dplyr::distinct() %>%
          dplyr::group_by(anchor_id) %>%
          dplyr::summarise(linked_loops = extract_ids(loop_ID), .groups = "drop")
        hit_df <- hit_df %>% dplyr::left_join(anchor_loop_agg, by = "anchor_id")
        summary_df <- hit_df %>%
          dplyr::group_by(qid) %>%
          dplyr::summarise(All_Loop_Connected_Genes = extract_genes(tgt_genes_pg), Regulated_promoter_genes = extract_genes(tgt_genes_p), Assigned_Target_Genes = extract_genes(tgt_genes_prio), Linked_Loop_IDs = extract_ids(linked_loops), .groups = "drop") %>%
          dplyr::mutate(join_id = paste0("Peak_", qid))
        bed_info <- dplyr::left_join(bed_info, summary_df, by = c("input_id" = "join_id")) %>% dplyr::select(-any_of(c("join_id", "qid")))
      } else {
        bed_info$All_Loop_Connected_Genes <- NA
        bed_info$Regulated_promoter_genes <- NA
        bed_info$Assigned_Target_Genes <- NA
        bed_info$Linked_Loop_IDs <- NA
      }

      fill_logic <- function(target, fallback) dplyr::case_when(!is.na(target) & target != "" ~ target, !is.na(fallback) & fallback != "" ~ fallback, TRUE ~ NA_character_)

      bed_info <- bed_info %>% dplyr::mutate(
        All_Loop_Connected_Genes_Filled = fill_logic(All_Loop_Connected_Genes, SYMBOL),
        Regulated_promoter_genes_Filled = fill_logic(Regulated_promoter_genes, SYMBOL),
        Assigned_Target_Genes_Filled = fill_logic(Assigned_Target_Genes, SYMBOL)
      )

      if ("Linked_Loop_IDs" %in% colnames(bed_info)) {
        target_col <- if ("Gene_description" %in% colnames(bed_info)) "Gene_description" else "SYMBOL"
        if (target_col %in% colnames(bed_info)) bed_info <- bed_info %>% dplyr::relocate(Linked_Loop_IDs, .after = all_of(target_col))
      }
    }
  }


  # Step 6: Visualization
  message("Step 6: Generating Visualizations (Returning plot objects)...")
  plot_df <- loop_annotation_final
  plot_df$loop_genes <- plot_df$All_Anchor_Genes

  plot_list <- build_annotation_plots(
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

  message("Step 7: Exporting to Excel...")
  loop_annotation_clean <- loop_annotation_final %>% dplyr::select(-any_of(c("a1_id", "a2_id")))
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Loop Annotation")
  openxlsx::writeData(wb, "Loop Annotation", loop_annotation_clean)
  openxlsx::addWorksheet(wb, "Anchor Annotation")
  openxlsx::writeData(wb, "Anchor Annotation", cluster_info)
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
  tryCatch(
    openxlsx::saveWorkbook(wb, file.path(out_dir, paste0(project_name, "_Basic_Results.xlsx")), overwrite = TRUE),
    error = function(e) warning("Failed to save Excel workbook: ", conditionMessage(e), call. = FALSE)
  )
  message("    Excel file saved.")

  message("Analysis Complete.")
  return(list(
    anchor_annotation = cluster_info,
    loop_annotation = loop_annotation_clean,
    promoter_centric_stats = promoter_centric_df,
    distal_element_stats = distal_element_df,
    target_annotation = bed_info,
    plots = plot_list
  ))
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
build_annotation_plots <- function(plot_df, bed_info, cluster_info,
  target_connected_loops, txdb_obj, org_db_pkg, species, project_name,
  color_palette, karyo_bin_size) {

  loop_types_sorted <- sort(unique(plot_df$loop_type))
  custom_colors <- get_colors(length(loop_types_sorted), color_palette)
  names(custom_colors) <- loop_types_sorted

  red_palette <- c("#FFFFFF", "#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C", "#FC4E2A", "#E31A1C", "#BD0026", "#800026", "#000000")
  blue_palette <- c("#FFFFFF", "#E1F5FE", "#B3E5FC", "#4FC3F7", "#039BE5", "#0277BD", "#01579B", "#000000")
  purple_palette <- c("#FFFFFF", "#F3E5F5", "#E1BEE7", "#BA68C8", "#9C27B0", "#7B1FA2", "#4A148C", "#000000")

  plot_list <- list()

  donut_data <- plot_df %>%
    dplyr::group_by(loop_type) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(
      prop = n / sum(n),
      legend_label = paste0(loop_type, " (n=", n, ", ", round(prop * 100, 1), "%)")
    ) %>%
    dplyr::arrange(dplyr::desc(n))

  donut_data$loop_type <- factor(donut_data$loop_type, levels = donut_data$loop_type)

  plot_list$Basic_Donut <- ggplot2::ggplot(donut_data, ggplot2::aes(x = 2, y = n, fill = loop_type)) +
    ggplot2::geom_bar(stat = "identity", color = "white") +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::xlim(0.5, 2.9) +
    ggplot2::geom_text(ggplot2::aes(x = 2.65, label = loop_type),
      position = ggplot2::position_stack(vjust = 0.5), size = 2.5) +
    ggplot2::scale_fill_manual(values = rev(custom_colors),
      labels = setNames(donut_data$legend_label, donut_data$loop_type)) +
    ggplot2::theme_void() +
    ggplot2::labs(title = paste0(project_name, ": Loop Type Distribution")) +
    ggplot2::theme(legend.position = "right",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

  plot_list$Basic_Circular <- draw_circular_bar_plot(plot_df, project_name, color_vec = custom_colors)

  if ("All_Anchor_Genes" %in% colnames(plot_df)) {
    genes_loop <- clean_gene_names(plot_df$All_Anchor_Genes, ";")
    if (length(genes_loop) > 0) {
      all_genes_gr <- GenomicFeatures::genes(txdb_obj)
      map <- AnnotationDbi::select(
        utils::getFromNamespace(org_db_pkg, org_db_pkg),
        keys = as.character(S4Vectors::mcols(all_genes_gr)$gene_id),
        columns = "SYMBOL", keytype = "ENTREZID")
      S4Vectors::mcols(all_genes_gr)$SYMBOL <- map$SYMBOL[match(
        S4Vectors::mcols(all_genes_gr)$gene_id, map$ENTREZID)]
      target_genes_gr <- all_genes_gr[S4Vectors::mcols(all_genes_gr)$SYMBOL %in% genes_loop]
      plot_list$Karyo_LoopGenes <- draw_karyo_heatmap_internal(
        target_genes_gr, "Loop Genes Distribution", karyo_bin_size, 0.99,
        txdb_obj, species, "Genes", custom_colors = red_palette)
    }
  }

  all_anchors <- dplyr::bind_rows(
    plot_df %>% dplyr::select(chr = chr1, start = start1, end = end1),
    plot_df %>% dplyr::select(chr = chr2, start = start2, end = end2)
  ) %>% dplyr::distinct()
  if (nrow(all_anchors) > 0) {
    plot_list$Karyo_Anchors <- draw_karyo_heatmap_internal(
      GenomicRanges::makeGRangesFromDataFrame(all_anchors),
      "Loop Anchor Load", karyo_bin_size, 0.99, txdb_obj, species, "Anchors",
      custom_colors = blue_palette)
  }

  temp_df_flower <- plot_df %>%
    dplyr::filter(!is.na(All_Anchor_Genes) & All_Anchor_Genes != "") %>%
    tidyr::separate_rows(All_Anchor_Genes, sep = ";") %>%
    dplyr::mutate(All_Anchor_Genes = trimws(All_Anchor_Genes)) %>%
    dplyr::filter(All_Anchor_Genes != "")
  gene_sets <- split(temp_df_flower$All_Anchor_Genes, temp_df_flower$loop_type)
  gene_sets <- lapply(gene_sets, unique)
  if (length(gene_sets) > 1) {
    plot_list$Basic_Flower <- draw_flower_simplified(gene_sets, project_name, custom_colors)
  }

  if (nrow(cluster_info) > 0) {
    plot_list$Anchor_Genomic_Distribution <- draw_pie_with_outside_labels(
      cluster_info, "annotation",
      paste0(project_name, ": All Anchors Genomic Distribution"), color_palette)
  }

  if (!is.null(bed_info)) {
    if ("Assigned_Target_Genes_Filled" %in% colnames(bed_info)) {
      genes_target <- clean_gene_names(bed_info$Assigned_Target_Genes_Filled, ";")
      if (length(genes_target) > 0) {
        all_genes_gr <- GenomicFeatures::genes(txdb_obj)
        map <- AnnotationDbi::select(
          utils::getFromNamespace(org_db_pkg, org_db_pkg),
          keys = as.character(S4Vectors::mcols(all_genes_gr)$gene_id),
          columns = "SYMBOL", keytype = "ENTREZID")
        S4Vectors::mcols(all_genes_gr)$SYMBOL <- map$SYMBOL[match(
          S4Vectors::mcols(all_genes_gr)$gene_id, map$ENTREZID)]
        target_genes_gr <- all_genes_gr[S4Vectors::mcols(all_genes_gr)$SYMBOL %in% genes_target]
        plot_list$Karyo_TargetGenes <- draw_karyo_heatmap_internal(
          target_genes_gr, "Target Genes (Assigned+Local)", karyo_bin_size,
          0.99, txdb_obj, species, "Genes", custom_colors = purple_palette)
      }
    }

    if (!is.null(target_connected_loops) && nrow(target_connected_loops) > 0) {
      target_rose_data <- target_connected_loops %>%
        dplyr::group_by(loop_type) %>%
        dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
        dplyr::mutate(prop = n / sum(n),
          legend_label = paste0(loop_type, " (n=", n, ", ", round(prop * 100, 1), "%)")) %>%
        dplyr::arrange(dplyr::desc(n))

      target_rose_data$loop_type <- factor(target_rose_data$loop_type,
        levels = target_rose_data$loop_type)

      plot_list$Target_Rose <- ggplot2::ggplot(
        target_rose_data, ggplot2::aes(x = loop_type, y = n, fill = loop_type)) +
        ggplot2::geom_bar(stat = "identity", width = 1, color = "white") +
        ggplot2::coord_polar(theta = "x") +
        ggplot2::scale_fill_manual(values = custom_colors,
          labels = setNames(target_rose_data$legend_label, target_rose_data$loop_type)) +
        ggplot2::theme_void() +
        ggplot2::theme(legend.position = "right",
          plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")) +
        ggplot2::labs(title = paste0(project_name, ": Target Connected Loops (Rose)"))
    }

    plot_list$Target_Genomic_Distribution <- draw_pie_with_outside_labels(
      bed_info, "annotation",
      paste0(project_name, ": All Targets Genomic Distribution"), color_palette)

    linked_bed <- bed_info %>% dplyr::filter(!is.na(Linked_Loop_IDs) & Linked_Loop_IDs != "")
    if (nrow(linked_bed) > 0) {
      plot_list$Target_Loop_Genomic_Distribution <- draw_pie_with_outside_labels(
        linked_bed, "annotation",
        paste0(project_name, ": Loop-Connected Targets Distribution"), color_palette)
    }
  }

  return(plot_list)
}

#' @title Expression-Aware refinement of loop anchors and target linkages
#'
#' @description
#' Integrates quantitative RNA-seq data (e.g., TPM/FPKM) with 3D structural data
#' to filter and reclassify regulatory elements, deriving a functionally active
#' regulatory network from physical chromatin contacts.
#'
#' @details
#' \strong{Algorithmic Framework:}
#' \itemize{
#'   \item \strong{Target Filtration:} Parses merged gene assignments (e.g., \code{"GeneA;GeneB"}), evaluates individual genes against a defined expression threshold, and retains only transcriptionally active targets.
#'   \item \strong{Biological Reclassification:} Reclassifies physically annotated promoters (\code{P}) and gene bodies (\code{G}) lacking active transcription as enhancer-like regulatory elements (\code{eP}, \code{eG}). This adjusts the regulatory syntax to reflect functional states (e.g., reannotating a silent \code{P-P} loop to an \code{eP-P} interaction).
#'   \item \strong{Structural Hub Preservation:} Inherits foundational 3D Hub classifications (e.g., \code{Is_High_Connectivity_Gene}) derived from the raw physical network. This decouples intrinsic structural network topology from tissue-specific transcriptional activation states.
#'   \item \strong{External Target Refinement:} Filters auxiliary target mapping columns (e.g., \code{Assigned_Target_Genes_Filled}) based on expression criteria, ensuring that mapped 1D genomic features are exclusively linked to active genes.
#' }
#'
#' @param annotation_res List. The raw foundational output object returned by \code{\link{annotate_peaks_and_loops}}.
#' @param expr_matrix_file Path to a normalised expression matrix (TPM/FPKM, genes × samples). Required for refinement. Default: \code{NULL}.
#' @param sample_columns Character vector or integer indices. Columns in \code{expr_matrix_file} to average. Default: \code{NULL}.
#' @param threshold Numeric. Minimum expression (e.g. TPM > 1) for a gene to be considered active. Default: \code{1}.
#' @param unit_type Character. Expression unit for plot labels (e.g., \code{"TPM"}). Default: \code{"TPM"}.
#' @param species Character. Genome assembly. One of \code{"hg38"}, \code{"hg19"}, \code{"mm10"}, \code{"mm9"}. Default: \code{"hg38"}.
#' @param out_dir Character. Output directory for the Excel results file. Default: \code{"./results/filtered"}.
#' @param project_name Character. Prefix for output files (automatically appends \code{"_Filtered"}). Default: \code{"HiChIP"}.
#' @param color_palette Character. RColorBrewer palette name. Default: \code{"Set2"}.
#' @param karyo_bin_size Integer. Bin width in bp for karyotype heatmaps. Default: \code{1e5}.
#' @param reclassify_by_expression Logical. If \code{TRUE} (default), silent promoters (P) and gene bodies (G) are reclassified as eP/eG.
#' @param hub_percentile Numeric (0–1). Node-degree quantile for hub detection. Default: \code{0.95}.
#'
#' @return An invisible named list:
#' \itemize{
#'   \item \code{loop_annotation} — Filtered 3D network with updated \code{loop_type} (e.g., eP-P).
#'   \item \code{anchor_annotation} — Cluster annotations with expressed targets.
#'   \item \code{promoter_centric_stats} — Gene-level connectivity statistics.
#'   \item \code{distal_element_stats} — Distal-element connectivity statistics.
#'   \item \code{target_annotation} — External features linked to active loop components.
#'   \item \code{plot_list} — Named list of ggplot objects (dumbbell, rose, karyotype).
#' }
#' Also writes \code{_Refined_Results.xlsx} to \code{out_dir}.
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
#' # Extract the first object found in the RData file
#' raw_annotation <- temp_env[[ls(temp_env)[1]]]
#'
#' # =========================================================================
#' # Example : Advanced filtering WITH Transcriptome-Guided Reclassification
#' # =========================================================================
#' res_reclassified <- refine_loop_anchors_by_expression(
#'   annotation_res = raw_annotation,
#'   expr_matrix_file = expr_path,
#'   sample_columns = c("con1", "con2"),
#'   threshold = 1.0,
#'   unit_type = "TPM",
#'   species = "hg38",
#'   out_dir = tempdir(),
#'   project_name = "Example_Reclassified",
#'   reclassify_by_expression = TRUE
#' )
#'
#' # View the biologically corrected loop types (e.g., transition from P-P to eP-P)
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
  color_palette = "Set2",
  karyo_bin_size = 1e5,
  reclassify_by_expression = TRUE,
  hub_percentile = 0.95
) {
  # --- 0. Setup ---
  if (!grepl("_Filtered$", project_name)) project_name <- paste0(project_name, "_Filtered")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  message(">>> [Refinement] Project Name: ", project_name)

  # --- 1. Load Data ---
  message(">>> [Step 1] Loading Data & Expression Matrix...")
  if (is.null(annotation_res$loop_annotation)) stop("'loop_annotation' missing.")

  original_loop_df <- annotation_res$loop_annotation
  loop_df <- annotation_res$loop_annotation
  clust_info <- annotation_res$anchor_annotation
  bed_info <- annotation_res$target_annotation

  # Check and reconstruct IDs if missing
  if (!"a1_id" %in% colnames(loop_df) || !"a2_id" %in% colnames(loop_df)) {
    message("    [Info] 'a1_id'/'a2_id' columns missing. Reconstructing from coordinates...")
    loop_df <- loop_df %>%
      dplyr::mutate(
        a1_id = paste(chr1, start1, end1, sep = "_"),
        a2_id = paste(chr2, start2, end2, sep = "_")
      )
  }

  upstream_promoter_stats <- annotation_res$promoter_centric_stats
  upstream_distal_stats <- annotation_res$distal_element_stats

  if (is.null(expr_matrix_file) || is.null(sample_columns)) {
    stop("Expression matrix file and sample columns are required for refinement.")
  }
  vals <- load_expression_matrix(expr_matrix_file, sample_columns)
  whitelist <- names(vals)[vals > threshold & !is.na(vals) & names(vals) != ""]
  message(sprintf("    >>> Active Genes (> %s %s): %d", threshold, unit_type, length(whitelist)))

  # Warn if whitelist has suspiciously low overlap with annotation genes
  anno_genes <- unique(c(trimws(unlist(strsplit(na.omit(loop_df$anchor1_gene), ";"))),
    trimws(unlist(strsplit(na.omit(loop_df$anchor2_gene), ";")))))
  anno_genes <- anno_genes[nzchar(anno_genes)]
  if (length(anno_genes) > 0) {
    overlap_rate <- length(intersect(whitelist, anno_genes)) / length(anno_genes)
    if (overlap_rate < 0.1)
      warning(sprintf("Only %.1f%% of annotation gene symbols match the expression matrix row names. ",
        overlap_rate * 100),
        "Check that expression matrix row names use the same gene identifier convention (e.g., SYMBOL).")
  }

  # --- 2. Update Anchors & Loops ---
  message(">>> [Step 2] Updating Anchors & Loops...")

  a1_res <- mapply(clean_anchor, loop_df$anchor1_gene, loop_df$anchor1_type,
    MoreArgs = list(allow = whitelist, down = reclassify_by_expression),
    SIMPLIFY = FALSE)
  a2_res <- mapply(clean_anchor, loop_df$anchor2_gene, loop_df$anchor2_type,
    MoreArgs = list(allow = whitelist, down = reclassify_by_expression),
    SIMPLIFY = FALSE)
  loop_df$anchor1_type <- vapply(a1_res, function(x) x$type, character(1))
  loop_df$anchor1_gene <- vapply(a1_res, function(x) x$gene, character(1))
  loop_df$anchor2_type <- vapply(a2_res, function(x) x$type, character(1))
  loop_df$anchor2_gene <- vapply(a2_res, function(x) x$gene, character(1))
  loop_df$loop_type <- mapply(function(t1, t2) paste(sort(c(t1, t2)), collapse = "-"), loop_df$anchor1_type, loop_df$anchor2_type)

  is_enh_like <- function(t) t %in% c("E", "eP", "eG")
  is_promoter <- function(t) t == "P"

  loop_df <- loop_df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      Putative_Target_Genes = dplyr::case_when(
        (is_promoter(anchor1_type) & is_enh_like(anchor2_type)) ~ extract_genes(anchor1_gene),
        (is_enh_like(anchor1_type) & is_promoter(anchor2_type)) ~ extract_genes(anchor2_gene),
        (is_promoter(anchor1_type) & is_promoter(anchor2_type)) ~ extract_genes(c(anchor1_gene, anchor2_gene)),
        TRUE ~ extract_genes(c(anchor1_gene, anchor2_gene))
      )
    ) %>%
    dplyr::ungroup()

  # --- 3. Stats Update ---
  message(">>> [Step 3] Updating Stats...")

  agg_cluster <- loop_df %>%
    dplyr::filter(!is.na(cluster_id)) %>%
    dplyr::group_by(cluster_id) %>%
    dplyr::summarise(Cluster_All_Genes = extract_genes(Putative_Target_Genes), .groups = "drop")
  loop_df <- loop_df %>%
    dplyr::select(-any_of("Cluster_All_Genes")) %>%
    dplyr::left_join(agg_cluster, by = "cluster_id")
  if (!is.null(clust_info)) {
    clust_info <- clust_info %>%
      dplyr::select(-any_of("Cluster_All_Genes")) %>%
      dplyr::left_join(agg_cluster, by = "cluster_id")
  }

  stats_res <- compute_refined_stats(
    loop_df = loop_df,
    upstream_promoter_stats = upstream_promoter_stats,
    upstream_distal_stats = upstream_distal_stats,
    vals = vals,
    threshold = threshold,
    hub_percentile = hub_percentile
  )
  promoter_centric_df <- stats_res$promoter_centric
  distal_element_df <- stats_res$distal_element

  message(">>> [Step 4] Refining Target Annotations...")
  if (!is.null(bed_info)) {
    cols_to_clean <- grep("Strict|Physical|Loop_Genes|promoter|Filled|Target_Genes|Assigned", colnames(bed_info), value = TRUE)
    raw_tgt_col <- "Assigned_Target_Genes_Filled"
    if (!raw_tgt_col %in% colnames(bed_info))
      raw_tgt_col <- grep("Filled", cols_to_clean, value = TRUE)[1]
    if (!is.na(raw_tgt_col) && raw_tgt_col %in% colnames(bed_info)) bed_info$SANKEY_RAW_GENES <- bed_info[[raw_tgt_col]]
    for (col in cols_to_clean) {
      if (col %in% colnames(bed_info)) {
        bed_info[[col]] <- vapply(as.character(bed_info[[col]]), function(x) {
          if (is.na(x) || x == "") return(NA_character_)
          gs <- unlist(strsplit(x, ";"))
          gs_active <- gs[trimws(gs) %in% whitelist]
          if (length(gs_active) == 0) return(NA_character_)
          return(paste(unique(sort(trimws(gs_active))), collapse = ";"))
        }, FUN.VALUE = character(1))
      }
    }
  }

  # --- 5. Visualization ---
  message(">>> [Step 5] Generating Visualizations (Returning plot objects)...")
  plot_list <- build_refinement_plots(
    original_loop_df = original_loop_df,
    loop_df = loop_df,
    bed_info = bed_info,
    whitelist = whitelist,
    project_name = project_name,
    karyo_bin_size = karyo_bin_size,
    species = species
  )

  # --- 6. Export ---
  message(">>> [Step 6] Exporting Refined Results...")
  wb <- openxlsx::createWorkbook()
  loop_export <- loop_df %>% dplyr::select(-any_of(c("loop_genes", "single_loop_genes", "proximate_loop_gene")))
  openxlsx::addWorksheet(wb, "Filtered Loop Annotation")
  openxlsx::writeData(wb, "Filtered Loop Annotation", loop_export)
  openxlsx::addWorksheet(wb, "Filtered Anchor Annotation")
  openxlsx::writeData(wb, "Filtered Anchor Annotation", clust_info)

  if (!is.null(promoter_centric_df)) {
    openxlsx::addWorksheet(wb, "Filtered Promoter Stats")
    openxlsx::writeData(wb, "Filtered Promoter Stats", promoter_centric_df)
  }
  if (!is.null(distal_element_df)) {
    openxlsx::addWorksheet(wb, "Filtered Distal Stats")
    openxlsx::writeData(wb, "Filtered Distal Stats", distal_element_df)
  }

  if (!is.null(bed_info)) {
    bed_export <- bed_info %>% dplyr::select(-any_of("SANKEY_RAW_GENES"))
    openxlsx::addWorksheet(wb, "Filtered Target Annotation")
    openxlsx::writeData(wb, "Filtered Target Annotation", bed_export)
  }

  tryCatch(
    openxlsx::saveWorkbook(wb, file.path(out_dir, paste0(project_name, "_Refined_Results.xlsx")), overwrite = TRUE),
    error = function(e) warning("Failed to save refined Excel workbook: ", conditionMessage(e), call. = FALSE)
  )
  message("    Excel saved.")

  message("Refinement Complete.")
  return(list(
    loop_annotation = loop_df,
    anchor_annotation = clust_info,
    promoter_centric_stats = promoter_centric_df,
    distal_element_stats = distal_element_df,
    target_annotation = bed_info,
    plots = plot_list
  ))
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
build_refinement_plots <- function(original_loop_df, loop_df, bed_info,
  whitelist, project_name, karyo_bin_size, species) {

  red_palette <- c("#FFFFFF", "#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C", "#FC4E2A", "#E31A1C", "#BD0026", "#800026", "#000000")
  purple_palette <- c("#FFFFFF", "#F3E5F5", "#E1BEE7", "#BA68C8", "#9C27B0", "#7B1FA2", "#4A148C", "#000000")

  # Assign Paired colours by descending loop-type frequency (not alphabetical)
  type_counts <- loop_df %>%
    dplyr::group_by(loop_type) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(n))
  all_types <- type_counts$loop_type
  custom_colors <- grDevices::colorRampPalette(
    RColorBrewer::brewer.pal(12, "Paired"))(length(all_types))
  names(custom_colors) <- all_types

  plot_list <- list()

  # Dumbbell comparison
  df_orig <- original_loop_df %>% dplyr::group_by(loop_type) %>%
    dplyr::summarise(Original = dplyr::n(), .groups = "drop")
  df_filt <- loop_df %>% dplyr::group_by(loop_type) %>%
    dplyr::summarise(Filtered = dplyr::n(), .groups = "drop")
  df_dumbbell <- dplyr::full_join(df_orig, df_filt, by = "loop_type") %>%
    dplyr::mutate(Original = ifelse(is.na(Original), 0, Original),
      Filtered = ifelse(is.na(Filtered), 0, Filtered),
      is_e_type = grepl("e", loop_type)) %>%
    dplyr::arrange(is_e_type, dplyr::desc(Original))
  df_dumbbell$loop_type <- factor(df_dumbbell$loop_type,
    levels = rev(df_dumbbell$loop_type))
  df_long <- df_dumbbell %>%
    tidyr::pivot_longer(cols = c("Original", "Filtered"),
      names_to = "Source", values_to = "Count")

  plot_list$Comparison_Dumbbell <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = df_dumbbell,
      ggplot2::aes(y = loop_type, yend = loop_type, x = Original, xend = Filtered),
      color = "#b2b2b2", linewidth = 0.8) +
    ggplot2::geom_point(data = df_long,
      ggplot2::aes(x = Count, y = loop_type, color = Source), size = 3) +
    ggplot2::scale_color_manual(values = c("Original" = "#999999", "Filtered" = "#E69F00")) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = paste0(project_name, ": Filtration Effect (Dumbbell)"),
      x = "Number of Loops", y = "Loop Type") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "top")

  # Donut (target-connected loops)
  if (!is.null(bed_info)) {
    gr_bed <- GenomicRanges::makeGRangesFromDataFrame(bed_info,
      keep.extra.columns = TRUE)
    active_anc <- dplyr::bind_rows(
      loop_df %>% dplyr::select(chr = chr1, start = start1, end = end1, cluster_id),
      loop_df %>% dplyr::select(chr = chr2, start = start2, end = end2, cluster_id)
    ) %>% dplyr::distinct()
    if (nrow(active_anc) > 0) {
      gr_anc <- GenomicRanges::makeGRangesFromDataFrame(active_anc,
        keep.extra.columns = TRUE)
      hits <- GenomicRanges::findOverlaps(gr_bed, gr_anc)
      if (length(hits) > 0) {
        hit_ids <- unique(gr_anc$cluster_id[S4Vectors::subjectHits(hits)])
        tgt_loops <- loop_df %>% dplyr::filter(cluster_id %in% hit_ids)
        if (nrow(tgt_loops) > 0) {
          donut_data <- tgt_loops %>%
            dplyr::group_by(loop_type) %>%
            dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
            dplyr::mutate(fraction = count / sum(count),
              legend_label = paste0(loop_type, " (n=", count, ", ",
                round(fraction * 100, 1), "%)"),
              plot_label = loop_type, is_lower_e = grepl("^e", loop_type)) %>%
            dplyr::arrange(is_lower_e, dplyr::desc(count)) %>%
            dplyr::mutate(loop_type = factor(loop_type, levels = loop_type))
          plot_list$Target_Loop_Donut <- ggplot2::ggplot(donut_data,
            ggplot2::aes(x = 2, y = count, fill = loop_type)) +
            ggplot2::geom_bar(stat = "identity", width = 1, color = "white") +
            ggplot2::coord_polar(theta = "y") +
            ggplot2::xlim(0.5, 2.9) +
            ggplot2::geom_text(ggplot2::aes(x = 2.8, label = plot_label),
              position = ggplot2::position_stack(vjust = 0.5), size = 3) +
            ggplot2::scale_fill_manual(values = custom_colors,
              labels = setNames(donut_data$legend_label,
                as.character(donut_data$loop_type)), name = "Loop Type") +
            ggplot2::theme_void() +
            ggplot2::labs(title = paste0(project_name, ": Loops Connected to Targets")) +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
              face = "bold", size = 14), legend.position = "right",
            legend.text = ggplot2::element_text(size = 10))
        }
      }
    }
  }

  # Sankey
  if (!is.null(bed_info) && "SANKEY_RAW_GENES" %in% colnames(bed_info) &&
    requireNamespace("networkD3", quietly = TRUE)) {
    raw_bed <- bed_info
    total_targets <- nrow(raw_bed)
    get_label_mapping <- function(vec) {
      tbl <- table(vec)
      tbl <- tbl[tbl > 0]
      if (length(tbl) == 0) return(character(0))
      labels <- paste0(names(tbl), " (n=", tbl, ", ",
        round(as.numeric(tbl) / total_targets * 100, 1), "%)")
      names(labels) <- names(tbl)
      return(labels)
    }
    check_status_strict <- function(g_str) {
      if (is.na(g_str) || g_str == "") return("No Gene Assigned")
      gs <- trimws(unlist(strsplit(as.character(g_str), ";")))
      if (any(gs %in% whitelist)) return("Active")
      return("Inactive")
    }
    sankey_data_raw <- raw_bed %>%
      dplyr::mutate(
        L1_Raw = dplyr::case_when(
          grepl("Distal", annotation, ignore.case = TRUE) ~ "Distal Intergenic",
          grepl("Exon", annotation, ignore.case = TRUE) ~ "Exon",
          grepl("Intron", annotation, ignore.case = TRUE) ~ "Intron",
          grepl("Promoter", annotation, ignore.case = TRUE) ~ "Promoter",
          TRUE ~ "Others"),
        L2_Raw = ifelse(!is.na(Linked_Loop_IDs) & Linked_Loop_IDs != "",
          "Connected", "Unconnected"),
        L3_Raw = vapply(SANKEY_RAW_GENES, check_status_strict,
          FUN.VALUE = character(1))
      ) %>%
      dplyr::filter(L3_Raw != "No Gene Assigned")
    if (nrow(sankey_data_raw) > 0) {
      l1_map <- get_label_mapping(sankey_data_raw$L1_Raw)
      l2_map <- get_label_mapping(sankey_data_raw$L2_Raw)
      l3_map <- get_label_mapping(sankey_data_raw$L3_Raw)
      sankey_data_ready <- sankey_data_raw %>%
        dplyr::mutate(Genomic_Distribution = l1_map[L1_Raw],
          Loop_Connection = l2_map[L2_Raw],
          Expression_Status = l3_map[L3_Raw]) %>%
        dplyr::filter(!is.na(Genomic_Distribution) &
          !is.na(Loop_Connection) & !is.na(Expression_Status))
      if (nrow(sankey_data_ready) > 0) {
        links <- dplyr::bind_rows(
          sankey_data_ready %>% dplyr::group_by(source = Genomic_Distribution,
            target = Loop_Connection) %>%
            dplyr::summarise(value = dplyr::n(), .groups = "drop"),
          sankey_data_ready %>% dplyr::group_by(source = Loop_Connection,
            target = Expression_Status) %>%
            dplyr::summarise(value = dplyr::n(), .groups = "drop"))
        nodes_vec <- unique(c(links$source, links$target))
        nodes <- data.frame(name = nodes_vec, stringsAsFactors = FALSE)
        links$IDsource <- match(links$source, nodes$name) - 1
        links$IDtarget <- match(links$target, nodes$name) - 1
        sankey_colors <- grDevices::colorRampPalette(
          RColorBrewer::brewer.pal(12, "Paired"))(length(nodes_vec))
        colourScale <- sprintf('d3.scaleOrdinal().range(["%s"])',
          paste(sankey_colors, collapse = '","'))
        sn <- networkD3::sankeyNetwork(Links = links, Nodes = nodes,
          Source = "IDsource", Target = "IDtarget", Value = "value",
          NodeID = "name", units = "Targets", fontSize = 14,
          fontFamily = "Arial", nodeWidth = 15, nodePadding = 15,
          iterations = 0, height = 600, width = 900,
          colourScale = colourScale, sinksRight = FALSE)
        plot_list$Target_Sankey <- htmlwidgets::onRender(sn,
          'function(el, x) { var svg = d3.select(el).select("svg"); function createValidID(name) { if (!name) return "unknown"; return name.replace(/[^a-zA-Z0-9-]/g, "_"); } svg.selectAll(".link").each(function(d) { var gradientID = "gradient-" + createValidID(d.source.name) + "-" + createValidID(d.target.name); if (svg.select("#" + gradientID).empty()) { var gradient = svg.append("defs").append("linearGradient").attr("id", gradientID).attr("gradientUnits", "userSpaceOnUse").attr("x1", d.source.x + d.source.dx / 2).attr("y1", d.source.y + d.source.dy / 2).attr("x2", d.target.x + d.target.dx / 2).attr("y2", d.target.y + d.target.dy / 2); var sourceColor = d3.select(el).selectAll(".node").filter(function(node) { return node.name === d.source.name; }).select("rect").style("fill"); var targetColor = d3.select(el).selectAll(".node").filter(function(node) { return node.name === d.target.name; }).select("rect").style("fill"); gradient.append("stop").attr("offset", "0%").attr("stop-color", sourceColor); gradient.append("stop").attr("offset", "100%").attr("stop-color", targetColor); } d3.select(this).style("stroke", "url(#" + gradientID + ")"); }); svg.selectAll(".node rect").style("stroke", "black").style("stroke-width", "1px"); }')
      }
    }
  }

  # Karyotype heatmaps
  txdb_pkg <- switch(species,
    "hg38" = "TxDb.Hsapiens.UCSC.hg38.knownGene",
    "hg19" = "TxDb.Hsapiens.UCSC.hg19.knownGene",
    "mm10" = "TxDb.Mmusculus.UCSC.mm10.knownGene",
    "mm9"  = "TxDb.Mmusculus.UCSC.mm9.knownGene", NULL)
  org_db <- switch(species,
    "hg38" = "org.Hs.eg.db",
    "hg19" = "org.Hs.eg.db",
    "mm10" = "org.Mm.eg.db",
    "mm9"  = "org.Mm.eg.db", NULL)
  if (!is.null(txdb_pkg) && !is.null(org_db) &&
      requireNamespace(txdb_pkg, quietly = TRUE) &&
      requireNamespace(org_db, quietly = TRUE)) {
    txdb_obj <- utils::getFromNamespace(txdb_pkg, txdb_pkg)
    all_genes <- GenomicFeatures::genes(txdb_obj)
    org_db_obj <- utils::getFromNamespace(org_db, org_db)
    map <- AnnotationDbi::select(org_db_obj,
      keys = as.character(all_genes$gene_id),
      columns = "SYMBOL", keytype = "ENTREZID")
    all_genes$SYMBOL <- map$SYMBOL[match(all_genes$gene_id, map$ENTREZID)]

    g_active <- clean_gene_names(loop_df$Putative_Target_Genes, ";")
    if (length(g_active) > 0) {
      plot_list$Refined_Karyo_Active <- draw_karyo_heatmap_internal(
        all_genes[all_genes$SYMBOL %in% g_active],
        "Refined Active Genes", karyo_bin_size, 0.99, txdb_obj, species,
        "Genes", custom_colors = red_palette)
    }

    clean_tgt_col <- if (!is.null(bed_info) &&
      "Assigned_Target_Genes_Filled" %in% colnames(bed_info))
      "Assigned_Target_Genes_Filled" else NULL
    if (!is.null(clean_tgt_col)) {
      g_tgt <- clean_gene_names(bed_info[[clean_tgt_col]], ";")
      if (length(g_tgt) > 0) {
        plot_list$Refined_Karyo_TargetGenes <- draw_karyo_heatmap_internal(
          all_genes[all_genes$SYMBOL %in% g_tgt],
          "Refined Target Genes", karyo_bin_size, 0.99, txdb_obj, species,
          "Genes", custom_colors = purple_palette)
      }
    }
  }

  # Rose
  rose_data <- loop_df %>%
    dplyr::group_by(loop_type) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(fraction = count / sum(count),
      legend_label = paste0(loop_type, " (n=", count, ", ",
        round(fraction * 100, 1), "%)"),
      is_lower_e = grepl("^e", loop_type)) %>%
    dplyr::arrange(dplyr::desc(count))
  plot_order <- rose_data$loop_type
  rose_data$loop_type <- factor(rose_data$loop_type, levels = plot_order)
  legend_order <- rose_data %>%
    dplyr::arrange(is_lower_e, dplyr::desc(count)) %>%
    dplyr::pull(loop_type)

  plot_list$Rose <- ggplot2::ggplot(rose_data,
    ggplot2::aes(x = loop_type, y = count, fill = loop_type)) +
    ggplot2::geom_bar(stat = "identity", width = 1, color = "white") +
    ggplot2::coord_polar(theta = "x") +
    ggplot2::scale_fill_manual(values = custom_colors,
      labels = setNames(rose_data$legend_label, as.character(rose_data$loop_type)),
      breaks = legend_order, name = "Loop Type") +
    ggplot2::theme_void() +
    ggplot2::labs(title = paste0(project_name, ": Loop Proportion (By Count)")) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
      face = "bold", size = 14), legend.position = "right",
    legend.text = ggplot2::element_text(size = 10))

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
format_annotation_columns <- function(df) {
  if ("annotation" %in% colnames(df)) {
    df <- df %>%
      dplyr::rename(detail_anno = annotation) %>%
      dplyr::mutate(annotation = gsub(" \\(.*", "", detail_anno)) %>%
      dplyr::relocate(annotation, .before = detail_anno)
  }
  return(df)
}
