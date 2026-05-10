#' Internal Package Imports
#'
#' @name looplook_imports
#' @noRd
#' @importFrom grDevices col2rgb rgb
#' @importFrom stats fisher.test median na.omit p.adjust quantile reorder runif setNames t.test var wilcox.test
#' @importFrom utils head read.table write.csv
NULL

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    ".", ".data", "All_Anchor_Genes", "All_Loop_Connected_Genes", "Assigned_Target_Genes",
    "CleanLoopType", "Conn_Group", "Count", "Degree", "Description",
    "Description_unique", "Distal_Anchor_ID", "Dominant_Interaction",
    "Dominant_Interaction_Filtered", "Expression", "Expression_Status", "FDR",
    "Feature", "Filtered", "Final_Label", "Fraction", "GENENAME", "Gene",
    "Genomic_Distribution", "Group", "High_Connectivity_Gene",
    "Is_Active_Gene", "Is_High_Connectivity_Distal_Element",
    "Is_High_Connectivity_Gene", "Is_High_Distal_Connectivity_Gene", "L1_Raw",
    "L2_Raw", "L3_Raw", "LFC", "Label", "LabelText", "Label_Text",
    "Linked_Loop_IDs", "Log10Degree", "LogFDR", "LogP", "Loop_Type",
    "Mean_Expression_Temp", "MotifLabel", "ONTOLOGY", "OddsRatio", "Original",
    "Percentage", "PlotFamily", "Putative_Target_Genes", "Rank",
    "Regulated_promoter_genes", "SANKEY_RAW_GENES", "SYMBOL", "SampleID",
    "Simplified", "Source", "Stage", "Status", "Target_Genes",
    "Target_Genes_Filtered", "Total_Loops", "Total_Loops_Filtered",
    "Unique_Gene_Count", "a1_id", "a2_id", "all_cluster_loop_genes", "all_of",
    "anchor1_gene", "anchor1_source", "anchor1_type", "anchor2_gene",
    "anchor2_source", "anchor2_type", "anchor_id", "annotation", "chr",
    "cluster_id", "col2rgb", "combined_score", "count", "deg", "detail_anno",
    "elementNROWS", "everything", "expansion", "final_color", "final_fill",
    "final_label", "final_symbols", "fisher.test", "fraction",
    "functional_anchor1_type", "functional_anchor2_type", "geneList",
    "gene_id", "gene_level", "geom_hline", "group", "has_active", "head",
    "hjust", "install.packages", "is_e_type", "is_lower_e", "label",
    "labelPosition", "label_text", "label_x", "len", "lfc", "linked_loops",
    "log_expr", "logP", "loop_ID", "loop_genes", "loop_genes_Total", "loop_i",
    "loop_type", "median", "mid1", "mid2", "n", "n_Linked_Distal",
    "n_Linked_Distal_Filtered", "n_Linked_Promoters",
    "n_Linked_Promoters_Filtered", "na.omit", "name", "p.adjust", "plot_label",
    "prop", "proximate_loop_gene", "pvalue", "qid", "quantile", "query_idx",
    "read.table", "reg_loop_genes", "reorder", "rgb", "runif", "runningScore",
    "scale_color_identity", "scale_fill_identity", "setNames",
    "single_loop_genes", "strand", "t.test", "t1", "t2", "tgt_genes_p",
    "tgt_genes_pg", "tgt_genes_prio", "topo_genes_p", "topo_genes_pg", "tpm",
    "tx_id", "type", "type_code", "type_rank", "valid_genes", "valid_tpms",
    "var", "width", "wilcox.test", "write.csv", "y_mid", "ymax", "ymin", "ypos", ":=",
    "Loop_Connection", "Neighbor_Gene", "Neighbor_Type", "s1", "s2", "x", "y"
  ))
}

#' Internal: Clean Gene Name Vector
#'
#' Removes empty strings, NA values, and duplicate entries from gene identifiers.
#' Optionally splits concatenated strings (e.g., "TP53;BRCA1") before cleaning.
#'
#' @param x Character vector of gene names, possibly containing delimiters.
#' @param split Character. If non-NULL, a regex passed to \code{\link{strsplit}}
#'   to split concatenated gene strings (e.g., \code{"[;,]"}). Set to \code{NULL}
#'   if \code{x} is already a clean character vector.
#' @return A unique, non-empty, non-NA character vector.
#' @keywords internal
clean_gene_names <- function(x, split = NULL) {
  if (is.null(x) || length(x) == 0) {
    return(character(0))
  }
  if (!is.null(split)) x <- unlist(strsplit(as.character(x), split))
  x <- unique(trimws(as.character(x)))
  x[x != "" & !is.na(x)]
}


#' Internal: Collapse Delimited Gene String
#'
#' Splits a semicolon-delimited gene string, removes duplicates and NAs,
#' and recollapses into a single string. Returns \code{NA_character_} if
#' no valid genes remain.
#'
#' @param genes_vec Character vector of delimited gene strings.
#' @return A single semicolon-delimited string, or \code{NA_character_}.
#' @keywords internal
extract_genes <- function(genes_vec) {
  res <- unique(na.omit(unlist(strsplit(as.character(genes_vec), ";"))))
  res <- res[nzchar(res)]
  if (length(res) == 0) {
    return(NA_character_)
  }
  paste(res, collapse = ";")
}


#' Internal: Resolve Gene Conflicts via Expression & Biotype
#'
#' For each genomic range in an annotation data frame, identifies all genes
#' whose promoters overlap the range, then selects the best candidate using
#' expression level and biotype priority (protein-coding > antisense > lncRNA >
#' pseudo). When multiple genes have similar expression, all are retained
#' (collapsed with ";").
#'
#' @param current_anno_df Data frame with columns suitable for
#'   \code{\link[GenomicRanges]{makeGRangesFromDataFrame}}.
#' @param txdb_obj A \code{TxDb} object for gene coordinate lookup.
#' @param org_db_pkg Character. Organism database package name.
#' @param tss_region Numeric vector of length 2. TSS region for promoter
#'   definition, e.g. \code{c(-2000, 2000)}.
#' @param gene_expr_map Named numeric vector of per-gene expression values,
#'   or \code{NULL} if unavailable.
#' @return The input data frame with \code{SYMBOL} and \code{annotation}
#'   columns resolved.
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps
#' @importFrom GenomicFeatures genes promoters
#' @importFrom S4Vectors queryHits subjectHits
#' @keywords internal
resolve_gene_conflicts <- function(
  current_anno_df, txdb_obj, org_db_pkg,
  tss_region, gene_expr_map
) {
  if (nrow(current_anno_df) == 0) {
    return(current_anno_df)
  }

  gr_input <- GenomicRanges::makeGRangesFromDataFrame(current_anno_df,
    keep.extra.columns = TRUE
  )
  all_genes <- GenomicFeatures::genes(txdb_obj)
  hits <- GenomicRanges::findOverlaps(
    gr_input,
    GenomicFeatures::promoters(all_genes,
      upstream = abs(tss_region[1]),
      downstream = abs(tss_region[2])
    )
  )

  if (length(hits) > 0) {
    candidates <- data.frame(
      query_idx = S4Vectors::queryHits(hits),
      gene_id = names(all_genes)[S4Vectors::subjectHits(hits)],
      stringsAsFactors = FALSE
    )

    org_db_obj <- utils::getFromNamespace(org_db_pkg, org_db_pkg)
    valid_keys <- AnnotationDbi::keytypes(org_db_obj)
    primary_key <- if ("ENTREZID" %in% valid_keys) "ENTREZID" else valid_keys[1]

    cols_to_get <- "SYMBOL"
    valid_cols <- AnnotationDbi::columns(org_db_obj)
    has_genetype <- "GENETYPE" %in% valid_cols
    if (has_genetype) cols_to_get <- c(cols_to_get, "GENETYPE")

    gene_map <- AnnotationDbi::select(
      org_db_obj,
      keys = unique(candidates$gene_id),
      columns = cols_to_get,
      keytype = primary_key
    )

    gene_map$tpm <- if (!is.null(gene_expr_map)) {
      ifelse(is.na(gene_expr_map[gene_map$SYMBOL]), 0,
        gene_expr_map[gene_map$SYMBOL]
      )
    } else {
      0
    }

    if (has_genetype) {
      gene_map <- gene_map %>%
        dplyr::mutate(
          type_rank = dplyr::case_when(
            grepl("protein", GENETYPE, ignore.case = TRUE) ~ 1,
            grepl("antisense", GENETYPE, ignore.case = TRUE) ~ 2,
            grepl("lncRNA|ncrna", GENETYPE, ignore.case = TRUE) ~ 3,
            grepl("pseudo", GENETYPE, ignore.case = TRUE) ~ 4,
            TRUE ~ 5
          )
        )
    } else {
      gene_map$type_rank <- 1 # All genes have equal rank if GENETYPE is missing
    }

    join_by_args <- setNames(primary_key, "gene_id")
    resolved_candidates <- candidates %>%
      dplyr::left_join(gene_map, by = join_by_args) %>%
      dplyr::group_by(query_idx) %>%
      dplyr::mutate(has_active = any(tpm > 0)) %>%
      dplyr::filter(!has_active | tpm > 0) %>%
      dplyr::filter(type_rank == min(type_rank, na.rm = TRUE)) %>%
      dplyr::summarise(
        valid_genes = list(SYMBOL[!is.na(SYMBOL) & SYMBOL != ""]),
        valid_tpms = list(tpm[!is.na(SYMBOL) & SYMBOL != ""]),
        .groups = "drop"
      ) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        final_symbols = {
          genes <- unlist(valid_genes)
          tpms <- unlist(valid_tpms)
          if (length(genes) == 0) {
            NA_character_
          } else if (length(genes) == 1) {
            genes[1]
          } else {
            max_tpm <- max(tpms, na.rm = TRUE)
            if (max_tpm <= 0) {
              paste(sort(unique(genes)), collapse = ";")
            } else {
              active_genes <- genes[tpms >= max_tpm * 0.1]
              paste(sort(unique(active_genes)), collapse = ";")
            }
          }
        }
      ) %>%
      dplyr::ungroup() %>%
      dplyr::filter(!is.na(final_symbols))

    if (!"SYMBOL" %in% colnames(current_anno_df)) {
      current_anno_df$SYMBOL <- NA_character_
    }
    if (!"annotation" %in% colnames(current_anno_df)) {
      current_anno_df$annotation <- NA_character_
    }

    match_idx <- match(
      resolved_candidates$query_idx,
      seq_len(nrow(current_anno_df))
    )
    valid_idx <- !is.na(match_idx)

    if (any(valid_idx)) {
      safe_match <- match_idx[valid_idx]
      current_anno_df$SYMBOL[safe_match] <-
        resolved_candidates$final_symbols[valid_idx]
      current_anno_df$annotation[safe_match] <- "Promoter"
    }
  }

  return(current_anno_df)
}


#' Internal: Reclassify Anchor by Expression
#'
#' Given an anchor's gene symbol and type, filters to active genes (present in
#' \code{allow}) and optionally reclassifies silent promoters/enhancers.
#'
#' @param g Character. Semicolon-delimited gene string.
#' @param t Character. Anchor type code (P, E, G, eP, eG).
#' @param allow Character vector. Whitelist of active gene symbols.
#' @param down Logical. If \code{TRUE}, reclassify silent P→eP and G→eG.
#' @return A list with \code{type} and \code{gene}.
#' @keywords internal
clean_anchor <- function(g, t, allow, down) {
  g_char <- as.character(g)
  t_char <- as.character(t)
  if (is.na(g_char) || g_char == "") {
    return(list(type = t_char, gene = NA_character_))
  }
  gs <- unlist(strsplit(g_char, ";"))
  active_gs <- gs[trimws(gs) %in% allow]
  if (length(active_gs) > 0) {
    return(list(type = t_char, gene = paste(unique(active_gs), collapse = ";")))
  }
  if (down) {
    new_type <- dplyr::case_when(
      t_char == "P" ~ "eP", t_char == "G" ~ "eG", TRUE ~ t_char
    )
    return(list(type = new_type, gene = NA_character_))
  }
  return(list(type = t_char, gene = NA_character_))
}


#' Internal: Compute Refined Network Statistics
#'
#' Recalculates promoter-centric and distal-element connectivity statistics
#' after expression-aware filtering, merging with upstream annotation stats
#' where available.
#'
#' @param loop_df Loop annotation data frame after expression filtering.
#' @param upstream_promoter_stats Upstream promoter-centric stats (or NULL).
#' @param upstream_distal_stats Upstream distal-element stats (or NULL).
#' @param vals Named numeric vector of per-gene mean expression.
#' @param threshold Numeric. Expression threshold for active gene classification.
#' @param hub_percentile Numeric. Quantile for hub cutoff.
#' @return A list with \code{promoter_centric} and \code{distal_element}
#'   data frames.
#' @importFrom stats quantile
#' @keywords internal
compute_refined_stats <- function(
  loop_df, upstream_promoter_stats,
  upstream_distal_stats, vals, threshold, hub_percentile
) {
  get_dom <- function(x) {
    if (length(x) == 0) {
      return(NA_character_)
    }
    names(which.max(table(x)))
  }
  get_expr <- function(g) {
    e <- vals[g]
    e[is.na(e)] <- 0
    return(e)
  }

  raw_stats_df <- dplyr::bind_rows(
    loop_df %>% dplyr::filter(anchor1_type == "P" & !is.na(anchor1_gene)) %>%
      dplyr::select(
        Gene = anchor1_gene, Neighbor_Type = anchor2_type,
        Loop_Type = loop_type
      ),
    loop_df %>% dplyr::filter(anchor2_type == "P" & !is.na(anchor2_gene)) %>%
      dplyr::select(
        Gene = anchor2_gene, Neighbor_Type = anchor1_type,
        Loop_Type = loop_type
      )
  ) %>%
    tidyr::separate_rows(Gene, sep = ";") %>%
    dplyr::mutate(Gene = trimws(Gene)) %>%
    dplyr::filter(Gene != "" & !is.na(Gene)) %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(
      Total_Loops_Filtered = dplyr::n(),
      n_Linked_Promoters_Filtered = sum(Neighbor_Type == "P", na.rm = TRUE),
      n_Linked_Distal_Filtered = sum(
        Neighbor_Type %in% c("E", "eP", "eG", "G"),
        na.rm = TRUE
      ),
      Dominant_Interaction_Filtered = get_dom(Loop_Type),
      .groups = "drop"
    )

  promoter_centric_df <- NULL
  if (nrow(raw_stats_df) > 0) {
    final_cutoff <- max(stats::quantile(
      raw_stats_df$Total_Loops_Filtered, hub_percentile,
      na.rm = TRUE
    ), 3)
    distal_cutoff <- max(stats::quantile(
      raw_stats_df$n_Linked_Distal_Filtered, hub_percentile,
      na.rm = TRUE
    ), 2)

    if (!is.null(upstream_promoter_stats)) {
      promoter_centric_df <- upstream_promoter_stats %>%
        dplyr::left_join(raw_stats_df, by = "Gene") %>%
        dplyr::mutate(
          Total_Loops = dplyr::coalesce(Total_Loops_Filtered, 0),
          n_Linked_Promoters = dplyr::coalesce(n_Linked_Promoters_Filtered, 0),
          n_Linked_Distal = dplyr::coalesce(n_Linked_Distal_Filtered, 0),
          Dominant_Interaction = dplyr::coalesce(
            Dominant_Interaction_Filtered, "None"
          ),
          Mean_Expression_Temp = get_expr(Gene),
          Is_Active_Gene = dplyr::if_else(
            Mean_Expression_Temp > threshold, "Yes", "No"
          ),
          Is_High_Connectivity_Gene = dplyr::if_else(
            Total_Loops >= final_cutoff, "Yes", "No"
          ),
          Is_High_Distal_Connectivity_Gene = dplyr::if_else(
            n_Linked_Distal >= distal_cutoff, "Yes", "No"
          )
        ) %>%
        dplyr::select(
          Gene, Total_Loops, n_Linked_Promoters, n_Linked_Distal,
          Dominant_Interaction, Is_High_Connectivity_Gene,
          Is_High_Distal_Connectivity_Gene, Is_Active_Gene,
          dplyr::everything()
        ) %>%
        dplyr::select(-any_of(c(
          "Total_Loops_Filtered",
          "n_Linked_Promoters_Filtered", "n_Linked_Distal_Filtered",
          "Dominant_Interaction_Filtered", "Is_Regulatory_Hub",
          "Mean_Expression_Temp", "n_Linked_Enhancers", "n_Linked_GeneBodies"
        )))
    } else {
      promoter_centric_df <- raw_stats_df %>%
        dplyr::rename(
          Total_Loops = Total_Loops_Filtered,
          n_Linked_Promoters = n_Linked_Promoters_Filtered,
          n_Linked_Distal = n_Linked_Distal_Filtered,
          Dominant_Interaction = Dominant_Interaction_Filtered
        ) %>%
        dplyr::mutate(
          Mean_Expression_Temp = get_expr(Gene),
          Is_Active_Gene = dplyr::if_else(
            Mean_Expression_Temp > threshold, "Yes", "No"
          ),
          Is_High_Connectivity_Gene = dplyr::if_else(
            Total_Loops >= final_cutoff, "Yes", "No"
          ),
          Is_High_Distal_Connectivity_Gene = dplyr::if_else(
            n_Linked_Distal >= distal_cutoff, "Yes", "No"
          )
        ) %>%
        dplyr::select(
          Gene, Total_Loops, n_Linked_Promoters, n_Linked_Distal,
          Dominant_Interaction, Is_High_Connectivity_Gene,
          Is_High_Distal_Connectivity_Gene, Is_Active_Gene,
          dplyr::everything()
        ) %>%
        dplyr::select(-any_of("Mean_Expression_Temp"))
    }
    promoter_centric_df <- promoter_centric_df %>%
      dplyr::arrange(dplyr::desc(n_Linked_Distal))
  }

  distal_element_df <- NULL
  if ("a1_id" %in% colnames(loop_df)) {
    distal_raw_df <- dplyr::bind_rows(
      loop_df %>% dplyr::filter(anchor1_type %in% c("E", "eP", "eG", "G")) %>%
        dplyr::select(
          Distal_Anchor_ID = a1_id, Neighbor_Type = anchor2_type,
          Loop_Type = loop_type, Neighbor_Gene = anchor2_gene
        ),
      loop_df %>% dplyr::filter(anchor2_type %in% c("E", "eP", "eG", "G")) %>%
        dplyr::select(
          Distal_Anchor_ID = a2_id, Neighbor_Type = anchor1_type,
          Loop_Type = loop_type, Neighbor_Gene = anchor1_gene
        )
    ) %>%
      dplyr::group_by(Distal_Anchor_ID) %>%
      dplyr::summarise(
        Total_Loops_Filtered = dplyr::n(),
        n_Linked_Distal_Filtered = sum(
          Neighbor_Type %in% c("E", "eP", "eG", "G"),
          na.rm = TRUE
        ),
        n_Linked_Promoters_Filtered = sum(Neighbor_Type == "P", na.rm = TRUE),
        Dominant_Interaction_Filtered = get_dom(Loop_Type),
        Target_Genes_Filtered = extract_genes(
          Neighbor_Gene[Neighbor_Type == "P"]
        ),
        .groups = "drop"
      )

    anchor_map <- dplyr::bind_rows(
      loop_df %>% dplyr::select(
        anchor_id = a1_id, chr = chr1,
        start = start1, end = end1, cluster_id
      ),
      loop_df %>% dplyr::select(
        anchor_id = a2_id, chr = chr2,
        start = start2, end = end2, cluster_id
      )
    ) %>% dplyr::distinct()

    if (nrow(distal_raw_df) > 0) {
      final_cutoff_dist <- max(stats::quantile(
        distal_raw_df$Total_Loops_Filtered, hub_percentile,
        na.rm = TRUE
      ), 3)
      if (!is.null(upstream_distal_stats) &&
        "Distal_Anchor_ID" %in% colnames(upstream_distal_stats)) {
        temp_df <- upstream_distal_stats %>%
          dplyr::left_join(distal_raw_df, by = "Distal_Anchor_ID") %>%
          dplyr::mutate(
            Total_Loops = dplyr::coalesce(Total_Loops_Filtered, 0),
            n_Linked_Distal = dplyr::coalesce(n_Linked_Distal_Filtered, 0),
            n_Linked_Promoters = dplyr::coalesce(
              n_Linked_Promoters_Filtered, 0
            ),
            Dominant_Interaction = dplyr::coalesce(
              Dominant_Interaction_Filtered, "None"
            ),
            Target_Genes = dplyr::coalesce(Target_Genes_Filtered, ""),
            Is_High_Connectivity_Distal_Element = dplyr::if_else(
              Total_Loops >= final_cutoff_dist, "Yes", "No"
            )
          )
      } else {
        temp_df <- distal_raw_df %>%
          dplyr::rename(
            Total_Loops = Total_Loops_Filtered,
            n_Linked_Distal = n_Linked_Distal_Filtered,
            n_Linked_Promoters = n_Linked_Promoters_Filtered,
            Dominant_Interaction = Dominant_Interaction_Filtered,
            Target_Genes = Target_Genes_Filtered
          ) %>%
          dplyr::mutate(
            Is_High_Connectivity_Distal_Element = dplyr::if_else(
              Total_Loops >= final_cutoff_dist, "Yes", "No"
            )
          )
      }
      temp_df <- temp_df %>%
        dplyr::select(-any_of(c(
          "chr", "start", "end", "cluster_id",
          "Distal_Type", "Distal_Type_Filtered", "Total_Loops_Filtered",
          "Target_Genes_Filtered", "n_Linked_Distal_Filtered",
          "n_Linked_Promoters_Filtered", "Dominant_Interaction_Filtered"
        )))
      distal_element_df <- temp_df %>%
        dplyr::left_join(anchor_map,
          by = c("Distal_Anchor_ID" = "anchor_id")
        ) %>%
        dplyr::select(
          chr, start, end, cluster_id, Total_Loops,
          n_Linked_Promoters, n_Linked_Distal, Dominant_Interaction,
          any_of("Is_High_Connectivity_Distal_Element"), Target_Genes
        ) %>%
        dplyr::filter(Total_Loops > 0) %>%
        dplyr::arrange(dplyr::desc(Total_Loops))
    }
  }

  list(
    promoter_centric = promoter_centric_df,
    distal_element = distal_element_df
  )
}


#' Internal: Load Expression Matrix
#'
#' Reads a normalized expression matrix (TPM/FPKM), sets gene identifiers as
#' row names, extracts the specified sample columns, and returns per-gene mean
#' expression values.
#'
#' @param expr_matrix_file Character. Path to the expression matrix file.
#' @param sample_columns Character or integer vector. Sample columns to average.
#' @return Named numeric vector of per-gene mean expression values.
#' @importFrom data.table fread
#' @keywords internal
load_expression_matrix <- function(expr_matrix_file, sample_columns = NULL) {
  if (!file.exists(expr_matrix_file)) {
    stop("Expression matrix file not found: ", expr_matrix_file)
  }
  d <- as.data.frame(data.table::fread(expr_matrix_file))
  rownames(d) <- d[[1]]
  d <- d[, -1, drop = FALSE]

  if (is.null(sample_columns)) {
    sub_mat <- d
  } else if (is.character(sample_columns)) {
    sub_mat <- d[, intersect(sample_columns, colnames(d)), drop = FALSE]
  } else {
    sub_mat <- d[, sample_columns, drop = FALSE]
  }
  if (ncol(sub_mat) == 0) stop("No valid sample columns found in expression matrix.")

  vals <- if (ncol(sub_mat) > 1) rowMeans(sub_mat, na.rm = TRUE) else sub_mat[, 1]
  names(vals) <- rownames(d)
  vals
}


#' Internal: Generate Colors
#'
#' Helper to generate a vector of n colors.
#'
#' @param n Integer. Number of colors.
#' @param palette_input Character. Palette name or custom colors.
#' @return Hex color codes.
#'
#' @importFrom RColorBrewer brewer.pal.info brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom scales hue_pal
#' @keywords internal
get_colors <- function(n, palette_input) {
  safe_n <- max(1, n)

  if (length(palette_input) == 1 && palette_input %in% row.names(RColorBrewer::brewer.pal.info)) {
    max_avail <- RColorBrewer::brewer.pal.info[palette_input, "maxcolors"]
    pal <- RColorBrewer::brewer.pal(min(safe_n, max_avail), palette_input)
    cols <- grDevices::colorRampPalette(pal)(safe_n)
    return(if (n == 0) character(0) else cols)
  } else if (length(palette_input) >= 1) {
    if (length(palette_input) < safe_n) {
      cols <- rep(palette_input, length.out = safe_n)
    } else {
      cols <- palette_input[seq_len(safe_n)]
    }
    return(if (n == 0) character(0) else cols)
  } else {
    cols <- scales::hue_pal()(safe_n)
    return(if (n == 0) character(0) else cols)
  }
}


#' Internal: Draw Karyotype Heatmap
#'
#' Creates a genome-wide heatmap of genomic feature density (e.g., loops) across chromosomes,
#' binned by a fixed window size, and rendered as a karyotype plot.
#'
#' @param gr_data (GRanges) Genomic ranges to visualize (e.g., loop anchors).
#' @param title_prefix (character) Subtitle descriptor (e.g., sample name).
#' @param bin_size (integer) Bin width in base pairs (e.g., 1e6 for 1 Mb).
#' @param sat_level (numeric) Quantile (0–1) for color saturation (e.g., 0.95).
#' @param ref_txdb (TxDb or similar) Reference genome annotation for chromosome lengths.
#' @param plot_species (character) Genome build/species code (e.g., "hg38", "mm10").
#' @param unit_label (character) Unit for load annotation (e.g., "loops").
#' @param custom_colors (character vector) Optional custom color palette.
#' @keywords internal
#' @importFrom GenomeInfoDb seqinfo seqlevelsStyle keepSeqlevels seqlengths seqlevels
#' @importFrom GenomicRanges GRanges tileGenome countOverlaps
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors mcols
#' @importFrom grDevices colorRampPalette
#' @importFrom karyoploteR plotKaryotype kpRect getDefaultPlotParams
#' @importFrom fields image.plot
#' @return A \code{looplook_karyo} object wrapping a rendered PNG payload. Use
#'   \code{print()} to display.
draw_karyo_heatmap_internal <- function(gr_data, title_prefix, bin_size, sat_level, ref_txdb, plot_species, unit_label, custom_colors = NULL) {
  standard_chroms <- paste0("chr", c(seq_len(22), "X", "Y"))
  if (grepl("mm", plot_species, fixed = TRUE)) standard_chroms <- paste0("chr", c(seq_len(19), "X", "Y"))

  std_seqinfo <- GenomeInfoDb::seqinfo(ref_txdb)
  try(
    {
      GenomeInfoDb::seqlevelsStyle(gr_data) <- "UCSC"
    },
    silent = TRUE
  )

  existing <- intersect(GenomeInfoDb::seqlevels(gr_data), standard_chroms)
  if (length(existing) == 0) {
    return(invisible(NULL))
  }

  gr_data <- GenomeInfoDb::keepSeqlevels(gr_data, existing, pruning.mode = "coarse")
  GenomeInfoDb::seqlevels(gr_data) <- standard_chroms
  common <- intersect(GenomeInfoDb::seqlevels(gr_data), GenomeInfoDb::seqlevels(std_seqinfo))
  GenomeInfoDb::seqlengths(gr_data)[common] <- GenomeInfoDb::seqlengths(std_seqinfo)[common]
  valid_chroms <- intersect(standard_chroms, GenomeInfoDb::seqlevels(std_seqinfo))

  if (length(valid_chroms) > 0) {
    full_genome_gr <- GenomicRanges::GRanges(seqnames = valid_chroms, ranges = IRanges::IRanges(start = 1, end = GenomeInfoDb::seqlengths(std_seqinfo)[valid_chroms]))
    GenomeInfoDb::seqinfo(full_genome_gr) <- std_seqinfo[valid_chroms]
    tiles <- GenomicRanges::tileGenome(GenomeInfoDb::seqinfo(full_genome_gr), tilewidth = bin_size, cut.last.tile.in.chrom = TRUE)
    hits <- GenomicRanges::countOverlaps(tiles, gr_data)

    bin_size_mb <- bin_size / 1e6
    median_val <- median(hits[hits > 0], na.rm = TRUE)
    if (is.na(median_val)) median_val <- 0

    heatmap_colors <- if (is.null(custom_colors)) c("#FFFFFF", "#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C", "#FC4E2A", "#E31A1C", "#BD0026", "#800026", "#000000") else custom_colors

    if (max(hits) == 0) {
      max_load <- 0
      S4Vectors::mcols(tiles)$color <- "white"
      cols <- c("white")
    } else {
      cutoff <- as.numeric(quantile(hits[hits > 0], probs = sat_level, names = FALSE))
      if (is.na(cutoff) || cutoff < 1) cutoff <- max(hits)
      max_load <- round(cutoff / bin_size_mb, 1)
      capped <- ifelse(hits > cutoff, cutoff, hits)

      col_func <- grDevices::colorRampPalette(heatmap_colors)
      cols <- col_func(100)

      idx <- ceiling((capped / cutoff) * 99) + 1
      idx[hits == 0] <- 1
      S4Vectors::mcols(tiles)$color <- cols[idx]
    }

    # Render once to PNG and keep the bytes in-memory so deferred report
    # rendering does not depend on a temp file surviving until a later chunk.
    f <- tempfile(fileext = ".png")
    grDevices::png(f, width = 10, height = 8, units = "in", res = 150)

    graphics::par(oma = c(2, 2, 6, 2))
    pp <- karyoploteR::getDefaultPlotParams(plot.type = 1)
    pp$leftmargin <- 0.08
    pp$rightmargin <- 0.08
    pp$data1height <- 100
    kp <- karyoploteR::plotKaryotype(genome = plot_species, plot.type = 1, chromosomes = valid_chroms, plot.params = pp, main = NULL)
    karyoploteR::kpRect(kp, data = tiles, y0 = 0, y1 = 1, col = S4Vectors::mcols(tiles)$color, border = NA)
    main_title <- paste0("Loop Analysis: ", title_prefix, "\n(Genomic Load: Median ~", round(median_val / bin_size_mb, 1), " ", unit_label, "/MB)")
    graphics::mtext(main_title, side = 3, line = 1, outer = TRUE, cex = 1.2, font = 2)
    fields::image.plot(
      legend.only = TRUE, zlim = c(0, max_load), col = cols,
      legend.lab = paste0("Load (", unit_label, "/MB)"), legend.mar = 4.5, smallplot = c(0.88, 0.91, 0.3, 0.7)
    )

    grDevices::dev.off()
    png_raw <- NULL
    if (file.exists(f)) {
      png_raw <- tryCatch(
        readBin(f, what = "raw", n = file.info(f)$size),
        error = function(e) NULL
      )
      if (!is.null(png_raw)) {
        unlink(f)
      }
    }

    payload <- list(type = "karyo_heatmap")
    if (!is.null(png_raw)) {
      payload$png_raw <- png_raw
    } else if (file.exists(f)) {
      payload$file <- f
    }

    structure(payload, class = "looplook_karyo")
  }
}

#' Print looplook karyogram
#'
#' Displays a previously captured karyotype heatmap. Uses embedded PNG bytes
#' when available, otherwise falls back to a file-backed image.
#'
#' @param x A \code{looplook_karyo} object.
#' @param ... Additional arguments (unused).
#' @export
print.looplook_karyo <- function(x, ...) {
  if (!is.null(x$png_raw)) {
    f <- tempfile(fileext = ".png")
    writeBin(x$png_raw, f)
    on.exit(unlink(f), add = TRUE)
    if (requireNamespace("png", quietly = TRUE)) {
      img <- png::readPNG(f)
      grid::grid.newpage()
      grid::grid.raster(img, interpolate = TRUE)
    } else {
      utils::browseURL(f)
    }
    return(invisible(x))
  }

  if (is.null(x$file) || !file.exists(x$file)) {
    stop("Karyo image is unavailable.", call. = FALSE)
  }
  if (requireNamespace("png", quietly = TRUE)) {
    img <- png::readPNG(x$file)
    grid::grid.newpage()
    grid::grid.raster(img, interpolate = TRUE)
  } else {
    utils::browseURL(x$file)
  }
  invisible(x)
}

#' Draw Expression Violin Plot
#'
#' Creates a violin + boxplot showing log2-transformed gene expression
#' grouped by loop type (e.g., promoter-enhancer, enhancer-enhancer).
#'
#' @param plot_data (data.frame) Must contain columns: `loop_type`, `expression_value`.
#' @param project_name (character) Project or sample name for plot title.
#' @param filename (character) Output file path (e.g., "expr_violin.pdf").
#' @param unit_type (character) Expression unit (e.g., "TPM", "FPKM").
#' @param group_colors (character vector) Named or ordered colors for each `loop_type`.
#' @keywords internal
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot scale_fill_manual theme_minimal theme element_text element_blank labs ggsave
#' @return A `ggplot` object of the violin + boxplot.
draw_expression_violin <- function(plot_data, project_name, filename = NULL, unit_type, group_colors) {
  # 1. Log transform
  plot_data$log_expr <- log2(plot_data$expression_value + 0.01)


  plot_data_unique <- plot_data %>%
    dplyr::select(loop_type, loop_genes, log_expr) %>%
    dplyr::distinct(loop_type, loop_genes, .keep_all = TRUE)


  count_df <- plot_data_unique %>%
    dplyr::group_by(loop_type) %>%
    dplyr::summarise(n = dplyr::n()) %>%
    dplyr::ungroup()


  new_labels <- stats::setNames(
    paste0(count_df$loop_type, "\n(n=", count_df$n, ")"),
    count_df$loop_type
  )


  p <- ggplot2::ggplot(plot_data_unique, ggplot2::aes(x = loop_type, y = log_expr, fill = loop_type)) +
    ggplot2::geom_violin(scale = "width", trim = FALSE, alpha = 0.7, color = "grey30", linewidth = 0.3) +
    ggplot2::geom_boxplot(width = 0.15, fill = "white", color = "black", outlier.shape = NA, alpha = 0.9) +
    ggplot2::scale_fill_manual(values = group_colors) +
    ggplot2::scale_x_discrete(labels = new_labels) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, face = "bold", color = "black", size = 11),
      panel.grid.major.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    ) +
    ggplot2::labs(
      title = paste0(project_name, ": Gene Expression by Loop Type"),
      x = NULL,
      y = paste0("Expression (Log2 ", unit_type, ")")
    )

  return(p)
}

#' Draw Enhancer Source Distribution
#'
#' Bar plot showing the origin of functional enhancer anchors
#' (e.g., native vs. promoter-derived).
#'
#' @param loop_data (data.frame) Must contain:
#'   - `functional_anchor1_type`, `anchor1_source`
#'   - `functional_anchor2_type`, `anchor2_source`
#' @param project_name (character) Project name for title.
#'   Only anchors where type == "E" are considered.
#' @param filename (character) Output file path (e.g., "enh_sources.pdf").
#' @keywords internal
#' @importFrom dplyr select filter group_by summarise ungroup mutate arrange bind_rows
#' @importFrom ggplot2 ggplot aes geom_bar geom_text scale_fill_manual scale_y_continuous theme_classic theme element_text element_blank labs ggsave
#' @importFrom scales percent
#' @return A `ggplot` object of the bar plot showing enhancer source distribution.
draw_enhancer_source_distribution <- function(loop_data, project_name, filename = NULL) {
  a1 <- loop_data %>% dplyr::select(type = functional_anchor1_type, source = anchor1_source)
  a2 <- loop_data %>% dplyr::select(type = functional_anchor2_type, source = anchor2_source)
  plot_data <- dplyr::bind_rows(a1, a2) %>%
    dplyr::filter(type == "E") %>%
    dplyr::group_by(source) %>%
    dplyr::summarise(Count = dplyr::n(), .groups = "drop") %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Percentage = Count / sum(Count)) %>%
    dplyr::arrange(dplyr::desc(Count))

  source_colors <- c(
    "Native" = "#E0E0E0",
    "Promoter-derived enhancer" = "#D95F02",
    "Gene_body-derived enhancer" = "#66A61E"
  )
  avail_sources <- unique(plot_data$source)
  final_colors <- source_colors[avail_sources]
  if (anyNA(final_colors)) {
    final_colors <- get_colors(length(avail_sources), "Set2")
    names(final_colors) <- avail_sources
  }

  plot_data$Label <- paste0(plot_data$Count, "\n(", scales::percent(plot_data$Percentage, accuracy = 0.1), ")")

  p <- ggplot(plot_data, aes(x = reorder(source, -Count), y = Count, fill = source)) +
    geom_bar(stat = "identity", width = 0.7, color = "black", alpha = 0.8) +
    geom_text(aes(label = Label), vjust = -0.5, fontface = "bold", size = 4) +
    scale_fill_manual(values = final_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1, face = "bold", color = "black"),
      axis.title.x = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    ) +
    labs(
      title = paste0(project_name, ": Origin of Functional Enhancers"),
      y = "Number of Anchors",
      subtitle = "Breakdown of anchors classified as Enhancers after filtering"
    )
  return(p)
}


#' Internal: Draw Rose Chart (Loop Counts)
#'
#' Visualizes the proportion of loop types based on LOOP COUNTS.
#'
#' @param data_df Data frame containing loop information.
#' @param project_name Character string for the project title.
#' @param filename Character string for the output filename.
#' @param color_vec Named character vector for colors.
#'
#' @return A ggplot object representing the rose plot.
#'
#' @keywords internal
draw_rose_plot <- function(data_df, project_name, filename = NULL, color_vec) {
  color_vec <- as.character(color_vec)
  rose_data <- data_df %>%
    dplyr::group_by(loop_type) %>%
    dplyr::summarise(Count = dplyr::n()) %>%
    dplyr::mutate(Percentage = Count / sum(Count)) %>%
    dplyr::arrange(dplyr::desc(Count))
  rose_data$Label <- paste0(rose_data$Count, "\n(", scales::percent(rose_data$Percentage, 0.1), ")")

  p <- ggplot2::ggplot(rose_data, ggplot2::aes(x = reorder(loop_type, -Count), y = Count, fill = loop_type)) +
    ggplot2::geom_bar(stat = "identity", width = 1, color = "white") +
    ggplot2::coord_polar(theta = "x", start = 0) +
    ggplot2::scale_fill_manual(values = color_vec) +
    ggplot2::geom_text(ggplot2::aes(y = Count, label = Label), size = 3.5, color = "black", vjust = -0.5) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(axis.title = ggplot2::element_blank(), axis.text = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(), legend.position = "right") +
    ggplot2::labs(title = paste0(project_name, ": Loop Proportion (By Count)"))
  return(p)
}

#' Internal: Draw Circular Bar Plot (Gene Counts)
#'
#' Optimized for strictly vertical label alignment using y=0 anchor.
#'
#' @param data_df Data frame containing loop and gene information.
#' @param project_name Character string for the project title.
#' @param filename Character string for the output filename.
#' @param color_vec Named character vector for colors.
#'
#' @return A ggplot object representing the circular bar plot.
#'
#' @keywords internal
draw_circular_bar_plot <- function(data_df, project_name, filename = NULL, color_vec) {
  color_vec <- as.character(color_vec)
  circ_data <- data_df %>%
    dplyr::filter(!is.na(loop_genes) & loop_genes != "") %>%
    tidyr::separate_rows(loop_genes, sep = ";") %>%
    dplyr::group_by(loop_type) %>%
    dplyr::summarise(Unique_Gene_Count = dplyr::n_distinct(trimws(loop_genes))) %>%
    dplyr::arrange(Unique_Gene_Count) %>%
    dplyr::mutate(Label_Text = paste0(loop_type, " : ", Unique_Gene_Count))

  circ_data$loop_type <- factor(circ_data$loop_type, levels = circ_data$loop_type)
  if (nrow(circ_data) == 0) {
    return(NULL)
  }
  max_gene_count <- max(circ_data$Unique_Gene_Count, na.rm = TRUE)

  p <- ggplot2::ggplot(circ_data, ggplot2::aes(x = loop_type, fill = loop_type)) +
    ggplot2::geom_col(ggplot2::aes(y = max_gene_count), width = 0.05, fill = "grey92", color = NA) +
    ggplot2::geom_col(ggplot2::aes(y = Unique_Gene_Count), width = 0.8, color = "white", linewidth = 0.2) +
    ggplot2::geom_text(ggplot2::aes(y = Unique_Gene_Count + max_gene_count * 0.02, label = Label_Text), hjust = 0, size = 3.5, fontface = "bold") +
    ggplot2::coord_polar(theta = "y", start = 0, clip = "off") +
    ggplot2::scale_fill_manual(values = color_vec, name = "Loop Type") +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(axis.title = ggplot2::element_blank(), axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(), legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"), plot.margin = ggplot2::margin(t = 20, r = 100, b = 20, l = 20, unit = "pt")) +
    ggplot2::scale_y_continuous(limits = c(-max_gene_count * 0.4, max_gene_count * 1.3)) +
    ggplot2::labs(title = paste0(project_name, ": Unique Target Genes (Ascending)"))

  return(p)
}

#' Internal: Draw Comparison Bar Chart
#'
#' Visualizes the comparison of loop counts between original and filtered datasets.
#'
#' @param original_df Data frame containing the original loop data.
#' @param filtered_df Data frame containing the filtered loop data.
#' @param filename Character string for the output filename.
#' @param color_vec Named character vector for colors.
#'
#' @return A ggplot object visualizing the comparison between groups.
#'
#' @keywords internal
draw_comparison_bar <- function(original_df, filtered_df, filename = NULL, color_vec) {
  fmt_type <- function(x) {
    return(x)
  }
  df_compare <- dplyr::bind_rows(
    original_df %>% dplyr::mutate(loop_type = fmt_type(loop_type)) %>% dplyr::group_by(loop_type) %>% dplyr::summarise(Count = dplyr::n()) %>% dplyr::mutate(Stage = "Original"),
    filtered_df %>% dplyr::group_by(loop_type) %>% dplyr::summarise(Count = dplyr::n()) %>% dplyr::mutate(Stage = "Filtered")
  ) %>% dplyr::mutate(Stage = factor(Stage, levels = c("Original", "Filtered")))

  p <- ggplot2::ggplot(df_compare, ggplot2::aes(x = loop_type, y = Count, fill = loop_type, alpha = Stage)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge", color = "black") +
    ggplot2::scale_fill_manual(values = color_vec) +
    ggplot2::scale_alpha_manual(values = c(0.4, 1.0)) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::labs(title = "Loop Classification Changes", y = "Number of Loops", x = "Loop Type")
  return(p)
}

#' Internal: Draw Target Genomic Annotation Distribution (Pie Chart)
#'
#' Optimized: Labels only show Count (%), and hides labels for small slices (<2%) to avoid overlap.
#'
#' @param bed_info Data frame containing annotation information.
#' @param project_name Character string for the project title.
#' @param filename Character string for the output filename.
#' @param color_palette Character string for the color palette (default: "Set2").
#'
#' @return A ggplot object representing the annotation pie chart.
#'
#' @keywords internal
draw_target_annotation_pie <- function(bed_info, project_name, filename = NULL, color_palette = "Set2") {
  plot_data <- bed_info %>%
    dplyr::mutate(Feature = simplify_annotation(annotation)) %>%
    dplyr::group_by(Feature) %>%
    dplyr::summarise(Count = dplyr::n()) %>%
    dplyr::mutate(Percentage = Count / sum(Count)) %>%
    dplyr::arrange(dplyr::desc(Count))

  plot_data <- plot_data %>%
    dplyr::mutate(
      Label_Text = paste0(Count, "\n(", scales::percent(Percentage, 0.1), ")"),
      Final_Label = ifelse(Percentage >= 0.02, Label_Text, "")
    )


  n_groups <- nrow(plot_data)
  if (!exists("get_colors")) custom_colors <- scales::hue_pal()(n_groups) else custom_colors <- get_colors(n_groups, color_palette)


  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = "", y = Count, fill = Feature)) +
    ggplot2::geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.5) +
    ggplot2::coord_polar("y", start = 0) +
    ggplot2::geom_text(
      ggplot2::aes(label = Final_Label),
      position = ggplot2::position_stack(vjust = 0.5),
      size = 3.5,
      fontface = "bold",
      color = "black"
    ) +
    ggplot2::scale_fill_manual(values = custom_colors, name = "Genomic Feature") +
    ggplot2::theme_void(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", margin = ggplot2::margin(b = 10)),
      legend.position = "right",
      legend.title = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::labs(title = paste0(project_name, ": Target Peak Genomic Distribution"))

  return(p)
}

#' Internal: Draw Target-Loop Connectivity
#'
#' Visualizes how many target peaks overlap with loops, and breakdown by Loop Type.
#'
#' @param bed_info Data frame containing peak information.
#' @param cluster_info Data frame containing cluster information (currently unused in plot but kept for consistency).
#' @param project_name Character string for the project title.
#' @param filename Character string for the output filename.
#' @param color_palette Character string for the color palette (default: "Set2").
#'
#' @return A ggplot object representing the connectivity bar chart.
#'
#' @keywords internal
draw_target_connectivity_bar <- function(bed_info, cluster_info, project_name, filename = NULL, color_palette = "Set2") {
  summary_df <- bed_info %>%
    dplyr::mutate(Status = ifelse(!is.na(loop_genes_Total) & loop_genes_Total != "", "Connected (Has Loop)", "Orphan (No Loop)")) %>%
    dplyr::group_by(Status) %>%
    dplyr::summarise(Count = dplyr::n()) %>%
    dplyr::mutate(Percentage = Count / sum(Count))

  p1 <- ggplot2::ggplot(summary_df, ggplot2::aes(x = Status, y = Count, fill = Status)) +
    ggplot2::geom_col(width = 0.6, color = "black") +
    ggplot2::geom_text(ggplot2::aes(label = Count), vjust = -0.5, fontface = "bold") +
    ggplot2::scale_fill_brewer(palette = "Set2") +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(legend.position = "none", axis.title.x = ggplot2::element_blank()) +
    ggplot2::labs(title = "Target Connectivity", y = "Number of Peaks")

  return(p1)
}


#' Internal: Draw Target-Associated Loop Distribution (Donut Chart)
#'
#' Visualizes the distribution of loop types that are connected to target peaks.
#'
#' @param loop_data Data frame containing loop information.
#' @param project_name Character string for the project title.
#' @param filename Character string for the output filename.
#' @param color_vec Named character vector for colors.
#'
#' @return A ggplot object representing the donut chart of loop type distribution.
#'
#' @keywords internal
draw_target_loop_donut <- function(loop_data, project_name, filename = NULL, color_vec) {
  plot_data <- loop_data %>%
    dplyr::group_by(loop_type) %>%
    dplyr::summarise(Count = dplyr::n()) %>%
    dplyr::mutate(Percentage = Count / sum(Count)) %>%
    dplyr::arrange(dplyr::desc(Count))

  plot_data$Label <- paste0(
    plot_data$loop_type, "\n",
    plot_data$Count, " (", scales::percent(plot_data$Percentage, 0.1), ")"
  )


  plot_data$ymax <- cumsum(plot_data$Percentage)
  plot_data$ymin <- c(0, head(plot_data$ymax, n = -1))
  plot_data$labelPosition <- (plot_data$ymax + plot_data$ymin) / 2


  p <- ggplot2::ggplot(plot_data, ggplot2::aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = loop_type)) +
    ggplot2::geom_rect(color = "white") +
    ggplot2::geom_text(x = 4.2, ggplot2::aes(y = labelPosition, label = Label), size = 3.5, color = "black") +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::scale_fill_manual(values = color_vec) +
    ggplot2::xlim(c(2, 4.5)) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none", plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")) +
    ggplot2::labs(title = paste0(project_name, ": Loops Connected to Targets"))

  return(p)
}


#' Internal: Simplify Genomic Annotation to Broad Categories
#'
#' Collapses detailed ChIPseeker annotation strings into five broad
#' categories: Promoter, Intron, Exon, Distal Intergenic, and Downstream.
#' Anything unrecognised is labelled "Others".
#'
#' @param x Character vector of annotation strings.
#' @return Character vector of simplified categories.
#' @keywords internal
simplify_annotation <- function(x) {
  vapply(x, function(s) {
    if (grepl("Promoter", s, ignore.case = TRUE)) {
      return("Promoter")
    }
    if (grepl("Intron", s, ignore.case = TRUE)) {
      return("Intron")
    }
    if (grepl("Exon", s, ignore.case = TRUE)) {
      return("Exon")
    }
    if (grepl("Intergenic", s, ignore.case = TRUE)) {
      return("Distal Intergenic")
    }
    if (grepl("Downstream", s, ignore.case = TRUE)) {
      return("Downstream")
    }
    return("Others")
  }, FUN.VALUE = character(1))
}


#' Internal: Draw Pie Chart with Outside Labels
#'
#' Simplified pie chart with labels placed outside the slices, using
#' RColorBrewer palettes for genomic annotation categories.
#'
#' @param data_df Data frame with an annotation column.
#' @param group_col Character. Column name for grouping.
#' @param title Character. Plot title.
#' @param palette Character. RColorBrewer palette name.
#' @return A ggplot object, or NULL if data is empty.
#' @importFrom ggplot2 ggplot aes geom_bar geom_segment geom_text coord_polar
#'   xlim scale_fill_brewer theme_void labs theme element_text
#' @keywords internal
draw_pie_with_outside_labels <- function(data_df, group_col, title, palette) {
  if (is.null(data_df) || nrow(data_df) == 0) {
    return(NULL)
  }

  plot_data <- data_df
  plot_data$Simplified <- simplify_annotation(plot_data[[group_col]])

  stats <- plot_data %>%
    dplyr::group_by(Simplified) %>%
    dplyr::summarise(Count = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(
      Fraction = Count / sum(Count),
      LabelText = ifelse(Fraction >= 0.01, paste0(Count, " (", round(Fraction * 100, 1), "%)"), "")
    ) %>%
    dplyr::arrange(dplyr::desc(Simplified))

  if (nrow(stats) == 0) {
    return(NULL)
  }

  stats <- stats %>%
    dplyr::mutate(
      ymax = cumsum(Fraction),
      ymin = c(0, head(ymax, n = -1)),
      ypos = (ymax + ymin) / 2,
      hjust = ifelse(ypos < 0.5, 0, 1)
    )

  ggplot2::ggplot(stats, ggplot2::aes(y = Fraction, fill = Simplified)) +
    ggplot2::geom_bar(ggplot2::aes(x = 1), width = 1, stat = "identity", color = "white") +
    ggplot2::geom_segment(
      data = subset(stats, LabelText != ""),
      ggplot2::aes(x = 1.51, xend = 1.62, y = ypos, yend = ypos),
      color = "grey50", linewidth = 0.5
    ) +
    ggplot2::geom_text(
      data = subset(stats, LabelText != ""),
      ggplot2::aes(x = 1.65, y = ypos, label = LabelText, hjust = hjust),
      size = 3.5, fontface = "bold", check_overlap = FALSE
    ) +
    ggplot2::coord_polar("y", start = 0) +
    ggplot2::xlim(0.5, 2.5) +
    ggplot2::scale_fill_brewer(palette = palette, name = "Genomic Feature") +
    ggplot2::theme_void() +
    ggplot2::labs(title = title) +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 12)
    )
}
