#' @title Integrative functional annotation and profiling of target genes
#'
#' @description
#' Integrates 3D genomic interaction data (e.g., Hi-C, HiChIP) with transcriptomic
#' profiles (RNA-seq) to evaluate the functional and regulatory landscape of target genes.
#' The workflow sequentially performs differential expression profiling, gene set enrichment,
#' transcription factor motif scanning, Gene Ontology (GO) enrichment, and
#' protein-protein interaction (PPI) network construction.
#'
#' @details
#' Two analysis steps use R's random number generator: GSEA target-gene down-sampling
#' (controlled by \code{gsea_nSample}) and motif background loop sampling (limited to
#' 2 000 background regions). For fully reproducible results, call \code{set.seed()}
#' before running this function. The \code{clusterProfiler::GSEA} call internally
#' uses \code{seed = TRUE} and is not affected by the global RNG state.
#'
#' @param annotation_res List. The result object returned by \code{\link{annotate_peaks_and_loops}}.
#' @param diff_file Character. Path to the differential expression file (CSV/TSV).
#' @param lfc_col Character. The column name in \code{diff_file} representing Log2 Fold Change.
#' @param expr_matrix_file Character. Path to the normalized expression matrix.
#' @param metadata_file Character. Path to the sample metadata file.
#' @param target_source Character vector. Source of target genes to analyze.
#' @param target_mapping_mode Character. Specifies the mapping strategy for target genes.
#' @param loop_types Character vector. The specific loop types to analyze.
#' @param include_Filled Logical. If \code{TRUE}, utilizes the comprehensively merged gene assignment.
#' @param use_nearest_gene Logical. If \code{TRUE}, bypasses 3D loop-based gene assignment.
#' @param group_order Character vector. Optional factor levels to sort sample groups.
#' @param project_name Character. Prefix for all output files and plot titles.
#' @param org_db Character. Organism annotation database (e.g., "org.Hs.eg.db").
#' @param run_motif Logical. Whether to perform Transcription Factor Binding Site motif analysis.
#' @param genome_id Character. Reference genome assembly for motif sequence extraction.
#' @param motif_p_thresh Numeric. P-value threshold for scanning.
#' @param motif_ntop Numeric. Number of top enriched motifs to output.
#' @param run_go Logical. Whether to perform Gene Ontology (GO) enrichment.
#' @param run_ppi Logical. Whether to construct Protein-Protein Interaction networks.
#' @param ppi_score Numeric. Minimum combined confidence score for STRING edges.
#' @param ppi_nSample Numeric. Maximum number of genes to include in PPI.
#' @param heatmap_nSample Numeric. Maximum number of genes to plot in heatmap.
#' @param gsea_nSample Numeric. Maximum number of target genes to sample for GSEA.
#' @param cnet_nSample Numeric. Number of top GO terms to display in cnetplot.
#' @param stat_test Character. Statistical test for LFC comparisons.
#' @param cor_method Character. Method for sample correlation matrices.
#'
#' @return An invisible nested list containing interactive plot objects, GO results, and gene sets.
#'
#' @examples
#' \donttest{
#' rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
#' diff_path <- system.file("extdata", "example_deg.txt", package = "looplook")
#' expr_path <- system.file("extdata", "example_tpm.txt", package = "looplook")
#' meta_path <- system.file("extdata", "example_coldata.txt", package = "looplook")
#' tmp <- new.env()
#' load(rdata_path, envir = tmp)
#' res <- tmp[[ls(tmp)[1]]]
#' profile_target_genes(
#'   annotation_res = res,
#'   diff_file  = diff_path,
#'   expr_matrix_file = expr_path,
#'   metadata_file = meta_path,
#'   run_go  = TRUE,
#'   run_ppi = FALSE
#' )
#' }
#'
#' @export
profile_target_genes <- function(
  annotation_res,
  diff_file,
  lfc_col = "log2FoldChange",
  expr_matrix_file,
  metadata_file,
  target_source = c("loops", "targets"),
  target_mapping_mode = c("all", "promoter"),
  loop_types = c("E-P", "P-P"),
  include_Filled = TRUE,
  use_nearest_gene = FALSE,
  group_order = NULL,
  project_name = "Analysis",
  org_db = "org.Hs.eg.db",
  run_motif = FALSE,
  genome_id = "hg38",
  motif_p_thresh = 1e-4,
  motif_ntop = 5,
  run_go = TRUE,
  run_ppi = FALSE,
  ppi_score = 400,
  ppi_nSample = 400,
  heatmap_nSample = 99999,
  gsea_nSample = 99999,
  cnet_nSample = 50,
  stat_test = "wilcox.test",
  cor_method = "pearson"
) {
  target_source <- match.arg(target_source, several.ok = TRUE)
  target_mapping_mode <- match.arg(target_mapping_mode)

  root_project_name <- project_name
  if (use_nearest_gene && "targets" %in% target_source) {
    root_project_name <- paste0(root_project_name, "_RefNearest")
  } else if ("targets" %in% target_source) {
    if (target_mapping_mode == "promoter") root_project_name <- paste0(root_project_name, "_Promoter")
    if (!include_Filled) root_project_name <- paste0(root_project_name, "_LoopOnly")
  }

  message(">>> Analysis Init | Root Project: ", root_project_name)

  # 仅保留极其核心的跨平台依赖检查，去除泛滥的 requireNamespace
  if (!requireNamespace(org_db, quietly = TRUE)) stop("Package ", org_db, " missing. Please install it.")

  if (run_motif) {
    bs_pkg <- switch(genome_id, "hg19" = "BSgenome.Hsapiens.UCSC.hg19", "hg38" = "BSgenome.Hsapiens.UCSC.hg38", "mm9" = "BSgenome.Mmusculus.UCSC.mm9", "mm10" = "BSgenome.Mmusculus.UCSC.mm10", NULL)
    if (is.null(bs_pkg) || !requireNamespace(bs_pkg, quietly = TRUE)) {
      warning("Unsupported genome or missing BSgenome package. Disabling Motif Analysis.")
      run_motif <- FALSE
    }
  }

  message("--- Reading files...")
  diff_df_raw <- read_robust_general(diff_file, header = TRUE, row_name = 1, desc = "Diff", min_cols = 1)
  tpm_mat_raw <- read_robust_general(expr_matrix_file, header = TRUE, row_name = 1, desc = "Expr", min_cols = 1)
  meta_raw <- read_robust_general(metadata_file, header = TRUE, row_name = NULL, desc = "Meta", min_cols = 1)

  if (ncol(meta_raw) >= 2) colnames(meta_raw)[c(1, 2)] <- c("SampleID", "Group")
  meta_raw$SampleID <- trimws(as.character(meta_raw$SampleID))
  if (!is.null(group_order)) meta_raw$Group <- factor(meta_raw$Group, levels = group_order)
  meta_raw <- meta_raw %>% dplyr::arrange(Group)

  if (!lfc_col %in% colnames(diff_df_raw)) stop("LFC column ", lfc_col, " not found")
  clean_diff <- diff_df_raw[!is.na(diff_df_raw[[lfc_col]]) & is.finite(diff_df_raw[[lfc_col]]), , drop = FALSE]
  global_glist <- sort(setNames(clean_diff[[lfc_col]], rownames(clean_diff)), decreasing = TRUE)

  loop_stats_df <- annotation_res$promoter_centric_stats
  final_master_list <- list()
  warn_env <- new.env(parent = emptyenv())
  warn_env$w <- character()

  safe_step <- function(label, expr) {
    tryCatch(
      expr,
      error = function(e) {
        msg <- paste0(label, " skipped: ", conditionMessage(e))
        warn_env$w <- c(warn_env$w, msg)
        warning(msg, call. = FALSE)
        NULL
      }
    )
  }

  for (src in target_source) {
    warn_env$w <- character()
    current_source_proj_name <- paste0(root_project_name, "_", src)
    message("\n================================================================")
    message(">>> Processing Source: [", src, "]")

    active_loop_types <- if (src == "loops") loop_types else NULL
    raw_gene_sets <- extract_target_gene_sets(annotation_res, src, active_loop_types, include_Filled, use_nearest_gene, target_mapping_mode)

    if (length(raw_gene_sets) == 0) {
      warning("No gene sets found. Skipping.")
      next
    }

    analysis_queue <- raw_gene_sets
    source_go_results <- list()
    source_plots <- list()

    for (task_name in names(analysis_queue)) {
      target_genes <- analysis_queue[[task_name]]
      current_proj_name <- paste0(current_source_proj_name, "_", task_name)
      target_genes <- intersect(target_genes, names(global_glist))

      message("\n--- Task: ", task_name, " (Valid Genes: ", length(target_genes), ") ---")
      if (length(target_genes) < 3) {
        message("  Too few genes (<3), skipping.")
        next
      }
      task_plots <- list()

      # 3.1 Violin
      p_vio <- safe_step("LFC violin", {
        run_lfc_violin(target_genes, global_glist, stat_test, current_proj_name)
      })
      if (!is.null(p_vio)) task_plots$LFC_Violin <- p_vio

      # 3.2 GSEA
      gsea_out <- safe_step("GSEA", {
        run_gsea_analysis(target_genes, global_glist, gsea_nSample, current_proj_name)
      })
      if (!is.null(gsea_out$plot)) task_plots$GSEA <- gsea_out$plot

      # 3.3 Visualizations (Total Loops)
      heat_plots <- safe_step("Expression/connectivity plots", {
        run_heatmap_and_connectivity(target_genes, tpm_mat_raw, meta_raw, loop_stats_df,
          global_glist, heatmap_nSample, cor_method, current_proj_name,
          source_type = src, target_col = NULL)
      })
      if (length(heat_plots) > 0) task_plots <- c(task_plots, heat_plots)

      # 3.3b Visualizations (Distal Loops)
      if (!is.null(loop_stats_df) && "n_Linked_Distal" %in% colnames(loop_stats_df)) {
        dist_plots <- safe_step("Distal connectivity plots", {
          run_heatmap_and_connectivity(target_genes, tpm_mat_raw, meta_raw, loop_stats_df,
            global_glist, heatmap_nSample, cor_method, current_proj_name,
            source_type = src, target_col = "n_Linked_Distal", skip_heatmap = TRUE)
        })
        if (length(dist_plots) > 0) task_plots <- c(task_plots, dist_plots)
      }

      # 3.4 Motif Analysis
      if (run_motif) {
        motif_plots <- safe_step("Motif analysis", {
          run_distal_motif_analysis(target_genes, annotation_res$loop_annotation,
            genome_id, motif_p_thresh, current_proj_name, motif_ntop)
        })
        if (length(motif_plots) > 0) task_plots <- c(task_plots, motif_plots)
      }

      # 3.5 GO
      if (run_go) {
        go_out <- safe_step("GO enrichment", {
          run_go_enrichment(target_genes, org_db, global_glist, cnet_nSample,
            current_proj_name)
        })
        if (!is.null(go_out$result) && nrow(go_out$result) > 0) {
          top_go <- if ("ONTOLOGY" %in% colnames(go_out$result)) {
            go_out$result %>% dplyr::group_by(ONTOLOGY) %>% dplyr::arrange(pvalue) %>% dplyr::slice_head(n = 5) %>% dplyr::ungroup()
          } else {
            head(go_out$result[order(go_out$result$pvalue), ], 15)
          }
          top_go$CleanLoopType <- task_name
          top_go$LoopType <- if ("ONTOLOGY" %in% colnames(go_out$result)) paste0(task_name, "\n(", top_go$ONTOLOGY, ")") else task_name
          top_go$Source <- src
          source_go_results[[length(source_go_results) + 1]] <- top_go
        }
        if (!is.null(go_out$plot)) task_plots$GO_Network <- go_out$plot
      }

      # 3.6 PPI
      if (run_ppi) {
        p_ppi <- safe_step("PPI network", {
          run_ppi_analysis(target_genes, global_glist, org_db, ppi_score,
            ppi_nSample, current_proj_name)
        })
        if (!is.null(p_ppi)) task_plots$PPI_Network <- p_ppi
      }

      source_plots[[task_name]] <- task_plots
    }

    if (run_go && length(source_go_results) > 0) {
      p_go_sum <- safe_step("GO summary plot", {
        plot_summary_go_lollipop(source_go_results, current_source_proj_name)
      })
      if (length(p_go_sum) > 0) source_plots$Summary_GO <- p_go_sum
    }
    final_master_list[[src]] <- list(go_results = source_go_results,
      target_gene_sets = analysis_queue, plots = source_plots,
      warnings = unique(warn_env$w))
  }
  message("\n All analysis complete.")
  return(invisible(final_master_list))
}

#' @title Robust Data Reader
#' @description Safely reads standard genomic formats utilizing `data.table::fread` for intelligent format inference.
#' @param f Character. Path to the input file.
#' @param header Logical. Whether the file contains a header row.
#' @param row_name Integer or NULL. Column index to be used as row names.
#' @param desc Character. Short description for error logging.
#' @param min_cols Integer. Minimum number of columns required.
#' @return A data frame.
#' @keywords internal
read_robust_general <- function(f, header = FALSE, row_name = NULL, desc = "file", min_cols = 3) {
  if (is.null(f) || length(f) == 0 || f == "") stop(desc, " path is empty.")
  if (!file.exists(f)) stop(desc, " not found: ", f)

  # 彻底听取审稿人建议，抛弃繁琐的 tryCatch 瞎猜，信任底层引擎
  d_dt <- data.table::fread(f, header = header, data.table = FALSE, showProgress = FALSE)

  if (!is.null(row_name) && ncol(d_dt) > 1) {
    rownames(d_dt) <- d_dt[, row_name]
    d_dt <- d_dt[, -row_name, drop = FALSE]
  }

  if (ncol(d_dt) < min_cols) {
    stop(desc, " has insufficient columns (found ", ncol(d_dt), ", required ", min_cols, ").")
  }
  return(d_dt)
}

#' @title Extract Target Gene Sets from Annotation Results
#' @description Parses loop and target annotations to extract valid gene lists.
#' @return A named list of character vectors, each containing target gene symbols.
#' @keywords internal
extract_target_gene_sets <- function(annotation_res, src, active_loop_types = NULL, include_Filled = TRUE, use_nearest_gene = FALSE, target_mapping_mode = "all") {
  raw_gene_sets <- list()

  if ("targets" %in% src && !is.null(annotation_res$target_annotation)) {
    bed_info <- annotation_res$target_annotation
    target_col <- NULL
    if (use_nearest_gene) {
      if ("SYMBOL" %in% colnames(bed_info)) target_col <- "SYMBOL"
      else if ("geneId" %in% colnames(bed_info)) target_col <- "geneId"
      if (is.null(target_col)) stop("Targets: 'SYMBOL' or 'geneId' required when use_nearest_gene is TRUE.")
    } else {
      base_col <- if (target_mapping_mode == "promoter") "Regulated_promoter_genes" else "Assigned_Target_Genes"
      desired_col <- if (include_Filled) paste0(base_col, "_Filled") else base_col
      if (desired_col %in% colnames(bed_info)) target_col <- desired_col
      else stop("Targets: Required column '", desired_col, "' not found.")
    }
    if (!is.null(target_col)) {
      gs <- clean_gene_names(bed_info[[target_col]], "[;,]")
      if (length(gs) > 0) raw_gene_sets[["Target_Genes"]] <- gs
    }
  }

  if ("loops" %in% src && !is.null(annotation_res$loop_annotation)) {
    loop_df <- annotation_res$loop_annotation
    gene_col <- "Putative_Target_Genes"
    if (!gene_col %in% colnames(loop_df)) stop("Loops: Required column '", gene_col, "' not found.")
    use_types <- if (is.null(active_loop_types)) unique(loop_df$loop_type) else intersect(active_loop_types, unique(loop_df$loop_type))
    if (length(use_types) > 0) {
      for (lt in use_types) {
        sub_df <- loop_df[loop_df$loop_type == lt, ]
        if (nrow(sub_df) > 0) {
          gs <- clean_gene_names(sub_df[[gene_col]], "[;,]")
          if (length(gs) > 0) {
            safe_name <- paste0(gsub("-", "", lt, fixed = TRUE), "_Genes")
            raw_gene_sets[[safe_name]] <- gs
          }
        }
      }
    }
  }
  return(raw_gene_sets)
}

#' @title Generate LFC Violin and Boxplot
#' @return A \code{ggplot} object, or \code{NULL} if fewer than 3 valid targets.
#' @keywords internal
run_lfc_violin <- function(target_genes, global_glist, stat_test = "wilcox.test", project_name) {
  valid_targets <- intersect(target_genes, names(global_glist))
  if (length(valid_targets) < 3) return(NULL)

  target_lfc <- global_glist[valid_targets]
  other_genes <- setdiff(names(global_glist), valid_targets)
  other_lfc <- global_glist[other_genes]

  plot_data <- data.frame(
    LFC = c(target_lfc, other_lfc),
    Group = factor(c(rep("Target", length(target_lfc)), rep("Background", length(other_lfc))), levels = c("Target", "Background"))
  )

  n_target <- length(target_lfc)
  n_back <- length(other_lfc)

  p_val <- if (stat_test == "wilcox.test") stats::wilcox.test(target_lfc, other_lfc)$p.value else stats::t.test(target_lfc, other_lfc)$p.value
  p_label <- formatC(p_val, format = "e", digits = 2)
  x_labels <- c("Target" = paste0("Target\n(n=", n_target, ")"), "Background" = paste0("Background\n(n=", n_back, ")"))

  y_min <- quantile(plot_data$LFC, 0.001, na.rm = TRUE)
  y_max <- quantile(plot_data$LFC, 0.999, na.rm = TRUE)
  y_pad <- (y_max - y_min) * 0.03
  cols <- c("Target" = "#E41A1C", "Background" = "#999999")

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Group, y = LFC, fill = Group)) +
    ggplot2::geom_violin(trim = TRUE, alpha = 0.5, color = NA) +
    ggplot2::geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.9, color = "black") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.5) +
    ggplot2::scale_x_discrete(labels = x_labels) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::coord_cartesian(ylim = c(y_min - y_pad, y_max + y_pad)) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::labs(title = project_name, subtitle = paste0("Stat: ", stat_test, ", P: ", p_label), y = "Log2 Fold Change", x = NULL) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 12), plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10), legend.position = "none")

  return(p)
}

#' @title Run Custom Gene Set Enrichment Analysis (GSEA)
#' @return A list with \code{result} (data frame) and \code{plot} (ggplot) elements.
#' @keywords internal
run_gsea_analysis <- function(target_genes, global_glist, gsea_ntop, current_proj_name) {
  curr_glist <- global_glist
  if (any(duplicated(curr_glist))) {
    curr_glist <- curr_glist + runif(length(curr_glist), 0, 1e-6)
    curr_glist <- sort(curr_glist, decreasing = TRUE)
  }
  names(curr_glist) <- toupper(names(curr_glist))
  curr_targets <- toupper(target_genes)

  if (!is.null(gsea_ntop) && length(curr_targets) > gsea_ntop) curr_targets <- sample(curr_targets, gsea_ntop)
  if (length(intersect(curr_targets, names(curr_glist))) < 2) return(list(result = NULL, plot = NULL))

  term_df <- data.frame(term = current_proj_name, gene = curr_targets, stringsAsFactors = FALSE)
  gsea_res <- clusterProfiler::GSEA(curr_glist, TERM2GENE = term_df, pvalueCutoff = 1.1, minGSSize = 2, maxGSSize = 50000, verbose = FALSE, seed = TRUE)

  p_out <- NULL
  if (!is.null(gsea_res) && nrow(as.data.frame(gsea_res)) > 0) {
    p_temp <- enrichplot::gseaplot2(gsea_res, geneSetID = 1, subplots = 1)
    d <- NULL
    if (inherits(p_temp, "ggplot")) d <- p_temp$data
    else if (inherits(p_temp, "aplot") || inherits(p_temp, "gglist") || is.list(p_temp)) {
      for (sub_p in p_temp) if (inherits(sub_p, "ggplot") && !is.null(sub_p$data) && "runningScore" %in% colnames(sub_p$data)) {
        d <- sub_p$data
        break
      }
    } else if (!is.null(p_temp$data)) d <- p_temp$data

    if (!is.null(d) && is.data.frame(d) && "runningScore" %in% colnames(d)) {
      if (!"geneList" %in% colnames(d)) d$geneList <- curr_glist[d$x]
      if (!"position" %in% colnames(d)) d$position <- as.numeric(names(curr_glist)[d$x] %in% curr_targets)
      max_rank <- max(d$x)
      res_df <- as.data.frame(gsea_res)
      nes_val <- res_df$NES[1]
      pval_val <- res_df$pvalue[1]
      main_col <- if (!is.na(nes_val) && nes_val >= 0) "#E41A1C" else "#377EB8"

      p1 <- ggplot2::ggplot(d, ggplot2::aes(x = x, y = runningScore)) +
        ggplot2::geom_line(color = main_col, linewidth = 1) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
        ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(0, max_rank)) +
        ggplot2::theme_bw() +
        ggplot2::labs(x = NULL, y = "ES", title = paste0(current_proj_name, "\nNES: ", round(nes_val, 3), "  P: ", formatC(pval_val, format = "e", digits = 2))) +
        ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(), plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 10))

      hit_data <- d[d$position == 1, ]
      p2 <- ggplot2::ggplot(hit_data, ggplot2::aes(x = x, y = 1)) +
        ggplot2::geom_segment(ggplot2::aes(xend = x, yend = 0), color = "black", alpha = 0.6) +
        ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(0, max_rank)) +
        ggplot2::scale_y_continuous(expand = c(0, 0)) +
        ggplot2::theme_void() +
        ggplot2::theme(panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5))

      p3 <- ggplot2::ggplot(d, ggplot2::aes(x = x, y = geneList)) +
        ggplot2::geom_segment(ggplot2::aes(xend = x, yend = 0, color = geneList)) +
        ggplot2::scale_color_gradient2(low = "#1B7837", mid = "white", high = "#762A83", midpoint = 0) +
        ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(0, max_rank)) +
        ggplot2::coord_cartesian(ylim = c(quantile(d$geneList, 0.005, na.rm = TRUE), quantile(d$geneList, 0.995, na.rm = TRUE))) +
        ggplot2::theme_classic() +
        ggplot2::labs(x = "Rank", y = "LFC") +
        ggplot2::theme(legend.position = "none", axis.text.y = ggplot2::element_text(size = 8))

      if (requireNamespace("aplot", quietly = TRUE)) p_out <- aplot::plot_list(p1, p2, p3, ncol = 1, heights = c(2, 0.5, 1.5))
      else p_out <- enrichplot::gseaplot2(gsea_res, geneSetID = 1, title = as.data.frame(gsea_res)$Description[1])
    }
  }
  return(list(result = as.data.frame(gsea_res), plot = p_out))
}

#' @title Perform GO Enrichment and Generate Network Plot
#' @importFrom methods slot<-
#' @return A list with \code{result} (data frame) and \code{plot} (ggplot) elements.
#' @keywords internal
run_go_enrichment <- function(genes, org_db, universe_genes, cnet_nSample = 50, project_name = "Analysis") {
  clean_genes <- clean_gene_names(genes)
  org_db_obj <- utils::getFromNamespace(org_db, org_db)

  gene_entrez <- AnnotationDbi::mapIds(org_db_obj, keys = clean_genes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  valid_entrez <- na.omit(gene_entrez)

  use_symbol_mode <- length(valid_entrez) < 5 || (length(valid_entrez) / length(clean_genes) < 0.1)

  if (use_symbol_mode) {
    final_genes <- clean_genes
    final_keytype <- "SYMBOL"
  } else {
    final_genes <- valid_entrez
    final_keytype <- "ENTREZID"
  }

  final_universe <- NULL
  if (!is.null(universe_genes)) {
    if (use_symbol_mode) {
      final_universe <- names(universe_genes)
    } else {
      univ_entrez <- AnnotationDbi::mapIds(org_db_obj, keys = names(universe_genes), column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
      final_universe <- na.omit(univ_entrez)
    }
  }

  ego <- clusterProfiler::enrichGO(gene = final_genes, universe = final_universe, OrgDb = org_db_obj, keyType = final_keytype, ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 1.0, minGSSize = 5, maxGSSize = 800, readable = (final_keytype == "ENTREZID"))

  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(list(result = NULL, plot = NULL))

  p_cnet <- NULL
  top_n <- if (!is.null(cnet_nSample)) {
    max(1L, min(as.integer(cnet_nSample), nrow(as.data.frame(ego))))
  } else {
    min(5L, nrow(as.data.frame(ego)))
  }
  fc_vec <- universe_genes

  genes_to_label <- c()
  top_df <- head(as.data.frame(ego), top_n)

  gene_to_pathways <- list()
  for (i in seq_len(nrow(top_df))) {
    gs <- unlist(strsplit(top_df$geneID[i], "/", fixed = TRUE))
    for (g in gs) gene_to_pathways[[g]] <- c(gene_to_pathways[[g]], top_df$ID[i])
  }

  for (i in seq_len(nrow(top_df))) {
    gs <- unlist(strsplit(top_df$geneID[i], "/", fixed = TRUE))
    valid_g <- intersect(gs, names(fc_vec))
    if (length(valid_g) > 0) genes_to_label <- c(genes_to_label, head(valid_g[order(abs(fc_vec[valid_g]), decreasing = TRUE)], 3))
  }

  hub_genes <- names(gene_to_pathways)[lengths(gene_to_pathways) >= 2]
  valid_hub <- intersect(hub_genes, names(fc_vec))
  if (length(valid_hub) > 0) genes_to_label <- c(genes_to_label, head(valid_hub[order(abs(fc_vec[valid_hub]), decreasing = TRUE)], 5))

  genes_to_label <- unique(genes_to_label)

  ego_df <- as.data.frame(ego)
  ego_df$Description <- vapply(ego_df$Description, function(x) paste(strwrap(x, width = 35), collapse = "\n"), FUN.VALUE = character(1))
  slot(ego, "result") <- ego_df

  options(ggrepel.max.overlaps = 100)
  p_cnet <- enrichplot::cnetplot(ego, foldChange = fc_vec, circular = FALSE, colorEdge = TRUE, showCategory = top_n, node_label = "category", repel = TRUE, cex_label_category = 1.3)

  if (length(genes_to_label) > 0 && requireNamespace("ggraph", quietly = TRUE)) {
    p_cnet <- p_cnet + ggraph::geom_node_text(ggplot2::aes(filter = name %in% genes_to_label, label = name), repel = TRUE, size = 3.5, fontface = "bold.italic", bg.color = "white", bg.r = 0.15, max.overlaps = Inf)
  }

  p_cnet <- p_cnet + ggplot2::scale_color_distiller(palette = "PuOr", name = "Log2FC") + ggplot2::labs(title = paste0("GO Network: ", project_name)) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

  return(list(result = as.data.frame(ego), plot = p_cnet))
}

#' @title Construct and Visualize STRING PPI Network
#' @importFrom utils capture.output
#' @return A \code{ggplot} object representing the PPI network, or \code{NULL} if no interactions found.
#' @keywords internal
run_ppi_analysis <- function(target_genes, global_glist, org_db, ppi_score, ppi_ntop, current_proj_name) {
  species_id <- if (grepl("Mm", org_db, ignore.case = TRUE)) 10090 else 9606
  string_db_obj <- STRINGdb::STRINGdb$new(version = "11.5", species = species_id, score_threshold = ppi_score, input_directory = "")

  ppi_genes <- target_genes
  if (!is.null(ppi_ntop) && length(ppi_genes) > ppi_ntop) {
    valid_in_lfc <- intersect(ppi_genes, names(global_glist))
    if (length(valid_in_lfc) > 0) ppi_genes <- head(valid_in_lfc[order(abs(global_glist[valid_in_lfc]), decreasing = TRUE)], ppi_ntop)
    else ppi_genes <- head(ppi_genes, ppi_ntop)
  }

  # Capture STRINGdb mapping output to prevent raw text in the report
  map_output <- capture.output(
    targets_mapped <- string_db_obj$map(data.frame(gene = ppi_genes), "gene", removeUnmappedRows = TRUE)
  )
  if (nrow(targets_mapped) == 0) return(NULL)

  hits <- targets_mapped$STRING_id
  if (length(hits) <= 1) return(NULL)

  g_string <- string_db_obj$get_subnetwork(hits)
  g_string <- igraph::delete_vertices(g_string, igraph::V(g_string)[igraph::degree(g_string) == 0])
  if (igraph::vcount(g_string) == 0) return(NULL)

  map_df <- targets_mapped[targets_mapped$STRING_id %in% igraph::V(g_string)$name, ]
  symbol_map <- setNames(map_df$gene, map_df$STRING_id)
  igraph::V(g_string)$symbol <- symbol_map[igraph::V(g_string)$name]

  lfc_vals <- setNames(as.numeric(global_glist), names(global_glist))[igraph::V(g_string)$symbol]
  lfc_vals[is.na(lfc_vals)] <- 0
  igraph::V(g_string)$lfc <- as.numeric(lfc_vals)
  igraph::V(g_string)$deg <- as.numeric(igraph::degree(g_string))

  if (is.null(igraph::E(g_string)$combined_score)) igraph::E(g_string)$combined_score <- ppi_score
  igraph::E(g_string)$combined_score <- as.numeric(igraph::E(g_string)$combined_score)

  top_n_labels <- 25
  num_nodes <- length(igraph::V(g_string)$deg)
  threshold_deg <- if (num_nodes > top_n_labels) sort(igraph::V(g_string)$deg, decreasing = TRUE)[top_n_labels] else 0
  igraph::V(g_string)$label_text <- ifelse(igraph::V(g_string)$deg >= threshold_deg, igraph::V(g_string)$symbol, NA)

  p_ppi <- ggraph::ggraph(g_string, layout = "fr") +
    ggraph::geom_edge_link(ggplot2::aes(alpha = combined_score), color = "grey60", edge_width = 0.5, show.legend = FALSE) +
    ggraph::geom_node_point(ggplot2::aes(color = lfc, size = deg), stroke = 0.5) +
    ggraph::geom_node_text(ggplot2::aes(label = label_text), repel = TRUE, size = 3.5, max.overlaps = Inf, fontface = "bold", bg.color = "white", bg.r = 0.1) +
    ggplot2::scale_color_distiller(palette = "PuOr", name = "LFC") +
    ggplot2::scale_size_continuous(range = c(2, 8), guide = "none") +
    ggraph::scale_edge_alpha_continuous(range = c(0.4, 0.9)) +
    ggraph::theme_graph(base_family = "sans", background = "white") +
    ggplot2::labs(title = paste0("PPI Network: ", current_proj_name), subtitle = paste0("Interacting Nodes: ", num_nodes, " | Score: ", min(igraph::E(g_string)$combined_score, na.rm = TRUE), "-", max(igraph::E(g_string)$combined_score, na.rm = TRUE))) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

  return(p_ppi)
}

#' @title Generate Summary GO Lollipop Facet Plot
#' @return A list of \code{ggplot} objects, one per ontology facet.
#' @keywords internal
#' @noRd
plot_summary_go_lollipop <- function(all_go_results, base_project_name) {
  plot_list <- list()
  valid_results <- Filter(function(x) !is.null(x) && is.data.frame(x) && nrow(x) > 0, all_go_results)
  if (length(valid_results) == 0) return(plot_list)

  final_go_df <- do.call(rbind, valid_results)
  if (is.null(final_go_df)) return(plot_list)

  if (!"CleanLoopType" %in% colnames(final_go_df)) final_go_df$CleanLoopType <- final_go_df$LoopType
  use_ggtext <- requireNamespace("ggtext", quietly = TRUE)

  for (ltype in unique(final_go_df$CleanLoopType)) {
    sub_df <- final_go_df %>% dplyr::filter(CleanLoopType == ltype)
    if (nrow(sub_df) == 0) next

    sub_df$logP <- -log10(sub_df$pvalue)
    scale_f <- max(sub_df$logP, na.rm = TRUE) / max(sub_df$Count, na.rm = TRUE) * 1.0

    sub_df <- sub_df %>%
      dplyr::group_by(ONTOLOGY) %>%
      dplyr::arrange(logP) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(Description_unique = factor(Description, levels = unique(Description)))

    onto_levels <- sort(unique(sub_df$ONTOLOGY))
    sub_df$ONTOLOGY <- factor(sub_df$ONTOLOGY, levels = onto_levels)
    onto_colors <- RColorBrewer::brewer.pal(max(3, length(onto_levels)), "Dark2")[seq_len(length(onto_levels))]

    if (use_ggtext) sub_df$ONTOLOGY_Plot <- factor(sub_df$ONTOLOGY, levels = onto_levels, labels = paste0("<span style='color:", onto_colors, "'>", onto_levels, "</span>"))
    else sub_df$ONTOLOGY_Plot <- sub_df$ONTOLOGY

    p_go <- ggplot2::ggplot(sub_df, ggplot2::aes(y = Description_unique)) +
      ggplot2::geom_segment(ggplot2::aes(x = 0, xend = logP, yend = Description_unique, color = ONTOLOGY), linewidth = 3) +
      ggplot2::geom_point(ggplot2::aes(x = logP, color = ONTOLOGY), size = 5) +
      ggplot2::geom_path(ggplot2::aes(x = Count * scale_f, group = 1), color = "grey60", linewidth = 1.5, linetype = "11") +
      ggplot2::geom_point(ggplot2::aes(x = Count * scale_f), color = "grey60", size = 4, shape = 17) +
      ggplot2::scale_x_continuous(name = expression(-log[10](p - value)), expand = ggplot2::expansion(mult = c(0, 0.6)), sec.axis = ggplot2::sec_axis(~ . / scale_f, name = "Gene Counts")) +
      ggplot2::facet_grid(ONTOLOGY_Plot ~ ., scales = "free_y", space = "free_y", switch = "y") +
      ggplot2::labs(y = NULL, title = paste0("GO Enrichment: ", ltype), subtitle = "Colored Dot: Significance | Grey Triangle: Gene Count") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10, color = "black"), strip.placement = "outside", strip.background = ggplot2::element_blank(), panel.grid.major.y = ggplot2::element_line(color = "grey95", linetype = "dashed"), legend.position = "none") +
      ggplot2::scale_color_brewer(palette = "Dark2")

    if (use_ggtext) p_go <- p_go + ggplot2::theme(strip.text.y.left = ggtext::element_markdown(angle = 0, face = "bold", size = 12))
    else p_go <- p_go + ggplot2::theme(strip.text.y.left = ggplot2::element_text(angle = 0, face = "bold", size = 12, color = onto_colors))

    plot_list[[ltype]] <- p_go
  }
  return(plot_list)
}

#' @title Generate Expression Heatmap and Connectivity Raincloud Plots
#' @return A named list of plot objects (Heatmap, Raincloud, Scatter).
#' @keywords internal
run_heatmap_and_connectivity <- function(target_genes, tpm_mat_raw, meta_raw, loop_stats_df, global_glist, heatmap_ntop, cor_method, current_proj_name, source_type, target_col = NULL, skip_heatmap = FALSE) {
  plots_list <- list()
  colnames(tpm_mat_raw) <- trimws(colnames(tpm_mat_raw))
  valid_s <- intersect(meta_raw$SampleID, colnames(tpm_mat_raw))
  if (length(valid_s) == 0) return(plots_list)

  curr_mat <- tpm_mat_raw[, valid_s, drop = FALSE]
  curr_meta <- meta_raw %>% dplyr::filter(SampleID %in% valid_s)

  if (!skip_heatmap) {
    expr_genes <- intersect(target_genes, rownames(curr_mat))
    if (length(expr_genes) >= 5) {
      mat_plot <- log2(curr_mat[expr_genes, , drop = FALSE] + 1)
      if (!is.null(heatmap_ntop) && nrow(mat_plot) > heatmap_ntop) {
        row_vars <- apply(mat_plot, 1, var, na.rm = TRUE)
        mat_plot <- mat_plot[head(names(sort(row_vars, decreasing = TRUE)), heatmap_ntop), , drop = FALSE]
      }
      mat_scaled <- t(scale(t(mat_plot)))
      mat_scaled[mat_scaled > 2] <- 2
      mat_scaled[mat_scaled < -2] <- -2
      mat_scaled[is.na(mat_scaled)] <- 0

      col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("#2BB2D1", "white", "#FF8181"))
      groups <- unique(curr_meta$Group)
      n_groups <- length(groups)
      cols <- if (n_groups <= 8) RColorBrewer::brewer.pal(max(3, n_groups), "Set2")[seq_len(n_groups)] else grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(n_groups)
      ha <- ComplexHeatmap::HeatmapAnnotation(Group = curr_meta$Group, col = list(Group = setNames(cols, groups)), simple_anno_size = unit(0.3, "cm"))
      plots_list$Heatmap <- ComplexHeatmap::Heatmap(mat_scaled, name = "Z-score", col = col_fun, cluster_columns = FALSE, show_row_names = (nrow(mat_scaled) <= 80), top_annotation = ha, border = TRUE, column_title = paste0("Expression Heatmap\n", current_proj_name), use_raster = FALSE)
    }
  }

  if (is.null(loop_stats_df)) return(plots_list)
  use_col <- NULL
  display_tag <- "Total Loops"

  if (!is.null(target_col)) {
    if (target_col %in% colnames(loop_stats_df)) {
      use_col <- target_col
      display_tag <- target_col
    } else return(plots_list)
  } else {
    if ("Total_Loops" %in% colnames(loop_stats_df)) use_col <- "Total_Loops"
    else if ("Loop_Degree" %in% colnames(loop_stats_df)) use_col <- "Loop_Degree"
    else if ("degree" %in% colnames(loop_stats_df)) use_col <- "degree"
  }

  if (is.null(use_col)) return(plots_list)

  gene_col_name <- colnames(loop_stats_df)[1]
  valid_targets <- intersect(target_genes, loop_stats_df[[gene_col_name]])
  if (length(valid_targets) < 5) return(plots_list)

  cols_to_extract <- unique(c(gene_col_name, use_col, intersect(colnames(loop_stats_df), c("Is_High_Connectivity_Gene", "Is_High_Distal_Connectivity_Gene", "High_Connectivity_Gene"))))
  stats_subset <- loop_stats_df[loop_stats_df[[gene_col_name]] %in% valid_targets, cols_to_extract]
  colnames(stats_subset)[which(colnames(stats_subset) == use_col)] <- "Degree"
  colnames(stats_subset)[1] <- "Gene"

  valid_expr_targets <- intersect(stats_subset$Gene, rownames(curr_mat))
  if (length(valid_expr_targets) < 5) return(plots_list)

  stats_subset <- stats_subset[stats_subset$Gene %in% valid_expr_targets, ]
  plot_df <- stats_subset %>%
    dplyr::mutate(Expression = as.numeric(log2(rowMeans(curr_mat[stats_subset$Gene, , drop = FALSE], na.rm = TRUE) + 1)), LFC = as.numeric(global_glist[stats_subset$Gene]), Log10Degree = log10(Degree)) %>%
    dplyr::filter(!is.na(Expression), !is.na(LFC), Degree >= 1)

  if (nrow(plot_df) < 5) return(plots_list)

  full_title_suffix <- paste0(if (source_type == "loops") "Looped (Specified Types)" else "Looped Targets", " | ", display_tag)

  plots_list$Scatter <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Log10Degree, y = Expression)) +
    ggpointdensity::geom_pointdensity(alpha = 0.6, size = 1.5) +
    viridis::scale_color_viridis(option = "D", name = "Density") +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, color = "black", se = TRUE, linewidth = 0.8) +
    ggpubr::stat_cor(method = cor_method, label.x.npc = "left", label.y.npc = "top", size = 4) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), panel.background = ggplot2::element_rect(color = "black", fill = "transparent"), legend.key = ggplot2::element_rect(fill = "transparent"), plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")) +
    ggplot2::labs(title = "Connectivity vs Expression (Scatter)", subtitle = paste0(full_title_suffix, "\nGenes: ", nrow(plot_df)), x = paste0("Log10 (", display_tag, ")"), y = "Log2(Mean Expression + 1)")

  if ("Is_High_Distal_Connectivity_Gene" %in% colnames(plot_df) && "Is_High_Connectivity_Gene" %in% colnames(plot_df)) {
    plot_df_rc <- plot_df %>%
      dplyr::mutate(
        high_distal = Is_High_Distal_Connectivity_Gene %in% c("Yes", "TRUE", TRUE, 1),
        high_total = Is_High_Connectivity_Gene %in% c("Yes", "TRUE", TRUE, 1)
      )
    plot_df_rc$Conn_Group <- ifelse(plot_df_rc$high_distal, "High Distal", ifelse(plot_df_rc$high_total, "High Total", "Others"))
    plot_df_rc$high_distal <- NULL
    plot_df_rc$high_total <- NULL
    plot_df_rc$Conn_Group <- factor(plot_df_rc$Conn_Group, levels = c("High Distal", "High Total", "Others"))
    custom_colors <- c("High Distal" = "#9BC985", "High Total" = "#ECB884", "Others" = "#82969D")
  } else {
    deg_thresh <- max(quantile(plot_df$Degree, 0.75, na.rm = TRUE), 2)
    plot_df_rc <- plot_df %>% dplyr::mutate(Conn_Group = factor(ifelse(Degree >= deg_thresh, "High Total", "Others"), levels = c("High Total", "Others")))
    custom_colors <- c("High Total" = "#ECB884", "Others" = "#82969D")
  }

  plot_df_rc <- plot_df_rc %>% dplyr::filter(!is.na(Conn_Group)) %>% droplevels()

  if (nlevels(plot_df_rc$Conn_Group) > 1) {
    clean_theme <- ggplot2::theme_classic() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14), plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 11, color = "black"), legend.position = "none", axis.text.x = ggplot2::element_text(angle = 20, hjust = 1, size = 10, color = "black"), axis.text.y = ggplot2::element_text(size = 10, color = "black"), axis.title.y = ggplot2::element_text(size = 12, face = "bold"), axis.line = ggplot2::element_line(color = "black", linewidth = 0.6), axis.ticks = ggplot2::element_line(color = "black"), panel.grid = ggplot2::element_blank())
    base_box <- function(y_var, y_lab) {
      ggplot2::ggplot(plot_df_rc, ggplot2::aes(fill = Conn_Group)) +
        ggplot2::geom_jitter(ggplot2::aes_string(x = "as.numeric(Conn_Group) - 0.12", y = y_var, color = "Conn_Group"), shape = 16, width = 0.03, height = 0, alpha = 0.6, size = 0.8) +
        ggplot2::stat_boxplot(ggplot2::aes_string(x = "as.numeric(Conn_Group)", y = y_var, color = "Conn_Group"), geom = "errorbar", width = 0.05, linewidth = 0.5) +
        ggplot2::geom_boxplot(ggplot2::aes_string(x = "as.numeric(Conn_Group)", y = y_var, color = "Conn_Group"), width = 0.12, notch = TRUE, outlier.shape = NA, alpha = 1, linewidth = 0.5) +
        ggplot2::stat_summary(ggplot2::aes_string(x = "as.numeric(Conn_Group)", y = y_var), fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.1, color = "black", linewidth = 0.4) +
        ggdist::stat_slab(ggplot2::aes_string(x = "as.numeric(Conn_Group) + 0.07", y = y_var, fill = "Conn_Group"), adjust = 0.5, width = 0.35, justification = 0, alpha = 0.3, color = NA) +
        ggdist::stat_slab(ggplot2::aes_string(x = "as.numeric(Conn_Group) + 0.07", y = y_var, color = "Conn_Group"), adjust = 0.5, width = 0.35, justification = 0, fill = NA, alpha = 0.5, linewidth = 0.4) +
        ggplot2::scale_x_continuous(breaks = seq_along(levels(plot_df_rc$Conn_Group)), labels = levels(plot_df_rc$Conn_Group)) +
        ggplot2::coord_cartesian(xlim = c(0.75, length(levels(plot_df_rc$Conn_Group)) + 0.6)) +
        ggplot2::scale_fill_manual(values = custom_colors) + ggplot2::scale_color_manual(values = custom_colors) +
        ggplot2::labs(title = "Regulation: High_connectivity vs Others", x = NULL, y = y_lab) + clean_theme
    }
    plots_list$Raincloud_LFC <- base_box("LFC", "Log2 Fold Change (LFC)") + ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.6)
    plots_list$Raincloud_Expr <- base_box("Expression", "Log2(Mean Expression + 1)")
  }
  return(plots_list)
}

.anchor_matches_targets <- function(gene_string, target_genes) {
  genes <- clean_gene_names(gene_string, ";")
  length(genes) > 0 && any(genes %in% target_genes)
}

#' @title Run Dual Motif Analysis for Loop Anchors
#' @param jaspar_db A JASPAR database object (e.g., \code{JASPAR2020::JASPAR2020} or \code{JASPAR2024::JASPAR2024}). Default: \code{JASPAR2020::JASPAR2020}.
#' @param jaspar_collection Character. JASPAR collection to query (e.g., \code{"CORE"}, \code{"CNE"}). Default: \code{"CORE"}.
#' @return A named list containing motif enrichment results and plot objects.
#' @keywords internal
run_distal_motif_analysis <- function(target_genes, loop_df, genome_id, pval_thresh,
  current_proj_name, top_n = 5, jaspar_db = JASPAR2020::JASPAR2020,
  jaspar_collection = "CORE") {
  bs_pkg <- switch(genome_id, "hg19" = "BSgenome.Hsapiens.UCSC.hg19", "hg38" = "BSgenome.Hsapiens.UCSC.hg38", "mm9" = "BSgenome.Mmusculus.UCSC.mm9", "mm10" = "BSgenome.Mmusculus.UCSC.mm10")
  species_id <- if (grepl("mm", genome_id)) 10090 else 9606
  genome_obj <- get0(bs_pkg, envir = asNamespace(bs_pkg))
  if (!is.data.frame(loop_df)) loop_df <- as.data.frame(loop_df)

  col_g1 <- intersect(c("anchor1_gene", "Anchor1_Gene", "gene_name_1", "Gene_Name_1", "Symbol_1", "nearest_gene_1", "gene1"), colnames(loop_df))[1]
  col_g2 <- intersect(c("anchor2_gene", "Anchor2_Gene", "gene_name_2", "Gene_Name_2", "Symbol_2", "nearest_gene_2", "gene2"), colnames(loop_df))[1]
  if (is.na(col_g1) || is.na(col_g2)) return(list())

  target_genes <- clean_gene_names(target_genes)
  a1_hits <- vapply(loop_df[[col_g1]], .anchor_matches_targets,
    target_genes = target_genes, FUN.VALUE = logical(1))
  a2_hits <- vapply(loop_df[[col_g2]], .anchor_matches_targets,
    target_genes = target_genes, FUN.VALUE = logical(1))
  target_indices <- which(a1_hits | a2_hits)
  if (length(target_indices) < 5) return(list())

  target_loops <- loop_df[target_indices, ]
  target_a1_hits <- a1_hits[target_indices]
  target_a2_hits <- a2_hits[target_indices]
  proximal_gr_list <- list()
  distal_gr_list <- list()

  for (i in seq_len(nrow(target_loops))) {
    r <- target_loops[i, ]
    if (target_a1_hits[i]) {
      proximal_gr_list[[length(proximal_gr_list) + 1]] <- GenomicRanges::GRanges(seqnames = r$chr1, ranges = IRanges::IRanges(start = r$start1, end = r$end1))
      distal_gr_list[[length(distal_gr_list) + 1]] <- GenomicRanges::GRanges(seqnames = r$chr2, ranges = IRanges::IRanges(start = r$start2, end = r$end2))
    }
    if (target_a2_hits[i]) {
      proximal_gr_list[[length(proximal_gr_list) + 1]] <- GenomicRanges::GRanges(seqnames = r$chr2, ranges = IRanges::IRanges(start = r$start2, end = r$end2))
      distal_gr_list[[length(distal_gr_list) + 1]] <- GenomicRanges::GRanges(seqnames = r$chr1, ranges = IRanges::IRanges(start = r$start1, end = r$end1))
    }
  }

  bg_indices <- setdiff(seq_len(nrow(loop_df)), target_indices)
  if (length(bg_indices) > 2000) bg_indices <- sample(bg_indices, 2000)
  bg_loops <- loop_df[bg_indices, ]
  bg_gr <- c(GenomicRanges::GRanges(seqnames = bg_loops$chr1, ranges = IRanges::IRanges(start = bg_loops$start1, end = bg_loops$end1)), GenomicRanges::GRanges(seqnames = bg_loops$chr2, ranges = IRanges::IRanges(start = bg_loops$start2, end = bg_loops$end2)))

  plots_list <- list()
  enrich_prox <- .calc_motif_enrichment(do.call(c, proximal_gr_list), bg_gr,
    genome_obj, pval_thresh, species_id, jaspar_db, jaspar_collection)
  res_prox <- .annotate_motif_families(enrich_prox, jaspar_db, jaspar_collection)
  plots_list$Proximal_Motif_Bar <- .plot_save_motif(res_prox, paste0(current_proj_name, "_Motif_Proximal"))
  plots_list$Proximal_Motif_Logos <- .plot_top_motif_logos(res_prox, top_n, jaspar_db)
  plots_list$Proximal_Motif_Rank <- .plot_motif_rank_scatter(res_prox, paste0(current_proj_name, "_Motif_Proximal"))

  enrich_dist <- .calc_motif_enrichment(do.call(c, distal_gr_list), bg_gr,
    genome_obj, pval_thresh, species_id, jaspar_db, jaspar_collection)
  res_dist <- .annotate_motif_families(enrich_dist, jaspar_db, jaspar_collection)
  plots_list$Distal_Motif_Bar <- .plot_save_motif(res_dist, paste0(current_proj_name, "_Motif_Distal"))
  plots_list$Distal_Motif_Logos <- .plot_top_motif_logos(res_dist, top_n, jaspar_db)
  plots_list$Distal_Motif_Rank <- .plot_motif_rank_scatter(res_dist, paste0(current_proj_name, "_Motif_Distal"))

  return(plots_list)
}

#' @title Plot Motif Rank Scatter
#' @return A \code{ggplot} object, or \code{NULL} if input is empty.
#' @keywords internal
.plot_motif_rank_scatter <- function(res_df, prefix, fdr_thresh = 0.05) {
  if (is.null(res_df) || nrow(res_df) == 0) return(NULL)
  if (!"Family" %in% colnames(res_df)) res_df$Family <- "Unknown"

  plot_df <- res_df %>%
    dplyr::mutate(FDR = ifelse(is.na(FDR), 1, FDR), LogFDR = -log10(FDR + 1e-300), Is_Sig = FDR < fdr_thresh, OddsRatio = pmax(0.8, pmin(1.6, as.numeric(OddsRatio)))) %>%
    dplyr::arrange(dplyr::desc(LogFDR)) %>%
    dplyr::mutate(Rank = dplyr::row_number())

  if (sum(plot_df$Is_Sig, na.rm = TRUE) == 0) plot_df$PlotFamily <- "Not Significant"
  else {
    fam_counts <- table(plot_df$Family[plot_df$Is_Sig & !plot_df$Family %in% c("Unknown", "", NA)])
    plot_df$PlotFamily <- ifelse(plot_df$Is_Sig, ifelse(plot_df$Family %in% names(sort(fam_counts, decreasing = TRUE))[seq_len(min(10, length(fam_counts)))], plot_df$Family, "Others"), "Not Significant")
  }

  plot_df$PlotFamily <- factor(plot_df$PlotFamily, levels = unique(c(setdiff(unique(plot_df$PlotFamily), c("Others", "Not Significant")), "Others", "Not Significant")))
  color_map <- setNames(c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A")[seq_len(length(setdiff(levels(plot_df$PlotFamily), c("Others", "Not Significant"))))], setdiff(levels(plot_df$PlotFamily), c("Others", "Not Significant")))
  color_map["Others"] <- "black"
  color_map["Not Significant"] <- "grey85"

  return(ggplot2::ggplot(plot_df, ggplot2::aes(x = Rank, y = LogFDR, size = OddsRatio, color = PlotFamily)) + ggplot2::geom_point(alpha = 0.8) + ggplot2::scale_color_manual(values = color_map, name = "TF Family (Top 10)") + ggplot2::scale_radius(name = "Odds Ratio", range = c(1, 7), breaks = c(0.8, 1.0, 1.2, 1.4, 1.6)) + ggplot2::geom_hline(yintercept = -log10(fdr_thresh), linetype = "dashed", color = "black", alpha = 0.5) + ggplot2::labs(title = paste0("Motif Enrichment Rank: ", basename(prefix)), x = "Rank", y = "-log10(FDR)") + ggplot2::theme_classic() + ggplot2::theme(legend.position = "right", legend.text = ggplot2::element_text(size = 8), legend.title = ggplot2::element_text(size = 9, face = "bold")))
}

#' @title Calculate Motif Enrichment via Fisher's Exact Test
#' @return A data frame of enrichment results, or \code{NULL} if input is empty.
#' @keywords internal
.calc_motif_enrichment <- function(fg_gr, bg_gr, genome_obj, pval_thresh, species_id,
  jaspar_db = JASPAR2020::JASPAR2020, jaspar_collection = "CORE") {
  fg_gr <- GenomicRanges::resize(fg_gr[GenomicRanges::start(fg_gr) > 0], width = 500, fix = "center")
  bg_gr <- GenomicRanges::resize(bg_gr[GenomicRanges::start(bg_gr) > 0], width = 500, fix = "center")
  if (length(fg_gr) == 0 || length(bg_gr) == 0) return(NULL)

  fg_seq <- BSgenome::getSeq(genome_obj, fg_gr)
  bg_seq <- BSgenome::getSeq(genome_obj, bg_gr)
  pfm_list <- TFBSTools::getMatrixSet(jaspar_db, list(species = species_id, collection = jaspar_collection))
  if (length(pfm_list) == 0) pfm_list <- TFBSTools::getMatrixSet(jaspar_db, list(collection = jaspar_collection))
  if (length(pfm_list) == 0) return(NULL)

  fg_counts <- colSums(as.matrix(motifmatchr::motifMatches(motifmatchr::matchMotifs(pfm_list, fg_seq, out = "matches", p.cutoff = pval_thresh))))
  bg_counts <- colSums(as.matrix(motifmatchr::motifMatches(motifmatchr::matchMotifs(pfm_list, bg_seq, out = "matches", p.cutoff = pval_thresh))))

  results_list <- lapply(union(names(fg_counts), names(bg_counts)), function(m) {
    a <- if (m %in% names(fg_counts)) fg_counts[[m]] else 0
    b <- if (m %in% names(bg_counts)) bg_counts[[m]] else 0
    if (a > 0) {
      ft <- fisher.test(matrix(c(a, b, length(fg_seq) - a, length(bg_seq) - b), nrow = 2), alternative = "greater")
      data.frame(MotifID = m, MotifName = TFBSTools::name(pfm_list[[m]]), Pvalue = ft$p.value, OddsRatio = ft$estimate, FG_Hits = a, FG_Total = length(fg_seq), BG_Hits = b, BG_Total = length(bg_seq))
    } else NULL
  })

  res_df <- do.call(rbind, Filter(Negate(is.null), results_list))
  if (!is.null(res_df) && nrow(res_df) > 0) {
    res_df$FDR <- p.adjust(res_df$Pvalue, method = "BH")
    res_df <- res_df[order(res_df$Pvalue), ]
  }
  return(res_df)
}

#' @title Plot and Save Motif Results (Barplot)
#' @return A \code{ggplot} barplot, or \code{NULL} if input is empty.
#' @keywords internal
.plot_save_motif <- function(res_df, prefix) {
  if (is.null(res_df) || nrow(res_df) == 0) return(NULL)
  top_df <- head(res_df, 15)
  top_df$MotifLabel <- factor(paste0(top_df$MotifName, " (", top_df$MotifID, ")"), levels = rev(paste0(top_df$MotifName, " (", top_df$MotifID, ")")))
  return(ggplot2::ggplot(top_df, ggplot2::aes(x = -log10(.data$Pvalue), y = MotifLabel)) + ggplot2::geom_col(fill = "#E7298A", width = 0.7) + ggplot2::labs(title = paste0("Motif Enrichment: ", basename(prefix)), x = "-log10(P-value)", y = NULL) + ggplot2::theme_classic())
}

#' @title Plot Top Motif Sequence Logos
#' @return A list of sequence logo plots, or \code{NULL} if input is empty.
#' @keywords internal
.plot_top_motif_logos <- function(res_df, top_n,
  jaspar_db = JASPAR2020::JASPAR2020) {
  if (is.null(res_df) || nrow(res_df) == 0) return(NULL)
  top_df <- head(res_df[order(res_df$Pvalue), ], top_n)
  pfm_list <- TFBSTools::getMatrixSet(jaspar_db, opts = list(ID = top_df$MotifID))
  if (length(pfm_list) == 0) return(NULL)

  plot_list <- list()
  for (i in seq_along(top_df$MotifID)) if (top_df$MotifID[i] %in% names(pfm_list)) plot_list[[paste0(top_df$MotifName[i], " (", top_df$MotifID[i], ")")]] <- TFBSTools::Matrix(pfm_list[[top_df$MotifID[i]]])
  if (length(plot_list) == 0) return(NULL)

  return(ggseqlogo::ggseqlogo(plot_list, ncol = 1) + ggplot2::theme_classic() + ggplot2::theme(axis.text.x = ggplot2::element_blank(), strip.text = ggplot2::element_text(size = 10, face = "bold", hjust = 0), strip.background = ggplot2::element_rect(fill = "grey95", color = NA)) + ggplot2::labs(y = "Bits", title = paste0("Top ", top_n, " Enriched Motifs (SeqLogo)")))
}

#' @title Annotate Motif Families
#' @return A data frame with added \code{Family} column, sorted by P-value.
#' @keywords internal
.annotate_motif_families <- function(res_df,
  jaspar_db = JASPAR2020::JASPAR2020, jaspar_collection = "CORE") {
  if (is.null(res_df) || nrow(res_df) == 0) return(res_df)
  meta_df <- do.call(rbind, lapply(TFBSTools::getMatrixSet(jaspar_db, list(collection = jaspar_collection)), function(x) data.frame(MotifID = TFBSTools::ID(x), Family = paste(if (is.null(TFBSTools::tags(x)$family)) "Unknown" else TFBSTools::tags(x)$family, collapse = "; "), stringsAsFactors = FALSE)))
  res_df <- merge(res_df[, !colnames(res_df) %in% c("Family", "Class")], meta_df, by = "MotifID", all.x = TRUE)
  return(res_df[order(res_df$Pvalue), c(c("MotifID", "MotifName", "Family", "Pvalue", "FDR", "OddsRatio", "FG_Hits", "FG_Total", "BG_Hits", "BG_Total"), setdiff(colnames(res_df), c("MotifID", "MotifName", "Family", "Pvalue", "FDR", "OddsRatio", "FG_Hits", "FG_Total", "BG_Hits", "BG_Total")))])
}


#' Render a Publication-Ready 3D Annotation Report
#'
#' One-click parameterised R Markdown report that executes the full looplook
#' pipeline (annotation → refinement → profiling) and renders an
#' interpretation-ready HTML document suitable for sharing with collaborators.
#'
#' @details
#' The profiling stage uses R's random number generator for GSEA gene-set
#' down-sampling (\code{gsea_nSample}) and motif background sampling. Call
#' \code{set.seed()} before \code{looplook_report()} for fully reproducible
#' results.
#'
#' @param bedpe_file Path to a BEDPE file of chromatin loops.
#' @param target_bed Optional path to a BED file of genomic features.
#' @param expr_matrix_file Optional path to a normalised expression matrix.
#' @param sample_columns Sample columns in the expression matrix to average.
#' @param species Genome assembly (\code{"hg38"}, \code{"hg19"}, \code{"mm10"}, \code{"mm9"}).
#' @param project_name Character. Project prefix for the report title.
#' @param out_dir Output directory. Created if missing.
#' @param threshold Numeric. Expression threshold for active gene classification.
#' @param reclassify_by_expression Logical. Reclassify silent promoters as eP/eG.
#' @param run_go Logical. Run GO enrichment (requires clusterProfiler).
#' @param diff_file Optional differential expression result file.
#' @param lfc_col Column name for log2 fold change in \code{diff_file}.
#' @param metadata_file Optional sample metadata file.
#' @param precomputed_res Optional. Either a \code{.RData} file path or an
#'   in-memory list object returned by \code{annotate_peaks_and_loops}.
#'   When provided, annotation is skipped and refinement starts from this object.
#' @param unit_type Character. Expression unit label for plot annotations. Default \code{"TPM"}.
#' @param tss_region Numeric vector of length 2. TSS flanking region in bp. Default \code{c(-2000, 2000)}.
#' @param neighbor_hop Integer. Neighbourhood extension for fragment overlap. Default \code{0}.
#' @param hub_percentile Numeric. Quantile threshold for hub classification. Default \code{0.95}.
#' @param color_palette Character. RColorBrewer qualitative palette name. Default \code{"Set2"}.
#' @param target_mapping_mode Character. Mapping strategy for target genes. Default \code{"all"}.
#' @param include_Filled Logical. Include comprehensively merged gene assignments. Default \code{TRUE}.
#' @param use_nearest_gene Logical. Bypass 3D loop-based assignment and use linear proximity. Default \code{FALSE}.
#' @param target_source Character vector. Source of target genes to profile. Default \code{"targets"}.
#' @param loop_types Character vector. Loop types to include in profiling. Default \code{c("E-P", "P-P")}.
#' @param stat_test Character. Statistical test for violin comparisons. Default \code{"wilcox.test"}.
#' @param run_ppi Logical. Run protein-protein interaction network analysis. Default \code{FALSE}.
#' @param run_motif Logical. Run transcription factor motif analysis. Default \code{FALSE}.
#' @param genome_id Character. Reference genome for motif scanning. Defaults to \code{species}.
#' @param motif_p_thresh Numeric. P-value threshold for motif enrichment. Default \code{1e-4}.
#' @param motif_ntop Numeric. Number of top motifs to display. Default \code{5}.
#' @param ppi_score Numeric. Minimum STRING combined score. Default \code{400}.
#' @param ppi_nSample Numeric. Maximum genes to include in PPI. Default \code{400}.
#' @param heatmap_nSample Numeric. Maximum genes in expression heatmap. Default \code{99999}.
#' @param gsea_nSample Numeric. Maximum genes sampled for GSEA. Default \code{99999}.
#' @param cnet_nSample Numeric. Number of GO terms in cnetplot. Default \code{50}.
#' @param karyo_bin_size Numeric. Bin size for karyotype heatmaps. Default \code{100000}.
#' @param output_file Character. Output HTML file name. \code{NULL} derives from \code{project_name}.
#' @param quiet Logical. Suppress rendering output. Default \code{FALSE}.
#' @param ... Additional arguments passed to \code{rmarkdown::render}.
#'
#' @return The path to the generated HTML report (invisibly).
#'
#' @export
#'
#' @examples
#' \donttest{
#' looplook_report(
#'   bedpe_file   = system.file("extdata", "example_loops_1.bedpe", package = "looplook"),
#'   target_bed   = system.file("extdata", "example_peaks.bed", package = "looplook"),
#'   project_name = "My Study",
#'   out_dir      = tempdir()
#' )
#' }
looplook_report <- function(
  bedpe_file,
  target_bed = NULL,
  expr_matrix_file = NULL,
  sample_columns = NULL,
  species = "hg38",
  project_name = "looplook Analysis",
  out_dir = "looplook_results",
  threshold = 1.0,
  unit_type = "TPM",
  reclassify_by_expression = TRUE,
  tss_region = c(-2000, 2000),
  neighbor_hop = 0,
  hub_percentile = 0.95,
  color_palette = "Set2",
  target_mapping_mode = "all",
  include_Filled = TRUE,
  use_nearest_gene = FALSE,
  target_source = "targets",
  loop_types = c("E-P", "P-P"),
  stat_test = "wilcox.test",
  run_go = TRUE,
  run_ppi = FALSE,
  run_motif = FALSE,
  genome_id = species,
  motif_p_thresh = 1e-4,
  motif_ntop = 5,
  ppi_score = 400,
  ppi_nSample = 400,
  heatmap_nSample = 99999,
  gsea_nSample = 99999,
  cnet_nSample = 50,
  karyo_bin_size = 1e5,
  diff_file = NULL,
  lfc_col = "log2FoldChange",
  metadata_file = NULL,
  precomputed_res = NULL,
  output_file = NULL,
  quiet = FALSE,
  ...
) {
  # Locate template
  template <- system.file("rmarkdown", "templates", "looplook-report",
    "skeleton", "skeleton.Rmd", package = "looplook")
  if (!nzchar(template)) stop("Report template not found. Reinstall looplook.")

  # Create output directory
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  }

  # Prepare output filename
  if (is.null(output_file)) {
    output_file <- paste0(gsub("[[:space:]]+", "_", project_name), "_Report.html")
  }

  if (!quiet) message(">> Generating looplook report: ", output_file)

  # Render report
  out_path <- rmarkdown::render(
    input = template,
    params = list(
      bedpe_file = bedpe_file,
      target_bed = target_bed,
      expr_matrix_file = expr_matrix_file,
      sample_columns = sample_columns,
      tss_region = tss_region,
      species = species,
      project_name = project_name,
      out_dir = out_dir,
      threshold = threshold,
      unit_type = unit_type,
      reclassify_by_expression = reclassify_by_expression,
      neighbor_hop = neighbor_hop,
      hub_percentile = hub_percentile,
      color_palette = color_palette,
      target_mapping_mode = target_mapping_mode,
      include_Filled = include_Filled,
      use_nearest_gene = use_nearest_gene,
      target_source = target_source,
      loop_types = loop_types,
      stat_test = stat_test,
      run_go = run_go,
      run_ppi = run_ppi,
      run_motif = run_motif,
      genome_id = genome_id,
      motif_p_thresh = motif_p_thresh,
      motif_ntop = motif_ntop,
      ppi_score = ppi_score,
      ppi_nSample = ppi_nSample,
      heatmap_nSample = heatmap_nSample,
      gsea_nSample = gsea_nSample,
      cnet_nSample = cnet_nSample,
      karyo_bin_size = karyo_bin_size,
      diff_file = diff_file,
      lfc_col = lfc_col,
      metadata_file = metadata_file,
      precomputed_res = precomputed_res
    ),
    output_dir = out_dir,
    output_file = output_file,
    quiet = quiet,
    ...
  )

  if (!quiet) message(">> Report saved to: ", out_path)
  return(invisible(out_path))
}
