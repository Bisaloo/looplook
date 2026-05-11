#' Internal: Prepare Track Data for IGV-Style Plot
#'
#' Reads BEDPE, optional BED, and gene annotation data; computes bezier curves,
#' alpha colours, and exon features for the integrative track plot.
#'
#' @param bedpe_file Path to BEDPE file.
#' @param target_bed Optional path to BED file.
#' @param chr Chromosome name.
#' @param from Start coordinate.
#' @param to End coordinate.
#' @param species Genome assembly.
#' @param max_levels Maximum visible arc stacking levels. Overlapping loops are
#'   assigned to distinct vertical bands up to this cap; denser regions are
#'   compressed into the available height while preserving relative layering.
#' @param base_anchor_height Anchor rectangle height.
#' @param loop_color Default arc colour.
#' @param anchor_color Default anchor colour.
#' @param score_to_alpha Map score to transparency.
#' @param min_score Optional score floor.
#' @return A named list of all data frames and plot parameters.
#' @keywords internal
#' @noRd
prepare_track_data <- function(
  bedpe_file, target_bed, chr, from, to, species,
  max_levels, base_anchor_height, loop_color, anchor_color, score_to_alpha,
  min_score
) {
  .scale_group_positions <- function(values, groups) {
    out <- numeric(length(values))
    idx_split <- split(seq_along(values), groups)
    for (idx in idx_split) {
      if (length(idx) == 1) {
        out[idx] <- 0.5
      } else {
        out[idx] <- scales::rescale(rank(values[idx], ties.method = "first"), to = c(0, 1))
      }
    }
    out
  }

  .assign_arc_levels <- function(left, right, max_levels) {
    n <- length(left)
    if (n == 0) {
      return(data.frame(
        raw_level = integer(0),
        display_level = integer(0),
        band_pos = numeric(0),
        stringsAsFactors = FALSE
      ))
    }

    left <- as.numeric(left)
    right <- as.numeric(right)
    span <- pmax(right - left, 1)
    max_levels <- max(1L, as.integer(max_levels[1]))

    raw_level <- integer(n)
    level_right_edge <- numeric(0)
    ord <- order(left, right)
    for (idx in ord) {
      reusable <- which(level_right_edge <= left[idx])
      if (length(reusable) == 0) {
        level_right_edge <- c(level_right_edge, right[idx])
        raw_level[idx] <- length(level_right_edge)
      } else {
        use_level <- reusable[1]
        raw_level[idx] <- use_level
        level_right_edge[use_level] <- right[idx]
      }
    }

    raw_max <- max(raw_level)
    compressed <- raw_max > max_levels
    display_level <- if (compressed) {
      ceiling(raw_level / raw_max * max_levels)
    } else {
      raw_level
    }
    band_values <- if (compressed) raw_level else span
    band_pos <- .scale_group_positions(band_values, display_level)

    data.frame(
      raw_level = raw_level,
      display_level = display_level,
      band_pos = band_pos,
      stringsAsFactors = FALSE
    )
  }

  # Species config
  species_map <- list(
    hg38 = c("TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db"),
    hg19 = c("TxDb.Hsapiens.UCSC.hg19.knownGene", "org.Hs.eg.db"),
    mm10 = c("TxDb.Mmusculus.UCSC.mm10.knownGene", "org.Mm.eg.db"),
    mm9  = c("TxDb.Mmusculus.UCSC.mm9.knownGene", "org.Mm.eg.db")
  )
  cfg <- species_map[[species]]
  if (is.null(cfg)) stop("Species not supported: ", species)
  txdb_pkg <- cfg[1]
  org_db_pkg <- cfg[2]

  if (!requireNamespace(txdb_pkg, quietly = TRUE)) stop("Please install ", txdb_pkg)
  if (!requireNamespace(org_db_pkg, quietly = TRUE)) stop("Please install ", org_db_pkg)
  txdb <- utils::getFromNamespace(txdb_pkg, txdb_pkg)

  # Read BEDPE (0-based) and convert to 1-based for internal use
  loops_raw <- as.data.frame(data.table::fread(bedpe_file, header = FALSE))
  colnames(loops_raw)[seq_len(6)] <- c(
    "chr1", "start1", "end1", "chr2",
    "start2", "end2"
  )
  loops_raw$start1 <- loops_raw$start1 + 1
  loops_raw$start2 <- loops_raw$start2 + 1

  has_score <- FALSE
  if (ncol(loops_raw) >= 7) {
    is_numeric_col <- function(x) {
      nums <- as.numeric(x)
      sum(!is.na(nums)) / length(nums) > 0.9 && !all(is.na(nums))
    }
    if (ncol(loops_raw) >= 8 && is_numeric_col(loops_raw[[8]])) {
      colnames(loops_raw)[8] <- "score"
      has_score <- TRUE
    } else if (is_numeric_col(loops_raw[[7]])) {
      colnames(loops_raw)[7] <- "score"
      has_score <- TRUE
    }
  }

  loops <- loops_raw %>%
    dplyr::mutate(
      dplyr::across(c(start1, end1, start2, end2), as.numeric),
      chr1 = as.character(chr1), chr2 = as.character(chr2)
    )
  if (has_score) {
    loops$score <- as.numeric(loops$score)
    if (!is.null(min_score)) {
      loops <- loops %>% dplyr::filter(score >= min_score)
      if (nrow(loops) == 0) stop("All loops filtered by min_score.")
    }
  }

  if (is.null(chr)) {
    all_chr <- c(loops$chr1, loops$chr2)
    chr <- names(sort(table(all_chr), decreasing = TRUE))[1]
  }
  loops_chr <- loops %>% dplyr::filter(chr1 == chr2 & chr1 == chr)
  if (is.null(from)) {
    from <- if (nrow(loops_chr) > 0) {
      min(c(loops_chr$start1, loops_chr$start2))
    } else {
      warning("No loops found on ", chr, "; using coordinate 1")
      1
    }
  }
  if (is.null(to)) {
    to <- if (nrow(loops_chr) > 0) {
      max(c(loops_chr$end1, loops_chr$end2))
    } else {
      warning("No loops found on ", chr, "; using coordinate 1")
      1
    }
  }

  loops_view <- loops_chr %>%
    dplyr::filter((end1 >= from & start1 <= to) |
      (end2 >= from & start2 <= to))

  text_indent <- from + ((to - from) * 0.005)
  anchor_ymax <- base_anchor_height

  bez_df <- data.frame()
  anchors <- data.frame()
  plot_ymax <- anchor_ymax + 0.5

  if (nrow(loops_view) > 0) {
    anchors <- dplyr::bind_rows(
      loops_view %>% dplyr::transmute(
        chr = chr1, start = start1, end = end1,
        loop_i = dplyr::row_number(), score = if (has_score) score else 1
      ),
      loops_view %>% dplyr::transmute(
        chr = chr2, start = start2, end = end2,
        loop_i = dplyr::row_number(), score = if (has_score) score else 1
      )
    ) %>% dplyr::mutate(ymin = 0, ymax = anchor_ymax)

    loops_calc <- loops_view %>%
      dplyr::mutate(
        mid1 = (start1 + end1) / 2, mid2 = (start2 + end2) / 2,
        left_mid = pmin(mid1, mid2), right_mid = pmax(mid1, mid2)
      ) %>%
      dplyr::mutate(
        span = pmax(.data$right_mid - .data$left_mid, 1),
        loop_i = dplyr::row_number(), center = (mid1 + mid2) / 2
      )
    level_df <- .assign_arc_levels(loops_calc$left_mid, loops_calc$right_mid, max_levels)
    loops_calc <- cbind(loops_calc, level_df)

    peak_base <- anchor_ymax + 0.08
    level_band <- 0.11
    band_height <- 0.045
    loops_calc$peak <- peak_base +
      (loops_calc$display_level - 1) * level_band +
      loops_calc$band_pos * band_height
    loops_calc <- loops_calc[order(loops_calc$peak, loops_calc$span), , drop = FALSE]

    bez_list <- lapply(seq_len(nrow(loops_calc)), function(i) {
      d <- loops_calc[i, ]
      data.frame(
        loop_i = d$loop_i,
        x = c(d$mid1, d$center, d$mid2),
        y = c(anchor_ymax, d$peak, anchor_ymax),
        arc_level = d$display_level,
        score = if (has_score) d$score else 1, stringsAsFactors = FALSE
      )
    })
    bez_df <- do.call(rbind, bez_list)
    plot_ymax <- max(bez_df$y) + 0.08

    calc_alpha_color <- function(scores, base_col, use_alpha) {
      if (anyNA(scores)) scores[is.na(scores)] <- min(scores, na.rm = TRUE)
      alphas <- if (use_alpha && max(scores) != min(scores)) {
        scales::rescale(scores, to = c(0.1, 1.0))
      } else {
        rep(0.8, length(scores))
      }
      rgb_val <- col2rgb(base_col)
      rgb(rgb_val[1], rgb_val[2], rgb_val[3],
        alpha = alphas * 255,
        maxColorValue = 255
      )
    }
    do_map <- has_score && score_to_alpha
    bez_df$final_color <- calc_alpha_color(bez_df$score, loop_color, do_map)
    anchors$final_fill <- calc_alpha_color(anchors$score, anchor_color, do_map)
  }

  # Overlap BED
  overlap_df_plot <- NULL
  if (!is.null(target_bed)) {
    ob <- as.data.frame(data.table::fread(target_bed, header = FALSE, select = seq_len(3)))
    colnames(ob) <- c("chr", "start", "end")
    ob$start <- ob$start + 1 # BED is 0-based; convert to 1-based
    overlap_df_plot <- ob %>%
      dplyr::filter(chr == !!chr, !(end < from | start > to)) %>%
      dplyr::mutate(ymin = -0.15, ymax = -0.10)
    if (nrow(overlap_df_plot) == 0) overlap_df_plot <- NULL
  }

  # Gene annotation
  genes_gr <- .with_known_upstream_noise_suppressed(GenomicFeatures::genes(txdb))
  genes_df <- as.data.frame(genes_gr) %>%
    dplyr::filter(as.character(seqnames) == chr, end > from, start < to)

  feature_df <- data.frame()
  if (nrow(genes_df) > 0) {
    try(
      {
        symbol_map <- .map_txdb_gene_ids(
          gene_ids = unique(as.character(genes_df$gene_id)),
          org_db = org_db_pkg,
          columns = "SYMBOL",
          context = "prepare_genomic_plot_data gene track",
          warn = TRUE
        )
        symbol_map <- symbol_map[!duplicated(symbol_map$gene_id), ]
        genes_df <- dplyr::left_join(genes_df, symbol_map, by = "gene_id")
      },
      silent = TRUE
    )
    if (!"SYMBOL" %in% colnames(genes_df)) {
      genes_df$SYMBOL <- genes_df$gene_id
    }
    genes_df <- genes_df %>%
      dplyr::mutate(
        final_label = ifelse(is.na(SYMBOL), gene_id, SYMBOL),
        label_x = pmax(from, pmin(to, ifelse(strand == "+", start, end)))
      )
    genes_df$gene_level <- IRanges::disjointBins(
      IRanges::IRanges(genes_df$start, genes_df$end)
    )

    tx_keytype <- if ("TXNAME" %in% AnnotationDbi::keytypes(txdb)) {
      "TXNAME"
    } else {
      "TXID"
    }
    tx2gene <- AnnotationDbi::select(txdb,
      keys = AnnotationDbi::keys(txdb, tx_keytype),
      columns = "GENEID", keytype = tx_keytype
    )
    colnames(tx2gene) <- c("tx_id", "gene_id")

    exons_list <- GenomicFeatures::exonsBy(txdb, "tx", use.names = TRUE)
    exons_gr <- unlist(exons_list)
    names(exons_gr) <- NULL
    exons_flat <- as.data.frame(exons_gr)
    exons_flat$tx_id <- rep(names(exons_list),
      times = S4Vectors::elementNROWS(exons_list)
    )
    exons_joined <- exons_flat %>%
      dplyr::left_join(tx2gene, by = "tx_id") %>%
      dplyr::filter(gene_id %in% genes_df$gene_id) %>%
      dplyr::filter(start < to & end > from)

    if (nrow(exons_joined) > 0) {
      longest_tx <- exons_joined %>%
        dplyr::group_by(gene_id, tx_id) %>%
        dplyr::summarise(len = sum(width), .groups = "drop") %>%
        dplyr::arrange(desc(len)) %>%
        dplyr::group_by(gene_id) %>%
        dplyr::slice(1)
      feature_df <- exons_joined %>%
        dplyr::filter(tx_id %in% longest_tx$tx_id) %>%
        dplyr::left_join(genes_df[, c("gene_id", "gene_level")],
          by = "gene_id"
        )
    }
  }

  list(
    chr = chr, from = from, to = to, text_indent = text_indent,
    bez_df = bez_df, anchors = anchors, overlap_df_plot = overlap_df_plot,
    genes_df = genes_df, feature_df = feature_df,
    plot_ymax = plot_ymax
  )
}

#' Integrative visualization of 3D chromatin loops and genomic features
#'
#' Generates an integrative genomic track plot displaying chromatin loops as arcs,
#' loop anchors as rectangles, optional overlapping features (e.g., ChIP-seq peaks),
#' and annotated genes. Loop arcs can be colored or sized by interaction score
#' (7th column in BEDPE).
#'
#' @param bedpe_file Character. Path to a BEDPE file (at least 6 columns; 7th
#'   column used as score if present).
#' @param target_bed Optional character. Path to a BED file (e.g., peaks) to
#'   overlay below the loop track.
#' @param chr Character. Chromosome name (e.g., "chr8"). If NULL, inferred from
#'   the most frequent chromosome in the BEDPE.
#' @param from Numeric. Start coordinate of the region to plot.
#' @param to Numeric. End coordinate of the region to plot.
#' @param species Character. Genome assembly: "hg38", "hg19", "mm10", or "mm9".
#' @param max_levels Integer. Maximum number of visible vertical levels for loop
#'   arc stacking (default: 10). Overlapping loops are separated into stacked
#'   bands up to this limit; denser regions are compressed into the available
#'   height while preserving relative layering.
#' @param base_anchor_height Numeric. Height of anchor rectangles (default: 0.05).
#' @param loop_color Character. Default color for arcs when no score is provided
#'   (default: "#5D6D7E").
#' @param anchor_color Character. Color for loop anchor rectangles
#'   (default: "#3498DB").
#' @param overlap_color Character. Color for overlap track (default: "#02ABB4").
#' @param exon_color Character. Gene exon fill color (default: "#2C3E50").
#' @param intron_color Character. Gene intron line color (default: "black").
#' @param score_to_alpha Logical. Whether to map interaction scores to arc
#'   transparency.
#' @param min_score Optional numeric. Floor value for score mapping.
#' @param save_file Character. Optional path to save the plot via
#'   \code{ggplot2::ggsave()}. When set, the plot is written to this file and
#'   the same \code{ggplot} object is still returned.
#' @return A \code{ggplot} object. If \code{save_file} is provided, the plot is
#'   also written to disk via \code{ggplot2::ggsave()}.
#' @importFrom dplyr %>%
#' @importFrom ggplot2 ggplot geom_rect geom_segment annotate coord_cartesian
#'   scale_x_continuous theme_classic labs ggsave arrow unit theme element_blank
#'   element_rect element_text margin
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggforce geom_bezier
#' @importFrom scales rescale comma
#' @importFrom GenomicFeatures genes exonsBy
#' @importFrom AnnotationDbi keys keytypes
#' @export
#' @examples
#' if (requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE) &&
#'   requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
#'   bedpe_path <- tempfile(fileext = ".bedpe")
#'   writeLines(
#'     "chr1\t11890000\t11890500\tchr1\t11905000\t11905500",
#'     bedpe_path
#'   )
#'   p <- plot_peaks_interactions(
#'     bedpe_file = bedpe_path,
#'     chr = "chr1",
#'     from = 11884299,
#'     to = 12106581,
#'     species = "hg38"
#'   )
#'   print(p)
#' }
plot_peaks_interactions <- function(
  bedpe_file,
  target_bed = NULL,
  chr = NULL,
  from = NULL,
  to = NULL,
  species = "hg38",
  max_levels = 10,
  base_anchor_height = 0.05,
  loop_color = "#5D6D7E",
  anchor_color = "#3498DB",
  overlap_color = "#02ABB4",
  exon_color = "#2C3E50",
  intron_color = "black",
  score_to_alpha = TRUE,
  min_score = NULL,
  save_file = NULL
) {
  d <- prepare_track_data(
    bedpe_file, target_bed, chr, from, to, species,
    max_levels, base_anchor_height, loop_color, anchor_color,
    score_to_alpha, min_score
  )

  chr <- d$chr
  from <- d$from
  to <- d$to

  gene_start_y <- -0.25
  row_height <- 0.12
  plot_ymin <- gene_start_y

  if (nrow(d$genes_df) > 0) {
    d$genes_df <- d$genes_df %>%
      dplyr::mutate(y_mid = gene_start_y - (gene_level - 1) * row_height)
    plot_ymin <- min(d$genes_df$y_mid) - 0.2
    if (nrow(d$feature_df) > 0) {
      d$feature_df <- d$feature_df %>%
        dplyr::mutate(
          y_mid = gene_start_y - (gene_level - 1) * row_height,
          ymin = y_mid - 0.025, ymax = y_mid + 0.025
        )
    }
  }

  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = -0.04, linetype = "dashed", color = "grey85", linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = -0.18, linetype = "dashed", color = "grey85", linewidth = 0.5) +
    ggplot2::annotate("text",
      x = d$text_indent, y = d$plot_ymax - 0.01,
      label = "loop track", hjust = 0, vjust = 1, size = 4, fontface = "bold", color = "black"
    )

  if (!is.null(d$overlap_df_plot)) {
    p <- p + ggplot2::annotate("text",
      x = d$text_indent, y = -0.08,
      label = "target track", hjust = 0, vjust = 0, size = 4, fontface = "bold", color = "black"
    )
  }
  p <- p + ggplot2::annotate("text",
    x = d$text_indent, y = -0.21,
    label = "gene track", hjust = 0, vjust = 0, size = 4, fontface = "bold", color = "black"
  )

  if (!is.null(d$overlap_df_plot)) {
    p <- p + ggplot2::geom_rect(
      data = d$overlap_df_plot,
      ggplot2::aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax),
      fill = overlap_color, alpha = 1
    )
  }

  if (nrow(d$genes_df) > 0) {
    p <- p +
      ggplot2::geom_segment(
        data = d$genes_df,
        ggplot2::aes(
          x = pmax(start, from), xend = pmin(end, to),
          y = y_mid, yend = y_mid
        ), color = intron_color, linewidth = 0.5
      ) +
      ggplot2::geom_segment(
        data = d$genes_df,
        ggplot2::aes(
          x = ifelse(strand == "+", pmin(end, to), pmax(start, from)),
          xend = ifelse(strand == "+", pmin(end, to), pmax(start, from)),
          y = y_mid, yend = y_mid
        ),
        arrow = ggplot2::arrow(length = ggplot2::unit(0.15, "cm"), type = "open"),
        color = intron_color, linewidth = 0.5
      )

    if (nrow(d$feature_df) > 0) {
      p <- p + ggplot2::geom_rect(
        data = d$feature_df,
        ggplot2::aes(
          xmin = pmax(start, from), xmax = pmin(end, to),
          ymin = ymin, ymax = ymax
        ), fill = exon_color, color = NA
      )
    } else {
      p <- p + ggplot2::geom_rect(
        data = d$genes_df,
        ggplot2::aes(
          xmin = pmax(start, from), xmax = pmin(end, to),
          ymin = y_mid - 0.025, ymax = y_mid + 0.025
        ), fill = exon_color
      )
    }
    p <- p + ggrepel::geom_text_repel(
      data = d$genes_df,
      ggplot2::aes(x = label_x, y = y_mid, label = final_label),
      nudge_y = -0.05, direction = "x", force = 1, size = 3,
      segment.size = 0.3, segment.color = "grey60",
      segment.linetype = "dashed", min.segment.length = 0
    )
  }

  if (nrow(d$bez_df) > 0) {
    p <- p +
      ggplot2::geom_rect(
        data = d$anchors,
        ggplot2::aes(
          xmin = start, xmax = end, ymin = ymin, ymax = ymax,
          fill = final_fill
        ), color = NA
      ) +
      ggplot2::scale_fill_identity() +
      ggforce::geom_bezier(
        data = d$bez_df,
        ggplot2::aes(x = x, y = y, group = loop_i, color = final_color),
        linewidth = 0.6
      ) +
      ggplot2::scale_color_identity()
  }

  p <- p +
    ggplot2::coord_cartesian(
      xlim = c(from, to), ylim = c(plot_ymin, d$plot_ymax),
      expand = FALSE
    ) +
    ggplot2::scale_x_continuous(labels = scales::comma) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.line.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
      plot.title = ggplot2::element_text(hjust = 0.5)
    ) +
    ggplot2::labs(
      x = paste0(chr, ":", from, "-", to),
      title = "Loops Integrative View"
    )

  if (!is.null(save_file)) {
    ggplot2::ggsave(save_file, plot = p, width = 10, height = 5)
  }
  return(p)
}

#' Draw Simplified Flower Plot for Core vs. Unique Genes
#'
#' Creates a circular "flower" diagram where each petal represents a gene set,
#' showing the number of genes unique to that set. The center displays the size
#' of the core intersection (genes shared by all sets). Designed for intuitive
#' comparison of shared vs. condition-specific genes across 2–6 groups.
#'
#' @param gene_lists A named list of character vectors containing gene identifiers.
#' @param project_name Character. Prefix for the plot title.
#' @param group_colors Named character vector for specific group mappings.
#' @return Invisibly returns the \code{ggplot} object.
#' @importFrom ggplot2 ggplot geom_polygon annotate coord_fixed theme_void labs theme element_text margin
#' @importFrom ggforce geom_circle
#' @importFrom scales hue_pal
#' @export
#' @examples
#' gene_sets <- list(
#'   Control = c("TP53", "BRCA1", "MYC", "EGFR"),
#'   Treated = c("BRCA1", "MYC", "EGFR", "KRAS"),
#'   Resistant = c("MYC", "EGFR", "KRAS", "BRAF")
#' )
#' draw_flower_simplified(
#'   gene_lists = gene_sets,
#'   project_name = "Drug Response",
#'   group_colors = c(Control = "#E41A1C", Treated = "#377EB8", Resistant = "#4DAF4A")
#' )
draw_flower_simplified <- function(gene_lists, project_name, group_colors) {
  gene_lists <- gene_lists[lengths(gene_lists) > 0]
  n_groups <- length(gene_lists)
  if (n_groups < 2) {
    message("Less than 2 non-empty gene lists; skipping flower plot.")
    return(invisible(NULL))
  }
  group_names <- names(gene_lists)

  final_colors <- if (!is.null(names(group_colors))) {
    group_colors[names(gene_lists)]
  } else {
    group_colors[seq_along(gene_lists)]
  }

  # Defensive scale extraction
  if (anyNA(final_colors)) final_colors <- scales::hue_pal()(max(1, n_groups))

  core_genes <- Reduce(intersect, gene_lists)
  core_count <- length(core_genes)

  petal_counts <- vapply(group_names, function(g) {
    others_union <- unique(unlist(gene_lists[setdiff(group_names, g)]))
    unique_to_g <- setdiff(gene_lists[[g]], others_union)
    length(unique_to_g)
  }, FUN.VALUE = integer(1))

  center_x <- 0
  center_y <- 0
  ellipse_a <- 3.8
  ellipse_b <- 1.6
  r_offset <- 1.1
  core_radius <- 1.4

  get_ellipse <- function(angle_deg, cx, cy, a, b, offset, group_lbl) {
    t <- seq(0, 2 * pi, length.out = 100)
    rad <- angle_deg * pi / 180
    x <- a * cos(t)
    y <- b * sin(t)
    x_rot <- x * cos(rad) - y * sin(rad)
    y_rot <- x * sin(rad) + y * cos(rad)
    data.frame(x = x_rot + offset * cos(rad) + cx, y = y_rot + offset * sin(rad) + cy, group = group_lbl)
  }

  plot_df <- data.frame()
  group_label_df <- data.frame()
  count_label_df <- data.frame()
  angles <- seq(90, 90 + 360 * (n_groups - 1) / n_groups, length.out = n_groups)

  for (i in seq_len(n_groups)) {
    nm <- group_names[i]
    ang <- angles[i]
    coords <- get_ellipse(ang, center_x, center_y, ellipse_a, ellipse_b, r_offset, nm)
    plot_df <- rbind(plot_df, coords)

    group_lab_r <- r_offset + ellipse_a + 0.6
    group_label_df <- rbind(group_label_df, data.frame(x = center_x + group_lab_r * cos(ang * pi / 180), y = center_y + group_lab_r * sin(ang * pi / 180), label = nm, group = nm))

    count_lab_r <- r_offset + ellipse_a * 0.65
    count_label_df <- rbind(count_label_df, data.frame(x = center_x + count_lab_r * cos(ang * pi / 180), y = center_y + count_lab_r * sin(ang * pi / 180), label = as.character(petal_counts[i]), group = nm))
  }

  p <- ggplot() +
    geom_polygon(data = plot_df, aes(x = x, y = y, fill = group), alpha = 0.6, color = "white", linewidth = 0.7) +
    ggforce::geom_circle(aes(x0 = 0, y0 = 0, r = core_radius), fill = "white", color = "grey70", linewidth = 1) +
    annotate("text", x = 0, y = 0.35, label = "Core", color = "grey50", size = 5.5, fontface = "bold") +
    annotate("text", x = 0, y = -0.3, label = core_count, color = "black", size = 9, fontface = "bold") +
    geom_text(data = count_label_df, aes(x = x, y = y, label = label), color = "black", size = 7, fontface = "bold") +
    geom_text(data = group_label_df, aes(x = x, y = y, label = label, color = group), size = 6, fontface = "bold") +
    scale_fill_manual(values = final_colors) +
    scale_color_manual(values = final_colors) +
    coord_fixed() +
    theme_void() +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = 20)), plot.margin = margin(30, 30, 30, 30)) +
    labs(title = paste0(project_name, ": Simplified Flower Plot\n(Core vs. Unique)"))

  return(invisible(p))
}


#' Generate UpSet Plot for Gene Set Intersections
#'
#' Visualizes intersections among multiple gene sets using the classic UpSetR package.
#' Uses grid graphics capture to ensure plot generation in all environments.
#'
#' @param gene_lists A named list of character vectors of gene identifiers.
#' @param project_name Character. Used for the file title.
#' @param group_colors Optional mapped color schemes.
#' @return Invisibly returns the \code{grob} object.
#' @importFrom UpSetR fromList upset
#' @importFrom grid grid.grabExpr grid.draw
#' @export
#' @examples
#' gene_sets <- list(
#'   Upregulated = c("TP53", "BRCA1", "MYC", "EGFR"),
#'   Downregulated = c("BRCA1", "MYC", "CDKN1A", "BAX"),
#'   Bound_by_TF = c("MYC", "EGFR", "CDKN1A", "KRAS")
#' )
#' draw_upset_intersections(
#'   gene_lists = gene_sets,
#'   project_name = "Transcriptional Regulation"
#' )
draw_upset_intersections <- function(gene_lists, project_name = "UpSet Plot", group_colors = NULL) {
  gene_lists <- gene_lists[lengths(gene_lists) > 0]
  if (length(gene_lists) < 2) {
    message("Less than 2 non-empty gene lists; skipping UpSet plot.")
    return(invisible(NULL))
  }

  input_data <- UpSetR::fromList(gene_lists)

  if (nrow(input_data) == 0 || ncol(input_data) == 0) {
    warning("Input data for UpSet plot is empty. Skipping.")
    return(invisible(NULL))
  }


  upset_grob <- tryCatch(
    {
      .with_known_upstream_noise_suppressed(
        grid::grid.grabExpr({
          UpSetR::upset(
            input_data,
            nsets = length(gene_lists),
            nintersects = 40,
            mb.ratio = c(0.55, 0.45),
            order.by = "freq",
            mainbar.y.label = "Gene Intersection Size",
            sets.x.label = "Set Size",
            text.scale = c(1.3, 1.3, 1, 1, 1.3, 1)
          )
        })
      )
    },
    error = function(e) {
      warning("UpSetR plotting failed: ", e$message)
      return(NULL)
    }
  )

  if (is.null(upset_grob)) {
    return(invisible(NULL))
  }

  return(invisible(upset_grob))
}
