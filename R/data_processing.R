#' Read BEDPE File into a GInteractions Object
#'
#' Reads a standard BEDPE file and converts it into a Bioconductor
#' \code{\link[InteractionSet]{GInteractions}} object.
#'
#' @details
#' \strong{Anchor Normalization:}
#' Anchor order is automatically normalized so that the first anchor is
#' lexicographically less than or equal to the second (e.g., chr1 < chr2),
#' ensuring compatibility with \code{GInteractions(mode = "strict")}.
#'
#' \strong{Score Detection:}
#' The function attempts to automatically detect a numeric score column.
#' It checks the 8th column first (standard for many tools); if not numeric,
#' it falls back to the 7th column. Non-numeric values are treated as 0.
#'
#' @param bedpe_file Character. Path to a BEDPE file. Must contain at least six columns:
#'   \code{chr1, start1, end1, chr2, start2, end2}.
#' @param score_col Integer or \code{NULL}. Column index to use as interaction
#'   score (e.g. PET count, -log10(p-value)). If \code{NULL} (default), column 8
#'   is tried first, then column 7; if neither is numeric, scores default to 0.
#'   Set explicitly when the score column position differs from the standard or
#'   when auto-detection picks the wrong column. Note: downstream filtering
#'   parameters such as \code{min_score} in
#'   \code{\link{consolidate_chromatin_loops}} assume \emph{higher scores = better
#'   interactions}. If your score column contains p-values or other metrics where
#'   \emph{smaller is better}, convert it first (e.g. \code{-log10(p)}) or filter
#'   manually after consolidation.
#' @return A \code{\link[InteractionSet]{GInteractions}} object with a \code{score} metadata column
#'   (defaulting to 0 if not provided).
#' @importFrom data.table fread
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom InteractionSet GInteractions
#' @importFrom S4Vectors mcols mcols<-
#' @export
#' @examples
#' # 1. Locate the example BEDPE file included in the package
#' # system.file finds the absolute path to 'inst/extdata/example_loops_1.bedpe'
#' bedpe_path <- system.file("extdata", "example_loops_1.bedpe", package = "looplook")
#'
#' # 2. Run the function (ensure file was found)
#' gi <- bedpe_to_gi(bedpe_path)
#'
#' # 3. Inspect the result
#' print(gi)
#'
#' # Check the imported score column
#' S4Vectors::mcols(gi)$score
bedpe_to_gi <- function(bedpe_file, score_col = NULL) {
  if (is.null(bedpe_file) || !file.exists(bedpe_file)) {
    stop("BEDPE file does not exist or path is invalid: ", bedpe_file)
  }
  df <- data.table::fread(bedpe_file, header = FALSE)

  if (ncol(df) < 6) {
    stop("BEDPE file must have at least 6 columns: ", bedpe_file)
  }

  # Convert coordinate columns to numeric before validation
  coord_cols <- c(2, 3, 5, 6)
  for (cc in coord_cols) {
    df[[cc]] <- suppressWarnings(as.numeric(df[[cc]]))
  }
  if (any(is.na(df[, ..coord_cols]))) {
    stop("BEDPE file contains non-numeric coordinate columns: ", bedpe_file)
  }
  if (any(df[[2]] > df[[3]], na.rm = TRUE) || any(df[[5]] > df[[6]], na.rm = TRUE)) {
    stop("BEDPE file contains rows with start > end coordinate: ", bedpe_file)
  }

  df <- as.data.frame(df)

  swap <- (df[, 1] > df[, 4]) | (df[, 1] == df[, 4] & df[, 2] > df[, 5])
  if (any(swap)) {
    df[swap, seq_len(6)] <- df[swap, c(4, 5, 6, 1, 2, 3)]
  }

  # BEDPE is 0-based half-open; GRanges expects 1-based closed
  gr1 <- .with_known_upstream_noise_suppressed(
    GenomicRanges::GRanges(df[, 1], IRanges::IRanges(df[, 2] + 1, df[, 3]))
  )
  gr2 <- .with_known_upstream_noise_suppressed(
    GenomicRanges::GRanges(df[, 4], IRanges::IRanges(df[, 5] + 1, df[, 6]))
  )

  gi <- .with_known_upstream_noise_suppressed(
    InteractionSet::GInteractions(gr1, gr2, mode = "strict")
  )

  .numeric_ratio <- function(x) {
    x <- as.character(x)
    x <- x[!is.na(x) & trimws(x) != ""]
    if (length(x) == 0) {
      return(0)
    }
    mean(!is.na(suppressWarnings(as.numeric(x))))
  }

  final_scores <- rep(0, nrow(df))
  found <- FALSE

  if (!is.null(score_col)) {
    if (!is.numeric(score_col) || length(score_col) != 1L ||
        is.na(score_col) || score_col < 1 ||
        score_col != as.integer(score_col)) {
      stop("score_col must be a single positive integer column index.", call. = FALSE)
    }
    if (score_col > ncol(df)) {
      stop("score_col = ", score_col, " exceeds file column count (", ncol(df), ")")
    }
    final_scores <- as.numeric(df[[score_col]])
    found <- TRUE
  } else {
    if (ncol(df) >= 8 && .numeric_ratio(df[[8]]) >= 0.5) {
      final_scores <- suppressWarnings(as.numeric(df[[8]]))
      found <- TRUE
    }
    if (!found && ncol(df) >= 7 && .numeric_ratio(df[[7]]) >= 0.5) {
      final_scores <- suppressWarnings(as.numeric(df[[7]]))
      found <- TRUE
    }
    if (!found && ncol(df) >= 7) {
      warning("Scores defaulted to 0: no column beyond 6 had >50% numeric values")
    }
  }

  final_scores[is.na(final_scores)] <- 0
  S4Vectors::mcols(gi)$score <- final_scores

  return(gi)
}


#' Spatial Clustering of GInteractions
#'
#' Merges spatially proximal chromatin loops into consensus interactions using graph-based clustering.
#' Loops are considered overlapping if both anchors are within \code{gap} bp of each other.
#' Each resulting cluster is represented by the \strong{union genomic range (min start to max end)} spanning all its members.
#'
#' @param gi A \code{\link[InteractionSet]{GInteractions}} object.
#' @param gap Numeric. Maximum distance (in base pairs) allowed between anchors to consider two loops overlapping. Default: 1000.
#' @return A list with two elements:
#' \describe{
#'   \item{\code{gi}}{Reduced \code{\link[InteractionSet]{GInteractions}} object, one per cluster.}
#'   \item{\code{membership}}{Integer vector indicating cluster assignment for each input loop.}
#' }
#' Metadata columns include \code{cluster_id}, \code{n_members}, and averaged \code{score}.
#' @importFrom igraph make_empty_graph add_edges components
#' @importFrom S4Vectors queryHits subjectHits mcols
#' @importFrom GenomicRanges GRangesList
#' @export
#' @examples
#' # 1. Load example data (loops that are close to each other)
#' bedpe_path <- system.file("extdata", "example_loops_1.bedpe", package = "looplook")
#' # Convert BEDPE to GInteractions object
#' gi_raw <- bedpe_to_gi(bedpe_path)
#'
#' # 2. Run clustering
#' res <- reduce_ginteractions(gi_raw, gap = 1000)
#'
#' # 3. Inspect results
#' print(res$gi)
#' head(res$membership)
#' table(res$membership)
reduce_ginteractions <- function(gi, gap = 1000) {
  if (length(gi) == 0) {
    return(list(gi = gi, membership = integer(0)))
  }

  dt <- gi_to_dt(gi)
  dt <- cluster_loops_dt(dt, gap)
  reduced_dt <- reduce_clusters_dt(dt)
  gi_red <- dt_to_gi(reduced_dt)

  list(gi = gi_red, membership = dt$cluster)
}

#' Read a Simple BED File into a GRanges Object
#'
#' Reads the first three columns of a BED file (chrom, start, end) and returns a
#' \code{\link[GenomicRanges]{GRanges}} object. Additional columns are ignored.
#'
#' @param bed_file Character. Path to a BED file (must be tab-delimited).
#' @return A \code{\link[GenomicRanges]{GRanges}} object, or \code{NULL} if \code{bed_file} is \code{NULL}.
#' @importFrom data.table fread
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @export
#' @examples
#' # 1. Locate the example BED file included in the package
#' bed_path <- system.file("extdata", "example_peaks.bed", package = "looplook")
#'
#' # 2. Run the function
#' gr <- read_simple_bed(bed_path)
#'
#' # 3. Inspect the result
#' print(gr)
#' length(gr)
read_simple_bed <- function(bed_file) {
  if (is.null(bed_file)) {
    return(NULL)
  }
  if (!file.exists(bed_file)) {
    stop("BED file does not exist: ", bed_file)
  }

  df <- data.table::fread(bed_file, header = FALSE, select = c(1, 2, 3))
  # Auto-detect header: only examine the first row to avoid misidentifying
  # a malformed data row as a header line.
  has_header <- nrow(df) > 0 && (
    is.na(suppressWarnings(as.numeric(df[[2]][1]))) ||
    is.na(suppressWarnings(as.numeric(df[[3]][1]))))
  if (has_header) {
    df <- data.table::fread(bed_file, header = FALSE, select = c(1, 2, 3), skip = 1)
  }
  # Validate all coordinate columns are numeric after header handling
  df[[2]] <- suppressWarnings(as.numeric(df[[2]]))
  df[[3]] <- suppressWarnings(as.numeric(df[[3]]))
  if (any(is.na(df[[2]]) | is.na(df[[3]]))) {
    stop("BED file contains non-numeric start/end coordinates after header handling.",
      call. = FALSE)
  }
  if (nrow(df) > 0 && any(df[[2]] > df[[3]], na.rm = TRUE)) {
    stop("BED file contains intervals with start > end.", call. = FALSE)
  }
  # BED is 0-based half-open; GRanges expects 1-based closed
  .with_known_upstream_noise_suppressed(
    GenomicRanges::GRanges(df[[1]], IRanges::IRanges(df[[2]] + 1, df[[3]]))
  )
}


#' Consolidate and Integrate Chromatin Loops from Replicates or Multiple Sources
#'
#' @description
#' This function consolidates chromatin loops from multiple BEDPE files. It is designed for two main purposes:
#' \enumerate{
#'   \item \strong{Replicate Consolidation}: Merging biological or technical replicates to identify high-confidence, reproducible loops (e.g., 3 replicates of H3K27ac HiChIP).
#'   \item \strong{Multi-Omics Integration}: The framework can be used to identify multi-source consensus by integrating datasets from various experimental designs, such as HiChIP assays targeting different factors (e.g., integrating \strong{H3K27ac} and \strong{H3K4me3}, or overlapping Hi-C with ChIA-PET).
#' }
#'
#' The function supports three modes:
#' \itemize{
#'   \item \code{"consensus"}: Implements graph-based connected component analysis to cluster spatially proximal anchors across samples. Only retains clusters detected in >= min_consensus biological replicates.
#'   \item \code{"intersect"}: Enforces strict reference-based filtering, retaining loops that show full genomic overlap with the reference file (File 1).
#'   \item \code{"union"}: Retains all chromatin interactions across the entire cohort, ideal for exploratory pan-tissue analyses.
#' }
#'
#' \strong{Connected-component chaining}: Graph-based clustering may
#' transitively chain loci (A–B and B–C merge, pulling A and C into the same
#' cluster even if they are far apart). Inspect \code{n_members} in the
#' output and reduce \code{gap} if clusters appear inflated.
#'
#' It also supports a \strong{two-stage filtering strategy} to maximize signal-to-noise ratio:
#' \itemize{
#'   \item \strong{Pre-filtering} (\code{min_raw_score}): Removes low-confidence noise (e.g., singleton reads) from raw files \emph{before} merging.
#'   \item \strong{Post-filtering} (\code{min_score}): Filters the final consensus loops based on their aggregated score.
#'   \item \strong{Replicate-Balanced Aggregation}: In \code{"consensus"} and \code{"union"} modes, each cluster score is computed as the mean of per-source mean scores, so one replicate with many fragmented calls cannot dominate the final representative score.
#' }
#'
#' @param files Character vector. Paths to BEDPE files (at least two).
#' @param gap Numeric. Distance (bp) to consider loops as overlapping. Default 1000.
#' @param mode Character. Choose one of the following: "consensus", "intersect", "union". Merge strategy:
#'   \itemize{
#'     \item \code{"intersect"}: Strict reference-based filtering (keeps loops in File 1 supported by ALL other files).
#'     \item \code{"union"}: Merges all detected loops into a comprehensive map.
#'     \item \code{"consensus"}: Graph-based clustering to find a consensus set supported by a majority of samples. (Formerly "reproducible").
#'   }
#' @param min_consensus Integer. Minimum number of replicates a loop must appear in
#'   (only effective when \code{mode = "consensus"}).
#'   If \code{NULL} (default), the threshold is automatically calculated:
#'   \itemize{
#'     \item For 2–3 replicates: Requires all (100\%).
#'     \item For \eqn{N \ge 4}: Requires \eqn{\lfloor 0.75N \rfloor + 1}
#'       (e.g., 3 for N=4, 4 for N=5, 7 for N=8).
#'   }
#' @param min_raw_score Numeric. \strong{Pre-filtering threshold}. Loops with a raw score (e.g., read count) below this value in individual files will be discarded \strong{before} any merging or intersection.
#'   \itemize{
#'     \item Recommended value: \code{2} (to remove singleton noise loops with count=1).
#'     \item Default: \code{NULL} (no pre-filtering).
#'   }
#' @param min_score Numeric. \strong{Post-filtering threshold}. Minimum score to keep a consolidated loop \strong{after} merging.
#'   \itemize{
#'     \item For \code{"consensus"} and \code{"union"} modes, this filters loops based on a replicate-balanced representative score (mean of per-source cluster means).
#'     \item For \code{"intersect"} mode, this filters the retained File 1 loops by their original score.
#'     \item Default: \code{NULL} (no post-filtering).
#'   }
#' @param score_col Integer or \code{NULL}. Column index to use as interaction
#'   score when reading BEDPE files. Passed to \code{\link{bedpe_to_gi}}.
#'   If \code{NULL} (default), auto-detection is used (see \code{\link{bedpe_to_gi}}).
#' @param blacklist_species Character. Species/build for built-in blacklist
#'   (e.g., "hg38", "hg19", "mm10", "mm9"), or a path to a custom BED file.
#' @param region_of_interest Character. Path to BED file defining regions of interest (ROI). Only loops overlapping these regions will be kept.
#' @param roi_mode Character. How loops must overlap \code{region_of_interest}.
#'   \code{"any"} (default): keep loops where \emph{either} anchor overlaps the ROI
#'   (suitable for promoter-centric or enhancer-gene queries).
#'   \code{"both"}: keep loops where \emph{both} anchors overlap the ROI
#'   (suitable for TAD confinement or domain-internal interaction queries).
#' @param out_file Character. The file name (including the file path) for saving results in the extended BEDPE format.
#' @param write_output Logical. If \code{TRUE} (default), write the consolidated BEDPE file when \code{out_file} is provided. If \code{FALSE}, return the \code{GInteractions} object without creating directories or files.
#' @param quiet Logical. If \code{TRUE}, suppress progress messages while preserving warnings. Default: \code{FALSE}.
#' @return A filtered \code{\link[InteractionSet]{GInteractions}} object with metadata columns:
#'   \describe{
#'     \item{\code{score}}{Replicate-balanced consensus score.}
#'     \item{\code{n_members}}{Number of raw loops merged into this entry
#'       (1 for intersect mode where no coordinate merging occurs).}
#'     \item{\code{n_reps}}{Number of input files that support this entry.}
#'   }
#'   When \code{write_output = TRUE} and \code{out_file} is provided, an extended
#'   BEDPE file is written with the additional columns \code{n_members} and
#'   \code{n_reps} appended after the standard BEDPE fields.
#' @importFrom data.table fread
#' @importFrom GenomicRanges GRanges seqnames start end
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors mcols mcols<- queryHits subjectHits
#' @importFrom igraph make_empty_graph add_edges components
#' @export
#' @examples
#' # 1. Get paths to example BEDPE files included in the package
#' f1 <- system.file("extdata", "example_loops_1.bedpe", package = "looplook")
#' f2 <- system.file("extdata", "example_loops_2.bedpe", package = "looplook")
#'
#' # 2. Run consolidation (ensure files exist)
#' # Example A: Intersect Mode
#' # Only keeps loops present in f1 that are also supported by f2
#' res_intersect <- consolidate_chromatin_loops(
#'   files = c(f1, f2),
#'   mode = "intersect",
#'   gap = 1000,
#'   out_file = tempfile(fileext = ".bedpe")
#' )
#'
#' # Example B: Consensus Mode (formerly Reproducible)
#' # Finds consensus loops supported by both replicates (default for N=2)
#' res_consensus <- consolidate_chromatin_loops(
#'   files = c(f1, f2),
#'   mode = "consensus",
#'   gap = 1000,
#'   out_file = tempfile(fileext = ".bedpe")
#' )
#'
#' # Example C: Union Mode
#' # Merges all loops into a single map
#' res_union <- consolidate_chromatin_loops(
#'   files = c(f1, f2),
#'   mode = "union",
#'   gap = 1000,
#'   out_file = tempfile(fileext = ".bedpe")
#' )
#'
#' # Example D: Dual Filtering Strategy (Recommended for HiChIP)
#' # 1. Pre-filter: Discard singletons (score < 2) to remove noise.
#' # 2. Merge: Find loops present in both replicates.
#' # 3. Post-filter: Keep only strong consensus loops (score > 5).
#' res_clean <- consolidate_chromatin_loops(
#'   files = c(f1, f2),
#'   mode = "consensus",
#'   min_raw_score = 2, # Pre-filter (remove noise)
#'   min_score = 5, # Post-filter (keep strong loops)
#'   gap = 1000,
#'   out_file = tempfile(fileext = ".bedpe")
#' )
#'
#' # Inspect results
#' length(res_intersect)
#' length(res_clean)
consolidate_chromatin_loops <- function(
  files = NULL,
  gap = 1000,
  mode = c("consensus", "intersect", "union"),
  min_consensus = NULL,
  score_col = NULL,
  min_raw_score = NULL,
  min_score = NULL,
  blacklist_species = NULL,
  region_of_interest = NULL,
  roi_mode = c("any", "both"),
  out_file = NULL,
  write_output = TRUE,
  quiet = FALSE
) {
  log_message <- function(...) {
    if (!quiet) message(...)
  }

  if (is.null(files) || length(files) < 2) {
    stop("`files` must contain at least two BEDPE file paths.", call. = FALSE)
  }
  missing_files <- files[!file.exists(files)]
  if (length(missing_files) > 0) {
    stop("BEDPE files not found: ", paste(missing_files, collapse = ", "), call. = FALSE)
  }
  mode <- match.arg(mode)
  roi_mode <- match.arg(roi_mode)
  n_reps <- length(files)

  log_message(">>> Reading BEDPE files")

  gi_list <- lapply(files, function(f) {
    gi <- bedpe_to_gi(f, score_col = score_col)

    if (!is.null(min_raw_score)) {
      if ("score" %in% colnames(S4Vectors::mcols(gi))) {
        keep_idx <- S4Vectors::mcols(gi)$score >= min_raw_score
        gi <- gi[keep_idx]
      }
    }
    return(gi)
  })

  for (i in seq_along(gi_list)) {
    S4Vectors::mcols(gi_list[[i]])$source <- i
    log_message("    File ", i, ": ", length(gi_list[[i]]), " loops")
  }

  result_gi <- NULL

  if (mode == "intersect") {
    log_message(">>> Intersect mode: Reference-based filtering (No Coordinate Merging)")
    log_message("    Base: File 1. Criterion: Must overlap with ALL other files.")

    current_gi <- gi_list[[1]]

    for (i in 2:n_reps) {
      if (length(current_gi) == 0) break
      log_message("    Intersecting with File ", i, "...")

      hits <- InteractionSet::findOverlaps(
        current_gi,
        gi_list[[i]],
        maxgap = gap,
        use.region = "both"
      )

      keep_idx <- unique(S4Vectors::queryHits(hits))
      current_gi <- current_gi[keep_idx]
    }

    result_gi <- current_gi
    S4Vectors::mcols(result_gi)$n_reps <- n_reps
    S4Vectors::mcols(result_gi)$n_members <- 1L
  } else {
    log_message(">>> Clustering mode (Union/Consensus): Merging coordinates via Graph")

    combined_dt <- data.table::rbindlist(lapply(gi_list, gi_to_dt))
    clustered <- cluster_loops_dt(combined_dt, gap)
    reduced_dt <- reduce_clusters_dt(clustered)


    if (mode == "consensus") {
      if (is.null(min_consensus)) {
        if (n_reps == 2) {
          min_consensus <- 2
        } else if (n_reps == 3) {
          min_consensus <- 3
        } else if (n_reps == 4) {
          min_consensus <- 3
        } else {
          min_consensus <- floor(0.75 * n_reps) + 1
        }
      }
      log_message(">>> Consensus mode: Keeping clusters in >= ", min_consensus, " replicates")
      reduced_dt <- reduced_dt[n_reps >= min_consensus]
    } else {
      log_message(">>> Union mode: Keeping all clusters")
    }

    result_gi <- dt_to_gi(reduced_dt)
  }


  if (!is.null(min_score)) {
    keep <- S4Vectors::mcols(result_gi)$score >= min_score
    result_gi <- result_gi[keep]
  }

  if (!is.null(blacklist_species)) {
    known_lists <- list(
      "hg38" = "hg38-blacklist.v2.bed",
      "hg19" = "hg19-blacklist.v2.bed",
      "mm10" = "mm10-blacklist.v2.bed",
      "mm9"  = "mm9-blacklist.v2.bed"
    )
    blacklist_path <- NULL
    if (blacklist_species %in% names(known_lists)) {
      blacklist_path <- system.file("extdata", known_lists[[blacklist_species]], package = "looplook")
    }
    if (is.null(blacklist_path) || blacklist_path == "") {
      blacklist_path <- blacklist_species
    }
    if (file.exists(blacklist_path)) {
      log_message(">>> Filtering blacklist: ", basename(blacklist_path))
      bl <- read_simple_bed(blacklist_path)
      bl <- .harmonize_seqlevels(bl, InteractionSet::anchors(result_gi, "first"), "blacklist")
      h1 <- InteractionSet::findOverlaps(InteractionSet::anchors(result_gi, "first"), bl)
      h2 <- InteractionSet::findOverlaps(InteractionSet::anchors(result_gi, "second"), bl)
      bad <- unique(c(S4Vectors::queryHits(h1), S4Vectors::queryHits(h2)))
      if (length(bad)) result_gi <- result_gi[-bad]
    } else {
      warning("Blacklist file not found: ", blacklist_species)
    }
  }

  if (!is.null(region_of_interest)) {
    log_message(">>> Filtering by region of interest (", roi_mode, "): ", basename(region_of_interest))
    if (file.exists(region_of_interest)) {
      tg <- read_simple_bed(region_of_interest)
      tg <- .harmonize_seqlevels(tg, InteractionSet::anchors(result_gi, "first"), "ROI")

      h1 <- InteractionSet::findOverlaps(InteractionSet::anchors(result_gi, "first"), tg)
      h2 <- InteractionSet::findOverlaps(InteractionSet::anchors(result_gi, "second"), tg)

      keep <- if (roi_mode == "any") {
        base::union(S4Vectors::queryHits(h1), S4Vectors::queryHits(h2))
      } else {
        base::intersect(S4Vectors::queryHits(h1), S4Vectors::queryHits(h2))
      }

      if (length(keep) > 0) {
        result_gi <- result_gi[keep]
        log_message("    Kept ", length(result_gi), " loops overlapping ROI.")
      } else {
        log_message("    No loops overlapped with the ROI. Returning empty set.")
        result_gi <- result_gi[0]
      }
    } else {
      warning("Region of interest file not found: ", region_of_interest)
    }
  }

  if (write_output && !is.null(out_file)) {
    a1 <- InteractionSet::anchors(result_gi, "first")
    a2 <- InteractionSet::anchors(result_gi, "second")

    # BEDPE export: convert back from 1-based closed to 0-based half-open
    out_df <- data.frame(
      chr1 = as.character(GenomicRanges::seqnames(a1)),
      start1 = GenomicRanges::start(a1) - 1L,
      end1 = GenomicRanges::end(a1),
      chr2 = as.character(GenomicRanges::seqnames(a2)),
      start2 = GenomicRanges::start(a2) - 1L,
      end2 = GenomicRanges::end(a2),
      name = ".",
      score = round(S4Vectors::mcols(result_gi)$score, 2),
      n_members = if (!is.null(S4Vectors::mcols(result_gi)$n_members)) {
        S4Vectors::mcols(result_gi)$n_members
      } else {
        rep(1L, length(result_gi))
      },
      n_reps = if (!is.null(S4Vectors::mcols(result_gi)$n_reps)) {
        S4Vectors::mcols(result_gi)$n_reps
      } else {
        rep(1L, length(result_gi))
      },
      stringsAsFactors = FALSE
    )

    dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
    utils::write.table(out_df, out_file,
      sep = "\t",
      row.names = FALSE, col.names = FALSE, quote = FALSE
    )
    log_message("Finished! Saved to ", out_file)
  }

  log_message("Finished! Final loops: ", length(result_gi))
  result_gi
}


# --- Helpers ---

gi_to_dt <- function(gi) {
  a1 <- InteractionSet::anchors(gi, "first")
  a2 <- InteractionSet::anchors(gi, "second")

  data.table::data.table(
    chr1 = as.character(GenomicRanges::seqnames(a1)),
    start1 = GenomicRanges::start(a1),
    end1 = GenomicRanges::end(a1),
    chr2 = as.character(GenomicRanges::seqnames(a2)),
    start2 = GenomicRanges::start(a2),
    end2 = GenomicRanges::end(a2),
    score = if (!is.null(S4Vectors::mcols(gi)$score)) {
      S4Vectors::mcols(gi)$score
    } else {
      0
    },
    source = if (!is.null(S4Vectors::mcols(gi)$source)) {
      S4Vectors::mcols(gi)$source
    } else {
      1L
    }
  )
}

reduce_clusters_dt <- function(dt) {
  cluster_coords <- dt[, list(
    chr1 = chr1[1],
    start1 = min(start1),
    end1 = max(end1),
    chr2 = chr2[1],
    start2 = min(start2),
    end2 = max(end2),
    n_members = .N
  ), by = cluster]

  cluster_scores <- dt[, .(
    score = .mean_or_na(score)
  ), by = .(cluster, source)][!is.na(score), .(
    score = .mean_or_na(score),
    n_reps = .N
  ), by = cluster]

  reduced_dt <- merge(
    cluster_coords,
    cluster_scores,
    by = "cluster",
    all.x = TRUE,
    sort = FALSE
  )
  data.table::setcolorder(reduced_dt, c(
    "cluster", "chr1", "start1", "end1", "chr2", "start2", "end2",
    "score", "n_members", "n_reps"
  ))
  reduced_dt
}

.mean_or_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  mean(x)
}

dt_to_gi <- function(dt) {
  gr1 <- .with_known_upstream_noise_suppressed(
    GenomicRanges::GRanges(
      dt$chr1, IRanges::IRanges(dt$start1, dt$end1)
    )
  )
  gr2 <- .with_known_upstream_noise_suppressed(
    GenomicRanges::GRanges(
      dt$chr2, IRanges::IRanges(dt$start2, dt$end2)
    )
  )

  gi <- .with_known_upstream_noise_suppressed(
    InteractionSet::GInteractions(gr1, gr2, mode = "strict")
  )
  S4Vectors::mcols(gi)$cluster_id <- dt$cluster
  S4Vectors::mcols(gi)$n_members <- dt$n_members
  S4Vectors::mcols(gi)$score <- dt$score
  S4Vectors::mcols(gi)$n_reps <- dt$n_reps
  gi
}

cluster_loops_dt <- function(dt, gap) {
  dt[, idx := .I]
  dt[, `:=`(
    a1_l = start1 - gap - 1L,
    a1_r = end1 + gap + 1L,
    a2_l = start2 - gap - 1L,
    a2_r = end2 + gap + 1L
  )]

  hits <- dt[dt, on = .(
    chr1 = chr1,
    a1_l <= end1,
    a1_r >= start1,
    chr2 = chr2,
    a2_l <= end2,
    a2_r >= start2
  ), nomatch = NULL, allow.cartesian = TRUE]

  edges <- hits[idx < i.idx, .(from = idx, to = i.idx)]

  g <- igraph::make_empty_graph(n = nrow(dt), directed = FALSE)
  if (nrow(edges) > 0) {
    edge_vec <- as.vector(t(as.matrix(edges)))
    g <- igraph::add_edges(g, edge_vec)
  }

  comp <- igraph::components(g)
  dt[, cluster := comp$membership]
  dt[, `:=`(idx = NULL, a1_l = NULL, a1_r = NULL, a2_l = NULL, a2_r = NULL)]
  dt
}

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "V1", "V2", "V3", "V4", "V5", "V6", "V7",
    "chr1", "start1", "end1", "chr2", "start2", "end2",
    "idx", "i.idx", "cluster", "score", "source", "n_members", "n_reps",
    "a1_l", "a1_r", "a2_l", "a2_r", ".N", ".I", "..coord_cols"
  ))
}
