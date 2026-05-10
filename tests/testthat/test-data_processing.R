# tests/testthat/test-data_processing.R

test_that("data_processing modules run successfully on example bedpe files", {
  loop1 <- system.file("extdata", "example_loops_1.bedpe", package = "looplook")
  loop2 <- system.file("extdata", "example_loops_2.bedpe", package = "looplook")
  h3k27ac_peaks <- system.file("extdata", "example_k27ac_peaks.bed", package = "looplook")
  skip_if(loop1 == "" || loop2 == "")

  res_clean <- suppressWarnings(suppressMessages(
    consolidate_chromatin_loops(
      files = c(loop1, loop2),
      mode = "consensus",
      min_raw_score = 2,
      min_score = 5,
      gap = 1000,
      blacklist_species = "hg38",
      region_of_interest = h3k27ac_peaks,
      out_file = tempfile(fileext = ".bedpe")
    )
  ))

  expect_true(!is.null(res_clean))
})

test_that("consolidate_chromatin_loops balances clustered scores across replicates", {
  f1 <- tempfile(fileext = ".bedpe")
  f2 <- tempfile(fileext = ".bedpe")

  writeLines(
    "chr1\t0\t100\tchr1\t200\t300\t100",
    f1
  )
  writeLines(c(
    "chr1\t5\t105\tchr1\t205\t305\t20",
    "chr1\t10\t110\tchr1\t210\t310\t25",
    "chr1\t15\t115\tchr1\t215\t315\t30"
  ), f2)

  gi <- suppressMessages(consolidate_chromatin_loops(
    files = c(f1, f2),
    mode = "consensus",
    gap = 25
  ))
  expect_equal(length(gi), 1)
  expect_equal(S4Vectors::mcols(gi)$n_members[[1]], 4L)
  expect_equal(S4Vectors::mcols(gi)$n_reps[[1]], 2L)
  expect_equal(S4Vectors::mcols(gi)$score[[1]], 62.5)

  gi_keep <- suppressMessages(consolidate_chromatin_loops(
    files = c(f1, f2),
    mode = "consensus",
    gap = 25,
    min_score = 50
  ))
  expect_equal(length(gi_keep), 1)

  gi_drop <- suppressMessages(consolidate_chromatin_loops(
    files = c(f1, f2),
    mode = "consensus",
    gap = 25,
    min_score = 63
  ))
  expect_equal(length(gi_drop), 0)

  unlink(c(f1, f2))
})

test_that("consolidate_chromatin_loops respects quiet and write_output flags", {
  loop1 <- system.file("extdata", "example_loops_1.bedpe", package = "looplook")
  loop2 <- system.file("extdata", "example_loops_2.bedpe", package = "looplook")
  skip_if(loop1 == "" || loop2 == "")

  out_dir <- tempfile(pattern = "looplook_cons_nowrite_")
  out_file <- file.path(out_dir, "loops.bedpe")
  unlink(out_dir, recursive = TRUE, force = TRUE)
  expect_false(dir.exists(out_dir))

  expect_no_message({
    gi <- consolidate_chromatin_loops(
      files = c(loop1, loop2),
      mode = "consensus",
      gap = 1000,
      out_file = out_file,
      write_output = FALSE,
      quiet = TRUE
    )
  })

  expect_s4_class(gi, "GInteractions")
  expect_false(dir.exists(out_dir))
  expect_false(file.exists(out_file))
})
