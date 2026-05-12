# Response to Reviewer Comments (Bisaloo, 2026-04-28)

We thank the reviewer for the careful and constructive first-round assessment. Below, we summarise the main blocking issues that have been addressed in the current revision and indicate the corresponding code, documentation, and testing updates. For broader stylistic points (for example, further modularisation of the longest pipeline functions), we indicate what has already been improved in this revision and what can be prioritised in a second pass if the reviewer considers that most helpful.

---

## General package development

### `R CMD check` timing

> The `annotate_peaks_and_loops()` example takes a long time. Would you be able to make it run faster by using a smaller dataset?

Yes. In response to this comment, we reduced the man-page example for `annotate_peaks_and_loops()` to a minimal tempfile-based BEDPE input together with the small sample `TxDb` SQLite database distributed with `GenomicFeatures`. The example writes to `tempdir()`, disables unnecessary file output, and therefore remains runnable while avoiding unnecessary work in the check directory.

On my local macOS arm64 development machine:

- `R CMD check --no-build-vignettes` for the current revision completed in about 3 minutes.
- A full `devtools::check()` completed in about 5.5 minutes.

This brings the `--no-build-vignettes` check below the 10-minute Bioconductor guideline.

### README file

> Why are code blocks set to `eval = FALSE` in the README?

Only two README chunks remain intentionally unevaluated:

1. The installation chunk. This is intentionally not executed during README rendering because it would attempt package installation.
2. The `looplook_report()` example. This call launches a nested R Markdown rendering step and creates report output on disk, which is not appropriate for README rendering.

All other README code chunks run normally and display output in the rendered README. The full executable workflow is documented in the vignette.

---

## The DESCRIPTION file

### `Depends`, `Imports`, `Suggests`, `Enhances`

> The package code includes a lot of `if (!requireNamespace(...))`. This makes for a very frustrating user experience [...] These dependencies should be listed in the `Imports` field.

We audited the dependency declarations and the runtime namespace checks.

- `ChIPseeker` has been moved to `Imports`, because it is used by the core annotation workflow.
- Missing packages used in optional workflows or documentation have been added to `DESCRIPTION`, including species-specific `TxDb` / `BSgenome` packages, `org.Mm.eg.db`, `DT`, `kableExtra`, and `htmltools`.
- We removed stale imports, including the unused `install.packages` import from `utils`.
- Conditional availability checks are now limited to genuinely optional downstream steps (for example, motif analysis, PPI analysis, or report/widget rendering), rather than the core annotation path.

We also appreciate the broader usability point raised here. In this revision, we focused on making the dependency graph explicit and removing the most problematic namespace inconsistencies. If the reviewer would prefer a stricter `Imports`-heavy design for additional downstream modules, we would be happy to prioritise that in the next revision.

### `SystemRequirements`

> It's uncommon to list recursive system requirements. You should only list system requirements that come from your own package.

This has been removed. `looplook` does not contain compiled code and no longer declares recursive system requirements from downstream dependencies.

---

## Documentation

### Vignette present and core functionality documented

> Vignette present and does describe the core functionality of the package.

The vignette now covers the full intended workflow of the package: installation, loop consolidation, 3D-guided annotation, expression-aware refinement, functional profiling, report generation, comparison with related tools, and `sessionInfo()`.

### Comparison with related tools

> When relevant, vignette provides review/comparison to other packages with similar functionality or scope.

We added and retained a dedicated "Comparison with Existing Tools" section in the vignette. It compares `looplook` against ChIPseeker, GREAT, GenomicInteractions, FUMA, and the ABC Model across the main capability classes relevant to this package.

### Package man page

> Please add a package man page. You can add it with `usethis::use_package_doc()`.

This has now been added. The package-level documentation is now provided via `R/looplook-package.R` / `man/looplook-package.Rd`.

### Runnable examples

> `if (bedpe_path != "")` is not necessary in examples.
>
> Same for `if (res_integrated != "")` in the vignette.

We removed these defensive guards from man-page examples and from the vignette. Examples now use direct `system.file()` paths where package data are needed. Where output files are expected, we use `tempdir()` rather than writing into the package root.

### All exported functions have man pages and runnable examples

All exported functions now have man pages. Runnable examples are provided for the user-facing exported functions:

- `annotate_peaks_and_loops()` uses a minimal temporary BEDPE input together with a sample `TxDb` SQLite file from `GenomicFeatures`.
- `plot_peaks_interactions()` uses a single temporary BEDPE example with the gene track disabled.
- `refine_loop_anchors_by_expression()` writes to `tempdir()`.
- `profile_target_genes()` now uses a lightweight runnable example based on the packaged `analysis_results.RData`, with expensive optional steps disabled.
- `looplook_report()` now has a runnable example using `precomputed_res`, `tempdir()`, and lightweight settings rather than `\donttest{}`.

Internal helpers that should not be part of the public man-page index are now documented with `@noRd`.

---

## Unit tests

> Unit tests present and covering large part of core functionality.

The placeholder test file has been replaced by actual `testthat` coverage across the core workflows.

We also added a dedicated regression test file covering, among other cases:

- BED / BEDPE 0-based to 1-based coordinate handling.
- Single-column expression-matrix gene-name preservation.
- `OrgDb` / `TxDb` object handling.
- `NULL target_annotation` in refinement.
- Semicolon-separated multi-gene handling.
- BEDPE export roundtrip consistency.
- Edge cases in conflict resolution and refined statistics.

We removed `skip_on_bioc()` from the test suite. The only strongly environment-dependent case is the network-dependent STRING integration test, which is now opt-in via `LOOPLOOK_RUN_NETWORK_TESTS=true` rather than being conditionally skipped by platform.

---

## R code

### Coding and syntax

> Coding and syntax:
> - `vapply` instead of `sapply`
> - `TRUE`, `FALSE` instead of `T`, `F`
> - no `<<-`
> - `message()`, `warning()`, `stop()` instead of `cat`

This has been addressed in the current revision:

- `sapply()` calls were replaced with `vapply()`.
- Long-form logicals are used consistently.
- `<<-` was removed from `profile_target_genes()`.
- User-facing output uses `message()`, `warning()`, or `stop()`.

We also updated deprecated ggplot2 line-width usage where the package itself was responsible (`size` -> `linewidth` in the affected ggplot2 calls, including theme elements such as `element_rect()`).

### Re-use of classes and functionality

> Re-use of classes and functionality (if appropriate).

The package now relies more explicitly on standard Bioconductor structures:

- `GRanges` and `GInteractions` remain the internal genomic / interaction representations.
- `annotate_peaks_and_loops()` now accepts `TxDb` and `OrgDb` objects directly.
- `org_db` lookup in GO analysis uses `getFromNamespace()` rather than requiring the user to attach the annotation package manually.

### Input validation and `tryCatch()`

> Do you need to wrap anything in `tryCatch()`? In R, it is more common to do robust input checking and validation.

We reviewed the remaining `tryCatch()` blocks and retained them only where graceful degradation is intentional:

- Optional downstream analyses (motif, GO, PPI, report/widget rendering).
- Workbook writing, so disk or permission problems do not discard in-memory results.
- Selected plotting/reporting steps where best-effort rendering is preferable to aborting the full workflow.

At the same time, we strengthened explicit input validation:

- Empty-input guards before `ChIPseeker::annotatePeak()`.
- Empty-region guards in `prepare_track_data()`.
- Stricter checks for missing columns and unsupported species.
- Stricter expression-matrix validation so missing or duplicated sample columns are rejected explicitly instead of being silently dropped.
- Overlap-rate warning when expression identifiers do not match annotation symbols well.

### `formatC()` / separator detection / ggplot2 deprecations

> `formatC()` handling of `NA` is not necessary.
>
> `fread()` should be able to identify the `sep` argument automatically.
>
> `size` aesthetic for lines was deprecated in ggplot2 3.4.0.

This has been addressed as follows:

- Removed the unnecessary pre-check around `formatC()`.
- Removed the manual separator-detection logic before `fread()`.
- Replaced the package-level deprecated line `size` usages with `linewidth`.

### Long functions / further refactoring

> Functional programming / no code repetition / no excessively long functions / vectorization / argument validation.

We improved several of these points in this revision, especially by:

- Reducing repeated helper logic.
- Replacing remaining `sapply()` use.
- Strengthening argument checks in user-facing functions.
- Expanding tests around edge cases and failure modes.

We did **not** perform a broad stylistic rewrite of the longest pipeline functions in this pass. We prioritised correctness, dependency cleanup, documentation, test coverage, and check compliance first. If the reviewer considers it useful, we would be happy to modularise the longest functions further in a subsequent revision.

---

## Usability

### Plot return values and report usability

> I find not displaying the plots in an interactive manner reduces the usability of the package. Letting the user manage the plots themselves would also allow them to customise the plots more easily.

We agree that returning plot objects improves both usability and downstream customisability. The plotting functions now return plot objects that can be displayed, modified, combined, or saved by the user, instead of enforcing a save-only workflow.

For users who prefer a shareable document, we also enhanced the `looplook_report()` interface and its R Markdown template under `inst/rmarkdown/templates/looplook-report/`. This optional one-command workflow wraps annotation, refinement, and profiling into a shareable HTML report, accessible from both the R console and the RStudio template menu. Template improvements in this revision include refined metric cards, colourblind-friendly default palettes, improved column-aligned tables, the ability to inject precomputed results via `precomputed_res` (as a file path or in-memory object), and companion Sankey widget rendering when that panel is present. This is a convenience wrapper, not a replacement for the core analysis functions.

### Removal of `NA` values

> Are you sure it's valid to remove `"NA"` values? Couldn't these be valid gene names?

The package filters missing values via `is.na()` / `na.omit()`, not by matching the literal character string `"NA"`. A literal gene symbol `"NA"` would be retained as a character value; only true missing values are removed.

### Colour palette provenance

> Where does this palette come from? It would be helpful to provide a reference.

We clarified the palette choices and their provenance. The package now relies on well-known colourblind-aware palettes:

- `Set2` / `Paired` from `RColorBrewer` for qualitative categories.
- `PuOr` from `RColorBrewer` for diverging fold-change scales.
- `viridis` for sequential heatmaps / density scales.

---

## Additional implementation corrections made during revision

During revision, we also corrected several implementation issues that were identified during internal validation:

- Consistent 0-based BED / BEDPE import and 1-based internal coordinate handling.
- Corrected BEDPE export back to 0-based format.
- Improved `OrgDb` object handling (`getFromNamespace` instead of `get`).
- Corrected multi-gene (semicolon-separated) splitting in refined promoter-centric statistics.
- Included expression-reclassified anchor types (`eP`, `eG`) in refined distal-element counting.
- Corrected the exported `n_members` field for consolidated BEDPE output.

These changes are covered by additional regression tests.

---

We thank the reviewer again for the careful and helpful assessment. We believe that the main blocking issues raised in this first round have now been addressed. If the reviewer considers it useful, we would be very happy to make a second pass focused specifically on dependency placement or further modularisation.

Best regards,

Ying Zhang  
on behalf of the looplook authors
