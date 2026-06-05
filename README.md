
# looplook <img src="man/figures/logo.png" align="right" width="160" alt="Looplook Logo" />

An integrative suite for expression-aware target assignment and
functional annotation of chromatin interactions.

[![R-CMD-check](https://github.com/zying106/looplook/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/zying106/looplook/actions/workflows/R-CMD-check.yaml)
[![License: GPL (\>=
3)](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Bioconductor
Ready](https://img.shields.io/badge/Bioconductor-Ready-success.svg)](#)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)

------------------------------------------------------------------------

## Introduction

Welcome to **`looplook`**, a versatile R/Bioconductor toolkit developed
to **integrate 3D chromatin architecture data** (e.g., HiChIP, ChIA-PET,
Hi-C) with **other tabular omics datasets**, including transcriptomics,
chromatin accessibility, protein-DNA interactions (derived from
ChIP-seq, CUT&Tag, or CUT&RUN), and genetic variants annotated by
genome-wide association studies.

Numerous studies have demonstrated that many distal regulatory elements
physically interact with target gene promoters via 3D chromatin
loopings, thereby regulating the expression of genes located tens of
kilobases to megabases away in the linear genome. However, conventional
annotations tend to assign putative elements (peaks) to their **nearest
genes in cis**, which often fails to reflect biological reality. Hence,
the accurate assignment of non-coding genetic variants or orphan peaks
to their cognate target genes remains a major bottleneck in the target
annotation of functional elements. To address this, `looplook`
systematically prioritizes **physical spatial chromatin contacts** to
batch-annotate thousands of regulatory elements at a **genome-wide**,
**high-throughput scale**, thereby identifying their candidate target
genes with high confidence and unprecedented efficiency.

Beyond its utility as a tool for integrative target annotation,
`looplook` can be used as a **standalone utility for loop analysis** per
se. Even in the absence of auxiliary omics data, it systematically
annotates the 3D chromatin interactome itself, **classifying complex
spatial topologies** (e.g., Enhancer-Promoter, Promoter-Promoter
interactions) and **quantifying node connectivity** to uncover **dense
regulatory hubs and enhancer cliques** (e.g., super-enhancers) that
drive cell-type-specific transcriptional programs.

------------------------------------------------------------------------

## Installation

`looplook` extensively leverages the Bioconductor ecosystem for **robust
genomic arithmetic and annotation**. To ensure optimal compatibility,
please ensure your **system environment** is fully up to date prior to
installation:

``` r
# Installation from GitHub
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db"))

if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("zying106/looplook")
```

------------------------------------------------------------------------

## Detailed Workflow & Core Modules

### Module 1: Data Consolidation & Preprocessing

In 3D genomics analyses, individual replicates typically exhibit
**certain degrees of inconsistency or noise**. The
`consolidate_chromatin_loops` function serves as the foundational
data-cleaning module, merging multiple replicates into a **standardized,
unified 3D chromatin interaction coordinate framework**.

**Key Parameters:**

- **`mode`**: Defines the overarching merging algorithm. The
  **`"consensus"`** mode (recommended) **employs graph-based connected
  component analysis** to cluster nearby chromatin loop anchors across
  biological/technical samples. The **`"intersect"`** mode applies
  strict reference-based filtering to retain only overlapping
  interactions, while the **`"union"`** mode retains all detected
  interactions for **exploratory pan-tissue analyses**.
- **`min_raw_score` & `min_score`** (The Dual-Filter): `min_raw_score`
  acts as a **pre-filter** applied to individual BEDPE files before
  clustering (e.g., removing singleton noise to substantially reduce
  computational memory overhead). `min_score` serves as a
  **post-filter** applied to the final merged chromatin interactome to
  ensure high-confidence interactions. In `"consensus"` and `"union"`
  modes, the representative score is **replicate-balanced**: the package
  first averages clustered loop scores within each replicate, then
  averages across replicates, so one replicate with denser loop calls
  does not dominate the final score.
- **`gap`**: Defines the **maximum spatial distance** (in base pairs)
  allowed between loop anchors for consideration as part of the same
  physical cluster.
- **`blacklist_species`**: Automatically excludes chromatin loops
  overlapping with high-variance, artifact-prone genomic regions (e.g.,
  centromeres, telomeres) by integrating the official ENCODE blacklist
  for specified species (e.g., `"hg38"`, `"mm10"`).
- **`region_of_interest`**: Accepts an **auxiliary BED file** (e.g., a
  specific disease-associated locus or ChIP-seq peak set) to exclude
  global background interactions, outputting only loops with physical
  connectivity to the target genomic region.

``` r
library(looplook)
out_dir <- tempdir()

f1 <- system.file("extdata", "example_loops_1.bedpe", package = "looplook")
f2 <- system.file("extdata", "example_loops_2.bedpe", package = "looplook")

consensus_global <- consolidate_chromatin_loops(
  files = c(f1, f2),
  mode = "consensus",
  gap = 1000,
  out_file = file.path(out_dir, "consensus_loops.bedpe")
)
```

### Module 2: 3D-Guided Peak Annotation & Mapping

This module serves as the **core target mapping engine**. It seamlessly
executes a rigorous hierarchical pipeline to resolve locus assignment
conflicts: Expression Pre-filter → Functional Biotype Prioritization →
Dominant Expression Tiebreaker.

**Key Parameters:**

- **`target_bed`**: pecifies the **auxiliary genomic features of
  interest** (e.g., GWAS SNPs, ATAC-seq peaks, or transcription factor
  binding sites) that require spatial target gene assignment.
- **`expr_matrix_file` & `sample_columns`**: Providing an RNA-seq matrix
  allows the engine to activate the Expression Pre-filter and Tiebreaker
  logic, drastically reducing false-positive gene assignments in
  **genomic regions harboring multiple genes**.
- **`neighbor_hop`**: An advanced topological parameter for 3D chromatin
  interactome network traversal. A value of 0 restricts annotation to
  direct physical chromatin interactions only; a value of 1 (Hub Mode)
  extends analysis to secondary network effects within enhancer cliques.
- **`tss_region`**: Defines the **spatial boundary of gene promoters**
  relative to the Transcription Start Site (TSS).

``` r
# Use pre-computed example annotation (or run annotate_peaks_and_loops with example data)
rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
load(rdata_path) # loads res_integrated
```

**Output Data Dictionary: The Comprehensive 3D Spatial Catalog**

The module systematically exports a multi-layered tabular catalog (e.g.,
`*_Basic_Results.xlsx`) detailing the spatial interactome:

- **Integrative Target Mapping (`target_annotation`)** Delineates the
  spatial coverage of genomic variants/features inputted by the user. It
  separates strict loop-derived columns (`Regulated_promoter_genes`,
  `Assigned_Target_Genes`) from `*_Filled` fallback columns, and records
  provenance through `Regulated_promoter_Evidence`,
  `Regulated_promoter_Fallback_Evidence`, and the long-format
  `target_gene_links` table. This keeps promoter-supported loop targets,
  historical assigned targets, and nearest-gene fallback choices
  traceable rather than conflated.

- **3D Network Architecture (`loop_annotation`)** Resolves the
  biological syntax of the structural interactome. It classifies
  topological interactions (loop_type, e.g., E-P, P-P) and distinguishes
  biologically relevant `Putative_Target_Genes` from the broader, raw
  physical interaction footprints (`All_Anchor_Genes`).

- **Topological Hub Detection** Quantifies structural node degrees
  (e.g., `n_Linked_Promoters`, `n_Linked_Distal`) to rigorously
  deconstruct the 3D interactome from two complementary perspectives:
  **`promoter_centric_stats`**: Identifies core target genes regulated
  by **complex regulatory architectures** (e.g., enhancer arrays or
  transcription factories)，while **`distal_element_stats`** highlights
  high-connectivity non-coding regions to facilitate the discovery of
  putative enhancer cliques.

<div align="center">

<img src="man/figures/g1_anno.jpg" width="800" style="border: 1px solid #ddd; border-radius: 4px; padding: 5px;" alt="Annotation Results" />
<p>

<em>Figure 1: <strong>Representative outputs of 3D-Guided
Annotation.</strong> This composite plot displays a curated subset of
the automated profiling suite, featuring macro-scale chromosomal
ideograms and topological overlap analysis of the annotated 3D
interactome.</em>
</p>

</div>

### Module 3: Expression-Aware Refinement

Physical proximity is a structural prerequisite, but not a direct proxy
for active **transcriptional regulation**. This module integrates
quantitative transcriptome data to annotate each loop with
expression-aware functional status. All structural loops are preserved;
the pipeline reclassifies silent anchors (P → eP, G → eG), flags which
loops belong to the high-confidence functional subset
(`Retained_In_Functional_Network`), and exposes `Refinement_Action` for
transparent interpretation.

**Key Parameters:**

- `threshold` & `unit_type`: Defines the quantitative cutoff (e.g.,
  `threshold = 1.0`, `unit_type = "TPM"`) required to consider a gene
  biologically active. This parameter enables downward compatibility
  with various normalization methods.
- `reclassify_by_expression`: When enabled (`TRUE`), **transcriptionally
  silent promoters** are not simply discarded; instead, they are
  biologically reclassified to enhancer-like regulatory elements. This
  correction refines the **regulatory topology** (e.g., reclassifying a
  functionally silent **P-P** loop into a curated **eP-P** loop).

``` r
expr_path <- system.file("extdata", "example_tpm.txt", package = "looplook")

refined_res <- refine_loop_anchors_by_expression(
  annotation_res = res_integrated,
  expr_matrix_file = expr_path,
  sample_columns = c("con1", "con2"),
  threshold = 1,
  unit_type = "TPM",
  reclassify_by_expression = TRUE,
  out_dir = out_dir,
  project_name = "Refined_Network"
)
```

**Output Data Dictionary: The Functionally Active Regulome** Following
transcriptome integration, the refined tabular outputs represent a
high-confidence, functionally active subset of the initial 3D chromatin
interactome. The **key features** of these output are as follows:

- **Expression-Aware Refinement**: All structural loops are preserved.
  The pipeline annotates each loop with expression-aware functional
  status (`Has_Active_Target`, `Retained_In_Functional_Network`,
  `Refinement_Action`) and provides a high-confidence functional subset
  for downstream analysis, without discarding structural evidence.
- **Dynamic Topological Reclassification**: The topological annotations
  within the `loop_type` are biologically recalibrated. By dynamically
  reclassifying transcriptionally silent promoters (**P**) as
  enhancer-like elements (**eP**), the catalog fundamentally corrects
  the spatial regulatory syntax (e.g., seamlessly transforming a
  transcriptionally silent **P-P** loop into a curated **eP-P**
  interaction axis).
- **Filtered Target Gene Links**: The refined provenance sheet retains
  only peak-gene links still used by the refined target columns and
  appends `Mean_Expression` plus `Passes_Expression_Filter`, so inactive
  Basic-stage links are not carried forward as active refined evidence.

<div align="center">

<img src="man/figures/g2_refine.jpg" width="800" style="border: 1px solid #ddd; border-radius: 4px; padding: 5px;" alt="Refinement Results" />
<p>

<em>Figure 2: <strong>Representative outputs of Expression-Aware
Refinement.</strong> As shown by the Multi-Omics Sankey Tracking, this
curated visualization dynamically traces the fate of peak through 3D
chromatin topological interactions to their corresponding
transcriptionally active target genes.</em>
</p>

</div>

### Module 4: Automated Functional Profiling

This module provides a fully automated, end-to-end multi-omics analysis
pipeline that integrates 3D genomic interactions with transcriptomic
data to unveil the regulatory mechanisms of targets of interest.

**Key Parameters:**

- **`target_source` (The Biological Scope)** Defines the biological
  scope of functional profiling.
  - **`"targets"`**: Focuses exclusively on the putative target genes
    regulated by inputted genomic features (**Peak-Centric mode**).
  - **`"loops"`**: Evaluates the entire 3D interactome independently of
    the inputted genomic features (**Global Network-Centric mode**).
- **`target_mapping_mode`** Controls the stringency of 3D target
  assignment. **`"all"`** accepts broad 3D target regulation, while
  **`"promoter"`** is highly stringent, requiring 3D loops to explicitly
  anchor at canonical promoter regions, excluding distal E-G
  connections.
- **`include_Filled` (The Stringency Toggle)** Adjusts the stringency of
  annotation integration.
  - **`TRUE` (Hybrid Mode)**: Utilizes the comprehensively merged
    annotation, prioritizing 3D loop-derived target genes while
    **rescuing unlooped genomic elements** by assigning them to their
    nearest linear genes.
  - **`FALSE` (Pure Spatial Mode)**: Strictly isolates and analyzes only
    the **3D interactome**.
- **`use_nearest_gene` (The Control)** Serves as a **classical baseline
  reference**. If set to `TRUE`, the engine bypasses 3D spatial topology
  and strictly assigns genomic features to their nearest linear genes,
  facilitating a direct comparison to demonstrate the **novel functional
  insights** gained from 3D-guided mapping.

``` r
diff_path <- system.file("extdata", "example_deg.txt", package = "looplook")
meta_path <- system.file("extdata", "example_coldata.txt", package = "looplook")

res_profile <- profile_target_genes(
  annotation_res = refined_res,
  diff_file = diff_path,
  lfc_col = "log2FoldChange",
  expr_matrix_file = expr_path,
  metadata_file = meta_path,
  target_source = c("loops", "targets"),
  target_mapping_mode = "all",
  include_Filled = TRUE,
  use_nearest_gene = FALSE,
  project_name = "Functional_Profiling",
  run_motif = FALSE,
  run_go = FALSE,
  run_ppi = FALSE
)
```

<div align="center">

<img src="man/figures/g3_profile.jpg" width="800" style="border: 1px solid #ddd; border-radius: 4px; padding: 5px;" alt="Profiling Results" />
<p>

<em>Figure 3: <strong>Representative outputs of Functional
Profiling.</strong> This curated composite highlights Divergent Concept
Networks and Asymmetric Motif Signatures, offering a partial glimpse
into the downstream visualizations that decode the trans-regulatory
logic of the spatial hubs.</em>
</p>

</div>

### Module 5: IGV-Style Track Visualization

This module precisely renders the **local** 3D chromatin spatial
interactome through a multi-tiered genomic browser-style visualization
interface, similar to the layout of the `Integrative Genomics Viewer`
(IGV).

**Key Parameters:**

- **`score_to_alpha`** Logical parameter. If `TRUE`, it maps
  quantitative **chromatin interaction scores** to the alpha
  (transparency) channel of the Bezier arcs, enabling visual
  differentiation of interaction strength.
- **`species`** Specifies the **organism of interest**, directing the
  function to automatically load the corresponding `TxDb` and `OrgDb`
  Bioconductor packages. This ensures **precise rendering of gene
  tracks**, including exon-intron structures and strand directionality.

``` r
bedpe_path <- system.file("extdata", "example_loops_1.bedpe", package = "looplook")
bed_path <- system.file("extdata", "example_peaks.bed", package = "looplook")

track_plot <- plot_peaks_interactions(
  bedpe_file = bedpe_path,
  target_bed = bed_path,
  chr = "chr1",
  from = 11884299,
  to = 12106581,
  species = "hg38"
)
```

<div align="center">

<img src="man/figures/plot1.jpg" width="800" style="border: 1px solid #ddd; border-radius: 4px; padding: 5px;" alt="Track Plot Results" />
<p>

<em>Figure 4: Integrative genomic browser view displaying 3D chromatin
loops, genomic regions, and directional gene models.</em>
</p>

</div>

------------------------------------------------------------------------

### One-Click Parameterised Report

To facilitate intuitive data exploration and result interpretation, an
integrated HTML report can be compiled to encapsulate the complete
analytical workflow (annotation, refinement, and profiling), presenting
the outputs in a structured and accessible format.

``` r
# One-click report (renders a standalone HTML via nested rmarkdown)
looplook::looplook_report(
  bedpe_file = system.file("extdata", "example_loops_1.bedpe", package = "looplook"),
  target_bed = system.file("extdata", "example_peaks.bed", package = "looplook"),
  expr_matrix_file = system.file("extdata", "example_tpm.txt", package = "looplook"),
  diff_file = system.file("extdata", "example_deg.txt", package = "looplook"),
  metadata_file = system.file("extdata", "example_coldata.txt", package = "looplook"),
  project_name = "My HiChIP Study"
)
```

The report is also available from the RStudio menu: **File → New File →
R Markdown → From Template → looplook Report**.

------------------------------------------------------------------------

## Contact

**Ying ZHANG** Zhejiang University Email: <12207129@zju.edu.cn>

For bug reports, feature requests, or questions regarding the package,
please open an issue at the [looplook GitHub
repository](https://github.com/zying106/looplook/issues).

------------------------------------------------------------------------

## Resources

- **Package website:** <https://zying106.github.io/looplook/> — full
  documentation, function reference, and vignettes
- **Preprint:** [bioRxiv
  10.64898/2026.04.03.715516](https://www.biorxiv.org/content/10.64898/2026.04.03.715516v1)
  — package manuscript and benchmarks

------------------------------------------------------------------------

## Citation

If you use `looplook` in your research, please cite the preprint:

> ZHANG Y, HUANG X, CHEN Y, XU L. **looplook: An integrative suite for
> expression-aware target assignment and functional annotation of
> chromatin interactions.** *bioRxiv*, 2026. DOI:
> [10.64898/2026.04.03.715516](https://www.biorxiv.org/content/10.64898/2026.04.03.715516v1)

``` bibtex
@article{zhang2026looplook,
  title = {looplook: An integrative suite for expression-aware target assignment and functional annotation of chromatin interactions},
  author = {Zhang, Ying and Huang, Xingze and Chen, Ye and Xu, Liang},
  year = {2026},
  journal = {bioRxiv},
  doi = {10.64898/2026.04.03.715516}
}
```

------------------------------------------------------------------------

## Session Information

For reproducibility, `looplook` is developed and tested under the
following environment:

    ## R version 4.5.1 (2025-06-13)
    ## Platform: aarch64-apple-darwin20
    ## Running under: macOS Sequoia 15.3
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: Asia/Singapore
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] looplook_0.99.13
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] splines_4.5.1                           
    ##   [2] BiocIO_1.18.0                           
    ##   [3] bitops_1.0-9                            
    ##   [4] ggplotify_0.1.3                         
    ##   [5] fields_17.1                             
    ##   [6] tibble_3.3.0                            
    ##   [7] R.oo_1.27.1                             
    ##   [8] polyclip_1.10-7                         
    ##   [9] XML_3.99-0.19                           
    ##  [10] rpart_4.1.24                            
    ##  [11] karyoploteR_1.34.2                      
    ##  [12] lifecycle_1.0.4                         
    ##  [13] rstatix_0.7.3                           
    ##  [14] doParallel_1.0.17                       
    ##  [15] lattice_0.22-7                          
    ##  [16] ensembldb_2.32.0                        
    ##  [17] MASS_7.3-65                             
    ##  [18] ggdist_3.3.3                            
    ##  [19] backports_1.5.0                         
    ##  [20] magrittr_2.0.4                          
    ##  [21] openxlsx_4.2.8.1                        
    ##  [22] Hmisc_5.2-4                             
    ##  [23] rmarkdown_2.30                          
    ##  [24] yaml_2.3.10                             
    ##  [25] ggtangle_0.0.7                          
    ##  [26] otel_0.2.0                              
    ##  [27] spam_2.11-1                             
    ##  [28] zip_2.3.3                               
    ##  [29] cowplot_1.2.0                           
    ##  [30] DBI_1.2.3                               
    ##  [31] RColorBrewer_1.1-3                      
    ##  [32] maps_3.4.3                              
    ##  [33] abind_1.4-8                             
    ##  [34] zlibbioc_1.54.0                         
    ##  [35] GenomicRanges_1.60.0                    
    ##  [36] purrr_1.2.1                             
    ##  [37] R.utils_2.13.0                          
    ##  [38] AnnotationFilter_1.32.0                 
    ##  [39] biovizBase_1.56.0                       
    ##  [40] BiocGenerics_0.54.1                     
    ##  [41] RCurl_1.98-1.17                         
    ##  [42] yulab.utils_0.2.1                       
    ##  [43] nnet_7.3-20                             
    ##  [44] VariantAnnotation_1.54.1                
    ##  [45] tweenr_2.0.3                            
    ##  [46] rappdirs_0.3.3                          
    ##  [47] circlize_0.4.16                         
    ##  [48] GenomeInfoDbData_1.2.14                 
    ##  [49] IRanges_2.42.0                          
    ##  [50] S4Vectors_0.48.0                        
    ##  [51] enrichplot_1.28.4                       
    ##  [52] data.tree_1.2.0                         
    ##  [53] ggrepel_0.9.8                           
    ##  [54] tidytree_0.4.6                          
    ##  [55] codetools_0.2-20                        
    ##  [56] DelayedArray_0.34.1                     
    ##  [57] DOSE_4.2.0                              
    ##  [58] ggforce_0.5.0                           
    ##  [59] shape_1.4.6.1                           
    ##  [60] tidyselect_1.2.1                        
    ##  [61] aplot_0.2.9                             
    ##  [62] UCSC.utils_1.4.0                        
    ##  [63] farver_2.1.2                            
    ##  [64] viridis_0.6.5                           
    ##  [65] matrixStats_1.5.0                       
    ##  [66] stats4_4.5.1                            
    ##  [67] base64enc_0.1-3                         
    ##  [68] bamsignals_1.40.0                       
    ##  [69] GenomicAlignments_1.44.0                
    ##  [70] jsonlite_2.0.0                          
    ##  [71] GetoptLong_1.1.0                        
    ##  [72] Formula_1.2-5                           
    ##  [73] iterators_1.0.14                        
    ##  [74] foreach_1.5.2                           
    ##  [75] tools_4.5.1                             
    ##  [76] treeio_1.32.0                           
    ##  [77] Rcpp_1.1.0                              
    ##  [78] glue_1.8.0                              
    ##  [79] gridExtra_2.3                           
    ##  [80] SparseArray_1.8.1                       
    ##  [81] xfun_0.55                               
    ##  [82] distributional_0.6.0                    
    ##  [83] qvalue_2.40.0                           
    ##  [84] MatrixGenerics_1.20.0                   
    ##  [85] GenomeInfoDb_1.44.3                     
    ##  [86] dplyr_1.1.4                             
    ##  [87] withr_3.0.2                             
    ##  [88] fastmap_1.2.0                           
    ##  [89] ggpointdensity_0.2.1                    
    ##  [90] digest_0.6.37                           
    ##  [91] R6_2.6.1                                
    ##  [92] gridGraphics_0.5-1                      
    ##  [93] colorspace_2.1-2                        
    ##  [94] networkD3_0.4.1                         
    ##  [95] GO.db_3.21.0                            
    ##  [96] dichromat_2.0-0.1                       
    ##  [97] RSQLite_2.4.3                           
    ##  [98] R.methodsS3_1.8.2                       
    ##  [99] UpSetR_1.4.0                            
    ## [100] tidyr_1.3.1                             
    ## [101] generics_0.1.4                          
    ## [102] data.table_1.17.8                       
    ## [103] rtracklayer_1.68.0                      
    ## [104] InteractionSet_1.36.1                   
    ## [105] httr_1.4.7                              
    ## [106] htmlwidgets_1.6.4                       
    ## [107] S4Arrays_1.8.1                          
    ## [108] regioneR_1.40.1                         
    ## [109] pkgconfig_2.0.3                         
    ## [110] gtable_0.3.6                            
    ## [111] blob_1.2.4                              
    ## [112] ComplexHeatmap_2.26.0                   
    ## [113] S7_0.2.0                                
    ## [114] XVector_0.48.0                          
    ## [115] clusterProfiler_4.16.0                  
    ## [116] htmltools_0.5.8.1                       
    ## [117] carData_3.0-5                           
    ## [118] dotCall64_1.2                           
    ## [119] fgsea_1.34.2                            
    ## [120] clue_0.3-66                             
    ## [121] ProtGenerics_1.40.0                     
    ## [122] scales_1.4.0                            
    ## [123] Biobase_2.68.0                          
    ## [124] TxDb.Hsapiens.UCSC.hg38.knownGene_3.21.0
    ## [125] png_0.1-9                               
    ## [126] ggfun_0.2.0                             
    ## [127] knitr_1.51                              
    ## [128] rstudioapi_0.17.1                       
    ## [129] reshape2_1.4.4                          
    ## [130] rjson_0.2.23                            
    ## [131] nlme_3.1-168                            
    ## [132] checkmate_2.3.3                         
    ## [133] curl_7.0.0                              
    ## [134] org.Hs.eg.db_3.21.0                     
    ## [135] GlobalOptions_0.1.2                     
    ## [136] cachem_1.1.0                            
    ## [137] stringr_1.5.2                           
    ## [138] parallel_4.5.1                          
    ## [139] foreign_0.8-90                          
    ## [140] AnnotationDbi_1.70.0                    
    ## [141] restfulr_0.0.16                         
    ## [142] pillar_1.11.1                           
    ## [143] grid_4.5.1                              
    ## [144] vctrs_0.7.2                             
    ## [145] ggpubr_0.6.2                            
    ## [146] car_3.1-3                               
    ## [147] cluster_2.1.8.1                         
    ## [148] htmlTable_2.4.3                         
    ## [149] evaluate_1.0.5                          
    ## [150] magick_2.9.0                            
    ## [151] GenomicFeatures_1.60.0                  
    ## [152] cli_3.6.5                               
    ## [153] compiler_4.5.1                          
    ## [154] bezier_1.1.2                            
    ## [155] Rsamtools_2.24.1                        
    ## [156] rlang_1.2.0                             
    ## [157] crayon_1.5.3                            
    ## [158] ggsignif_0.6.4                          
    ## [159] plyr_1.8.9                              
    ## [160] fs_1.6.6                                
    ## [161] stringi_1.8.7                           
    ## [162] viridisLite_0.4.2                       
    ## [163] BiocParallel_1.42.2                     
    ## [164] Biostrings_2.76.0                       
    ## [165] lazyeval_0.2.2                          
    ## [166] GOSemSim_2.34.0                         
    ## [167] Matrix_1.7-3                            
    ## [168] BSgenome_1.76.0                         
    ## [169] patchwork_1.3.2                         
    ## [170] bit64_4.6.0-1                           
    ## [171] ggplot2_4.0.0                           
    ## [172] KEGGREST_1.48.1                         
    ## [173] SummarizedExperiment_1.38.1             
    ## [174] broom_1.0.11                            
    ## [175] igraph_2.2.0                            
    ## [176] memoise_2.0.1                           
    ## [177] ggtree_3.16.3                           
    ## [178] fastmatch_1.1-6                         
    ## [179] bit_4.6.0                               
    ## [180] gson_0.1.0                              
    ## [181] ape_5.8-1
