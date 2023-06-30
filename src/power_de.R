#' ---
#' title: Power analysis for RNAseq differential expression analysis
#' subtitle: 'SAR: NA , Study: NA'
#' author:  Siying Huang (E0482362), Biomarker statistics team
#' date: 'created: 2023-06-30 , updated (`r Sys.Date()`)'
#' always_allow_html: true
#' output:
#'   html_document:
#'     css: sanofi.css
#'     code_folding: "hide"
#'     toc: yes
#'     toc_float:
#'       collapse: no
#' ---
#+ setup, include = FALSE
knitr::opts_knit$set(root_dir = "/mnt/c/Users/e0482362/Work/pathway_analysis/src")
knitr::opts_chunk$set(echo = T, comment = "", message = F, warning = F, error = F)
options(width = 100)
#+ libs
library(here)

#' # Suggested papers
#' - [Statistical Power Analysis for Designing Bulk, Single-Cell, and Spatial Transcriptomics Experiments: Review, Tutorial, and Perspectives](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9952882/)
#' - [Power analysis and sample size estimation for RNA-Seq differential expression](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4201821/)
#' - [Power analysis for RNA-Seq differential expression studies]()
#' - [powsimR: power analysis for bulk and single cell RNA-seq experiments](https://academic.oup.com/bioinformatics/article/33/21/3486/3952669)
#' - [Dedicated transcriptomics combined with power analysis lead to functional understanding of genes with weak phenotypic changes in knockout lines](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008354)
#' - [ssizeRNA](https://cran.r-project.org/web/packages/ssizeRNA/vignettes/ssizeRNA.pdf)

#' <details><summary>Session Info</summary>
sessionInfo()
#' </details>
#+ echo = F, eval = F
# Markdown --------------------------------------------------------
# rmarkdown::render('src/power_de.R', output_dir = 'output')
# knitr::spin('src/power_de.R', format = 'Rmd', knit = FALSE)