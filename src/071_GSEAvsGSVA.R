#' ---
#' title: GSEA vs. GSVA
#' subtitle: 'SAR: sar , Study: study'
#' author:  Siying Huang (E0482362), Biomarker statistics team
#' date: 'created: 2023-07-07 , updated (`r Sys.Date()`)'
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
knitr::opts_knit$set(root_dir='/mnt/c/Users/e0482362/Work/pathway_analysis/src')
knitr::opts_chunk$set(echo = T, comment = '',message = F, warning = F, error=F)
options(width = 100)
#+ libs
library(here)

#' # Understanding the difference of GSEA vs. GSVA
#'
# find a dataset with batch effect
# GSVA variability can come from different sources, i.e. biological variation, techinical variation.

/*#' <details><summary>Session Info</summary>
sessionInfo()
#' </details>
#+ echo = F, eval = F
# Markdown --------------------------------------------------------
# rmarkdown::render('src/071_GSEAvsGSVA.R', output_dir = 'output')
# knitr::spin('src/071_GSEAvsGSVA.R', format = 'Rmd', knit = FALSE)
*/