#' ---
#' title: Gene set variation analysis
#' subtitle: 'SAR: NA , Study: NA'
#' author:  Siying Huang (E0482362), Biomarker statistics team
#' date: 'created: 2023-06-06 , updated (`r Sys.Date()`)'
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

#' <details><summary>Session Info</summary>
sessionInfo()
#' </details>
#+ echo = F, eval = F
# Markdown --------------------------------------------------------
# rmarkdown::render('src/07_GSVA.R', output_dir = 'output')
# knitr::spin('src/07_GSVA.R', format = 'Rmd', knit = FALSE)