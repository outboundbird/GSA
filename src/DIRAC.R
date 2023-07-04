#' ---
#' title:  Differential Rank Conservation (DIRAC)
#' subtitle: 'SAR: NA , Study: NA'
#' author:  Siying Huang (E0482362), Biomarker statistics team
#' date: 'created: 2023-07-04 , updated (`r Sys.Date()`)'
#' always_allow_html: true
#' output:
#'   html_document:
#'     css: sanofi.css
#'     code_folding: "hide"
#'     toc: yes
#'     toc_float:
#'       collapse: no
#' bibliography: references.bib
#' ---
#+ setup, include = FALSE
knitr::opts_knit$set(root_dir='/mnt/c/Users/e0482362/Work/pathway_analysis/src')
knitr::opts_chunk$set(echo = T, comment = '',message = F, warning = F, error=F)
options(width = 100)
#+ libs
library(here)

#' #  Differential Rank Conservation (DIRAC)
#' DIRAC was proposed by [Eddy et al ](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2877722/).
#'
#' > Differential Rank Conservation (DIRAC), which permits one to assess these combinatorial interactions to quantify various biological pathways or networks in a comparative sense, and to determine how they change in different individuals experiencing the same disease process.
#' > -- @eddyIdentifyingTightlyRegulated2010
#'
#' ## Overview of DIRAC method
#'
#' ![](images/dirac.jpg)
#'
#' ## DIRAC in R
#'
library(GSReg)

#' ## Reference
/*
#' <details><summary>Session Info</summary>
sessionInfo()
#' </details>
#+ echo = F, eval = F
# Markdown --------------------------------------------------------
# rmarkdown::render('src/DIRAC.R', output_dir = 'output')
# knitr::spin('src/DIRAC.R', format = 'Rmd', knit = FALSE)
*/