#' ---
#' title: Intro
#' subtitle: 'SAR: sar , Study: study'
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

knitr::opts_chunk$set(echo = T, comment = '',message = F, warning = F, error=F)
options(width = 100)
#+ echo = F

#' # Gene sets and pathways
#'
#' ::: {.definition #label name="definition name"}
#' gene set , biological pathway
#' :::
#'
#'
#' ## Biological process data bases

#' <details><summary>Session Info</summary>
sessionInfo()
#' </details>
#+ echo = F, eval = F
# Markdown --------------------------------------------------------
# rmarkdown::render('src/01_intro.R', output_dir = 'output')
# knitr::spin('src/01_intro.R', format = 'Rmd', knit = FALSE)