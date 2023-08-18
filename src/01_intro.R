/*#' ---
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
#+ setup, include = FALSE*/
knitr::opts_chunk$set(echo = T, comment = '',message = F, warning = F, error=F)
options(width = 100)

#' # Gene sets and pathways
#'
#' ::: {.definition #label name="definition name"}
#' gene set , biological pathway
#' :::
#' The purpose of gene set analysis is to circumvent the limiations from single-gene
#' analysis, such as difficulty of interpreting the biological meaning of single genes,
#' understanding the multiple-comparison adjusted results, poor reproducibility
#' from independent studies.
#'
#' ## Overal of the gene set analysis
#' @deleeuwStatisticalPropertiesGeneset2016
#' @goemanAnalyzingGeneExpression2007
#'
#' ## Biological vs. statistical categorizations of gene set analysis
#'
#' ### Satatistical categorization of GSA
#'
#' The statistical categorization of GSA methods is based on the type of null
#' hypothesis: self-contained vs. competetive null hypothesis.
#'
#' 
#'
#'
#' ## Organization of this e-book
#' I will illustrate the gene pathway analysis using a public dataset for smoking study.
#' Chapter 2 contains information on this study and how to obtain the dataset
#' from NCBI website to an R environment.
#' Chapter 3 quickly illustrate the differential gene expression analysis. This
#' result will be used in some of the parametrical pathway analysis method.
#' Chapter 4 introduce over-representation analysis
#' Chapter 5 introduce the classic gene set enrichment analysis
#' Chapter 6 talks abou the gene set variation analysis
#' Chapter 7 discusses other methods for gene set enrichment analysis
#' Chapter 8 talks about advanced topic regarding gene set analysis
/*
#' <details><summary>Session Info</summary>
sessionInfo()
#' </details>
#+ echo = F, eval = F
# Markdown --------------------------------------------------------
# rmarkdown::render('src/01_intro.R', output_dir = 'output')
# knitr::spin('src/01_intro.R', format = 'Rmd', knit = FALSE)
*/