#' ---
#' title: Differential coexpression analysis
#' subtitle: 'SAR: NA , Study: tutorial'
#' author:  Siying Huang (E0482362), Biomarker statistics team
#' date: 'created: 2023-06-26 , updated (`r Sys.Date()`)'
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
knitr::opts_knit$set(root_dir = "/mnt/c/Users/e0482362/Work/pathway_analysis/src")
knitr::opts_chunk$set(echo = T, comment = "", message = F, warning = F, error = F)
options(width = 100)
#+ libs
library(here)
#' Gene set coexpression analysis and gene set net correlation analysis are classified as
#' topology-based methods. This type of analysis addresses the issue that not all
#' genes play equally important role in a biological process. Therefore it considers
#' the correlations among genes.
#'
#' # Gene Set Coexpression Analysis (GSCA)
#' Gene set coexpression analysis was first proposed by @choiStatisticalMethodsGene2009a.
#'
#'  ## Scheme
#' ![](images/gsca.png)
#'
#' # Gene Sets Net Correlation Analysis (GSNCA)
#' GSNCA method was proposed by @rahmatallahGeneSetsNet2014.
#'
#' ## Scheme
#' ![](images/gsnca.jpg)
#'
#' > in the gene coexpression network between two conditions. Importantly, we do not infer ‘gene coexpression networks’ explicitly, but, instead, we estimate net correlation changes by introducing for each gene a weight factor that characterizes its cross-correlations in the coexpression networks. Weight vectors in both conditions are found as eigenvectors of correlation matrices with zero diagonal elements. The Gene Sets Net Correlations Analysis (GSNCA) tests the hypothesis that for a gene set there is no difference in the gene weight vectors between two conditions.
#'
#' # Differential gene set coexpression analysis in R
#'
library(GSAR)

#' ## Reference
/*
#' <details><summary>Session Info</summary>
sessionInfo()
#' </details>
#+ echo = F, eval = F
# Markdown --------------------------------------------------------
# rmarkdown::render('src/DCA.R', output_dir = 'output')
# knitr::spin('src/DCA.R', format = 'Rmd', knit = FALSE)
*/