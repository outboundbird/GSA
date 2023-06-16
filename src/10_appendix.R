/*
#' ---
#' title: Appendix
#' subtitle: 'SAR: sar , Study: study'
#' author:  Siying Huang (E0482362), Biomarker statistics team
#' date: 'created: 2023-06-16 , updated (`r Sys.Date()`)'
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
*/

#' # Appendix {-}
#+ echo = F, cache = T
library(readxl)
library(DT)
library(here)
ora <- read_xlsx(file.path(here(), "data/pa_method.xlsx"), sheet = 1)
fcs <- read_xlsx(file.path(here(), "data/pa_method.xlsx"), sheet = 2)
ptb <- read_xlsx(file.path(here(), "data/pa_method.xlsx"), sheet = 3)
multi <- read_xlsx(file.path(here(), "data/pa_method.xlsx"), sheet = 4)
path_db <- read.csv(file.path(here::here(), "data/public_pathway_db.csv"))

#' ## Biological process data bases
#'
#+ echo = F

datatable(path_db,
  rownames = F,
  options = list(scrollX = "400px")
)
#' ## Summary of currently available pathway analysis methods
#' ### Over-representation analysis
#+ echo = F

datatable(ora)
#' ### Gene set enrichment analysis
#+ echo = F
datatable(fcs)

/*#' <details><summary>Session Info</summary>
sessionInfo()
#' </details>
#+ echo = F, eval = F
# Markdown --------------------------------------------------------
# rmarkdown::render('src/10_appendix.R', output_dir = 'output')
# knitr::spin('src/10_appendix.R', format = 'Rmd', knit = FALSE)
*/