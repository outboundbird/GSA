#' ---
#' title: Over-representation Analysis
#' subtitle: 'SAR: NA , Study: NA'
#' author:  Siying Huang (E0482362), Biomarker statistics team
#' date: 'created: 2023-06-01 , updated (`r Sys.Date()`)'
#' always_allow_html: true
#' output:
#'   html_document:
#'     css: sanofi.css
#'     code_folding: "show"
#'     toc: yes
#'     toc_float:
#'       collapse: no
#' ---
#+ setup, include = FALSE
knitr::opts_chunk$set(echo = T, comment = "", message = F, warning = F, error = F)
options(width = 100)
#+ libs
library(here)
#' # Experiment Data
#'
#' This pathway analysis will be using the smoking study data GDS3713 (GSE18723).
#' The results from DEA will be used for ORA.
#+ io
rst <- readRDS(file.path(here(), "data/de_rst.rds"))

#' In this micro-array data due to that there are multiple probes match to
#' the same gene. For illustration of the PA purpose I simply removed the
#' duplicated gene probes.In pratice these duplicates should be removed with cautions.

rm_idx <- duplicated(rst$ID)
rst <- rst[!rm_idx, ]
#' This leaves `r nrow(rst)` genes in the dataset.
#'
#' # Gene set selection
#'  Normally these gene sets should represent major players in biological
#' pathways. This can be achieved through the query from the biological
#' databases or from other analytical approach such as co-expression analysis.
#' For illustration purpose, here I randomly select 10 genes for each gene set.

set.seed(123)
gs1 <- sample(rst$ID, 10)
gs2 <- sample(rst$ID, 10)

#' Here I choose five most significantly differentially expressed genes,
#' five genes not differentially expressed.
gs3 <- c(head(rst$ID, 5), tail(rst$ID, 5))


#' ## define differentially expressed gene
#'
#' Here I define differentially expressed genes with following criterion:
#' - adjusted p-value < 0.001
#' - log fold change > 5
#'
#' These thresholds are not fixed can be varied depending on the context.
#' Or it can be obtained from other analytical methods.
#'
rst <- rst %>%
  mutate(
    de = adj.P.Val < 0.001 & abs(logFC) > 5,
    gs1 = ID %in% gs1,
    gs2 = ID %in% gs2,
    gs3 = ID %in% gs3
  )

table(rst[, c("gs1", "de")])
table(rst[, c("gs2", "de")])
table(rst[, c("gs3", "de")])
#' ## Over-representation test
#' **The null hypothesis of ORA test:** being differentally expressed and
#' belong to a particular pathway (gene set) are independent.
#'
#' **Common statistical tests: ** chi-square test, Fisher's exact test. Note
#' that when counts in a cell of a contigency table is small, chi-square
#' approximation may not be accurate.
#'
#' ### Gene set 1

conttab_gs1 <- as.table(table(rst[, c("gs1", "de")]))
chisq.test(conttab_gs1)
fisher.test(conttab_gs1)
#' ### Gene set 3
conttab_gs3 <- as.table(table(rst[, c("gs3", "de")]))
chisq.test(conttab_gs3)
fisher.test(conttab_gs3)

#' <details><summary>Session Info</summary>
sessionInfo()
#' </details>
#+ echo = F, eval = F
# Markdown --------------------------------------------------------
# rmarkdown::render('src/ORA.R', output_dir = 'output')