#' ---
#' title: title
#' subtitle: 'SAR: sar , Study: study'
#' author:  Siying Huang (E0482362), Biomarker statistics team
#' date: 'created: 2023-06-08 , updated (`r Sys.Date()`)'
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
library(dplyr)
library(ggplot2)
rst <- readRDS(file.path(here(), "data/de_rst.rds"))
gene_ids <- readRDS(file.path(here(), "data/gene_id.rds"))
#' # Other methods
#' - [p value combination](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3134237/)
#' - [maxmean statsitics, restandardization](https://arxiv.org/pdf/math/0610667.pdf)
#' @efronTestingSignificanceSets GSA package
#' - [quantile method]()ttps://www.pnas.org/doi/10.1073/pnas.0506577102)
#' - PAGE
set.seed(123)
gs1 <- sample(rst$ID, 30)
gs2 <- c(head(rst$ID, 15), tail(rst$ID, 15))
gs3 <- head(rst$ID, 30)

rst <- rst %>%
  mutate(
    gs1 = ID %in% gs1,
    gs2 = ID %in% gs2,
    gs3 = ID %in% gs3
  )

# compare density plot

plot(density(rst$t),
  ylim = c(0, 1),
  lwd = 4,
  col = "grey"
)
lines(density(rst[rst$gs1, "t"]),
  ylim = c(0, 1), col = "red"
)
lines(density(rst[rst$gs2, "t"]),
  ylim = c(0, 1), col = "darkgreen"
)
lines(density(rst[rst$gs3, "t"]),
  ylim = c(0, 1),
  col = "blue"
)
legend("topleft",
  c("genome", "gs1", "gs2", "gs3"),
  col = c("grey", "red", "darkgreen", "blue"),
  lty = 1,
  bty = "n"
)

#' # Comparison of gene set analysis methods
#' https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9882404/
#' ## Reviews
#' https://www.frontiersin.org/articles/10.3389/fphys.2015.00383/full
#'
#' <details><summary>Session Info</summary>
sessionInfo()
#' </details>
#+ echo = F, eval = F
# Markdown --------------------------------------------------------
# rmarkdown::render('src/07_others.R', output_dir = 'output')
# knitr::spin('src/07_others.R', format = 'Rmd', knit = FALSE)