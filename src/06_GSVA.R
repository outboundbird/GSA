#' ---
#' title: GSVA
#' subtitle: 'SAR: sar , Study: study'
#' author:  Siying Huang (E0482362), Biomarker statistics team
#' date: 'created: 2023-06-09 , updated (`r Sys.Date()`)'
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
knitr::opts_chunk$set(echo = T, comment = "", message = F, warning = F, error = F)
options(width = 100)
#+ libs
library(here)
library(dplyr)
library(pheatmap)
library(GSVA)
#' # Gene set variation analysis (GSVA)
#' The gene set variation analysis (GSVA) was developed by Hänzelmann et al in 2013, which is a a non-parametric, unsupervised method. @hanzelmannGSVAGeneSet2013
#' 
#' > GSVA calculates sample-wise gene set enrichment scores as a function of
#' > genes inside and outside the gene set, analogously to a competitive gene
#' > set test. Further, it estimates variation of gene set enrichment over
#' >  the samples independently of any class label.
#'
#' ![](/mnt/c/Users/e0482362/Work/pathway_analysis/figure/gsva.png)
#'

#+ cache = T
expr <- readRDS(file.path(here(), "data/expr.rds"))
rst <- readRDS(file.path(here(), "data/de_rst.rds"))
# pheno <- readRDS(file.path(here(),'data/pheno.rds'))
# data management

# pick 3 smokers , 3 non-smokers
ctrl_ids <- c("GSM464699", "GSM464705", "GSM464715")
case_ids <- c("GSM464678", "GSM464683", "GSM464692")

# remove dup probes, select case, control samples, log-transform
expr <- expr[!duplicated(expr$IDENTIFIER), ] %>%
  data.frame() %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames("IDENTIFIER") %>%
  dplyr::select_at(c(ctrl_ids, case_ids)) %>%
  mutate_if(is.numeric, log2)

# same gene sets
set.seed(123)
gs1 <- sample(rst$ID, 10)
gs2 <- c(head(rst$ID, 5), tail(rst$ID, 5))
gs3 <- head(rst$ID, 10)
cat("set1:", gs1, "\nset2:", gs2, "\nset3:", gs3)
n_bg = 50
# to simplify, I randomly choose another `r n_bg` genes in the expression dataset as background gene.
bg <- sample(rst$ID, n_bg)
gs <- Reduce(union, list(gs1, gs2, gs3, bg))
length(gs)
# GSVA
expr_g <- expr[gs, ]

# pheatmap::pheatmap(expr_g)
# esimate cumulative density function for each gene
gd <- apply(expr_g, 1, function(x) {
  f <- ecdf(x)
  f(x)
}) %>%
  t() %>%
  data.frame() %>%
  setNames(colnames(expr_g))
head(gd)
gene_density <- log(gd / (1 - gd))
num_genes <- length(gs)

sort_sgn_idxs <- apply(gene_density, 2, order, decreasing = TRUE)
rownames(sort_sgn_idxs) <- rownames(expr_g)
head(sort_sgn_idxs)
# for sample GSM464699:
sample1 <- sort_sgn_idxs[, "GSM464699", drop = F] %>%
  data.frame() %>%
  dplyr::arrange(desc(GSM464699)) %>%
  tibble::rownames_to_column('ID') %>%
  mutate(gs1 = ID %in% gs1) %>%
    mutate(
      step_gs1 = if_else(gs1, GSM464699, 0L),
      es_gs1 = cumsum(step_gs1),
      step = if_else(gs1, 0L, GSM464699),
      es = cumsum(step)
    )

head(sample1)
ymax <- max(sample1[, c("es_gs1", "es")])
plot(sample1$es_gs1,
  type = "s",
  ylim = c(0, ymax),
  col = "red"
)
lines(sample1$es,
  type = "s"
)
# org -----------------------------------------
compute_rank_score <- function(sort_idx_vec) {
  tmp <- rep(0, num_genes)
  tmp[sort_idx_vec] <- abs(seq(from = num_genes, to = 1) - num_genes / 2)
  return(tmp)
}

rank.scores <- rep(0, num_genes)
rank.scores <- apply(sort_sgn_idxs, 2, compute_rank_score)
head(rank.scores)
#' ## GSVA using R package

/*
#' <details><summary>Session Info</summary>
sessionInfo()
#' </details>
#+ echo = F, eval = F
# Markdown --------------------------------------------------------
# rmarkdown::render('src/06_GSVA.R', output_dir = 'output')
# knitr::spin('src/06_GSVA.R', format = 'Rmd', knit = FALSE)
*/