#' ---
#' title: GSVA
#' subtitle: 'SAR: sar , Study: study'
#' author:  Siying Huang (E0482362), Biomarker statistics team
#' date: 'created: 2023-06-09 , updated (`r Sys.Date()`)'
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
ctrl_ids <- c("GSM464697", "GSM464699", "GSM464705", "GSM464702", "GSM464715")
case_ids <- c("GSM464678", "GSM464688", "GSM464683", "GSM464692", "GSM464696")

# remove dup probes, select case, control samples, log-transform
expr <- expr[!duplicated(expr$IDENTIFIER), ] %>%
  data.frame() %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames("IDENTIFIER") %>%
  dplyr::select_at(c(ctrl_ids, case_ids)) %>%
  mutate_if(is.numeric, log2)
#' ## Gene set selection
#' Here I select the genes of each set1 in the same manner as before. While for
#' geneset 2 and 3, I choose up and down-regulated genes only for each set.
# same gene sets
set.seed(123)
gs1 <- sample(rst$ID, 10)
gs2 <- head(rst[which(rst$logFC > 20 & rst$adj.P.Val < 0.0001), "ID"], 10)
gs3 <- head(rst[which(rst$logFC < -20 & rst$adj.P.Val < 0.0001), "ID"], 10)

cat("set1:", gs1, "\nset2:", gs2, "\nset3:", gs3)
n_bg <- 50
#' For the purpose of illustration and simplification, I randomly choose another `r n_bg` genes in the expression dataset as background gene.
bg <- sample(rst$ID, n_bg)
gs <- Reduce(union, list(gs1, gs2, gs3, bg))
length(gs)

#' ## Gene set variation calculation
#+ hmap, fig.dim = c(4, 8)
expr_g <- expr[gs, ]
labels <- data.frame(status = rep(c("control", "case"), each = 5))
rownames(labels) <- c(ctrl_ids, case_ids)
grid::grid.newpage()
h <- pheatmap::pheatmap(expr_g,
  annotation_col = labels,
  border_color = NA,
  # breaks = seq(-6, 6, length.out = 100),
  legend_labels = "",
  main = "",
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  clusstering_method = "complete",
  cluster_rows = F, show_rownames = FALSE,
  cluster_cols = F, show_colnames = FALSE,
  silent = TRUE
)
grid::grid.draw(h)

#' ### esimate cumulative density function for each gene in a gene set
#' Here I use the example of gene set 1.
#+ fig.dim = c(8, 4)
par(mfrow = c(1,2), cex = 0.8)
plot(density(as.numeric(expr_g[gs1[1], ])),
  ylim = c(0, 2.5), xlim = c(0, 12),
  main = "Gene expression levels across samples"
)
invisible(lapply(seq(2, length(gs1)), function(i) {
  g <- gs1[i]
  lines(density(as.numeric(expr_g[g, ])),
    col = i,
    lty = i
  )
}))
legend("topleft", gs1, col = 1:10, lty = 1:10, cex = 0.6, bty = "n")

# cdf before normaliization
plot(ecdf(expr_g[gs1[1], ]),
  verticals = T, do.points = F, xlim = c(0, 12),
  main = "Culmulative density of gene expression levels"
)
invisible(lapply(seq(2, length(gs1)), function(i) {
  g <- gs1[i]
  plot(ecdf(expr_g[g, ]),
    verticals = T, do.points = F,
    col = i,
    lty = i,
    add = T
  )
}))

legend("topleft", gs1, col = 1:10, lty = 1:10, cex = 0.6, bty = "n")

#' ## normalizing expr levels
#' normalizing expression levels with culmulative density function Fx
gd <- apply(expr_g, 1, function(x) {
  f <- ecdf(x)
  f(x)
}) %>%
  t() %>%
  data.frame() %>%
  setNames(colnames(expr_g))

# gd <- apply(gd, 1, function(x) ifelse(x == 1, 0.99, x)) %>% t()
# normalize gene density
gene_density <- log(gd / (1 - gd))

#' Here I plot the emprical CDF instead of the explicit CDF mentioned in the paper for the purpose of illustration.
#'
#+ cdf , echo = F, fig.dim = c(8, 4)
par(mfrow = c(1, 2))
plot(density(as.numeric(gene_density[1, ])))
invisible(lapply(seq(2, length(gs1)), function(i) {
  g <- gs1[i]
  lines(density(as.numeric(gene_density[g, ])),
    col = i,
    lty = i
  )
}))
legend("topleft", gs1, col = 1:10, lty = 1:10, cex = 0.6, bty = "n")

plot(ecdf(as.numeric(gene_density[1, ])),
  verticals = T, do.points = F,
  xlim = c(-10, 10),
  ylim = c(0,1),
  main = "Culmulative density of gene expression levels after normalization"
)

invisible(lapply(seq(2, length(gs1)), function(i) {
  plot(ecdf(as.numeric(gene_density[i, ])),
    verticals = T, do.points = F,
    col = i,
    lty = i,
    add = T
  )
}))

legend("bottomright", gs1, col = 1:10, lty = 1:10, cex = 0.6, bty = "n")

#' ## Sample-wise enrichment scores
#' ### compute rank score
#' To reduce the influence of potential outliers, we first convert $z_{ij}$ to ranks $z(i)j$ for each sample j and normalize further $r_{ij}=|p/2−z_{(i)j}|$ to make the ranks symmetric around zero.
#'
sort_sgn_idxs <- apply(gene_density, 2, order, decreasing = TRUE)
rownames(sort_sgn_idxs) <- rownames(expr_g)
head(sort_sgn_idxs)

rank_norm <- apply(sort_sgn_idxs, 2, compute_rank_score)
rownames(rank_norm) <- rownames(sort_sgn_idxs)

#' ## enrichemnent statistics -GSVA score
#' here i set tau = 1
head(rank_norm[, 1:3])
n_tot = nrow(expr_g)
n_gs = length(gs1)

sample <- rank_norm[, 1, drop = F] %>%
  data.frame() %>%
  tibble::rownames_to_column("ID") %>%
  mutate(gs1 = ID %in% gs1) %>%
  rename(rscore = GSM464697)
sum_gs = sum(sample[sample$gs1, "rscore"])
dnom = 1/ (n_tot - n_gs)

sample <- sample %>%
  mutate(score = if_else(gs1, abs(rscore) / sum_gs - as.numeric(gs1) * dnom, 0),
  score_all =  if_else(!gs1, abs(rscore) / sum_gs - as.numeric(gs1) * dnom, 0)) %>%
  mutate(
    score = case_when(row_number() == 1 ~ 0, TRUE ~ score),
    score_all = case_when(row_number() == 1 ~ 0, TRUE ~ score_all),
    es = cumsum(score),
    es_all = cumsum(score_all)
  )

plot(sample$es, type = "s")
lines(sample$es_all, type = "s")

# calculate random walk stat

gset_idxs <- which(rownames(gene_density) %in% gs1)
gset_idxs <- gs1
ks_test_Rcode(data.frame(gene_density[, "GSM464697", drop = F]), gset_idxs)

gene.density <- gene_density[, 1, drop = F]
ks_test_Rcode <- function(gene.density, gset_idxs, tau = 1, make.plot = FALSE) {
  n.genes <- length(gene.density)
  n.gset <- length(gset_idxs)

  sum.gset <- sum(abs(gene.density[gset_idxs])^tau)

  dec <- 1 / (n.genes - n.gset)

  sort.idxs <- order(gene.density, decreasing = T)
  offsets <- sort(match(gset_idxs, sort.idxs))

  last.idx <- 0
  values <- rep(NaN, length(gset_idxs))
  current <- 0
  for (i in seq_along(offsets)) {
    current <- current + abs(gene.density[sort.idxs[offsets[i]]])^tau / sum.gset - dec * (offsets[i] - last.idx - 1)

    values[i] <- current
    last.idx <- offsets[i]
  }
  check_zero <- current - dec * (n.genes - last.idx)
  # if(check_zero > 10^-15){
  # 	stop(paste=c("Expected zero sum for ks:", check_zero))
  # }
  if (make.plot) {
    plot(offsets, values, type = "l")
  }

  max.idx <- order(abs(values), decreasing = T)[1]
  mx.value <- values[max.idx]

  return(mx.value)
}

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