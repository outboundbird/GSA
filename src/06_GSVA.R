#' ---
#' title: GSVA
#' subtitle: 'SAR: NA , Study: NA'
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
source(file.path(here(), "src/utils.R"))
#' # Gene set variation analysis (GSVA)
#' The gene set variation analysis (GSVA) was developed by [*Hänzelmann et al*](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3618321/#sec-4title) in 2013, which is a non-parametric, unsupervised method for gene set enrichement analysis. @hanzelmannGSVAGeneSet2013 ^[All quoted text are from original paper.]
#'
#' > GSVA calculates sample-wise gene set enrichment scores as a function of
#' > genes inside and outside the gene set, analogously to a competitive gene
#' > set test. Further, it estimates variation of gene set enrichment over
#' >  the samples independently of any class label.
#'
#' ![](/mnt/c/Users/e0482362/Work/pathway_analysis/figure/gsva.png)
#' The major steps for GSVA method:
#'
#' - rank and normalize gene expression levels across samples with CDF function.
#' - rank the gene expression levels within a sample
#' - calculate enrichment score (ES) and KS test
#' - calculate enrichment statistics (GSVA score)
#'
#+ cache = T, echo = F
expr <- readRDS(file.path(here(), "data/expr.rds"))
rst <- readRDS(file.path(here(), "data/de_rst.rds"))
n_bg <- 100
#' ## Smoking study
#'
#' For simlipcity, I choose 5 smokers and 5 non-smokers out of the original
#' dataset. And for the background gene set, I randomly choose `r n_bg` genes.
#'
#+ echo = F
# data management
ctrl_ids <- c("GSM464697", "GSM464699", "GSM464705", "GSM464702", "GSM464715")
case_ids <- c("GSM464678", "GSM464688", "GSM464683", "GSM464692", "GSM464696")

# remove dup probes, select case, control samples, log-transform
expr <- expr[!duplicated(expr$IDENTIFIER), ] %>%
  data.frame() %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames("IDENTIFIER") %>%
  dplyr::select_at(c(ctrl_ids, case_ids)) %>%
  mutate_if(is.numeric, log2)
#' ### Gene set selection
#' Here I select the genes of each set1 in the same manner as before. While for
#' geneset 2 and 3, I choose up and down-regulated genes only for each set.
#+ cache = T
set.seed(123)
gs1 <- sample(rst$ID, 10)
gs2 <- sample(rst[which(rst$logFC > 20 & rst$adj.P.Val < 0.0001), "ID"], 10)
gs3 <- sample(rst[which(rst$logFC < -20 & rst$adj.P.Val < 0.0001), "ID"], 10)

cat("set1:", gs1, "\nset2:", gs2, "\nset3:", gs3)
#' For the purpose of illustration and simplification, I randomly choose another `r n_bg` genes in the expression dataset as background gene.
#+ echo = F, cache = T
bg <- sample(rst$ID, n_bg)
gs <- Reduce(union, list(gs1, gs2, gs3, bg))
expr_g <- expr[gs, ]
#' ## Gene set variation calculation
#' Here I use the geneset 1 (randomly selected genes) to illustration the GSVA process.
#'
#+ hmap, fig.dim = c(4, 8), eval = F, echo = F
labels <- data.frame(status = rep(c("non-smokers", "smokers"), each = 5))
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
  cluster_rows = T, show_rownames = FALSE,
  cluster_cols = T, show_colnames = FALSE,
  silent = TRUE
)
grid::grid.draw(h)

#' ### Density function of gene expression levels

#+ fig.dim = c(8, 4), echo = F
par(mfrow = c(1, 2), cex = 0.8, font.main = 1)
plot(density(as.numeric(expr_g[gs1[1], ])),
  ylim = c(0, 2.5), xlim = c(0, 12),
  main <- "PDF of gene expression levels"
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
plot(ecdf(as.numeric(expr_g[gs1[1], ])),
  verticals = T, do.points = F, xlim = c(0, 12),
  main = "Empricial CDF of gene expression levels"
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

#' ## Normalizing gene expression levels
#' In order to be able to compare the gene expression levels on the same scale for the genes in a gene set. One will need to normalize the expression levels.This is done using culmulative density function $F_x$.
gd <- apply(expr_g, 1, function(x) {
  f <- ecdf(x)
  f(x)
}) %>%
  t() %>%
  data.frame() %>%
  setNames(colnames(expr_g))

gd <- apply(gd, 2, function(x) ifelse(x == 1, 0.99, x))
# normalize gene density
gene_density <- data.frame(log(gd / (1 - gd)))

#' Here I plot the **emprical CDF** instead of the explicit CDF mentioned in the paper for the purpose of illustration. For microarray data or continous data, the kernel for estimating CDF is Gaussian function, whereas in the count data, the kernel for CDF estimation is Poisson function `r emo::ji('fish')`.
#'
#' ![Gaussian kernel](images/gaus_kernel.png){width="300"}
#'
#' ![Poisson kernel](images/pois_kernel.png){width="300" height="86"}
#'
#+ cdf , echo = F, fig.dim = c(8, 4)
par(mfrow = c(1, 2), cex = 0.8, font.main = 1)
plot(density(as.numeric(gene_density[1, ])),
  main = "PDF after normalization"
)
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
  ylim = c(0, 1),
  main = "CDF after normalization"
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
#' Next, within each sample, the gene expression levels are reordered. To normalize further, a rank score is computated within each sample.
#'
#' ### Compute rank score
#' > To reduce the influence of potential outliers, we first convert $z_{ij}$ to ranks $z(i)j$ for each sample j and normalize further $r_{ij}=|p/2−z_{(i)j}|$ to make the ranks symmetric around zero.
#'
sort_sgn_idxs <- apply(gene_density, 2, order, decreasing = TRUE)
rownames(sort_sgn_idxs) <- rownames(expr_g)

rank_norm <- apply(sort_sgn_idxs, 2, compute_rank_score)
rownames(rank_norm) <- rownames(sort_sgn_idxs)

#' ### Enrichement score (ES)
#' The ES in GSVA is similar to GSVA using [two-sample Kolmogorov-Smirnov (KS)](https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test#Two-sample_Kolmogorov%E2%80%93Smirnov_test) like random walk statistic. Here I illustrate calculating ES in one sample. I set $\tau = 1$.
#+ cache = T
n_tot = nrow(expr_g)
n_gs = length(gs1)

sample <- rank_norm[, 1, drop = F] %>%
  data.frame() %>%
  tibble::rownames_to_column("ID") %>%
  mutate(gs1 = ID %in% gs1) %>%
  rename(rscore = GSM464697) %>%
    arrange(rscore)
sum_gs = sum(sample[sample$gs1, "rscore"])
dnom = 1/ (n_tot - n_gs)

sample <- sample %>%
  mutate(
    score = if_else(gs1, abs(rscore) / sum_gs, -as.numeric(!gs1) * dnom)
  ) %>%
  mutate(
    score = case_when(row_number() == 1 ~ 0, TRUE ~ score),
    es = cumsum(score)
  )
plot(sample$es,
  type = "l",
  xlab = "Gene rank",
  ylab = "Enrichment score",
  main = "Random walk in sample 1"
)
axis(side = 3, at = which(sample$gs1), label = F, col.ticks = "red")

#' #### KS statistics
#' Here I simplify the test to a classic KS test. Notice the difference in weight at each step from the random walk plot above.
#+ cache = T
ks_test_plot(sample, "ID", "rscore",
  gs1, "gs1 enrichment in sample 1",
  xlab = "rank score",
  exact = F,
  mexclude = T
)
#' Similarly one can look at the other two pathways
#+ fig.dim = c(8, 3), cache = T, echo = F
p2 <- ks_test_plot(sample, "ID", "rscore",
  gs2, "gs2 enrichment in sample 1",
  xlab = "rank score",
  exact = F,
  mexclude = T
)
p3 <- ks_test_plot(sample, "ID", "rscore",
  gs3, "gs3 enrichment in sample 1",
  xlab = "rank score",
  exact = F,
  mexclude = T
)
ggpubr::ggarrange(p2, p3, nrow = 1, ncol = 2)

#' ### Enrichemnent statistics
#' In GSVA, two approaches were prposed to turn KS random walk statistics into enrichment statistics,AKA, GSVA score.
#' > We offer two approaches for turning the KS like random walk statistic into an enrichment statistic (ES) (also called GSVA score), the classical maximum deviation method and a normalized ES. The first ES is the maximum deviation from zero of the random walk of the j-th sample with respect to the k-th gene set:
#'
#' $$ES_{jk}^{max}= \nu_{jk}[arg max_{l = 1..p}(|\nu_{jk}{l}|)]$$
#'
max(abs(sample[, "es", drop = F]))
# plot(sample$es, type = "l")
# abline(v = which.max(abs(sample$es)), col = "red", lty = 3)

#' #### Normalized enrichment score
#' $$ES_{jk}^{diff} = |ES_{jk}^+|+|ES_{jk}^-|$$
#'
#' > where $|ES_{jk}^+|$ and $|ES_{jk}^-|$ are the largest positive and negative random walk deviations from zero, respectively, for sample j and gene set k.
max(sample[, "es", drop = F], 0) - min(sample[, "es", drop = F], 0)
#' ####  Distribution of GSVA score
#'
#' Like the distribution in Eenrichment Score in GSEA, the null distribution of first option Eenrichment Statistics in GSVA  is also bimodal.
#'
#' > ![](images/gsva_null.png){width="520"}
#' > -- Supplemental material
#'
#' ## GSVA using R package
#'
#' In this section, I demonstrate GSVA using `GSVA` library in R.
#' The `gsva`function mainly takes expression data, gene set list and calculate enrichment statistics (GSVA score) with two options mentioned above. By default `gsva` and `gaussian` kernel are selected for the GSVA calculation. Other methods are implemented in this function as well (details in the next section).
library(GSVA)
gene_sets <- list(g1 = gs1, gs2 = gs2, gs3 = gs3)
expr_mtx <- as.matrix(expr_g)
es_gsva <- gsva(
  expr_mtx,
  gset.idx.list = gene_sets,
  mx.diff = F
)

#' ### Visualize pathway level enrichment in each sample
#+ fig.dim = c(6, 2), echo = F, cache = T
labels <- data.frame(status = rep(c("non-smokers", "smokers"), each = 5))
rownames(labels) <- c(ctrl_ids, case_ids)
grid::grid.newpage()
p <- pheatmap::pheatmap(es_gsva,
  annotation_col = labels,
  border_color = NA,
  legend_labels = "",
  main = "",
  cluster_rows = F, show_rownames = T,
  cluster_cols = F, show_colnames = T,
  fontsize = 8,
  fontsize_col = 5,
  fontsize_row = 10,
  angle_col = 45,
  silent = TRUE
)
grid::grid.draw(p)

#' ## Other GSVA-like methods
#' It is possible to use three other methods in the `GSVA` package. These are PLAGE, Z-scores and ssGSEA. I briefly describe these methods below.
#'
#' ### [PLAGE](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1261155/) @tomfohrPathwayLevelAnalysis2005
#' > standardizes each gene expression profile over the samples and then estimates the pathway activity profiles for each gene set as the coefficients of the first right-singular vector of the singular value decomposition of the gene set. It assume joint normal distributions in gene expression profiles.
#'
#' ### [Z-scores](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000217) @leeInferringPathwayActivity2008
#' > The combined z-score method standardizes first, as PLAGE, each gene expression profile into z-scores but the pathway activity profile is then obtained by combining the individual gene z-scores per sample. It assumes the gene expression profile are jointly normally distributioned. The combined z-score additionally assumes that genes act **independently** within each gene set.
#' ![](images/zscore.png)
#'
#' ### [ssGSEA](https://www.nature.com/articles/nature08460) @barbieSystematicRNAInterference2009
#' > It uses the difference in *empirical cumulative distribution functions* of gene expression ranks inside and outside the gene set to calculate an enrichment statistic per sample which is further normalized by the range of values taken throughout all gene sets and samples.
#' (The demonstration above uses eCDF.)
#'
#' ## Implication of GSVA score
#' The most frequencly used output from GSVA is the GSVA score, which is usually in a format of matrix. This matrix can be used for downstream analyses, such as hypothesis testing, classification, etc. There are some questions raise from using GSVA score in the downstream analysis.
#'
/*
#' <details><summary>Session Info</summary>
sessionInfo()
#' </details>
#+ echo = F, eval = F
# Markdown --------------------------------------------------------
# rmarkdown::render('src/06_GSVA.R', output_dir = 'output')
# knitr::spin('src/06_GSVA.R', format = 'Rmd', knit = FALSE)
*/
