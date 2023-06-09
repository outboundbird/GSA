#' ---
#' title: title
#' subtitle: 'SAR: sar , Study: study'
#' author:  Siying Huang (E0482362), Biomarker statistics team
#' date: 'created: 2023-06-06 , updated (`r Sys.Date()`)'
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
#+
library(here)
library(dplyr)
source(file.path(here(), "src/utils.R"))
#+ io
rst <- readRDS(file.path(here(), "data/de_rst.rds"))
gene_ids <- readRDS(file.path(here(), "data/gene_id.rds"))

#' # Gene Set Enrichment analysis
#' ![](https://www.gsea-msigdb.org/gsea/images/GSEA-logo.gif)
#' The first Gene set enrichment proposed
#' @moothaPGC1aresponsiveGenesInvolved2003
#' ![](https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fng1180/MediaObjects/41588_2003_Article_BFng1180_Fig1_HTML.gif?as=webp)
#' @subramanianGeneSetEnrichment2005

#' ## Enrichment test
#'
#' ![](https://www.ncbi.nlm.nih.gov/core/lw/2.0/html/tileshop_pmc/tileshop_pmc_inline.html?title=Click%20on%20image%20to%20zoom&p=PMC3&id=1266115_zpq0370595180001.jpg)
#'
#' I demonstrate the principle of GSEA using the DGEA results.
#' And I use the same gene sets (S) in the ORA session
#'

set.seed(123)
gs1 <- sample(rst$ID, 10)
gs2 <- c(head(rst$ID, 5), tail(rst$ID, 5))
gs3 <- head(rst$ID, 10)
cat("set1:", gs1, "\nset2:", gs2, "\nset3:", gs3)
#' The whole ranked list (L)
#' the whole  list can be ranked based on different parameters here i use t-statistics. it can also be fc, signal to noise ratio.
#' Sorting the DGEA result by t-statistics

rst <- dplyr::arrange(rst, desc(t)) %>%
  mutate(gs1 = ID %in% gs1)

t_stat <- setNames(rst$t, rst$ID)

barplot(unname(t_stat),
  ylab = "T-statistics",
  xlab = 'Rank',
  col = 'grey',
  border = NA,
  axes  = T
)
axis(
  side = 1, at = seq(1, length(t_stat), 1000),
  labels = F
)
axis(
  side = 3, at = which(rst$gs1),
  col.ticks = "red",
  labels = F
)


#' ### Calculating Enrichment Score
#'
#' > The score is calculated by walking down the list L, increasing a running-sum statistic when we encounter a gene in S and decreasing it when we encounter genes not in S. The magnitude of the increment depends on the correlation of the gene with the phenotype. The enrichment score is the maximum deviation from zero encountered in the random walk; [cite gsea]
#'
#' Here I use *gene set 1* as an example to illustrate the calculation of enrichment score (ES).
#'
#' **Method from *Mootha et al 2003* ** [cite]
#' $$P_{hit}(S,i)= \sum_{g_{j}{\in}S}{\sqrt{\frac{(N-N_h)}{N_h}}}$$
#' $$P_{miss}(S,i) =\sum_{g_{j}{\not\in}S}{ \sqrt{\frac{N_h}{(N-N_h)}}}$$
n_h <- length(gs1)
n_tot <- length(t_stat)
# method1
es_neg <- -sqrt(n_h / (n_tot - n_h))
es_pos <- sqrt((n_tot - n_h) / n_h)
rst <- rst %>%
  mutate(score_gs1 = if_else(gs1, es_pos, es_neg)) %>%
  mutate(
    score_gs1=case_when(row_number()==1 ~0 , TRUE ~ score_gs1),
  es_gs1 = cumsum(score_gs1))

#' The max value of the enrichment score is `r round(max(abs(rst$es_gs1)),2)`
#'
#' **Method from *Subramanian et al 2005* **[cite]
#' $$\begin{equation*}\;P_{{\mathrm{hit}}}(S,i)={{\sum_{\begin{matrix}\scriptstyle{g_{j}{\in}S}\\ \scriptstyle{j{\leq}i}\end{matrix}}}}\frac{|r_{j}|^{p}}{N_{R}},\hspace{1em}{\mathrm{where}}{\;}N_{R}={{\sum_{g_{j}{\in}S}}}|r_{j}|^{p} \end{equation*}$$
#' $$\begin{equation*}\;P_{{\mathrm{miss}}}(S,i)={{\sum_{\begin{matrix}\scriptstyle{g_{j}{\not\in}S}\\ \scriptstyle{j{\leq}i}\end{matrix}}}}\frac{1}{(N-N_{H})}.\;\end{equation*}$$
p_miss <- -1 / (n_tot - n_h)
nr <- sum(abs(t_stat[gs1]))
rst <- rst %>%
  mutate(score2_gs1 = if_else(gs1, abs(t) / nr, p_miss)) %>%
  mutate(
    score2_gs1 = case_when(row_number() == 1 ~ 0, TRUE ~ score2_gs1),
    es2_gs1 = cumsum(score2_gs1)
  )
  #+ walking plot, fig.dim = c(8, 4)
  par(mfrow = c(1, 2), cex = 0.8)
  plot(rst$es_gs1, type = "l", ylab = "Enrichment score", main = "ES (Mootha et al)")
  abline(h = 0, lty = 2, col = "red")
  axis(
    side = 3, at = which(rst$gs1),
    col.ticks = "red",
    labels = F
  )
  plot(rst$es2_gs1, type = "l", ylab = "Enrichment Score", main = "ES score (Subramanian et al 2005)")
  abline(h = 0, lty = 2, col = "red")
  axis(
    side = 3, at = which(rst$gs1),
    col.ticks = "red",
    labels = F
  )

#' Similarily, I can calculate the ES for gene set 2 and gene set 3 (Figure below). And one can clearly see the extreme pattern in the walking plot due to the way that genes were selected in the sets.
#'
#+ fig.dim = c(8,4)
par(mfrow = c(1, 2), cex = 0.8)
calc_enrich_score(rst, gs2, "ID", "t", "mootha")
calc_enrich_score(rst, gs3, "ID", "t", "mootha")

#' ## Estimation of Significance Level of ES
#' This random walk feature correspond to Kolmogorov–Smirnov-like statistic. In the Subramanian et al method, each step the augmentation is differentially weighted given the correlation with the profile or phenotype.
#'
#' ### Kolmogorov-Smirnov Tests
#'
#' [K-S test](https://fr.wikipedia.org/wiki/Test_de_Kolmogorov-Smirnov) is designed to find general difference in the cumulative distribution. Here I plot simplied the K-S test of the three gene sets.
#'

p1 <- ks_test_plot(rst, "ID", "t", gs1, "gs1")
p2 <- ks_test_plot(rst, "ID", "t", gs2, "gs2")
p3 <- ks_test_plot(rst, "ID", "t", gs3, "gs3")
ggpubr::ggarrange(p1, p2, p3, ncol = 3)

#' limitations:
#' ### Permutation test
#'
#' ###  Adjustment for Multiple Hypothesis Testing.
#'
#' ## The Leading-Edge Subset.



#' ## GSEA using R packages
/*
#' <details><summary>Session Info</summary>
sessionInfo()
#' </details>
#+ echo = F, eval = F
# Markdown --------------------------------------------------------
# rmarkdown::render('src/05_GSEA.R', output_dir = 'output')
# knitr::spin('src/05_GSEA.R', format = 'Rmd', knit = F )
*/