#' ---
#' title: GSEA
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
#' bibliography: references.bib
#' ---
#+ setup, include = FALSE
knitr::opts_chunk$set(echo = T, comment = "", message = F, warning = F, error = F)
options(width = 100)
#+
library(here)
library(dplyr)
source(file.path(here(), "src/utils.R"))
#+ io, cache = T
rst <- readRDS(file.path(here(), "data/de_rst.rds"))
gene_ids <- readRDS(file.path(here(), "data/gene_id.rds"))

#' # Gene Set Enrichment analysis
#' The idea of gene set enrichment analysis was proposed first by @moothaPGC1alpharesponsiveGenesInvolved2003.
#' Later on, [*@subramanianGeneSetEnrichment2005 *](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1239896/) made the modifictions based on their idea. This method now is widely used in the genetic study nowadays. These methods are considered as univariate funcitonal class scoring methods. It uses *gene scores* to summarize the differences between the comparison groups then uses gene set scores to summarize the expression levels of genes in a set and use a single statstic for signficance assessment.
#'
#' I'd like to point out that this enrichment analysis can be applied in other omic data, not soley in the transcriptomic field.
#'
#' **Disclaimer:  All quoted text are from original @subramanianGeneSetEnrichment2005 paper.**
#'
#' ## Scheme of gene set enrichment analysis
#' ![](https://www.pnas.org/cms/10.1073/pnas.0506580102/asset/c5e213a9-4247-4506-bae4-908054152f97/assets/graphic/zpq0370595180001.jpeg)
#'
#' ## Enrichment test
#' I demonstrate the principle of GSEA using the DGEA results from the [smoking study](https://link.springer.com/article/10.1007/s00251-010-0431-6).
#' And I use the same gene sets $(S)$ in the ORA session for illustration of GSEA.
#+ cache = T
set.seed(123)
gs1 <- sample(rst$ID, 10)
gs2 <- c(head(rst$ID, 5), tail(rst$ID, 5))
gs3 <- head(rst$ID, 10)
cat("set1:", gs1, "\nset2:", gs2, "\nset3:", gs3)
#' The whole ranked list $(L)$
#' There are many choices for ranking parameter. I choose t-statistics. Other options are fold changes, signal to noise ratio, etc.
#'
#' ### Sorting the DGEA result by t-statistics
rst <- dplyr::arrange(rst, desc(t)) %>%
  mutate(gs1 = ID %in% gs1)

t_stat <- setNames(rst$t, rst$ID)
#+ cache = T, fig.cap = "Gene rank by t-statistics. The red ticks represents the position of genes insided gene set1."
plot(unname(t_stat),
  ylab = "T-statistics",
  xlab = 'Rank',
  col = 'grey',
  type ='h',
  axes  = T
)

axis(
  side = 3, at = which(rst$gs1),
  col.ticks = "red",
  labels = F
)

#' ### Calculating Enrichment Score
#'
#' > The score is calculated by walking down the list L, increasing a running-sum statistic when we encounter a gene in S and decreasing it when we encounter genes not in S. The magnitude of the increment depends on the correlation of the gene with the phenotype. The enrichment score is the maximum deviation from zero encountered in the random walk;
#'
#' Here I use *gene set 1* as an example to illustrate the calculation of enrichment score (ES).
#'
#' **Method from @moothaPGC1alpharesponsiveGenesInvolved2003 **
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
#' **Method from @subramanianGeneSetEnrichment2005 **
#' $$\begin{equation*}\;P_{{\mathrm{hit}}}(S,i)={{\sum_{\begin{matrix}\scriptstyle{g_{j}{\in}S}\\ \scriptstyle{j{\leq}i}\end{matrix}}}}\frac{|r_{j}|^{p}}{N_{R}},\hspace{1em}{\mathrm{where}}{\;}N_{R}={{\sum_{g_{j}{\in}S}}}|r_{j}|^{p} \end{equation*}$$
#' $$\begin{equation*}\;P_{{\mathrm{miss}}}(S,i)={{\sum_{\begin{matrix}\scriptstyle{g_{j}{\not\in}S}\\ \scriptstyle{j{\leq}i}\end{matrix}}}}\frac{1}{(N-N_{H})}.\;\end{equation*}$$
#+ fig.cap="Comparison between two methods", cache = T
p_miss <- -1 / (n_tot - n_h)
nr <- sum(abs(t_stat[gs1]))
rst <- rst %>%
  mutate(score2_gs1 = if_else(gs1, abs(t) / nr, p_miss)) %>%
  mutate(
    score2_gs1 = case_when(row_number() == 1 ~ 0, TRUE ~ score2_gs1),
    es2_gs1 = cumsum(score2_gs1)
  )
  #+  plot, fig.dim = c(8, 4), cache = T
  par(mfrow = c(1, 2), cex = 0.8)
  plot(rst$es_gs1, type = "l", ylab = "Enrichment score", main = "ES (Mootha et al)")
  abline(h = 0, lty = 2, col = "red")
  axis(
    side = 3, at = which(rst$gs1),
    col.ticks = "red",
    labels = F
  )
  plot(rst$es2_gs1, type = "l", ylab = "Enrichment Score", main = "ES (Subramanian et al)")
  abline(h = 0, lty = 2, col = "red")
  axis(
    side = 3, at = which(rst$gs1),
    col.ticks = "red",
    labels = F
  )

#' Similarily, I can calculate the ES for gene set 2 and gene set 3 (Figure \@ref(fig:gs23) below).
#' And one can clearly see the extreme pattern in the random walk plot due to
#' the way that genes were selected in the sets. The Subramanian panelizes the
#' genes in the middle of rank by assigning less weight on the ES.
#'
#+ gs23, fig.dim = c(8,6), fig.cap ="Random walk of ES for gene set 2 and 3. Upper panel with Mootha method, lower panel with Subramanian method", cache = T
par(mfrow = c(2, 2), cex = 0.8, mar = c(4, 4, 1, 1))
calc_enrich_score(rst, gs2, "ID", "t", "mootha")
calc_enrich_score(rst, gs3, "ID", "t", "mootha")
calc_enrich_score(rst, gs2, "ID", "t")
calc_enrich_score(rst, gs3, "ID", "t")

#' #### Normalized Enrichment Score (NES)
#' > when adjusting for variation in gene set size, we normalize the $ES(S, \pi)$ for a
#' >  given S, separately rescaling the positive and negative scores by dividing by
#' >  their mean value, yielding $NES(S, \pi)$ and $NES(S)$ (normalized scores, NES).
#' >
#' > -- Supplement @subramanianGeneSetEnrichment2005
#'
#' Here I use gene set 1 as an example and I only illustrate the normalization
#' for the observed ES.
rst %>%
  group_by(gs1) %>%
  summarise(mean_es = mean(es2_gs1))

rst <- rst %>%
  mutate(nes = if_else(gs1, es2_gs1 / 0.16, es2_gs1 / 0.0853))
#' Ehe enrichment score  `r max(abs(rst[, "es2_gs1"]))` and the normalized
#' enrichement score is `r max(abs(rst[, 'nes']))`.
#'
#' ### Estimation of Significance Level of ES
#' This random walk feature correspond to [one-sample Kolmogorov–Smirnov](https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test#One-sample_Kolmogorov%E2%80%93Smirnov_statistic)-like statistic. In the Subramanian et al method, each step the augmentation is differentially weighted given the correlation with the profile or phenotype.
#'
#' #### Kolmogorov-Smirnov Tests
#' > **Null hypothei**s: No gene set is associated with the class distinction (@moothaPGC1alpharesponsiveGenesInvolved2003), where the rank ordering is used as measure of association.
#'
#' [K-S test](https://fr.wikipedia.org/wiki/Test_de_Kolmogorov-Smirnov) is designed to find general difference in the cumulative distribution. Here I plot simplied the K-S test of the three gene sets.
#'
#+ fig.dim = c(8, 3.5), cache = T
p1 <- ks_test_plot(rst, "ID", "t", gs1, "gs1")
p2 <- ks_test_plot(rst, "ID", "t", gs2, "gs2")
p3 <- ks_test_plot(rst, "ID", "t", gs3, "gs3")
ggpubr::ggarrange(p1, p2, p3, ncol = 3)
#' #### Significance assessment
#'
#' - Phenotype permutation
#'
#' In Subramanian et al method, the test significance is achieved by permutation:
#'
#' > We estimate the statistical significance (nominal P value) of the ES by using
#' > an empirical phenotype-based permutation test procedure that
#' >  **preserves the complex correlation structure of the gene expression data**.
#'
#' The phenotype is permutated and the ES scores are derived to achieve the
#' null distribution, whereas the genetic
#' structure is perserved. This method is wildly accepted in univariate genetic sudies.
#' However, it turned out to be controversial and received the criticism of
#' generating "false positive results" later on. Several methods were proposed to
#' alleviate the intercorrelation of the genes.
#'
#' - Gene sampling
#'
#' In gene sampling method, a group of genes are randomly selected from the
#' reference or background gene set to form the null gene sets. The the ES of
#' a given gene set $G_i$ is compared against the ES formed from the null
#' gene sets. The advantage of this method is that it does not depends on the
#' sample size therefore can be applied to the small sample size datasets.
#' However this method unrealistically assumes that all the genes in the
#' dataset are independent. Another draw back of this method is that
#' this method is computationally demanding.
#'
#' ###  Adjustment for Multiple Hypothesis Testing.
#' When testing biological pathway enrichment, several gene sets are usually queried.
#' Therefore GSEA usually involves adjustment for multiple comparisons.
#' Classic methods are availble in the R package for addressing this issue.
#' However, the pathways themselves can be correlated thus violating the
#' independence assumption in many of the methods for adjusting multiple testing,
#' such as Bonferroni, Benjamini-Hochberg.
#' To properly adjust for multiplicity issue, requires advanced techniques.
#' I will not discuss this subject here. A detailed discussion should be proposed
#' elsewhere.
#'
#' ## The Leading-Edge Subset analysis.
#'
#' > We define the leading-edge subset to be those genes in the gene set S
#' > that appear in the ranked list L at, or before, the point where the
#' > running sum reaches its maximum deviation from zero
#' > The leading-edge subset can be interpreted as the core of a gene set
#' > that accounts for the enrichment signal.
#' > ...
#' > High scoring gene sets can be grouped on the basis of leading-edge
#' > subsets of genes that they share. Such groupings can reveal which of
#' > those gene sets correspond to the same biological processes and
#' > which represent distinct processes.
#' > -- @subramanianGeneSetEnrichment2005
#'
#' I illustrate this idea with an example using gene set 1. Here I find the max
#' ES in the random walk plot. I can see the ES reaches to max at gene
#' `207479_at` then decreases. According to the definition by Subramanian et al,
#' The leading edge genes are EIF3I, DYRK3, and 207479_at.
#'
#+ cache = T
calc_enrich_score(rst, gs1, "ID", "t")
abline(v = which.max(rst$es2_gs1), lty = 3, col = "purple")
rst[rst$gs1, c("ID", "es2_gs1", "nes")]

#' ## Distribution of GSEA enrichment score
#'
#' Unlike classic gaussian Null distribution, the null distribution of GSEA
#' is **bimodal**.
#' ![GSEA paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1239896/bin/pnas_0506580102v2_06580Fig4A.jpg)
#'
#' This character was addressed in the supplemental material in their paper in the section
#' of computing statistical significance:
#'
#' > **Computing Significance By Using Positive or Negative Sides of the Observed and
#' > Null Bimodal ES Distributions**
#' >
#' >  ... constructing the null by
#' > using random phenotype assignments tends to produce a more symmetric distribution that
#' > may not exactly coincide with the bulk, nonextreme part of the distribution of the
#' > observed values. To address this issue, we determine significance and adjust for multiple
#' > hypotheses testing by independently using the positive and negative sides of the observed
#' > and null bimodal ES distributions. In this way, the significance tests [nominal P value,
#' > familywise-error rate (FWER), and false discovery rate (FDR)] are single tail tests on the
#' > appropriate (positive/negative) side of the null distribution.
#' >
#' > -- Supllemental material @subramanianGeneSetEnrichment2005
#'
#' ## GSEA using R packages
#'
#' ## Pros and Cons
#'
#' ## Reference
#'
/*
#' <details><summary>Session Info</summary>
sessionInfo()
#' </details>
#+ echo = F, eval = F
# Markdown --------------------------------------------------------
# rmarkdown::render('src/05_GSEA.R', output_dir = 'output')
# knitr::spin('src/05_GSEA.R', format = 'Rmd', knit = F )
*/