#' ---
#' title: Other methods in gene set analysis
#' subtitle: 'SAR: sar , Study: study'
#' author:  Siying Huang (E0482362), Biomarker statistics team
#' date: 'created: 2023-06-08 , updated (`r Sys.Date()`)'
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
#+ libs
library(here)
library(dplyr)
library(ggplot2)
source(file.path(here(), "src/utils.R"))
#+
rst <- readRDS(file.path(here(), "data/de_rst.rds"))
rownames(rst) <- rst$ID
# gene_ids <- readRDS(file.path(here(), "data/gene_id.rds"))

#' # Other methods
#' These methods are challengers of the GSEA methods that first proposed by
#' @moothaPGC1alpharesponsiveGenesInvolved2003. Many of them critize the Kolmogorov–Smirnov
#' test used for enrichment score are not sensitive. And all of these alternatives
#' claims that they have better power than K-S test in GSEA. The others aim to
#' solve the issue of intercorrelation among genes within the same or close
#' biological process.
#'
#' - [maxmean statsitics, restandardization](https://arxiv.org/pdf/math/0610667.pdf) @efronTestingSignificanceSets2007
#' - [mean shift test](https://www.pnas.org/doi/10.1073/pnas.0506577102) @tianDiscoveringStatisticallySignificant2005 , @irizarryGeneSetEnrichment2009a
#' - [PAGE](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-6-144) @kimPAGEParametricAnalysis2005
#'  - [SAFE](https://academic.oup.com/bioinformatics/article/21/9/1943/408983?login=true) @barrySignificanceAnalysisFunctional2005
#' - [Camera](https://academic.oup.com/nar/article/40/17/e133/2411151?login=true) @wuCameraCompetitiveGene2012
#' - [ROAST](https://academic.oup.com/bioinformatics/article/26/17/2176/200022?login=true) @wuROASTRotationGene2010
#'
#' Here I use the smoking study to illustrate the ideas of these methods. Instead of choosing a particular biological pathway, I choose these pathway with intention to evaluate the property of each method.
#'
#' # Gene set selection
#' The results for the DEG analysis were sorted based on p values.
#' Here I generate three gene sets:
#'
#' - geneset 1: randomly select 10 genes
#' - geneset 2: select 50 genes that are not differentially expressed
#' - geneset 3: select top 30 differentially expressed genes
set.seed(123)
gs1 <- sample(rst$ID, 10)
gs2 <- sample(rst[rst$P.Value > 0.5, "ID"], 50)
gs3 <- head(rst$ID, 30)
#' Of note, in gene set 3, 50% of the genes are up-regulated, 50% of genes are
#' down regulated.
#+ echo = F
gs3set <- rst[gs3, ] %>%
  mutate(reg = if_else(logFC > 0, "UP", "DOWN"))

table(gs3set[, "reg"])
#' ## Comparison of the distributions of gene sets
#'
#+ fig.dim= c(9, 4), echo = F, cache = T
rst <- rst %>%
  mutate(
    gs1 = ID %in% gs1,
    gs2 = ID %in% gs2,
    gs3 = ID %in% gs3
  )
par(mfrow = c(1, 3), font.main = 1, cex.main = 1)
plot(density(rst$t),
  ylim = c(0, 1),
  lwd = 4,
  col = "grey",
  main = "t-statistic"
)
lines(density(rst[rst$gs1, "t"]), ylim = c(0, 1), col = "red")
lines(density(rst[rst$gs2, "t"]), ylim = c(0, 1), col = "darkgreen")
lines(density(rst[rst$gs3, "t"]), col = "blue")
abline(v = 0, lty = 3, col = "gray")
legend("topleft",
  c("genome", "gs1", "gs2", "gs3"),
  col = c("grey", "red", "darkgreen", "blue"),
  lty = 1,
  bty = "n",
  cex = 0.8
)

plot(density(rst$logFC),
  lwd = 4,
  col = "grey",
  main = "Log(FC)",
  xlim = c(-250, 250)
)
lines(density(rst[rst$gs1, "logFC"]), ylim = c(0, 1), col = "red")
lines(density(rst[rst$gs2, "logFC"]), ylim = c(0, 1), col = "darkgreen")
lines(density(rst[rst$gs3, "logFC"]),  ylim = c(0, 1),  col = "blue")

plot(density(rst$logFC),
  lwd = 4,
  col = "grey",
  main = "Log(FC)"
)
lines(density(rst[rst$gs1, "logFC"]), ylim = c(0, 1), col = "red")
lines(density(rst[rst$gs2, "logFC"]), ylim = c(0, 1), col = "darkgreen")
lines(density(rst[rst$gs3, "logFC"]), ylim = c(0, 1), col = "blue")

mtext("Results from DEG analysis",
  cex = 1, side = 3,
  line = 2, outer = F, adj = 0
)

#' ## Parametric Analysis of Gene Set Enrichment (PAGE)
#' > PAGE calculates a Z score for a given gene set from a parameter such as **fold change** value between two experimental groups and infers statistical significance of the Z score against standard normal distribution.
#' > -- @kimPAGEParametricAnalysis2005
#'
#' **Null hypothesis:**
#'  the distribution of the parameters for a given geneset follows normal distribution.
#'
m <- mean(rst$logFC)
s <- sd(rst$logFC)
# calc z score
rst <- rst %>%
  mutate(z_fc = (logFC - m) / s)
#+ echo = F, cache = T, fig.dim = c(9,4)
# KS test
par(mfrow = c(1, 3), font.main = 1, cex.main = 1)
qqnorm(rst$z_fc, main = "Genomewide Z scores")
qqline(rst$z_fc)

qqnorm(rst[gs3, "z_fc"],
  xlab = "Z score",
  main = "Gene set 3"
)
qqline(rst[gs3, "z_fc"])

rst_ks <- ks.test(rst[gs3, "z_fc"], "pnorm")
plot(ecdf(rst[gs3, "z_fc"]),
  do.points = F, verticals = T,
  main = "Gene set 3"
)
plot(ecdf(rnorm(length(gs3))),
  do.points = F, verticals = T,
  lty = 3, lwd = 2,
  add = T, col = "red"
)
mtext(
  c(paste(
    "D =", round(rst_ks$statistic, 2),
    ", p =", round(rst_ks$p.value, 4)
  )),
  cex = 0.7
)
legend("topleft", c("Theoretical", "Observed"),
  col = c("red", "black"),
  lty = c(2, 1),
  bty = "n",
  cex = 0.8
)

#' The other two gene sets. Notice that gene set 2 is composed of 50 non-significantly
#' differentially expressed genes. However the distribution of the fold change
#' does not follows unimodal gaussian distribution.
#+ fig.dim = c(5,7), echo = F
par(mfrow = c(3, 2), mar = c(4, 3, 2, 1))
qqnorm(rst[gs1, "z_fc"],
  xlab = "Fold change",
  main = "Gene set1"
)
qqline(rst[gs1, "z_fc"])
ref_ks_test(rst, gs1, "z_fc", title = "gs1")

qqnorm(rst[gs2, "z_fc"],
  xlab = "Fold change",
  main = "Gene set2"
)
qqline(rst[gs2, "z_fc"])
ref_ks_test(rst, gs2, "z_fc", title = "gs2")

#' ## Mean shift methods
#' These are the challengers of K-S test in the GSEA.
#' The first test: t-test or wilcoxon test, as gene set gets large it approximate
#' Z distribution.
#' $E_g^z= \sqrt{N_g}\bar{t}$ with $\bar{t}= 1/N_g\sum_{i \in A_g} t_i$
#'
#' @tianDiscoveringStatisticallySignificant2005 propose to estimate the significance
#' of a particular gene set by constructiong null distribution using permutations.
#'
#' ![](images/tian.jpeg){width=400}
m_gs1 <- mean(rst[gs1, "t"])
z_gs1 <- sqrt(length(gs1)) * m_gs1
pnorm(z_gs1, lower.tail = F)

z_test(rst, gs2, "t")
z_test(rst, gs3, "t")

#' Chi square test was proposed by @irizarryGeneSetEnrichment2009a
#' to circumvent that the one-sample z test cannot detect the pathway where
#' it contain half up-regulated, half-downregulated genes.
#'
#' > A possible advantage of GSEA, i.e. the K-S test, over the one sample z-test is that the latter is specifically designed to identify gene sets with mean shifts and the K-S test is designed to find general difference in the cumulative distribution. In principle, we want to be able to detect gene sets for which some members are up-regulated and others are down-regulated. The z-test is not sensitive to this change as there is no shift in mean. We therefore, propose the use of another standard statistical test useful for detecting changes in scale: the χ2 test. -- @irizarryGeneSetEnrichment2009a
#'
#'
#' $E_g^{X^2} = \frac{\sum_{i\in A_g}(t_i - \bar{t})^2 -(N_g -1) }{2(N_g -1)}$
#'
x2_gs1 <- (sum((rst[gs1, "t"] - m_gs1)^2) - (length(gs1) - 1)) / (2 * (length(gs1) - 1))

1 - pchisq(x2_gs1, df = length(gs1) - 1)
pchisq(x2_gs1, df = 1)

chisq_test <- function(t_vec, df = NULL) {
  if (is.null(df)) df <- length(t_vec) - 1
  t_bar <- mean(t_vec, na.rm = T)
  m <- length(t_vec) - 1
  numr <- sum((t_vec - t_bar)^2) - m
  x2 <- numr / (2 * m)
  print(x2)
  1 - pchisq(x2, df = df)
}

t_vec <- rst[gs3, "t"]
chisq_test(rst[gs1, "t"])
chisq_test(rst[gs2, "t"])
chisq_test(rst[gs3, "t"])

#' ## Maxmean statistics and restarndardization
#' This method was proposed by @efronTestingSignificanceSets2007 and is available in the R package `GSA`.
library(GSA)
expr <- readRDS(file.path(here(), "data/expr.rds"))
# I simplily removed the duplicated genes here. In practice removing
# duplicated genes should be given careful consideration.
expr <- expr[!duplicated(expr$IDENTIFIER), ] %>%
  `rownames<-`(NULL) %>%
  tibble::column_to_rownames("IDENTIFIER")%>%
  dplyr::select(-ID_REF)

pheno <- readRDS(file.path(here(), "data/pheno.rds")) %>%
  mutate(status = if_else(stress == "control", 1, 2)) %>%
  dplyr::select(status)
class(pheno)
GSA.func(
  as.matrix(log(expr)),
  pheno[, 1],
  genesets = list(gs1),
  genenames= rownames(expr), resp.type = "Two class unpaired"
)

summary(expr)

expr[gs2,]

x <- matrix(rnorm(1000 * 20), ncol = 20)
summary(x)
dd <- sample(1:1000, size = 100)
u <- matrix(2 * rnorm(100), ncol = 10, nrow = 100)
x[dd, 11:20] <- x[dd, 11:20] + u
y <- c(rep(1, 10), rep(2, 10))
genenames <- paste("g", 1:1000, sep = "")
# create some random gene sets
genesets <- vector("list", 50)
for (i in 1:50) {
  genesets[[i]] <- paste("g", sample(1:1000, size = 30), sep = "")
}
geneset.names <- paste("set", as.character(1:50), sep = "")
GSA.func.obj <- GSA.func(x, y, genenames = genenames, genesets = genesets, resp.type = "Two class unpaired")

str(GSA.func.obj)

#' ## SAFE
#'
#' ## Camera
#'
#' ## Roast
#'
#'
#' # Comparison of gene set analysis methods
#'
#' [Comparison of gene set scoring methods for reproducible evaluation of multiple tuberculosis gene signatures](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9882404/)
#'
#' ## Reviews
#' - [Pathway Analysis: State of the Art](https://www.frontiersin.org/articles/10.3389/fphys.2015.00383/full) @garcia-camposPathwayAnalysisState2015
#' - [Analyzing gene expression data in terms of gene sets: methodological issues](https://academic.oup.com/bioinformatics/article/23/8/980/198511?login=true)
#' - [Gene set analysis methods: statistical models and methodological differences](https://academic.oup.com/bib/article/15/4/504/407653?login=true) @maciejewskiGeneSetAnalysis2014
#' - The Statistical Properties of Gene-Set Analysis @deleeuwStatisticalPropertiesGeneset2016
#' - @khatriTenYearsPathway2012
#' - @malekiGeneSetAnalysis2020
#' - @malekiMethodChoiceGene2019
#' - @mathurGeneSetAnalysis2018a
#'
#'
#' ## Reference
#'
/*#' <details><summary>Session Info</summary>
sessionInfo()
#' </details>
#+ echo = F, eval = F
# Markdown --------------------------------------------------------
# rmarkdown::render('src/07_others.R', output_dir = 'output')
# knitr::spin('src/07_others.R', format = 'Rmd', knit = FALSE)*/