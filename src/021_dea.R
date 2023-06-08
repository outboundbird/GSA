#' ---
#' title: Differential expression analysis
#' subtitle: 'SAR: NA , Study: NA'
#' author:  Siying Huang (E0482362), Biomarker statistics team
#' date: 'created: 2023-06-01 , updated (`r Sys.Date()`)'
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
library(limma)
library(dplyr)
source(file.path(here(), "src/utils.R"))
#' # Differential gene expression analysis (DGEA)
#' > Comparison of differential gene expression patterns in these conditions has enabled the identification of common elements that are significantly enriched in gene classes with particular functions such as protein synthesis, hormone delivery, and morphological plasticity.
#' >      - [Hormones, Brain and Behavior (Third Edition), 2017](https://www.sciencedirect.com/topics/biochemistry-genetics-and-molecular-biology/differential-gene-expression)
#'
#' Differential expression analysis result is a prerequisite for some of the
#' pathway analysis methods. This chapter does not aim to demonstrate how to
#' properly peform DGEA. Details and explainations of DEGA are omitted.
#+
expr <- readRDS(file.path(here(), "data/expr.rds"))
pheno <- readRDS(file.path(here(), "data/pheno.rds"))
# data management
pheno <- pheno %>%
  mutate(status = as.factor(as.integer(stress)))
#' In this micro-array data due to that there are multiple probes match to
#' the same gene (isoforms).
#' For simplicity, I removed duplicated genes of which
#' have multiple probes correspond to the same gene . In practice, the
#' removal of these isoforms should be done with examination of coefficient
#' variation (CV). Before removing the isforms, there are {{nrow(expr)}} rows.
expr <- expr[!duplicated(expr$IDENTIFIER), ]
#' Removing the multiple isoforms leaves {{nrow(expr)}} in the dataset.

/*# match(names(expr[, 3:ncol(expr)]), pheno$sample)
*/

#' ## Comparison between smokers and non-smokers
#' A linear model is used with the constrast made between smokers and controls.
#' In this example, I use `limma` package and derived test statistics with
#' emperical bayes mderation on the stand errors.
#'
#+ cache = T
y <- as.matrix(expr[, 3:ncol(expr)])
rownames(y) <- expr$IDENTIFIER
group <- pheno[["status"]]
mod <- model.matrix(~ 0 + group)
fit <- lmFit(y, mod)
contr <- makeContrasts(
  contrasts = "group1-group2",
  levels = colnames(coef(fit))
)
tmp <- contrasts.fit(fit, contr)
ebfit <- eBayes(tmp, robust = T)
rst_eb <- topTable(ebfit, number = Inf) %>%
  tibble::rownames_to_column('ID')
rst <- setNames(rst_eb, c("ID", "estimate", "AveExp", "t", "p.value", "FDR", "B"))
#' ## Results visualization
#' The marked genes pass the FDR significant threshold at adjusted p-value <= 0.01.
#' This DEA result will be used for the following pathway analysis.
#'
#+ fig.dim = c(9,4), cache = T
par(mfrow = c(1, 3))
hist(rst_eb$logFC, breaks = 50, main = 'log (Fold change)')
hist(rst_eb$P.Value, breaks = 50, main = 'P values')
qqplot(
  -log10(runif(nrow(rst_eb))),
  -log10(rst_eb$P.Value),
  xlab = 'Theoretical uniform distribution',
  ylab = '-log10(P values)'
)
qqline(-log10(rst_eb$P.Value), distribution = qunif, col = "red")

#+ fig.dim = c(9,8), cache = T
volc_plot(rst, lab_col = "ID", lab_method = "fdr", fdr_thresh = 0.01)

/*#+ echo = F, eval = F
saveRDS(rst_eb, "data/de_rst.rds")
*/
/*
#' <details><summary>Session Info</summary>
sessionInfo()
#' </details>
#+ echo = F, eval = F
# Markdown --------------------------------------------------------
# rmarkdown::render('src/021_dea.R', output_dir = 'output')
# knitr::spin('src/021_dea.R', format = 'Rmd', knit = FALSE)*/