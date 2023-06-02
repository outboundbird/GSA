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
#' # Experiment data
#' In this example, I use the public dataset GDS3713 (GSE18723) downloaded from NCBI GEO database. A short description of this experiment:
#'
#' > Analysis of peripheral circulating B cells from smoking and non-smoking
#' > healthy US white females. B cells are directly associated with the onset
#' > and development of many smoking-induced diseases. Results provide insight
#' > into the molecular basis of B cell involvement in smoking-related pathogenesis.
#'
#' There are 79 samples contained in this dataset, 39 smokers and 40 controls.
#' A total of 22283 genes of whom the expession levels are measured on
#' Affymetrix HG-133A GeneChip microarray.
#'
#+
expr <- readRDS(file.path(here(), "data/expr.rds"))
pheno <- readRDS(file.path(here(), "data/pheno.rds"))
# data management
pheno <- pheno %>%
  mutate(status = as.factor(as.integer(stress)))

# match(names(expr[, 3:ncol(expr)]), pheno$sample)

#' # Differential expression analysis
#' A linear model is used with the constrast made between smokers and controls.
#' In this example, I use `limma` package and derived test statistics with
#' emperical bayes mderation on the stand errors.
#'
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
rst_eb <- topTable(ebfit, number = Inf)
rst <- setNames(rst_eb, c("ID", "estimate", "AveExp", "t", "p.value", "FDR", "B"))

#' The marked genes pass the FDR significant threshold at adjusted p-value <= 0.01.
#' This DEA result will be used for the following pathway analysis.
#'
#+ fig.dim = c(9,4)
par(mfrow = c(1, 3))
hist(rst_eb$logFC)
hist(rst_eb$P.Value)
qqplot(-log10(runif(nrow(rst_eb))), -log10(rst_eb$P.Value))
qqline(-log10(rst_eb$P.Value), distribution = qunif, col = "red")

#+ fig.dim = c(9,8)
volc_plot(rst, lab_col = "ID", lab_method = "fdr", fdr_thresh = 0.01)
#+ echo = F, eval = F
saveRDS(rst_eb, "data/de_rst.rds")

#' <details><summary>Session Info</summary>
sessionInfo()
#' </details>
#+ echo = F, eval = F
# Markdown --------------------------------------------------------
# rmarkdown::render('src/dea.R', output_dir = 'output')