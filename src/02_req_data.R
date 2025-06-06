#' ---
#' title: Study dataset
#' subtitle: 'SAR: NA , Study: pathway analysis'
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
/*knitr::opts_chunk$set(echo = T, comment = '',message = F, warning = F, error=F)
options(width = 100)*/

/*
# require data from NCBI Gene Expression Omnibus (GEO) [database](https://www.ncbi.nlm.nih.gov/gds)
# datasets:
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS5022 AMI
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS5074 AMI 48 hrs
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS4500 AML train set
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS4501 AML test set
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS4181 AML
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS4844 CF
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS3115 CHF
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS3383 chronic stress
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS838 CML
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS3503 exercise
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS4345 GI cancer muscle tissue
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS1412 hormone therapy
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS3005 IL1 and IL6
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS5260 Psoriasis
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS3628 RA
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS5403 RA and OA
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS4266 renal transplant
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS5186 RIA
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS3898 social isolation
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS4193 SLE
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS3713 smoking
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS1436 smoking
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS1269 smoking
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS534 smoking
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS3318 scikle cell disease
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS4719 SLE
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS3920 MS
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS3886 MS
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS4270 UC
# publication of smoking study: https://link.springer.com/article/10.1007/s00251-010-0431-6
*/

#' # Study dataset
#'
#' In this tutorial, I still use the pulic smoking study dataset from NCBI
#' Gene Expression Omnibus (GEO) [database](https://www.ncbi.nlm.nih.gov/gds)
#' - GDS3713 (GSE18723).
#'
#' ### Summary of the study
#' This study aims to analyze peripheral circulating B cells from smoking and non-smoking healthy US white females. B cells are directly associated with the onset and development of many smoking-induced diseases. Results provide insight into the molecular basis of B cell involvement in smoking-related pathogenesis. It contains 39 smoker and 40 non-smokers. A total of 22283 genes of whom the expession levels are measured on Affymetrix HG-133A GeneChip microarray.
#'
#' ## Extracting data from NCBI Gene Expression Omnibus databse
#'
#' Using `GEOquery` library to retrieve data from NCBI database. The GDS dataset
#' is stored in the `gds` R4 class object.
#+ cache = T
library(here)
library(GEOquery)
library(dplyr)
gds <- getGEO("GDS3713")
# gds@header # meta data of the experiment
# gds@dataTable@table # expression data
# gds@dataTable@columns # phenotype data
#' structure of the `gds` object
str(gds, max.level = 3)
/*show(Meta(gds))*/
#' ### Assigning the expression and phenotype data
#' The expression data is stored in the `gds@dataTable@tabl` with gene in
#' The rows and subjects in the columns.
#' The phenotype data is store at `gds@dataTable@columns`.
expr <- Table(gds)
summary(expr[, 1:4])
pheno <- Columns(gds)
summary(pheno)

#' ## convert to GDS to expression set
#' Converting the gds object to expression set will generate a series of sub documents, including an annotation file of the genes in the expression data.

eset <- GDS2eSet(gds)
dim(pData(eset)) # pheno data, data.frame
str(eset, max.level =4)
# gene annotations
annot_cols <- eset@featureData@varMetadata
annot_data <- eset@featureData@data

/*# experiment info
# eset@experimentData
# assay <- eset@assayData
*/
/*#+ save, eval= F, echo = F
saveRDS(expr, "data/expr.rds")
saveRDS(pheno, "data/pheno.rds")
saveRDS(annot_data, "data/annot.rds")*/

/*
#' <details><summary>Session Info</summary>
sessionInfo()
#' </details>
# Markdown --------------------------------------------------------
# rmarkdown::render('src/02_req_data.R', output_dir = 'output')
# knitr::spin('src/02_req_data.R', format = 'Rmd', knit = FALSE)
*/
