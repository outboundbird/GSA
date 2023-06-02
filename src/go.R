#' ---
#' title: GO analysis
#' subtitle: 'Original input: Rémi BRAZEILLES (E0425705) '
#' author:  compiled by Siying Huang (E0482362), Biomarker statistics team
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
knitr::opts_chunk$set(echo = T, comment = '',message = F, warning = F, error=F)
options(width = 100)
################################################################################
# Program information
#
# Program name       : GO.analysis.R
# Description        : How to perform a GO enrichment analysis with clusterProfiler package
# Author             : Rémi BRAZEILLES (E0425705)
# Date completed     : 2020-06-02
# Input files        : none
# Children programs  : none
# Outputs created    : ego
# Packages needed    : clusterprofilter, enrichplot, DOSE, org.Hs.eg.db
# R version          : 3.6.1
# Modification status :
################################################################################
#+ libs
library(here)
require(clusterProfiler)
require(enrichplot)
require(DOSE)
require(org.Hs.eg.db)

#' ## Data loading
data(geneList, package="DOSE")
head(geneList)
#' ## Preparation of two vectors containing information on genes:
#' - the first one require the list of significant genes (for ex genes which are modualted by a treatment)
#' - the second contains the whole list of genes on which the analysis was done: called "the universe".
#'
gene <- names(geneList)[abs(geneList) > 2]

#' If needed, possibility to match gene name to get the entrezid, or vice-versa with the following line Or use other type of gene ID like ENSEMBLE ID.
#'
gene.df <- bitr(gene, fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)
head(gene.df)

#' ## Inputs needed for GO enrichment analysis
#' - gene: list of significant genes of interest
#' - universe: list of all genes
#' - OrgDb: annotation DB
#' - keytype: which gene ID is used in gene and universe: choice between all genes ID like "ENTREZID","SYMBOL","ENSEMBL"
#' - ont: on which ontologies the analysis has to be done: choice between "CC" (cellular component), "BP" (Biological process), "MF" (Molecular Function) and "ALL" for the 3
#' - pAdjustMethod: type of multiplicity control method: choice between "BH", "bonferroni", etc
#' - readable: if TRUE, results will be displayed by gene SYMBOL. if FALSE, results will  be displayed by gene keytype ID

ego <- enrichGO(gene          = gene,
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

#' Several functions of the enrichplot package are available to visualize GO or other enrichment analysis results
#+ fig.dim = c(8,5)
barplot(ego, showCategory=10)
dotplot(ego, showCategory=10)

#' Others visualizations are available. for more information please refer to [`clusterProfiler`](https://yulab-smu.github.io/clusterProfiler-book/chapter12.html).

#' <details><summary>Session Info</summary>
sessionInfo()
#' </details>
#+ echo = F
# Markdown --------------------------------------------------------
# rmarkdown::render('src/go.R', output_dir = 'output')