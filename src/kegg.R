#' ---
#' title: KEGG
#' subtitle: 'Original input: Rémi BRAZEILLES (E0425705)'
#' author:  Complied by Siying Huang (E0482362), Biomarker statistics team
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
knitr::opts_chunk$set(echo = T, comment = "", message = F, warning = T, error = T)
options(width = 100)

################################################################################
# Program information
#
# Program name        : KEGG_analysis.R
# Description         : How to perform a KEGG enrichment analysis with clusterProfiler package
# Author              : Rémi BRAZEILLES (E0425705)
# Date completed      : 2020-06-02
# Input files         : none
# Children programs   : none
# Outputs created     : kegg
# Packages needed     : clusterprofilter, enrichplot, DOSE, org.Hs.eg.db
# R version           : 3.6.1
# Modification status :
################################################################################
#+ libs
library(here)
require(clusterProfiler)
require(enrichplot)
require(DOSE)
require(org.Hs.eg.db)


#' ## Data loading
data(geneList, package = "DOSE")

#' ## Preparation of two vectors containing information on genes:
#' - the first one require the list of significant genes (for ex genes which are modualted by a treatment)

gene <- names(geneList)[abs(geneList) > 2]

#' ## Inputs needed for KEGG enrichment analysis
#' - gene : list of significant genes of interest
#' - organism : supported organism listed in 'http://www.genome.jp/kegg/catalog/org_list.html'
#' - pAdjustMethod : type of multiplicity control method: choice between "BH", "bonferroni", etc
#' - readable : if TRUE, results will be displayed by gene SYMBOL. if FALSE, results will  be displayed by gene keytype ID

kegg <- enrichKEGG(
  gene = gene,
  keyType = "ncbi-geneid",
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

#' Several functions of the enrichplot package are available to visualize GO or other enrichment analysis results
#'
#+ error = T
barplot(kegg, showCategory = nrow(kegg))
dotplot(kegg, showCategory = nrow(kegg))

#' Others visualizations are available. for more information please refer to [`clusterProfiler`](https://yulab-smu.github.io/clusterProfiler-book/chapter12.html)

#' <details><summary>Session Info</summary>
sessionInfo()
#' </details>
#+ echo = F
# Markdown --------------------------------------------------------
# rmarkdown::render('src/kegg.R', output_dir = 'output')