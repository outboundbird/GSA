#' ---
#' title: query bioloical pathway DataBase
#' subtitle: 'SAR: NA , Study: NA'
#' author:  Siying Huang (E0482362), Biomarker statistics team
#' date: 'created: 2023-06-02 , updated (`r Sys.Date()`)'
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
knitr::opts_knit$set(root_dir = "/mnt/c/Users/e0482362/Work/pathway_analysis/src")
knitr::opts_chunk$set(echo = T, comment = "", message = F, warning = F, error = F)
options(width = 100)
#+ libs
library(here)
require(clusterProfiler)
require(enrichplot)
require(DOSE)
require(org.Hs.eg.db)
library(biomaRt)
library(DT)

#' Study dataset
#' In this tutorial, I still use the pulic smoking study dataset from NCBI GEO
#' database - GDS3713 (GSE18723).
#'
rst <- readRDS(file.path(here(), "data/de_rst.rds"))
annot <- readRDS(file.path(here(), "data/annot.rds")) %>%
  data.frame()

#' # Annotation
#' The GDS3713 dataset comes with an annotation file that has information
#' for genes in the array:
cat(paste0(names(annot), collapse = '\t'))
#' One can directly query a specific gene from this file. For example, we can
#' look up the ICOSLG information.
datatable(annot[which(annot$Gene.symbol == "ICOSLG"), ],
  rownames = F
)
#' However, most of time, we don't have such a file that comes with the
#' trancriptomic data. We will then need to annotate the genes.
#'
#' # Convert gene symbol to gene ID using `biomaRt`.
#' Since many pathway database requires entrez gene IDs as entries and most of
#' sequencing/ microarray output provides gene symbols or ensembl IDs, we'll
#' need to convert gene symbol to entrez IDs for the pathway analysis.
#'
#' We can query from the [ensembl](http://www.ensembl.org/index.html) site manually.
#' Or we can do better using`biomaRt`package.
#'  For a full tutorial of `biomaRt` package, see [here](https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html). It covers
#' some more advanced use case. Here I only cover basic usage of this package.
#'
#' ## 1. Connect to the selected *BioMart* database using `useEnsembl` function.
#' Before connecting to a specific database, we can first check available databases.
#'
#' ### list of available BioMart databases hosted by Ensembl
# available database
datatable(listEnsembl())
# available biological species
mart <- useEnsembl("genes")
datatable(listDatasets(mart))
#' Since we want to convert gene symbol to Entrez gene IDs in human subjects,
#' we select `biomart = 'genes'` and `dataset= 'hsapiens_gene_ensembl'`.
#' Use `useEnsembl` function to setup the connection with ensembl host.
# convert from gene symbol to gene ID

ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  mirror = "www"
)

#' ## 2. Annotate the gene set
#' Here instead of selecting genes in a biological pathway, I simply randomly
#' choose some genes in the dataset.
#'
set.seed(123)
gs1 <- sample(rst$ID, 10)
gs2 <- sample(rst$ID, 10)

#' Here I choose five most significantly differentially expressed genes,
#' five genes not differentially expressed.
gs3 <- c(head(rst$ID, 5), tail(rst$ID, 5))
print(glue::glue('geneset 1: {paste0(gs1, collapse = ", ")}
geneset2: {paste0(gs2, collapse =", ")}\ngeneset3: {paste0(gs3, collapse =", ")}'))

#' Retreive the annotated file with `getBM` function. Since we use gene symbols
#' to match aginst gene IDs, we'll need to specify the type of input data in
#' the `filter` argument (filters = "hgnc_symbol"). We would like to output
#' three annotation info `hgnc_symbol`,`entrezgene_id`.
#'
annotgs1 <- getBM(
  attributes = c(
    "hgnc_symbol", "entrezgene_id"
  ),
  filters = "hgnc_symbol",
  values = gs1,
  mart = ensembl
)
print(annotgs1)

#+ echo = F, cache = T
annot_list <- lapply(list(gs1, gs2, gs3), function(gs) {
  getBM(
    attributes = c(
      "hgnc_symbol", "entrezgene_id"
    ),
    filters = "hgnc_symbol",
    values = gs1,
    mart = ensembl
  )
}) %>% setNames(c("gs1", "gs2", "gs3"))

#' Note that not all gene symbols can be perfectly match against gene ID.
#' Alternatively, you can annotate all the gene in the dataset then extract
#' the ones in the genesets.

#' <details><summary>Session Info</summary>
sessionInfo()
#' </details>
#+ echo = F, eval = F
# Markdown --------------------------------------------------------
# rmarkdown::render('src/db_query.R', output_dir = 'output')