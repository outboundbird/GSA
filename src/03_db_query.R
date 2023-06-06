#' ---
#' title: Query bioloical pathway databases
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
knitr::opts_chunk$set(echo = T, comment = "", message = F, warning = T, error = F)
options(width = 100)
#+ libs
library(here)
require(clusterProfiler)
require(enrichplot)
require(DOSE)
require(org.Hs.eg.db)
library(biomaRt)
library(DT)

#' # Study dataset
#'
#' In this tutorial, I still use the pulic smoking study dataset from NCBI GEO
#' database - GDS3713 (GSE18723).
#'
rst <- readRDS(file.path(here(), "data/de_rst.rds"))
annot <- readRDS(file.path(here(), "data/annot.rds")) %>%
  data.frame()

#' # Annotation of genes
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
#' Before connecting to a specific database, we can first check available
#' databases.
#'
#' ### List of available BioMart databases hosted by Ensembl
# available database
datatable(listEnsembl())
# available biological species for the gene databases
mart <- useEnsembl("genes", mirror = 'www')
datatable(listDatasets(mart))
#' Since we want to convert gene symbol to Entrez gene IDs in human subjects,
#' we select `biomart = 'genes'` and `dataset= 'hsapiens_gene_ensembl'`.
#' Use `useEnsembl` function to setup the connection with ensembl site host.
# convert from gene symbol to gene ID

ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  mirror = "www"
)

#' ## 2. Annotate genes
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
geneset 2: {paste0(gs2, collapse =", ")}\ngeneset 3: {paste0(gs3, collapse =", ")}'))

#' Retreive the annotated file with `getBM` function. Since we use gene symbols
#' to match aginst gene IDs, we'll need to specify the type of input data in
#' the `filter` argument (`filters = "hgnc_symbol"`). We would like to output
#' the `hgnc_symbol`,`entrezgene_id` reference.
#'
annotgs3 <- getBM(
  attributes = c("hgnc_symbol", "entrezgene_id"),
  filters = "hgnc_symbol",
  values = gs3,
  mart = ensembl
)
print(annotgs3)

#+ echo = F, cache = T, eval = F
annot_list <- lapply(list(gs1, gs2, gs3), function(gs) {
  getBM(
    attributes = c("hgnc_symbol", "entrezgene_id"),
    filters = "hgnc_symbol",
    values = gs,
    mart = ensembl
  )
}) %>% setNames(c("gs1", "gs2", "gs3"))

#' **Note that not all gene symbols can be perfectly match against gene ID.**
#' Alternatively, you can annotate all the gene in the dataset then extract
#' the ones in the genesets. For simplicity, I removed duplicated genes of which
#' have multiple probes correspond to the same gene (isoforms). In practice, the
#' removal of these isoforms should be done with examination of coefficient
#' variation (CV).

all_genes <- rst$ID[!duplicated(rst$ID)]
#' Upon removal of duplicates, there are `r length(all_genes)` genes left in
#' the dataset. We annotate these genes.
#+ cache = T, eval = F
annot_all <- getBM(
  attributes = c("hgnc_symbol", "entrezgene_id"),
  filters = "hgnc_symbol",
  values = all_genes,
  mart = ensembl
)
#+ eval = F, echo = F
saveRDS(annot_all, file.path(here(), "data/gene_id.rds"))

#' ## Alternative methods for annotation
#'
#' Alternatively one can use `bitr` function from `clusterProfiler` package.
#' Here we need to specify the database for curation. This requires the pacakge
#' `org.Hs.eg.db`. If you decide to go with method, make sure to update this
#' package often to obtain the most current information. The `org.Hs.eg.db`
#' package is updated biannually.

bitr(gs1,
  fromType = "SYMBOL",
  toType = c("ENTREZID", "SYMBOL"),
  OrgDb = org.Hs.eg.db
)

#' # Curate biological functions
#'
#' Upon selection genes of interest, we can start to curate the biological
#' functions of these genes. The most popular databases are GO, KEGG, Reactome,
#' etc. For full understanding the property of each of these databases, see
#' these references [xxx].
#'
#' For example, I'd like to query the biological functions of genes in the set 3
#' . One can use the `groupGO` function to look up biological functions:
#'
#' >  - MF: Molecular Function
#' >     - molecular activities of gene products
#' >  - CC: Cellular Component
#' >     - where gene products are active
#' >  - BP: Biological Process
#' >     - pathways and larger processes made up of the activities of multiple
#' gene products
#'
#' Additionally, one can look up onto [KEGG](https://www.genome.jp/kegg/ko.html),
#' [WikiPathway](https://www.wikipathways.org/), [Reactome](https://reactome.org/),
#' [Diease Ontology (DO)](https://disease-ontology.org/),
#' [Medical Subject Headings (MeSH)](https://www.nlm.nih.gov/mesh/meshhome.html),
#' [Molecular Signatures Database (MSigDb)](https://www.gsea-msigdb.org/gsea/msigdb).
#' I didn't find any ready to use functions from `clusterProfile` package,except
#' for `GO` query. However, personally, I don't think statiticians are experts to
#' interpret the finding in biological sense. Since this is less of interest for
#' analysts, the `clusterProfiler` simplifed this curation process and combined
#' with the tests (ORA, GSEA) in their function construction. For details see
#' separate files.
#+ cache = T
ggo <- groupGO(
  gene = as.character(annotgs3$entrezgene_id),
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  level = 3,
  readable = TRUE
)
tail(ggo)
#' # Reference
#'
#' <details><summary>Session Info</summary>
sessionInfo()
#' </details>
#+ echo = F, eval = F
# Markdown --------------------------------------------------------
# rmarkdown::render('src/03_db_query.R', output_dir = 'output')
# knitr::spin('src/03_db_query.R', format = 'Rmd', knit = FALSE)