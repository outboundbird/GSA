#' ---
#' title: Over-representation Analysis
#' subtitle: 'SAR: NA , Study: NA'
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
#' bibliography: references.bib
#' ---
#+ setup, include = FALSE
knitr::opts_chunk$set(echo = T, comment = "", message = F, warning = F, error = F)
options(width = 100)
#+ echo = F
library(here)
library(DT)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
#' # Over-representation Analysis
#' In this section, I will illustrate the principle of overpresentation analysis.
#' There exist a variety of tools and methods for ORA.
#' A few review papers for details:
#'
#'  - @huangBioinformaticsEnrichmentTools2009
#'
#'
#' ## Data prepartion
#' I will use the microarray data from smoking study as an example in this illustration.
#' The details of this study is documented in this [paper](https://link.springer.com/article/10.1007/s00251-010-0431-6).
#' The results (`rst` object) from differential gene expression analysis (DGEA)
#' will be used here.
#' Additionally, annotation file for gene IDs will be needed for R package based
#' functions.
#+ echo = F
rst <- readRDS(file.path(here(), "data/de_rst.rds")) %>%
  arrange(P.Value)
gene_ids <- readRDS(file.path(here(), "data/gene_id.rds"))
#'
#' ## Gene set selection
#'  Normally these gene sets should represent major players in biological
#' pathways. This can be achieved through the query from the biological
#' databases or from other analytical approach such as co-expression analysis.
#' For illustration purpose, here I randomly select 10 genes for gene set 1.

set.seed(123)
gs1 <- sample(rst$ID, 10)

#' For gene set 2, I choose five most significantly differentially expressed genes,
#' five genes not differentially expressed.
gs2 <- c(head(rst$ID, 5), tail(rst$ID, 5))

print(glue::glue('geneset 1: {paste0(gs1, collapse = ", ")}
geneset 2: {paste0(gs2, collapse =", ")}'))

#' ### Define differentially expressed gene
#'
#' Here I define differentially expressed genes with following criterion:
#'
#' - adjusted p-value < 0.001
#' - log fold change > 5
#'
#' These threshold criterion are not fixed can be varied depending on the context.
#' Or it can be obtained from other analytical methods. The definition of these
#' "over-represented" genes remains a controversial topic. However, it is out of
#' the scope of this illustration.
#'
rst <- rst %>%
  mutate(
    de = adj.P.Val < 0.001 & abs(logFC) > 5,
    gs1 = ID %in% gs1,
    gs2 = ID %in% gs2,
  )
table(rst[, c("de", "gs1")])
table(rst[, c("de", "gs2")])

#' ## Over-representation analysis
#'
#' > Over-representation analysis (ORA) is used to determine which a priori defined
#' > gene sets are more present (over-represented) in a subset of “interesting”
#' > genes than what would be expected by chance.
#' > -- @huangBioinformaticsEnrichmentTools2009.
#'
#' **The null hypothesis of ORA test:** being differentally expressed and
#' belong to a particular pathway (gene set) are independent.
#'
#' **Common statistical tests: ** chi-square test, Fisher's exact test. Note
#' that when counts in a cell of a contigency table is small, chi-square
#' approximation may not be accurate.
#'
#' The following examples tests the hypothesis for the gene sets described above
#' , e.g. if the genes in gene set 1 are differentially expressed between smokers
#' and non-smokers.
#'
#' **Gene set 1**

conttab_gs1 <- as.table(table(rst[, c("gs1", "de")]))
print(conttab_gs1)
chisq.test(conttab_gs1)
fisher.test(conttab_gs1, alternative = "greater")
#' **Gene set 2**
conttab_gs2 <- as.table(table(rst[, c("gs2", "de")]))
print(conttab_gs2)
chisq.test(conttab_gs2)
fisher.test(conttab_gs2, alternative = "greater")

#' ## Over-respensentation analysis from R package
#' The advantage of performing ORA using R packages `clusterProfile`
#' (which is the fusion of a few packages) is
#' that the functions combine the step of database query and the step of testing
#' in one go. The `clusterProfile` package currently provide
#' biological database indicated in the *Over-representation analysis* box in
#' the graph below.
#' ![curtesy of clusterProfiler](images/clusterProfiler-diagram.png)
#'
#' As shown in the graph, the ORA tests with different db queries are prefixed
#' with `enrich`, such as `enrichGO`, `enrichDO`. Note that these functions
#' are different than the test above. They first query the biological process
#' databases and find the genes lie in the pathway(s). Then it will perform
#' the test to see if the genes in the pathways are differentially expressed compared
#' to the background genes. For example, I'm interested in
#' the top 10 most significantly differentially expressed genes. I would like
#' to query the biological pathways that these genes involve in then test if the
#' gene in these pathways are significantly differentially expressed between
#' smokers and non-smokers.

gs3 <- head(rst$ID, 10)
gs3_ids <- gene_ids[which(gene_ids$hgnc_symbol %in% gs3), "entrezgene_id"] %>%
  as.character()
#' Gene IDs of geneset 3: `r gs3_ids`
#' **Note that not all gene symbols can be mapped to the gene IDs.**
#'
#' ## Test set up and result interpretations
#' To perform ORA in `clusterProfiler` one needs to provide following input.
#'
#' **Setup parameters:**
#'
#' - gene: a **character** vector of gene entrez ID in a gene set to be tested
#' - universe: a **character** vector of gene entrez ID in the background gene set. This background gene is defined by user. It could be all genes in the dataset or housekeeping genes or negative controls, depending on the purpose of comparison.
#' - ontology type: biological ontology type. This varies from different reference database. Check specific database for details.
#' - reference database: some function requires the input of database. one can
#' either use the connection to online db or download to local, takes big
#' storage space.
#'
#' **Interpretation of results:**
#'
#' I use the result from the following GO enrichment analysis for illustration. This applies
#' to the analytical results queried from other databases. From the result object,
#' one can see that the these top genes were mapped onto 62 gene sets (GO terms).
#' Out of these sets, 10 were tested significant. The significance cutoff set
#' up is illustrated in the result object.
#' To interpret the result table, take the examplle of GO term `GO:0033116` as
#' an example. This GO term contains 60 genes, out of these genes, only one
#' gene is found in our queried list, which contains 5 genes. So the `GeneRatio`
#' is 1/5. The pathway. In the background gene set (BgRatio), there are 60
#' out of 11599 genes. From this we can construct a contingency table below and
#' perform the Fisher test:
#'
#' || GO term| non-Go term|
#' | --- | --- | --- |
#' |enriched|1|4|
#' |background| 59 |11535|
tab <- as.table(matrix(c(1, 59, 4, 11535), ncol = 2))
fisher.test(tab, alternative = "greater")
#' ## Examples {.tabset}
#' I list the example queries in the mostly used databases.
#' If not specified, all the functions are from the `clusterProfile` package.
#'
#' ### GO
#' [Gene Ontology](http://geneontology.org/)

#'
#+ ora, cache = T
rst_go <- enrichGO(gs3_ids,
  universe = as.character(gene_ids$entrezgene_id),
  OrgDb = org.Hs.eg.db,
  ont = "CC",
  pAdjustMethod = "BH",
  readable = F
)

str(rst_go, max.level = 2)
datatable(rst_go@result,
  rownames = F,
  options = list(scrollX = "400px")
)
rst_go@geneSets["GO:0033116"]

#' ### KEGG
#' [Kyoto Encyclopedia of Genes and Genomes](https://www.genome.jp/kegg/)
#' **Note**
#'
#' I ran into bug from the older version of `clusterProfiler`.
#' In order to get this function work, update R version and upgrade `DOSE`
#'  package if you use online reference. Alternatively download KEGG
#' reference and set `use_internal_data = TRUE`.
#' I can't upgrade R in my local enviroment and don't have space to download
#' massive data. So I won't be able to illustrate KEGG here.
#+ cache = T, warning = T, error = T
rst_kegg <- enrichKEGG(gs3_ids,
  universe = as.character(gene_ids$entrezgene_id),
  minGSSize = 1,
  pvalueCutoff = 0.3,
  organism = "hsa",
  use_internal_data = F
)
#' KEGG module
#+ cache = T, error = T
rst_mkegg <- enrichMKEGG(gs3_ids,
  universe = as.character(gene_ids$entrezgene_id),
  organism = "hsa",
  keyType = 'ncbi-geneid'
)

str(rst_mkegg, max.levels = 1)
datatable(rst_mkegg@result,
  rownames = F,
  options = list(scrollX = "400px")
)

#' ### WikiPathways
#' [WikiPathways](https://www.wikipathways.org/)
#+ cache = T, error = T
enrichWP(gs3_ids,
  organism = "Homo sapiens",
  universe = as.character(gene_ids$entrezgene_id)
)

#' ### Reactome
#' [Reactome](https://reactome.org/)
#+ cache = T, error = T
rst_react <- ReactomePA::enrichPathway(gs3_ids,
  universe = as.character(gene_ids$entrezgene_id),
  qvalueCutoff = 0.1
)

str(rst_react, max.level = 2)
datatable(rst_react@result,
  rownames = F,
  options = list(scrollX = "400px")
)

#' ### DO
#' [Disease Ontology](https://disease-ontology.org/)
#+ cache = T, error = T
DOSE::enrichDO(gs3_ids,
  universe = as.character(gene_ids$entrezgene_id)
)

#' ### DisGeNET
#' [DisGeNET](https://www.disgenet.org/)
#+ cache = T, error = T
rst_dgn <- DOSE::enrichDGN(gs3_ids,
  universe = as.character(gene_ids$entrezgene_id),
  readable = F
)
str(rst_dgn, max.level = 2)
datatable(rst_dgn@result,
  rownames = F,
  options = list(scrollX = "400px")
)

#' ### MeSH
#' [MeSH](https://www.ncbi.nlm.nih.gov/mesh/)
#'
#' The enrichement test function requires the feed the reference argument. One
#' can either download the whole reference or connect on online db.
#'  The following code  will create a local database in the cache directory.
#' The location of the local cache can be found (and updated)
#' with `getAnnotationHubCache` and `setAnnotationHubCache`;
#' `removeCache` removes all cache resources.
#+  eval = F, error = T
annot <- AnnotationHub::AnnotationHub(localHub = F)
print(annot)
# query annotation hub for MeSH database
hsa <- AnnotationHub::query(annot, c("MeSHDb", "Homo sapiens"))
print(str(hsa))
# download rosources to local cache directory
ref_hsa <- hsa[[1]]
MeSHDB <- MeSHDbi::MeSHDb(ref_hsa)
#' Otherwise one can install `MeSH.Hsa.eg.db` package in the local environment.
#+ cache = T, error = T
# require installation of `MeSH.Hsa.eg.db` package.
# replace `MeSHDb = MeSHDB` if using local reference database.
rst_mesh <- meshes::enrichMeSH(gs3_ids,
  MeSHDb = "MeSH.Hsa.eg.db"
)
str(rst_mesh, max.level = 2)
datatable(rst_mesh@result,
  rownames = F,
  options = list(scrollX = "400px")
)


# library(RDAVIDWebService)
# enrichDAVID(gs3_ids,
# universe = as.character(gene_ids$entrezgene_id)
# )
#' ## Pros and Cons of ORA
#'
#' - Simple and straitforward
#' - can be used for exploratory analysis when there's no prior hypothesis about
#' any biological mechanism.
#' - However, can not infer the directly of differentially expressed pathway
#' e.g. up-regulated or down-regulated. This is not so helpful when particularly
#' interested in the pathways that a drug modulates.
#' - assumes the genes in the test are iid.
#'
#' ## Reference
/*
#' <details><summary>Session Info</summary>
sessionInfo()
#' </details>
#+ echo = F, eval = F
# Markdown --------------------------------------------------------
# rmarkdown::render('src/04_ORA.R', output_dir = 'output')
# knitr::spin("src/04_ORA.R", format = "Rmd", knit = F)
*/