#' ---
#' title: Differential coexpression analysis
#' subtitle: 'SAR: NA , Study: tutorial'
#' author:  Siying Huang (E0482362), Biomarker statistics team
#' date: 'created: 2023-06-26 , updated (`r Sys.Date()`)'
#' always_allow_html: true
#' output:
#'   html_document:
#'     css: sanofi.css
#'     code_folding: "hide"
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
source(file.path(here(), "src/utils.R"))

#+ io, cache = T
rst <- readRDS(file.path(here(), "data/de_rst.rds"))
expr <- readRDS(file.path(here(), "data/expr.rds"))
gene_ids <- readRDS(file.path(here(), "data/gene_id.rds"))

ctrl_ids <- c("GSM464697", "GSM464699", "GSM464705", "GSM464702", "GSM464715")
case_ids <- c("GSM464678", "GSM464688", "GSM464683", "GSM464692", "GSM464696")

expr <- expr[!duplicated(expr$IDENTIFIER), ] %>%
  data.frame() %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames("IDENTIFIER") %>%
  dplyr::select_at(c(ctrl_ids, case_ids)) %>%
  mutate_if(is.numeric, log2)

set.seed(123)
gs1 <- sample(rst$ID, 10)
gs3 <- head(rst$ID, 10)
n_bg <- 40
bg <- sample(rst$ID, n_bg)
gs <- Reduce(union, list(gs1, gs3, bg))
expr_case <- expr[gs, case_ids]
expr_ctrl <- expr[gs, ctrl_ids]
expr_g <- expr[gs, c(case_ids, ctrl_ids)]
# dim(expr_case)
# dim(expr_ctrl)

#' Gene set coexpression analysis and gene set net correlation analysis are classified as
#' topology-based methods. This type of analysis addresses the issue that not all
#' genes play equally important role in a biological process. Therefore it considers
#' the correlations among genes.
#'
#' # Gene Set Coexpression Analysis (GSCA)
#' [Gene set coexpression analysis](https://www.biostat.wisc.edu/~kendzior/GSCA/) was first proposed by @choiStatisticalMethodsGene2009a.
#'
#' ## Scheme
#'
#' ![](images/gsca.png){width="344"}
#' T1 and T2 represent two biological conditions, Nk represents the number of
#' arrays in condition k, k= 1,2. The dispersion index is given by the Edulidean distance, adjusted for the size of the gene set considered:
#'
#' $$Ds(\rho_c^{T1}, \rho_c^{T2})= \sqrt{\frac{1}{P_c} \sum_{p=1}^{P_c}(\tilde{\rho}_p^{T_1, T_2})^2} $$
#' 
#+ cache = T
library(GSCA)
#' # Gene Sets Net Correlation Analysis (GSNCA)
#' GSNCA method was proposed by @rahmatallahGeneSetsNet2014.
#'
#' ## Scheme
#' ![](images/gsnca.jpg){width="394"}
#' > in the gene coexpression network between two conditions. Importantly, we do not infer ‘gene coexpression networks’ explicitly, but, instead, we estimate net correlation changes by introducing for each gene a weight factor that characterizes its cross-correlations in the coexpression networks. Weight vectors in both conditions are found as eigenvectors of correlation matrices with zero diagonal elements. The Gene Sets Net Correlations Analysis (GSNCA) tests the hypothesis that for a gene set there is no difference in the gene weight vectors between two conditions.
#'

library(GSAR)
library(igraph)

B1 <- cor(t(expr_case[gs3, ]))
rst1_igen <- eigen(B1)
l11 <- sum(abs(rst1_igen$values))
w1 <- rst1_igen$values / l11

B0 <- cor(t(expr_ctrl[gs3, ]))
rst0_igen <- eigen(B0)
l10 <- sum(abs(rst0_igen$values))
w0 <- rst0_igen$values / l10

sum(abs(w1 - w0))

# minimum span tree on genes
B <- cor(t(expr_g))
dst <- 1 - abs(B)
rst_mst <- ape::mst(dst)
plot(rst_mst, graph = "nsca")
str(rst_mst)
which(rst_mst["MYDGF", ] == 1)

mst_g <- findMST2(as.matrix(expr_g))
str(mst_g)

plotMST2.pathway(as.matrix(expr_g[gs3, ]),
  group = rep(c(1, 2), each = 5)
)
dim(expr_g)

# minimum spanning tree on samples
B <- cor(expr_g)
dst <- 1 - abs(B)
rst_mst <- ape::mst(dst)
head(rst_mst)
rst_pca <- prcomp(expr_g)
plot(rst_mst,
  x1 = rst_pca$x[, 1],
  x2 = rst_pca$x[, 2],
  # graph = "circle"
)

gr <- igraph::graph.adjacency(dst)
mst_ <- igraph::minimum.spanning.tree(gr)

#' # Differential gene set coexpression analysis in R
#'
library(GSAR)
labs <- rep(c(1, 2), each = 5)
WWtest(as.matrix(expr_g[gs1, ]), labs)
WWtest(as.matrix(expr_g[gs3, ]), labs)

GSNCAtest(as.matrix(expr_g[gs1, ]), labs)
GSNCAtest(as.matrix(expr_g[gs3, ]), labs)

plotMST2.pathway(as.matrix(expr_g[gs1, ]), labs)
plotMST2.pathway(as.matrix(expr_g[gs3, ]), labs)
#' ## Reference
/*
#' <details><summary>Session Info</summary>
sessionInfo()
#' </details>
#+ echo = F, eval = F
# Markdown --------------------------------------------------------
# rmarkdown::render('src/DCA.R', output_dir = 'output')
# knitr::spin('src/DCA.R', format = 'Rmd', knit = FALSE)
*/
