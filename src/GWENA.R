#' ---
#' title: Gene co-Expression Network Analysis pipeline
#' subtitle: 'SAR: sar , Study: study'
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
#' ---
#+ setup, include = FALSE
# knitr::opts_knit$set(root_dir='/mnt/c/Users/e0482362/Work/pathway_analysis/src')
knitr::opts_chunk$set(echo = T, comment = "", message = F, warning = F, error = F)
options(width = 100)
#+ libs
library(here)
library(dplyr)
library(GWENA)
#+ io, cache = T
# io -----------------------------------------

expr <- readRDS(here("data/expr.rds"))

expr <- expr[!duplicated(expr$IDENTIFIER), ] %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames("IDENTIFIER") %>%
  select(-ID_REF)

# expr <- expr[1:200, ]

pheno <- readRDS(here("data/pheno.rds")) %>%
  mutate(
    status = as.integer(stress),
    stress = as.character(stress)
  ) %>%
  tibble::column_to_rownames("sample")
str(pheno)

#' # Gene co-expression analysis pipeline
#' ![](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8152313/bin/12859_2021_4179_Fig1_HTML.jpg)
#'
#' [vignette](https://bioconductor.org/packages/release/bioc/vignettes/GWENA/inst/doc/GWENA_guide.html)
#'
#' [pipeline paper, Lemoine et al, 2021](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8152313/)
#' [here insert description of the co-expression analysis summary]
#' Snapshot of the expression data.
SummarizedExperiment::SummarizedExperiment(
  assays = list(expr = expr),
  colData = S4Vectors::DataFrame(pheno)
)
#' ## Gene filtering
#'   use `filter_low_var` or `filter_RNA_seq` functions.
expr_filtered <- filter_low_var(t(expr))
dim(expr_filtered)
#'
#' ## Network building
#'
#' > Gene co-expression networks are an ensemble of genes (nodes) linked to each other (edges) according to the strength of their relation. In GWENA, this strength is estimated by the computation of a **(dis)similarity score** which can start with a distance (euclidian, minkowski, ...) but is usually a correlation. Among these, Pearson's one is the most popular, however in GWENA we use Spearman correlation by default. It is less sensitive to outliers which are frequent in transcriptomics datasets and does not assume that the data follows the normal distribution.
#' > The co-expression network is built according to the following sub-steps :
#' >
#' > 1. A **correlation** (or distance) between each pair of genes is computed.
#' > 2. A **power law** is fitted on the correlation matrix. This step can be performed by itself through the function `get_fit.expr` if needed.
#' > 3. An **adjacency score** is computed by adjusting previous correlations by the fitted power law.
#' > 4. A **topological overlap score** is computed by accounting for the network's topology.
#'
#' **Caution**: this process takes time!
#+ build-network, cache = T
net <- build_net(expr_filtered,
  cor_func = "spearman",
  n_threads = 6
)
# Power selected :
net$metadata$power
# Fit of the power law to data ($R^2$) :
fit_power_table <- net$metadata$fit_power_table
fit_power_table[fit_power_table$Power == net$metadata$power, "SFT.R.sq"]

#' ## Module detection
#' > At this point, the network is a [complete graph](https://en.wikipedia.org/wiki/Complete_graph): all nodes are connected to all other nodes with different strengths. Because gene co-expression networks have a scale free property, groups of genes are strongly linked with one another. In co-expression networks these groups are called **modules** and assumed to be representative of genes working together to a common set of functions.
#' >
#' > Such modules can be detected using unsupervised learning or modeling. GWENA use the hierarchical clustering but other methods can be used (kmeans, Gaussian mixture models, etc.).

modules <- detect_modules(expr_filtered,
  net$network,
  detailled_result = TRUE,
  merge_threshold = 0.25
)

#' **Important**: Module 0 contains all genes that did not fit into any modules.
#' > Since this operation tends to create multiple smaller modules with highly similar expression profile (based on the [eigengene](FAQ.html/#what-is-an-eigengene) of each), they are usually merged into one.
#'
# Number of modules before merging :
length(unique(modules$modules_premerge))
# Number of modules after merging:
length(unique(modules$modules))

layout_mod_merge <- plot_modules_merge(
  modules_premerge = modules$modules_premerge,
  modules_merged = modules$modules
)

ggplot2::ggplot(
  data.frame(modules$modules %>% stack()),
  ggplot2::aes(x = ind)
) +
  ggplot2::stat_count() +
  ggplot2::ylab("Number of genes") +
  ggplot2::xlab("Module")

#' > Each of the modules presents a distinct profile, which can be plotted in two figures to separate the positive (+ facet) and negative (- facet) correlations profile. As a summary of this profile, the eigengene (red line) is displayed to act as a signature.
#' >
#+ plot, cache = F
plot_expression_profiles(expr_filtered, modules$modules)
#'
#' ## Biological integration
#'
#' ### Functional enrichment
#' > A popular way to explore the modules consists of linking them with a known biological function by using currated gene sets. Among the available ones, [Gene Ontology (GO)](https://www.nature.com/articles/srep18871#ref-CR1), [Kyoto Encyclopedia of Genes and Genomes (KEGG)](https://www.nature.com/articles/srep18871#ref-CR2), [WikiPathways](https://doi.org/10.1093/nar/gkx1064), [Reactome](https://www.nature.com/articles/srep18871#ref-CR4), [Human Phenotype Ontology (HPO)](https://doi.org/10.1093%2Fnar%2Fgkt1026) put modules into a broader systemic perspective.
#' >
#' > In oppositions, databases references like [TRANSFAC](https://doi.org/10.1093%2Fbib%2Fbbn016), [miRTarBase](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3013699), [Human Protein Atlas (HPA)](https://dx.doi.org/10.1126%2Fscience.1260419), and [CORUM](http://www.ncbi.nlm.nih.gov/pubmed/30357367) give more details about tissue/cell/condition information.
#' >
#' > Using the [over-representation analysis (ORA)](http://doi.org/10.1371/journal.pcbi.1002375) tool [GOSt](https://biit.cs.ut.ee/gprofiler/gost) from g:Profiler, we can retrieve the biological association for each module and plot it as follows.
#' >
#' ### Phenotypic association
#' > If phenotypic information is available about the samples provided, an association test can help to determine if a module is specifically linked to a trait. In this case, module 1 seems to be strongly linked to `Age`.
#' >
phenotype_association <- associate_phenotype(
  modules$modules_eigengenes,
  pheno[, c("stress", "status")]
)

plot_modules_phenotype(phenotype_association)

#' ## Graph visualization and topological analysis
#'
#' > Information can be retrieved from the network topology itself. For example, hub genes are highly connected genes known to be associated with key biological functions. They can be detected by different methods :
#' >
#' > * `get_hub_high_co`: Highest connectivity, select the top n (n depending on parameter given) highest connected genes. Similar to WGCNA's selection of hub genes
#' > * `get_hub_degree`: Superior degree, select genes whose degree is greater than the average connection degree of the network. Definition from network theory.
#' > * `get_hub_kleinberg`: Kleinberg's score, select genes whose Kleinberg's score is superior to the provided threshold.
#' >
#' > Manipulation of graph objects can be quite demanding in memory and CPU usage. Caution is advised when choosing to plot networks larger than 100 genes.
#' > Since co-expression networks are complete graphs, readability is hard because all genes are connected with each other. In order to clarity visualization, edges with a similarity score below a threshold are removed.
#+ net-vis
module_example <- modules$modules$`0`
graph <- build_graph_from_sq_mat(net$network[module_example, module_example])

layout_mod_2 <- plot_module(
  graph,
  # upper_weight_th = 0.999,
  # vertex.label.cex = 0,
  # node_scaling_max = 7,
  # legend_cex = 1,
  # zoom = 1.5
)

#' > As modules also follow a modular topology inside, it may be interesting to detect the sub clusters inside them to find genes working toward the same function through enrichment. The sub cluster can then be plotted on the graph to see their interaction.
#+ sub-clust, cache =F
net_mod_1 <- net$network[modules$modules$`1`, modules$modules$`1`]
sub_clusters <- get_sub_clusters(net_mod_1)
layout_mod_2_sub_clust <- plot_module(
  graph,
  # upper_weight_th = 0.999995,
  # groups = sub_clusters,
  # vertex.label.cex = 0,
  # node_scaling_max = 7,
  # legend_cex = 1
)
#' ## Networks comparison
#' > A co-expression network can be built for each of the experimental conditions studied (e.g. control/test) and then be compared with each other to detect differences of patterns in co-expression. These may indicate breaks of inhibition, inefficiency of a factor of transcription, etc. These analyses can focus on preserved modules between conditions (e.g. to detect housekeeping genes), or unpreserved modules (e.g. to detect genes contributing to a disease).
#' >
#' > GWENA uses a comparison test based on random re-assignment of gene names inside modules to see whether patterns inside modules change (from [NetRep](https://cran.r-project.org/web/packages/NetRep/vignettes/NetRep.html) package). This permutation test is repeated a large number of times to evaluate the significance of the result obtained.
#' > To perform the comparison, all previous steps leading to modules detection need to be done for each condition. To save CPU, memory and time, the parameter `keep_cor_mat` from the `build_net` function can be switched to TRUE so the similarity matrix is kept and can be passed to `compare_conditions`. If not, the matrix is re-computed in `compare_conditions`.
#+ net-compare
# Expression by condition with data.frame/matrix
#+ net_compare, eval = F
samples_by_cond <- lapply(levels(as.factor(pheno$stress)), function(cond) {
  df <- pheno %>%
    dplyr::filter(stress == cond) %>%
    dplyr::select(TRTLES1)
  apply(df, 1, paste, collapse = "_")
}) %>% setNames(meta$LESION %>% unique())

expr_by_cond <- lapply(samples_by_cond %>% names(), function(cond) {
  samples <- samples_by_cond[[cond]]
  t(expr_cut[, which(colnames(expr_cut) %in% names(samples))])
}) %>% setNames(samples_by_cond %>% names())

expr_by_cond[[2]] %>%
  is.na() %>%
  any()
#+ build-net, cache = F, eval = F
# Network building and modules detection by condition
bound_spearman <- function(x) {
  c <- cor(x, method = "spearman", use = "pairwise")
  c[c > 1] <- 1
  c[c <- 1] <- 1
  c
}

net_by_cond <- lapply(expr_by_cond,
  build_net,
  your_func = bound_spearman,
  cor_func = "other",
  n_threads = 4,
  stop_if_fit_pb = F,
  keep_matrices = "both"
)


mod_by_cond <- mapply(detect_modules, expr_by_cond,
  lapply(net_by_cond, `[[`, "network"),
  MoreArgs = list(detailled_result = TRUE),
  SIMPLIFY = FALSE
)

comparison <- compare_conditions(expr_by_cond,
  lapply(net_by_cond, `[[`, "adja_mat"),
  lapply(net_by_cond, `[[`, "cor_mat"),
  lapply(mod_by_cond, `[[`, "modules"),
  pvalue_th = 0.05
)
str(comparison)
comparison$result$NORMAL$LESION$comparison
plot_comparison_stats(comparison$result$NORMAL$LESION$p.values)

enrichment <- bio_enrich(modules$modules)
plot_enrichment(enrichment)


#' <details><summary>Session Info</summary>
sessionInfo()
#' </details>
#+ echo = F, eval = F
# Markdown --------------------------------------------------------
# rmarkdown::render('src/GWENA.R', output_dir = 'output')
# knitr::spin('src/GWENA.R', format = 'Rmd', knit = FALSE)