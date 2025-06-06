#' ---
#' title: Power analysis for RNAseq differential expression analysis
#' subtitle: 'SAR: NA , Study: NA'
#' author:  Siying Huang (E0482362), Biomarker statistics team
#' date: 'created: 2023-06-30 , updated (`r Sys.Date()`)'
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
knitr::opts_knit$set(root_dir = "/mnt/c/Users/e0482362/Work/pathway_analysis/src")
knitr::opts_chunk$set(echo = T, comment = "", message = F, warning = F, error = F)
options(width = 100)
#+ libs
library(here)

#' # Suggested papers
#' - [Statistical Power Analysis for Designing Bulk, Single-Cell, and Spatial Transcriptomics Experiments: Review, Tutorial, and Perspectives](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9952882/)
#' - [Power analysis and sample size estimation for RNA-Seq differential expression](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4201821/)
#' - [Power analysis for RNA-Seq differential expression studies]()
#' - [powsimR: power analysis for bulk and single cell RNA-seq experiments](https://academic.oup.com/bioinformatics/article/33/21/3486/3952669)
#' - [Dedicated transcriptomics combined with power analysis lead to functional understanding of genes with weak phenotypic changes in knockout lines](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008354)
#' - [ssizeRNA](https://cran.r-project.org/web/packages/ssizeRNA/vignettes/ssizeRNA.pdf)
#' -[Tong T, Zhao H. Practical guidelines for assessing power and false discovery rate for a fixed sample size in microarray experiments. Stat Med. 2008 May 20;27(11):1960-72. doi: 10.1002/sim.3237. PMID: 18338314; PMCID: PMC3157366.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3157366/)
#' - [Yudi Pawitan et al., “False Discovery Rate, Sensitivity and Sample Size for Microarray Studies,” Bioinformatics 21, no. 13 (July 1, 2005): 3017–24, https://doi.org/10.1093/bioinformatics/bti448.](https://academic.oup.com/bioinformatics/article/21/13/3017/196891)
#' - [sizepower](https://www.bioconductor.org/packages/devel/bioc/vignettes/sizepower/inst/doc/sizepower.pdf)
#' - [RNASeqPower](http://bioconductor.org/packages/release/bioc/html/RNASeqPower.html)
#' - [MOPower](https://github.com/HSyed91/MOPower)
#' - [Sonia Tarazona et al., “Harmonization of Quality Metrics and Power Calculation in Multi-Omic Studies,” Nature Communications 11, no. 1 (June 18, 2020): 3092, https://doi.org/10.1038/s41467-020-16937-8.](https://www.nature.com/articles/s41467-020-16937-8)
#'      - [POWSC](https://github.com/suke18/POWSC)
#'
#' ## scRNAseq
#' - [Kenong Su, Zhijin Wu, and Hao Wu, “Simulation, Power Evaluation and Sample Size Recommendation for Single-Cell RNA-Seq,” Bioinformatics 36, no. 19 (July 2, 2020): 4860–68, https://doi.org/10.1093/bioinformatics/btaa607.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7824866/)
#' # Compare power calcualtion results from different R packages
#'
#' # Sample size consideration for multi-omic studies
#'
#' - [David A. A. Baranger et al., “Multi-Omics Cannot Replace Sample Size in Genome-Wide Association Studies,” Genes, Brain and Behavior n/a, no. n/a: e12846, accessed July 5, 2023, https://doi.org/10.1111/gbb.12846.](https://onlinelibrary.wiley.com/doi/full/10.1111/gbb.12846)
#' - [Alexander Kirpich et al., “Variable Selection in Omics Data: A Practical Evaluation of Small Sample Sizes,” PLOS ONE 13, no. 6 (June 21, 2018): e0197910, https://doi.org/10.1371/journal.pone.0197910.](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0197910)
#' <details><summary>Session Info</summary>
sessionInfo()
#' </details>
#+ echo = F, eval = F
# Markdown --------------------------------------------------------
# rmarkdown::render('src/power_de.R', output_dir = 'output')
# knitr::spin('src/power_de.R', format = 'Rmd', knit = FALSE)