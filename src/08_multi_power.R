#' ---
#' title: Multiplicity issue and power analysis
#' subtitle: 'SAR: NA , Study: NA'
#' author:  Siying Huang (E0482362), Biomarker statistics team
#' date: 'created: 2023-06-19 , updated (`r Sys.Date()`)'
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

#' In this section, I talk about the null hypothesis, multiple comparison issue and statistical power for pathway analysis.
#'
#' # Null hypothesis
#' Here I cited the concepts of different types of null hypothesis from @deleeuwStatisticalPropertiesGeneset2016, @goemanAnalyzingGeneExpression2007 [Goeman et al](https://academic.oup.com/bioinformatics/article/23/8/980/198511?login=true#401164258)
#'
#'
#' || Differentially expressed gene| non-differentially expressed gene| Total|
#' | --- | --- | --- | --- |
#' |In gene set| $m_{gD}$| $m_{g\bar{D}}$| $m_G$|
#' |Not in gene set| $m_{\bar{g}D}$| $m_{\bar{G}\bar{D}}$| $m_{\bar{G}}$|
#' |Total| $m_D$| $m_{\bar{D}}$| $m$|
#'
#' ## Competitive null hypothesis
#' > competitive GSA considers all genes in the data, testing the null hypothesis that the genes in the gene set are no more strongly associated with the phenotype than other genes.
#'
#' ## Self-contained null hypothesis
#' > Self-contained GSA considers only the genes in the gene set and tests the null hypothesis that none of those genes are associated with the phenotype. @deleeuwStatisticalPropertiesGeneset2016
#'
#' ## Hybrid null hypothesis
#'
#' # Multiplicity issue
#' Gene set analysis or pathway analysis usually contains multiple pathways or gene sets. This poses a question of adjustment of multiple testing. Conventially, plenty of methods have been developed to address this issue. However the unique character in the pathway analysis is that these pathways are often interrelated to each other. Thus each of the pathway cannot be considered as "independent" entity. This makes us to reconsiderate the way of setting up threshold for controlling false positive results.
#'
#' # Statistical power estimation
#'
#' ## Reference
#' [Selective inference](https://academic.oup.com/jrsssb/article/76/1/297/7075946)
#'
/*#' <details><summary>Session Info</summary>
sessionInfo()
#' </details>
#+ echo = F, eval = F
# Markdown --------------------------------------------------------
# rmarkdown::render('src/08_multi_power.R', output_dir = 'output')
# knitr::spin('src/08_multi.R', format = 'Rmd', knit = FALSE)*/