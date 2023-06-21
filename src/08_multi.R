#' ---
#' title: Multiplicity issue
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
#' ---
#+ setup, include = FALSE
knitr::opts_chunk$set(echo = T, comment = "", message = F, warning = F, error = F)
options(width = 100)

#' # Multiplicity issue
#' Gene set analysis or pathway analysis usually contains multiple pathways or gene sets. This poses a question of adjustment of multiple testing. Conventially, plenty of methods have been developed to address this issue. However the unique character in the pathway analysis is that these pathways are often interrelated to each other. Thus each of the pathway cannot be considered as "independent" entity. This makes us to reconsiderate the way of setting up threshold for controlling false positive results.
#'
#' ## Reference
#' [@benjaminiSelectiveInferenceMultiple](https://academic.oup.com/jrsssb/article/76/1/297/7075946)
#' @panEffectsThresholdChoice2005
#' 1.Pan KH, Lih CJ, Cohen SN. Effects of threshold choice on biological conclusions reached during analysis of gene expression by DNA microarrays. Proceedings of the National Academy of Sciences [Internet]. 2005 Jun 21 [cited 2023 Jun 16];102(25):8961–5. Available from: https://www.pnas.org/doi/full/10.1073/pnas.0502674102




#' <details><summary>Session Info</summary>
sessionInfo()
#' </details>
#+ echo = F, eval = F
# Markdown --------------------------------------------------------
# rmarkdown::render('src/08_multi.R', output_dir = 'output')
# knitr::spin('src/08_multi.R', format = 'Rmd', knit = FALSE)