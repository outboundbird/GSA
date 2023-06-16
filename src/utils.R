library(ggplot2)
library(dplyr)
volc_plot <- function(data,
           lab_method = "p",
           lab_col = "x",
           p_thresh = 0.05,
           fdr_thresh = 0.1,
           log_cf = FALSE,
           show_text= TRUE,
           ...) {
    switch(lab_method,
      p = {
        idx <- data[, "p.value"] < p_thresh
      },
      fdr = {
        idx <- data[, "FDR"] < fdr_thresh
      }
    )

    pch_idx <- ifelse(idx, 16, 1)
    col_idx <- ifelse(idx, 2, "gray")

    if (lab_col == "rownames") {
      labels <- ifelse(idx, rownames(data), "")
    } else {
      labels <- ifelse(idx, as.character(data[, lab_col]), "")
    }

    if (log_cf) {
      est <- log2(data[, "estimate"])
      xlab_text <- "log_2(Effect estimates)"
    } else {
      est <- data[, "estimate"]
      xlab_text <- "Effect estimates"
    }

    plot(
      est,
      -log10(data[, "p.value"]),
      xlab = xlab_text,
      ylab = "-log10(P value)",
      col = col_idx,
      pch = pch_idx,
      ...
    )
    abline(
      h = -log10(p_thresh),
      col = rgb(0, 0, 1, 0.5),
      lty = 2
    )
    abline(
      v = 0,
      col = rgb(0, 0, 0, 0.5),
      lty = 2
    )
    if (show_text) {
          text(data[, "estimate"], -log10(data[, "p.value"]),
            labels,
            pos = 2,
            cex = 0.5
          )

    }

  }

get_args <- function(args, key, default) {
  if (is.null(args[[key]])) {
    return(default)
  }
  args[[key]]
}

#' ks_test_plot
#'
#' Perform Kolmogorov-Smirnov non-parametric test. A warapper of the stat:ks.test(), produce visualization of the tests
#' "" """"
#' @param data data frame that contains a column of genesymbol name, fold change.
#' @param genename str, colname of gene symbol.
#' @param fc str, colname of fold change, should not be on the log scale.
#' @param pathway_genes a character vector indicating genes of interest, which should be contained in the data set.
#' @param pathway_name str, the name to put on the plot.
#' @param xlab, str, x-axis label.
#' @param ... other keyword argument pass to ks.test
#' @mexclude boolean, logical; indicate if the gene set and background gene set test should be mutually exclusive. TRUE, correspond to GSEA KS test, FALSE (default) correspond to GSVA KS test.
#' @param plot logical, if TRUE (default) will return ggplot, otherwise, will return the result of ks.test.
#'
#' @return RETURN_DESCRIPTION
ks_test_plot <- function(data, genename, fc, pathway_genes,
                         pathway_name, ..., xlab = "Fold Change", mexclude = F,
                         plot = TRUE) {
  if ("tbl" %in% class(data[, 1])) {
    stop("The data cannot be tbl_df, tbl objects")
  }

args <- list(...)
ks_exact <- get_args(args, "exact", T)

  rank_fc <- arrange(data, .data[[fc]]) %>%
    select_at(c(genename, fc))
  idx <- rank_fc[, genename] %in% pathway_genes
  genelist <- rank_fc[idx, ]
  if (mexclude) {
    rank_fc <- rank_fc %>%
      filter(!ID %in% pathway_genes)
  }
  ks_rst <- ks.test(genelist[, fc], rank_fc[, fc], exact = ks_exact)

  if (!plot) {
    return(ks_rst)
  }

  ggplot(rank_fc, aes_string(x = fc)) +
    stat_ecdf(geom = "step", linetype = "dashed") +
    stat_ecdf(
      data = data.frame(genelist), aes_string(x = fc),
      color = "#6f00ff"
    ) +
    geom_hline(
      yintercept = 0.5,
      color = "#6f00ff57", linetype = "dashed"
    ) +
    labs(
      x = xlab, y = "CDF",
      title = "Kolmogorov-Smirnov test",
      subtitle = pathway_name,
      caption = glue::glue(
        "P value: ", format(ks_rst$p.value, digits = 3, scientific = T),
        " D= ", round(unname(ks_rst$statistic), 3)
      )
    ) +
    theme_bw()
}


#' Calculate enrichment score
#'
#' Calculate enrichment score for GSEA
#'
#' @param data data.frame object with DEGA result.
#' @param gs a string vector of gene set.
#' @param gid_col string or integer to index the the column in the data frame to match against gene set content. Can be gene symbol, gene ID, etc
#' @param corr_col string or integer index to indicate the correlation column to be used in the calculation.
#' @param method methods for calculation, index by the first authors: mootha - Mooth et al, 2003; 'sub' for Subramanian et al 2005
#' @param plot boolean, if TRUE plot the walking plot.
#'
#' @return NULL
calc_enrich_score <- function(data, gs, gid_col, corr_col,
method ='sub', plot= TRUE) {
  n_tot <- nrow(data)
  n_s <- length(gs)
  nr <- sum(abs(data[which(data[, gid_col] %in% gs), corr_col]))

  p_hit <- switch(method,
    mootha = sqrt((n_tot - n_s) / n_s),
    sub = 1 / nr
  )

  p_miss <- switch(method,
    mootha = -sqrt(n_s / (n_tot - n_s)),
    sub = -1 / (n_tot - n_s)
  )

  data$gs_idx <- data[, gid_col] %in% gs

  data <- switch(method,
    mootha = {
      data %>%
        mutate(score = if_else(gs_idx, p_hit, p_miss)) %>%
        mutate(
          score = case_when(row_number() == 1 ~ 0, TRUE ~ score),
          es = cumsum(score)
        )
    },
    sub = {
      data$abs_r <- abs(data[, corr_col])
      data %>%
        mutate(score = if_else(gs_idx, abs_r * p_hit, p_miss)) %>%
        mutate(
          score = case_when(row_number() == 1 ~ 0, TRUE ~ score),
          es = cumsum(score)
        )
    }
  )
  if (plot) {
    plot(data$es, type = "l", ylab = "Enrichment score")
    abline(h = 0, lty = 2, col = "red")
    axis(
      side = 3, at = which(data$gs_idx),
      col.ticks = "red",
      labels = F
    )
  }
}

# GSVA
compute_rank_score <- function(sort_idx_vec ) {
  num_genes <- length(sort_idx_vec)
  tmp <- rep(0, num_genes)
  tmp[sort_idx_vec] <- abs(seq(from = num_genes, to = 1) - num_genes / 2)
  return(tmp)
}

sort_idx_vec <- c(1, 3, 4, 2)
compute_rank_score(sort_idx_vec)
