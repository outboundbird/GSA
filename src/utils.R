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
