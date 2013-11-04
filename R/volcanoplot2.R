volcanoplot2 <-
function (fit, coef = 1, highlight = 0, names = fit$genes$ID,
                          xlab = "Log Fold Change", ylab = "Log Odds", pch = 16, cex = 0.35,
                          ...){
###A small modification of the limma volcanoplot to add scatterSmooth 
  if (!is(fit, "MArrayLM"))
    stop("fit must be an MArrayLM")
  if (is.null(fit$lods))
    stop("No B-statistics found, perhaps eBayes() not yet run")
  x <- as.matrix(fit$coef)[, coef]
  y <- as.matrix(fit$lods)[, coef]
  smoothScatter(x, y, xlab = xlab, ylab = ylab, pch = pch, cex = cex,...)
  if (highlight > 0) {
    if (is.null(names))
      names <- 1:length(x)
    names <- as.character(names)
    o <- order(y, decreasing = TRUE)
    i <- o[1:highlight]
    text(x[i], y[i], labels = substring(names[i], 1, 8), cex = 0.8, col = "blue")
  }
  invisible()
}
