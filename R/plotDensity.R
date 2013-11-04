plotDensity <-
function (mat, ylab = "density", xlab = "log2 Intensity", type = "l", col = 1:6, 
na.rm = TRUE, ...) 
{
   x.density <- apply(mat, 2, density, na.rm = na.rm)
    all.x <- do.call("cbind", lapply(x.density, function(x) x$x))
    all.y <- do.call("cbind", lapply(x.density, function(x) x$y))
    matplot(all.x, all.y, ylab = ylab, xlab = xlab, type = type, 
			col = col, ...)
    invisible(list(all.x = all.x, all.y = all.y))
}
