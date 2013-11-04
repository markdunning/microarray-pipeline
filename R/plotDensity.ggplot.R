plotDensity.ggplot <-
function (mat, ylab = "density", xlab = "log2 Intensity", type = "l", col = 1:6, 
na.rm = TRUE, ...)
{

df <- melt(data.frame(mat))

p <- ggplot(df, aes(x = value, fill=factor(variable)),alpha=0.2) + geom_density()

p

}
