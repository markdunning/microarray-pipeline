\name{geneHeatmap}
\alias{geneHeatmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make heatmap for a set of genes
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
geneHeatmap(dbcon, genes, annoName = NULL, outfile = "Results.html", outfile_path = ".")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dbcon}{
%%     ~~Describe \code{dbcon} here~~
}
  \item{genes}{
%%     ~~Describe \code{genes} here~~
}
  \item{annoName}{
%%     ~~Describe \code{annoName} here~~
}
  \item{outfile}{
%%     ~~Describe \code{outfile} here~~
}
  \item{outfile_path}{
%%     ~~Describe \code{outfile_path} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (dbcon, genes, annoName = NULL, outfile = "Results.html", 
    outfile_path = ".") 
{
    dir.create(outfile_path)
    library(hwriter)
    outPage = openPage(filename = outfile)
    if (is.null(annoName)) {
        meta <- dbGetQuery(dbcon, "SELECT * FROM Meta")
        annoName <- meta$value[which(meta$id == "BioCAnnotation")]
    }
    allIDs <- list()
    mappedSymbols <- list()
    for (i in 1:length(genes)) {
        queryIDs <- mapToProbes(genes[i], annoName = annoName)
        if (!is.na(queryIDs)) {
            allIDs[[i]] <- as.character(queryIDs)
            mappedSymbols[[i]] <- rep(as.character(genes[i]), 
                length(queryIDs))
        }
    }
    names(mappedSymbols) <- genes
    allIDs <- unlist(allIDs)
    mappedSymbols <- unlist(mappedSymbols)
    combMat <- attachCovarsToExpression(dbcon, queryIDs = allIDs)
    combMat <- data.frame(combMat, Symbol = mappedSymbols[match(allIDs, 
        combMat$ProbeID)])
    plotMat <- data.frame(ProbeID = combMat$ProbeID, ArrayID = combMat$ArrayID, 
        Expression = combMat$Expression)
    plotMat <- cast(plotMat, ProbeID ~ ArrayID)
    rownames(plotMat) <- plotMat[, 1]
    plotMat <- plotMat[, -1]
    newRows <- as.character(combMat$Symbol[match(rownames(plotMat), 
        combMat$ProbeID)])
    plotMat <- plotMat[, -1]
    colmat <- matrix(nrow = ncol(plotMat), ncol = 2)
    colmat[, 2] <- "white"
    grouping <- as.factor(combMat$Sample_Group)
    groupcols <- grouping
    levels(groupcols) <- brewer.pal(length(levels(groupcols)), 
        "Paired")
    sampsize <- table(combMat$Sample_Group[match(colnames(plotMat), 
        combMat$ArrayID)])
    hwrite(data.frame(Group = levels(grouping), n = sampsize)[, 
        -2], bgcolor = levels(groupcols), page = outPage, border = 0)
    colmat[, 1] <- as.character(groupcols[match(colnames(plotMat), 
        combMat$ArrayID)])
    hmCol = rev(colorRampPalette(brewer.pal(10, "RdBu"))(64))
    pngfile <- "heatmap.png"
    png(paste(outfile_path, pngfile, sep = "/"), width = 1200, 
        height = 800)
    heatmap.plus(as.matrix(plotMat), ColSideColors = colmat, 
        labRow = newRows, col = hmCol)
    dev.off()
    hwriteImage(pngfile, outPage, br = TRUE)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
