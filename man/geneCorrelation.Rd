\name{geneCorrelation}
\alias{geneCorrelation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot correlation of two genes
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
geneCorrelation(dbcon, gene1, gene2, factor = "Sample_Group", selectBestProbe = FALSE, annoName = NULL, outfile = "Results.html", outfile_path = ".")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dbcon}{
%%     ~~Describe \code{dbcon} here~~
}
  \item{gene1}{
%%     ~~Describe \code{gene1} here~~
}
  \item{gene2}{
%%     ~~Describe \code{gene2} here~~
}
  \item{factor}{
%%     ~~Describe \code{factor} here~~
}
  \item{selectBestProbe}{
%%     ~~Describe \code{selectBestProbe} here~~
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
function (dbcon, gene1, gene2, factor = "Sample_Group", selectBestProbe = FALSE, 
    annoName = NULL, outfile = "Results.html", outfile_path = ".") 
{
    library(hwriter)
    dir.create(outfile_path)
    outPage = openPage(filename = outfile)
    if (is.null(annoName)) {
        meta <- dbGetQuery(dbcon, "SELECT * FROM Meta")
        annoName <- meta$value[which(meta$id == "annotation")]
    }
    gene1IDs <- mapToProbes(gene1, annoName = annoName)
    plotCount <- 1
    if (!is.na(gene1IDs)) {
        locsPic <- paste(outfile_path, "/", gene1, "-locations.png", 
            sep = "")
        png(locsPic, width = 800, height = 400)
        plotIlluminaLocation(gene1, annoName = annoName)
        dev.off()
        heading = paste("Positions of probes for ", gene1, " on the Illumina ", 
            annoName, " chip", sep = "")
        hwrite(heading, outPage, heading = 2, br = TRUE)
        hwriteImage(locsPic, outPage)
        fields <- c("PROBEQUALITY", "CODINGZONE", "GENOMICLOCATION", 
            "REPEATMASK", "SECONDMATCHES", "ENTREZID")
        probeSummary <- matrix(nrow = length(gene1IDs), ncol = length(fields))
        colnames(probeSummary) <- fields
        rownames(probeSummary) <- gene1IDs
        for (k in 1:length(gene1IDs)) {
            probeSummary[k, ] <- unlist(sapply(fields, function(x) mget(gene1IDs[j], 
                eval(as.name(paste("illumina", annoName, x, sep = ""))))))
        }
        heading = paste("Summary of annotation of probes for ", 
            gene1, " on the Illumina ", annoName, " chip", sep = "")
        hwrite(heading, outPage)
        hwrite(probeSummary, outPage)
        hwrite("Correlation plots", outPage)
        hwrite(gene1IDs, link = paste("#", gene1IDs, sep = ""), 
            page = outPage)
        for (j in 1:length(gene2)) {
            queryIDs <- mapToProbes(gene2[j], annoName = annoName)
            if (!is.na(queryIDs)) {
                heading = paste("Correltion of ", gene1, " with ", 
                  gene2[j], sep = "")
                hwrite(heading, outPage, heading = 2, br = TRUE)
                locsPic <- paste(outfile_path, "/", gene2[j], 
                  "-locations.png", sep = "")
                png(locsPic, width = 800, height = 400)
                plotIlluminaLocation(gene2[j], annoName = annoName)
                dev.off()
                heading = paste("Positions of probes for ", gene2[j], 
                  " on the Illumina ", annoName, " chip", sep = "")
                hwrite(heading, outPage, heading = 2, br = TRUE)
                hwriteImage(locsPic, outPage, br = TRUE)
                fields <- c("PROBEQUALITY", "CODINGZONE", "GENOMICLOCATION", 
                  "REPEATMASK", "SECONDMATCHES", "ENTREZID")
                probeSummary <- matrix(nrow = length(queryIDs), 
                  ncol = length(fields))
                colnames(probeSummary) <- fields
                rownames(probeSummary) <- queryIDs
                for (k in 1:length(queryIDs)) {
                  probeSummary[k, ] <- unlist(sapply(fields, 
                    function(x) mget(queryIDs[j], eval(as.name(paste("illumina", 
                      annoName, x, sep = ""))))))
                }
                heading = paste("Summary of annotation of probes for ", 
                  gene2[j], " on the Illumina ", annoName, " chip", 
                  sep = "")
                hwrite(heading, outPage, br = TRUE)
                hwrite(probeSummary, outPage)
                allIDs <- c(gene1IDs, queryIDs)
                symbols <- c(rep(gene1, length(gene1IDs)), rep(gene2[j], 
                  length(queryIDs)))
                names(symbols) <- allIDs
                eMat <- getExpressionForIDs(dbcon, queryIDs = allIDs)
                corMat <- cor(t(cast(eMat, ProbeID ~ ArrayID)))
                corMat <- corMat[allIDs, allIDs]
                colMat <- matrix(nrow = length(gene1IDs) + length(queryIDs), 
                  ncol = length(gene1IDs) + length(queryIDs), 
                  "#ccffdd")
                colMat[1:length(gene1IDs), 1:length(gene1IDs)] <- "#ffffaa"
                colMat[(length(gene1IDs) + 1):nrow(colMat), (length(gene1IDs) + 
                  1):nrow(colMat)] <- "#ffffaa"
                hwrite("Summary of correlation", page = outPage, 
                  heading = 2)
                hwrite(corMat, outPage, br = TRUE, bgcolor = colMat)
                hwrite("Intensity Summary", page = outPage, heading = 2)
                pngfile <- paste(gene2[j], "-comparison.png", 
                  sep = "")
                eMat <- data.frame(eMat, Symbol = symbols[eMat$ProbeID])
                plot <- ggplot(eMat, aes(x = Symbol, y = Expression, 
                  fill = ProbeID)) + geom_boxplot() + coord_flip()
                ggsave(plot, filename = paste(outfile_path, pngfile, 
                  sep = "/"), width = 6, height = 6, dpi = 100)
                hwriteImage(pngfile, outPage, br = TRUE)
                for (id1 in gene1IDs) {
                  heading = paste("Correlations of ", id1, " with  ", 
                    gene2[j], sep = "")
                  hwrite(heading, outPage, heading = 2, br = TRUE, 
                    name = id1)
                  for (id2 in queryIDs) {
                    combMat <- attachCovarsToExpression(dbcon, 
                      queryIDs = id1)
                    exp2 <- getExpressionForIDs(dbcon, queryIDs = id2)
                    factorcol <- which(colnames(combMat) == factor)
                    combMat <- data.frame(combMat, Gene2 = exp2$Expression)
                    if (any(is.na(combMat[, factorcol]))) {
                      combMat2 <- combMat[-which(is.na(combMat[, 
                        factorcol])), ]
                    }
                    else combMat2 <- combMat
                    pngfile1 <- paste(plotCount, "-all.png", 
                      sep = "")
                    colnames(combMat2)[factorcol] <- "FactorToPlot"
                    myplot <- ggplot(combMat2, aes(x = Expression, 
                      y = Gene2, col = FactorToPlot)) + geom_point() + 
                      xlab(paste(gene1, id1, sep = ":")) + ylab(paste(gene2[j], 
                      id2, sep = ":"))
                    ggsave(myplot, filename = paste(outfile_path, 
                      pngfile1, sep = "/"), width = 6, height = 3, 
                      dpi = 100)
                    hwriteImage(pngfile1, outPage)
                    pngfile2 <- paste(plotCount, "-groups.png", 
                      sep = "")
                    myplot2 <- ggplot(combMat2, aes(x = Expression, 
                      y = Gene2)) + geom_point() + geom_smooth() + 
                      facet_wrap(~FactorToPlot) + xlab(paste(gene1, 
                      id1, sep = ":")) + ylab(paste(gene2[j], 
                      id2, sep = ":"))
                    ggsave(myplot2, filename = paste(outfile_path, 
                      pngfile2, sep = "/"), width = 6, height = 3, 
                      dpi = 100)
                    hwriteImage(pngfile2, outPage, br = TRUE)
                    plotCount <- plotCount + 1
                  }
                }
            }
        }
        closePage(outPage)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
