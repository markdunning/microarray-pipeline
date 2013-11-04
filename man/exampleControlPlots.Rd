\name{exampleControlPlots}
\alias{exampleControlPlots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make example control plots}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
exampleControlPlots(BSData, config)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BSData}{
%%     ~~Describe \code{BSData} here~~
}
  \item{config}{
%%     ~~Describe \code{config} here~~
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
function (BSData, config) 
{
    snr <- qcData(BSData)[, "P95Grn"]/qcData(BSData)[, "P05Grn"]
    worstArray <- sampleNames(BSData)[which.min(snr)]
    bld <- readIllumina(sectionNames = worstArray, illuminaAnnotation = annotation(BSData))
    p1 <- combinedControlPlot(bld) + opts(title = paste("Control plot for array with signal-to-noise ratio of", 
        round(min(snr), 3)))
    bestArray <- sampleNames(BSData)[which.max(snr)]
    bld <- readIllumina(sectionNames = bestArray, illuminaAnnotation = annotation(BSData))
    p2 <- combinedControlPlot(bld) + opts(title = paste("Control plot for array with signal-to-noise ratio of", 
        round(max(snr), 3)))
    png("QA/controlPlotsComparison.png", width = 1200, height = 1200)
    gridExtra::grid.arrange(p1, p2, ncol = 1)
    dev.off()
    if (config["bash", ] == 1) {
        BASHMetrics <- qcData(BSData)[, "BeadsMasked"]
        worstArray <- sampleNames(BSData)[which.max(BASHMetrics)]
        load(paste("Robjects/", worstArray, "-BLData.Rda", sep = ""))
        pctMasked <- round(max(BASHMetrics)/numBeads(bld), 3) * 
            100
        p1 <- beadarray:::imageplot(BLData.bashed) + opts(title = paste("Array most affected by spatial artefacts, with", 
            pctMasked, "\% of beads masked by BASH"))
        p2 <- beadarray:::showArrayMask2(BLData.bashed, array = 1)
        bestArray <- sampleNames(BSData)[which.min(BASHMetrics)]
        load(paste("Robjects/", bestArray, "-BLData.Rda", sep = ""))
        pctMasked <- round(min(BASHMetrics)/numBeads(bld), 3) * 
            100
        p3 <- beadarray:::imageplot(BLData.bashed) + opts(title = paste("Array least affected by spatial artefacts, with", 
            pctMasked, "\% of beads masked by BASH"))
        p4 <- beadarray:::showArrayMask2(BLData.bashed, array = 1)
        png("QA/spatialArtefactComparison.png", width = 1200, 
            height = 800)
        gridExtra::grid.arrange(p1, p3, p2, p4, ncol = 2)
        dev.off()
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
