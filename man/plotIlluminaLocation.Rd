\name{plotIlluminaLocation}
\alias{plotIlluminaLocation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot locations of Illumina probes in relation to genes
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotIlluminaLocation(Symbol, ensembl = NULL, annoName = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Symbol}{
%%     ~~Describe \code{Symbol} here~~
}
  \item{ensembl}{
%%     ~~Describe \code{ensembl} here~~
}
  \item{annoName}{
%%     ~~Describe \code{annoName} here~~
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
function (Symbol, ensembl = NULL, annoName = NULL) 
{
    require("GenomeGraphs")
    if (is.null(ensembl)) {
        ensembl = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
    }
    annoPkg <- paste(annoName, ".db", sep = "")
    annoLoaded <- require(annoPkg, character.only = TRUE)
    if (annoLoaded) {
        mapEnv <- as.name(paste(annoName, "SYMBOL", sep = ""))
        t <- try(eval(mapEnv), silent = TRUE)
        if (class(t) == "try-error") {
            message(paste("Could not find a ", from, "mapping in annotation package ", 
                annoPkg, ". Perhaps it needs updating?", sep = ""))
        }
        else {
            pIDs <- unlist(mget(Symbol, revmap(eval(mapEnv)), 
                ifnotfound = NA))
        }
    }
    else message("Could not load package ", annoPkg)
    gn <- try(makeGene(Symbol, type = "hgnc_symbol", biomart = ensembl))
    chrName <- paste("chr", getBM(attributes = "chromosome_name", 
        filters = "hgnc_symbol", value = Symbol, ensembl)[[1]], 
        sep = "")
    uExon <- unique(gn@ens[, 3])
    exStart <- unlist(lapply(split(gn@ens[, 5], gn@ens[, 3]), 
        unique))
    exEnd <- unlist(lapply(split(gn@ens[, 4], gn@ens[, 3]), unique))
    mapEnv <- as.name(paste(annoName, "GENOMICLOCATION", sep = ""))
    probeLocs <- mget(pIDs, eval(mapEnv))
    pStarts <- as.numeric(unlist(lapply(probeLocs, function(x) strsplit(x, 
        ":")[[1]][2])))
    pEnds <- as.numeric(unlist(lapply(probeLocs, function(x) strsplit(x, 
        ":")[[1]][3])))
    pStrands <- unlist(lapply(probeLocs, function(x) strsplit(x, 
        ":")[[1]][4]))
    pOverlays = list()
    for (i in 1:length(pIDs)) {
        if (pStrands[i] == "+") {
            pOverlays = append(pOverlays, makeRectangleOverlay(start = as.numeric(pStarts[i]), 
                end = as.numeric(pEnds[i]), region = c(1, 2)))
            pOverlays = append(pOverlays, makeTextOverlay(pIDs[i], 
                xpos = as.numeric(pStarts[i]), ypos = 0.95))
        }
        else {
            pOverlays = append(pOverlays, makeRectangleOverlay(start = as.numeric(pStarts[i]), 
                end = as.numeric(pEnds[i]), region = c(1, 2)))
            pOverlays = append(pOverlays, makeTextOverlay(pIDs[i], 
                xpos = as.numeric(pStarts[i]), ypos = 0.95))
        }
    }
    minX = min(min(exStart), min(pStarts)) - 1000
    maxX = max(max(exEnd), max(pEnds)) + 1000
    gdPlot(list(makeGenomeAxis(add35 = TRUE, add53 = TRUE), gn), 
        overlays = pOverlays, minBase = minX, maxBase = maxX)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
