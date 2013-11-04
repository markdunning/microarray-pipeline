\name{getProbesInRegion}
\alias{getProbesInRegion}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find probes within specified genomic region
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getProbesInRegion(chr, start, end, annoName = "illuminaHumanv3")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chr}{
%%     ~~Describe \code{chr} here~~
}
  \item{start}{
%%     ~~Describe \code{start} here~~
}
  \item{end}{
%%     ~~Describe \code{end} here~~
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
function (chr, start, end, annoName = "illuminaHumanv3") 
{
    annoPkg <- paste(annoName, ".db", sep = "")
    annoLoaded <- require(annoPkg, character.only = TRUE)
    if (annoLoaded) {
        mapEnv <- as.name(paste(annoName, "GENOMICLOCATION", 
            sep = ""))
        t <- try(eval(mapEnv), silent = TRUE)
        if (class(t) == "try-error") {
            message(paste("Could not find a GENOMICLOCATION mapping in annotation package ", 
                annoPkg, ". Perhaps it needs updating?", sep = ""))
        }
        else {
            allLocs <- as.list(eval(mapEnv))
            cat("Get probes in Region", chr, ":", start, ":", 
                end, "\n")
            chrs <- unlist(lapply(allLocs, function(x) strsplit(as.character(x), 
                ":")[[1]][1]))
            spos <- unlist(lapply(allLocs, function(x) strsplit(as.character(x), 
                ":")[[1]][2]))
            epos <- unlist(lapply(allLocs, function(x) strsplit(as.character(x), 
                ":")[[1]][3]))
            chrRangs <- RangedData(seqnames = chr, ranges = IRanges(start = as.numeric(spos[which(chrs == 
                chr)]), end = as.numeric(epos[which(chrs == chr)])), 
                names = names(chrs)[which(chrs == chr)])
            query <- RangedData(IRanges(start = as.numeric(start), 
                end = as.numeric(end)))
            matches <- as.matrix(findOverlaps(query, chrRangs))[, 
                2]
            pIDs <- chrRangs$names[matches]
        }
    }
    pIDs
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
