\name{boxplotFromPublicData}
\alias{boxplotFromPublicData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
make a boxplot from a public dataset
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
boxplotFromPublicData(dbcon, gene, factor = "Sample_Group", annoName = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dbcon}{
%%     ~~Describe \code{dbcon} here~~
}
  \item{gene}{
%%     ~~Describe \code{gene} here~~
}
  \item{factor}{
%%     ~~Describe \code{factor} here~~
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
function (dbcon, gene, factor = "Sample_Group", annoName = NULL) 
{
    if ("fData" \%in\% dbListTables(dbcon)) {
        queryIDs <- dbGetQuery(dbcon, paste("Select probe FROM fData where Symbol = '", 
            gene, "'", sep = ""))
    }
    else {
        if (is.null(annoName)) {
            meta <- dbGetQuery(dbcon, "SELECT * FROM Meta")
            annoName <- meta$value[which(meta$id == "BiocAnnotation")]
        }
        queryIDs <- mapToProbes(gene, annoName = annoName)
    }
    if (!is.na(queryIDs)) {
        combMat <- attachCovarsToExpression(dbcon, queryIDs = queryIDs)
        factorcol <- which(colnames(combMat) == factor)
        combMat <- attachCovarsToExpression(dbcon, queryIDs = queryIDs)
        if (nrow(combMat) > 0) {
            if (any(is.na(combMat[, factorcol]))) {
                combMat2 <- combMat[-which(is.na(combMat[, factorcol])), 
                  ]
            }
            else combMat2 <- combMat
            colnames(combMat2)[factorcol] <- "FactorToPlot"
            myplot <- ggplot(combMat2, aes(x = FactorToPlot, 
                y = Expression, fill = FactorToPlot)) + geom_boxplot() + 
                xlab(factor) + scale_fill_discrete(name = factor) + 
                facet_wrap(~ProbeID, ncol = 4)
            myplot
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
