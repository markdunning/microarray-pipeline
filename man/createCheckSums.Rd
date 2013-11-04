\name{createCheckSums}
\alias{createCheckSums}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
make md5 sums
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
createCheckSums(sampleSheet)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sampleSheet}{
%%     ~~Describe \code{sampleSheet} here~~
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
function (sampleSheet) 
{
    library(tools)
    checksum <- NULL
    checked <- NULL
    for (i in 1:nrow(sampleSheet)) {
        chip <- sampleSheet$Sentrix_ID[i]
        section <- sampleSheet$Sentrix_Position[i]
        dirFiles <- list.files(as.character(chip))
        toCheck <- paste(chip, dirFiles[grep(paste(chip, section, 
            sep = "_"), dirFiles)], sep = "/")
        checksum <- c(checksum, sapply(toCheck, md5sum))
        checked <- c(checked, toCheck)
    }
    data.frame(checked, checksum)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
