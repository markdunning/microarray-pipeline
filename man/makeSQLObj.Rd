\name{makeSQLObj}
\alias{makeSQLObj}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create new database object from an R object or matrix
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
makeSQLObj(Robj, name = "mydata", covarTable = NULL, annotation = "Humanv3", users = c("mark.dunning@cancer.org.uk"), deTable = NULL, fData = NULL, GSE = NULL, StudyTitle = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Robj}{
%%     ~~Describe \code{Robj} here~~
}
  \item{name}{
%%     ~~Describe \code{name} here~~
}
  \item{covarTable}{
%%     ~~Describe \code{covarTable} here~~
}
  \item{annotation}{
%%     ~~Describe \code{annotation} here~~
}
  \item{users}{
%%     ~~Describe \code{users} here~~
}
  \item{deTable}{
%%     ~~Describe \code{deTable} here~~
}
  \item{fData}{
%%     ~~Describe \code{fData} here~~
}
  \item{GSE}{
%%     ~~Describe \code{GSE} here~~
}
  \item{StudyTitle}{
%%     ~~Describe \code{StudyTitle} here~~
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
library(beadarrayExampleData)
data(exampleSummaryData)
maqc <- channel(exampleSummaryData, "G")
makeSQLObj(Robj=maqc,annotation="illuminaHumanv3", GSE="GSEXXX", name="MAQC", StudyTitle="MAQC")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
