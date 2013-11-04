\name{validUser}
\alias{validUser}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Check if a user is allowed to access database
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
validUser(dbcon = NULL, email, table = "users", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dbcon}{
%%     ~~Describe \code{dbcon} here~~
}
  \item{email}{
%%     ~~Describe \code{email} here~~
}
  \item{table}{
%%     ~~Describe \code{table} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (dbcon = NULL, email, table = "users", ...) 
{
    if (is.null(dbcon)) {
        dbcon <- connectToAnalysisDb(, ...)
    }
    userList <- dbGetQuery(dbcon, paste("select * from", table))
    return(tolower(email) \%in\% tolower(userList[, 2]))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
