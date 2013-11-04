\name{generateColourMatrix}
\alias{generateColourMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate a colour matrix
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
generateColourMatrix(targets)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{targets}{
%%     ~~Describe \code{targets} here~~
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
function (targets) 
{
    col_matrix <- matrix(nrow = nrow(targets), ncol = ncol(targets), 
        dimnames = list(row.names(targets), colnames(targets)))
    reject <- vector()
    for (c in 1:length(colnames(targets))) {
        types <- levels(as.factor(targets[, c]))
        if (colnames(col_matrix)[c] == "Sentrix_ID" || length(types) > 
            1 && length(types) != nrow(targets) && colnames(col_matrix)[c] != 
            "RIN" && colnames(col_matrix)[c] != "registrationScore" && 
            colnames(col_matrix)[c] != "focusScore" && colnames(col_matrix)[c] != 
            "Controls") {
            if (length(types) == 2) {
                cols <- brewer.pal(3, config["colScheme", ])[c(1, 
                  2)]
            }
            else {
                cols <- colorRampPalette(brewer.pal(8, config["colScheme", 
                  ]))(length(types))
            }
            col_matrix[, c] <- cols[match(targets[, c], types)]
        }
        else {
            reject[length(reject) + 1] <- c
        }
    }
    if (length(reject) != 0) {
        col_matrix <- matrix(col_matrix[, -reject], nrow = nrow(col_matrix))
        row.names(col_matrix) <- row.names(targets)
        colnames(col_matrix) <- colnames(targets)[-reject]
    }
    col_matrix
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
