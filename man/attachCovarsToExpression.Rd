\name{attachCovarsToExpression}
\alias{attachCovarsToExpression}
\title{
Combine expression and meta-data for probes
}
\description{
Queries both the assay data and meta data for a set of IDs of interest, and combine into a data frame that is ameanable to plotting with ggplot2
}
\usage{
attachCovarsToExpression(DB, assayDataTable = "assayData", covarTable = "Covars", sampleCol = "ArrayID", idCol = "ProbeID", queryIDs = NULL, covars = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DB}{
A database connection
}
  \item{assayDataTable}{
Name of table within the database that contains assaydata information

}
  \item{covarTable}{
Name of table within the database that contains covariate (meta) data
}
  \item{sampleCol}{
Column name that indicates the array identifier. Must be consistent between assay and covariate tables
}
  \item{idCol}{
Column name that indicates the probe identifier. Must be consistent between assay and covariate tables
}
  \item{queryIDs}{
Probe IDs of interest
}
  \item{covars}{
names of covariates to return, the default being all
}
}
\details{

}
\value{

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Mark Dunning
}


\examples{

maqc <- connectToLocalDb(paste(system.file("data", package="crimp"), "MAQC.sqlite",sep="/"))

ids <- as.character(mapToProbes("PTEN", annoName="illuminaHumanv3"))
df <- attachCovarsToExpression(maqc, queryIDs = ids)
head(df)


}

