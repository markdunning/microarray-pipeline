\name{combineSummaryDataFromDisk}
\alias{combineSummaryDataFromDisk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Combine summary data from disk
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
combineSummaryDataFromDisk(fileNames, config)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fileNames}{
%%     ~~Describe \code{fileNames} here~~
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
function (fileNames, config) 
{
    annoName <- config["platform", ]
    annoLoaded <- require(paste("illumina", annoName, ".db", 
        sep = ""), character.only = TRUE)
    if (all(sapply(fileNames, function(x) file.exists(x)))) {
        if (annoLoaded) {
            mapEnv <- as.name(paste("illumina", annoName, "ARRAYADDRESS", 
                sep = ""))
            IlluminaIDs = as.character(unlist(mappedkeys(x = eval(mapEnv))))
            eMat <- seMat <- nObsMat <- matrix(nrow = length(IlluminaIDs), 
                ncol = length(fileNames))
            rownames(eMat) <- rownames(seMat) <- rownames(nObsMat) <- IlluminaIDs
        }
        load(fileNames[1])
        pData <- phenoData(bsd)@data
        qcData <- bsd@QC@data
        eMat[, 1] <- exprs(bsd)[match(IlluminaIDs, featureNames(bsd)), 
            ]
        seMat[, 1] <- se.exprs(bsd)[match(IlluminaIDs, featureNames(bsd)), 
            ]
        nObsMat[, 1] <- nObservations(bsd)[match(IlluminaIDs, 
            featureNames(bsd)), ]
        sampleNames <- sampleNames(bsd)
        for (i in 2:length(fileNames)) {
            load(fileNames[i])
            sampleNames = c(sampleNames, sampleNames(bsd))
            eMat[, i] <- exprs(bsd)[match(IlluminaIDs, featureNames(bsd)), 
                ]
            seMat[, i] <- se.exprs(bsd)[match(IlluminaIDs, featureNames(bsd)), 
                ]
            nObsMat[, i] <- nObservations(bsd)[match(IlluminaIDs, 
                featureNames(bsd)), ]
            pData <- rbind(pData, pData(bsd))
            qcData <- merge(qcData, qcData(bsd), all = TRUE, 
                sort = FALSE)
        }
        colnames(eMat) <- colnames(seMat) <- colnames(nObsMat) <- sampleNames
        BSData <- new("ExpressionSetIllumina")
        assayData(BSData) = assayDataNew(exprs = eMat, se.exprs = seMat, 
            nObservations = nObsMat, storage.mode = "list")
        p = new("AnnotatedDataFrame", data.frame(pData, row.names = sampleNames))
        phenoData(BSData) = p
        QC = new("AnnotatedDataFrame", data.frame(qcData, row.names = sampleNames))
        annotation(BSData) <- annoName
        BSData@QC = QC
        BSData@channelData <- list()
        BSData@channelData[[1]] <- rep("G", length(sampleNames))
        status = rep("Unknown", length(IlluminaIDs))
        mapEnv <- as.name(paste("illumina", annoName, "REPORTERGROUPNAME", 
            sep = ""))
        t <- try(eval(mapEnv), silent = TRUE)
        if (class(t) == "try-error") {
            message(paste("Could not find a REPORTERGROUPNAME mapping in annotation package ", 
                annoPkg, ". Perhaps it needs updating?", sep = ""))
        }
        else {
            status[which(!is.na(IlluminaIDs))] = unlist(mget(IlluminaIDs[which(!is.na(IlluminaIDs))], 
                eval(mapEnv), ifnotfound = NA))
            status[which(is.na(status))] = "regular"
        }
        featureData(BSData) = new("AnnotatedDataFrame", data = data.frame(IlluminaID = IlluminaIDs, 
            Status = status, row.names = IlluminaIDs))
        BSData@annotation = annoName
        save(BSData, file = "Robjects/Combined_BSData.Rda")
    }
    else (message("Fatal error; Not all summary files could be found\n"))
    BSData
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
