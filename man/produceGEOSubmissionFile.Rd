\name{produceGEOSubmissionFile}
\alias{produceGEOSubmissionFile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Produce GEO submission file
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
produceGEOSubmissionFile(lumiNormalized, lumiRaw, lib.mapping = NULL, idType = "Probe", sampleInfo = NULL, fileName = "GEOSubmissionFile.txt", supplementaryRdata = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lumiNormalized}{
%%     ~~Describe \code{lumiNormalized} here~~
}
  \item{lumiRaw}{
%%     ~~Describe \code{lumiRaw} here~~
}
  \item{lib.mapping}{
%%     ~~Describe \code{lib.mapping} here~~
}
  \item{idType}{
%%     ~~Describe \code{idType} here~~
}
  \item{sampleInfo}{
%%     ~~Describe \code{sampleInfo} here~~
}
  \item{fileName}{
%%     ~~Describe \code{fileName} here~~
}
  \item{supplementaryRdata}{
%%     ~~Describe \code{supplementaryRdata} here~~
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
function (lumiNormalized, lumiRaw, lib.mapping = NULL, idType = "Probe", 
    sampleInfo = NULL, fileName = "GEOSubmissionFile.txt", supplementaryRdata = FALSE, 
    ...) 
{
    if (missing(lumiNormalized)) 
        stop("Please provide all required input parameters!\n")
    if (is.matrix(lumiNormalized)) {
        expr.norm <- lumiNormalized
        detect <- NULL
        se.expr <- NULL
        expr <- NULL
        beadNum <- NULL
        if (!missing(lumiRaw)) 
            expr <- lumiRaw
    }
    else {
        expr.norm <- exprs(lumiNormalized)
        detect <- Detection(lumiRaw)
        se.expr <- se.exprs(lumiRaw)
        expr <- exprs(lumiRaw)
        if (!is.null(nObservations(lumiRaw))) 
            beadNum <- nObservations(lumiRaw)
        else beadNum <- BSData@assayData$NoBeads
    }
    if (is.null(sampleInfo)) {
        sampleInfo <- produceGEOSampleInfoTemplate(lumiNormalized, 
            lib.mapping = lib.mapping, fileName = NULL)
    }
    else if (length(sampleInfo) == 1 && is.character(sampleInfo)) {
        sampleInfo <- read.table(sampleInfo, sep = "\t", colClasses = "character", 
            skip = 1, head = TRUE, strip.white = TRUE, quote = "")
    }
    else if (is.null(nrow(sampleInfo))) {
        stop("Please provide correct sample information (a data.frame, matrix, or sampleInfo file)!\n")
    }
    sampleInfoTitle <- colnames(sampleInfo)
    if (supplementaryRdata) 
        sampleInfo[, "Sample_supplementary_file"] <- "supplementaryData.Rdata"
    if (is.matrix(lumiNormalized)) {
        nuID <- rownames(lumiNormalized)
    }
    else {
        nuID <- featureNames(lumiNormalized)
    }
    chipVersion <- annotation(lumiNormalized)
    annoPkg <- paste("illumina", chipVersion, ".db", sep = "")
    probeId <- nuID
    annoLoaded <- require(annoPkg, character.only = TRUE)
    if (annoLoaded) {
        mapEnv <- as.name(paste("illumina", chipVersion, "ARRAYADDRESS", 
            sep = ""))
        probeId = as.character(unlist(mget(as.character(probeId), 
            revmap(eval(mapEnv)), ifnotfound = NA)))
    }
    sampleID <- sampleInfo[, "sampleID"]
    sampleTitle <- sampleInfo[, "Sample_title"]
    for (i in 1:length(sampleID)) {
        if (i == 1) {
            cat("^SAMPLE =", sampleTitle[i], "\n", sep = "", 
                file = fileName, append = FALSE)
        }
        else {
            cat("^SAMPLE =", sampleTitle[i], "\n", sep = "", 
                file = fileName, append = TRUE)
        }
        sampleInfo.i <- paste("!", sampleInfoTitle[-1], " = ", 
            sampleInfo[i, -1], "\n", sep = "", collapse = "")
        sampleInfo.i <- gsub("'", "\\'", sampleInfo.i)
        cat(sampleInfo.i, file = fileName, append = TRUE, sep = "")
        tableHead <- "ID_REF"
        cat("#ID_REF = \n", file = fileName, append = TRUE)
        if (!is.null(nuID)) {
            cat(paste("#ArrayAddressID = Illumina bead-level identifier. See Bioconductor package,", 
                annoPkg, "for more details.\n"), file = fileName, 
                append = TRUE)
            tableHead <- c(tableHead, "ArrayAddress")
        }
        cat("#VALUE = normalized signal intensity\n", file = fileName, 
            append = TRUE)
        tableHead <- c(tableHead, "VALUE")
        if (!is.null(expr)) {
            cat("#RAW_VALUE = raw signal intensity\n", file = fileName, 
                append = TRUE)
            tableHead <- c(tableHead, "RAW_VALUE")
        }
        if (!is.null(se.expr)) {
            cat("#BEAD_STDERR = the standard error of the probe measurements\n", 
                file = fileName, append = TRUE)
            tableHead <- c(tableHead, "BEAD_STDERR")
        }
        if (!is.null(detect)) {
            cat("#Detection_Pval = the detection p-value of the probe\n", 
                file = fileName, append = TRUE)
            tableHead <- c(tableHead, "Detection_Pval")
        }
        if (!is.null(beadNum)) {
            cat("#Avg_NBEADS = Number of beads for the probe\n", 
                file = fileName, append = TRUE)
            tableHead <- c(tableHead, "Avg_NBEADS")
        }
        sampleTable.i <- probeId
        if (!is.null(nuID)) 
            sampleTable.i <- cbind(sampleTable.i, nuID)
        sampleTable.i <- cbind(sampleTable.i, expr.norm[, i])
        if (!is.null(expr)) 
            sampleTable.i <- cbind(sampleTable.i, expr[, i])
        if (!is.null(se.expr)) 
            sampleTable.i <- cbind(sampleTable.i, se.expr[, i])
        if (!is.null(detect)) 
            sampleTable.i <- cbind(sampleTable.i, detect[, i])
        if (!is.null(beadNum)) 
            sampleTable.i <- cbind(sampleTable.i, beadNum[, i])
        sampleTable.i <- rbind(tableHead, sampleTable.i)
        cat("!sample_table_begin\n", file = fileName, append = TRUE)
        write.table(sampleTable.i, sep = "\t", quote = FALSE, 
            file = fileName, append = TRUE, col.names = FALSE, 
            row.names = FALSE)
        cat("!sample_table_end\n", file = fileName, append = TRUE)
    }
    if (supplementaryRdata) {
        lumiNormalized <- lumiNormalized[, sampleID]
        if (!missing(lumiRaw)) {
            lumiRaw <- lumiRaw[, sampleID]
            save(lumiNormalized, lumiRaw, sampleInfo, file = "supplementaryData.Rdata")
        }
        else {
            save(lumiNormalized, sampleInfo, file = "supplementaryData.Rdata")
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
