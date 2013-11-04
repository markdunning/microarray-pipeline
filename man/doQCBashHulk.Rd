\name{doQCBashHulk}
\alias{doQCBashHulk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Do BASH and HULK
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
doQCBashHulk(x, config, sampleSheet)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{config}{
%%     ~~Describe \code{config} here~~
}
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
function (x, config, sampleSheet) 
{
    library(beadarray)
    mouse = FALSE
    sID <- strsplit(x, "_")[[1]][1]
    secNames <- x
    if (mouse) {
        secNames <- paste(x, c(1, 2), sep = "_")
    }
    message("Reading the specified section...")
    if (config["baseDir", ] == "WORKINGDIR") {
        BLData <- readIllumina(useImages = as.logical(config["useImages", 
            ]), illuminaAnnotation = config["platform", ], sectionNames = secNames, 
            sampleSheet = sampleSheet)
    }
    else {
        if (file.exists(config["baseDir", ])) {
            BLData <- readIllumina(dir = config["baseDir", ], 
                useImages = as.logical(config["useImages", ]), 
                illuminaAnnotation = config["platform", ], sectionNames = secNames, 
                sampleSheet = sampleSheet)
        }
        else {
            message("base directory defined in config file does not exist")
            BLData <- readIllumina(useImages = as.logical(config["useImages", 
                ]), illuminaAnnotation = config["platform", ], 
                sectionNames = secNames, sampleSheet = sampleSheet)
        }
    }
    if (config["doqc", ] == 1) {
        message("Creating QC plots...")
        tmp <- try(expressionQCPipeline(BLData, qcDir = config["qc_path", 
            ], plotType = ".png", overWrite = F, zlim = NULL, 
            tagsToDetect = NULL))
    }
    probableAnnotation <- suggestAnnotation(BLData)
    if (probableAnnotation != annotation(BLData)) 
        message(paste("This data was read in as ", annotation(BLData), 
            " but beadarray thinks it is ", probableAnnotation, 
            ". Please check your config file\n", sep = ""))
    if (config["beadRegistration", ]) {
        message("Making bead registration checks")
        dir.create(paste(config["qc_path", ], "/registrationChecks/", 
            sep = ""), recursive = TRUE, showWarnings = FALSE)
        regCheck <- try(checkRegistration(BLData))
        if (class(regCheck) != "try-error") {
            png(paste(config["qc_path", ], "/registrationChecks/", 
                x, ".png", sep = ""), width = as.numeric(config["png_w", 
                ]), height = as.numeric(config["png_h", ]), type = "cairo")
            boxplot(regCheck, main = x, plotP95 = TRUE)
            dev.off()
        }
        else message("Bead registration check was not successful, proceeding anyway..")
    }
    if (config["bash", ] == 1) {
        weights <- BASH(BLData, array = 1)
        BLData.bashed <- setWeights(BLData, wts = weights$wts, 
            array = 1)
        QCmat <- weights$QC
        BLData.bashed = insertSectionData(BLData.bashed, what = "BASHQC", 
            data = QCmat)
    }
    else {
        BLData.bashed <- BLData
    }
    if (as.numeric(config["bash", ])) {
        dir.create("QA", showWarnings = FALSE)
        dir.create("beadarray_QA/imagePlots", showWarnings = FALSE, 
            recursive = T)
        ipHt <- 1000
        png(file = paste("beadarray_QA/imagePlots/", x, ".png", 
            sep = ""), width = 1000, height = ipHt)
        plots <- 2 + as.numeric((config["hulk", ]))
        p1 <- beadarray:::showArrayMask2(BLData.bashed, array = 1)
        logData <- logGreenChannelTransform(BLData.bashed, array = 1)
        BLData.logged <- insertBeadData(BLData.bashed, array = 1, 
            data = logData, what = "Grn")
        p2 <- imageplot(BLData.logged, array = 1, low = "green", 
            high = "red", main = paste(x, "Before", sep = "-"))
        if (config["hulk", ] == 1) {
            hulk <- HULK(BLData.bashed, array = 1)
            BLData.bashed <- insertBeadData(BLData.bashed, array = 1, 
                data = hulk, what = "Grn")
            p3 <- imageplot(BLData.bashed, array = 1, low = "green", 
                high = "red", main = paste(x, "After", sep = "-"))
            gridExtra:::grid.arrange(p1, p2, p3, ncol = 1)
        }
        else gridExtra:::grid.arrange(p1, p2, ncol = 1)
        dev.off()
    }
    if (config["makeBabFile", ]) {
        txtFile <- as.character(BLData@sectionData$Targets$textFile)
        locsGreen <- as.character(BLData@sectionData$Targets$greenLocs)
        path <- as.character(BLData@sectionData$Targets$directory)
        if (length(grep(".bab", txtFile)) == 0) {
            message("Creating .bab compressed data")
            txtFile <- as.character(BLData@sectionData$Targets$textFile)
            locsGreen <- as.character(BLData@sectionData$Targets$greenLocs)
            path <- as.character(BLData@sectionData$Targets$directory)
            BeadDataPackR:::compressBeadData(txtFile = txtFile, 
                locsGrn = locsGreen, path = path)
            if (config["removeTextAndLocs", ]) {
                file.remove(c(paste(path, "/", txtFile, sep = ""), 
                  paste(path, "/", locsGreen, sep = "")))
            }
        }
        else message("Compressed bab file already exists...")
    }
    outfile <- paste("Robjects/", x, "-BLData.Rda", sep = "")
    message("saving bead-level data and weights as ", outfile, 
        sep = "")
    save(weights, BLData.bashed, file = outfile)
    myMean = function(x) mean(x, na.rm = TRUE)
    mySd = function(x) sd(x, na.rm = TRUE)
    fooFun = function(BLData, array) {
        x = getBeadData(BLData, array = array, what = "Grn")
    }
    greenChannel = new("illuminaChannel", fooFun, illuminaOutlierMethod, 
        myMean, mySd, "G")
    bsd <- beadarray::summarize(BLData.bashed, list(greenChannel), 
        useSampleFac = F)
    outfile <- paste("Robjects/", x, "-summary.Rda", sep = "")
    save(bsd, file = outfile)
    outfile
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
