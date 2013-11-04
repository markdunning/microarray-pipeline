\name{geneSummary}
\alias{geneSummary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarise expression of set of genes in a dataset
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
geneSummary(dbcon, genes, factor = "Sample_Group", selectBestProbe = FALSE, annoName = NULL, outfile = "Results.html", outfile_path = ".", externalData = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dbcon}{
%%     ~~Describe \code{dbcon} here~~
}
  \item{genes}{
%%     ~~Describe \code{genes} here~~
}
  \item{factor}{
%%     ~~Describe \code{factor} here~~
}
  \item{selectBestProbe}{
%%     ~~Describe \code{selectBestProbe} here~~
}
  \item{annoName}{
%%     ~~Describe \code{annoName} here~~
}
  \item{outfile}{
%%     ~~Describe \code{outfile} here~~
}
  \item{outfile_path}{
%%     ~~Describe \code{outfile_path} here~~
}
  \item{externalData}{
%%     ~~Describe \code{externalData} here~~
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
function (dbcon, genes, factor = "Sample_Group", selectBestProbe = FALSE, 
    annoName = NULL, outfile = "Results.html", outfile_path = ".", 
    externalData = NULL) 
{
    library(hwriter)
    library(reshape)
    dir.create(outfile_path)
    outPage = openPage(filename = outfile)
    if (is.null(annoName)) {
        meta <- dbGetQuery(dbcon, "SELECT * FROM Meta")
        annoName <- meta$value[which(meta$id == "BiocAnnotation")]
    }
    allIDs <- list()
    mappedSymbols <- list()
    for (i in 1:length(genes)) {
        queryIDs <- mapToProbes(genes[i], annoName = annoName)
        if (!is.na(queryIDs)) {
            allIDs[[i]] <- as.character(queryIDs)
            mappedSymbols[[i]] <- rep(as.character(genes[i]), 
                length(queryIDs))
        }
    }
    names(mappedSymbols) <- genes
    linksToPlots <- paste("#", genes, "Boxplot", sep = "")
    linksToAnnotation <- paste("#", genes, "Annotation", sep = "")
    linksToDE <- paste("#", genes, "DE", sep = "")
    NumberOfProbes = unlist(lapply(mappedSymbols, length))
    linksToPlots[which(NumberOfProbes == 0)] <- ""
    linksToAnnotation[which(NumberOfProbes == 0)] <- ""
    linksToDE[which(NumberOfProbes == 0)] <- ""
    summary <- data.frame(Symbol = genes, NumberOfProbes = unlist(lapply(mappedSymbols, 
        length)), Annotation = "Annotation", Boxplots = "Boxplots", 
        DifferentialExpression = "DifferentialExpression")
    hwrite("Report Summary", outPage, heading = 1, br = TRUE, 
        name = "Summary")
    hwrite(summary, outPage, br = TRUE, col.link = list(Annotation = linksToAnnotation, 
        Boxplots = linksToPlots, DifferentialExpression = linksToDE))
    hwrite("Diferential expression summary", heading = 2, br = TRUE, 
        outPage, link = "#DESummary")
    allIDs <- unlist(allIDs)
    mappedSymbols <- unlist(mappedSymbols)
    fullDE <- dbGetQuery(dbcon, "SELECT * FROM DeTable")
    nprobes <- length(unique(fullDE$ProbeID))
    Rank = rep(NA, nrow(fullDE))
    rnk <- lapply(split(as.numeric(fullDE$LogOdds), as.character(fullDE$Contrasts)), 
        rank)
    for (i in 1:length(rnk)) {
        Rank[grep(names(rnk)[i], fullDE$Contrasts)] <- rnk[[i]]
    }
    fullDE <- data.frame(fullDE, Significant = as.numeric(fullDE$Adjusted) < 
        0.05, Rank = nprobes - Rank + 1)
    genes <- genes[which(NumberOfProbes > 0)]
    for (i in 1:length(genes)) {
        hwrite("<----------------------------------------------------------------------------------------------->", 
            outPage, br = TRUE)
        hwrite("", outPage, br = TRUE)
        hwrite("", outPage, br = TRUE)
        hwrite("", outPage, br = TRUE)
        hwrite("", outPage, br = TRUE)
        locsPic <- paste(outfile_path, "/", genes[i], "-locations.png", 
            sep = "")
        png(locsPic, width = 800, height = 400)
        plotIlluminaLocation(genes[i], annoName = annoName)
        dev.off()
        hwrite(genes[i], outPage, heading = 2, br = TRUE, name = genes[i])
        hwrite("Navigation", outPage, br = TRUE)
        hwrite(genes, link = paste("#", genes, sep = ""), page = outPage, 
            br = TRUE, border = 0)
        hwrite("Diferential expression summary", br = TRUE, outPage, 
            link = "#DESummary")
        hwrite("Report Summary", outPage, br = TRUE, link = "#Summary")
        heading = paste("Positions of probes for ", genes[i], 
            " on the Illumina ", annoName, " chip", sep = "")
        hwrite(heading, outPage, heading = 2, br = TRUE, name = paste(genes[i], 
            "Annotation", sep = ""))
        hwriteImage(locsPic, outPage, br = TRUE)
        fields <- c("PROBEQUALITY", "CODINGZONE", "GENOMICLOCATION", 
            "REPEATMASK", "SECONDMATCHES", "ENTREZREANNOTATED", 
            "PROBESEQUENCE")
        queryIDs <- mapToProbes(genes[i], annoName = annoName)
        qualCols <- matrix(nrow = 4, ncol = 2)
        qualCols[, 1] <- c("Perfect", "Good", "Bad", "No match")
        qualCols[, 2] <- c("#aaffaa", "#66ffff", "#ffbbaa", "#ff5500")
        probeSummary <- matrix(nrow = length(queryIDs), ncol = length(fields))
        colnames(probeSummary) <- fields
        rownames(probeSummary) <- queryIDs
        for (j in 1:length(queryIDs)) {
            probeSummary[j, ] <- unlist(sapply(fields, function(x) mget(queryIDs[j], 
                eval(as.name(paste(annoName, x, sep = ""))))))
        }
        Entrez = paste("<a href=http://www.ncbi.nlm.nih.gov/gene?term=", 
            probeSummary[, 6], ">", probeSummary[, 6], "</a>", 
            sep = "")
        probeSummary[, 6] <- Entrez
        heading = paste("Summary of annotation of probes for ", 
            genes[i], " on the ", annoName, " chip", sep = "")
        hwrite(heading, outPage, heading = 2, br = TRUE)
        hwrite(probeSummary, outPage, br = TRUE, row.bgcolor = c(NA, 
            qualCols[match(probeSummary[, 1], qualCols[, 1]), 
                2]))
        combMat <- attachCovarsToExpression(dbcon, queryIDs = queryIDs)
        factorcol <- which(colnames(combMat) == factor)
        heading <- paste("showing the expression of ", genes[i], 
            " against ", factor, sep = "")
        hwrite("Boxplots ", link = "http://en.wikipedia.org/wiki/Box_plot", 
            outPage, heading = 2)
        hwrite(heading, outPage, heading = 2, br = TRUE, name = paste(genes[i], 
            "Boxplot", sep = ""))
        message("Plotting gene ", genes[i])
        if (any(is.na(combMat[, factorcol]))) {
            combMat2 <- combMat[-which(is.na(combMat[, factorcol])), 
                ]
        }
        else combMat2 <- combMat
        pngfile <- paste(genes[i], ".png", sep = "")
        colnames(combMat2)[factorcol] <- "FactorToPlot"
        myplot <- ggplot(combMat2, aes(x = FactorToPlot, y = Expression, 
            fill = FactorToPlot)) + geom_boxplot() + xlab(factor) + 
            scale_fill_discrete(name = factor) + facet_wrap(~ProbeID, 
            ncol = 4)
        ggsave(myplot, filename = paste(outfile_path, pngfile, 
            sep = "/"), width = 6, height = 4 * ceiling(length(queryIDs)/4), 
            dpi = 100)
        hwriteImage(pngfile, outPage, br = TRUE)
        if (!is.null(externalData)) {
            for (l in 1:length(externalData)) {
                myplot <- try(boxplotFromPublicData(externalData[[l]], 
                  genes[i]), silent = TRUE)
                if (any(class(myplot) == "ggplot")) {
                  heading <- paste("Boxplot using the ", names(externalData)[l], 
                    "Dataset")
                  GEO <- as.character(dbGetQuery(cai, "SELECT value FROM Meta WHERE id = 'GSE'"))
                  GEOurl <- paste("http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", 
                    GEO, sep = "")
                  hwrite(heading, outPage, heading = 2, br = TRUE)
                  hwrite(GEO, outPage, link = GEOurl, br = TRUE)
                  ext_file <- paste(genes[i], "-", names(externalData)[l], 
                    ".png", sep = "")
                  ggsave(myplot, filename = paste(outfile_path, 
                    ext_file, sep = "/"), width = 6, height = 4, 
                    dpi = 100)
                  hwriteImage(ext_file, outPage, br = TRUE)
                }
            }
        }
        heading <- paste("Differential expression results for ", 
            genes[i], sep = "")
        hwrite(heading, outPage, heading = 2, br = TRUE, name = paste(genes[i], 
            "DE", sep = ""))
        plotFrame <- fullDE
        probeStats <- plotFrame[which(fullDE$ProbeID \%in\% queryIDs), 
            ]
        rownames(probeStats) <- 1:nrow(probeStats)
        hwrite(probeStats, outPage, br = TRUE, col.link = list(Contrasts = paste("#Contrasts:", 
            probeStats$Contrast, sep = "")))
        hwrite("Volcano plots", outPage, heading = 2, br = TRUE, 
            link = "http://en.wikipedia.org/wiki/Volcano_plot_(statistics)")
        pngfile <- paste(genes[i], "-DE.png", sep = "")
        myplot <- ggplot(fullDE, aes(x = as.numeric(LogFC), y = as.numeric(LogOdds))) + 
            geom_point(colour = "blue", alpha = 1/10) + geom_text(data = probeStats, 
            aes(x = as.numeric(LogFC), y = as.numeric(LogOdds), 
                label = ProbeID, colour = ProbeID)) + ylab("Log Odds") + 
            xlab("Log FC") + facet_wrap(~Contrasts)
        ggsave(myplot, filename = paste(outfile_path, pngfile, 
            sep = "/"), width = 8, height = 6, dpi = 100)
        hwriteImage(pngfile, outPage, br = TRUE)
    }
    deTable <- fullDE[which(fullDE$ProbeID \%in\% allIDs), ]
    deTable$LogFC = as.numeric(deTable$LogFC)
    deTable$Pvalue = as.numeric(deTable$Pvalue)
    deTable$Adjusted = as.numeric(deTable$Adjusted)
    deTable$LogOdds = as.numeric(deTable$LogOdds)
    deTable$Rank = as.numeric(deTable$Rank)
    myTable <- data.frame(ProbeID = deTable[, 1], Symbol = paste("<a href=#", 
        mappedSymbols, "DE>", mappedSymbols, "</a>", sep = ""), 
        deTable[, -1])
    hwrite("Summary of differential expression", outPage, heading = 2, 
        br = TRUE, name = "DESummary")
    hwrite("Selected Genes", outPage, heading = 2, br = TRUE)
    gTab <- gvisTable(myTable, options = list(width = 1200))
    cat(createGoogleGadget(gTab), file = outPage)
    hwrite("Details of specific contrasts", outPage, heading = 2, 
        br = TRUE)
    splitContrast <- split(fullDE, fullDE$Contrasts)
    for (nc in 1:length(splitContrast)) {
        topHits <- order(as.numeric(splitContrast[[nc]]$LogOdds), 
            decreasing = TRUE)[1:50]
        topMat <- splitContrast[[nc]][topHits, ]
        annoToGet <- c("SYMBOLREANNOTATED", "GENOMICLOCATION", 
            "ENTREZREANNOTATED")
        annoSummary <- matrix(nrow = nrow(topMat), ncol = length(annoToGet))
        colnames(annoSummary) <- annoToGet
        for (j in 1:nrow(topMat)) {
            annoSummary[j, ] <- unlist(sapply(annoToGet, function(x) mget(topMat[j, 
                1], eval(as.name(paste(annoName, x, sep = ""))))))
        }
        Symbol = paste("<a href=http://www.genecards.org/index.php?path=/Search/keyword/", 
            annoSummary[, 1], ">", annoSummary[, 1], "</a>", 
            sep = "")
        Entrez = paste("<a href=http://www.ncbi.nlm.nih.gov/gene?term=", 
            annoSummary[, 3], ">", annoSummary[, 3], "</a>", 
            sep = "")
        annoSummary[, 1] <- Symbol
        annoSummary[, 3] <- Entrez
        totalDEGenes <- sum(splitContrast[[nc]]$Adjusted < 0.05)
        cName <- unique(splitContrast[[nc]]$Contrasts)
        hwrite(cName, outPage, heading = 2, br = TRUE, name = paste("Contrasts", 
            cName, sep = ":"))
        hwrite(paste("Total Number of DE genes", totalDEGenes, 
            sep = ":"), outPage, heading = 2, br = TRUE)
        newMat <- data.frame(ProbeID = topMat$ProbeID, annoSummary, 
            LogFC = as.numeric(topMat$LogFC), Pvalue = as.numeric(topMat$Pvalue), 
            Adjusted = as.numeric(topMat$Adjusted), LogOdds = as.numeric(topMat$LogOdds))
        allDE <- splitContrast[[nc]]$ProbeID[which(splitContrast[[nc]]$Adjusted < 
            0.05)]
        mapEnv <- as.name(paste(annoName, "GENOMICLOCATION", 
            sep = ""))
        genLoc <- unlist(mget(as.character(allDE), eval(mapEnv), 
            ifnotfound = NA))
        chr <- sapply(genLoc, function(x) strsplit(x, ":", fixed = TRUE)[[1]][1])
        start <- sapply(genLoc, function(x) strsplit(x, ":", 
            fixed = TRUE)[[1]][2])
        end <- sapply(genLoc, function(x) strsplit(x, ":", fixed = TRUE)[[1]][3])
        if (any(is.na(chr))) {
            start <- start[-which(is.na(start))]
            end <- end[-which(is.na(chr))]
            chr <- chr[-which(is.na(chr))]
        }
        deRngs <- GRanges(seqnames = chr, ranges = IRanges(start = as.numeric(start), 
            end = as.numeric(end), names = names(start)))
        library(ggbio)
        data("hg19Ideogram", package = "biovizBase")
        hg19Ideo <- hg19Ideogram
        chr.sub <- paste("chr", 1:22, sep = "")
        hg19Ideo <- keepSeqlevels(hg19Ideogram, chr.sub)
        head(hg19Ideo)
        p <- ggplot() + layout_circle(hg19Ideo, geom = "ideo", 
            fill = "gray70", radius = 30, trackwidth = 4) + layout_circle(hg19Ideo, 
            geom = "text", aes(label = seqnames), vjust = 0, 
            radius = 38, trackwidth = 7) + layout_circle(deRngs, 
            geom = "rect", color = "steelblue", radius = 23, 
            trackwidth = 6)
        file <- paste(cName, "DElocations.png", sep = "")
        hwrite("Location of all DE probes", outPage, heading = 2, 
            br = TRUE)
        ggsave(p, filename = paste(outfile_path, file, sep = "/"), 
            width = 8, height = 8, dpi = 100)
        hwriteImage(file, outPage, br = TRUE)
        hwrite("Top 50 probes", outPage, heading = 2, br = TRUE)
        gTab <- gvisTable(newMat, options = list(width = 1200))
        cat(createGoogleGadget(gTab), file = outPage)
    }
    closePage(outPage)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
