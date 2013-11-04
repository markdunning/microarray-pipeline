summarizeFromDisk <-
function(arrNames){

  for(i in 1:length(arrNames))
  {
    load(paste("Robjects/", arrNames[i], "-BLData.Rda", sep=""))
    if(config["bash",]==1){
      if(mouse){
        BASHMetrics = matrix(nrow=1, ncol=4)
        rownames(BASHMetrics) = arrNames[i]
        BASHMetrics[,1] = unlist(lapply(weights, function(x) length(which(x$wts==0))))[1]
        BASHMetrics[,2] = unlist(lapply(weights, function(x) length(which(x$wts==0))))[2]
        BASHMetrics[,3] = weights[[1]]$ext[1]
        BASHMetrics[,4] = weights[[2]]$ext[2]
        colnames(BASHMetrics) = c("Number Masked 1", "Number Masked 2", "Extended Score 1", "Extended Score 2")
      } else {
        BASHMetrics = matrix(nrow=1, ncol=2)
        rownames(BASHMetrics) = arrNames[i]
        BASHMetrics[,1] = length(which(weights$wts==0)) ##$wts is now a vector and not a list
        BASHMetrics[,2] = weights$ext
        colnames(BASHMetrics) = c("Number Masked", "Extended Score")
      }

    }
    if(i==1){
      BLData <- BLData.bashed
      #unlink(paste("Robjects/", arrayNames[i], "-BLData.Rda", sep="")) if $rmbl;
      ### UPDATE - names/location metrics changed, BLData_at_sectiondata 
      if(config["bash",]==1){
      write.table(BASHMetrics, file=paste(config["qc_path",],"BASHMetrics.csv",sep=""),
                   row.names=TRUE,col.names=TRUE,sep=",",append=F)
      }
     } else {
      ### gdata includes a function called combine which takes precedence. However, Biobase:::combine somehow works out that its the beadarray method that we want
      BLData <- Biobase:::combine(BLData,BLData.bashed)
      # "unlink(\"$_\")\n" if $rmbl;
      if(config["bash",]==1){
      write.table(BASHMetrics, file=paste(config["qc_path",],"BASHMetrics.csv",sep=""),
                    row.names=TRUE,col.names=FALSE,sep=",",append=T)
      }
    }
    print(paste(i, "/", length(arrNames), "....done"))
  }


  ##Do plots of the best and worst examples 


  BASHMetrics <- read.csv("beadarray_QA/BASHMetrics.csv")

  worstArray <- which.max(BASHMetrics[,1])

  pctMasked <- round(BASHMetrics[worstArray,1] / numBeads(BLData)[worstArray],3)*100
  p1 <- beadarray:::imageplot(BLData, worstArray) +opts(title=paste("Array most affected by spatial artefacts, with", pctMasked, "% of beads masked by BASH"))
  p2 <- beadarray:::showArrayMask2(BLData,array=worstArray)

  bestArray <- which.min(BASHMetrics[,1])

  pctMasked <- round(BASHMetrics[bestArray,1] / numBeads(BLData)[bestArray] ,3)*100
  p3 <- beadarray:::imageplot(BLData, bestArray) +opts(title=paste("Array least affected by spatial artefacts, with", pctMasked, "% of beads masked by BASH"))
  p4 <- beadarray:::showArrayMask2(BLData,array=bestArray)


  png("QA/spatialArtefactComparison.png",width=1200,height=800)
  gridExtra::grid.arrange(p1,p3,p2,p4,ncol=2)
  dev.off()

  metrics <- read.csv("beadarray_QA/scanMetrics.csv")

  snr <- metrics$P95Grn / metrics$P05Grn

  worstArray <- which.min(snr)
  
  p1 <- combinedControlPlot(BLData, array = worstArray) + opts(title=paste("Control plot for array with signal-to-noise ratio of", round(snr[worstArray],3)))

  bestArray <- which.max(snr)
  
  p2 <- combinedControlPlot(BLData, array = bestArray) + opts(title=paste("Control plot for array with signal-to-noise ratio of", round(snr[bestArray],3)))
  

  png("QA/controlPlotsComparison.png",width=1200,height=1200)
  gridExtra::grid.arrange(p1,p2,ncol=1)
  dev.off()
  
  if(config["blmerge",]==1){
  save(BLData,file="Robjects/Combined_BLData.Rda")
  }
  
  ### UPDATE expressionQCpipeline replaces QC
  if(config["doqc",]==1){
    forDel <- paste(config["qc_path",],c("detectionMetrics.csv", "outlierMetrics.csv", "probeMetrics.csv", "scanMetrics.csv", "Boxplot.png", "Summary.htm"),sep="")
    unlink(forDel)
    ##\problem of Mousev2 not having 'housekeeping' controls
    if(mouse) expressionQCPipeline(BLData,qcDir=config["qc_path",], plotType=".png", overWrite=F,tagsToDetect = list(Biotin = "biotin"))
    else expressionQCPipeline(BLData,qcDir=config["qc_path",], plotType=".png", overWrite=F, tagsToDetect=NULL)
#    expressionQCPipeline(BLData,qcDir=config["qc_path",], plotType=".png", overWrite=F)
 
 }
  
  ### createBeadSummaryData now caled summarize - needs a function to merge and log data
  if(config["bsmerge",]==1){
    myMean = function(x) mean(x, na.rm = TRUE)
    mySd = function(x) sd(x, na.rm = TRUE)
    fooFun = function (BLData, array){
      x = getBeadData(BLData, array = array, what = "Grn")
    }
    greenChannel = new("illuminaChannel", fooFun, illuminaOutlierMethod, myMean, mySd, "G")
    
    ##Taking out the bit that uses a sample factor when summarizing
    #if(config["platform",]=="Mousev2"| config["platform",]=="Mousev3"){
    #  BSData <- beadarray::summarize(BLData, list(greenChannel), useSampleFac=T, sampleFac=as.vector(strtrim(sectionNames(BLData),12)))
   # }
    #else{
      BSData <- beadarray::summarize(BLData, list(greenChannel), useSampleFac=F)
   # }
    save(BSData,file="Robjects/Combined_BSData.Rda")
  }
  #source("metricsPlots.R") ### RUN Basic Metrics Plots


}
