exampleControlPlots <-
function(BSData, config){

    snr <- qcData(BSData)[,"P95Grn"] / qcData(BSData)[,"P05Grn"]
    
    worstArray <- sampleNames(BSData)[which.min(snr)]

    bld <- readIllumina(sectionNames = worstArray, illuminaAnnotation=annotation(BSData))
       
    p1 <- combinedControlPlot(bld) + opts(title=paste("Control plot for array with signal-to-noise ratio of", round(min(snr),3)))

    bestArray <- sampleNames(BSData)[which.max(snr)]

    bld <- readIllumina(sectionNames = bestArray, illuminaAnnotation=annotation(BSData))
  
    p2 <- combinedControlPlot(bld) + opts(title=paste("Control plot for array with signal-to-noise ratio of", round(max(snr),3)))
  
    png("QA/controlPlotsComparison.png",width=1200,height=1200)
    gridExtra::grid.arrange(p1,p2,ncol=1)
    dev.off()


   if(config["bash",]==1){

    BASHMetrics <- qcData(BSData)[,"BeadsMasked"]

    worstArray <- sampleNames(BSData)[which.max(BASHMetrics)]

    load(paste("Robjects/", worstArray, "-BLData.Rda",sep=""))

    pctMasked <- round(max(BASHMetrics) / numBeads(bld),3)*100

    p1 <- beadarray:::imageplot(BLData.bashed) +opts(title=paste("Array most affected by spatial artefacts, with", pctMasked, "% of beads masked by BASH"))
    p2 <- beadarray:::showArrayMask2(BLData.bashed,array=1)

    bestArray <- sampleNames(BSData)[which.min(BASHMetrics)]
  
    load(paste("Robjects/", bestArray, "-BLData.Rda",sep=""))

    pctMasked <- round(min(BASHMetrics) / numBeads(bld),3)*100
    p3 <- beadarray:::imageplot(BLData.bashed) +opts(title=paste("Array least affected by spatial artefacts, with", pctMasked, "% of beads masked by BASH"))
    p4 <- beadarray:::showArrayMask2(BLData.bashed,array=1)


    png("QA/spatialArtefactComparison.png",width=1200,height=800)
    gridExtra::grid.arrange(p1,p3,p2,p4,ncol=2)
    dev.off()	
  
  }
    

}
