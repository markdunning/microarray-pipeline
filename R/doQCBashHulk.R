doQCBashHulk <-
function(x, config,sampleSheet){
  library(beadarray)
  mouse=FALSE

#  if(config["platform",]=="Mousev2"| config["platform",]=="Mousev3"){
 #   mouse<- TRUE
 # }


  sID <- strsplit(x, "_")[[1]][1]
  secNames <- x

  if(mouse){
    secNames <- paste(x, c(1,2), sep="_")
  }
  message("Reading the specified section...")

  if(config["baseDir",] == "WORKINGDIR"){
  ##the directories are assumed to be sub-directories of the current working directory

    ##Check that the TIFF exists before trying to read
    BLData <- readIllumina(useImages=as.logical(config["useImages",]), illuminaAnnotation= config["platform",],
                sectionNames=secNames,sampleSheet=sampleSheet)

  }

  else {
    if(file.exists(config["baseDir",])){
  
    BLData <- readIllumina(dir = config["baseDir",],useImages=as.logical(config["useImages",]), illuminaAnnotation= config["platform",],
                sectionNames=secNames,sampleSheet=sampleSheet)
    

    }
    else {
      message("base directory defined in config file does not exist")
        BLData <- readIllumina(useImages=as.logical(config["useImages",]), illuminaAnnotation= config["platform",],
                sectionNames=secNames,sampleSheet=sampleSheet)
     
    }

  }


  ### UPDATE expressionQCpipeline replactes QC

  if(config["doqc",]==1){
    message("Creating QC plots...")

   tmp <- try(expressionQCPipeline(BLData,qcDir=config["qc_path",], plotType=".png", overWrite=F, zlim=NULL, tagsToDetect=NULL))
  
  }
  

  probableAnnotation <- suggestAnnotation(BLData)

  if(probableAnnotation != annotation(BLData)) message(paste("This data was read in as ", annotation(BLData), " but beadarray thinks it is ", probableAnnotation, ". Please check your config file\n",sep=""))


  
  if(config["beadRegistration",]){

    message("Making bead registration checks")
    
    dir.create(paste(config["qc_path",], "/registrationChecks/",sep=""),recursive=TRUE,showWarnings=FALSE)

    regCheck <- try(checkRegistration(BLData))

    if(class(regCheck) != "try-error"){

      png(paste(config["qc_path",],"/registrationChecks/",x,".png",sep=""),width=as.numeric(config["png_w",]),height=as.numeric(config["png_h",]),type="cairo")
      boxplot(regCheck, main=x, plotP95=TRUE)
      dev.off()
    }
    else message("Bead registration check was not successful, proceeding anyway..")

  }

### RUN BASH
  ## UPDATE - changes to BASH output type since previous pipeline

#check log parameter, will assume check how to log
  if(config["bash",]==1){
    weights <- BASH(BLData,array=1)
    BLData.bashed <- setWeights(BLData,wts=weights$wts,array=1)
    ##Bash $wts item is now a vector and not a list
    QCmat <- weights$QC
   
    BLData.bashed = insertSectionData(BLData.bashed, what = "BASHQC", data = QCmat)
  } else {
    BLData.bashed <- BLData
  }


  if(as.numeric(config["bash",])){

  #Generate Array Masks and image plots for each array
  dir.create("QA",showWarnings=FALSE)
  dir.create("beadarray_QA/imagePlots",showWarnings = FALSE,recursive=T)
  ipHt <- 1000

  png(file=paste("beadarray_QA/imagePlots/",x,".png",sep=""),width=1000,height=ipHt)


  ### Add an extra plot if we do HULK

  plots <- 2 +  as.numeric((config["hulk",]))
#  par(mfrow=c(plots,1))

  p1 <- beadarray:::showArrayMask2(BLData.bashed,array=1)
  
  logData <- logGreenChannelTransform(BLData.bashed, array=1)

  BLData.logged <-  insertBeadData(BLData.bashed, array=1, data=logData, what="Grn") 

  p2<- imageplot(BLData.logged,array=1,low="green",high="red",main=paste(x,"Before",sep="-"))

    if(config["hulk",]==1){
      hulk <- HULK(BLData.bashed, array=1)
      ##hulk return type is now a vector and not a list
      BLData.bashed <- insertBeadData(BLData.bashed, array=1, data=hulk, what="Grn") 
      p3 <- imageplot(BLData.bashed,array=1,low="green",high="red",main=paste(x,"After",sep="-"))


      gridExtra:::grid.arrange(p1,p2,p3,ncol=1)
    }
    else gridExtra:::grid.arrange(p1,p2,ncol=1)
  
  dev.off()

  }
  

  if(config["makeBabFile",]){


    txtFile <- as.character(BLData@sectionData$Targets$textFile)
    locsGreen <- as.character(BLData@sectionData$Targets$greenLocs)
    path <- as.character(BLData@sectionData$Targets$directory)

    if(length(grep(".bab", txtFile)) == 0){
  
    message("Creating .bab compressed data")

   

    txtFile <- as.character(BLData@sectionData$Targets$textFile)
    locsGreen <- as.character(BLData@sectionData$Targets$greenLocs)
    path <- as.character(BLData@sectionData$Targets$directory)
    BeadDataPackR:::compressBeadData(txtFile=txtFile, locsGrn=locsGreen, path = path)
    
    ##if we know it worked, remove the text and locs files?

    if(config["removeTextAndLocs",]){
      file.remove(c(paste(path, "/", txtFile,sep=""), paste(path, "/", locsGreen,sep="")))
    }


    }
    else message("Compressed bab file already exists...")
  
  }
  
  outfile <- paste("Robjects/", x, "-BLData.Rda", sep="")

  message("saving bead-level data and weights as ", outfile,sep="")

  ### Save the BLData object along with weights
  save(weights,BLData.bashed,file=outfile)

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
      bsd <- beadarray::summarize(BLData.bashed, list(greenChannel), useSampleFac=F)

  
  outfile <- paste("Robjects/", x, "-summary.Rda", sep="")

  save(bsd, file=outfile)
    
  outfile    


}
