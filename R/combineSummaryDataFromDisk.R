combineSummaryDataFromDisk <-
function(fileNames, config){

annoName <- config["platform",]

annoLoaded <- require(paste("illumina", annoName, ".db",sep=""), character.only=TRUE)


if(all(sapply(fileNames, function(x) file.exists(x)))){

  if(annoLoaded){


    mapEnv <-  as.name(paste("illumina", annoName, "ARRAYADDRESS",sep=""))

    IlluminaIDs = as.character(unlist(mappedkeys(x=eval(mapEnv))))

    
    eMat <- seMat <- nObsMat <- matrix(nrow = length(IlluminaIDs), ncol = length(fileNames))

    rownames(eMat) <- rownames(seMat) <- rownames(nObsMat) <- IlluminaIDs

  }


  load(fileNames[1])

#  pData <- data.frame(nrow = length(fileNames), ncol=ncol(pData(bsd)))
#  qcData <- data.frame(nrow = length(fileNames), ncol=ncol(qcData(bsd)))

  pData <- phenoData(bsd)@data

  qcData <- bsd@QC@data

      eMat[,1] <- exprs(bsd)[match(IlluminaIDs, featureNames(bsd)),]

      seMat[,1] <- se.exprs(bsd)[match(IlluminaIDs, featureNames(bsd)),]

      nObsMat[,1] <- nObservations(bsd)[match(IlluminaIDs, featureNames(bsd)),]
    

	      
  ###Create template matrices to store as many arrays as the number of files


    sampleNames <- sampleNames(bsd)

    for(i in 2:length(fileNames)){

      load(fileNames[i])
    
      sampleNames = c(sampleNames, sampleNames(bsd))

      eMat[,i] <- exprs(bsd)[match(IlluminaIDs, featureNames(bsd)),]

      seMat[,i] <- se.exprs(bsd)[match(IlluminaIDs, featureNames(bsd)),]

      nObsMat[,i] <-  nObservations(bsd)[match(IlluminaIDs, featureNames(bsd)),]
    
      pData <- rbind(pData, pData(bsd))
  
      qcData <- merge(qcData, qcData(bsd), all=TRUE,sort=FALSE)


      }

  colnames(eMat) <- colnames(seMat) <- colnames(nObsMat) <- sampleNames

  BSData <- new("ExpressionSetIllumina")
  assayData(BSData)=assayDataNew(exprs = eMat, se.exprs = seMat, nObservations = nObsMat,storage.mode="list")

  p = new("AnnotatedDataFrame", data.frame(pData, row.names=sampleNames))
 
  phenoData(BSData) = p    


  QC = new("AnnotatedDataFrame", data.frame(qcData, row.names = sampleNames))

  annotation(BSData) <- annoName

  BSData@QC = QC

  BSData@channelData <- list()
  BSData@channelData[[1]] <- rep("G", length(sampleNames))

  status = rep("Unknown", length(IlluminaIDs))

  mapEnv <-  as.name(paste("illumina", annoName, "REPORTERGROUPNAME",sep=""))

  t <- try(eval(mapEnv),silent=TRUE)

  if(class(t) == "try-error"){
    message(paste("Could not find a REPORTERGROUPNAME mapping in annotation package ", annoPkg,". Perhaps it needs updating?", sep=""))

  }
  
  else{

  status[which(!is.na(IlluminaIDs))] = unlist(mget(IlluminaIDs[which(!is.na(IlluminaIDs))], eval(mapEnv), ifnotfound=NA))	

  status[which(is.na(status))] = "regular"

  }

  featureData(BSData) = new("AnnotatedDataFrame", data=data.frame(IlluminaID  = IlluminaIDs, Status = status, row.names=IlluminaIDs))

  BSData@annotation = annoName

  save(BSData,file="Robjects/Combined_BSData.Rda")
    
  }

else (message("Fatal error; Not all summary files could be found\n"))

BSData

}
