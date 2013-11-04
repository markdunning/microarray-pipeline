
createCheckSums <- function(sampleSheet){
library(tools)

checksum <- NULL
checked <- NULL
for(i in 1:nrow(sampleSheet)){


	chip <- sampleSheet$Sentrix_ID[i]
	section <- sampleSheet$Sentrix_Position[i]
	dirFiles <- list.files(as.character(chip))
	toCheck <- paste(chip, dirFiles[grep(paste(chip, section,sep="_"),dirFiles)],sep="/")	
	checksum <- c(checksum, sapply(toCheck, md5sum))
	checked <- c(checked, toCheck)

}

data.frame(checked, checksum)

}
