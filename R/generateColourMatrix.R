generateColourMatrix <-
function(targets){

col_matrix <- matrix(nrow=nrow(targets),ncol=ncol(targets),dimnames=list(row.names(targets),colnames(targets)))
reject<-vector()
for (c in 1:length(colnames(targets))){
	types <- levels(as.factor(targets[,c]))
if (colnames(col_matrix)[c] =="Sentrix_ID" ||length(types) >1 && length(types)!=nrow(targets) && colnames(col_matrix)[c] !="RIN" && colnames(col_matrix)[c] !="registrationScore" && colnames(col_matrix)[c] !="focusScore" && colnames(col_matrix)[c] !="Controls"){
          if (length(types)==2){
            cols <- brewer.pal(3,config["colScheme",])[c(1,2)]
          }else{
            cols <- colorRampPalette(brewer.pal(8,config["colScheme",]))(length(types))
          }
	col_matrix[,c] <- cols[match(targets[,c],types)]
	}else {
		reject[length(reject)+1] <- c
	}	
	}
if (length(reject) !=0){
	col_matrix <- matrix(col_matrix[,-reject],nrow=nrow(col_matrix))
	row.names(col_matrix) <- row.names(targets)
	colnames(col_matrix) <- colnames(targets)[-reject]
}

#if (length(xmlFiles)!=0){
#	col_matrix <- cbind(col_matrix, registrationScore = rgb(colorRamp(c("Red","Green"))(targets$registrationScore),maxColorValue=255))
#	#col_matrix <- cbind(col_matrix, focusScore = rgb(colorRamp(c("Red","Green"))(targets$focusScore),maxColorValue=255))
#	targets$focusScore <- round(as.numeric(targets$focusScore)*100,0)
#	col_matrix <- cbind(col_matrix, focusScore = redgreen(diff(range(targets$focusScore))+1)[targets$focusScore-(min(targets$focusScore)-1)])
#}

col_matrix 


}
