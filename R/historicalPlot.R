historicalPlot <-
function(MetricsData, metrics, addToHist = FALSE){

if (file.exists(MetricsData)){
cat("MetricsData File Found\n")
oldMetrics <- read.delim(MetricsData)
cat(dim(oldMetrics),"\n")
oldMetrics$Date <- sub(x=oldMetrics$Date,pattern="/09",replacement="/2009")
tmp1 <- strsplit(sapply(strsplit(as.character(oldMetrics$Date)," "),"[",1),"/")

HalfYear <- paste(sapply(tmp1,"[",3),cut(as.numeric(sapply(tmp1,"[",1)),breaks=2,labels=F),sep="-")
SentrixID <- as.character(metrics$Matrix)

oldMetrics <- data.frame(oldMetrics, HalfYear)
metrics <- data.frame(metrics, SentrixID)

p <- ggplot(oldMetrics,aes(x=P95Grn,y=P05Grn))
p <- p + geom_point(aes(color=HalfYear))
p <- p + scale_color_brewer(palette="Paired")

p <- p + geom_point(data=metrics,aes(x=P95Grn,y=P05Grn,shape=SentrixID))
p <- p + geom_abline(intercept=0,slope=0.1,color="Red",size=1)

make_png("QA/P95/P95historical.png");show(p);dev.off()

if (as.logical(addToHist)){
write.table(rbind(oldMetrics,metrics),file=MetricsData,row.names=T,col.names=T,quote=F,sep="\t")
}

}else {
###If we  can not find the historical metrics file then carry on, the P95historical.png will not be made though
  cat("MetricsData File not found\n")
}


}
