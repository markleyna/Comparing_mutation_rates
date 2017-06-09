install.packages("readr")
install.packages("dplyr")
library("dplyr")
library("readr", lib.loc = "/home/french15/Rlibs")
SalCB <- read.csv("SALCBKBLAST.csv")
SalStart <- SalCB$X3568909
SalEnd<-SalCB$X3604712
CBStart<-SalCB$X4324058
CBEnd<-SalCB$X4359892
#7 and 8 are query start/end, 9 and 10 are sample start/end
SalmonellaDNA<-read_file("SalmonellaPureDNA.txt")
SalmonellaDNA<-gsub("\n","",SalmonellaDNA)
CitroBacDNA<-read_file("CitroBacKPureDNA.txt")
CitroBacDNA<-gsub("\n","",CitroBacDNA)
names<-c("gapnum","salgenestart","salgeneend","salgeneseq","cbgenestart","cbgeneend","cbgeneseq")
ExportDF<-data.frame(-1,-1,-1,"q",-1,-1,"q", stringsAsFactors = FALSE)
colnames(ExportDF)<-names
i = 0
while (i<length(CBStart)){

 i = i + 1
  row=c(i,SalStart[i],SalEnd[i],substring(SalmonellaDNA,SalStart[i],SalEnd[i]),CBStart[i],CBEnd[i],substring(CitroBacDNA,CBStart[i],CBEnd[i]))
  ExportDF<-rbind(ExportDF,row)
}
SalmonellaCBExportDF <- ExportDF[-1,]
write.csv(SalmonellaCBExportDF, "RVersionCBSalGenes.csv")

ECSAL <- read.csv("ECSALBLAST.csv")
EColiDNA <- read_file("EColiPureDNA.txt")
EColiDNA<-gsub("\n","",EColiDNA)
SalStart <- ECSAL$X3570417
SalEnd<-ECSAL$X3596716
ECStart<-ECSAL$X3427172
ECEnd<-ECSAL$X3453457
names<-c("gapnum","salgenestart","salgeneend","salgeneseq","ecgenestart","ecgeneend","ecgeneseq")
ExportDF<-data.frame(-1,-1,-1,"q",-1,-1,"q", stringsAsFactors = FALSE)
i = 0
while (i<length(ECStart)){
  
  i = i + 1
  row=c(i,SalStart[i],SalEnd[i],substring(SalmonellaDNA,SalStart[i],SalEnd[i]),ECStart[i],ECEnd[i],substring(EColiDNA,ECStart[i],ECEnd[i]))
  ExportDF<-rbind(ExportDF,row)
}
SalEC<-ExportDF[-1,]
write.csv(SalEC,"RVersionSalECGenes.csv")

CBEC <- read.csv("ECCBKBLAST.csv")
EColiDNA <- read_file("EColiPureDNA.txt")
EColiDNA<-gsub("\n","",EColiDNA)
CBstart <- CBEC$X3427172
CBEnd<-CBEC$X3453454
ECStart<-CBEC$X4325524
ECEnd<-CBEC$X4351821
names<-c("gapnum","cbgenestart","cbgeneend","cbgeneseq","ecgenestart","ecgeneend","ecgeneseq")
ExportDF<-data.frame(-1,-1,-1,"q",-1,-1,"q", stringsAsFactors = FALSE)
colnames(ExportDF)<-names
i = 1
while (i<length(ECEnd)){
  
  i = i + 1
  row=c(i,CBstart[i],CBEnd[i],substring(CitroBacDNA,CBstart[i],CBEnd[i]),ECStart[i],ECEnd[i],substring(EColiDNA,ECStart[i],ECEnd[i]))
  ExportDF<-rbind(ExportDF,row)
}
CBEC<-ExportDF[-1,]
CBEC<- subset(CBEC, (nchar(cbgeneseq)>1000))
CBEC<- subset(CBEC, (nchar(ecgeneseq)>1000))
write.csv(CBEC,"RVersionCBECGenes.csv")

#SalECGaps

#very ugly search/sort idea
#sweep through one end and find the shortest distance to a start and check the other pair and if both <500, then tag as a gap

GappingCBEC <- unique(CBEC[order(as.numeric(CBEC$cbgenestart)),])
q= 1
#gapgeneration <- GappingCBEC[1,]
# while (q<length(GappingCBEC$gapnum)){
#   r = 1
#   distance = 9999999999
#   while ( r<length(GappingCBEC$gapnum)){
#   
#   newdistance= sqrt(((as.numeric(GappingCBEC$ecgenestart[r])-as.numeric(GappingCBEC$ecgeneend[q]))^2) + ((as.numeric(GappingCBEC$cbgenestart[r])-as.numeric(GappingCBEC$cbgeneend[q]))^2))
#   if ((newdistance > 0) && (as.numeric(GappingCBEC$cbgenestart[r])>as.numeric(GappingCBEC$cbgeneend[q])) && (as.numeric(GappingCBEC$ecgenestart[r])>as.numeric(GappingCBEC$ecgeneend[q]))){
#     if(newdistance<distance){
#       distance = newdistance
#       bestrow = r
#       print(bestrow)
#     }
#     
#   }
#   r=r+1
#   
#   }
#   newrow = GappingCBEC[bestrow,]
#   #GappingCBEC<-GappingCBEC[-bestrow,]
#   gapgeneration = rbind(gapgeneration, newrow)
# #   q = q+1
# # } 
# gapgeneration<-GappingCBEC
# gapgeneration<-unique(gapgeneration)
# s = 1
# 
# outputgaps<-gapgeneration[1,]
# while (s<length(gapgeneration$gapnum)-1){
#   if (as.numeric(gapgeneration$cbgeneend[s])-as.numeric(gapgeneration$cbgenestart[s+1])<2000){
#     if(as.numeric(gapgeneration$cbgeneend[s])-as.numeric(gapgeneration$cbgenestart[s+1])>0){
#       if(as.numeric(gapgeneration$ecgeneend[s])-as.numeric(gapgeneration$ecgenestart[s+1])<2000){
#         if(as.numeric(gapgeneration$ecgeneend[s])-as.numeric(gapgeneration$ecgenestart[s+1])>0){
#            gap = c(gapgeneration$cbgeneend[s+1],gapgeneration$cbgenestart[s],substring(CitroBacDNA,gapgeneration$cbgeneend[s+1],gapgeneration$cbgenestart[s]),gapgeneration$ecgeneend[s+1],gapgeneration$ecgenestart[s],substring(EColiDNA,gapgeneration$ecgeneend[s+1],gapgeneration$ecgenestart[s]))
#            outputgaps<-rbind(outputgaps,gap)  
#         }
#       }
#     }
#     
#   }
#   s=s+1    
# }
myDF
factor <- 1

gapIden <- function(cbend, cbstart, ecend, ecstart) {
  diff <- cbstart - cbend
  if ((ecstart - ecend) <= (factor * diff)) {
    tempDF <- data.frame(cbend, cbstart, ecend, ecstart)
    return(tempDF)
    #return(TRUE)
  }
  else {
    return(NA)
    #return(FALSE)
  }
}

gapList <- mapply(gapIden, myDF$cbgeneend[-1], myDF$cbgenestart[-608], myDF$ecgeneend[-1], myDF$ecgenestart[-608])
gapDF <-do.call(rbind, gapList) 
gapDF <- gapDF[complete.cases(gapDF),]
write.csv(gapDF, "UpdatedCBECgaps.csv")
class(myDF[1,])
colnames(myDF[1,])

