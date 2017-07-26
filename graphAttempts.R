#install.packages("plotly")
#install.packages("RColorBrewer")
library("RColorBrewer")
library("plotly")
library(readr)
Genes <- read.csv("/home/nmarkle/Comparing_mutation_rates/CBSalGapData.csv", stringsAsFactors = FALSE)
print(head(Genes$X))

try <- strsplit(Genes$X, " ")
mutation1 <- unlist(try)[2*(1:nrow(Genes))-1]
mutation2 <- unlist(try)[2*(1:nrow(Genes))]

sanityC1 <- data.frame(paste(substr(mutation1,0,1), substr(mutation2,0,1),sep=""))
sanityC2 <- data.frame(paste(substr(mutation1,2,2), substr(mutation2,2,2),sep=""))
sanity<-data.frame(sanityC1,sanityC2,Genes$x)
colnames(sanity)<-c("first","second",'Freq')

workingSet <- subset(sanity, substr(first,0,1)=="T" & substr(first,2,2)!="T")
workingSet <- subset(workingSet, substr(second,2,2)!=substr(second,0,1))

workingSet[workingSet$first=="T-","Freq"] <- (6.2-as.numeric(as.character((workingSet$Freq[workingSet$first=="T-"]))))/6.2
workingSet[workingSet$first=="TA","Freq"] <- (14.5-as.numeric(as.character((workingSet$Freq[workingSet$first=="TA"]))))/14.5
workingSet[workingSet$first=="TC","Freq"] <- (28.2-as.numeric(as.character((workingSet$Freq[workingSet$first=="TC"]))))/28.2
workingSet[workingSet$first=="TG","Freq"] <- (17.05-as.numeric(as.character((workingSet$Freq[workingSet$first=="TG"]))))/17.05
workingSet$Freq <- workingSet$Freq*10

attemptSet<-workingSet[order(workingSet$second),]
workingSet <- attemptSet

xlab<-c("T-","TA","TC","TG")
i<-1
rownames(attemptSet)<-1:nrow(attemptSet)
ylab<-as.vector(attemptSet$second[seq(1,length(attemptSet$second),4)])
m<-matrix(workingSet$Freq,nrow=4, ncol = 20)
p<-plot_ly(x=ylab,y=xlab,z=m,colors=colorRamp(c("blue","red")), type="heatmap", colorbar=list(title="Percent Away\nfrom Uniform")) %>%
  layout(title="Percentage Away From Uniform For Two Mutations in Salmonella and Citrobacter",
         xaxis=list(title = "Second Mutation"),
         yaxis=list(title="First Mutation"))
ggplotly(p) 

#,marker=list(color=colorRampPalette(brewer.pal(9,"OrRd"))(100))
