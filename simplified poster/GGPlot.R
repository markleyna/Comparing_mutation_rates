#install.packages("ggplot2")
library(ggplot2)
CBSalGap<-read.csv("CBSalGapData.csv")
CBSalGene<-read.csv("CBSalGeneData.csv")
CBDF<-rbind(CBSalGap,CBSalGene)
ggplot(CBDF)+geom_histogram(binwidth = 1)
