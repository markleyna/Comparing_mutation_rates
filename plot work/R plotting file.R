CBECGene<-read.csv("CBECGeneData.csv")
CBECGap<-read.csv("CBECGapData.csv")
CBSalGene<-read.csv("CBSalGeneData.csv")
CBSalGap<-read.csv("CBSalGapData.csv")
SalECGene<-read.csv("SalECGeneData.csv")
SalECGap<-read.csv("ECSalGapData.csv")

#filters to apply to data
drops<-c("Adash","dashA","Tdash","dashT","Cdash","dashC","Gdash","dashG")
unif<-c(80,80,80,80,80,80,80,80,80,80,80,80)


#Testing CBEC
mx1<-as.vector(data.matrix(CBECGene))
mx2<-as.vector(data.matrix(CBECGap))
mx<-rbind(mx1,mx2)
chisq.test(mx)
#with drops
mx1d<-as.vector(data.matrix(CBECGene[ , !(names(CBECGene) %in% drops)]))
mx2d<-as.vector(data.matrix(CBECGap[ , !(names(CBECGap) %in% drops)]))
mxd<-rbind(mx1d,mx2d)
chisq.test(mxd)
#against unif
mx1du<-rbind(mx1d,unif)
chisq.test(mx1du)
mx2du<-rbind(mx2d,unif)
chisq.test(mx2du)
#pie charts
pie(as.matrix(CBECGap), labels=names(CBECGap))
pie(as.matrix(CBECGene), labels=names(CBECGene))

#gap/gene comparison plots
test<-as.matrix(rbind(CBECGap,CBECGene))
percentagetest<-test/rowSums(test)
barplot(test,beside=T)
barplot(percentagetest,beside=TRUE)
#with drops
droptest<-as.matrix(rbind(CBECGap[ , !(names(CBECGene) %in% drops)],CBECGene[ , !(names(CBECGap) %in% drops)]))
barplot(droptest,beside=T)
percentagedroptest<-droptest/rowSums(droptest)
jpeg("CBECPercentage.jpeg", width = 1920, height=1080)
barplot(percentagedroptest, beside=T, main="Citrobacter/EColi comparison", cex.main=7, cex.axis = 4, cex.names = 4)
dev.off()




#Testing CBSal
mx1<-as.vector(data.matrix(CBSalGene))
mx2<-as.vector(data.matrix(CBSalGap))
mx<-rbind(mx1,mx2)
chisq.test(mx)
#with drops
mx1d<-as.vector(data.matrix(CBSalGene[ , !(names(CBSalGene) %in% drops)]))
mx2d<-as.vector(data.matrix(CBSalGap[ , !(names(CBSalGap) %in% drops)]))
mxd<-rbind(mx1d,mx2d)
chisq.test(mxd)
#against unif
mx1du<-rbind(mx1d,unif)
chisq.test(mx1du)
mx2du<-rbind(mx2d,unif)
chisq.test(mx2du)
#pie charts
pie(as.matrix(CBSalGap), labels=names(CBSalGap))
pie(as.matrix(CBSalGene), labels=names(CBSalGene))

#gap/gene comparison plots
test<-as.matrix(rbind(CBSalGap,CBSalGene))
percentagetest<-test/rowSums(test)
barplot(test,beside=T)
barplot(percentagetest,beside=TRUE)
#with drops
droptest<-as.matrix(rbind(CBSalGap[ , !(names(CBSalGene) %in% drops)],CBSalGene[ , !(names(CBSalGap) %in% drops)]))
barplot(droptest,beside=T)
percentagedroptest<-droptest/rowSums(droptest)
jpeg("CBSalPercentage.jpeg", width = 1920, height=1080)
barplot(percentagedroptest, beside=T, main="Salmonella/Citrobacter comparison", cex.main=7, cex.axis = 4, cex.names = 4)
dev.off()





#Testing SalEC
mx1<-as.vector(data.matrix(SalECGene))
mx2<-as.vector(data.matrix(SalECGap))
mx<-rbind(mx1,mx2)
chisq.test(mx)
#with drops
mx1d<-as.vector(data.matrix(SalECGene[ , !(names(SalECGene) %in% drops)]))
mx2d<-as.vector(data.matrix(SalECGap[ , !(names(SalECGap) %in% drops)]))
mxd<-rbind(mx1d,mx2d)
chisq.test(mxd)
pie(as.matrix(SalECGap), labels=names(CBSalGap))
pie(as.matrix(SalECGene), labels=names(CBSalGene))


#gap/gene comparison plots
test<-as.matrix(rbind(SalECGap,SalECGene))
percentagetest<-test/rowSums(test)
barplot(test,beside=T)
barplot(percentagetest,beside=TRUE)
#with drops
droptest<-as.matrix(rbind(SalECGap[ , !(names(SalECGene) %in% drops)],SalECGene[ , !(names(SalECGap) %in% drops)]))
barplot(droptest,beside=T)
percentagedroptest<-droptest/rowSums(droptest)
jpeg("SalECPercentage.jpeg", width = 1920, height=1080)
barplot(percentagedroptest, beside=T, main="Salmonella/EColi comparison", cex.main=7, cex.axis = 4, cex.names = 4)
legend(1700,y=9000,legend=c("Salmonella","EColi"))
dev.off()

