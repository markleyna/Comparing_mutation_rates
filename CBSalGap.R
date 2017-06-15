#Import Packages
#install.packages("stringr")
#install.packages("readr")
library(stringr)
library(parallel)
library(lattice)
library(MASS)
#library('readr')
library('readr', lib.loc="/home/nmarkle/Rlibs/")

#Variable Declaration Block
numInRow <- 2
Genes <- read.csv("/home/nmarkle/Comparing_mutation_rates/UpdatedSalCBgaps.csv", stringsAsFactors = FALSE)
CitroBacDNA<-read_file("/home/nmarkle/Comparing_mutation_rates/CitroBacKPureDNA.txt")
CitroBacDNA<-gsub("\n","",CitroBacDNA)
SalDNA <- read_file("/home/nmarkle/Comparing_mutation_rates/SalmonellaPureDNA.txt")
SalDNA<-gsub("\n","",SalDNA)
Genes$salgeneseq<-substring(SalDNA,first = Genes$cbend,last = Genes$cbstart)
Genes$cbgeneseq<-substring(CitroBacDNA,first = Genes$ecstart, last = Genes$ecend)
#Genes <- data.frame(Genes$X, Genes$cbgeneseq, Genes$ecgeneseq, stringsAsFactors = FALSE)
#testDF <- subset(Genes, Genes$cbstart > Genes$cbend)
Genes<-subset(Genes, nchar(Genes$salgeneseq) <= 800)
Genes<-subset(Genes, nchar(Genes$cbgeneseq) <= 800)
Genes$X<-seq(1,nrow(Genes),1)
#-----------Testing Line----------------
#Genes <- head(Genes, n=1)
#---------------------------------------
G1<-as.vector(Genes$X)
G2<-as.vector(Genes$cbgeneseq)
G3<-as.vector(Genes$salgeneseq)

buildchar = "Q" #done to prevent things from crashing later on, don't ask
Gene1Base <- c(buildchar,buildchar)
Gene2Base<-c(buildchar,buildchar) 
PreCount<-c(-1,-1)
PostCount<-c(-1,-1)
DF = data.frame(Gene1Base,Gene2Base,PreCount,PostCount)
DF <- data.frame(lapply(DF, as.character), stringsAsFactors=FALSE) #done to stop shenanigans with strings later on
Gene1RunTitle <- paste(">","ECGap_number","k", sep = "") #must be 2 digit representation for code to work
Gene2RunTitle <- paste(">","SalGap_number","k", sep = "") #must be 3 digit representation for code to work

#CLUSAL Run on Genes
clusal_run <- function(testNum, Gene1, Gene2) {
  Gene1Base <- c(buildchar,buildchar)
  Gene2Base<-c(buildchar,buildchar) 
  PreCount<-c(-1,-1)
  PostCount<-c(-1,-1)
  Gene1break = ""
  Gene2break = ""
  tempDF = data.frame(Gene1Base,Gene2Base,PreCount,PostCount)
  tempDF <- data.frame(lapply(tempDF, as.character), stringsAsFactors=FALSE) #done to stop shenanigans with strings later on
  if (Gene1 == "" || Gene2 == "") {
    return(tempDF)
  }

  Gene1Test <- toString(Gene1)
  Gene2Test<-toString(Gene2)
  Gene1RunTitle <- paste(">","ECGap_number","k", sep = "") #I don't feel like rewriting this to hold generality, so these need to stay as is to keep the parsing working correctly
  Gene2RunTitle <- paste(">","SalGap_number","k", sep = "")
  
  #write a text file that is in proper FASTA format
  fileName <- paste("FastaInGCS", testNum, ".txt", sep="")
  fileConn<-file(fileName) #reset the file to being blank text again
  writeLines("\n",fileConn)
  close(fileConn)
  sink(fileName)
  cat(Gene1RunTitle)
  cat("\n")
  cat(Gene1Test)
  cat("\n")
  cat(Gene2RunTitle)
  cat("\n")
  cat(Gene2Test)
  sink()
  #Call the clustalw function
  systemCall = paste("clustalw2 -infile=", fileName, " -type=DNA", sep="")
  system(systemCall)
  #pull the gapped strings out of the .aln file
  alnlines<-readLines(paste("FastaInGCS", testNum, ".aln", sep=""))
  #Test Line for llc server
  alnlines<-readLines("/home/nmarkle/Comparing_mutation_rates/Tests/FastaTest.aln")
  
  alnlines1<-grep(substring(Gene1RunTitle,2), alnlines, value=TRUE)
  alnlines2<-grep(substring(Gene2RunTitle,2), alnlines, value=TRUE)
  alnlines1<-gsub(substring(Gene1RunTitle,2), "", alnlines1)
  alnlines1<-gsub(" ", "", alnlines1)
  alnlines2<-gsub(substring(Gene2RunTitle,2), "", alnlines2)
  alnlines2<-gsub(" ", "", alnlines2)
  alnlines1<-paste(alnlines1, collapse='')
  alnlines2<-paste(alnlines2, collapse='')
  Gene1gapstring <- alnlines1
  Gene2gapstring <- alnlines2
  
  s = 1
  Flag1 = FALSE
  Flag2 = FALSE
  PreCounter = 0
  PostCounter =0
  #start of intra-gene break isolation
  while (s <= nchar(Gene1gapstring)){
    Gene1workingcharacter = substring(Gene1gapstring,s,s)
    Gene2workingcharacter = substring(Gene2gapstring,s,s)
    s = s+1
    if(Gene1workingcharacter != Gene2workingcharacter && Flag2 == TRUE && PostCounter > 0)
    {
      if (nchar(Gene1break) == numInRow) {
        Addrow = c(Gene1break,Gene2break,PreCounter,PostCounter)
        tempDF = rbind(tempDF,Addrow)
      }
      Flag2=TRUE
      Flag1=FALSE
      PreCounter=PostCounter
      PostCounter=0
      Gene1break = Gene1workingcharacter
      Gene2break = Gene2workingcharacter
    }
    else if (Gene1workingcharacter != Gene2workingcharacter && Flag2 == TRUE && PostCounter == 0) {
      if (nchar(Gene1break) < numInRow) {
        Gene1break <- paste(Gene1break, Gene1workingcharacter, sep="")
        Gene2break <- paste(Gene2break, Gene2workingcharacter, sep="")
      }
      else {
        PreCounter = 0
        Flag2 = FALSE
        Flag1 = FALSE
        Gene1break = ""
        Gene2break = ""
      }
    }
    else if(Gene1workingcharacter == Gene2workingcharacter && Flag2==TRUE)
    {
      PostCounter = PostCounter +1
    }
    else if(Gene1workingcharacter != Gene2workingcharacter && Flag1==TRUE)
    {
      Flag1 = FALSE
      Flag2 = TRUE
      Gene1break = Gene1workingcharacter
      Gene2break = Gene2workingcharacter
    }
    else if(Gene1workingcharacter == Gene2workingcharacter && Flag1==TRUE)
    {
      PreCounter = PreCounter+1
    }
    else if(Gene1workingcharacter== Gene2workingcharacter && Flag1==FALSE)
    {
      PreCounter = PreCounter + 1
      Flag1 = TRUE
    }
    
    
    
  }
  #end of intra-gene break identification
  return(tempDF)
}

DF <- mcmapply(clusal_run, G1, G2, G3, SIMPLIFY = FALSE, mc.cores = 8)
DF <- do.call(rbind,DF)

DF$PostCount <- as.integer(DF$PostCount)
DF$PreCount <- as.integer((DF$PreCount))
ModDF = DF #DF is the entries from tallying up the mismatches
ModDF = subset(ModDF, PostCount > 0) #clears up things where an error got through
PlotData = subset(ModDF, PostCount>=3 & PreCount>=3) #this is the arbitry definition of what we considered conserved, this can change
#PlotData = DF

#hist(as.numeric(PlotData$PreCount)) #plotting up everything
#hist(as.numeric(PlotData$PostCount))
#hist(as.numeric(ModDF$PreCount))
#hist(as.numeric(ModDF$PostCount))
#plot(jitter(as.numeric(DF$PreCount)),jitter(as.numeric(DF$PostCount)))
#plot(as.numeric(PlotData$PreCount),as.numeric((PlotData$PostCount)))

make_graphs <- function(finalDF, firstM) {
GeneAdash = subset(DF, (Gene1Base == "A" & Gene2Base == "-")) #subsetting out all the relevant data to chart Transition and PAM later
GenedashA = subset(DF, (Gene1Base == "-" & Gene2Base == "A"))
GeneAC = subset(DF, (Gene1Base == "A" & Gene2Base == "C"))
GeneCA = subset(DF, (Gene1Base == "C" & Gene2Base == "A"))
GeneAG = subset(DF, (Gene1Base == "A" & Gene2Base == "G"))
GeneGA = subset(DF, (Gene1Base == "G" & Gene2Base == "A"))
GeneAT = subset(DF, (Gene1Base == "A" & Gene2Base == "T"))
GeneTA = subset(DF, (Gene1Base == "T" & Gene2Base == "A"))
GeneTC = subset(DF, (Gene1Base == "T" & Gene2Base == "C"))
GeneCT = subset(DF, (Gene1Base == "C" & Gene2Base == "T"))
GeneTG = subset(DF, (Gene1Base == "T" & Gene2Base == "G"))
GeneGT = subset(DF, (Gene1Base == "G" & Gene2Base == "T"))
GeneTdash =subset(DF, (Gene1Base == "T" & Gene2Base == "-"))
GenedashT = subset(DF, (Gene1Base == "-" & Gene2Base == "T"))
GeneCG = subset(DF, (Gene1Base == "C" & Gene2Base == "G"))
GeneGC= subset(DF, (Gene1Base == "G" & Gene2Base == "C"))
GeneCdash = subset(DF, (Gene1Base == "C" & Gene2Base == "-"))
GenedashC = subset(DF, (Gene1Base == "-" & Gene2Base == "C"))
GeneGdash = subset(DF, (Gene1Base == "G" & Gene2Base == "-"))
GenedashG = subset(DF, (Gene1Base == "-" & Gene2Base == "G"))

Gene_vector_of_lengths = c(nrow(GeneAdash),nrow(GenedashA),nrow(GeneAC),nrow(GeneCA),nrow(GeneAG),nrow(GeneGA),nrow(GeneAT),nrow(GeneTA),nrow(GeneTC),nrow(GeneCT),nrow(GeneTG),nrow(GeneGT),nrow(GeneTdash),nrow(GenedashT),nrow(GeneCG),nrow(GeneGC),nrow(GeneCdash),nrow(GenedashC),nrow(GeneGdash),nrow(GenedashG))
Gene_vector_of_names = c("Adash","dashA","AC","CA","AG","GA","AT","TA","TC","CT","TG","GT","Tdash","dashT","CG","GC","Cdash","dashC","Gdash","dashG")
barplot(Gene_vector_of_lengths,names.arg = Gene_vector_of_names, main=firstM) #constructing barplot of mutation frequencies
}

counts <- data.frame(table(PlotData$Gene1Base, PlotData$Gene2Base))
gene_vector_of_names <- paste(counts$Var1, counts$Var2)
#write.csv(Gene_vector_of_lengths, file = "CBSalGap2Data.csv")

pdf('CBsalgapplot.pdf')

#-------Lattice and Chi-------------------------------------------
lDf <- data.frame(paste(substring(DF$Gene1Base,1,1),substring(DF$Gene2Base,1,1)),paste(substring(DF$Gene1Base,2),substring(DF$Gene2Base,2)),DF$PreCount,DF$PostCount)
names(lDf) <- c('Mutation1', 'Mutation2', 'PreCount', 'PostCount')
ctbl <- table(lDf$Mutation1, lDf$Mutation2)
lcounts <- data.frame(ctbl)
gap_vector_of_names <- paste(lcounts$Var1, lcounts$Var2)
barchart(lcounts$Freq~lcounts$Var2|lcounts$Var1,ylab="Mutation Frequencies",xlab="First Mutation",main="Second Mutation By First",layout=c(2,10))
print('Chi Test')
chisq.test(ctbl)
chisq.test(ctbl, simulate.p.value = TRUE, B=10000)
chisq.test(ctbl, simulate.p.value = TRUE, B=10000)
#-----------------------------------------------------------------

# Plot for each beginning mutation
# TC
graphDF <- subset(PlotData, substring(Gene1Base,1,1) == 'T' & substring(Gene2Base,1,1) == 'C')
if (nrow(graphDF) != 0) {
  make_graphs(graphDF, 'TC')
  #old way
  #graphCounts <- data.frame(table(substring(graphDF$Gene1Base,2), substring(graphDF$Gene2Base,2)))
  #graph_vector_of_names <- paste(graphCounts$Var1, graphCounts$Var2)
  #barplot(graphCounts$Freq,names.arg = graph_vector_of_names,main="TC")
}
#TG
graphDF <- subset(PlotData, substring(Gene1Base,1,1) == 'T' & substring(Gene2Base,1,1) == 'G')
if (nrow(graphDF) != 0) {
  make_graphs(graphDF, 'TG')
}
#TA
graphDF <- subset(PlotData, substring(Gene1Base,1,1) == 'T' & substring(Gene2Base,1,1) == 'A')
if (nrow(graphDF) != 0) {
  make_graphs(graphDF, 'TA')
}
#T-
graphDF <- subset(PlotData, substring(Gene1Base,1,1) == 'T' & substring(Gene2Base,1,1) == '-')
if (nrow(graphDF) != 0) {
  make_graphs(graphDF, 'T-')
}
#-C
graphDF <- subset(PlotData, substring(Gene1Base,1,1) == '-' & substring(Gene2Base,1,1) == 'C')
if (nrow(graphDF) != 0) {
  make_graphs(graphDF, '-C')
}
#-G
graphDF <- subset(PlotData, substring(Gene1Base,1,1) == '-' & substring(Gene2Base,1,1) == 'G')
if (nrow(graphDF) != 0) {
  make_graphs(graphDF, '-G')
}
#-A
graphDF <- subset(PlotData, substring(Gene1Base,1,1) == '-' & substring(Gene2Base,1,1) == 'A')
if (nrow(graphDF) != 0) {
  make_graphs(graphDF, '-A')
}
#-T
graphDF <- subset(PlotData, substring(Gene1Base,1,1) == '-' & substring(Gene2Base,1,1) == 'T')
if (nrow(graphDF) != 0) {
  make_graphs(graphDF, '-T')
}
#AC
graphDF <- subset(PlotData, substring(Gene1Base,1,1) == 'A' & substring(Gene2Base,1,1) == 'C')
if (nrow(graphDF) != 0) {
  make_graphs(graphDF, 'AC')
}
#AG
graphDF <- subset(PlotData, substring(Gene1Base,1,1) == 'A' & substring(Gene2Base,1,1) == 'G')
if (nrow(graphDF) != 0) {
  make_graphs(graphDF, 'AG')
}
#AT
graphDF <- subset(PlotData, substring(Gene1Base,1,1) == 'A' & substring(Gene2Base,1,1) == 'T')
if (nrow(graphDF) != 0) {
  make_graphs(graphDF, 'AT')
}
#A-
graphDF <- subset(PlotData, substring(Gene1Base,1,1) == 'A' & substring(Gene2Base,1,1) == '-')
if (nrow(graphDF) != 0) {
  make_graphs(graphDF, 'A-')
}
#GC
graphDF <- subset(PlotData, substring(Gene1Base,1,1) == 'G' & substring(Gene2Base,1,1) == 'C')
if (nrow(graphDF) != 0) {
  make_graphs(graphDF, 'GC')
}
#GA
graphDF <- subset(PlotData, substring(Gene1Base,1,1) == 'G' & substring(Gene2Base,1,1) == 'A')
if (nrow(graphDF) != 0) {
  make_graphs(graphDF, 'GA')
}
#GT
graphDF <- subset(PlotData, substring(Gene1Base,1,1) == 'G' & substring(Gene2Base,1,1) == 'T')
if (nrow(graphDF) != 0) {
  make_graphs(graphDF, 'GT')
}
#G-
graphDF <- subset(PlotData, substring(Gene1Base,1,1) == 'G' & substring(Gene2Base,1,1) == '-')
if (nrow(graphDF) != 0) {
  make_graphs(graphDF, 'G-')
}
#CG
graphDF <- subset(PlotData, substring(Gene1Base,1,1) == 'C' & substring(Gene2Base,1,1) == 'G')
if (nrow(graphDF) != 0) {
  make_graphs(graphDF, 'CG')
}
#CA
graphDF <- subset(PlotData, substring(Gene1Base,1,1) == 'C' & substring(Gene2Base,1,1) == 'A')
if (nrow(graphDF) != 0) {
  make_graphs(graphDF, 'CA')
}
#CT
graphDF <- subset(PlotData, substring(Gene1Base,1,1) == 'C' & substring(Gene2Base,1,1) == 'T')
if (nrow(graphDF) != 0) {
  make_graphs(graphDF, 'CT')
}
#C-
graphDF <- subset(PlotData, substring(Gene1Base,1,1) == 'C' & substring(Gene2Base,1,1) == '-')
if (nrow(graphDF) != 0) {
  make_graphs(graphDF, 'C-')
}

barplot(counts$Freq, names.arg = gene_vector_of_names)

system("rm /home/nmarkle/FastaInGCS*")

dev.off()

