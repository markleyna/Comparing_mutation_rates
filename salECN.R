#Import Packages
#http://stackoverflow.com/questions/27893230/installation-of-package-file-path-had-non-zero-exit-status-in-r
#install.packages("stringr", "~/Rlibs", "https://cran.cnr.berkeley.edu/")
library(stringr)
library(parallel)
#library('readr')
library('readr', lib.loc="/home/nmarkle/Rlibs/")

#Variable Declaration Block
numInRow <- 2
Genes <- read.csv("/home/nmarkle/Comparing_mutation_rates/RVersionSalECGenes.csv", stringsAsFactors = FALSE)
SalDNA<-read_file("/home/nmarkle/Comparing_mutation_rates/SalmonellaPureDNA.txt")
SalDNA<-gsub("\n","",SalDNA)
EColiDNA<-read_file("/home/nmarkle/Comparing_mutation_rates/EColiPureDNA.txt")
EColiDNA<-gsub("\n","",EColiDNA)
Genes <- data.frame(Genes$gapnum, Genes$ecgeneseq, Genes$salgeneseq, stringsAsFactors = FALSE)
#-------Testing line!-------
#Genes <- head(Genes, n=1)
#---------------------------
G1<-as.vector(Genes$Genes.gapnum)
G2<-as.vector(Genes$Genes.ecgeneseq)
G3<-as.vector(Genes$Genes.salgeneseq)

buildchar = "Q"
Gene1Base <- c(buildchar,buildchar)
Gene2Base<-c(buildchar,buildchar)
PreCount<-c(-1,-1)
PostCount<-c(-1,-1)
DF = data.frame(Gene1Base,Gene2Base,PreCount,PostCount)
DF <- data.frame(lapply(DF, as.character), stringsAsFactors=FALSE)
Gene1RunTitle <- paste(">", "ECGap_number","k",sep = "")
Gene2RunTitle <- paste(">", "SalGap_number", "k", sep = "")

#CLUSAL Run on Genes
clusal_run <- function(testNum, Gene1, Gene2) {
  Gene1Base <- c(buildchar,buildchar)
  Gene2Base<-c(buildchar,buildchar)
  PreCount<-c(-1,-1)
  PostCount<-c(-1,-1)
  tempDF = data.frame(Gene1Base,Gene2Base,PreCount,PostCount)
  tempDF <- data.frame(lapply(DF, as.character), stringsAsFactors=FALSE)
  Gene1Test = toString(Gene1)
  Gene2Test = toString(Gene2)
  Gene1break = ""
  Gene2break = ""
  Gene1RunTitle <- paste(">", "ECGap_number","k",sep = "")
  Gene2RunTitle <- paste(">", "SalGap_number", "k", sep = "")
 
  #write text file
  fileName <- paste("FastaInSE", testNum, ".txt", sep = "")
  fileConn<-file(fileName)
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
  systemCall = paste("clustalw2 -infile=", fileName, " -type=DNA", sep="")
  system(systemCall)
  alnlines<-readLines(paste("FastaInSE", testNum, ".aln", sep= ""))
  #Testing line for llc server
  #alnlines<-readLines("/home/nmarkle/Comparing_mutation_rates/FastaTest.aln")
  #alnlines<-readLines(paste("/home/nmarkle/Comparing_mutation_rates/Fastas/FastaIn", testNum, ".aln", sep=""))

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
  while (s <= nchar(Gene1gapstring)) {
    Gene1workingcharacter = substring(Gene1gapstring,s,s)
    Gene2workingcharacter = substring(Gene2gapstring,s,s)
    s = s + 1
    if (Gene1workingcharacter != Gene2workingcharacter && Flag2 == TRUE && PostCounter > 0) {
      #print(Gene1break)
      #print(Gene2break)
      #print(nchar(Gene1break))
      if (nchar(Gene1break) == numInRow) {
        Addrow = c(Gene1break,Gene2break,PreCounter,PostCounter)
        tempDF = rbind(tempDF, Addrow)
      }
      Flag2 = TRUE
      Flag1 = FALSE
      PreCounter = PostCounter
      PostCounter = 0
      Gene1break = Gene1workingcharacter
      Gene2break = Gene2workingcharacter
    }
    else if (Gene1workingcharacter != Gene2workingcharacter && Flag2 == TRUE && PostCounter == 0) {
      if (nchar(Gene1break) < numInRow) {
        Gene1break <- paste(Gene1break, Gene1workingcharacter, sep = '')
        Gene2break <- paste(Gene2break, Gene2workingcharacter,sep = '')
      }
      else {
	      PreCounter = 0
	      Flag2 = FALSE
	      Flag1 = FALSE
	      Gene1break = ""
	      Gene2break = ""
      }
    }
    else if (Gene1workingcharacter == Gene2workingcharacter && Flag2 == TRUE) {
      PostCounter = PostCounter + 1
    }
    else if (Gene1workingcharacter != Gene2workingcharacter && Flag1 == TRUE) {
      Flag1 = FALSE
      Flag2 = TRUE
      Gene1break = Gene1workingcharacter
      Gene2break = Gene2workingcharacter
    }
    else if (Gene1workingcharacter == Gene2workingcharacter && Flag1 == TRUE) {
      PreCounter = PreCounter + 1
    }
    else if (Gene1workingcharacter == Gene2workingcharacter && Flag1 == FALSE) {
      PreCounter = PreCounter + 1
      Flag1 = TRUE
    }
    
    # end of intra-gene break identification
  }
  return(tempDF)
}

#DF <- mapply(clusal_run, G1, G2, G3, SIMPLIFY = FALSE)
DF <- mcmapply(clusal_run, G1, G2, G3, SIMPLIFY = FALSE, mc.cores = 8)
DF <- do.call(rbind,DF)

# now looking at results for n mismatches--but isn't depedent on n.
DF$PostCount <- as.integer(DF$PostCount)
DF$PreCount <- as.integer((DF$PreCount))
ModDF = DF
ModDF = subset(ModDF, PostCount > 0)
PlotData = subset(ModDF, PostCount>=3 & PreCount>=3)
#PlotData = ModDF

make_graphs <- function(finalDF, firstM) {
  GeneAdash = subset(finalDF, (substring(graphDF$Gene1Base,2) == "A" & substring(graphDF$Gene2Base,2) == "-"))
  GenedashA = subset(finalDF, (substring(graphDF$Gene1Base,2) == "-" & substring(graphDF$Gene2Base,2) == "A"))
  GeneAC = subset(finalDF, (substring(graphDF$Gene1Base,2) == "A" & substring(graphDF$Gene2Base,2) == "C"))
  GeneCA = subset(finalDF, (substring(graphDF$Gene1Base,2) == "C" & substring(graphDF$Gene2Base,2) == "A"))
  GeneAG = subset(finalDF, (substring(graphDF$Gene1Base,2) == "A" & substring(graphDF$Gene2Base,2) == "G"))
  GeneGA = subset(finalDF, (substring(graphDF$Gene1Base,2) == "G" & substring(graphDF$Gene2Base,2) == "A"))
  GeneAT = subset(finalDF, (substring(graphDF$Gene1Base,2) == "A" & substring(graphDF$Gene2Base,2) == "T"))
  GeneTA = subset(finalDF, (substring(graphDF$Gene1Base,2) == "T" & substring(graphDF$Gene2Base,2) == "A"))
  GeneTC = subset(finalDF, (substring(graphDF$Gene1Base,2) == "T" & substring(graphDF$Gene2Base,2) == "C"))
  GeneCT = subset(finalDF, (substring(graphDF$Gene1Base,2) == "C" & substring(graphDF$Gene2Base,2) == "T"))
  GeneTG = subset(finalDF, (substring(graphDF$Gene1Base,2) == "T" & substring(graphDF$Gene2Base,2) == "G"))
  GeneGT = subset(finalDF, (substring(graphDF$Gene1Base,2) == "G" & substring(graphDF$Gene2Base,2) == "T"))
  GeneTdash = subset(finalDF, (substring(graphDF$Gene1Base,2) == "T" & substring(graphDF$Gene2Base,2) == "-"))
  GenedashT = subset(finalDF, (substring(graphDF$Gene1Base,2) == "-" & substring(graphDF$Gene2Base,2) == "T"))
  GeneCG = subset(finalDF, (substring(graphDF$Gene1Base,2) == "C" & substring(graphDF$Gene2Base,2) == "G"))
  GeneGC = subset(finalDF, (substring(graphDF$Gene1Base,2) == "G" & substring(graphDF$Gene2Base,2) == "C"))
  GeneCdash = subset(finalDF, (substring(graphDF$Gene1Base,2) == "C" & substring(graphDF$Gene2Base,2) == "-"))
  GenedashC = subset(finalDF, (substring(graphDF$Gene1Base,2) == "-" & substring(graphDF$Gene2Base,2) == "C"))
  GeneGdash = subset(finalDF, (substring(graphDF$Gene1Base,2) == "G" & substring(graphDF$Gene2Base,2) == "-"))
  GenedashG = subset(finalDF, (substring(graphDF$Gene1Base,2) == "-" & substring(graphDF$Gene2Base,2) == "G"))
  
  Gene_vector_of_Lengths = c(nrow(GeneAdash),nrow(GenedashA),nrow(GeneAC),nrow(GeneCA),nrow(GeneAG),nrow(GeneGA),nrow(GeneAT),nrow(GeneTA),nrow(GeneTC),nrow(GeneCT),nrow(GeneTG),nrow(GeneGT),nrow(GeneTdash),nrow(GenedashT),nrow(GeneCG),nrow(GeneGC),nrow(GeneCdash),nrow(GenedashC),nrow(GeneGdash),nrow(GenedashG))
  Gene_vector_of_names = c("Adash","dashA","AC","CA","AG","GA","AT","TA","TC","CT","TG","GT","Tdash","dashT","CG","GC","Cdash","dashC","Gdash","dashG")
  barplot(Gene_vector_of_Lengths,names.arg = Gene_vector_of_names, main=firstM)
}

#creates a data frame where it finds the count for each unique combination of ECBase and SalBase
counts <- data.frame(table(PlotData$Gene1Base, PlotData$Gene2Base)) 
#for future reference data.frame(table(PlotData$ECBase, PlotData$SalBase)[,]) is another way
#for this method each row represents a unique value for ECBase
#each of the columns then is the count for each unique value of SalBase

gap_vector_of_names <- paste(counts$Var1, counts$Var2)

pdf("SalECGene Plots.pdf")
hist(as.numeric(PlotData$PreCount))
hist(as.numeric(PlotData$PostCount))
hist(as.numeric(ModDF$PreCount))
hist(as.numeric(ModDF$PostCount))
plot(jitter(as.numeric(DF$PreCount)),jitter(as.numeric(DF$PostCount)))
plot(as.numeric(PlotData$PreCount),as.numeric((PlotData$PostCount)))
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

barplot(counts$Freq,names.arg = gap_vector_of_names) #constructing barplot of mutation frequencies

system("rm /home/nmarkle/FastaInSE*")

dev.off()
