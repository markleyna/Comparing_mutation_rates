#Import Packages
#install.packages("stringr")
#install.packages('readr',"~/Rlibs","https://cran.cnr.berkeley.edu/")
library(stringr)
library('readr', lib.loc = "/home/nmarkle/Rlibs/")
#library('readr')
# Additional package
#install.packages("parallel", "~/Rlibs", "https://cran.cnr.berkeley.edu/")
library("parallel")

#Variable Declaration Block
# Where you enter the files you'll need
Genes <- read.csv("/home/nmarkle/Comparing_mutation_rates/RVersionCBECGenes.csv", stringsAsFactors = FALSE)
CitroBacDNA<-read_file("/home/nmarkle/Comparing_mutation_rates/CitroBacKPureDNA.txt")
CitroBacDNA<-gsub("\n","",CitroBacDNA)
EColiDNA <- read_file("/home/nmarkle/Comparing_mutation_rates/EColiPureDNA.txt")
EColiDNA<-gsub("\n","",EColiDNA)
Genes <- data.frame(Genes$gapnum, Genes$cbgeneseq, Genes$ecgeneseq, stringsAsFactors = FALSE)
#-------Testing line!--------
#Genes <- head(Genes, n=20)
#----------------------------
G1<-as.vector(Genes$Genes.gapnum)
G2<-as.vector(Genes$Genes.cbgeneseq)
G3<-as.vector(Genes$Genes.ecgeneseq)
buildchar = "Q" #done to prevent things from crashing later on, don't ask
Gene1Base <- c(buildchar,buildchar)
Gene2Base<-c(buildchar,buildchar) 
PreCount<-c(-1,-1)
PostCount<-c(-1,-1)
DF = data.frame(Gene1Base,Gene2Base,PreCount,PostCount)
DF <- data.frame(lapply(DF, as.character), stringsAsFactors=FALSE) #done to stop shenanigans with strings later on
Gene1RunTitle <- paste(">","ECGap_number","k", sep = "") #must be 2 digit representation for code to work
Gene2RunTitle <- paste(">","SalGap_number","k", sep = "") #must be 3 digit representation for code to work

clusal_run <- function(testNum, Gene1, Gene2) {
  # TESTING LINES
  #print(paste("testNum: ", testNum, sep=""))
  #return(data.frame(Gene1,Gene2))
  # Temp DF to be returned at the end
  Gene1Base <- c(buildchar,buildchar)
  Gene2Base<-c(buildchar,buildchar) 
  PreCount<-c(-1,-1)
  PostCount<-c(-1,-1)
  tempDF <- data.frame(Gene1Base,Gene2Base,PreCount,PostCount)
  tempDF <- data.frame(lapply(DF, as.character), stringsAsFactors=FALSE)
  Gene1Test = toString(Gene1)
  Gene2Test = toString(Gene2)
  # TESTING LINES
  #print(paste("Gene1Test:", strsplit(Gene1Test,1,5,'')[[1]], sep=""))
  #return(data.frame(Gene1Test,Gene2Test))
  Gene1break = ""
  Gene2break = ""
  
  Gene1RunTitle <- paste(">","ECGap_number","k", sep = "") #I don't feel like rewriting this to hold generality, so these need to stay as is to keep the parsing working correctly
  Gene2RunTitle <- paste(">","SalGap_number","k", sep = "")
  
  fileName <- paste("FastaIn_", testNum, ".txt", sep = "") # want the files to be unique when going parallel
  #Where to start commenting out for llc.stat.purdue.edu
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
  systemCall = paste("clustalw2 -infile=", fileName, " -type=DNA", sep = "")
  system(systemCall)
  #pull the gapped strings out of the .aln file
  alnlines <- readLines(paste("FastaIn_", testNum, ".aln", sep = ""))
  #Stop commenting and un-comment the next line for llc that is
  #alnlines <- readLines(paste("Comparing_mutation_rates/Fastas/FastaIn", testNum, ".aln", sep = ""))
  
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
  # TESTING LINE
  #return(data.frame(Gene1gapstring,Gene2gapstring))
   
  #intra-gene break isolation
  s = 1
  Flag1 = FALSE
  Flag2 = FALSE
  PreCounter = 0
  PostCounter = 0
  while (s <= nchar(Gene1gapstring)) {
    Gene1workingcharacter = substring(Gene1gapstring,s,s)
    Gene2workingcharacter = substring(Gene2gapstring,s,s)
    s = s + 1
    if (Gene1workingcharacter != Gene2workingcharacter && Flag2 == TRUE) {
      Addrow = c(Gene1break,Gene2break,PreCounter,PostCounter)
      tempDF = rbind(tempDF, Addrow)
      Flag2 = TRUE
      Flag1 = FALSE
      PreCounter = PostCounter
      PostCounter = 0
      Gene1break = Gene1workingcharacter
      Gene2break = Gene2workingcharacter
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

# here's the parallel version:
# mcmapply(clusal_run, G1, G2, G3, SIMPLIFY = FALSE, mc.cores=parallel::detectCores()-1)
# need to reevaluate the number of cores to use--probably too many

DF <- mcmapply(clusal_run, G1, G2, G3, SIMPLIFY = FALSE, mc.cores=8)
finalDF <- do.call(rbind, DF)

# Graphs
print(names(finalDF))
finalDF$PostCount <- as.integer(finalDF$PostCount)
finalDF$PreCount <- as.integer(finalDF$PreCount)
Modfinal = finalDF
Modfinal = subset(Modfinal, PostCount > 0)
PlotData = subset(Modfinal, PostCount >= 3 & PreCount >= 3)
hist(as.numeric(PlotData$PreCount))
hist(as.numeric(PlotData$PostCount))
hist(as.numeric(Modfinal$PreCount))
hist(as.numeric(Modfinal$PostCount))
plot(jitter(as.numeric(finalDF$PreCount)),jitter(as.numeric(finalDF$PostCount)))
plot(as.numeric(PlotData$PreCount),as.numeric((PlotData$PostCount)))

GeneAdash = subset(finalDF, (Gene1Base == "A" & Gene2Base == "-"))
GenedashA = subset(finalDF, (Gene1Base == "-" & Gene2Base == "A"))
GeneAC = subset(finalDF, (Gene1Base == "A" & Gene2Base == "C"))
GeneCA = subset(finalDF, (Gene1Base == "C" & Gene2Base == "A"))
GeneAG = subset(finalDF, (Gene1Base == "A" & Gene2Base == "G"))
GeneGA = subset(finalDF, (Gene1Base == "G" & Gene2Base == "A"))
GeneAT = subset(finalDF, (Gene1Base == "A" & Gene2Base == "T"))
GeneTA = subset(finalDF, (Gene1Base == "T" & Gene2Base == "A"))
GeneTC = subset(finalDF, (Gene1Base == "T" & Gene2Base == "C"))
GeneCT = subset(finalDF, (Gene1Base == "C" & Gene2Base == "T"))
GeneTG = subset(finalDF, (Gene1Base == "T" & Gene2Base == "G"))
GeneGT = subset(finalDF, (Gene1Base == "G" & Gene2Base == "T"))
GeneTdash = subset(finalDF, (Gene1Base == "T" & Gene2Base == "-"))
GenedashT = subset(finalDF, (Gene1Base == "-" & Gene2Base == "T"))
GeneCG = subset(finalDF, (Gene1Base == "C" & Gene2Base == "G"))
GeneGC = subset(finalDF, (Gene1Base == "G" & Gene2Base == "C"))
GeneCdash = subset(finalDF, (Gene1Base == "C" & Gene2Base == "-"))
GenedashC = subset(finalDF, (Gene1Base == "-" & Gene2Base == "C"))
GeneGdash = subset(finalDF, (Gene1Base == "G" & Gene2Base == "-"))
GenedashG = subset(finalDF, (Gene1Base == "-" & Gene2Base == "G"))

Gene_vector_of_Lengths = c(nrow(GeneAdash),nrow(GenedashA),nrow(GeneAC),nrow(GeneCA),nrow(GeneAG),nrow(GeneGA),nrow(GeneAT),nrow(GeneTA),nrow(GeneTC),nrow(GeneCT),nrow(GeneTG),nrow(GeneGT),nrow(GeneTdash),nrow(GenedashT),nrow(GeneCG),nrow(GeneGC),nrow(GeneCdash),nrow(GenedashC),nrow(GeneGdash),nrow(GenedashG))
Gene_vector_of_names = c("Adash","dashA","AC","CA","AG","GA","AT","TA","TC","CT","TG","GT","Tdash","dashT","CG","GC","Cdash","dashC","Gdash","dashG")

pdf('CBECgeneplotMapply.pdf')
barplot(Gene_vector_of_Lengths,names.arg = Gene_vector_of_names)
  

system("rm /home/nmarkle/FastaIn_*")

dev.off()
