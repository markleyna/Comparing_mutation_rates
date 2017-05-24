#Import Packages
#http://stackoverflow.com/questions/27893230/installation-of-package-file-path-had-non-zero-exit-status-in-r
#install.packages("stringr", "~/Rlibs", "https://cran.cnr.berkeley.edu/")
library(stringr)
#library('readr')
library('readr', lib.loc="/home/nmarkle/Rlibs/")

#Variable Declaration Block
numInRow <- 2
Genes <- read.csv("/home/nmarkle/Comparing_mutation_rates/RVersionCBECGenes.csv", stringsAsFactors = FALSE)
CitroBacDNA<-read_file("/home/nmarkle/Comparing_mutation_rates/CitroBacKPureDNA.txt")
CitroBacDNA<-gsub("\n","",CitroBacDNA)
EColiDNA<-read_file("/home/nmarkle/Comparing_mutation_rates/EColiPureDNA.txt")
EColiDNA<-gsub("\n","",EColiDNA)
Genes <- data.frame(Genes$gapnum, Genes$cbgeneseq, Genes$ecgeneseq, stringsAsFactors = FALSE)
G1<-as.vector(Genes$Genes.gapnum)
G2<-data.frame(Genes$Genes.cbgeneseq)
G3<-data.frame(Genes$Genes.ecgeneseq)

#redeclaring everything resets the variables since I didn't bother to rename when I adapted code, and this is probably was a mistake
buildchar = "Q"
Gene1Base <- c(buildchar,buildchar)
Gene2Base<-c(buildchar,buildchar)
PreCount<-c(-1,-1)
PostCount<-c(-1,-1)
Gene1break = ""
Gene2break = ""
DF = data.frame(Gene1Base,Gene2Base,PreCount,PostCount)
DF <- data.frame(lapply(DF, as.character), stringsAsFactors=FALSE)
Gene1RunTitle <- paste(">", "ECGap_number","k",sep = "")
Gene2RunTitle <- paste(">", "SalGap_number", "k", sep = "")
n = 1
#CLUSAL Run on Genes
while (n <= nrow(Genes)){
#while(n<=1){ #testing line when not running full version of code
  testnum<-G1[n]
  Gene1Test <- toString(G2[n,1])
  Gene2Test <- toString(G3[n,1])
  
  #write text file
  fileConn<-file("FastaIn.txt")
  writeLines("\n",fileConn)
  close(fileConn)
  sink("FastaIn.txt")
  cat(Gene1RunTitle)
  cat("\n")
  cat(Gene1Test)
  cat("\n")
  cat(Gene2RunTitle)
  cat("\n")
  cat(Gene2Test)
  sink()
  system("clustalw2 -infile=FastaIn.txt -type=DNA")
  alnlines<-readLines("/home/nmarkle/FastaIn.aln")
  #Testing line for llc server
  #alnlines<-readLines("/home/nmarkle/Comparing_mutation_rates/FastaTest.aln")
  
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
   
  n=n+1
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
        DF = rbind(DF, Addrow)
      }
      Flag2 = TRUE
      Flag1 = FALSE
      PreCounter = PostCounter
      PostCounter = 0
      Gene1break = Gene1workingcharacter
      Gene2break = Gene2workingcharacter
    }
    else if (Gene1workingcharacter != Gene2workingcharacter && Flag2 == TRUE && PostCounter == 0) {
      #Gene1break <- c(Gene1break, Gene1workingcharacter)
      Gene1break <- paste(Gene1break, Gene1workingcharacter, sep = '')
      #Gene2break <- c(Gene2break, Gene2workingcharacter)
      Gene2break <- paste(Gene2break, Gene2workingcharacter,sep = '')
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
}

# now looking at results for n mismatches--but isn't depedent on n.
DF$PostCount <- as.integer(DF$PostCount)
DF$PreCount <- as.integer((DF$PreCount))
ModDF = DF
ModDF = subset(ModDF, PostCount > 0)
PlotData = subset(ModDF, PostCount>=3 & PreCount>=3)
#PlotData = ModDF

#creates a data frame where it finds the count for each unique combination of ECBase and SalBase
counts <- data.frame(table(PlotData$Gene1Base, PlotData$Gene2Base)) 
#for future reference data.frame(table(PlotData$ECBase, PlotData$SalBase)[,]) is another way
#for this method each row represents a unique value for ECBase
#each of the columns then is the count for each unique value of SalBase

gap_vector_of_names <- paste(counts$Var1, counts$Var2)

pdf("SalCBGap Plots")
hist(as.numeric(PlotData$PreCount))
hist(as.numeric(PlotData$PostCount))
hist(as.numeric(ModDF$PreCount))
hist(as.numeric(ModDF$PostCount))
plot(jitter(as.numeric(DF$PreCount)),jitter(as.numeric(DF$PostCount)))
plot(as.numeric(PlotData$PreCount),as.numeric((PlotData$PostCount)))
barplot(counts$Freq,names.arg = gap_vector_of_names) #constructing barplot of mutation frequencies
#^-- don't actually know what the name of this column will be
dev.off()
