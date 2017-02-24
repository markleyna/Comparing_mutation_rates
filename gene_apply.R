# Additional package
install.packages("parallel")
library("parallel")

#Variable Declaration Block
Genes <- read.csv("RVersionCBECGenes.csv", stringsAsFactors = FALSE)
CitroBacDNA<-read_file("CitroBacKPureDNA.txt")
CitroBacDNA<-gsub("\n","",CitroBacDNA)
EColiDNA <- read_file("EColiPureDNA.txt")
EColiDNA<-gsub("\n","",EColiDNA)
Genes <- data.frame(Genes$gapnum, Genes$cbgeneseq, Genes$ecgeneseq, stringsAsFactors = FALSE)

G1<-as.vector(Genes$Genes.gapnum)
G2<-data.frame(Genes$Genes.cbgeneseq)
G3<-data.frame(Genes$Genes.ecgeneseq)
n=1
buildchar = "Q" #done to prevent things from crashing later on, don't ask
Gene1Base <- c(buildchar,buildchar)
Gene2Base<-c(buildchar,buildchar) 
PreCount<-c(-1,-1)
PostCount<-c(-1,-1)
Gene1break = ""
Gene2break = ""
DF = data.frame(Gene1Base,Gene2Base,PreCount,PostCount)
DF <- data.frame(lapply(DF, as.character), stringsAsFactors=FALSE) #done to stop shenanigans with strings later on
Gene1RunTitle <- paste(">","ECGap_number","k", sep = "") #must be 2 digit representation for code to work
Gene2RunTitle <- paste(">","SalGap_number","k", sep = "") #must be 3 digit representation for code to work

clusal_run <- function(testNum, Gene1, Gene2) {
  Gene1Test = toString(Gene1[1])
  Gene2Test = toString(Gene2[1])
  
  Gene1RunTitle <- paste(">","ECGap_number","k", sep = "") #I don't feel like rewriting this to hold generality, so these need to stay as is to keep the parsing working correctly
  Gene2RunTitle <- paste(">","SalGap_number","k", sep = "")
  
  fileName <- paste("FastaIn", testNum, sep = "") # want the files to be unique when going parallel
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
  
  #pull the gapped strings out of the .aln file
  
  # Will either steal from Brian or help him look at it more
}

# need to test to make sure bread and butter mapply works first but here's the parallel version:
# mcmapply(clusal_run, G1, G2, G3, SIMPLIFY = FALSE, mc.cores=parallel::detectCores()-1)
# need to reevaluate the number of cores to use--probably too many

DF <- mapply(clusal_run, G1, G2, G3, SIMPLIFY = FALSE)
finalDF <- do.call(rbind, DF)

#CLUSAL Run on Genes
#while (n <= nrow(Genes)){
while(n<=20){ #testing line when not running full version of code
  testnum<-G1[n]
  Gene1Test <- toString(G2[n,1])
  Gene2Test<-toString(G3[n,1])
  Gene1RunTitle <- paste(">","ECGap_number","k", sep = "") #I don't feel like rewriting this to hold generality, so these need to stay as is to keep the parsing working correctly
  Gene2RunTitle <- paste(">","SalGap_number","k", sep = "")
  
  #write a text file that is in proper FASTA format
  fileConn<-file("FastaIn.txt") #reset the file to being blank text again
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
  #Call the clustalw function 
  
  #pull the gapped strings out of the .aln file
  
  p= 0
  flag = FALSE
  tempstring=""
  Gene1Bit <- vector(mode="character", length=10)
  counter = 0
  while (p<=nchar(pagecode))
  {
    if (substring(pagecode,p,p)=='E'){
      if(substring(pagecode,p,p+4)=='ECGap'){
        flag = TRUE
      }
    }
    if(flag){
      if(substring(pagecode,p,p)=='S')
      {
        flag = FALSE
        Gene1Bit[counter] = tempstring
        tempstring = ""
        counter = counter +1 
      }
      if(substring(pagecode,p,p)!='S')
      {
        tempstring = paste(tempstring, substring(pagecode,p,p),sep="")
      }
      
    }
    p= p+1
  }
  Gene1gapstring = str_c(substring(Gene1Bit, 1+nchar(Gene1RunTitle)),collapse = "")
  Gene1gapstring = substring(Gene1gapstring,0,nchar(Gene1gapstring)-9)
  nchar(Gene1gapstring)
  ##replacingblock
  p= 0
  flag = FALSE
  tempstring=""
  Gene2Bit <- vector(mode="character", length=10)
  counter = 0
  while (p<=nchar(pagecode))
  {
    if (substring(pagecode,p,p)=='S'){
      if(substring(pagecode,p,p+5)=='SalGap'){
        flag = TRUE
      }
    }
    if(flag){
      if(substring(pagecode,p,p)=='*'||substring(pagecode,p,p)=='<'||substring(pagecode,p,p)=='E')
      {
        flag = FALSE
        Gene2Bit[counter] = tempstring
        tempstring = ""
        counter = counter +1 
      }
      if(substring(pagecode,p,p)!='*')
      {
        tempstring = paste(tempstring, substring(pagecode,p,p),sep="")
      }
      
    }
    p= p+1
  }
  Gene2gapstring = str_c(substring(Gene2Bit, 1+nchar(Gene2RunTitle)),collapse = "")
  Gene2gapstring = substring(Gene2gapstring,0,nchar(Gene2gapstring)-10)
  nchar(Gene2gapstring)
  
  r=1
  matchstring = ""
  while (r<=nchar(Gene1gapstring))
  {
    if (substring(Gene1gapstring,r,r) == substring(Gene2gapstring,r,r))
    {
      matchstring = paste(matchstring, 'T', sep="")
    }
    if(substring(Gene1gapstring,r,r) != substring(Gene2gapstring,r,r))
    {
      matchstring = paste(matchstring, 'F', sep="")
    }
    r=r+1
  }
  
  n=n+1
  Sys.sleep(3) # done to preven getting kicked off the CLUSTALw website, probably can bring this down that I'm doing more analysis, but w/e
  s = 2
  Flag1 = FALSE
  Flag2 = FALSE
  PreCounter = 0
  PostCounter =0
  #start of intra-gene break isolation
  while (s <= nchar(Gene1gapstring)){
    
    Gene1workingcharacter = substring(Gene1gapstring,s,s)
    Gene2workingcharacter = substring(Gene2gapstring,s,s)
    s = s+1
    if(Gene1workingcharacter != Gene2workingcharacter && Flag2 == TRUE)
    {
      Addrow = c(Gene1break,Gene2break,PreCounter,PostCounter)
      DF = rbind(DF,Addrow)
      print("ASDF")
      Flag2=FALSE
      PostCounter=0
      PreCounter=0
    }
    else if(Gene1workingcharacter == Gene2workingcharacter && Flag2==TRUE)
    {
      PostCounter = PostCounter +1
      #print(PostCounter) causes too much lag
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
      Flag1 = TRUE
    }
    
    
    
  }
  #end of intra-gene break identification
}
