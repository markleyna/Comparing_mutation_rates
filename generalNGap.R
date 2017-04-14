#Import Packages
#http://stackoverflow.com/questions/27893230/installation-of-package-file-path-had-non-zero-exit-status-in-r
install.packages("RSelenium")
install.packages("stringr")
install.packages("binman")
install.packages("/tmp/RtmpZOJwv6/downloaded_packages/RSelenium_1.7.0.tar.gz", repos = NULL, type="source")
library("RSelenium")
library(lib.loc = "/tmp/RtmpZOJwv6/downloaded_packages")
library(stringr)
find.package("stringr")
library("XML")

#Variable Declaration Block
Gaps <- read.csv("RVersionCBECGaps.csv")
colnames(Gaps) <- c("id", "gapnum", "ecstart", "ecend", "ecgapseq", "salstart", "salend", "salgapseq")
Gaps <- data.frame(Gaps$gapnum, Gaps$ecgapseq, Gaps$salgapseq)
G1<-as.vector(Gaps$Gaps.gapnum)
G2<-data.frame(Gaps$Gaps.ecgapseq)
G3<-data.frame(Gaps$Gaps.salgapseq)
n=1

#redeclaring everything resets the variables since I didn't bother to rename when I adapted code, and this is probably was a mistake
buildchar = "Q"
ECBase <- c(buildchar,buildchar)
SalBase<-c(buildchar,buildchar)
PreCount<-c(-1,-1)
PostCount<-c(-1,-1)
ECbreak = ""
Salbreak = ""
DF = data.frame(ECBase,SalBase,PreCount,PostCount)
DF <- data.frame(lapply(DF, as.character), stringsAsFactors=FALSE)

#Selenium Server Initialization
file.path(find.package("RSelenium"))
RSelenium::checkForServer()   
myDriver <- remoteDriver(remoteServerAddr = "localhost" 
                         , port = 4006
                         , browserName = "firefox"
)
myDriver$open()
myDriver$close()
help.search("remoteDriver")
#CLUSAL Run on Genes
#while (n <= nrow(Genes)){
while(n<=20){ #testing line when not running full version of code
  myDriver$navigate("http://www.genome.jp/tools/clustalw")
  DNACheckBox <- myDriver$findElement(using = 'xpath', "//*/input[@value = 'dna']")
  DNACheckBox$clickElement()
  testnum<-G1[n]
  Gap1Test <- toString(G2[n,1])
  Gap2Test<-toString(G3[n,1])
  TextBox <- myDriver$findElement(using = 'xpath', "//textarea[@name = 'sequence']")
  TextBox$clickElement()
  Gap1RunTitle <- paste(">","ECGap_number","k", sep = "") #I don't feel like rewriting this to hold generality, so these need to stay as is to keep the parsing working correctly
  Gap2RunTitle <- paste(">","SalGap_number","k", sep = "")
  TextBox$sendKeysToElement(list(Gap1RunTitle))
  TextBox$sendKeysToElement(list("\uE007")) #line break
  TextBox$sendKeysToElement(list(Gap1Test))
  TextBox$sendKeysToElement(list("\uE007"))
  TextBox$sendKeysToElement(list(Gap2RunTitle))
  TextBox$sendKeysToElement(list("\uE007"))
  TextBox$sendKeysToElement(list(Gap2Test))
  GoButton <- myDriver$findElement(using = 'xpath', "//*/input[@value = 'Execute Multiple Alignment']")
  GoButton$clickElement()
  pagecode <- myDriver$getPageSource()
  pagecode<-as.character(pagecode)#this is basically wizardry to pull the page code and parse everything useful out of it
  pagecode<-gsub("\n","",pagecode)# this was build on lots of trial and error, don't break it please
  Gap1RunTitle <- substring(Gap1RunTitle, 2)
  Gap2RunTitle <- substring(Gap2RunTitle,2)
  pagecode <- str_replace_all(string=pagecode, pattern=" ", repl="")
  head(pagecode)
  
  p= 0
  flag = FALSE
  tempstring=""
  Gap1Bit <- vector(mode="character", length=10)
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
        Gap1Bit[counter] = tempstring
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
  Gap1gapstring = str_c(substring(Gap1Bit, 1+nchar(Gap1RunTitle)),collapse = "")
  Gap1gapstring = substring(Gap1gapstring,0,nchar(Gap1gapstring)-9)
  nchar(Gap1gapstring)
  
  q=0
  Gap2Flag = FALSE
  tempstring=""
  Gap2Bit <- vector(mode="character", length=10)
  counter = 0
  while (q<=nchar(pagecode))
  {
    if (substring(pagecode,q,q)=='S'){
      if(substring(pagecode,q,q+6)=='SalGap_'){
        Gap2Flag = TRUE
      }
    }
    if(Gap2Flag){
      if(substring(pagecode,q,q)=='*'||substring(pagecode,q,q)=='<'||substring(pagecode,q,q)=='E')
      {
        Gap2Flag = FALSE
        Gap2Bit[counter] = tempstring
        tempstring = ""
        counter = counter + 1 
      }
      if(substring(pagecode,q,q)!='*')
      {
        tempstring = paste(tempstring, substring(pagecode,q,q),sep="")
      }
      
    }
    q= q+1
  }
  Gap2gapstring = str_c(substring(Gap2Bit, 1+nchar(Gap2RunTitle)),collapse = "")
  Gap2gapstring = substring(Gap2gapstring,0,nchar(Gap2gapstring)-10)
  Gap2gapstring = gsub("k","",Gap2gapstring)
  nchar(Gap2gapstring)
  
numInRow = 2
#findMismatches <- function(numInRow, Gap1gapstring, Gap2gapstring) {
n=n+1
Sys.sleep(3) # done to preven getting kicked off the CLUSTALw website, probably can bring this down that I'm doing more analysis, but w/e
s = 1
Flag1 = FALSE
Flag2 = FALSE
PreCounter = 0
PostCounter =0

while (s <= nchar(Gap1gapstring)){
  
  Gap1workingcharacter = substring(Gap1gapstring,s,s)
  Gap2workingcharacter = substring(Gap2gapstring,s,s)
  s = s+1
  if(Gap1workingcharacter != Gap2workingcharacter && Flag2 == TRUE)
  {
    Addrow = c(Gap1break,Gap2break,PreCounter,PostCounter)
    DF = rbind(DF,Addrow)
    Flag2=FALSE
    PostCounter=0
    PreCounter=0
  }
  else if(Gap1workingcharacter == Gap2workingcharacter && Flag2==TRUE)
  {
    PostCounter = PostCounter +1
    #print(PostCounter) causes too much lag
  }
  else if(Gap1workingcharacter != Gap2workingcharacter && Flag1==TRUE)
  {
    Flag3 = TRUE
    for (i in 1:numInRow) {
     if((i+s > nchar(Gap1gapstring)) || (substring(Gap1gapstring,s+i,s+i) == substring(Gap2gapstring,s+i,s+i))) {
       Flag3 = FALSE
       break
     }
    }
    if (Flag3 == TRUE) {
     Flag1 = FALSE
     Flag2 = TRUE
     Gap1break = substring(Gap1gapstring, s, s+numInRow+1)
     Gap2break = substring(Gap2gapstring, s, s+numInRow+1)
   }
  }
  else if(Gap1workingcharacter == Gap2workingcharacter && Flag1==TRUE)
  {
    PreCounter = PreCounter+1
  }
  else if(Gap1workingcharacter == Gap2workingcharacter && Flag1==FALSE)
  {
    Flag1 = TRUE
  }
  
}
}

# now looking at results for n mismatches--but isn't depedent on n.
DF$PostCount <- as.integer(DF$PostCount)
DF$PreCount <- as.integer((DF$PreCount))
ModDF = DF
ModDF = subset(ModDF, PostCount > 0)
PlotData = subset(ModDF, PostCount>=3 & PreCount>=3)

#creates a data frame where it finds the count for each unique combination of ECBase and SalBase
counts <- data.frame(table(PlotData$ECBase, PlotData$SalBase)) 
#for future reference data.frame(table(PlotData$ECBase, PlotData$SalBase)[,]) is another way
#for this method each row represents a unique value for ECBase
#each of the columns then is the count for each unique value of SalBase

gap_vector_of_names <- paste(counts$ECBase, counts$SalBase)

pdf("SalCBGap Plots")
hist(as.numeric(PlotData$PreCount))
hist(as.numeric(PlotData$PostCount))
hist(as.numeric(ModDF$PreCount))
hist(as.numeric(ModDF$PostCount))
plot(jitter(as.numeric(DF$PreCount)),jitter(as.numeric(DF$PostCount)))
plot(as.numeric(PlotData$PreCount),as.numeric((PlotData$PostCount)))
barplot(counts$count,names.arg = gap_vector_of_names) #constructing barplot of mutation frequencies
#^-- don't actually know what the name of this column will be
dev.off()
