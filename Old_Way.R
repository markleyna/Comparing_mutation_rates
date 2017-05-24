while (s <= nchar(Gene1gapstring)){
  
  Gene1workingcharacter = substring(Gene1gapstring,s,s)
  Gene2workingcharacter = substring(Gene2gapstring,s,s)
  
  if(Gene1workingcharacter != Gene2workingcharacter && Flag2 == TRUE)
  {
    Addrow = c(Gene1break,Gene2break,PreCounter,PostCounter)
    DF = rbind(DF,Addrow)
    Flag2=FALSE
    PostCounter=0
    PreCounter=0
  }
  else if(Gene1workingcharacter == Gene2workingcharacter && Flag2==TRUE)
  {
    PostCounter = PostCounter +1
  }
  else if(Gene1workingcharacter != Gene2workingcharacter && Flag1==TRUE)
  {
    Flag3 = TRUE
    for (i in 1:(numInRow-1)) {
      if((i+s > nchar(Gene1gapstring)) || (substring(Gene1gapstring,s+i,s+i) == substring(Gene2gapstring,s+i,s+i))) {
        s = s + i - 1
        Flag3 = FALSE
        break
      }
    }
    if (Flag3 == TRUE) {
      Flag1 = FALSE
      Flag2 = TRUE
      Gene1break = substring(Gene1gapstring, s, s+numInRow-1)
      Gene2break = substring(Gene2gapstring, s, s+numInRow-1)
      s = s + numInRow - 1
    }
  }
  else if(Gene1workingcharacter == Gene2workingcharacter && Flag1==TRUE)
  {
    PreCounter = PreCounter+1
  }
  else if(Gene1workingcharacter == Gene2workingcharacter && Flag1==FALSE)
  {
    Flag1 = TRUE
  }
  s = s+1
}