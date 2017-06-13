#Import Packages
print("entered the script")
setwd("/home/french15")
#install.packages("RSelenium")
#install.packages("stringr")
#install.packages("readr")
library("RSelenium")
library(stringr)
library('readr')
print("loaded everything successfully")
#Predeclaration of Functions
WhatAmino<-function(threebase){
  #converts a 3 base region of DNA to the corresponding amino acid
  if (threebase == 'TAG' ||threebase ==  'TAA' ||threebase ==  'TGA'){
    return ("Stop")
  }
  if (threebase == "ATT"||threebase == "ATC"||threebase == "ATA"){
    return ("Isoleucine")
  }
  if(threebase == "CTT"||threebase == "CTC"||threebase == "CTA"||threebase == "CTG"||threebase == "TAA"||threebase == "TTG"||threebase=="TTA"){
    return("Leucine")
  }
  if(threebase == "GTC"||threebase == "GTT"||threebase == "GTA"||threebase == "GTG")
  {
    return("Valine")
  }
  if(threebase == "TTT"||threebase == "TTC")
  {
    return("Phenylalanine")
  }
  if(threebase == "ATG"){
    return("Methionine")
  }
  if(threebase == "TGT"||threebase == "TGC"){
    return("Cysteine")
  }
  if(threebase == "GCT"||threebase == "GCC"||threebase == "GCA"||threebase=="GCG"){
    return("Alanine")
  }
  if(threebase == "GCT"||threebase == "GGC"||threebase == "GGA"||threebase == "GGG"||threebase=="GGT"){
    return("Glycine")
  }
  if(threebase == "CCT"||threebase == "CCC"||threebase == "CCA"||threebase == "CCG"){
    return("Proline")
  }
  if(threebase == "ACT"||threebase == "ACC"||threebase == "ACA"||threebase == "ACG"){
    return("Threonine")
  }
  if(threebase == "TCT"||threebase == "TCC"||threebase == "TCA"||threebase == "TCG"||threebase == "AGT"||threebase == "AGC"){
    return("Serine")
  }
  if(threebase == "TAT"||threebase == "TAC"){
    return("Tyrosine")
  }
  if(threebase == "TGG"){
    return("Tryptophan")
  }
  if(threebase == "CAA"||threebase == "CAG"){
    return("Glutamine")
  }
  if(threebase == "AAT"||threebase == "AAC"){
    return("Asparagine")
  }
  if(threebase == "CAT"||threebase == "CAC"){
    return("Histidine")
  }
  if(threebase == "GAA"||threebase == "GAC"||threebase=="GAG"){
    return("Glutamic Acid")
  }
  if(threebase == "GAT"||threebase == "GAC"){
    return("Aspartic Acid")
  }
  if(threebase == "AAA"||threebase == "AAG"){
    return("Lysine")
  }
  if(threebase == "CGT"||threebase == "CGC"||threebase == "CGA"||threebase == "CGG"||threebase == "AGA"||threebase == "AGG"){
    return("Arginine")
  }
  else{
    return(threebase)
  }
}
Transition <- function(Amino1,Amino2,MatrixMethod){
  #converts amion acid pair to the corresponding BLOSUM62 score
  #will need to rename, don't want to change that until I know where it's called from
  #MatrixMethod corresponds to either Transition (B) or PAM (P)
    #Transition is always the first return, PAM is the else statement
  if ((Amino1=="Stop")||(Amino2=="Stop")){
    
  }
  if((Amino1=="Alanine"&&Amino2=="Alanine")||(Amino2=="Alanine"&&Amino1=="Alanine")){
    if (MatrixMethod=="B") return(4)
    else return (2)
  }
  if((Amino1=="Arginine"&&Amino2=="Alanine")||(Amino2=="Arginine"&&Amino1=="Alanine")){
    if(MatrixMethod=="B") return(-1)
    else return (-2)
  }
  if((Amino1=="Asparagine"&&Amino2=="Alanine")||(Amino2=="Asparagine"&&Amino1=="Alanine")){
    if(MatrixMethod=="B") return(-2)
    else return (0)  
  }
  if((Amino1=="Aspartic Acid"&&Amino2=="Alanine")||(Amino2=="Aspartic Acid"&&Amino1=="Alanine")){
    if(MatrixMethod=="B") return(-2)
    else return (0)
  }
  if((Amino1=="Cysteine"&&Amino2=="Alanine")||(Amino2=="Cysteine"&&Amino1=="Alanine")){
    if(MatrixMethod=="B") return(0)
    else return (-2)
  }
  if((Amino1=="Glutamine"&&Amino2=="Alanine")||(Amino2=="Glutamine"&&Amino1=="Alanine")){
    if(MatrixMethod=="B") return(-1)
    else return (0)
  }
  if((Amino1=="Glutamic Acid"&&Amino2=="Alanine")||(Amino2=="Glutamic Acid"&&Amino1=="Alanine")){
    if(MatrixMethod=="B") return(-1)
    else return (0)
  }
  if((Amino1=="Glycine"&&Amino2=="Alanine")||(Amino2=="Glycine"&&Amino1=="Alanine")){
    if(MatrixMethod=="B") return(0)
    else return (1)
  }
  if((Amino1=="Histidine"&&Amino2=="Alanine")||(Amino2=="Histidine"&&Amino1=="Alanine")){
    if(MatrixMethod=="B") return(-2)
    else return (-1)
  }
  if((Amino1=="Isoleucine"&&Amino2=="Alanine")||(Amino2=="Isoleucine"&&Amino1=="Alanine")){
    if(MatrixMethod=="B") return(-1)
    else return (-1)
  }
  if((Amino1=="Leucine"&&Amino2=="Alanine")||(Amino2=="Leucine"&&Amino1=="Alanine")){
    #didn't see a return option for Leucine and Alanine in PAM method
    return(-1)
  }
  if((Amino1=="Lysine"&&Amino2=="Alanine")||(Amino2=="Lysine"&&Amino1=="Alanine")){
    if(MatrixMethod=="B") return(-1)
    else return (-2)
  }
  if((Amino1=="Methionine"&&Amino2=="Alanine")||(Amino2=="Methionine"&&Amino1=="Alanine")){
    if(MatrixMethod=="B") return(-1)
    else return (-1)
  }
  if((Amino1=="Phenylalanine"&&Amino2=="Alanine")||(Amino2=="Phenylalanine"&&Amino1=="Alanine")){
    if(MatrixMethod=="B") return(-2)
    else return (-4)
  }
  if((Amino1=="Proline"&&Amino2=="Alanine")||(Amino2=="Proline"&&Amino1=="Alanine")){
    if(MatrixMethod=="B") return(-1)
    else return (1)
  }
  if((Amino1=="Serine"&&Amino2=="Alanine")||(Amino2=="Serine"&&Amino1=="Alanine")){
    if(MatrixMethod=="B") return(1)
    else return (1)
  }
  if((Amino1=="Threonine"&&Amino2=="Alanine")||(Amino2=="Threonine"&&Amino1=="Alanine")){
    if(MatrixMethod=="B") return(0)
    else return (1)
  }
  if((Amino1=="Tryptophan"&&Amino2=="Alanine")||(Amino2=="Tryptophan"&&Amino1=="Alanine")){
    if(MatrixMethod=="B") return(-3)
    else return (-6)
  }
  if((Amino1=="Tyrosine"&&Amino2=="Alanine")||(Amino2=="Tyrosine"&&Amino1=="Alanine")){
    if(MatrixMethod=="B") return(-2)
    else return (-3)
  }
  if((Amino1=="Valine"&&Amino2=="Alanine")||(Amino2=="Valine"&&Amino1=="Alanine")){
    if(MatrixMethod=="B") return(0)
    else return (0)
  }
  #Arganine Block
  if((Amino1=="Arginine"&&Amino2=="Arginine")||(Amino2=="Arginine"&&Amino1=="Arginine")){
    if(MatrixMethod=="B") return(5)
    else return (6)
  }
  if((Amino1=="Asparagine"&&Amino2=="Arginine")||(Amino2=="Asparagine"&&Amino1=="Arginine")){
    if(MatrixMethod=="B") return(0)
    else return (0)
  }
  if((Amino1=="Aspartic Acid"&&Amino2=="Arginine")||(Amino2=="Aspartic Acid"&&Amino1=="Arginine")){
    if(MatrixMethod=="B") return(-2)
    else return (-1)
  }
  if((Amino1=="Cysteine"&&Amino2=="Arginine")||(Amino2=="Cysteine"&&Amino1=="Arginine")){
    if(MatrixMethod=="B") return(-3)
    else return (-4)
  }
  if((Amino1=="Glutamine"&&Amino2=="Arginine")||(Amino2=="Glutamine"&&Amino1=="Arginine")){
    if(MatrixMethod=="B") return(1)
    else return (1)
  }
  if((Amino1=="Glutamic Acid"&&Amino2=="Arginine")||(Amino2=="Glutamic Acid"&&Amino1=="Arginine")){
    if(MatrixMethod=="B") return(0)
    else return (-1)
  }
  if((Amino1=="Glycine"&&Amino2=="Arginine")||(Amino2=="Glycine"&&Amino1=="Arginine")){
    if(MatrixMethod=="B") return(-2)
    else return (-3)
  }
  if((Amino1=="Histidine"&&Amino2=="Arginine")||(Amino2=="Histidine"&&Amino1=="Arginine")){
    if(MatrixMethod=="B") return(0)
    else return(2)
  }
  if((Amino1=="Isoleucine"&&Amino2=="Arginine")||(Amino2=="Isoleucine"&&Amino1=="Arginine")){
    if(MatrixMethod=="B") return(-3)
    else return(-2)
  }
  if((Amino1=="Leucine"&&Amino2=="Arginine")||(Amino2=="Leucine"&&Amino1=="Arginine")){
    if(MatrixMethod=="B") return(-2)
    else return (-3)
  }
  if((Amino1=="Lysine"&&Amino2=="Arginine")||(Amino2=="Lysine"&&Amino1=="Arginine")){
    if(MatrixMethod=="B") return(2)
    else return (3)
  }
  if((Amino1=="Methionine"&&Amino2=="Arginine")||(Amino2=="Methionine"&&Amino1=="Arginine")){
    if(MatrixMethod=="B") return(-1)
    else return (0)
  }
  if((Amino1=="Phenylalanine"&&Amino2=="Arginine")||(Amino2=="Phenylalanine"&&Amino1=="Arginine")){
    if(MatrixMethod=="B") return(-3)
    else return (-4)
  }
  if((Amino1=="Proline"&&Amino2=="Arginine")||(Amino2=="Proline"&&Amino1=="Arginine")){
    if(MatrixMethod=="B") return(-2)
    else return (0)
  }
  if((Amino1=="Serine"&&Amino2=="Arginine")||(Amino2=="Serine"&&Amino1=="Arginine")){
    if(MatrixMethod=="B") return(-1)
    else return (0)
  }
  if((Amino1=="Threonine"&&Amino2=="Arginine")||(Amino2=="Threonine"&&Amino1=="Arginine")){
    if(MatrixMethod=="B") return(-1)
    else return (-1)
  }
  if((Amino1=="Tryptophan"&&Amino2=="Arginine")||(Amino2=="Tryptophan"&&Amino1=="Arginine")){
    if(MatrixMethod=="B") return(-3)
    else return (2)
  }
  if((Amino1=="Tyrosine"&&Amino2=="Arginine")||(Amino2=="Tyrosine"&&Amino1=="Arginine")){
    if(MatrixMethod=="B") return(-2)
    else return (-4)
  }
  if((Amino1=="Valine"&&Amino2=="Arginine")||(Amino2=="Valine"&&Amino1=="Arginine")){
    if(MatrixMethod=="B") return(-2)
    else return (-2)
  }
  #aspargine Block
  if((Amino1=="Asparagine"&&Amino2=="Asparagine")||(Amino2=="Asparagine"&&Amino1=="Asparagine")){
    if(MatrixMethod=="B") return(6)
    else return (2)
  }
  if((Amino1=="Aspartic Acid"&&Amino2=="Asparagine")||(Amino2=="Aspartic Acid"&&Amino1=="Asparagine")){
    if(MatrixMethod=="B") return(1)
    else return (2)
  }
  if((Amino1=="Cysteine"&&Amino2=="Asparagine")||(Amino2=="Cysteine"&&Amino1=="Asparagine")){
    if(MatrixMethod=="B") return(-3)
    else return (-4)
  }
  if((Amino1=="Glutamine"&&Amino2=="Asparagine")||(Amino2=="Glutamine"&&Amino1=="Asparagine")){
    if(MatrixMethod=="B") return(0)
    else return (1)
  }
  if((Amino1=="Glutamic Acid"&&Amino2=="Asparagine")||(Amino2=="Glutamic Acid"&&Amino1=="Asparagine")){
    if(MatrixMethod=="B") return(0)
    else return (1)
  }
  if((Amino1=="Glycine"&&Amino2=="Asparagine")||(Amino2=="Glycine"&&Amino1=="Asparagine")){
    if(MatrixMethod=="B") return(0)
    else return (0)
  }
  if((Amino1=="Histidine"&&Amino2=="Asparagine")||(Amino2=="Histidine"&&Amino1=="Asparagine")){
    if(MatrixMethod=="B") return(1)
    else return (2)
  }
  if((Amino1=="Isoleucine"&&Amino2=="Asparagine")||(Amino2=="Isoleucine"&&Amino1=="Asparagine")){
    if(MatrixMethod=="B") return(-3)
    else return (-2)
  }
  if((Amino1=="Leucine"&&Amino2=="Asparagine")||(Amino2=="Leucine"&&Amino1=="Asparagine")){
    if(MatrixMethod=="B") return(-3)
    else return (-3)
  }
  if((Amino1=="Lysine"&&Amino2=="Asparagine")||(Amino2=="Lysine"&&Amino1=="Asparagine")){
    if(MatrixMethod=="B") return(-1)
    else return (1)
  }
  if((Amino1=="Methionine"&&Amino2=="Asparagine")||(Amino2=="Methionine"&&Amino1=="Asparagine")){
    if(MatrixMethod=="B") return(-3)
    else return (-2)
  }
  if((Amino1=="Phenylalanine"&&Amino2=="Asparagine")||(Amino2=="Phenylalanine"&&Amino1=="Asparagine")){
    if(MatrixMethod=="B") return(-3)
    else return (-4)
  }
  if((Amino1=="Proline"&&Amino2=="Asparagine")||(Amino2=="Proline"&&Amino1=="Asparagine")){
    if(MatrixMethod=="B") return(-2)
    else return (-1)
  }
  if((Amino1=="Serine"&&Amino2=="Asparagine")||(Amino2=="Serine"&&Amino1=="Asparagine")){
    if(MatrixMethd=="B") return(1)
    else return (1)
  }
  if((Amino1=="Threonine"&&Amino2=="Asparagine")||(Amino2=="Threonine"&&Amino1=="Asparagine")){
    if(MatrixMethod=="B") return(0)
    else return (0)
  }
  if((Amino1=="Tryptophan"&&Amino2=="Asparagine")||(Amino2=="Tryptophan"&&Amino1=="Asparagine")){
    if(MatrixMethod=="B") return(-4)
    else return (-4)
  }
  if((Amino1=="Tyrosine"&&Amino2=="Asparagine")||(Amino2=="Tyrosine"&&Amino1=="Asparagine")){
    if(MatrixMethod=="B") return(-2)
    else return (-2)
  }
  if((Amino1=="Valine"&&Amino2=="Asparagine")||(Amino2=="Valine"&&Amino1=="Asparagine")){
    if(MatrixMethod=="B") return(-3)
    else return (-2)
  }
  #Aspartic Acid Block
  if((Amino1=="Aspartic Acid"&&Amino2=="Aspartic Acid")||(Amino2=="Aspartic Acid"&&Amino1=="Aspartic Acid")){
    if(MatrixMethod=="B") return(6)
    else return (4)
  }
  if((Amino1=="Cysteine"&&Amino2=="Aspartic Acid")||(Amino2=="Cysteine"&&Amino1=="Aspartic Acid")){
    if(MatrixMethod=="B") return(-3)
    else return (-5)
  }
  if((Amino1=="Glutamine"&&Amino2=="Aspartic Acid")||(Amino2=="Glutamine"&&Amino1=="Aspartic Acid")){
    if(MatrixMethod=="B") return(0)
    else return (2)
  }
  if((Amino1=="Glutamic Acid"&&Amino2=="Aspartic Acid")||(Amino2=="Glutamic Acid"&&Amino1=="Aspartic Acid")){
    if(MatrixMethod=="B") return(2)
    else return (3)
  }
  if((Amino1=="Glycine"&&Amino2=="Aspartic Acid")||(Amino2=="Glycine"&&Amino1=="Aspartic Acid")){
    if(MatrixMethod=="B") return(-1)
    else return (1)
  }
  if((Amino1=="Histidine"&&Amino2=="Aspartic Acid")||(Amino2=="Histidine"&&Amino1=="Aspartic Acid")){
    if(MatrixMethod=="B") return(-1)
    else return(1)
  }
  if((Amino1=="Isoleucine"&&Amino2=="Aspartic Acid")||(Amino2=="Isoleucine"&&Amino1=="Aspartic Acid")){
    if(MatrixMethod=="B") return(-3)
    else return(-2)
  }
  if((Amino1=="Leucine"&&Amino2=="Aspartic Acid")||(Amino2=="Leucine"&&Amino1=="Aspartic Acid")){
    if(MatrixMethod=="B") return(-4)
    else return(-4)
  }
  if((Amino1=="Lysine"&&Amino2=="Aspartic Acid")||(Amino2=="Lysine"&&Amino1=="Aspartic Acid")){
    if(MatrixMethod=="B") return(-1)
    else return(0)
  }
  if((Amino1=="Methionine"&&Amino2=="Aspartic Acid")||(Amino2=="Methionine"&&Amino1=="Aspartic Acid")){
    if(MatrixMethod=="B") return(-3)
    else return(-3)
  }
  #Should Amino1 be Phenylalanine?
  if((Amino1=="Proline"&&Amino2=="Aspartic Acid")||(Amino2=="Phenylalanine"&&Amino1=="Aspartic Acid")){
    if(MatrixMethod=="B") return(-3)
    else return(-6)
  }
  if((Amino1=="Proline"&&Amino2=="Aspartic Acid")||(Amino2=="Proline"&&Amino1=="Aspartic Acid")){
    if(MatrixMethod=="B") return(-1)
    else return(-1)
  }
  if((Amino1=="Serine"&&Amino2=="Aspartic Acid")||(Amino2=="Serine"&&Amino1=="Aspartic Acid")){
    if(MatrixMethod=="B") return(0)
    else return(0)
  }
  if((Amino1=="Threonine"&&Amino2=="Aspartic Acid")||(Amino2=="Threonine"&&Amino1=="Aspartic Acid")){
    if(MatrixMethod=="B") return(-1)
    else return(0)
  }
  if((Amino1=="Tryptophan"&&Amino2=="Aspartic Acid")||(Amino2=="Tryptophan"&&Amino1=="Aspartic Acid")){
    if(MatrixMethod=="B") return(-4)
    else return(-7)
  }
  if((Amino1=="Tyrosine"&&Amino2=="Aspartic Acid")||(Amino2=="Tyrosine"&&Amino1=="Aspartic Acid")){
    if(MatrixMethod=="B") return(-3)
    else return(-4)
  }
  if((Amino1=="Valine"&&Amino2=="Aspartic Acid")||(Amino2=="Valine"&&Amino1=="Aspartic Acid")){
    if(MatrixMethod=="B") return(-3)
    else return(-2)
  }
  #Cysteine Block
  if((Amino1=="Cysteine"&&Amino2=="Cysteine")||(Amino2=="Cysteine"&&Amino1=="Cysteine")){
    if(MatrixMethod=="B") return(9)
    else return(4)
  }
  if((Amino1=="Glutamine"&&Amino2=="Cysteine")||(Amino2=="Glutamine"&&Amino1=="Cysteine")){
    if(MatrixMethod=="B") return(-3)
    else return(-5)
  }
  if((Amino1=="Glutamic Acid"&&Amino2=="Cysteine")||(Amino2=="Glutamic Acid"&&Amino1=="Cysteine")){
    if(MatrixMethod=="B") return(-4)
    else return(-5)
  }
  if((Amino1=="Glycine"&&Amino2=="Cysteine")||(Amino2=="Glycine"&&Amino1=="Cysteine")){
    if(MatrixMethod=="B") return(-3)
    else return(-3)
  }
  if((Amino1=="Histidine"&&Amino2=="Cysteine")||(Amino2=="Histidine"&&Amino1=="Cysteine")){
    if(MatrixMethod=="B") return(-3)
    else return(-3)
  }
  if((Amino1=="Isoleucine"&&Amino2=="Cysteine")||(Amino2=="Isoleucine"&&Amino1=="Cysteine")){
    if(MatrixMethod=="B") return(-1)
    else return(-2)
  }
  if((Amino1=="Leucine"&&Amino2=="Cysteine")||(Amino2=="Leucine"&&Amino1=="Cysteine")){
    if(MatrixMethod=="B") return(-1)
    else return(-6)
  }
  if((Amino1=="Lysine"&&Amino2=="Cysteine")||(Amino2=="Lysine"&&Amino1=="Cysteine")){
    if(MatrixMethod=="B") return(-3)
    else return(-5)
  }
  if((Amino1=="Methionine"&&Amino2=="Cysteine")||(Amino2=="Methionine"&&Amino1=="Cysteine")){
    if(MatrixMethod=="B") return(-1)
    else return(-5)
  }
  if((Amino1=="Phenylalanine"&&Amino2=="Cysteine")||(Amino2=="Phenylalanine"&&Amino1=="Cysteine")){
    if(MatrixMethod=="B") return(-2)
    else return(-4)
  }
  if((Amino1=="Proline"&&Amino2=="Cysteine")||(Amino2=="Proline"&&Amino1=="Cysteine")){
    if(MatrixMethod=="B") return(-3)
    else return(-3)
  }
  if((Amino1=="Serine"&&Amino2=="Cysteine")||(Amino2=="Serine"&&Amino1=="Cysteine")){
    if(MatrixMethod=="B") return(-1)
    else return(0)
  }
  if((Amino1=="Threonine"&&Amino2=="Cysteine")||(Amino2=="Threonine"&&Amino1=="Cysteine")){
    if(MatrixMethod=="B") return(-1)
    else return(-2)
  }
  if((Amino1=="Tryptophan"&&Amino2=="Cysteine")||(Amino2=="Tryptophan"&&Amino1=="Cysteine")){
    if(MatrixMethod=="B") return(-2)
    else return(-8)
  }
  if((Amino1=="Tyrosine"&&Amino2=="Cysteine")||(Amino2=="Tyrosine"&&Amino1=="Cysteine")){
    if(MatrixMethod=="B") return(-2)
    else return(0)
  }
  if((Amino1=="Valine"&&Amino2=="Cysteine")||(Amino2=="Valine"&&Amino1=="Cysteine")){
    if(MatrixMethod=="B") return(-1)
    else return(-2)
  }
  #Glu Block
  if((Amino1=="Glutamine"&&Amino2=="Glutamine")||(Amino2=="Glutamine"&&Amino1=="Glutamine")){
    if(MatrixMethod=="B") return(5)
    else return(4)
  }
  if((Amino1=="Glutamic Acid"&&Amino2=="Glutamine")||(Amino2=="Glutamic Acid"&&Amino1=="Glutamine")){
    if(MatrixMethod=="B") return(2)
    else return(2)
  }
  if((Amino1=="Glycine"&&Amino2=="Glutamine")||(Amino2=="Glycine"&&Amino1=="Glutamine")){
    if(MatrixMethod=="B") return(-2)
    else return(-1)
  }
  if((Amino1=="Histidine"&&Amino2=="Glutamine")||(Amino2=="Histidine"&&Amino1=="Glutamine")){
    if(MatrixMethod=="B") return(0)
    else return(3)
  }
  if((Amino1=="Isoleucine"&&Amino2=="Glutamine")||(Amino2=="Isoleucine"&&Amino1=="Glutamine")){
    if(MatrixMethod=="B") return(-3)
    else return(-2)
  }
  if((Amino1=="Leucine"&&Amino2=="Glutamine")||(Amino2=="Leucine"&&Amino1=="Glutamine")){
    if(MatrixMethod=="B") return(-2)
    else return(-2)
  }
  if((Amino1=="Lysine"&&Amino2=="Glutamine")||(Amino2=="Lysine"&&Amino1=="Glutamine")){
    if(MatrixMethod=="B") return(1)
    else return(1)
  }
  if((Amino1=="Methionine"&&Amino2=="Glutamine")||(Amino2=="Methionine"&&Amino1=="Glutamine")){
    if(MatrixMethod=="B") return(0)
    else return(-1)
  }
  if((Amino1=="Phenylalanine"&&Amino2=="Glutamine")||(Amino2=="Phenylalanine"&&Amino1=="Glutamine")){
    if(MatrixMethod=="B") return(-3)
    else return(-5)
  }
  if((Amino1=="Proline"&&Amino2=="Glutamine")||(Amino2=="Proline"&&Amino1=="Glutamine")){
    if(MatrixMethod=="B") return(-1)
    else return(0)
  }
  if((Amino1=="Serine"&&Amino2=="Glutamine")||(Amino2=="Serine"&&Amino1=="Glutamine")){
    if(MatrixMethod=="B") return(0)
    else return(-1)
  }
  if((Amino1=="Threonine"&&Amino2=="Glutamine")||(Amino2=="Threonine"&&Amino1=="Glutamine")){
    if(MatrixMethod=="B") return(-1)
    else return(-1)
  }
  if((Amino1=="Tryptophan"&&Amino2=="Glutamine")||(Amino2=="Tryptophan"&&Amino1=="Glutamine")){
    if(MatrixMethod=="B") return(-2)
    else return(-5)
  }
  if((Amino1=="Tyrosine"&&Amino2=="Glutamine")||(Amino2=="Tyrosine"&&Amino1=="Glutamine")){
    if(MatrixMethod=="B") return(-1)
    else return(-4)
  }
  if((Amino1=="Valine"&&Amino2=="Glutamine")||(Amino2=="Valine"&&Amino1=="Glutamine")){
    if(MatrixMethod=="B") return(-2)
    else return(-2)
  }
  
  #Glutamic Acid Block
  if((Amino1=="Glutamic Acid"&&Amino2=="Glutamic Acid")||(Amino2=="Glutamic Acid"&&Amino1=="Glutamic Acid")){
    if(MatrixMethod=="B") return(5)
    else return(4)
  }
  if((Amino1=="Glycine"&&Amino2=="Glutamic Acid")||(Amino2=="Glycine"&&Amino1=="Glutamic Acid")){
    if(MatrixMethod=="B") return(-2)
    else return(0)
  }
  if((Amino1=="Histidine"&&Amino2=="Glutamic Acid")||(Amino2=="Histidine"&&Amino1=="Glutamic Acid")){
    if(MatrixMethod=="B") return(0)
    else return(1)
  }
  if((Amino1=="Isoleucine"&&Amino2=="Glutamic Acid")||(Amino2=="Isoleucine"&&Amino1=="Glutamic Acid")){
    if(MatrixMethod=="B") return(-3)
    else return(-2)
  }
  if((Amino1=="Leucine"&&Amino2=="Glutamic Acid")||(Amino2=="Leucine"&&Amino1=="Glutamic Acid")){
    if(MatrixMethod=="B") return(-3)
    else return(-3)
  }
  if((Amino1=="Lysine"&&Amino2=="Glutamic Acid")||(Amino2=="Lysine"&&Amino1=="Glutamic Acid")){
    if(MatrixMethod=="B") return(1)
    else return(0)
  }
  if((Amino1=="Methionine"&&Amino2=="Glutamic Acid")||(Amino2=="Methionine"&&Amino1=="Glutamic Acid")){
    if(MatrixMethod=="B") return(-2)
    else return(-2)
  }
  if((Amino1=="Phenylalanine"&&Amino2=="Glutamic Acid")||(Amino2=="Phenylalanine"&&Amino1=="Glutamic Acid")){
    if(MatrixMethod=="B") return(-3)
    else return(-5)
  }
  if((Amino1=="Proline"&&Amino2=="Glutamic Acid")||(Amino2=="Proline"&&Amino1=="Glutamic Acid")){
    if(MatrixMethod=="B") return(-1)
    else return(-1)
  }
  if((Amino1=="Serine"&&Amino2=="Glutamic Acid")||(Amino2=="Serine"&&Amino1=="Glutamic Acid")){
    if(MatrixMethod=="B") return(0)
    else return(0)
  }
  if((Amino1=="Threonine"&&Amino2=="Glutamic Acid")||(Amino2=="Threonine"&&Amino1=="Glutamic Acid")){
    if(MatrixMethod=="B") return(-1)
    else return(0)
  }
  if((Amino1=="Tryptophan"&&Amino2=="Glutamic Acid")||(Amino2=="Tryptophan"&&Amino1=="Glutamic Acid")){
    if(MatrixMethod=="B") return(-3)
    else return(-7)
  }
  if((Amino1=="Tyrosine"&&Amino2=="Glutamic Acid")||(Amino2=="Tyrosine"&&Amino1=="Glutamic Acid")){
    if(MatrixMethod=="B") return(-2)
    else return(-4)
  }
  if((Amino1=="Valine"&&Amino2=="Glutamic Acid")||(Amino2=="Valine"&&Amino1=="Glutamic Acid")){
    if(MatrixMethod=="B") return(-2)
    else return(-2)
  }
  #Glycine Block
  if((Amino1=="Glycine"&&Amino2=="Glycine")||(Amino2=="Glycine"&&Amino1=="Glycine")){
    if(MatrixMethod=="B") return(6)
    else return(5)
  }
  if((Amino1=="Histidine"&&Amino2=="Glycine")||(Amino2=="Histidine"&&Amino1=="Glycine")){
    if(MatrixMethod=="B") return(-2)
    else return(-2)
  }
  if((Amino1=="Isoleucine"&&Amino2=="Glycine")||(Amino2=="Isoleucine"&&Amino1=="Glycine")){
    if(MatrixMethod=="B") return(-4)
    else return(-3)
  }
  if((Amino1=="Leucine"&&Amino2=="Glycine")||(Amino2=="Leucine"&&Amino1=="Glycine")){
    if(MatrixMethod=="B") return(-4)
    else return(-4)
  }
  if((Amino1=="Lysine"&&Amino2=="Glycine")||(Amino2=="Lysine"&&Amino1=="Glycine")){
    if(MatrixMethod=="B") return(-2)
    else return(-2)
  }
  if((Amino1=="Methionine"&&Amino2=="Glycine")||(Amino2=="Methionine"&&Amino1=="Glycine")){
    if(MatrixMethod=="B") return(-3)
    else return(-3)
  }
  if((Amino1=="Phenylalanine"&&Amino2=="Glycine")||(Amino2=="Phenylalanine"&&Amino1=="Glycine")){
    if(MatrixMethod=="B") return(-3)
    else return(-5)
  }
  if((Amino1=="Proline"&&Amino2=="Glycine")||(Amino2=="Proline"&&Amino1=="Glycine")){
    if(MatrixMethod=="B") return(-2)
    else return(-1)
  }
  if((Amino1=="Serine"&&Amino2=="Glycine")||(Amino2=="Serine"&&Amino1=="Glycine")){
    if(MatrixMethod=="B") return(0)
    else return(1)
  }
  if((Amino1=="Threonine"&&Amino2=="Glycine")||(Amino2=="Threonine"&&Amino1=="Glycine")){
    if(MatrixMethod=="B") return(-2)
    else return(0)
  }
  if((Amino1=="Tryptophan"&&Amino2=="Glycine")||(Amino2=="Tryptophan"&&Amino1=="Glycine")){
    if(MatrixMethod=="B") return(-2)
    else return(-7)
  }
  if((Amino1=="Tyrosine"&&Amino2=="Glycine")||(Amino2=="Tyrosine"&&Amino1=="Glycine")){
    if(MatrixMethod=="B") return(-3)
    else return(-5)
  }
  if((Amino1=="Valine"&&Amino2=="Glycine")||(Amino2=="Valine"&&Amino1=="Glycine")){
    if(MatrixMethod=="B") return(-3)
    else return(-1)
  }
  #Histidine Block
  if((Amino1=="Histidine"&&Amino2=="Histidine")||(Amino2=="Histidine"&&Amino1=="Histidine")){
    if(MatrixMethod=="B") return(8)
    else return(6)
  }
  if((Amino1=="Isoleucine"&&Amino2=="Histidine")||(Amino2=="Isoleucine"&&Amino1=="Histidine")){
    if(MatrixMethod=="B") return(-3)
    else return(-2)
  }
  if((Amino1=="Leucine"&&Amino2=="Histidine")||(Amino2=="Leucine"&&Amino1=="Histidine")){
    if(MatrixMethod=="B") return(-3)
    else return(-2)
  }
  if((Amino1=="Lysine"&&Amino2=="Histidine")||(Amino2=="Lysine"&&Amino1=="Histidine")){
    if(MatrixMethod=="B") return(-1)
    else return(0)
  }
  if((Amino1=="Methionine"&&Amino2=="Histidine")||(Amino2=="Methionine"&&Amino1=="Histidine")){
    if(MatrixMethod=="B") return(2)
    else return(-2)
  }
  if((Amino1=="Phenylalanine"&&Amino2=="Histidine")||(Amino2=="Phenylalanine"&&Amino1=="Histidine")){
    if(MatrixMethod=="B") return(-1)
    else return(-2)
  }
  if((Amino1=="Proline"&&Amino2=="Histidine")||(Amino2=="Proline"&&Amino1=="Histidine")){
    if(MatrixMethod=="B") return(-2)
    else return(0)
  }
  if((Amino1=="Serine"&&Amino2=="Histidine")||(Amino2=="Serine"&&Amino1=="Histidine")){
    if(MatrixMethod=="B") return(-1)
    else return(1)
  }
  if((Amino1=="Threonine"&&Amino2=="Histidine")||(Amino2=="Threonine"&&Amino1=="Histidine")){
    if(MatrixMethod=="B") return(-2)
    else return(-1)
  }
  if((Amino1=="Tryptophan"&&Amino2=="Histidine")||(Amino2=="Tryptophan"&&Amino1=="Histidine")){
    if(MatrixMethod=="B") return(-2)
    else return(-3)
  }
  if((Amino1=="Tyrosine"&&Amino2=="Histidine")||(Amino2=="Tyrosine"&&Amino1=="Histidine")){
    if(MatrixMethod=="B") return(2)
    else return(0)
  }
  if((Amino1=="Valine"&&Amino2=="Histidine")||(Amino2=="Valine"&&Amino1=="Histidine")){
    if(MatrixMethod=="B") return(-3)
    else return(-2)
  }
  
  #Isoleucine
  if((Amino1=="Isoleucine"&&Amino2=="Isoleucine")||(Amino2=="Isoleucine"&&Amino1=="Isoleucine")){
    if(MatrixMethod=="B") return(4)
    else return(5)
  }
  if((Amino1=="Leucine"&&Amino2=="Isoleucine")||(Amino2=="Leucine"&&Amino1=="Isoleucine")){
    if(MatrixMethod=="B") return(2)
    else return(2)
  }
  if((Amino1=="Lysine"&&Amino2=="Isoleucine")||(Amino2=="Lysine"&&Amino1=="Isoleucine")){
    if(MatrixMethod=="B") return(-3)
    else return(-2)
  }
  if((Amino1=="Methionine"&&Amino2=="Isoleucine")||(Amino2=="Methionine"&&Amino1=="Isoleucine")){
    if(MatrixMethod=="B") return(1)
    else return(2)
  }
  if((Amino1=="Phenylalanine"&&Amino2=="Isoleucine")||(Amino2=="Phenylalanine"&&Amino1=="Isoleucine")){
    if(MatrixMethod=="B") return(0)
    else return(1)
  }
  if((Amino1=="Proline"&&Amino2=="Isoleucine")||(Amino2=="Proline"&&Amino1=="Isoleucine")){
    if(MatrixMethod=="B") return(-3)
    else return(-2)
  }
  if((Amino1=="Serine"&&Amino2=="Isoleucine")||(Amino2=="Serine"&&Amino1=="Isoleucine")){
    if(MatrixMethod=="B") return(-2)
    else return(-1)
  }
  if((Amino1=="Threonine"&&Amino2=="Isoleucine")||(Amino2=="Threonine"&&Amino1=="Isoleucine")){
    if(MatrixMethod=="B") return(-1)
    else return(0)
  }
  if((Amino1=="Tryptophan"&&Amino2=="Isoleucine")||(Amino2=="Tryptophan"&&Amino1=="Isoleucine")){
    if(MatrixMethod=="B") return(-3)
    else return(-5)
  }
  if((Amino1=="Tyrosine"&&Amino2=="Isoleucine")||(Amino2=="Tyrosine"&&Amino1=="Isoleucine")){
    if(MatrixMethod=="B") return(-1)
    else return(-1)
  }
  if((Amino1=="Valine"&&Amino2=="Isoleucine")||(Amino2=="Valine"&&Amino1=="Isoleucine")){
    if(MatrixMethod=="B") return(3)
    else return(4)
  }
  #Leucine Block
  if((Amino1=="Leucine"&&Amino2=="Leucine")||(Amino2=="Leucine"&&Amino1=="Leucine")){
    if(MatrixMethod=="B") return(4)
    else return(6)
  }
  if((Amino1=="Lysine"&&Amino2=="Leucine")||(Amino2=="Lysine"&&Amino1=="Leucine")){
    if(MatrixMethod=="B") return(-2)
    else return(-3)
  }
  if((Amino1=="Methionine"&&Amino2=="Leucine")||(Amino2=="Methionine"&&Amino1=="Leucine")){
    if(MatrixMethod=="B") return(2)
    else return(4)
  }
  if((Amino1=="Phenylalanine"&&Amino2=="Leucine")||(Amino2=="Phenylalanine"&&Amino1=="Leucine")){
    if(MatrixMethod=="B") return(0)
    else return(2)
  }
  if((Amino1=="Proline"&&Amino2=="Leucine")||(Amino2=="Proline"&&Amino1=="Leucine")){
    if(MatrixMethod=="B") return(-3)
    else return(-3)
  }
  if((Amino1=="Serine"&&Amino2=="Leucine")||(Amino2=="Serine"&&Amino1=="Leucine")){
    if(MatrixMethod=="B") return(-2)
    else return(-3)
  }
  if((Amino1=="Threonine"&&Amino2=="Leucine")||(Amino2=="Threonine"&&Amino1=="Leucine")){
    if(MatrixMethod=="B") return(-1)
    else return(-2)
  }
  if((Amino1=="Tryptophan"&&Amino2=="Leucine")||(Amino2=="Tryptophan"&&Amino1=="Leucine")){
    if(MatrixMethod=="B") return(-2)
    else return(-2)
  }
  if((Amino1=="Tyrosine"&&Amino2=="Leucine")||(Amino2=="Tyrosine"&&Amino1=="Leucine")){
    if(MatrixMethod=="B") return(-1)
    else return(-1)
  }
  if((Amino1=="Valine"&&Amino2=="Leucine")||(Amino2=="Valine"&&Amino1=="Leucine")){
    if(MatrixMethod=="B") return(1)
    else return(2)
  }
  #Lysine
  if((Amino1=="Lysine"&&Amino2=="Lysine")||(Amino2=="Lysine"&&Amino1=="Lysine")){
    if(MatrixMethod=="B") return(5)
    else return(5)
  }
  if((Amino1=="Methionine"&&Amino2=="Lysine")||(Amino2=="Methionine"&&Amino1=="Lysine")){
    if(MatrixMethod=="B") return(-1)
    else return(0)
  }
  if((Amino1=="Phenylalanine"&&Amino2=="Lysine")||(Amino2=="Phenylalanine"&&Amino1=="Lysine")){
    if(MatrixMethod=="B") return(-3)
    else return(-5)
  }
  if((Amino1=="Proline"&&Amino2=="Lysine")||(Amino2=="Proline"&&Amino1=="Lysine")){
    if(MatrixMethod=="B") return(-1)
    else return(-1)
  }
  if((Amino1=="Serine"&&Amino2=="Lysine")||(Amino2=="Serine"&&Amino1=="Lysine")){
    if(MatrixMethod=="B") return(0)
    else return(0)
  }
  if((Amino1=="Threonine"&&Amino2=="Lysine")||(Amino2=="Threonine"&&Amino1=="Lysine")){
    if(MatrixMethod=="B") return(-1)
    else return(0)
  }
  if((Amino1=="Tryptophan"&&Amino2=="Lysine")||(Amino2=="Tryptophan"&&Amino1=="Lysine")){
    if(MatrixMethod=="B") return(-3)
    else return(-3)
  }
  if((Amino1=="Tyrosine"&&Amino2=="Lysine")||(Amino2=="Tyrosine"&&Amino1=="Lysine")){
    if(MatrixMethod=="B") return(-2)
    else return(-4)
  }
  if((Amino1=="Valine"&&Amino2=="Lysine")||(Amino2=="Valine"&&Amino1=="Lysine")){
    if(MatrixMethod=="B") return(-2)
    else return(-2)
  }
  
  #methionine
  if((Amino1=="Methionine"&&Amino2=="Methionine")||(Amino2=="Methionine"&&Amino1=="Methionine")){
    if(MatrixMethod=="B") return(5)
    else return(6)
  }
  if((Amino1=="Phenylalanine"&&Amino2=="Methionine")||(Amino2=="Phenylalanine"&&Amino1=="Methionine")){
    if(MatrixMethod=="B") return(0)
    else return(0)
  }
  if((Amino1=="Proline"&&Amino2=="Methionine")||(Amino2=="Proline"&&Amino1=="Methionine")){
    if(MatrixMethod=="B") return(-2)
    else return(-2)
  }
  if((Amino1=="Serine"&&Amino2=="Methionine")||(Amino2=="Serine"&&Amino1=="Methionine")){
    if(MatrixMethod=="B") return(-1)
    else return(-2)
  }
  if((Amino1=="Threonine"&&Amino2=="Methionine")||(Amino2=="Threonine"&&Amino1=="Methionine")){
    if(MatrixMethod=="B") return(-1)
    else return(-1)
  }
  if((Amino1=="Tryptophan"&&Amino2=="Methionine")||(Amino2=="Tryptophan"&&Amino1=="Methionine")){
    if(MatrixMethod=="B") return(-1)
    else return(-4)
  }
  if((Amino1=="Tyrosine"&&Amino2=="Methionine")||(Amino2=="Tyrosine"&&Amino1=="Methionine")){
    if(MatrixMethod=="B") return(-1)
    else return(-2)
  }
  if((Amino1=="Valine"&&Amino2=="Methionine")||(Amino2=="Valine"&&Amino1=="Methionine")){
    if(MatrixMethod=="B") return(1)
    else return(2)
  }
  #phenylalanine
  if((Amino1=="Phenylalanine"&&Amino2=="Phenylalanine")||(Amino2=="Phenylalanine"&&Amino1=="Phenylalanine")){
    if(MatrixMethod=="B") return(6)
    else return(9)
  }
  if((Amino1=="Proline"&&Amino2=="Phenylalanine")||(Amino2=="Proline"&&Amino1=="Phenylalanine")){
    if(MatrixMethod=="B") return(-4)
    else return(-5)
  }
  if((Amino1=="Serine"&&Amino2=="Phenylalanine")||(Amino2=="Serine"&&Amino1=="Phenylalanine")){
    if(MatrixMethod=="B") return(-2)
    else return(-3)
  }
  if((Amino1=="Threonine"&&Amino2=="Phenylalanine")||(Amino2=="Threonine"&&Amino1=="Phenylalanine")){
    if(MatrixMethod=="B") return(-2)
    else return(-2)
  }
  if((Amino1=="Tryptophan"&&Amino2=="Phenylalanine")||(Amino2=="Tryptophan"&&Amino1=="Phenylalanine")){
    if(MatrixMethod=="B") return(1)
    else return(0)
  }
  if((Amino1=="Tyrosine"&&Amino2=="Phenylalanine")||(Amino2=="Tyrosine"&&Amino1=="Phenylalanine")){
    if(MatrixMethod=="B") return(3)
    else return(7)
  }
  if((Amino1=="Valine"&&Amino2=="Phenylalanine")||(Amino2=="Valine"&&Amino1=="Phenylalanine")){
    if(MatrixMethod=="B") return(-1)
    else return(-1)
  }
  #Proline
  if((Amino1=="Proline"&&Amino2=="Proline")||(Amino2=="Proline"&&Amino1=="Proline")){
    if(MatrixMethod=="B") return(7)
    else return(6)
  }
  if((Amino1=="Serine"&&Amino2=="Proline")||(Amino2=="Serine"&&Amino1=="Proline")){
    if(MatrixMethod=="B") return(-1)
    else return(1)
  }
  if((Amino1=="Threonine"&&Amino2=="Proline")||(Amino2=="Threonine"&&Amino1=="Proline")){
    if(MatrixMethod=="B") return(-1)
    else return(1)
  }
  if((Amino1=="Tryptophan"&&Amino2=="Proline")||(Amino2=="Tryptophan"&&Amino1=="Proline")){
    if(MatrixMethod=="B") return(-4)
    else return(-6)
  }
  if((Amino1=="Tyrosine"&&Amino2=="Proline")||(Amino2=="Tyrosine"&&Amino1=="Proline")){
    if(MatrixMethod=="B") return(-3)
    else return(-5)
  }
  if((Amino1=="Valine"&&Amino2=="Proline")||(Amino2=="Valine"&&Amino1=="Proline")){
    if(MatrixMethod=="B") return(-2)
    else return(-1)
  }
  #Serine
  if((Amino1=="Serine"&&Amino2=="Serine")||(Amino2=="Serine"&&Amino1=="Serine")){
    if(MatrixMethod=="B") return(4)
    else return(3)
  }
  if((Amino1=="Threonine"&&Amino2=="Serine")||(Amino2=="Threonine"&&Amino1=="Serine")){
    if(MatrixMethod=="B") return(1)
    else return(1)
  }
  if((Amino1=="Tryptophan"&&Amino2=="Serine")||(Amino2=="Tryptophan"&&Amino1=="Serine")){
    if(MatrixMethod=="B") return(-3)
    else return(-2)
  }
  if((Amino1=="Tyrosine"&&Amino2=="Serine")||(Amino2=="Tyrosine"&&Amino1=="Serine")){
    if(MatrixMethod=="B") return(-2)
    else return(-3)
  }
  if((Amino1=="Valine"&&Amino2=="Serine")||(Amino2=="Valine"&&Amino1=="Serine")){
    if(MatrixMethod=="B") return(-2)
    else return(-1)
  }
  #Threonine
  if((Amino1=="Threonine"&&Amino2=="Threonine")||(Amino2=="Threonine"&&Amino1=="Threonine")){
    if(MatrixMethod=="B") return(5)
    else return(3)
  }
  if((Amino1=="Tryptophan"&&Amino2=="Threonine")||(Amino2=="Tryptophan"&&Amino1=="Threonine")){
    if(MatrixMethod=="B") return(-2)
    else return(-5)
  }
  if((Amino1=="Tyrosine"&&Amino2=="Threonine")||(Amino2=="Tyrosine"&&Amino1=="Threonine")){
    if(MatrixMethod=="B") return(-2)
    else return(-3)
  }
  if((Amino1=="Valine"&&Amino2=="Threonine")||(Amino2=="Valine"&&Amino1=="Threonine")){
    if(MatrixMethod=="B") return(0)
    else return(0)
  }
  #Tryptophan
  if((Amino1=="Tryptophan"&&Amino2=="Tryptophan")||(Amino2=="Tryptophan"&&Amino1=="Tryptophan")){
    if(MatrixMethod=="B") return(11)
    else return(17)
  }
  if((Amino1=="Tyrosine"&&Amino2=="Tryptophan")||(Amino2=="Tyrosine"&&Amino1=="Tryptophan")){
    if(MatrixMethod=="B") return(2)
    else return(0)
  }
  if((Amino1=="Valine"&&Amino2=="Tryptophan")||(Amino2=="Valine"&&Amino1=="Tryptophan")){
    if(MatrixMethod=="B") return(-3)
    else return(-6)
  }
  #tyrosine
  if((Amino1=="Tyrosine"&&Amino2=="Tyrosine")||(Amino2=="Tyrosine"&&Amino1=="Tyrosine")){
    if(MatrixMethod=="B") return(7)
    else return(10)
  }
  if((Amino1=="Valine"&&Amino2=="Tyrosine")||(Amino2=="Valine"&&Amino1=="Tyrosine")){
    if(MatrixMethod=="B") return(-1)
    else return(-2)
  }
  #valine
  if((Amino1=="Valine"&&Amino2=="Valine")||(Amino2=="Valine"&&Amino1=="Valine")){
    if(MatrixMethod=="B") return(4)
    else return(4)
  }
  #else{
  # print(Amino1)
  #print(Amino2)
  #return(-20)
  #}
} 
#Should be able to delete, fully implemented in BLOSOM method
PAM<- function(Amino1, Amino2){
  #converts amion acid pair to the corresponding PAM250 score
  if ((Amino1=="Stop")||(Amino2=="Stop")){
    
  }
  if((Amino1=="Alanine"&&Amino2=="Alanine")||(Amino2=="Alanine"&&Amino1=="Alanine")){
    return(2)
  }
  if((Amino1=="Arginine"&&Amino2=="Alanine")||(Amino2=="Arginine"&&Amino1=="Alanine")){
    return(-2)
  }
  if((Amino1=="Asparagine"&&Amino2=="Alanine")||(Amino2=="Asparagine"&&Amino1=="Alanine")){
    return(0)
  }
  if((Amino1=="Aspartic Acid"&&Amino2=="Alanine")||(Amino2=="Aspartic Acid"&&Amino1=="Alanine")){
    return(0)
  }
  if((Amino1=="Cysteine"&&Amino2=="Alanine")||(Amino2=="Cysteine"&&Amino1=="Alanine")){
    return(-2)
  }
  if((Amino1=="Glutamine"&&Amino2=="Alanine")||(Amino2=="Glutamine"&&Amino1=="Alanine")){
    return(0)
  }
  if((Amino1=="Glutamic Acid"&&Amino2=="Alanine")||(Amino2=="Glutamic Acid"&&Amino1=="Alanine")){
    return(0)
  }
  if((Amino1=="Glycine"&&Amino2=="Alanine")||(Amino2=="Glycine"&&Amino1=="Alanine")){
    return(1)
  }
  if((Amino1=="Histidine"&&Amino2=="Alanine")||(Amino2=="Histidine"&&Amino1=="Alanine")){
    return(-1)
  }
  if((Amino1=="Isoleucine"&&Amino2=="Alanine")||(Amino2=="Isoleucine"&&Amino1=="Alanine")){
    return(-1)
  }
  if((Amino1=="Lysine"&&Amino2=="Alanine")||(Amino2=="Lysine"&&Amino1=="Alanine")){
    return(-2)
  }
  if((Amino1=="Methionine"&&Amino2=="Alanine")||(Amino2=="Methionine"&&Amino1=="Alanine")){
    return(-1)
  }
  if((Amino1=="Phenylalanine"&&Amino2=="Alanine")||(Amino2=="Phenylalanine"&&Amino1=="Alanine")){
    return(-4)
  }
  if((Amino1=="Proline"&&Amino2=="Alanine")||(Amino2=="Proline"&&Amino1=="Alanine")){
    return(1)
  }
  if((Amino1=="Serine"&&Amino2=="Alanine")||(Amino2=="Serine"&&Amino1=="Alanine")){
    return(1)
  }
  if((Amino1=="Threonine"&&Amino2=="Alanine")||(Amino2=="Threonine"&&Amino1=="Alanine")){
    return(1)
  }
  if((Amino1=="Tryptophan"&&Amino2=="Alanine")||(Amino2=="Tryptophan"&&Amino1=="Alanine")){
    return(-6)
  }
  if((Amino1=="Tyrosine"&&Amino2=="Alanine")||(Amino2=="Tyrosine"&&Amino1=="Alanine")){
    return(-3)
  }
  if((Amino1=="Valine"&&Amino2=="Alanine")||(Amino2=="Valine"&&Amino1=="Alanine")){
    return(0)
  }
  #Arganine Block
  if((Amino1=="Arginine"&& Amino2=="Arginine")||(Amino2=="Arginine" && Amino1=="Arginine")){
    return(6)
  }
  if((Amino1=="Asparagine"&&Amino2=="Arginine")||(Amino2=="Asparagine"&&Amino1=="Arginine")){
    return(0)
  }
  if((Amino1=="Aspartic Acid"&&Amino2=="Arginine")||(Amino2=="Aspartic Acid"&&Amino1=="Arginine")){
    return(-1)
  }
  if((Amino1=="Cysteine"&&Amino2=="Arginine")||(Amino2=="Cysteine"&&Amino1=="Arginine")){
    return(-4)
  }
  if((Amino1=="Glutamine"&&Amino2=="Arginine")||(Amino2=="Glutamine"&&Amino1=="Arginine")){
    return(1)
  }
  if((Amino1=="Glutamic Acid"&&Amino2=="Arginine")||(Amino2=="Glutamic Acid"&&Amino1=="Arginine")){
    return(-1)
  }
  if((Amino1=="Glycine"&&Amino2=="Arginine")||(Amino2=="Glycine"&&Amino1=="Arginine")){
    return(-3)
  }
  if((Amino1=="Histidine"&&Amino2=="Arginine")||(Amino2=="Histidine"&&Amino1=="Arginine")){
    return(2)
  }
  if((Amino1=="Isoleucine"&&Amino2=="Arginine")||(Amino2=="Isoleucine"&&Amino1=="Arginine")){
    return(-2)
  }
  if((Amino1=="Leucine"&&Amino2=="Arginine")||(Amino2=="Leucine"&&Amino1=="Arginine")){
    return(-3)
  }
  if((Amino1=="Lysine"&&Amino2=="Arginine")||(Amino2=="Lysine"&&Amino1=="Arginine")){
    return(3)
  }
  if((Amino1=="Methionine"&&Amino2=="Arginine")||(Amino2=="Methionine"&&Amino1=="Arginine")){
    return(0)
  }
  if((Amino1=="Phenylalanine"&&Amino2=="Arginine")||(Amino2=="Phenylalanine"&&Amino1=="Arginine")){
    return(-4)
  }
  if((Amino1=="Proline"&&Amino2=="Arginine")||(Amino2=="Proline"&&Amino1=="Arginine")){
    return(0)
  }
  if((Amino1=="Serine"&&Amino2=="Arginine")||(Amino2=="Serine"&&Amino1=="Arginine")){
    return(0)
  }
  if((Amino1=="Threonine"&&Amino2=="Arginine")||(Amino2=="Threonine"&&Amino1=="Arginine")){
    return(-1)
  }
  if((Amino1=="Tryptophan"&&Amino2=="Arginine")||(Amino2=="Tryptophan"&&Amino1=="Arginine")){
    return(2)
  }
  if((Amino1=="Tyrosine"&&Amino2=="Arginine")||(Amino2=="Tyrosine"&&Amino1=="Arginine")){
    return(-4)
  }
  if((Amino1=="Valine"&&Amino2=="Arginine")||(Amino2=="Valine"&&Amino1=="Arginine")){
    return(-2)
  }
  #aspargine Block
  if((Amino1=="Asparagine"&&Amino2=="Asparagine")||(Amino2=="Asparagine"&&Amino1=="Asparagine")){
    return(2)
  }
  if((Amino1=="Aspartic Acid"&&Amino2=="Asparagine")||(Amino2=="Aspartic Acid"&&Amino1=="Asparagine")){
    return(2)
  }
  if((Amino1=="Cysteine"&&Amino2=="Asparagine")||(Amino2=="Cysteine"&&Amino1=="Asparagine")){
    return(-4)
  }
  if((Amino1=="Glutamine"&&Amino2=="Asparagine")||(Amino2=="Glutamine"&&Amino1=="Asparagine")){
    return(1)
  }
  if((Amino1=="Glutamic Acid"&&Amino2=="Asparagine")||(Amino2=="Glutamic Acid"&&Amino1=="Asparagine")){
    return(1)
  }
  if((Amino1=="Glycine"&&Amino2=="Asparagine")||(Amino2=="Glycine"&&Amino1=="Asparagine")){
    return(0)
  }
  if((Amino1=="Histidine"&&Amino2=="Asparagine")||(Amino2=="Histidine"&&Amino1=="Asparagine")){
    return(2)
  }
  if((Amino1=="Isoleucine"&&Amino2=="Asparagine")||(Amino2=="Isoleucine"&&Amino1=="Asparagine")){
    return(-2)
  }
  if((Amino1=="Leucine"&&Amino2=="Asparagine")||(Amino2=="Leucine"&&Amino1=="Asparagine")){
    return(-3)
  }
  if((Amino1=="Lysine"&&Amino2=="Asparagine")||(Amino2=="Lysine"&&Amino1=="Asparagine")){
    return(1)
  }
  if((Amino1=="Methionine"&&Amino2=="Asparagine")||(Amino2=="Methionine"&&Amino1=="Asparagine")){
    return(-2)
  }
  if((Amino1=="Phenylalanine"&&Amino2=="Asparagine")||(Amino2=="Phenylalanine"&&Amino1=="Asparagine")){
    return(-4)
  }
  if((Amino1=="Proline"&&Amino2=="Asparagine")||(Amino2=="Proline"&&Amino1=="Asparagine")){
    return(-1)
  }
  if((Amino1=="Serine"&&Amino2=="Asparagine")||(Amino2=="Serine"&&Amino1=="Asparagine")){
    return(1)
  }
  if((Amino1=="Threonine"&&Amino2=="Asparagine")||(Amino2=="Threonine"&&Amino1=="Asparagine")){
    return(0)
  }
  if((Amino1=="Tryptophan"&&Amino2=="Asparagine")||(Amino2=="Tryptophan"&&Amino1=="Asparagine")){
    return(-4)
  }
  if((Amino1=="Tyrosine"&&Amino2=="Asparagine")||(Amino2=="Tyrosine"&&Amino1=="Asparagine")){
    return(-2)
  }
  if((Amino1=="Valine"&&Amino2=="Asparagine")||(Amino2=="Valine"&&Amino1=="Asparagine")){
    return(-2)
  }
  #Aspartic Acid Block
  if((Amino1=="Aspartic Acid"&&Amino2=="Aspartic Acid")||(Amino2=="Aspartic Acid"&&Amino1=="Aspartic Acid")){
    return(4)
  }
  if((Amino1=="Cysteine"&&Amino2=="Aspartic Acid")||(Amino2=="Cysteine"&&Amino1=="Aspartic Acid")){
    return(-5)
  }
  if((Amino1=="Glutamine"&&Amino2=="Aspartic Acid")||(Amino2=="Glutamine"&&Amino1=="Aspartic Acid")){
    return(2)
  }
  if((Amino1=="Glutamic Acid"&&Amino2=="Aspartic Acid")||(Amino2=="Glutamic Acid"&&Amino1=="Aspartic Acid")){
    return(3)
  }
  if((Amino1=="Glycine"&&Amino2=="Aspartic Acid")||(Amino2=="Glycine"&&Amino1=="Aspartic Acid")){
    return(1)
  }
  if((Amino1=="Histidine"&&Amino2=="Aspartic Acid")||(Amino2=="Histidine"&&Amino1=="Aspartic Acid")){
    return(1)
  }
  if((Amino1=="Isoleucine"&&Amino2=="Aspartic Acid")||(Amino2=="Isoleucine"&&Amino1=="Aspartic Acid")){
    return(-2)
  }
  if((Amino1=="Leucine"&&Amino2=="Aspartic Acid")||(Amino2=="Leucine"&&Amino1=="Aspartic Acid")){
    return(-4)
  }
  if((Amino1=="Lysine"&&Amino2=="Aspartic Acid")||(Amino2=="Lysine"&&Amino1=="Aspartic Acid")){
    return(0)
  }
  if((Amino1=="Methionine"&&Amino2=="Aspartic Acid")||(Amino2=="Methionine"&&Amino1=="Aspartic Acid")){
    return(-3)
  }
  if((Amino1=="Proline"&&Amino2=="Aspartic Acid")||(Amino2=="Phenylalanine"&&Amino1=="Aspartic Acid")){
    return(-6)
  }
  if((Amino1=="Proline"&&Amino2=="Aspartic Acid")||(Amino2=="Proline"&&Amino1=="Aspartic Acid")){
    return(-1)
  }
  if((Amino1=="Serine"&&Amino2=="Aspartic Acid")||(Amino2=="Serine"&&Amino1=="Aspartic Acid")){
    return(0)
  }
  if((Amino1=="Threonine"&&Amino2=="Aspartic Acid")||(Amino2=="Threonine"&&Amino1=="Aspartic Acid")){
    return(0)
  }
  if((Amino1=="Tryptophan"&&Amino2=="Aspartic Acid")||(Amino2=="Tryptophan"&&Amino1=="Aspartic Acid")){
    return(-7)
  }
  if((Amino1=="Tyrosine"&&Amino2=="Aspartic Acid")||(Amino2=="Tyrosine"&&Amino1=="Aspartic Acid")){
    return(-4)
  }
  if((Amino1=="Valine"&&Amino2=="Aspartic Acid")||(Amino2=="Valine"&&Amino1=="Aspartic Acid")){
    return(-2)
  }
  #Cysteine Block
  if((Amino1=="Cysteine"&&Amino2=="Cysteine")||(Amino2=="Cysteine"&&Amino1=="Cysteine")){
    return(4)
  }
  if((Amino1=="Glutamine"&&Amino2=="Cysteine")||(Amino2=="Glutamine"&&Amino1=="Cysteine")){
    return(-5)
  }
  if((Amino1=="Glutamic Acid"&&Amino2=="Cysteine")||(Amino2=="Glutamic Acid"&&Amino1=="Cysteine")){
    return(-5)
  }
  if((Amino1=="Glycine"&&Amino2=="Cysteine")||(Amino2=="Glycine"&&Amino1=="Cysteine")){
    return(-3)
  }
  if((Amino1=="Histidine"&&Amino2=="Cysteine")||(Amino2=="Histidine"&&Amino1=="Cysteine")){
    return(-3)
  }
  if((Amino1=="Isoleucine"&&Amino2=="Cysteine")||(Amino2=="Isoleucine"&&Amino1=="Cysteine")){
    return(-2)
  }
  if((Amino1=="Leucine"&&Amino2=="Cysteine")||(Amino2=="Leucine"&&Amino1=="Cysteine")){
    return(-6)
  }
  if((Amino1=="Lysine"&&Amino2=="Cysteine")||(Amino2=="Lysine"&&Amino1=="Cysteine")){
    return(-5)
  }
  if((Amino1=="Methionine"&&Amino2=="Cysteine")||(Amino2=="Methionine"&&Amino1=="Cysteine")){
    return(-5)
  }
  if((Amino1=="Phenylalanine"&&Amino2=="Cysteine")||(Amino2=="Phenylalanine"&&Amino1=="Cysteine")){
    return(-4)
  }
  if((Amino1=="Proline"&&Amino2=="Cysteine")||(Amino2=="Proline"&&Amino1=="Cysteine")){
    return(-3)
  }
  if((Amino1=="Serine"&&Amino2=="Cysteine")||(Amino2=="Serine"&&Amino1=="Cysteine")){
    return(0)
  }
  if((Amino1=="Threonine"&&Amino2=="Cysteine")||(Amino2=="Threonine"&&Amino1=="Cysteine")){
    return(-2)
  }
  if((Amino1=="Tryptophan"&&Amino2=="Cysteine")||(Amino2=="Tryptophan"&&Amino1=="Cysteine")){
    return(-8)
  }
  if((Amino1=="Tyrosine"&&Amino2=="Cysteine")||(Amino2=="Tyrosine"&&Amino1=="Cysteine")){
    return(0)
  }
  if((Amino1=="Valine"&&Amino2=="Cysteine")||(Amino2=="Valine"&&Amino1=="Cysteine")){
    return(-2)
  }
  #Glu Block
  if((Amino1=="Glutamine"&&Amino2=="Glutamine")||(Amino2=="Glutamine"&&Amino1=="Glutamine")){
    return(4)
  }
  if((Amino1=="Glutamic Acid"&&Amino2=="Glutamine")||(Amino2=="Glutamic Acid"&&Amino1=="Glutamine")){
    return(2)
  }
  if((Amino1=="Glycine"&&Amino2=="Glutamine")||(Amino2=="Glycine"&&Amino1=="Glutamine")){
    return(-1)
  }
  if((Amino1=="Histidine"&&Amino2=="Glutamine")||(Amino2=="Histidine"&&Amino1=="Glutamine")){
    return(3)
  }
  if((Amino1=="Isoleucine"&&Amino2=="Glutamine")||(Amino2=="Isoleucine"&&Amino1=="Glutamine")){
    return(-2)
  }
  if((Amino1=="Leucine"&&Amino2=="Glutamine")||(Amino2=="Leucine"&&Amino1=="Glutamine")){
    return(-2)
  }
  if((Amino1=="Lysine"&&Amino2=="Glutamine")||(Amino2=="Lysine"&&Amino1=="Glutamine")){
    return(1)
  }
  if((Amino1=="Methionine"&&Amino2=="Glutamine")||(Amino2=="Methionine"&&Amino1=="Glutamine")){
    return(-1)
  }
  if((Amino1=="Phenylalanine"&&Amino2=="Glutamine")||(Amino2=="Phenylalanine"&&Amino1=="Glutamine")){
    return(-5)
  }
  if((Amino1=="Proline"&&Amino2=="Glutamine")||(Amino2=="Proline"&&Amino1=="Glutamine")){
    return(0)
  }
  if((Amino1=="Serine"&&Amino2=="Glutamine")||(Amino2=="Serine"&&Amino1=="Glutamine")){
    return(-1)
  }
  if((Amino1=="Threonine"&&Amino2=="Glutamine")||(Amino2=="Threonine"&&Amino1=="Glutamine")){
    return(-1)
  }
  if((Amino1=="Tryptophan"&&Amino2=="Glutamine")||(Amino2=="Tryptophan"&&Amino1=="Glutamine")){
    return(-5)
  }
  if((Amino1=="Tyrosine"&&Amino2=="Glutamine")||(Amino2=="Tyrosine"&&Amino1=="Glutamine")){
    return(-4)
  }
  if((Amino1=="Valine"&&Amino2=="Glutamine")||(Amino2=="Valine"&&Amino1=="Glutamine")){
    return(-2)
  }
  
  #Glutamic Acid Block
  if((Amino1=="Glutamic Acid"&&Amino2=="Glutamic Acid")||(Amino2=="Glutamic Acid"&&Amino1=="Glutamic Acid")){
    return(4)
  }
  if((Amino1=="Glycine"&&Amino2=="Glutamic Acid")||(Amino2=="Glycine"&&Amino1=="Glutamic Acid")){
    return(0)
  }
  if((Amino1=="Histidine"&&Amino2=="Glutamic Acid")||(Amino2=="Histidine"&&Amino1=="Glutamic Acid")){
    return(1)
  }
  if((Amino1=="Isoleucine"&&Amino2=="Glutamic Acid")||(Amino2=="Isoleucine"&&Amino1=="Glutamic Acid")){
    return(-2)
  }
  if((Amino1=="Leucine"&&Amino2=="Glutamic Acid")||(Amino2=="Leucine"&&Amino1=="Glutamic Acid")){
    return(-3)
  }
  if((Amino1=="Lysine"&&Amino2=="Glutamic Acid")||(Amino2=="Lysine"&&Amino1=="Glutamic Acid")){
    return(0)
  }
  if((Amino1=="Methionine"&&Amino2=="Glutamic Acid")||(Amino2=="Methionine"&&Amino1=="Glutamic Acid")){
    return(-2)
  }
  if((Amino1=="Phenylalanine"&&Amino2=="Glutamic Acid")||(Amino2=="Phenylalanine"&&Amino1=="Glutamic Acid")){
    return(-5)
  }
  if((Amino1=="Proline"&&Amino2=="Glutamic Acid")||(Amino2=="Proline"&&Amino1=="Glutamic Acid")){
    return(-1)
  }
  if((Amino1=="Serine"&&Amino2=="Glutamic Acid")||(Amino2=="Serine"&&Amino1=="Glutamic Acid")){
    return(0)
  }
  if((Amino1=="Threonine"&&Amino2=="Glutamic Acid")||(Amino2=="Threonine"&&Amino1=="Glutamic Acid")){
    return(0)
  }
  if((Amino1=="Tryptophan"&&Amino2=="Glutamic Acid")||(Amino2=="Tryptophan"&&Amino1=="Glutamic Acid")){
    return(-7)
  }
  if((Amino1=="Tyrosine"&&Amino2=="Glutamic Acid")||(Amino2=="Tyrosine"&&Amino1=="Glutamic Acid")){
    return(-4)
  }
  if((Amino1=="Valine"&&Amino2=="Glutamic Acid")||(Amino2=="Valine"&&Amino1=="Glutamic Acid")){
    return(-2)
  }
  #Glycine Block
  if((Amino1=="Glycine"&&Amino2=="Glycine")||(Amino2=="Glycine"&&Amino1=="Glycine")){
    return(5)
  }
  if((Amino1=="Histidine"&&Amino2=="Glycine")||(Amino2=="Histidine"&&Amino1=="Glycine")){
    return(-2)
  }
  if((Amino1=="Isoleucine"&&Amino2=="Glycine")||(Amino2=="Isoleucine"&&Amino1=="Glycine")){
    return(-3)
  }
  if((Amino1=="Leucine"&&Amino2=="Glycine")||(Amino2=="Leucine"&&Amino1=="Glycine")){
    return(-4)
  }
  if((Amino1=="Lysine"&&Amino2=="Glycine")||(Amino2=="Lysine"&&Amino1=="Glycine")){
    return(-2)
  }
  if((Amino1=="Methionine"&&Amino2=="Glycine")||(Amino2=="Methionine"&&Amino1=="Glycine")){
    return(-3)
  }
  if((Amino1=="Phenylalanine"&&Amino2=="Glycine")||(Amino2=="Phenylalanine"&&Amino1=="Glycine")){
    return(-5)
  }
  if((Amino1=="Proline"&&Amino2=="Glycine")||(Amino2=="Proline"&&Amino1=="Glycine")){
    return(-1)
  }
  if((Amino1=="Serine"&&Amino2=="Glycine")||(Amino2=="Serine"&&Amino1=="Glycine")){
    return(1)
  }
  if((Amino1=="Threonine"&&Amino2=="Glycine")||(Amino2=="Threonine"&&Amino1=="Glycine")){
    return(0)
  }
  if((Amino1=="Tryptophan"&&Amino2=="Glycine")||(Amino2=="Tryptophan"&&Amino1=="Glycine")){
    return(-7)
  }
  if((Amino1=="Tyrosine"&&Amino2=="Glycine")||(Amino2=="Tyrosine"&&Amino1=="Glycine")){
    return(-5)
  }
  if((Amino1=="Valine"&&Amino2=="Glycine")||(Amino2=="Valine"&&Amino1=="Glycine")){
    return(-1)
  }
  #Histidine Block
  if((Amino1=="Histidine"&&Amino2=="Histidine")||(Amino2=="Histidine"&&Amino1=="Histidine")){
    return(6)
  }
  if((Amino1=="Isoleucine"&&Amino2=="Histidine")||(Amino2=="Isoleucine"&&Amino1=="Histidine")){
    return(-2)
  }
  if((Amino1=="Leucine"&&Amino2=="Histidine")||(Amino2=="Leucine"&&Amino1=="Histidine")){
    return(-2)
  }
  if((Amino1=="Lysine"&&Amino2=="Histidine")||(Amino2=="Lysine"&&Amino1=="Histidine")){
    return(0)
  }
  if((Amino1=="Methionine"&&Amino2=="Histidine")||(Amino2=="Methionine"&&Amino1=="Histidine")){
    return(-2)
  }
  if((Amino1=="Phenylalanine"&&Amino2=="Histidine")||(Amino2=="Phenylalanine"&&Amino1=="Histidine")){
    return(-2)
  }
  if((Amino1=="Proline"&&Amino2=="Histidine")||(Amino2=="Proline"&&Amino1=="Histidine")){
    return(0)
  }
  if((Amino1=="Serine"&&Amino2=="Histidine")||(Amino2=="Serine"&&Amino1=="Histidine")){
    return(1)
  }
  if((Amino1=="Threonine"&&Amino2=="Histidine")||(Amino2=="Threonine"&&Amino1=="Histidine")){
    return(-1)
  }
  if((Amino1=="Tryptophan"&&Amino2=="Histidine")||(Amino2=="Tryptophan"&&Amino1=="Histidine")){
    return(-3)
  }
  if((Amino1=="Tyrosine"&&Amino2=="Histidine")||(Amino2=="Tyrosine"&&Amino1=="Histidine")){
    return(0)
  }
  if((Amino1=="Valine"&&Amino2=="Histidine")||(Amino2=="Valine"&&Amino1=="Histidine")){
    return(-2)
  }
  
  #Isoleucine
  if((Amino1=="Isoleucine"&&Amino2=="Isoleucine")||(Amino2=="Isoleucine"&&Amino1=="Isoleucine")){
    return(5)
  }
  if((Amino1=="Leucine"&&Amino2=="Isoleucine")||(Amino2=="Leucine"&&Amino1=="Isoleucine")){
    return(2)
  }
  if((Amino1=="Lysine"&&Amino2=="Isoleucine")||(Amino2=="Lysine"&&Amino1=="Isoleucine")){
    return(-2)
  }
  if((Amino1=="Methionine"&&Amino2=="Isoleucine")||(Amino2=="Methionine"&&Amino1=="Isoleucine")){
    return(2)
  }
  if((Amino1=="Phenylalanine"&&Amino2=="Isoleucine")||(Amino2=="Phenylalanine"&&Amino1=="Isoleucine")){
    return(1)
  }
  if((Amino1=="Proline"&&Amino2=="Isoleucine")||(Amino2=="Proline"&&Amino1=="Isoleucine")){
    return(-2)
  }
  if((Amino1=="Serine"&&Amino2=="Isoleucine")||(Amino2=="Serine"&&Amino1=="Isoleucine")){
    return(-1)
  }
  if((Amino1=="Threonine"&&Amino2=="Isoleucine")||(Amino2=="Threonine"&&Amino1=="Isoleucine")){
    return(0)
  }
  if((Amino1=="Tryptophan"&&Amino2=="Isoleucine")||(Amino2=="Tryptophan"&&Amino1=="Isoleucine")){
    return(-5)
  }
  if((Amino1=="Tyrosine"&&Amino2=="Isoleucine")||(Amino2=="Tyrosine"&&Amino1=="Isoleucine")){
    return(-1)
  }
  if((Amino1=="Valine"&&Amino2=="Isoleucine")||(Amino2=="Valine"&&Amino1=="Isoleucine")){
    return(4)
  }
  #Leucine Block
  if((Amino1=="Leucine"&&Amino2=="Leucine")||(Amino2=="Leucine"&&Amino1=="Leucine")){
    return(6)
  }
  if((Amino1=="Lysine"&&Amino2=="Leucine")||(Amino2=="Lysine"&&Amino1=="Leucine")){
    return(-3)
  }
  if((Amino1=="Methionine"&&Amino2=="Leucine")||(Amino2=="Methionine"&&Amino1=="Leucine")){
    return(4)
  }
  if((Amino1=="Phenylalanine"&&Amino2=="Leucine")||(Amino2=="Phenylalanine"&&Amino1=="Leucine")){
    return(2)
  }
  if((Amino1=="Proline"&&Amino2=="Leucine")||(Amino2=="Proline"&&Amino1=="Leucine")){
    return(-3)
  }
  if((Amino1=="Serine"&&Amino2=="Leucine")||(Amino2=="Serine"&&Amino1=="Leucine")){
    return(-3)
  }
  if((Amino1=="Threonine"&&Amino2=="Leucine")||(Amino2=="Threonine"&&Amino1=="Leucine")){
    return(-2)
  }
  if((Amino1=="Tryptophan"&&Amino2=="Leucine")||(Amino2=="Tryptophan"&&Amino1=="Leucine")){
    return(-2)
  }
  if((Amino1=="Tyrosine"&&Amino2=="Leucine")||(Amino2=="Tyrosine"&&Amino1=="Leucine")){
    return(-1)
  }
  if((Amino1=="Valine"&&Amino2=="Leucine")||(Amino2=="Valine"&&Amino1=="Leucine")){
    return(2)
  }
  #Lysine
  if((Amino1=="Lysine"&&Amino2=="Lysine")||(Amino2=="Lysine"&&Amino1=="Lysine")){
    return(5)
  }
  if((Amino1=="Methionine"&&Amino2=="Lysine")||(Amino2=="Methionine"&&Amino1=="Lysine")){
    return(0)
  }
  if((Amino1=="Phenylalanine"&&Amino2=="Lysine")||(Amino2=="Phenylalanine"&&Amino1=="Lysine")){
    return(-5)
  }
  if((Amino1=="Proline"&&Amino2=="Lysine")||(Amino2=="Proline"&&Amino1=="Lysine")){
    return(-1)
  }
  if((Amino1=="Serine"&&Amino2=="Lysine")||(Amino2=="Serine"&&Amino1=="Lysine")){
    return(0)
  }
  if((Amino1=="Threonine"&&Amino2=="Lysine")||(Amino2=="Threonine"&&Amino1=="Lysine")){
    return(0)
  }
  if((Amino1=="Tryptophan"&&Amino2=="Lysine")||(Amino2=="Tryptophan"&&Amino1=="Lysine")){
    return(-3)
  }
  if((Amino1=="Tyrosine"&&Amino2=="Lysine")||(Amino2=="Tyrosine"&&Amino1=="Lysine")){
    return(-4)
  }
  if((Amino1=="Valine"&&Amino2=="Lysine")||(Amino2=="Valine"&&Amino1=="Lysine")){
    return(-2)
  }
  
  #methionine
  if((Amino1=="Methionine"&&Amino2=="Methionine")||(Amino2=="Methionine"&&Amino1=="Methionine")){
    return(6)
  }
  if((Amino1=="Phenylalanine"&&Amino2=="Methionine")||(Amino2=="Phenylalanine"&&Amino1=="Methionine")){
    return(0)
  }
  if((Amino1=="Proline"&&Amino2=="Methionine")||(Amino2=="Proline"&&Amino1=="Methionine")){
    return(-2)
  }
  if((Amino1=="Serine"&&Amino2=="Methionine")||(Amino2=="Serine"&&Amino1=="Methionine")){
    return(-2)
  }
  if((Amino1=="Threonine"&&Amino2=="Methionine")||(Amino2=="Threonine"&&Amino1=="Methionine")){
    return(-1)
  }
  if((Amino1=="Tryptophan"&&Amino2=="Methionine")||(Amino2=="Tryptophan"&&Amino1=="Methionine")){
    return(-4)
  }
  if((Amino1=="Tyrosine"&&Amino2=="Methionine")||(Amino2=="Tyrosine"&&Amino1=="Methionine")){
    return(-2)
  }
  if((Amino1=="Valine"&&Amino2=="Methionine")||(Amino2=="Valine"&&Amino1=="Methionine")){
    return(2)
  }
  #phenylalanine
  if((Amino1=="Phenylalanine"&&Amino2=="Phenylalanine")||(Amino2=="Phenylalanine"&&Amino1=="Phenylalanine")){
    return(9)
  }
  if((Amino1=="Proline"&&Amino2=="Phenylalanine")||(Amino2=="Proline"&&Amino1=="Phenylalanine")){
    return(-5)
  }
  if((Amino1=="Serine"&&Amino2=="Phenylalanine")||(Amino2=="Serine"&&Amino1=="Phenylalanine")){
    return(-3)
  }
  if((Amino1=="Threonine"&&Amino2=="Phenylalanine")||(Amino2=="Threonine"&&Amino1=="Phenylalanine")){
    return(-2)
  }
  if((Amino1=="Tryptophan"&&Amino2=="Phenylalanine")||(Amino2=="Tryptophan"&&Amino1=="Phenylalanine")){
    return(0)
  }
  if((Amino1=="Tyrosine"&&Amino2=="Phenylalanine")||(Amino2=="Tyrosine"&&Amino1=="Phenylalanine")){
    return(7)
  }
  if((Amino1=="Valine"&&Amino2=="Phenylalanine")||(Amino2=="Valine"&&Amino1=="Phenylalanine")){
    return(-1)
  }
  #Proline
  if((Amino1=="Proline"&&Amino2=="Proline")||(Amino2=="Proline"&&Amino1=="Proline")){
    return(6)
  }
  if((Amino1=="Serine"&&Amino2=="Proline")||(Amino2=="Serine"&&Amino1=="Proline")){
    return(1)
  }
  if((Amino1=="Threonine"&&Amino2=="Proline")||(Amino2=="Threonine"&&Amino1=="Proline")){
    return(1)
  }
  if((Amino1=="Tryptophan"&&Amino2=="Proline")||(Amino2=="Tryptophan"&&Amino1=="Proline")){
    return(-6)
  }
  if((Amino1=="Tyrosine"&&Amino2=="Proline")||(Amino2=="Tyrosine"&&Amino1=="Proline")){
    return(-5)
  }
  if((Amino1=="Valine"&&Amino2=="Proline")||(Amino2=="Valine"&&Amino1=="Proline")){
    return(-1)
  }
  #Serine
  if((Amino1=="Serine"&&Amino2=="Serine")||(Amino2=="Serine"&&Amino1=="Serine")){
    return(3)
  }
  if((Amino1=="Threonine"&&Amino2=="Serine")||(Amino2=="Threonine"&&Amino1=="Serine")){
    return(1)
  }
  if((Amino1=="Tryptophan"&&Amino2=="Serine")||(Amino2=="Tryptophan"&&Amino1=="Serine")){
    return(-2)
  }
  if((Amino1=="Tyrosine"&&Amino2=="Serine")||(Amino2=="Tyrosine"&&Amino1=="Serine")){
    return(-3)
  }
  if((Amino1=="Valine"&&Amino2=="Serine")||(Amino2=="Valine"&&Amino1=="Serine")){
    return(-1)
  }
  #Threonine
  if((Amino1=="Threonine"&&Amino2=="Threonine")||(Amino2=="Threonine"&&Amino1=="Threonine")){
    return(3)
  }
  if((Amino1=="Tryptophan"&&Amino2=="Threonine")||(Amino2=="Tryptophan"&&Amino1=="Threonine")){
    return(-5)
  }
  if((Amino1=="Tyrosine"&&Amino2=="Threonine")||(Amino2=="Tyrosine"&&Amino1=="Threonine")){
    return(-3)
  }
  if((Amino1=="Valine"&&Amino2=="Threonine")||(Amino2=="Valine"&&Amino1=="Threonine")){
    return(0)
  }
  #Tryptophan
  if((Amino1=="Tryptophan"&&Amino2=="Tryptophan")||(Amino2=="Tryptophan"&&Amino1=="Tryptophan")){
    return(17)
  }
  if((Amino1=="Tyrosine"&&Amino2=="Tryptophan")||(Amino2=="Tyrosine"&&Amino1=="Tryptophan")){
    return(0)
  }
  if((Amino1=="Valine"&&Amino2=="Tryptophan")||(Amino2=="Valine"&&Amino1=="Tryptophan")){
    return(-6)
  }
  #tyrosine
  if((Amino1=="Tyrosine"&&Amino2=="Tyrosine")||(Amino2=="Tyrosine"&&Amino1=="Tyrosine")){
    return(10)
  }
  if((Amino1=="Valine"&&Amino2=="Tyrosine")||(Amino2=="Valine"&&Amino1=="Tyrosine")){
    return(-2)
  }
  #valine
  if((Amino1=="Valine"&&Amino2=="Valine")||(Amino2=="Valine"&&Amino1=="Valine")){
    return(4)
  }
  else{
    print(Amino1)
    print(Amino2)
    return(-20)
  }
}
#Variable Declaration Block
Genes <- read.csv("RVersionCBSalGenes.csv", stringsAsFactors = FALSE)
CitroBacDNA<-read_file("CitroBacKPureDNA.txt")
CitroBacDNA<-gsub("\n","",CitroBacDNA)
SalDNA <- read_file("SalmonellaPureDNA.txt")
SalDNA<-gsub("\n","",SalDNA)
Genes <- data.frame(Genes$gapnum, Genes$cbgeneseq, Genes$salgeneseq, stringsAsFactors = FALSE)
G1<-as.vector(Genes$Genes.gapnum)
G2<-data.frame(Genes$Genes.cbgeneseq)
G3<-data.frame(Genes$Genes.salgeneseq)
length(G3)
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



#CLUSAL Run on Genes
while (n <= nrow(G2)){
#  while(n<=4){ #testing line when not running full version of code
  testnum<-G1[n]
  Gene1Test <- toString(G2[n,1])
  Gene2Test<-toString(G3[n,1])
  Gene1RunTitle <- paste(">","ECGap_number","k", sep = "") #I don't feel like rewriting this to hold generality, so these need to stay as is to keep the parsing working correctly
  Gene2RunTitle <- paste(">","SalGap_number","k", sep = "")
  
  #write a text file that is in proper FASTA format
  teststring <-"FastaInCBSal.txt"
  fileConn<-file(teststring) #reset the file to being blank text again
  writeLines("\n",fileConn)
  close(fileConn)
  sink(teststring)
  cat(Gene1RunTitle)
  cat("\n")
  cat(Gene1Test)
  cat("\n")
  cat(Gene2RunTitle)
  cat("\n")
  cat(Gene2Test)
  sink()
  if(is.null(Gene1Test) || is.null(Gene2Test)){
    file.remove("FastaInCBSal.txt")
    sink(teststring)
    cat(Gene1RunTitle)
    cat("\n")
    cat("AAA")
    cat("\n")
    cat(Gene2RunTitle)
    cat("\n")
    cat("AAA")
    sink()
  }
  #Call the clustalw function
  print("Before")
  systemstring <- "clustalw2 -infile=FastaInCBSal.txt -type=DNA"
  system(systemstring)
  #pull the gapped strings out of the .aln file
  alnlines<-readLines("FastaInCBSal.aln")
  print("After")
  print(Gene1RunTitle)
  print(Gene2RunTitle)
  alnlines1<-grep(substring(Gene1RunTitle,2), alnlines, value=TRUE)
  alnlines2<-grep(substring(Gene2RunTitle,2), alnlines, value=TRUE)
  alnlines1<-gsub(substring(Gene1RunTitle,2), "", alnlines1)
  alnlines1<-gsub(" ", "", alnlines1)
  alnlines2<-gsub(substring(Gene2RunTitle,2), "", alnlines2)
  alnlines2<-gsub(" ", "", alnlines2)
  alnlines1<-paste(alnlines1, collapse='')
  alnlines2<-paste(alnlines2, collapse='')
  print(alnlines2)
  Gene1gapstring <- alnlines1
  Gene2gapstring <- alnlines2
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
#      print("ASDF")
      Flag2=FALSE
      PostCounter=0
      if(Flag1==FALSE){
      PreCounter=0
      }
    }
    if(Gene1workingcharacter == Gene2workingcharacter && Flag2==TRUE)
    {
      PostCounter = PostCounter +1
      #print(PostCounter) causes too much lag
    }
    if(Gene1workingcharacter != Gene2workingcharacter && Flag1==TRUE)
    {
      Flag1 = FALSE
      Flag2 = TRUE
      Gene1break = Gene1workingcharacter
      Gene2break = Gene2workingcharacter
    }
    if(Gene1workingcharacter == Gene2workingcharacter && Flag1==TRUE)
    {
      PreCounter = PreCounter+1
    }
    if(Gene1workingcharacter== Gene2workingcharacter && Flag1==FALSE)
    {
      Flag1 = TRUE
      PreCounter=PreCounter+1
    }
    
    
    
  }
  #end of intra-gene break identification
}
print(names(DF))
DF$PostCount <- as.integer(DF$PostCount)
DF$PreCount <- as.integer((DF$PreCount))
ModDF = DF #DF is the entries from tallying up the mismatches
ModDF = subset(ModDF, PostCount > 0) #clears up things where an error got through
PlotData = subset(ModDF, PostCount>=3 & PreCount>=3) #this is the arbitry definition of what we considered conserved, this can change
hist(as.numeric(PlotData$PreCount)) #plotting up everything
hist(as.numeric(PlotData$PostCount))
hist(as.numeric(ModDF$PreCount))
hist(as.numeric(ModDF$PostCount))
plot(jitter(as.numeric(DF$PreCount)),jitter(as.numeric(DF$PostCount)))
plot(as.numeric(PlotData$PreCount),as.numeric((PlotData$PostCount)))

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
names(Gene_vector_of_lengths)=Gene_vector_of_names
write.csv(Gene_vector_of_lengths,file= "CBSalGeneData.csv")
pdf('CBSalplot.pdf')
barplot(Gene_vector_of_lengths,names.arg = Gene_vector_of_names) #constructing barplot of mutation frequencies
dev.off()

