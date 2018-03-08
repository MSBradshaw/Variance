library(tibble)
library(readr)
library(dplyr)
library(plotrix)
#biocLite('Rsamtools') 
library(Rsamtools)

#pileup <- readPileup('8802/pileup.pileup')

setwd('C:/Users/Michael/Documents/Variance')
data <- read_tsv('8802/pileup.pileup')
colnames(data) <- c('name','position','bp','depth','info','quality')
for( i in seq(1,nrow(data))){
  #get the current base
  base = data$bp[i]
  string <- strsplit(data$info[i], "")[[1]]
  As <- c()
  Cs <- c()
  Ts <- c()
  Gs <- c()
  for( char in string){
    if(char %in% c('A','a')) {
      
    }else if(char %in% c('C','c')){
      
    }else if(char %in% c('T','t')){
      
    }else if(char %in% c('G','g')){
      
    }else if(char %in% c(.))
  }
}