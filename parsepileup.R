library(tibble)
library(readr)
library(dplyr)
library(plotrix)
#biocLite('Rsamtools') 

#pileup <- readPileup('8802/pileup.pileup')

setwd('C:/Users/Michael/Documents/Variance')
data <- read_tsv('8802/pileup.pileup',col_names = FALSE)
colnames(data) <- c('name','position','bp','depth','info','quality')
bp_info <- matrix(, nrow = 0, ncol = 4)
for( i in seq(1,nrow(data))){
  #get the current base
  base = data$bp[i]
  string <- strsplit(data$info[i], "")[[1]]
  As <- 0
  Cs <- 0
  Ts <- 0
  Gs <- 0
  for( char in string){
    if(char %in% c('A','a') | data$bp[i] %in% c('A','a')) {
      As <- As + 1
    }else if(char %in% c('C','c') | (data$bp[i] %in% c('C','c')  & char %in% c('.',',','') )){
      Cs <- Cs + 1
    }else if(char %in% c('T','t') | (data$bp[i] %in% c('T','t')  & char %in% c('.',',','') )){
      Ts <- Ts + 1
    }else if(char %in% c('G','g') | (data$bp[i] %in% c('G','g')  & char %in% c('.',',','') )  ){
      Gs <- Gs + 1
    }
  }
  bp_info <- rbind(bp_info, c(As,Cs,Gs,Ts))
}