library(tibble)
library(readr)
library(dplyr)
library(plotrix)
#biocLite('Rsamtools') 

#Thymine	6	YELLOW
#Adenine	6	BLUE
#Cytosine	6	RED
#Guanine	6	GREEN

#pileup <- readPileup('8802/pileup.pileup')

setwd('C:/Users/Michael/Documents/Variance')
data <- read_tsv('8802/pileup.pileup',col_names = FALSE)
colnames(data) <- c('name','position','bp','depth','info','quality')

#return a matrix of the bp
get_bp_info <- function(data){
  bp_info <- matrix(, nrow = 0, ncol = 4)
  for( i in seq(1,nrow(data))){
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
  return(bp_info)
}

#generate the sequence
create_svg_string <- function(bp_info){
  x_count <- 10
  output <- '<!DOCTYPE html>\n<html>\n<head>\n<meta charset="UTF-8">\n<title>title</title>\n</head>\n
  <body>
  <style>
  .svg_item{
    width: 100%;
  }
  .hover_me:hover{
  fill: green;
  }
.variance .variance_none{
  fill: black;
}
.variance .variance_01{
  fill: grey;
}
.variance .variance_10{
  fill: #E6E6E6;
}
.variance .variance_high{
  fill: #FFF;
}
  </style><svg class="svg_item">'
  for( i in seq(1,nrow(bp_info))){
    #get the most common bp
    index <- which.max(bp_info[i,])
    color <- c('blue','red','green','yellow')
    variance <- 1 - (max(bp_info[i,]) / sum(bp_info[i,]))
    extra_class <- ''
    
    if(variance == 0){
      extra_class <- paste(extra_class, 'variance_none')
    }else if(variance < .01){
      #illumina error rate
      extra_class <- paste(extra_class, 'variance_01')
    }else if(variance < .10){
      extra_class <- paste(extra_class, 'variance_10')
    }else{
      extra_class <- paste(extra_class, 'variance_high')
    }
    
    line <- paste('<rect class="hover_me ', extra_class ,' " x="', x_count ,'" y="10" width="10" height="20" stroke="black" fill= "', color[index], '"/>\n')
    output <- paste(output, line)
    x_count <- x_count + 10
  }  
  output <- paste(output, '</svg>\n</body>\n</html>')
  return(output)
}

string <- get_bp_info(data)
svg_item <- create_svg_string(string)

fileConn<-file("output.html")
writeLines(svg_item, fileConn)
close(fileConn)



