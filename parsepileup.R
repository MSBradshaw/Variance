library(tibble)
library(readr)
library(dplyr)
library(plotrix)
library(stringr)

#Thymine	6	YELLOW
#Adenine	6	BLUE
#Cytosine	6	RED
#Guanine	6	GREEN

#pileup <- readPileup('8802/pileup.pileup')

setwd('C:/Users/Michael/Documents/Variance')
data <- read_tsv('8802/pileup.pileup',col_names = FALSE)
data2 <- read_tsv('final_dirs/8805/pileup.pileup',col_names = FALSE)
data3 <- read_tsv('final_dirs/8800/pileup.pileup',col_names = FALSE)
data4 <- read_tsv('final_dirs/8810/pileup.pileup',col_names = FALSE)

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

#number of bps to be removed from the front to of the read
front_trim <- 600

#number of bps to be removed from the end to of the read
back_trim <- 600

#generate the sequence string with bp coloring or variance coloring
create_conscensus_svg_string <- function(bp_info,conscensus,variance){
  extra_class <- ''
  variance_string <- ''
  if(variance){
    extra_class <- ' variance '
  }
  x_count <- 10
  cons_vec <- strsplit(conscensus,'')

  output <- paste('
  <svg class="svg_item ', extra_class ,'">',sep='') 
  
  for( i in seq(1,nrow(bp_info))){
    class <- 'match'
    if(!bp_info[i,1] %in% cons_vec[i]){
      class <- bp_info[i,1]
    }
    
    variance <- 1 - (max(bp_info[i,]) / sum(bp_info[i,]))
    
    if(variance == 0){
      class <- paste(class, 'variance_none')
    }else if(variance < .01){
      #illumina error rate
      class <- paste(class, 'variance_01')
    }else if(variance < .10){
      class <- paste(class, 'variance_10')
    }else{
      class <- paste(class, 'variance_high')
    }
    
    line <- paste('<rect class="base ', class ,' " x="', x_count ,'" y="10" width="10" height="20" stroke="black"/>\n')
    output <- paste(output, line)
    x_count <- x_count + 10
  }
}  

#generate the sequence string with bp coloring or variance coloring
create_svg_string <- function(bp_info,front_trim,back_trim,variance,y_count){
  variance_string <- ''
  if(variance){
    variance_string <- ' variance '
  }
  x_count <- 10
  output <- paste('
  <svg class="svg_item ', variance_string ,'">',sep='') #this is inside a paste function so I can define the width of the svg

  for( i in seq(1,nrow(bp_info))){
    #get the most common bp
    index <- which.max(bp_info[i,])
    color <- c('A','C','G','T')
    variance <- 1 - (max(bp_info[i,]) / sum(bp_info[i,]))
    extra_class <- color[index]
    
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
    
    line <- paste('<rect class="base ', extra_class ,' " x="', x_count ,'" y="', y_count ,'" width="10" height="20" stroke="black"/>\n')
    output <- paste(output, line)
    x_count <- x_count + 10
  }  

  output <- paste(output, '</svg>')
  return(output)
}

#creates the style tag
#bp_rows is the number of rows in the bp_info matrix
get_style_string <- function(bp_rows){
style <- paste('  <style>
    .A{
      fill: blue;
    }
  .C{
    fill: red;
  }
  .G{
    fill:green;
  }
  .T{
    fill:yellow;
  }
  base:hover{
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
  .match{
  fill: #E6E6E6;
  }
  .svg_item {
    width:',(bp_rows * 10),'px;
    height: 30px;
  }
  </style>',sep='')
  return(style)
}

#wrap a given string in an html skeleton so it can be viewed in the browser
htlm_wrap <- function(string){
  paste('<!DOCTYPE html>\n<html>\n<head>\n<meta charset="UTF-8">\n<title>title</title>\n</head>\n<body>',string,'</body>\n</html>')
}

#given a vector or strings, get the consesous sequence
get_consensus <- function(strings){
  seq <- ''
  range <- seq(1,length(strings))
  first_string <- strsplit(strings[1],'')[[1]]
  for( i in seq(1:length(first_string))){
    bps <- str_sub(strings,i,i)
    seq <- paste(seq,names(sort(table(bps),decreasing = TRUE))[1],sep='')
  }
  return(seq)
}

#trim the data
data_small <- data[front_trim:(nrow(data) - back_trim),]
bp_data <- get_bp_info(data_small)

style <- get_style_string(nrow(bp_data))
svg_item <- create_svg_string(bp_data,front_trim,back_trim,TRUE,10)
svg_item2 <- create_svg_string(bp_data,front_trim,back_trim,FALSE,10)
final_string <- paste(style,svg_item,svg_item2)
html <- htlm_wrap(final_string)

fileConn<-file("output.html")
writeLines(html, fileConn)
close(fileConn)



