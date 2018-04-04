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

data <- read_tsv('final_dirs/8664_2/pileup.pileup',col_names = FALSE)
data2 <- read_tsv('final_dirs/8664_3/pileup.pileup',col_names = FALSE)
data3 <- read_tsv('final_dirs/8664_4/pileup.pileup',col_names = FALSE)
data4 <- read_tsv('final_dirs/8664_5/pileup.pileup',col_names = FALSE)
data5 <- read_tsv('final_dirs/8664_6/pileup.pileup',col_names = FALSE)

colnames(data) <- c('name','position','bp','depth','info','quality')
colnames(data2) <- c('name','position','bp','depth','info','quality')
colnames(data3) <- c('name','position','bp','depth','info','quality')
colnames(data4) <- c('name','position','bp','depth','info','quality')
colnames(data5) <- c('name','position','bp','depth','info','quality')


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
  color <- c('A','C','G','T')
  
  for( i in seq(1,nrow(bp_info))){
    class <- 'match'
    index <- which.max(bp_info[i,])
    letter <- color[index]
    if(!letter %in% cons_vec[[1]][i]){
      class <- letter
    }
    
    variance <- 1 - (max(as.numeric(bp_info[i,])) / sum(as.numeric(bp_info[i,])))
    
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
  return(output)
}  

#generate the sequence string with bp coloring or variance coloring
create_svg_string <- function(bp_info,front_trim,back_trim,variance,y_count,bright){
  variance_string <- ''
  if(variance){
    variance_string <- ' variance '
  }
  if(bright){
    variance_string <- paste(variance_string,' bright ')
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
    fill: #333333;
  }
  .variance .variance_10{
    fill: #999999;
  }
  .variance .variance_high{
    fill: #FFF;
  }
  .bright .variance_none{
    fill: black;
}
.bright .variance_01{
fill: #662900;
}
.bright .variance_10{
fill: #e65c00;
}
.bright .variance_high{
fill: #ffe0cc;
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
  lengths <- nchar(strings)
  index <- which.min(lengths)
  first_string <- strsplit(strings[index],'')[[1]]
  for( i in seq(1:length(first_string))){
    bps <- str_sub(strings,i,i)
    seq <- paste(seq,names(sort(table(bps),decreasing = TRUE))[1],sep='')
  }
  return(seq)
}

get_bp_string <- function(bp_info){
  bases <- apply(bp_info,1,function(row){
    bases <- c('A','C','G','T')
    index <- which.max(row)
    return(bases[index])
  })
  final <- paste(bases,collapse='')
  return(final)
}
#trim the data
data_small <- data[front_trim:(nrow(data) - back_trim),]
data_small2 <- data2[front_trim:(nrow(data2) - back_trim),]
data_small3 <- data3[front_trim:(nrow(data3) - back_trim),]
data_small4 <- data4[front_trim:(nrow(data4) - back_trim),]
data_small5 <- data5[front_trim:(nrow(data5) - back_trim),]


bp_data <- get_bp_info(data_small)
bp_data2 <- get_bp_info(data_small2)
bp_data3 <- get_bp_info(data_small3)
bp_data4 <- get_bp_info(data_small4)
bp_data5 <- get_bp_info(data_small5)


style <- get_style_string(nrow(bp_data))
svg_item <- create_svg_string(bp_data,front_trim,back_trim,TRUE,10,TRUE)
svg_item2 <- create_svg_string(bp_data2,front_trim,back_trim,TRUE,10,TRUE)
svg_item3 <- create_svg_string(bp_data3,front_trim,back_trim,TRUE,10,TRUE)
svg_item4 <- create_svg_string(bp_data4,front_trim,back_trim,TRUE,10,TRUE)
svg_item5 <- create_svg_string(bp_data5,front_trim,back_trim,TRUE,10,TRUE)


strings <- c(get_bp_string(bp_data),get_bp_string(bp_data2),get_bp_string(bp_data3),get_bp_string(bp_data4))
con <- get_consensus(strings)

#final_string <- create_conscensus_svg_string(bp_data,con,FALSE)

final_string <- paste(style,svg_item,svg_item2,svg_item3,svg_item4,svg_item5)

#final_string <- paste(style,final_string)
html <- htlm_wrap(final_string)

fileConn<-file("output.html")
writeLines(html, fileConn)
close(fileConn)

build_image <- function(files,image_option,color_option,front_trim,back_trim){
  #read in the files
  read_tsv('final_dirs/8664_2/pileup.pileup',col_names = FALSE)
  files <- as.tibble(files)
  file_infos <- apply(files,1,function(path){
    data <- read_tsv(path,col_names = FALSE)
    colnames(data) <- c('name','position','bp','depth','info','quality')
    return(data)
  })
  #trim the data
  for(i in seq(1,length(file_infos))){
    file_infos[[i]] <- file_infos[[i]][front_trim:(nrow(file_infos[[i]]) - back_trim),]
  }
  return(file_infos)
}

paths <- c('final_dirs/8664_2/pileup.pileup','final_dirs/8664_3/pileup.pileup','final_dirs/8664_4/pileup.pileup')

thing2 <- build_image(paths,'','',600,600)

