library(stringr)
library(readxl)
library(reshape2)
library(plyr)
library(readr)
library(dplyr)
library(tidyr)
library(XML)
library(data.table)
library(magrittr)

IndexCal <- function(Index){
  Index <- readLines(Index)
  Index <- str_replace_all(Index, "<a href=\"", "")
  Index <- str_replace_all(Index, "\">Peak List</a>", "")
  Index <- htmlParse(Index)
  Index <- readHTMLTable(Index)
  Index <- as.data.frame(Index[["NULL"]])
  Index <- data.frame(lapply(Index, as.character), stringsAsFactors = FALSE)
  colnames(Index) <- str_replace_all(Index[1,], "Â", "")
  Index <- subset.data.frame(Index, Index$Confidence == "High")
  #Index <- Index[,seq(2, 34, by = 2)]
  #when Target Decoy Validator for FDR
  Index <- Index[,seq(2, 40, by = 2)]
  #when Percolator for FDR
  Index$File <- Index$`Peak List`
  return(Index)
}


SubFrag <- function(Address){
  File1 <- read.delim(Address, stringsAsFactors=FALSE)
  File1$File <- Address
  FragList <- rbind(FragList, File1)
}
