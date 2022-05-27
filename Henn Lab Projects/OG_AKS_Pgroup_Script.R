---
  title: "Omixon_GenDx"
output: html_document
---
  
  setwd("/Users/aks/Henn Lab Projects") # change for your specific wd

library(tidyverse)
library(stats)
library(stringi)
library(stringr)
require(tidyr)
require(ggplot2)
library(plyr)
library(pipeR)
library(rlist)
require(dplyr)
require(data.table)
library(janitor)
library(rapport)

GenDX <- read.delim(file="GenDX_updated_7_26_Himba_resolved_reformat.txt")
Omixon <- read.delim(file="14June2021_Himba_Omixon_HLA_calls.txt")

#Getting rid of excess stuff in Sample_ID and making sure the IDs all in the same column
Omixon <- separate(Omixon, Sample, c("Sample_ID", "extra"), sep="-")
Omixon <- separate(Omixon, Sample_ID, c("Sample_ID", "extra"), sep="_")
Omixon <- Omixon[,-2]
Omixon <- separate(Omixon, Sample_ID, c("HMB", "Sample_ID"), sep="HMB")
GenDX <- separate(GenDX, Sample, c("Sample_ID", "extra"), sep="-")
GenDX <- separate(GenDX, Sample_ID, c("Sample_ID", "extra"), sep="_")
GenDX <- GenDX[,-2]
GenDX <- separate(GenDX, Sample_ID, c("HMB", "Sample_ID"), sep="HMB")
# run code up to here and see what's different about row 923 and above for Omixon, and 174 above for GenDX, etc.
# for GenDX and Omixon remove HMB and D at the beginning of all sample IDs
# also reference NCTB data (all samples start with 'NC')
# because that is what will be used in script (except with the peoples' genetic data)
# so can get rid of below while loops using specific numbers

k<-1
while (k<=nrow(GenDX)) {
  if (is.na(GenDX[k,2])==TRUE){
    GenDX[k,2] <- GenDX[k,1]};
  k<-k+1
}


k<-1
while (k<=nrow(Omixon)) {
  if (is.na(Omixon[k,2])==TRUE){
    Omixon[k,2] <- Omixon[k,1]};
  k<-k+1
}


GenDX <- select(GenDX, Sample_ID, HLA_A, HLA_B, HLA_C, DPA1, DPB1, DQA1, DQB1, DRB1) # remove empty HMB columns
Omixon <- select(Omixon, -HMB)  #fix for omixon file

GenDX$Sample_ID <- gsub('[^0-9]', '', GenDX$Sample_ID) # removes all letters in Sample ID's
Omixon$Sample_ID <- gsub('[^0-9]', '', Omixon$Sample_ID) 

# another method:
# GenDX$Sample_ID <- gsub('[A-Z]*', '', GenDX$Sample_ID) 
# GenDX$Sample_ID <- gsub('[a-z]*', '', GenDX$Sample_ID)

# should I add 'NC' to the start of all sample ID's to match the NCTB data? or remove the 'NC' in the NCTB data sample ID's?


GenDX <- GenDX[1:20,]
Omixon <- Omixon[1:20,]
# subset data to only keep rows 1-50 to make testing script faster


GenDX <- separate(GenDX, HLA_A, c("A_1", "A_2"), sep=",")
GenDX <- separate(GenDX, HLA_B, c("B_1", "B_2"), sep=",")
GenDX <- separate(GenDX, HLA_C, c("C_1", "C_2"), sep=",")
GenDX <- separate(GenDX, DPA1, c("DPA1_1", "DPA1_2"), sep=",")
GenDX <- separate(GenDX, DPB1, c("DPB1_1", "DPB1_2"), sep=",")
GenDX <- separate(GenDX, DQA1, c("DQA1_1", "DQA1_2"), sep=",")
GenDX <- separate(GenDX, DQB1, c("DQB1_1", "DQB1_2"), sep=",")
GenDX <- separate(GenDX, DRB1, c("DRB1_1", "DRB1_2"), sep=",")


k<-1 # start at row 1
while (k<=nrow(GenDX)) {
  GenDX[k,] <- str_remove_all(GenDX[k,], ".*\\*");
  k<-k+1 # go through each row
}


# gsub(".*\\*", "", GenDX[k,])

Omixon <- Omixon[,c(1:13)]
Omixon <- Omixon[,-c(7,9,12)]

#take out the () at the end of some of the calls
Omixon$HLA.A <- word(Omixon$HLA.A, 1, sep = "\\(")
Omixon$HLA.B <-  word(Omixon$HLA.B, 1, sep = "\\(")
Omixon$HLA.C <-  word(Omixon$HLA.C, 1, sep = "\\(")
Omixon$HLA.DPA1 <- word(Omixon$HLA.DPA1, 1, sep = "\\(")
Omixon$HLA.DPB1 <- word(Omixon$HLA.DPB1, 1, sep = "\\(")
Omixon$HLA.DQA1 <- word(Omixon$HLA.DQA1, 1, sep = "\\(")
Omixon$HLA.DQB1 <- word(Omixon$HLA.DQB1, 1, sep = "\\(")
Omixon$HLA.DRB1 <- word(Omixon$HLA.DRB1, 1, sep = "\\(")


# remove the HMB column at the end?
library(tidyr)
library(stringr)

# to make sure we're always using the most updated version of the p groups file
# files are updated every 3 months
pgroup.url <- "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/wmda/hla_nom_p.txt"
destfile <- "/Users/aks/Henn Lab Projects/hla_nom_p.txt" # update with your own wd
download.file(pgroup.url, destfile)
P_group <- read.delim(file="hla_nom_p.txt")
P_group <- separate(P_group, X..file..hla_nom_p.txt, c("Locus_2", "Alleles", "P_group"), sep=";")
P_group$P_group <- gsub('[A-Z]*', '', P_group$P_group) # remove all letters in P_group column

#good to go
P_group$Locus <- "null"
k<-1
while(k<=nrow(P_group)){
  if(P_group[k,1] == "A*"){
    P_group[k,4] <- "A"
  }
  else if(P_group[k,1] == "B*"){
    P_group[k,4] <- "B"
  }
  else if(P_group[k,1] == "C*"){
    P_group[k,4] <- "C"
  }
  else if(P_group[k,1] == "DPA1*"){
    P_group[k,4] <- "DPA1"
  }
  else if(P_group[k,1] == "DPB1*"){
    P_group[k,4] <- "DPB1"
  }
  else if(P_group[k,1] == "DQA1*"){
    P_group[k,4] <- "DQA1"
  }
  else if(P_group[k,1] == "DQB1*"){
    P_group[k,4] <- "DQB1"
  }
  else if(P_group[k,1] == "DRB1*"){
    P_group[k,4] <- "DRB1"
  }
  k <- k+1
}

P_group$P_group[P_group$P_group==""] <- NA
k<-1
while (k<=nrow(P_group)) {
  if(is.na(P_group[k,3])){P_group[k,3]<-P_group[k,2]};
  k <- k+1
}

P_group <- P_group[-c(1:5),]

#sometimes for class2 the typing calls only have 3 field but the conversion sheet has 4, so need to drop the 4th field on the conversion sheet
my_fun_2<- function(i) {
  word(i, start=1, end=2, sep="\\:")
}
my_fun_3<- function(i) {
  word(i, start=1, end=3, sep="\\:")
}

#GenDX P_group conversion 

GenDX[is.na(GenDX)]<-""


k<-1
while(k<=nrow(GenDX))
{
  i <- 2;
  while(i<=17){
    if(GenDX[k,i]!=""){
      locus <- str_split(GenDX[k,i], "\\*") [[1]] [[1]];
      P<-P_group[P_group$Locus==locus,];
      if(locus=="A"|locus=="B"|locus=="C"){
        j<-1;
        while(j <= nrow(P)){
          print(j);
          if((paste(unlist(str_split(GenDX[k,i], ""))[-1:-2], collapse ="")) %in% sapply(list(str_split(P[j, 2], "/")) [[1]] [[1]], my_fun_3) ||
             (paste(unlist(str_split(GenDX[k,i], ""))[-1:-2], collapse ="")) %in% list(str_split(P[j, 2], "/")) [[1]] [[1]] ||
             (my_fun_3(paste(unlist(str_split(trimws(GenDX[k,i]), ""))[-1:-2], collapse =""))) %in% list(str_split(P[j, 2], "/")) [[1]] [[1]]){
            GenDX[k,i] <- P[j,3];print(P[j,3]);break}
          else{j<-j+1}};
        i <- i + 1
      };
      if(locus =="DPA1" | locus=="DPB1" | locus=="DQA1" | locus=="DQB1" | locus=="DRB1"){
        j<-1;
        while(j <= nrow(P)){
          print(j);
          if((paste(unlist(str_split(GenDX[k,i], ""))[-1:-5], collapse ="")) %in% sapply(list(str_split(P[j, 2], "/")) [[1]] [[1]], my_fun_3) ||
             (paste(unlist(str_split(GenDX[k,i], ""))[-1:-5], collapse ="")) %in% list(str_split(P[j, 2], "/")) [[1]] [[1]] ||
             (my_fun_3(paste(unlist(str_split(trimws(GenDX[k,i]), ""))[-1:-5], collapse =""))) %in% list(str_split(P[j, 2], "/")) [[1]] [[1]] || 
             (my_fun_3(paste(unlist(str_split(trimws(GenDX[k,i]), ""))[-1:-5], collapse =""))) %in% sapply(list(str_split(P[j, 2], "/")) [[1]] [[1]], my_fun_3)){
            GenDX[k,i]<-P[j,3];print(P[j,3]);break}
          else{j<-j+1}};
        i <- i + 1
      }}
    else{i <- i + 1}};
  k <- k + 1
}

write.table(GenDX, file="GenDX_Pgroup")


#Omixon P group conversion
k<-1
while(k<=nrow(Omixon))
{
  i <- 3;
  while(i<=10){
    if(Omixon[k,i]!=""){
      locus <- str_split(str_split(Omixon[k,i], "\\*") [[1]] [1], "-") [[1]] [2]; # error
      print(locus);
      P<-P_group[P_group$Locus==locus,];
      if(locus=="A"|locus=="B"|locus=="C"){
        j<-1;
        print(j);
        while(j <= nrow(P)){
          if((paste(unlist(str_split(trimws(Omixon[k,i]), ""))[-1:-6], collapse ="")) %in% sapply(list(str_split(P[j, 2], "/")) [[1]] [[1]], my_fun_3) ||
             (paste(unlist(str_split(trimws(Omixon[k,i]), ""))[-1:-6], collapse ="")) %in% list(str_split(P[j, 2], "/")) [[1]] [[1]] ||
             (my_fun_3(paste(unlist(str_split(trimws(Omixon[k,i]), ""))[-1:-6], collapse =""))) %in% list(str_split(P[j, 2], "/")) [[1]] [[1]]){
            Omixon[k,i] <- P[j,3];print(P[j,3]);break}
          else{j<-j+1}};
        i <- i + 1
      };
      if(locus =="DPA1" | locus=="DPB1" | locus=="DQA1" | locus=="DQB1" | locus=="DRB1"){
        j<-1;
        while(j <= nrow(P)){
          if((paste(unlist(str_split(trimws(Omixon[k,i]), ""))[-1:-9], collapse ="")) %in% sapply(list(str_split(P[j, 2], "/")) [[1]] [[1]], my_fun_3) ||
             (paste(unlist(str_split(trimws(Omixon[k,i]), ""))[-1:-9], collapse ="")) %in% list(str_split(P[j, 2], "/")) [[1]] [[1]] ||
             (my_fun_3(paste(unlist(str_split(trimws(Omixon[k,i]), ""))[-1:-9], collapse =""))) %in% list(str_split(P[j, 2], "/")) [[1]] [[1]]){
            Omixon[k,i] <- P[j,3];print(P[j,3]);break}
          else{j<-j+1}};
        i <- i + 1
      }}
    else{i<- i+1}};
  k <- k + 1
}

write.table(Omixon, file="Omixon_Pgroup")
write.table(GenDX, file="GenDX_Pgroup")


