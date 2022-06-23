---
title: "Omixon_GenDx"
output: html_document
---
 
setwd("/Users/aks/Henn Lab Projects") # change to your wd

library(tidyverse)
library(stats)
library(stringi)
library(stringr)
library(tidyr)
require(ggplot2)
library(plyr)
library(pipeR)
library(rlist)
require(dplyr)
require(data.table)
library(janitor)
library(rapport)

# install.packages("devtools")
# devtools::install_github("gadenbuie/regexplain")


GenDX <- read.delim(file="GenDX_updated_7_26_Himba_resolved_reformat.txt")
Omixon <- read.delim(file="14June2021_Himba_Omixon_HLA_calls.txt")

#Getting rid of excess stuff in the Sample ID's
GenDX <- separate(GenDX, Sample, c("Sample_ID", "extra"), sep="-")
GenDX <- separate(GenDX, Sample_ID, c("Sample_ID", "extra"), sep="_")
GenDX <- GenDX[,-2]
GenDX$Sample_ID <- gsub('[^0-9]', '', GenDX$Sample_ID) # removes all letters in Sample ID's
GenDX <- select(GenDX, Sample_ID, HLA_A, HLA_B, HLA_C, DPA1, DPB1, DQA1, DQB1, DRB1) # only keep relevant columns

Omixon <- separate(Omixon, Sample, c("Sample_ID", "extra"), sep="-")
Omixon <- separate(Omixon, Sample_ID, c("Sample_ID", "extra"), sep="_")
Omixon <- Omixon[,-2]
Omixon$Sample_ID <- gsub('[^0-9]', '', Omixon$Sample_ID) 
Omixon <- select(Omixon, Sample_ID, Allele, HLA.A, HLA.B, HLA.C, HLA.DPA1, HLA.DPB1, HLA.DQA1, HLA.DQB1, HLA.DRB1) 
colnames(Omixon) <- gsub("HLA.", "", colnames(Omixon))
Omixon$Allele <- str_remove_all(Omixon$Allele, "ALLELE_")

GenDX <- separate(GenDX, HLA_A, c("A_1", "A_2"), sep=",")
GenDX <- separate(GenDX, HLA_B, c("B_1", "B_2"), sep=",")
GenDX <- separate(GenDX, HLA_C, c("C_1", "C_2"), sep=",")
GenDX <- separate(GenDX, DPA1, c("DPA1_1", "DPA1_2"), sep=",")
GenDX <- separate(GenDX, DPB1, c("DPB1_1", "DPB1_2"), sep=",")
GenDX <- separate(GenDX, DQA1, c("DQA1_1", "DQA1_2"), sep=",")
GenDX <- separate(GenDX, DQB1, c("DQB1_1", "DQB1_2"), sep=",")
GenDX <- separate(GenDX, DRB1, c("DRB1_1", "DRB1_2"), sep=",")

# take out the () and spaces at the end of all data points in Omixon
k<-1
while (k<=nrow(Omixon)) {
  Omixon[k,] <- str_remove_all(Omixon[k,], " \\(.*");
  k<-k+1
}


# downloads file via url to make sure we're always using the most updated version of the p groups file
# NOTE: files are updated every 3 months
pgroup.url <- "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/wmda/hla_nom_p.txt"
destfile <- "/Users/aks/Henn Lab Projects/hla_nom_p.txt" # change to your wd
download.file(pgroup.url, destfile)
P_group <- read.delim(file="hla_nom_p.txt")
P_group <- separate(P_group, X..file..hla_nom_p.txt, c("Locus", "Alleles", "P_group"), sep=";")


P_group$P_group[P_group$P_group==""] <- NA
k<-1
while (k<=nrow(P_group)) {
  if(is.na(P_group[k,3])){P_group[k,3]<-P_group[k,2]};
  k <- k+1
}


P_group$P_group <- gsub('[A-Z]*', '', P_group$P_group) # remove all letters in P_group column
P_group$P_group <- str_replace(P_group$P_group, "(.*)", "\\1:P") # put the letter P at the end of each P group
P_group$Alleles <- str_replace(P_group$Alleles, "(.*)", "/\\1")


# cleans up Locus column
k<-1
while(k<=nrow(P_group)){
  if(P_group[k,1] == "A*"){
    P_group[k,1] <- "A"
  }
  else if(P_group[k,1] == "B*"){
    P_group[k,1] <- "B"
  }
  else if(P_group[k,1] == "C*"){
    P_group[k,1] <- "C"
  }
  else if(P_group[k,1] == "DPA1*"){
    P_group[k,1] <- "DPA1"
  }
  else if(P_group[k,1] == "DPB1*"){
    P_group[k,1] <- "DPB1"
  }
  else if(P_group[k,1] == "DQA1*"){
    P_group[k,1] <- "DQA1"
  }
  else if(P_group[k,1] == "DQB1*"){
    P_group[k,1] <- "DQB1"
  }
  else if(P_group[k,1] == "DRB1*"){
    P_group[k,1] <- "DRB1"
  }
  else{P_group[k,1] <- "null"}
  k <- k+1
}


P_group <- filter(P_group, Locus!="null") # removes irrelevant rows


# GenDX P group conversion

GenDX[is.na(GenDX)]<-""
# empty cells in GenDX represent NA's (meaning there was no data/sample for that patient at that loci)

# only comparing the first two fields of the data points to the P_group column to make the search faster
# make the loci match so we only search through those p groups
k<-1
while(k<=nrow(GenDX)) {
  i<-2;
  while (i<=ncol(GenDX)) {
    if(GenDX[k,i]!=""){
      locus <- str_replace(GenDX[k,i], "(.*)\\*.*:.*", "\\1");
      allele <- str_replace(GenDX[k,i], ".*\\*([0-9]{2,3}:[0-9]{2,3}).*", "\\1:P"); # add the P to enable search
      r <- 1;
      while (r<=nrow(P_group)) {
        if(P_group$Locus[r]==locus & P_group$P_group[r]==allele){
          GenDX[k,i] <- P_group$P_group[r]}; # will be labeled with P at the end to denote that P group has been identified
        r<-r+1};
      }
    else{GenDX[k,i]<-GenDX[k,i]};
    i<-i+1
  };
  k<-k+1
}

# to find any exceptions (i.e. data points where the first 2 fields do not match the p group name)
# only looking at the data points without 'P' at the end and run them through all the P_group$Alleles to find their p group
k<-1
while (k<=nrow(GenDX)) {
  i<-2;
  while (i<=ncol(GenDX)) {
    if(str_detect(GenDX[k,i], ":P$")==FALSE){
      locus <- str_replace(GenDX[k,i], "(.*)\\*.*:.*", "\\1");
      allele <- str_replace(GenDX[k,i], ".*\\*(.*)", "/\\1");
      r<-1;
      while (r<=nrow(P_group)) {
        if(P_group$Locus[r]==locus & str_detect(P_group$Alleles[r], allele)==TRUE){
          GenDX[k,i] <- P_group$P_group[r]}; # will be labeled with P at the end to denote that P group has been identified
        r<-r+1};
      }
    else{GenDX[k,i]<-GenDX[k,i]};
    i<-i+1
  };
  k<-k+1
}


# Omixon P group conversion

# only comparing the first 2 fields of the data to the p groups to find a match (quicker method)
k<-1
while(k<=nrow(Omixon)) {
  i<-3;
  while (i<=ncol(Omixon)) {
    if(Omixon[k,i]!=""){
      locus <- str_replace(Omixon[k,i], "HLA-(.*)\\*.*:.*", "\\1");
      allele <- str_replace(Omixon[k,i], ".*\\*([0-9]{2,3}:[0-9]{2,3}).*", "\\1:P"); # add the P to enable search
      r <- 1;
      while (r<=nrow(P_group)) {
        if(P_group$Locus[r]==locus & P_group$P_group[r]==allele){
          Omixon[k,i] <- P_group$P_group[r]}; # will be labeled with P at the end to denote that P group has been identified
        r<-r+1};
    }
    else{Omixon[k,i]<-Omixon[k,i]};
    i<-i+1
  };
  k<-k+1
}


# to find exceptions
k<-1
while (k<=nrow(Omixon)) {
  i<-3;
  while (i<=ncol(Omixon)) {
    if(str_detect(Omixon[k,i], ":P$")==FALSE && Omixon[k,i]!=""){
      locus <- str_replace(Omixon[k,i], "HLA-(.*)\\*.*:.*", "\\1");
      allele <- str_replace(Omixon[k,i], ".*\\*([0-9]{2,3}:[0-9]{2,3}).*", "/\\1");
      r<-1;
      while (r<=nrow(P_group)) {
        if(P_group$Locus[r]==locus & str_detect(P_group$Alleles[r], allele)==TRUE){
          Omixon[k,i] <- P_group$P_group[r]}; # will be labeled with P at the end to denote that P group has been identified
        r<-r+1};
    }
    else if(Omixon[k,i]==""){Omixon[k,i] <- NA}
    else{Omixon[k,i]<-Omixon[k,i]};
    i<-i+1
  };
  k<-k+1
}


# reformat Omixon data frame so it matches the GenDX data frame
Omixon <- Omixon[rowSums(is.na(Omixon)) != 8,] # remove empty rows from Omixon

Omixon <- pivot_wider(Omixon, names_from = "Allele", values_from = c("A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"))

Omixon2 <- data.frame(matrix(ncol = 17, nrow=0))
colnames(Omixon2) <- colnames(Omixon)

k<-1
while (k<=nrow(Omixon)) {
  i<-1
  while (i<=ncol(Omixon)) {
    x <- unlist(Omixon[k,i]);
    x <- x[!is.na(x)];
    if (is_empty(x) == FALSE){
      if (length(unique(x)) == 1){Omixon2[k,i] <- unique(x)}
      else if(length(unique(x)) > 1){Omixon2[k,i] <- toString(unique(x))}}
    else if(is_empty(x) == TRUE){Omixon2[k,i] <- ""};
    i<-i+1
  };
  k<-k+1
}




write.table(Omixon2, file="Omixon_Pgroup")
write.table(GenDX, file="GenDX_Pgroup")


