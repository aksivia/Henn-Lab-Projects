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

# tinkering with a simplified idea of isolating sample IDs
Omixon$Sample
Omixon$Sample <- str_replace(Omixon$Sample, "(\d*-|\d*_).*", "\\1")
w <- str_replace(Omixon[923,1], "([0-9]*)[:punct:].*", "\\1")

# GenDX$Sample <- str_replace(GenDX$Sample, "(\d*-|\d*_).*", "\\1")

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

GenDX <- GenDX[1:10,]
Omixon <- Omixon[1:10,]
# subset data to make testing script faster

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
# works

GenDX[is.na(GenDX)]<-""

# only comparing the first two fields of the data points to the P_group column to make the search faster
# make the loci match -- keep letter and make sure it matches P_group$Locus and only search through those p groups
# locus <- str_replace(GenDX[k,i], ".*\\*([0-9]{2,3}:[^:]{2,3}).*", "\\1:P") # also works


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

# make sure first field is matching first field in the exceptions and not the second or third field (especially if call only has 2-3 fields)

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
      print(locus);
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
      allele <- str_replace(Omixon[k,i], ".*\\*(.*)", "/\\1");
      r<-1;
      while (r<=nrow(P_group)) {
        if(P_group$Locus[r]==locus & str_detect(P_group$Alleles[r], allele)==TRUE){
          Omixon[k,i] <- P_group$P_group[r]}; # will be labeled with P at the end to denote that P group has been identified
        r<-r+1};
    }
    else{Omixon[k,i]<-Omixon[k,i]};
    i<-i+1
  };
  k<-k+1
}


# talk to Gillian about how to handle multiple alleles for same sample ID and lots of blank spaces
# maybe find all p groups for all data points, then compare between both allele 1 versions for a patient and combine?
# or compare allele 1 and allele 2 data for the same patient and combine?


# at the end, after all the P groups have been found, can go through and remove all the ":P" at the end of all the data points
# keep all :P for now, so at the end we can identify any data points that did not match to a P group at all


write.table(Omixon, file="Omixon_Pgroup")
write.table(GenDX, file="GenDX_Pgroup")


