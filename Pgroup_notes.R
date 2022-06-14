
# another method to remove all letters from Sample_ID column
# GenDX$Sample_ID <- gsub('[A-Z]*', '', GenDX$Sample_ID) 
# GenDX$Sample_ID <- gsub('[a-z]*', '', GenDX$Sample_ID)



# tinkering with a simplified idea of isolating sample IDs
Omixon$Sample
Omixon$Sample <- str_replace(Omixon$Sample, "(\d*-|\d*_).*", "\\1")
w <- str_replace(Omixon[923,1], "([0-9]*)[:punct:].*", "\\1")

# GenDX$Sample <- str_replace(GenDX$Sample, "(\d*-|\d*_).*", "\\1")

# locus <- str_replace(GenDX[k,i], ".*\\*([0-9]{2,3}:[^:]{2,3}).*", "\\1:P") # also works
# another way to parse locus info using regex



# if the Sample ID is in the HMB column and not the Sample_ID column 
# then this moves the Sample ID into the Sample_ID
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



#GenDX P_group conversion 

# works
k<-1
while(k<=nrow(GenDX))
{
  i<-2;
  while (i<=ncol(GenDX)) {
    if(GenDX[k,i]!=""){
      GenDX[k,i] <- str_replace(GenDX[k,i], "(.{2,3}:..):.*", "\\1")}
    else{GenDX[k,i]<-GenDX[k,i]};
    print(GenDX[k,i]);
    i<-i+1
  };
  k<-k+1
}

GenDX[is.na(GenDX)]<-""


k<-1
while(k<=nrow(GenDX))
{
  i<-2;
  while (i<=ncol(GenDX)) {
    if(GenDX[k,i]!=""){
      locus <- str_replace(GenDX[k,i], "(..:..):.*", "\\1");
      r <- 1;
      while (r<=nrow(P_group)) {
        if(P_group$P_group[r]==locus){
          GenDX[k,i] <- P_group$P_group[r]};
        r<-r+1
      };
      i<-i+1
    }
  };
  k<-k+1
}




#duplicate
k<-1
while(k<=nrow(GenDX))
{
  i<-2;
  while (i<=ncol(GenDX)) {
    if(GenDX[k,i]!=""){
      GenDX[k,i] <- str_replace(GenDX[k,i], "(..:..):.*", "\\1");
      #compare to P_group$P_group
      #then replace field info with P_group
      i<-i+1
    }
  };
  k<-k+1
}

#one idea, doesn't work
k<-1
while(k<=nrow(GenDX))
{
  i<-2;
  while (i<=ncol(GenDX)) {
    if(GenDX[k,i]!=""){
      locus <- str_replace(GenDX[k,i], "(..:..):.*", "\\1");
      for pgroup in P_group$P_group
      do 
      if(pgroup==locus) {
        GenDX[k,i] <- locus}
      done;
      i<-i+1
    }
  };
  k<-k+1
}

# regular expressions to parse gene information and only compare first two fields to first two fields in p group file
# can remove letter and asterisk 
# "(.*)(..:..):..:.."" remove \\1 and compare \\2 to P_group$P_group

A*02:01:01:01