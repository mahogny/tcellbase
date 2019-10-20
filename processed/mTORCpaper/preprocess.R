library(stringr)
library(Biobase)

dat <- read.csv("41590_2019_495_MOESM4_ESM.csv", sep="\t", stringsAsFactors = FALSE, check.names = FALSE)

dat <- dat[,c(
  "Gene.names",
  "CD4 naïve copies average",
  "CD4 TCR copies average",         
  "TH1 copies average",
  "CD8 TCR copies average",
  "CD8 naïve copies average")]

colnames(dat) <- c(
  "GeneSym",
  "CD4 naive",
  "CD4 tcr",
  "CD4 Th1",
  "CD8 tcr",
  "CD8 naive")

colnames(dat) 

tosplit <- grep(";",dat$GeneSym)

dat_done <- dat[-tosplit,]

dat_fix <- dat[tosplit,]

ensconvert <- read.csv("../../processed/ensembl_mouse.csv", stringsAsFactors = FALSE)
head(ensconvert)

for(currow in 1:nrow(dat_fix)){
  genesyms <- str_split(dat_fix$GeneSym[currow],";")[[1]]
  
  if(any(genesyms %in% ensconvert$Associated.Gene.Name)){
    genesym <- intersect(genesyms,ensconvert$Associated.Gene.Name)[1]
    dat_fix$GeneSym[currow] <- genesym
  } else {
    print(sprintf("missing %s",dat_fix$GeneSym[currow]))
    dat_fix$GeneSym[currow] <- NA
  }
}
dat_fix <- dat_fix[!is.na(dat_fix$GeneSym),]

dat_all <- rbind(dat_done,dat_fix)
dim(dat_all)

all(isUnique(ensconvert$Associated.Gene.Name)) ## no
all(isUnique(ensconvert$Ensembl.Gene.ID))   ### yes

colnames(dat_all)[1] <- "Associated.Gene.Name"

out <- merge(ensconvert,dat_all)

head(out)
out$Associated.Gene.Name[duplicated(out$Ensembl.Gene.ID)]
out$Ensembl.Gene.ID[duplicated(out$Ensembl.Gene.ID)]

out[out$Ensembl.Gene.ID=="ENSMUSG00000026319",]
out[out$Ensembl.Gene.ID=="ENSMUSG00000018697",]

### Filter out rows with NAs - can we do better?
keep <- apply(!is.na(out[,-(1:2)]),1,all)  ### Fairly brutal! Can we save more?
out <- out[keep,]


out$Associated.Gene.Name[duplicated(out$Ensembl.Gene.ID)]
out$Ensembl.Gene.ID[duplicated(out$Ensembl.Gene.ID)]


out[out$Ensembl.Gene.ID=="ENSMUSG00000036880",]
out[out$Ensembl.Gene.ID=="ENSMUSG00000025791",]

library("sqldf")

out <- sqldf("select `Associated.Gene.Name`,`Ensembl.Gene.ID`, avg(`CD4 naive`) as `CD4 naive`, avg(`CD4 tcr`) as `CD4 tcr`, avg(`CD4 Th1`) as `CD4 Th1`, avg(`CD8 tcr`) as `CD8 tcr`, avg(`CD8 naive`) as `CD8 naive` from out group by `Ensembl.Gene.ID`")

out$Associated.Gene.Name[duplicated(out$Ensembl.Gene.ID)]
out$Ensembl.Gene.ID[duplicated(out$Ensembl.Gene.ID)]

dim(out)

colnames(out)[1:2] <- c("Associated Gene Name","Ensembl Gene ID")
write.table(out,"processed_raw.csv",row.names = FALSE,sep="\t",quote = FALSE)

max(out[,4])

hist(log10(out[,5]))


#CD4 NAÏVE REP *
#[7] "CD4 TCR REP 1"                             "CD4 TCR REP 2"                             "CD4 TCR REP 3"                            
#[13] "TH1 REP 1"                                 "TH1 REP 2"                                 "TH1 REP 3"                                
#[19] "CD8 NAÏVE REP 1"                           "CD8 NAÏVE REP 2"                           "CD8 NAÏVE REP 3"                          [22] "CD8 NAÏVE REP 4"                           "CD8 NAÏVE REP 5"                           "CD8 NAÏVE REP 6"                          
#[25] "CD8 TCR REP 1"                             "CD8 TCR REP 2"                             "CD8 TCR REP 3"                            
#[31] "CTL REP 1"                                 "CTL REP 2"                                 "CTL REP 3"                                #what is this?

