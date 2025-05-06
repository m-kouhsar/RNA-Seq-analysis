message("#########################################################")
message("kallisto pipeline version 1.0.0")
message("Wrote by m.kouhsar@exeter.ac.uk")
message("R script for merging results")
message("#########################################################")

Arguments <- commandArgs(T)

ResultsDir <- trimws(Arguments[1])
OutPrefix <- trimws(Arguments[2])
if(is.na(OutPrefix)){
  OutPrefix <- ""
}

suppressMessages(library(stringr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))

message("Reading results directory...")
abu_files <- list.files(path = ResultsDir , pattern = "abundance.tsv" , full.names = T , recursive = T)
abu_files <- str_split(abu_files , pattern = "/",simplify = T)
message(nrow(abu_files) , " samples are detected.")

message("Merging...")
message("     Sample 1 ...")

length <- read.delim(file = paste(abu_files[1,], collapse = "/"),header = T,sep = "\t")[,c(1,2)]
colnames(length)[2] <- abu_files[1,ncol(abu_files)-1]

eff_length <- read.delim(file = paste(abu_files[1,], collapse = "/"),header = T,sep = "\t")[,c(1,3)]
colnames(eff_length)[2] <- abu_files[1,ncol(abu_files)-1]

est_counts <- read.delim(file = paste(abu_files[1,], collapse = "/"),header = T,sep = "\t")[,c(1,4)]
colnames(est_counts)[2] <- abu_files[1,ncol(abu_files)-1]

tpm <- read.delim(file = paste(abu_files[1,], collapse = "/"),header = T,sep = "\t")[,c(1,5)]
colnames(tpm)[2] <- abu_files[1,ncol(abu_files)-1]

for (i in 2:nrow(abu_files)) {
  
  message("     Sample ",i," ...")
  
  length.temp <- read.delim(file = paste(abu_files[i,], collapse = "/"),header = T,sep = "\t")[,c(1,2)]
  colnames(length.temp)[2] <- abu_files[i,ncol(abu_files)-1]
  
  eff_length.temp <- read.delim(file = paste(abu_files[i,], collapse = "/"),header = T,sep = "\t")[,c(1,3)]
  colnames(eff_length.temp)[2] <- abu_files[i,ncol(abu_files)-1]
  
  est_counts.temp <- read.delim(file = paste(abu_files[i,], collapse = "/"),header = T,sep = "\t")[,c(1,4)]
  colnames(est_counts.temp)[2] <- abu_files[i,ncol(abu_files)-1]
  
  tpm.temp <- read.delim(file = paste(abu_files[i,], collapse = "/"),header = T,sep = "\t")[,c(1,5)]
  colnames(tpm.temp)[2] <- abu_files[i,ncol(abu_files)-1]
  
  length <- full_join(length , length.temp , by= "target_id") 
  eff_length <- full_join(eff_length , eff_length.temp , by= "target_id") 
  est_counts <- full_join(est_counts , est_counts.temp , by= "target_id") 
  tpm <- full_join(tpm , tpm.temp , by= "target_id") 
  
}

length <- mutate_all(length,~replace_na(., 0))
eff_length <- mutate_all(eff_length,~replace_na(., 0))
est_counts <- mutate_all(est_counts,~replace_na(., 0))
tpm <- mutate_all(tpm,~replace_na(., 0))

message("Saving merged files...")
write.table(length , file = paste(OutPrefix , "kallisto.length.txt" ,sep = ifelse(OutPrefix == "" , "",".")),
            row.names = F,col.names = T,sep = '\t',quote = F)
write.table(eff_length , file = paste(OutPrefix , "kallisto.eff_length.txt" ,sep = ifelse(OutPrefix == "" , "",".")),
            row.names = F,col.names = T,sep = '\t',quote = F)
write.table(est_counts , file = paste(OutPrefix , "kallisto.est_counts.txt" ,sep = ifelse(OutPrefix == "" , "",".")),
            row.names = F,col.names = T,sep = '\t',quote = F)
write.table(tpm , file = paste(OutPrefix , "kallisto.TPM.txt" ,sep = ifelse(OutPrefix == "" , "",".")),
            row.names = F,col.names = T,sep = '\t',quote = F)

message("All done!")



