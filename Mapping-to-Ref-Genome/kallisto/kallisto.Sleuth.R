message("#########################################################")
message("kallisto pipeline version 1.0.0")
message("Wrote by m.kouhsar@exeter.ac.uk")
message("R script for running DEG analysis using the Sleuth package")
message("#########################################################")


suppressMessages(library(stringr))
suppressMessages(library(sleuth))
suppressMessages(library("biomaRt"))

#################### Preparing input arguments #################

Arguments <- commandArgs(T)

kallisto_res_dir <- trimws(Arguments[1])
pheno_file <- trimws(Arguments[2]) 
# Target mapping file contains information (eg. gene ID and type) about the transcripts
target_map_file <- trimws(Arguments[3]) 
# pheno_file is a csv file which must contains the following columns:
#        sample: Samples ID
#        path: path of the kallisto resuls (inside kallisto_res_dir) folder for each sample
#        All the variables that are included in lm_model must be represented by a column with the same name

lm_model <- as.formula(Arguments[4])
# Factor variable include condition variable in lm_model
var_factor <- trimws(str_split_1(Arguments[5],pattern = ","))    
# Numerical varaibles in lm_model
var_numeric <- trimws(str_split_1(Arguments[6],pattern = ","))   
OutPrefix <- trimws(Arguments[6])

if(is.na(OutPrefix)){
  OutPrefix <- ""
}

message("Input arguments:")
message("        kallisto results directory: ",kallisto_res_dir)
message("        Phenotype file: ",pheno_file)
message("        Target mapping file: ",target_map_file)
message("        Linear regression model: ",lm_model)
message("        Factor variables in the model: ",var_factor)
message("        Numeric variables in the model: ",var_numeric)
message("        Output files prefix: ",OutPrefix)


if(!all(tolower(all.vars(lm_model)) %in% tolower(c(var_factor,var_numeric)))){
  stop("The following variables in the regression model are not specified as factor or numeric variables:\n",
       setdiff(all.vars(lm_model) , c(var_factor,var_numeric)))
}

pheno <- read.csv(pheno_file , stringsAsFactors = F)

if(!all(c("sample","path",tolower(all.vars(lm_model))) %in% tolower(colnames(pheno)))){
  
  stop("The following variables can't be find in phenotype file: \n",
       setdiff(c("sample","path",all.vars(lm_model)), colnames(pheno)))
}

for (v in var_factor) {
  pheno[,v] <- as.factor(pheno[,v])
}
for (v in var_numeric) {
  pheno[,v] <- as.numeric(pheno[,v])
}
pheno$path <- paste(kallisto_res_dir , pheno$path , sep = "/")

################### Reading kallisto results #################

SO_file <- paste0(OutPrefix , ".SO.rdat")

if(file.exists(SO_file)){
  message("Loading Sleuth Object from ",SO_file, " ...")
  load(SO_file)
}else{
  target_map <- read.delim(target_map_file, stringsAsFactors = F , header = T , sep = "\t")
  so <- sleuth_prep(pheno,target_mapping = target_map,extra_bootstrap_summary = TRUE)
  
  message("Saving Sleuth Object in ",SO_file," ...")
  save(so,file = SO_file)
}

################## General QC ############################




################### Finding DEGs using linear regression ###############















