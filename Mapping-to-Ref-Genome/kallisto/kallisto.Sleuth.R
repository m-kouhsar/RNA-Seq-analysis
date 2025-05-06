message("#########################################################")
message("kallisto pipeline version 1.0.0")
message("Wrote by m.kouhsar@exeter.ac.uk")
message("R script for running DEG analysis using the Sleuth package")
message("#########################################################")

Arguments <- commandArgs(T)

kallisto_res_dir <- trimws(Arguments[1])
pheno_file <- trimws(Arguments[2]) 

# pheno_file is a csv file which must contains the following columns:
#        sample: Samples ID
#        path: path of the kallisto resuls (inside kallisto_res_dir) folder for each sample
#        All the variables that are included in lm_model must be represented by a column with the same name

lm_model <- trimws(Arguments[3])
var_factor <- trimws(Arguments[3])    #Factor variable in lm_model
var_numeric <- trimws(Arguments[3])   #Numerical varaibles in lm_model
OutPrefix <- trimws(Arguments[4])
if(is.na(OutPrefix)){
  OutPrefix <- ""
}

suppressMessages(library(stringr))
suppressMessages(library(sleuth))


