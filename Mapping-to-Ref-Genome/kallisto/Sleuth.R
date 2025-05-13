message("#########################################################")
message("kallisto pipeline version 1.0.0")
message("Wrote by m.kouhsar@exeter.ac.uk")
message("R script for running DEG analysis using the Sleuth package")
message("#########################################################")

suppressMessages(library(stringr))
suppressMessages(library(sleuth))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(tibble))

#################### Preparing input arguments #################

Arguments <- commandArgs(T)

kallisto_res_dir <- trimws(Arguments[1])
# pheno_file is a csv file which must contains the following columns:
#        sample: Samples ID
#        path: path of the kallisto resuls (inside kallisto_res_dir) folder for each sample
#        All the variables that are included in lm_model must be represented by a column with the same name
pheno_file <- trimws(Arguments[2]) 
# Target mapping file contains information (eg. gene ID and type) about the transcripts
target_map_file <- trimws(Arguments[3]) 
lm_model <- trimws(Arguments[4])
# Factor variable include condition variable in lm_model
var_factor <- trimws(Arguments[5])    
# Numerical varaibles in lm_model
var_numeric <- trimws(Arguments[6])  
PCs <- as.numeric(trimws(Arguments[7]))
OutPrefix <- trimws(Arguments[8])
ScriptDir <- trimws(Arguments[9])

if(is.na(OutPrefix)){
  OutPrefix <- ""
}
source(paste0(ScriptDir , "/mahalanobis.outlier.R"))

lm_model <- as.formula(lm_model)
var_factor <- str_split_1(var_factor , pattern = ",")
var_numeric <- str_split_1(var_numeric , pattern = ",")

message("Input arguments:")
message("        kallisto results directory: ",kallisto_res_dir)
message("        Phenotype file: ",pheno_file)
message("        Target mapping file: ",target_map_file)
message("        Linear regression model: ",lm_model)
message("        Factor variables in the model: ",paste(var_factor, collapse = ", "))
message("        Numeric variables in the model: ",paste(var_numeric,collapse=", "))
message("        Number of PCs to add to the model: ",PCs)
message("        Output files prefix: ",OutPrefix)

message("#########################################################")

if(!all(tolower(all.vars(lm_model)) %in% tolower(c(var_factor,var_numeric)))){
  stop("The following variables in the regression model are not specified as factor or numeric variables:\n",
       paste(setdiff(all.vars(lm_model) , c(var_factor,var_numeric)),collapse = ", "))
}

pheno <- read.csv(pheno_file , stringsAsFactors = F)

if(!all(c("sample","path",tolower(all.vars(lm_model))) %in% tolower(colnames(pheno)))){
  
  stop("The following variables can't be find in phenotype file: \n",
       paste(setdiff(c("sample","path",all.vars(lm_model)), colnames(pheno)),collapse = ", "))
}

for (v in var_factor) {
  pheno[,v] <- as.factor(pheno[,v])
}
for (v in var_numeric) {
  pheno[,v] <- as.numeric(pheno[,v])
}
pheno$path <- paste(kallisto_res_dir , pheno$path , sep = "/")

################### Reading kallisto results #################

target_map <- read.delim(target_map_file, stringsAsFactors = F , header = T , sep = "\t")

SO_file <- paste0(OutPrefix , ".SO.rdat")

if(file.exists(SO_file)){
  
  message("Loading Sleuth Object from ",SO_file, " ...")
  load(SO_file)
  
}else{
  
  
  so <- sleuth_prep(pheno,target_mapping = target_map,extra_bootstrap_summary = TRUE)
  
  message("Saving Sleuth Object in ",SO_file," ...")
  save(so,file = SO_file)
}

################### Removing Outliers ##################################

message("Finding Outlier samples using Mahalanobis distance and chi-squared distribution...")

tpm.norm <- so$obs_norm[, c(1,2,4)]
tpm.norm <- pivot_wider(tpm.norm, id_cols = "target_id", names_from = "sample", values_from = "tpm")
tpm.norm <- column_to_rownames(tpm.norm, var = "target_id")
tpm.norm <- as.matrix(tpm.norm)

Outliers <- mahalanobis.outlier(Data=tpm.norm , method="pca",pca.scale=F)

Outliers.id <- Outliers$Data.2D[Outliers$Data.2D$Outlier=="Yes" , ]

if(nrow(Outliers.id) > 0){
  
  message(nrow(Outliers.id) , " outliers found:\n" , paste(rownames(Outliers.id), collapse = ", "))
  
  message("Saving Outlier plots..")
  tiff(filename = paste0(OutPrefix , ".OutlierPlot.2D.tif") , res = 300 , units = "in" , height = 8 , width = 8)
  print(Outliers$Plot.2D)
  graphics.off()
  
  tiff(filename = paste0(OutPrefix , ".OutlierPlot.QQ.tif") , res = 300 , units = "in" , height = 8 , width = 8)
  print(Outliers$Plot.Dist)
  graphics.off()
  
  tiff(filename = paste0(OutPrefix , ".OutlierPlot.Dist.tif") , res = 300 , units = "in" , height = 8 , width = 8)
  print(Outliers$Plot.QQ)
  graphics.off()
  
  message("Removing outliers and updating SO object...")
  pheno.so <- so$sample_to_covariates
  pheno.so <- pheno.so[!(pheno.so$sample %in% rownames(Outliers.id)),]
  so <- sleuth_prep(pheno.so,target_mapping = target_map,extra_bootstrap_summary = TRUE)
  message("Saving updated Sleuth Object in ",SO_file," ...")
  save(so, Outliers,file = SO_file)
  
}else{
  message("There is no outlier sample in the data.")
}


################### Finding DEGs using linear regression ###############
message("DEG analysis..")


if(PCs > 0){
  tpm.norm <- so$obs_norm[, c(1,2,4)]
  tpm.norm <- pivot_wider(tpm.norm, id_cols = "target_id", names_from = "sample", values_from = "tpm")
  tpm.norm <- column_to_rownames(tpm.norm, var = "target_id")
  tpm.norm <- as.matrix(tpm.norm)
  tpm.norm <- t(tpm.norm)
  pca <- prcomp(tpm.norm, rank. = PCs)
  pca <- pca$x
  pca <- as.data.frame(scale(pca))
  pca <- rownames_to_column(pca)
  index <- match(pheno$sample , pca$rowname)
  pheno <- cbind.data.frame(pheno , pca[index,-1])
  so$sample_to_covariates <- pheno
}

so_fit <- sleuth_fit(so, lm_model)

d_matrix.group <- colnames(so_fit$fits$full$design_matrix)[2]
so_fit.wt <- sleuth_wt(so_fit,d_matrix.group)
results_table <- sleuth_results(so_fit.wt, d_matrix.group, test_type = 'wt')
print(head(results_table))

message("Saving results in ",paste0(OutPrefix , ".sleuth.DEG.tsv")," ...")
write.table(results_table , file = paste0(OutPrefix , ".sleuth.DEG.tsv"), quote = F , sep = "\t" , row.names = F)

message("All done!")











