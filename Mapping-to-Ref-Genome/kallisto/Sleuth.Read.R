message("#########################################################")
message("kallisto pipeline version 1.0.0")
message("Wrote by m.kouhsar@exeter.ac.uk")
message("R script for reading kallisto results using the Sleuth package.")
message("A Sleuth object will be generated and saved.")
message("#########################################################")

suppressMessages(library(stringr))
suppressMessages(library(sleuth))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(tibble))
suppressMessages(library(qqman))

#################### Input Arguments ####################

# kallisto_res_dir: Directory of kallisto results to read by sleuth
# pheno_file: is a csv file which must contains the following columns:
#             sample: Samples ID
#             path: path of the kallisto resuls (inside kallisto_res_dir) folder for each sample
#             All factor and numeric variables must be represented by a column with the same name
# target_map_file: Target mapping file contains information (eg. gene ID and type) about the transcripts
# var_factor: Factor variable in correlation plot
# var_numeric : Numerical varaibles in correlation plot
# RemoveOutliers: Do you want to remove outlier samples from the data? (set it to 'yes' or 'no')
# OutPrefix: Results files/images prefix (can contains a directory)
# ScriptDir: Directory of all Scripts related to this analysis 

####################################################################
Arguments <- commandArgs(T)

kallisto_res_dir <- trimws(Arguments[1])
pheno_file <- trimws(Arguments[2]) 
target_map_file <- trimws(Arguments[3]) 
var_factor <- trimws(Arguments[4])    
var_numeric <- trimws(Arguments[5]) 
RemoveOutliers <- tolower(trimws(Arguments[6]))
OutPrefix <- trimws(Arguments[7])
ScriptDir <- trimws(Arguments[8])

if(is.na(OutPrefix)){
  OutPrefix <- ""
}


source(paste0(ScriptDir , "/mahalanobis.outlier.R"))
source(paste0(ScriptDir , "/CovariatePlot.R"))

var_factor <- str_split_1(var_factor , pattern = ",")
var_numeric <- str_split_1(var_numeric , pattern = ",")

message("Input arguments:")
message("        kallisto results directory: ",kallisto_res_dir)
message("        Phenotype file: ",pheno_file)
message("        Target mapping file: ",target_map_file)
message("        Factor variables in the model: ",paste(var_factor, collapse = ", "))
message("        Numeric variables in the model: ",paste(var_numeric,collapse=", "))
message("        Removing Outliers: ",RemoveOutliers)
message("        Output files prefix: ",OutPrefix)
message("        Script directory: ",ScriptDir)

message("#########################################################")

dir.create(path = dirname(OutPrefix),recursive = T,showWarnings = F)

pheno <- read.csv(pheno_file , stringsAsFactors = F)

if(!all(c("sample","path",var_factor , var_numeric) %in% colnames(pheno))){
  
  stop("The following variables can't be find in phenotype file: \n",
       paste(setdiff(c("sample","path",var_factor , var_numeric), colnames(pheno)),collapse = ", "))
}

for (v in var_factor) {
  pheno[,v] <- as.factor(pheno[,v])
}
for (v in var_numeric) {
  pheno[,v] <- as.numeric(pheno[,v])
}
pheno$path <- paste(kallisto_res_dir , pheno$path , sep = "/")

################### Reading kallisto results #################

if(target_map_file=="" | is.na(target_map_file)){
  message("No target mapping file provided.")
  target_map <- NULL
}else{
  target_map <- read.delim(target_map_file, stringsAsFactors = F , header = T , sep = "\t")
}


SO_file <- paste0(OutPrefix , ".SO.rdat")

if(file.exists(SO_file)){
  
  message("Loading Sleuth Object from ",SO_file, " ...")
  load(SO_file)
  
}else{
  
  
  so <- sleuth_prep(pheno,target_mapping = target_map,extra_bootstrap_summary = TRUE)
  
  message("Saving Sleuth Object in ",SO_file," ...")
  save(so,file = SO_file)
}

tpm.norm <- so$obs_norm[, c(1,2,4)]
tpm.norm <- pivot_wider(tpm.norm, id_cols = "target_id", names_from = "sample", values_from = "tpm")
tpm.norm <- column_to_rownames(tpm.norm, var = "target_id")
pheno.so <- so$sample_to_covariates
index <- match(pheno.so$sample , colnames(tpm.norm))
tpm.norm <- tpm.norm[,index]
rownames(pheno.so) <- pheno.so$sample


P <- CovariatePlot(Data = tpm.norm , Phenotype = pheno.so , Factor_var = var_factor , 
                   Numeric_var = var_numeric,Plot_titel = "Correlation Plot" )
tiff(filename = paste0(OutPrefix , ".PC.Corr.tif") , res = 300 , units = "in" , height = 8 , width = 8)
print(P)
graphics.off()

if(RemoveOutliers=="yes"){
  
  ################### Removing Outliers ##################################
  
  message("Finding Outlier samples using Mahalanobis distance and chi-squared distribution...")
  
  Outliers <- mahalanobis.outlier(Data=tpm.norm , method="pca",pca.scale=F)
  
  Outliers.id <- Outliers$Data.2D[Outliers$Data.2D$Outlier=="Yes" , ]
  
  if(nrow(Outliers.id) > 0){
    
    message(nrow(Outliers.id) , " outliers found:\n" , paste(rownames(Outliers.id), collapse = ", "))
    
    message("Saving Outlier plots..")
    tiff(filename = paste0(OutPrefix , ".OutlierPlot.2D.tif") , res = 300 , units = "in" , height = 8 , width = 8)
    print(Outliers$Plot.2D)
    graphics.off()
    
    tiff(filename = paste0(OutPrefix , ".OutlierPlot.QQ.tif") , res = 300 , units = "in" , height = 8 , width = 8)
    print(Outliers$Plot.QQ)
    graphics.off()
    
    tiff(filename = paste0(OutPrefix , ".OutlierPlot.Dist.tif") , res = 300 , units = "in" , height = 8 , width = 8)
    print(Outliers$Plot.Dist)
    graphics.off()
    
    message("Removing outliers and updating SO object...")
    
    pheno.so <- pheno.so[!(pheno.so$sample %in% rownames(Outliers.id)),]
    index <- match(pheno.so$sample , pheno$sample)
    pheno.so$path <- pheno$path[index]
    
    so <- sleuth_prep(pheno.so,target_mapping = target_map,extra_bootstrap_summary = TRUE)
    message("Saving updated Sleuth Object in ",SO_file," ...")
    save(so, Outliers,file = SO_file)
    
  }else{
    message("There is no outlier sample in the data.")
  }
  
}

message("All done!")











