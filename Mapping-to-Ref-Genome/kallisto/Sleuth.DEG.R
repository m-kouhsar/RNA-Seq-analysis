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
suppressMessages(library(qqman))

#################### Input Arguments ####################

# SO_file: Sleuth Object file (created by 1.Sleuth.Read.R)
# lm_model: linear regression model to run the test
# var_factor: Factor variable include condition variable in lm_model
# var_numeric : Numerical varaibles in lm_model
# PCs: Number of principal components you want to add to the analysis as covariates
# OutPrefix: Results files/images prefix (can contains a directory)

####################################################################
Arguments <- commandArgs(T)

SO_file <- trimws(Arguments[1])
lm_model <- trimws(Arguments[2])
var_factor <- trimws(Arguments[3])    
var_numeric <- trimws(Arguments[4])  
PCs <- as.numeric(trimws(Arguments[5]))
save_fit <- trimws(Arguments[6])
OutPrefix <- trimws(Arguments[7])

if(is.na(OutPrefix)){
  OutPrefix <- ""
}
var_factor <- trimws(str_split_1(var_factor , pattern = ","))
var_numeric <- trimws(str_split_1(var_numeric , pattern = ","))
if(tolower(save_fit) == "yes"){
  save_fit = T
}else{
  save_fit=F
}

message("Input arguments:")
message("        Sleuth Object file: ",SO_file)
message("        Linear regression model: ",lm_model)
message("        Factor variables in the model: ",paste(var_factor, collapse = ", "))
message("        Numeric variables in the model: ",paste(var_numeric,collapse=", "))
message("        Number of PCs to add to the model: ",PCs)
message("        Save fit object: ",save_fit)
message("        Output files prefix: ",OutPrefix)

message("#########################################################")

dir.create(path = dirname(OutPrefix),recursive = T,showWarnings = F)

if(!all(all.vars(lm_model) %in% c(var_factor,var_numeric))){
  stop("The following variables in the regression model are not specified as factor or numeric variables:\n",
       paste(setdiff(all.vars(lm_model) , c(var_factor,var_numeric)),collapse = ", "))
}

message("Loading Sleuth Object from ",SO_file, " ...")
load(SO_file)

pheno.so <- so$sample_to_covariates

if(!all(c("sample",all.vars(lm_model)) %in% colnames(pheno.so))){
  
  stop("The following variables can't be find in phenotype data in SO file: \n",
       paste(setdiff(c("sample",all.vars(lm_model)), colnames(pheno.so)),collapse = ", "))
}

message("DEG analysis..")

if(PCs > 0){
  tpm.norm <- so$obs_norm[, c(1,2,4)]
  tpm.norm <- pivot_wider(tpm.norm, id_cols = "target_id", names_from = "sample", values_from = "tpm")
  tpm.norm <- column_to_rownames(tpm.norm, var = "target_id")
  pca <- prcomp(t(tpm.norm), rank. = PCs)
  pca <- pca$x
  pca <- as.data.frame(scale(pca))
  pca <- rownames_to_column(pca)
  index <- match(pheno.so$sample , pca$rowname)
  pheno.so <- cbind.data.frame(pheno.so , pca[index,-1])
  so$sample_to_covariates <- pheno.so
  lm_model <- paste0(lm_model ,"+", paste0("PC",1:PCs, collapse = "+"))
}

lm_model <- as.formula(lm_model)
so_fit <- sleuth_fit(so, lm_model)

if(save_fit){
  message("Saving sleuth fit object to ",paste0(OutPrefix , ".sleuth.DEG.tsv"),"...")
  save(so_fit ,file = paste0(OutPrefix , ".sleuth.Fit.rdat") )
}

d_matrix.group <- colnames(so_fit$fits$full$design_matrix)[2]
message("which_beta argument in sleuth_wt function is set to ",d_matrix.group)
so_fit.wt <- sleuth_wt(so_fit,d_matrix.group)
results_table <- sleuth_results(so_fit.wt, d_matrix.group, test_type = 'wt')

################################# Checking the results #############################
print(head(results_table))

chisq <- qchisq(1-results_table$pval,1)
inflation <- median(chisq, na.rm = T)/qchisq(0.5,1)
message("Inflation index: ",round(inflation, digits = 2))

tiff(filename = paste0(OutPrefix , ".sleuth.DEG.QQ.tif") , res = 300 , units = "in" , height = 8 , width = 8)
qq(results_table$pval, main="QQ plot")
text(x = 0.5,y = (par("usr")[4]-0.2),
     label = bquote(lambda == .(round(inflation, 2))),
     adj = c(0, 1),
     cex = 1)
graphics.off()

message("Saving results in ",paste0(OutPrefix , ".sleuth.DEG.tsv")," ...")
write.table(results_table , file = paste0(OutPrefix , ".sleuth.DEG.tsv"), quote = F , sep = "\t" , row.names = F)

message("All done!")
