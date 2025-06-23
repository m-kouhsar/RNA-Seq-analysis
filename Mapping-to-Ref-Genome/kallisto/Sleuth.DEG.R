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
suppressMessages(library(sva))

#################### Input Arguments ####################

# SO_file: Sleuth Object file (created by 1.Sleuth.Read.R)
# lm_model: linear regression model to run the test
# PCs: Number of principal components you want to add to the model
# SVs: Number of surrogate variables you want to add to the model
# save_updated_so: Would you like to save SO (and SVA) object after running test? (yes/no)
# OutPrefix: Results files/images prefix (can contains a directory)

####################################################################
Arguments <- commandArgs(T)

SO_file <- trimws(Arguments[1])
lm_model <- trimws(Arguments[2])
PCs <- as.numeric(trimws(Arguments[3]))
SVs <- as.numeric(trimws(Arguments[4]))
save_updated_so <- trimws(Arguments[5])
OutPrefix <- trimws(Arguments[6])

if(is.na(OutPrefix)){
  OutPrefix <- ""
}

message("Input arguments:")
message("        Sleuth Object file: ",SO_file)
message("        Linear regression model: ",lm_model)
message("        Number of PCs to add to the model: ",PCs)
message("        Number of Surrogate Variable to add to the model: ",SVs)
message("        Save updated SO object after test: ",save_updated_so)
message("        Output files prefix: ",OutPrefix)

message("#########################################################")

dir.create(path = dirname(OutPrefix),recursive = T,showWarnings = F)

message("Loading Sleuth Object from ",SO_file, " ...")
load(SO_file)

pheno.so <- so$sample_to_covariates

if(!all(c("sample",all.vars(lm_model)) %in% colnames(pheno.so))){
  
  stop("The following variables can't be find in phenotype data in SO file: \n",
       paste(setdiff(c("sample",all.vars(lm_model)), colnames(pheno.so)),collapse = ", "))
}

tpm.norm <- so$obs_norm[, c(1,2,4)]
tpm.norm <- pivot_wider(tpm.norm, id_cols = "target_id", names_from = "sample", values_from = "tpm")
tpm.norm <- column_to_rownames(tpm.norm, var = "target_id")

if(PCs > 0){
  message("Calculate the Principal Components...")
  pca <- prcomp(t(tpm.norm), rank. = PCs)
  pca <- pca$x
  pca <- as.data.frame(scale(pca))
  pca <- rownames_to_column(pca)
  message("Adding PCs to the model...")
  index <- match(pheno.so$sample , pca$rowname)
  pheno.so <- cbind.data.frame(pheno.so , pca[index,-1])
  so$sample_to_covariates <- pheno.so
  lm_model <- paste0(lm_model ,"+", paste0("PC",1:PCs, collapse = "+"))
}

if(SVs > 0){
  if(!exists("sva_df")){
    
    print("Calculating Surrogate Variables...")
    mod1 = model.matrix(as.formula(lm_model), data=pheno.so)
    mod0 = cbind(mod1[,1])
    index = apply(tpm.norm , 1 , function(x){all(x==0)})
    tpm.norm <- tpm.norm[!index , ]
    set.seed(12345)
    sva_obj= sva::sva(as.matrix(tpm.norm),mod1,mod0)
  }
  message("Adding SVs to the model...")
  sva_df = as.data.frame(sva_obj$sv)
  names(sva_df) <- paste0("SV",1:dim(sva_df)[2])
  rownames(sva_df) <- colnames(tpm.norm)
  index <- match(pheno.so$sample , rownames(sva_df))
  sva_df = sva_df[index , , drop = FALSE]
  pheno.so <- cbind.data.frame(pheno.so , sva_df)
  so$sample_to_covariates <- pheno.so
  if(SVs > ncol(sva_df))
    SVs = ncol(sva_df)
  
  lm_model <- paste0(lm_model ,"+", paste0("SV",1:SVs, collapse = "+"))
}

lm_model <- as.formula(lm_model)

d_matrix <- model.matrix(lm_model , data = pheno.so)
print(head(d_matrix))
so <- sleuth_fit(so, d_matrix)

d_matrix.group <- colnames(so$fits$full$design_matrix)[2]
message("which_beta argument in sleuth_wt function is set to ",d_matrix.group)
so <- sleuth_wt(so,d_matrix.group)
results_table <- sleuth_results(so, d_matrix.group, test_type = 'wt')

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

message("Saving DEG summary statistics in ",paste0(OutPrefix , ".sleuth.DEG.tsv")," ...")
write.table(results_table , file = paste0(OutPrefix , ".sleuth.DEG.tsv"), quote = F , sep = "\t" , row.names = F)

if(tolower(save_updated_so)=="yes"){
  message("Saving updated SO object in ",SO_file," ...")
  
  possible_objects <- c("so", "sva_obj", "Outliers")
  possible_objects <- possible_objects[sapply(possible_objects, exists)]
  
  save(list = possible_objects , file = SO_file)
}

message("All done!")
