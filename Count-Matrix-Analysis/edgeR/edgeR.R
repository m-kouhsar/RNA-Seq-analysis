library(DESeq2)
library(edgeR)
library(stringr)
library(sva)
set.seed(12345)

########################################################################
#
#          Input parameters
#
########################################################################

setwd("C:/Users/mk693/OneDrive - University of Exeter/Desktop/2021/NIH/Data/miRNA-Seq/Sep2024")
counts.file <- "Raw/BDR.miRNA.AD.C.Psy2.Count.txt"
pheno.file <- "Raw/BDR.miRNA.AD.C.Psy2.Pheno.csv"
var.trait <- "Trait"
var.batch.num <- "Age,RIN"
var.batch.fact <- "Sex,Plate"
outliers = "BBN10205"
logFC = round(log2(1.2) , digits = 2)
Pvalue = 0.05
P.adjust = 0.05
p.adjust.method = "bonferroni"
runSVA = T
n.SV = 3
OutPrefix <- "Results/BDR/BDR.miRNA.AD.C.Psy2.edgeR"

########################################################################
#
#          Reading the data
#
########################################################################

counts <- read.table(counts.file , header = T , row.names = 1 , sep = "\t", stringsAsFactors = F)
pheno <- read.csv(pheno.file , row.names = 1 , stringsAsFactors = F)
var.batch.num <- trimws(str_split_1(var.batch.num , pattern = ","))
var.batch.fact <- trimws(str_split_1(var.batch.fact , pattern = ","))
outliers <- trimws(str_split_1(outliers , pattern = ","))

if(!identical(colnames(counts) , rownames(pheno))){
  warning("Row names in Phenotype data are not matched with column names in count data!")
  index <- intersect(rownames(pheno) , colnames(counts))
  if(length(index)<2){
    stop("Phenotype data and counts data cannot be matched!")
  }else{
    pheno <- pheno[index , ]
    counts <- counts[, index]
  }
}

for (var_ in var.batch.num) {
  pheno[,var_] <- as.numeric(pheno[,var_])
}

########################################################################
#
#          Removing outliers
#
########################################################################

if(all(outliers != "")){
  message("Removing outliers...")
  if(! all (outliers %in% colnames(counts))){
    warning("The followin outlier IDs are not exist in the data: ")
    paste(outliers[!(outliers %in% colnames(counts))] , collapse = ";")
  }
  counts <- counts[,!(colnames(counts) %in% outliers)]
  pheno <- pheno[!(rownames(pheno) %in% outliers),]
  
  paste("Is count and phenotype data are matched?",ifelse(identical(colnames(counts) , rownames(pheno)),"Yes","NO"))
}

for (var_ in var.batch.fact) {
  pheno[,var_] <- as.factor(pheno[,var_])
}

########################################################################
#
#          Filtering low counts
#
########################################################################
message("Filtering low count genes...")
keep <- edgeR::filterByExpr(counts,group = pheno[,var.trait],min.count = 10)
message(sum(keep), " genes remained.")
counts <- counts[keep,]

########################################################################
#
#          DEG analysis
#
########################################################################
dge.list <- DGEList(counts = counts,samples = pheno)

dge.list <- calcNormFactors(dge.list)

design.formula = as.formula(paste0("~0+",var.trait,"+",paste(c(var.batch.fact , var.batch.num ) , collapse = "+")))
design.matrix = model.matrix(design.formula , data = pheno)
groups <- levels(as.factor(pheno[,var.trait]))
colnames(design.matrix)[1:length(groups)] <- groups

dge.list <- estimateDisp(dge.list , design = design.matrix , robust = T)

fit <- edgeR::glmFit(dge.list , design = design.matrix)

contrasts_  <- t(combn(groups, 2))

for (i in 1:nrow(contrasts_)){
  
  if(runSVA){
    pheno.1 <- pheno[pheno$Trait %in% contrasts_[i,],]
    counts.1 <- counts[,pheno$Trait %in% contrasts_[i,]]
    mod0 <- model.matrix(~1,data=pheno.1)
    design.sva <- as.formula(paste0("~",paste(c(var.batch.fact , var.batch.num ) , collapse = "+"),"+",var.trait))
    mod1 <- model.matrix(design.sva , data = pheno.1)
    
    svs = sva(dat = as.matrix(counts.1),mod = mod1 , mod0 = mod0)$sv
    
    colnames(svs) <- paste0("SV", c(1:ncol(svs)))
    
    pheno.1 <- cbind.data.frame(pheno.1 , svs)
    
    if(ncol(svs) < n.SV){
      n.SV = ncol(svs)
    }
    var.batch.all <- c(var.batch.fact , var.batch.num , paste0("SV", c(1:n.SV)))
    
    dge.list <- DGEList(counts = counts.1,samples = pheno.1)
    
    dge.list <- calcNormFactors(dge.list)
    
    design.formula = as.formula(paste0("~0+",var.trait,"+",paste(var.batch.all, collapse = "+")))
    design.matrix = model.matrix(design.formula , data = pheno.1)
    groups <- levels(as.factor(pheno.1[,var.trait]))
    colnames(design.matrix)[1:length(groups)] <- groups
    
    dge.list <- estimateDisp(dge.list , design = design.matrix , robust = T)
    
    fit <- edgeR::glmFit(dge.list , design = design.matrix)
    
    lrt <- edgeR::glmLRT(fit, contrast = makeContrasts(contrasts = paste(contrasts_[i,] , collapse = "-") , levels = design.matrix))
    result = topTags(lrt,adjust.method = p.adjust.method , p.value = P.adjust,n = nrow(dge.list))
    result = result$table
    if(is.null(result)){
      result = as.data.frame(matrix(data = NA , nrow = 1 , ncol = 5))
      colnames(result) = c("logFC" , "logCPM", "LR"   ,  "PValue" ,"FWER")
    }else{
      result.filter <- result[(abs(result$logCPM) > logFC ) & (result$PValue < Pvalue),]
    }
  
    }else{
    
    
    lrt <- edgeR::glmLRT(fit, contrast = makeContrasts(contrasts = paste(contrasts_[i,] , collapse = "-") , levels = design.matrix))
    result = topTags(lrt,adjust.method = p.adjust.method , p.value = P.adjust,n = nrow(dge.list))
    result = result$table
    if(is.null(result)){
      result = as.data.frame(matrix(data = NA , nrow = 1 , ncol = 5))
      colnames(result) = c("logFC" , "logCPM", "LR"   ,  "PValue" ,"FWER")
    }else{
      result.filter <- result[(abs(result$logCPM) > logFC ) & (result$PValue < Pvalue),]
    }
    
    
  }
  
  write.csv(result.filter , file = paste0(OutPrefix,".",paste(contrasts_[i,], collapse = "."),".logFC.",logFC,
                                          ".Pval.",Pvalue,".",p.adjust.method,".",P.adjust,".csv"))
}
