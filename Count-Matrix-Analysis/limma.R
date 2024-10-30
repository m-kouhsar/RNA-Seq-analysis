library(DESeq2)
library(edgeR)
library(limma)
library(stringr)
library(sva)
set.seed(12345)
########################################################################
#
#          Input parameters
#
########################################################################

setwd("C:/Users/mk693/OneDrive - University of Exeter/Desktop/2021/NIH/Data/miRNA-Seq/Sep2024")
counts.file <- "Raw/Project.11008.PITT.miRNA.AD.C.Psy2.Count.txt"
pheno.file <- "Raw/Project.11008.PITT.miRNA.AD.C.Psy2.Pheno.csv"
var.trait <- "Trait"
var.batch.num <- "Age,RIN"
var.batch.fact <- "Sex,Plate"
outliers = "CW12006,CW03001,CW98019,CW99011" #"BBN00226278,BBN20018,BBN00228673"
logFC = round(log2(1.2) , digits = 2)
Pvalue = 1.1
P.adjust = 1.1
p.adjust.method = "bonferroni"
runSVA = T
n.SV = 4
OutPrefix <- "Results/PITT-ADRC/PITT.miRNA.AD.C.Psy2.limma"

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

for (var_ in var.batch.num) {
  pheno[,var_] <- as.numeric(pheno[,var_])
}

########################################################################
#
#          Balancing the number of samples in groups
#
########################################################################
pheno <- pheno[c(rownames(pheno)[pheno$Trait == "C"],
                 rownames(pheno)[pheno$Trait == "AD"],
                 sample(rownames(pheno)[pheno$Trait=="ADP2"],size = 36,replace = F)),]
table(pheno$Trait)
index = match(rownames(pheno) , colnames(counts))
counts = counts[,index]
identical(rownames(pheno) , colnames(counts))

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

fit <- voom(dge.list , design = design.matrix , plot = F)
fit <- lmFit(fit , design = design.matrix)

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
    
    fit <- voom(dge.list , design = design.matrix , plot = F)
    fit <- lmFit(fit , design = design.matrix)
    contrast = makeContrasts(contrasts = paste(contrasts_[i,] , collapse = "-") , levels = design.matrix)
    fit <- contrasts.fit(fit , contrast)
    fit <- eBayes(fit)
    result = topTable(fit,adjust.method = p.adjust.method ,n = nrow(dge.list))
    SE = sqrt(fit$s2.post)*fit$stdev.unscaled
    index = match(rownames(result) , rownames(SE))
    result$SE = SE[index]
    if(is.null(result)){
      result = as.data.frame(matrix(data = NA , nrow = 1 , ncol = 6))
      colnames(result) = c("logFC" , "AveExpr", "t"   ,  "PValue" ,"adj.P.Val" , "B")
    }
    
  }else{
    
    
    contrast = makeContrasts(contrasts = paste(contrasts_[i,] , collapse = "-") , levels = design.matrix)
    fit <- contrasts.fit(fit , contrast)
    fit <- eBayes(fit)
    result = topTable(fit,adjust.method = p.adjust.method ,n = nrow(dge.list))
    SE = sqrt(fit$s2.post)*fit$stdev.unscaled
    index = match(rownames(result) , rownames(SE))
    result$SE = SE[index]
    if(is.null(result)){
      result = as.data.frame(matrix(data = NA , nrow = 1 , ncol = 6))
      colnames(result) = c("logFC" , "AveExpr", "t"   ,  "P.Value" ,"adj.P.Val" , "B")
    }
    
    
  }
  
  result.filter <- result[(abs(result$logFC) > logFC ) & (result$P.Value < Pvalue) & (result$adj.P.Val < P.adjust),]
  write.csv(result.filter , file = paste0(OutPrefix,".",paste(contrasts_[i,], collapse = "."),".logFC.",logFC,
                                          ".Pval.",Pvalue,".",p.adjust.method,".",P.adjust,".csv"))
}

########################################################################
#
#          Visualization
#
########################################################################


