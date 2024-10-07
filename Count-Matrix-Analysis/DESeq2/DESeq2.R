library(DESeq2)
library(edgeR)
library(stringr)
library(sva)
set.seed(12345)

####################################################################################################
DEG.DESeq2 <- function(count.data , phenotype.data , trait , batches, p.adjust.method = "BH"){
  
  result.table <- NA
  
  design_ <- as.formula(paste0("~",paste(batches , collapse = "+"),"+",trait))
  
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = count.data , colData = phenotype.data , design = design_ )
  
  dds <- DESeq(dds)
  
  result.table <- as.data.frame(results(dds,contrast = c(trait , unique(phenotype.data[,trait])[1:2])))
  result.table$padj <- p.adjust(result.table$pvalue , method = p.adjust.method)
  result.table <- result.table[order(result.table$pvalue, decreasing = F),]
  
  return(result.table)
}

########################################################################
#
#          Input parameters
#
########################################################################

setwd("./")
counts.file <- "Raw/Count.txt"
pheno.file <- "Raw/Pheno.csv"
var.trait <- "Trait"
var.batch.num <- "Age,RIN"
var.batch.fact <- "Sex,Plate"
outliers = "Sample1,Sample2"
logFC = 1
Pvalue = 0.05
P.adjust = 0.05
p.adjust.method = "bonferroni"
runSVA = T
n.SV = 3
OutPrefix <- "Results/Project01"

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

for (var_ in var.batch.fact) {
  pheno[,var_] <- as.factor(pheno[,var_])
}

for (var_ in var.batch.num) {
  pheno[,var_] <- as.numeric(pheno[,var_])
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

########################################################################
#
#          DEG analysis
#
########################################################################

groups <- unique(pheno[,var.trait])
contrasts_  <- t(combn(groups, 2))

for (i in 1:nrow(contrasts_)){
  
  pheno.1 <- pheno[pheno$Trait %in% contrasts_[i,],]
  counts.1 <- counts[,pheno$Trait %in% contrasts_[i,]]
  
  if(runSVA){
    mod0 <- model.matrix(~1,data=pheno.1)
    design_ <- as.formula(paste0("~",paste(c(var.batch.fact , var.batch.num ) , collapse = "+"),"+",var.trait))
    mod1 <- model.matrix(design_ , data = pheno.1)
    
    svs = sva(dat = as.matrix(counts.1),mod = mod1 , mod0 = mod0)$sv
    
    colnames(svs) <- paste0("SV", c(1:ncol(svs)))
    
    pheno.1 <- cbind.data.frame(pheno.1 , svs)
    
    if(ncol(svs) < n.SV){
      n.SV = ncol(svs)
    }
    
    var.batch.all <- c(var.batch.fact , var.batch.num , paste0("SV", c(1:n.SV)))
    result <- as.data.frame(DEG.DESeq2(count.data = counts.1 , phenotype.data = pheno.1,trait = var.trait , 
                         batches = var.batch.all,p.adjust.method = p.adjust.method))
    
  }else{
    
    var.batch.all <- c(var.batch.fact , var.batch.num )
    result <- DEG.DESeq2(count.data = counts.1 , phenotype.data = pheno.1,trait = var.trait , 
                         batches = var.batch.all,p.adjust.method = p.adjust.method)
    
  }
  
  result.filter <- result[(abs(result$log2FoldChange) > logFC ) & (result$pvalue < Pvalue) & (result$padj < P.adjust),]
  write.csv(result.filter , file = paste0(OutPrefix,".",paste(contrasts_[i,], collapse = "."),".logFC.",logFC,
                                   ".Pval.",Pvalue,".",p.adjust.method,".",P.adjust,".csv"))
}








































