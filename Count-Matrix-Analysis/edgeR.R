message("loading requireed libraries...")
suppressMessages(suppressWarnings(library(DESeq2)))
suppressMessages(suppressWarnings(library(edgeR)))
suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(sva)))
set.seed(12345)
########################################################################
calculate_lambda <- function(pvals) {
  pvals <- pvals[!is.na(pvals) & pvals > 0 & pvals < 1]
  chisq <- qchisq(1 - pvals, df = 1)
  lambda <- median(chisq) / qchisq(0.5, df = 1)
  return(lambda)
}
########################################################################
#
#          Input parameters
#
########################################################################

args = commandArgs(T)

counts.file <- trimws(args[1])
pheno.file <- trimws(args[2])
var.trait <- trimws(args[3])
var.batch.num <- args[4]
var.batch.fact <- args[5]
outliers <- args[6]
gFilter.min.count <- as.numeric(trimws(args[7]))
gFilter.min.prop <- as.numeric(trimws(args[8]))
n.SV <- tolower(trimws(args[9]))
n.PC <- as.numeric(trimws(args[10]))
OutPrefix <- trimws(args[11])

cat("##########################################################################\n")
message("Input arguments:")
message("        Count matrix file: ",counts.file)
message("        Phenotype file: ",pheno.file)
message("        Trait variable: ",var.trait)
message("        Numeric batches: ",var.batch.num)
message("        Categorical batches: ",var.batch.fact)
message("        Outlie samples: ",outliers)
message("        Minimum count threshold for gene filtering: ", gFilter.min.count)
message("        Minimum proportion of the samples for gene filtering: ", gFilter.min.prop)
message("        Number of Sorrogate Variables added to the model: ",n.SV)
message("        Number of Principal Components added to the model: ",n.PC)
message("        Output files prefix: ",OutPrefix)
cat("##########################################################################\n")

########################################################################
#
#          Reading the data
#
########################################################################

message("Reading input data...")

counts <- read.table(counts.file , header = T , row.names = 1 , sep = "\t", stringsAsFactors = F, check.names = F)
counts <- round(counts)
#Assume that the count data are integer. If it is not, you should use textimport. If you can't use teximport. you may have to round the values or use limma instead. 

pheno <- read.csv(pheno.file , row.names = 1 , stringsAsFactors = F)

var.batch.num <- trimws(str_split_1(var.batch.num , pattern = ","))
var.batch.num <- var.batch.num[nchar(var.batch.num) > 0]

var.batch.fact <- trimws(str_split_1(var.batch.fact , pattern = ","))
var.batch.fact <- var.batch.fact[nchar(var.batch.fact) > 0]

outliers <- trimws(str_split_1(outliers , pattern = ","))

if(!identical(colnames(counts) , rownames(pheno))){
  warning("Row names in Phenotype data are not matched with column names in count data. Shared IDs will be considered.")
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
message("filtering low count genes...")
message("Genes that don't have minimum count of ",gFilter.min.count, " in at least ",(gFilter.min.prop*100) , "% of the samples will be removed.")
keep <- edgeR::filterByExpr(counts,group = pheno[,var.trait],min.count = gFilter.min.count, min.prop = gFilter.min.prop)
message(sum(!keep),"/",nrow(counts)," genes removed. Remaining genes:", sum(keep))
counts <- counts[keep,]

########################################################################
#
#          Removing outlier samples
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
  
  #paste("Is count and phenotype data are matched?",ifelse(identical(colnames(counts) , rownames(pheno)),"Yes","NO"))
}

var.all = var.trait
pheno[,var.trait] <- as.factor(pheno[,var.trait])
if(length(var.batch.fact) > 0){
  for (var_ in var.batch.fact) {
    pheno[,var_] <- as.factor(pheno[,var_])
  }
  var.all = c(var.all , var.batch.fact)
}

if(length(var.batch.num) > 0){
  for (var_ in var.batch.num) {
    pheno[,var_] <- as.numeric(pheno[,var_])
  }
  var.all = c(var.all , var.batch.num)
}

########################################################################
#
#          DEG analysis
#
########################################################################

OutPrefix <- paste0(OutPrefix , ".edgeR")

if(n.SV == "all" | (as.numeric(n.SV) > 0)){
  message("Calculating surrogate variables...")
  mod0 <- model.matrix(as.formula(paste0("~",paste(c(var.batch.num,var.batch.fact),collapse = "+"))),data=pheno)
  design.sva <- as.formula(paste0("~",paste(var.all , collapse = "+")))
  message("SVA mod0:\n",mod0 , "\nSVA main model:\n",design.sva)
  mod1 <- model.matrix(design.sva , data = pheno)
  
  counts.norm <- edgeR::cpm(counts , log = T)
  svs = sva(dat = as.matrix(counts.norm),mod = mod1 , mod0 = mod0)$sv
  
  cat("\n")
  colnames(svs) <- paste0("SV", c(1:ncol(svs)))
  print(head(svs))
  
  pheno <- cbind.data.frame(pheno , svs)
  
  if((ncol(svs) < as.numeric(n.SV)) | (n.SV == "all")){
    n.SV = ncol(svs)
  }
  
  var.all <- c(var.all , paste0("SV", c(1:n.SV)))
  
  OutPrefix = paste0(OutPrefix , ".SV",n.SV)
  
}

if(n.PC > 0){
  message("Calculate the Principal Components...")
  counts.norm <- edgeR::cpm(counts , log = T)
  
  pca <- prcomp(t(counts.norm), rank. = n.PC)
  PCs <- pca$x
  PCs <- as.data.frame(scale(PCs))
  pheno <- cbind.data.frame(pheno , PCs)
  
  var.all <- c(var.all , paste0("PC", c(1:n.PC)))
  
  OutPrefix = paste0(OutPrefix , ".PC",n.PC)
}

message("Running DEG analysis using glmFit function in edgeR...")

design.formula = as.formula(paste0("~0+",paste(var.all , collapse = "+")))

message("Linea regression model:\n",design.formula)

design.matrix = model.matrix(design.formula , data = pheno)
groups <- levels(as.factor(pheno[,var.trait]))
colnames(design.matrix)[1:length(groups)] <- groups

dge.list <- DGEList(counts = counts,samples = pheno , group = pheno[,var.trait])

dge.list <- calcNormFactors(dge.list)

dge.list <- estimateDisp(dge.list , design = design.matrix , robust = T)

fit <- edgeR::glmFit(dge.list , design = design.matrix)

contrasts_  <- t(combn(groups, 2))

for (i in 1:nrow(contrasts_)){
  
  out_name = paste0(OutPrefix,".",paste(contrasts_[i,], collapse = "."))
  
  lrt <- edgeR::glmLRT(fit, contrast = makeContrasts(contrasts = paste(contrasts_[i,] , collapse = "-") , levels = design.matrix))
  result = topTags(lrt,adjust.method = "BH" ,n = nrow(dge.list))
  result = result$table
  if(is.null(result)){
    result = as.data.frame(matrix(data = NA , nrow = 1 , ncol = 5))
  }
  
  colnames(result) = c("logFC" , "logCPM", "LR"   ,  "PValue" ,"adjPValue.BH")
  result$adjPvalue.bnf <- p.adjust(result$PValue , method = "bonferroni")
  result <- result[order(result$PValue, decreasing = F),]
  result <- cbind.data.frame(Gene = rownames(result) , result)
  
  inflation = calculate_lambda(result$PValue)
  message("The inflation index of the results (Lambda) is ", inflation)
  
  pdf(file = paste0(out_name , ".QQ.pdf") , height = 8 , width = 8)
  qq(result$PValue, main=paste0("QQ Plot - " , paste(contrasts_[i,], collapse = " vs ")))
  text(x = 0.5,y = (par("usr")[4]-0.2),
       label = bquote(lambda == .(round(inflation, 2))),
       adj = c(0, 1),
       cex = 1)
  graphics.off()
  
  write.csv(result , file = paste0(out_name , ".csv"), row.names = F)
}
