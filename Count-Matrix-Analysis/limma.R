message("loading requireed libraries...")
suppressMessages(suppressWarnings(library(DESeq2)))
suppressMessages(suppressWarnings(library(edgeR)))
suppressMessages(suppressWarnings(library(limma)))
suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(sva)))
set.seed(12345)

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
n.SV <- as.numeric(trimws(args[9]))
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
pheno <- read.csv(pheno.file , row.names = 1 , stringsAsFactors = F)
var.batch.num <- trimws(str_split_1(var.batch.num , pattern = ","))
var.batch.fact <- trimws(str_split_1(var.batch.fact , pattern = ","))
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
keep <- edgeR::filterByExpr(counts,group = pheno[,var.trait],min.count = 5, min.prop = 0.75)
message(sum(!keep),"/",nrow(counts)," genes removed. Remaining genes:", sum(keep))
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
  
  #paste("Is count and phenotype data are matched?",ifelse(identical(colnames(counts) , rownames(pheno)),"Yes","NO"))
}

for (var_ in var.batch.fact) {
  pheno[,var_] <- as.factor(pheno[,var_])
}

for (var_ in var.batch.num) {
  pheno[,var_] <- as.numeric(pheno[,var_])
}

########################################################################
#
#          DEG analysis
#
########################################################################

OutPrefix <- paste0(OutPrefix , ".limma")
var.batch.all <- c(var.batch.fact , var.batch.num)
if(n.SV > 0){
  message("Calculating sorrogate variables...")
  mod0 <- model.matrix(~1,data=pheno)
  design.sva <- as.formula(paste0("~",var.trait,"+",paste(c(var.batch.fact , var.batch.num ) , collapse = "+")))
  message("SVA model:\n",design.sva)
  mod1 <- model.matrix(design.sva , data = pheno)
  counts.norm <- edgeR::cpm(counts , log = T)
  svs = sva(dat = as.matrix(counts.norm),mod = mod1 , mod0 = mod0)$sv
  
  colnames(svs) <- paste0("SV", c(1:ncol(svs)))
  
  pheno <- cbind.data.frame(pheno , svs)
  
  if(ncol(svs) < n.SV){
    n.SV = ncol(svs)
  }
  
  var.batch.all <- c(var.batch.all , paste0("SV", c(1:n.SV)))
  head(pheno)
  OutPrefix = paste0(OutPrefix , ".SV",n.SV)
  
}

if(n.PC > 0){
  message("Calculate the Principal Components...")
  counts.norm <- edgeR::cpm(counts , log = T)
  
  pca <- prcomp(t(tpm.norm), rank. = PCs)
  PCs <- pca$x
  PCs <- as.data.frame(scale(PCs))
  pheno <- cbind.data.frame(pheno , PCs)
  
  var.batch.all <- c(var.batch.all , paste0("PC", c(1:n.SV)))
  head(pheno)
  OutPrefix = paste0(OutPrefix , ".PC",n.SV)
}

message("Running DEG analysis using voom and lmfit functions in limma...")

design.formula = as.formula(paste0("~0+",var.trait,"+",paste(var.batch.all , collapse = "+")))
message("Linea regression model:\n",design.formula)

design.matrix = model.matrix(design.formula , data = pheno)
groups <- levels(as.factor(pheno[,var.trait]))
colnames(design.matrix)[1:length(groups)] <- groups

dge.list <- DGEList(counts = counts,samples = pheno , group = pheno[,var.trait])
dge.list <- calcNormFactors(dge.list)

fit.voom <- voom(dge.list , design = design.matrix , plot = F)
fit.lm <- lmFit(fit.voom , design = design.matrix)

contrasts_  <- t(combn(groups, 2))

for (i in 1:nrow(contrasts_)){
  out_name = paste0(OutPrefix,".",paste(contrasts_[i,], collapse = "."))

  contrast = makeContrasts(contrasts = paste(contrasts_[i,] , collapse = "-") , levels = design.matrix)
  fit.contrast <- contrasts.fit(fit.lm , contrast)
  fit.ebays <- eBayes(fit.contrast)
  result = topTable(fit.ebays ,n = nrow(dge.list))
  SE = sqrt(fit.ebays$s2.post)*fit.ebays$stdev.unscaled
  index = match(rownames(result) , rownames(SE))
  result$SE = SE[index]
  if(is.null(result)){
    result = as.data.frame(matrix(data = NA , nrow = 1 , ncol = 6))
  }
  
  colnames(result) = c("logFC" , "AveExpr", "t"   ,  "PValue" ,"adjPvalue.BH" , "B","SE")
  result$adjPvalue.bnf <- p.adjust(result$PValue , method = "bonferroni")
  result <- result[order(result$PValue, decreasing = F),]
  result <- cbind.data.frame(Gene = rownames(result) , result)
  write.csv(result , file = paste0(out_name , ".csv"), row.names = F)
}
