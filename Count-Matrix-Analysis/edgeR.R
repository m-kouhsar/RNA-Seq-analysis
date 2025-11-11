message("loading requireed libraries...")
suppressMessages(suppressWarnings(library(DESeq2)))
suppressMessages(suppressWarnings(library(edgeR)))
suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(sva)))
set.seed(12345)

########################################################################
#
#          Input parameters
#
########################################################################

args = commandArgs(T)

counts.file <- args[1]
pheno.file <- args[2]
var.trait <- args[3]
var.batch.num <- args[4]
var.batch.fact <- args[5]
outliers <- args[6]
gFilter.min.count <- as.numeric(trimws(args[7]))
gFilter.min.prop <- as.numeric(trimws(args[8]))
runSVA <- args[9]
n.SV <- args[10]
OutPrefix <- args[11]

cat("##########################################################################\n")
message("Input arguments:")
message("        Count matrix file: ",counts.file)
message("        Phenotype file: ",pheno.file)
message("        Trait variable: ",var.trait)
message("        Numeric batches: ",var.batch.num)
message("        Categorical batches: ",var.batch.fact)
message("        Outlie samples: ",outliers)
message("        Do you want to add sorrogate variables to the model? ",runSVA)
message("        Number of sorrogate variables: ",n.SV)
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
runSVA = ifelse(trimws(tolower(runSVA))=="yes",T ,F)

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

OutPrefix <- paste0(OutPrefix , ".edgeR")

if(runSVA){
  message("Calculating surrogate variables...")
  mod0 <- model.matrix(~1,data=pheno)
  design.sva <- as.formula(paste0("~",var.trait,"+",paste(c(var.batch.fact , var.batch.num ) , collapse = "+")))
  message("SVA model:\n",design.sva)
  mod1 <- model.matrix(design.sva , data = pheno)
  
  svs = sva(dat = as.matrix(counts),mod = mod1 , mod0 = mod0)$sv
  
  colnames(svs) <- paste0("SV", c(1:ncol(svs)))
  
  pheno <- cbind.data.frame(pheno , svs)
  
  if(ncol(svs) < n.SV){
    n.SV = ncol(svs)
  }
  
  var.batch.all <- c(var.batch.fact , var.batch.num , paste0("SV", c(1:n.SV)))
  
  OutPrefix = paste0(OutPrefix , ".SV",n.SV)
  
}else{
  var.batch.all <- c(var.batch.fact , var.batch.num )
}

message("Running DEG analysis using glmFit function in edgeR...")

design.formula = as.formula(paste0("~0+",var.trait,"+",paste(var.batch.all , collapse = "+")))

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
  write.csv(result , file = paste0(out_name , ".csv"), row.names = F)
}
