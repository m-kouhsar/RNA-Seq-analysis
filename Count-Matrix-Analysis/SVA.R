########################################################################################################
#                                                                                                      #
# Running SVA analysis to find Sorrugate Variables and deal with hidden batch effect with RNA-Seq data #
#
#
#
#



########################################################################################################
suppressMessages(library(sva))
suppressMessages(library(edgeR))
suppressMessages(library(stringr))

args <- commandArgs(T)

counts.file <- trimws(args[1])
pheno.file <- trimws(args[2])
gFilter.min.count <- as.numeric(trimws(args[3]))
gFilter.min.prop <- as.numeric(trimws(args[4]))
model.main <- as.formula(trimws(args[5]))
model.null <- as.formula(trimws(args[6]))
var.interest <- trimws(args[7])
outliers <- trimws(args[8])
OutPrefix <- trimws(args[9])

cat("##########################################################################\n")
message("Input arguments:")
message("        Count matrix file: ",counts.file)
message("        Phenotype file: ",pheno.file)
message("        Main model (mod in sva): ",model.main)
message("        Null model (mod0 in sva): ",model.null)
message("        Variable of Interest: ",var.interest)
message("        Minimum count threshold for gene filtering: ", gFilter.min.count)
message("        Minimum proportion of the samples for gene filtering: ", gFilter.min.prop)
message("        Outlie samples: ",outliers)
message("        Output files prefix: ",OutPrefix)
cat("##########################################################################\n")

dir.create(dirname(OutPrefix) , recursive = T)

message("Reading the input data...")

counts <- read.table(counts.file , header = T , row.names = 1 , sep = "\t", stringsAsFactors = F, check.names = F)
pheno <- read.csv(pheno.file , row.names = 1 , stringsAsFactors = F)

vars <- unique(var.interest,c(all.vars(model.main) , all.vars(model.null)))
if(!all(vars %in% colnames(pheno))){
  stop("The following variables can't find in the phenotype data:\n",paste(setdiff(vars , colnames(pheno)) , collapse = ";"))
}

if(!identical(colnames(counts) , rownames(pheno))){
  warning("Row names in Phenotype data are not matched with column names in count data. Shared IDs will be considered.")
  index <- intersect(rownames(pheno) , colnames(counts))
  message("Number of shared Sample IDs: ",length(index))
  if(length(index)<2){
    stop("Phenotype data and counts data cannot be matched!")
  }else{
    pheno <- pheno[index , ]
    counts <- counts[, index]
  }
}

#######################################################################################################
outliers <- trimws(str_split_1(outliers , pattern = ","))
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

########################################################################################################
message("filtering low count genes...")
message("Genes that don't have minimum count of ",gFilter.min.count, " in at least ",(gFilter.min.prop*100) , "% of the samples will be removed.")
keep <- edgeR::filterByExpr(round(counts),group = pheno[,var.interest],min.count = gFilter.min.count, min.prop = gFilter.min.prop)
message(sum(!keep),"/",nrow(counts)," genes removed. Remaining genes:", sum(keep))
counts <- counts[keep,]

########################################################################################################
message("Running SVA...")
model.main <- model.matrix(model.main , data = pheno)
model.null <- model.matrix(model.null , data = pheno)
counts.norm <- edgeR::cpm(counts , log = T)
SVs <- sva(dat = as.matrix(counts.norm),mod = model.main,mod0 = model.null)$sv
colnames(SVs) <- paste0("SV",c(1:ncol(SVs)))
pheno.new <- cbind.data.frame(pheno , SVs)

write.csv(pheno.new ,file = paste0(OutPrefix , ".SVS.csv") )

message("Detected SVs are saved in \n",paste0(OutPrefix , ".SVS.csv"))

























