
args = commandArgs(T)

script_dir <- trimws(args[1])
counts.file <- trimws(args[2])
pheno.file <- trimws(args[3])
var.trait <- trimws(args[4])
var.fact <- trimws(args[5])
var.num <- trimws(args[6])
normalize.method <- trimws(args[7])
lib.size.threshold <- as.numeric(trimws(args[8]))
gFilter.min.count <- as.numeric(trimws(args[9]))
gFilter.min.prop <- as.numeric(trimws(args[10]))
OutPrefix <- trimws(args[11])

message("Input arguments:")
message("        Count matrix: ", counts.file)
message("        Phenotype file: ", pheno.file)
message("        Trait variable: ", var.trait)
message("        Numeric variables: ", var.num)
message("        Categorical variables: ", var.fact)
message("        Normalization method: ", normalize.method)
message("        Library size threshold: ", lib.size.threshold)
message("        Minimum count threshold for gene filtering: ", gFilter.min.count)
message("        Minimum proportion of the samples for gene filtering: ", gFilter.min.prop)
message("        Output files prefix: ", OutPrefix)
cat("\n")
########################################################################
#
#          Reading the data
#
########################################################################
source(paste0(script_dir , "/DataChecking.Functions.R"))
message("Reading the data ...")
counts <- read.table(file = counts.file, stringsAsFactors = F,header = T, row.names = 1,check.names=F)
pheno <- read.csv(pheno.file , row.names = 1 , stringsAsFactors = F)

if(!identical(colnames(counts) , rownames(pheno))){
  message("Warning message:\nColnames in the count matrix are not equal to the rownames in the phenotype file!\nShared names will be considered.")
  index <- intersect(colnames(counts), rownames(pheno))
  message("Number of shared names: " , length(index))
  if(length(index) < 1){
    stop("There is no shared ID between Phneotype and Count data!")
  }
  counts <- counts[,index]
  pheno <- pheno[index , ]
}

var.fact <- trimws(str_split_1(var.fact , pattern = ","))
var.num <- trimws(str_split_1(var.num , pattern = ","))

pheno <- pheno[,c(var.fact,var.num,var.trait)]

pheno[,var.trait] <- as.factor(pheno[,var.trait])

for (i in 1:length(var.fact)) {
  pheno[,var.fact[i]] <- as.factor(pheno[,var.fact[i]])
}
for (i in 1:length(var.num)) {
  pheno[,var.num[i]] <- as.numeric(pheno[,var.num[i]])
}

########################################################################
#
#          Samples read depth
#
########################################################################
dge <- DGEList(counts)
plot_data <- cbind.data.frame(sample = colnames(dge) ,pheno , lib.size = dge$samples$lib.size)
p <- ggplot(data = plot_data , aes(x=sample , y=lib.size))+
  geom_bar(stat = "identity")+
  ylab("Library size")+
  ggtitle("Library Sizes")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),plot.title = element_text(hjust = 0.5))+
  geom_hline(yintercept = lib.size.threshold, colour ="red")

message(sum(dge$samples$lib.size < lib.size.threshold) ,"/" , ncol(dge), " samples have library size < " , lib.size.threshold)
plot_data$lib.size.pass <- (dge$samples$lib.size >= lib.size.threshold)

pdf(file = paste0(OutPrefix , ".LibSize.pdf"),width = 20,height = 10)
print(p)
graphics.off()
########################################################################
#
#          Filtering low count genes
#
########################################################################
message("filtering low count genes...")
message("Genes that do not have a minimum count of ",gFilter.min.count, " in at least ",(gFilter.min.prop*100) , "% of samples in the smallest group will be removed.")
keep <- edgeR::filterByExpr(counts,group = pheno[,var.trait],min.count = gFilter.min.count, min.prop = gFilter.min.prop)
message(sum(!keep),"/",nrow(counts)," genes removed. Remaining genes:", sum(keep))
counts <- counts[keep,]

########################################################################
#
#          Checking for possible outliers and batch effects
#
########################################################################

if(normalize.method == "cpm"){
  message("Normalizing gene counts uisng cpm function in edgeR...")
  counts.norm <- edgeR::cpm(counts , log = TRUE)
  
}else{
  if(normalize.method == "vst"){
    message("Normalizing gene counts uisng vst method in DESeq2...")
    counts.norm <- DESeq2::varianceStabilizingTransformation(as.matrix(counts) , blind = T, fitType = "parametric")
  }
}

pdf(file = paste0(OutPrefix , ".Count.BoxPlot.pdf"),width = 16,height = 10)
boxplot(as.matrix(counts),xlab = "Samples" , ylab= "Counts")
title("Raw count")
boxplot(as.matrix(counts.norm),xlab = "Samples" , ylab= paste0("Normalized counts (" , ifelse(normalize.method == "cpm" , "CPM" , "VST"),")"))
title(paste0("Normalized counts (" , ifelse(normalize.method == "cpm" , "Count Per Million" , "variance stabilizing transformation"),")"))
graphics.off()

######################################################################
cor_plot1 <- generate_corr_analysis(data = subset(plot_data , select =c(-lib.size.pass,-sample)),
                                   categorical_vars =c(var.trait , var.fact) ,
                                   plot.title = "Correlation test on all covariates" )
message("PCA analysis...")
pcs_obj <- prcomp(t(counts.norm), scale. = T , center = T , rank. = 10)
PCs <- as.data.frame(pcs_obj$x)

plot_data <- cbind.data.frame(plot_data , PCs)

plot_data$PC1_z <- scale(PCs$PC1)
plot_data$PC2_z <- scale(PCs$PC2)
plot_data$Outliers.PC.ZScore <- (abs(plot_data$PC1_z) > 3 | abs(plot_data$PC2_z) > 3)

plot_data$AveExpr <- colMeans(counts.norm , na.rm = T)
PC.PVar <- round(summary(pcs_obj)$importance[2,1:10]*100,digits = 2)

all_var = c(var.trait , var.fact , var.num)

cor_plot2 <- generate_corr_analysis(data = subset(plot_data , select =c(-lib.size.pass,-sample,-Outliers.PC.ZScore,-PC1_z , -PC2_z)),
                                   categorical_vars =c(var.trait , var.fact) ,
                                   plot.title = "Correlation test on all covariates" )

cor_ <- cor_plot2$cor.mat[paste0("PC",c(1:10)),all_var]
cor_pval <- cor_plot2$p.mat[paste0("PC",c(1:10)),all_var]
impact_val <- matrix(data = NA,nrow =10,ncol = length(all_var) )
colnames(impact_val) <- all_var
rownames(impact_val) <- colnames(PCs)[1:10]

for (i in 1:10) {
  for (j in 1:length(all_var)) {
    var_ <- all_var[j]
    
    if(cor_pval[i,var_] < 0.05){
      impact_val[i , var_] <- cor_[i,var_]**2 * PC.PVar[i]
    }else{
      impact_val[i , var_] = 0
    }
  }
}

textMatrix_cor = paste(signif(cor_, 2), "\n(",signif(cor_pval, 1), ")", sep = "")
dim(textMatrix_cor) = dim(cor_)

pdf(file = paste0(OutPrefix , ".PCA.pdf"),width = 10,height = 10)

print(cor_plot1$plot)

par(mar = c(12, 8,3, 3))
labeledHeatmap(Matrix = cor_,
               xLabels = colnames(cor_),
               yLabels = paste0(rownames(cor_),"(",PC.PVar[1:10],"%)"),
               ySymbols = rownames(cor_),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix_cor,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = "PCA Analysis (Correlation and P-value)")

textMatrix_imp = signif(impact_val, 2)
dim(textMatrix_imp) = dim(cor_)

par(mar = c(12, 8,3, 3))
labeledHeatmap(Matrix = cor_,
               xLabels = colnames(cor_),
               yLabels = paste0(rownames(cor_),"(",PC.PVar[1:10],"%)"),
               ySymbols = rownames(cor_),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix_imp,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               legendLabel = "Cor",
               main = expression("Significant PC Correlations" ~ (Impact ~ value == r^2 %*% "%Var"[PC])))


PCA.mahal <- mahalanobis.outlier(Data = counts.norm , method = "pca" , plot.title = "Outliers by Mahalanobis Distance")
plot_data$mdist <- PCA.mahal$Data.2D$mdist
plot_data$pchisq <- PCA.mahal$Data.2D$pchisq
plot_data$qchisq <- PCA.mahal$Data.2D$qchisq
plot_data$Outliers.Mahalanobis <- PCA.mahal$Data.2D$Outlier

print(PCA.mahal$Plot.2D)
print(PCA.mahal$Plot.QQ)

for (var_ in c(var.trait , var.fact , var.num)) {
  p <- ggplot(data = plot_data, aes_string(x = "PC1_z" , y = "PC2_z" , colour = var_))  + geom_point() + 
    # Add the +/- 3 z-score lines
    geom_vline(xintercept = 3, linetype = "dashed", color = "red") +
    geom_vline(xintercept = -3, linetype = "dashed", color = "red") +
    geom_hline(yintercept = 3, linetype = "dashed", color = "red") +
    geom_hline(yintercept = -3, linetype = "dashed", color = "red") +
    ggtitle("PCA Plot (Z-Scores) with Outlier Thresholds")+
    xlab("PC1 (Z Score)")+
    ylab("PC2 (Z Score)")+
    theme_bw()
  print(p)
  if(var_ %in% c(var.fact,var.trait)){
    p <- ggplot() + geom_boxplot(data = plot_data , aes_string(x = var_ ,  y = "AveExpr" , fill = var_))+
      ylab("Average gene expression")
    print(p)
  }
}
graphics.off()
######################################################################
message("Hirarchical clustering...")

distance <- dist(t(counts.norm) , method = "euclidean")
hc = hclust(distance, method = "average")

pdf(file = paste0(OutPrefix , ".hClust.pdf"),width = 18,height = 10)
plot(hc,xlab = "", sub = "",cex=0.5)
graphics.off()

#################################################
write.csv(plot_data , file = paste0(OutPrefix , ".DataChecking.csv") , row.names = F)

