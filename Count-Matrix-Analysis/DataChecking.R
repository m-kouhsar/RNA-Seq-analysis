message("loading requiered packages...")
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(DESeq2)))
suppressWarnings(suppressMessages(library(edgeR)))
suppressWarnings(suppressMessages(library(WGCNA)))
suppressWarnings(suppressMessages(library(stringr)))

args = commandArgs(T)

counts.file <- trimws(args[1])
pheno.file <- trimws(args[2])
var.trait <- trimws(args[3])
var.fact <- trimws(args[4])
var.num <- trimws(args[5])
normalize.method <- trimws(args[6])
lib.size.threshold <- as.numeric(trimws(args[7]))
gFilter.min.count <- as.numeric(trimws(args[8]))
gFilter.min.prop <- as.numeric(trimws(args[9]))
OutPrefix <- trimws(args[10])


message("Input arguments:")
message("        Count matrix: ", counts.file)
message("        Phenotype file: ", pheno.file)
message("        Trait variable: ", var.trait)
message("        Numeric variables: ", var.num)
message("        Categorical variables: ", var.fact)
message("        Output files prefix: ", OutPrefix)
cat("\n")
########################################################################
#
#          Reading the data
#
########################################################################
message("Reading the data ...")
counts <- read.table(file = counts.file, stringsAsFactors = F,header = T, row.names = 1,check.names=F)
pheno <- read.csv(pheno.file , row.names = 1 , stringsAsFactors = F)

if(!identical(colnames(counts) , rownames(pheno))){
  stop("Colnames in the count matrix are not equal to the rownames in the phenotype file!")
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
plot_data <- cbind.data.frame(pheno , sample = colnames(dge) , lib.size = dge$samples$lib.size)
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
keep <- edgeR::filterByExpr(counts,group = pheno[,var.trait],min.count = 5, min.prop = 0.75)
message("Low count genes: ",sum(!keep) , " out of ",nrow(counts))
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
message("PCA analysis...")
pcs_obj <- prcomp(t(counts.norm), scale. = T , center = T , rank. = 10)
PCs <- as.data.frame(pcs_obj$x)

plot_data <- cbind.data.frame(plot_data , PCs)

plot_data$PC1_z <- scale(PCs$PC1)
plot_data$PC2_z <- scale(PCs$PC2)
plot_data$outliers <- (abs(plot_data$PC1_z) > 3 | abs(plot_data$PC2_z) > 3)

plot_data$AveExpr <- colMeans(counts.norm , na.rm = T)
PC.PVar <- round(summary(pcs_obj)$importance[2,1:2]*100,digits = 2)

all_var = c(var.trait , var.fact , var.num)
cor_ <- matrix(data = NA,nrow =10,ncol = length(all_var) )
colnames(cor_) <- all_var
rownames(cor_) <- colnames(PCs)[1:10]

cor_pval <- cor_

for (i in 1:10) {
  for (j in 1:length(all_var)) {
    var_ <- all_var[j]
    
    if(var_ %in% var.fact){
      res1<-cor.test(as.numeric(PCs[,i]),as.numeric(as.factor(pheno[,var_])), method="spearman",exact = FALSE)
      cor_[i,var_]<-as.numeric(res1$estimate)
      cor_pval[i,var_]<-as.numeric(res1$p.value)
    }else{
      res1<-cor.test(as.numeric(PCs[,i]),as.numeric(pheno[,var_]), method="pearson",exact = FALSE)
      cor_[i,var_]<-as.numeric(res1$estimate)
      cor_pval[i,var_]<-as.numeric(res1$p.value)
    }
  }
}

textMatrix = paste(signif(cor_, 2), "\n(",signif(cor_pval, 1), ")", sep = "")
dim(textMatrix) = dim(cor_)

pdf(file = paste0(OutPrefix , ".PCA.pdf"),width = 10,height = 10)
par(mar = c(12, 8,3, 3))
labeledHeatmap(Matrix = cor_,
               xLabels = colnames(cor_),
               yLabels = paste0(rownames(cor_),"(",plot_data$PC.PVar[1:10],"%)"),
               ySymbols = rownames(cor_),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("PCA Analysis"))

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
pdf(file = paste0(OutPrefix , ".hClust.pdf"),width = 16,height = 10)
hc = hclust(distance, method = "average")
plot(hc,xlab = "", sub = "")
graphics.off()

#################################################
write.csv(plot_data , file = paste0(OutPrefix , ".DataChecking.csv") , row.names = F)

