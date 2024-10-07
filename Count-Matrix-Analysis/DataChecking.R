
library(ggplot2)
library(DESeq2)
library(WGCNA)
library(stringr)

setwd("./")
counts.file <- "Raw/Count.txt"
pheno.file <- "Raw/Pheno.csv"
var.trait <- "Trait"
var.fact <- "Sex,Plate"
var.num <- "Age,RIN"
OutPrefix <- "Results/Project01"

########################################################################
#
#          Reading the data
#
########################################################################

counts <- read.table(file = counts.file, stringsAsFactors = F,header = T, row.names = 1)
pheno <- read.csv(pheno.file , row.names = 1 , stringsAsFactors = F)
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


identical(colnames(counts) , rownames(pheno))

########################################################################
#
#          Filtering low count genes
#
########################################################################

keep <- edgeR::filterByExpr(counts,group = pheno$Trait,min.count = 5)
table(keep)
counts <- counts[keep,]

########################################################################
#
#          Checking for possible outliers and batch effects
#
########################################################################
#counts.logCPM <- edgeR::cpm(counts , log = T)
counts.vst <- DESeq2::varianceStabilizingTransformation(as.matrix(counts) , blind = T, fitType = "parametric")

pcs <- prcomp(t(counts.vst), scale. = T , center = T)
plot.data <- cbind.data.frame(pheno , pcs$x[,c(1,2)],MeanExpr = colMeans(counts.vst , na.rm = T))
plot.data$PC.PVar <- round(summary(pcs)$importance[2,]*100,digits = 2)

all_var = c(var.trait , var.fact , var.num)
cor_ <- matrix(data = NA,nrow =10,ncol = length(all_var) )
colnames(cor_) <- all_var
rownames(cor_) <- colnames(pcs$x)[1:10]

cor_pval <- cor_

for (i in 1:10) {
  for (j in 1:length(all_var)) {
    var_ <- all_var[j]
    
    if(var_ %in% var.fact){
      res1<-cor.test(as.numeric(pcs$x[,i]),as.numeric(as.factor(pheno[,var_])), method="spearman",exact = FALSE)
      cor_[i,var_]<-as.numeric(res1$estimate)
      cor_pval[i,var_]<-as.numeric(res1$p.value)
    }else{
      res1<-cor.test(as.numeric(pcs$x[,i]),as.numeric(pheno[,var_]), method="pearson",exact = FALSE)
      cor_[i,var_]<-as.numeric(res1$estimate)
      cor_pval[i,var_]<-as.numeric(res1$p.value)
    }
  }
}

textMatrix = paste(signif(cor_, 2), "\n(",signif(cor_pval, 1), ")", sep = "")
dim(textMatrix) = dim(cor_)

pdf(file = paste0(OutPrefix , ".Count.BoxPlot.pdf"),width = 16,height = 10)
boxplot(as.matrix(counts),xlab = "Samples" , ylab= "Counts")
title("Raw count")
boxplot(as.matrix(counts.vst),xlab = "Samples" , ylab= "Normalized counts (VST)")
title("Normalized count (variance stabilizing transformation)")
graphics.off()

pdf(file = paste0(OutPrefix , ".hClust.pdf"),width = 16,height = 10)
distance <- dist(t(counts.vst) , method = "euclidean")
hc = hclust(distance, method = "average")
plot(hc,xlab = "", sub = "")
graphics.off()

pdf(file = paste0(OutPrefix , ".PCA.pdf"),width = 10,height = 10)
par(mar = c(12, 8,3, 3))
labeledHeatmap(Matrix = cor_,
               xLabels = colnames(cor_),
               yLabels = paste0(rownames(cor_),"(",plot.data$PC.PVar[1:10],"%)"),
               ySymbols = rownames(cor_),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("PCA Analysis"))

for (var_ in c(var.trait , var.fact , var.num)) {
  p <- ggplot(data = plot.data, aes_string(x = "PC1" , y = "PC2" , colour = var_))  + geom_point() + theme_bw()
  print(p)
  if(var_ %in% c(var.fact,var.trait)){
    p <- ggplot() + geom_boxplot(data = plot.data , aes_string(x = var_ ,  y = "MeanExpr" , fill = var_))+
      ylab("Average gene expression")
    print(p)
  }
}

graphics.off()

#################################################
#clusters = cbind.data.frame(sample = rownames(pheno),cluster = cutree(hc , k = 3))


