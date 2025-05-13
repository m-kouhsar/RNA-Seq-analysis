CovariatePlot <- function(Data, Phenotype , Factor_var , Numeric_var , PCs = 10, Plot_titel=""){
  
  ##############################################################################################################
  # Data: normalized epression/methylation matrix or data frame with genes/CpG in rows and samples in column
  # Phenotype: A data frame contains information about samples
  # Factor_var: Column names of Categorical information in Phenotype 
  # Numeric_var: Column names of numerical information in Phenotype 
  # PCs: Number of principal components in the analysis
  ##############################################################################################################
  # This function calculate principal components in Data and then create a correlation plot to show the correlation 
  #      (pearson for numerical and spearman for categorical information) between the PCs and the phenotype information about the samples
  #############################################################################################################
  
  suppressMessages(library(ggplot2))
  suppressMessages(library(reshape2))
  
  
  if(!all(c(Factor_var , Numeric_var) %in% colnames(Phenotype))){
    stop("The following variables can't be find in Phenotype data. Check Phenotype columns names.\n",
         paste(setdiff(c(Factor_var , Numeric_var) , colnames(Phenotype)), collapse = ", "))
  }
  
  if(!identical(colnames(Data) , rownames(Phenotype))){
    stop("Column names and row names in Data and Phenotype must be equal!")
  }
  
  if(PCs < 1){stop("Number of PC (PCs) must be >= 1!")}
  
  pca <- prcomp(t(Data) , rank. =PCs)
  pca.pc <- pca$x
  pca.importance <- summary(pca)$importance[2,1:PCs]
  
  corr_mat <- matrix(data = NA , nrow = PCs , ncol = length(c(Factor_var , Numeric_var)))
  rownames(corr_mat) <- paste0("PC",1:PCs,"(",round(pca.importance*100 , digits = 2),"%)")
  colnames(pca.pc) <- paste0("PC",1:PCs,"(",round(pca.importance*100 , digits = 2),"%)")
  colnames(corr_mat) <- c(Factor_var , Numeric_var)
  
  corr_pval <- corr_mat
  
  for (i in rownames(corr_mat)){
    for (j in Factor_var) {
      corr_ <- cor.test(as.numeric(as.factor(Phenotype[,j])) , pca.pc[,i], method="spearman" , exact=F)
      corr_mat[i,j] <- corr_$estimate
      corr_pval[i,j] <- corr_$p.value
    }
  }
  
  for (i in rownames(corr_mat)){
    for (j in Numeric_var) {
      corr_ <- cor.test(as.numeric(Phenotype[,j]) , pca.pc[,i], method="pearson")
      corr_mat[i,j] <- corr_$estimate
      corr_pval[i,j] <- corr_$p.value
    }
    
  }
  
  PlotData <- melt(corr_mat)
  PlotData$value2 <- melt(corr_pval)[,3]
  
  p <- ggplot(data = PlotData, aes(x=Var2, y=Var1, fill=value)) + 
        geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Correlation")+
    theme_minimal()+
    geom_text(aes(x=Var2, y=Var1, label = paste0(round(value, digits = 2),"\n(",formatC(value2, format = "e",digits = 1),")")), 
              color = "black", size = 3)+
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
    xlab("")+ylab("")+
    ggtitle(Plot_titel)
  
  return(p)
  
}