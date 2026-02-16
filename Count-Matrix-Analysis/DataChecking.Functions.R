message("loading requiered packages...")
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(DESeq2)))
suppressWarnings(suppressMessages(library(edgeR)))
suppressWarnings(suppressMessages(library(WGCNA)))
suppressWarnings(suppressMessages(library(stringr)))
suppressWarnings(suppressMessages(library(ggcorrplot)))
suppressWarnings(suppressMessages(library(reshape2)))

########################################################################
mahalanobis.outlier <- function(Data , method = "pca", plot.title=NA , tsne.seed = NA, pca.scale=T , pca.center=T){
  
  suppressMessages(library(car))
  suppressMessages(library(dplyr))
  suppressMessages(library(ggplot2))
  
  if(!is.na(tsne.seed)){
    set.seed(seed = tsne.seed) 
  }
  method = match.arg(arg = method , choices = c("pca" , "tsne") , several.ok = F)
  
  if(method == "tsne"){
    suppressMessages(library(Rtsne))
    tsne_out <- Rtsne(t(Data),dims = 2,)
    Data.2D <- data.frame(D1 = tsne_out$Y[,1], 
                          D2 = tsne_out$Y[,2])
    rownames(Data.2D) <- colnames(Data)
  }
  if(method == "pca"){
    pc <- prcomp(t(Data),scale. = pca.scale,center = pca.center , rank. =2)
    pc.importance <- summary(pc)$importance[2,]
    Data.2D <- as.data.frame(pc$x)
  }
  Data.2D[,1] <- scale(Data.2D[,1] , center = T , scale = T)
  Data.2D[,2] <- scale(Data.2D[,2] , center = T , scale = T)
  center_ <- colMeans(Data.2D)
  cov_ <- cov(Data.2D)
  
  # Calculating the squared Mahalanobis distance
  Data.2D$mdist <- mahalanobis(
    x = Data.2D,
    center = center_,
    cov = cov_
  )
  
  cutoff <- qchisq(p = 0.95, df = 2)
  R <- sqrt(cutoff)
  
  ellipse_ <- car::ellipse(
    center = center_[1:2],
    shape = cov_[1:2,1:2],
    radius = R,
    segments = 150,
    draw = FALSE
  )
  ellipse_ <- as.data.frame(ellipse_)
  colnames(ellipse_) <- colnames(Data.2D)[1:2]
  
  Data.2D$pchisq <- pchisq(Data.2D$mdist, df = 2, lower.tail = FALSE)
  
  Data.2D <- Data.2D %>%
    mutate(Outlier = ifelse(mdist > cutoff, 'Yes', 'No'))
  if(method == "pca"){
    p1 <- ggplot(Data.2D, aes(x = PC1 , y = PC2, color = Outlier))+
      xlab(paste0("PC1 ",round(pc.importance[1],digits = 2)*100,"%"))+
      ylab(paste0("PC2 ",round(pc.importance[2],digits = 2)*100,"%"))
    plot.subtitle = "Dimentinality reduction method: PCA"
  }
  if(method == "tsne"){
    p1 <- ggplot(Data.2D, aes(x = D1 , y = D2, color = Outlier))
    plot.subtitle = "Dimentinality reduction method: tSNE"
  }
  if(is.na(plot.title)){
    plot.title = ""
  }
  p1 <- p1 +
    geom_point(size = 3) +
    geom_point(aes(center_[1], center_[2]) , size = 5 , color = 'blue') +
    geom_polygon(data = ellipse_, fill = 'white', color = 'black', alpha = 0.3) +
    scale_color_manual(values = c('gray44', 'red')) +
    labs(title =plot.title,subtitle = paste0("Outliers in 2D Plot, ",plot.subtitle)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0))
  
  p2 <- ggplot() + geom_point(data = Data.2D , aes(x=mdist , y= pchisq , color = Outlier)) + theme_bw() +
    scale_color_manual(values = c('black', 'red'))+
    geom_vline(xintercept = cutoff , color="red")+
    ylab("Chi-Square probability")+
    xlab("Square Mahalanobis distance")+
    labs(title =plot.title,subtitle = paste0("Outliers in Chi-Square Plot, ",plot.subtitle))
  
  Data.2D$qchisq =qchisq(ppoints(length(Data.2D$mdist)), df = 2)
  p3 <- ggplot() + geom_point(data = Data.2D, aes(x=sort(qchisq), y=sort(mdist))) + theme_bw() +
    geom_abline(aes(slope = 1, intercept = 0),color="red")+
    xlab("Chi-Square quantiles")+
    ylab("Square Mahalanobis distance quantiles")+
    labs(title =plot.title,subtitle = paste0("QQ Plot, ",plot.subtitle))
  
  return(list(Data.2D = Data.2D , Plot.2D=p1 , Plot.Dist=p2 , Plot.QQ = p3))
}

########################################################################

generate_corr_analysis <- function(data, 
                                   categorical_vars = NULL, 
                                   plot.title = "Correlation Heatmap") {
  
  # 1. Data Preparation
  # Create a working copy of the data
  df_work <- data
  
  # Convert specified categorical variables to numeric ranks (1, 2, 3...)
  # This is required because correlation functions only accept numeric inputs.
  if (!is.null(categorical_vars)) {
    for (cat_col in categorical_vars) {
      if (cat_col %in% names(df_work)) {
        # Convert to factor first to ensure levels, then to integer
        df_work[[cat_col]] <- as.numeric(as.factor(df_work[[cat_col]]))
      }
    }
  }
  
  # Select only numeric columns (now including the converted categorical ones)
  df_num <- df_work[, sapply(df_work, is.numeric)]
  
  # 2. Initialize Matrices
  n <- ncol(df_num)
  col_names <- colnames(df_num)
  
  M <- matrix(NA, n, n)      # Correlation Matrix
  p_mat <- matrix(NA, n, n)  # P-value Matrix
  
  rownames(M) <- colnames(M) <- col_names
  rownames(p_mat) <- colnames(p_mat) <- col_names
  
  # 3. Calculate Mixed Correlations (Double Loop)
  for (i in 1:n) {
    for (j in 1:n) {
      
      # Diagonal is always 1, p-value 0
      if (i == j) {
        M[i, j] <- 1
        p_mat[i, j] <- 0
        next
      }
      
      # Determine Method:
      # If EITHER variable is in the categorical list -> Spearman
      var1_name <- col_names[i]
      var2_name <- col_names[j]
      
      is_cat <- (var1_name %in% categorical_vars) || (var2_name %in% categorical_vars)
      current_method <- if (is_cat) "spearman" else "pearson"
      
      # Run Test
      # suppressWarnings to handle ties in Spearman rank
      test_res <- suppressWarnings(
        cor.test(df_num[, i], df_num[, j], method = current_method)
      )
      
      M[i, j] <- test_res$estimate
      p_mat[i, j] <- test_res$p.value
    }
  }
  
  # 4. Prepare Data for Plotting (Long Format)
  # Convert matrix to long format for ggplot
  # We use expand.grid logic to ensure order matches the matrix
  plot_data <- expand.grid(Var1 = rownames(M), Var2 = colnames(M))
  plot_data$Corr <- as.vector(M)
  plot_data$P_val <- as.vector(p_mat)
  
  # Create Label: "0.85 \n (p=1.2e-03)"
  plot_data$Label <- paste0(
    format(round(plot_data$Corr, 2), nsmall = 2), 
    "\n(p=", 
    formatC(plot_data$P_val, format = "e", digits = 2), 
    ")"
  )
  
  # 5. Generate Plot
  gg_plot <- ggcorrplot(M, 
                        method = "square", 
                        type = "full", 
                        lab = FALSE, 
                        title = plot.title, 
                        outline.color = "black",
                        show.diag = F,
                        ggtheme = ggplot2::theme_minimal()) +
    geom_text(data = plot_data, 
              aes(x = Var1, y = Var2, label = Label), 
              size = 3, 
              inherit.aes = FALSE) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  
  # 6. Create Output Dataframe (Non-Redundant)
  ut_index <- upper.tri(M, diag = FALSE)
  
  stats_df <- data.frame(
    Variable_1 = rownames(M)[row(M)[ut_index]],
    Variable_2 = colnames(M)[col(M)[ut_index]],
    Correlation = M[ut_index],
    P_value = p_mat[ut_index]
  )
  
  # Add a column indicating which method was likely used
  stats_df$Method <- apply(stats_df, 1, function(row) {
    if (row['Variable_1'] %in% categorical_vars || row['Variable_2'] %in% categorical_vars) {
      return("Spearman")
    } else {
      return("Pearson")
    }
  })
  
  stats_df$Correlation_Rounded <- round(stats_df$Correlation, 2)
  stats_df$P_value_Sci <- formatC(stats_df$P_value, format = "e", digits = 2)
  
  return(list(plot = gg_plot, data = stats_df , cor.mat = M , p.mat = p_mat))
}

###################################################################################

