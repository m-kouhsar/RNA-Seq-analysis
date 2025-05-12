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
