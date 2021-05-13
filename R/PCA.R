#' Performs PCA
#'
#' Performs PCA using the base function prcomp with scaling. Outputs ggplot graphs
#' and a summary table
#'
#'
#' @param data csv file with OTU identifiers in the first column, followed by morphological characters, each in in a separate column.
#' @return Returns a summary table called PCA_summary.csv that includes eigenvalues and PCA loadings in addition
#'  to PCA plots of the first three PC's
#' @import ggplot2 gridExtra
#' @export

GS_pca <- function(data){
  pca <- prcomp(log10(data[,2:ncol(data)]), scale=TRUE)
  scores <- data.frame(Species=data[,1], pca$x[,1:3])
  p1 <- ggplot(scores, aes(PC1, PC2, group=data[,1])) + theme_bw() +
    geom_point(aes(color=data[,1]), size=5) +
    theme(legend.key = element_blank()) +
    stat_ellipse(aes(x=PC1, y=PC2,color=data[,1], group=data[,1]), level=0.67) +
    theme(legend.position="bottom", legend.title = element_blank(), legend.text=element_text(size=12), axis.text=element_text(size=10), axis.title=element_text(size=12))

  p2 <- ggplot(scores, aes(PC2, PC3, group=data[,1])) + theme_bw() +
    geom_point(aes(color=data[,1]), size=5) +
    theme(legend.key = element_blank()) +
    stat_ellipse(aes(x=PC2, y=PC3,color=data[,1], group=data[,1]), level=0.67) +
    theme(legend.position="bottom", legend.title = element_blank(), legend.text=element_text(size=12), axis.text=element_text(size=10), axis.title=element_text(size=12))

  grid.arrange(p1, p2, ncol=2)

  loadings <- pca$rotation
  temp <- summary(pca)
  sum <- temp$importance #Extract important variables
  eigen <- pca$sdev^2 #Extract eigenvalues
  sum <- rbind(sum, eigen, loadings) #Combine eigenvalues, loadings, and summary stats into single table
  write.csv(sum, "PCA_summary.csv")
  print(scores)
}
