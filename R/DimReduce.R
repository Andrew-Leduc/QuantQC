#' Add two numbers.
#'
#' This function takes two numeric inputs and returns their sum.
#'
#' @param x A numeric value.
#' @param y A numeric value.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add_numbers(2, 3)
#' @export
RunPCA <- function(nPOP_obj){
  sc.batch_cor <- nPOP_obj@protein.imputed

  # Correlation matrix for PCA
  cor_mat <- cor(sc.batch_cor,use = 'pairwise.complete.obs')

  # Eigen values of correlation matrix for PCA
  sc.pca <- eigen(cor_mat)
  scx<-as.data.frame(sc.pca$vectors)
  colnames(scx)<-paste0("PC",1:ncol(scx))

  # Calculate and plot % variance of PCs
  percent_var <- sc.pca$values/sum(sc.pca$values)*100
  plot(1:length(percent_var), percent_var, xlab="PC", ylab="% of variance explained")


  # Map meta data for plotting PCAs
  scx$ID <-colnames(sc.batch_cor)
  scx <- scx %>% left_join(batch_label,by = c('ID'))


}


PlotPCA <- function(nPOP_obj, by = condition){

}
