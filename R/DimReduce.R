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
ComputePCA <- function(QQC){
  sc.batch_cor <- QQC@matricies@protein.imputed
  cellenONE_meta <- QQC@meta.data

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
  scx <- scx %>% left_join(cellenONE_meta,by = c('ID'))

  QQC@reductions <- list(PCA = scx)

  return(QQC)

}


PlotPCA <- function(QQC, by = 'Condition'){
  PCA_plot <- QQC@reductions[['PCA']]

  if(by == 'Condition'){
    pca_plot <- ggplot(PCA_plot, aes(x = PC1, y = PC2, color = sample)) + geom_point()+
      dot_plot + ggtitle('PCA by Condition')
  }
  if(by == 'Total protein'){
    pca_plot <- ggplot(PCA_plot, aes(x = PC1, y = PC2, color = prot_total)) + geom_point()+ ggtitle('PCA by Total Cell Intensity')+
      dot_plot + scale_color_gradient2(midpoint = median(PCA_plot$prot_total), low = 'red',mid = 'white', high = 'blue')
  }

  if(by == 'Label'){
    pca_plot <- ggplot(PCA_plot, aes(x = PC1, y = PC2, color = label)) + geom_point()+
      dot_plot + ggtitle('PCA by Label')
  }

  if(by == 'Run order'){
    pca_plot <- ggplot(PCA_plot, aes(x = PC1, y = PC2, color = Order)) + geom_point() + ggtitle('PCA by Run Order')+
      dot_plot + scale_color_gradient2(midpoint = median(PCA_plot$Order,na.rm = T), low = 'red',mid = 'white', high = 'blue')
  }

  pca_plot

}


ComputeUMAP <- function(QQC){
  protein_Data <- QQC@matricies@protein.imputed
  scx <- QQC@reductions[['PCA']]

  prot_umap <- CreateSeuratObject(counts = protein_Data, project = "prot_mat")
  prot_umap <- NormalizeData(prot_umap, normalization.method = "LogNormalize", scale.factor = 10000)
  prot_umap@assays$RNA@scale.data <- protein_Data

  all.genes <- rownames(protein_Data)
  prot_umap <- ScaleData(prot_umap, features = all.genes)
  prot_umap@assays$RNA@scale.data <- protein_Data
  prot_umap <- Seurat::RunPCA(prot_umap, features = all.genes)
  prot_umap <- FindNeighbors(prot_umap, dims = 1:6)
  prot_umap <- FindClusters(prot_umap, resolution = 0.5)
  prot_umap <- RunUMAP(prot_umap, dims = 1:6)

  um_plot <- as.data.frame(prot_umap@reductions[["umap"]]@cell.embeddings)

  um_plot$sample <- scx$sample
  um_plot$cluster <- prot_umap@meta.data[["RNA_snn_res.0.5"]]
  um_plot$lab <- scx$label
  um_plot$prot_total <- scx$prot_total




  # ggplot(um_plot, aes(x = UMAP_1,y = UMAP_2, color = cluster)) + geom_point(size = 3) + theme_classic()+
  #   ggtitle('UMAP colored by Sample')

  QQC@reductions[['UMAP']] <- um_plot

  return(QQC)
}


PlotUMAP <- function(QQC, by = 'Cluster'){
  UMAP_plot <- QQC@reductions[['UMAP']]

  if(by == 'Cluster'){
    umap_plot <- ggplot(UMAP_plot, aes(x = UMAP_1, y = UMAP_2, color = cluster)) + geom_point()+
      dot_plot
  }

  if(by == 'Condition'){
    umap_plot <- ggplot(UMAP_plot, aes(x = UMAP_1, y = UMAP_2, color = sample)) + geom_point()+
      dot_plot
  }
  if(by == 'Total protein'){
    umap_plot <- ggplot(UMAP_plot, aes(x = UMAP_1, y = UMAP_2, color = prot_total)) + geom_point()+
      dot_plot + scale_color_gradient2(midpoint = median(UMAP_plot$prot_total), low = 'red',mid = 'white', high = 'blue')
  }

  if(by == 'Label'){
    umap_plot <- ggplot(UMAP_plot, aes(x = UMAP_1, y = UMAP_2, color = label)) + geom_point()+
      dot_plot
  }

  if(by == 'Run order'){

  }

  umap_plot

}


FeatureUMAP <- function(QQC, prot = NA, imputed = T){
  UMAP_plot <- QQC@reductions[['UMAP']]

  if(imputed == T){
    prot_mat <- QQC@matricies@protein.imputed
  }
  if(imputed == F){
    prot_mat <- QQC@matricies@protein
  }

  UMAP_plot$protein <- prot_mat[prot,]

  umap_plot <- ggplot(UMAP_plot, aes(x = UMAP_1, y = UMAP_2, color = protein)) + geom_point()+
      dot_plot + scale_color_gradient2(midpoint = 0, low = 'blue',mid = 'white', high = 'red')



  umap_plot

}

