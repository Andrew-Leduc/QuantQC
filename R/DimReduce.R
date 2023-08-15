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
ComputePCA <- function(nPOP_obj){
  sc.batch_cor <- nPOP_obj@protein.imputed
  cellenONE_meta <- nPOP_obj@meta.data

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

  nPOP_obj@reductions <- list(PCA = scx)

  return(nPOP_obj)

}


PlotPCA <- function(nPOP_obj, by = 'Condition'){
  PCA_plot <- nPOP_obj@reductions[['PCA']]

  if(by == 'Condition'){
    pca_plot <- ggplot(PCA_plot, aes(x = PC1, y = PC2, color = sample)) + geom_point()+
      dot_plot
  }
  if(by == 'Total protein'){
    pca_plot <- ggplot(PCA_plot, aes(x = PC1, y = PC2, color = prot_total)) + geom_point()+
      dot_plot + scale_color_gradient2(midpoint = median(PCA_plot$prot_total), low = 'red',mid = 'white', high = 'blue')
  }

  if(by == 'Label'){
    pca_plot <- ggplot(PCA_plot, aes(x = PC1, y = PC2, color = label)) + geom_point()+
      dot_plot
  }

  if(by == 'Run order'){
    pca_plot <- ggplot(PCA_plot, aes(x = PC1, y = PC2, color = Order)) + geom_point()+
      dot_plot + scale_color_gradient2(midpoint = median(PCA_plot$Order,na.rm = T), low = 'red',mid = 'white', high = 'blue')
  }

  pca_plot

}


ComputeUMAP <- function(nPOP_obj){
  protein_Data <- nPOP_obj@protein.imputed
  scx <- nPOP_obj@reductions[['PCA']]

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

  nPOP_obj@reductions[['UMAP']] <- um_plot

  return(nPOP_obj)
}


PlotUMAP <- function(nPOP_obj, by = 'Cluster'){
  UMAP_plot <- nPOP_obj@reductions[['UMAP']]

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


FeatureUMAP <- function(nPOP_obj, prot = NA, imputed = T){
  UMAP_plot <- nPOP_obj@reductions[['UMAP']]

  if(imputed == T){
    prot_mat <- nPOP_obj@protein.imputed
  }
  if(imputed == F){
    prot_mat <- nPOP_obj@protein
  }

  UMAP_plot$protein <- prot_mat[prot,]

  umap_plot <- ggplot(UMAP_plot, aes(x = UMAP_1, y = UMAP_2, color = protein)) + geom_point()+
      dot_plot + scale_color_gradient2(midpoint = 0, low = 'blue',mid = 'white', high = 'red')



  umap_plot

}

