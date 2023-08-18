

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
EvaluateNegativeControls <- function(nPOP_obj,CV_thresh){

  if(nrow(nPOP_obj@meta.data)==0){
    return("Must map sample identites to data first")
  }


  if(nPOP_obj@ms_type == 'DDA'){
    nPOP_obj <- EvaluateNegativeControls_DDA(nPOP_obj,CV_thresh)
  }

  if(nPOP_obj@ms_type == 'DIA' | nPOP_obj@ms_type ==  'DIA_C'){
    nPOP_obj <- EvaluateNegativeControls_DIA(nPOP_obj)
  }


  return(nPOP_obj)


}

EvaluateNegativeControls_DDA <- function(nPOP_obj,CV_thresh){



  # Compute CVs of cells and negative controls, function outputs CV plot and list of cells with CVs
  CVm <- CVs(nPOP_obj,CV_thresh)


  # Get IDs of cells with median protein CVs (good_cells is a df with cell ID and CV)
  good_cells <- CVm[[2]] %>% filter(cvq < CV_thresh)


  if(length(good_cells) < 3){

    return('less than three good cells, try increasing CV filter')

  }


  # Count the number of peptides in negative controls and good single cells
  Peptide_counts_by_sample <- Count_peptides_per_cell(nPOP_obj@matricies@peptide,nPOP_obj@meta.data,good_cells)


  # Plot distributions
  Numb_data_points <- ggplot(Peptide_counts_by_sample, aes(x = Number_precursors, fill = type)) + geom_histogram(position = 'identity', alpha = .5) + ggtitle(paste0('# precursors per sample')) + ylab('# of samples')+dot_plot


  #neg_ctrl_data <- left_join(CVm[[2]],Peptide_counts_by_sample)

  # Filter for only good cells
  #cols_to_keep <- c('Protein.Group','seqcharge',as.character(good_cells$variable))

  #Ref_norm_data_filtered <- Ref_norm_data[,colnames(Ref_norm_data) %in% cols_to_keep]
  Peptide_counts_by_sample$variable <- rownames(Peptide_counts_by_sample)
  neg_meta <-  CVm[[2]] %>% left_join(Peptide_counts_by_sample, by = c('variable'))

  nPOP_obj@neg_ctrl.info <- neg_meta

  return(nPOP_obj)


}

EvaluateNegativeControls_DIA <- function(nPOP_obj){


  # Count the number of peptides in negative controls and good single cells
  Peptide_counts_by_sample <- Count_peptides_per_cell(nPOP_obj@matricies@peptide,nPOP_obj@meta.data)

  #Ref_norm_data_filtered <- Ref_norm_data[,colnames(Ref_norm_data) %in% cols_to_keep]
  Peptide_counts_by_sample$variable <- rownames(Peptide_counts_by_sample)

  nPOP_obj@neg_ctrl.info <- Peptide_counts_by_sample

  return(nPOP_obj)


}

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
PlotNegCtrl <- function(nPOP_obj,CV_thresh){


  plot_data <- nPOP_obj@neg_ctrl.info

  if(nPOP_obj@ms_type == 'DDA'){

    peps <- ggplot(plot_data, aes(x = Number_precursors, fill = type)) + geom_histogram(position = 'identity', alpha = .5) + ggtitle(paste0('# precursors per sample')) + ylab('# of samples')+dot_plot

    CV_mat_pos <- plot_data %>% filter(value == 'cell')
    CV_mat_neg <- plot_data %>% filter(value == 'neg')

    cvs <- ggplot(data=plot_data, aes(x=cvq,fill=value)) + geom_density( alpha=0.5,adjust=1.5) + theme_pubr() +
      scale_fill_manual(values=my_col3) +
      xlab("CV of peptides mapping to a protein") + ylab("Fraction of cells") + rremove("y.ticks") + rremove("y.text") +
      font("xylab", size=17) +
      font("x.text", size=15) +
      font('title',size=12)+
      coord_cartesian(xlim=c(.1,.65))+
      annotate("text", x=0.2, y= 14, label=paste0(sum(CV_mat_pos$cvq < CV_thresh)," cells"), size=10, color=my_col3[c(1)])+
      annotate("text", x=0.64, y= 12, label=paste0(sum(CV_mat_neg$cvq > CV_thresh,na.rm = T)," Ctr -"), size=10, color=my_col3[c(2)])+
      annotate("text", x=0.63, y= 14, label=paste0(sum(CV_mat_pos$cvq > CV_thresh)," cells"), size=10, color=my_col3[c(1)])+
      annotate("text", x=0.2, y= 12, label=paste0((sum(CV_mat_neg$cvq < CV_thresh,na.rm = T)-1)," Ctr -"), size=10, color=my_col3[c(2)])+
      ggtitle('Cells need atleast 3 proteins with multiple peptides')+
      rremove("legend") + geom_vline(xintercept=CV_thresh, lty=2, size=2, color="gray50") + theme(plot.margin = margin(1, 1, 0, 1, "cm"))


    return(peps+cvs)
  }


  if(nPOP_obj@ms_type == 'DIA' | nPOP_obj@ms_type ==  'DIA_C'){
    peps <- ggplot(plot_data, aes(x = Number_precursors, fill = type)) + geom_histogram(position = 'identity', alpha = .5) + ggtitle(paste0('# precursors per sample')) + ylab('# of samples')+dot_plot

    return(peps)
  }

}

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
FilterBadCells <- function(nPOP_obj, CV_thresh = NA, min_pep = NA){

  neg_filter <- nPOP_obj@neg_ctrl.info
  peptide_data <- nPOP_obj@matricies@peptide


  if(nPOP_obj@ms_type == 'DIA' | nPOP_obj@ms_type == 'DIA_C'){
    neg_filter <- neg_filter %>% dplyr::filter(type != 'negative ctrl')
    neg_filter <- neg_filter %>% dplyr::filter(Number_precursors > min_pep)

  }

  if(nPOP_obj@ms_type == 'DDA'){
    neg_filter <- neg_filter %>% dplyr::filter(value != 'neg')
    if(is.na(min_pep)==F){
      neg_filter <- neg_filter %>% dplyr::filter(Number_precursors > min_pep)
    }
    neg_filter <- neg_filter %>% dplyr::filter(cvq < CV_thresh)


  }

  peptide_data <- peptide_data[,colnames(peptide_data) %in% neg_filter$variable]
  nPOP_obj@matricies@peptide <- peptide_data
  return(nPOP_obj)

}



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
filt.mat.cr<-function(mat, pct.r,pct.c){

  kc<-c()
  for(k in 1:ncol(mat)){

    pct.na<-length(which(is.na(mat[,k]))) / length(mat[,k])
    if(pct.na <= pct.c){ kc<-c(kc,k)}
    #print(pct.na)


  }

  mat<-mat[,kc]

  kr<-c()
  for(k in 1:nrow(mat)){

    pct.na<-length(which(is.na(mat[k,]))) / length(mat[k,])
    if(pct.na <= pct.r){ kr<-c(kr,k)}
    #print(pct.na)


  }

  mat<-mat[kr,]



  return(mat)

}
