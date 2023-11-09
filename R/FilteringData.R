

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
EvaluateNegativeControls <- function(QQC,CV_thresh){

  if(nrow(QQC@meta.data)==0){
    return("Must map sample identites to data first")
  }


  if(QQC@ms_type == 'DDA'){
    QQC <- EvaluateNegativeControls_DDA(QQC,CV_thresh)
  }

  if(QQC@ms_type == 'DIA' | QQC@ms_type ==  'DIA_C'){
    QQC <- EvaluateNegativeControls_DIA(QQC)
  }


  return(QQC)


}

EvaluateNegativeControls_DDA <- function(QQC,CV_thresh){


  # Compute CVs of cells and negative controls, function outputs CV plot and list of cells with CVs
  CVm <- CVs(QQC,CV_thresh)


  # Get IDs of cells with median protein CVs (good_cells is a df with cell ID and CV)
  good_cells <- CVm %>% filter(cvq < CV_thresh)


  if(length(good_cells) < 3){

    return('less than three good cells, try increasing CV filter')

  }


  # Count the number of peptides in negative controls and good single cells
  Peptide_counts_by_sample <- Count_peptides_per_cell(QQC@matricies@peptide,QQC@meta.data,good_cells)


  # Plot distributions
  Numb_data_points <- ggplot(Peptide_counts_by_sample, aes(x = Number_precursors, fill = type)) + geom_histogram(position = 'identity', alpha = .5) + ggtitle(paste0('# precursors per sample')) + ylab('# of samples')+dot_plot



  # Filter for only good cells
  #cols_to_keep <- c('Protein.Group','seqcharge',as.character(good_cells$variable))

  Peptide_counts_by_sample$variable <- rownames(Peptide_counts_by_sample)
  neg_meta <-  CVm %>% left_join(Peptide_counts_by_sample, by = c('variable'))

  QQC@neg_ctrl.info <- neg_meta

  return(QQC)


}

EvaluateNegativeControls_DIA <- function(QQC){


  # Count the number of peptides in negative controls and good single cells
  Peptide_counts_by_sample <- Count_peptides_per_cell(QQC@matricies@peptide,QQC@meta.data)

  #Ref_norm_data_filtered <- Ref_norm_data[,colnames(Ref_norm_data) %in% cols_to_keep]
  Peptide_counts_by_sample$variable <- rownames(Peptide_counts_by_sample)

  QQC@neg_ctrl.info <- Peptide_counts_by_sample

  return(QQC)


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
PlotNegCtrl <- function(QQC,CV_thresh){


  plot_data <- QQC@neg_ctrl.info

  if(QQC@ms_type == 'DDA'){

    peps <- ggplot(plot_data, aes(x = Number_precursors, fill = type)) + geom_histogram(position = 'identity', alpha = .5) + ggtitle(paste0('precursors per cell')) + ylab('# of samples')+dot_plot

    CV_mat_pos <- plot_data %>% filter(value == 'cell')
    CV_mat_neg <- plot_data %>% filter(value == 'neg')

    cvs <- ggplot(data=plot_data, aes(x=cvq,fill=value)) + geom_density( alpha=0.5,adjust=1.5) + dot_plot+
      scale_fill_manual(values=my_col3) +
      xlab("CV, peptides from same protein") + ylab("Fraction of cells") + rremove("y.ticks") + rremove("y.text") +
      coord_cartesian(xlim=c(.1,.65))+
      annotate("text", x=0.2, y= 14, label=paste0(sum(CV_mat_pos$cvq < CV_thresh)," cells"), size=10, color=my_col3[c(1)])+
      annotate("text", x=0.64, y= 12, label=paste0(sum(CV_mat_neg$cvq > CV_thresh,na.rm = T)," Ctr -"), size=10, color=my_col3[c(2)])+
      annotate("text", x=0.63, y= 14, label=paste0(sum(CV_mat_pos$cvq > CV_thresh)," cells"), size=10, color=my_col3[c(1)])+
      annotate("text", x=0.2, y= 12, label=paste0((sum(CV_mat_neg$cvq < CV_thresh,na.rm = T)-1)," Ctr -"), size=10, color=my_col3[c(2)])+
      ggtitle('Min 3 proteins w/ many peps')+
      rremove("legend") + geom_vline(xintercept=CV_thresh, lty=2, size=2, color="gray50") + theme(plot.margin = margin(1, 1, 0, 1, "cm"))


    return(peps+cvs)
  }


  if(QQC@ms_type == 'DIA' | QQC@ms_type ==  'DIA_C'){
    peps <- ggplot(plot_data, aes( x = log10(intense), fill = type)) + geom_histogram(position = 'identity',alpha = .5) + ggtitle(paste0('Negative ctrl Vs Single cells')) + ylab('# of samples')+dot_plot+
      xlab('log10(Intensity)')

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
FilterBadCells <- function(QQC, CV_thresh = NA, min_intens = NA){

  neg_filter <- QQC@neg_ctrl.info
  peptide_data <- QQC@matricies@peptide
  peptide_mask <- QQC@matricies@peptide_mask

  if(QQC@ms_type == 'DIA' | QQC@ms_type == 'DIA_C'){
    neg_filter <- neg_filter %>% dplyr::filter(type != 'negative ctrl')
    neg_filter <- neg_filter %>% dplyr::filter(log10(intense) > min_intens)

  }

  if(QQC@ms_type == 'DDA'){
    neg_filter <- neg_filter %>% dplyr::filter(value != 'neg')
    if(is.na(min_pep)==F){
      neg_filter <- neg_filter %>% dplyr::filter(log10(intense) > min_intens)
    }
    neg_filter <- neg_filter %>% dplyr::filter(cvq < CV_thresh)


  }

  peptide_data <- peptide_data[,colnames(peptide_data) %in% neg_filter$variable]
  peptide_mask <- peptide_mask[,colnames(peptide_mask) %in% neg_filter$variable]
  QQC@matricies@peptide <- peptide_data
  QQC@matricies@peptide_mask <- peptide_mask

  return(QQC)

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
