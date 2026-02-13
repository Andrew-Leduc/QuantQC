

#' Evaluate negative controls in a QQC object.
#'
#' Computes quality metrics comparing single cells to negative controls.
#' Dispatches to DDA or DIA-specific evaluation depending on the MS type
#' stored in the QQC object. For DDA data, coefficient of variation (CV)
#' analysis is performed; for DIA data, peptide counts are compared between
#' single cells and negative controls.
#'
#' @param QQC A QQC object with populated \code{meta.data} and \code{raw_data} or \code{matricies} slots.
#' @return The QQC object with the \code{neg_ctrl.info} slot populated with negative control evaluation results.
#' @export
EvaluateNegativeControls <- function(QQC){

  if(nrow(QQC@meta.data)==0){
    return("Must map sample identites to data first")
  }


  if(QQC@ms_type == 'DDA'){
    QQC <- EvaluateNegativeControls_DDA(QQC)
  }

  if(QQC@ms_type == 'DIA' | QQC@ms_type ==  'DIA_C'){
    QQC <- EvaluateNegativeControls_DIA(QQC)
  }


  return(QQC)


}




CVs <- function(QQC){
  cell_id <- QQC@meta.data
  # Normalize peptide data
  mat_norm <- Normalize_reference_vector(QQC@matricies@peptide)

  mat_norm <- as.data.frame(mat_norm)
  mat_norm$Protein <- QQC@matricies@peptide_protein_map$Protein
  mat_norm$pep <- QQC@matricies@peptide_protein_map$seqcharge

  # convert to data.table for fast computation
  mat_norm <- as.data.table(mat_norm)

  # Melt data frame
  mat_norm.melt <- data.table::melt(mat_norm, id = c('Protein','pep'))


  # Call the function on the data.table
  CV_mat <- fast_cv(mat_norm.melt)

  # Revert back to DF
  CV_mat <- as.data.frame(CV_mat)
  CV_mat$value <- NULL

  # count number for protein with multiple peptides in each cell
  CV_mat_count <- CV_mat %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(counter = sum(is.na(cvq)==F))


  # only look at cells with atleast 20 proteins with multiple peptides
  CV_mat_count <- CV_mat_count %>% filter(counter > 20)


  # Take median protein CV per cell
  CV_mat <- CV_mat %>%
    dplyr::group_by( variable) %>%
    dplyr::summarise(cvq = median(cvq,na.rm = T))

  CV_mat <- CV_mat %>% filter(variable %in% CV_mat_count$variable)


  # store values from cells, not negative ctrls
  pos <- cell_id %>% filter(sample != 'neg')

  CV_mat$value <- NA
  CV_mat$value[CV_mat$variable %in% pos$ID] <- 'cell'
  CV_mat$value[!CV_mat$variable %in% pos$ID] <- 'neg'

  CV_mat_pos <- CV_mat %>% filter(value == 'cell')
  CV_mat_neg <- CV_mat %>% filter(value == 'neg')

  return(CV_mat)

}



Count_peptides_per_cell <- function(sc.data,cellenONE_meta,good_cells = NULL){



  # Get IDs for negative controls
  negative_IDs <- cellenONE_meta$ID[cellenONE_meta$sample == 'neg']


  # If negative control has 0 peptides detected, it will be left out of sc.data,
  # find negative controls with 0 peptides detected so 0s can be included on plot

  numb_neg_controls <- length(intersect(negative_IDs,colnames(sc.data)))


  #If negative controls are not all 0 peptides and there are more than 2, count number peptides for non 0 negative controls
  if(numb_neg_controls > 1){

    # Data matrix for negative controls

    neg_mat <- sc.data[,negative_IDs]

    # Count number peptides in each negative control
    sum_neg <- colSums(is.na(neg_mat)==F)
    sum_neg_int <- colSums(neg_mat,na.rm = T)
    # Make data frame for plotting
    neg_df <- as.data.frame(sum_neg)
    colnames(neg_df) <- 'Number_precursors'
    neg_df$intense <- sum_neg_int
    neg_df$type <- 'negative ctrl'



    # # Add in 0s for any negative controls with 0 peptides
    # if(zero_peptide_negs != 0 ){
    #   for(i in 1:zero_peptide_negs){
    #     neg_df[nrow(neg_df) + 1,] = c(0,0,"negative ctrl")
    #   }
    # }


  }else if(numb_neg_controls == 1){
    # If the negative controls all have 0 peptides measured
    neg_vect <- sc.data[,negative_IDs]

    # Make data frame for plotting
    neg_df <- as.data.frame(matrix(data = 0,ncol = 1,nrow = length(negative_IDs)))
    colnames(neg_df) <- 'Number_precursors'
    neg_df$intense <- sum(is.na(neg_vect)==F)
    neg_df$type <- 'negative ctrl'

  }else{
    # Make data frame for plotting
    neg_df <- as.data.frame(matrix(data = 0,ncol = 1,nrow = 1))
    colnames(neg_df) <- 'Number_precursors'
    neg_df$intense <- 0
    neg_df$type <- 'negative ctrl'
  }



  if(is.null(good_cells)==F){
    # Filter for single cells that passed CV test (for pSCoPE data ONLY)
    sum_cell <- sc.data[,colnames(sc.data) %in% good_cells$variable]

  }else{

    # for DIA its all non zero single cells
    good_cells <- cellenONE_meta %>% filter(sample != 'neg')
    sum_cell <- sc.data[,colnames(sc.data) %in% (good_cells$ID)]
  }

  # Sum number peptides for all real cells
  sum_cell_id <- colSums(is.na(sum_cell)==F)
  sum_cell_intense <- colSums(sum_cell,na.rm = T)

  # Make cell count data frame for plotting
  pos_df <- as.data.frame(sum_cell_id)
  colnames(pos_df) <- 'Number_precursors'
  pos_df$intense <-  sum_cell_intense
  pos_df$type <- 'single cells'

  # Combind positive and negative control data frames
  neg_vs_pos_DF <- rbind(pos_df,neg_df)

  return(neg_vs_pos_DF)
}



EvaluateNegativeControls_DDA <- function(QQC){

  # Compute CVs of cells and negative controls, function outputs CV plot and list of cells with CVs
  CVm <- CVs(QQC)


  # Get IDs of cells with median protein CVs (good_cells is a df with cell ID and CV)
  good_cells <- CVm %>% filter(cvq < .9)

  good_cells <- good_cells %>% filter(value != 'neg')

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



#' Plot negative control evaluation results.
#'
#' Generates diagnostic plots comparing single cells to negative controls.
#' For DDA data, produces a combined plot of intensity distributions and a
#' density plot of median protein CVs with an annotated CV threshold. For DIA
#' data, produces a histogram of log10 intensity distributions.
#'
#' @param QQC A QQC object with a populated \code{neg_ctrl.info} slot (from \code{EvaluateNegativeControls}).
#' @param CV_thresh Numeric CV threshold used to separate good cells from failed cells (DDA only).
#' @return A ggplot object displaying the negative control diagnostic plots.
#' @export
PlotNegCtrl <- function(QQC,CV_thresh){


  plot_data <- QQC@neg_ctrl.info

  if(QQC@ms_type == 'DDA'){

    #peps <- ggplot(plot_data, aes(x = Number_precursors, fill = type)) + geom_histogram(position = 'identity', alpha = .5) + ggtitle(paste0('precursors per cell')) + ylab('# of samples')+dot_plot

    peps <- ggplot(plot_data, aes( x = log10(intense), fill = type)) + geom_histogram(position = 'identity',alpha = .5)  + ylab('# of samples')+dot_plot+
      xlab('sum(log10(Intensity))')

    CV_mat_pos <- plot_data %>% filter(value == 'cell')
    CV_mat_neg <- plot_data %>% filter(value == 'neg')

    cvs <- ggplot(data=plot_data, aes(x=cvq,fill=value)) + geom_density( alpha=0.5,adjust=1.5) + dot_plot+
      scale_fill_manual(values=my_col3) +
      xlab("CV, peptides from same protein") + ylab("Fraction of cells") + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
      coord_cartesian(xlim=c(.1,.65))+
      annotate("text", x=0.2, y= 14, label=paste0(sum(CV_mat_pos$cvq < CV_thresh)," cells"), size=10, color=my_col3[c(1)])+
      annotate("text", x=0.64, y= 12, label=paste0(sum(CV_mat_neg$cvq > CV_thresh,na.rm = T)," Ctr -"), size=10, color=my_col3[c(2)])+
      annotate("text", x=0.63, y= 14, label=paste0(sum(CV_mat_pos$cvq > CV_thresh)," cells"), size=10, color=my_col3[c(1)])+
      annotate("text", x=0.2, y= 12, label=paste0((sum(CV_mat_neg$cvq < CV_thresh,na.rm = T)-1)," Ctr -"), size=10, color=my_col3[c(2)])+
      ggtitle('Min 3 proteins w/ many peps')+
      theme(legend.position = "none") + geom_vline(xintercept=CV_thresh, lty=2, size=2, color="gray50") + theme(plot.margin = margin(1, 1, 0, 1, "cm"))


    return(peps+cvs)
  }


  if(QQC@ms_type == 'DIA' | QQC@ms_type ==  'DIA_C'){
    peps <- ggplot(plot_data, aes( x = log10(intense), fill = type)) + geom_histogram(position = 'identity',alpha = .5) + ggtitle(paste0('Negative ctrl Vs Single cells')) + ylab('# of samples')+dot_plot+
      xlab('log10(Intensity)')

    return(peps)
  }

}

#' Filter out low-quality cells from a QQC object.
#'
#' Removes cells that fail quality thresholds based on negative control evaluation.
#' For DIA data, cells are filtered by a minimum log10 intensity threshold.
#' For DDA data, cells are filtered by both a CV threshold and an optional
#' minimum intensity threshold. Negative controls are always removed.
#'
#' @param QQC A QQC object with a populated \code{neg_ctrl.info} slot.
#' @param CV_thresh Numeric CV threshold; cells with CVs above this value are removed (DDA only). Default \code{NA}.
#' @param min_intens Numeric minimum log10 intensity threshold; cells below this are removed. Default \code{NA}.
#' @return The QQC object with the peptide matrix (and peptide mask for DIA) filtered to retain only passing cells.
#' @export
FilterBadCells <- function(QQC, CV_thresh = NA, min_intens = NA){

  neg_filter <- QQC@neg_ctrl.info
  peptide_data <- QQC@matricies@peptide


  if(QQC@ms_type == 'DIA' | QQC@ms_type == 'DIA_C'){
    peptide_mask <- QQC@matricies@peptide_mask
    neg_filter <- neg_filter %>% dplyr::filter(type != 'negative ctrl')
    neg_filter <- neg_filter %>% dplyr::filter(log10(intense) > min_intens)
    peptide_mask <- peptide_mask[,colnames(peptide_mask) %in% neg_filter$variable]
    QQC@matricies@peptide_mask <- peptide_mask

  }

  if(QQC@ms_type == 'DDA'){
    neg_filter <- neg_filter %>% dplyr::filter(value != 'neg')
    if(is.na(min_intens)==F){
      neg_filter <- neg_filter %>% dplyr::filter(log10(intense) > min_intens)
    }
    neg_filter <- neg_filter %>% dplyr::filter(cvq < CV_thresh)


  }

  peptide_data <- peptide_data[,colnames(peptide_data) %in% neg_filter$variable]

  QQC@matricies@peptide <- peptide_data


  return(QQC)

}


#' Trim excess peptides per protein to the top five.
#'
#' For each protein with more than five mapped peptides, retains up to five
#' peptides selected by the intersection of those ranked highest by median
#' intensity and those ranked highest by observation count. This reduces
#' redundancy while keeping the most informative peptides per protein.
#'
#' @param QQC A QQC object whose \code{matricies} slot contains the peptide matrix and peptide-protein map.
#' @return The QQC object with the peptide matrix and peptide-protein map trimmed to the selected peptides.
#' @export
Trim_extra_peptides <- function(QQC){
  peptide_data <- QQC@matricies@peptide
  #protein_dat <- r3_10day_male@matricies@protein
  peptide_protein_map <- QQC@matricies@peptide_protein_map

  upl <- unique(peptide_protein_map$Protein)

  # List describing which protein each peptide comes from
  prot_list <- peptide_protein_map$Protein
  rownames(peptide_protein_map) <- peptide_protein_map$seqcharge

  count <- 0

  # Loop over each protein, calculate correlations between peptides mapping to a protein

  new_pep_list <- c()
  for(i in upl){

    # Matrix for a single protein
    mat_p1 <- peptide_data[which(prot_list == i),]
    if(is.null(nrow(mat_p1)) == FALSE ){

      if(nrow(mat_p1) >5){
        abs_peps1 <- rowMedians(as.matrix(mat_p1[rownames(mat_p1),]),na.rm=T)
        abs_peps2 <- rowSums(is.na(mat_p1[rownames(mat_p1),])==F)
        abs_peps1 <- abs_peps1[order(-abs_peps1)]
        abs_peps2 <- abs_peps2[order(-abs_peps2)]
        sect <- intersect(names(abs_peps1)[1:5],names(abs_peps2)[1:5])

        if(length(sect) < 4){
          sect <- names(abs_peps1)[1:5]
        }

        new_pep_list <- c(new_pep_list,sect)
      }else{
        new_pep_list <- c(new_pep_list,rownames(mat_p1))
      }
    }else{
      new_pep_list <- c(new_pep_list,rownames(peptide_protein_map)[peptide_protein_map$Protein == i])
    }
  }

  peptide_data <- peptide_data[new_pep_list,]

  peptide_protein_map <- peptide_protein_map[new_pep_list,]

  QQC@matricies@peptide <- peptide_data
  QQC@matricies@peptide_protein_map <- peptide_protein_map

  return(QQC)

}


#' Filter a matrix by maximum allowed missingness in rows and columns.
#'
#' Removes columns (cells) and then rows (features) that exceed the specified
#' fraction of missing (NA) values. Column filtering is applied first, followed
#' by row filtering on the reduced matrix.
#'
#' @param mat A numeric matrix to filter.
#' @param pct.r Maximum allowed fraction of NA values per row (0 to 1).
#' @param pct.c Maximum allowed fraction of NA values per column (0 to 1).
#' @return The filtered matrix with rows and columns exceeding the missingness thresholds removed.
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
