### CV and cell filtering


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

CVs <- function(nPOP_obj,thresh){
  cell_id <- nPOP_obj@meta.data
  # Normalize peptide data
  mat_norm <- Normalize_reference_vector(nPOP_obj@matricies@peptide)

  mat_norm <- as.data.frame(mat_norm)
  mat_norm$Protein <- nPOP_obj@matricies@peptide_protein_map$Protein
  mat_norm$pep <- nPOP_obj@matricies@peptide_protein_map$seqcharge

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

  # Plot distribution
  CV_plot <- ggplot(data=CV_mat, aes(x=cvq,fill=value)) + geom_density( alpha=0.5,adjust=1.5) + theme_pubr() +
    scale_fill_manual(values=my_col3) +
    xlab("CV of peptides mapping to a protein") + ylab("Fraction of cells") + rremove("y.ticks") + rremove("y.text") +
    font("xylab", size=17) +
    font("x.text", size=15) +
    font('title',size=12)+
    coord_cartesian(xlim=c(.1,.65))+
    annotate("text", x=0.2, y= 14, label=paste0(sum(CV_mat_pos$cvq < thresh)," cells"), size=10, color=my_col3[c(1)])+
    annotate("text", x=0.64, y= 12, label=paste0(sum(CV_mat_neg$cvq > thresh,na.rm = T)," Ctr -"), size=10, color=my_col3[c(2)])+
    annotate("text", x=0.63, y= 14, label=paste0(sum(CV_mat_pos$cvq > thresh)," cells"), size=10, color=my_col3[c(1)])+
    annotate("text", x=0.2, y= 12, label=paste0((sum(CV_mat_neg$cvq < thresh,na.rm = T)-1)," Ctr -"), size=10, color=my_col3[c(2)])+
    ggtitle('Cells need atleast 3 proteins with multiple peptides')+
    rremove("legend") + geom_vline(xintercept=thresh, lty=2, size=2, color="gray50") + theme(plot.margin = margin(1, 1, 0, 1, "cm"))



  return(list(CV_plot,CV_mat))

}



### Peptide correlations


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
SharedPeptideCor <- function(nPOP_obj, res = 'sc'){
  peptide_data <- nPOP_obj@matricies@peptide
  protein_dat <- nPOP_obj@matricies@protein
  peptide_protein_map <- nPOP_obj@matricies@peptide_protein_map


  peptide_data <- Normalize_reference_vector(peptide_data,log = T)



  # Initialized empty matricies to store correlations between peptides
  mat_stor = matrix(data = NA,nrow = 10000,ncol = 3)
  mat_stor_fake = matrix(data = NA,nrow = 10000,ncol = 2)

  # List of all the unique proteins in data
  upl <- unique(peptide_protein_map$Protein)

  # List describing which protein each peptide comes from
  prot_list <- peptide_protein_map$Protein


  count <- 0

  # Loop over each protein, calculate correlations between peptides mapping to a protein
  for(i in upl){

    # Matrix for a single protein
    mat_p1 <- peptide_data[which(prot_list == i),]
    if(is.null(nrow(mat_p1)) == FALSE ){

      # calculate pairwise observations (how many times peptides are observed in same cells)
      obs_mat <- psych::pairwiseCount(t(mat_p1))

      # calculate correlations between peptides
      cor_mat <- cor(t(mat_p1),use = 'pairwise.complete.obs')

      # Require atleast 5 pairwise observations to plot the correlations
      cor_mat[obs_mat < 4] <- NA

      # Turn correlation matrix to a vector for storage, store correlations and protein identity
      cor_mat <- cor_mat[lower.tri(cor_mat)]
      if(is.na(median(cor_mat,na.rm = T)) == F){
        count = count + 1
        mat_stor[count,1] <- i
        mat_stor[count,2] <- median(cor_mat,na.rm = T)
        mat_stor[count,3] <- nrow(mat_p1)
      }
    }
  }


  # Shuffle around the protein list so its correlating peptides from different proteins to eachother
  # This should form our "NULL" distribution
  prot_list_mixed <- sample(peptide_protein_map$Protein)

  # Same loop as above but on missmatched peptides
  for(i in upl){
    mat_p1 <- peptide_data[which(prot_list_mixed == i),]
    if(is.null(nrow(mat_p1)) == FALSE ){
      obs_mat <- psych::pairwiseCount(t(mat_p1))
      cor_mat <- cor(t(mat_p1),use = 'pairwise.complete.obs')
      cor_mat[obs_mat < 4] <- NA
      cor_mat <- cor_mat[lower.tri(cor_mat)]
      if(is.na(median(cor_mat,na.rm = T)) == F){
        count = count + 1
        mat_stor_fake[count,1] <- i
        mat_stor_fake[count,2] <- median(cor_mat,na.rm = T)
      }
    }
  }

  # Format storage matrix
  mat_stor <- as.data.frame(mat_stor)
  colnames(mat_stor) <- c('Protein','Cor','Obs')
  mat_stor <- mat_stor %>% filter(is.na(Cor) == F)
  mat_stor$Cor <- as.numeric(mat_stor$Cor)

  # Store a protein matrix with of protein quant for proteins with multiple peptides
  mat_stor <- mat_stor %>% filter(Protein %in% rownames(protein_dat))
  rownames(mat_stor) <- mat_stor$Protein
  protein_mat_mult_pep <- protein_dat[mat_stor$Protein,]
  pep_cor <- mat_stor[rownames(protein_mat_mult_pep),]

  # Binning protein abundance by fold change for faceting correlations
  # Peptides from proteins that change more should correlate better (more signal = better signal/noise)

  pep_cor$FC <- rowMeans(abs(protein_mat_mult_pep),na.rm = T)
  pep_cor$FC[pep_cor$FC <= .4] <- .4
  pep_cor$FC[pep_cor$FC > .4 & pep_cor$FC <= .8] <- .8
  pep_cor$FC[pep_cor$FC > .8 & pep_cor$FC <= 1.2] <- 1.2
  pep_cor$FC[pep_cor$FC > 1.2 & pep_cor$FC <= 2] <- 2
  pep_cor$FC[pep_cor$FC > 2] <- 3
  pep_cor$FC <- as.character(pep_cor$FC)


  nPOP_obj@pep.cor <- list(pep_cor,median(as.numeric(mat_stor_fake[,2]),na.rm = T))

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
PlotPepCor <- function(nPOP_obj){

  pep_cor <- nPOP_obj@pep.cor[[1]]
  null_dist <- nPOP_obj@pep.cor[[2]]

  ggplot(pep_cor, aes(y = Cor, x = FC)) + geom_boxplot(color="black", fill = 'gray') + xlab('Mean abs(protein fold change)') +
    ylab('Correlation between peptides mapping to a protein')+
    stat_summary(fun.data=f, geom="text", vjust=-0.5, col="blue")+
    geom_hline(yintercept = null_dist, col = "red")+theme_linedraw()




}






#####
# Counting Peptide and protein numbers
####

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
inSet_completness <- function(Raw_data, cellenOne_meta){
  count = 0
  for(i in unique(paste0(cellenONE_meta$injectWell,cellenONE_meta$plate))){
    i
    Data_inset <- as.data.frame(Raw_cell_mat)[,which(grepl(i,colnames(Raw_cell_mat)))]

    Data_inset <- Data_inset[which(is.na(rowMeans(Data_inset,na.rm = T))==F),]

    One_set <- Data_completeness(Data_inset)

    if(count == 0){
      Cell_complete <- One_set[[1]]
      #Protein_complete <- One_set[[2]]
    }else{
      Cell_complete <- rbind(Cell_complete,One_set[[1]])
      #Protein_complete <- rbind(Cell_complete,One_set[[2]])
    }

    count = count+1
  }

  return(Cell_complete)
}



### General Data Analysis functions



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
Count_peptides_per_cell <- function(sc.data,cellenONE_meta,good_cells = NULL){
  #sc.data <- Raw_cell_mat
  # Get IDs for negative controls
  negative_IDs <- cellenONE_meta$ID[cellenONE_meta$sample == 'neg']


  # If negative control has 0 peptides detected, it will be left out of sc.data,
  # find negative controls with 0 peptides detected so 0s can be included on plot
  zero_peptide_negs <- length(negative_IDs) - length(intersect(negative_IDs,colnames(sc.data)))


  #If negative controls are not all 0 peptides, count number peptides for non 0 negative controls
  if(zero_peptide_negs != length(negative_IDs)){

    # Data matrix for negative controls
    neg_mat <- sc.data[,negative_IDs]

    # Count number peptides in each negative control
    sum_neg <- colSums(is.na(neg_mat)==F)

    # Make data frame for plotting
    neg_df <- as.data.frame(sum_neg)
    colnames(neg_df) <- 'Number_precursors'
    neg_df$type <- 'negative ctrl'

    # Add in 0s for any negative controls with 0 peptides
    if(zero_peptide_negs != 0 ){
      for(i in 1:zero_peptide_negs){
        neg_df[nrow(neg_df) + 1,] = c(0,"negative ctrl")
      }
    }


  }else{
    # If the negative controls all have 0 peptides measured

    # Make data frame for plotting
    neg_df <- as.data.frame(matrix(data = 0,ncol = 1,nrow = length(negative_IDs)))
    colnames(neg_df) <- 'Number_precursors'
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
  sum_cell <- colSums(is.na(sum_cell)==F)


  # Make cell count data frame for plotting
  pos_df <- as.data.frame(sum_cell)
  colnames(pos_df) <- 'Number_precursors'
  pos_df$type <- 'single cells'

  # Combind positive and negative control data frames
  neg_vs_pos_DF <- rbind(pos_df,neg_df)

  return(neg_vs_pos_DF)
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
PlotProtAndPep <- function(nPOP_obj){
  nPOP_obj <- Trachea_3_7_prot
  # count peptide numbers
  numb_pep <- colSums(is.na(nPOP_obj@matricies@peptide)==F)
  numb_pep <- as.data.frame(numb_pep)
  colnames(numb_pep) <- 'Number_precursors'

  # count protein numbers
  numb_prot <- colSums(is.na(nPOP_obj@matricies@protein)==F)
  numb_prot <- as.data.frame(numb_prot)
  colnames(numb_prot) <- 'Number_proteins'

  # Plot peptide and protein numbers
  if(nPOP_obj@ms_type == 'DDA'){
    pep_number <- ggplot(numb_pep, aes(x = Number_precursors)) + geom_histogram(bins = 30,position = 'identity',alpha = .5) + ggtitle('# precursors per sample') + rremove('legend')+ylab('# of single cells')+dot_plot

    prot_number<- ggplot(numb_prot, aes(x = Number_proteins)) + geom_histogram(bins = 30,position = 'identity',alpha = .5) + ggtitle('# proteins per sample') + ylab('# of single cells')+dot_plot

    return(pep_number+prot_number)
  }

  if(nPOP_obj@ms_type == 'DIA' | nPOP_obj@ms_type == 'DIA_C'){

    numb_pep_NF <- colSums(nPOP_obj@matricies@peptide_mask==T)
    numb_pep_NF <- as.data.frame(numb_pep_NF)
    colnames(numb_pep_NF) <- 'Number_precursors'
    numb_pep_NF$Filter <- 'No Ch Qvalue Filter'
    numb_pep$Filter <- paste0(nPOP_obj@misc[['ChQ']], ' Channel Q Value')

    numb_pep <- rbind(numb_pep,numb_pep_NF)

    # count protein numbers
    numb_prot_NF <- colSums(nPOP_obj@matricies@protein_mask==T)
    numb_prot_NF <- as.data.frame(numb_prot_NF)
    colnames(numb_prot_NF) <- 'Number_proteins'
    numb_prot_NF$Filter <- 'No Ch Qvalue Filter'
    numb_prot$Filter <- paste0(nPOP_obj@misc[['ChQ']], ' Channel Q Value')

    numb_prot <- rbind(numb_prot,numb_prot_NF)

    pep_number <- ggplot(numb_pep, aes(x = Number_precursors, fill = Filter)) + geom_histogram(bins = 30,position = 'identity',alpha = .5) + ggtitle('# precursors per sample') + rremove('legend')+ylab('# of single cells')+dot_plot

    prot_number<- ggplot(numb_prot, aes(x = Number_proteins, fill = Filter)) + geom_histogram(bins = 30,position = 'identity',alpha = .5) + ggtitle('# proteins per sample') + ylab('# of single cells')+dot_plot


    return(pep_number+prot_number)
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
PlotDataComplete <- function(nPOP_obj){
  nPOP_obj <- Trachea_3_7_prot
  data <- nPOP_obj@matricies@protein

  missingness_cell_filt <- colSums(is.na(data) == F)/nrow(data)
  missingness_cell_mat <- reshape2::melt(missingness_cell_filt)

  missingness_prot_filt <- rowSums(is.na(data)==F)/ncol(data)
  missingness_prot_mat <- reshape2::melt(missingness_prot_filt)

  if(nPOP_obj@ms_type == 'DDA'){
    mp <- ggplot(missingness_prot_mat, aes(x = value)) +
      geom_histogram(bins = 20,position = 'identity',alpha = .5) + ggtitle(paste0('Protein completness, ', nrow(data) ,' proteins'))+rremove('legend') +ylab('# of proteins')+xlab('fraction values present')+dot_plot

    mc <- ggplot(missingness_cell_mat, aes(x = value)) +
      geom_histogram(bins = 20,position = 'identity',alpha = .5) + ggtitle('Cell completness') + ylab('# of single cells')+
      xlab('fraction values present')+dot_plot


    return(mp+mc)
  }


  if(nPOP_obj@ms_type == 'DIA' | nPOP_obj@ms_type == 'DIA_C'){

    missingness_cell_NF <- colSums(nPOP_obj@matricies@protein_mask)/nrow(data)
    missingness_cell_mat_NF <- reshape2::melt(missingness_cell_NF)
    colnames(missingness_cell_mat_NF) <- 'Cell_miss'
    missingness_cell_mat_NF$Filter <- 'No Ch Qvalue Filter'
    colnames(missingness_cell_mat) <- 'Cell_miss'
    missingness_cell_mat$Filter <- paste0(nPOP_obj@misc[['ChQ']], ' Channel Q Value')


    missingness_cell_mat <- rbind(missingness_cell_mat,missingness_cell_mat_NF)

    # count protein numbers
    missingness_prot_NF <- rowSums(nPOP_obj@matricies@protein_mask)/ncol(data)
    missingness_prot_mat_NF <- reshape2::melt(missingness_prot_NF)
    colnames(missingness_prot_mat_NF) <- 'Miss_proteins'
    missingness_prot_mat_NF$Filter <- 'No Ch Qvalue Filter'
    colnames(missingness_prot_mat) <- 'Miss_proteins'
    missingness_prot_mat$Filter <- paste0(nPOP_obj@misc[['ChQ']], ' Channel Q Value')

    missingness_prot_mat <- rbind(missingness_prot_mat,missingness_prot_mat_NF)

    mp <- ggplot(missingness_cell_mat, aes(x = Cell_miss, fill = Filter)) +
      geom_histogram(bins = 20,position = 'identity',alpha = .5) + ggtitle('Cell completness') + ylab('# of single cells')+
      xlab('fraction values present')+dot_plot
    mc <-  ggplot(missingness_prot_mat, aes(x = Miss_proteins,fill = Filter)) +
      geom_histogram(bins = 20,position = 'identity',alpha = .5) + ggtitle(paste0('Protein completness, ', nrow(data) ,' proteins'))+rremove('legend') +ylab('# of proteins')+xlab('fraction values present')+dot_plot



    return(mp+mc)
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
PlotSCtoCarrierRatio <- function(nPOP_obj){

  do_plot <- "No carrier used"

  if(nPOP_obj@ms_type == "DDA"){

    sc.data <- nPOP_obj@raw_data
    good_cells <- colnames(nPOP_obj@matricies@peptide)

    cellenOne_meta <- nPOP_obj@meta.data
    good_cells_p_negs <- c(good_cells,cellenOne_meta$ID[cellenOne_meta$sample == 'neg'])

    ri.index<-which(colnames(sc.data)%in%paste0("Reporter.intensity.",1:18))

    sc.data[sc.data==0] <- NA
    sc.data[, ri.index] <- sc.data[, ri.index] / sc.data[, ri.index[1]]


    sc.data <- sc.data[,c('Raw.file','Well','plate',paste0("Reporter.intensity.",5:18))]
    sc.data <- reshape2::melt(sc.data, id = c('Raw.file','Well','plate'))

    sc.data$ID <- paste0(sc.data$Well,sc.data$plate,sc.data$variable)
    sc.data <- sc.data %>% filter(ID %in% good_cells_p_negs)
    sc.data <- sc.data %>% group_by(ID) %>% dplyr::summarise(value = median(value,na.rm = T))

    sc.data <- sc.data %>% dplyr::left_join(cellenOne_meta, by = c('ID'))

    sc.data$value <- 1/sc.data$value

    do_plot <- ggplot(sc.data, aes(x = sample, y = value)) + geom_boxplot() + ylab('# cells / carrier amount') +dot_plot

  }

  if(nPOP_obj@ms_type == "DIA_C"){

    cellenOne_meta <- nPOP_obj@meta.data
    Filt_ratio <- nPOP_obj@matricies@peptide #%>% dplyr::select(good_cells$variable)
    Filt_ratio <- reshape2::melt(Filt_ratio)
    Filt_ratio_cell <- Filt_ratio %>% group_by(Var2) %>% dplyr::summarise(value_new = median(value,na.rm = T))
    Filt_ratio_cell <- Filt_ratio_cell %>% left_join(cellenOne_meta,by = c('Var2'='ID'))


    Filt_rat_plot <- ggplot(Filt_ratio_cell, aes(y = 1/(value_new), x =  sample)) + geom_boxplot() +
      ylab('# single cells in carrier') + ggtitle( paste0('Carrier is median ', round(median(1/Filt_ratio$value,na.rm = T),2),' single cells'))

    do_plot <- Filt_rat_plot


  }

  do_plot

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
PlotCellSizeVsIntensity <- function(nPOP_obj, type = 'sample'){
  meta <- nPOP_obj@meta.data

  good_cells <- colnames(nPOP_obj@matricies@peptide)
  meta <- meta %>% dplyr::filter(ID %in% good_cells)

  title_text <- paste0('Correlation = ', round(cor((meta$diameter/2)^3,meta$prot_total,use = 'pairwise.complete.obs'),2))

  if(type == 'sample'){
    plot_ <-  ggplot(meta, aes(x = (diameter/2)^3,y = prot_total,color = sample)) + geom_point() +
      ggtitle(title_text)

    return(plot_)
  }

  if(type == 'Run order'){
    plot_ <-    ggplot(meta, aes(x = (diameter/2)^3,y = prot_total,color = Order)) +
      geom_point() + ggtitle(title_text)+
      scale_color_gradient2(midpoint = median(meta$Order,na.rm = T), low = 'blue',mid = 'white', high = 'red')+
      dot_plot

    return(plot_)
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
ProteinClustConsistency <- function(nPOP_obj, prot = NA, type = 'line'){

  if(is.null(nPOP_obj@reductions[['UMAP']])==F){

    #nPOP_obj <- Trachea_3_7_prot
    #prot = "P33267"

    if(nPOP_obj@ms_type =='DDA'){
      raw_pep_intense <- nPOP_obj@raw_data
      raw_pep_intense <- raw_pep_intense %>% filter(Leading.razor.protein == prot)
      raw_pep_intense <- raw_pep_intense %>% group_by(seqcharge) %>% summarise(pep_raw = median(Reporter.intensity.1,na.rm=T))
      raw_pep_intense <- raw_pep_intense[order(-raw_pep_intense$pep_raw),]

    }
    if(nPOP_obj@ms_type =='DIA' | nPOP_obj@ms_type =='DIA_C'){
      raw_pep_intense <- nPOP_obj@raw_data
      raw_pep_intense <- raw_pep_intense %>% filter(Protein.Group == prot)
      raw_pep_intense <- raw_pep_intense %>% group_by(seqcharge) %>% summarise(pep_raw = median(Ms1.Area,na.rm=T))
      raw_pep_intense <- raw_pep_intense[order(-raw_pep_intense$pep_raw),]

    }





    # Get 4 best peptides


    clusters <- nPOP_obj@reductions[['UMAP']]
    prot_map <- nPOP_obj@matricies@peptide_protein_map %>% filter(Protein == prot)

    raw_pep_intense_prot <- raw_pep_intense %>% filter(seqcharge %in% prot_map$seqcharge)
    if(nrow(raw_pep_intense_prot)>4){
      raw_pep_intense_prot <- raw_pep_intense_prot[1:4,]
    }
    if(nrow(raw_pep_intense_prot)==1){
      return('Only 1 peptide')
    }

    prot_map <- prot_map %>% filter(seqcharge %in% raw_pep_intense_prot$seqcharge)
    pep_mat <- (Normalize_reference_vector(nPOP_obj@matricies@peptide,log = T))
    rownames(pep_mat) <- nPOP_obj@matricies@peptide_protein_map$seqcharge

    pep_mat <- pep_mat[prot_map$seqcharge,]

    pep_mat <- reshape2::melt(pep_mat)
    clusters$ID <- rownames(clusters)
    pep_mat <- pep_mat %>%  dplyr::left_join(clusters, by = c('Var2' = 'ID'))

    pep_mat_values <-  pep_mat %>% dplyr::group_by(Var1,cluster) %>%
      dplyr::summarise(med_abs = median(value,na.rm = T))

    pep_mat_SD<- pep_mat %>%  dplyr::group_by(Var1,cluster) %>%
      dplyr::summarise(sd_error = sd(value,na.rm = T)/sum(is.na(value)==F))

    pep_mat_count <- pep_mat %>%  dplyr::group_by(Var1,cluster) %>%
      dplyr::summarise(numb_data = sum(is.na(value)==F)/(sum(is.na(value)==F)+sum(is.na(value)==T)))

    prot_values <- pep_mat %>%  dplyr::group_by(cluster) %>%
      dplyr::summarise(med_abs = median(value,na.rm = T))



    for(i in 1:nrow(raw_pep_intense_prot)){
      prot_values_pep1 <- prot_values
      prot_values_pep1$Var1 <- unique(pep_mat_values$Var1)[i]
      if(i == 1){
        prot_values <- prot_values_pep1
      }else{
        prot_values <- rbind(prot_values,prot_values_pep1)
      }

    }

    prot_values$sd_error <- 0
    prot_values$numb_dp <- 0
    prot_values$numb_dp_z <- 0
    prot_values$FractionCellsExpress <- " Less 33%"

    pep_mat_values$sd_error <- pep_mat_SD$sd_error
    pep_mat_values$numb_dp <- pep_mat_count$numb_data
    pep_mat_values$numb_dp_z <- scale(pep_mat_values$numb_dp)[,1]

    pep_mat_values$Var3 <- pep_mat_values$Var1
    prot_values$Var3 <- 'Protein'


    # Transform z-scores into quantiles
    pep_mat_values$FractionCellsExpress <- cut(pep_mat_values$numb_dp, breaks = 3,
                                 labels = c( " Less 33%", "Greater 33%", "Greater 66%"))



    all <- rbind(pep_mat_values,prot_values)

    #pep_mat_values_test <- pep_mat_values

    #pep_mat_values_test$prote <- prot_values$med_abs


    all$cluster <- (as.character(all$cluster))

    if(type == 'line'){
      if(length(all$Var1) > 0){


        plot_ <- ggplot(all, aes(x = cluster, y = med_abs, color = Var3, group = Var3)) +
          geom_line( aes(color = Var3)) +
          geom_point(aes(size = FractionCellsExpress,color = Var3)) +
          scale_size_manual(values = c(" Less 33%" = 2, "Greater 33%" = 4, "Greater 66%" = 6)) +
          geom_errorbar(aes(ymin = med_abs - sd_error, ymax = med_abs + sd_error, color = Var3), width = 0.2) +
          labs(x = "X-axis", y = "Y-axis") +
          ggtitle(paste0(prot,", peptide agreemeent between clusters")) +
          theme_bw() + scale_colour_manual(values = c("red", "blue", "green","purple","black"))+
          xlab('Clusters') +
          ggplot2::facet_wrap(~Var1, scales = "free_y") +
          theme(strip.background = element_blank())+ylab('Log2(Abs)')

        return(plot_)
      }
    }

    if(type == 'bubble'){

      all$var4 <- NA
      all$var4[all$Var3 == 'Protein'] <- 'Protein'
      all$var4[all$Var3 != 'Protein'] <- 'Peptide'

      all$FractionCellsExpress[all$Var3 == 'Protein'] <- 'Greater 66%'
      all$FractionCellsExpress[all$FractionCellsExpress == 'Less 33%'] <- ' Less 33%'

      plot_ <- ggplot(all, aes(x = cluster, y = Var3, fill = med_abs)) +
        geom_point(aes(size = FractionCellsExpress), shape = 21) + theme_bw() +
        scale_size_manual(values = c(" Less 33%" = 3, "Greater 33%" = 5, "Greater 66%" = 8))+
        scale_fill_gradient2(midpoint = 0, low = 'blue',mid = 'white', high = 'red')+
        facet_grid(vars(var4),scales = "free_y", space = "free_y") + ylab('')

      return(plot_)



    }



  }


}


ImputationComparison <- function(nPOP_obj, cluster = 1){
  cluster = 1
  clusters <- nPOP_obj@reductions[['UMAP']]
  clusters$ID <- rownames(clusters)



  prot_imp <- nPOP_obj@protein.imputed
  prot_noimp <- nPOP_obj@protein


  prot_imp <- reshape2::melt(prot_imp)
  prot_imp <- prot_imp %>% left_join(clusters, by = c('Var2' = 'ID'))
  prot_imp <- prot_imp %>% group_by(cluster,Var1) %>% summarise(prot_score = median(value,na.rm=T))
  prot_imp <- prot_imp %>% filter(cluster == 2)

  prot_noimp <- reshape2::melt(prot_noimp)
  prot_noimp <- prot_noimp %>% left_join(clusters, by = c('Var2' = 'ID'))
  prot_noimp <- prot_noimp %>% group_by(cluster,Var1) %>% summarise(prot_score = median(value,na.rm=T))
  prot_noimp <- prot_noimp %>% filter(cluster == 2)


  prot_noimp$imp <- prot_imp$prot_score

  ggplot(prot_noimp, aes(x = prot_score, y = imp)) + geom_point()+ theme_bw()+
    xlab('No Imputation') + ylab('With Imputation') +
      ggtitle(paste0('Cluster ',cluster, ' Protein Averages, Cor = ',
              round(cor(prot_noimp$prot_score,prot_imp$prot_score, use = 'pairwise.complete.obs'),3)))






}



PlotMS1vMS2 <- function(nPOP_obj){

  peps <- nPOP_obj@matricies@peptide
  raw_data <- nPOP_obj@raw_data
  raw_data <- raw_data %>% filter(Channel.Q.Value < nPOP_obj@misc[['ChQ']])

  rownames(peps) <- nPOP_obj@matricies@peptide_protein_map$seqcharge

  raw_data <- raw_data %>% filter(ID %in% colnames(peps))

  MS2_mat <- reshape2::dcast(raw_data, seqcharge ~ ID, value.var = 'Precursor.Quantity')

  rownames(MS2_mat)<- MS2_mat$seqcharge
  MS2_mat$seqcharge <- NULL
  MS2_mat <- as.matrix(MS2_mat)

  sect_pep <- intersect(rownames(MS2_mat),  rownames(peps) )
  sect_col <- intersect(colnames(MS2_mat),  colnames(peps) )

  MS2_mat <- MS2_mat[sect_pep,sect_col]
  peps <- peps[sect_pep,sect_col]

  MS2_mat <- Normalize_reference_vector(MS2_mat, log = T)
  peps <- Normalize_reference_vector(peps, log = T)

  cors <- c()
  pep <- c()
  for(i in 1:nrow(MS2_mat)){
    same <- psych::pairwiseCount(MS2_mat[i,],peps[i,])
    if(same > 10){
      pep <- c(pep, rownames(peps)[i])
      cors <- c(cors, cor(MS2_mat[i,],peps[i,], use = 'pairwise.complete.obs'))
    }

  }


  df <- as.data.frame(cors)
  df$FC <- rowMeans(abs(peps[pep,]),na.rm = T)
  df$FC[df$FC < .4] <- .4
  df$FC[df$FC > .4 & df$FC < .8] <- .8
  df$FC[df$FC > .8 & df$FC < 1.2] <- 1.2
  df$FC[df$FC > 1.2 & df$FC < 1.6] <- 1.6
  df$FC[df$FC > 1.6] <- 2

  df$FC <- as.character(df$FC)


  ggplot(df, aes(x = FC, y = cors))+ geom_boxplot(color="black", fill = 'gray') +
    ggtitle('Ms1, Ms2 Peptide Correlation') + ylab('Correlations') +
    xlab('Average abs(log2 Fold Change)')+ theme_linedraw()+
    stat_summary(fun.data=f, geom="text", vjust=-0.5, col="blue")





}

