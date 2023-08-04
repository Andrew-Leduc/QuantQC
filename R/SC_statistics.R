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
  mat_norm <- Normalize_reference_vector(nPOP_obj@peptide)

  mat_norm <- as.data.frame(mat_norm)
  mat_norm$Protein <- nPOP_obj@peptide_protein_map$Protein
  mat_norm$pep <- nPOP_obj@peptide_protein_map$seqcharge

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
SharedPeptideCor <- function(nPOP_obj){
  peptide_data <- nPOP_obj@peptide
  protein_dat <- nPOP_obj@protein
  peptide_protein_map <- nPOP_obj@peptide_protein_map

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

  ggplot(pep_cor, aes(y = Cor, x = FC)) + geom_boxplot() + xlab('Mean abs(protein fold change)') +
    ylab('Correlation between peptides mapping to a protein')+
    stat_summary(fun.data=f, geom="text", vjust=-0.5, col="blue")+
    geom_hline(yintercept = null_dist, col = "red")




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

  # count peptide numbers
  numb_pep <- colSums(is.na(nPOP_obj@peptide)==F)
  numb_pep <- as.data.frame(numb_pep)
  colnames(numb_pep) <- 'Number_precursors'

  # count protein numbers
  numb_prot <- colSums(is.na(nPOP_obj@protein)==F)
  numb_prot <- as.data.frame(numb_prot)
  colnames(numb_prot) <- 'Number_proteins'

  # Plot peptide and protein numbers
  pep_number <- ggplot(numb_pep, aes(x = Number_precursors)) + geom_histogram(bins = 30,position = 'identity',alpha = .5) + ggtitle('# precursors per sample') + rremove('legend')+ylab('# of single cells')+dot_plot

  prot_number<- ggplot(numb_prot, aes(x = Number_proteins)) + geom_histogram(bins = 30,position = 'identity',alpha = .5) + ggtitle('# proteins per sample') + ylab('# of single cells')+dot_plot

  pep_number+prot_number
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

  data <- nPOP_obj@protein

  missingness_cell_filt <- colSums(is.na(data) == F)/nrow(data)
  missingness_prot_filt <- rowSums(is.na(data)==F)/ncol(data)

  missingness_prot_mat <- (missingness_prot_filt)
  missingness_prot_mat <- reshape2::melt(missingness_prot_mat)

  missingness_cell_mat <- missingness_cell_filt
  missingness_cell_mat <- reshape2::melt(missingness_cell_mat)

  mp <- ggplot(missingness_prot_mat, aes(x = value)) +
    geom_histogram(bins = 20,position = 'identity',alpha = .5) + ggtitle(paste0('Protein completness, ', nrow(protein_mat) ,' proteins'))+rremove('legend') +ylab('# of proteins')+xlab('fraction values present')+dot_plot

  mc <- ggplot(missingness_cell_mat, aes(x = value)) +
    geom_histogram(bins = 20,position = 'identity',alpha = .5) + ggtitle('Cell completness') + ylab('# of single cells')+
    xlab('fraction values present')+dot_plot


  mp+mc


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

  if(nPOP_obj@type = "DDA"){

    sc.data <- nPOP_obj@raw_data
    cellenOne_meta <- nPOP_obj@meta.data

    ri.index<-which(colnames(sc.data)%in%paste0("Reporter.intensity.",1:18))

    sc.data[sc.data==0] <- NA
    sc.data[, ri.index] <- sc.data[, ri.index] / sc.data[, ri.index[1]]


    sc.data <- sc.data[,c('Raw.file','Well',paste0("Reporter.intensity.",5:18))]
    sc.data <- reshape2::melt(sc.data, id = c('Raw.file','Well'))

    sc.data <- sc.data %>% filter(value < 2)
    sc.data$ID <- paste0(sc.data$Well,sc.data$variable)

    sc.data <- sc.data %>% group_by(ID) %>% dplyr::summarise(value = median(value,na.rm = T))

    sc.data <- sc.data %>% dplyr::left_join(cellenOne_meta, by = c('ID'))

    sc.data$value <- 1/sc.data$value

    ggplot(sc.data, aes(x = sample, y = value)) + geom_boxplot() + ylab('# cells / carrier amount') +dot_plot

  }

  if(nPOP_obj@type == "DIA_C"){


  }
  if(nPOP_obj@type == "DIA"){
    return("No carrier used")
  }


}
