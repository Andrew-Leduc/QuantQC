#####
# TMT Processing
#####

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
TMT_Reference_channel_norm <- function(nPOP_obj){

  sc.data <- nPOP_obj@raw_data

  ri.index<-which(colnames(sc.data)%in%paste0("Reporter.intensity.",2:18))


  sc.data[, ri.index] <- sc.data[, ri.index] / sc.data[, ri.index[1]]

  sc.data <- as.data.table(sc.data)

  sc.data <- sc.data[,c('seqcharge','Leading.razor.protein','Raw.file','Well',paste0("Reporter.intensity.",5:18))]
  sc.data <- data.table::melt(sc.data, id = c('seqcharge','Leading.razor.protein','Raw.file','Well'))

  sc.data <- sc.data[sc.data$value < 2.5,]
  sc.data$ID <- paste0(sc.data$Well,sc.data$variable)

  sc.data <- data.table::dcast(sc.data,Leading.razor.protein+seqcharge~ID,value.var = 'value')


  # Add in 0 peptides for negative controls that were totally filtered out
  sc.data[sc.data==0] <- NA

  sc.data <- as.data.frame(sc.data)

  sc.data <- sc.data %>% distinct(seqcharge,.keep_all = T)


  prot_pep_map <- as.data.frame(cbind(sc.data$Leading.razor.protein,sc.data$seqcharge))
  colnames(prot_pep_map) <- c('Protein','seqcharge')
  nPOP_obj@peptide_protein_map <- prot_pep_map

  nPOP_obj@peptide <- as.matrix(sc.data[,3:ncol(sc.data)])

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
inSet_norm <- function(Raw_data, cellenOne_meta){
  count = 0

  for(i in unique(paste0(cellenONE_meta$injectWell))){

    set <- Raw_data[,which(grepl(i,colnames(Raw_data)))]

    Raw_data[,which(grepl(i,colnames(Raw_data)))] <- set - rowMeans(set,na.rm = T)

  }

  return(Raw_data)
}


#####
# mTRAQ Processing
#####


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
remove_sets_carrier_lessthan_SCs <- function(df,files){
  files_remove <- c()
  #df[is.nan(df)] <- NA
  df[df ==0] <- NA
  df[df == Inf] <- NA
  df[df == -Inf] <- NA


  for(i in unique(files)){
    if(length(intersect(colnames(df),(paste0(i,'8')))) > 0){
      count = FALSE
      count1 = FALSE
      count2 = FALSE
      if(length(intersect(colnames(df),(paste0(i,'0')))) > 0){
        count1 = TRUE
        dif1 <- median(as.matrix((df %>% select(paste0(i,'8')))/(df %>% select(paste0(i,'0'))))[,1],na.rm = T)
        if(dif1 > 3){
          count = TRUE
          count1 = TRUE
        }
      }
      if(length(intersect(colnames(df),(paste0(i,'4')))) > 0){
        count2 = TRUE
        dif2 <- median(as.matrix((df %>% select(paste0(i,'8')))/(df %>% select(paste0(i,'4'))))[,1],na.rm = T)
        if(dif2 > 3){
          count = TRUE

        }else{
          count = FALSE
        }
      }

      if(count == FALSE){
        if(count1 == TRUE){
          df <- df %>% select(-paste0(i,'0'))
        }
        if(count2 == TRUE){
          df <- df %>% select(-paste0(i,'4'))
        }

        df <- df %>% select(-paste0(i,'8'))
        files_remove <- c(files_remove,i)
      }

    }else{
      df <- df %>% select(-paste0(i,'4'))
      df <- df %>% select(-paste0(i,'0'))
      files_remove <- c(files_remove,i)
    }
  }
  files <- files[!files %in% files_remove]
  return(list(files,df))

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
carrier_norm <- function(Raw_data_lim.d,carrier_CH,plex_used){

  Carrier_list <- colnames(Raw_data_lim.d)[str_sub(colnames(Raw_data_lim.d), start= -1) == carrier_CH]

  Remove_cells <- c()
  for(i in Carrier_list){
    Flag = FALSE
    i <- str_sub(i, 1, -2)
    for(j in 1:length(plex_used)){
      if(paste0(i,plex_used[j]) %in% colnames(Raw_data_lim.d)){
        if(median(as.matrix((Raw_data_lim.d %>% dplyr::select(paste0(i,plex_used[j])))/(Raw_data_lim.d %>% dplyr::select(paste0(i,carrier_CH)))),na.rm = T) > 1){
          Flag = TRUE
        }
      }
    }
    if(Flag == T){
      Remove_cells <- c(Remove_cells,i)
    }
  }

  non_carrier <- colnames(Raw_data_lim.d)[!colnames(Raw_data_lim.d) %in% Carrier_list]
  for(i in Remove_cells){
    Carrier_list <- Carrier_list[Carrier_list != paste0(i,carrier_CH)]
  }
  for(i in 1:length(plex_used)){
    non_carrier <-non_carrier[non_carrier != paste0(Remove_cells,plex_used[i])]
  }
  all_cols <- c('Protein.Group','seqcharge',non_carrier)
  Carrier_list <- str_sub(Carrier_list, 1, -2)



  for(i in Carrier_list){
    for(j in 1:length(plex_used)){
      if(paste0(i,plex_used[j]) %in% colnames(Raw_data_lim.d)){
        Raw_data_lim.d[,paste0(i,plex_used[j])] <- (Raw_data_lim.d %>% dplyr::select(paste0(i,plex_used[j])))/(Raw_data_lim.d %>% dplyr::select(paste0(i,carrier_CH)))
      }
    }
  }


  Raw_data_lim.d <- Raw_data_lim.d %>% dplyr::select(intersect(all_cols,colnames(Raw_data_lim.d)))

  Raw_data_lim.d[Raw_data_lim.d == Inf] <- NA
  Raw_data_lim.d[Raw_data_lim.d == -Inf] <- NA
  Raw_data_lim.d[Raw_data_lim.d ==0] <- NA

  return(Raw_data_lim.d)
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
cellXgene <- function(Raw_data, carrier, plex_used ,quant){

  Raw_data <- Raw_data %>% filter(plex %in% plex_used)

  Raw_data_lim <- Raw_data %>% dplyr::select(Protein.Group,seqcharge,Ms1.Area,File.Name)
  Raw_data_lim.d <- dcast(Raw_data_lim,Protein.Group+seqcharge~File.Name,value.var = 'Ms1.Area')

  # Normalize data by carrier
  ## This code also removes all sets where a single cell is larger in mean intensity than the carrier
  if(carrier == T){

    Raw_data_lim.d <- carrier_norm(Raw_data_lim.d,carrier_CH,plex_used)

  }

  return(Raw_data_lim.d)

}


#####
#General Proc
#####

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
KNN_impute<-function(nPOP_obj, k = 3){

  sc.data <- nPOP_obj@protein

  # Create a copy of the data, NA values to be filled in later
  sc.data.imp<-sc.data

  # Calculate similarity metrics for all column pairs (default is Euclidean distance)
  dist.mat<-as.matrix( dist(t(sc.data)) )
  #dist.mat<- 1-as.matrix(cor((dat), use="pairwise.complete.obs"))

  #dist.mat<-as.matrix(as.dist( dist.cosine(t(dat)) ))

  # Column names of the similarity matrix, same as data matrix
  cnames<-colnames(dist.mat)

  # For each column in the data...
  for(X in cnames){

    # Find the distances of all other columns to that column
    distances<-dist.mat[, X]

    # Reorder the distances, smallest to largest (this will reorder the column names as well)
    distances.ordered<-distances[order(distances, decreasing = F)]

    # Reorder the data matrix columns, smallest distance to largest from the column of interest
    # Obviously, first column will be the column of interest, column X
    dat.reordered<-sc.data[ , names(distances.ordered ) ]

    # Take the values in the column of interest
    vec<-sc.data[, X]

    # Which entries are missing and need to be imputed...
    na.index<-which( is.na(vec) )

    # For each of the missing entries (rows) in column X...
    for(i in na.index){

      # Find the most similar columns that have a non-NA value in this row
      closest.columns<-names( which( !is.na(dat.reordered[i, ])  ) )

      #print(length(closest.columns))

      # If there are more than k such columns, take the first k most similar
      if( length(closest.columns)>k ){

        # Replace NA in column X with the mean the same row in k of the most similar columns
        vec[i]<-mean( sc.data[ i, closest.columns[1:k] ] )

      }


      # If there are less that or equal to k columns, take all the columns
      if( length(closest.columns)<=k ){

        # Replace NA in column X with the mean the same row in all of the most similar columns
        vec[i]<-mean( sc.data[ i, closest.columns ] )

      }


    }

    # Populate a the matrix with the new, imputed values
    sc.data.imp[,X]<-vec

  }

  # Normalize imputed data
  sc.data.imp <- Normalize_reference_vector_log(sc.data.imp)


  nPOP_obj@protein.imputed <- sc.data.imp

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
Normalize_reference_vector <- function(dat,log = F){
  dat <- as.matrix(dat)
  dat[dat==0] <- NA

  refVec <- matrixStats::rowMedians(x = dat, cols = 1:ncol(dat), na.rm = T)
  for(k in 1:ncol(dat)){
    dat[,k]<-dat[,k] * median(refVec/dat[,k], na.rm = T)
  }
  for(k in 1:nrow(dat)){
    dat[k,]<-dat[k,]/mean(dat[k,], na.rm = T)
  }
  if(log == T){
    dat <- log2(dat)
  }
  return(dat)
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
Normalize_reference_vector_log <- function(dat){
  dat <- as.matrix(dat)

  refVec <- matrixStats::rowMedians(x = dat, cols = 1:ncol(dat), na.rm = T)
  for(k in 1:ncol(dat)){
    dat[,k]<-dat[,k] + median(refVec - dat[,k], na.rm = T)
  }
  for(k in 1:nrow(dat)){
    dat[k,]<-dat[k,]-mean(dat[k,], na.rm = T)
  }
  return(dat)
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
CollapseToProtein <- function(nPOP_obj, opt){
  sc.data <- nPOP_obj@peptide

  # This function colapses peptide level data to the protein level
  # There are different ways to collapse peptides mapping from the same
  # protein to a single data point. The simplest way is to take the median
  # of the peptide values after normalizing. There are also more sophisticated
  # ways that rely on the raw peptide intensities like MaxLFQ


  if(opt == 1){

    # Remove peptide protein name columns
    Normalize_peptide_data <- as.matrix(sc.data)

    # Normalize peptide data for cell size and then to relative abundances and log transform
    Normalize_peptide_data <- Normalize_reference_vector(Normalize_peptide_data,log = T)


    # Remove unwanted values
    Normalize_peptide_data[Normalize_peptide_data == Inf] <- NA
    Normalize_peptide_data[Normalize_peptide_data == -Inf] <- NA

    #Re-Join data
    Normalize_peptide_data <- as.data.table(cbind(nPOP_obj@peptide_protein_map,Normalize_peptide_data))

    # Remove peptides observed less than 10 times
    Normalize_peptide_data <- Normalize_peptide_data[rowSums(is.na(Normalize_peptide_data)==F) > 9,]


    # Collapse peptide levels to median protein level, first melt, then collapse then expand back out
    Normalize_peptide_data <- data.table::melt(Normalize_peptide_data,id.vars = c('Protein','seqcharge'))
    Normalize_peptide_data$seqcharge <- NULL
    Normalize_protein_data <- Normalize_peptide_data[, lapply(.SD, median,na.rm = TRUE), by = c('Protein','variable')]

    # Create Protein x Cell matrix

    Normalize_protein_data <- data.table::dcast(Normalize_protein_data, Protein ~ variable, value.var = 'value')

    Normalize_protein_data <- as.data.frame(Normalize_protein_data)
    rownames(Normalize_protein_data) <- Normalize_protein_data$Protein
    Normalize_protein_data$Protein <- NULL

    # Re-column and row normalize:
    Normalize_protein_data<-Normalize_reference_vector_log(Normalize_protein_data)

    nPOP_obj@protein <- Normalize_protein_data
    return(nPOP_obj)

  }




  if(opt == 2){

    #Max LFQ protein level
    protein_mat <- diann_maxlfq(sc.data,sample.header = "File.Name",group.header = "Protein.Group",id.header = "seqcharge",quantity.header = "Ms1.Area")

    #Normalize protein level data and log transform
    Normalize_protein_data <- Normalize_reference_vector(protein_mat, log = 'yes')

    return(Normalize_protein_data)
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
BatchCorrect <- function(nPOP_obj){

  cellenONE_meta <- nPOP_obj@meta.data
  protein_mat_imputed <- nPOP_obj@protein.imputed
  protein_mat <- nPOP_obj@protein

  # Get meta data for batch correction
  batch_label  <- cellenONE_meta %>% dplyr::filter(ID %in% colnames(protein_mat_imputed))
  batch_label <- batch_label[order(match(batch_label$ID,colnames(protein_mat_imputed))),]

  #linker <- linker %>% dplyr::select(Well,Order)
  #batch_label <- batch_label %>% dplyr::left_join(linker, by =c('injectWell' = 'Well'))


  # Perform batch corrections, possible sources label bias, Every LC/MS runs or groups of LC/MS runs
  #sc.batch_cor <- ComBat(protein_mat_imputed, batch=factor(batch_label$label))
  sc.batch_cor <- limma::removeBatchEffect(protein_mat_imputed,batch = batch_label$injectWell, batch2 = batch_label$label)

  # Re normalize and NA out imputed values
  sc.batch_cor <- Normalize_reference_vector_log(sc.batch_cor)

  # Store unimputed matrix
  sc.batch_cor_noimp <- sc.batch_cor
  sc.batch_cor_noimp[is.na(protein_mat)==T] <- NA

  nPOP_obj@protein.imputed <- sc.batch_cor
  nPOP_obj@protein <- sc.batch_cor_noimp


  return(nPOP_obj)

}
