### Statistics for LC/MS QC in order samples run

#' test2.
#'
#' This function takes two numeric inputs and returns their sum.
#'
#' @param x A numeric value.
#' @param y A numeric value.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add_numbers(2, 3)
#' @export
Calculate_run_order_statistics <- function(nPOP_obj){

  if(nPOP_obj@ms_type == 'DDA'){
    nPOP_obj <- Calculate_run_order_statistics_DDA(nPOP_obj)
  }
  if(nPOP_obj@ms_type == 'DIA' | nPOP_obj@ms_type == 'DIA_C'){
    nPOP_obj <- Calculate_run_order_statistics_DIA(nPOP_obj)
  }

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
Calculate_run_order_statistics_DDA <- function(nPOP_obj){


  Raw_data <- nPOP_obj@raw_data

  Raw_data_IDs <- reshape2::dcast(Raw_data,seqcharge ~ Order, value.var = 'Intensity')
  NumbRunIDs <- colSums(is.na(Raw_data_IDs[,2:ncol(Raw_data_IDs)])==F)

  IDs_df <- as.data.frame(NumbRunIDs)
  IDs_df$Run <- as.numeric(names(NumbRunIDs))

  #IDs_df <- IDs_df[match(linker_intersect$Run, IDs_df$Run),]
  #IDs_df$Order <- 1:nrow(IDs_df)



  Raw_data_MS1 <- reshape2::dcast(Raw_data,seqcharge ~ Order, value.var = 'Intensity')
  Raw_data_MS1[,2:ncol(Raw_data_MS1)] <- log2(Raw_data_MS1[,2:ncol(Raw_data_MS1)]/Raw_data_MS1[,2])
  MS1_means <- colMeans(Raw_data_MS1[,2:ncol(Raw_data_MS1)],na.rm = T)

  MS1_means_df <- as.data.frame(MS1_means)
  MS1_means_df$Run <- as.numeric(names(MS1_means))



  Raw_data_MS2 <- reshape2::dcast(Raw_data,seqcharge ~ Order, value.var = 'Reporter.intensity.1')
  Raw_data_MS2[,2:ncol(Raw_data_MS2)] <- log2(Raw_data_MS2[,2:ncol(Raw_data_MS2)]/Raw_data_MS2[,2])
  MS2_means <- colMeans(Raw_data_MS2[,2:ncol(Raw_data_MS2)],na.rm = T)
  MS2_means_df <- as.data.frame(MS2_means)
  MS2_means_df$Run <- as.numeric(names(MS2_means))


  Raw_data_RT <- reshape2::dcast(Raw_data,seqcharge ~ Order, value.var = 'Retention.time')
  Raw_data_RT[,2:ncol(Raw_data_RT)] <- Raw_data_RT[,2:ncol(Raw_data_RT)] - Raw_data_RT[,2]
  RT_means <- colMeans(Raw_data_RT[,2:ncol(Raw_data_RT)],na.rm = T)
  RT_sds <- colSds(as.matrix(Raw_data_RT[,2:ncol(Raw_data_RT)]),na.rm = T)

  RT_df <- as.data.frame(RT_means)
  RT_df$RT_sds <- RT_sds
  RT_df$Run <- as.numeric(names(RT_means))


  nPOP_obj@run_order.statistics <- list(IDs_df,MS1_means_df,MS2_means_df,RT_df)

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
Calculate_run_order_statistics_DIA <- function(nPOP_obj){


  Raw_data <- nPOP_obj@raw_data


  Raw_data_IDs <- Raw_data %>% group_by(Order) %>% dplyr::summarise(NumbRunIDs = sum(is.na(Ms1.Area)==F))
  colnames(Raw_data_IDs) <- c('Run','NumbRunIDs')


  Raw_data_MS1 <- Raw_data %>% group_by(Order) %>% dplyr::summarise(MS1_means = sum(Ms1.Area,na.rm = T))
  min_val <- min(Raw_data_MS1$Order)
  Raw_data_MS1$MS1_means <- Raw_data_MS1$MS1_means / Raw_data_MS1$MS1_means[Raw_data_MS1$Order == min_val]
  colnames(Raw_data_MS1) <- c('Run','MS1_means')


  Raw_data_MS2 <- Raw_data %>% group_by(Order) %>% dplyr::summarise(MS2_means = sum(Precursor.Quantity,na.rm = T))
  min_val <- min(Raw_data_MS2$Order)
  Raw_data_MS2$MS2_means <- Raw_data_MS2$MS2_means / Raw_data_MS2$MS2_means[Raw_data_MS2$Order == min_val]
  colnames(Raw_data_MS2) <- c('Run','MS2_means')

  RT_df <- Raw_data %>% group_by(Order) %>% dplyr::summarise(RT_sds = sd(RT,na.rm = T))
  colnames(RT_df) <- c('Run','RT_sds')

  RT_df_mean <- Raw_data %>% group_by(Order) %>% dplyr::summarise(RT_means = mean(RT,na.rm = T))

  RT_df$RT_means <- RT_df_mean$RT_means


  nPOP_obj@run_order.statistics <- list(Raw_data_IDs,Raw_data_MS1,Raw_data_MS2,RT_df)

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
PlotIntensityDrift <- function(nPOP_obj){

  IDs <- ggplot(nPOP_obj@run_order.statistics[[1]], aes(x = Run,y = NumbRunIDs)) + geom_point()+dot_plot +
    ggtitle('Run IDs as runs progress') + ylab('# Precursor IDs')


  MS1 <- ggplot(nPOP_obj@run_order.statistics[[2]], aes(x = Run,y = MS1_means)) + geom_point()+dot_plot +
    ggtitle('Intersected MS1 Intensity as runs progress') + ylab('Log2( Normalized to run 1)')


  MS2 <- ggplot(nPOP_obj@run_order.statistics[[3]], aes(x = Run,y = MS2_means)) + geom_point()+dot_plot+
    ggtitle('Carrier RI Intensity as runs progress')+ylab('Log2( Normalized to run 1)')


  IDs/MS1/MS2
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
PlotRTDrift <- function(nPOP_obj){

  Mean <-ggplot(nPOP_obj@run_order.statistics[[4]], aes(x = Run,y = RT_means)) + geom_point()+dot_plot +
    ggtitle('Mean RT Drift as runs progress (seconds)')

  SDs <- ggplot(nPOP_obj@run_order.statistics[[4]], aes(x = Run,y = RT_sds)) + geom_point() +dot_plot+
    ggtitle('RT Standard Deviation as runs progress')

  Mean/SDs

}


