### Statistics for LC/MS QC in order samples run

#' Calculate run-order statistics for LC/MS performance monitoring.
#'
#' Computes per-run quality metrics (precursor IDs, intensity drift, retention
#' time drift) across the injection order. Dispatches to DDA- or DIA-specific
#' implementations based on the MS type stored in the QQC object.
#'
#' @param QQC A QQC object with populated \code{raw_data} and \code{ms_type} slots.
#' @return The QQC object with the \code{run_order.statistics} slot populated.
#' @export
Calculate_run_order_statistics <- function(QQC){

  if(QQC@ms_type == 'DDA'){
    QQC <- Calculate_run_order_statistics_DDA(QQC)
  }
  if(QQC@ms_type == 'DIA' | QQC@ms_type == 'DIA_C'){
    QQC <- Calculate_run_order_statistics_DIA(QQC)
  }

  return(QQC)

}

#' Calculate run-order statistics for DDA data.
#'
#' Computes per-run metrics across LC/MS injection order for DDA experiments,
#' including the number of precursor IDs, mean MS1 intensity drift (log2
#' normalized to the first run), mean MS2 reporter intensity drift, and
#' retention time mean and standard deviation drift.
#'
#' @param QQC A QQC object with DDA \code{raw_data} containing columns such as
#'   \code{seqcharge}, \code{Order}, \code{Intensity}, \code{Reporter.intensity.1},
#'   and \code{Retention.time}.
#' @return The QQC object with the \code{run_order.statistics} slot populated as a list
#'   of four data frames: precursor IDs, MS1 means, MS2 means, and RT statistics.
#' @export
Calculate_run_order_statistics_DDA <- function(QQC){


  Raw_data <- QQC@raw_data

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


  QQC@run_order.statistics <- list(IDs_df,MS1_means_df,MS2_means_df,RT_df)

  return(QQC)
}

#' Calculate run-order statistics for DIA data.
#'
#' Computes per-run metrics across LC/MS injection order for DIA experiments,
#' including the number of precursor IDs, total MS1 area (normalized to the
#' first run), total precursor quantity (normalized to the first run), and
#' retention time mean and standard deviation.
#'
#' @param QQC A QQC object with DIA \code{raw_data} containing columns such as
#'   \code{Order}, \code{Ms1.Area}, \code{Precursor.Quantity}, and \code{RT}.
#' @return The QQC object with the \code{run_order.statistics} slot populated as a list
#'   of four data frames: precursor IDs, MS1 means, MS2 means, and RT statistics.
#' @export
Calculate_run_order_statistics_DIA <- function(QQC){


  Raw_data <- QQC@raw_data


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


  QQC@run_order.statistics <- list(Raw_data_IDs,Raw_data_MS1,Raw_data_MS2,RT_df)

  return(QQC)
}

#' Plot intensity drift across LC/MS run order.
#'
#' Generates a three-panel plot showing how precursor identification counts,
#' MS1 intensity, and MS2 intensity change across consecutive LC/MS runs.
#' Useful for diagnosing instrument performance degradation over time.
#'
#' @param QQC A QQC object with a populated \code{run_order.statistics} slot
#'   (from \code{Calculate_run_order_statistics}).
#' @return A patchwork composite of three ggplot panels: precursor IDs, MS1 intensity,
#'   and MS2 intensity versus run order.
#' @export
PlotIntensityDrift <- function(QQC){

  IDs <- ggplot(QQC@run_order.statistics[[1]], aes(x = Run,y = NumbRunIDs)) + geom_point()+dot_plot +
    ggtitle('Run IDs as runs progress') + ylab('# Precursor IDs')


  MS1 <- ggplot(QQC@run_order.statistics[[2]], aes(x = Run,y = MS1_means)) + geom_point()+dot_plot +
    ggtitle('Intersected MS1 Intensity as runs progress') + ylab('Log2( Normalized to run 1)')


  MS2 <- ggplot(QQC@run_order.statistics[[3]], aes(x = Run,y = MS2_means)) + geom_point()+dot_plot+
    ggtitle('MS2 Intensity as runs progress')+ylab('Log2(Normalized to run 1)')


  IDs/MS1/MS2
}

#' Plot retention time drift across LC/MS run order.
#'
#' Generates a two-panel plot showing how mean retention time and retention
#' time standard deviation change across consecutive LC/MS runs. Useful for
#' monitoring chromatographic stability over the course of an experiment.
#'
#' @param QQC A QQC object with a populated \code{run_order.statistics} slot
#'   (from \code{Calculate_run_order_statistics}).
#' @return A patchwork composite of two ggplot panels: mean RT drift and RT
#'   standard deviation versus run order.
#' @export
PlotRTDrift <- function(QQC){

  Mean <-ggplot(QQC@run_order.statistics[[4]], aes(x = Run,y = RT_means)) + geom_point()+dot_plot +
    ggtitle('Mean RT Drift as runs progress (seconds)')

  SDs <- ggplot(QQC@run_order.statistics[[4]], aes(x = Run,y = RT_sds)) + geom_point() +dot_plot+
    ggtitle('RT Standard Deviation as runs progress')

  Mean/SDs

}


