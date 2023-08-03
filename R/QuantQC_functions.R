# Welcome to the QuantQC functions,
  # Code is carefully annotated to explain the data analysis process

  # For easy navigation option+command+0 to fold all functions on Mac, Alt+0 Windows

# nPOP Object

nPOP <- setClass(
  Class = 'nPOP',
  slots = c(
    ms_type = 'character',
    plex = 'numeric',
    raw_data = 'data.frame',
    peptide = 'matrix',
    protein = 'matrix',
    protein.imputed = 'matrix',
    peptide_protein_map = 'data.frame',
    cellenONE.meta = 'data.frame',
    meta.data = 'data.frame',
    run_order.statistics = 'list',
    pep.cor = 'list',
    neg_ctrl.info = 'data.frame',
    graphs = 'list',
    neighbors = 'list',
    reductions = 'list',
    images = 'list',
    project.name = 'character',
    misc = 'list'
  )
)



#Fasta Paths for gene name

Mouse <- 'mouse_convert.tsv'
Human <- 'human_convert.tsv'



#themes

my_col3 <- c("purple2","black")

dot_plot <-  theme_bw()+theme(plot.title = element_text(hjust = .5,size = 24),
                              axis.title.x = element_text(size = 20),
                              axis.title.y = element_text(size = 20),
                              axis.text.x = element_text(size = 12),
                              axis.text.y = element_text(size = 12))






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
  if(nPOP_obj@ms_type == 'DIA'){
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
  #MS1_means_df <- MS1_means_df[match(linker_intersect$Run, MS1_means_df$Run),]
  #MS1_means_df$Order <- 1:nrow(MS1_means_df)


  Raw_data_MS2 <- reshape2::dcast(Raw_data,seqcharge ~ Order, value.var = 'Reporter.intensity.1')
  Raw_data_MS2[,2:ncol(Raw_data_MS2)] <- log2(Raw_data_MS2[,2:ncol(Raw_data_MS2)]/Raw_data_MS2[,2])
  MS2_means <- colMeans(Raw_data_MS2[,2:ncol(Raw_data_MS2)],na.rm = T)
  MS2_means_df <- as.data.frame(MS2_means)
  MS2_means_df$Run <- as.numeric(names(MS2_means))
  #MS2_means_df <- MS2_means_df[match(linker_intersect$Run, MS2_means_df$Run),]
  #MS2_means_df$Order <- 1:nrow(MS2_means_df)

  Raw_data_RT <- reshape2::dcast(Raw_data,seqcharge ~ Order, value.var = 'Retention.time')
  Raw_data_RT[,2:ncol(Raw_data_RT)] <- Raw_data_RT[,2:ncol(Raw_data_RT)] - Raw_data_RT[,2]
  RT_means <- colMeans(Raw_data_RT[,2:ncol(Raw_data_RT)],na.rm = T)
  RT_sds <- colSds(as.matrix(Raw_data_RT[,2:ncol(Raw_data_RT)]),na.rm = T)

  RT_df <- as.data.frame(RT_means)
  RT_df$RT_sds <- RT_sds
  RT_df$Run <- as.numeric(names(RT_means))
  #RT_df <- RT_df[match(linker_intersect$Run, RT_df$Run),]
  #RT_df$Order <- 1:nrow(RT_df)

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

  Raw_data_IDs <- reshape2::dcast(Raw_data,seqcharge ~ Order, value.var = 'Intensity')
  NumbRunIDs <- colSums(is.na(Raw_data_IDs[,2:ncol(Raw_data_IDs)])==F)

  IDs_df <- as.data.frame(NumbRunIDs)
  #IDs_df$Run <- names(NumbRunIDs)
  #IDs_df <- IDs_df[match(linker_intersect$Run, IDs_df$Run),]
  #IDs_df$Order <- 1:nrow(IDs_df)



  Raw_data_MS1 <- reshape2::dcast(Raw_data,seqcharge ~ Raw.file, value.var = 'Intensity')
  Raw_data_MS1[,2:ncol(Raw_data_MS1)] <- log2(Raw_data_MS1[,2:ncol(Raw_data_MS1)]/Raw_data_MS1[,2])
  MS1_means <- colMeans(Raw_data_MS1[,2:ncol(Raw_data_MS1)],na.rm = T)

  MS1_means_df <- as.data.frame(MS1_means)
  #MS1_means_df$Run <- names(MS1_means)
  #MS1_means_df <- MS1_means_df[match(linker_intersect$Run, MS1_means_df$Run),]
  #MS1_means_df$Order <- 1:nrow(MS1_means_df)


  Raw_data_MS2 <- reshape2::dcast(Raw_data,seqcharge ~ Raw.file, value.var = 'Reporter.intensity.1')
  Raw_data_MS2[,2:ncol(Raw_data_MS2)] <- log2(Raw_data_MS2[,2:ncol(Raw_data_MS2)]/Raw_data_MS2[,2])
  MS2_means <- colMeans(Raw_data_MS2[,2:ncol(Raw_data_MS2)],na.rm = T)
  MS2_means_df <- as.data.frame(MS2_means)
  #MS2_means_df$Run <- names(MS2_means)
  #MS2_means_df <- MS2_means_df[match(linker_intersect$Run, MS2_means_df$Run),]
  #MS2_means_df$Order <- 1:nrow(MS2_means_df)


  nPOP_obj@run_order.statistics <- list(IDs_df,MS1_means_df,MS2_means_df)

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


### TMT only Functions

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
MQ_to_nPOP <- function(path1,linker,PIF_in,PEP_in){


  columns_to_read <- c('Modified sequence','Intensity','Retention time','Charge','Raw file','PEP','PIF', 'Leading razor protein',
                       'Potential contaminant','Reverse', paste0("Reporter intensity ",1:18))

  # Read in file
  data1 <- data.table::fread(path1,select = columns_to_read)
  #data2 <- data.table::fread(path2,select = columns_to_read)
  data <- data1

  # Remove spaces from column names
  colnames(data) <- str_replace_all(colnames(data),' ','.')

  # Convert from data.table to data.frame
  data <- as.data.frame(data)

  # Filter for only raw files in the linker
  data <- data %>% filter(Raw.file %in% linker$Run)

  # Denote run order
  linker$Order <- 1:nrow(linker)

  #Link data to inject wells
  data <- data %>% left_join(linker, by = c('Raw.file'= 'Run'))

  # Unique precursor ID
  data$seqcharge <- paste0(data$Modified.sequence,data$Charge)
  data$seqRun <- paste0(data$seqcharge, data$Raw.file)
  data <- data %>% distinct(seqRun,.keep_all = T)


  # Clean up Leading razor protein strings
  parse_row<-grep("|",data$Leading.razor.protein, fixed=T)
  if(length(parse_row)>0){
    split_prot<-str_split(data$Leading.razor.protein[parse_row], pattern = fixed("|"))
    split_prot2<-unlist(split_prot)[seq(2,3*length(split_prot),3)]
    data$Leading.razor.protein[parse_row]<-split_prot2
  }
  parse_row<-grep("-",data$Leading.razor.protein, fixed=T)
  if(length(parse_row)>0){
    split_prot<-str_split(data$Leading.razor.protein[parse_row], pattern = fixed("-"))
    split_prot2<-unlist(split_prot)[seq(1,2*length(split_prot),2)]
    data$Leading.razor.protein[parse_row]<-split_prot2
  }


  # Filter data
  data<-data %>% filter(PEP < PEP_in)
  data<-data %>% filter(PIF > PIF_in)
  data<-data %>% filter(Potential.contaminant != '+')
  data<-data %>% filter(Reverse != '+')



  nPOP_obj <- new('nPOP',raw_data = data,ms_type = 'DDA')


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
analyzeCellenONE_TMT <- function(allDays){

  # Code to parse cellenONE files and map cell diameters, a mess and not too important,
  # dont feel obligeted to read


  #file_paths
  labelPath <- system.file("extdata", "14plex_files/Labels.fld", package = "QuantQC")
  pickupPath1 <-  system.file("extdata", "14plex_files/Pickup_mock.fld", package = "QuantQC")


  allDays[grepl("Transmission",allDays$X),]$X <- NA
  allDays <- allDays %>% fill(2:7, .direction = "up") %>% drop_na(XPos)

  #### Labelling and Pickup Field Files
  ## labelling file

  con_lab <-file(labelPath)
  lines_lab<- readLines(con_lab)
  close(con_lab)
  slines_lab <- strsplit(lines_lab,"\t")
  colCount_lab <- max(unlist(lapply(slines_lab, length)))

  label <- read.table(labelPath, sep="\t",fill=TRUE,header = F,col.names=c("position", "well","volume","field"), quote = "", skip = 22)

  fieldRaw_label <- label[grepl("\\[\\d",label$position),]$position
  fieldBrack_label <- gsub("\\[|\\]", "", fieldRaw_label)
  fieldComma_label <- strsplit(fieldBrack_label, "," )
  fieldNum_label <- unlist(lapply(fieldComma_label, `[[`, 1))
  label$field[grepl("\\[\\d",label$position)] <- fieldNum_label
  label <- label %>% fill(field, .direction = "down")
  label <-  label[!label$well=="",]
  label$field <- as.numeric(label$field) + 1

  labelxyPos <- strsplit(label$position, "\\/")
  label$yPos <- unlist(lapply(labelxyPos, '[[', 1))
  label$xPos <- unlist(lapply(labelxyPos, '[[', 2))

  ## sample pickup file

  con_pickup <-file(pickupPath1)
  lines_pickup<- readLines(con_pickup)
  close(con_pickup)
  slines_pickup <- strsplit(lines_pickup,"\t")
  colCount_pickup <- max(unlist(lapply(slines_pickup, length)))

  pickup <- read.table(pickupPath1, sep="\t",fill=TRUE,header = F,col.names=c("position", "well","volume","field"), quote = "", skip = 22)


  fieldRaw_pickup <- pickup[grepl("\\[\\d",pickup$position),]$position
  fieldBrack_pickup <- gsub("\\[|\\]", "",fieldRaw_pickup)
  fieldComma_pickup <- strsplit(fieldBrack_pickup, "," )
  fieldNum_pickup <- unlist(lapply(fieldComma_pickup, `[[`, 1))
  pickup$field[grepl("\\[\\d",pickup$position)] <- fieldNum_pickup
  pickup <- pickup %>% fill(field, .direction = "down")
  pickup <-  pickup[!pickup$well=="",]
  pickup$field <- as.numeric(pickup$field) + 1

  ## trying to get pickup working in the same way
  # fixing pickup file
  pickup$target <- 1

  pickup$well <- substring(pickup$well, 2)
  pickup$well <- gsub(",","",pickup$well)



  pickupxyPos <- strsplit(pickup$position, "\\/")
  pickup$yPos <- unlist(lapply(pickupxyPos, '[[', 1))
  pickup$xPos <- unlist(lapply(pickupxyPos, '[[', 2))
  pickup$xPos <- as.numeric(pickup$xPos)
  pickup$yPos <- as.numeric(pickup$yPos)

  ### making sure that label "slide" refers to fields or actual slides.
  ## There is a sequence of 108 x,y and there are 108 spots per slide
  # order is y/x positions

  label$yPos <- as.numeric(label$yPos)
  label$xPos <- as.numeric(label$xPos)
  label$well <- substring(label$well, 3)
  label$well <- as.numeric(gsub(",","",label$well))
  matchTMTSCP <- paste0("TMT", 1:length(unique(label$well)))

  for (i in 1:length(matchTMTSCP)) {

    label[which(label$well == i),]$well <- matchTMTSCP[i]

  }

  label$well <- substring(label$well, 4)


  label <- label %>% filter(field %in% unique(allDays$Field))
  pickup <- pickup %>% filter(field %in% unique(allDays$Field))

  #### Trying to map label to cell
  ###  lets try and keep all three important in one
  allDays <- transform(allDays, xyf = paste0(XPos, YPos, Field))
  label <- transform(label, xyf = paste0(xPos, yPos, field))


  ### Isolation and Label merged

  isoLab <- allDays %>% group_by(Target) %>% left_join(label, by = 'xyf')
  labelCount <- isoLab %>% group_by(well) %>% dplyr::summarize(count=n())


  ### clustering together pickup points with sample points using ANN
  ## super convolutedly/idiotically

  isoLab$ann <- NA
  isoLab$pickupX <- NA
  isoLab$pickupY <- NA

  ann_ <- ann(ref = as.matrix(unique(pickup[, c("xPos","yPos")])),  target = as.matrix(isoLab[ , c("xPos","yPos")]), k=1)
  isoLab$ann <-  ann_$knnIndexDist[,1]
  isoLab_new <- unique(pickup[, c("xPos","yPos")])


  isoLab$pickupX <- isoLab_new[isoLab$ann,]$xPos
  isoLab$pickupY <- isoLab_new[isoLab$ann,]$yPos
  isoLab_bound <-isoLab



  ### Merge pickup and isoLab
  isoLab_bound <- transform(isoLab_bound, xyft = paste0(pickupX, pickupY, Field, Target))
  pickup <-  transform(pickup, xyft = paste0(xPos, yPos, field, target))


  isoLab_final <- isoLab_bound %>% left_join(pickup, by = 'xyft')
  wellCount <- isoLab_final %>% group_by(well.y) %>% dplyr::summarize(count=n())


  ### Clean up to yield final dataframe
  cellenOne_data <- data.frame(sample = isoLab_final$condition, isoTime = isoLab_final$Time, diameter = isoLab_final$Diameter, elongation = isoLab_final$Elongation, slide = isoLab_final$Target, field = isoLab_final$Field, dropXPos = isoLab_final$XPos, dropYPos = isoLab_final$YPos, label = isoLab_final$well.x, pickupXPos = isoLab_final$pickupX, pickupYPos = isoLab_final$pickupY, injectWell = isoLab_final$well.y)


  ##*sigh* not done yet

  cellenOne_data$wellAlph <- substring(cellenOne_data$injectWell, 1, 1)
  cellenOne_data$wellNum <- substring(cellenOne_data$injectWell, 2)

  cellenOne_data <- cellenOne_data %>% arrange(wellAlph, as.numeric(wellNum), as.numeric(label))
  #cellenOne_data <- cellenOne_data %>% dplyr::select(dropYPos,dropXPos,sample,diameter,elongation,field,label,injectWell)

  cellenOne_data$pickupXPos_numb <- as.numeric(cellenOne_data$pickupXPos)
  cellenOne_data$pickupYPos_numb  <- as.numeric(cellenOne_data$pickupYPos)


  # Assigning TMT tags to wells where label was picked up out of during prep
  # Each tag was dispensed multiple times so maps to multiple wells of plate
  cellenOne_data$label <- paste0('Reporter.intensity.' , (as.numeric(cellenOne_data$label) + 4))


  # Create cell ID to match MaxQuant report
  cellenOne_data$ID <- paste0(cellenOne_data$injectWell,cellenOne_data$label)




  return(cellenOne_data)


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
analyzeCellenONE_mTRAQ <- function(allDays,plex){
  #plex = 3
  #allDays = all_cells


  # Code to parse cellenONE files and map cell diameters, a mess and not too important,
  # dont feel obligeted to read

  if(plex == 2){
    #File paths to pickup/label files

    # 2plex
    labelPath <- system.file("extdata", "2plex_files/Labels.fld", package = "QuantQC")
    pickupPath1 <- system.file("extdata", "2plex_files/Pickup_1_mock.fld", package = "QuantQC")
    pickupPath2 <- system.file("extdata", "2plex_files/Pickup_2_mock.fld", package = "QuantQC")

  }
  if(plex == 3){

    # 3plex
    labelPath <- system.file("extdata", "3plex_files/Labels.fld", package = "QuantQC")
    pickupPath1 <- system.file("extdata", "3plex_files/Pickup_1_mock.fld", package = "QuantQC")
    pickupPath2 <- system.file("extdata", "3plex_files/Pickup_2_mock.fld", package = "QuantQC")
  }


  allDays[grepl("Transmission",allDays$X),]$X <- NA
  allDays <- allDays %>% fill(2:7, .direction = "up") %>% drop_na(XPos)

  #### Labelling and Pickup Field Files
  ## labelling file

  con_lab <-file(labelPath)
  lines_lab<- readLines(con_lab)
  close(con_lab)
  slines_lab <- strsplit(lines_lab,"\t")
  colCount_lab <- max(unlist(lapply(slines_lab, length)))

  label <- read.table(labelPath, sep="\t",fill=TRUE,header = F,col.names=c("position", "well","volume","field"), quote = "", skip = 22)

  fieldRaw_label <- label[grepl("\\[\\d",label$position),]$position
  fieldBrack_label <- gsub("\\[|\\]", "", fieldRaw_label)
  fieldComma_label <- strsplit(fieldBrack_label, "," )
  fieldNum_label <- unlist(lapply(fieldComma_label, `[[`, 1))
  label$field[grepl("\\[\\d",label$position)] <- fieldNum_label
  label <- label %>% fill(field, .direction = "down")
  label <-  label[!label$well=="",]
  label$field <- as.numeric(label$field) + 1

  labelxyPos <- strsplit(label$position, "\\/")
  label$yPos <- unlist(lapply(labelxyPos, '[[', 1))
  label$xPos <- unlist(lapply(labelxyPos, '[[', 2))

  ## sample pickup file

  con_pickup <-file(pickupPath1)
  lines_pickup<- readLines(con_pickup)
  close(con_pickup)
  slines_pickup <- strsplit(lines_pickup,"\t")
  colCount_pickup <- max(unlist(lapply(slines_pickup, length)))

  pickup <- read.table(pickupPath1, sep="\t",fill=TRUE,header = F,col.names=c("position", "well","volume","field"), quote = "", skip = 22)



  con_pickup <-file(pickupPath2)
  lines_pickup<- readLines(con_pickup)
  close(con_pickup)
  slines_pickup <- strsplit(lines_pickup,"\t")
  colCount_pickup <- max(unlist(lapply(slines_pickup, length)))
  pickup2 <- read.table(pickupPath2, sep="\t",fill=TRUE,header = F,col.names=c("position", "well","volume","field"), quote = "", skip = 27)

  pickup <- rbind(pickup,pickup2)

  fieldRaw_pickup <- pickup[grepl("\\[\\d",pickup$position),]$position
  fieldBrack_pickup <- gsub("\\[|\\]", "",fieldRaw_pickup)
  fieldComma_pickup <- strsplit(fieldBrack_pickup, "," )
  fieldNum_pickup <- unlist(lapply(fieldComma_pickup, `[[`, 1))
  pickup$field[grepl("\\[\\d",pickup$position)] <- fieldNum_pickup
  pickup <- pickup %>% fill(field, .direction = "down")
  pickup <-  pickup[!pickup$well=="",]
  pickup$field <- as.numeric(pickup$field) + 1



  pickupxyPos <- strsplit(pickup$position, "\\/")
  pickup$yPos <- unlist(lapply(pickupxyPos, '[[', 1))
  pickup$xPos <- unlist(lapply(pickupxyPos, '[[', 2))


  ### making sure that label "slide" refers to fields or actual slides.
  ## There is a sequence of 108 x,y and there are 108 spots per slide
  # order is y/x positions
  LabelxyPos <- strsplit(label$position, "\\/")
  label$yPos <- unlist(lapply(LabelxyPos, '[[', 1))
  label$xPos <- unlist(lapply(LabelxyPos, '[[', 2))
  label$well <- substring(label$well, 3)
  label$well <- as.numeric(gsub(",","",label$well))
  matchTMTSCP <- paste0("TMT", 1:length(unique(label$well)))

  for (i in 1:length(matchTMTSCP)) {

    label[which(label$well == i),]$well <- matchTMTSCP[i]

  }

  label$well <- substring(label$well, 4)


  label <- label %>% filter(field %in% unique(allDays$Field))
  pickup <- pickup %>% filter(field %in% unique(allDays$Field))


  #### Trying to map label to cell
  ###  lets try and keep all three important in one
  allDays <- transform(allDays, xyf = paste0(XPos, YPos, Field))
  label <- transform(label, xyf = paste0(xPos, yPos, field))

  ### Isolation and Label merged

  isoLab <- allDays %>% group_by(Target) %>% left_join(label, by = 'xyf')
  labelCount <- isoLab %>% group_by(well) %>% dplyr::summarize(count=n())

  ## trying to get pickup working in the same way
  # fixing pickup file
  pickup$target <- 1

  pickup$well <- substring(pickup$well, 2)
  pickup$well <- gsub(",","",pickup$well)

  #pickup <- pickup %>% filter(!well %in% No_use$Wells_no )


  slidesUsedForPrep <- unique(isoLab$Target)

  #pickup <- pickup[which(pickup$target %in% slidesUsedForPrep),]

  ### clustering together pickup points with sample points using ANN
  ## super convolutedly/idiotically

  isoLab$ann <- NA
  isoLab$pickupX <- NA
  isoLab$pickupY <- NA

  #ann_123 <- ann(ref = as.matrix(unique(pickup[(-which(pickup$field == 4)) , c("xPos","yPos")])),  target = as.matrix(isoLab[(-which(isoLab$Field == 4)) , c("xPos","yPos")]), k=1)

  ann_ <- ann(ref = as.matrix(unique(pickup[, c("xPos","yPos")])),  target = as.matrix(isoLab[ , c("xPos","yPos")]), k=1)

  #ann_4 <- ann(ref = as.matrix(unique(pickup[(which(pickup$field == 4)) , c("xPos","yPos")])),  target = as.matrix(isoLab[(which(isoLab$Field == 4)) , c("xPos","yPos")]), k=1)

  #isoLab[(-which(isoLab$Field == 4)),]$ann <-  ann_123$knnIndexDist[,1]
  #isoLab[(which(isoLab$Field == 4)),]$ann <-  ann_4$knnIndexDist[,1]

  isoLab$ann <-  ann_$knnIndexDist[,1]


  ## split - combine
  #isoLab_123 <- isoLab[-which(isoLab$Field == 4),]
  #isoLab_4 <- isoLab[which(isoLab$Field == 4),]

  isoLab_new <- unique(pickup[, c("xPos","yPos")])
  #notFieldFourPickUnique <- unique(pickup[-(which(pickup$field == 4)) , c("xPos","yPos")])
  #fieldFourPickUnique <- unique(pickup[(which(pickup$field == 4)) , c("xPos","yPos")])

  #isoLab_123$pickupX <- notFieldFourPickUnique[isoLab_123$ann,]$xPos
  #isoLab_123$pickupY <- notFieldFourPickUnique[isoLab_123$ann,]$yPos

  #isoLab_4$pickupX <- fieldFourPickUnique[isoLab_4$ann,]$xPos
  #isoLab_4$pickupY <- fieldFourPickUnique[isoLab_4$ann,]$yPos

  isoLab$pickupX <- isoLab_new[isoLab$ann,]$xPos
  isoLab$pickupY <- isoLab_new[isoLab$ann,]$yPos

  isoLab_bound <-isoLab
  #isoLab_bound <- rbind(isoLab_123, isoLab_4)



  ### Merge pickup and isoLab
  isoLab_bound <- transform(isoLab_bound, xyft = paste0(pickupX, pickupY, Field, Target))
  pickup <-  transform(pickup, xyft = paste0(xPos, yPos, field, target))

  #intersect(pickup$xyft,isoLab_bound$xyft)

  isoLab_final <- isoLab_bound %>% left_join(pickup, by = 'xyft')
  wellCount <- isoLab_final %>% group_by(well.y) %>% dplyr::summarize(count=n())



  ### Clean up to yield final dataframe
  cellenOne_data <- data.frame (sample = isoLab_final$condition, isoTime = isoLab_final$Time, diameter = isoLab_final$Diameter, elongation = isoLab_final$Elongation, slide = isoLab_final$Target, field = isoLab_final$Field, dropXPos = isoLab_final$XPos, dropYPos = isoLab_final$YPos, label = isoLab_final$well.x, pickupXPos = isoLab_final$pickupX, pickupYPos = isoLab_final$pickupY, injectWell = isoLab_final$well.y)


  ##*sigh* not done yet

  cellenOne_data$wellAlph <- substring(cellenOne_data$injectWell, 1, 1)
  cellenOne_data$wellNum <- substring(cellenOne_data$injectWell, 2)

  cellenOne_data <- cellenOne_data %>% arrange(wellAlph, as.numeric(wellNum), as.numeric(label))
  #cellenOne_data <- cellenOne_data %>% dplyr::select(dropYPos,dropXPos,sample,diameter,elongation,field,label,injectWell)
  cellenOne_data$pickupXPos_numb <- as.numeric(cellenOne_data$pickupXPos)
  cellenOne_data$pickupYPos_numb  <- as.numeric(cellenOne_data$pickupYPos)


  # if using plex/scopeDIA samples are injected into 2 different plates,
  # there are 16 fields and samples from 1-8 go into plate 1, 9-16 go into plate 2
  cellenOne_data$plate <- NA
  cellenOne_data$plate[cellenOne_data$field > 8] <- 2
  cellenOne_data$plate[cellenOne_data$field <= 8] <- 1


  # Assigning mTRAQ tags to wells where label was picked up out of during prep
  # Each tag was dispensed multiple times so maps to multiple wells of plate
  cellenOne_data$tag <- NA

  if(plex == 2){
    cellenOne_data$tag[cellenOne_data$label %in% c('1','3','5','7')] <- '0'
    cellenOne_data$tag[cellenOne_data$label %in% c('2','4','6','8')] <- '4'
  }
  if(plex == 3){
    cellenOne_data$tag[cellenOne_data$label %in% c('1','4','7','10')] <- '0'
    cellenOne_data$tag[cellenOne_data$label %in% c('2','5','8','11')] <- '4'
    cellenOne_data$tag[cellenOne_data$label %in% c('3','6','9','12')] <- '8'
  }

  cellenOne_data$label <- cellenOne_data$tag
  # Create cell ID to match DIA-NN report
  cellenOne_data$ID <- paste0(cellenOne_data$injectWell,cellenOne_data$plate,'.',cellenOne_data$tag)

  return(cellenOne_data)
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
link_cellenONE_Raw <- function(nPOP_obj,allDays){

  if(nPOP_obj@ms_type == 'DDA'){
    cellenOne_data <- analyzeCellenONE_TMT(allDays)
  }
  if(nPOP_obj@ms_type == 'DIA'){
    cellenOne_data <- analyzeCellenONE_TMT(allDays,nPOP_obj@plex)
  }



  peptide_data <- nPOP_obj@peptide
  # Get list of unique cell IDs
  cellID <- colnames(peptide_data)#[1:ncol(peptide_data)]
  cellID <- as.data.frame(cellID)
  colnames(cellID) <- 'ID'

  cellenOne_data_small <- cellenOne_data %>% dplyr::select(any_of(c('ID','diameter','sample','label','injectWell','plate')))
  cellenOne_data_small <- as.data.frame(cellenOne_data_small)


  cellID <- cellID %>% left_join(cellenOne_data_small,by = c('ID'))

  cellID$sample[is.na(cellID$sample)==T] <- 'neg'

  cellID$prot_total <- log2(colSums(peptide_data[,1:ncol(peptide_data)],na.rm = T))

  nPOP_obj@cellenONE.meta <- cellenOne_data

  nPOP_obj@meta.data <- cellID

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
PlotSCtoCarrierRatio <- function(nPOP_obj){
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


### DIA only functions

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
DIANN_to_nPOP <- function(path,linker,plex){

  columns_to_read <-c('Genes','Run','Lib.PG.Q.Value','Precursor.Id','Stripped.Sequence',
                  'Precursor.Charge','Ms1.Area','Protein.Group','Channel.Q.Value')

  Raw_data <- data.table::fread(path,select = columns_to_read)

  Raw_data <- as.data.frame(Raw_data)

  #Raw_data <- Raw_data %>% filter(Lib.PG.Q.Value < .01)
  Raw_data <- Raw_data %>% filter(Run %in% linker$Run)
  Raw_data <- Raw_data %>% left_join(linker, by = c('Run'))

  #Unique precursor ID
  Raw_data$seqcharge <- paste0(Raw_data$Stripped.Sequence,Raw_data$Precursor.Charge)
  Raw_data <- Raw_data %>% filter(Protein.Group != '')


  # this grabs the mTRAQ tag used (may need to be adjusted if using different multiplexing tag)
  Raw_data$plex <- substr(Raw_data$Precursor.Id[1:nrow(Raw_data)], 10, 10)

  # Unique cell ID
  Raw_data$ID <- paste0(Raw_data$Well,Raw_data$Plate,'.',Raw_data$plex)
  Raw_data$File.Name <- Raw_data$ID

  #Remove redundant data points
  Raw_data$uq <- paste0(Raw_data$File.Name,Raw_data$Protein.Group,Raw_data$seqcharge)
  Raw_data <- Raw_data %>% distinct(uq,.keep_all = T)
  Raw_data$uq <- NULL

  nPOP_obj <- new('nPOP',raw_data = Raw_data,type = 'DIA')

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
EvaluateNegativeControls <- function(nPOP_obj,CV_thresh){

  # Compute CVs of cells and negative controls, function outputs CV plot and list of cells with CVs
  CVm <- CVs(nPOP_obj,CV_thresh)


  # Get IDs of cells with median protein CVs (good_cells is a df with cell ID and CV)
  good_cells <- CVm[[2]] %>% filter(cvq < CV_thresh)


  if(length(good_cells) < 3){

    return('less than three good cells, try increasing CV filter')

  }


  # Count the number of peptides in negative controls and good single cells
  Peptide_counts_by_sample <- Count_peptides_per_cell(nPOP_obj@peptide,nPOP_obj@meta.data,good_cells)


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
PlotNegCtrl <- function(nPOP_obj,thresh){

  plot_data <- nPOP_obj@neg_ctrl.info

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
    annotate("text", x=0.2, y= 14, label=paste0(sum(CV_mat_pos$cvq < thresh)," cells"), size=10, color=my_col3[c(1)])+
    annotate("text", x=0.64, y= 12, label=paste0(sum(CV_mat_neg$cvq > thresh,na.rm = T)," Ctr -"), size=10, color=my_col3[c(2)])+
    annotate("text", x=0.63, y= 14, label=paste0(sum(CV_mat_pos$cvq > thresh)," cells"), size=10, color=my_col3[c(1)])+
    annotate("text", x=0.2, y= 12, label=paste0((sum(CV_mat_neg$cvq < thresh,na.rm = T)-1)," Ctr -"), size=10, color=my_col3[c(2)])+
    ggtitle('Cells need atleast 3 proteins with multiple peptides')+
    rremove("legend") + geom_vline(xintercept=thresh, lty=2, size=2, color="gray50") + theme(plot.margin = margin(1, 1, 0, 1, "cm"))


  peps+cvs

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
FilterBadCells <- function(nPOP_obj,CV_thresh){
  neg_filter <- nPOP_obj@neg_ctrl.info
  peptide_data <- nPOP_obj@peptide
  neg_filter <- neg_filter %>% dplyr::filter(cvq < CV_thresh)
  neg_filter <- neg_filter %>% dplyr::filter(value != 'neg')
  peptide_data <- peptide_data[,colnames(peptide_data) %in% neg_filter$variable]

  nPOP_obj@peptide <- peptide_data
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
RunPCA <- function(nPOP_obj){
  sc.batch_cor <- nPOP_obj@protein.imputed

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
  scx <- scx %>% left_join(batch_label,by = c('ID'))


}

# Other small less complex functions

f <- function(y) {c(label=length(y), y=median(y))}


cv<-function(x){

  sd(x, na.rm=T) / mean(x, na.rm=T)

}


fast_cv <- function(dt) {
  summary <- dt[!is.na(value), .(sd_value = sd(value), mean_value = mean(value)), by = .(variable, Protein)]
  dt[!is.na(value), cvq := summary[.SD, on = .(variable, Protein), sd_value / mean_value]]
}


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


normalize <- function(evnew,log = F){
  evnew <- as.matrix(evnew)
  evnew[evnew==0] <- NA
  for(i in 1:ncol(evnew)){
    evnew[,i] <- evnew[,i]/median(evnew[,i],na.rm = T)
  }
  for(i in 1:nrow(evnew)){
    evnew[i,] <- evnew[i,]/mean(evnew[i,],na.rm = T)
  }
  if(log == T){
    evnew <- log2(evnew)
  }

  return(evnew)
}






### RNA seq integration functions




### pSCoPE

Gen_Tlist_DIANN <- function(DIA,start_time){
  DIA_filt <- DIA %>%
    dplyr::select(Stripped.Sequence,Modified.Sequence,RT,PEP,Precursor.Charge,Precursor.Mz,Ms1.Area,Genes,Run,Protein.Group)

  DIA_filt$Charge_double <- as.double(DIA_filt$Precursor.Charge)
  DIA_filt$Mass <- ((DIA_filt$Precursor.Mz)*(DIA_filt$Charge_double))-(1.0072826748*((DIA_filt$Charge_double)))
  DIA_filt$MassRound <- round(DIA_filt$Mass,7)

  DIA_filt$SeqCharge <- paste0(DIA_filt$Stripped.Sequence, DIA_filt$Precursor.Charge)
  DIA_filt_single <- DIA_filt %>% dplyr::group_by(SeqCharge) %>% dplyr::top_n(-1,PEP) %>% dplyr::slice(1)

  DIA_filt_single$RT_sec <- round(DIA_filt_single$RT*60)
  DIA_filt_single$ind <- rownames(DIA_filt_single)
  DIA_filt_single <- DIA_filt_single %>% ungroup()

  filtDIA <- DIA_filt_single %>% group_by(SeqCharge) %>% slice(1)
  #unfiltDIA <- DIA_unfilt_single %>% group_by(SeqCharge) %>% slice(1)

  filtDIA$TargBool <- TRUE
  #unfiltDIA$TargBool <- FALSE

  DIA_input <- filtDIA

  DIA_input$Prot <- gsub(";.*","",DIA_input$Protein.Group)
  DIA_input$Prot <- gsub("-.*","",DIA_input$Prot)


  DIA_input$RTBool <- TRUE
  DIA_input$Masses <- "376.27"
  DIA_input$Leading.razor.protein <- DIA_input$Protein.Group
  DIA_input$Leading.razor.protein <- sub(';.*$','', DIA_input$Leading.razor.protein)


  colsToSelect <- c("Stripped.Sequence","Modified.Sequence","RT","Precursor.Charge","Ms1.Area","MassRound","TargBool","RTBool","Masses", "Leading.razor.protein","SeqCharge")

  FinalDF <- DIA_input %>% dplyr::ungroup() %>% dplyr::select(all_of(colsToSelect))

  #FinalDF$newMod <- paste0(FinalDF$PEPStrippedSequence,"0",rownames(FinalDF))

  colnames(FinalDF) <-c("Sequence","ModSeq","Retention.time","Charge","Apex.intensity","Mass","TargBoolean","RTC_Boolean","Fragments.mz","Leading.razor.protein","Modified")


  FinalDF$Retention.time <- round(FinalDF$Retention.time,3)
  FinalDF$Apex.intensity <- round(FinalDF$Apex.intensity,6)

  ###
  # Uncomment this if you're using a delayed acquisition method, otherwise not
  ###
  FinalDF$Retention.time <- FinalDF$Retention.time - start_time




  # Defining intensity tertiles for intensity-based priority tiers
  IntQuant <-quantile(FinalDF$Apex.intensity, c(.33,.66))

  # To define priority based on precursor intensity
  FinalDF$Priority <- ifelse((FinalDF$TargBoolean == TRUE) & (FinalDF$Apex.intensity <= IntQuant[[1]]),1,ifelse((FinalDF$TargBoolean == TRUE) & (FinalDF$Apex.intensity > IntQuant[[1]])&(FinalDF$Apex.intensity <= IntQuant[[2]]),2,ifelse((FinalDF$TargBoolean == TRUE)& (FinalDF$Apex.intensity > IntQuant[[2]]), 3,0)))

  FinalDF$ModSeq <- NULL

  FinalDF_differential <- FinalDF %>% filter(Priority == 3)

  prot <- c()
  numb <- c()

  for(i in unique(FinalDF_differential$Leading.razor.protein)){
    numb <- c(numb,sum(FinalDF_differential$Leading.razor.protein == i))
    prot <- c(prot,i)
  }

  df_lim <- as.data.frame(prot)
  df_lim$numn <- numb

  df_lim <- df_lim %>% filter(numb > 4)

  FinalDF_differential <- FinalDF_differential[order(-FinalDF_differential$Apex.intensity),]



  downgrade <- c()
  for(i in df_lim$prot){
    df_prot_lim <- FinalDF_differential %>% filter(Leading.razor.protein == i)
    downgrade <- c(downgrade,df_prot_lim$Modified[5:nrow(df_prot_lim)])

  }

  FinalDF$Priority[FinalDF$Modified %in% downgrade] <- 1


  return(FinalDF)
}

reformat <- function(FinalDF,scout){
  scout$seqcharge <- paste0(scout$Sequence, scout$Charge)
  scout <- scout %>% filter(PEP < .04)
  scout <- scout %>% filter(PIF > .75)

  FinalDF$Priority[FinalDF$Priority ==3] <- 2
  FinalDF$Priority[FinalDF$Modified %in% scout$seqcharge ] <- 3


  FinalDF_differential <- FinalDF %>% filter(Priority == 3)

  prot <- c()
  numb <- c()

  for(i in unique(FinalDF_differential$Leading.razor.protein)){
    numb <- c(numb,sum(FinalDF_differential$Leading.razor.protein == i))
    prot <- c(prot,i)
  }

  df_lim <- as.data.frame(prot)
  df_lim$numn <- numb

  df_lim <- df_lim %>% filter(numb > 4)

  FinalDF_differential <- FinalDF_differential[order(-FinalDF_differential$Apex.intensity),]



  downgrade <- c()
  for(i in df_lim$prot){
    df_prot_lim <- FinalDF_differential %>% filter(Leading.razor.protein == i)
    downgrade <- c(downgrade,df_prot_lim$Modified[5:nrow(df_prot_lim)])

  }

  FinalDF$Priority[FinalDF$Modified %in% downgrade] <- 1


  return(FinalDF)

}



# if(carrier == F){
#   # Compile  number of proteins per cell with and without channel q value filter
#   #Max LFQ protein level
#   Raw_data_lim <- Raw_data_lim %>% filter(File.Name %in% good_cells$variable)
#   protein_mat <- diann_maxlfq(Raw_data_lim,sample.header = "File.Name",group.header = "Protein.Group",id.header = "seqcharge",quantity.header = "Ms1.Area")
#
#   #Max LFQ protein level
#   Raw_data_lim_NF <- Raw_data_lim_NF %>% filter(File.Name %in% good_cells$variable)
#   protein_mat_NF <- diann_maxlfq(Raw_data_lim_NF,sample.header = "File.Name",group.header = "Protein.Group",id.header = "seqcharge",quantity.header = "Ms1.Area")
#
# }else{# Collapse to median peptide abundance because we normalized by carrier
#   protein_mat <- matrix(data = NA, nrow = length(unique(protein_list)),ncol = nrow(good_cells))
#   rownames(protein_mat) <- 1:nrow(protein_mat)
#   colnames(protein_mat) <- colnames(cells_filt)
#   for(i in 1:length(unique(protein_list))){
#     if(length(which(protein_list==unique(protein_list)[i])) > 1){
#       protein_mat[i,] <- colMedians(as.matrix(cells_filt[which(protein_list==unique(protein_list)[i]),]),na.rm = T)
#     }else{
#       protein_mat[i,] <- cells_filt[which(protein_list==unique(protein_list)[i]),]
#     }
#
#     rownames(protein_mat)[i] <- unique(protein_list)[i]
#   }
#   protein_mat[protein_mat==0] <- NA
#   protein_mat[is.nan(protein_mat)==T] <- NA
#   protein_mat_NF <-  matrix(data = NA, nrow = length(unique(protein_listNF)),ncol = nrow(good_cells))
#   rownames(protein_mat_NF) <- 1:nrow(protein_mat_NF)
#   colnames(protein_mat_NF) <- colnames(cells_NF_filt)
#   for(i in 1:length(unique(protein_listNF))){
#     if(length(which(protein_listNF==unique(protein_listNF)[i])) > 1){
#       protein_mat_NF[i,]<-colMedians(as.matrix(cells_NF_filt[which(protein_listNF==unique(protein_listNF)[i]),]),na.rm = T)
#     }else{
#       protein_mat_NF[i,] <- cells_NF_filt[which(protein_listNF==unique(protein_listNF)[i]),]
#     }
#
#     rownames(protein_mat_NF)[i] <- unique(protein_listNF)[i]
#   }
#   protein_mat_NF[protein_mat_NF==0] <- NA
#   protein_mat_NF[is.nan(protein_mat_NF)==T] <- NA
#
# }

# Raw_cell_mat <- hjj
# cells_filt <- Normalize_reference_vector(cells_filt)
# cells_NF_filt <- Raw_data_lim_NF.d %>% dplyr::select(good_cells$variable)
# cells_NF_filt <- Normalize_reference_vector(cells_NF_filt)
# numb_pep_filt <- colSums(cells_filt > 0,na.rm = T)
# numb_pep_NF <- colSums(cells_NF_filt > 0,na.rm = T)
# numb_pep <- as.data.frame(numb_pep_filt)
# colnames(numb_pep) <- 'number_peptides'
# numb_pep$filt <- paste0(ChQval,' Channel Qval filter')
# numb_pepnf <- as.data.frame(numb_pep_NF)
# colnames(numb_pepnf) <- 'number_peptides'
# numb_pepnf$filt <- 'No Channel Qval filter'
# numb_pep <- rbind(numb_pep,numb_pepnf)

# nPOP_obj <- list(Raw_data,Ref_norm_data)
# names(nPOP_obj) <- c('RAW','RefNorm')
#
# nPOP <- setClass(
#   Class = 'nPOP',
#   slots = c(
#     assays = 'list',
#     meta.data = 'data.frame',
#     active.assay = 'character',
#     active.ident = 'factor',
#     graphs = 'list',
#     neighbors = 'list',
#     reductions = 'list',
#     images = 'list',
#     project.name = 'character',
#     misc = 'list',
#     version = 'package_version',
#     commands = 'list',
#     tools = 'list'
#   )
# )
