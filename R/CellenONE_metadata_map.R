
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
  meta <- nPOP_obj@meta.data
#  for(i in 1:length(allDays)){
#
#  }

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

  #cellID <- cellID %>% left_join(meta, by = c('InjectWell' = 'Well'))

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




#' test1.
#'
#' This function takes two numeric inputs and returns their sum.
#'
#' @param x A numeric value.
#' @param y A numeric value.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add_numbers(2, 3)
#' @export

PlotSlideLayout_celltype <- function(nPOP_obj){

  ggplot(nPOP_obj@cellenONE.meta) +
    geom_point(aes(x = dropXPos,y = dropYPos,color = sample)) +
    geom_text(aes(x = pickupXPos_numb,y = pickupYPos_numb,label = injectWell,size = 5),hjust= .5, vjust=-.6) +
    facet_wrap(~field,ncol = 4)+
    scale_y_reverse()

}


#' test1.
#'
#' This function takes two numeric inputs and returns their sum.
#'
#' @param x A numeric value.
#' @param y A numeric value.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add_numbers(2, 3)
#' @export
PlotSlideLayout_label <- function(nPOP_obj){

  # print the mTRAQ labels overlayed on the positions of the slide
  ggplot(nPOP_obj@cellenONE.meta, aes(x = dropXPos,y = dropYPos,color = label)) +
    geom_point() +scale_y_reverse()+ facet_wrap(~field,ncol = 4)

}
