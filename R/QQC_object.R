# QQC Object

QQC <- setClass(
  Class = 'QQC',
  slots = c(
    ms_type = 'character',
    plex = 'numeric',
    raw_data = 'data.frame',

    matricies = 'ANY',

    cellenONE.meta = 'data.frame',
    meta.data = 'data.frame',
    run_order.statistics = 'list',
    pep.cor = 'list',
    neg_ctrl.info = 'data.frame',
    LC_batch_deviations = 'list',
    reductions = 'list',
    misc = 'list',

    miceotopes = 'ANY'


  )
)

matricies_DDA <- setClass(
  Class = 'matricies_DDA',
  slots = c(

    peptide = 'matrix',
    protein = 'matrix',
    protein_abs = 'matrix',
    protein.imputed = 'matrix',
    peptide_protein_map = 'data.frame'


  )
)

matricies_DIA <- setClass(
  Class = 'matricies_DIA',
  slots = c(

    peptide = 'matrix',
    peptide_mask = 'matrix',
    #peptide_filt.MS2 = 'matrix',
    protein = 'matrix',
    protein_mask = 'matrix',
    #protein_filt.MS2 = 'matrix',
    protein_abs = 'matrix',
    protein.imputed = 'matrix',
    peptide_protein_map = 'data.frame'


  )
)

matricies_Miceotopes <- setClass(
  Class = 'matricies_Miceotopes',
  slots = c(

    Raw_H = 'matrix',
    Raw_L = 'matrix',
    HovL_pep = 'matrix',
    Beta_pep = 'matrix',
    Alpha_pep = 'matrix',
    HovL_prot = 'matrix',
    Beta_prot = 'matrix',
    Alpha_prot = 'matrix',

    peptide_protein_map = 'data.frame'


  )
)


#' Import MaxQuant evidence file and create a QQC object for DDA data.
#'
#' Reads a MaxQuant evidence.txt file (or a directory of such files), joins it
#' with a sample-to-well linker file, filters by PEP and PIF thresholds,
#' removes reverse hits, parses protein accessions, and constructs a QQC
#' object ready for downstream processing.
#'
#' @param data_path Path to an evidence.txt file or a directory containing multiple evidence files.
#' @param linker Path to a CSV linker file mapping \code{Run} names to well and plate metadata.
#' @param plex Integer TMT plex size (e.g., \code{14}, \code{29}, or \code{32}).
#' @param PIF_in Numeric minimum precursor ion fraction (PIF) threshold.
#' @param PEP_in Numeric maximum posterior error probability (PEP) threshold.
#' @return A QQC object with \code{raw_data}, \code{meta.data}, \code{ms_type = "DDA"}, and plex information.
#' @export
MQ_to_QQC <- function(data_path,linker,plex ,PIF_in,PEP_in){


  linker <- read.csv(linker)

  if(plex == 14){
    RI_numb = 18
  }
  if(plex == 29){
    RI_numb = 32
  }

  if(plex == 32){
    RI_numb = 32
  }

  columns_to_read <- c('Modified sequence','Intensity','Retention time','Charge','Raw file','PEP','PIF', 'Leading razor protein',
                       'Potential contaminant','Reverse', paste0("Reporter intensity ",1:RI_numb))

  # Read in file

  if(dir.exists(data_path) == F){
    data <- data.table::fread(data_path,select = columns_to_read)
  }

  if(dir.exists(data_path)){

    if(str_sub(data_path,-1) != '/'){
      data_path <- paste0(data_path,'/')
    }

    raw_files <- list.files(data_path)
    for(i in 1:length(raw_files)){
      if(i == 1){
        data <- data.table::fread(paste0(data_path,raw_files[i]),select = columns_to_read)
      }
      if(i > 1){
        data_next <- data.table::fread(paste0(data_path,raw_files[i]),select = columns_to_read)
        data_list <- list(data, data_next)
        data <- data.table::rbindlist(data_list)
      }
    }
  }




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
  #data<-data %>% dplyr::filter(Potential.contaminant != '+')
  data<-data %>% filter(Reverse != '+')



  QQC <- new('QQC',raw_data = data, meta.data = linker ,ms_type = 'DDA', misc = list(plex = plex))


  return(QQC)

}


#' Import DIA-NN output and create a QQC object for DIA data.
#'
#' Reads a DIA-NN report file (TSV or Parquet, single file or directory),
#' joins it with a sample-to-well linker, filters by library protein group
#' Q-value, extracts the mTRAQ plex channel from precursor IDs, and
#' constructs a QQC object for downstream processing.
#'
#' @param data_path Path to a DIA-NN report file or a directory of report files (TSV or Parquet).
#' @param linker_path Path to a CSV linker file mapping \code{Run} names to well and plate metadata.
#' @param plex Integer plex size (e.g., \code{2} or \code{3} for mTRAQ channels).
#' @param carrier Logical; if \code{TRUE}, sets the MS type to \code{"DIA_C"} (carrier-assisted DIA).
#'   Default \code{FALSE}.
#' @return A QQC object with \code{raw_data}, \code{meta.data}, \code{ms_type}, and plex information.
#' @export
DIANN_to_QQC <- function(data_path,linker_path,plex,carrier = F){

  linker <- read.csv(linker_path)
  linker$Order <- 1:nrow(linker)



  columns_to_read <-c('Genes','Run','Lib.PG.Q.Value','RT','Precursor.Id','Stripped.Sequence','Precursor.Mz',
                      'Precursor.Charge','Precursor.Quantity','Ms1.Area','Protein.Group','Translated.Q.Value','Channel.Q.Value')

  if(dir.exists(data_path) == F){
    if(grepl('\\.parquet$', data_path, ignore.case = TRUE)){
      Raw_data <- arrow::read_parquet(data_path)
      Raw_data <- Raw_data %>% dplyr::select(dplyr::all_of(columns_to_read))
    } else {
      Raw_data <- data.table::fread(data_path,select = columns_to_read)
    }
  }

  if(dir.exists(data_path)){

    if(str_sub(data_path,-1) != '/'){
      data_path <- paste0(data_path,'/')
    }

    raw_files <- list.files(data_path)
    for(i in 1:length(raw_files)){
      file_path <- paste0(data_path, raw_files[i])
      if(grepl('\\.parquet$', raw_files[i], ignore.case = TRUE)){
        data_next <- arrow::read_parquet(file_path)
        data_next <- data_next %>% dplyr::select(dplyr::all_of(columns_to_read))
      } else {
        data_next <- data.table::fread(file_path, select = columns_to_read)
      }
      if(i == 1){
        Raw_data <- data_next
      } else {
        data_list <- list(Raw_data, data_next)
        Raw_data <- data.table::rbindlist(data_list)
      }
    }
  }

  Raw_data <- as.data.frame(Raw_data)

  Raw_data <- Raw_data %>% filter(Lib.PG.Q.Value < .01)
  Raw_data <- Raw_data %>% filter(Run %in% linker$Run)
  Raw_data <- Raw_data %>% left_join(linker, by = c('Run'))

  #Unique precursor ID
  Raw_data$seqcharge <- paste0(Raw_data$Stripped.Sequence,Raw_data$Precursor.Charge)
  Raw_data <- Raw_data %>% filter(Protein.Group != '')


  # this grabs the mTRAQ tag used (may need to be adjusted if using different multiplexing tag)
  Raw_data$plex <- substr(Raw_data$Precursor.Id[1:nrow(Raw_data)], 10, 10)

  # Unique cell ID
  Raw_data$ID <- paste0(Raw_data$Well,Raw_data$plate,'.',Raw_data$plex)
  Raw_data$File.Name <- Raw_data$ID

  #Remove redundant data points
  Raw_data$uq <- paste0(Raw_data$File.Name,Raw_data$Protein.Group,Raw_data$seqcharge)
  Raw_data <- Raw_data %>% distinct(uq,.keep_all = T)
  Raw_data$uq <- NULL

  if(carrier == F){
    QQC <- new('QQC',raw_data = Raw_data,ms_type = 'DIA',meta.data = linker ,misc = list(plex = plex))

  }
  if(carrier == T){
    QQC <- new('QQC',raw_data = Raw_data, meta.data = linker, ms_type = 'DIA_C', misc = list(plex = plex))
  }



  return(QQC)

}



#' Import JMOD Parquet output and create a QQC object for DIA data.
#'
#' Reads a JMOD-format Parquet file, renames columns to the standard DIA-NN
#' schema, joins with a sample-to-well linker, and constructs a QQC object.
#' This importer is intended for data processed with the JMOD pipeline as
#' an alternative to DIA-NN.
#'
#' @param data_path Path to a JMOD Parquet file.
#' @param linker_path Path to a CSV linker file mapping \code{Run} names to well and plate metadata.
#' @param plex Integer plex size (e.g., \code{2} or \code{3} for mTRAQ channels).
#' @param carrier Logical; if \code{TRUE}, sets the MS type to \code{"DIA_C"} (carrier-assisted DIA).
#'   Default \code{FALSE}.
#' @return A QQC object with \code{raw_data}, \code{meta.data}, \code{ms_type}, and plex information.
#' @export
JMOD_to_QQC <- function(data_path,linker_path,plex,carrier = F){

  linker <- read.csv(linker_path)
  linker$Order <- 1:nrow(linker)

  # {'Run':"file_name",
  #   'Modified.Sequence':"seq",
  #   'Stripped.Sequence':"stripped_seq",
  #   'Protein.Group':"protein",
  #   'Channel':"channel",
  #   'Channel.Q.Value':None,
  #   'Translated.Q.Value':"Qvalue",
  #   'Ms1.Area':"plex_Area",
  #   'Ms1.Normalised':None,
  #   'Precursor.Quantity':"coeff", # would not recommend using
  #   'Proteotypic':None,
  #   'Precursor.Normalised':None,
  #   'PEP':None
  # }

  columns_to_read <-c('file_name','seq','stripped_seq',
                      'z','plex_Area','protein','channel')

  Raw_data <- arrow::read_parquet(
    data_path
  )

  Raw_data <- Raw_data %>% dplyr::select(dplyr::all_of(columns_to_read))

  # setnames(
  #   Raw_data,
  #   c('Run','RT','Precursor.Id','Stripped.Sequence',
  #     'Precursor.Mz','Precursor.Charge','Precursor.Quantity',
  #     'Ms1.Area','Protein.Group','plex')
  # )



  colnames(Raw_data) <- c('Run','Precursor.Id','Stripped.Sequence',
                          'Precursor.Charge',
                          'Ms1.Area','Protein.Group','plex')

  Raw_data$RT <- 0
  Raw_data$Precursor.Quantity <- 0
  Raw_data$Precursor.Mz <- 0
  Raw_data$Channel.Q.Value <- 0



  #Raw_data <- data.table::fread(data_path,select = columns_to_read)

  Raw_data <- as.data.frame(Raw_data)


  Raw_data <- Raw_data %>% filter(Run %in% linker$Run)
  Raw_data <- Raw_data %>% left_join(linker, by = c('Run'))

  #Unique precursor ID
  Raw_data$seqcharge <- paste0(Raw_data$Stripped.Sequence,Raw_data$Precursor.Charge)
  Raw_data <- Raw_data %>% filter(Protein.Group != '')



  # Unique cell ID
  Raw_data$ID <- paste0(Raw_data$Well,Raw_data$plate,'.',Raw_data$plex)
  Raw_data$File.Name <- Raw_data$ID

  #Remove redundant data points
  #Raw_data$uq <- paste0(Raw_data$File.Name,Raw_data$Protein.Group,Raw_data$seqcharge)
  #Raw_data <- Raw_data %>% distinct(uq,.keep_all = T)
  #Raw_data$uq <- NULL

  if(carrier == F){
    QQC <- new('QQC',raw_data = Raw_data,ms_type = 'DIA',meta.data = linker ,misc = list(plex = plex))

  }
  if(carrier == T){
    QQC <- new('QQC',raw_data = Raw_data, meta.data = linker, ms_type = 'DIA_C', misc = list(plex = plex))
  }



  return(QQC)

}









#' Import FragPipe PSM output and create a QQC object for DDA data.
#'
#' Reads PSM files from FragPipe output directories, renames columns to match
#' the MaxQuant-style schema used internally, maps TMT reporter ion channels,
#' joins with a sample-to-well linker, filters by purity (PIF), and
#' constructs a QQC object ready for downstream processing.
#'
#' @param data_path Path to a FragPipe results directory containing per-experiment subdirectories with \code{psm.tsv} files.
#' @param linker Path to a CSV linker file mapping \code{Run} names to well and plate metadata.
#' @param plex Integer TMT plex size (currently supports \code{32}).
#' @param PIF_in Numeric minimum purity (PIF) threshold.
#' @return A QQC object with \code{raw_data}, \code{meta.data}, \code{ms_type = "DDA"}, and plex information.
#' @export
FragPipe_to_QQC <- function(data_path,linker,plex ,PIF_in){


  folders <- list.dirs(data_path, recursive = FALSE)

  columns_to_use <- c('Peptide','Intensity','Apex.Retention.Time','Charge','Spectrum.File','Purity','Protein.ID')
  psm_list <- lapply(folders, function(folder) {
    psm_file <- file.path(folder, "psm.tsv")
    if (file.exists(psm_file)) {
      read.delim(psm_file)
    } else {
      NULL
    }
  })

  count = 0
  for(i in psm_list){
    if(is.null(i)==F){
      if(count == 0){
        data <- i
        TMT_grab <- colnames(data)[grep('Intensity.',colnames(data))] #<- str_extract(colnames(data)[grep('Intensity.',colnames(data))], "(?<=_)[^_]+$")
        columns_to_use <- c('Peptide','Intensity','Apex.Retention.Time','Charge','Spectrum.File','Purity','Protein.ID',TMT_grab)

        data <- data %>% dplyr::select(all_of(columns_to_use))


        colnames(data)[grep('Intensity.',colnames(data))] <- str_extract(colnames(data)[grep('Intensity.',colnames(data))], "(?<=_)[^_]+$")

      }else{
        data1 <- i

        TMT_grab <- colnames(data1)[grep('Intensity.',colnames(data1))] #<- str_extract(colnames(data)[grep('Intensity.',colnames(data))], "(?<=_)[^_]+$")
        columns_to_use <- c('Peptide','Intensity','Apex.Retention.Time','Charge','Spectrum.File','Purity','Protein.ID',TMT_grab)
        data1 <- data1 %>% dplyr::select(all_of(columns_to_use))
        colnames(data1)[grep('Intensity.',colnames(data1))] <- str_extract(colnames(data1)[grep('Intensity.',colnames(data1))], "(?<=_)[^_]+$")

        data <-rbind(data,data1)
      }
      count <- 1
    }
  }




  parts <- str_split(data$Spectrum.File, "\\\\")
  data$Spectrum.File <- sapply(parts, function(x) x[length(x) - 1])






  linker <- read.csv(linker)




  if(plex == 32){


    map_FP <- c("126", "127N", "127C","128N",
             "128C","129N","129C",
             "130N","130C","131N",
             "131C","132N","132C",
             "133N","133C","134N",
             "134C","135N","127D",
             "128ND","128CD","129ND",
             "129CD","130ND","130CD",
             "131ND","131CD","132ND",
             "132CD","133ND","133CD",
             "134ND","134CD","135ND","135CD")

    map_MQ <-  c("1", "2", "3","5",
                 "7","9","11",
                 "13","15","17",
                 "19","21","23",
                 "25","27","29",
                 "31","33","4",
                 "6","8","10",
                 "12","14","16",
                 "18","20","22",
                 "24","26","28",
                 "30","32","34","35")

    map_MQ <- paste0('Reporter.ion.', map_MQ)

    columns_renamed <- c('Modified.sequence','Intensity','Retention.time','Charge','Raw.file','PIF','Leading.razor.protein',map_MQ)

    colnames(data) <- columns_renamed

  }




  # Filter for only raw files in the linker
  data <- data %>% dplyr::filter(Raw.file %in% linker$Run)

  # Denote run order
  linker$Order <- 1:nrow(linker)

  #Link data to inject wells
  data <- data %>% left_join(linker, by = c('Raw.file'= 'Run'))

  # Unique precursor ID
  data$seqcharge <- paste0(data$Modified.sequence,data$Charge)
  data$seqRun <- paste0(data$seqcharge, data$Raw.file)
  data <- data %>% distinct(seqRun,.keep_all = T)





  # Filter data
  #data<-data %>% filter(PEP < PEP_in)
  data<-data %>% dplyr::filter(PIF > PIF_in)
  #data<-data %>% dplyr::filter(Potential.contaminant != '+')
  #data<-data %>% filter(Reverse != '+')



  QQC <- new('QQC',raw_data = data, meta.data = linker ,ms_type = 'DDA', misc = list(plex = plex))


  return(QQC)

}







