

#formatting

my_col3 <- c("purple2","black")

dot_plot <-  theme_bw()+theme(plot.title = element_text(hjust = .5,size = 24),
                              axis.title.x = element_text(size = 20),
                              axis.title.y = element_text(size = 20),
                              axis.text.x = element_text(size = 12),
                              axis.text.y = element_text(size = 12))




# Other small less complex functions


inSet_norm <- function(Raw_data, cellenOne_meta){
  count = 0

  for(i in unique(paste0(cellenONE_meta$injectWell))){

    set <- Raw_data[,which(grepl(i,colnames(Raw_data)))]

    Raw_data[,which(grepl(i,colnames(Raw_data)))] <- set - rowMeans(set,na.rm = T)

  }

  return(Raw_data)
}


f <- function(y) {c(label=length(y), y=median(y))}


cv<-function(x){

  sd(x, na.rm=T) / mean(x, na.rm=T)

}


fast_cv <- function(dt) {
  summary <- dt[!is.na(value), .(sd_value = sd(value), mean_value = mean(value)), by = .(variable, Protein)]
  dt[!is.na(value), cvq := summary[.SD, on = .(variable, Protein), sd_value / mean_value]]
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


