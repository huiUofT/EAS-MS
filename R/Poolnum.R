#' This function is used to identify the # of pools and samples
#'
#' @param folder
#' @param Control
#' @param STD
#'
#' @return
#' @export
#'
#' @examples

Poolnum<-function(folder,STD,Control){
  
  #' the path to save results
  path<-getwd()
  
  #' the path to save raw data
  path.data<-paste0(path,"/",folder)
  
  #' ----------------------------------------------
  #' detect peaks from mass spec raw files
  #' ---------------------------------------------
  setwd(path.data)
  msfiles<-list.files()
  
  Anno <- data.frame(
    File = msfiles,
    Protein = sapply(strsplit(msfiles, "_"), function(x) x[1]),
    Sample = sapply(strsplit(msfiles, "_"), function(x) x[2]),
    Replicate = sapply(strsplit(msfiles, "_"), function(x) x[3])
  )
  Sample<-Anno$Sample[Anno$Sample!=Control]
  setwd(path)
  return(unique(Sample))
}


