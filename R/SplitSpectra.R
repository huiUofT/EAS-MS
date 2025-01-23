#' function to splitting data to negative and positve spectra
#'
#' @param input
#' @export
SplitSpectra<-function(input){
  #' the path to save results
  path<-getwd()

  #' the path to save raw data
  path.data<-paste0(path,"/data")
  setwd(path.data)
  msfiles<-list.files()


  path_results_neg<-paste0(path,"/negative")
  path_results_pos<-paste0(path,"/positive")
  setwd(path.data)
  msfiles<-list.files()
  msfiles<-mixedsort(msfiles)###sort the file name according to number

  ##########split files to negative and positive results according to signals######
  for (i in 1:length(msfiles)){
    setwd(path.data)
    xraw<-xcmsRaw(msfiles[i])
    group<-factor(rep(c(1,2),floor(0.5*length(xraw@scanindex))))
    if (round(0.5*length(xraw@scanindex))<0.5*length(xraw@scanindex)){group<-factor(c(group,1))}
    xraw<-split(xraw,f=group)
    filename<-strsplit(msfiles[i],"[.]")
    filename<-unlist(filename[[1]])[1]

    if (input=='pos'){##if the first scan is positve
    setwd(path_results_pos)
    write.mzdata(xraw$'1',filename)
    setwd(path_results_neg)
    write.mzdata(xraw$'2',filename)}
    if (input=='neg'){##if the first scan is negative
      setwd(path_results_pos)
      write.mzdata(xraw$'2',filename)
      setwd(path_results_neg)
      write.mzdata(xraw$'1',filename)}
  }
}
