#' This function is used for initial searching
#'
#' @param mydata
#' @param Control
#'
#' @return
#' @export
#'
#' @examples

Diffpeak_STD<-function(peaks.raw,STD,Control,Fold){
    header <- names(peaks.raw)
  ControlID<-grep(Control,header)
  STDID<-grep(STD,header)
  
  indexsave<-NULL
  for ( i in 1:nrow(peaks.raw)){
    test<-peaks.raw[i,STDID]#data of standard
    ctrl<-peaks.raw[i,ControlID]#data of control
    fold<-(sum(test)*length(ctrl))/(sum(ctrl)*length(test))
    if (fold<Fold){next}
    
    if (sd(test) > 0){
      ttest <- t.test(test, ctrl)
      pvalue <- ttest$p.value
      if (pvalue<0.05){
        indexsave<-c(indexsave,i)
      }
    }}
  
  
  setwd(path)
  return(peaks.raw[indexsave,])
}


