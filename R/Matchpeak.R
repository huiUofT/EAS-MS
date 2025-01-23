#' Match peaks from the library
#'
#' @param mylib
#' @param ppm
#' @param RT
#' @param btw
#' @param Pool
#' @param path
#' @param folder
#'
#' @return a complete peak feature matrix
#'
Matchpeak<-function(mylib,ppm,btw,Pool,path,folder){
  ppm<-ppm/10^6
  
  #' the path to save raw data
  path.data<-paste0(path,"/",folder)
  
  #' ----------------------------------------------
  #' detect peaks from mass spec raw files
  #' ---------------------------------------------
  setwd(path.data)
  msfiles<-list.files()
  index1<-grep(Control,msfiles)##the msfiles from controls
  
  index2<-grep(Pool,msfiles)##the msfiles from the same chemical pools
  msfiles<-msfiles[c(index1,index2)]#the data of STD and controls from the same pool
  
  #---------------------------------------
  #the MS and RT range of MS documents
  #---------------------------------------
  xset<-xcmsSet(msfiles[1],method='centWave',ppm=2.5,peakwidth=c(10,30),snthresh=10,nSlaves=1)
  minmz<-min(xset@peaks[,1])
  maxmz<-max(xset@peaks[,1])
  minrt<-min(xset@peaks[,4])
  maxrt<-max(xset@peaks[,4])
  
  
  #-------------------------------------
  #match peaks
  #-------------------------------------
  newpeak<-matrix(rep(0,nrow(mylib)*(2+length(msfiles))),nrow=nrow(mylib),ncol=2+length(msfiles))
  newpeak[,1]<-mylib[,1]
  newpeak[,2]<-mylib[,2]
  for (k in 1:length(msfiles)){
      print(c('matching peaks across...',Pool))
      xraw<-xcmsRaw(msfiles[k])
      for (n in 1:nrow(newpeak)){
        mz.value<-mylib[n,1]
        mzmin<-max(minmz,mz.value-mz.value*ppm)
        mzmax<-min(maxmz,mz.value+mz.value*ppm)
        if (mzmin>mz.value||mzmax<mz.value){next}
        rt.value<-mylib[n,2]*60
        rtmin<-max(minrt,rt.value-btw)
        rtmax<-min(maxrt,rt.value+btw)
        if (rtmin>rt.value||rtmax<rt.value){next}
        
        #' extract peak intensities
        peak<-rawEIC(xraw,mzrange=cbind(mzmin,mzmax),rtrange=cbind(rtmin,rtmax))
        intensity<-max(peak$intensity)
        newpeak[n,k+2]<-intensity
      }
  }
  
  colnames(newpeak)<-c('mz','rt',msfiles)
  newpeak[which(newpeak==0)]<-100
  
  setwd(path)
  return(data.frame(newpeak))
}