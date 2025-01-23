peaks<-peaks.raw
setwd('C:/Rprogram/SGC/ASMS_drug/positive')
msfiles<-list.files()

##keep highly abundant ones
index<-NULL
for (i in 1:nrow(peaks)){
  intensity<-mean(as.numeric(peaks.raw[i,3:5]))
  if (intensity>1000000){
    index<-c(index,i)
  }
}

peaks<-peaks[index,]
mydata<-matrix(rep(0,length(index)*102),nrow=length(index),ncol=102)
xrawdata<-xcmsRaw(msfiles[1])
ppmwin<-5*10^(-6)
rtwin<-60
#get the scans
for (j in 1:nrow(peaks)){
  RT<-peaks$rt[j]
  
  if (RT<1||RT>7){##skip the peaks with weird retention time
    next
  }
  
  mz<-peaks$mz[j]
  
  mydata[j,1]<-mz
  mydata[j,2]<-RT
  
  #extract intensity
  rt<-RT*60
  
  #'smooth the exact mass and retention time
  mz.min<-mz*(1-ppmwin)
  mz.min<-max(mz.min, xrawdata@mzrange[1])
  mz.max<-mz*(1+ppmwin)
  mz.max<-min(mz.max,xrawdata@mzrange[2])
  
  #'range of analysis is the entire chromatograph
  rt.min<-max(min(xrawdata@scantime),rt-rtwin)
  
  #' the first and last 3 minutes are skipped to avoid sections with no analytes
  rt.max<-min(max(xrawdata@scantime),rt+rtwin)
  peaks.scan<-rawEIC(xrawdata,mzrange=cbind(mz.min,mz.max),rtrange=cbind(rt.min, rt.max))
  
  if (length(peaks.scan$intensity)>100){##cut if too long
    peaks.scan$intensity<-peaks.scan$intensity[1:100]
  }
  mydata[j,3:(length(peaks.scan$intensity)+2)]<-peaks.scan$intensity
  }