#' This function is used to for ASMS processing according to ion mode and chemical pools
#'
#' @param Pool
#' @param Control
#' @param STD
#' @param Ionmode
#'
#' @return
#' @export
#'
#' @examples

ASMSrun<-function(Ionmode,STD,Control){
  
  #------------------------
  #find out the chemical pools
  #------------------------
  Chempool<-Poolnum(Ionmode,STD,Control)# the chemical pool name
  
  #-****************************************************************
  #-****************************************************************
  #data processing according to each pool
  #-****************************************************************
  #-****************************************************************
  for (i in 1:length(Chempool)){##extracting data according to each pool
  Pool<-Chempool[i]
  
  #----------------------------------------------------
  #extracting peaks from standards, intensity cutoff is 10^6, ppm = 2.5
  #----------------------------------------------------
  peaks.raw<-Peakextract(5*10^5,2.5,Ionmode,Pool,Control,STD)##extracting data of standards and methanol from the same pool
  
  #-----------------------------------------------------
  #'Extracting differential peaks between standard and control
  #-----------------------------------------------------
  peaks.STD<-Diffpeak_STD(peaks.raw,STD,Control,10)#This is the function extracting differential peaks compared to controL
 
  #------------------------------------------------------------
  #matching peaks across samples from the same pool of standards
  #------------------------------------------------------------
  peaks.sample<-Matchpeak(peaks.STD,5,30,Pool,path,Ionmode)
 
  
  #---------------------------------------------------------
  #'Extracting differential peaks across proteins
  #'data has to be organized like NR_XX_X, WT_XX_X, NR_DMSO_X
  #----------------------------------------------------------
  peaks.Diff<-LigandFeatureID(peaks.sample,Ionmode,Control,Fold)#This is the function extracting differential peaks compared to controL
  write.table(peaks.Diff,file=paste0('AllDiff_',Pool,'_',Ionmode,'.csv'),sep=',',row.names = FALSE)
  
  
  #-----------------------------
  #removing isotopic peaks
  #-----------------------------
  #2 ppm is the mass tolerance to find isotopic peaks
  #'0.80 is the correlation coefficient to extract isotopic peaks
  #'3 is the intensity ratio cutoff between primary and other isotopic peaks
  peaks.iso<-FindIsotope(peaks.Diff,2,0.80,2,Ionmode)
  
  
  #----------------------------
  #searching database
  #----------------------------
  setwd(path)
  peaks.iso$mz<-peaks.iso$mz*(1-ppmshift*10^(-6))##simple mass calibration
  if (Ionmode=='positive'){
  Library<-InitialSearch(peaks.iso,ppm,1,Database)}
  if (Ionmode=='negative'){
    Library<-InitialSearch(peaks.iso,ppm,-1,Database)}
  
  #delete redundant peaks
  mz<-NULL
  rt<-NULL
  indexsave<-NULL
  for (i in 1:nrow(Library)){
    index<-which(mz==Library$mz[i])
    if (length(index)==0){#new mz
      mz<-c(mz,Library$mz[i])
      rt<-c(rt,Library$rt[i])
      next
    }
    if (Library$rt[i]==rt[index[1]]){#retention time is identical
      indexsave<-c(indexsave,i)
      next
    }
    mz<-c(mz,Library$mz[i])
    rt<-c(rt,Library$rt[i])
  }
  if (length(indexsave)>0){
    Library<-Library[-indexsave,]}
  
  #-------------------------------------
  #outputs of the ASMS data
  #------------------------------------
  if (length(Library)==0){
    print('No hit found lOl')
  }
  if (nrow(Library)>0){
  write.table(Library,file=paste0('SGC_',Pool,'_',Ionmode,'.csv'),sep=',',row.names = FALSE)
  }
  
  #------------------------------------
  #searching all STD data, not only the difference one
  #------------------------------------
  peaks.iso<-FindIsotope(peaks.sample,2,0.80,2,Ionmode)
  setwd(path)
  peaks.iso$mz<-peaks.iso$mz*(1-ppmshift*10^(-6))##simple mass calibration
  peaks.iso$ChemicalID<-rep(Pool,nrow(peaks.iso))
  if (Ionmode=='positive'){
    Library<-InitialSearch(peaks.iso,ppm,1,Database)}
  if (Ionmode=='negative'){
    Library<-InitialSearch(peaks.iso,ppm,-1,Database)}
  
  #delete redundant peaks
  mz<-NULL
  rt<-NULL
  indexsave<-NULL
  for (i in 1:nrow(Library)){
    index<-which(mz==Library$mz[i])
    if (length(index)==0){#new mz
      mz<-c(mz,Library$mz[i])
      rt<-c(rt,Library$rt[i])
      next
    }
    if (Library$rt[i]==rt[index[1]]){#retention time is identical
      indexsave<-c(indexsave,i)
      next
    }
    mz<-c(mz,Library$mz[i])
    rt<-c(rt,Library$rt[i])
  }
  if (length(indexsave)>0){
    Library<-Library[-indexsave,]}
  
  #-------------------------------------
  #outputs of the ASMS data
  #------------------------------------
  if (nrow(Library)>0){
    write.table(Library,file=paste0('SGCall_',Pool,'_',Ionmode,'.csv'),sep=',',row.names = FALSE)
  }
  }
}


