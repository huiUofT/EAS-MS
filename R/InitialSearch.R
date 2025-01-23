#' This function is used for initial searching
#'
#' @param Library
#' @param ppm
#' @param polarity
#' @param Database
#'
#' @return
#' @export
#'
#' @examples
InitialSearch<-function(Library,ppm,polarity,Database){
  if (polarity==1){
    #Adducts<-c('[M+H]+','[M+H-H2O]+','[M+NH4]+','[M+CH3OH+H]+')
    #MW.adducts<-c(1.007825,-17.00274,18.03437,33.03404)-0.00054
    Adducts<-c('[M+H]+')
    MW.adducts<-c(1.007825)-0.00054
  }

  if (polarity==-1){
    #Adducts<-c('[M]-','[M-H]-','[M-Br+O]-','[M-H-H2O]-','[M+Cl]-','[M+CH2O2-H]-')
    #MW.adducts<-c(0,-1.007825,-62.92342,-19.01839,34.96885,44.99765)+0.00054
    Adducts<-c('[M-H]-')
    MW.adducts<-c(-1.007825)+0.00054
  }

  indexresults<-NULL
  Library$ID<-rep(0,nrow(Library))
  Library$SMILES<-rep(0,nrow(Library))
  Library$Adducts<-rep(0,nrow(Library))
  Library$LOGP<-rep(0,nrow(Library))
  Library$formula<-rep(0,nrow(Library))
  Library$mserror<-rep(0,nrow(Library))
  Library$`SGC ID for Component`<-rep(0,nrow(Library))
  Library$`SGC ID for Pool`<-rep(0,nrow(Library))
  for (i in 1:nrow(Library)){
    mz<-Library$mz[i]
    index_pool<-which(Database$`Pooled Well`==Library$ChemicalID[i])##find the chemical from the same pool
    if (length(index_pool)==0){next}
    for (j in 1:length(MW.adducts)){
      mserror<-(mz-Database$MONOISOTOPIC_MASS[index_pool]-MW.adducts[j])/mz
      index1<-which(abs(mserror)<(ppm*10^(-6)))
      if (length(index1)==0){next}
      index1<-index_pool[index1]


      Hbondsave<-0
      if (length(index1)>0){#save the id with match

        ##check Hbond acceptor and donor
        cpd<-Database$SMILES[index1]
        for (kk in 1:length(cpd)){
          mols<-parse.smiles(cpd[kk])
          dc <- get.desc.categories()
          dn <- get.desc.names(dc[4])
          decs<-eval.desc(mols,dn)
          HBDonor<-as.numeric(decs[3])
          HBAccept<-as.numeric(decs[4])

          if (polarity==1&&HBAccept>0){Hbondsave<-c(Hbondsave,kk)}##positive ion mode with >1 hbond acceptor
          if (polarity==-1&&HBDonor>0){Hbondsave<-c(Hbondsave,kk)}##negative ion mode with >1 hbond donor
        }
        if (length(Hbondsave)==0){next}
        index1<-index1[Hbondsave]##only keep the cpds with H bonds



        indexresults<-c(indexresults,i)
        if (Library$ID[i]==0){
          Library$ID[i]<-paste(index1,collapse = ';')
        }else{
          Library$ID[i]<-paste(c(Library$ID[i],index1),collapse = ';')
        }
        if (Library$SMILES[i]==0){
          Library$SMILES[i]<-paste(c(Database$SMILES[index1]),collapse = ';')
          Library$`SGC ID for Component`[i]<-paste(c(Database$`SGC ID for Component`[index1]),collapse = ';')
          Library$`SGC ID for Pool`[i]<-paste(c(Database$`SGC ID for Pool`[index1]),collapse = ';')
        }else{
          Library$SMILES[i]<-paste(c(Library$SMILES[i],Database$SMILES[index1]),collapse = ';')
          Library$`SGC ID for Component`[i]<-paste(c(Library$`SGC ID for Component`[i],Database$`SGC ID for Component`[index1]),collapse = ';')
          Library$`SGC ID for Pool`[i]<-paste(c(Library$`SGC ID for Pool`[i],Database$`SGC ID for Pool`[index1]),collapse = ';')
        }
        Adduct.paste<-rep(Adducts[j],length(index1))
        if (Library$Adducts[i]==0){
          Library$Adducts[i]<-paste(Adduct.paste,collapse = ';')
        }else{
          Library$Adducts[i]<-paste(c(Library$Adducts[i],Adduct.paste),collapse = ';')
        }
        if (Library$formula[i]==0){
          Library$formula[i]<-paste(Database$formula[index1],collapse = ';')
          Library$mserror[i]<-paste(mserror[index1]*10^6,collapse = ';')
        }else{
          Library$formula[i]<-paste(c(Library$formula[i],Database$formula[index1]),collapse = ';')
          Library$mserror[i]<-paste(c(Library$mserror[i],mserror[index1]*10^6),collapse = ';')
        }
        for (k in 1:length(index1)){
          Library$LOGP[i]<-paste(c(Library$LOGP[i],Database$LOGP[index1[k]]),collapse = ';')
        }
      }
    }##findout matching, ppm
  }
  Library<-Library[indexresults,]
  Library<-Library%>%dplyr::filter(formula!='')
  return(Library)
}
