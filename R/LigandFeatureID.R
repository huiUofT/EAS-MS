#' This function is used for initial searching
#'
#' @param mydata
#' @param Control
#'
#' @return
#' @export
#'
#' @examples

LigandFeatureID<-function(mydata,polarity,control,Fold){


  header <- names(mydata)
  mydata$ID<-1:nrow(mydata)

  header_list <- header[3:length(header)]

  Anno <- data.frame(
    File = header_list,
    Protein = sapply(strsplit(header_list, "_"), function(x) x[1]),
    Sample = sapply(strsplit(header_list, "_"), function(x) x[2]),
    Replicate = sapply(strsplit(header_list, "_"), function(x) x[3])
  )

  Anno$Replicate <- as.numeric(Anno$Replicate)


  ###### 5.Determine numbers of protein target

  prname <- unique(Anno$Protein)

  prname <- prname[prname != STD]

  prnum <- length(prname)

  ###### 6.Statistics

  mydata <- as_tibble(mydata)

  datafinal <- mydata[, 0]
  peaks<-NULL

  for (i in 1:prnum) {

    ###### Select protein target

    prni <- prname[i]

    Annopri <- subset(Anno, Protein == prni)

    sampname <- unique(Annopri$Sample)

    sampname <- sampname[sampname != control]

    sampnum <- length(sampname)
    if (sampnum==0){next}

    for (j in 1: sampnum){

      ###### Select environmental sample type

      data <- mydata

      new_column_fc1 <- rep(1, nrow(data))

      new_column_p1 <- rep(1, nrow(data))

      data <- cbind(data,
                    new_column_fc1,
                    new_column_p1)

      sampnj <- sampname[j]

      test_anno <- Anno[which(Anno$Protein == prni), ]
      control1 <- Anno[which(Anno$Sample == sampnj), ]##the chemical pool of other proteins, as controls
      index1<-which(control1$Protein!=STD)##exclude standards
      index2<-which(control1$Protein[index1]!=prni)##exclude itself
      index1<-index1[index2]
      control1<-control1[index1,]##the other proteins
      test_anno <- test_anno[which(test_anno$Sample == sampnj), ]
      control_solvent<-Anno[which(Anno$Sample==control),]##the other control of methanol

      FileID_test <- test_anno$File

      FileID_control1 <- control1$File
      FileID_control_solvent<-control_solvent$File

      testdata <- data[, colnames(data)%in%FileID_test]

      data.crtl1 <- data[, colnames(data)%in%FileID_control1]
      data.solvent <- data[, colnames(data)%in%FileID_control_solvent]

      ###### Calculate fc and pvalue

      for (g in 1:nrow(data)){
        
        
        #------------------------------------------
        #skip peak features detected in methanol
        #-------------------------------------------
        test <- testdata[g,]
        solvent<-data.solvent[g,]
        fold.solvent<-(sum(test)*length(solvent))/(sum(solvent)*length(test))
        if (fold.solvent<5){next}##skip the chemical if it is detected in methanol at 1/5 times
        
        

        test <- testdata[g,]
        ctrl1 <- data.crtl1[g,]

        fold1 <- (sum(test)*length(ctrl1))/(sum(ctrl1)*length(test))

        data$new_column_fc1[g] <- fold1

        if (sd(test) > 0){
          ttest1 <- t.test(test, ctrl1)
          pvalue1 <- ttest1$p.value
          data$new_column_p1[g] <- pvalue1
        } else {data$new_column_p1[g] <- 1}
      }


      ###### Find significant peaks

      sig <- data %>% dplyr::filter(new_column_fc1 > Fold & new_column_p1 < 0.05)
      sig$ChemicalID<-rep(sampname[j],nrow(sig))
      data$ChemicalID<-rep(sampname[j],nrow(data))
      sig$ProtID<-rep(paste0(prname[i],'_',sampname[j]),nrow(sig))
      data$ProtID<-rep(paste0(prname[i],'_',sampname[j]),nrow(data))

      if (nrow(sig)==0){next}#no significant peak

      ###### Calculate average intensity for each sig peak

      sig$avg <- rep(1, nrow(sig))

      data_avg <- sig[, colnames(sig)%in%FileID_test]

      for (k in 1:nrow(sig)){

        print(prni)

        print(sampnj)

        print(k)

        sam <- data_avg[k,]

        avgsam <- sum(sam)/ncol(sam)

        sig$avg[k] <- avgsam

      }

      ###### Prepare data for volcano plot

      dataset <- data

      dataset$avg <- rep(100, nrow(dataset))

      dataset$change <- rep("stable", nrow(dataset))

      sig$change <- rep("up", nrow(sig))


      ###### Remove sig peaks since the avg has been calculated somewhere else

      dataset <- dataset[-which(dataset$new_column_p1 < 0.05 &
                                  dataset$new_column_fc1 > Fold),]

      dataset <- bind_rows(dataset, sig)

      ###### Plotting

      cut_off_pvalue <- -log10(0.05)

      cut_off_logFC <- log10(Fold)

      temp_p <- arrange(sig, desc(avg))

      tempn <- paste(prni, sampnj, polarity, sep = "_")

      p <- ggplot(
        dataset, aes(x = log10(new_column_fc1), y = -log10(new_column_p1), colour = change, size = avg)) +
        scale_color_manual(values=c("#231f20", "#be1e2d")) +
        geom_point(alpha = 0.7, shape = 16) +
        geom_hline(yintercept = cut_off_pvalue, col="#be1e2d", linetype = 2) +
        geom_vline(xintercept = cut_off_logFC, col="#be1e2d",linetype = 2) +
        theme_bw() +
        theme(panel.grid.major=element_line(colour=NA)) +
        theme(panel.grid.minor=element_line(colour=NA)) +
        geom_text_repel(
          data=temp_p[1:5,],
          aes(label=round(mz,digits=4)), #use mzmed for xcms peak picking
          size=4,
          color="#be1e2d",
          segment.color="black",
          show.legend = FALSE) +
        theme(legend.position = "none") +
        xlab("Log(Fold change)") +
        ylab("-Log P value") +
        theme(axis.text.x = element_text(size = rel(1.5), color = "black")) +
        theme(axis.text.y = element_text(size = rel(1.5), color = "black")) +
        theme(axis.title.x = element_text(size = rel(1.2), color = "black")) +
        theme(axis.title.y = element_text(size = rel(1.2), color = "black")) +
        ggtitle(tempn)

      p

      tempn <- paste(prni, sampnj, polarity, sep = "_")

      tempn <- paste(tempn, ".tiff", sep = "")

      ggsave(p, filename = tempn, width = 6, height = 5)

      ###### Save data

      postfix_fc1 <- "fc1"

      postfix_p1 <- "pvalue1"

      colname1 <- paste(prni, sampnj, postfix_fc1)

      colname2 <- paste(prni, sampnj, postfix_p1)

      colnames(sig)[colnames(sig) == "new_column_fc1"] <- colname1

      colnames(sig)[colnames(sig) == "new_column_p1"] <- colname2

      tempn <- paste(prni, sampnj, polarity, sep = "_")

      tempn <- paste(tempn, ".csv", sep = "")

      write.table(sig, file = tempn, sep=',', row.names = FALSE)

      if (length(peaks)==0){
        peaks<-data[sig$ID,]
      }else{
      peaks<-rbind(peaks,data[sig$ID,])}

      datafinal <- cbind(datafinal, data)

      datafinal <- datafinal[, !duplicated(colnames(datafinal))]

      sig <- NULL

    }

  }

  tempn <- paste(prname, collapse = "_")

  tempn <- paste(tempn, polarity, sep = "_")

  tempn <- paste(tempn, ".csv", sep = "")
  return(peaks)
}


