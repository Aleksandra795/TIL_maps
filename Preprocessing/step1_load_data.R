rm(list=ls())
gc

library(ggplot2)

setwd("./TIL_maps")

# data path (path to folders downloaded from https://stonybrookmedicine.app.box.com/v/til-results-new-model)
dpath <- ""

# list all cancers
cancer <- list.dirs(dpath, full.names=F, recursive=F)

# set colors
col_bin <- c("#377eb8","#e41a1c")
col_prob <- scale_fill_gradientn(colours=c("#377eb8","#081D58","#49006A","#e41a1c"), na.value="white")

for(a in 1:length(cancer)){
  cat("Cancer: ",cancer[a]," \n")
  
  # create folders for saving
  dir.create(paste0("res/plots/img_binary/",cancer[a]), showWarnings=F)
  dir.create(paste0("res/plots/img_prob/",cancer[a]), showWarnings=F)
  dir.create(paste0("res/txt_files/",cancer[a]), showWarnings=F)
  
  # list all files and remove low_resolution data
  d <- gsub("prediction-","",list.files(path=paste0(dpath,cancer[a],"/heatmap_txt_binary"), pattern="prediction"))
  ind <- grepl("low_res",d)
  d <- d[!ind]
  
  # create data with samples
  samples <- data.frame(d)
  colnames(samples) <- "FileName"
  tmp <- sapply(samples$FileName,function(x){strsplit(x,"[.]")[[1]][1]})
  samples$Barcode <- tmp
  samples$Patient <- substring(tmp, 1,12)
  samples$Type <- substring(tmp, 14,15)
  samples$Cancer <- cancer[a]
  samples$TILs_perc <- samples$TILs <- samples$Tissue <- samples$Background <- samples$N <- samples$y <- samples$x <-  NA

  cat("No. of prediction files: ",length(d)," \n")
  
  if (!file.exists(paste0("data/Data_",cancer[a],".RData"))){
    data_all <- list()
    for (b in 1:length(d)){
      
      # check if data exist
      if(file.exists(paste0(dpath,cancer[a],"/heatmap_txt_binary/prediction-", samples$FileName[b])) &
         file.exists(paste0(dpath,cancer[a],"/heatmap_txt_binary/color-",samples$FileName[b])) &
         file.exists(paste0(dpath,cancer[a],"/heatmap_txt/prediction-", samples$FileName[b]))){
        
        # load data
        data1 <- read.delim(file=paste0(dpath,cancer[a],"/heatmap_txt_binary/prediction-", samples$FileName[b]), sep=" ", header=F)
        data2 <- read.delim(file=paste0(dpath,cancer[a],"/heatmap_txt_binary/color-",samples$FileName[b]), sep=" ", header=F)
        data3 <- read.delim(file=paste0(dpath,cancer[a],"/heatmap_txt/prediction-", samples$FileName[b]), sep=" ", header=F)
        # data4 <- read.delim(file=paste0(dpath,cancer[a],"/heatmap_txt/color-",samples$FileName[b]), sep=" ", header=F)
        
        n <- nrow(data1)
        if (nrow(data2) > n) data2 <- data2[1:n,]
        if (nrow(data3) > n) data3 <- data3[1:n,]
        
        # make x,y coordinates
        if ((sum(data1[,1] == data2[,1]) == nrow(data1)) & (sum(data1[,1] == data3[,1]) == nrow(data1)) &
            (sum(data1[,2] == data2[,2]) == nrow(data1)) & (sum(data1[,2] == data3[,2]) == nrow(data1))){
          # cat("Matching OK\n")
        } else{
          cat("Data not matched !!!\n")
        }
        
        # create data frame with all data
        data <- data.frame(unclass(factor(data1[,1])),unclass(factor(data1[,2])))
        colnames(data) <- c("x","y")
        data$Background <- data2[,3] < 12
        data$TILs <- data1[,3]
        data$TILs_prob <- data3[,3]
        
        # remove background before saving
        data_save <- data
        data_save[data$Background,"TILs"] <- NA
        data_save[data$Background,"TILs_prob"] <- NA
        p <- ggplot(data_save, aes(x,y)) + geom_tile(aes(fill=TILs), show.legend = F) +
          theme_void() + theme(panel.background = element_rect(fill="white", colour="white")) +
          scale_fill_gradientn(colours=c("#377eb8","#e41a1c"), values=c(0,1), na.value="white") +
          scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand= c(0,0))
        png(file=paste0("res/plots/img_binary/",cancer[a],"/",samples$Barcode[b],".png"),width=length(unique(data_save$x)),height=length(unique(data_save$y)))
        print(p)
        dev.off()
        
        p <- ggplot(data_save, aes(x,y)) + geom_tile(aes(fill=TILs_prob), show.legend = F) +
          theme_void() + theme(panel.background = element_rect(fill="white", colour="white")) +
          col_prob +
          scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand= c(0,0))
        png(file=paste0("res/plots/img_prob/",cancer[a],"/",samples$Barcode[b],".png"),width=length(unique(data_save$x)),height=length(unique(data_save$y)))
        print(p)
        dev.off()
        
        data_all[[samples$Barcode[b]]] <- data
        write.table(data, file=paste0("res/txt_files/",cancer[a],"/",samples$Barcode[b],".txt"), sep="\t", row.names=F, quote=F)
        
        # save some info about sample
        samples$x[b] <- length(unique(data$x))
        samples$y[b] <- length(unique(data$y))
        samples$N[b] <- nrow(data)
        samples$Background[b] <- sum(data$Background)
        samples$Tissue[b] <- sum(data$TILs==0 & !data$Background)
        samples$TILs[b] <- sum(data$TILs==1 & !data$Background)
        samples$TILs_perc[b] <- 100 * samples$TILs[b]/(samples$Tissue[b] + samples$TILs[b])   
      }
      
    }
    # remove not complete cases
    samples <- samples[!is.na(samples$N),]
    save(data_all, samples, file=paste0("data/Data_",cancer[a],".RData"))

    cat("No. of complete cases: ",nrow(samples)," \n")
  } 
  
}

# get list of samples
samples_all <- list()
cancer_stats_all <- data.frame(matrix(nrow=length(cancer),ncol=4))
colnames(cancer_stats_all) <- c("Cancer","N.predictions","N.Barcodes","N.Patients")
for(a in 1:length(cancer)){
  cat("Cancer: ",cancer[a]," \n")
  cancer_stats_all$Cancer[a] <- cancer[a]
  
  d <- gsub("prediction-","",list.files(path=paste0(dpath,cancer[a],"/heatmap_txt_binary"), pattern="prediction"))
  ind <- grepl("low_res",d)
  d <- d[!ind]
  tmp <- sapply(d,function(x){strsplit(x,"[.]")[[1]][1]})
  cancer_stats_all$N.predictions[a] <- length(tmp)
  
  load(file=paste0("data/Data_",cancer[a],".RData"))
  samples_all[[cancer[a]]] <- samples[,1:5]
  cancer_stats_all$N.Barcodes[a] <- nrow(samples)
  cancer_stats_all$N.Patients[a] <- length(unique(samples$Patient))
}
samples_all <- do.call(rbind,samples_all)

write.table(samples_all, file="res/Samples_summary.txt", sep="\t", row.names=F, quote=F)
write.table(cancer_stats_all, file="res/Cancer_summary.txt", sep="\t", row.names=F, quote=F)