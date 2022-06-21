rm(list=ls())
gc

library(EBImage)
library(ggplot2)
library(doParallel)
library(foreach)

setwd("./TIL_maps")
registerDoParallel(12)

# list all cancers
cancers_all <- read.table("res/Cancer_summary.txt", header=T)
cancer <- cancers_all$Cancer
# cancer <- "brca"

area_thr <- 10
line_thr <- 0.2

area_perc_all <- list()
for(a in seq(1,length(cancer),1)){
  cat("Cancer: ",cancer[a]," \n")
  
  # load samples description
  load(file=paste0("data/Data_",cancer[a],".RData"))
  
  # create folders for saving
  dir.create(paste0("res/plots_regions/img_binary/",cancer[a]), showWarnings=F)
  dir.create(paste0("res/plots_regions/img_prob/",cancer[a]), showWarnings=F)
  dir.create(paste0("res/txt_files_regions/",cancer[a]), showWarnings=F)
  
  if (file.exists(paste0("data/Data_",cancer[a],".RData"))){
    data_all <- list()
    for (b in 1:nrow(samples)){
      cat("Sample: ",samples$Barcode[b]," \n")
      
      
      
      if (!length(list.files(path=paste0("res/txt_files_regions/",cancer[a]),pattern=samples$Barcode[b]))){
        # load data
        data <- data_all[[samples$Barcode[b]]]
        
        # check if image is not too small
        if(nrow(data) > 10000){
          # load images
          img_bin <- readImage(paste0("res/plots/img_binary/",cancer[a],"/",samples$Barcode[b],".png"))
          img_prob <- readImage(paste0("res/plots/img_prob/",cancer[a],"/",samples$Barcode[b],".png"))
          
       
                  
          # remove Background lines from plot
          img_bw <- imageData(img_bin < 1)[,,1]
          row_back <- rowSums(img_bw)
          diff_row_back <- diff(diff(row_back))
          # plot(diff_row_back)
          ind <- row_back >0
          diff_row_back[ind[2:(length(ind)-1)]] <- 0
          ind_row <- sort(which(diff_row_back > line_thr * length(row_back)) + 1)
          if (length(ind_row)){
            img_bin <- img_bin[-ind_row,,]
            img_prob <- img_prob[-ind_row,,]
            data <- data[!(data$x %in% ind_row),]
            for(r in seq(length(ind_row),1,-1)) data$x[data$x > ind_row[r]] <- data$x[data$x > ind_row[r]] - 1  
          } 
          col_back <- colSums(img_bw)
          diff_col_back <- diff(diff(col_back))
          ind <- col_back >0
          diff_col_back[ind[2:(length(ind)-1)]] <- 0
          ind_col <- sort(which(diff_col_back > line_thr * length(col_back)) + 1)
          if (length(ind_col)){
            img_bin <- img_bin[,-ind_col,]
            img_prob <- img_prob[,-ind_col,]
            ind_col <- sort(max(data$y)-ind_col+1)
            data <- data[!(data$y %in% ind_col),]
            for(r in seq(length(ind_col),1,-1)) data$y[data$y > ind_col[r]] <- data$y[data$y > ind_col[r]] - 1  
          } 
          
          if(dim(img_bin)[1] != length(unique(data$x)) | dim(img_bin)[2] != length(unique(data$y))){
            cat("Data not matched !!!\n")
          }
          
          # separate regions and calculate the area
          labels_bin <- bwlabel(img_bin<1)
          labels_area <- table(labels_bin)
          labels_area <- labels_area[2:length(labels_area)]
          labels_area_perc <- as.numeric(100 * labels_area/sum(labels_area))
          area_perc_all[[paste0(cancer[a], "_",samples$Barcode[b])]] <- labels_area_perc
          
          # remove regions lower than a threshold
          reg_del <- rownames(labels_area)[labels_area_perc < area_thr]
          mask <- labels_bin[,,1]
          mask[mask %in% reg_del] <- 0
          
          # separate regions, crop and save
          un_reg <- rownames(table(mask))
          if(length(un_reg) > 1){
            for (r in 2:length(un_reg)){
              img_bin_tmp <- img_bin
              img_prob_tmp <- img_prob
              img_bin_tmp[mask != un_reg[r]] <- 1
              img_prob_tmp[mask != un_reg[r]] <- 1
              ind_row <- rowSums(mask == un_reg[r]) > 0
              ind_col <- colSums(mask == un_reg[r]) > 0
              img_bin_tmp <- img_bin_tmp[ind_row,ind_col,]
              img_prob_tmp <- img_prob_tmp[ind_row,ind_col,]
              # plot(img_bin_tmp)
              
              # write cropped region
              writeImage(img_bin_tmp,files=paste0("res/plots_regions/img_binary/",cancer[a],"/",samples$Barcode[b],"_R",r-1,".png"))
              writeImage(img_prob_tmp,files=paste0("res/plots_regions/img_prob/",cancer[a],"/",samples$Barcode[b],"_R",r-1,".png"))
              
              
              # remove small regions from txt file
              ind <- which(imageData(mask) == un_reg[r], arr.ind=T)
              ind[,2] <- max(data$y) - ind[,2] + 1
              ind2 <- foreach (ii=1:nrow(ind), combine=c) %dopar%{which(data$x == ind[ii,1] & data$y == ind[ii,2])}
              ind2 <- unlist(ind2) 
              data_tmp <- data[ind2,]        
              write.table(data_tmp, file=paste0("res/txt_files_regions/",cancer[a],"/",samples$Barcode[b],"_R",r-1,".txt"), sep="\t", row.names=F, quote=F)
            } 
          }
          
        } else {
          cat("Less than 10 000 pixels in file!!!\n")
        }
      } 
      save(data_all, samples, file=paste0("data/Data_",cancer[a],".RData"))
      cat("No. of complete cases: ",nrow(samples)," \n")
    } 
  }
}
