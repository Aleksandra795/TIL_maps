rm(list=ls())
gc()

library(apcluster)
library(clusterCrit)
library(foreach)
library(doParallel)

setwd("/mnt/data/Projects/HE/TIL_maps")
source("scripts/calc_spat_meas.R")

# list all cancers
cancers_all <- read.table("res/Cancer_summary.txt", header=T)
cancer <- cancers_all$Cancer
# cancer <- "brca"

# for(a in seq(length(cancer),1,-1)){
for(a in seq(2,8,2)){
  cat("Cancer: ",cancer[a]," \n")
  
  # create folders for saving
  dir.create(paste0("res/ap_plots/img_binary/",cancer[a]), showWarnings=F)
  dir.create(paste0("res/ap_plots/img_prob/",cancer[a]), showWarnings=F)
  
  # load samples description
  load(file=paste0("data/Data_",cancer[a],".RData"))
  
  # list all files in the folder
  files_all <-list.files(path=paste0("res/txt_files_regions/",cancer[a]), pattern="*.txt")
  
  # analyze all regions in parallel
  
  start_time <- Sys.time()
  spatial_meas <- mclapply(files_all, calc_spat_meas, mc.cores=12, mc.preschedule=F, mc.silent=T)
  cat(Sys.time() - start_time, "\n")
  # spatial_meas <- list()
  # for (file_name in files_all){
  #   cat(file_name,"\n")
  #   spatial_meas[[file_name]] <- calc_spat_meas(file_name)
  # }
  
  spatial_meas <- do.call(rbind,spatial_meas)
   
  # save results per cancer
  saveRDS(spatial_meas,file=paste0("res/Spatial_measures/",cancer[a],".rds"))
  
}

