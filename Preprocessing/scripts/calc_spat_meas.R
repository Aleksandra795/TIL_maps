calc_spat_meas <- function(file_name){
  name <- substring(file_name, 1, nchar(file_name)-4)
  cat("Sample: ",name," \n")
  spatial_meas_tmp <- data.frame(matrix(NA, nrow=1, ncol=9)) 
  colnames(spatial_meas_tmp) <- c("TILs_perc","Ball_Hall_discr","Banfeld_Raftery_discr","C_index_discr","Det_Ratio_discr",
                                  "Ball_Hall_prob","Banfeld_Raftery_prob","C_index_prob","Det_Ratio_prob")
  
  
  if (!file.exists(paste0("res/Spatial_measures/tmp/",file_name,".rds"))){
    # load data and select only TILs pixels
    data_raw <- read.delim(file=paste0("res/txt_files_regions/",cancer[a],"/",file_name))
    data_raw <- data_raw[!data_raw$Background,]
    data <- data_raw[data_raw$TILs==1,]
    
    # too big data cause memory error during clustering
    if (nrow(data) > 40000){
      ind <- sample.int(nrow(data),40000)
      data <- data[ind,]
    }
    
    if (nrow(data) < 10){
      cat("No TILs found!! \n")
    } else {
      # Scale data to range 0-1
      data$x <- (data$x - min(data$x))/diff(range(data$x))
      data$y <- (data$y - min(data$y))/diff(range(data$y))
      
      # Calculate affinity propagation
      cat("No. of pixels:",nrow(data),"\n")
      if (nrow(data) > 2000){
        frac <- 1000/nrow(data)
        cat("Short version with ",frac*nrow(data),"pixels.\n")
        apres1 <- apclusterL(negDistMat(r=2), data[,c("x","y")], q=0, seed=1,
                             frac=frac, sweeps=min(round(1/frac),5), maxits=100, convits=10) 
        apres2 <- apclusterL(negDistMat(r=2), data[,c("x","y","TILs_prob")], q=0, seed=1,
                             frac=frac, sweeps=min(round(1/frac),5), maxits=100, convits=10)
        
      } else {
        apres1 <- apcluster(negDistMat(r=2), data[,c("x","y")], q=0, seed=1, maxits=100, convits=10)
        apres2 <- apcluster(negDistMat(r=2), data[,c("x","y","TILs_prob")], q=0, seed=1, maxits=100, convits=10)
      }
      
      # Save the apcluster plots
      png(paste0("res/ap_plots/img_binary/",cancer[a],"/",name,".png"))
      plot(apres1,data[,c("x","y")])
      dev.off()
      png(paste0("res/ap_plots/img_prob/",cancer[a],"/",name,".png"))
      plot(apres2,data[,c("x","y")])
      dev.off()
      
      # calculate cluster metrics
      metrics1 <- intCriteria(as.matrix(data[,c("x","y")]),
                              as.integer(apres1@idx),
                              c("Ball_Hall","Banfeld_Raftery","C_index","Det_Ratio"))
      metrics2 <- intCriteria(as.matrix(data[,c("x","y","TILs_prob")]),
                              as.integer(apres2@idx),
                              c("Ball_Hall","Banfeld_Raftery","C_index","Det_Ratio"))
      
      # save results
      spatial_meas_tmp$TILs_perc <- 100 * sum(data_raw$TILs)/nrow(data_raw)
      spatial_meas_tmp$Ball_Hall_discr <- metrics1$ball_hall
      spatial_meas_tmp$Banfeld_Raftery_discr <- metrics1$banfeld_raftery
      spatial_meas_tmp$C_index_discr <- metrics1$c_index
      spatial_meas_tmp$Det_Ratio_discr <- metrics1$det_ratio
      spatial_meas_tmp$Ball_Hall_prob <- metrics2$ball_hall
      spatial_meas_tmp$Banfeld_Raftery_prob <- metrics2$banfeld_raftery
      spatial_meas_tmp$C_index_prob <- metrics2$c_index
      spatial_meas_tmp$Det_Ratio_prob <- metrics2$det_ratio   
      rownames(spatial_meas_tmp) <- name
      saveRDS(spatial_meas_tmp,file=paste0("res/Spatial_measures/tmp/",file_name,".rds"))
    }
  } else {
    cat("Result file exist.\n")
    spatial_meas_tmp <- readRDS(paste0("res/Spatial_measures/tmp/",file_name,".rds"))
  }
 
  return(spatial_meas_tmp)
}
