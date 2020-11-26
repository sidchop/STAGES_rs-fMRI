

get_stages_fmri_data_a1a2 <- function() {
  #setwd("~/Dropbox/Sid/R_files/STAGES_fmri/data/ts_p300_n7/fmriprep_aroma_2p_hpfdt_gmr_NoSmooth/")
  
  #if working on cluster:
  #setwd("/home/scho0011/kg98/Sid/STAGES/STAGES_fmriprep/timeseries/schaefer_300parc_7net/fmriprep_aroma_2p_hpfdt_gmr_NoSmooth")
  
  
  filenames <- list.files(path = "~/Dropbox/Sid/R_files/STAGES_fmri/data/ts_p300_n7/fmriprep_aroma_2p_hpfdt_gmr_NoSmooth/", 
                          pattern = "ts", full.names = TRUE)  
  
  
  All_ts <- lapply(filenames, function(i){
    read.delim2(i, header = FALSE, sep = "", stringsAsFactors=FALSE)
  })
  
  corr_mat = list()
  
  ############ Creating FC mats for Analysis 1 Baseline FEPvHC (a1)
  #exclude subs
  #Headmotion, artefact exclude
  a2.exclude.qc <- c("011_1", "019_1", "070c_1", "036c_1", "049_2")
  #placebo group who were exposed to med exclude
  
  #bad carpet plots:35c_1 40c_1 61_1 70c_1 15c_2 32c_2 40_2, 49c_2
  # bbad conmats 70c_4
  
  a2.exclude.medexp <- c("035_2", "040_2","061_2", "062_2")
  
  
  filenames <- list.files(path = "~/Dropbox/Sid/R_files/STAGES_fmri/data/ts_p300_n7/fmriprep_aroma_2p_hpfdt_gmr_NoSmooth/", 
                          pattern = c("ts"), full.names = TRUE)  
  filenames.ses12 <- grep(filenames, pattern='_4', inv=T, value=T)
  filenames.ses12.qc <- grep(filenames.ses12, 
                             pattern='ts_sub-011_1.txt|ts_sub-019_1*.txt|ts_sub-070c_1.txt|ts_sub-036c_1.txt|ts_sub-035_2.txt|ts_sub-040_2.txt|ts_sub-061_2.txt|ts_sub-062_2.txt|ts_sub-049_2.txt', inv=T, value=T)
  
  
  a2_ts <- lapply(filenames.ses12.qc, function(i){
    read.delim2(i, header = FALSE, sep = "", stringsAsFactors=FALSE)
  })
  
  corr_mat <- list()
  
  for (i in 1:length(a2_ts)) {
    temp <- as.data.frame(a2_ts[i])
    #exclude low-signal ROIs (determined using the elbow method)
    temp <- temp[-c(87, 88, 89, 90, 91, 92, 94, 114, 238, 240, 241, 242, 243, 244, 245, 246), ] 
    indx <- sapply(temp, is.character)
    temp[indx] <- lapply(temp[indx], function(x) as.numeric(as.character(x)))
    corr_mat[[i]]<-cor(t(temp))
  }
  
  # Data laod complete
  #--------------------------------
  #Prepare vectors for the analysis
  # and declate number of BS
  #----------------------------------
  
#number of bootstraps
  parc <- dim(corr_mat[[1]])[1] #number of parcs
  edge.count <- (parc*(parc-1)/2)   #number of edges (N(N-1)/2)  #leave edgecount as is for creating con_mat_vec
  subs <- length(corr_mat)[1] #N subs
  corr_mat_vec <- matrix(nrow = subs, ncol = edge.count)
  
  
  
  #vectorise each subjects corr mat as a row in a matrix
  for (i in 1:subs) {
    corr_mat_vec[i,] <- t(corr_mat[[i]][upper.tri(corr_mat[[i]], diag = FALSE)])
  }
  return(corr_mat_vec)
}
