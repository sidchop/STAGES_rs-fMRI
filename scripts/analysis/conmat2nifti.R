# ==============================================
# This script converts conectivity matricies into 
# .nii files, to they can be read into the SwE toolbox
# ==============================================
packages <- c("readxl", "sandwich", "clubSandwich", "multiwayvcov", "clusterSEs",
              "contrast", "amen", "lmerTest", "voxel", "oro.nifti", "neurobase",
              "foreach", "doParallel", "MASS", "multcomp", "corrplot", "R.matlab", 
              "qgraph", "biclust", "igraph", "ggplot2", "pheatmap", "jtools","ggstance", "reshape2",
              "sjPlot", "sna", "ggm", "SDMTools", "Rfast", "tidyverse")

lapply(packages, require, character.only = TRUE)


#===============================
#convert conmats from .txt to .nii
#==============================
setwd("~/path/to/conmats.txt")
library(neurobase)

fname = system.file(
  file.path("nifti", "mniRL.nii.gz"),
  package = "oro.nifti")
eve = neurobase::readnii(fname)
zeroes = niftiarr(eve, 0)

for (s in 1:dim(corr_mat_vec)[1]){
  a <- c(corr_mat_vec[s,], rep(0,852859))
  fimg <- neurobase::remake_img(a, zeroes)
  writeNIfTI(fimg, paste0("scan_",s), gzipped = FALSE, verbose = TRUE)
}


#================
#write out a mask
#=================

mask <- rep(1, 49770)
mask <- c(mask, rep(0,852859))
mask <- neurobase::remake_img(mask, zeroes)
writeNIfTI(mask, "mask", gzipped = FALSE, verbose = TRUE)
