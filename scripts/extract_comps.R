#Make null distribution of max connectomes for specified contrasts 
packages <- c("readxl", "sandwich", "clubSandwich", "multiwayvcov", "clusterSEs",
              "contrast", "amen", "lmerTest", "voxel", "oro.nifti", "neurobase",
              "foreach", "doParallel", "MASS", "multcomp", "corrplot", "R.matlab", 
              "qgraph", "biclust", "igraph", "ggplot2", "pheatmap", "jtools","ggstance", "reshape2",
              "sjPlot", "sna", "ggm", "SDMTools", "Rfast", "tidyverse", "network")

lapply(packages, require, character.only = TRUE)
library("rhdf5")
library("neurobase")
################################################
#obs Fstat as matrix for matlab TFNBS analysis
###############################################

setwd("~/kg98/Sid/STAGES/STAGES_fmriprep/analyses/swe_validation/baseline_diff/")
fmat <- c(readNIfTI2("swe_vox_Fstat_c01.nii"))
fmat <- fmat[!is.na(fmat)]
write.table(fmat, 'obs_fstat.txt', col.names = F, row.names = F)

x <- read.csv('tfce_0.5_3.csv')



##################
#Illness effect 
##################


setwd("~/kg98/Sid/STAGES/STAGES_fmriprep/analyses/swe_validation/a3/med_v_placebo/")
pmat <- as.data.frame(1-t(as.data.frame(h5read("wb_pvals.mat", "wb_pvals"))))
compute_nulls_max(pmat)



##### Observerd data
pstat <- c(readNIfTI2("swe_vox_Fstat_lp-WB_c01.nii"))
pstat <- pstat[!is.na(pstat)]

compute_obs_max(mat = pstat, write.comps = T)

##################
#conjunction 1 - medication effect
###################
setwd("~/kg98/Sid/STAGES/STAGES_fmriprep/analyses/swe_validation/a3/med_v_hc/")
pmat1 <- as.data.frame(1-t(as.data.frame(h5read("wb_pvals.mat", "wb_pvals"))))

setwd("~/kg98/Sid/STAGES/STAGES_fmriprep/analyses/swe_validation/a3/med_v_placebo/")
pmat2 <- as.data.frame(1-t(as.data.frame(h5read("wb_pvals.mat", "wb_pvals"))))

setwd("~/kg98/Sid/STAGES/STAGES_fmriprep/analyses/swe_validation/a3/medication_effect/")

compute_nulls_max_conj(mat1 = pmat1, mat2 = pmat2)

##Compute observerd conjunction1
setwd("~/kg98/Sid/STAGES/STAGES_fmriprep/analyses/swe_validation/a3/med_v_hc/")
pstat1 <- c(readNIfTI2("swe_vox_Fstat_lp-WB_c01.nii"))
pstat1 <- pstat1[!is.na(pstat1)]

setwd("~/kg98/Sid/STAGES/STAGES_fmriprep/analyses/swe_validation/a3/med_v_placebo/")
pstat2 <- c(readNIfTI2("swe_vox_Fstat_lp-WB_c01.nii"))
pstat2 <- pstat2[!is.na(pstat2)]

setwd("/projects/kg98/Sid/STAGES/STAGES_fmriprep/analyses/swe_validation/a3/medication_effect/")
compute_obs_max_conj(mat1 = pstat1, mat2 = pstat2, write.comps = T)



##################
#conjunction 2 - medication protective/normalisation
###################
setwd("~/kg98/Sid/STAGES/STAGES_fmriprep/analyses/swe_validation/illness_effect/")
pmat1 <- as.data.frame(1-t(as.data.frame(h5read("wb_pvals.mat", "wb_pvals"))))

setwd("~/kg98/Sid/STAGES/STAGES_fmriprep/analyses/swe_validation/change_in_controls/")
pmat2 <- as.data.frame(1-t(as.data.frame(h5read("wb_pvals.mat", "wb_pvals"))))

setwd("/projects/kg98/Sid/STAGES/STAGES_fmriprep/analyses/swe_validation/illness_effect_conj_with_chang_in_controls/")
compute_nulls_max_conj2(mat1 = pmat1, mat2 = pmat2)

##Compute observerd conjunction2
setwd("~/kg98/Sid/STAGES/STAGES_fmriprep/analyses/swe_validation/illness_effect/")
pstat1 <- c(readNIfTI2("swe_vox_Fstat_lp-WB_c01.nii"))
pstat1 <- pstat1[!is.na(pstat1)]

setwd("~/kg98/Sid/STAGES/STAGES_fmriprep/analyses/swe_validation/change_in_controls/")
pstat2 <- c(readNIfTI2("swe_vox_Fstat_lp-WB_c01.nii"))
pstat2 <- pstat2[!is.na(pstat2)]

compute_obs_max_conj2(mat1 = pstat1, mat2 = pstat2, write.comps = TRUE)









##
setwd("/projects/kg98/Sid/STAGES/STAGES_fmriprep/analyses/swe_validation/medication_effect_2_conjunction/med_vs_hc/")
###
null.dist_0.05 <- scan("max.comp.list_0.95")
null.dist_0.01 <- scan("max.comp.list_0.99")
null.dist_0.001 <- scan("max.comp.list_0.999")


observed_0.999 = 1
observed_0.99 = 30
observed_0.95 = 402


#0.999
thresh <- quantile(null.dist_0.001, c(.95))
hist(null.dist_0.001, breaks = 100)
abline(v = thresh, col = "red", lwd = 2, lty = 2)
abline(v = observed_0.999, col = "blue", lwd = 2, lty = 1)


#0.99
thresh <- quantile(null.dist_0.01, c(.95))
hist(null.dist_0.01, breaks = 100)
abline(v = thresh, col = "red", lwd = 2, lty = 2)
abline(v = observed_0.99, col = "blue", lwd = 2, lty = 1)

#0.95
thresh <- quantile(null.dist_0.05, c(.95))
hist(null.dist_0.05, breaks = 100)
abline(v = thresh, col = "red", lwd = 2, lty = 2)
abline(v = observed_0.95, col = "blue", lwd = 2, lty = 1)



