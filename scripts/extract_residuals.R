#Extract residuals
library("neurobase")

setwd("/projects/kg98/Sid/STAGES/STAGES_fmriprep/analyses/swe_validation/residuals/a3_covar_only_model/")
resid_files <- list.files(pattern = "resid")


baseline_resid_mat <- matrix(ncol=length(resid_files), nrow=49770)
for (i in 1:length(resid_files)){
  resid_vect <- c(readNIfTI2(paste(resid_files[i])))
  baseline_resid_mat[,i] <- resid_vect[!is.na(resid_vect)]
}
dim(baseline_resid_mat)
write.table(baseline_resid_mat, "a3_covars_only_resid_matrix.txt")
