#========================================================
# This script extracts and writes out a binary adj-matrix 
#for the largest (max extant) component in the observed data
#=========================================================

#Read in -log10 wild-bootstrapped p-values for the observed data (from SwE toolbox in matlab) 
pstat <- readNIfTI2("swe_vox_Fstat_lp-WB_c01.nii")
pstat <- c(pstat)
pstat <- pstat[!is.na(pstat)]

#converting the vector of non-parametric -log10(p-values) back into a square matrix
p.vect <- pstat
temp <- matrix(nrow=316, ncol=316) 
temp[upper.tri(temp)] <- as.matrix(p.vect)
temp <- forceSymmetric(temp, uplo = "U")
diag(temp) <- 0

thresh <- c(0.05,0.01, 0.001) 
thresh_nlog10 <- -log10(thresh)
max_observerd_comp <- NA
for (s in 1:length(thresh_nlog10)) {
  p.mat.bin <- binarize(temp, threshold=thresh_nlog10[s]) 
  p.mat.bin.net <- as.network(as.matrix(temp_bin), directed = FALSE)
  plot(p.mat.bin.net)
  max.comps <- component.largest(p.mat.bin.net, return.as.edgelist = TRUE)
  p.mat.bin.maxcomp <- as.matrix(p.mat.bin)
  r <- rep(0, 316)
  for (i in 1:316) {
    ifelse(max.comps[i]==FALSE, p.mat.bin.maxcomp[i, ] <- as.vector(r), p.mat.bin.maxcomp[i, ] <- p.mat.bin.maxcomp[i, ])
  }
  diag(p.mat.bin.maxcomp) <- 0
  p.mat.bin.maxcomp <- as.matrix(forceSymmetric(p.mat.bin.maxcomp))
  #qgraph(p.mat.bin.maxcomp) [optiomal], visualise graph
  cluster_sig_edges <- which(p.mat.bin.maxcomp[upper.tri(p.mat.bin.maxcomp)] > 0)
  #write out largest cluster
  write.csv(x = as.matrix(p.mat.bin.maxcomp), file = paste0("max_bin_mat_",thresh_nlog10[s],".txt"), quote = FALSE, row.names = FALSE, col.names = FALSE )
}