write_largest_comp <- function(matx, filename) {
  #stopifnot(sum(mat==0) + sum(mat==1) == length(mat))
  p.mat.bin.net <- as.network(as.matrix(matx), directed = FALSE)

  max.comps <- sna::component.largest(p.mat.bin.net, return.as.edgelist = FALSE)
  
  p.mat.bin.maxcomp <- as.matrix(matx)

  r <- rep(0, dim(matx)[1])
  
  for (i in 1:dim(matx)[1]) {
    ifelse(max.comps[i]==FALSE, p.mat.bin.maxcomp[i, ] <- as.vector(r), p.mat.bin.maxcomp[i, ] <- p.mat.bin.maxcomp[i, ])
  }
  diag(p.mat.bin.maxcomp) <- 0
  p.mat.bin.maxcomp <- as.matrix(forceSymmetric(p.mat.bin.maxcomp, uplo = "U"))
  plot(as.network(as.matrix(p.mat.bin.maxcomp), directed = FALSE))
  #cluster_sig_edges <- which(p.mat.bin.maxcomp[upper.tri(p.mat.bin.maxcomp)] > 0)
  #setwd("~/Dropbox/Sid/R_files/STAGES_fmri/data/swe_validation_contrasts/gmr/a2_illness_effect/")
  #write out largest cluster
  write.csv(x = as.matrix(p.mat.bin.maxcomp), file = filename, quote = FALSE, row.names = FALSE, col.names = FALSE  )
  
}

