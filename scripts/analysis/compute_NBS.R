# ==============================================
# This script computes a Network Based Statistic (NBS) family-wise error corrected
# p-value, using brain-wide wald-statistics from the observed and null (wild bootstrapped) data.
# The wald-statistic used was generated using the SwE toolbox, but can also be done within
# R, using the clubSandwich toolbox, which has identical results (see X)
# This script was run separately for each of the 3 primary contrasts. 
# ==============================================


# Read in -log10 wild-bootstrapped p-values for the observed data (from SwE toolbox in matlab) 
pstat <- readNIfTI2("swe_vox_Fstat_lp-WB_c01.nii")
pstat <- c(pstat)
pstat <- pstat[!is.na(pstat)]

# Read in wald statistic for the null (wild-bootstrapped) data (from SwE toolbox in matlab) 
h5f = rhdf5::h5read("wb_wstat.mat", "wb_wstat")

# Read in -log10 p-values for the null (wild-bootstrapped) data (from SwE toolbox in matlab) 
h5f = rhdf5::h5read("wb_pvals.mat", "wb_pvals")
pmat <- as.data.frame(h5f)

#read in wmat our put from a2
#wmat <- as.data.frame(h5f)

#wmat <- rmatio::read.mat("wb_wstat.mat")
#wmat <- as.data.frame(wmat[[1]])
t_wmat <- t(wmat)
dim(t_wmat)



#==========
# Compute FWE null distribution
# For each primary 'component-forming' threshold, compute
# the size (extant) of the largest componant, for each bootstraped null.
# And write this out as a .txt file
#=========
alpha <- c(0.999, 0.99, 0.95) #primary 'component-forming' alpha
fwe.thresh <- 0.05

max.comp.list  <- vector(length=dim(t_pmat_inv)[2])

for (a in 1:length(alpha)){
  bs <- dim(t_pmat_inv)[2]
  for (k in 1:bs) {
    temp <- matrix(nrow=316, ncol=316)
    
    temp[upper.tri(temp)] <- as.matrix(t_pmat_inv[,k])
    
    temp <- forceSymmetric(temp, uplo = "U")
    diag(temp) <- 0
    # calculate nulls 
    temp_bin<- binarize(temp, threshold=alpha[a]) 
    temp.i <- graph_from_adjacency_matrix(as.matrix(temp_bin), 
                                          weighted = NULL,
                                          mode = c("undirected"))
    
    decomposed_comps <-  decompose.graph(temp.i)
    max.comp.list[k] <- max(unlist((lapply(decomposed_comps, gsize))))
    print(paste(k, ": max = ", max.comp.list[k]))
  }
  write.table(max.comp.list, paste("max.comp.list_",alpha[a], sep = ""), col.names = F, row.names = F)
}

max.comp.list.extant  <- vector(length=dim(t_pmat)[2])

#load in the null distributions created above
null.dist_0.05 <- scan("max.comp.list_0.95")
null.dist_0.01 <- scan("max.comp.list_0.99")
null.dist_0.001 <- scan("max.comp.list_0.999")

#===============================================================
#Calculate size (extant) of the largest component in observed data
#===============================================================
#converting the vector of non-parametric -log10(p-values) back into a square matrix
p.vect <- pstat
temp <- matrix(nrow=316, ncol=316) 
temp[upper.tri(temp)] <- as.matrix(p.vect)
temp <- forceSymmetric(temp, uplo = "U")
diag(temp) <- 0

thresh <- c(0.05,0.01, 0.001) 
thresh_nlog10 <- -log10(thresh)
max_observerd_comp <- NA
for(s in 1:length(thresh_nlog10)){
  temp_bin<- binarize(temp, threshold=thresh_nlog10[s]) 
  #View(igraph::degree(temp.i)) #[optional] - visualize graph
  temp.i <- graph_from_adjacency_matrix(as.matrix(temp_bin), 
                                        weighted = NULL,
                                        mode = c("undirected"))
  decomposed_comps <-  decompose.graph(temp.i)
  max_observerd_comp[s] <- max(unlist((lapply(decomposed_comps, gsize))))
}

#=============================================================================
# Visualize the fwe-null distribution, constructed from the size of the largest 
# component from each bootstrap, and where on that distribution the observed data
# falls (i.e. observed > 95th percentile/p_fwe < 0.016?)
#===========================================================================

#0.999
thresh <- quantile(null.dist_0.001, c(.984))
hist(null.dist_0.001, breaks = 100)
abline(v = thresh, col = "red", lwd = 2, lty = 2)
abline(v =   max_observerd_comp[3], col = "blue", lwd = 2, lty = 1)

#0.99
thresh <- quantile(null.dist_0.01, c(.984))
hist(null.dist_0.01, breaks = 100)
abline(v = thresh, col = "red", lwd = 2, lty = 2)
abline(v =   max_observerd_comp[2], col = "blue", lwd = 2, lty = 1)

#0.95
thresh <- quantile(null.dist_0.05, c(.984))
hist(null.dist_0.05, breaks = 100)
abline(v = thresh, col = "red", lwd = 2, lty = 2)
abline(v =   max_observerd_comp[1], col = "blue", lwd = 2, lty = 1)

