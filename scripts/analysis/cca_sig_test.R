## ===========================
## Script name: cca_sig_test
##
## Purpose of script: Non-parametric significance testing between a composite variable (e.g. a variate from a cca or pls analysis)
## and a matrix of original metrics used to compute the composite variable (e.g. functional connectivity)
##
## Author: Sid Chopra (sid.chopra@monash.edu)
## Date Created: 24/07/2020
## ---------------------------
## Inputs: 
##
## variate = a numeric matrix of dim n x 1.
## edgeMetric = a numeric matrix of dimensions n x e (where e is the number of edges)
## b = an integer for home many bootstraps (method1) or permutations (method2) to compute. 
## method = a integer of value 0, 1, or 2
## Method 0: Returns uncorrected p-values, (same as method 1, but no FDR correction applied)
## Method 1: Returns FDR-corrected p-values by: 
## 1) Using bootstrapping to compute a standard error at each edge (SE)
## 2) Convert to Z-scores (r/se)
## 3) Convert Z-scores to two-tailed p-values and apply FDR corrected to these p-values. 
## Method 2 (Experimental): returns FWE-corrected p-values by:
## 1) Use permutation to compute a null distribution of r values at each edge (r) 
## 2) Create a fwe-null distribution be saving the max permuted r-value from each edge. 
## 3) compute a fwe-p-value using the original r value from step one and the distribution of max r values from step two. 
##
## ---------------------------
## Outputs: a vector of p-values
##
## ----------------------------
## Notes on using FDR/FWE: Alex: Method 1 and 2 seem to be answering different questions: 
## (1) is trying to get a robust estimate of the correlation that accounts for the sampling error of the statistics. 
## (2) is asking whether the correlation is larger than null expectation.
## FDR method is usually preferred and used in the literature. 
## ---------------------------
## Example usage
## mat <- rmatio::read.mat("a2_cca.mat")
## variate <- as.matrix(mat[["U"]][,1])
## edgeMetric <- read.table(patient_FC.txt)
## p <- cca_sig_test(variate, edgeMetric, method = 1, b=500)
##
## ---------------------------
## Required packages: Base r only
## ===========================

cca_sig_test <- function(variate, edgeMetric, method=1, b=500) {
  
  variate <- as.matrix(variate)
  edgeMetric <- as.matrix(edgeMetric)
  
  observerd_cor <- cor(variate, edgeMetric)
  
  if(method == 0 |  method == 1) {
    boot.se <- NULL
    for (e in 1:dim(edgeMetric)[2]){
      n <- dim(variate)[1]
      boot.cor.all<-NULL
      
      for (i in 1:b){
        index<-sample(1:n, replace=T)
        boot.patients.a2.resid<-as.matrix(edgeMetric[index,e])
        boot.brain_var<-as.matrix(variate[index,1])
        boot.cor<-cor(boot.patients.a2.resid, boot.brain_var)
        boot.cor.all<-c(boot.cor.all, boot.cor)
      }
      
      ## Bootstrap standard error
      boot.se[e] <- sd(boot.cor.all)
      print(paste("Edge:", e, b, "bootstraps."))
    }
    #hist(boot.se) 
    Z_score <- observerd_cor/boot.se
    #compute 2 sided pvalues
    pvalue2sided <- 2*pnorm(-abs(Z_score))
    
    #if method = 0, return uncorrected p values
    if(method==0) {return(t(pvalue2sided))}
    
    #do a FDR correction
    pvalue2sided_fdr <- p.adjust(pvalue2sided, method = "fdr", n = length(pvalue2sided))
    
    if(method==1) {return(pvalue2sided_fdr)}
  
    }
  
  

  
  if(method== 2) {
    max.perm.p <- NULL
    #set.seed(419)    # Make the results reproducible
    
    nperm <- b
    variate.perm <- numeric(nperm)
    n <- nrow(edgeMetric)
    max.vec <- NULL
    for (i in 1:dim(edgeMetric)[2]){
      for (ii in 1:nperm) {
        ind <- sample(n, replace = F)
        variate.perm[ii] <- cor(edgeMetric[ind,i], variate[,1])
      }
      max.vec[i] <- max(variate.perm[i])
      print(paste("Edge:", i, nperm, "permutations."))
    }
    #fwe_corrected pvalue
    p_vec_fwe <- NULL
    for (p in 1:dim(observerd_cor)[2]) {
      p_vec_fwe[p] <- sum(abs(max.vec) >= abs(observerd_cor[p])) / length(max.vec)
    }
    return(as.numeric(p_vec_fwe))
  
  }
}
