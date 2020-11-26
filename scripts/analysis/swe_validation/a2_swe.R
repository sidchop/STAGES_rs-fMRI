# Analysis 2 - longitudinal differences - using swe+wb+nbs

#---------------------------------
# if running on a cluster and splitting chage the values below or pass on split num after call Rscript 
# 
#for split in {1..100} ; do export s=${split} ; sbatch batch_computeNBS.sh  ; done

#args <- commandArgs(TRUE)
N <- 1
#N  <- as.numeric(args[1]) #<---------------------- turn on for cluster 
splits <- 1 #change to one if not splitting
#
split.num <- N
#---------------------------------

packages <- c("readxl", "sandwich", "clubSandwich", "multiwayvcov", "clusterSEs",
              "contrast", "amen", "lmerTest", "voxel", "oro.nifti", "neurobase",
              "foreach", "doParallel", "MASS", "multcomp", "corrplot", "R.matlab", 
              "qgraph", "biclust", "igraph", "ggplot2", "pheatmap", "jtools","ggstance", "reshape2",
              "sjPlot", "sna", "ggm", "SDMTools", "Rfast", "tidyverse")

lapply(packages, require, character.only = TRUE)
#Load ts file for each sub (time series)
#If working locally:
setwd("~/Dropbox/Sid/R_files/STAGES_fmri/data/ts_p300_n7/fmriprep_aroma_2p_hpfdt_gmr_NoSmooth/")

#if working on cluster:
#setwd("/home/scho0011/kg98/Sid/STAGES/STAGES_fmriprep/timeseries/schaefer_300parc_7net/fmriprep_aroma_2p_hpfdt_gmr_NoSmooth")


filenames <- list.files(pattern = "ts", full.names = TRUE)  


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


filenames <- list.files(pattern = c("ts"), full.names = TRUE)  
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

b <- 999 #number of bootstraps
parc <- dim(corr_mat[[1]])[1] #number of parcs
edge.count <- (parc*(parc-1)/2)   #number of edges (N(N-1)/2)  #leave edgecount as is for creating con_mat_vec
subs <- length(corr_mat)[1] #N subs
corr_mat_vec <- matrix(nrow = subs, ncol = edge.count)



#vectorise each subjects corr mat as a row in a matrix
for (i in 1:subs) {
  corr_mat_vec[i,] <- t(corr_mat[[i]][upper.tri(corr_mat[[i]], diag = FALSE)])
}

#-----------------------------------
#load covars: Group, age_c, sex, fd_c
#------------------------------------
#if working on cluster 
#setwd("/home/scho0011/kg98/Sid/STAGES/STAGES_fmriprep/analyses/analysis1_baseline")
##if working locally
setwd("~/Dropbox/Sid/R_files/STAGES_fmri/data")

covs <- read.csv("covs_a2.csv") 





#run across multiple cores
#cores <- detectCores() #Turn off for cluster 
#registerDoParallel(cores-4) #Turn off for cluster 


#FUNCTIONS
wild <- function(n) sample(c(-1, 1), n, replace = TRUE)


#wild_b <- function(fitted, residual, b, beta) {
#  w.vect.temp  <- vector(length=b) 
##  mod.0 <- m
#  cli <- factor(covs$subject, levels = unique(covs$subject))
#   w.vect.temp <-  foreach(ii=1:b, .combine=rbind) %dopar% {
#    j <- wild(nlevels(cli))  
#    dv_star <<- as.vector((fitted + (residual*j[cli])))   ##muliply each subs residuals by -1/1 and add back to predictevalues
#    
#    dv_star_fit = lm(formula=dv_star ~ g1_s1 + g1_s2 + g2_s1 + g2_s2 + g3_s1 + g3_s2 + age_c + sex + fd_c + 0,  data=tempData)
#    
#    cov_star <- vcovCL(dv_star_fit, type = "HC2", cluster = ~covs$subject, fix = TRUE, cadjust = FALSE)     
#    
#    wald <- t(K %*% (dv_star_fit$coef-beta)) %*% solve((K %*% cov_star %*% t(K))) %*% (K %*% (dv_star_fit$coef-beta))/2 #####Divide this by rank of contrast!!
#    wald
#  } 
#  
#  return(w.vect.temp)
#}


#try apply-afying wb funcion
#generate 5000 samples from the radenmacher dist
randmat <- matrix(nrow = dim(covs)[[1]], ncol = b)
cli <- factor(covs$subject, levels = unique(covs$subject))

for (bb in 1:b) {
j <- wild(nlevels(cli))  
randmat[,bb] <- j[cli]
}

wild_app <- function(randvec, fitted, residual, b, beta){
  dv_star <<- as.vector((fitted + (residual*randvec)))   ##muliply each subs residuals by -1/1 and add back to predictevalues
  dv_star_fit = lm(formula=dv_star ~ g1_s1 + g1_s2 + g2_s1 + g2_s2 + g3_s1 + g3_s2 + age_c + sex + fd_c + 0,  data=tempData)
  cov_star <- vcovCL(dv_star_fit, type = "HC2", cluster = ~covs$subject, fix = TRUE, cadjust = FALSE)     
  wald <- t(K %*% (dv_star_fit$coef-beta)) %*% solve((K %*% cov_star %*% t(K))) %*% (K %*% (dv_star_fit$coef-beta))/1 #####Divide this by rank of contrast!!
  wald
}












adj_resid <- function(x, cluster = NULL, type = NULL, sandwich = TRUE, fix = FALSE, ...)
{
  ## compute meat of sandwich
  resid <- calc_resid(x, cluster = cluster, type = type, ...)
  
  return(resid)
}

calc_resid <- function(x, cluster = NULL, type = NULL, cadjust = TRUE, multi0 = FALSE, ...)
{
  ## extract estimating functions / aka scores
  if (is.list(x) && !is.null(x$na.action)) class(x$na.action) <- "omit"
  ef <- estfun(x)
  k <- NCOL(ef)
  n <- NROW(ef)
  
  ## set up return value with correct dimension and names
  rval <- matrix(0, nrow = k, ncol = k,
                 dimnames = list(colnames(ef), colnames(ef)))
  
  ## cluster can either be supplied explicitly or
  ## be an attribute of the model...FIXME: other specifications?
  if (is.null(cluster)) cluster <- attr(x, "cluster")
  
  ## resort to cross-section if no clusters are supplied
  if (is.null(cluster)) cluster <- 1L:n
  
  ## collect 'cluster' variables in a data frame
  if(inherits(cluster, "formula")) {
    cluster_tmp <- expand.model.frame(x, cluster, na.expand = FALSE)
    cluster <- model.frame(cluster, cluster_tmp, na.action = na.pass)
  } else {
    cluster <- as.data.frame(cluster)
  }
  
  ## handle omitted or excluded observations
  if((n != NROW(cluster)) && !is.null(x$na.action) && (class(x$na.action) %in% c("exclude", "omit"))) {
    cluster <- cluster[-x$na.action, , drop = FALSE]
  }
  
  if(NROW(cluster) != n) stop("number of observations in 'cluster' and 'estfun()' do not match")
  
  ## for multi-way clustering: set up interaction patterns
  p <- NCOL(cluster)
  if (p > 1L) {
    cl <- lapply(1L:p, function(i) combn(1L:p, i, simplify = FALSE))
    cl <- unlist(cl, recursive = FALSE)
    sign <- sapply(cl, function(i) (-1L)^(length(i) + 1L))    
    paste_ <- function(...) paste(..., sep = "_")
    for (i in (p + 1L):length(cl)) {
      cluster <- cbind(cluster, Reduce(paste_, unclass(cluster[, cl[[i]] ]))) ## faster than: interaction()
    }
    if(multi0) cluster[[length(cl)]] <- 1L:n
  } else {
    cl <- list(1)
    sign <- 1
  }
  
  ## number of clusters (and cluster interactions)
  g <- sapply(1L:length(cl), function(i) {
    if(is.factor(cluster[[i]])) {
      length(levels(cluster[[i]]))
    } else {
      length(unique(cluster[[i]]))
    }
  })
  #gmin <- min(g[1L:p])
  ## FIXME: additional argument for optionally using only smallest number of clusters?
  ## See also Cameron, Gelbach and Miller (2011, page 241)
  #if(FALSE) g[] <- gmin
  
  ## type of bias correction
  if(is.null(type)) {
    type <- if(class(x)[1L] == "lm") "HC1" else "HC0"
  }
  type <- match.arg(type, c("HC", "HC0", "HC1", "HC2", "HC3"))
  if(type == "HC") type <- "HC0"
  
  ## building blocks for HC2/HC3
  if(type %in% c("HC2", "HC3"))
  {
    if(any(g == n)) h <- hatvalues(x)
    
    if(!all(g == n)) {
      if(!(class(x)[1L] %in% c("lm", "glm"))) warning("clustered HC2/HC3 are only applicable to (generalized) linear regression models")
      
      ## regressor matrix
      X <- model.matrix(x)
      if(any(alias <- is.na(coef(x)))) X <- X[, !alias, drop = FALSE]
      attr(X, "assign") <- NULL
      
      ## working weights
      w <- weights(x, "working")
      
      ## (X'X)^(-1)
      XX1 <- if(is.null(w)) chol2inv(qr.R(qr(X))) else chol2inv(qr.R(qr(X * sqrt(w))))
      
      ## working residuals
      res <- rowMeans(ef/X, na.rm = TRUE)
      res[apply(abs(ef) < .Machine$double.eps, 1L, all)] <- 0
      
      ## matrix square root
      matpower <- function(X, p) {
        if((ncol(X) == 1L) && (nrow(X) == 1L)) return(X^p)
        Xeig <- eigen(X, symmetric = TRUE)
        if(any(Xeig$values < 0)) stop("matrix is not positive semidefinite")
        sqomega <- diag(Xeig$values^p)
        return(Xeig$vectors %*% sqomega %*% t(Xeig$vectors))
      }
    }
  }
  
  ## add OPG for each cluster-aggregated estfun
  for (i in 1L:length(cl))
  {
    ## estimating functions for aggregation by i-th clustering variable
    efi <- ef
    
    ## add cluster adjustment g/(g - 1) or not?
    ## only exception: HC0 adjustment for multiway clustering at "last interaction"
    adj <- if(multi0 & (i == length(cl))) {
      if(type == "HC1") (n - k)/(n - 1L) else 1
    } else {
      if(cadjust) g[i]/(g[i] - 1L) else 1
    }  
    
    ## HC2/HC3
    if(type %in% c("HC2", "HC3")) {
      if(g[i] == n) {
        efi <- if(type == "HC2") {
          efi/sqrt(1 - h)
        } else {
          efi/(1 - hatvalues(x))
        }
      } else {
        for(j in unique(cluster[[i]])) {
          ij <- which(cluster[[i]] == j)
          Hij <- if(is.null(w)) {
            X[ij, , drop = FALSE] %*% XX1 %*% t(X[ij, , drop = FALSE])
          } else {
            X[ij, , drop = FALSE] %*% XX1 %*% t(X[ij, , drop = FALSE]) %*% diag(w[ij], nrow = length(ij), ncol = length(ij))
          }
          Hij <- if(type == "HC2") {
            matpower(diag(length(ij)) - Hij, -0.5)
          } else {
            solve(diag(length(ij)) - Hij)
          }
          efi[ij, ] <- drop(Hij %*% res[ij]) * X[ij, , drop = FALSE]
        }
      }
      
      ## "inverse" cluster adjustment that Bell & McCaffrey (2002) and hence also
      ## Cameron & Miller (2005, Eq. 25) recommend for HC3 (but not HC2)
      ## -> canceled out again if cadjust = TRUE
      efi <- sqrt((g[i] - 1L)/g[i]) * efi
    }
    
    ## aggregate within cluster levels      
    #efi <- if(g[i] < n) apply(efi, 2L, rowsum, cluster[[i]]) else efi
    
    efi_res <- rowMeans(efi/X, na.rm = TRUE) #added to extract adj res
    ## aggregate across cluster variables
    #rval <- rval + sign[i] * adj * crossprod(efi)/n
  }
  
  ## HC1 adjustment with residual degrees of freedom: (n - 1)/(n - k)
  if(type == "HC1") rval <- (n - 1L)/(n - k) * rval
  
  return(efi_res)
}






#########end of functions #########3

#use this if doing two layer bs
#wild_b <- function(m, bb, beta) {
#  w.vect.temp_2  <- vector(length=bb) 
#  mod.0 <- m
#  cli <- factor(covs$subject, levels = unique(covs$subject))
#for(iii in 1:bb) { #chnage to dopar for paralell
#    j <- wild(nlevels(cli))  
#    dv_star <- as.vector((mod.0$fit + (mod.0$residuals*j[cli])))   ##muliply each subs residuals by -1/1 and add back to predictevalues
#    
#    dv_star_fit = lm(formula=dv_star ~ g1_s1 + g1_s2 + g2_s1 + g2_s2 + g3_s1 + g3_s2 + age_c + sex + fd_c + 0,  data=tempData)
#    
#    cov_star <- vcovCL(dv_star_fit, type = "HC2", cluster = ~covs$subject, fix = TRUE, cadjust = FALSE)     
#    
#    w.vect.temp_2[iii]  <- t(K %*% (dv_star_fit$coef-beta)) %*% solve((K %*% cov_star %*% t(K))) %*% (K %*% (dv_star_fit$coef-beta))/2 #####Divide this by rank of contrast!!
#    
#  } 
#  
#  return(as.vector(w.vect.temp_2))
#}
#

##permute edges by wild bootstrap
# use this if doing 2layred bootstrap 
#wild_b_2 <- function(m, b, beta) {
#  p.vect.temp  <- vector(length=b) 
#  w.vect.temp  <- vector(length=b) 
#  w.vect.temp_2 <- 10*b
#  mod.0 <- m
#  cli <- factor(covs$subject, levels = unique(covs$subject))
#for(ii in 1:b){  
#    j <- wild(nlevels(cli))  
#    dv_star <- as.vector((mod.0$fit + (mod.0$residuals*j[cli])))   ##muliply each subs residuals by -1/1 and add back to predictevalues
#    
#    dv_star_fit = lm(formula=dv_star ~ g1_s1 + g1_s2 + g2_s1 + g2_s2 + g3_s1 + g3_s2 + age_c + sex + fd_c + 0,  data=tempData)
#    
#    cov_star <- vcovCL(dv_star_fit, type = "HC2", cluster = ~covs$subject, fix = TRUE, cadjust = FALSE)     
#    
#    wald <- t(K %*% (dv_star_fit$coef-beta)) %*% solve((K %*% cov_star %*% t(K))) %*% (K %*% (dv_star_fit$coef-beta))/2 #####Divide this by rank of contrast!!
#    
#    w.vect.temp[ii] <- wald
#    
#    coef_2 <- dv_star_fit$coefficients
#    mod2 <- dv_star_fit
#    b2 <- 10
#    temp <- wild_b(m=mod2, bb=b2, beta=coef_2)
#    ifelse(ii==1,  w.vect.temp_2 <- temp,  w.vect.temp_2 <- c(w.vect.temp_2, temp))
#  print(ii)
#} 
#  percentile.t <<- ecdf(w.vect.temp_2) ##options(scipen=999) (disable sci notation) ---- calulating WB -pvalue (minimum p value =  1/(BS+1))
#  for (p in 1:length(p.vect.temp)) {
#    p.vect.temp[p] <- (1 - percentile.t(w.vect.temp[p])) 
#  }
#  return(list(V1=p.vect.temp,V2=w.vect.temp))
#}


#-----------------------------------
# If splitting file to run on cluster change the starting edge.number, do not need to change anything here, 
# just working out the first and last edge number of each split
# only change splits, and split.num at the top of the script
first.edge <- 1 
edge.count <-  (parc*(parc-1)/2)/splits # potentially redundant 
if(splits==1) {
  edge.count <- edge.count
  first.edge <- 1 
} else {
  x <- vector(length = edge.count*splits)
  n <- splits
  x[1:length(x)] <- 1:length(x)
  split.x <- split(x, sort(x%%n))
  first.edge <- split.x[[split.num]][1] 
  edge.count <- split.x[[split.num]][edge.count] 
}

#----------------------------------

#use this if doing 2layer bs
#p.vect.null <- matrix(nrow=b, ncol = (parc*(parc-1)/2)/splits)

w.vect.null <- matrix(nrow=b, ncol = (parc*(parc-1)/2)/splits)
w.vect.max <- w.vect <- p.vect <-vector(length=(parc*(parc-1)/2)/splits) 


K <- matrix(c(1, -1, 0, 0, -1, 1, 0, 0, 0, 0, 0, 1, -1, -1, 1, 0, 0, 0), 
            ncol = 9, nrow = 2, byrow = TRUE)

null.index <- 1 # so the vector of nulls in not written to the nth*split location of the vector
for (i in first.edge:edge.count) {   #for each edge
  tempData <- mod0 <- mod1 <- stats <- NULL
  tempData <- cbind(covs[,-10], corr_mat_vec[,i]) #add covars
  names(tempData)[10] <- "connectivity"
  mod1 <- lm(connectivity ~ g1_s1 + g1_s2 + g2_s1 + g2_s2 + g3_s1 + g3_s2 + age_c + sex + fd_c + 0,  data=tempData)
  cov <- vcovCL(x = mod1, type = "HC2", cluster = ~covs$subject, fix = TRUE, cadjust = FALSE)     

  #####Divide this by rank of contrast !!######## also change in bs function
  wald_1 <- t(K %*% mod1$coef) %*% solve((K %*% cov %*% t(K))) %*% (K %*% mod1$coef)/2
  #####Divide above by rank of contrast !!######## also change in BS function
  w.vect[null.index] <- wald_1 
  coef <- mod1$coef
  residual <- adj_resid(x = mod1, type = "HC2", cluster = ~covs$subject, fix = TRUE, cadjust = FALSE)
  w.vect.temp <- apply(randmat, 2, wild_app, fitted=mod1$fit,residual=residual,
             beta = coef)
#w.vect.temp <- wild_b(fitted=mod1$fit, residual=residual, b = b, beta = coef) 
# 
#  #use this for 2 layer bs
# #p.w.list <- wild_b_2(m = mod1, b = b, beta = coef) #change number of bs here 
# #w.vect.null <-  p.w.list[[2]]
# #p.vect.temp <-  p.w.list[[1]]
# 
  w.vect.max[null.index] <- max(w.vect.temp)
  w.vect.null[,null.index] <- w.vect.temp
  percentile.t <<- ecdf(w.vect.temp) ##options(scipen=999) (disable sci notation) ---- calulating WB -pvalue (minimum p value =  1/(BS+1))
  p.vect[null.index] <- (1 - percentile.t(wald_1)) 
  null.index <- null.index + 1
  print(paste0("edge number: ", i))
}

##Change table here to where you want the files written 
setwd("/home/scho0011/kg98/Sid/STAGES/STAGES_fmriprep/analyses/analysis1_baseline/gmr_noSmooth")


write.table(round(w.vect.null, digits = 4), paste("w.vect.null_", split.num, sep = ""), col.names = F, row.names = F)
write.table(round(p.vect, digits = 4), paste("p.vect_", split.num, sep = ""), col.names = F, row.names = F)
write.table(round(w.vect.max, digits = 4), paste("w.vect.max_", split.num, sep = ""), col.names = F, row.names = F)





# how to read the null vector back to concat with other nulls 
#x.t.vect.null <- as.matrix(read.table("test", header = F))

#-----------------------------------------------------------------------------------------
# Convert each row on the t.vect.null in to adj matrix to compute the size (no. of edges) 
# of the largest componant at differnt alpha thresholds.  
#-----------------------------------------------------------------------------------------

##Creater a vector with the critical values of the thresholds you want applied to the nulls

#use this to calcluate critical thresholds at differnt dfs and alphas qt(.999, 82)
alpha <- c(9.9558, 6.1756, 3.852) #derived from quantiles of the w.vect dist
fwe.thresh <- 0.05



w.vect.null.t <- t(w.vect.null)# collums are now null images  
max.comp.list  <- vector(length=dim(w.vect.null.t)[2])



for (a in 1:length(alpha)){
  bs <- dim(w.vect.null.t)[2]
  for (k in 1:bs) {
    temp <- matrix(nrow=316, ncol=316)
  
    temp[upper.tri(temp)] <- as.matrix(w.vect.null.t[,k])
    
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

# removed as we are now using a 2 layer WB
#for (a in 1:length(alpha)){
#  t.vect.null.t <- t(t.vect.null)# collums are now null images  
#  max.comp.list <- vector(length=dim(t.vect.null.t)[2])
#  bs <- dim(t.vect.null.t)[2]
#  for (k in 1:bs) {
#    temp <- matrix(nrow=dim(corr_mat[[1]])[1], ncol=dim(corr_mat[[1]])[1])
#    temp[upper.tri(temp)] <- as.matrix(t.vect.null.t[,k])
#    temp <- forceSymmetric(temp, uplo = "U")
#    diag(temp) <- 1
#    
#    # calculate nulls 
#    temp_bin<- binarize(temp, threshold=critical.thresh[a]) # need to find proper threshold. Change to match threshold above 
#    temp.i <- as.igraph(qgraph(temp_bin))
#    decomposed_comps <-  decompose.graph(temp.i)
#    max.index <- which.max(clusters(temp.i)$csize)
#    ifelse(length(max.index)!=0, max.comp.list[k] <- gsize(decomposed_comps[[max.index]]), max.comp.list[k] <- 0)
#    print(k)
#  }
#  write.table(max.comp.list, paste("max.comp.list_",alpha[a],"_", split.num, sep = ""), col.names = F, row.names = F)
#}




#===================================================================================
#The above (NBS computation) should be probably run on the cluster if bs/permutations are >999
#====================================================================================
setwd("~/Dropbox/Sid/R_files/STAGES_fmri/data/a2/res_fixed/")
#Load in p.vect
#Load and concat observerd p.vales
p.vect <- vector()
for (p in 1:630) {
  temp <- scan(file = as.character(paste0("p.vect_",p)))
  p.vect <- c(p.vect, temp)
  remove(temp)
}


#Load in w.vect.max
#Load and concat observerd p.vales
w.vect.max <- vector()
for (p in 1:630) {
  temp <- scan(file = as.character(paste0("w.vect.max_",p)))
  w.vect.max <- c(w.vect.max, temp)
  remove(temp)
}


#Load in w.vect.max
#Load and concat observerd p.vales ####fix
#p.vect.null <- matrix(nrow = 49770, ncol = 9999)
#Load and concat observerd p.vales
w.vect.null <- matrix()
for (p in 1:630) {
  temp <- read.table(file = as.character(paste0("w.vect.null_",p)))
  w.vect.null <- cbind(w.vect.null, temp)
  remove(temp)
}
w.vect.null <- w.vect.null[,-1]


#----------------------------------

##Calculate largest comp in observerd data
temp <- matrix(nrow=316, ncol=316) #next 3 lines of code just converting the pvector back into a con mat
temp[upper.tri(temp)] <- as.matrix(1-p.vect)
temp <- forceSymmetric(temp, uplo = "U")
diag(temp) <- 0

temp_bin<- binarize(temp, threshold=0.999) # this is what to change when looking at sig comps
#View(igraph::degree(temp.i))
temp.i <- graph_from_adjacency_matrix(as.matrix(temp_bin), 
                                      weighted = NULL,
                                      mode = c("undirected")) #maybe need to add as.matrix

# splits into comps (finds largest comp and gives size edges)
decomposed_comps <-  decompose.graph(temp.i)
max(unlist((lapply(decomposed_comps, gsize))))

####edge-wise FC
hist(w.vect.max)
quantile(w.vect.max, c(.95)) 
which(w.vect > 14.9919)
max(w.vect)


#critical.thresh <- c(8.1958, 5.2034, 3.27) ###Calculate these using a mega distribution


####load null vectors at each thresh
setwd("~/Dropbox/Sid/R_files/STAGES_fmri/data/a2/")
null.dist.f_0.05 <- scan("res_fixed/max.comp.list_3.852")
null.dist.f_0.01 <- scan("res_fixed/max.comp.list_6.1756")
null.dist.f_0.001 <- scan("res_fixed/max.comp.list_9.9558")



observed_0.999 = 19
observed_0.99 = 491
observed_0.95 = 2489

#0.999
thresh <- quantile(null.dist.f_0.001, c(.95))
hist(null.dist.f_0.001)
abline(v = thresh, col = "red", lwd = 2, lty = 2)
abline(v = observed_0.999, col = "blue", lwd = 2, lty = 1)


#0.99
thresh <- quantile(null.dist.f_0.01, c(.95))
hist(null.dist.f_0.01)
abline(v = thresh, col = "red", lwd = 2, lty = 2)
abline(v = observed_0.99, col = "blue", lwd = 2, lty = 1)

#0.95
thresh <- quantile(null.dist.f_0.05, c(.95))
hist(null.dist.f_0.05)
abline(v = thresh, col = "red", lwd = 2, lty = 2)
abline(v = observed_0.95, col = "blue", lwd = 2, lty = 1)








###Now that we have a comp (p<0.05 FWE - (cft p<0.001))
# What is the effect? 



##extract binary adj-matrix for the largest comp 
p.mat.bin <- temp_bin
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

qgraph(p.mat.bin.maxcomp)
cluster_sig_edges <- which(p.mat.bin.maxcomp[upper.tri(p.mat.bin.maxcomp)] > 0)
setwd("~/Dropbox/Sid/R_files/STAGES_fmri/data/a2")
#write out largest cluster
write.csv(x = as.matrix(p.mat.bin.maxcomp), file = "test_clust_0.999_a2.txt", quote = FALSE, row.names = FALSE, col.names = FALSE  )

x <- as.matrix(read.csv("test_clust_0.999_a2.txt"))

library(ggconn)
setwd("~/Dropbox/Sid/R_files/ggconn")
ggconn(atlas = "Stages_melbBrain", conmat = x, node.size = 1, labels = T)



#search a index to find the vector/edge  number
index_vect <- vector()
index_vect <- 1:49770
t.mat <- matrix(nrow=316, nco=316)
setwd("~/Dropbox/Sid/R_files/STAGES_fmri/data/a1")
labs <- read.csv(file = "labs_melbbrain.csv", header = F) 
rownames(t.mat) <- labs$V1
colnames(t.mat) <- labs$V1
t.mat[upper.tri(t.mat)] <- index_vect #mod for 1-p, log for -log(10)
t.mat <- forceSymmetric(t.mat, uplo = "U")
diag(t.mat) <- 1
t.mat[107,4]



## categotise edges as illness relates, medication related, or other
long.plots <- list()
for (a in 1:length(cluster_sig_edges)){
 # index <- t.mat[cluster_sig_edges[a,1], cluster_sig_edges[a,2]]
  index <- cluster_sig_edges[a] 
  tempData <- cbind(covs[,-10], corr_mat_vec[, index]) #add covars
  names(tempData)[10] <- "connectivity"
  mod1 <- lm(connectivity ~ g1_s1 + g1_s2 + g2_s1 + g2_s2 + g3_s1 + g3_s2 + age_c + sex + fd_c + 0,  data=tempData)
  cov <- vcovCL(mod1, type = "HC3", cluster = ~covs$subject, fix = TRUE, cadjust = FALSE)     
  #####Divide this by rank of contrast !!######## also change in bs function
  wald_1 <- t(K %*% mod1$coef) %*% solve((K %*% cov %*% t(K))) %*% (K %*% mod1$coef)/2
  #####Divide above by rank of contrast !!######## also change in BS function
  SEmat <- sqrt(diag(cov[1:6,1:6]))
  plotData <- matrix(nrow=6, ncol=4)
  colnames(plotData) <- c("coef", "ses", "group", "se")
  plotData[,1] <- t(mod1$coef[1:6])
  plotData[,2] <- t(c(1,2,1,2,1,2))
  plotData[,3] <- t(c(1,1,2,2,3,3))
  plotData[,4] <- t(SEmat)
  plotData[,3] <- factor(plotData[,3],
                         levels = c(1,2,3),
                         labels = c("placebo", "medication", "healthy"))
  pd <- position_dodge(0.1)
  p <-  ggplot(data=as.data.frame(plotData), aes(x=as.factor(ses), y=coef, group=group)) + 
    geom_point(aes(color=as.factor(group)), position=pd) + 
    geom_line(aes(group=group, color=as.factor(group)), position=pd) +
    geom_errorbar(aes(ymin=coef-se, ymax=coef+se, color=as.factor(group)), width=.1, position=pd) +
    ggtitle(paste(rownames(which(t.mat==cluster_sig_edges[a], arr.ind = T))[1], "<-->",
                  rownames(which(t.mat==cluster_sig_edges[a], arr.ind = T))[2]))
  #p <- p + theme(legend.position="none")
  long.plots[[a]] <- p
  
 # Kmed <- matrix(c(0,1,0,0,0,-1,0,0,0), nrow=1, ncol=9, byrow = TRUE)
 # t(Kmed %*% mod1$coef) %*% solve((Kmed %*% cov %*% t(Kmed))) %*% (Kmed %*% mod1$coef)
 # t <- glht(mod1, linfct = Kmed)
 # summary(t)
}
long.plots[[5]]

library("patchwork")
long.plots[[1]] + long.plots[[2]] + long.plots[[3]] +
  long.plots[[4]] + long.plots[[5]] + long.plots[[6]] +
  long.plots[[7]] + long.plots[[8]] + long.plots[[9]] +
  long.plots[[10]] + long.plots[[11]] + plot_layout(ncol=3)




w.vect.temp <- wild_b(m = mod1, b = b, beta = coef) #change number of bs here 
w.vect.max[null.index] <- max(w.vect.temp)
w.vect.null[,null.index] <- w.vect.temp
percentile.t <<- ecdf(w.vect.temp) ##options(scipen=999) (disable sci notation) ---- calulating WB -pvalue (minimum p value =  1/(BS+1))
p.vect[null.index] <- (1 - percentile.t(wald_1)) 




##visualise using ggconnectome
source("~/Dropbox/Sid/R_files/ggconnectome/functions/add_geoms.R")
source("~/Dropbox/Sid/R_files/ggconnectome/functions/ggConnectome.R")
source("~/Dropbox/Sid/R_files/ggconnectome/functions/ggConnectome3D.R")

ggConnectome(atlas = "Stages_melbBrain", conmat = x)
ggConnectome3D(atlas = "", conmat = x)


##################################################
# TFNBS TFNBS TFNBS TFNBS TFNBS TFNBS TFNBS TFNBS 
#TFNBS TFNBS TFNBS TFNBS TFNBS TFNBS TFNBS TFNBS 
##################################################
setwd("~/Dropbox/Sid/R_files/STAGES_fmri/data/a2/tfnbs/")
obs_fstat <- read.csv("a2_wstat.txt", header = F)

hist(c(unlist(obs_fstat)), breaks = 10000)

obs_tfnbs <- tfnbs(E = 0.5, H=2.25, mat = obs_fstat)

pheatmap::pheatmap(temp, cluster_rows = F, cluster_cols = F)

as.igraph(qgraph(temp_mat))
hist(obs_tfnbs, breaks=10000)

View(obs_tfnbs)

nulls1 <- readRDS("MYLIST_1000.Rds")
nulls2 <- readRDS("MYLIST_2000.Rds")
nulls3 <- readRDS("MYLIST_3000.Rds")
nulls2 <- plyr::compact(nulls2)
nulls3 <- plyr::compact(nulls3)

listonulls <- append(nulls1, nulls2)
listonulls <- append(listonulls, nulls3)

max_vect <- vector()
for (i in 1:length(listonulls)){
  max_vect[i] <- max(listonulls[[i]])
}
hist(max_vect, breaks = 1000)
abline(v = quantile(max_vect, 0.95), col = "red", lwd = 2, lty = 2)




#Visualising the wild bs ----
write.table(x = tempData)

