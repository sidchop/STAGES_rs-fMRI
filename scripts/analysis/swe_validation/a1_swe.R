# Analysis 1 - Baseline differences - using standard t-test

#---------------------------------
# if running on a cluster and splitting chage the values below or pass on split num after call Rscript 
# 
#for split in {1..10} ; do export s=${split} ; sbatch batch_computeNBS.sh  ; done

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
a1.exclude <- c("011", "019", "070c", "036c")
filenames <- list.files(pattern = c("ts"), full.names = TRUE)  
filenames.ses1 <- grep(filenames, pattern='_1', inv=F, value=T)
filenames.ses1.qc <- grep(filenames.ses1, 
                          pattern='*011_*|*019_*|*070c_*|*036c_*', inv=T, value=T)


a1_ts <- lapply(filenames.ses1.qc, function(i){
  read.delim2(i, header = FALSE, sep = "", stringsAsFactors=FALSE)
})



corr_mat <- list()

for (i in 1:length(a1_ts)) {
  temp <- as.data.frame(a1_ts[i])
  #exclude low-signal ROIs (determined using the elbow method)
  temp <- temp[-c(87, 88, 89, 90, 91, 92, 94, 114, 238, 240, 241, 242, 243, 244, 245, 246), ] 
 rownames(temp) <- temp[,1] 
  temp[indx] <- lapply(temp[indx], function(x) as.numeric(as.character(x)))
  corr_mat[[i]]<-cor(t(temp))
}

#if you want to write out the conmats
for (k in 1:length(corr_mat)) {
  write.table(corr_mat[[k]], paste0(k,".txt"), row.names = F, col.names = F)
}

View(Stages_melbBrain
     )
x <- Stages_melbBrain[,1]
write.table(x,"labels.txt", row.names = F, col.names = F)

#
# Data laod complete
#--------------------------------
#Prepare vectors for the analysis
# and declate number of BS
#----------------------------------

b <- 9999 #number of bootstraps
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
setwd("~/Dropbox/Sid/R_files/STAGES_fmri/data/")

covs <- read.csv("covs_a1.csv") 

#run across multiple cores
cores <- detectCores() #Turn off for cluster 
registerDoParallel(cores=cores-2) #Turn off for cluster 


#FUNCTION
wild <- function(n) sample(c(-1, 1), n, replace = TRUE)
##permute edges by wild bootstrap
wild_b <- function(m, b, beta) {
  w.vect.temp  <- vector(length=b) 
  mod.0 <- m
  cli <- factor(covs$subject, levels = unique(covs$subject))
  w.vect.temp <-  foreach(ii=1:b, .combine=rbind) %dopar% {
    j <- wild(nlevels(cli))  
    dv_star <<- as.vector((mod.0$fit + (mod.0$residuals*j[cli])))   ##muliply each subs residuals by -1/1 and add back to predictevalues
    
    dv_star_fit = lm(formula=dv_star ~  group_1 + group_2 + age_c + sex + fd_c + 0,  data=tempData)
    
    cov_star <- vcovCL(dv_star_fit, type = "HC2", cluster = ~covs$subject, fix = TRUE, cadjust = FALSE)     
    
    wald <- t(K %*% (dv_star_fit$coef-beta)) %*% solve((K %*% cov_star %*% t(K))) %*% (K %*% (dv_star_fit$coef-beta)) #####Divide this by rank of contrast!!
    wald
    #here we would could also calculate a p value for cluster thresholding
    #  install_github("bozenne/lavaSearch2")
    #  library(lavaSearch2)
    #  dfSigmaRobust(K, vcov(dv_star_fit), cov_star,  dv_star_fit$model )
  } 
  
  return(w.vect.temp)
}



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

w.vect.null <- matrix(nrow=b, ncol = edge.count) 
w.vect.max <- w.vect <- p.vect <-vector(length=edge.count) 

K <-  matrix(c(-1, 1, 0, 0, 0), ncol = 5, nrow = 1, byrow = TRUE) # set contrast 


for (i in first.edge:edge.count) {   #for each edge
  tempData <- mod0 <- mod1 <- stats <- NULL
  tempData <- cbind(covs, corr_mat_vec[,i]) #add covars
  names(tempData)[7] <- "connectivity"
  mod1 <- lm(connectivity ~ group_1 + group_2 + age_c + sex + fd_c + 0,  data=tempData)
  cov <- vcovCL(mod1, type = "HC3", cluster = ~covs$subject, fix = TRUE, cadjust = FALSE)     
  #####Divide this by rank of contrast !!######## also change in bs function
  wald_1 <- t(K %*% mod1$coef) %*% solve((K %*% cov %*% t(K))) %*% (K %*% mod1$coef)
  #####Divide above by rank of contrast !!######## also change in BS function
  w.vect[i] <- wald_1 
  coef <- mod1$coef
  w.vect.temp <- wild_b(m = mod1, b = b, beta = coef) #change number of bs here 
  w.vect.max[i] <- max(w.vect.temp)
  w.vect.null[,i] <- w.vect.temp
  percentile.t <<- ecdf(w.vect.temp) ##options(scipen=999) (disable sci notation) ---- calulating WB -pvalue (minimum p value =  1/(BS+1))
  p.vect[i] <- (1 - percentile.t(wald_1)) 
  print(i)
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
alpha <- c(0.999, 0.99, 0.95)
critical.thresh <- c(3.193, 2.373, 1.664)
fwe.thresh <- 0.05


for (a in 1:length(alpha)){
  t.vect.null.t <- t(t.vect.null)# collums are now null images  
  max.comp.list <- vector(length=dim(t.vect.null.t)[2])
  bs <- dim(t.vect.null.t)[2]
  for (k in 1:bs) {
    temp <- matrix(nrow=dim(corr_mat[[1]])[1], ncol=dim(corr_mat[[1]])[1])
    temp[lower.tri(temp)] <- as.matrix(t.vect.null.t[,k])
    temp <- forceSymmetric(temp, uplo = L)
    diag(temp) <- 1
    
    # calculate nulls 
    temp_bin<- binarize(temp, threshold=critical.thresh[a]) # need to find proper threshold. Change to match threshold above 
    temp.i <- as.igraph(qgraph(temp_bin))
    decomposed_comps <-  decompose.graph(temp.i)
    max.index <- which.max(clusters(temp.i)$csize)
    ifelse(length(max.index)!=0, max.comp.list[k] <- gsize(decomposed_comps[[max.index]]), max.comp.list[k] <- 0)
    print(k)
  }
  write.table(max.comp.list, paste("max.comp.list_",alpha[a],"_", split.num, sep = ""), col.names = F, row.names = F)
}




#===================================================================================
#The above (NBS computation) should be probably run on the cluster if bs/permutations are >999
#====================================================================================


#-------------------------------------------------
# concatnate p.vect, t.vect.max, max.comp.list_$alpha, 
# if split and run on cluster
#-------------------------------------------------




#----------------------------------

##Calculate largest comp in observerd data
temp <- matrix(nrow=316, ncol=316) #next 3 lines of code just converting the pvector back into a con mat
temp[lower.tri(temp)] <- as.matrix(1-p.vect)
temp <- forceSymmetric(temp, uplo = "L")
diag(temp) <- 1

temp_bin<- binarize(temp, threshold=3.6) # this is what to change when looking at sig comps


temp.i <- as.igraph(qgraph(temp_bin)) #maybe need to add as.matrix

# splits into comps (finds largest comp and gives size edges)
decomposed_comps <-  decompose.graph(temp.i)
max.index <- which.max(clusters(temp.i)$csize)
ifelse(length(max.index)!=0, max.comp <- gsize(decomposed_comps[[max.index]]), max.comp <- 0)




#0.001 #18
thr_3.2_t #<- max.comp.list
hist(thr_3.2_t)
quantile(thr_3.32_t, c(.95)) 



#0.01 #max_comp #13
thr_2.60 <- #max.comp.list
  hist(thr_2.60)
quantile(thr_2.60, c(.95)) 

#0.05 #max_comp 146
thr_1.97 <- max.comp.list
hist(thr_1.97)
quantile(thr_1.97, c(.95)) 


####edge-wise FC
hist(r.vect.max)
quantile(r.vect.max, c(.95)) 
which(r.vect > 0.242)
max(r.vect)


#---------------------------------------
# below has been writtin for aspreeFC project (spearmans rho)
# Need to be adapted for t-value rather than rho
#---------------------------------------


###Now that we have a comp (p<0.05 FWE - (cft p<0.001))
# What is the effect? 

#search a index to find the vector/edge  number
index_vect <- vector()
index_vect <- 1:3321
t.mat <- matrix(nrow=82, ncol=82)
labs <- read_xlsx("labs.xlsx") 
rownames(t.mat) <- labs$Label
colnames(t.mat) <- labs$Label
t.mat[lower.tri(t.mat)] <- index_vect #mod for 1-p, log for -log(10)
t.mat <- forceSymmetric(t.mat, uplo = L)
diag(t.mat) <- 1
t.mat[8,5]



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

qgraph(p.mat.bin.maxcomp)
cluster_sig_edges <- which(p.mat.bin.maxcomp > 0)

#write out largest cluster
write.csv(x = as.matrix(p.mat.bin.maxcomp), file = "test_clust_a1.txt", quote = FALSE, row.names = FALSE, col.names = FALSE  )

x <- as.matrix(read.csv("test_clust_a1.txt"))






##visualise using ggconnectome
source("~/Dropbox/Sid/R_files/ggconnectome/functions/add_geoms.R")
source("~/Dropbox/Sid/R_files/ggconnectome/functions/ggConnectome.R")
source("~/Dropbox/Sid/R_files/ggconnectome/functions/ggConnectome3D.R")

ggConnectome(atlas = "Stages_melbBrain", conmat = x)
ggConnectome3D(atlas = "", conmat = x)
