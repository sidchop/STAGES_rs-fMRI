##make circle plots
####### Visualisng networks in STAGES fmri
library(igraph)
library(qgraph)
library(reshape2)
library(neurobase)
library(circlize)
library(tidyverse)
library(viridis)
library(patchwork)
library(hrbrthemes)
library(circlize)


source("~/Dropbox/Sid/R_files/functions/plotClassifiedEdges.R")
source("~/Dropbox/Sid/R_files/functions/vec_2_mat.R")
load("~/Dropbox/Sid/R_files/brainconn/data/stages_study.rda") 
setwd("~/Dropbox/Sid/R_files/STAGES_fmri/output/NetworkPlots")





#load all comps as a list ----
comps <- list()
comps[[1]] <- read.table("t_pos_a1_baseline.txt")
comps[[2]] <- read.table("t_neg_a1_baseline.txt")
comps[[3]] <- read.table("t_pos_a2_illness.txt")
comps[[4]] <- read.table("t_neg_a2_illness.txt")
comps[[5]] <- read.table("t_pos_a2_med.txt")
comps[[6]] <- read.table("t_neg_a2_med.txt")
comps[[7]] <- read.table("t_pos_a3_illness.txt")
comps[[8]] <- read.table("t_neg_a3_illness.txt")
comps[[9]] <- read.table("t_pos_a3_med.txt")
comps[[10]] <- read.table("t_neg_a3_med.txt")



network_plot_titles <- c("Patients lower than controls - baseline", 
                         "Patients higher than controls - baseline",
                         "PIPT decreasing compared to controls - baseline to 3m", 
                         "PIPT increasing compared to controls - baseline to 3m",
                         "MIPT increasing compared to controls & PIPT - baseline to 3m", 
                         "MIPT decreasing compared to controls & PIPT - baseline to 3m", 
                         "PIPT decreasing compared to controls - baseline to 12m - p<0.001 effect only",
                         "PIPT increasing compared to controls - baseline to 12m - p<0.001 effect only",
                         "MIPT increasing compared to controls & PIPT - baseline to 12m", 
                         "MIPT decreasing compared to controls & PIPT - baseline to 12m")

comps_t_stat <- list()
#load in tmaps if making circple plots weighted by t-value
setwd("~/Dropbox/Sid/R_files/STAGES_fmri/data/swe_validation_contrasts/gmr/fwe_contrast_corrected/t_stats/")
tstat_baseline <- c(readNIfTI2("t_a1_baseline.nii"))
tstat_baseline <- tstat_baseline[!is.na(tstat_baseline)]
tmat_baseline <- vec_2_mat(tstat_baseline, 316, 0)

tstat_illness_3m <- c(readNIfTI2("t_a2_illness.nii"))
tstat_illness_3m <- tstat_illness_3m[!is.na(tstat_illness_3m)]
tmat_illness_3m <- vec_2_mat(tstat_illness_3m, 316, 0)

tstat_illness_12m <- c(readNIfTI2("t_a3_illness.nii"))
tstat_illness_12m <- tstat_illness_12m[!is.na(tstat_illness_12m)]
tmat_illness_12m <- vec_2_mat(tstat_illness_12m, 316, 0)
#medication effect tscors are the maximum of two t-scores, already made using the visualise t-score script
tmat_med_3m <- as.matrix(read.csv("med_3m_max_t.csv", header = F))
tmat_med_12m <- as.matrix(read.csv("med_12m_max_t.csv", header = F))

comps_t_stat[[1]] <- tmat_baseline 
comps_t_stat[[2]] <- tmat_baseline 
comps_t_stat[[3]] <- tstat_illness_3m 
comps_t_stat[[4]] <- tstat_illness_3m 
comps_t_stat[[5]] <- tmat_med_3m
comps_t_stat[[6]] <- tmat_med_3m
comps_t_stat[[7]] <- tmat_illness_12m
comps_t_stat[[8]] <- tmat_illness_12m
comps_t_stat[[9]] <- tmat_med_12m
comps_t_stat[[10]] <- tmat_med_12m


setwd("~/Dropbox/Sid/R_files/STAGES_fmri/output/NetworkPlots")

#==================
####circle plot
#====================


##schafer to dk - optional conversion
#labels <- read.csv("dk_short_labs", header = F)
#networks <- read.csv("labels_dk.txt", header = F)
## rsn only 
#labels <- c("Visual", "Somatomotor", "Dorsal Attention", "Ventral Attention" , "Limbic",
#            "Frontoparietal", "Default", "Amyg & Hippo", "Thalamus", "Striatum")
#networks <- read.table("labels_mtl_str_thal.txt")
#===========
#run pce
#===========
#thr=4 #thresholding 
thresh_level <- 1-c(.99, .95, .90, .80)
setwd("~/Dropbox/Sid/R_files/STAGES_fmri/output/circle_plots")
for (c in 1:length(comps)) {
  for (t in 1:4) { #this loop is making cps at diff thresh at thr roi level (ie when using tscores for weights)
    
    #get tscores
    comps_thr <- as.matrix(comps[[c]]*comps_t_stat[[c]])
    #threshold to top 10% of scores
    t_thr <- quantile(comps_thr[comps_thr != 0], c(.99, .95, .90, .80))
    
    comps_thr[comps_thr < t_thr[t]] <- 0.1 #& comps_thr > 0] <- 0.1#switch sign 

    
    
    
    #pce <- plotClassifiedEdges(adj = comps[[c]], ids = networks, labels = labels)
    
    
    #network_plot_list[[c]] <- makeNetworkMatrix(pce[[1]])
    
    #==========
    # reorder back into resting state networks 
    #========
    #cols <- t((colnames(pce[[1]])))
    #write.table(cols, "labels_aal_to_sch.txt", col.names = F, row.names = F, quote = F)
    #pce_thr <- pce[[3]]
    pce_thr <- comps_thr 
    
    #load in short names
  short_labs <- read.table("stages_schaefer_short_labs", header = T)
    #colnames(pce_thr) <- rownames(pce_thr) <- stages_study$ROI.Name
  colnames(pce_thr) <- rownames(pce_thr) <-   t(short_labs)
  
  #pce_thr[pce_thr==0] <- 0.5 #only do this for sparce graphs 
    ##reorder to networks and ant -> post. 
    order.file <- read.table("~/Dropbox/Sid/R_files/STAGES_fmri/output/circle_plots/stages_reorderd_by_network.txt",
                             header = T)
    
#reorder by the feature, then use the orig_ordering
order <- order.file[order(order.file$network_ordering),]
    
pce_thr_reordered <-  pce_thr[c(order$orig_ordering), c(order$orig_ordering)]

    
    
    #reorder matrix based on networks
    #ordering <- c(rep(2, 10), rep(6,20), rep(3, 5), rep(4,5))
    #pce_thr_reordered <- graph4lg::reorder_mat(pce_thr, ordering)
    #quantiles <- quantile(pce_thr, c(.10, .20, .50, .80, .90, .95))
    #pce_thr[abs(pce_thr) < quantiles[4]] <- 0
    
    # Make the circular plot
    #make thin lines faded like : https://medium.com/@akkalbist55/multiple-group-chord-diagram-5e25666cb0d
    #col_mat = rand_color(length(pce_thr), transparency = 0.5)
    #
    
    ###below is for dk atlas onlt
  # col_mat = rand_color(length(pce_thr), transparency = 0.5)
  # dim(col_mat) = dim(pce_thr)
  # col_mat[1:47, ] <-"#9C755F80"
  # col_mat[,1:47] <- "#9C755F80"
  # col_mat[48:104, ] <-"#F28E2B80"
  # col_mat[,48:104] <- "#F28E2B80"
  # col_mat[105:138, ] <-"#E1575980"
  # col_mat[,105:138] <- "#E1575980"
  # col_mat[139:172, ] <-"#76B7B280"
  # col_mat[,139:172] <- "#76B7B280"
  # col_mat[173:177, ] <-"#59A14F80"
  # col_mat[,173:177] <- "#59A14F80"
  # col_mat[178:217, ] <-"#EDC94880"
  # col_mat[,178:217] <- "#EDC94880"
  # col_mat[218:284, ] <-"#B07AA180"
  # col_mat[,218:284] <- "#B07AA180"
  # col_mat[285:316, ] <-"#4E79A780"
  # col_mat[,285:316] <- "#4E79A780"
    #
   # col_mat[pce_thr < 3] = "#00000000" #leace links but make transparent DONT DO THIS FOR SPARCE GRAPHS
    #col_mat[pce_thr < 1] = "#00000000" #leace links but make transparent
    #dim(col_mat) = dim(pce_thr)
    ### above is for dk atlas only
    
    #grid.col = c(Visual = "#9C755FFF", Somatomotor = "#F28E2BFF", "Dorsal Attention" = "#E15759FF",
    #             "Ventral Attention"= "#76B7B2FF" , "Limbic" = "#59A14FFF"
    #             , "Frontoparietal" = "#EDC948FF", "Default" = "#B07AA1FF", 
    #             "Amyg & Hippo" = "#4E79A7FF", "Thalamus" = "#BAB0ACFF", "Striatum" = "#FF9DA7FF")
   grid.col = c(rep("#9C755FFF", 47), rep("#F28E2BFF", 57), rep("#E15759FF",34),
                rep("#76B7B2FF",34), rep("#59A14FFF",5), rep("#EDC948FF",40),
                rep("#B07AA1FF",67),rep("#4E79A7FF",32))
    
    circos.clear()
    
    #if(c==5) {
    #  circos.par(gap.after = c(rep(2,26),15,rep(2,6),15,rep(2,6)), start.degree = 165) 
    #}
    #if(c==6) {
    #  circos.par(gap.after = c(rep(2,26),15,rep(2,6),15,rep(2,6)), start.degree = 165, clock.wise = FALSE) 
    #}
 
#png(paste0("cp2_roi_", c, "thresh_", thresh_level[t],".png"),res=300, bg=NA, width = 1500, height = 1500)

    
    
     circos.par(gap.after = rep(0,316))

      chordDiagram(pce_thr_reordered, symmetric = F, annotationTrack = "grid",
 preAllocateTracks = list(track.height = uh(7, 'mm'),
                                                        track.height = uh(10, 'mm')), 
                        #col = grid.col, 
                 grid.col = as.vector(grid.col),self.link = 0, 
 link.visible = pce_thr_reordered > 0.1)

 
    
    
    
    
    circos.track(track.index = 1, panel.fun = function(x, y) {
      circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                  facing = "clockwise", niceFacing = TRUE, adj = c(0, 1), cex = 0.3)
    }, bg.border = NA)
    

    
    
    #only if using dk atlas
    
    #highlight.sector(row.names(pce_thr)[1:6], 
    #                 track.index = 2, 
    #                 col = "#9C755FFF", 
    #  #                text = 'Visual', 
    #                  cex = 0.8, 
    #                 text.col = 'white', niceFacing = TRUE)
    #
    #highlight.sector(row.names(pce_thr)[7:11], 
    #                 track.index = 2, 
    #                 col = "#F28E2BFF", 
    ##                 text = 'Somatomotor', 
    #                 cex = 0.8, 
    #                 text.col = 'white', niceFacing = TRUE)
    #highlight.sector(row.names(pce_thr)[12:14], 
    #                 track.index = 2, 
    #                 col = "#E15759FF", 
    #  #               text = 'Dorsal Attention', 
    #                 cex = 0.8, 
    #                 text.col = 'white', niceFacing = TRUE)
    #highlight.sector(row.names(pce_thr)[15:18], 
    #                 track.index = 2, 
    #                 col = "#76B7B2FF", 
    #      #           text = 'Ventral Attention', 
    #                 cex = 0.8, 
    #                 text.col = 'white', niceFacing = TRUE)
    #highlight.sector(row.names(pce_thr)[19:20], 
    #                 track.index = 2, 
    #                 col = "#59A14FFF", 
    #       #          text = 'Limbic', 
    #                 cex = 0.8, 
    #                 text.col = 'white', niceFacing = TRUE)
    #highlight.sector(row.names(pce_thr)[21:24], 
    #                 track.index = 2, 
    #                 col = '#EDC948FF', 
    #    #             text = 'Frontoparietal', 
    #                 cex = 0.8, 
    #                 text.col = 'white', niceFacing = TRUE)
    #highlight.sector(row.names(pce_thr)[25:33], 
    #                 track.index = 2, 
    #                 col = '#B07AA1FF', 
    #       #          text = 'Default', 
    #                 cex = 0.8, 
    #                 text.col = 'white', niceFacing = TRUE)
    #highlight.sector(row.names(pce_thr)[34:40], 
    #                 track.index = 2, 
    #                 col = "#4E79A7FF", 
    #            #     text = 'Subcortical', 
    #                 cex = 0.8, 
    #                 text.col = 'white', niceFacing = TRUE)
    #
    
    
   # title(network_plot_titles[c], cex = 0.8)
    
    
    dev.off()
    
  }
}


 

