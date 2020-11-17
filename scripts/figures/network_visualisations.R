
####### Visualisng networks in STAGES fmri
library(igraph)
library(qgraph)
library(reshape2)
library(ggplot2)

#plot classified edges
source("~/Dropbox/Sid/R_files/functions/plotClassifiedEdges.R")
source("~/Dropbox/Sid/R_files/STAGES_fmri/scripts/makeNetworkMatrix2.R")
load("~/Dropbox/Sid/R_files/brainconn/data/stages_study.rda") #add colnames
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

#cca illness primary outcomes 3m
comps[[11]] <- read.table("~/Dropbox/Sid/R_files/STAGES_fmri/data/cca/a2/cor_to_brain_FDR_pos.txt")
comps[[12]] <- read.table("~/Dropbox/Sid/R_files/STAGES_fmri/data/cca/a2/cor_to_brain_FDR_neg.txt")

#cca med exploratory outcomes 12m
comps[[13]] <- read.table("~/Dropbox/Sid/R_files/STAGES_fmri/data/cca/a3/cor_to_brain_FDR_pos.txt")
comps[[14]] <- read.table("~/Dropbox/Sid/R_files/STAGES_fmri/data/cca/a3/cor_to_brain_FDR_neg.txt")


## Sup results
#baseline 
comps[[15]] <- read.table("~/Dropbox/Sid/R_files/STAGES_fmri/output/NetworkPlots/t_pos_a1_baseline_01.txt")
comps[[16]] <- read.table("~/Dropbox/Sid/R_files/STAGES_fmri/output/NetworkPlots/t_neg_a1_baseline_01.txt")
comps[[17]] <- read.table("~/Dropbox/Sid/R_files/STAGES_fmri/output/NetworkPlots/t_pos_a1_baseline_001.txt")
comps[[18]] <- read.table("~/Dropbox/Sid/R_files/STAGES_fmri/output/NetworkPlots/t_neg_a1_baseline_001.txt")

#illness 3m
comps[[19]] <- read.table("~/Dropbox/Sid/R_files/STAGES_fmri/output/NetworkPlots/t_pos_a2_illness_01.txt")
comps[[20]] <- read.table("~/Dropbox/Sid/R_files/STAGES_fmri/output/NetworkPlots/t_neg_a2_illness_01.txt")

comps[[21]] <- read.table("~/Dropbox/Sid/R_files/STAGES_fmri/output/NetworkPlots/t_pos_a2_illness_001.txt")
comps[[22]] <- read.table("~/Dropbox/Sid/R_files/STAGES_fmri/output/NetworkPlots/t_neg_a2_illness_001.txt")

#med 3m
comps[[23]] <- read.table("~/Dropbox/Sid/R_files/STAGES_fmri/output/NetworkPlots/t_pos_a2_med_01.txt")
comps[[24]] <- read.table("~/Dropbox/Sid/R_files/STAGES_fmri/output/NetworkPlots/t_neg_a2_med_01.txt")

comps[[25]] <- read.table("~/Dropbox/Sid/R_files/STAGES_fmri/output/NetworkPlots/t_pos_a2_med_001.txt")
comps[[26]] <- read.table("~/Dropbox/Sid/R_files/STAGES_fmri/output/NetworkPlots/t_neg_a2_med_001.txt")

#med 12m
comps[[27]] <- read.table("~/Dropbox/Sid/R_files/STAGES_fmri/output/NetworkPlots/t_pos_a3_med_01.txt")
comps[[28]] <- read.table("~/Dropbox/Sid/R_files/STAGES_fmri/output/NetworkPlots/t_neg_a3_med_01.txt")

comps[[29]] <- read.table("~/Dropbox/Sid/R_files/STAGES_fmri/output/NetworkPlots/t_pos_a3_med_001.txt")
comps[[30]] <- read.table("~/Dropbox/Sid/R_files/STAGES_fmri/output/NetworkPlots/t_neg_a3_med_001.txt")


#split subcortex in to mtl, stri, thal
labels <- c("Visual", "Somatomotor", "Dorsal Attention", "Ventral Attention" , "Limbic",
                        "Frontoparietal", "Default", "Amyg & Hippo", "Thalamus", "Striatum")


networks <- read.table("labels_mtl_str_thal.txt")




#split subcortex in to stri, thal and class hipp/amy as limbic
#labels <- c("Visual", "Somatomotor", "Dorsal Attention", "Ventral Attention" , "Limbic",
#            "Frontoparietal", "Default", "Thalamus", "Striatum")
#networks <- read.table("labels_str_thal.txt")


#network_plot_titles <- c("Patients lower than controls - baseline", 
#                       "Patients higher than controls - baseline",
#                       "PIPT decreasing compared to controls - baseline to 3m", 
#                       "PIPT increasing compared to controls - baseline to 3m",
#                       "MIPT increasing compared to controls & PIPT - baseline to 3m", 
#                       "MIPT decreasing compared to controls & PIPT - baseline to 3m", 
#                       "PIPT decreasing compared to controls - baseline to 12m - p<0.001 effect only",
#                       "PIPT increasing compared to controls - baseline to 12m - p<0.001 effect only",
#                       "MIPT increasing compared to controls & PIPT - baseline to 12m", 
#                       "MIPT decreasing compared to controls & PIPT - baseline to 12m")

#max and min for colourbars

network_plot_list_a <- list()
network_plot_list_b <- list()
for (c in 1:(length(comps)/2)) {

  if(c==1) {comp1 <- 1
  comp2 <- 2}
  if(c==2) {comp1 <- 3
  comp2 <- 4}
  if(c==3) {comp1 <- 5
  comp2 <- 6}
  if(c==4) {comp1 <- 7
  comp2 <- 8}
  if(c==5) {comp1 <- 9
  comp2 <- 10}
  if(c==6) {comp1 <- 11
  comp2 <- 12}
  if(c==7) {comp1 <- 13
  comp2 <- 14}
  if(c==8) {comp1 <- 15
  comp2 <- 16}
  if(c==9) {comp1 <- 17
  comp2 <- 18}
  if(c==10) {comp1 <- 19
  comp2 <- 20}
  if(c==11) {comp1 <- 21
  comp2 <- 22}
  if(c==12) {comp1 <- 23
  comp2 <- 24}
  if(c==13) {comp1 <- 25
  comp2 <- 26}
  if(c==14) {comp1 <- 27
  comp2 <- 28}
  if(c==15) {comp1 <- 29
  comp2 <- 30}
  
  #make input into pce binary
  matrix1 <- as.matrix(comps[[comp1]])
  matrix1[matrix1 != 0] <- 1
  
  matrix2 <- as.matrix(comps[[comp2]])
  matrix2[matrix2 != 0] <- 1
  
  ##run pce
  pce1 <- plotClassifiedEdges(adj =   matrix1, ids = networks, labels = labels)
  pce2 <- plotClassifiedEdges(adj =   matrix2, ids = networks, labels = labels)
  
  #round
  pce1[[3]] <- round(pce1[[3]]*10, 3)
  pce2[[3]] <- round(pce2[[3]]*10, 3) #make it a more interp % and round
  pce1[[1]]<- round(pce1[[1]], 3) #round
  pce2[[1]]<- round(pce2[[1]], 3) 
  
  #add short labs 
  short_labs <- c("Vis", "SomMot", "DorsAttn", 
                  "VentAttn", "Lim", "FPN", "DMN", 
                  "MTL", "Thal", "Stri")

  
  rownames(pce1[[1]]) <- rownames(pce1[[3]]) <-   rownames(pce2[[1]]) <- rownames(pce2[[3]]) <-   short_labs 
 colnames(pce1[[1]]) <- colnames(pce1[[3]]) <-   colnames(pce2[[1]]) <- colnames(pce2[[3]]) <-   short_labs 
  
  #set min and max of colour bar
  ifelse(min(pce1[[3]]) < min(pce2[[3]]), min1 <- min(pce1[[3]]), min1 <- min(pce2[[3]]))
  ifelse(min(pce1[[1]]) < min(pce2[[1]]), min2 <- min(pce1[[1]]), min2 <- min(pce2[[1]]))
  
  ifelse(max(pce1[[3]]) > max(pce2[[3]]), max1 <- max(pce1[[3]]), max1 <- max(pce2[[3]]))
  ifelse(max(pce1[[1]]) > max(pce2[[1]]), max2 <- max(pce1[[1]]), max2 <- max(pce2[[1]]))
  #network matrix =====================
  network_plot_list_a[[c]] <- makeNetworkMatrix2(matrix1 = pce1[[3]], matrix2 = pce1[[1]], 
                                               min1 = min1, max1 = max1, 
                                               min2 = min2, max2 = max2,
                                               pal = "Reds") #multiply by 10 for pce 3
  
  network_plot_list_b[[c]] <- makeNetworkMatrix2(matrix1 = pce2[[3]], matrix2 = pce2[[1]], 
                                               min1 = min1, max1 = max1, min2 = min2,max2= max2,
                                               pal = "Reds")
  #make cirecle plot
  #png(paste0("network2_pce3_pce1_", c,".png"), res = 1200)
network_plot_list_a[[c]] + theme(panel.grid.major = element_blank(), 
                               panel.grid.minor = element_blank(),
                               panel.background = element_rect(fill = "transparent",colour = NA),
                               plot.background = element_rect(fill = "transparent",colour = NA))
ggsave(filename = paste0("network3_pce3_pce1_", comp1,".png"), 
       plot = last_plot(), bg = "transparent")

network_plot_list_b[[c]] + theme( panel.grid.major = element_blank(), 
                                 panel.grid.minor = element_blank(),
                                 panel.background = element_rect(fill = "transparent",colour = NA),
                                 plot.background = element_rect(fill = "transparent",colour = NA))


ggsave(filename = paste0("network3_pce3_pce1_", comp2,".png"), 
         plot = last_plot(), bg = "transparent")
   # dev.off()
  }



##illness related plots
#patients lower
#cowplot::plot_grid(network_plot_list[[1]], network_plot_list[[3]], network_plot_list[[7]],
#                   rows = 1, cols = 3)
#patients higher
#cowplot::plot_grid(network_plot_list[[2]], network_plot_list[[4]], network_plot_list[[8]],
#                   rows = 1, cols = 3)

#medication effect
#medication increasing
#cowplot::plot_grid(network_plot_list[[5]], network_plot_list[[9]],
#                   rows = 1, cols = 2)
#medication decreasing
#cowplot::plot_grid(network_plot_list[[6]], network_plot_list[[10]],
#                   rows = 1, cols = 2)

















###=============hitogram plot of degree
#load("~/Dropbox/Sid/R_files/brainconn/data/stages_study.rda") #add colnames
##rownames(pc_degree) <- stages_study$ROI.Name
##View(pc_degree)
#pc_degree <- as.data.frame(pc_degree)
#pc_degree_bind <- cbind(pc_degree, as.data.frame(stages_study$ROI.Name))
#pc_degree_order <- pc_degree_bind[order(pc_degree_bind$`igraph::degree(igraph::graph_from_adjacency_matrix(as.matrix(comp), mode = c("undirected"), weighted = NULL))`) , ]
##pc_degree <- sort(pc_degree$`unlist(pc_degree)`, decreasing = T)
#
#colnames(pc_degree_order) <- c("degree", "roi")
#pc_degree_order$roi <- factor(pc_degree_order$roi, levels = pc_degree_order$roi[order(pc_degree_order$degree)])
#pc_degree_order <- pc_degree_order[250:316,]
#library(viridis)
#ggplot(pc_degree_order,aes(roi, degree))+geom_bar(stat="identity", position = position_dodge(width = .5)) + 
#  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
#                                   size = 7, hjust = 1)) +theme_minimal()
#
#
#
#
#
#
#
#
#
#
#pheatmap::pheatmap(pce[[3]], treeheight_row = 0, treeheight_col = 0, 
#                   cluster_rows=FALSE, cluster_cols=FALSE, display_numbers = round(pce[[1]], 3))
#
#