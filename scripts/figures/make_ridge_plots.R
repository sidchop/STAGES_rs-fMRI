##Make ridge plot
library(dplyr)
library(forcats)
library(ggridges)

net_names <- c("Visual", "Somatomotor", "Dorsal Attention", 
               "Ventral Attention", "Limbic", "Frontoparietal", 
               "Default", "Amyg & Hippo", "Thalamus", "Striatum")
network_loadings <- list()

loading_matrix <- as.matrix(read.table("~/Dropbox/Sid/R_files/STAGES_fmri/data/cca/a2/cor_to_brain_pos.txt", header = F))



network_loadings[[1]] <- rbind(loading_matrix[1:24,], loading_matrix[143:165,])
network_loadings[[2]] <- rbind(loading_matrix[25:53,], loading_matrix[166:193,])
network_loadings[[3]] <- rbind(loading_matrix[54:67,], loading_matrix[194:211,])
network_loadings[[4]] <- rbind(loading_matrix[70:85,], loading_matrix[212:229])
network_loadings[[5]] <- rbind(loading_matrix[86:88,], loading_matrix[230:231])
network_loadings[[6]] <- rbind(loading_matrix[89:105,], loading_matrix[232:254])
network_loadings[[7]] <- rbind(loading_matrix[106:142,], loading_matrix[255:284])
network_loadings[[8]] <- rbind(loading_matrix[285:288,], loading_matrix[301:304])
network_loadings[[9]] <- rbind(loading_matrix[289:292,], loading_matrix[305:308])
network_loadings[[10]] <- rbind(loading_matrix[293:300,], loading_matrix[309:316])

#make ridgeplot dataframe
network_loadings_ridges <- list()
for (n in 1:length(network_loadings)) {
  network_loadings_ridges[[n]] <- cbind(c(network_loadings[[n]]), rep(net_names[n], length(c(network_loadings[[n]]))))
}

loadings_ridges_vec <- rbind(network_loadings_ridges[[1]], network_loadings_ridges[[2]],network_loadings_ridges[[3]],network_loadings_ridges[[4]],
                             network_loadings_ridges[[5]],network_loadings_ridges[[6]],network_loadings_ridges[[7]],network_loadings_ridges[[8]],
                             network_loadings_ridges[[9]],network_loadings_ridges[[10]])

#remove zeros and split by sign
#== choose neg or pos or both for loadings
loadings_ridges_vec_pos <- as.data.frame(loadings_ridges_vec[which(loadings_ridges_vec[,1]!=0),])



colnames(loadings_ridges_vec_pos) <- c("loading", "network")
loadings_ridges_vec_pos$network <- as.factor(loadings_ridges_vec_pos$network)
loadings_ridges_vec_pos$loading <- as.numeric(as.character(loadings_ridges_vec_pos$loading))



#quantiles <- quantile(loadings_ridges_vec_pos$loading, c(.10, .20, .50, .80, .90, .95))
#loadings_ridges_vec_pos_thresh <- loadings_ridges_vec_pos[which(loadings_ridges_vec_pos$loading<quantiles[1]),]


#for neg corrlations 
ridge_plot <- loadings_ridges_vec_pos %>%
  mutate(network = fct_reorder(network, loading, .fun='mean')) %>%
  ggplot(aes(y=reorder(network, loading), x=loading)) + 
  geom_density_ridges(aes(fill=network), alpha=0.7) +
  theme_minimal() + 
  theme(legend.position = "none", axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17),
        axis.title=element_text(size=18, vjust = 2)) +
  scale_fill_manual(values = c("Visual" = "#9C755FFF", "Somatomotor" = "#F28E2BFF",
                               "Dorsal Attention"= "#E15759FF","Ventral Attention" = "#76B7B2FF", 
                               "Limbic" = "#59A14FFF", "Frontoparietal" = "#EDC948FF", 
                               "Default" = "#B07AA1FF", "Amyg & Hippo" = "#5575b3FF", 
                               "Thalamus" = "#496499FF", "Striatum" = "#4E79A7FF")) + 
  scale_x_reverse() + 
  xlab("Correlation to Brain Variate") + 
  ylab("") + xlim(c(-0.1,0.8)) 
ridge_plot
setwd("~/Dropbox/Sid/R_files/STAGES_fmri/output/figures")
ggsave(plot = ridge_plot, 
       filename = "clin_pos_ridge_plot.tiff", 
       device = "tiff")



#for pos corrlations 
# ridge_plot <- loadings_ridges_vec_pos %>%
#  mutate(network = fct_reorder(network, loading, .fun='mean')) %>%
#  ggplot(aes(y=reorder(network, desc(loading)), x=loading)) + 
#  geom_density_ridges(aes(fill=network), alpha=0.7) +
#  theme_minimal() + 
#  theme(legend.position = "none", axis.text.x = element_text(size = 15),
#        axis.text.y = element_text(size = 15),
#        axis.title=element_text(size=18, vjust = 3)) +
#  scale_fill_manual(values = c("Visual" = "#9C755FFF", "Somatomotor" = "#F28E2BFF",
#                               "Dorsal Attention"= "#E15759FF","Ventral Attention" = "#76B7B2FF", 
#                               "Limbic" = "#59A14FFF", "Frontoparietal" = "#EDC948FF", 
#                               "Default" = "#B07AA1FF", "Amyg & Hippo" = "#5575b3FF", 
#                               "Thalamus" = "#496499FF", "Striatum" = "#4E79A7FF")) + 
#  xlab("Correlation to Brain Variate") + 
#  ylab("")
#ridge_plot
#
#setwd("~/Dropbox/Sid/R_files/STAGES_fmri/output/figures")
#ggsave(plot = ridge_plot, filename = "dose_pos_redge_plot.tiff", device = "tiff")
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
###########################################
##load in top 5 pct matrix
#
#loading_matrix_5pct <- as.matrix(read.table("~/Dropbox/Sid/R_files/STAGES_fmri/data/cca/a2/cor_to_brain_pos_5_pct.txt"))
#source("~/Dropbox/Sid/R_files/functions/plotClassifiedEdges.R")
##split subcortex in to mtl, stri, thal
#labels <- c("Visual", "Somatomotor", "Dorsal Attention", "Ventral Attention" , "Limbic",
#            "Frontoparietal", "Default", "Amyg & Hippo", "Thalamus", "Striatum")
#networks <- read.table("~/Dropbox/Sid/R_files/STAGES_fmri/output/NetworkPlots/labels_mtl_str_thal.txt")
#
#matrix <- as.matrix(loading_matrix_5pct)
#matrix[matrix != 0] <- 1
###run pce
#pce <- plotClassifiedEdges(adj = matrix, ids = networks, labels = labels)
#
#
#library(reshape2)
#library(scales)
#matrix1 <- pce[[3]]
#
#matrix1 <- round(matrix1*10, 3)
##reorder to match ridgeplot ordering 
#ranks <- 11-(rank(tapply(loadings_ridges_vec_pos$loading, loadings_ridges_vec_pos$network, mean)))
#ranks_sorted <- rev(sort(ranks))
#matrix1_reordered <- matrix1[c(names(ranks_sorted )), c(names(ranks_sorted ))]
#
## Get lower triangle of the correlation matrix
#get_upper_tri<-function(cormat){
#  cormat[lower.tri(cormat)] <- NA
#  return(cormat)
#}
#
#
#lower_tri <- get_upper_tri(matrix1_reordered) 
#library(reshape2)
#melted_cormat <- melt(lower_tri, na.rm = TRUE)
#
#head(melted_cormat)
#
#
#ggheatmap <- ggplot(data = melted_cormat, aes(x=Var2, y=Var1, fill=value)) + 
#  geom_tile(color = "black")  +
#  scale_y_discrete(position = "right") + 
#  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
#                                   size = 18, hjust = 1), 
#        axis.text.y = element_text(vjust = 0.5, 
#                                   size = 18, hjust = 2),
#        axis.title.x = element_blank(),
#        axis.title.y = element_blank(),
#        panel.border = element_blank(),
#        axis.ticks = element_blank(),
#        legend.text=element_text(size=14),
#        legend.title = element_text(size=15),
#        legend.position = c(0.3, 0.8),
#        plot.margin=grid::unit(c(0,0,0,0), "mm"),
#        panel.grid.major = element_blank(), 
#        panel.grid.minor = element_blank(),
#        panel.background = element_rect(fill = "transparent",colour = NA),
#        plot.background = element_rect(fill = "transparent",colour = NA)) +
#  coord_fixed() + 
#  scale_fill_gradientn(colours=c("#F5F5F5","light yellow", "yellow", "orange", "red", "dark red")) +
#  labs(fill = "Normalised\n proportion")
#
#
#ggsave(plot = ggheatmap, filename = "clin_pos_matrix_5pct.tiff", device = "tiff")
