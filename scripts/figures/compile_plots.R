##Compile figures 
library(ggplot2)
library(cowplot)
library(magick)
library(cowplot)
library(patchwork)


#illness 
idi <- c(1,3,7,2,4,8)
#medication
mdi <- c(6,10, 5,9)

##compile plots main plots
panel_list <- list()
for (p in 1:length(idi)) {
  setwd("~/Dropbox/Sid/R_files/STAGES_fmri/output/brain_figs/")
  brain <- ggdraw() + draw_image(paste0(idi[p],".png"), 
                                 scale = 0.9) 
  
  setwd("~/Dropbox/Sid/R_files/STAGES_fmri/output/surfaces")
  surf_rl <- ggdraw() + draw_image(paste0("t_analysis_rh_",idi[p],"_inflated_transparent_reds_lat.jpg"), scale = 1) 
  surf_rm <- ggdraw() + draw_image(paste0("t_analysis_rh_",idi[p],"_inflated_transparent_reds_med.jpg"), scale = 1)
  surf_ll <- ggdraw() + draw_image(paste0("t_analysis_lh_",idi[p],"_inflated_transparent_reds_lat.jpg"), scale = 1)
  surf_lm <- ggdraw() + draw_image(paste0("t_analysis_lh_",idi[p],"_inflated_transparent_reds_med.jpg"), scale = 1)
  setwd("~/Dropbox/Sid/R_files/STAGES_fmri/output/surfaces/colourbars")
  cb <- ggdraw() + draw_image(paste0("colourbar_",idi[p],".tiff"))
  setwd("~/Dropbox/Sid/R_files/STAGES_fmri/output/subcortex/imgs/")
  surf_subl <- ggdraw() + draw_image(paste0("f_t_sub_",idi[p],"_l.png"))
  surf_subr <- ggdraw() + draw_image(paste0("f_t_sub_",idi[p],"_r.png"))
  
  setwd("~/Dropbox/Sid/R_files/STAGES_fmri/output/NetworkPlots/")
  net <- ggdraw() + draw_image(paste0("network3_pce3_pce1_", idi[p],".png"), scale = 1.1)
  
  panel_list[[p]] <- compile_figures(brain = brain, 
                                     networkplot = net, 
                                     lh_l = surf_ll, 
                                     lh_m = surf_lm, 
                                     rh_l =  surf_rl, 
                                     rh_m =  surf_rm,  
                                     subcortl = surf_subl,
                                     subcortr = surf_subr,
                                     colorbar = cb)
  
}

#ii <- cowplot::plot_grid(panel_list[[1]],  panel_list[[2]],  panel_list[[3]], 
#                         rows = 3, cols = 1)#, labels = "auto", label_size = 30)
#id <- cowplot::plot_grid(panel_list[[1]],  panel_list[[2]],  panel_list[[3]], 
#                   rows = 3, cols = 1)#, labels = "auto", label_size = 30)

#plot <- cowplot::plot_grid(panel_list[[1]],panel_list[[2]],
#                   panel_list[[3]],panel_list[[4]],
#                   panel_list[[5]],panel_list[[6]],rows = 3, cols = 2)


#library(patchwork)
#panel_list[[1]]+panel_list[[2]]+panel_list[[3]]

#Write out panels 
setwd("~/Dropbox/Sid/R_files/STAGES_fmri/output/figures/")
for (a in 1:length(idi)) {
  ggsave(panel_list[[a]], filename = paste0("fig1_panel_",idi[a],".tiff"), device = "tiff")
}







panel_list <- list()
for (p in 1:length(mdi)) {
  setwd("~/Dropbox/Sid/R_files/STAGES_fmri/output/brain_figs/")
  brain <- ggdraw() + draw_image(paste0(mdi[p],".png"), 
                                 scale = 0.9) 
  
  setwd("~/Dropbox/Sid/R_files/STAGES_fmri/output/surfaces")
  surf_rl <- ggdraw() + draw_image(paste0("t_analysis_rh_",mdi[p],"_inflated_transparent_reds_lat.jpg"), scale = 1)
  surf_rm <- ggdraw() + draw_image(paste0("t_analysis_rh_",mdi[p],"_inflated_transparent_reds_med.jpg"), scale = 1)
  surf_ll <- ggdraw() + draw_image(paste0("t_analysis_lh_",mdi[p],"_inflated_transparent_reds_lat.jpg"), scale = 1)
  surf_lm <- ggdraw() + draw_image(paste0("t_analysis_lh_",mdi[p],"_inflated_transparent_reds_med.jpg"), scale = 1)
  setwd("~/Dropbox/Sid/R_files/STAGES_fmri/output/surfaces/colourbars")
  cb <- ggdraw() + draw_image(paste0("colourbar_",mdi[p],".tiff"))
  setwd("~/Dropbox/Sid/R_files/STAGES_fmri/output/subcortex/imgs/")
  #flip
  
  surf_subl <- ggdraw() + draw_image(paste0("f_t_sub_",mdi[p],"_l.png"))
  surf_subr <- ggdraw() + draw_image(paste0("f_t_sub_",mdi[p],"_r.png"))
  
  setwd("~/Dropbox/Sid/R_files/STAGES_fmri/output/NetworkPlots/")
  net <- ggdraw() + draw_image(paste0("network3_pce3_pce1_", mdi[p],".png"), scale = 1.1)
  
  panel_list[[p]] <- compile_figures(brain = brain, 
                                     networkplot = net, 
                                     lh_l = surf_ll, 
                                     lh_m = surf_lm, 
                                     rh_l =  surf_rl, 
                                     rh_m =  surf_rm,  
                                     subcortl = surf_subl,
                                     subcortr = surf_subr,
                                     colorbar = cb)
  
}

#ii <- cowplot::plot_grid(panel_list[[1]],  panel_list[[2]],  panel_list[[3]], 
#                         rows = 3, cols = 1)#, labels = "auto", label_size = 30)
#id <- cowplot::plot_grid(panel_list[[1]],  panel_list[[2]],  panel_list[[3]], 
#                   rows = 3, cols = 1)#, labels = "auto", label_size = 30)

#plot <- cowplot::plot_grid(panel_list[[1]],panel_list[[2]],
#                   panel_list[[3]],panel_list[[4]],
#                   panel_list[[5]],panel_list[[6]],rows = 3, cols = 2)


#library(patchwork)
#panel_list[[1]]+panel_list[[2]]+panel_list[[3]]

#Write out panels 
setwd("~/Dropbox/Sid/R_files/STAGES_fmri/output/figures/")
for (a in 1:length(mdi)) {
  ggsave(panel_list[[a]], filename = paste0("fig2_panel_",mdi[a],".tiff"), device = "tiff")
}






###################################
##plot symptom correlation figure #
##################################

#if i also want to plot ridge plots, i need to make unthresholded maps 
library(cowplot)
library(ggplot2)
#illness a2
setwd("~/Dropbox/Sid/R_files/STAGES_fmri/data/cca/a2")
a1_mat <- rmatio::read.mat("a2_illness_1.mat")

## correlation for PC
pcs_raw <- read.csv('a2_illness_cca_pcs.csv')[,-1]
variates <- as.data.frame(cbind(a1_mat[["U"]][,1], a1_mat[["V"]][,1]))
cor_pc <- as.data.frame(cbind(cor(pcs_raw, variates[,1]), 1:dim(a1_mat[["A"]])[1]))

pc_cor_plot <- ggplot(data=cor_pc, aes(x=as.factor(V2), y=V1)) + 
  geom_bar(stat="identity", fill = "#CD5C5C") + xlab("\u0394 Brain Principal Components") + 
  ylab("Canonical Loading") + ylim(-1, 1) +  
  theme_classic() + theme(axis.title.x = element_text(size = 16), 
                          axis.title.y = element_text(size = 16),
                          axis.text.x = element_text(size = 10), 
                          axis.text.y = element_text(size = 10))

pc_cor_plot

#add astrix for sig loadings
#source("~/Dropbox/Sid/R_files/functions/cca_sig_test.R")
#p_brain_fdr <- cca_sig_test(variates[,1], pcs_raw, b = 5000, method = "1")
#label.df <- data.frame(V2 = c("5", "6", "11"),
#                       V1 = c(-0.5315, 0.406,   0.5175)) #added. .1 for space


#pc_cor_plot  <- pc_cor_plot   + geom_text(data = label.df, label = "*",
#                                          size=6)





## correlation for clin
clin_raw <- read.csv('a2_illness_cca_clin1.csv')
variates <- as.data.frame(cbind(a1_mat[["U"]][,1], a1_mat[["V"]][,1]))
cor_clin <- as.data.frame(cbind(cor(clin_raw, variates[,2]), 1:dim(a1_mat[["B"]])[1]))
cor_clin$V2 <- factor(cor_clin$V2, labels = c("SOFAS", "BPRS"))

clin_cor_plot <- ggplot(data=cor_clin, aes(x=as.factor(V2), y=V1)) + 
  geom_bar(stat="identity", fill = "#CD5C5C") + xlab("\u0394 Behavioural Scales") + 
  ylab("Canonical Loading") + 
  ylim(-1, 1) + theme_classic() + theme(axis.title.x = element_text(size = 16), 
                                        axis.title.y = element_text(size = 16),
                                        axis.text.x = element_text(size = 10), 
                                        axis.text.y = element_text(size = 10))

clin_cor_plot
#p_brain_fdr <- cca_sig_test(variates[,2], clin_raw, b = 5000, method = "1")
#label.df <- data.frame(V2 = c("SOFAS", "BPRS"),
#                      V1 = c(-1, 0.591)) #added. .1 for space


#clin_cor_plot <- clin_cor_plot  + geom_text(data = label.df, label = "*",
#                                          size=6)



#jamovi colours ""#56B4E9","#E69F00", "#999999"","#E69F00", "#999999"
cca_scatter <- as.data.frame(cbind(a1_mat[["U"]][,1], a1_mat[["V"]][,1]))
colnames(cca_scatter) <- c("\u0394 Brain Canonical Variate", "\u0394 Behaviour Canonical Variate")

scatter_plot <- ggplot(cca_scatter, aes(x=`Δ Brain Canonical Variate`, y=`Δ Behaviour Canonical Variate`)) +
  geom_point(size=4, color = "#B22222") +  geom_smooth(method='lm',linetype = 2, colour = "grey", fill = "#CD5C5C") +
#  geom_text(x=1.6, y=-1.6, label=paste0("R = ", round(a1_mat[["r"]][1],3)))  +
#  geom_text(x=1.6, y=-1.8, label=paste("p-fwe = ",a1_mat[["pfwer"]][1])) +
  theme_classic() + theme(axis.title.x = element_text(size = 16), 
                          axis.title.y = element_text(size = 16),
                          axis.text.x = element_text(size = 10), 
                          axis.text.y = element_text(size = 10))


setwd("~/Dropbox/Sid/R_files/STAGES_fmri/output/figures/")
cca_plot <- cowplot::plot_grid(pc_cor_plot, scatter_plot,clin_cor_plot, nrow = 1, ncol = 3,rel_widths = c(1,1,0.5))

ggsave(cca_plot, filename = paste0("fig3_cca.tiff"), device = "tiff")


brain_pos_ax <- ggdraw() + 
  draw_image("~/Dropbox/Sid/R_files/STAGES_fmri/output/brain_figs/11_fdr_axial.png")
#brain_pos_sag <- ggdraw() + 
#  draw_image("~/Dropbox/Sid/R_files/STAGES_fmri/output/brain_figs/11_fdr_sag.tif")

brain_neg_ax <- ggdraw() + 
  draw_image("~/Dropbox/Sid/R_files/STAGES_fmri/output/brain_figs/12_fdr_axial.png")
#brain_neg_sag <- ggdraw() + 
#  draw_image("~/Dropbox/Sid/R_files/STAGES_fmri/output/brain_figs/12_fdr_sag.tif")


net_pos <- ggdraw() + 
  draw_image("~/Dropbox/Sid/R_files/STAGES_fmri/output/NetworkPlots/network3_pce3_pce1_11.png", scale = 1)
net_neg <- ggdraw() +
  draw_image("~/Dropbox/Sid/R_files/STAGES_fmri/output/NetworkPlots/network3_pce3_pce1_12.png", scale = 1)

rig_pos <- ggdraw() + draw_image("~/Dropbox/Sid/R_files/STAGES_fmri/output/figures/clin_pos_ridge_plot.tiff")
rig_neg <- ggdraw() + draw_image("~/Dropbox/Sid/R_files/STAGES_fmri/output/figures/clin_neg_ridge_plot.tiff")


setwd("~/Dropbox/Sid/R_files/STAGES_fmri/output/surfaces")
rh_l <- ggdraw() + draw_image(paste0("t_analysis_rh_",11,"_inflated_transparent_reds_lat.jpg"), scale = 1)
rh_m <- ggdraw() + draw_image(paste0("t_analysis_rh_",11,"_inflated_transparent_reds_med.jpg"), scale = 1)
lh_l <- ggdraw() + draw_image(paste0("t_analysis_lh_",11,"_inflated_transparent_reds_lat.jpg"), scale = 1)
lh_m <- ggdraw() + draw_image(paste0("t_analysis_lh_",11,"_inflated_transparent_reds_med.jpg"), scale = 1)
setwd("~/Dropbox/Sid/R_files/STAGES_fmri/output/surfaces/colourbars")
cb <- ggdraw() + draw_image(paste0("colourbar_",11,".tiff"))
setwd("~/Dropbox/Sid/R_files/STAGES_fmri/output/subcortex/imgs/")
subcortl <- ggdraw() + draw_image(paste0("t_sub_",11,"_l.png"))
subcortr <- ggdraw() + draw_image(paste0("t_sub_",11,"_r.png"))

surf <- cowplot::plot_grid(rh_l, lh_l,
                           rh_m, lh_m, 
                           nrow=2, ncol = 2)

#brainsl <- cowplot::ggdraw(surf) + cowplot::draw_plot(subcortl, x = -.01, y = 0.02 , scale = 0.5)
brainsl <- cowplot::ggdraw(surf) + cowplot::draw_plot(subcortl, x = 0.28, y = -0.35 , scale = 0.258)

brains <- cowplot::ggdraw(brainsl) + cowplot::draw_plot(subcortr, x = -0.28, y = -0.35 , scale = 0.258)

brains
#add colour bar
brains_cb <- cowplot::plot_grid(brains, cb, rel_widths = c(1, 0.1), rel_heights = c(1,0.3), 
                                ncol = 2)



ggsave(plot = brains_cb, filename = "temp", device = "png")

brains_x <- ggdraw() + draw_image("temp", scale = 1)



p1 <- cowplot::plot_grid(brain_pos_ax, net_pos, rig_pos,brains_x,
                         nrow = 1, ncol = 4, rel_widths = c(0.7,1.5,1.2,1.2), rel_heights = c(1,0.8,1.4,1))
ggsave(p1, filename = paste0("~/Dropbox/Sid/R_files/STAGES_fmri/output/figures/fig3_panel1.tiff"), device = "tiff")






setwd("~/Dropbox/Sid/R_files/STAGES_fmri/output/surfaces")
rh_l <- ggdraw() + draw_image(paste0("t_analysis_rh_",12,"_inflated_transparent_reds_lat.jpg"), scale = 1)
rh_m <- ggdraw() + draw_image(paste0("t_analysis_rh_",12,"_inflated_transparent_reds_med.jpg"), scale = 1)
lh_l <- ggdraw() + draw_image(paste0("t_analysis_lh_",12,"_inflated_transparent_reds_lat.jpg"), scale = 1)
lh_m <- ggdraw() + draw_image(paste0("t_analysis_lh_",12,"_inflated_transparent_reds_med.jpg"), scale = 1)
setwd("~/Dropbox/Sid/R_files/STAGES_fmri/output/surfaces/colourbars")
cb <- ggdraw() + draw_image(paste0("colourbar_",12,".tiff"))
setwd("~/Dropbox/Sid/R_files/STAGES_fmri/output/subcortex/imgs/")
subcortl <- ggdraw() + draw_image(paste0("t_sub_",12,"_l.png"))
subcortr <- ggdraw() + draw_image(paste0("t_sub_",12,"_r.png"))

surf <- cowplot::plot_grid(rh_l, lh_l,
                           rh_m, lh_m, 
                           nrow=2, ncol = 2)

#brainsl <- cowplot::ggdraw(surf) + cowplot::draw_plot(subcortl, x = -.01, y = 0.02 , scale = 0.5)
brainsl <- cowplot::ggdraw(surf) + cowplot::draw_plot(subcortl, x = 0.28, y = -0.35 , scale = 0.258)

brains <- cowplot::ggdraw(brainsl) + cowplot::draw_plot(subcortr, x = -0.28, y = -0.35 , scale = 0.258)

brains
#add colour bar
brains_cb <- cowplot::plot_grid(brains, cb, rel_widths = c(1, 0.1), rel_heights = c(1,0.3), 
                                ncol = 2)


ggsave(plot = brains_cb, filename = "temp", device = "png")

brains_x <- ggdraw() + draw_image("temp", scale = 1)


p2 <- cowplot::plot_grid(brain_neg_ax, net_neg, rig_neg,brains_x,
                         nrow = 1, ncol = 4, rel_widths = c(0.7,1.5,1.2,1.2), rel_heights = c(1,0.8,1.4,1))
ggsave(p2, filename = paste0("~/Dropbox/Sid/R_files/STAGES_fmri/output/figures/fig3_panel2.tiff"), device = "tiff")



#p1 <- p1 + draw_text("Increased connectivity associated with worse functional outcome", x = 0.4, y =0.95) +
#  draw_text("Decreased connectivity associated with worse functional outcome", x = 0.4, y =0.45)
  

cowplot::plot_grid(cca_plot, p1, nrow = 2, ncol = 1, rel_heights  =c(0.7,1.4))
ggsave2(filename = "test.tiff", device = "tiff", width = 16, height = 13)


#================================================================
# 12-month medication_brain to behaviour cca (exploratory)


setwd("~/Dropbox/Sid/R_files/STAGES_fmri/data/cca/a3")
a1_mat <- rmatio::read.mat("a3_med_2.mat")

## correlation for PC
pcs_raw <- a1_mat[["pcs"]][["data"]][[1]]
variates <- as.data.frame(cbind(a1_mat[["U"]][,1], a1_mat[["V"]][,1]))
cor_pc <- as.data.frame(cbind(cor(pcs_raw, variates[,1]), 1:dim(a1_mat[["A"]])[1]))

pc_cor_plot <- ggplot(data=cor_pc, aes(x=as.factor(V2), y=V1)) + 
  geom_bar(stat="identity") + xlab("\u0394 Brain Principal Components") + 
  ylab("Canonical Loading") +
  theme_classic() + theme(axis.title.x = element_text(size = 16), 
                          axis.title.y = element_text(size = 16),
                          axis.text.x = element_text(size = 10), 
                          axis.text.y = element_text(size = 10))

pc_cor_plot
## correlation for PC - add astrix for sig
source("~/Dropbox/Sid/R_files/functions/cca_sig_test.R")
#p_brain_fdr <- cca_sig_test(variates[,1], pcs_raw, b = 5000, method = "1")
label.df <- data.frame(V2 = c("5"),
                       V1 = c(-0.5408)) #added. .1 for space

pc_cor_plot  <- pc_cor_plot   + geom_text(data = label.df, label = "*",
                                            size=6)

## correlation for clin
clin_raw <- a1_mat[["clin"]][["data"]][[1]]
variates <- as.data.frame(cbind(a1_mat[["U"]][,1], a1_mat[["V"]][,1]))
cor_clin <- as.data.frame(cbind(cor(clin_raw, variates[,2]), 1:dim(a1_mat[["B"]])[1]))

cor_clin$V2 <- factor(cor_clin$V2, labels = c("BPRS-psychotic", "BPRS-affect", "BPRS-activation", "BPRS-negative",
                                              "BPRS-disorganization","SANS-Total", "HAM-Depression", "HAM-Anxiety",
                                              "WHO-QoL", "QLS"))
clin_cor_plot <- ggplot(data=cor_clin, aes(x=as.factor(V2), y=V1)) + 
  geom_bar(stat="identity") + xlab("\u0394 Behavioural Scales") + 
  ylab("Canonical Loading")  + theme_classic() + theme(axis.title.x = element_text(size = 16), 
                                        axis.title.y = element_text(size = 16),
                                        axis.text.x = element_text(size = 10, angle = 45, hjust = 1), 
                                        axis.text.y = element_text(size = 10))

#add astrix to sig scales
source("~/Dropbox/Sid/R_files/functions/cca_sig_test.R")
#p_clin_fdr <- cca_sig_test(variates[,2], clin_raw, b = 5000, method = "1")
label.df <- data.frame(V2 = c("BPRS-affect","BPRS-negative","SANS-Total", "QLS"),
                       V1 = c(0.5616, 0.8181, 0.6576,-0.5848)) #added. .1 for space

clin_cor_plot <- clin_cor_plot  + geom_text(data = label.df, label = "*",
                                            size=6)




cca_scatter <- as.data.frame(cbind(a1_mat[["U"]][,1], a1_mat[["V"]][,1]))
colnames(cca_scatter) <- c("\u0394 Brain Canonical Variate", "\u0394 Behaviour Canonical Variate")

scatter_plot <- ggplot(cca_scatter, aes(x=`Δ Brain Canonical Variate`, y=`Δ Behaviour Canonical Variate`)) +
  geom_point(size=4, color = "#56B4E9") +  geom_smooth(method='lm',linetype = 2, colour = "grey", fill = "light grey") +
  geom_text(x=1.6, y=-1.6, label=paste0("R = ", round(a1_mat[["r"]][1],3)))  +
  geom_text(x=1.6, y=-1.8, label=paste("p-fwe = ",a1_mat[["pfwer"]][1])) +
  theme_classic() + theme(axis.title.x = element_text(size = 16), 
                          axis.title.y = element_text(size = 16),
                          axis.text.x = element_text(size = 10), 
                          axis.text.y = element_text(size = 10))




setwd("~/Dropbox/Sid/R_files/STAGES_fmri/output/figures/")
cca_plot <- cowplot::plot_grid(pc_cor_plot, scatter_plot,clin_cor_plot, nrow = 1, ncol = 3, scale = 0.9)
ggsave(cca_plot, filename = paste0("fig4_cca.tiff"), device = "tiff", width = 14, height = 6)

#


brain_pos_ax <- ggdraw() + 
  draw_image("~/Dropbox/Sid/R_files/STAGES_fmri/output/brain_figs/13_12med_cca_axial.tif")
brain_pos_sag <- ggdraw() + 
  draw_image("~/Dropbox/Sid/R_files/STAGES_fmri/output/brain_figs/13_12med_cca_sag.tif")

brain_neg_ax <- ggdraw() + 
  draw_image("~/Dropbox/Sid/R_files/STAGES_fmri/output/brain_figs/14_12med_cca_axial.tif")
brain_neg_sag <- ggdraw() + 
  draw_image("~/Dropbox/Sid/R_files/STAGES_fmri/output/brain_figs/14_12med_cca_sag.tif")


net_pos <- ggdraw() + 
  draw_image("~/Dropbox/Sid/R_files/STAGES_fmri/output/NetworkPlots/network3_pce3_pce1_13.png", scale = 1)
net_neg <- ggdraw() +
  draw_image("~/Dropbox/Sid/R_files/STAGES_fmri/output/NetworkPlots/network3_pce3_pce1_14.png", scale = 1)

#rig_pos <- ggdraw() + draw_image("~/Dropbox/Sid/R_files/STAGES_fmri/output/figures/old/clin_pos_redge_plot.tiff")
#rig_neg <- ggdraw() + draw_image("~/Dropbox/Sid/R_files/STAGES_fmri/output/figures/old/clin_neg_redge_plot.tiff")


setwd("~/Dropbox/Sid/R_files/STAGES_fmri/output/surfaces")
rh_l <- ggdraw() + draw_image(paste0("t_analysis_rh_",13,"_inflated_transparent_reds_lat.jpg"), scale = 1)
rh_m <- ggdraw() + draw_image(paste0("t_analysis_rh_",13,"_inflated_transparent_reds_med.jpg"), scale = 1)
lh_l <- ggdraw() + draw_image(paste0("t_analysis_lh_",13,"_inflated_transparent_reds_lat.jpg"), scale = 1)
lh_m <- ggdraw() + draw_image(paste0("t_analysis_lh_",13,"_inflated_transparent_reds_med.jpg"), scale = 1)
setwd("~/Dropbox/Sid/R_files/STAGES_fmri/output/surfaces/colourbars")
cb <- ggdraw() + draw_image(paste0("colourbar_",13,".tiff"))
setwd("~/Dropbox/Sid/R_files/STAGES_fmri/output/subcortex/imgs/")
subcortl <- ggdraw() + draw_image(paste0("t_sub_",13,"_l.png"))
subcortr <- ggdraw() + draw_image(paste0("t_sub_",13,"_r.png"))

surf <- cowplot::plot_grid(rh_l, lh_l,
                           rh_m, lh_m, 
                           nrow=2, ncol = 2)

#brainsl <- cowplot::ggdraw(surf) + cowplot::draw_plot(subcortl, x = -.01, y = 0.02 , scale = 0.5)
brainsl <- cowplot::ggdraw(surf) + cowplot::draw_plot(subcortl, x = 0.27, y = -0.33 , scale = 0.4)

brains <- cowplot::ggdraw(brainsl) + cowplot::draw_plot(subcortr, x = -0.27, y = -0.33 , scale = 0.4)
#add colour bar
brains_cb <- cowplot::plot_grid(brains, cb, rel_widths = c(1, 0.1), rel_heights = c(1,0.3), 
                                ncol = 2)
ggsave(plot = brains_cb, filename = "temp", device = "png")

brains_x <- ggdraw() + draw_image("temp", scale = 1)






p1 <- cowplot::plot_grid(brain_pos_ax, brain_pos_sag,net_pos, brains_x,
                         nrow = 1, ncol = 4, rel_widths = c(1,0.9,1.4,1), rel_heights = c(1,0.9,1.4,1))
ggsave(p1, filename = paste0("fig4_panel1.tiff"), device = "tiff")





setwd("~/Dropbox/Sid/R_files/STAGES_fmri/output/surfaces")
rh_l <- ggdraw() + draw_image(paste0("t_analysis_rh_",14,"_inflated_transparent_reds_lat.jpg"), scale = 1)
rh_m <- ggdraw() + draw_image(paste0("t_analysis_rh_",14,"_inflated_transparent_reds_med.jpg"), scale = 1)
lh_l <- ggdraw() + draw_image(paste0("t_analysis_lh_",14,"_inflated_transparent_reds_lat.jpg"), scale = 1)
lh_m <- ggdraw() + draw_image(paste0("t_analysis_lh_",14,"_inflated_transparent_reds_med.jpg"), scale = 1)
setwd("~/Dropbox/Sid/R_files/STAGES_fmri/output/surfaces/colourbars")
cb <- ggdraw() + draw_image(paste0("colourbar_",14,".tiff"))
setwd("~/Dropbox/Sid/R_files/STAGES_fmri/output/subcortex/imgs/")
subcortl <- ggdraw() + draw_image(paste0("t_sub_",14,"_l.png"))
subcortr <- ggdraw() + draw_image(paste0("t_sub_",14,"_r.png"))

surf <- cowplot::plot_grid(rh_l, lh_l,
                           rh_m, lh_m, 
                           nrow=2, ncol = 2)

#brainsl <- cowplot::ggdraw(surf) + cowplot::draw_plot(subcortl, x = -.01, y = 0.02 , scale = 0.5)
brainsl <- cowplot::ggdraw(surf) + cowplot::draw_plot(subcortl, x = 0.27, y = -0.33 , scale = 0.4)

brains <- cowplot::ggdraw(brainsl) + cowplot::draw_plot(subcortr, x = -0.27, y = -0.33 , scale = 0.4)
#add colour bar
brains_cb <- cowplot::plot_grid(brains, cb, rel_widths = c(1, 0.1), rel_heights = c(1,0.3), 
                                ncol = 2)
ggsave(plot = brains_cb, filename = "temp", device = "png")

brains_x <- ggdraw() + draw_image("temp", scale = 1)




p2 <- cowplot::plot_grid(brain_neg_ax, brain_neg_sag, net_neg, brains_x,
                         nrow = 1, ncol = 4, rel_widths = c(1,0.9,1.4,1,1), rel_heights = c(1,0.9,1.4,1,1))
ggsave(p2, filename = paste0("fig4_panel2.tiff"), device = "tiff")

















###################################
##plot dose correlation   figure #   #### Not reporting this analysis for now
##################################

######==================
# Medication olz correlation
#####==================

library(ggplot2)

#med a3
setwd("~/Dropbox/Sid/R_files/STAGES_fmri/data/cca/dose_analysis")
a1_mat <- rmatio::read.mat("12m_pcs_lt0.mat")

## correlation for PC
pcs_raw <- read.csv('pcs_gt0_12m.csv')#[,-1]
variates <- as.data.frame(cbind(a1_mat[["U"]], a1_mat[["V"]]))
cor_pc <- as.data.frame(cbind(cor(pcs_raw, variates[,1]), 1:length(a1_mat[["A"]])))

pc_cor_plot <- ggplot(data=cor_pc, aes(x=as.factor(V2), y=V1)) + 
  geom_bar(stat="identity") + xlab("Brain Principal Components") + 
  ylab("Canonical Loading") + ylim(-1, 1)+ 
  theme_classic() + theme(axis.title.x = element_text(size = 16), 
                                                                   axis.title.y = element_text(size = 16),
                                                                   axis.text.x = element_text(size = 10), 
                                                                   axis.text.y = element_text(size = 10))

##correlation for clin
## correlation for clin
#clin_raw <- read.csv('olz_dose_gt0_12m.csv')
#variates <- as.data.frame(cbind(a1_mat[["U"]], a1_mat[["V"]]))
#cor_clin <- as.data.frame(cbind(cor(clin_raw, variates[,2]), 1:length(a1_mat[["B"]])))
#cor_clin$V2 <- factor(cor_clin$V2, labels = c("olz_eq_dose"))
#clin_cor_plot <- ggplot(data=cor_clin, aes(x=as.factor(V2), y=V1)) + 
#  geom_bar(stat="identity") + xlab("PC") + 
#  ylab("correlation to behavior variate") + 
#  ylim(-1, 1) + theme(axis.text.x = element_text(angle = 45, hjust = 1))



#scatter
cca_scatter <- as.data.frame(cbind(a1_mat[["U"]], a1_mat[["V"]]))
colnames(cca_scatter) <- c("Brain-change Canonical Variate", "Antipsychotic Dose Canonical Variate")
scatter_plot <- ggplot(cca_scatter, aes(y=`Brain-change Canonical Variate`, x=`Antipsychotic Dose Canonical Variate`)) +
  geom_point(size=4, color = "#56B4E9") +  geom_smooth(method='lm',linetype = 2, colour = "grey", fill = "light grey") +
  geom_text(x=1.6, y=-1.5, label=paste0("R = ", round(a1_mat[["r"]][1],3)))  +
  geom_text(x=1.6, y=-1.4, label=paste("p-fwe = ",a1_mat[["pfwer"]][1])) +
  theme_classic() + theme(axis.title.x = element_text(size = 16), 
                          axis.title.y = element_text(size = 16),
                          axis.text.x = element_text(size = 10), 
                          axis.text.y = element_text(size = 10))

p1 <- cowplot::plot_grid(pc_cor_plot, scatter_plot,nrow = 1)
#read in brain figs
net_pos <- ggdraw() + draw_image("~/Dropbox/Sid/R_files/STAGES_fmri/output/NetworkPlots/network2_pce3_pce1_9.png", scale = 1.3)
net_neg <- ggdraw() + draw_image("~/Dropbox/Sid/R_files/STAGES_fmri/output/NetworkPlots/network2_pce3_pce1_10.png", scale = 1.3)

brain_pos <- ggdraw() + draw_image("~/Dropbox/Sid/R_files/STAGES_fmri/output/brain_figs/13.png")
brain_neg <- ggdraw() + draw_image("~/Dropbox/Sid/R_files/STAGES_fmri/output/brain_figs/14.png")


rig_pos <- ggdraw() + draw_image("~/Dropbox/Sid/R_files/STAGES_fmri/output/figures/dose_pos_redge_plot.tiff")
rig_neg <- ggdraw() + draw_image("~/Dropbox/Sid/R_files/STAGES_fmri/output/figures/dose_neg_redge_plot.tiff")
  
p2 <- cowplot::plot_grid(brain_pos, net_pos, rig_pos, nrow=1, ncol=3)
p3 <- cowplot::plot_grid(brain_neg, net_neg, rig_neg, nrow=1, ncol=3)

cowplot::plot_grid(p1, p2, p3, nrow = 3, ncol = 1)







#####Explore netplots
net_plot_list <- list()
  setwd("~/Dropbox/Sid/R_files/STAGES_fmri/output/NetworkPlots/")
for (p in 1:14) {
  net_plot_list[[p]] <- ggdraw() + draw_image(paste0("network3_pce3_pce1_",p,".png"), scale = 1.3)
}
  id <- c(1,3,7)
  ii <- c(2,4,8)
  #medication
  md <- c(6,10)
  mi <- c(5,9)
  
  cowplot::plot_grid(net_plot_list[[1]], net_plot_list[[5]], 
                     net_plot_list[[9]], net_plot_list[[2]], 
                     net_plot_list[[6]], net_plot_list[[10]], nrow = 2, ncol = 3)
  
  
  cowplot::plot_grid(net_plot_list[[1]], net_plot_list[[3]], 
                     net_plot_list[[7]], net_plot_list[[2]], 
                     net_plot_list[[4]], net_plot_list[[8]], nrow = 2, ncol = 3)
  