#Generate .node file for BrainNet-Viewer (matlab)


#load in MNI coords 

load("~/Dropbox/Sid/R_files/brainconn/data/stages_study.rda") 
mni <- stages_study[,2:4]
#read in network labs
network_labs <- c(rep(1,24), rep(2,29), rep(3, 16), rep(4, 16), rep(5, 3),
                  rep(6, 17), rep(7, 37), rep(1,23), rep(2,28), rep(3,18), 
                  rep(4, 18), rep(5, 2),  rep(6, 23), rep(7, 30), rep(8, 32))


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

comps[[11]] <- read.table("~/Dropbox/Sid/R_files/STAGES_fmri/data/cca/a2/cor_to_brain_FDR_pos.txt")
comps[[12]] <- read.table("~/Dropbox/Sid/R_files/STAGES_fmri/data/cca/a2/cor_to_brain_FDR_neg.txt")

comps[[13]] <- read.table("~/Dropbox/Sid/R_files/STAGES_fmri/data/cca/a3/cor_to_brain_FDR_pos.txt")
comps[[14]] <- read.table("~/Dropbox/Sid/R_files/STAGES_fmri/data/cca/a3/cor_to_brain_FDR_neg.txt")
comps[[15]] <- read.table("~/Dropbox/Sid/R_files/STAGES_fmri/data/cca/a3/cor_to_brain_FDR.txt")

# Not reporting the dose-respose analysis 
#comps[[13]] <- read.table("~/Dropbox/Sid/R_files/STAGES_fmri/data/cca/dose_analysis/cor_to_brain_pos_pct.txt")
#comps[[14]] <- read.table("~/Dropbox/Sid/R_files/STAGES_fmri/data/cca/dose_analysis/cor_to_brain_neg_pct.txt")

setwd("~/Dropbox/Sid/R_files/STAGES_fmri/output/brain_figs")
for (c in 1:length(comps)) {
network.graph <- igraph::graph_from_adjacency_matrix(as.matrix(comps[[c]]), mode = "undirected") 
degree <- igraph::degree(network.graph, mode = "total")
node.file <- cbind(mni, network_labs, degree, stages_study$ROI.Name)
write.table(node.file, paste0("analysis_", c, ".node"), row.names = F, col.names = F, quote = F)
write.table(as.matrix(comps[[c]]), paste0("analysis_", c, ".edge"), row.names = F, col.names = F, quote = F)

}


