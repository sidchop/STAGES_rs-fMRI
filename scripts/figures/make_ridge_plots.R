#' ===========================
#' Script name: make_ridge_plots 
#'
#' Purpose of script: Take the input of a weighted adjacency matrix and network labels, and
#' plot distributions of values for each network as a ordered, stacked ridge plot 
#' Author: Sid Chopra (sid.chopra@monash.edu)
#' Date Created: 5/08/2020
#' ---------------------------
#' Inputs: 
#' loading_matrix = weighted and symmetric matrix of size n x n.
#' network_labels =  n x 1 vector of network ids. Each network should be represented
#' as a unique number. Each row is a different node., with label for each network 
#' network_names = A string vector of networks names corresponding to the network_labels,
#' of size unique(n)
#' loadings(optional) = "pos", "neg" or "both" (show only positive/negative or all non-zero loadings)
#' quantile(optional) = Numeric value between 0-1, to threshold the matrix (i.e. .95 = top 5% of loadings only)
#' ---------------------------
#' Outputs: 
#' ggplot object
#' ----------------------------
#' Notes:
#' ---------------------------
#' Example usage:
#' loading_matrix <- as.matrix(read.table("~/Dropbox/Sid/R_files/STAGES_fmri/data/cca/a2/cor_to_brain_pos.txt", header = F))
#' network_labels <- read.table("~/Dropbox/Sid/R_files/STAGES_fmri/output/NetworkPlots/labels_mtl_str_thal.txt")
#' network_names <- c("Visual", "Somatomotor", "Dorsal Attention",  "Ventral Attention", "Limbic", "Frontoparietal", 
#' "Default", "Amyg & Hippo", "Thalamus", "Striatum") 
#' 
#'  make_ridge_plots(loading_matrix, network_labels, network_names, loadings = "both", quantile=0.50)
#' ---------------------------
#' Required packages:
#' library(dplyr)
#' library(forcats)
#' library(ggridges)
#' library(tidyr)
#'
#'=====================


make_ridge_plots <- function(loading_matrix, 
                             network_labels, 
                             network_names, 
                             loadings = "both",
                             quantile=NULL) {
 
  #clump netowrks together in a list 
  
  loading_matrix <- as.matrix(loading_matrix)
  diag(loading_matrix) <- 0
  network.count <- unique(network_labels)
  n.net <- dim(network.count)[1]
  network_loadings <- list()
  for (n in 1:n.net) {
    network_loadings[[n]] <- loading_matrix[which(network_labels==network.count[n,]), ]
  }
  
  #make ridgeplot dataframe
  network_loadings_ridges <- list()
  for (n in 1:length(network_loadings)) {
    network_loadings_ridges[[n]] <- cbind(c(network_loadings[[n]]), rep(network_names[n], length(c(network_loadings[[n]]))))
  }
  
  loadings_ridges_vec <- do.call("rbind", network_loadings_ridges)
  
  #select both, pos or neg loadings only
  if(loadings=="both"){
    loadings_ridges_vec_fil <- as.data.frame(loadings_ridges_vec[which(loadings_ridges_vec[,1]!=0),])
  }
  if(loadings=="pos"){
    loadings_ridges_vec_fil <- as.data.frame(loadings_ridges_vec[which(loadings_ridges_vec[,1]>0),])
  }
  if(loadings=="neg"){
    loadings_ridges_vec_fil <- as.data.frame(loadings_ridges_vec[which(loadings_ridges_vec[,1]<0),])
  }
  
  
  colnames(loadings_ridges_vec_fil) <- c("loading", "network")
  loadings_ridges_vec_fil$network <- as.factor(loadings_ridges_vec_fil$network)
  loadings_ridges_vec_fil$loading <- as.numeric(as.character(loadings_ridges_vec_fil$loading))

  #select top x% of edges (i.e. .95 = top 5%)
  if(!is.null(quantile)) {
    quantiles <- quantile(loadings_ridges_vec_fil$loading, c(quantile))
    loadings_ridges_vec_fil  <- loadings_ridges_vec_fil[which(loadings_ridges_vec_fil$loading > quantiles),]
    }

  

  
  #make ridge plot
  ridge_plot <- loadings_ridges_vec_fil %>%
    mutate(network = fct_reorder(network, loading, .fun='median')) %>%
    ggplot(aes(y=reorder(network, loading), x=loading)) + 
    geom_density_ridges(aes(fill=network), alpha=0.7) +
    theme_minimal() + 
    theme(legend.position = "none", axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title=element_text(size=18, vjust = 3)) +
    #scale_fill_manual(values = c("Visual" = "#9C755FFF", "Somatomotor" = "#F28E2BFF",             #optional, custom colours
    #                             "Dorsal Attention"= "#E15759FF","Ventral Attention" = "#76B7B2FF", 
    #                             "Limbic" = "#59A14FFF", "Frontoparietal" = "#EDC948FF", 
    #                             "Default" = "#B07AA1FF", "Amyg & Hippo" = "#5575b3FF", 
    #                             "Thalamus" = "#496499FF", "Striatum" = "#4E79A7FF")) + 
    xlab("Correlation to Brain Variate") + 
    ylab("Network")
  
  if(loadings=="neg") {ridge_plot <-   ridge_plot + scale_x_reverse()} #flip x axis with neg loadings, looks nicer
  #ggsave(plot = ridge_plot, filename = "dose_neg_redge_plot.tiff", device = "tiff")
  
  return(ridge_plot)
}
