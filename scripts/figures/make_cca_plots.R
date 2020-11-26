#' ===========================
#' Script name: make_cca_plots
#'
#' Purpose of script: Visualizes and explore output of permcca.m script
#' Author: Sid Chopra (sid.chopra@monash.edu)
#' Date Created: 2/8/2020
#' ---------------------------
#' Inputs:
#' matfile = a list matching the .mat file saved output from permcca.m script (can
#' be read into r using rmatio::read.mat()
#' use.coef: Boolean operator to indicate whether canonical coefficients or canonical
#' loading should be plotted. Default is loadings.
#' ---------------------------
#' Outputs:
#' cowplot/ggplot object
#' ----------------------------
#' Notes:
#' ---------------------------
#' Example usage:
#' mat <- rmatio::read.mat("permcca_output1.mat")
#' make_cca_plots(mat, use.coef = F)
#' ---------------------------
#' Required packages:
#' library(rmatio)
#' library(ggplot2)
#' library(cowplot)
#'
#' ===========================



make_cca_plots <- function(mat, use.coef = FALSE) {
  if(use.coef == FALSE){
    ## compute and plot left side loadings
    pcs_raw <- mat["pc"][[1]] #will have to change this index to match your own data
    variates <- as.data.frame(cbind(mat[["U"]][,1], mat[["V"]][,1]))
    cor_pc <- as.data.frame(cbind(cor(pcs_raw, variates[,1]), 1:dim(mat["pc"][[1]])[2]))
    pc_plot <- ggplot(data=cor_pc, aes(x=as.factor(V2), y=V1)) +
      geom_bar(stat="identity") + xlab("Brain Principal Components") +
      ylab("Canonical Loading") + ylim(-1, 1) +
      theme_classic()
    pc_plot
    ## compute and plot right side loadings
    clin_raw <- mat["scales"][[1]] #will have to change this index to match your own data
    variates <- as.data.frame(cbind(mat[["U"]][,1], mat[["V"]][,1]))
    cor_clin <- as.data.frame(cbind(cor(clin_raw, variates[,2]), 1:dim(mat["B"][[1]])[1]))
    cor_clin$V2 <- factor(cor_clin$V2, labels = c(1:11))
    clin_plot <- ggplot(data=cor_clin, aes(x=as.factor(V2), y=V1)) +
      geom_bar(stat="identity") + xlab("Behavioural Scales") +
      ylab("Canonical Loading") +
      ylim(-1, 1) + theme_classic()
    clin_plot
  }


  if(use.coef == TRUE) {
    ## plot left side coefs
    pcs_coef <- as.data.frame(as.matrix(cbind(mat[["A"]], 1:dim(mat[["A"]])[1])))
    pcs_coef$V1 <- round(as.numeric(as.character(pcs_coef$V1)), 3)
    pc_plot <- ggplot(data=pcs_coef, aes(x=as.factor(V3), y=V1)) +
      geom_bar(stat="identity") + xlab("Brain Principal Components") +
      ylab("Brain Cannonical Coeffeicents") +
      ylim(-1, 1) +  theme_classic()
    pc_plot

    # plot right side coefs
    pcs_clin <- as.data.frame(as.matrix(cbind(mat[["B"]], 1:dim(mat[["B"]])[1])))
    pcs_clin$V1 <- round(as.numeric(as.character(pcs_clin$V1)), 3)
    #pcs_clin$V3 <- factor(pcs_clin$V3, labels = c("SOFAS", "BPRS"))
    pcs_clin$V2 <- factor(pcs_clin$V2, labels = c(1:length(pcs_clin$V2)))
    clin_plot <- ggplot(data=pcs_clin, aes(x=as.factor(V2), y=V1)) +
      geom_bar(stat="identity") +
      xlab("Clinical and Functional Scales") +
      ylab("Behaviour Cannonical Coeffeicents") +
      ylim(-1.06, 1)   + theme_classic()
    clin_plot
  }

  #make scatter plots
  cca_scatter <- as.data.frame(cbind(mat[["U"]][,1], mat[["V"]][,1]))
  colnames(cca_scatter) <- c("Brain-change Canonical Variate", "Behaviour-change Canonical Variate")
  scatter_plot <- ggplot(cca_scatter,
                         aes(x=`Behaviour-change Canonical Variate`,
                             y=`Brain-change Canonical Variate`)) +
    geom_point(size=4, color = "#56B4E9") +
    geom_smooth(method='lm',linetype = 2, colour = "grey", fill = "light grey") +
    geom_text(x=1.6, y=-1.6, label=paste0("R = ", round(mat[["r"]][1],3)))  +
    geom_text(x=1.6, y=-1.8, label=paste("p-fwe = ",mat[["pfwer"]][1])) +
    theme_classic()
  scatter_plot

  #combine plots
  combined_plot <- cowplot::plot_grid(pc_plot, scatter_plot,clin_plot,
                                      nrow = 1,rel_widths = c(0.5,1,1), rel_heights = c(0.5,1,1)) #change the rel_wid/hei to make it pretty


  return(combined_plot)

}
