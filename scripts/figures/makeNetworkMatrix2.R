## ===========================
## Script name:  makeNetworkMatrix2
##
## Purpose of script: Create a tile plot with independent upper and lower triangle inputs. 
## Author: Sid Chopra (sid.chopra@monash.edu)
## Date Created: 24/07/2020
## ---------------------------
## Inputs: 
## matrix1 =  n x n matrix
## matrix2 =  n x n matrix
## min/max functions (optional) allow you to control the min/max of colour bar independently for each matrix
## pal (optional) = scale_fill_gradientn input. e.g. vector of colour/hex names. 
## ---------------------------
## Outputs: 
## ggplot object
## ----------------------------
## Notes: Wrote this to visualize output from plot plotClassifiedEdges.R scripts. 
## ---------------------------
## Example usage:
## adj_mat <- read.table("t_pos_a1_baseline.txt")
## labels <- c(""Vis", "SomMot", "DorsAttn", "VentAttn", "Lim", "FPN", "DMN", "MTL", "Thal", "Stri")
## networks <- read.table("labels_mtl_str_thal.txt")
## pce <- plotClassifiedEdges(adj =  adj_mat, ids = networks, labels = labels)
## makeNetworkMatrix2(matrix1 = pce[[3]], matrix2 = pce[[1]]
## ---------------------------
## Required packages:  
## library(reshape2)
## library(scales)
## library(gtable)
## library(patchwork)
## library(ggplot2)
## ===========================



makeNetworkMatrix2 <- function(matrix1,
                               matrix2,
                               min1=NULL,
                               max1=NULL,
                               min2=NULL,
                               max2=NULL,
                               pal= c("light yellow", "yellow", "orange", "red", "dark red"),
                               title="") {
  library(reshape2)
  library(scales)
  library(gtable)
  library(patchwork)
  library(ggplot2)


  
  # Get lower triangle of the correlation matrix
  get_upper_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  # Get upper triangle of the correlation matrix
  get_lower_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  
  
  lower_tri <- get_lower_tri(matrix1) 
  
  melted_cormat <- melt(lower_tri, na.rm = TRUE)
  melted_cormat[melted_cormat==0] <- NA #sent zero values as NA so they come our grey
  head(melted_cormat)
  
  ggheatmap1 <- ggplot(data = melted_cormat, aes(x=Var2, y=Var1, fill=value)) + 
    geom_tile(color = "black")  +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 18, hjust = 1), 
          axis.text.y = element_text(vjust = 1, 
                                     size = 18, hjust = 1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          legend.text=element_text(size=14),
          legend.title = element_text(size=15),
          plot.margin=grid::unit(c(0,0,0,0), "mm"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA)) +
    coord_fixed() + scale_y_discrete(position = "right") + 
    scale_fill_gradientn(colours=pal,
                         limits=c(min1, max1), na.value  = "#eeeeee") +
    labs(fill = "Normalised\n proportion") + 
    guides(fill = guide_colourbar(nbin = 100))
  
  
  
  
  upper_tri <- get_upper_tri(matrix2) 
  melted_cormat <- melt(upper_tri, na.rm = TRUE) #set zero values as NA so they come out grey/"#eeeeee"
  melted_cormat[melted_cormat==0] <- NA 
  
  
  ggheatmap2 <- ggplot(data = melted_cormat, aes(x=Var2, y=Var1, fill=value)) + 
    geom_tile(color = "black")  +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 18, hjust = 0), 
          axis.text.y = element_text(vjust = 1, 
                                     size = 18, hjust = 1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          legend.position="left",
          legend.text=element_text(size=14),
          legend.title = element_text(size=15),
          plot.margin=grid::unit(c(0,0,0,0), "mm"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA)) +
    coord_fixed() + scale_y_discrete(position = "left") +
    scale_x_discrete(position = "top") +
    scale_fill_gradientn(colours=pal,
                         limits=c(min2, max2), na.value = "#eeeeee") +
    labs(fill = "No. of\n edges") +
    guides(fill = guide_colourbar(nbin = 100))
  
  
  ggheatmap1 <- ggheatmap1 +
    theme(
      rect = element_rect(fill = "transparent") # all rectangles
    )
  
  ggheatmap2 <- ggheatmap2 +
    theme(
      rect = element_rect(fill = "transparent") # all rectangles
    )
  
  layout <- c(area(1,1, 5,4), area(1,2,5,5))
  
  ggheatmap_fin <- ggheatmap2 +  ggheatmap1 +  plot_layout(design = layout)
  
  #Add title 
  ggheatmap_fin <- ggheatmap_fin + ggtitle(title)   
  
  return(ggheatmap_fin)
  
}

