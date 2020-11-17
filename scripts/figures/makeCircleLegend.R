#make legend for circle plota

makeCircleLegend <- function(){
  ####Make circos legend
  leg <- matrix(1, 8, 8)
  rownames(leg) <- colnames(leg) <- c("Visual", "Somatomotor", "Dorsal Attention", "Ventral Attention" , "Limbic",
                                      "Frontoparietal", "Default", "Subcortical")
  grid.col = c("Visual" = "#9C755FFF", "Somatomotor" = "#F28E2BFF",
               "Dorsal Attention" = "#E15759FF", 'Ventral Attention' = "#76B7B2FF",
               'Limbic' = "#59A14FFF", "Frontoparietal" = '#EDC948FF',
               'Default' = "#B07AA1FF", "Subcortical" = "#4E79A7FF")
  col_mat[leg !=0] = "#00000000"
  
 plot <-  chordDiagram(leg, annotationTrack = "grid", 
               annotationTrackHeight = c(0.4), grid.col = grid.col, col = col_mat) +  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                facing = "bending.inside", niceFacing = TRUE, adj = c(0.5,-2), cex = 0.9)
  }, bg.border = NA)
  
  return(plot)
}
