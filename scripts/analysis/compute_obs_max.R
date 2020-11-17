compute_obs_max <- function(mat, write.comps=FALSE) {
  
  thresh <- c(1.301, 2, 3)
  vals <- vector()
  for (t in 1:length(thresh)) {
    temp <- matrix(nrow=316, ncol=316) #next 3 lines of code just converting the pvector back into a con mat
    temp[upper.tri(temp)] <- as.matrix(mat)
    temp <- forceSymmetric(temp, uplo = "U")
    diag(temp) <- 0
    
    temp_bin<- binarize(temp, threshold=thresh[t]) # this is what to change when looking at sig comps
    #View(igraph::degree(temp.i))
    temp.i <- graph_from_adjacency_matrix(as.matrix(temp_bin), 
                                          weighted = NULL,
                                          mode = c("undirected")) #maybe need to add as.matrix
    # splits into comps (finds largest comp and gives size edges)
    decomposed_comps <-  decompose.graph(temp.i)
    vals[t] <- max(unlist((lapply(decomposed_comps, gsize))))
    if(write.comps == TRUE) { 
      filename <- paste0(thresh[t], "_obs_comp.csv")
      write_largest_comp(matx=temp_bin, filename = filename)
      
    }
    
  }
  return(vals)
  
  
}
