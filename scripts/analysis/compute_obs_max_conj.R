compute_obs_max_conj <- function(mat1, mat2, write.comps = FALSE) {
  
  thresh <- c(1.301, 2, 3)
 # thresh <- c(0.651, 1, 1.5) #-log10(sqrt(p))
  vals <- vector()
  for (t in 1:length(thresh)) {
    temp1 <- vec_2_mat(mat1, 316, 0)
    temp2 <- vec_2_mat(mat2, 316, 0)
    
    #binarise
    temp1_bin<- binarize(temp1, threshold=thresh[t]) 
    temp2_bin<- binarize(temp2, threshold=thresh[t]) 
    
    temp <- temp1_bin + temp2_bin
    temp_bin <- binarize(temp, threshold = 1.9) 
    temp.i <- graph_from_adjacency_matrix(as.matrix(temp_bin), 
                                          weighted = NULL,
                                          mode = c("undirected"))
    
    decomposed_comps <-  decompose.graph(temp.i)
    
    vals[t] <- max(unlist((lapply(decomposed_comps, gsize))))
    if(write.comps == TRUE) { 
      filename <- paste0(thresh[t], "_obs_comp.csv")
      write_largest_comp(matx=temp_bin, filename = filename)
      
    }
    
  }
  return(vals)
  
  
}
