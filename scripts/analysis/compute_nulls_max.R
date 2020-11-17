compute_nulls_max <- function(mat) {
  alpha <- c(0.999, 0.99, 0.95) #derived from quantiles of the w.vect dist
  fwe.thresh <- 0.05
  
  max.comp.list  <- vector(length=dim(mat)[2])
  
  
  for (a in 1:length(alpha)){
    bs <- dim(mat)[2]
    for (k in 1:bs) {
      temp <- matrix(nrow=316, ncol=316)
      
      temp[upper.tri(temp)] <- as.matrix(mat[,k])
      
      temp <- forceSymmetric(temp, uplo = "U")
      diag(temp) <- 0
      # calculate nulls 
      temp_bin<- binarize(temp, threshold=alpha[a]) 
      
      temp.i <- graph_from_adjacency_matrix(as.matrix(temp_bin), 
                                            weighted = NULL,
                                            mode = c("undirected"))
      
      decomposed_comps <-  decompose.graph(temp.i)
      
      max.comp.list[k] <- max(unlist((lapply(decomposed_comps, gsize))))
      print(paste(k, ": max = ", max.comp.list[k]))
    }
    write.table(max.comp.list, paste("max.comp.list_",alpha[a], sep = ""), col.names = F, row.names = F)
  }
  
  
  
}
