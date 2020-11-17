#This function takins in 2 natricies (collums are images) and computes a cunjunction contrast (the intersection)
#max size of componant in edges that are significant in BOTH contrasts, 


compute_nulls_max_conj <- function(mat1, mat2) {
 # alpha <- c(0.968, 0.90, 0.776) #1-sqrt(p)
  alpha <- c(0.999, 0.99, 0.95) 
  max.comp.list  <- vector(length=dim(mat1)[2])
  
  
  for (a in 1:length(alpha)){
    bs <- dim(mat1)[2]
    for (k in 1:bs) {
      temp1 <- vec_2_mat(mat1[,k], 316, 0)
      temp2 <- vec_2_mat(mat2[,k], 316, 0)
      #binarise
      temp1_bin<- binarize(temp1, threshold=alpha[a]) 
      temp2_bin<- binarize(temp2, threshold=alpha[a]) 
      
      temp <- temp1_bin + temp2_bin
      temp_bin <- binarize(temp, threshold = 1.9) 
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
