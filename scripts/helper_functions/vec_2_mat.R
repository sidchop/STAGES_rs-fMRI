vec_2_mat <- function(vec, length, diag) {
  temp <- matrix(nrow=length, ncol=length)
  temp[upper.tri(temp)] <- as.matrix(vec)
  temp <- Matrix::forceSymmetric(temp, uplo = "U")
  diag(temp) <- diag
return(as.matrix(temp))
}
