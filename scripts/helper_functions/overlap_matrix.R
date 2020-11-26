###This funtion computes overlap between 2 binary marricies. Useful for conjunction inference in fwe corrected matricies. 
overlap_matrix <- function(x, y)
  if(dim(x) != dim(y)) {print("Matrix dimentions do not match")} 
xy <- matric(n)