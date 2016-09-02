FPCA_2D_score_fast <-
function(X){
  ##read the dimension of the array
  m <- dim(X)[1]
  n <- dim(X)[2]
  k <- dim(X)[3]
  ####calculate the fourier series for each image and combine then into a big matrix
  FC <- matrix(rep(0,(2*m-1)*(2*n-1)*k),ncol=k)
  for (i in 1:k){
    FC[,i] <-  as.vector(FFT2FS_2D(X[,,i]))
  }
  ####calculate the svd of the C matrix
  C_svd <- fast.svd(FC/sqrt(k))
  C_value <- (C_svd$d)^2
  C_vector <- C_svd$u
  C_score <- t(FC) %*% C_vector
  FPCA <- list("eigen_value"=C_value,"FPC_score"=C_score,"Eigen_vector"=C_vector)
  return(FPCA)
}
