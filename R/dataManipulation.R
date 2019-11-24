#' @export
#' @import Matrix
tf_idf_transform<-function(matrix){
  idf <- log( ncol(matrix) / ( 1 + Matrix::rowSums(matrix != 0) ) )
  idf<-.sparseDiagonal(x=idf)
  tf_idf <- crossprod(matrix, idf)
  colnames(tf_idf) <- rownames(matrix)
  tf_idf_out<-t(tf_idf / sqrt( rowSums( tf_idf^2 ) ))
  return(tf_idf_out)
}
