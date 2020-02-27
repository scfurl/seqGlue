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


#' Removes data (collapses) with duplicated names
#'
#' @description Will collapse a matrix based on duplicated labels that have the greatest variance
#' @param matrix Input object with unique rownames
#' @param labels the non-unique labels (i.e. Gene Symbol)
#' @export
#' 

remove_duplicated_rows<-function(matrix, labels){
  rn<-labels
  if (length(make.unique(labels)) != length(unique(labels))) {
    unique_rn<-rownames(matrix)[!rn %in% rn[duplicated(rn)]]
    duplicated_rn<-rownames(matrix)[rn %in% rn[duplicated(rn)]]
    duplicated_lab<-rn[rn %in% rn[duplicated(rn)]]
    duplicated<-matrix[rn %in% rn[duplicated(rn)],]
    ts = data.table(rowV = matrixStats::rowVars(as.matrix(duplicated)), rn_D=duplicated_lab, rn_U=duplicated_rn)
    tokeep<-c(ts[ , .SD[which.max(rowV)], by = rn_D]$rn_U, unique_rn)
    matrix_out<-matrix[rownames(matrix) %in% tokeep,]
    rownames(matrix_out)<-labels[rownames(matrix) %in% tokeep]
    return(matrix_out)
  }
}
