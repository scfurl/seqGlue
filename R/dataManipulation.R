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

#' Converts to long data format
#'
#' @description Will convert a data object (either DESeq2 object or a matrix) to a data table for further work.  
#' @param matrix Input object with unique rownames
#' @param coldata a factor or dataframe of factors used to add information to long data table, the string "all" can be used if the class of the input
#' object is a DESeq2 dataset and all coldata information is to be copied to long data table
#' @export
#' 
convert_to_long<-function(object, coldata="all", values_to="NormalizedExpression", names_to = "sample"){
  if(class(object) %in% c("DESeqTransform", "DESeqDataSet")){
    if(coldata=="all"){
      coldata=colData(object)
    }else{
      coldata<-data.frame(coldata, row.names = colnames(object))
    }
  }
  if(class(object) %in% c("matrix")){
    if(coldata=="all"){
      stop("The 'all' feature is not supported for matrix input")
    }else{
      coldata<-data.frame(coldata, row.names = colnames(object))
    }
  }
  dt<-data.table::as.data.table(assay(object))
  dt$Gene<-rownames(object)
  dtl<-dt %>% tidyr::pivot_longer(-Gene, values_to = values_to, names_to = names_to)
  dtout<-cbind(dtl, data.table::as.data.table(coldata[match(dtl[[names_to]], rownames(coldata)),]))
  colnames(dtout)<-make.unique(colnames(dtout))
  return(dtout)
}

