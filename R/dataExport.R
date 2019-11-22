#export data

#'Export data to Excel file
#'
#' @description
#' This function takes data in various forms (DESeq, ExpressionSet and exports them to an Excel File)
#' @param object the object from which function will extract data.  Currently supported objects include:
#' 1) DESeq2 object and 2) Bioconductor expression set.
#' @param filename filename of Excel file
#' @param data_type can be one of normalized, raw, transformed.
#' @references Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and dispersion for
#' RNA-seq data with DESeq2 Genome Biology 15(12):550 (2014)
#' @references Orchestrating high-throughput genomic analysis with Bioconductor. W. Huber, V.J. Carey, 
#' R. Gentleman, ..., M. Morgan Nature Methods, 2015:12, 115.
#' @importFrom openxlsx createWorkbook
#' @importFrom openxlsx addWorksheet
#' @importFrom openxlsx writeData
#' @importFrom openxlsx saveWorkbook
#' @export


dataExport<-function (object, filename, data_type=c("normalized_counts", "raw_counts", "vst_transformed"), unique_id_column="id", symbol_name="gene_short_name") 
{
  if(length(data_type)>1 & all(data_type %in% c("normalized_counts", "raw_counts", "vst_transformed"))){
    data_type="normalized_counts"
  }
  if(class(object)=="DESeqDataSet"){
    pData<-colData(object)
    fData<-mcols(object)
    if(data_type=="normalized_counts"){
      exprMat<-counts(object, normalized=TRUE)
    }
    if(data_type=="raw_counts"){
      exprMat<-counts(object, normalized=FALSE)
    }
    if(data_type=="vst_transformed"){
      exprMat<-vst(object)@assays$data[[1]]
    }
    symbol<-fData[[symbol_name]]
    uid<-fData[[unique_id_column]]
    df<-data.frame(symbol=symbol, uid=uid)
    sexprMat<-format(exprMat, scientific = FALSE)
    bounddf<-cbind(df, sexprMat)
    chdf<-apply(bounddf, 2, as.character)
    pData_df<-data.frame(t(as.matrix(pData)))
    colnames(pData_df)<-colnames(exprMat)
    symbol_p<-rep("", nrow(pData_df))
    uid_p<-rep("", nrow(pData_df))
    empty_df<-data.frame(symbol=symbol_p, uid=uid_p)
    pData_r2m<-cbind(empty_df, pData_df)
    final_data<-rbind(pData_r2m, chdf)
    #filename<-"~/test.xlsx"
    wb<-openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, sheetName = "dataExport")
    openxlsx::writeData(x = final_data, wb = wb, sheet = "dataExport")
    openxlsx::saveWorkbook(wb, file = filename)
  }
  if(class(object)=="ExpressionSet"){
    pData<-pData(object)
    fData<-fData(object)
    if(data_type=="normalized_counts"){
      stop("normalized_counts not currently supported for ExpressionSet; use raw_counts")
    }
    if(data_type=="raw_counts"){
      exprMat<-exprs(object)
    }
    if(data_type=="vst_transformed"){
      stop("vst_transformed not currently supported for ExpressionSet; use raw_counts")
    }
    symbol<-fData[[symbol_name]]
    uid<-fData[[unique_id_column]]
    df<-data.frame(symbol=symbol, uid=uid)
    sexprMat<-format(exprMat, scientific = FALSE)
    bounddf<-cbind(df, sexprMat)
    chdf<-apply(bounddf, 2, as.character)
    pData_df<-data.frame(t(as.matrix(pData)))
    colnames(pData_df)<-colnames(exprMat)
    symbol_p<-rep("", nrow(pData_df))
    uid_p<-rep("", nrow(pData_df))
    empty_df<-data.frame(symbol=symbol_p, uid=uid_p)
    pData_r2m<-cbind(empty_df, pData_df)
    final_data<-rbind(pData_r2m, chdf)
    #filename<-"~/test.xlsx"
    wb<-openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, sheetName = "dataExport")
    openxlsx::writeData(x = final_data, wb = wb, sheet = "dataExport")
    openxlsx::saveWorkbook(wb, file = filename)
    }
}
