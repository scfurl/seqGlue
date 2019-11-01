#'Write GSEA files
#'
#' @description
#' GSEA formats are cumbersome.  This function writes working .GCT and .CLS files for use in a variety of GSEA
#' software including but not limited to the desktop (Java) application and the java cli.
#' @param object the object from which function will extract data.  Currently supported objects include:
#' 1) DESeq2 object and 2) Bioconductor expression set.
#' @param class_labels the meta data column name (as a string) whose data will be written into the .CLS file.
#' @param gene_labels the column name of the feature data that contains the desired name of .GCT rownames (ie genes; default is "gene_short_name").
#' @param gct_filename filename of the .GCT file.
#' @param cls_filename filename of the .CLS file.
#' @references Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and dispersion for
#' RNA-seq data with DESeq2 Genome Biology 15(12):550 (2014)
#' @references Orchestrating high-throughput genomic analysis with Bioconductor. W. Huber, V.J. Carey, 
#' R. Gentleman, ..., M. Morgan Nature Methods, 2015:12, 115.
#' @importFrom tools file_path_as_absolute
#' @importFrom matrixStats rowVars
#' @export


writeGSEAfiles<-function (object, class_labels, gene_labels="gene_short_name", gct_filename="data.gct", cls_filename="data.cls") 
{
  gct_filename<-file_path_as_absolute(gct_filename)
  cls_filename<-file_path_as_absolute(cls_filename)
  if(class(object)=="DESeqDataSet"){
    object<-object[,order(colData(object)[[class_labels]])]
    exprMat<-counts(object, normalized=TRUE)
    class_labels=colData(object)[[class_labels]]
    rn<-mcols(object)[[gene_labels]]
    
    if (length(make.unique(rn)) != length(unique(rn))) {
      unique_rn<-rownames(exprMat)[!rn %in% rn[duplicated(rn)]]
      duplicated_rn<-rownames(exprMat)[rn %in% rn[duplicated(rn)]]
      duplicated_lab<-rn[rn %in% rn[duplicated(rn)]]
      warning("Non unique gene_labels found in data; these have been collapsed keeping the label that has the greatest variance of normalized data")
      duplicated<-exprMat[rn %in% rn[duplicated(rn)],]
      ts = data.table(rowV = rowVars(duplicated), rn_D=duplicated_lab, rn_U=duplicated_rn)
      tokeep<-c(ts[ , .SD[which.max(rowV)], by = rn_D]$rn_U, unique_rn)
      rn<-rn[rownames(exprMat) %in% tokeep]
      exprMat<-exprMat[rownames(exprMat) %in% tokeep,]
    }
  }
  if(class(object)=="ExpressionSet"){
    object<-object[,order(pData(object)[[class_labels]])]
    exprMat<-exprs(object)
    class_labels=pData(object)[[class_labels]]
    rn<-fData(object)[[gene_labels]]
    if (length(make.unique(rn)) != length(unique(rn))) {
      unique_rn<-rownames(exprMat)[!rn %in% rn[duplicated(rn)]]
      duplicated_rn<-rownames(exprMat)[rn %in% rn[duplicated(rn)]]
      duplicated_lab<-rn[rn %in% rn[duplicated(rn)]]
      warning("Non unique gene_labels found in data; these have been collapsed keeping the label that has the greatest variance of normalized data")
      duplicated<-exprMat[rn %in% rn[duplicated(rn)],]
      ts = data.table(rowV = rowVars(duplicated), rn_D=duplicated_lab, rn_U=duplicated_rn)
      tokeep<-c(ts[ , .SD[which.max(rowV)], by = rn_D]$rn_U, unique_rn)
      rn<-rn[rownames(exprMat) %in% tokeep]
    }
  }
  nGenes = nrow(exprMat)
  nConds = ncol(exprMat)
  write("#1.2", file = gct_filename, append = F)
  write(paste(nGenes, nConds, sep = "\t"), file = gct_filename, append = T)
  if (length(make.unique(colnames(exprMat))) != length(unique(colnames(exprMat)))) {
    colnames(exprMat) <- make.unique(colnames(exprMat))
  }
  write(paste("Name", "Description", paste(colnames(exprMat), 
                                           collapse = "\t"), sep = "\t"), file = gct_filename, append = T)
  rownames(exprMat) = paste(rn, "na", sep = "\t")
  write.table(exprMat, file = gct_filename, append = T, quote = F, 
              sep = "\t", na = "", row.names = T, col.names = F)
  message(paste0("GCT File written:\n", gct_filename, "\n"))
  nConds = length(class_labels)
  uniqueLabels = levels(class_labels)
  write(paste(nConds, length(uniqueLabels), "1", sep = " "), 
        file = cls_filename, append = F)
  write(paste("#", paste(uniqueLabels, collapse = " "), sep = " "), 
        file = cls_filename, append = T)
  write(as.numeric(class_labels)-1, file = cls_filename, append = T, 
        sep = " ", ncolumns = nConds)
  message(paste0("CLS File written:\n", cls_filename, "\n"))
}
