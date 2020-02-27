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
#' @importFrom DelayedMatrixStats rowVars
#' @import data.table
#' @export


writeGSEAfiles<-function (object, class_labels, gene_labels="gene_short_name", gct_filename="data.gct", cls_filename="data.cls", normalized=TRUE) 
{
  #gct_filename<-file_path_as_absolute(gct_filename)
  #cls_filename<-file_path_as_absolute(cls_filename)
  if(class(object)=="DESeqDataSet"){
    object<-object[,order(colData(object)[[class_labels]])]
    exprMat<-counts(object, normalized=normalized)
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
  if(class(object)=="cell_data_set"){
    object<-object[,order(pData(object)[[class_labels]])]
    exprMat<-as.matrix(counts(object, normalized=normalized))
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
      exprMat<-exprMat[rownames(exprMat) %in% tokeep]
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
  message(paste0("GCT File written:\n", file_path_as_absolute(gct_filename), "\n"))
  nConds = length(class_labels)
  uniqueLabels = levels(class_labels)
  write(paste(nConds, length(uniqueLabels), "1", sep = " "), 
        file = cls_filename, append = F)
  write(paste("#", paste(uniqueLabels, collapse = " "), sep = " "), 
        file = cls_filename, append = T)
  write(as.numeric(class_labels)-1, file = cls_filename, append = T, 
        sep = " ", ncolumns = nConds)
  message(paste0("CLS File written:\n", file_path_as_absolute(cls_filename), "\n"))
}



#' GSEA and plot
#' @description Performs GSEA (using fgsea) and returns a GSEA-style enrichment plot.  Optionally create an \
#' enrichment plot with multiple layers that correspong to enrichnment of a given "pathway" in a list of ranked genes \
#' such that each entry in a list corresponds to a different enrichment plot.  Alternatively, one may supply a list of \
#' pathways to overlay the enrichment of multiple pathways on one ranked list.  Can't do both though....
#' @param pathway = vector (or list) of genes.  Names of list will define plot coloring.
#' @param stats = vector (or list) of ranked gene stats (usually fold change or SNR) with names \
#' that contain the gene name.  Names of list will define plot coloring.
#' @param rug = whether to make a rug plot(s).
#' @param dot_enhance character string denoting a color that enhances the dot appearance \
#'  with another color
#' @import ggplot2
#' @importFrom fgsea calcGseaStat
#' @param all_the_rest_of_them Should be self explanatory
#' @return Performs GSEA of "pathway" genes on stats'
#' @references fgsea package
#' @export

enrichmentPlot<-function (pathway, stats, 
                          gseaParam = 1,  rug=T, dot_enhance="darkgrey",
                          dot_enhance_size=2, dot_shape=21, 
                          dot_enhance_alpha=0.7, dot_size=1,
                          return_data=FALSE,
                          print_plot=FALSE,
                          return_plot=TRUE) 
{
  list_amenable_args<-c("stats", "pathway")
  for(arg in list_amenable_args){
    if(class(get(arg))!="list") {assign(arg, list(get(arg)))}
  }
  if(is.null(names(stats))){names(stats)=1:length(stats)}
  if(is.null(names(pathway))){names(pathway)=1:length(pathway)}
  stats_list_process<-function(stats, pathway){
    datalist<-lapply(1:length(stats), function(n){
      rnk <- rank(-stats[[n]])
      ord <- order(rnk)
      statsAdj <- stats[[n]][ord]
      statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
      statsAdj <- statsAdj/max(abs(statsAdj), na.rm = T)
      pw <- unname(as.vector(na.omit(match(pathway[[1]], names(statsAdj)))))
      pw <- sort(pw)
      gseaRes <- calcGseaStat(statsAdj, selectedStats = pw, 
                              returnAllExtremes = TRUE, returnLeadingEdge = TRUE)
      tops <- gseaRes$tops
      i <- length(statsAdj)
      xs <- as.vector(pw)
      ys <- as.vector(rbind(tops))
      le <- c(rep(1, length(xs)))
      le[xs %in% gseaRes$le]<-5
      le_bool<-le==5
      toPlot <- data.frame(x = c(0, xs, i + 1), y = c(0, ys, 0), le=c(0,le,0))
      #diff <- (max(ys) - min(ys))/6
      df_out<-data.frame(Rank = c(xs), ES = c(ys), LE=le_bool, row.names = names(tops))
      df_out$group<-names(stats)[n]
      toPlot$group<-names(stats)[n]
      list(df_out=df_out, toPlot=toPlot, gseaRes=gseaRes)
    })
    return(datalist)
  }
  pathway_list_process<-function(stats, pathway){
    datalist<-lapply(1:length(pathway), function(n){
      rnk <- rank(-stats[[1]])
      ord <- order(rnk)
      statsAdj <- stats[[1]][ord]
      statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
      statsAdj <- statsAdj/max(abs(statsAdj), na.rm = T)
      pw <- unname(as.vector(na.omit(match(pathway[[n]], names(statsAdj)))))
      pw <- sort(pw)
      gseaRes <- calcGseaStat(statsAdj, selectedStats = pw, 
                              returnAllExtremes = TRUE, returnLeadingEdge = TRUE)
      tops <- gseaRes$tops
      i <- length(statsAdj)
      xs <- as.vector(pw)
      ys <- as.vector(rbind(tops))
      le <- c(rep(1, length(xs)))
      le[xs %in% gseaRes$le]<-5
      le_bool<-le==5
      toPlot <- data.frame(x = c(0, xs, i + 1), y = c(0, ys, 0), le=c(0,le,0))
      #diff <- (max(ys) - min(ys))/6
      df_out<-data.frame(Rank = c(xs), ES = c(ys), LE=le_bool, row.names = names(tops))
      df_out$group<-names(pathway)[n]
      toPlot$group<-names(pathway)[n]
      list(df_out=df_out, toPlot=toPlot, gseaRes=gseaRes)
    })
    return(datalist)
  }
  
  if(length(stats)==1 & length(pathway) > 1){
    datalist<-pathway_list_process(stats, pathway)
  }else{
    datalist<-stats_list_process(stats, pathway)
  }
  plotList<-lapply(datalist, "[[", 2)
  maxs<-max(sapply(plotList, function(df) max(df$y)))
  mins<-min(sapply(plotList, function(df) min(df$y)))
  y_rug_counter<-mins - (maxs-mins)/6
  y_rug_diff<-(maxs-mins)/10
  for(i in 1:length(plotList)){
    plotList[[i]]$y_rug<-y_rug_counter
    y_rug_counter <- y_rug_counter - y_rug_diff
  }
  #find allgroup_min
  
  dataList<-lapply(datalist, "[[", 1)
  gseaRes<-lapply(datalist, "[[", 3)
  dataList<-lapply(dataList, function(df){
    df$marker<-rownames(df)
    rownames(df)<-NULL
    df
  })
  mergedPlot<-do.call(rbind, plotList)
  dataList<-do.call(rbind, dataList)
  
  x = y = NULL
  g <- ggplot(mergedPlot, aes(x = x, y = y, color=group))+ 
    geom_point(size = dot_size)+ 
    geom_hline(yintercept = maxs, linetype = "dashed")+ 
    geom_hline(yintercept = mins, linetype = "dashed")+ 
    geom_hline(yintercept = 0, colour = "black")+ 
    geom_line()+ 
    theme_bw()
  # if(rug) {g<-g+geom_segment(data = data.frame(x = pathway), mapping = aes(x = x, 
  #                                                                          y = min(ys)-diff/2, xend = x, yend = min(ys)-diff), size = 0.2, color=rug_color)}  
  if(!is.null(dot_enhance)) {g<-g+
    geom_point(color = dot_enhance, size = dot_enhance_size, shape=dot_shape, alpha=dot_enhance_alpha)}
  
  if(rug){ g<-g+ geom_point(data = mergedPlot, aes(x = x,
                                                   y = y_rug, size = le, color=group), shape = 124)}
  g<-g+theme(panel.border = element_blank(), panel.grid.minor = element_blank()) + 
    labs(x = "rank", y = "enrichment score")+ guides( size = FALSE)
  if(print_plot) print(g)
  if(return_data) return(list(plot=g, gseaRes=gseaRes, df_out=dataList))
  if(return_plot) return(g)
}


#' @keywords internal
compileStats<-function(gsa, gs=NULL){
  #this function collapses adjusted stats from a gsa (piano package) generated using gsea or fgsea
  stats<-data.frame(up=as.vector(gsa$pAdjDistinctDirUp), down=as.vector(gsa$pAdjDistinctDirDn))
  if(is.null(gs)){
    gs<-names(gsa$gsc)
  }
  stats[is.na(stats)]<-1
  stats<-1-stats
  dir<-colnames(stats)[max.col(stats, ties.method = "first")]
  stats<-abs(stats-1)
  stats$dir<-dir
  stats$name<-names(gsa$gsc)
  return( stats[stats$name %in% gs, ])
}

#' @export
returnFDR<-function(gsa, gs=NULL){
  #this function returns a stat line describing the FDR and direction of a gsa generated using gsea or fsea
  if(class(gsa)=="GSAres" & length(gs)==1){
    dat<-compileStats(gsa=gsa, gs=gs)
    return(paste0("FDR = ", round(dat[[dat$dir]], 6), " in the ",dat$dir,  " direction"))
  }
  if(class(gsa)=="list"){
    ru<-lapply(1:length(gsa), function(num){
      dat<-compileStats(gsa=gsa[[num]], gs=gs)
      lineout<-paste0(names(gsa)[num], " <= ",round(dat[[dat$dir]], 4))
    })
    return(paste0("FDR: ", paste0(unlist(ru), collapse="; ")))
  }
  if(length(gs)>1){
    ru<-lapply(1:length(gs), function(num){
      dat<-compileStats(gsa=gsa, gs=gs[num])
      lineout<-paste0(gs[num], " <= ",round(dat[[dat$dir]], 4))
    })
    return(paste0("FDR: ", paste0(unlist(ru), collapse="; ")))
  }
}
