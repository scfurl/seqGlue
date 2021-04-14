#' Plot genes from a matrix, DESeq2 object or limma object
#'
#' @param obj input object
#' @param gene gene input
#' @param intgroup factor to be used for grouping,
#' @param normalized normalize data (default is TRUE) 
#' @param transform transform data (default is FALSE)
#' @param main title
#' @param xlab label for x axis of plot
#' @param returnData return data instead of plot 
#' @param replaced = FALSE
#' @param theme_opts use theme_opts (default is T)
#' @param annotate_data default F
#' @param gene_short_name column containin gene symbol; default = "gene_short_name"
#' @param ENSID column containing gene id (ENSEMBL) default = "GENEID"
#' @param pc pseudocount
#' @param annotation_col the column of meta data to use for annotation (default ="sample")
#' @param point_size default 3
#' @param colors colors to use
#' @export
#' 
plotGene<-function(obj, ...){
  if(class(obj)=="DESeqDataSet"){
    return(plotGeneDESeq2(obj=obj, ...))
  }
  if(class(obj)=="EList"){
    return(plotGeneLimma(obj=obj, ...))
  }
  if(class(obj)=="matrix"){
    return(plotGeneMatrix(obj=obj, ...))
  }
}

#' @import ggplot2
#' @import DESeq2
#' @export
#' 
plotGeneMatrix<-function(obj=NULL, gene=NULL, intgroup = NULL, normalized = TRUE, 
                         transform = FALSE, main, xlab = "group", returnData = FALSE, 
                         replaced = FALSE, theme_opts=T,
                         annotate_data=F, annotation_col="sample", point_size=3, colors=NULL) {
  logxy <- if (transform)
    "y"
  else ""
  if (missing(main)) {
    main <- if (is.numeric(gene)) {
      rownames(obj)[gene]
    }
    else {
      gene
    }
  }
  ylab <- ifelse(normalized, "normalized count", "count")
  if (returnData) 
    #return(data.frame(count = data$count, colData(obj)[intgroup]))
    return(data)
  
  # plot(data$group + runif(ncol(obj), -0.05, 0.05), data$count,
  #      xlim = c(0.5, max(data$group) + 0.5), log = logxy, xaxt = "n",
  #      xlab = xlab, ylab = ylab, main = main)
  # axis(1, at = seq_along(levels(group)), levels(group))
  g<-ggplot(data, aes(x=group, y=count, fill=group)) +
    geom_boxplot(width=0.4) + #geom_dotplot(binaxis="y", stackdir="center", method="histodot", binwidth=1/nrow(data)) +
    geom_point(shape=21, size=point_size)  +
    xlab(xlab)+
    ylab(ylab)+
    ggtitle(main)
  if(transform){g<-g+scale_y_log10() }
  if(theme_opts){g<-g+theme_opts()}
  if(!is.null(colors)){g<-g+scale_fill_manual(values=colors)}
  if(annotate_data){g<-g+
    geom_label_repel(aes(x=group, y=count, label =sample), fill="white",
                     fontface = 'bold', size=4,
                     box.padding = unit(0.35, "lines"),
                     point.padding = unit(0.3, "lines"))
  }
  g
}

#' @import ggplot2
#' @importFrom edgeR cpm
#' @export
plotGeneLimma<-function(obj=NULL, gene=NULL, intgroup = NULL, normalized = TRUE, 
                        transform = FALSE, main, xlab = "group", returnData = FALSE, 
                        replaced = FALSE, theme_opts=T,
                        annotate_data=F, annotation_col="sample", point_size=3, colors=NULL) {
  if(class(obj)=="EList"){
    cnts<-cpm(obj)[gene,]
  }
  if(is.null(intgroup))stop("Must specify a factor in 'intgroup'")
  
  data <- data.frame(count = cnts, group = intgroup)
  if(annotate_data){
    data[["sample"]]<-colData(obj)[[annotation_col]]
  }
  logxy <- if (transform)
    "y"
  else ""
  if (missing(main)) {
    main <- if (is.numeric(gene)) {
      rownames(obj)[gene]
    }
    else {
      gene
    }
  }
  ylab <- ifelse(normalized, "normalized count", "count")
  if (returnData) 
    #return(data.frame(count = data$count, colData(obj)[intgroup]))
    return(data)
  
  # plot(data$group + runif(ncol(obj), -0.05, 0.05), data$count,
  #      xlim = c(0.5, max(data$group) + 0.5), log = logxy, xaxt = "n",
  #      xlab = xlab, ylab = ylab, main = main)
  # axis(1, at = seq_along(levels(group)), levels(group))
  g<-ggplot(data, aes(x=group, y=count, fill=group)) +
    geom_boxplot(width=0.4) + #geom_dotplot(binaxis="y", stackdir="center", method="histodot", binwidth=1/nrow(data)) +
    geom_point(shape=21, size=point_size)  +
    xlab(xlab)+
    ylab(ylab)+
    ggtitle(main)
  if(transform){g<-g+scale_y_log10() }
  if(theme_opts){g<-g+theme_opts()}
  if(!is.null(colors)){g<-g+scale_fill_manual(values=colors)}
  if(annotate_data){g<-g+
    geom_label_repel(aes(x=group, y=count, label =sample), fill="white",
                     fontface = 'bold', size=4,
                     box.padding = unit(0.35, "lines"),
                     point.padding = unit(0.3, "lines"))
  }
  g
}

#' @import ggplot2
#' @import DESeq2
#' @export
plotGeneDESeq2<-function(obj, gene, intgroup = "condition", normalized = TRUE, 
                         transform = TRUE, main, xlab = "group", returnData = FALSE, 
                         replaced = FALSE, gene_short_name="gene_short_name", ENSID="GENEID", pc, theme_opts=F,
                         annotate_data=F, annotation_col="sample", point_size=6, colors=NULL) {
  if(!is.null(gene_short_name)){
    main<-gene
    gene<-mcols(obj)[[ENSID]][mcols(obj)[[gene_short_name]] %in% gene][1]
  }
  stopifnot(length(gene) == 1 & (is.character(gene) | (is.numeric(gene) & 
                                                         (gene >= 1 & gene <= nrow(obj)))))
  if (!all(intgroup %in% names(colData(obj))))
    stop("all variables in 'intgroup' must be columns of colData")
  if (!returnData) {
    if (!all(sapply(intgroup, function(v) is(colData(obj)[[v]],
                                             "factor")))) {
      stop("all variables in 'intgroup' should be factors, or choose returnData=TRUE and plot manually")
    }
  }
  if (missing(pc)) {
    pc <- if (transform)
      0.5
    else 0
  }
  if (is.null(sizeFactors(obj)) & is.null(normalizationFactors(obj))) {
    obj <- estimateSizeFactors(obj)
  }
  cnts <- counts(obj, normalized = normalized, replaced = replaced)[gene,]
  data <- cbind( data.frame(count = cnts + pc),data.frame(colData(obj)))
  if(annotate_data){
    data[["sample"]]<-colData(obj)[[annotation_col]]
  }
  logxy <- if (transform)
    "y"
  else ""
  if (missing(main)) {
    main <- if (is.numeric(gene)) {
      rownames(obj)[gene]
    }
    else {
      gene
    }
  }
  ylab <- ifelse(normalized, "normalized count", "count")
  if (returnData) 
    #return(data.frame(count = data$count, colData(obj)[intgroup]))
    return(data)
  
  # plot(data$group + runif(ncol(obj), -0.05, 0.05), data$count,
  #      xlim = c(0.5, max(data$group) + 0.5), log = logxy, xaxt = "n",
  #      xlab = xlab, ylab = ylab, main = main)
  # axis(1, at = seq_along(levels(group)), levels(group))
  g<-ggplot(data, aes(x=group, y=count, fill=group)) +
    geom_boxplot(width=0.4) + #geom_dotplot(binaxis="y", stackdir="center", method="histodot", binwidth=1/nrow(data)) +
    geom_point(shape=21, size=point_size)  +
    xlab(xlab)+
    ylab(ylab)+
    ggtitle(main)
  if(transform){g<-g+scale_y_log10() }
  if(theme_opts){g<-g+theme_opts()}
  if(!is.null(colors)){g<-g+scale_fill_manual(values=colors)}
  if(annotate_data){g<-g+
    geom_label_repel(aes(x=group, y=count, label =sample), fill="white",
                     fontface = 'bold', size=4,
                     box.padding = unit(0.35, "lines"),
                     point.padding = unit(0.3, "lines"))
  }
  g
}



#' Plot genes from a DESeq2 object using a violin plot
#'
#' @param obj input object
#' @param gene gene input
#' @param intgroup factor to be used for grouping,
#' @param normalized normalize data (default is TRUE) 
#' @param transform transform data (default is FALSE)
#' @param main title
#' @param xlab label for x axis of plot
#' @param returnData return data instead of plot 
#' @param replaced = FALSE
#' @param theme_opts use theme_opts (default is T)
#' @param annotate_data default F
#' @param gene_short_name column containin gene symbol; default = "gene_short_name"
#' @param ENSID column containing gene id (ENSEMBL) default = "GENEID"
#' @param pc pseudocount
#' @param annotation_col the column of meta data to use for annotation (default ="sample")
#' @param point_size default 3
#' @param colors colors to use
#' @export
#' 
#' 
plotGeneViolin<-function (obj, gene, intgroup = "condition", normalized = TRUE, text_size=4,
               transform = TRUE, main, xlab = "group", returnData = FALSE, 
               replaced = FALSE, gene_short_name = "gene_short_name", ENSID = "gene_id", 
               pc, theme_opts = F, annotate_data = F, annotation_col = "sample", 
               point_size = 6, colors = NULL) 
{
  if (!is.null(gene_short_name)) {
    main <- gene
    gene <- mcols(obj)[[ENSID]][mcols(obj)[[gene_short_name]] %in% 
                                  gene][1]
  }
  stopifnot(length(gene) == 1 & (is.character(gene) | (is.numeric(gene) & 
                                                         (gene >= 1 & gene <= nrow(obj)))))
  if (!all(intgroup %in% names(colData(obj)))) 
    stop("all variables in 'intgroup' must be columns of colData")
  if (!returnData) {
    if (!all(sapply(intgroup, function(v) is(colData(obj)[[v]], 
                                             "factor")))) {
      stop("all variables in 'intgroup' should be factors, or choose returnData=TRUE and plot manually")
    }
  }
  if (missing(pc)) {
    pc <- if (transform) 
      0.5
    else 0
  }
  if (is.null(sizeFactors(obj)) & is.null(normalizationFactors(obj))) {
    obj <- estimateSizeFactors(obj)
  }
  cnts <- counts(obj, normalized = normalized, replaced = replaced)[gene, 
  ]
  data <- cbind(data.frame(count = cnts + pc), data.frame(colData(obj)))
  if (annotate_data) {
    data[["sample"]] <- colData(obj)[[annotation_col]]
  }
  logxy <- if (transform) 
    "y"
  else ""
  if (missing(main)) {
    main <- if (is.numeric(gene)) {
      rownames(obj)[gene]
    }
    else {
      gene
    }
  }
  ylab <- ifelse(normalized, "normalized count", "count")
  if (returnData) 
    return(data)
  g <- ggplot(data, aes_string(x = intgroup, y = "count", fill = intgroup)) + 
    geom_violin(scale = "width", trim = F) + geom_jitter(shape = 21, width = .1, size = point_size) + 
    xlab(xlab) + ylab(ylab) + ggtitle(main)
  if (transform) {
    g <- g + scale_y_log10()
  }
  if (theme_opts) {
    g <- g + theme_opts()
  }
  if (!is.null(colors)) {
    g <- g + scale_fill_manual(values = colors)
  }
  if (annotate_data) {
    g <- g + geom_label_repel(aes_string(x = intgroup, y = "count", label = annotation_col), 
                              fill = "white", fontface = "bold", size = text_size, box.padding = unit(0.35, 
                                                                                                      "lines"), point.padding = unit(0.3, "lines"))
  }
  g
}



