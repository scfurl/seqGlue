#' Plot dispersion
#'
#' @description Many genes typically co-vary with one another, and so the dimensionality of the
#' data can be reduced with a wide variety of different algorithms. After calculating dispersion for
#' a using the \code{calc_dispersion} function, the 
#' \code{select_genes} function allows the user to identify a set of genes
#' that will be used in downstream dimensionality reduction methods.
#' @param obj the output of the \code{calc_dispersion} function.
#' @export
plot_gene_dispersion<-function(obj){
    prd<-obj$disp_table
    prd$fit<-log(obj$disp_func(prd$mu))
    prd$mu<-log(prd$mu)
    prd$disp<-log(prd$disp)
    colnames(prd)<-c("log_mean", "log_dispersion", "gene_id", "fit")
    g<-ggplot2::ggplot(prd, ggplot2::aes(x = log_mean, y = fit)) 
    if("selected_features" %in% names(obj)){
      prd$selected_features = obj$selected_features
      g <- g + ggplot2::geom_point(data=prd, ggplot2::aes(x=log_mean, y=log_dispersion, color=selected_features), alpha=0.4)
    }else{
      g <- g + ggplot2::geom_point(data=prd, ggplot2::aes( x=log_mean, y=log_dispersion, color="grey"), alpha=0.4)
    }
    g<-g+
      ggplot2::theme_bw() +
      ggplot2::geom_line(data=prd, ggplot2::aes( x=log_mean, y=fit)) +
      ggplot2::geom_smooth(data=prd, formula = fit ~ log_mean, stat = "identity") +
      labs(y="log_disp")
    return(g)
}

#' Select features based on mean expression and variance
#'
#' @description Many genes typically co-vary with one another, and so the dimensionality of the
#' data can be reduced with a wide variety of different algorithms. After calculating dispersion for
#' a using the \code{calc_dispersion} function, the 
#' \code{select_genes} function allows the user to identify a set of genes
#' that will be used in downstream dimensionality reduction methods.
#'
#'
#' @param obj the output of the \code{calc_dispersion} function.
#' @param fit_min the minimum multiple of the dispersion fit calculation; default = 1
#' @param fit_max the maximum multiple of the dispersion fit calculation; default = Inf
#' @param logmean_ul the maximum multiple of the dispersion fit calculation; default = Inf
#' @param logmean_ll the maximum multiple of the dispersion fit calculation; default = Inf
#' @param top top_n if specified, will override the fit_min and fit_max to select the top n most 
#' variant features.  logmena_ul and logmean_ll can still be used.
#' @return an updated object (list) that records the selected features
#' @export

select_genes<-function(obj, fit_min=1, fit_max=Inf, logmean_ul=Inf, logmean_ll=-Inf, top_n=NULL){
    df<-obj$disp_table
    df$fit<-obj$disp_func(df$mu)
    df$ratio<-df$disp/df$fit
    df$log_disp=log(df$disp)
    df$log_mean<-log(df$mu)
    df$index<-1:nrow(df)
    if(!is.null(top_n)){
      in_range<-df[which(df$log_mean > logmean_ll & df$log_mean < logmean_ul),]
      obj$selected_features <- df$index %in% in_range[order(-in_range$ratio),][1:top_n,]$index
    }else{
      obj$selected_features <- df$ratio > fit_min & df$ratio < fit_max & df$log_mean > logmean_ll & df$log_mean < logmean_ul
    }
    return(obj)
}


#' @export
get_selected_genes<-function(obj, gene_column="id"){
    as.character(obj$disp_table[[gene_column]][obj$selected_features])
}

#' Calculate gene dispersion of a matrix
#'
#' @description Many genes typically co-vary with one another, and so the dimensionality of the
#' data can be reduced with a wide variety of different algorithms. After calculating dispersion for
#' a using the \code{calc_dispersion} function, the 
#' \code{select_genes} function allows the user to identify a set of genes
#' that will be used in downstream dimensionality reduction methods.
#' @return a list containing 1) a table of dispersion values and 2) the function to fit
#' @export
#' @references Dispersion fitting taken from monocle2
#' @importFrom Matrix rowSums
#' @importFrom DelayedMatrixStats rowMeans2
#' @importFrom DelayedMatrixStats rowVars
#' @importFrom DESeq2 estimateSizeFactorsForMatrix
#' 
calc_dispersion<-function (obj, min_cells_detected=1, min_exprs = 1, id_tag="id", removeOutliers=TRUE) 
{
  if(class(obj) %in% "matrix"){
    sf<-DESeq2::estimateSizeFactorsForMatrix(obj)
    x <-DelayedArray(t(t(obj)/sf))
    f_expression_mean <- DelayedMatrixStats::rowMeans2(obj)
    f_expression_var <- DelayedMatrixStats::rowVars(obj)
    xim <- mean(1/sf)
    disp_guess_meth_moments <- f_expression_var - xim * f_expression_mean
    disp_guess_meth_moments <- disp_guess_meth_moments/(f_expression_mean^2)
    res <- data.frame(mu = as.vector(f_expression_mean), disp = as.vector(disp_guess_meth_moments))
    # res[res$mu == 0,]$mu = NA
    # res[res$mu == 0,]$disp = NA
    res$disp[res$disp < 0] <- 0
    res[[id_tag]] <- row.names(obj)
    disp_table <- subset(res, is.na(mu) == FALSE)
    res <- parametricDispersionFit(disp_table, verbose = T)
    fit <- res[[1]]
    coefs <- res[[2]]
    if (removeOutliers) {
      CD <- cooks.distance(fit)
      cooksCutoff <- 4/nrow(disp_table)
      message(paste("Removing", length(CD[CD > cooksCutoff]), 
                    "outliers"))
      outliers <- union(names(CD[CD > cooksCutoff]), setdiff(row.names(disp_table), 
                                                             names(CD)))
      res <- parametricDispersionFit(disp_table[row.names(disp_table) %in% 
                                                            outliers == FALSE, ], verbose=T)
      fit <- res[[1]]
      coefs <- res[[2]]
      names(coefs) <- c("asymptDisp", "extraPois")
      ans <- function(q) coefs[1] + coefs[2]/q
      attr(ans, "coefficients") <- coefs
    }
    res <- list(disp_table = disp_table, disp_func = ans)
    return(res)
  }
}

parametricDispersionFit<-function (disp_table, verbose = FALSE, initial_coefs = c(1e-06, 
                                                         1)) 
{
  coefs <- initial_coefs
  iter <- 0
  while (TRUE) {
    residuals <- disp_table$disp/(coefs[1] + coefs[2]/disp_table$mu)
    good <- disp_table[which(disp_table$disp > 0 & (residuals < 
                                                      10000)), ]
    if (verbose) 
      fit <- glm(disp ~ I(1/mu), data = good, family = Gamma(link = "identity"), 
                 start = coefs)
    else suppressWarnings(fit <- glm(disp ~ I(1/mu), data = good, 
                                     family = Gamma(link = "identity"), start = coefs))
    oldcoefs <- coefs
    coefs <- coefficients(fit)
    if (coefs[1] < initial_coefs[1]) {
      coefs[1] <- initial_coefs[1]
    }
    if (coefs[2] < 0) {
      stop("Parametric dispersion fit failed. Try a local fit and/or a pooled estimation. (See '?estimateDispersions')")
    }
    if (sum(log(coefs/oldcoefs)^2) < coefs[1]) 
      break
    iter <- iter + 1
    if (iter > 10) {
      warning("Dispersion fit did not converge.")
      break
    }
  }
  if (!all(coefs > 0)) {
    stop("Parametric dispersion fit failed. Try a local fit and/or a pooled estimation. (See '?estimateDispersions')")
  }
  list(fit, coefs)
}


