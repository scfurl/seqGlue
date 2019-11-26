#' @export
plot3d<-function (obj, color_palette = NULL, size=10, alpha=1) {
  p <- plotly::plot_ly(obj, x = ~dim_1, 
                       y = ~dim_2, z = ~dim_3, type = "scatter3d", 
                       size = I(size), color = ~color, colors = color_palette, 
                       mode = "markers", alpha = I(alpha))
  p
}