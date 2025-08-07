
#' @title Plot SVM Classification Boundaries
#' @description
#' Generates a 2D visualization of a Support Vector Machine (SVM) classification model, including predicted regions, data points, and support vectors.
#' If the model is linear with two classes, it also plots the decision boundary and margins.
#'
#' @param modelo_svm A fitted SVM model object (e.g., from `e1071::svm`) with at least two variables.
#' @param datasc A data frame containing the PCA or input data used for classification.
#' @param varnames A character vector with the names of the variables to use in the plot (only the first two will be used).
#' @param opt A list-like object containing output options; must include `opt$out`, a directory where plots will be saved.
#' @param name A character string to identify the output file name.
#'
#' @return A `ggplot` object representing the SVM classification plot. The plot is also saved to a PDF.
#'
#' @details
#' This function was adapted from an example published at \url{https://rpubs.com/Joaquin_AR/267926}.
#' It performs interpolation across the feature space, applies the SVM model to the grid, and plots the decision boundaries.
#' If the kernel is linear and the problem is binary, the function also draws the decision hyperplane and margin lines.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   model <- e1071::svm(class ~ ., data = training_data, kernel = "linear")
#'   plotSVM(model, training_data, c("PC1", "PC2"), opt = list(out = "./results/"), name = "example")
#' }
#' }
#'
#' @seealso
#'  \code{\link[e1071]{svm}},
#'  \code{\link[stats]{predict}}
#'
#' @seealso
#'  \code{\link[dplyr]{select}}
#' @rdname plotSVM
#' @export
#' @importFrom dplyr select
plotSVM<-function(modelo_svm, datasc, varnames, opt, name){
  # Select and rename the predictor variables
  datos <- datasc[,varnames[1:2]]
  names(datos) <- paste("X", 1:ncol(datos), sep="")
  datos$y <- datasc$class

  # Compute ranges and generate new grid points for prediction
  rangos <- datos %>% dplyr::select(matches("^X[0-9]+")) %>% map(range)
  new_xs <- map(rangos, \(x)seq(from = x[1], to = x[2], length = 75))

  # Interpolate grid of points
  nuevos_puntos <- expand.grid(new_xs)

  # Predict class on the grid
  predicciones <- predict(object = modelo_svm, newdata = nuevos_puntos)

  # Store predicted points for coloring decision regions
  color_regiones <- data.frame(nuevos_puntos, y = predicciones)

  # # Compute hyperplane and margin (only for binary classification with linear kernel)
  beta <- drop(t(modelo_svm$coefs) %*% as.matrix(datos[, c("X1", "X2")])[modelo_svm$index,])
  beta0 <- modelo_svm$rho


  g1 <- ggplot() +
    geom_point(data = color_regiones, aes(x = X1, y = X2, color = as.factor(y)),
               size = 0.2, alpha=0.5) +
    geom_point(data = datos, aes(x = X1, y = X2, color = as.factor(y)),
               size = 2) +
    geom_point(data = datos[modelo_svm$index, ],
               aes(x = X1, y = X2, color = as.factor(y)),
               shape = 21, colour = "black",
               size = 2) +
    theme_bw() #+theme(legend.position = "none")

  if(modelo_svm$kernel == 0 & length(unique(datasc$class))==2){
    g1<-g1 + geom_abline(intercept = beta0/beta[2], slope = -beta[1]/beta[2]) +
      geom_abline(intercept = (beta0 - 1)/beta[2], slope = -beta[1]/beta[2],
                  linetype = "dashed") +
      geom_abline(intercept = (beta0 + 1)/beta[2], slope = -beta[1]/beta[2],
                  linetype = "dashed")
  }

  ggsave(filename = paste0(opt$out, name, "_SVMplot.pdf"), plot = g1, width = 8, height = 5)
  return(g1)
}
