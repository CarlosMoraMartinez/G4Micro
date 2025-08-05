
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param modelo_svm PARAM_DESCRIPTION
#' @param datasc PARAM_DESCRIPTION
#' @param varnames PARAM_DESCRIPTION
#' @param opt PARAM_DESCRIPTION
#' @param name PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[dplyr]{select}}
#' @rdname plotSVM
#' @export 
#' @importFrom dplyr select
plotSVM<-function(modelo_svm, datasc, varnames, opt, name){
  #Sacado de: https://rpubs.com/Joaquin_AR/267926
  datos <- datasc[,varnames[1:2]]
  names(datos) <- paste("X", 1:ncol(datos), sep="")
  datos$y <- datasc$class
  rangos <- datos %>% dplyr::select(matches("^X[0-9]+")) %>% map(range)
  new_xs <- map(rangos, \(x)seq(from = x[1], to = x[2], length = 75))

  # Interpolación de puntos
  nuevos_puntos <- expand.grid(new_xs)

  # Predicción según el modelo
  predicciones <- predict(object = modelo_svm, newdata = nuevos_puntos)

  # Se almacenan los puntos predichos para dar color a las regiones
  color_regiones <- data.frame(nuevos_puntos, y = predicciones)

  # Para extraer la ecuación del hiperplano y del margen es necesario aplicar
  # algebra lineal.
  beta <- drop(t(modelo_svm$coefs) %*% as.matrix(datos[, c("X1", "X2")])[modelo_svm$index,])
  beta0 <- modelo_svm$rho


  g1 <- ggplot() +
    # Representación de las 2 regiones empleando los puntos y coloreándolos
    # según la clase predicha por el modelo
    geom_point(data = color_regiones, aes(x = X1, y = X2, color = as.factor(y)),
               size = 0.2, alpha=0.5) +
    # Se añaden las observaciones
    geom_point(data = datos, aes(x = X1, y = X2, color = as.factor(y)),
               size = 2) +
    # Se identifican aquellas observaciones que son vectores soporte del modelo
    geom_point(data = datos[modelo_svm$index, ],
               aes(x = X1, y = X2, color = as.factor(y)),
               shape = 21, colour = "black",
               size = 2) +
    #scale_color_lancet()+
    #scale_fill_lancet() +
    theme_bw() #+theme(legend.position = "none")

  if(modelo_svm$kernel == 0 & length(unique(datasc$class))==2){
    # Se añaden las rectas del hiperplano y los márgenes
    g1<-g1 + geom_abline(intercept = beta0/beta[2], slope = -beta[1]/beta[2]) +
      geom_abline(intercept = (beta0 - 1)/beta[2], slope = -beta[1]/beta[2],
                  linetype = "dashed") +
      geom_abline(intercept = (beta0 + 1)/beta[2], slope = -beta[1]/beta[2],
                  linetype = "dashed")
  }

  ggsave(filename = paste0(opt$out, name, "_SVMplot.pdf"), plot = g1, width = 8, height = 5)
  return(g1)
}
