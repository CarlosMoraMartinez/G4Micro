#' @title Linear Regression of Diversity Variables Against Quantitative Predictors
#' @description
#' Performs linear regression analyses of alpha diversity (or other numeric response variables)
#' against one or more quantitative predictor variables. Both response and predictor variables
#' are scaled prior to fitting. Summarizes regression statistics including intercepts, slopes,
#' standard errors, and overall model fit metrics (R-squared, p-values, AIC, etc.).
#' Saves results as a TSV file and returns a summary data frame.
#'
#' @param divtab A data frame containing numeric response variables and quantitative predictor variables.
#' @param vars A character vector specifying the names of the response variables (e.g., diversity indices) to regress.
#' @param qvars A character vector specifying the names of quantitative predictor variables.
#' @param outdir A character string specifying the directory path to save the output results file.
#' @param name A character string prefix for the output filename (default: \code{"AlphaDiv_quant_vars_regressions"}).
#'
#' @return A data frame summarizing regression results for each combination of response variable and predictor.
#' Includes estimates for intercept and slope, their standard errors, and overall model statistics
#' (e.g., R-squared, F-statistic, p-value, AIC).
#'
#' @details
#' The function scales both response and predictor variables before fitting linear models to
#' ensure comparability of coefficients. For each pair of response and predictor, a linear
#' model is fit, and tidy model statistics are extracted using the \code{broom} package.
#'
#' The output file will be saved as \code{<name>.tsv} in the specified \code{outdir}.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   results <- testDiversityWithQuantVars(
#'     divtab = diversity_data,
#'     vars = c("Shannon", "Simpson"),
#'     qvars = c("Age", "BMI"),
#'     outdir = "results"
#'   )
#' }
#' }
#'
#' @seealso
#' \code{\link[dplyr]{mutate_at}},
#' \code{\link[dplyr]{mutate}},
#' \code{\link[broom]{glance}},
#' \code{\link[stats]{lm}}
#'
#' @rdname testDiversityWithQuantVars
#' @export
#' @importFrom dplyr mutate_at mutate_all mutate
#' @importFrom broom glance
testDiversityWithQuantVars <- function(divtab, vars, qvars,
                                       outdir,
                                       name = "AlphaDiv_quant_vars_regressions"){
  library(broom)
  res <- data.frame()

  getLm <- function(v, g){
    aux <- divtab %>% dplyr::mutate_at(c(v, g), scale )

    form <- as.formula(paste0(v, " ~ ", g))
    tryCatch(mod <- lm(form, data=aux), error=function(x)cat("Error fitting ", as.character(form)))
    return(mod)
  }

  df <- expand.grid(variable=vars, predictor=qvars) %>%
    dplyr::mutate_all(as.character) %>%
    dplyr::mutate(model = map2(variable, predictor,getLm)) %>%
    dplyr::mutate(model2 = map(model, broom::glance)) %>%
    dplyr::mutate(intercept = sapply(model, function(mod) summary(mod)$coefficients[1] ),
                  slope = sapply(model, function(mod) summary(mod)$coefficients[2] ),
                  std_error_intercept = sapply(model, function(mod) summary(mod)$coefficients[3] ),
                  std_error_slope = sapply(model, function(mod) summary(mod)$coefficients[4] )
    ) %>%
    unnest(model2)
  fname <- paste0(outdir, "/", name, ".tsv")
  write_tsv(df, file=fname)
  return(df)
}
