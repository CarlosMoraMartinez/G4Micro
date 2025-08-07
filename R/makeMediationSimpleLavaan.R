#' @title Run a Mediation Analysis with lavaan
#' @description
#' Performs a simple mediation analysis using the `lavaan` package. Estimates the direct, indirect, and total effects
#' of a predictor variable on an outcome variable via a mediator, using structural equation modeling (SEM).
#'
#' @param df Data frame. The dataset containing the variables used in the mediation model.
#' @param xname Character. Name of the independent (predictor) variable.
#' @param yname Character. Name of the dependent (outcome) variable.
#' @param medname Character. Name of the mediator variable.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{\code{sem}}{A fitted \code{lavaan} SEM object.}
#'   \item{\code{estimates}}{A data frame with the estimated path coefficients, standard errors, z-scores, and p-values
#'         for the direct effect (cp), indirect effect (a*b), and total effect (cp + a*b).}
#' }
#'
#' @details
#' The mediation model estimated is defined as:
#' \itemize{
#'   \item \code{medname ~ a * xname} (effect of independent variable on mediator)
#'   \item \code{yname ~ b * medname + cp * xname} (effect of mediator and independent variable on outcome)
#'   \item Indirect effect: \code{a*b}
#'   \item Total effect: \code{cp + a*b}
#' }
#' If the dependent variable (\code{yname}) has only one unique value in the data, the function returns a zero-filled
#' result with \code{p.value} set slightly above 0.05 to indicate non-significance.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   library(dplyr)
#'   df <- tibble(
#'     x = rnorm(100),
#'     m = 0.5 * x + rnorm(100),
#'     y = 0.6 * m + 0.3 * x + rnorm(100)
#'   )
#'   result <- makeMediationSimpleLavaan(df, "x", "y", "m")
#'   print(result$estimates)
#'   summary(result$sem)
#'  }
#' }
#' @seealso
#'  \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{select}}, \code{\link[dplyr]{filter}}
#'  \code{\link[lavaan]{sem}}
#' @rdname makeMediationSimpleLavaan
#' @export
#' @importFrom dplyr mutate select filter
#' @importFrom lavaan sem
makeMediationSimpleLavaan <- function(df, xname, yname, medname){
  library(lavaan)

  eq1 = paste0(medname, " ~ a*", xname)
  eq2 = paste0(yname, " ~ b*", medname, " + cp*", xname)
  eq3 = "ab := a*b"
  eq4 = "total := cp + (a*b)"

  if(length(unique(df[, yname])) < 2){
    ss <- data.frame(
      Estimate = numeric(5),
      S.E. = numeric(5),
      `z-score` = numeric(5),
      p.value = rep(0.05001, 5)
    )
    names(ss)[3] <- "z-score"
    rownames(ss) <- c("a", "b", "cp", "a*b", "cp+a*b")
    return(list(sem=list(), estimates=ss))
  }

  df <- df %>% dplyr::mutate(!!yname := ordered(!!sym(yname)))


  eqs <- paste(eq1, eq2, eq3, eq4, sep="\n", collapse="\n")

  fit <- lavaan::sem(eqs, data = df)
  ss <- summary(fit)$pe %>%
    dplyr::select(label, est, se, z, pvalue) %>%
    dplyr::filter(label %in% c("a", "b", "cp", "ab", "total")) %>%
    dplyr::mutate(label = ifelse(label == "ab", "a*b", label),
                  label = ifelse(label == "total", "cp+a*b", label)) %>%
    column_to_rownames("label")
  names(ss) <- c("Estimate", "S.E.", "z-score", "p.value")

  return(list(sem=fit, estimates=ss))

}
