#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param df PARAM_DESCRIPTION
#' @param xname PARAM_DESCRIPTION
#' @param yname PARAM_DESCRIPTION
#' @param medname PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
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
