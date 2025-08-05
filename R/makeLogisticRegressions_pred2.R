
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param depvar PARAM_DESCRIPTION
#' @param posval PARAM_DESCRIPTION
#' @param df PARAM_DESCRIPTION
#' @param vars2test PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[stats]{glm}}
#'  \code{\link[broom]{reexports}}
#'  \code{\link[dplyr]{desc}}
#' @rdname makeLogisticRegressions_pred2
#' @export 
#' @importFrom stats glm
#' @importFrom broom glance
#' @importFrom dplyr desc
makeLogisticRegressions_pred2 <- function(depvar, posval, df, vars2test){
  df[, depvar] <- ifelse(df[, depvar]==posval, 1, 0)
  names(df) <- gsub("-", ".", names(df))
  vars2test <- gsub("-", ".", vars2test)
  mods_alone <- lapply(vars2test, FUN=function(var){
    form <- paste0(depvar, "~ ", var) %>% as.formula
    fm1 <- glm(form, data = df, family = binomial(link = "logit"))
    return(fm1)
  })
  mods <- lapply(vars2test, FUN=function(var){
    form <- paste0(depvar, "~ ", var) %>% as.formula
    fm1 <- stats::glm(form, data = df, family = binomial(link = "logit"))
    ss <- summary(fm1)
    pred <- ifelse(predict(fm1, type="response") > 0.5, 1, 0)
    confmat <- confusionMatrix(as.factor(pred), as.factor(df$Psoriasis))
    tab <- cbind(
      variable = var,
      broom::glance(fm1),
      coef_intercept=ss$coefficients[1, 1],
      coef_slope=ss$coefficients[2, 1],
      stderr_intercept=ss$coefficients[1, 2],
      stderr_slope=ss$coefficients[2, 2],
      zval_intercept=ss$coefficients[1, 3],
      zval_slope=ss$coefficients[2, 3],
      pval_intercept=ss$coefficients[1, 4],
      pval_slope=ss$coefficients[2, 4],
      Sensitivity=confmat$byClass["Sensitivity"],
      Specificity=confmat$byClass["Specificity"],
      PosPredValue=confmat$byClass["Pos Pred Value"],
      NegPredValue=confmat$byClass["Neg Pred Value"],
      Accuracy=confmat$overall["Accuracy"]
    )
    return(tab)
  }) %>% bind_rows()
  rownames(mods) <- names_cyt
  mods <- mods %>% arrange(dplyr::desc(Accuracy))
  return(list(mods_alone=mods_alone, mods=mods))
}
