#' @title Identify Significant Variables Using Logistic Regression
#' @description Performs univariate logistic regressions for each variable in a dataset to determine their association with a binary outcome. Returns performance metrics for each variable.
#' @param datasc A data frame with a column `class` (factor) for binary classification and other numeric predictor variables. Must also contain a column `sample` which will be excluded.
#' @param levs Character vector of length 2 indicating the two levels of the binary outcome variable. The second level is considered the positive class.
#' @return A data frame with one row per variable including the p-value from logistic regression and metrics: Accuracy, Sensitivity, Specificity, Positive Predictive Value (PPV), and Negative Predictive Value (NPV).
#' @details For each predictor variable, a logistic regression is fit using `glm`. Predictions are binarized with a threshold of 0.5 and compared to the true class to generate confusion matrices. Variables are ranked by p-value.
#' @examples
#' \dontrun{
#' if(interactive()){
#'   # Example usage:
#'   data <- data.frame(
#'     sample = paste0("S", 1:100),
#'     var1 = rnorm(100),
#'     var2 = rnorm(100),
#'     class = factor(sample(c("Control", "Case"), 100, replace = TRUE))
#'   )
#'   get_signif_components(data, levs = c("Control", "Case"))
#' }
#' }
#' @seealso
#'  \code{\link[dplyr]{select}},
#'  \code{\link[dplyr]{arrange}},
#'  \code{\link[stats]{glm}},
#'  \code{\link[stats]{predict}},
#'  \code{\link[caret]{confusionMatrix}}
#' @rdname get_signif_components
#' @export
#' @importFrom dplyr select arrange
#' @importFrom stats glm predict
#' @importFrom caret confusionMatrix
get_signif_components <- function(datasc, levs){
  df <- datasc %>% dplyr::select(-sample)
  varnames <- names(df)[names(df)!="class"]
  res <- data.frame()
  ps <- c()
  for(i in varnames){
    df$aux <- df[, i]
    mod_glm <- glm(class ~ aux, data=df, family = binomial)
    ps <- c(ps, summary(mod_glm)$coefficients[2, 4])
    predict_glm1 <- predict(mod_glm, df)
    predict_glm1 <- ifelse(predict_glm1 > 0.5, levs[2], levs[1]) %>% factor(levels=levs)
    confmat <- caret::confusionMatrix(predict_glm1, datasc$class, positive = levs[2])
    auxdf <- data.frame(var=i,
                        pval=summary(mod_glm)$coefficients[2, 4],
                        Accuracy=confmat$overall["Accuracy"],
                        Sensitivity = confmat$byClass["Sensitivity"],
                        Specificity = confmat$byClass["Specificity"],
                        PPV = confmat$byClass["Pos Pred Value"],
                        NPV = confmat$byClass["Neg Pred Value"]
    )
    res <- rbind(res, auxdf)
  }
  return(res %>% dplyr::arrange(pval))
}
