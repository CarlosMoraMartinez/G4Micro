#'  @title Identify Significant Variables Across Multiple Classes
#' @description Evaluates each numeric variable for significant differences across multiple classes using linear models.
#' Returns model statistics and significance flags.
#' @param datasc A data frame with a column `class` (factor with 2 or more levels) and other numeric predictor variables. Must include a `sample` column which will be excluded.
#' @param levs Character vector indicating the levels of the `class` factor, in the desired order.
#' @param plim P-value threshold for determining statistical significance. Default: 0.05
#' @return A data frame with one row per variable including model p-value (`pval`), whether any class comparisons were significant (`any_sig`), and which coefficients were significant (`which_sig`), along with summary statistics from `broom::glance()`.
#' @details For each variable, a linear model is fitted with the formula `value ~ class`. If any class coefficient (vs baseline) is below the `plim` threshold, it is flagged as significant. The result includes model-level statistics (R-squared, F-statistic, etc.) and names of significant contrasts.
#' @examples
#' \dontrun{
#' if(interactive()){
#'   data <- data.frame(
#'     sample = paste0("S", 1:100),
#'     var1 = rnorm(100),
#'     var2 = rnorm(100),
#'     class = factor(sample(c("A", "B", "C"), 100, replace = TRUE))
#'   )
#'   get_signif_components_multiclass(data, levs = c("A", "B", "C"))
#' }
#' }
#' @seealso
#'  \code{\link[dplyr]{select}},
#'  \code{\link[dplyr]{filter}},
#'  \code{\link[dplyr]{rename}},
#'  \code{\link[dplyr]{arrange}},
#'  \code{\link[stats]{lm}},
#'  \code{\link[broom]{glance}}
#' @seealso
#'  \code{\link[dplyr]{select}}, \code{\link[dplyr]{filter}}, \code{\link[dplyr]{rename}}, \code{\link[dplyr]{arrange}}
#'  \code{\link[broom]{reexports}}
#' @rdname get_signif_components_multiclass
#' @export
#' @importFrom dplyr select filter rename arrange
#' @importFrom broom glance
get_signif_components_multiclass <- function(datasc, levs, plim=0.05){
  df <- datasc %>% dplyr::select(-sample) %>% dplyr::filter(!is.na(class))
  varnames <- names(df)[names(df)!="class"]
  res <- data.frame()
  ps <- c()
  for(i in varnames){
    df$aux <- df[, i]
    mod <- lm(aux ~ class, data=df)
    modsum <- summary(mod)

    any_sig <- any(modsum$coefficients[2:length(levs), 4] < plim)
    which_sig = which(modsum$coefficients[2:length(levs), 4] < 0.05) %>% names %>% paste(collapse="_")
    auxdf <- data.frame(var=i,
                        any_sig = any_sig,
                        which_sig = which_sig
    ) %>% cbind(broom::glance(mod))
    res <- rbind(res, auxdf)
  }
  return(res %>% dplyr::rename(pval = p.value) %>% dplyr::arrange(pval))
}
