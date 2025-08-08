#'  @title Detect Outliers in a Numeric Variable Using Grubbs' Test
#' @description
#' Iteratively identifies outliers in a numeric variable within a data frame
#' using Grubbs' test, returning the detected outlier rows.
#'
#' @param s_meta A data frame containing the variable to test.
#' @param vname Character string, name of the numeric variable/column to test for outliers.
#' @param p_lim Numeric, significance level threshold for Grubbs' test. Default is 0.05.
#'
#' @return A data frame with rows identified as outliers. Includes \code{sampleID} and the tested variable.
#'
#' @details
#' The function removes NA values for the tested variable, then repeatedly applies
#' Grubbs' test to detect the largest outlier until no outliers are detected or fewer
#' than 3 observations remain. Detected outliers are returned in a data frame.
#'
#' Requires the \code{outliers} package for \code{grubbs.test}.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   library(dplyr)
#'   library(outliers)
#'   df <- data.frame(sampleID = 1:10, value = c(1,2,2,3,3,4,4,100,5,6))
#'   detectOutliers(df, "value")
#' }
#' }
#'
#' @rdname detectOutliers
#' @export
#' @importFrom dplyr filter select
#' @importFrom rlang sym
#' @importFrom outliers grubbs.test
detectOutliers <- function(s_meta, vname, p_lim=0.05){
  library(outliers)
  aux <- s_meta %>%
    filter(!is.na(!!sym(vname)))
  outs <- data.frame()
  while(grubbs.test(aux[, vname])$p.value < p_lim & nrow(aux) > 2){
    outs <- rbind(outs,
                  aux %>% filter(!!sym(vname) == max(!!sym(vname) )) %>%
                    select(sampleID, all_of(vname))
    )
    aux <- aux %>% filter(!!sym(vname) < max(!!sym(vname) ))
  }
  return(outs)
}
