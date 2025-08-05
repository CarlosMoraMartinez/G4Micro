#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param s_meta PARAM_DESCRIPTION
#' @param vname PARAM_DESCRIPTION
#' @param p_lim PARAM_DESCRIPTION, Default: 0.05
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname detectOutliers
#' @export 
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
