#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param mat1 PARAM_DESCRIPTION
#' @param mat2 PARAM_DESCRIPTION
#' @param group1 PARAM_DESCRIPTION
#' @param group2 PARAM_DESCRIPTION
#' @param outdir PARAM_DESCRIPTION
#' @param name PARAM_DESCRIPTION
#' @param levnames PARAM_DESCRIPTION, Default: c("PSO", "CTRL")
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
#' @rdname getCorsByGroup
#' @export 
#' @importFrom dplyr select
getCorsByGroup <- function(mat1, mat2, group1, group2, outdir, name, levnames=c("PSO", "CTRL")){

  cors <- getCorrelations(mat1, mat2)
  cors_pso <- getCorrelations(mat1[group1,], mat2[group1,])
  cors_ct <- getCorrelations(mat1[group2,], mat2[group2,])

  cors_pso_df <- cors_pso[[1]] %>% dplyr::select(-var1, -var2)
  names(cors_pso_df) <- paste(levnames[1], names(cors_pso_df), sep="_")

  cors_ct_df <- cors_ct[[1]] %>% dplyr::select(-var1, -var2)
  names(cors_ct_df) <- paste(levnames[2], names(cors_ct_df), sep="_")


  cordf <- cbind(cors[[1]], cors_pso_df, cors_ct_df)

  oname <- paste0(outdir, name)
  write_tsv(cordf, file = oname)
  return(list(alldf=cordf, all=cors, gr1=cors_pso, gr2=cors_ct))
}
