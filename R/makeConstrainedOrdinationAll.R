#' @title Generate Multiple Constrained Ordination Plots
#' @description
#' Performs constrained ordination on a phyloseq object for multiple sets of variables.
#' It calls either `makeConstrainedOrdinationSingleVar` for single variables or
#' `makeConstrainedOrdination` for multiple variables, generating ggplot objects for each.
#' The plots are saved to disk using the provided output directory and naming options.
#'
#' @param phobj A phyloseq object containing microbiome data and sample metadata.
#' @param opt List of options used by the `WriteManyPlots` function (e.g., graphical parameters).
#' @param name Character string for the base name used when saving plots. Default is `"CCA"`.
#' @param outdir Directory path to save the output plots. Default is `"CCA"`.
#' @param variableList A list of character vectors. Each vector specifies one or more variables
#' to include in the ordination formula. Default is `list(c("Condition"))`.
#' @param dist_type Distance metric to use for ordination (passed to `phyloseq::distance`). Default is `"bray"`.
#'
#' @return A named list of ggplot objects corresponding to each variable set ordination.
#'
#' @details
#' This function iterates over sets of variables provided in `variableList` and performs constrained ordination
#' using the appropriate function depending on the number of variables in each set.
#' It then writes all plots to files using the external `WriteManyPlots` function.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   # Assume physeq_obj is your phyloseq object and opt is your plotting options list
#'   variable_sets <- list(c("Condition"), c("Treatment", "TimePoint"))
#'   plots <- makeConstrainedOrdinationAll(physeq_obj, opt, name="MyOrdination", outdir="results", variableList=variable_sets)
#'   print(plots[[1]])
#' }
#' }
#'
#' @rdname makeConstrainedOrdinationAll
#' @export
makeConstrainedOrdinationAll <- function(phobj, opt, name="CCA", outdir= "CCA", variableList = list(c("Condition")),  dist_type="bray"){
  plots <- list()
  i = 1
  for(vars in variableList){
    if(length(vars) ==1 ){
      plots[[i]] <- makeConstrainedOrdinationSingleVar(phobj, variables =vars,  dist_type=dist_type)
    }else{
      plots[[i]] <- makeConstrainedOrdination(phobj, variables =vars,  dist_type=dist_type)
    }
    i <- i+1
  }
  names(plots)<- sapply(variableList, paste, sep="-", collapse="-")
  WriteManyPlots(plots, name, outdir, opt=opt)
  return(plots)
}
