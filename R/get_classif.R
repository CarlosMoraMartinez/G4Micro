
#' @title Parse Taxonomic Classification String into a Data Frame
#' @description
#' Parses a taxonomic classification string formatted with hierarchical ranks separated by
#' pipe (`|`) and double underscore (`__`), and returns a one-row data frame with named taxonomic ranks.
#'
#' @param classtring Character. Taxonomic classification string, e.g. "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria".
#' @param classnames Named list or vector. Optional. Mapping of taxonomic rank prefixes to rank names.
#'        Default is a list with typical taxonomic levels:
#'        \itemize{
#'          \item d = "Kingdom"
#'          \item k = "Kingdom2"
#'          \item p = "Phylum"
#'          \item c = "Class"
#'          \item o = "Order"
#'          \item f = "Family"
#'          \item g = "Genus"
#'          \item s = "Species"
#'          \item xx = "Strain"
#'        }
#'
#' @return A one-row \code{data.frame} with columns named after taxonomic ranks and values parsed from the input string.
#' Missing ranks will be \code{NA}.
#'
#' @details
#' The function splits the input classification string by pipe and double underscore to extract rank prefixes and
#' associated taxon names. It handles the special case where multiple "k" ranks appear by converting the first "k" to "d".
#' The output data frame columns correspond to the provided or default rank names.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   clas <- "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__Escherichia|s__coli"
#'   get_classif(clas)
#' }
#' }
#'
#' @rdname get_classif
#' @export
get_classif <- function(classtring, classnames=NULL){

  if(is.null(classnames)){
    classnames <- list(d="Kingdom",
                       k="Kingdom2",
                       p="Phylum",
                       c="Class",
                       o="Order",
                       f="Family",
                       g="Genus",
                       s="Species",
                       xx="Strain")
  }

  classvec <- strsplit(classtring, "\\|")[[1]] %>%
    strsplit("__")
  classlist <- map(classvec, \(x)x[2]) %>% unlist
  class_init <- map(classvec, \(x)x[1]) %>% unlist
  if(length(which(class_init == "k")) > 1){
    class_init[which(class_init == "k")[1]] <- "d"
  }
  names(classlist) <- class_init
  #cat(classtring, "\n")

  aux <- data.frame(matrix(ncol=length(classnames), nrow=1, dimnames = list(NULL, unlist(classnames))))
  aux[1, unlist(classnames[class_init])] <- classlist[class_init]
  return(aux)
}
