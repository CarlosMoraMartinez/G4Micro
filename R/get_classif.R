
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param classtring PARAM_DESCRIPTION
#' @param classnames PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
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
