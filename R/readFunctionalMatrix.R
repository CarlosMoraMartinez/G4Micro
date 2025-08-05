
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param opt PARAM_DESCRIPTION
#' @param fname PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname readFunctionalMatrix
#' @export 
readFunctionalMatrix <- function(opt, fname){
  #ftab <- read.table(paste0(opt$input_funcional, fname), comment.char = "", sep="\t", head=T)
  ftab <- read_delim(paste0(opt$input_funcional, fname), delim="\t")
  names(ftab) <- sapply(names(ftab), function(x)strsplit(x, "_")[[1]][1]) %>% gsub("G4M0", "G4M", .)
  names(ftab)[1] <- "Pathway"
  return(ftab)
}
