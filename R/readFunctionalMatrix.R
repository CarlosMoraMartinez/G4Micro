
#' @title Read Functional Matrix from File
#' @description
#' Reads a tab-delimited functional matrix file from a specified directory,
#' processes the column names to remove suffixes after underscores,
#' and standardizes specific patterns in column names.
#'
#' @param opt A list containing at least the element \code{input_funcional},
#'   which specifies the directory path where the file is located.
#' @param fname A character string specifying the filename of the functional matrix to read.
#'
#' @return A tibble/data frame containing the functional matrix with cleaned column names.
#' The first column is renamed to "Pathway".
#'
#' @details
#' This function reads a tab-separated file using \code{read_delim} from the \code{readr} package.
#' Column names are truncated at the first underscore (_) and occurrences of "G4M0" are replaced with "G4M".
#' The first column is renamed to "Pathway" to indicate functional pathway names.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   opt <- list(input_funcional = "data/functional/")
#'   fname <- "functional_matrix.tsv"
#'   ftab <- readFunctionalMatrix(opt, fname)
#'   head(ftab)
#' }
#' }
#'
#' @rdname readFunctionalMatrix
#' @export
#' @importFrom readr read_delim
readFunctionalMatrix <- function(opt, fname){
  #ftab <- read.table(paste0(opt$input_funcional, fname), comment.char = "", sep="\t", head=T)
  ftab <- read_delim(paste0(opt$input_funcional, fname), delim="\t")
  names(ftab) <- sapply(names(ftab), function(x)strsplit(x, "_")[[1]][1]) %>% gsub("G4M0", "G4M", .)
  names(ftab)[1] <- "Pathway"
  return(ftab)
}
