restauraropt_mk <- function(opt){
  output <- opt$out
  restaurar <- function(optbad){
    optbad$out <- output
    return(opt)
  }
}
