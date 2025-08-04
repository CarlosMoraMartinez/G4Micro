updatePsWithLogs <- function(phobj, vars = c("nreads")){

  for(v in vars){
    newname <- paste0(v, "_log")
    sample_data(phobj)[[newname]] <- log(sample_data(phobj)[[v]] +1  )
  }
  return(phobj)
}
