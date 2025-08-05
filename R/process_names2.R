process_names2 <- function(namelist){

  res <- gsub("_", " ", namelist) %>%
    sapply(\(x){
      if(x=="") return("")
      x <- strsplit(x, " ")[[1]]
      y <- toupper(strsplit(x[1], "")[[1]][1])
      y <- paste0(y, ". ", x[2])
      return(y)
    })
}
