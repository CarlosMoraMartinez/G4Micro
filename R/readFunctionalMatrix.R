
readFunctionalMatrix <- function(opt, fname){
  #ftab <- read.table(paste0(opt$input_funcional, fname), comment.char = "", sep="\t", head=T)
  ftab <- read_delim(paste0(opt$input_funcional, fname), delim="\t")
  names(ftab) <- sapply(names(ftab), function(x)strsplit(x, "_")[[1]][1]) %>% gsub("G4M0", "G4M", .)
  names(ftab)[1] <- "Pathway"
  return(ftab)
}
