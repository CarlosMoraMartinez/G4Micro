defWriteMatAsDF <- function(mat, opt, name){
  fname <-paste(opt$out, name, sep="/", collapse="/")
  df <- mat %>% as.data.frame(row.names = rownames(.)) %>%
    rownames_to_column("gene")
  write_tsv(df, fname)
  return(df)
}
