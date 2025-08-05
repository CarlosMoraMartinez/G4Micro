process_names <- function(namelist){
  gsub("-", ".", namelist) %>%
    gsub("\\[|\\]|\\(|\\)", "", ., perl=T)
}
