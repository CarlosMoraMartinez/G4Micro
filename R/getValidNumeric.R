getValidNumeric <- function(df, max_nas=0.0, max_classes=4){
  is_number <- function(x){
    (is.numeric(x) &
       (length(unique(x)) >= max_classes))
  }
  meet_nas <- function(x){
    sum(is.na(x) | (x== "-"))/length(x) <= max_nas
  }

  res <- data.frame(var = names(df),
                    is_number = sapply(df, is_number),
                    meet_nas = sapply(df, meet_nas)
  ) %>% dplyr::mutate(is_valid = is_number & meet_nas)
  return(res)
}
