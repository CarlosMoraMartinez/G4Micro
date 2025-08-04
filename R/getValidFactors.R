getValidFactors <- function(df, max_nas=0.2, min_pct=0.1, max_classes=5, exclude=c()){
  is_factor <- function(x){
    (is.factor(x) |
       is.character(x) |
       (is.numeric(x) & length(unique(x)) <= max_classes))
  }
  meet_nas <- function(x){
    sum(is.na(x) | (x== "-"))/length(x) <= max_nas
  }
  is_balanced <- function(x){
    tt <- table(x)
    tt <- tt/sum(tt)
    min(tt) >= min_pct & (length(tt) > 1)
  }


  res <- data.frame(var = names(df),
                    is_factor = sapply(df, is_factor),
                    meet_nas = sapply(df, meet_nas),
                    is_balanced = sapply(df, is_balanced),
                    not_excluded = (!names(df) %in% exclude)
  ) %>% dplyr::mutate(is_valid = is_factor & meet_nas & is_balanced & not_excluded)
  return(res)
}
