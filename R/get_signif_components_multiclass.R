get_signif_components_multiclass <- function(datasc, levs, plim=0.05){
  df <- datasc %>% dplyr::select(-sample) %>% dplyr::filter(!is.na(class))
  varnames <- names(df)[names(df)!="class"]
  res <- data.frame()
  ps <- c()
  for(i in varnames){
    df$aux <- df[, i]
    mod <- mod <- lm(aux ~ class, data=df)
    modsum <- summary(mod)

    any_sig <- any(modsum$coefficients[2:length(levs), 4] < plim)
    which_sig = which(modsum$coefficients[2:length(levs), 4] < 0.05) %>% names %>% paste(collapse="_")
    auxdf <- data.frame(var=i,
                        any_sig = any_sig,
                        which_sig = which_sig
    ) %>% cbind(broom::glance(mod))
    res <- rbind(res, auxdf)
  }
  return(res %>% dplyr::rename(pval = p.value) %>% dplyr::arrange(pval))
}
