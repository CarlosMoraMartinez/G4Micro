testDiversityWithQuantVars <- function(divtab, vars, qvars,
                                       outdir,
                                       name = "AlphaDiv_quant_vars_regressions"){
  library(broom)
  res <- data.frame()

  getLm <- function(v, g){
    aux <- divtab %>% dplyr::mutate_at(c(v, g), scale )

    form <- as.formula(paste0(v, " ~ ", g))
    tryCatch(mod <- lm(form, data=aux), error=function(x)cat("Error fitting ", as.character(form)))
    return(mod)
  }

  df <- expand.grid(variable=vars, predictor=qvars) %>%
    dplyr::mutate_all(as.character) %>%
    dplyr::mutate(model = map2(variable, predictor,getLm)) %>%
    dplyr::mutate(model2 = map(model, broom::glance)) %>%
    dplyr::mutate(intercept = sapply(model, function(mod) summary(mod)$coefficients[1] ),
                  slope = sapply(model, function(mod) summary(mod)$coefficients[2] ),
                  std_error_intercept = sapply(model, function(mod) summary(mod)$coefficients[3] ),
                  std_error_slope = sapply(model, function(mod) summary(mod)$coefficients[4] )
    ) %>%
    unnest(model2)
  fname <- paste0(outdir, "/", name, ".tsv")
  write_tsv(df, file=fname)
  return(df)
}
