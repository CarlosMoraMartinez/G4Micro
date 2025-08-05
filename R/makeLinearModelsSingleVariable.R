#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param divtab PARAM_DESCRIPTION
#' @param interestvar PARAM_DESCRIPTION
#' @param extravars PARAM_DESCRIPTION
#' @param alphaindices PARAM_DESCRIPTION, Default: c("Observed", "Chao1", "Shannon", "InvSimpson")
#' @param combos PARAM_DESCRIPTION, Default: 1:3
#' @param outdir PARAM_DESCRIPTION, Default: ''
#' @param name PARAM_DESCRIPTION, Default: 'linearmodels'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[base]{subset}}
#'  \code{\link[dplyr]{mutate}}
#' @rdname makeLinearModelsSingleVariable
#' @export 
#' @importFrom base subset
#' @importFrom dplyr mutate
makeLinearModelsSingleVariable <- function(divtab,
                                           interestvar,
                                           extravars,
                                           alphaindices =c("Observed", "Chao1", "Shannon", "InvSimpson"),
                                           combos=1:3,
                                           outdir = "", name = "linearmodels" ){
  if(length(unique(divtab[, interestvar])) < 2) return(list()) #Error: only 1 level, not possible to fit model
  extravars <-  map_vec(divtab[, extravars], \(x)length(unique(x[!is.na(x)]))) %>%
    base::subset(. > 1) %>% names # remove variables without 2 or more levels

  models <- list()
  anovas_singlevar <- data.frame()
  anovas <- data.frame()
  for(aind in alpha_indices){
    formulanull_ch <- paste0(aind, " ~ 1")
    formulanull <- as.formula(formulanull_ch)
    models[[formulanull_ch]] <- lm(formula = formulanull, data = divtab )
    # Variable por variable
    combolist <- c(interestvar, extravars)
    for(varg in combolist){
      formula_ch <- paste0(aind, " ~ ", varg)
      formula <- as.formula(formula_ch)
      models[[formula_ch]] <-  lm(formula = formula, data = divtab )
      auxanova <- anova(models[[formula_ch]]) %>% as.data.frame
      m1sum <- summary(models[[formula_ch]])$coefficients
      m2sum <- summary(models[[formulanull_ch]])$coefficients
      aux <- cbind(data.frame(nvars = 0,
                              Index=aind,
                              model=formula_ch,
                              reduced_model = formulanull_ch,
                              mod1= paste(rownames(m1sum), "=", round(m1sum[, "Estimate"], 2), "(p=", round(m1sum[, "Pr(>|t|)"], 4), ")", sep="", collapse="|"),
                              mod2=paste(rownames(m2sum), "=", round(m2sum[, "Estimate"], 2), "(p=", round(m2sum[, "Pr(>|t|)"], 4), ")", sep="", collapse="|")
      ),
      auxanova[1, ])
      anovas_singlevar <- rbind(anovas_singlevar, aux)
    }
    for(n_comb in combos){
      if(n_comb == 0)next
      combolist <- combn(extravars, n_comb, simplify = F)
      for(varg in combolist){
        formulared_ch <- paste0(aind, " ~ ", paste(varg, sep=" + ", collapse=" + "))
        formulared <- as.formula(formulared_ch)
        models[[formulared_ch]] <-  lm(formula = formulared, data = divtab )
        formula_ch <- paste0(aind, " ~ ", paste(varg, sep=" + ", collapse=" + "), " + ", interestvar)
        formula <- as.formula(formula_ch)
        models[[formula_ch]] <-  lm(formula = formula, data = divtab )

        auxanova <- anova(models[[formula_ch]],
                          models[[formulared_ch]]) %>% as.data.frame
        m1sum <- summary(models[[formula_ch]])$coefficients
        m2sum <- summary(models[[formulared_ch]])$coefficients
        aux <- cbind(data.frame(nvars = n_comb,
                                Index=aind,
                                model=formula_ch,
                                reduced_model = formulared_ch,
                                mod1= paste(rownames(m1sum), "=", round(m1sum[, "Estimate"], 2), "(p=", round(m1sum[, "Pr(>|t|)"], 4), ")", sep="", collapse="|"),
                                mod2=paste(rownames(m2sum), "=", round(m2sum[, "Estimate"], 2), "(p=", round(m2sum[, "Pr(>|t|)"], 4), ")", sep="", collapse="|")
        ),
        auxanova[2, ])
        anovas <- rbind(anovas, aux)
      }#combos
    }#combo length
  }#alpha dif measures
  anovas <- anovas %>% dplyr::mutate(
    padj_all = p.adjust(`Pr(>F)`, method="BH"),
    varsonly = sapply(strsplit(model, "~ "), \(x)x[2])
  ) %>% group_by(varsonly) %>%
    dplyr::mutate(padj_bymodel = p.adjust(`Pr(>F)`, method="BH")) %>%
    ungroup() %>%
    select(-varsonly)

  anovas_singlevar <- anovas_singlevar %>% dplyr::mutate(
    padj_all = p.adjust(`Pr(>F)`, method="BH"),
    varsonly = sapply(strsplit(model, "~ "), \(x)x[2])
  ) %>%
    group_by(varsonly) %>%
    dplyr::mutate(padj_bymodel = p.adjust(`Pr(>F)`, method="BH")) %>%
    ungroup() %>%
    select(-varsonly)

  results <- list(single_anovas =anovas_singlevar, anovas=anovas,models=models)
  save(results, file=paste0(outdir, "/", name, "linearmodels.RData"))
  write_tsv(anovas_singlevar, file=paste0(outdir, "/", name, "linearmodels_anovas_singlevar.tsv"))
  write_tsv(anovas, file=paste0(outdir, "/", name, "linearmodels_anovas_allvars.tsv"))
  return(results)
}
