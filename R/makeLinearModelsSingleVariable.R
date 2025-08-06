#' @title makeLinearModelsSingleVariable
#' @description Linear models to test differences in alpha diversity indices.
#' @param divtab DataFrame containing metadata and alpha diversity indices, as returned by \code{calculateAlphaDiversityTable()}
#' @param interestvar Main variable to test for differences in alpha diversity indices (for instance, sex or disease status).
#' @param extravars Vector of covariates to include in the analysis. Also, differences in alpha diversity indices will be tested for these covariates.
#' @param alphaindices Alpha diversity indices present in the divtab DataFrame, Default: c("Observed", "Chao1", "Shannon", "InvSimpson")
#' @param combos Vector of integers. Linear models will be built with all the groups of covariates (with and without the main variable) of sizes included in this argument, Default: 1:3
#' @param outdir ath to the output directory where the results table will be saved, Default: ''
#' @param name Prefix for the output file names (TSV tables, .RData), Default: 'linearmodels'
#' @return A list with three elements:
#' `single_anovas`: a DataFrame with the result of testing differences for individual variables (both interestvar and extravars)
#' `anovas`: a DataFrame with the results of linear models including interstvar + groups of covariates, or only groups of covariates.
#' `models`
#' @details
#' This function performs linear modeling to assess differences in alpha diversity indices
#' across levels of a main variable of interest (e.g., disease status), optionally adjusting
#' for one or more covariates (e.g., age, sex). It fits two sets of models:
#'
#' 1. **Single-variable testing**: Each variable (both the main variable and covariates) is tested individually
#'    against each specified alpha diversity index using a simple linear model.
#'
#' 2. **Multiple-variable modeling**: For each alpha diversity index, the function builds multivariable linear models
#'    combining the main variable and groups of covariates. All combinations of covariates are generated based on the
#'    sizes specified in the `combos` argument (e.g., combinations of 1 to 3 covariates). For each model, the function
#'    compares a full model (main variable + covariates) against a reduced model (covariates only) using ANOVA.
#'
#' The output includes:
#' - A summary table with statistics from individual variable models (`single_anovas`).
#' - A summary table with statistics from multivariable models including covariates and the variable of interest (`anovas`).
#'
#' All models are returned in a named list (`models`), and results are saved to disk in both RData and TSV format
#' in the output directory.
#'
#' Variables with fewer than two levels are automatically excluded from the analysis.
#' P-values are adjusted using the Benjamini-Hochberg (BH) method across all models (`padj_all`)
#' and within models sharing the same variable structure (`padj_bymodel`).
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
#' @importFrom dplyr mutate
makeLinearModelsSingleVariable <- function(divtab,
                                           interestvar,
                                           extravars,
                                           alphaindices =c("Observed", "Chao1", "Shannon", "InvSimpson"),
                                           combos=1:3,
                                           outdir = "", name = "linearmodels" ){
  if(length(unique(divtab[, interestvar])) < 2) return(list()) #Error: only 1 level, not possible to fit model
  extravars <-  map_vec(divtab[, extravars], \(x)length(unique(x[!is.na(x)]))) %>%
    subset(. > 1) %>% names # remove variables without 2 or more levels

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
