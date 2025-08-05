#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param divtab PARAM_DESCRIPTION
#' @param vars PARAM_DESCRIPTION
#' @param groupvars PARAM_DESCRIPTION
#' @param pairvar PARAM_DESCRIPTION
#' @param outdir PARAM_DESCRIPTION
#' @param name PARAM_DESCRIPTION, Default: 'alpha_diversity'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname testDiversityDifferences_paired
#' @export 
testDiversityDifferences_paired <- function(divtab, vars, groupvars, pairvar ,outdir, name="alpha_diversity"){
  library(car)
  divtab <- divtab[order(divtab[,pairvar]) , ]
  res <- data.frame()
  for (v in vars){
    for(g in groupvars){
      divtab[, g] <- as.factor(divtab[, g])
      divtab[, v] <- as.numeric(divtab[, v])
      gr_a <-divtab[divtab[, g] == levels(divtab[, g])[1] ,v]
      gr_b <-divtab[divtab[, g] == levels(divtab[, g])[2] ,v]
      num_groups <- length(unique(divtab[!is.na(divtab[, g]), g]))
      form2<- as.formula(paste0(v, " ~ ", g, " + Error(", pairvar, "/", g, ")"))
      form <- as.formula(paste0(v, " ~ ", g))
      aovres <- aov(form2, divtab)
      sum_aovres <- summary(aovres)[[2]] %>% unlist
      tres <- tryCatch(ifelse(num_groups==2, t.test(gr_a, gr_b, paired = T)$p.value, NA), error=function(x)return(NA))
      wres <- tryCatch(ifelse(num_groups==2, wilcox.test(gr_a, gr_b, paired = T)$p.value, NA), error=function(x)return(NA))
      swres <- shapiro.test(divtab[, v])$p.value
      bt <-  bartlett.test(form, divtab)[[3]]
      levt <-  leveneTest(form, divtab)[1, 3]
      aux <- data.frame(
        variable = v,
        groups = g,
        comparison = "all",
        anova_F =  sum_aovres["F value1"],
        anova_p = sum_aovres["Pr(>F)1"],
        t_test = tres,
        wilcox_test = wres,
        shapiro_normality_test = swres,
        bartlett_test = bt,
        levene_test = levt
      )
      res <- rbind(res, aux)
      if(num_groups > 2){
        aux_t <- getTestsForAllCombinations(divtab[, v], divtab[, g])
        aux2 <- data.frame(
          variable = v,
          groups = g,
          comparison = aux_t$groups_compared,
          anova_F = NA,
          anova_p = NA,
          t_test = aux_t$t_pval,
          wilcox_test = aux_t$wilcox_pval,
          shapiro_normality_test = aux_t$shapiro_test,
          bartlett_test = aux_t$bartlett_test,
          levene_test = aux_t$levene_test
        )
        res <- rbind(res, aux2)
      } ## if num_groups > 2
    } # For grouping variabbles
  } # For problem variables
  res$t_corrected <- p.adjust(res$t_test, method="BH")
  res$wilcox_corrected <- p.adjust(res$wilcox_test, method="BH")
  write_tsv(res, paste0(outdir, "/", name, "_QualitVarsTestsPaired.tsv"))
  return(res)
}
