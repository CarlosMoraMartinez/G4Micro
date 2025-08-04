testDiversityDifferences <- function(divtab, vars, groupvars, outdir, name="alpha_diversity"){
  library(car)
  res <- data.frame()
  for (v in vars){
    if(length(unique(divtab[, v])) < 2 ) next
    for(g in groupvars){
      if(length(unique(divtab[, g])) < 2 | min(table(divtab[!is.na(divtab[, v]), g]))<2) next
      meantab <- tapply(divtab[, v], divtab[, g], mean, na.rm=T)
      if(any(is.na(meantab))) next
      divtab[, g] <- as.factor(divtab[, g])
      divtab[, v] <- as.numeric(divtab[, v])
      num_groups <- length(unique(divtab[!is.na(divtab[, g]), g]))
      form <- as.formula(paste0(v, " ~ ", g))
      aovres <- aov(form, divtab)
      tres <- tryCatch(ifelse(num_groups==2, t.test(form, divtab)$p.value, NA), error=function(x)return(NA))
      wres <- tryCatch(ifelse(num_groups==2, wilcox.test(form, divtab)$p.value, NA), error=function(x)return(NA))
      swres <- shapiro.test(divtab[, v])$p.value
      bt <-  bartlett.test(form, divtab)[[3]]
      levt <-  leveneTest(form, divtab)[1, 3]
      aux <- data.frame(
        variable = v,
        groups = g,
        comparison = "all",
        anova_F = summary(aovres)[[1]]["F value"][1, 1],
        anova_p = summary(aovres)[[1]]["Pr(>F)"][1, 1],
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
  write_tsv(res, paste0(outdir, "/", name, "_QualitVarsTests.tsv"))
  return(res)
}
