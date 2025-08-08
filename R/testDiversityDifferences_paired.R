#' @title Test Differences in Alpha Diversity Indices with Paired Samples
#' @description
#' Tests for differences in alpha diversity indices (or other continuous variables) across
#' grouping variables while accounting for paired samples defined by a pairing variable.
#' Performs repeated measures ANOVA, paired t-tests, and paired Wilcoxon tests where applicable.
#' Also conducts diagnostic tests including Shapiro-Wilk normality, Bartlett’s and Levene’s
#' tests for homogeneity of variance. Results are saved as a TSV file and returned as a data frame.
#'
#' @param divtab A data frame containing alpha diversity indices (or continuous numeric variables) as columns, along with grouping and pairing variables.
#' @param vars A character vector specifying the column names in \code{divtab} that contain the continuous variables to test (e.g., diversity indices).
#' @param groupvars A character vector specifying the column names in \code{divtab} that contain the grouping variables (factors) across which differences will be tested.
#' @param pairvar A single character string specifying the column name in \code{divtab} that identifies paired samples (e.g., subject ID).
#' @param outdir A character string specifying the path to the output directory where the results TSV file will be saved.
#' @param name A character string prefix for the output filename (default: \code{"alpha_diversity"}).
#'
#' @return A data frame containing the results of the statistical tests for each combination
#' of variable and grouping variable. Includes repeated measures ANOVA F-statistics and p-values,
#' paired t-test and paired Wilcoxon p-values (when applicable), Shapiro-Wilk normality test,
#' Bartlett’s and Levene’s tests for homogeneity of variance. Also includes Benjamini-Hochberg
#' corrected p-values for t-tests and Wilcoxon tests.
#'
#' @details
#' For each continuous variable in \code{vars} and each grouping variable in \code{groupvars}, the function:
#'
#' 1. Orders the data frame by the pairing variable to ensure paired samples align.
#' 2. Performs repeated measures ANOVA using the formula \code{variable ~ group + Error(pair/group)}.
#' 3. If the grouping variable has exactly two levels:
#'    - Performs a paired Student's t-test.
#'    - Performs a paired Wilcoxon signed-rank test.
#' 4. Applies diagnostic tests:
#'    - Shapiro-Wilk test for normality of the continuous variable.
#'    - Bartlett’s test for homogeneity of variances.
#'    - Levene’s test for homogeneity of variances.
#' 5. If the grouping variable has more than two levels:
#'    - Performs all pairwise comparisons between groups using t-tests and Wilcoxon tests (via the helper function \code{getTestsForAllCombinations}).
#'
#' The results are saved as a TSV file named \code{<name>_QualitVarsTestsPaired.tsv} in the specified output directory, and also returned as a data frame.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   res <- testDiversityDifferences_paired(
#'     divtab = my_diversity_data,
#'     vars = c("Shannon", "Chao1"),
#'     groupvars = c("Treatment"),
#'     pairvar = "SubjectID",
#'     outdir = "results/",
#'     name = "paired_diversity_stats"
#'   )
#' }
#' }
#'
#' @seealso
#' \code{\link[stats]{aov}},
#' \code{\link[stats]{t.test}},
#' \code{\link[stats]{wilcox.test}},
#' \code{\link[stats]{shapiro.test}},
#' \code{\link[stats]{bartlett.test}},
#' \code{\link[car]{leveneTest}},
#' \code{\link[stats]{p.adjust}}
#'
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
