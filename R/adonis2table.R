#' @title Format PERMANOVA and Dispersion Results into a Summary Table
#' @description
#' This function takes the output of \code{\link[vegan]{adonis2}}, and optionally
#' the results from \code{\link[vegan]{permutest}} on \code{\link[vegan]{betadisper}}
#' and \code{\link[vegan]{anova.cca}} on \code{\link[vegan]{capscale}}, to return a
#' summary table with relevant statistics for each tested variable.
#'
#' It is designed to support automated PERMANOVA workflows across multiple variables.
#'
#' @param mod1 A result object from \code{\link[vegan]{adonis2}}.
#' @param perm_disp (Optional) A result object from \code{\link[vegan]{permutest}} applied to a \code{\link[vegan]{betadisper}} model. Default: \code{NULL}
#' @param ancap (Optional) A result object from \code{\link[vegan]{anova.cca}} applied to a \code{\link[vegan]{capscale}} model. Default: \code{NULL}
#' @param var A string representing the name of the variable tested. Default: \code{""}
#' @param adonisby String passed to \code{\link[vegan]{adonis2}} in the \code{by} argument. Default: "terms"
#'
#' @return A \code{data.frame} with the following columns:
#' \itemize{
#'   \item \code{variable}: Name of the tested variable
#'   \item \code{DF_var}, \code{DF_Residual}, \code{DF_Total}: Degrees of freedom
#'   \item \code{SumOfSQs_var}, \code{SumOfSQs_Residual}, \code{SumOfSQs_Total}: Sum of squares
#'   \item \code{R2_var}, \code{R2_Residual}, \code{R2_Total}: R-squared values
#'   \item \code{F_statistic}, \code{P}: PERMANOVA F and p-values
#'   \item \code{perm_disp_F}, \code{perm_disp_P}, \code{perm_disp_npermuts}: Dispersion test results (if provided)
#'   \item \code{capscaleanopva_F}, \code{capscaleanova_P}: CAP scale ANOVA results (if provided)
#' }
#'
#' @details
#' This function is typically used inside a loop that runs PERMANOVA tests
#' across multiple metadata variables, to format and consolidate results for reporting.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   data("GlobalPatterns", package = "phyloseq")
#'   GP <- GlobalPatterns
#'   dist <- phyloseq::distance(GP, method = "bray")
#'   meta <- as.data.frame(sample_data(GP))
#'   mod <- vegan::adonis2(dist ~ SampleType, data = meta)
#'   adonis2table(mod, var = "SampleType")
#'  }
#' }
#' @rdname adonis2table
#' @export
adonis2table <- function(mod1, perm_disp=NULL, ancap=NULL, var="", adonisby=NULL){

  if(is.null(adonisby)){
    aux <- data.frame(
      variable = ifelse(var == "", rownames(mod1)[1], var),
      DF_var=mod1$Df[1],
      DF_Residual = mod1$Df[2],
      DF_Total = mod1$Df[3],
      SumOfSQs_var = mod1$SumOfSqs[1],
      SumOfSQs_Residual =mod1$SumOfSqs[2],
      SumOfSQs_Total = mod1$SumOfSqs[3],
      R2_var = mod1$R2[1],
      R2_Residual = mod1$R2[2],
      R2_Total = mod1$R2[3],
      F_statistic = mod1$F[1],
      P = mod1$`Pr(>F)`[1]
    )
  }else{
    aux <- data.frame(
      variable = rownames(mod1),
      DF_var=mod1$Df,
      DF_Residual = mod1$Df,
      DF_Total = mod1$Df,
      SumOfSQs_var = mod1$SumOfSqs,
      SumOfSQs_Residual =mod1$SumOfSqs,
      SumOfSQs_Total = mod1$SumOfSqs,
      R2_var = mod1$R2,
      R2_Residual = mod1$R2,
      R2_Total = mod1$R2,
      F_statistic = mod1$F,
      P = mod1$`Pr(>F)`
    )
  }
  if(!is.null(perm_disp)){
    aux$perm_disp_P = perm_disp$tab[1, "Pr(>F)"]
    aux$perm_disp_F = perm_disp$tab[1, "F"]
    aux$perm_disp_npermuts = perm_disp$tab[1, "N.Perm"]
  }else{
    aux$perm_disp_P = NA
    aux$perm_disp_F = NA
    aux$perm_disp_npermuts = NA
  }

  if(!is.null(ancap)){
    aux$capscaleanova_P = ancap[1, "Pr(>F)"]
    aux$capscaleanopva_F = ancap[1, "F"]
  }else{
    aux$capscaleanova_P = NA
    aux$capscaleanopva_F = NA
  }

  return(aux)
}
