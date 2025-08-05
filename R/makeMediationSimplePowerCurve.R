#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param df PARAM_DESCRIPTION
#' @param xname PARAM_DESCRIPTION
#' @param yname PARAM_DESCRIPTION
#' @param medname PARAM_DESCRIPTION
#' @param param_estimates PARAM_DESCRIPTION, Default: c(a = 0.1, b = 0.1, cp = 0.1)
#' @param nboot PARAM_DESCRIPTION, Default: 1000
#' @param nrep PARAM_DESCRIPTION, Default: 1000
#' @param nboot_curve PARAM_DESCRIPTION, Default: 100
#' @param nrep_curve PARAM_DESCRIPTION, Default: 100
#' @param min_n PARAM_DESCRIPTION, Default: 100
#' @param max_n PARAM_DESCRIPTION, Default: 500
#' @param inc_n PARAM_DESCRIPTION, Default: 20
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[bmemLavaan]{power.bmem}}
#'  \code{\link[dplyr]{mutate}}
#' @rdname makeMediationSimplePowerCurve
#' @export 
#' @importFrom bmemLavaan power.bmem
#' @importFrom dplyr mutate
makeMediationSimplePowerCurve <- function(df, xname, yname, medname, param_estimates = c(a=0.1, b=0.1, cp=0.1),
                                          nboot=1000, nrep=1000, nboot_curve=100, nrep_curve=100,
                                          min_n=100, max_n=500, inc_n=20){
  library(bmemLavaan)

  eq1 = paste0(yname, " ~ b*start(", param_estimates["b"], ")*", medname, " + cp*start(", param_estimates["a"], ")*", xname)
  eq2 = paste0(medname, " ~ a*start(", param_estimates["a"], ")*", xname)
  eq3 = "ab := a*b"
  #eq4 = "total := cp + (a*b)"
  eqs <- paste(eq1, eq2, eq3, sep="\n", collapse="\n")

  eq1 = paste0(medname, " ~ a*", xname)
  eq2 = paste0(yname, " ~ b*", medname, " + cp*", xname)
  eq3 = "ab:=a*b"
  eq4 = "total := cp + (a*b)"
  eqs <- paste(eq1, eq2, eq3, eq4, sep="\n", collapse="\n")

  pw <- power.basic(eqs, indirect = eq3, nobs = seq(100, 1000, 100), nrep=100)
  #eqs<-'
  #Condition_bin ~ b*start(5)*IMC_log + cp*start(0.1)*Actinomyces_naeslundii
  #IMC_log ~ a*start(0.2)*Actinomyces_naeslundii
  #ab := a*b

  power_res <- bmemLavaan::power.bmem(model = eqs,
                                      nobs = nrow(df),
                                      nrep = nrep,
                                      nboot = nboot,
                                      alpha=0.95, verbose = F)
  powsum <- summary(power_res)

  check_ssize <- seq(min_n, max_n, inc_n)
  pwc <- power.curve(model = eqs,
                     nobs=check_ssize,
                     method='normal', nrep=nrep_curve,
                     nboot=nboot_curve, alpha=.95,interactive=FALSE)
  pwcnames <- colnames(pwc)
  pwcmeanings <- c("b", "cp", "a", "Ycov", "Mcov", "Xcov", "Indir")
  names(pwcmeanings) <- pwcnames
  pwcdf <- pwc %>% as.data.frame %>%
    dplyr::mutate(N = check_ssize) %>%
    gather("param_name", "power", -N) %>%
    dplyr::mutate(param = pwcmeanings[param_name],
                  Taxon = xname)

  return(list(pwres = power_res, pwsum=powsum, pwcurve = pwcdf, estimates=param_estimates))

}
