#' @title Power Analysis for Mediation Models Using bmemLavaan
#' @description Performs power analysis for mediation models over a range of sample sizes, returning a power curve and summary statistics.
#' @param df A data frame containing the variables used in the mediation analysis.
#' @param xname Character. Name of the independent variable (X).
#' @param yname Character. Name of the dependent variable (Y).
#' @param medname Character. Name of the mediator variable (M).
#' @param param_estimates Named numeric vector with starting values for parameters a, b, and cp. Default: \code{c(a = 0.1, b = 0.1, cp = 0.1)}.
#' @param nboot Integer. Number of bootstrap samples for initial power estimation. Default: 1000.
#' @param nrep Integer. Number of repetitions for initial power estimation. Default: 1000.
#' @param nboot_curve Integer. Number of bootstrap samples for power curve estimation. Default: 100.
#' @param nrep_curve Integer. Number of repetitions for power curve estimation. Default: 100.
#' @param min_n Integer. Minimum sample size to consider in power curve. Default: 100.
#' @param max_n Integer. Maximum sample size to consider in power curve. Default: 500.
#' @param inc_n Integer. Step size between sample sizes in power curve. Default: 20.
#' @return A list containing:
#' \item{pwres}{Raw result from \code{power.bmem}.}
#' \item{pwsum}{Summary of power results.}
#' \item{pwcurve}{A data frame with power estimates for each parameter across sample sizes.}
#' \item{estimates}{Input parameter estimates used in the simulation.}
#' @details This function uses the `bmemLavaan` package to estimate statistical power for the indirect effect (\code{a*b}) and other parameters in a simple mediation model. It returns both summary statistics and a full power curve over a specified range of sample sizes.
#' @examples
#' \dontrun{
#' if(interactive()){
#'   makeMediationSimplePowerCurve(mydata, xname = "X", yname = "Y", medname = "M")
#' }
#' }
#' @seealso
#'  \code{\link[bmemLavaan]{power.bmem}},
#'  \code{\link[dplyr]{mutate}}
#' @rdname makeMediationSimplePowerCurve
#' @export
#' @importFrom bmemLavaan power.bmem power.curve
#' @importFrom dplyr mutate
makeMediationSimplePowerCurve <- function(df, xname, yname, medname, param_estimates = c(a=0.1, b=0.1, cp=0.1),
                                          nboot=1000, nrep=1000, nboot_curve=100, nrep_curve=100,
                                          min_n=100, max_n=500, inc_n=20){

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

  #pw <- bmem::power.basic(eqs, indirect = eq3, nobs = seq(100, 1000, 100), nrep=100)

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
  pwc <- bmemLavaan::power.curve(model = eqs,
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
