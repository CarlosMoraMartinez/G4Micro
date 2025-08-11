
#' @title Estimate Mediation Model Power Curve Using Lavaan
#' @description Simulates data and evaluates the statistical power of a mediation model using increasing sample sizes and the `lavaan` package.
#'
#' @param df Data frame with original observed data to extract distributions.
#' @param xname Name (string) of the independent variable (X).
#' @param yname Name (string) of the dependent variable (Y), must be binary.
#' @param medname Name (string) of the mediator variable (M).
#' @param med_res A data frame with estimated mediation coefficients from `makeMediationSimpleLavaan()`.
#' @param nrep Number of simulations to perform per sample size. Default is 1000.
#' @param min_n Minimum sample size for power simulation. Default is 100.
#' @param max_n Maximum sample size for power simulation. Default is 500.
#' @param inc_n Increment step of sample sizes from min_n to max_n. Default is 20.
#' @param error Noise level used in mediator simulation. (Not currently used.) Default is 0.1.
#' @param plim_pow P-value threshold to consider the result as significant.
#' @param ncores Number of cores to use for parallel computing. Default is 12.
#'
#' @return A data frame summarizing the mean estimates, standard errors, z-scores, and power for each parameter at each sample size.
#'
#' @details This function simulates new datasets by sampling the independent variable and generating mediator and binary outcome values based on a logistic model using coefficients from a previously fitted mediation model. For each simulated dataset, a mediation analysis is repeated and the proportion of significant results (p < `plim_pow`) is used to estimate statistical power across a range of sample sizes.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{summarise}}
#' @rdname makeMediationSimplePowerCurveLavaan
#' @export
#' @importFrom dplyr mutate summarise group_by pull filter bind_rows
#' @importFrom tidyr pivot_longer
#' @importFrom tibble rownames_to_column
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
makeMediationSimplePowerCurveLavaan <- function(df, xname, yname, medname, med_res,
                                                nrep=1000,
                                                min_n=100,
                                                max_n=500,
                                                inc_n=20,
                                                error=0.1,
                                                plim_pow=0.05,
                                                ncores=10){

  cat("Doing Power Test for ", xname, "\n")
  cl <- makeCluster(ncores)
  registerDoParallel(cl)

  valsMed <- df %>% pull(!!sym(medname))
  valsX <- df %>% pull(!!sym(xname))

  ssizes <- seq(min_n, max_n, inc_n)

  eq1_lm <- paste0(medname, " ~ ", xname)
  mod1 = lm(as.formula(eq1_lm), data=df)
  mod1res <- resid(mod1)

  ressum_all <-  foreach(ssize = ssizes, .combine = bind_rows,
                         .export = c("makeMediationSimpleLavaan"),
                         .packages = c("dplyr", "tidyr", "magrittr", "tibble", "bmem", "sem", "lavaan")) %dopar%{

    cat("simulating ", xname, ", sample size=", ssize, "\n")

    res <- data.frame()
    for(numrep in 1:nrep){
       X <- sample(valsX, ssize, replace=TRUE)
       #M = med_res["a", "Estimate"]*X + rnorm(n = ssize , mean = 0, sd = sd(mod1res))
       M = med_res["a", "Estimate"]*X + sample(mod1res, ssize, replace=TRUE)
       logit_Y <- med_res["cp", "Estimate"] * X + med_res["b", "Estimate"] * M
       p_Y <- 1 / (1 + exp(-logit_Y))
       Y <- rbinom(ssize, 1, p_Y)

       simdf <- data.frame(
         X = X,
         M = M,
         Y = Y
       )
       colnames(simdf) <- c(xname, medname, yname)

       sim_sobelres <- tryCatch({

       makeMediationSimpleLavaan(simdf, xname = xname, yname = yname, medname = medname)$estimates %>%
         rownames_to_column("param") %>%
         filter(param %in% c("a", "b", "cp", "a*b", "cp+a*b ")) %>%
         dplyr::mutate(numrep = numrep)

        }, error = function(e) {
              warning(e)
              xx <- data.frame(param=character(0),
                               Estimate=numeric(0),
                               S.E.=numeric(0),
                               `z-score`=numeric(0),
                               p.value=numeric(0),
                               numrep=numeric(0))
                               names(xx)[4] <- "z-score"
              })

      res <- res %>% bind_rows(sim_sobelres)
    } ## for loop

    ressum <- res %>% group_by(param) %>%
      dplyr::mutate(
        p.value = ifelse(is.na(p.value), 1, p.value)
      ) %>%
      dplyr::summarise(
        ssize = ssize,
        mean_est = mean(Estimate, na.rm=T),
        mean_se = mean(`S.E.`, na.rm=T),
        mean_zscore = mean(`z-score`, na.rm=T),
        n = n(),
        power = sum(p.value <= plim_pow)/n()
      )
      ressum
  } ## foreach
  stopCluster(cl)
  return(ressum_all)

}
