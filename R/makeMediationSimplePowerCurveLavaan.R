
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param df PARAM_DESCRIPTION
#' @param xname PARAM_DESCRIPTION
#' @param yname PARAM_DESCRIPTION
#' @param medname PARAM_DESCRIPTION
#' @param med_res PARAM_DESCRIPTION
#' @param nrep PARAM_DESCRIPTION, Default: 1000
#' @param min_n PARAM_DESCRIPTION, Default: 100
#' @param max_n PARAM_DESCRIPTION, Default: 500
#' @param inc_n PARAM_DESCRIPTION, Default: 20
#' @param error PARAM_DESCRIPTION, Default: 0.1
#' @param ncores PARAM_DESCRIPTION, Default: 12
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
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
#' @importFrom dplyr mutate summarise
makeMediationSimplePowerCurveLavaan <- function(df, xname, yname, medname, med_res,
                                                nrep=1000,
                                                min_n=100,
                                                max_n=500,
                                                inc_n=20,
                                                error=0.1,
                                                ncores=12){
  library(parallel)
  library(doParallel)
  library(foreach)

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

                           #res <- foreach(numrep = 1:nrep, .combine = bind_rows,
                           #               .packages = c("dplyr", "tidyr", "magrittr", "tibble", "bmem", "sem")) %dopar% {

                           res <- data.frame()
                           for(numrep in 1:nrep){
                             X <- sample(valsX, ssize, replace=TRUE)
                             #M = med_res["a", "Estimate"]*X + rnorm(n = ssize , mean = 0, sd = sd(mod1res))
                             M = med_res["a", "Estimate"]*X + sample(mod1res, ssize, replace=TRUE)

                             logit_Y <- med_res["cp", "Estimate"] * X + med_res["b", "Estimate"] * M
                             p_Y <- 1 / (1 + exp(-logit_Y))
                             Y <- rbinom(ssize, 1, p_Y)

                             #simdf <- data.frame(
                             #  !!xname := X,
                             #  !!medname := M,
                             #  !!yname := Y
                             #)

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
                           }

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
                               power = sum(p.value < 0.05)/n()
                             )
                           ressum
                         }
  stopCluster(cl)
  return(ressum_all)

}
