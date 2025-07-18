################################################################################
#' @title RD Heterogeneous Treatment Effects.  Linear combinations of parameters
#'
#' @description \code{rdhte_lincom} computes point estimates, p-values, and
#' robust bias-corrected confidence intervals for linear combinations of
#' parameters after any estimation using \code{\link{rdhte}}
#' (Calonico, Cattaneo, Farrell, Palomba and Titiunik, 2025a).
#' Inference is implemented using robust bias-correction methods
#' (Calonico, Cattaneo, and Titiunik, 2014). It is based on the \code{R} function 
#' \code{\link[multcomp]{glht}}.
#'
#' Companion commands: \code{\link{rdhte}} for estimation and inference of RD-HTE
#' and \code{\link{rdbwhte}} for data-driven bandwidth selection.
#'
#' Related software packages for analysis and interpretation of RD designs and
#' related methods are available in: \url{https://rdpackages.github.io/}.
#'
#'For background methodology, see Calonico, Cattaneo, Farrell, and Titiunik
#'(2019), Calonico, Cattaneo and Farrell (2020), Cattaneo and Titiunik (2022).
#'
#'
#' @param model a fitted model returned by \code{\link{rdhte}}.
#' @param linfct a specification of the linear hypotheses to be tested. 
#' Linear functions can be specified by either the matrix of coefficients or by 
#' symbolic descriptions of one or more linear hypotheses.
#' @param level Confidence level for intervals; default is \code{level = 95}.
#' @param digits Number of decimal places to format numeric outputs (default 3).
#'
#' @return A list with two data frames: 'individual' and 'joint', with rounded values.
#'
#'
#' @author
#' Sebastian Calonico, University of California, Davis \email{scalonico@ucdavis.edu}.
#'
#' Matias D. Cattaneo, Princeton University  \email{cattaneo@princeton.edu}.
#'
#' Max H. Farrell, University of California, Santa Barbara \email{maxhfarrell@ucsb.edu}.
#'
#' Filippo Palomba, Princeton University \email{fpalomba@princeton.edu}.
#'
#' Rocio Titiunik, Princeton University \email{titiunik@princeton.edu}.
#'
#' @references
#' Calonico, Cattaneo, Farrell, Palomba and Titiunik (2025): rdhte: Learning Conditional Average Treatment Effects in RD Designs. \emph{Working paper}.
#'
#' Calonico, Cattaneo, Farrell, Palomba and Titiunik (2025): Treatment Effect Heterogeneity in Regression Discontinuity Designs. \emph{Working paper}
#'
#' @seealso \code{\link{rdhte}}, \code{\link{rdbwhte}}
#'
#' @examples
#' set.seed(123)
#' n <- 1000
#' X <- runif(n, -1, 1)
#' W <- rbinom(n, 1, 0.5)
#' Y <- 3 + 2*X + 1.5*X^2 + 0.5*X^3 + sin(2*X) + 3*W*(X>=0) + rnorm(n)
#' m1 <- rdhte(y = Y, x = X, covs.hte = factor(W))
#' linfct <- c("`factor(W)0` - `factor(W)1` = 0")
#' rdhte_lincom(model = m1, linfct = linfct)
#'
#' @export

rdhte_lincom <- function(model,
                                 linfct,
                                 level = 95,
                                 digits     = 3) {
  if (!requireNamespace("multcomp", quietly = TRUE)) {
    stop("Package 'multcomp' is required to use this function.")
  }
  
  coef    <- model$coef
  coef.bc <- model$coef.bc
  vcov    <- model$vcov
  
  # Create the glht object
  gh_est <- glht(model = model, linfct = linfct, coef = coef,    vcov = vcov)
  gh_inf <- glht(model = model, linfct = linfct, coef = coef.bc, vcov = vcov)

  # Summary for individual tests (no adjustment)
  sgh_est <- summary(gh_est)
  sgh_inf <- summary(gh_inf)
  
  # Extract individual results
  coefs  <- sgh_est$test$coefficients
  ses    <- sgh_inf$test$sigma
  tstats <- sgh_inf$test$tstat
  pvals  <- sgh_inf$test$pvalues
  
  # Confidence intervals
  #ci <- confint(gh_inf, level = conf.level)$confint
  qz <- -qnorm(abs((1-(level/100))/2))
  ci <- cbind(sgh_inf$test$coefficients- qz*sgh_inf$test$sigma, sgh_inf$test$coefficients + qz*sgh_inf$test$sigma)
  
  individual <- data.frame(
    hypothesis = names(coefs),
    estimate   = as.numeric(coefs),
    #se         = as.numeric(ses),
    t_stat     = as.numeric(tstats),
    p_value    = as.numeric(pvals),
    #conf.low   = ci[, "lwr"],
    #conf.high  = ci[, "upr"],
    conf.low   = ci[, 1],
    conf.high  = ci[, 2],
    row.names  = NULL
  )
  
  # Joint Wald chi-squared test
  #wtest <- multcomp::waldtest(gh)
  wtest <- summary(gh_inf, test = Chisqtest())
  joint <- data.frame(
    statistic.chi2 = wtest$test$SSH,
    df             = wtest$test$df[1],
    p_value        = wtest$test$pvalue,
    row.names      = NULL
  )
  
  # Round numeric columns
  individual[] <- lapply(individual, function(x) if (is.numeric(x)) round(x, digits) else x)
  joint[]      <- lapply(joint,      function(x) if (is.numeric(x)) round(x, digits) else x)
  
  list(individual = individual, joint = joint)
}
