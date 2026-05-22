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
#' A detailed introduction to the software is given in Calonico, Cattaneo,
#' Farrell, Palomba and Titiunik (2025b).
#'
#' Related software packages for analysis and interpretation of RD designs and
#' related methods are available in: \url{https://rdpackages.github.io/}.
#'
#' For background methodology, see Calonico, Cattaneo, Farrell, and Titiunik
#' (2019), Calonico, Cattaneo and Farrell (2020), and Cattaneo and Titiunik
#' (2022).
#'
#'
#' @param model a fitted model returned by \code{\link{rdhte}}.
#' @param linfct a specification of the linear hypotheses to be tested. 
#' Linear functions can be specified by either the matrix of coefficients or by 
#' symbolic descriptions of one or more linear hypotheses.
#' @param level Confidence level for intervals (percentage form, in \code{[1, 100)}); default is \code{level = 95}. Passing the fraction form (e.g. \code{0.95}) is rejected with a clear error.
#' @param digits Number of decimal places to format numeric outputs (default 3).
#'
#' @return A list with two data frames:
#' \describe{
#'   \item{\code{individual}}{One row per hypothesis. Columns: \code{hypothesis}, \code{estimate} (conventional point estimate of the linear combination), \code{z_stat} (asymptotic z-statistic from the bias-corrected fit), \code{p_value} (two-sided p-value from the standard normal), \code{conf.low}, \code{conf.high} (robust bias-corrected CI bounds at the requested confidence level).}
#'   \item{\code{joint}}{One row. Columns: \code{statistic} (Wald chi-squared from the bias-corrected fit), \code{df} (number of restrictions), \code{p_value}.}
#' }
#' Numeric columns are rounded to \code{digits} decimal places.
#'
#'
#' @author
#' Sebastian Calonico, University of California, Davis \email{scalonico@ucdavis.edu}.
#'
#' Matias D. Cattaneo, Princeton University  \email{matias.d.cattaneo@gmail.com}.
#'
#' Max H. Farrell, University of California, Santa Barbara \email{mhfarrell@gmail.com}.
#'
#' Filippo Palomba, Princeton University \email{filippo.palomba19@gmail.com}.
#'
#' Rocio Titiunik, Princeton University \email{rocio.titiunik@gmail.com}.
#'
#' @references
#' Calonico, Cattaneo, Farrell, Palomba and Titiunik (2025): \href{https://arxiv.org/abs/2507.01128}{rdhte: Conditional Average Treatment Effects in RD Designs.} \emph{Working paper}.
#'
#' Calonico, Cattaneo, Farrell, Palomba and Titiunik (2025): \href{https://arxiv.org/abs/2503.13696}{Treatment Effect Heterogeneity in Regression Discontinuity Designs.} \emph{Working paper}.
#'
#' Cattaneo and Titiunik. 2022. \href{https://rdpackages.github.io/references/Cattaneo-Titiunik_2022_ARE.pdf}{Regression Discontinuity Designs.} \emph{Annual Review of Economics},  14: 821-851.
#'
#' Calonico, Cattaneo, and Farrell. 2020. \href{https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell_2020_ECTJ.pdf}{Optimal Bandwidth Choice for Robust Bias Corrected Inference in Regression Discontinuity Designs.} \emph{Econometrics Journal}, 23(2): 192-210.
#'
#' Calonico, Cattaneo, Farrell, and Titiunik. 2019. \href{https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell-Titiunik_2019_RESTAT.pdf}{Regression Discontinuity Designs using Covariates.} \emph{Review of Economics and Statistics}, 101(3): 442-451.
#'
#' Calonico, Cattaneo, and Titiunik. 2014a. \href{https://rdpackages.github.io/references/Calonico-Cattaneo-Titiunik_2014_ECMA.pdf}{Robust Nonparametric Confidence Intervals for Regression-Discontinuity Designs.} \emph{Econometrica} 82(6): 2295-2326.
#'
#' Granzier, Pons, and Tricaud. 2023. \href{https://www.aeaweb.org/articles?id=10.1257/app.20210840}{Coordination and Bandwagon Effects: How Past Rankings Shape the Behavior of Voters and Candidates.} \emph{American Economic Journal: Applied Economics}, 15(4): 177-217.
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
#' \dontrun{
#' data(rdhte_dataset)
#' with(rdhte_dataset, {
#'   rd_ideology <- rdhte(y = y, x = x, covs.hte = factor(w_ideology),
#'                        cluster = cluster_var)
#'   rdhte_lincom(rd_ideology,
#'                linfct = c("`factor(w_ideology)4` - `factor(w_ideology)3` = 0",
#'                           "`factor(w_ideology)4` = 0"))
#' })
#' }
#'
#' @export

rdhte_lincom <- function(model,
                         linfct,
                         level  = 95,
                         digits = 3) {
  if (!requireNamespace("multcomp", quietly = TRUE)) {
    stop("Package 'multcomp' is required to use this function.")
  }
  ## Validate level. rdhte uses 0-100 convention; reject `level = 0.95`
  ## etc. that would silently give a tiny CI. Lower bound 1 catches the
  ## common fraction-instead-of-percent mistake.
  if (!is.numeric(level) || length(level) != 1 || level < 1 || level >= 100) {
    stop("`level` must be a single number in [1, 100); got ", deparse(level),
         ". Use the percentage form (e.g. 95), not the fraction form (0.95).",
         call. = FALSE)
  }

  coef    <- model$coef
  coef.bc <- model$coef.bc
  vcov    <- model$vcov

  # Create the glht object
  gh_est <- glht(model = model, linfct = linfct, coef = coef,    vcov = vcov)
  gh_inf <- glht(model = model, linfct = linfct, coef = coef.bc, vcov = vcov)

  ## When the model was fit with target.contrast, the bandwidth is optimal
  ## only for that contrast direction. Warn (once per call) if the user is
  ## now testing a hypothesis whose direction differs from the target.
  ## Same-direction calls (up to a positive scalar) stay silent.
  if (!is.null(model$target.contrast)) {
    Lmat <- gh_est$linfct
    target_K <- as.numeric(model$target.contrast)
    target_norm <- target_K / sqrt(sum(target_K^2))
    different <- TRUE
    if (is.matrix(Lmat) && ncol(Lmat) == length(target_K)) {
      ## Test if every row of Lmat is parallel to target_K (up to a scalar).
      different <- FALSE
      for (i in seq_len(nrow(Lmat))) {
        r <- as.numeric(Lmat[i, ])
        rn <- sqrt(sum(r^2))
        if (rn < 1e-12) { different <- TRUE; break }
        if (abs(abs(sum((r / rn) * target_norm)) - 1) > 1e-8) {
          different <- TRUE; break
        }
      }
    }
    if (different) {
      warning(paste0(
        "Model was fit with target.contrast = (",
        paste(sprintf("%g", target_K), collapse = ", "),
        "); the requested hypothesis tests a different direction. ",
        "Inference uses the bandwidth optimized for the target contrast; ",
        "consider rdhte_contrast() to re-optimize the bandwidth for each ",
        "contrast separately."), call. = FALSE)
    }
  }

  # Per-hypothesis tests with raw (univariate) p-values from the standard
  # normal, not multcomp's default single-step multiplicity adjustment.
  # Matches the docstring contract and Stata's `lincom` per-test reporting.
  sgh_est <- summary(gh_est, test = adjusted(type = "none"))
  sgh_inf <- summary(gh_inf, test = adjusted(type = "none"))

  # Extract individual results. multcomp's `tstat` field is the asymptotic
  # standardized statistic; with the default Gaussian linfct it's a z-stat.
  coefs  <- sgh_est$test$coefficients
  zstats <- sgh_inf$test$tstat
  pvals  <- sgh_inf$test$pvalues

  # Confidence intervals (BC + Gaussian)
  qz <- qnorm(1 - (1 - level/100) / 2)
  ci <- cbind(sgh_inf$test$coefficients - qz*sgh_inf$test$sigma,
              sgh_inf$test$coefficients + qz*sgh_inf$test$sigma)

  individual <- data.frame(
    hypothesis = names(coefs),
    estimate   = as.numeric(coefs),
    z_stat     = as.numeric(zstats),
    p_value    = as.numeric(pvals),
    conf.low   = ci[, 1],
    conf.high  = ci[, 2],
    row.names  = NULL
  )
  
  # Joint Wald chi-squared test
  wtest <- summary(gh_inf, test = Chisqtest())
  joint <- data.frame(
    statistic = wtest$test$SSH,
    df        = wtest$test$df[1],
    p_value   = wtest$test$pvalue,
    row.names = NULL
  )
  
  # Round numeric columns
  individual[] <- lapply(individual, function(x) if (is.numeric(x)) round(x, digits) else x)
  joint[]      <- lapply(joint,      function(x) if (is.numeric(x)) round(x, digits) else x)
  
  list(individual = individual, joint = joint)
}
