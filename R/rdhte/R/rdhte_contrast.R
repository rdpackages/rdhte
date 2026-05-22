################################################################################
#' @title RD Heterogeneous Treatment Effects: contrast-targeted bandwidth.
#'
#' @description \code{rdhte_contrast} refits an \code{\link{rdhte}} model at the
#' MSE-optimal bandwidth for a single user-chosen linear combination of the
#' Estimate vector and returns the resulting point estimate, robust
#' bias-corrected estimate, standard error, p-value, and confidence interval.
#'
#' The bandwidth is computed from the vector-first MSE expansion: starting
#' from the (V, B) pair already extracted by the pilot \code{\link{rdhte}}
#' fit, the optimal bandwidth for a contrast \code{contrast} is
#' \deqn{h^\star_\iota \;=\; \big(\,\iota' V \iota \,/\, (4\,(\iota' B)^2\,n)\,\big)^{1/5}}{h*_iota = (iota' V iota / (4 (iota' B)^2 n))^(1/5)}
#' where \eqn{\iota = Q' \cdot \mathtt{contrast}}{iota = Q' contrast} maps
#' the user-facing Estimate-row vector to a selector on
#' \eqn{\boldsymbol{\varsigma}=(\hat\theta, \hat\boldsymbol{\xi}')'}{varsigma}.
#' For continuous \code{covs.hte}, \eqn{Q}{Q} is the identity and
#' \code{contrast} acts directly on \eqn{\boldsymbol{\varsigma}}{varsigma}.
#'
#' Use \code{\link{rdhte_lincom}} instead when you want joint Wald inference
#' on several linear restrictions at one bandwidth.
#'
#' @param model A fitted model returned by \code{\link{rdhte}}.
#' @param contrast Numeric vector with length equal to \code{nrow(model$Estimate)}.
#'   The user-facing linear combination of Estimate rows whose MSE is optimized.
#' @param refit Logical. \code{TRUE} (default) refits the model at
#'   \eqn{h^\star_\iota}{h*_iota}. \code{FALSE} projects the existing fit
#'   (equivalent to \code{\link{rdhte_lincom}}) and returns the resulting
#'   estimate and SE under the model's current bandwidth.
#' @param level Confidence level (percentage form, in [1, 100)); default 95.
#' @param tol Tolerance below which |iota' B| is treated as zero and
#'   \eqn{h^\star_\iota}{h*_iota} is declared undefined. Default \code{1e-10}.
#'
#' @return A list with elements
#' \describe{
#'   \item{\code{estimate}}{Conventional point estimate \eqn{\iota' \hat\boldsymbol{\varsigma}}{iota' varsigma_hat}.}
#'   \item{\code{estimate.bc}}{Robust bias-corrected point estimate.}
#'   \item{\code{se.rb}}{Robust bias-corrected standard error.}
#'   \item{\code{z}}{Asymptotic z-statistic (\code{estimate.bc / se.rb}).}
#'   \item{\code{pv.rb}}{Two-sided p-value from the standard normal.}
#'   \item{\code{ci.rb}}{Confidence interval at the requested \code{level}.}
#'   \item{\code{h.star}}{The bandwidth used for the refit (or the current
#'     model bandwidth, if \code{refit = FALSE}).}
#'   \item{\code{h.pilot}}{Pilot bandwidth used to derive (V, B).}
#'   \item{\code{V.iota}, \code{B.iota}}{The scalar quadratic forms
#'     \eqn{\iota' V \iota}{iota' V iota} and \eqn{\iota' B}{iota' B}.}
#'   \item{\code{target.contrast}, \code{target.iota}}{The contrast in
#'     Estimate-row form and the implied varsigma-space selector.}
#'   \item{\code{refit}}{Logical indicating whether the model was refit.}
#'   \item{\code{model}}{When \code{refit = TRUE}, the refit \code{rdhte}
#'     object; \code{NULL} otherwise.}
#' }
#'
#' @seealso \code{\link{rdhte}}, \code{\link{rdhte_lincom}}.
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' n <- 1000
#' X <- runif(n, -1, 1)
#' W <- rbinom(n, 1, 0.5)
#' Y <- 3 + 2*X + 1.5*X^2 + 3*W*(X >= 0) + rnorm(n)
#' m <- rdhte(y = Y, x = X, covs.hte = factor(W))
#' # MSE-optimal bandwidth for the (W=1 minus W=0) effect difference:
#' rdhte_contrast(m, contrast = c(-1, 1))
#' }
#'
#' @export
rdhte_contrast <- function(model, contrast, refit = TRUE, level = 95, tol = 1e-10) {
  if (!inherits(model, "rdhte"))
    stop("'model' must be an rdhte fit.", call. = FALSE)
  if (is.null(model$W.lev) || length(model$Estimate) < 2L)
    stop("rdhte_contrast requires a heterogeneous fit (covs.hte with >= 2 Estimate rows).",
         call. = FALSE)
  K <- length(model$Estimate)
  if (!is.numeric(contrast) || any(!is.finite(contrast)) || length(contrast) != K)
    stop(sprintf("'contrast' must be a finite numeric vector of length %d.", K), call. = FALSE)
  if (all(contrast == 0))
    stop("'contrast' cannot be the zero vector.", call. = FALSE)
  if (!is.numeric(level) || length(level) != 1L || level < 1 || level >= 100)
    stop("'level' must be a single number in [1, 100); use percentage form (e.g. 95).",
         call. = FALSE)

  ## All downstream projections happen on Estimate-row space, where
  ## model$coef and model$vcov live. `contrast` is already in that
  ## coordinate system: bc' * Estimate = bc' * Q * varsigma = iota' * varsigma,
  ## so the two presentations are numerically identical.
  if (isTRUE(refit)) {
    ## Re-evaluate the original rdhte() call with target.contrast set.
    cl <- model$call
    if (is.null(cl))
      stop("Cannot refit: the model object does not retain its original call.", call. = FALSE)
    cl$target.contrast <- contrast
    cl$h <- NULL; cl$h.l <- NULL; cl$h.r <- NULL
    cl$level <- level
    new_model <- eval(cl, envir = parent.frame())

    iota_target <- new_model$target.iota
    proj_vec    <- contrast
    coef_used   <- new_model$coef
    coefbc_used <- new_model$coef.bc
    vcov_used   <- new_model$vcov
    h.star.used  <- new_model$h.star
    h.pilot.used <- new_model$h.pilot
    V.iota       <- as.numeric(crossprod(iota_target,
                                         new_model$V %*% iota_target))
    B.iota       <- as.numeric(crossprod(iota_target, new_model$B))
  } else {
    ## Pure projection on the model's current fit: same machinery, no
    ## refit. proj_vec is the Estimate-row contrast itself.
    proj_vec    <- contrast
    coef_used   <- model$coef
    coefbc_used <- model$coef.bc
    vcov_used   <- model$vcov
    h.star.used  <- if (!is.null(model$h.star)) model$h.star else NA_real_
    h.pilot.used <- if (!is.null(model$h.pilot)) model$h.pilot else NA_real_
    Q_proj <- if (isTRUE(model$covs.cont)) diag(K) else .rdhte_Q_factor(K)
    iota_target <- as.numeric(crossprod(Q_proj, contrast))
    V.iota <- as.numeric(crossprod(proj_vec, vcov_used %*% proj_vec))
    B.iota <- if (!is.null(model$B))
                as.numeric(crossprod(iota_target, model$B))
              else NA_real_
  }

  estimate     <- as.numeric(crossprod(proj_vec, coef_used))
  estimate.bc  <- as.numeric(crossprod(proj_vec, coefbc_used))
  se.rb        <- sqrt(as.numeric(crossprod(proj_vec, vcov_used %*% proj_vec)))
  z            <- estimate.bc / se.rb
  pv.rb        <- 2 * pnorm(-abs(z))
  qz           <- -qnorm(abs((1 - level / 100) / 2))
  ci.rb        <- c(estimate.bc - qz * se.rb, estimate.bc + qz * se.rb)

  out <- list(estimate        = estimate,
              estimate.bc     = estimate.bc,
              se.rb           = se.rb,
              z               = z,
              pv.rb           = pv.rb,
              ci.rb           = ci.rb,
              h.star          = h.star.used,
              h.pilot         = h.pilot.used,
              V.iota          = V.iota,
              B.iota          = B.iota,
              target.contrast = contrast,
              target.iota     = iota_target,
              level           = level,
              refit           = isTRUE(refit),
              model           = if (isTRUE(refit)) new_model else NULL)
  class(out) <- "rdhte_contrast"
  out
}

#' @keywords internal
#' @export
print.rdhte_contrast <- function(x, ...) {
  cat("rdhte contrast estimate\n")
  cat(sprintf("  contrast (Estimate rows) : %s\n",
              paste(sprintf("%.3f", x$target.contrast), collapse = ", ")))
  cat(sprintf("  refit                    : %s\n", ifelse(x$refit, "yes", "no")))
  if (!is.na(x$h.star))  cat(sprintf("  bandwidth (h*_iota)      : %.4f\n", x$h.star))
  if (!is.na(x$h.pilot) && x$refit)
    cat(sprintf("  pilot bandwidth          : %.4f\n", x$h.pilot))
  cat("\n")
  cat(sprintf("  estimate                 : %10.5f\n", x$estimate))
  cat(sprintf("  estimate.bc (RBC)        : %10.5f\n", x$estimate.bc))
  cat(sprintf("  se.rb                    : %10.5f\n", x$se.rb))
  cat(sprintf("  z-stat                   : %10.3f\n", x$z))
  cat(sprintf("  p-value                  : %10.4f\n", x$pv.rb))
  cat(sprintf("  %d%% CI                   : [%.5f, %.5f]\n",
              round(x$level), x$ci.rb[1], x$ci.rb[2]))
  invisible(x)
}
