# Internal helpers for rdhte. None of these are exported.

# Resolve a set of unevaluated argument expressions against `data`,
# returning a named list. `mc` is `match.call()` from the calling
# function (so the bare expressions are visible); `caller` is the
# parent.frame() of the public-API call (for symbols not present in
# `data` -- e.g. user-defined helpers in their workspace).
#
# Mirrors rdrobust 4.0.0's .rdrobust_resolve_data.
.rdhte_resolve_data <- function(mc, data, caller, args) {
  out <- setNames(vector("list", length(args)), args)
  for (a in args) {
    if (a %in% names(mc) && !is.null(mc[[a]])) {
      out[[a]] <- eval(mc[[a]], envir = data, enclos = caller)
    }
  }
  out
}

# Pure-kernel constant k0 = e1' M_K^{-1} m (one-sided), where
#   M_K = [[mu_0, mu_1], [mu_1, mu_2]],  m = (mu_2, mu_3)',
# and mu_k := int_0^1 K(u) u^k du. This is the kernel-and-design
# constant that converts the auxiliary-quadratic-regression jump block
# into the asymptotic bias vector B of the rdhte estimator under the
# partially-linear-in-W structure of Assumption 2.
#
# Closed forms: triangular -> -1/10, uniform -> -1/6, Epanechnikov -> -11/95.
.rdhte_kernel_k0 <- function(kernel) {
  k <- tolower(kernel)
  mu_tri <- function(j) 1 / ((j + 1) * (j + 2))                  # K(u) = (1-|u|)
  mu_uni <- function(j) 1 / (j + 1)                              # K(u) = 1
  mu_epa <- function(j) (3/4) * (1/(j + 1) - 1/(j + 3))          # K(u) = (3/4)(1-u^2)
  mu_fun <- switch(k,
                   "tri"          = mu_tri,
                   "triangular"   = mu_tri,
                   "uni"          = mu_uni,
                   "uniform"      = mu_uni,
                   "epa"          = mu_epa,
                   "epanechnikov" = mu_epa,
                   stop("Unknown kernel: ", kernel, call. = FALSE))
  M  <- matrix(c(mu_fun(0), mu_fun(1), mu_fun(1), mu_fun(2)), 2L, 2L)
  mv <- c(mu_fun(2), mu_fun(3))
  as.numeric(solve(M, mv)[1L])
}

# Q maps Estimate rows to varsigma = (theta, xi_1, ..., xi_{d}).
# Continuous covs.hte: Estimate rows ARE varsigma entries, so Q = I_{1+d}.
# Factor covs.hte with K levels and the baseline level first:
#   Estimate[1] = theta,  Estimate[k] = theta + xi_{k-1}  for k = 2..K.
#   Q[1, ] = (1, 0_{K-1}); Q[k, ] = (1, e_{k-1}') for k = 2..K.
# Then iota = Q' c_tau converts a user-facing Estimate-row contrast c_tau
# into a varsigma-space selector iota.
.rdhte_Q_factor <- function(K) {
  if (K < 1L) stop("Q dimension K must be >= 1.", call. = FALSE)
  if (K == 1L) return(matrix(1, 1L, 1L))
  cbind(1, rbind(rep(0, K - 1L), diag(K - 1L)))
}

# Pull (V_hat, B_hat) from the two existing rdhte regression fits.
#  - rd_est : local-polynomial estimation fit (order p)
#  - rd_inf : augmented bias-correction fit (order q = p+1)
#  - vce.hte: type passed to sandwich::vcovCL ("HC0"/"HC1"/"HC2"/"HC3")
#  - cluster_sub: cluster vector restricted to the fit's subset (or NULL)
#  - varsigma_names: names of (theta, xi_*) in coef(rd_est) (length 1+d)
#  - bias_names    : names of the T*X^q jump block in coef(rd_inf) (length 1+d)
#  - n             : full pre-bandwidth N (denominator in MSE expansion)
#  - h_pilot       : the pilot bandwidth used for both fits (scalar)
#  - k0            : kernel constant from .rdhte_kernel_k0
.rdhte_extract_VB_from_fits <- function(rd_est, rd_inf, vce.hte, cluster_sub,
                                        varsigma_names, bias_names,
                                        n, h_pilot, k0) {
  Omega_full <- sandwich::vcovCL(rd_est, cluster = cluster_sub, type = vce.hte)
  if (!all(varsigma_names %in% rownames(Omega_full)))
    stop("Some varsigma coefficient names not found in fit; rdhte target-contrast helper aborted.",
         call. = FALSE)
  V_sample <- Omega_full[varsigma_names, varsigma_names, drop = FALSE]
  # V (asymptotic constant) = n * h * V_sample
  V_hat <- n * h_pilot * V_sample
  gamma_q <- coef(rd_inf)[bias_names]
  if (anyNA(gamma_q))
    stop("Bias-block coefficients contain NA; the augmented regression is rank-deficient at the pilot bandwidth.",
         call. = FALSE)
  B_hat <- k0 * as.numeric(gamma_q)
  list(V = V_hat, B = B_hat, V_sample = V_sample, gamma_q = as.numeric(gamma_q))
}

# Closed-form contrast-targeted MSE-optimal bandwidth.
#   h* = (iota' V iota / (4 * (iota' B)^2 * n))^{1/5}
# Errors out if |iota' B| is below tolerance (the formula diverges and
# the leading-order argument no longer applies).
.rdhte_hstar_iota <- function(V, B, iota, n, tol = 1e-10) {
  V_iota <- as.numeric(crossprod(iota, V %*% iota))
  B_iota <- as.numeric(crossprod(iota, B))
  if (V_iota <= 0)
    stop("Estimated iota' V iota <= 0; cannot compute contrast-targeted bandwidth.",
         call. = FALSE)
  if (abs(B_iota) < tol)
    stop(sprintf(paste0("Estimated bias for the requested contrast is below tolerance ",
                        "(|iota' B| = %.2e < %.0e); h*_iota is undefined. ",
                        "Pass `h` explicitly or pick a different contrast."),
                 abs(B_iota), tol), call. = FALSE)
  list(h_star  = (V_iota / (4 * B_iota^2 * n))^(1/5),
       V_iota  = V_iota,
       B_iota  = B_iota)
}

# Normalize the user-supplied `vce` argument to the canonical rdhte
# vocabulary, applying the documented cluster<->no-cluster remappings:
#   with cluster -> cr1 (default), cr2, cr3; hc0/hc1+cluster warn->cr1,
#                   hc2+cluster -> cr2, hc3+cluster -> cr3.
#   without cluster -> hc0/hc1/hc2/hc3 (default hc3); cr*+no-cluster
#                      warn->hc*; cr0+no-cluster warn->hc0.
# Returns list(vce = canonical_name, vce_label = display label).
#
# Caller is responsible for passing user_vce = !missing(vce), since
# missing() resolves in the immediate caller frame.
#
# Used by rdhte() and rdbwhte() to keep their vce parsing in one place.
.rdhte_normalize_vce <- function(vce, cluster, user_vce) {
  if (!is.character(vce) || length(vce) != 1L)
    stop("`vce` must be a single character string.", call. = FALSE)
  vce <- tolower(vce)
  if (!vce %in% c("hc0", "hc1", "hc2", "hc3", "cr0", "cr1", "cr2", "cr3"))
    stop(sprintf(
      "`vce` must be one of \"hc0\", \"hc1\", \"hc2\", \"hc3\", \"cr1\", \"cr2\", \"cr3\" (received: %s).",
      toString(vce)), call. = FALSE)

  if (!is.null(cluster)) {
    if (!isTRUE(user_vce)) {
      vce <- "cr1"  # silent default when cluster supplied
    } else if (vce %in% c("hc0", "hc1", "cr0")) {
      warning(paste0("vce='", vce, "' is not a cluster option. Switching to vce='cr1'."), call. = FALSE)
      vce <- "cr1"
    } else if (vce == "hc2") {
      warning("vce='hc2' is not a cluster option. Switching to vce='cr2'.", call. = FALSE)
      vce <- "cr2"
    } else if (vce == "hc3") {
      warning("vce='hc3' is not a cluster option. Switching to vce='cr3'.", call. = FALSE)
      vce <- "cr3"
    }
  } else {
    if (vce == "cr0") {
      warning("vce='cr0' requires a cluster variable. Falling back to vce='hc0'.", call. = FALSE)
      vce <- "hc0"
    } else if (vce == "cr1") {
      warning("vce='cr1' requires a cluster variable. Falling back to vce='hc1'.", call. = FALSE)
      vce <- "hc1"
    } else if (vce == "cr2") {
      warning("vce='cr2' requires a cluster variable. Falling back to vce='hc2'.", call. = FALSE)
      vce <- "hc2"
    } else if (vce == "cr3") {
      warning("vce='cr3' requires a cluster variable. Falling back to vce='hc3'.", call. = FALSE)
      vce <- "hc3"
    }
  }

  vce_label <- switch(vce,
                      "hc0" = "HC0", "hc1" = "HC1", "hc2" = "HC2", "hc3" = "HC3",
                      "cr1" = "CR1", "cr2" = "CR2", "cr3" = "CR3", toupper(vce))

  list(vce = vce, vce_label = vce_label)
}

# Fill the effective-sample-size matrix Nh.lev (n.lev x 2) given a
# populated h.lev (n.lev x 2 matrix of L/R bandwidths per group),
# the centered running variable Xc, the treatment indicator T, and
# (for factor covs.hte) the per-observation group label and level
# vector. Used by rdhte() to consolidate the six near-duplicate
# count loops in its bandwidth-selection block.
#
# Continuous covs.hte (or no covs.hte) corresponds to n.lev == 1 OR
# is.null(W.lev): a single row, no group restriction.
#
# Returns the Nh.lev matrix.
.rdhte_fill_Nh <- function(Xc, T, h.lev, W.covs = NULL, W.lev = NULL) {
  n.lev <- nrow(h.lev)
  Nh.lev <- matrix(NA_integer_, n.lev, 2)
  if (n.lev == 1L || is.null(W.lev)) {
    Nh.lev[1, 1] <- sum(abs(Xc) <= h.lev[1, 1] & T == 0)
    Nh.lev[1, 2] <- sum(abs(Xc) <= h.lev[1, 2] & T == 1)
  } else {
    for (l in seq_len(n.lev)) {
      ind <- W.covs == W.lev[l]
      Nh.lev[l, 1] <- sum(abs(Xc[ind]) <= h.lev[l, 1] & T[ind] == 0)
      Nh.lev[l, 2] <- sum(abs(Xc[ind]) <= h.lev[l, 2] & T[ind] == 1)
    }
  }
  Nh.lev
}
