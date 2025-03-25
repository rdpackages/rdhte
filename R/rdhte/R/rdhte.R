################################################################################
#' @title RD Heterogeneous Treatment Effects Estimation and Inference
#'
#' @description \code{rdhte} provides estimation and inference for heterogeneous
#' treatment effects in RDD using local polynomial regression
#' allowing for interactions with pretreatment covariates.
#' Inference is implemented using robust bias-correction methods.
#'
#' Companion commands: \code{\link{rdbwhte}} for data-driven bandwidth selection.
#'
#' Related Stata and R packages useful for inference in RD designs are described
#' in the website: \url{https://rdpackages.github.io/}.
#'
#' @param y Outcome variable.
#' @param x Running variable.
#' @param c Cutoff value (default = 0).
#' @param covs.hte Covariate(s) for heterogeneous treatment effects (required).
#' @param covs.eff Additional covariates for efficiency (optional).
#' @param p Polynomial order (default = 1).
#' @param kernel Kernel type (default = "tri").
#' @param h Choice of bandwidth (optional).
#' @param vce Variance estimator (default = "hc3").
#' @param cluster Optional cluster variable.
#' @param level Confidence level (default = 95).
#' @param bw.joint Logical, use joint bandwidth selection (default = FALSE).
#'
#' @return A list with selected RD HTE effects and model information.
#' \item{Estimate}{vector of conventional local-polynomial RD estimates.}
#' \item{Estimate_bc}{vector of bias-corrected local-polynomial RD estimates.}
#' \item{se_rb}{vector containing robust bias corrected standard errors of the local-polynomial RD estimates.}
#' \item{ci_rb}{matrix containing robust bias corrected confidence intervals.}
#' \item{t_rb}{vector containing the t-statistics associated with robust local-polynomial RD estimates.}
#' \item{pv_rb}{vector containing the p-values associated with robust local-polynomial RD estimates.}
#' \item{coefs}{vector containing the coefficients for the jointly estimated p-th order local polynomial model.}
#' \item{vcov}{estimated variance-covariance matrix.}
#' \item{W_lev}{vector of group level identifiers.}
#' \item{kernel}{kernel type used.}
#' \item{vce}{variance estimator used.}
#' \item{c}{cutoff value.}
#' \item{h}{vector containing the bandwidths used.}
#' \item{p}{order of the polynomial used for estimation of the regression function.}
#' \item{N}{vector with the original number of observations for each group.}
#' \item{Nh}{vector with the effective number of observations for each group.}
#' \item{coef_report}{internal value.}
#' \item{level}{confidence level used.}
#' \item{rdmodel}{rd model.}
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
#' Calonico, Cattaneo, Farrell, Palomba and Titiunik (2025): Heterogenous Treatment Effects in Regression Discontinuity Designs. \emph{Working paper}
#'
#' @seealso \code{\link{rdbwhte}}
#'
#'
#' @examples
#' set.seed(123)
#' n <- 5000
#' X <- runif(n, -1, 1)
#' W <- rbinom(n, 1, 0.5)
#' Y <- 3 + 2*X + 1.5*X^2 + 0.5*X^3 + sin(2*X) + 3*W*(X>=0) + rnorm(n)
#' rdhte.1 = rdhte(y=Y, x=X, covs.hte=factor(W))
#' summary(rdhte.1)
#' @export
rdhte <- function(y, x, c = 0, covs.hte = NULL, covs.eff = NULL, p = 1, kernel = "tri", h = NULL,
                  vce = "hc3", cluster = NULL, level = 95, bw.joint = FALSE) {

  if (is.null(covs.hte)) {
    warning("rdhte requires specifying heterogenity variables via covs.hte. RD ATE reported")

    # Create initial data frame and remove missing values
    dd <- data.frame(Y = y, X = x)
    if (!is.null(cluster))  dd$C    <- cluster
    if (!is.null(covs.eff)) dd$covs <- covs.eff
    dd <- na.omit(dd)
    N  <- nrow(dd)

    # Extract relevant variables
    Xc  <- dd$X - c
    Y   <- dd$Y
    covs    <- if (!is.null(covs.eff)) dd$covs else NULL
    cluster <- if (!is.null(cluster))  dd$C    else NULL
    T <- as.integer(Xc >= 0)

    if (!is.null(covs)) covs <- model.matrix(~covs-1)


    # Polynomial transformation of running variable (Xc)
    Xp0 <- poly(Xc, raw = TRUE, degree = p)
    Xp1 <- poly(Xc, raw = TRUE, degree = p + 1)

    # Standardize kernel and variance estimator case
    kernel  <- tolower(kernel)
    vce     <- tolower(vce)
    vce.hte <- toupper(vce)

    # Factor or continuous covariate handling for W
    coef_report <- FALSE
    W_lev <- "Overall"
    n_lev <- 1
    N_lev <- N


    # Bandwidth is not provided

    if (is.null(h)) {
        rd.bw  <- suppressWarnings(rdbwselect(y = Y, x = Xc, covs = covs, p = p, vce = vce, cluster = cluster, kernel = kernel))
        h      <- rd.bw$bws[1]
    }

    Nh_lev <- sum(abs(Xc)<=h)

    # Weighted Regression
    r.bw <- abs(Xc) <= h
    kweight = NULL
    Kernel = "Uniform"
    if (kernel=="tri") {
      kweight = 1-abs(Xc)/h
      Kernel = "Triangular"
    }
    if (kernel=="epa") {
      kweight = 1-(Xc/h)^2
      Kernel = "Epanechnikov"
    }
    if (is.null(covs.eff)) {
      reg_formula_est <- Y ~ T * Xp0
      reg_formula_inf <- Y ~ T * Xp1
    } else {
      reg_formula_est <- Y ~ T * Xp0 * W + covs
      reg_formula_inf <- Y ~ T * Xp1 * W + covs
    }

    # Run a single regression model with weights
    rd_est <- lm(reg_formula_est, weights = kweight, subset = r.bw)
    rd_inf <- lm(reg_formula_inf, weights = kweight, subset = r.bw)

    # Compute robust variance-covariance matrix
    rd_vcov <- vcovCL(rd_inf, cluster = cluster[r.bw], type = vce.hte)

    # Extract treatment effect estimate
    tau.hat    = rd_est$coefficients["T"]
    tau.hat.bc = rd_inf$coefficients["T"]
    tau.hat.se = sqrt(rd_vcov["T","T"])


    #tau.hat.pv = 2*pnorm(-abs(tau.hat.bc / tau.hat.se))
    #results_df = NULL
  } else {




  # Create initial data frame and remove missing values
  dd <- data.frame(Y = y, X = x, covs.hte = covs.hte)
  if (!is.null(cluster))  dd$C    <- cluster
  if (!is.null(covs.eff)) dd$covs <- covs.eff
  dd <- na.omit(dd)
  N  <- nrow(dd)

  # Extract relevant variables
  Xc  <- dd$X - c
  Y   <- dd$Y
  idx <- grep("^covs.hte", names(dd))
  W   <- dd[, idx]
  covs    <- if (!is.null(covs.eff)) dd$covs else NULL
  cluster <- if (!is.null(cluster))  dd$C    else NULL
  T <- as.integer(Xc >= 0)

  if (!is.null(covs))  covs <- model.matrix(~covs-1)

  # Polynomial transformation of running variable (Xc)
  Xp0 <- poly(Xc, raw = TRUE, degree = p)
  Xp1 <- poly(Xc, raw = TRUE, degree = p + 1)

  # Standardize kernel and variance estimator case
  kernel  <- tolower(kernel)
  vce     <- tolower(vce)
  vce.hte <- toupper(vce)

  # Factor or continuous covariate handling for W
  coef_report <- TRUE
  W_lev <- NULL
  n_lev <- 1
  N_lev <- 0

  if (is.factor(W)) {
    W_lev <- levels(W)
    n_lev <- nlevels(W)
    N_lev <- table(W)
    coef_report <- FALSE
  } else {
    if (length(idx) == 1) {
      cond <- mean((W == 0) | (W == 1))
        if (cond == 1 ) {
          W <- factor(W)
          W_lev <- levels(W)
          n_lev <- nlevels(W)
          N_lev <- table(W)
          coef_report <- FALSE
        }
    } else {
    cond1 <- mean((rowSums(W) == 1) | (ncol(W) == 1))
    cond2 <- mean((W == 0) | (W == 1))
    if (cond1 == 1 & cond2 == 1) {
      col_names <- names(W)
      W     <- factor(apply(W[, col_names], 1, function(x) which(x == 1)))
      W_lev <- levels(W)
      n_lev <- nlevels(W)
      N_lev <- table(W)
      coef_report <- FALSE
    }
  }
}

  # Initialize bandwidth storage
  h_vec <- rep(0, N)
  Nh_lev <- numeric(n_lev)

  # Bandwidth is provided (1 or many)
  if (!is.null(h)) {
    if (length(h) == 1) {
      h_vec <- rep(h, N)
      h     <- rep(h, n_lev)
    }
    if (n_lev>1){
      hdif = abs(n_lev - length(h))
      if (hdif > 0) stop("check the number of bandwidths provided")
      for (l in 1:n_lev) {
        ind        <- W == W_lev[l]
        h_vec[ind] <- h[l]
        Nh_lev[l]  <- sum(abs(Xc[ind]) <= h[l])
      }
    }
  }


  # Bandwidth is not provided

  if (is.null(h)) {

    if (bw.joint==TRUE | is.null(W_lev)) {
      # Joint bandwidth estimation
      rd.bw  <- suppressWarnings(rdbwselect(y = Y, x = Xc, covs = covs, p = p, vce = vce, cluster = cluster, kernel = kernel))
      h      <- rep(rd.bw$bws[1], n_lev)
      h_vec  <- rep(rd.bw$bws[1], N)
      if (n_lev>1) {
        for (l in 1:n_lev) {
          ind       <- which(W==W_lev[l])
          Nh_lev[l] <- sum(abs(Xc[ind])<=h[l])
        }
      }
    } else {
        for (l in 1:n_lev) {
          rd.bw      <- suppressWarnings(rdbwselect(y = Y, x = Xc, covs = covs, p = p, vce=vce, cluster = cluster, kernel = kernel, subset=W==W_lev[l]))
          h[l]       <- rd.bw$bws[1]
          ind        <- which(W==W_lev[l])
          h_vec[ind] <- h[l]
          Nh_lev[l]  <- sum( abs(Xc[ind])<=h[l]  )
        }
      }
    }

  if (is.null(W_lev)) {
    W      <- as.matrix(W)
    N_lev  <- N
    Nh_lev <- sum(abs(Xc)<=h)
  }


  # Weighted Regression
    r.bw <- abs(Xc) <= h_vec
    kweight = NULL
    Kernel = "Uniform"
    if (kernel=="tri") {
      kweight = 1-abs(Xc)/h_vec
      Kernel = "Triangular"
      }
    if (kernel=="epa") {
      kweight = 1-(Xc/h_vec)^2
      Kernel = "Epanechnikov"
    }
    if (is.null(covs.eff)) {
      reg_formula_est <- Y ~ T * Xp0 * W
      reg_formula_inf <- Y ~ T * Xp1 * W
    } else {
      reg_formula_est <- Y ~ T * Xp0 * W + covs * W
      reg_formula_inf <- Y ~ T * Xp1 * W + covs * W
    }

    # Run a single regression model with weights
    rd_est <- lm(reg_formula_est, weights = kweight, subset = r.bw)
    rd_inf <- lm(reg_formula_inf, weights = kweight, subset = r.bw)

    # Compute robust variance-covariance matrix
    rd_vcov <- vcovCL(rd_inf, cluster = cluster[r.bw], type = vce.hte)

    # Extract treatment effect estimate
    tau.hat    = rd_est$coefficients["T"]
    tau.hat.bc = rd_inf$coefficients["T"]
    tau.hat.se = sqrt(rd_vcov["T","T"])

    if (!is.factor(W)) {
      ncoeff = names(rd_est$coefficients[!is.na(rd_est$coefficients)])
      W_lev = c("T", ncoeff[grepl("T:W", ncoeff)])
      n_lev <- length(W_lev)
      for (j in 2:n_lev) {
        lev = W_lev[j]
        tau.hat[j]    = rd_est$coefficients[lev]
        tau.hat.bc[j] = rd_inf$coefficients[lev]
        tau.hat.se[j] = sqrt(rd_vcov[lev,lev])
      }
    } else {
    for (j in 2:n_lev) {
        lev = paste("T:W", W_lev[j], sep="")
        tau.hat[j]    = rd_est$coefficients["T"] + rd_est$coefficients[lev]
        tau.hat.bc[j] = rd_inf$coefficients["T"] + rd_inf$coefficients[lev]
        tau.hat.se[j] = sqrt(rd_vcov["T","T"] + rd_vcov[lev,lev] + 2*rd_vcov["T",lev])
    }
  }

    # Construct results table
    #results_df <- data.frame(
    #  Group       = W_lev,
    #  Estimate    = tau.hat,
    #  Estimate.bc = tau.hat.bc,
    #  SE          = tau.hat.se,
    #  p.value     = tau.hat.pv
    #)

  }

  qz <- -qnorm(abs((1-(level/100))/2))
  t_rb = tau.hat.bc/tau.hat.se
  pv_rb = 2*pnorm(-abs(t_rb))


    rdmodel <- "RD Heterogeneous Treatment Effects Estimation"
    if (!is.null(covs.eff)) rdmodel <-paste(rdmodel, ", Covariate-Adjusted", sep="")
    if (!is.null(cluster))  rdmodel <-paste(rdmodel, ", Cluster-Adjusted", sep="")


 # Store model information
  out = list(  Estimate    = tau.hat,
               Estimate_bc = tau.hat.bc,
               se_rb       = tau.hat.se,
               ci_rb       = cbind(tau.hat.bc - qz*tau.hat.se, tau.hat.bc + qz*tau.hat.se),
               t_rb        = t_rb,
               pv_rb       = pv_rb,
               coefs       = rd_est$coefficients,
               vcov        = rd_vcov,
               W_lev       = W_lev,
               kernel      = Kernel,
               vce         = vce.hte,
               c = c,
               h = h,
               p = p,
               N = N_lev,
               Nh = Nh_lev,
               coef_report = coef_report,
               level = level,
               rdmodel = rdmodel)

    out$call <- match.call()
    class(out) <- "rdhte"
    return(out)
}


#' Internal function.
#'
#' @param x Class \code{rdhte} objects.
#'
#' @keywords internal
#' @return No return value, called for side effects.
#' @export
print.rdhte <- function(x,...){
  cat("RD Heterogeneous Treatment Effects Estimation\n")
  cat(paste("Number of Obs.           ",  format(sum(x$N), width=10, justify="right"), "\n", sep=""))
  cat(paste("Kernel                   ",  format(x$kernel, width=10, justify="right"), "\n", sep=""))
  cat(paste("VCE method               ",  format(x$vce,    width=10, justify="right"), "\n", sep=""))
  cat(paste("Poly. Order (p)          ",  format(x$p,      width=10, justify="right"), "\n", sep=""))
  cat(paste("Bandwidth (h)            ",  format(sprintf("%10.3f",x$h)),               "\n", sep=""))
  cat("\n")
}



################################################################################
#' Internal function.
#'
#' @param object Class \code{rdhte} objects.
#'
#' @keywords internal
#' @return No return value, called for side effects.
#' @export
summary.rdhte <- function(object,...) {
  x    <- object
  args <- list(...)

  cat(paste(x$rdmodel,"\n", sep=""))
  cat(paste("","\n", sep=""))

  cat(paste("Number of Obs.             ",  format(sum(x$N), width=10, justify="right"),"\n", sep=""))
  if (x$coef_report==TRUE) {
  cat(paste("Bandwidth (h)              ",  format(sprintf("%10.3f",x$h)),     "\n", sep=""))
  cat(paste("Eff. Number of Obs.        ",  format(x$Nh, width=10, justify="right"),       "\n", sep=""))
  }
  cat(paste("Kernel                     ",  format(x$kernel,   width=10, justify="right"),"\n", sep=""))
  cat(paste("VCE method                 ",  format(x$vce,      width=10, justify="right"),"\n", sep=""))
  cat(paste("Poly. Order (p)            ",  format(x$p,      width=10, justify="right"),      "\n", sep=""))


  cat("\n")

  ### print output
  cat(paste(rep("=", 15 + 15 + 25 + 15 + 15 + 10), collapse="")); cat("\n")

  if (x$coef_report==FALSE) {
    cat(format("Group"            , width=15, justify="right"))
    cat(format("RD Effect"        , width=15, justify="right"))
  } else {
    cat(format("Coeff"            , width=15, justify="right"))
    cat(format("Estimate"         , width=15, justify="right"))
  }

  cat(format(paste("[", x$level, "% ", "Robust C.I.]", sep=""), width=25, justify="centre"))
  cat(format("Pr(>|t|)"            , width=15, justify="right"))

  if (x$coef_report==FALSE) {
    cat(format("Group Size"       , width=15, justify="right"))
    cat(format("h"                , width=10, justify="right"))
  }
  cat("\n")

  cat(paste(rep("=", 15 + 15 + 25 + 15 + 15 + 10), collapse="")); cat("\n")

  for (i in 1:length(x$W_lev)) {
    cat(format(x$W_lev[i],                                                width=15, justify="right"))
    cat(format(sprintf("%3.3f", x$Estimate[i]),                           width=15, justify="right"))
    cat(format(paste("[", sprintf("%3.3f", x$ci_rb[i,1]), " , ", sep=""), width=15, justify="right"))
    cat(format(paste(     sprintf("%3.3f", x$ci_rb[i,2]), "]",   sep=""), width=10, justify="left"))
    cat(format(sprintf("%3.3f", x$pv_rb[i]),                              width=15, justify="right"))

    if (x$coef_report==FALSE) {
      cat(format(x$Nh[i],                 width=15, justify="right"))
      cat(format(sprintf("%3.3f",x$h[i]), width=10, justify="right"))
    }
    cat("\n")
  }

  cat(paste(rep("=", 15 + 15 + 25 + 15 + 15 + 10), collapse="")); cat("\n")
}



