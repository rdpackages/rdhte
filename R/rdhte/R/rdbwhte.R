################################################################################
#' @title MSE-Optimal Bandwidth Selection for RD Heterogeneous Treatment
#' Effects Estimation
#'
#' @description \code{rdbwhte} computes MSE-optimal bandwidths for estimating
#' RD heterogeneous treatment effects based on covariates.
#'
#' Companion commands: \code{\link{rdhte}} for RD HTE estimation and inference.
#'
#' Related Stata and R packages useful for inference in RD designs are
#' described in the website: \url{https://rdpackages.github.io/}.
#'
#' @param y Outcome variable.
#' @param x Running variable.
#' @param c Cutoff value (default = 0).
#' @param covs.hte Covariate(s) for heterogeneous treatment effects (required).
#' @param covs.eff Additional covariates for efficiency (optional).
#' @param p Polynomial order (default = 1).
#' @param kernel Kernel type (default = "tri").
#' @param vce Variance estimator (default = "hc3").
#' @param cluster Optional cluster variable.
#' @param bw.joint Logical, use joint bandwidth selection (default = FALSE).
#'
#' @return A list with selected bandwidths and model information.
#' \item{W_lev}{vector of group level identifiers.}
#' \item{kernel}{kernel type used.}
#' \item{vce}{variance estimator used.}
#' \item{c}{cutoff value.}
#' \item{h}{vector containing the bandwidths used.}
#' \item{p}{order of the polynomial used for estimation of the regression function.}
#' \item{N}{vector with the original number of observations for each group.}
#' \item{Nh}{vector with the effective number of observations for each group.}
#' \item{coef_report}{internal value.}
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
#' @examples
#' set.seed(123)
#' n <- 5000
#' X <- runif(n, -1, 1)
#' W <- rbinom(n, 1, 0.5)
#' Y <- 3 + 2*X + 1.5*X^2 + 0.5*X^3 + sin(2*X) + 3*W*(X>=0) + rnorm(n)
#' rdbwhte.1 = rdbwhte(y=Y, x=X, covs.hte=factor(W))
#' summary(rdbwhte.1)
#' @export
rdbwhte <- function(y, x, c = 0, covs.hte = NULL, covs.eff = NULL, p = 1,
                    kernel = "tri", vce = "hc3", cluster = NULL, bw.joint = FALSE) {

  if (is.null(covs.hte)) {
    warning("rdbwhte requires specifying heterogenity variables via covs.hte.
            MSE-Optimal Bandwidth for Overall RD ATE reported")

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

      rd.bw  <- suppressWarnings(rdbwselect(y = Y, x = Xc, covs = covs, p = p, vce = vce, cluster = cluster, kernel = kernel))
      h      <- rd.bw$bws[1]
      Nh_lev <- sum(abs(Xc)<=h)


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
  if (!is.null(covs)) covs <- model.matrix(~covs-1)

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

  # Bandwidth Selection

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
      h = NULL
      for (l in 1:n_lev) {
        rd.bw      <- suppressWarnings(rdbwselect(y = Y, x = Xc, covs = covs, p = p, vce=vce, cluster = cluster, kernel = kernel, subset=W==W_lev[l]))
        h[l]       <- rd.bw$bws[1]
        ind        <- which(W==W_lev[l])
        h_vec[ind] <- h[l]
        Nh_lev[l]  <- sum( abs(Xc[ind])<=h[l]  )
      }
    }


  if (is.null(W_lev)) {
    W      <- as.matrix(W)
    N_lev  <- N
    Nh_lev <- sum(abs(Xc)<=h)
  }

}


  rdmodel <- "MSE-Optimal Bandwidth Selection for RD Heterogeneous Treatment Effects Estimation"
  if (!is.null(covs.eff)) rdmodel <-paste(rdmodel, ", Covariate-Adjusted", sep="")
  if (!is.null(cluster))  rdmodel <-paste(rdmodel, ", Cluster-Adjusted", sep="")


  # Store model information
  out = list(  W_lev       = W_lev,
               kernel      = rd.bw$kernel,
               vce         = vce.hte,
               c = c,
               h = h,
               p = p,
               N = N_lev,
               Nh = Nh_lev,
               coef_report = coef_report,
               rdmodel = rdmodel)

  out$call <- match.call()
  class(out) <- "rdbwhte"
  return(out)
}


#' Internal function.
#'
#' @param x Class \code{rdhte} objects.
#'
#' @keywords internal
#' @return No return value, called for side effects.
#' @export
print.rdbwhte <- function(x,...){
  cat("MSE-Optimal Bandwidth Selection for RD Heterogeneous Treatment Effects Estimation\n")
  cat(paste("Number of Obs.           ",  format(sum(x$N), width=10, justify="right"), "\n", sep=""))
  cat(paste("Kernel                   ",  format(x$kernel, width=10, justify="right"), "\n", sep=""))
  cat(paste("VCE method               ",  format(x$vce,    width=10, justify="right"), "\n", sep=""))
  cat(paste("Poly. Order (p)          ",  format(x$p,      width=10, justify="right"), "\n", sep=""))
 # cat(paste("Bandwidth (h)            ",  format(sprintf("%10.3f",x$h)),               "\n", sep=""))
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
summary.rdbwhte <- function(object,...) {
  x    <- object
  args <- list(...)

  cat(paste(x$rdmodel,"\n", sep=""))
  cat(paste("","\n", sep=""))

  cat(paste("Number of Obs.             ",  format(sum(x$N), width=10, justify="right"),"\n", sep=""))
  cat(paste("Kernel                     ",  format(x$kernel,   width=10, justify="right"),"\n", sep=""))
  cat(paste("VCE method                 ",  format(x$vce,      width=10, justify="right"),"\n", sep=""))
  cat(paste("Poly. Order (p)            ",  format(x$p,      width=10, justify="right"),      "\n", sep=""))
  cat("\n")

  ### print output
  cat(paste(rep("=", 15 + 15 ), collapse="")); cat("\n")

  if (x$coef_report==FALSE) {
    cat(format("Group"            , width=15, justify="right"))
    cat(format("Bandwidth"        , width=15, justify="right"))
  } else {
    cat(format(""                 , width=15, justify="right"))
    cat(format("Bandwidth"        , width=15, justify="right"))
  }

  cat("\n")

  cat(paste(rep("=", 15 + 15 ), collapse="")); cat("\n")

  if (x$coef_report==FALSE) {
  for (i in 1:length(x$W_lev)) {
    cat(format(x$W_lev[i],              width=15, justify="right"))
    cat(format(sprintf("%3.3f",x$h[i]), width=15, justify="right"))
    cat("\n")
  }
  } else {
    cat(format("Overall"              , width=15, justify="right"))
    cat(format(sprintf("%3.3f",x$h[1]), width=15, justify="right"))
    cat("\n")
  }

  cat(paste(rep("=", 15 + 15 ), collapse="")); cat("\n")

}



