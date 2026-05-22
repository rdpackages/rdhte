 ################################################################################
#' @title Data-Driven Optimal Bandwidth Selection for RD Heterogeneous Treatment
#' Effects Estimation
#'
#' @description \code{rdbwhte} computes MSE- and CER-optimal bandwidths for
#' estimating RD heterogeneous treatment effects based on covariates
#' (Calonico, Cattaneo, Farrell, Palomba and Titiunik, 2025a).
#'
#' Companion commands: \code{\link{rdhte}} for RD HTE estimation and inference,
#' and \code{\link{rdhte_lincom}} for testing linear restrictions of parameters.
#'
#' A detailed introduction to the software is given in Calonico, Cattaneo,
#' Farrell, Palomba and Titiunik (2025b). Related software packages for
#' analysis and interpretation of RD designs and related methods are available
#' in: \url{https://rdpackages.github.io/}.
#'
#' For background methodology, see Calonico, Cattaneo, Farrell, and Titiunik
#' (2019), Calonico, Cattaneo and Farrell (2020), and Cattaneo and Titiunik
#' (2022).
#'
#' @param y Outcome variable.
#' @param x Running variable.
#' @param c RD cutoff in \code{x}; default is \code{c = 0}.
#' @param covs.hte covariates for heterogeneous treatment effects. Factor variables can be used to distinguish between continuous and categorical variables, select reference categories, specify interactions between variables, and include polynomials of continuous variables.
#' @param covs.eff additional covariates to be used for efficiency improvements.
#' @param p order of the local polynomial used to construct the point estimator (default = 1).
#' @param q order of the local polynomial used to construct the bias correction.
#'   If \code{NULL} (default), \code{q} is set to \code{p + 1}.
#' @param kernel kernel function used to construct the RD estimators. Options are \code{triangular} (default option), \code{epanechnikov} and \code{uniform}.
#' @param weights variable used for optional weighting of the bandwidth-selection procedure. The unit-specific weights multiply the kernel function.
#' @param vce character string specifying the variance-covariance matrix
#'   estimator type. Without \code{cluster}: \code{"hc0"}, \code{"hc1"},
#'   \code{"hc2"}, \code{"hc3"} (default \code{"hc3"}). With
#'   \code{cluster}: \code{"cr1"} (default), \code{"cr2"}, \code{"cr3"}.
#'   Legacy aliases: \code{"hc0"}/\code{"hc1"} + \code{cluster} are
#'   remapped to \code{"cr1"} with a warning; \code{"hc2"} -> \code{"cr2"}
#'   and \code{"hc3"} -> \code{"cr3"} similarly. \code{"cr1"},
#'   \code{"cr2"}, \code{"cr3"} without \code{cluster} fall back to
#'   \code{"hc1"}, \code{"hc2"}, \code{"hc3"} with a warning.
#' @param cluster variable indicating the clustering of observations.
#' @param subset optional vector specifying a subset of observations to be used.
#' @param data optional data frame. When supplied, \code{y}, \code{x}, \code{covs.hte}, \code{covs.eff}, \code{weights}, \code{cluster}, and \code{subset} may be given as bare variable names referring to columns of \code{data}.
#' @param bwselect bandwidth selection procedure to be used.
#' Options are:
#' \code{mserd} one common MSE-optimal bandwidth selector for the RD treatment effect estimator.
#' \code{msetwo} two different MSE-optimal bandwidth selectors (below and above the cutoff) for the RD treatment effect estimator.
#' \code{msesum} one common MSE-optimal bandwidth selector for the sum of regression estimates (as opposed to difference thereof).
#' \code{msecomb1} for min(\code{mserd},\code{msesum}).
#' \code{msecomb2} for median(\code{msetwo},\code{mserd},\code{msesum}), for each side of the cutoff separately.
#' \code{cerrd} one common CER-optimal bandwidth selector for the RD treatment effect estimator.
#' \code{certwo} two different CER-optimal bandwidth selectors (below and above the cutoff) for the RD treatment effect estimator.
#' \code{cersum} one common CER-optimal bandwidth selector for the sum of regression estimates (as opposed to difference thereof).
#' \code{cercomb1} for min(\code{cerrd},\code{cersum}).
#' \code{cercomb2} for median(\code{certwo},\code{cerrd},\code{cersum}), for each side of the cutoff separately.
#' Note: MSE = Mean Square Error; CER = Coverage Error Rate. Default is \code{bwselect=mserd}.
#' @param bw.joint logical. If \code{TRUE}, forces all bandwidths to be the same across groups (default is \code{bw.joint = FALSE}). When \code{covs.hte} is continuous (rather than a factor or 0/1 indicator), a single joint bandwidth is always used regardless of this argument.
#'
#' @return A list with the following named elements:
#' \item{W.lev}{Group-level identifiers, or \code{NULL} for continuous \code{covs.hte}.}
#' \item{W.names}{Display labels for the rows of \code{h} (prefixed with the
#'   \code{covs.hte} expression for categorical \code{covs.hte}).}
#' \item{covs.hte_chr}{Character representation of the \code{covs.hte} argument.}
#' \item{kernel}{Kernel type used.}
#' \item{vce}{Variance estimator display label.}
#' \item{vce_select}{Canonical lowercase variance-estimator name.}
#' \item{c}{Cutoff value.}
#' \item{h}{An \code{n.lev x 2} matrix of left/right bandwidths, one row per group.}
#' \item{p}{Order of the polynomial used for estimation.}
#' \item{q}{Order of the polynomial used for bias correction.}
#' \item{bwselect}{Bandwidth selection procedure used.}
#' \item{N}{Length-2 vector \code{c(N_left, N_right)} of pre-bandwidth sample sizes.}
#' \item{Nh}{An \code{n.lev x 2} matrix of effective sample sizes (per group, per side).}
#' \item{covs.cont}{Logical; \code{TRUE} for continuous \code{covs.hte} (or no \code{covs.hte}), \code{FALSE} for factor.}
#' \item{rdmodel}{Human-readable model description string.}
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
#' @seealso \code{\link{rdhte}}, \code{\link{rdhte_lincom}}
#'
#' @examples
#' set.seed(123)
#' n <- 5000
#' X <- runif(n, -1, 1)
#' W <- rbinom(n, 1, 0.5)
#' Y <- 3 + 2*X + 1.5*X^2 + 0.5*X^3 + sin(2*X) + 3*W*(X>=0) + rnorm(n)
#' rdbwhte.1 = rdbwhte(y=Y, x=X, covs.hte=factor(W))
#' summary(rdbwhte.1)
#'
#' \dontrun{
#' data(rdhte_dataset)
#' with(rdhte_dataset, {
#'   summary(rdbwhte(y = y, x = x, covs.hte = factor(w_ideology),
#'                   cluster = cluster_var))
#'   summary(rdbwhte(y = y, x = x, covs.hte = factor(w_ideology),
#'                   cluster = cluster_var, bw.joint = TRUE))
#'   summary(rdbwhte(y = y, x = x, covs.hte = w_strength,
#'                   cluster = cluster_var))
#' })
#' }
#' @export
rdbwhte <- function(y, x, c = 0, covs.hte = NULL, covs.eff = NULL, p = 1, q = NULL,
                    kernel = "tri", weights = NULL, vce = "hc3", cluster = NULL,
                    bwselect = "mserd", bw.joint = FALSE, subset = NULL,
                    data = NULL) {

  ## Capture the literal covs.hte expression for downstream slot reporting.
  ## Mirrors rdhte() so R rdbwhte and Python rdbwhte agree on this slot.
  covs.hte_chr <- if (missing(covs.hte) || is.null(substitute(covs.hte))) NULL
                  else paste0(deparse(substitute(covs.hte)), collapse = " ")
  if (!is.null(covs.hte_chr) && nchar(covs.hte_chr) > 100) covs.hte_chr <- "covs.hte"

  if (!is.null(data)) {
    .vars <- .rdhte_resolve_data(match.call(), data, parent.frame(),
              c("y", "x", "covs.hte", "covs.eff", "weights", "cluster", "subset"))
    y        <- .vars$y
    x        <- .vars$x
    covs.hte <- .vars$covs.hte
    covs.eff <- .vars$covs.eff
    weights  <- .vars$weights
    cluster  <- .vars$cluster
    subset   <- .vars$subset
  }

  ## Input validation (mirror rdhte()).
  if (!is.numeric(c) || length(c) != 1L || !is.finite(c))
    stop(sprintf("Cutoff 'c' must be a single finite numeric value (received: %s).", toString(c)), call. = FALSE)
  if (!is.numeric(p) || length(p) != 1L || p < 0 || p != round(p))
    stop(sprintf("Polynomial order 'p' must be a single non-negative integer (received: %s).", toString(p)), call. = FALSE)
  if (is.null(q)) q <- p + 1L
  if (!is.numeric(q) || length(q) != 1L || q <= p || q != round(q))
    stop(sprintf("Polynomial order 'q' must be a single integer with q > p (received q=%s, p=%s).", toString(q), toString(p)), call. = FALSE)
  .n_orig <- length(x)
  if (length(y) != .n_orig)
    stop(sprintf("'y' and 'x' must have equal length (got y=%d, x=%d).", length(y), .n_orig), call. = FALSE)
  if (!is.null(weights) && length(weights) != .n_orig)
    stop(sprintf("'weights' must have length equal to length(x) (got %d, expected %d).", length(weights), .n_orig), call. = FALSE)
  if (!is.null(cluster) && length(cluster) != .n_orig)
    stop(sprintf("'cluster' must have length equal to length(x) (got %d, expected %d).", length(cluster), .n_orig), call. = FALSE)
  if (!is.null(covs.eff)) {
    .nc_eff <- if (is.matrix(covs.eff) || is.data.frame(covs.eff)) nrow(covs.eff) else length(covs.eff)
    if (.nc_eff != .n_orig)
      stop(sprintf("'covs.eff' must have nrow equal to length(x) (got %d, expected %d).", .nc_eff, .n_orig), call. = FALSE)
  }
  if (!is.null(covs.hte) && !inherits(covs.hte, "formula") && !is.character(covs.hte)) {
    .nc_hte <- if (is.matrix(covs.hte) || is.data.frame(covs.hte)) nrow(covs.hte) else length(covs.hte)
    if (.nc_hte != .n_orig)
      stop(sprintf("'covs.hte' must have nrow equal to length(x) (got %d, expected %d).", .nc_hte, .n_orig), call. = FALSE)
  }
  if (!is.null(subset)) {
    if (is.logical(subset)) {
      if (length(subset) != .n_orig)
        stop(sprintf("Logical 'subset' must have length equal to length(x) (got %d, expected %d).", length(subset), .n_orig), call. = FALSE)
    } else if (is.numeric(subset)) {
      if (any(!is.finite(subset)) || any(subset < 1) || any(subset > .n_orig) || any(subset != round(subset)))
        stop(sprintf("Numeric 'subset' must contain integer indices in 1..%d.", .n_orig), call. = FALSE)
    } else {
      stop("'subset' must be logical or integer.", call. = FALSE)
    }
  }

  ## vce normalization + label mapping. Shared with rdhte() via
  ## .rdhte_normalize_vce() in helpers.R. See rdhte.R for documentation.
  .user_vce <- !missing(vce)
  .vce_norm <- .rdhte_normalize_vce(vce, cluster, .user_vce)
  vce       <- .vce_norm$vce
  vce_label <- .vce_norm$vce_label

  if (is.null(covs.hte)) {
    # Create initial data frame and remove missing values
    dd <- data.frame(Y = y, X = x)
    if (!is.null(cluster))  dd$C    <- cluster
    if (!is.null(covs.eff)) dd$covs <- covs.eff
    if (!is.null(weights))  dd$weights <- weights
    if (!is.null(subset)) {
      dd <- dd[subset, , drop = FALSE]
    }
    dd <- na.omit(dd)
    N  <- nrow(dd)

    # Extract relevant variables
    Xc  <- dd$X - c
    Y   <- dd$Y
    covs    <- if (!is.null(covs.eff)) dd$covs    else NULL
    cluster <- if (!is.null(cluster))  dd$C       else NULL
    weights <- if (!is.null(weights))  dd$weights else NULL
    T <- as.integer(Xc >= 0)

    if (!is.null(covs)) covs <- model.matrix(~covs-1)

    # Standardize kernel and variance estimator case
    kernel  <- tolower(kernel)
    vce     <- tolower(vce)

    # Factor or continuous covariate handling for W
    covs.cont <- TRUE
    W.lev <- "Overall"
    n.lev <- 1

    rd.bw  <- suppressWarnings(rdbwselect(y = Y, x = Xc, covs = covs, p = p, q = q,
                                          bwselect = bwselect, vce = vce, cluster = cluster,
                                          kernel = kernel, weights = weights))
    h.lev      <- matrix(c(rd.bw$bws[1], rd.bw$bws[2]),1,2)
    N.lev <- c(sum(T==0), sum(T==1))
    Nh.lev.l <- sum( (abs(Xc) <= h.lev[1,1]) & T==0)
    Nh.lev.r <- sum( (abs(Xc) <= h.lev[1,2]) & T==1)
    Nh.lev <- matrix(c(Nh.lev.l, Nh.lev.r),1,2)

  } else {

  ### Check if covs.hte is a formula expression (passed as a string):
  is_char <- is.character(covs.hte)
  cls     <- if (is_char) "formula" else "operation"
  if (cls == "formula") {
    txt <- covs.hte[1]
    if (startsWith(trimws(txt), "~")) {
      fml.txt <- paste0(txt, "-1")
    } else {
      fml.txt <- paste0("~", txt, "-1")
    }
    f <- as.formula(fml.txt, env = parent.frame())
    if (!is.null(data)) {
      mf <- stats::model.frame(f, data = data, na.action = stats::na.pass)
      covs.hte.mm <- stats::model.matrix(f, data = mf)
    } else {
      covs.hte.mm <- stats::model.matrix(f)
    }
    ## See rdhte.R for the rationale on this prefix dance.
    covs.hte.df <- as.data.frame(covs.hte.mm)
    names(covs.hte.df) <- paste0("covs.hte.", names(covs.hte.df))
    dd <- cbind(data.frame(Y = y, X = x), covs.hte.df)
  } else {
    dd <- data.frame(Y = y, X = x, covs.hte = covs.hte)
  }

  if (!is.null(cluster))  dd$C    <- cluster
  if (!is.null(covs.eff)) dd$covs <- covs.eff
  if (!is.null(weights))  dd$weights <- weights

  if (!is.null(subset)) {
    dd <- dd[subset, , drop = FALSE]
  }
  dd <- na.omit(dd)
  N  <- nrow(dd)

  # Extract relevant variables
  Xc  <- dd$X - c
  Y   <- dd$Y
  idx <- grep("^covs.hte", names(dd))
  W   <- dd[, idx]
  covs    <- if (!is.null(covs.eff)) dd$covs    else NULL
  cluster <- if (!is.null(cluster))  dd$C       else NULL
  weights <- if (!is.null(weights))  dd$weights else NULL
  T <- as.integer(Xc >= 0)
  if (!is.null(covs)) covs <- model.matrix(~covs-1)

  # Standardize kernel and variance estimator case
  kernel  <- tolower(kernel)
  vce     <- tolower(vce)

  # Factor or continuous covariate handling for W
  covs.cont <- TRUE
  W.lev <- NULL
  n.lev <- 1
  N.lev <- c(sum(T == 0), sum(T == 1))

  if (is.factor(W)) {
    W.lev <- levels(W)
    n.lev <- nlevels(W)
    covs.cont <- FALSE
  } else {
    if (length(idx) == 1) {
      cond <- mean((W == 0) | (W == 1))
      if (cond == 1 ) {
        W <- factor(W)
        W.lev <- levels(W)
        n.lev <- nlevels(W)
        covs.cont <- FALSE
      }
    } else {
      cond1 <- mean((rowSums(W) == 1) | (ncol(W) == 1))
      cond2 <- mean((W == 0) | (W == 1))
      if (cond1 == 1 & cond2 == 1) {
        col_names <- names(W)
        W     <- factor(apply(W[, col_names], 1, function(x) which(x == 1)))
        W.lev <- levels(W)
        n.lev <- nlevels(W)
        covs.cont <- FALSE
      }
    }
  }

  # Initialize bandwidth storage
  #h.vec <- rep(0, N)
  #Nh.lev <- numeric(n.lev)

  h.vec  <- rep(0, N)
  h.lev  <- matrix(0, n.lev, 2)
  Nh.lev <- matrix(NA, n.lev, 2)

  # Bandwidth Selection

    if (isTRUE(bw.joint) | is.null(W.lev)) {
      # Joint bandwidth estimation
      rd.bw  <- suppressWarnings(rdbwselect(y = Y, x = Xc, covs = covs, p = p, q = q,
                                            bwselect = bwselect, vce = vce, cluster = cluster,
                                            kernel = kernel, weights = weights))

      h.lev[,1]  <- rep(rd.bw$bws[1], n.lev)
      h.lev[,2]  <- rep(rd.bw$bws[2], n.lev)

      ind.l <- Xc<0
      ind.r <- Xc>=0
      h.vec[ind.l] <- rd.bw$bws[1]
      h.vec[ind.r] <- rd.bw$bws[2]

      if (n.lev>1) {
        for (l in 1:n.lev) {
          ind       <- which(W==W.lev[l])
          Nh.lev[l,1]  <- sum( (abs(Xc[ind]) <= h.lev[l,1]) & (T[ind]==0))
          Nh.lev[l,2]  <- sum( (abs(Xc[ind]) <= h.lev[l,2]) & (T[ind]==1))
        }
      }
    } else {
      for (l in 1:n.lev) {
        rd.bw      <- suppressWarnings(rdbwselect(y = Y, x = Xc, covs = covs, p = p,  q = q,
                                                  bwselect = bwselect, vce=vce, cluster = cluster,
                                                  kernel = kernel, weights = weights,
                                                  subset = W==W.lev[l]))

        h.lev[l,1]       <- rd.bw$bws[1]
        h.lev[l,2]       <- rd.bw$bws[2]

        ind.l     <- W == W.lev[l] & T==0
        ind.r     <- W == W.lev[l] & T==1
        h.vec[ind.l] <- rd.bw$bws[1]
        h.vec[ind.r] <- rd.bw$bws[2]

        Nh.lev[l,1]  <- sum( (abs(Xc[ind.l]) <= rd.bw$bws[1] ))
        Nh.lev[l,2]  <- sum( (abs(Xc[ind.r]) <= rd.bw$bws[2] ))

      }
    }


  if (is.null(W.lev)) {
    W      <- as.matrix(W)

    Nh.lev.l <- sum( (abs(Xc) <= h.lev[1,1]) & T==0)
    Nh.lev.r <- sum( (abs(Xc) <= h.lev[1,2]) & T==1)
    ## Keep a 1 x 2 matrix shape (consistent with the factor branch and
    ## rdhte's $Nh slot). Avoids downstream code branching on shape.
    Nh.lev <- matrix(c(Nh.lev.l, Nh.lev.r), 1, 2)
  }

}


  rdmodel <- "MSE-Optimal Bandwidth Selection for RD Heterogeneous Treatment Effects Estimation"
  if (!is.null(covs.eff)) rdmodel <-paste(rdmodel, ", Covariate-Adjusted", sep="")
  if (!is.null(cluster))  rdmodel <-paste(rdmodel, ", Cluster-Adjusted", sep="")


  ## W.names: prefix each level with the original covs.hte expression so it
  ## reads sensibly in tidy()/glance() (parity with Python rdbwhte).
  W.names <- if (!is.null(W.lev) && !identical(W.lev, "Overall"))
               paste(covs.hte_chr, W.lev, sep = "")
             else W.lev

  # Store model information
  out = list(  W.lev       = W.lev,
               W.names     = W.names,
               covs.hte_chr = covs.hte_chr,
               kernel      = rd.bw$kernel,
               vce         = vce_label,
               vce_select  = vce,
               c = c,
               h = h.lev,
               p = p,
               q = q,
               bwselect = bwselect,
               N = N.lev,
               Nh = Nh.lev,
               covs.cont = covs.cont,
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
print.rdbwhte <- function(x, ...){
  cat("MSE-Optimal Bandwidth Selection for RD Heterogeneous Treatment Effects Estimation\n")
  cat(paste("Number of Obs.           ",  format(sum(x$N), width=10, justify="right"), "\n", sep=""))
  cat(paste("Kernel                   ",  format(x$kernel, width=10, justify="right"), "\n", sep=""))
  cat(paste("VCE method               ",  format(x$vce,    width=10, justify="right"), "\n", sep=""))
  cat(paste("Poly. Order (p)          ",  format(x$p,      width=10, justify="right"), "\n", sep=""))
  cat("\n")

  invisible(x)
}


################################################################################
#' Summary method for \code{rdbwhte} objects.
#'
#' @param object Class \code{rdbwhte} object returned by \code{\link{rdbwhte}}.
#' @param ... ignored.
#'
#' @keywords internal
#' @return An S3 object inheriting from \code{rdbwhte} with class
#'   \code{summary.rdbwhte}; printed by \code{print.summary.rdbwhte}.
#' @export
summary.rdbwhte <- function(object, ...) {
  out <- object
  out$.summary <- TRUE
  class(out) <- c("summary.rdbwhte", class(object))
  out
}

#' Internal function.
#'
#' @param x Class \code{summary.rdbwhte} object.
#' @param ... ignored.
#'
#' @keywords internal
#' @return The input \code{x} returned invisibly; called for side effects.
#' @export
print.summary.rdbwhte <- function(x, ...) {
  cat(paste(x$rdmodel,"\n", sep=""))
  cat(paste("","\n", sep=""))

  cat(paste("Number of Obs.           ",  format(sum(x$N), width=10, justify="right"),"\n", sep=""))
  cat(paste("BW type                  ",  format(x$bwselect, width=10, justify="right"),"\n", sep=""))
  cat(paste("Kernel                   ",  format(x$kernel,   width=10, justify="right"),"\n", sep=""))
  cat(paste("VCE method               ",  format(x$vce,      width=10, justify="right"),"\n", sep=""))
  cat("\n")

  cat(paste("Number of Obs.           ",  format(x$N[1], width=10, justify="right"),  "   ", format(x$N[2], width=10, justify="right"), "\n", sep=""))
  if (isTRUE(x$covs.cont)) {
    cat(paste("Eff. Number of Obs.      ",  format(x$Nh[1], width=10, justify="right"),  "   ", format(x$Nh[2], width=10, justify="right"), "\n", sep=""))
  }
  cat(paste("Order est. (p)           ",  format(x$p, width=10, justify="right"),  "   ", format(x$p, width=10, justify="right"), "\n", sep=""))
  cat(paste("Order bias  (q)          ",  format(x$q, width=10, justify="right"),  "   ", format(x$q, width=10, justify="right"), "\n", sep=""))

  cat("\n")

  cat(paste(rep("=", 15 + 15 + 15), collapse="")); cat("\n")

  if (isFALSE(x$covs.cont)) {
    cat(format("Group", width=15, justify="right"))
  } else {
    cat(format("",      width=15, justify="right"))
  }
  cat(format("h-", width=15, justify="right"))
  cat(format("h+", width=15, justify="right"))
  cat("\n")

  cat(paste(rep("-", 15 + 15 + 15), collapse="")); cat("\n")

  if (isFALSE(x$covs.cont)) {
    for (i in 1:length(x$W.lev)) {
      cat(format(x$W.lev[i],                width=15, justify="right"))
      cat(format(sprintf("%3.3f",x$h[i,1]), width=15, justify="right"))
      cat(format(sprintf("%3.3f",x$h[i,2]), width=15, justify="right"))
      cat("\n")
    }
  } else {
    cat(format("Overall",                 width=15, justify="right"))
    cat(format(sprintf("%3.3f",x$h[1]),   width=15, justify="right"))
    cat(format(sprintf("%3.3f",x$h[2]),   width=15, justify="right"))
    cat("\n")
  }

  cat(paste(rep("=", 15 + 15 + 15), collapse="")); cat("\n")

  invisible(x)
}

#' Tidy a \code{rdbwhte} object
#'
#' broom-compatible \code{tidy()} method. Returns one row per heterogeneity
#' subgroup (factor case) or one row labelled \code{"Overall"} (continuous).
#'
#' @param x A \code{rdbwhte} object.
#' @param ... ignored.
#' @return A data frame with columns \code{term}, \code{h.left}, \code{h.right}.
#' @keywords internal
#' @exportS3Method broom::tidy
tidy.rdbwhte <- function(x, ...) {
  is_cont <- isTRUE(x$covs.cont)
  if (is_cont) {
    h_l <- x$h[1]; h_r <- x$h[2]
    terms <- "Overall"
  } else {
    h_l <- x$h[, 1]; h_r <- x$h[, 2]
    terms <- as.character(x$W.lev)
  }
  data.frame(
    term    = terms,
    h.left  = as.numeric(h_l),
    h.right = as.numeric(h_r),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

#' Glance at a \code{rdbwhte} object
#'
#' broom-compatible \code{glance()} method. Returns a one-row summary of the
#' bandwidth-selection call.
#'
#' @param x A \code{rdbwhte} object.
#' @param ... ignored.
#' @return A one-row data frame.
#' @keywords internal
#' @exportS3Method broom::glance
glance.rdbwhte <- function(x, ...) {
  data.frame(
    n        = sum(x$N),
    n.left   = x$N[1],
    n.right  = x$N[2],
    cutoff   = x$c,
    p        = x$p,
    q        = x$q,
    kernel   = x$kernel,
    vce      = x$vce,
    vce_select = x$vce_select,
    bwselect = x$bwselect,
    n.terms  = length(x$W.lev),
    covs.continuous = isTRUE(x$covs.cont),
    stringsAsFactors = FALSE
  )
}
