 ################################################################################
#' @title Data-Driven Optimal Bandwidth Selection for RD Heterogeneous Treatment
#' Effects Estimation
#'
#' @description \code{rdbwhte} computes MSE- and CER-optimal bandwidths for estimating
#' RD heterogeneous treatment effects based on covariates.
#'
#' Companion commands: \code{\link{rdhte}} for RD HTE estimation and inference,
#' and \code{\link{rdhte_lincom}} for testing linear restrictions of parameters.
#'
#' Related Stata and R packages useful for inference in RD designs are
#' described in the website: \url{https://rdpackages.github.io/}.
#'
#' @param y Outcome variable.
#' @param x Running variable.
#' @param c RD cutoff in \code{x}; default is \code{c = 0}.
#' @param covs.hte covariates for heterogeneous treatment effects. Factor variables can be used to distinguish between continuous and categorical variables, select reference categories, specify interactions between variables, and include polynomials of continuous variables.
#' @param covs.eff additional covariates to be used for efficiency improvements.
#' @param p order of the local polynomial used to construct the point estimator (default = 1).
#' @param q order of the local polynomial used to construct the bias correction  (default = 2).
#' @param kernel kernel function used to construct the RD estimators. Options are \code{triangular} (default option), \code{epanechnikov} and \code{uniform}.
#' @param vce character string specifying the variance-covariance matrix estimator type (hc0–hc3) (default = "hc3").
#' @param cluster variable indicating the clustering of observations.
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
#' @param bw.joint logical. If \code{TRUE}, forces all bandwidths to be the same across groups (default is \code{bw.joint = FALSE}).
#'
#' @return A list with selected bandwidths and model information.
#' \item{W.lev}{vector of group level identifiers.}
#' \item{kernel}{kernel type used.}
#' \item{vce}{variance estimator used.}
#' \item{c}{cutoff value.}
#' \item{h}{vector containing the bandwidths used.}
#' \item{p}{order of the polynomial used for estimation.}
#' \item{q}{order of the polynomial used for inference.}
#' \item{N}{vector with the original number of observations for each group.}
#' \item{Nh}{vector with the effective number of observations for each group.}
#' \item{covs.cont}{internal value.}
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
#' Calonico, Cattaneo, Farrell, Palomba and Titiunik (2025): Treatment Effect Heterogeneity in Regression Discontinuity Designs. \emph{Working paper}
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
#' @export
rdbwhte <- function(y, x, c = 0, covs.hte = NULL, covs.eff = NULL, p = 1, q = 2,
                    kernel = "tri", vce = "hc3", cluster = NULL, bwselect = "mserd", bw.joint = FALSE) {

  if (is.null(covs.hte)) {
   # warning("rdbwhte requires specifying heterogenity variables via covs.hte.
   #         MSE-Optimal Bandwidth for Overall RD ATE reported")

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
    #Xp <- poly(Xc, raw = TRUE, degree = p)
    #Xq <- poly(Xc, raw = TRUE, degree = q)

    # Standardize kernel and variance estimator case
    kernel  <- tolower(kernel)
    vce     <- tolower(vce)
    vce.hte <- toupper(vce)

    # Factor or continuous covariate handling for W
    covs.cont <- FALSE
    W.lev <- "Overall"
    n.lev <- 1
    N.lev <- N




      rd.bw  <- suppressWarnings(rdbwselect(y = Y, x = Xc, covs = covs, p = p, q = q, bwselect = bwselect, vce = vce, cluster = cluster, kernel = kernel))
      h.lev      <- matrix(c(rd.bw$bws[1], rd.bw$bws[2]),1,2)
      #Nh.lev <- sum(abs(Xc)<=h)
      N.lev.l <- sum(T==0)
      N.lev.r <- sum(T==1)
      N.lev <- c(N.lev.l, N.lev.r)
      Nh.lev.l <- sum( (abs(Xc)<=h.lev[1]) & T==0)
      Nh.lev.r <- sum( (abs(Xc)<=h.lev[2]) & T==1)
      Nh.lev <- matrix(c(Nh.lev.l, Nh.lev.r),1,2)

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
  if (!is.null(covs)) covs <- model.matrix(~covs-1)

  # Polynomial transformation of running variable (Xc)
  #Xp <- poly(Xc, raw = TRUE, degree = p)
  #Xq <- poly(Xc, raw = TRUE, degree = q)

  # Standardize kernel and variance estimator case
  kernel  <- tolower(kernel)
  vce     <- tolower(vce)
  vce.hte <- toupper(vce)

  # Factor or continuous covariate handling for W
  covs.cont <- TRUE
  W.lev <- NULL
  n.lev <- 1
  N.lev <- 0

  if (is.factor(W)) {
    W.lev <- levels(W)
    n.lev <- nlevels(W)
    N.lev <- table(W)
    covs.cont <- FALSE
  } else {
    if (length(idx) == 1) {
      cond <- mean((W == 0) | (W == 1))
      if (cond == 1 ) {
        W <- factor(W)
        W.lev <- levels(W)
        n.lev <- nlevels(W)
        N.lev <- table(W)
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
        N.lev <- table(W)
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

    if (bw.joint==TRUE | is.null(W.lev)) {
      # Joint bandwidth estimation
      rd.bw  <- suppressWarnings(rdbwselect(y = Y, x = Xc, covs = covs, p = p, q = q, bwselect = bwselect, vce = vce, cluster = cluster, kernel = kernel))
      #h      <- rep(rd.bw$bws[1], n.lev)
      #h.vec  <- rep(rd.bw$bws[1], N)

      h.lev[,1]  <- rep(rd.bw$bws[1], n.lev)
      h.lev[,2]  <- rep(rd.bw$bws[2], n.lev)

      ind.l <- Xc<0
      ind.r <- Xc>=0
      h.vec[ind.l] <- rd.bw$bws[1]
      h.vec[ind.r] <- rd.bw$bws[2]



      h.lev[,1]  <- rep(rd.bw$bws[1], n.lev)
      h.lev[,2]  <- rep(rd.bw$bws[2], n.lev)

      if (n.lev>1) {
        for (l in 1:n.lev) {
          ind       <- which(W==W.lev[l])
          #Nh.lev[l] <- sum(abs(Xc[ind])<=h[l])
          Nh.lev[l,1]  <- sum( (abs(Xc[ind]) <= h.lev[l]) & (T[ind]==0))
          Nh.lev[l,2]  <- sum( (abs(Xc[ind]) <= h.lev[l]) & (T[ind]==1))
        }
      }
    } else {
      #h = NULL
      for (l in 1:n.lev) {
        rd.bw      <- suppressWarnings(rdbwselect(y = Y, x = Xc, covs = covs, p = p,  q = q, bwselect = bwselect, vce=vce, cluster = cluster, kernel = kernel, subset=W==W.lev[l]))
        #h[l]       <- rd.bw$bws[1]
        #ind        <- which(W==W.lev[l])
        #h.vec[ind] <- h[l]
        #Nh.lev[l]  <- sum( abs(Xc[ind])<=h[l]  )

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
    #N.lev  <- N
    #Nh.lev <- sum(abs(Xc)<=h)

    N.lev.l <- sum(T==0)
    N.lev.r <- sum(T==1)
    N.lev <- c(N.lev.l, N.lev.r)

    Nh.lev.l <- sum( (abs(Xc)<=h.lev[1]) & T==0)
    Nh.lev.r <- sum( (abs(Xc)<=h.lev[2]) & T==1)
    Nh.lev <- c(Nh.lev.l, Nh.lev.r)


  }

}


  rdmodel <- "MSE-Optimal Bandwidth Selection for RD Heterogeneous Treatment Effects Estimation"
  if (!is.null(covs.eff)) rdmodel <-paste(rdmodel, ", Covariate-Adjusted", sep="")
  if (!is.null(cluster))  rdmodel <-paste(rdmodel, ", Cluster-Adjusted", sep="")


  # Store model information
  out = list(  W.lev       = W.lev,
               kernel      = rd.bw$kernel,
               vce         = vce.hte,
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

  #cat("Call: rdrobust\n\n")
  cat(paste("Number of Obs.           ",  format(sum(x$N), width=10, justify="right"),"\n", sep=""))
  cat(paste("BW type                  ",  format(x$bwselect, width=10, justify="right"),"\n", sep=""))
  cat(paste("Kernel                   ",  format(x$kernel,   width=10, justify="right"),"\n", sep=""))
  cat(paste("VCE method               ",  format(x$vce,      width=10, justify="right"),"\n", sep=""))
  cat("\n")

  cat(paste("Number of Obs.           ",  format(x$N[1],   width=10, justify="right"),  "   ", format(x$N[2],   width=10, justify="right"),       "\n", sep=""))
  if (x$covs.cont==TRUE) {
    cat(paste("Eff. Number of Obs.      ",  format(x$Nh[1], width=10, justify="right"),  "   ", format(x$Nh[2], width=10, justify="right"),       "\n", sep=""))
    llenght <- 15 + 15 + 15 + 15 + 30
  }
  cat(paste("Order est. (p)           ",  format(x$p,      width=10, justify="right"),  "   ", format(x$p,      width=10, justify="right"),       "\n", sep=""))
  cat(paste("Order bias  (q)          ",  format(x$q,      width=10, justify="right"),  "   ", format(x$q,      width=10, justify="right"),       "\n", sep=""))
#  if (x$covs.cont==TRUE) {
 #   cat(paste("BW est. (h)              ",  format(sprintf("%10.3f",x$h[1])),  "   ", format(sprintf("%10.3f",x$h[2])),      "\n", sep=""))
  #}

  cat("\n")




  ### print output
  cat(paste(rep("=", 15 + 15 +15), collapse="")); cat("\n")

  if (x$covs.cont==FALSE) {
    cat(format("Group"            , width=15, justify="right"))
  } else {
    cat(format(""                 , width=15, justify="right"))
  }

  cat(format("h-"        , width=15, justify="right"))
  cat(format("h+"        , width=15, justify="right"))
  cat("\n")

  cat(paste(rep("-", 15 + 15 + 15 ), collapse="")); cat("\n")

  if (x$covs.cont==FALSE) {
  for (i in 1:length(x$W.lev)) {
    cat(format(x$W.lev[i],              width=15, justify="right"))
    cat(format(sprintf("%3.3f",x$h[i,1]), width=15, justify="right"))
    cat(format(sprintf("%3.3f",x$h[i,2]), width=15, justify="right"))
    cat("\n")
  }
  } else {
    cat(format("Overall"              , width=15, justify="right"))
    cat(format(sprintf("%3.3f",x$h[1]), width=15, justify="right"))
    cat(format(sprintf("%3.3f",x$h[2]), width=15, justify="right"))
    cat("\n")
  }

  cat(paste(rep("=", 15 + 15 +15), collapse="")); cat("\n")

}



