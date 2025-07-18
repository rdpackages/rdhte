################################################################################
#' @title RD Heterogeneous Treatment Effects Estimation and Inference
#'
#' @description \code{rdhte} provides estimation and inference for heterogeneous
#' treatment effects in RD designs using local polynomial regressions,
#' allowing for interactions with pretreatment covariates
#' (Calonico, Cattaneo, Farrell, Palomba and Titiunik, 2025a).
#' Inference is implemented using robust bias-correction methods
#' (Calonico, Cattaneo, and Titiunik, 2014)
#'
#' Companion commands: \code{\link{rdbwhte}} for data-driven bandwidth selection
#' and \code{\link{rdhte_lincom}} for testing linear restrictions of parameters.
#'
#' Related software packages for analysis and interpretation of RD designs and
#' related methods are available in: \url{https://rdpackages.github.io/}.
#'
#'For background methodology, see Calonico, Cattaneo, Farrell, and Titiunik
#'(2019), Calonico, Cattaneo and Farrell (2020), Cattaneo and Titiunik (2022).
#'
#' @param y Outcome variable.
#' @param x Running variable.
#' @param c RD cutoff in \code{x}; default is \code{c = 0}.
#' @param covs.hte covariates for heterogeneous treatment effects. Factor variables can be used to distinguish between continuous and categorical variables, select reference categories, specify interactions between variables, and include polynomials of continuous variables.
#' If not specified, the RD Average Treatment Effect is computed.
#' @param covs.eff additional covariates to be used for efficiency improvements.
#' @param p order of the local polynomial used to construct the point estimator (default = 1).
#' @param q order of the local polynomial used to construct the bias correction  (default = 2).
#' @param kernel kernel function used to construct the RD estimators. Options are \code{triangular} (default option), \code{epanechnikov} and \code{uniform}.
#' @param weights variable used for optional weighting of the estimation procedure. The unit-specific weights multiply the kernel function.
#' @param h main bandwidth used to construct the RD estimator. If not specified, bandwidth \code{h} is computed by the companion command \code{\link{rdbwhte}}. More than one bandwidth can be specified for categorical covariates.
#' @param h.l same as \code{h}, but only used for observations left of the cutoff \code{c}.
#' @param h.r same as \code{h}, but only used for observations right of the cutoff \code{c}.
#' @param vce character string specifying the variance-covariance matrix estimator type (hc0–hc3) (default = "hc3").
#' It is based on the \code{R} function \code{\link[sandwich]{vcovCL}}.
#' @param cluster variable indicating the clustering of observations.
#' @param level confidence level for confidence intervals; default is \code{level = 95}.
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
#' @param subset optional vector specifying a subset of observations to be used.

#'
#' @return A list with selected RD HTE effects and model information.
#' \item{Estimate}{vector of conventional local-polynomial RD estimates.}
#' \item{Estimate.bc}{vector of bias-corrected local-polynomial RD estimates.}
#' \item{se.rb}{vector containing robust bias corrected standard errors of the local-polynomial RD estimates.}
#' \item{ci.rb}{matrix containing robust bias corrected confidence intervals.}
#' \item{t.rb}{vector containing the t-statistics associated with robust local-polynomial RD estimates.}
#' \item{pv.rb}{vector containing the p-values associated with robust local-polynomial RD estimates.}
#' \item{coefs}{vector containing the coefficients for the jointly estimated p-th order local polynomial model.}
#' \item{vcov}{estimated variance-covariance matrix.}
#' \item{W.lev}{vector of group level identifiers.}
#' \item{kernel}{kernel type used.}
#' \item{vce}{variance estimator used.}
#' \item{c}{cutoff value.}
#' \item{h}{vector containing the bandwidths used.}
#' \item{p}{order of the polynomial used for estimation of the regression function.}
#' \item{q}{order of the polynomial used for inference on the regression function.}
#' \item{N}{vector with the original number of observations for each group.}
#' \item{Nh}{vector with the effective number of observations for each group.}
#' \item{covs.cont}{internal value.}
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
#' Calonico, Cattaneo, Farrell, Palomba and Titiunik (2025): \href{https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell-Palomba-Titiunik_2025_Stata.pdf}{rdhte: Learning Conditional Average Treatment Effects in RD Designs.} \emph{Working paper}.
#'
#' Calonico, Cattaneo, Farrell, Palomba and Titiunik (2025): \href{https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell-Palomba-Titiunik_2025_HTERD.pdf}{Treatment Effect Heterogeneity in Regression Discontinuity Designs.} \emph{Working paper}.
#'
#' Cattaneo, Farrell, and Titiunik. 2022. \href{https://rdpackages.github.io/references/Cattaneo-Titiunik_2022_ARE.pdf}{Regression Discontinuity Designs.} \emph{Annual Review of Economics},  14: 821-851.
#'
#' Calonico, Cattaneo, and Farrell. 2020. \href{https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell_2020_ECTJ.pdf}{Optimal Bandwidth Choice for Robust Bias Corrected Inference in Regression Discontinuity Designs.} \emph{Econometrics Journal}, 23(2): 192-210.
#'
#' Calonico, Cattaneo, Farrell, and Titiunik. 2019. \href{https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell-Titiunik_2019_RESTAT.pdf}{Regression Discontinuity Designs using Covariates.} \emph{Review of Economics and Statistics}, 101(3): 442-451.
#'
#' Calonico, Cattaneo, and Titiunik. 2014a. \href{https://rdpackages.github.io/references/Calonico-Cattaneo-Titiunik_2014_ECMA.pdf}{Robust Nonparametric Confidence Intervals for Regression-Discontinuity Designs.} \emph{Econometrica} 82(6): 2295-2326.
#'
#' @seealso \code{\link{rdbwhte}}, \code{\link{rdhte_lincom}}
#'
#' @examples
#' set.seed(123)
#' n <- 1000
#' X <- runif(n, -1, 1)
#' W <- rbinom(n, 1, 0.5)
#' Y <- 3 + 2*X + 1.5*X^2 + 0.5*X^3 + sin(2*X) + 3*W*(X>=0) + rnorm(n)
#' m1 = rdhte(y = Y, x = X, covs.hte = factor(W))
#' summary(m1)
#' @export
rdhte <- function(y, x, c = 0, covs.hte = NULL, covs.eff = NULL, 
                  p = 1, q = 2, kernel = "tri", weights = NULL, 
                  h = NULL, h.l = NULL, h.r = NULL,
                  vce = "hc3", cluster = NULL, level = 95, 
                  bwselect = NULL, bw.joint = FALSE, subset = NULL) {

  if (is.null(covs.hte)) {
    #warning("rdhte requires specifying heterogenity variables via covs.hte. RD ATE reported")

    # Create initial data frame and remove missing values
    dd <- data.frame(Y = y, X = x)
    if (!is.null(cluster))  dd$C    <- cluster
    if (!is.null(covs.eff)) dd$covs <- covs.eff
    if (!is.null(weights))  dd$weights <- weights
    if (!is.null(subset)) {
      dd$subset <- subset
      dd <- subset(dd, subset)
    }
    dd <- na.omit(dd)
    N  <- nrow(dd)

    # Extract relevant variables
    Xc  <- dd$X - c
    Y   <- dd$Y
    covs    <- if (!is.null(covs.eff)) dd$covs else NULL
    cluster <- if (!is.null(cluster))  dd$C    else NULL
    weights <- if (!is.null(weights))  dd$weights    else NULL
    
    T <- as.integer(Xc >= 0)

    if (!is.null(covs)) covs <- model.matrix(~covs-1)

    # Polynomial transformation of running variable (Xc)
    Xp <- poly(Xc, raw = TRUE, degree = p)
    Xq <- poly(Xc, raw = TRUE, degree = q)

    # Standardize kernel and variance estimator case
    kernel  <- tolower(kernel)
    vce     <- tolower(vce)
    vce.hte <- toupper(vce)

    # Factor or continuous covariate handling for W
    covs.cont <- TRUE
    W.lev <- "Overall"
    n.lev <- 1
    N.lev <- N
    cls <- NULL
    covs.hte_chr <- NULL
    W.names <-"T"
    
    # Bandwidth not provided
    if (is.null(h) & is.null(h.l) & is.null(h.r)) {
        if (is.null(bwselect)) bwselect="mserd"
        rd.bw  <- suppressWarnings(rdbwselect(y = Y, x = Xc, covs = covs, p = p, q = q, bwselect = bwselect, vce = vce, cluster = cluster, kernel = kernel))
        h.lev      <- c(rd.bw$bws[1], rd.bw$bws[2])
    } else if (is.null(h.r) & is.null(h.r))  {
        h.lev <- rep(h, 2)
        bwselect = "Manual"
    } else {
        h.lev <- c(h.l, h.r)
        bwselect = "Manual"
    }

    N.lev.l <- sum(T==0)
    N.lev.r <- sum(T==1)
    N.lev <- c(N.lev.l, N.lev.r)

    Nh.lev.l <- sum( (abs(Xc)<=h.lev[1]) & T==0)
    Nh.lev.r <- sum( (abs(Xc)<=h.lev[2]) & T==1)
    Nh.lev <- matrix(c(Nh.lev.l, Nh.lev.r),1,2)


    # Weighted Regression
    #r.bw <- abs(Xc) <= h
    r.bw <- (abs(Xc) <= h.lev[1] & T==0) |  (abs(Xc) <= h.lev[2] & T==1)
    k.weights = NULL
    Kernel = "Uniform"
    if (kernel=="tri") {
      k.weights = (1-abs(Xc)/h.lev[1])*(T==0) + (1-abs(Xc)/h.lev[2])*(T==1)
      Kernel = "Triangular"
    }
    if (kernel=="epa") {
      k.weights = (1-(Xc/h.lev[1])^2)*(T==0) + (1-(Xc/h.lev[2])^2)*(T==1)
      Kernel = "Epanechnikov"
    }
    if (is.null(covs.eff)) {
      reg_formula_est <- Y ~ T * Xp
      reg_formula_inf <- Y ~ T * Xq
    } else {
      reg_formula_est <- Y ~ T * Xp * W + covs * W
      reg_formula_inf <- Y ~ T * Xq * W + covs * W
    }
    
    if (!is.null(weights))  k.weights <- k.weights*weights    

    # Run a single regression model with weights
    rd_est <- lm(reg_formula_est, weights = k.weights, subset = r.bw)
    rd_inf <- lm(reg_formula_inf, weights = k.weights, subset = r.bw)

    # Compute robust variance-covariance matrix
    rd.vcov.full <- vcovCL(rd_inf, cluster = cluster[r.bw], type = vce.hte)

    # Extract treatment effect estimate
    tau.hat    <- rd_est$coefficients["T"]
    tau.hat.bc <- rd_inf$coefficients["T"]
    tau.hat.se <- sqrt(rd.vcov.full["T","T"])


  ###############################################################################################

  } else {
    
    y_chr        <- deparse(substitute(y))
    x_chr        <- deparse(substitute(x))
    covs.hte_chr <- deparse(substitute(covs.hte))
    
    #print(covs.hte_chr)
     
    ### Chec if covs.hte is a formula:
    is_char <- is.character(covs.hte)
    cls     <- if (is_char) "formula" else "operation"
    # validate formulas
    if (cls == "formula") {
      txt <- covs.hte[1]
      if (!startsWith(trimws(txt), "~"))
        fml.txt <- paste0("~", txt, "-1")             
        f <- as.formula(fml.txt)
        covs.hte.mm <- model.matrix(f)
        dd <- data.frame(Y = y, X = x, covs.hte = covs.hte.mm)
    } else {
      dd <- data.frame(Y = y, X = x, covs.hte = covs.hte)
    }
    
  if (!is.null(cluster))  dd$C    <- cluster
  if (!is.null(covs.eff)) dd$covs <- covs.eff
  if (!is.null(weights))  dd$weights <- weights
    
  if (!is.null(subset)) {
      dd$subset <- subset
      dd <- subset(dd, subset)
  }
  dd <- na.omit(dd)
  N  <- nrow(dd)

  # Extract relevant variables
  Xc  <- dd$X - c
  Y   <- dd$Y
  idx <- grep("^covs.hte", names(dd))
  names(dd) <- gsub("covs.hte.", "", names(dd))
  
  #print(head(dd))
  
  W.covs   <- dd[, idx]
  if (cls == "formula") W.names <- names(dd[, idx])
  covs    <- if (!is.null(covs.eff)) dd$covs else NULL
  cluster <- if (!is.null(cluster))  dd$C    else NULL
  weights <- if (!is.null(weights))  dd$weights    else NULL
  T <- as.integer(Xc >= 0)

  if (!is.null(covs))  covs <- model.matrix(~covs-1)

  # Polynomial transformation of running variable (Xc)
  Xp <- poly(Xc, raw = TRUE, degree = p)
  Xq <- poly(Xc, raw = TRUE, degree = q)

  # Standardize kernel and variance estimator case
  kernel  <- tolower(kernel)
  vce     <- tolower(vce)
  vce.hte <- toupper(vce)

  # Factor or continuous covariate handling for W
  covs.cont <- TRUE
  W.lev <- NULL
  n.lev <- 1
  N.lev <- 0

  if (is.factor(W.covs)) {
    W.lev <- levels(W.covs)
    n.lev <- nlevels(W.covs)
    N.lev <- table(W.covs)
    covs.cont <- FALSE
  } else {
    if (length(idx) == 1) {
      cond <- mean((W.covs == 0) | (W.covs == 1))
        if (cond == 1 ) {
          W.covs <- factor(W.covs)
          W.lev <- levels(W.covs)
          n.lev <- nlevels(W.covs)
          N.lev <- table(W.covs)
          covs.cont <- FALSE
        }
    } else {
    cond1 <- mean((rowSums(W.covs) == 1) | (ncol(W.covs) == 1))
    cond2 <- mean((W.covs == 0) | (W.covs == 1))
    if (cond1 == 1 & cond2 == 1) {
      col_names <- names(W.covs)
      W.covs     <- factor(apply(W.covs[, col_names], 1, function(x) which(x == 1)))
      W.lev <- levels(W.covs)
      n.lev <- nlevels(W.covs)
      N.lev <- table(W.covs)
      covs.cont <- FALSE
    }
  }
}




  # Initialize bandwidth storage
  h.vec  <- rep(0, N)
  h.lev  <- matrix(0, n.lev, 2)
  Nh.lev <- matrix(NA, n.lev, 2)

  # Bandwidth is provided (1 or many)
  if (!is.null(h)) {
    bwselect = "Manual"

    if (n.lev>1){

      if (length(h) == 1) {
        h.vec <- matrix(h, N, 1)
        h.lev <- matrix(h, n.lev, 2)

        for (l in 1:n.lev) {
          ind        <- W.covs == W.lev[l]
          Nh.lev[l,1]  <- sum( (abs(Xc[ind]) <= h) & (T[ind]==0))
          Nh.lev[l,2]  <- sum( (abs(Xc[ind]) <= h) & (T[ind]==1))
        }

      } else {

        hdif = abs(n.lev - length(h))
        if (hdif > 0) stop("check the number of bandwidths provided")

        for (l in 1:n.lev) {
          h.lev[l,1] <- h[l]
          h.lev[l,2] <- h[l]

          ind        <- W.covs == W.lev[l]
          h.vec[ind] <- h[l]
          Nh.lev[l,1]  <- sum( (abs(Xc[ind]) <= h.lev[l]) & (T[ind]==0))
          Nh.lev[l,2]  <- sum( (abs(Xc[ind]) <= h.lev[l]) & (T[ind]==1))
        }
      }
    } else {
      h.vec <- matrix(h, N, 1)
      h.lev <- matrix(h, n.lev, 2)
      Nh.lev[1,1]  <- sum( (abs(Xc) <= h & T==0) )
      Nh.lev[1,2]  <- sum( (abs(Xc) <= h & T==1) )
      }

  } else if (!is.null(h.l) & !is.null(h.r))  {

    bwselect = "Manual"

    if (n.lev>1){

    if (length(h.l) == 1 & length(h.r) == 1) {

      ind.l <- Xc<0
      ind.r <- Xc>=0
      h.vec[ind.l] <- h.l
      h.vec[ind.r] <- h.r
      h.l.lev <- matrix(h.l, n.lev, 1)
      h.r.lev <- matrix(h.r, n.lev, 1)
      h.lev <- cbind(h.l.lev, h.r.lev)

      for (l in 1:n.lev) {
        h.lev[l,1] <- h.l
        h.lev[l,2] <- h.r

        ind.l     <- W.covs == W.lev[l] & T==0
        ind.r     <- W.covs == W.lev[l] & T==1
        h.vec[ind.l] <- h.l
        h.vec[ind.r] <- h.r

        Nh.lev[l,1]  <- sum( (abs(Xc[ind.l]) <= h.l) )
        Nh.lev[l,2]  <- sum( (abs(Xc[ind.r]) <= h.r) )
      }

    } else {

      hdif.l = abs(n.lev - length(h.l))
      hdif.r = abs(n.lev - length(h.r))
      if ((hdif.l > 0 | hdif.r > 0) & (length(h.l)>1 & length(h.r)))     stop("check the number of bandwidths provided")

      for (l in 1:n.lev) {
        h.lev[l,1] <- h.l[l]
        h.lev[l,2] <- h.r[l]

        ind.l     <- W.covs == W.lev[l] & T==0
        ind.r     <- W.covs == W.lev[l] & T==1
        h.vec[ind.l] <- h.l[l]
        h.vec[ind.r] <- h.r[l]

        Nh.lev[l,1]  <- sum( (abs(Xc[ind.l]) <= h.l[l]) )
        Nh.lev[l,2]  <- sum( (abs(Xc[ind.r]) <= h.r[l]) )
      }
    }

    } else {

      ind.l <- Xc<0
      ind.r <- Xc>=0
      h.vec[ind.l] <- h.l
      h.vec[ind.r] <- h.r
      Nh.lev[1,1]  <- sum( (abs(Xc[ind.l]) <= h.l) )
      Nh.lev[1,2]  <- sum( (abs(Xc[ind.r]) <= h.r) )

    }
  }



  # Bandwidth not provided

  if (is.null(h) & is.null(h.l) & is.null(h.r)) {

    if (is.null(bwselect)) bwselect="mserd"

    if (bw.joint==TRUE | is.null(W.lev)) {
      # Joint bandwidth estimation
      rd.bw  <- suppressWarnings(rdbwselect(y = Y, x = Xc, covs = covs, p = p, q = q, bwselect = bwselect, vce = vce, cluster = cluster, kernel = kernel))

      h.lev[,1]  <- rep(rd.bw$bws[1], n.lev)
      h.lev[,2]  <- rep(rd.bw$bws[2], n.lev)

      ind.l <- Xc<0
      ind.r <- Xc>=0
      h.vec[ind.l] <- rd.bw$bws[1]
      h.vec[ind.r] <- rd.bw$bws[2]

      if (n.lev>1) {
        for (l in 1:n.lev) {
          ind       <- which(W.covs==W.lev[l])
          #Nh.lev[l] <- sum(abs(Xc[ind])<=h[l])
          Nh.lev[l,1]  <- sum( (abs(Xc[ind]) <= h.lev[l]) & (T[ind]==0))
          Nh.lev[l,2]  <- sum( (abs(Xc[ind]) <= h.lev[l]) & (T[ind]==1))
        }
      }
    } else {
        for (l in 1:n.lev) {
          rd.bw      <- suppressWarnings(rdbwselect(y = Y, x = Xc, covs = covs, p = p, q = q, bwselect = bwselect, vce=vce, cluster = cluster, kernel = kernel, subset = W.covs==W.lev[l]))
          h.lev[l,1]       <- rd.bw$bws[1]
          h.lev[l,2]       <- rd.bw$bws[2]

          ind.l     <- W.covs == W.lev[l] & T==0
          ind.r     <- W.covs == W.lev[l] & T==1
          h.vec[ind.l] <- rd.bw$bws[1]
          h.vec[ind.r] <- rd.bw$bws[2]

          Nh.lev[l,1]  <- sum( (abs(Xc[ind.l]) <= rd.bw$bws[1] ))
          Nh.lev[l,2]  <- sum( (abs(Xc[ind.r]) <= rd.bw$bws[2] ))

        }
      }
    }

  if (is.null(W.lev)) {
    W.covs      <- as.matrix(W.covs)
    N.lev.l <- sum(T==0)
    N.lev.r <- sum(T==1)
    N.lev <- c(N.lev.l, N.lev.r)

    Nh.lev.l <- sum( (abs(Xc)<=h[1]) & T==0)
    Nh.lev.r <- sum( (abs(Xc)<=h[2]) & T==1)
    Nh.lev <- c(Nh.lev.l, Nh.lev.r)
  }


  # Weighted Regression
    r.bw <- abs(Xc) <= h.vec
    k.weights = NULL
    Kernel = "Uniform"
    if (kernel=="tri") {
      k.weights = 1-abs(Xc)/h.vec
      Kernel = "Triangular"
      }
    if (kernel=="epa") {
      k.weights = 1-(Xc/h.vec)^2
      Kernel = "Epanechnikov"
    }
    if (is.null(covs.eff)) {
      reg_formula_est <- Y ~ T * Xp * W.covs
      reg_formula_inf <- Y ~ T * Xq * W.covs
    } else {
      reg_formula_est <- Y ~ T * Xp * W.covs + covs * W.covs
      reg_formula_inf <- Y ~ T * Xq * W.covs + covs * W.covs
    }

    if (!is.null(weights))  k.weights <- k.weights*weights    
    
    
    # Run a single regression model with weights
    rd_est <- lm(reg_formula_est, weights = k.weights, subset = r.bw)
    rd_inf <- lm(reg_formula_inf, weights = k.weights, subset = r.bw)

    # Compute robust variance-covariance matrix
    rd.vcov.full <- vcovCL(rd_inf, cluster = cluster[r.bw], type = vce.hte)

    # Extract treatment effect estimate
    tau.hat    = rd_est$coefficients["T"]
    tau.hat.bc = rd_inf$coefficients["T"]
    tau.hat.se = sqrt(rd.vcov.full["T","T"])

    
    
    
    if (!is.factor(W.covs)) {
      ncoeff = names(rd_est$coefficients[!is.na(rd_est$coefficients)])
      W.lev = c("T", ncoeff[grepl("T:W.covs", ncoeff)])
      n.lev <- length(W.lev)
      W.names <- c("T", NA) 
      for (j in 2:n.lev) {
        lev <- W.lev[j]
        W.names[j] <- lev
        tau.hat[j]    <- rd_est$coefficients[lev]
        tau.hat.bc[j] <- rd_inf$coefficients[lev]
        tau.hat.se[j] <- sqrt(rd.vcov.full[lev,lev])
      }
      
  
      ### Covariances
      vcov     = diag(tau.hat.se^2)
      
      for (i in 1:n.lev) {
        for (j in i:n.lev) {
					lev_i <- W.lev[i]
					lev_j <- W.lev[j]
          tmp <- rd.vcov.full[lev_i, lev_j]
          vcov[i, j] = vcov[j, i] = tmp
        }				
      }
  
  
  
    } else {
      vcov     = matrix(NA, n.lev, n.lev)
      vcov[1,1]     = tau.hat.se[1]^2
      
    for (j in 2:n.lev) {
      
        lev = paste("T:W.covs", W.lev[j], sep="")
        tau.hat[j]    = rd_est$coefficients["T"] + rd_est$coefficients[lev]
        tau.hat.bc[j] = rd_inf$coefficients["T"] + rd_inf$coefficients[lev]
        tau.hat.se[j] = sqrt(rd.vcov.full["T","T"] + rd.vcov.full[lev,lev] + 2*rd.vcov.full["T",lev])
        
        vcov[j,j]     = tau.hat.se[j]^2
        
        
        tmp = rd.vcov.full["T","T"] + rd.vcov.full["T", lev]
        if (j>1 & j <= n.lev) {
          vcov[1, j] = vcov[j, 1] = tmp
        }
        
    }
      
			### Covariances
			if (n.lev > 2) {
			  for (i in 2:n.lev) {
				  k = i + 1 
				  if (k <= n.lev) {
				    for (j in k:n.lev) {
    				  lev_i <- paste("T:W.covs", W.lev[i], sep="")
    					lev_j <- paste("T:W.covs", W.lev[j], sep="")		
    					
    				  tmp <- rd.vcov.full["T","T"] + 
    				    rd.vcov.full["T", lev_i] +  
    				    rd.vcov.full["T", lev_j] + 
    				    rd.vcov.full[lev_i, lev_j] 
    				  
    				    vcov[i, j] = vcov[j, i] = tmp		
				    }
				}				
      }
    }	
  }
    
  

    # Construct results table
    #results_df <- data.frame(
    #  Group       = W.lev,
    #  Estimate    = tau.hat,
    #  Estimate.bc = tau.hat.bc,
    #  SE          = tau.hat.se,
    #  p.value     = tau.hat.pv
    #)

    
    if (cls != "formula") {
      if (covs.cont==TRUE) {
        W.names = c("T", paste("T:",covs.hte_chr, sep=""))
        #W.names <- paste(covs.hte_chr,W.lev, sep="")  
      } else {
        W.names <- paste(covs.hte_chr,W.lev, sep="")  
      }
    } else {
      W.names = W.lev = gsub("W.covs", "", W.lev)
    }
    
    #print(W.lev)
    #print(W.names)
    
    names(tau.hat) = names(tau.hat.bc) = names(tau.hat.se) = W.names
    if (length(tau.hat.se)>1) colnames(vcov) = rownames(vcov) = W.names
    names(tau.hat) <- gsub("W.covs.", "", names(tau.hat))
    
  }

  qz <- -qnorm(abs((1-(level/100))/2))
  t.rb = tau.hat.bc/tau.hat.se
  pv.rb = 2*pnorm(-abs(t.rb))

    rdmodel <- "Sharp RD Average Treatment Effect."
    if (!is.null(covs.hte)) rdmodel <- "Sharp RD Heterogeneous Treatment Effects: Subgroups."
    if (!is.null(covs.hte) & covs.cont=="TRUE") rdmodel <- "Sharp RD Heterogeneous Treatment Effects: Continuous."
    if (!is.null(covs.eff)) rdmodel <-paste(rdmodel, ", Covariate-Adjusted", sep="")
    if (!is.null(cluster))  rdmodel <-paste(rdmodel, ", Cluster-Adjusted", sep="")

  
    #if (cls == "formula") {
      #W.names = c("T", paste("T:",W.names, sep=""))
    #} else {
    #  W.names <- paste(covs.hte_chr,W.lev, sep="")  
    #}
    


  # Store model information
  out = list(  Estimate    = tau.hat,
               Estimate.bc = tau.hat.bc,
               se.rb       = tau.hat.se,
               coef        = tau.hat,
               coef.bc     = tau.hat.bc,
               vcov        = vcov,
               ci.rb       = cbind(tau.hat.bc - qz*tau.hat.se, tau.hat.bc + qz*tau.hat.se),
               t.rb        = t.rb,
               pv.rb       = pv.rb,
               coef.full   = rd_est$coefficients,
               vcov.full   = rd.vcov.full,
               W.lev       = W.lev,
               W.names.    = W.names,
               covs.hte_chr = covs.hte_chr,
               kernel      = Kernel,
               bwselect    = bwselect,
               vce         = vce.hte,
               c = c,
               h = h.lev,
               p = p,
               q = q,
               N = N.lev,
               Nh = Nh.lev,
               covs.cont = covs.cont,
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

  name.long = max(12,nchar(x$covs.hte_chr))
  
  #llenght <- 15 + 15 + 15 + 15 + 25 +  10 + 10 + 10 + 10
  llenght <- name.long + 12 + 12 + 12 + 24 +  10 + 10 + 10 + 10
  

  #cat("Call: rdrobust\n\n")
  cat(paste("Number of Obs.           ",  format(sum(x$N),   width=10, justify="right"),"\n", sep=""))
  cat(paste("BW type                  ",  format(x$bwselect, width=10, justify="right"),"\n", sep=""))
  cat(paste("Kernel                   ",  format(x$kernel,   width=10, justify="right"),"\n", sep=""))
  cat(paste("VCE method               ",  format(x$vce,      width=10, justify="right"),"\n", sep=""))
  cat("\n")

  cat(paste("Number of Obs.           ",  format(x$N[1],      width=10, justify="right"),  "   ", format(x$N[2],   width=10, justify="right"),       "\n", sep=""))
  if (x$covs.cont==TRUE) {
    cat(paste("Eff. Number of Obs.      ",  format(x$Nh[1], width=10, justify="right"),  "   ", format(x$Nh[2], width=10, justify="right"),       "\n", sep=""))
    llenght <- name.long + 12 + 12 + 12 + 24
  }
  cat(paste("Order est. (p)           ",  format(x$p,      width=10, justify="right"),  "   ", format(x$p,      width=10, justify="right"),       "\n", sep=""))
  cat(paste("Order bias  (q)          ",  format(x$q,      width=10, justify="right"),  "   ", format(x$q,      width=10, justify="right"),       "\n", sep=""))
  if (x$covs.cont==TRUE) {
    cat(paste("BW est. (h)              ",  format(sprintf("%10.3f",x$h[1])),  "   ", format(sprintf("%10.3f",x$h[2])),      "\n", sep=""))
  }




  cat("\n")

  ### print output
  cat(paste(rep("=", llenght ), collapse="")); cat("\n")
  cat(format(""                    , width=name.long, justify="right"))
  cat(format("Point"               , width=12, justify="right"))
  cat(format("Robust Inference"    , width=24, justify="right"))

  cat("\n")

  if (x$covs.cont==FALSE) {
    cat(format(x$covs.hte_chr,                                              width=name.long, justify="right"))
    #cat(format("Group"            , width=12, justify="right"))
    cat(format("Estimate"         , width=12, justify="right"))
  } else {
    cat(format(""                 , width=name.long, justify="right"))
    cat(format("Estimate"         , width=12, justify="right"))
  }

  cat(format("z"              , width=12, justify="right"))
  cat(format("Pr(>|z|)"       , width=12, justify="right"))
  cat(format(paste("[ ", x$level, "% ", "C.I. ]", sep=""), width=24, justify="centre"))

  if (x$covs.cont==FALSE) {
    cat(format("Nh-"       , width=10, justify="right"))
    cat(format("Nh+"       , width=10, justify="right"))
    cat(format("h-"        , width=10, justify="right"))
    cat(format("h+"        , width=10, justify="right"))

  }
  cat("\n")

  cat(paste(rep("-", llenght ), collapse="")); cat("\n")

  for (i in 1:length(x$W.lev)) {
    if (x$covs.cont==TRUE) cat(format(x$W.names[i],                              width=name.long, justify="right"))
    else   cat(format(x$W.lev[i],                                                width=name.long, justify="right"))
    cat(format(sprintf("%3.3f", x$Estimate[i]),                           width=12, justify="right"))
    cat(format(sprintf("%3.3f", x$t.rb[i]),                               width=12, justify="right"))
    cat(format(sprintf("%3.3f", x$pv.rb[i]),                              width=12, justify="right"))
    cat(format(paste("[", sprintf("%3.3f", x$ci.rb[i,1]), " , ", sep=""), width=12, justify="right"))
    cat(format(paste(     sprintf("%3.3f", x$ci.rb[i,2]), "]",   sep=""), width=12, justify="left"))

    if (x$covs.cont==FALSE) {
      cat(format(x$Nh[i,1],                   width=10, justify="right"))
      cat(format(x$Nh[i,2],                   width=10, justify="right"))
      cat(format(sprintf("%3.3f",x$h[i,1]),   width=10, justify="right"))
      cat(format(sprintf("%3.3f",x$h[i,2]),   width=10, justify="right"))
    }
    cat("\n")
  }

  cat(paste(rep("=", llenght ), collapse="")); cat("\n")
}



