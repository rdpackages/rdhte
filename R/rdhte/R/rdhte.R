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
#' @param y Outcome variable.
#' @param x Running variable.
#' @param c RD cutoff in \code{x}; default is \code{c = 0}.
#' @param covs.hte covariates for heterogeneous treatment effects. Factor variables can be used to distinguish between continuous and categorical variables, select reference categories, specify interactions between variables, and include polynomials of continuous variables.
#' If not specified, the RD Average Treatment Effect is computed.
#' @param covs.eff additional covariates to be used for efficiency improvements.
#' @param p order of the local polynomial used to construct the point estimator (default = 1).
#' @param q order of the local polynomial used to construct the bias correction.
#'   If \code{NULL} (default), \code{q} is set to \code{p + 1}.
#' @param kernel kernel function used to construct the RD estimators. Options are \code{triangular} (default option), \code{epanechnikov} and \code{uniform}.
#' @param weights variable used for optional weighting of the estimation procedure. The unit-specific weights multiply the kernel function.
#' @param h main bandwidth used to construct the RD estimator. If not specified, bandwidth \code{h} is computed by the companion command \code{\link{rdbwhte}}. More than one bandwidth can be specified for categorical covariates.
#' @param h.l same as \code{h}, but only used for observations left of the cutoff \code{c}.
#' @param h.r same as \code{h}, but only used for observations right of the cutoff \code{c}.
#' @param vce character string specifying the variance-covariance matrix
#'   estimator type. Without \code{cluster}: \code{"hc0"}, \code{"hc1"},
#'   \code{"hc2"}, \code{"hc3"} (default \code{"hc3"}). With
#'   \code{cluster}: \code{"cr1"} (default; standard cluster-robust
#'   sandwich with small-sample correction), \code{"cr2"}
#'   (Bell-McCaffrey leverage-adjusted), \code{"cr3"} (block-jackknife).
#'   Legacy aliases:
#'   \code{"hc0"}/\code{"hc1"} + \code{cluster} are remapped to
#'   \code{"cr1"} with a warning; \code{"hc2"} -> \code{"cr2"} and
#'   \code{"hc3"} -> \code{"cr3"} similarly. \code{"cr1"}, \code{"cr2"},
#'   \code{"cr3"} without \code{cluster} fall back to \code{"hc1"},
#'   \code{"hc2"}, \code{"hc3"} with a warning.
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
#' @param bw.joint logical. If \code{TRUE}, forces all bandwidths to be the same across groups (default is \code{bw.joint = FALSE}). When \code{covs.hte} is continuous (rather than a factor or 0/1 indicator), a single joint bandwidth is always used regardless of this argument.
#' @param subset optional vector specifying a subset of observations to be used.
#' @param data optional data frame. When supplied, \code{y}, \code{x}, \code{covs.hte}, \code{covs.eff}, \code{weights}, \code{cluster}, and \code{subset} may be given as bare variable names referring to columns of \code{data}. \code{covs.hte} additionally accepts a one-sided formula (e.g. \code{"~ z1 + z2"}) whose variables are looked up in \code{data} first.
#' @param target.contrast (experimental, in-flight) optional contrast vector or
#'   matrix used to refit the bandwidth to be MSE-optimal for a particular
#'   contrast of the CATE vector. \code{NULL} (default) keeps the standard
#'   per-cell MSE-optimal bandwidth. API may change.
#'
#' @return A list with the following named elements:
#' \item{Estimate}{Vector of conventional local-polynomial RD estimates, one per group level (or per slope-coefficient for continuous \code{covs.hte}). Also available as \code{coef}.}
#' \item{Estimate.bc}{Vector of bias-corrected estimates. Also available as \code{coef.bc}.}
#' \item{se.rb}{Vector of robust bias-corrected standard errors.}
#' \item{ci.rb}{Matrix (\code{n.lev x 2}) of robust bias-corrected confidence-interval bounds.}
#' \item{t.rb}{Vector of asymptotic z-statistics (named \code{t.rb} for legacy reasons; the underlying inference is Gaussian).}
#' \item{pv.rb}{Vector of two-sided p-values from the standard normal.}
#' \item{vcov}{Group-level variance-covariance matrix of \code{Estimate.bc}.}
#' \item{coef.full}{Full coefficient vector from the underlying joint local-polynomial regression (used by \code{\link{rdhte_lincom}}).}
#' \item{vcov.full}{Full variance-covariance matrix of \code{coef.full}.}
#' \item{W.lev}{Group-level identifiers (or coefficient names for continuous \code{covs.hte}).}
#' \item{W.names}{Display labels for the rows of \code{Estimate}.}
#' \item{kernel}{Kernel type used (e.g. \code{"Triangular"}).}
#' \item{bwselect}{Bandwidth selection procedure used.}
#' \item{vce}{Variance estimator display label (e.g. \code{"CR1"}, \code{"HC3"}).}
#' \item{vce_select}{Canonical lowercase variance-estimator name
#'   (e.g. \code{"cr1"}, \code{"hc3"}).}
#' \item{c}{Cutoff value.}
#' \item{h}{An \code{n.lev x 2} matrix of left/right bandwidths, one row per group.}
#' \item{p}{Order of the polynomial used for estimation.}
#' \item{q}{Order of the polynomial used for bias correction.}
#' \item{N}{Length-2 vector \code{c(N_left, N_right)} of pre-bandwidth sample sizes.}
#' \item{Nh}{An \code{n.lev x 2} matrix of effective sample sizes (per group, per side).}
#' \item{covs.cont}{Logical; \code{TRUE} for continuous \code{covs.hte} (or no \code{covs.hte}), \code{FALSE} for factor.}
#' \item{level}{Confidence level used.}
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
#'
#' \dontrun{
#' # Empirical examples using the bundled Granzier, Pons, and Tricaud data.
#' data(rdhte_dataset)
#' with(rdhte_dataset, {
#'   rd_left <- rdhte(y = y, x = x, covs.hte = factor(w_left),
#'                    cluster = cluster_var)
#'   summary(rd_left)
#'   rdhte_lincom(rd_left,
#'                linfct = "`factor(w_left)1` - `factor(w_left)0` = 0")
#'
#'   summary(rdhte(y = y, x = x, covs.hte = factor(w_left),
#'                 cluster = cluster_var, bw.joint = TRUE))
#'   summary(rdhte(y = y, x = x, covs.hte = factor(w_left):factor(w_strong),
#'                 cluster = cluster_var))
#'   summary(rdhte(y = y, x = x, covs.hte = w_strength,
#'                 kernel = "uni", cluster = cluster_var))
#' })
#' }
#' @export
rdhte <- function(y, x, c = 0, covs.hte = NULL, covs.eff = NULL,
                  p = 1, q = NULL, kernel = "tri", weights = NULL,
                  h = NULL, h.l = NULL, h.r = NULL,
                  vce = "hc3", cluster = NULL, level = 95,
                  bwselect = "mserd", bw.joint = FALSE, subset = NULL,
                  data = NULL, target.contrast = NULL) {

  ## Resolve bare-name args against `data`, when supplied. Mirrors the
  ## rdrobust 4.0.0 NSE pattern: y = vote, x = margin, covs.hte = factor(class),
  ## covs.eff = z1, cluster = state, weights = w, subset = state != "Alaska"
  ## all become column references to `data`. Strings (formula syntax,
  ## column-name vectors) pass through to the existing covs.hte parser.
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

  ## Input validation. covs.hte may be a formula/string; skip length check then.
  if (!is.numeric(c) || length(c) != 1L || !is.finite(c))
    stop(sprintf("Cutoff 'c' must be a single finite numeric value (received: %s).", toString(c)), call. = FALSE)
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 100)
    stop("`level` must be a single number in (0, 100).", call. = FALSE)
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
  ## covs.hte: only length-check when supplied as a vector/factor (not formula/string)
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
  if (!is.null(h)   && (!is.numeric(h)   || any(!is.finite(h))   || any(h   <= 0) || length(h)   > 2))
    stop("'h' must be a positive scalar or a length-2 positive vector.", call. = FALSE)
  if (!is.null(h.l) && (!is.numeric(h.l) || any(!is.finite(h.l)) || any(h.l <= 0)))
    stop("'h.l' must be a positive numeric value.", call. = FALSE)
  if (!is.null(h.r) && (!is.numeric(h.r) || any(!is.finite(h.r)) || any(h.r <= 0)))
    stop("'h.r' must be a positive numeric value.", call. = FALSE)
  if (!is.null(target.contrast)) {
    if (is.null(covs.hte))
      stop("'target.contrast' requires 'covs.hte' (needs at least 2 Estimate rows).", call. = FALSE)
    if (!is.numeric(target.contrast) || any(!is.finite(target.contrast)))
      stop("'target.contrast' must be a finite numeric vector.", call. = FALSE)
    if (all(target.contrast == 0))
      stop("'target.contrast' cannot be the zero vector.", call. = FALSE)
    if (!is.null(h) || !is.null(h.l) || !is.null(h.r))
      stop("'target.contrast' selects the bandwidth internally; do not pass 'h', 'h.l', or 'h.r' together with it.", call. = FALSE)
  }

  ## vce normalization + label mapping. Shared with rdbwhte() via
  ## .rdhte_normalize_vce() in helpers.R. Apply the documented
  ## cluster<->no-cluster remappings; canonical names are hc0/.../hc3
  ## and cr1/.../cr3. Internally cr1/cr2/cr3 are forwarded to
  ## sandwich::vcovCL as type = "HC1"/"HC2"/"HC3" (which is exactly
  ## the CR1/CR2/CR3 sandwich when a cluster vector is supplied).
  .user_vce <- !missing(vce)
  .vce_norm <- .rdhte_normalize_vce(vce, cluster, .user_vce)
  vce       <- .vce_norm$vce
  vce_label <- .vce_norm$vce_label

  if (is.null(covs.hte)) {
    #warning("rdhte requires specifying heterogenity variables via covs.hte. RD ATE reported")

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
    ## sandwich::vcovCL only accepts HC0..HC3 as `type`. cr1/cr2/cr3 map
    ## to HC1/HC2/HC3 (which is the cluster-robust sandwich when a
    ## cluster vector is passed via the `cluster` arg).
    vce.hte <- switch(vce,
                      "cr1" = "HC1", "cr2" = "HC2", "cr3" = "HC3",
                      toupper(vce))

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
        rd.bw  <- suppressWarnings(rdbwselect(y = Y, x = Xc, covs = covs, p = p, q = q, bwselect = bwselect, vce = vce, cluster = cluster, kernel = kernel, weights = weights))
        h.lev      <- matrix(c(rd.bw$bws[1], rd.bw$bws[2]), 1, 2)
    } else if (is.null(h.l) & is.null(h.r))  {
        h.lev <- matrix(rep(h, 2), 1, 2)
        bwselect = "Manual"
    } else {
        h.lev <- matrix(c(h.l, h.r), 1, 2)
        bwselect = "Manual"
    }

    N.lev.l <- sum(T==0)
    N.lev.r <- sum(T==1)
    N.lev <- c(N.lev.l, N.lev.r)

    Nh.lev.l <- sum( (abs(Xc) <= h.lev[1,1]) & T==0)
    Nh.lev.r <- sum( (abs(Xc) <= h.lev[1,2]) & T==1)
    Nh.lev <- matrix(c(Nh.lev.l, Nh.lev.r),1,2)


    # Weighted Regression
    r.bw <- (abs(Xc) <= h.lev[1,1] & T==0) |  (abs(Xc) <= h.lev[1,2] & T==1)
    k.weights = NULL
    Kernel = "Uniform"
    if (kernel=="tri") {
      k.weights = (1-abs(Xc)/h.lev[1,1])*(T==0) + (1-abs(Xc)/h.lev[1,2])*(T==1)
      Kernel = "Triangular"
    }
    if (kernel=="epa") {
      k.weights = (1-(Xc/h.lev[1,1])^2)*(T==0) + (1-(Xc/h.lev[1,2])^2)*(T==1)
      Kernel = "Epanechnikov"
    }
    if (is.null(covs.eff)) {
      reg_formula_est <- Y ~ T * Xp
      reg_formula_inf <- Y ~ T * Xq
    } else {
      reg_formula_est <- Y ~ T * Xp + covs
      reg_formula_inf <- Y ~ T * Xq + covs
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

    ## deparse(substitute(...)) is used for nice display labels in summary
    ## ("factor(Wb)0", "factor(Wb)1", ...). When rdhte() is called via
    ## do.call(rdhte, list(..., covs.hte = factor(Wb))), substitute() returns
    ## the EVALUATED factor object (not the unevaluated symbol), and deparse
    ## of a length-n factor produces an n-line representation. Collapse to
    ## one line and fall back to a generic name when the expression is too
    ## long to be a real variable reference.
    covs.hte_chr <- paste0(deparse(substitute(covs.hte)), collapse = " ")
    if (nchar(covs.hte_chr) > 100) covs.hte_chr <- "covs.hte"

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
      ## When `data=` is supplied, look up formula vars in the data frame
      ## first; otherwise the formula's environment (= parent.frame above)
      ## handles lookup as before.
      if (!is.null(data)) {
        mf <- stats::model.frame(f, data = data, na.action = stats::na.pass)
        covs.hte.mm <- stats::model.matrix(f, data = mf)
      } else {
        covs.hte.mm <- stats::model.matrix(f)
      }
      ## Prefix matrix column names so the downstream `grep("^covs.hte", ...)`
      ## finds them. data.frame(covs.hte = <matrix>) uses the matrix's own
      ## column names verbatim, dropping the `covs.hte` prefix; we restore
      ## it manually before binding.
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
  ## sandwich::vcovCL only accepts HC0..HC3 as `type`. cr1/cr2/cr3 map
  ## to HC1/HC2/HC3 (which is the cluster-robust sandwich when a cluster
  ## vector is passed via the `cluster` arg).
  vce.hte <- switch(vce,
                    "cr1" = "HC1", "cr2" = "HC2", "cr3" = "HC3",
                    toupper(vce))

  # Factor or continuous covariate handling for W
  covs.cont <- TRUE
  W.lev <- NULL
  n.lev <- 1
  # N.lev is the (left, right) full-sample split, regardless of whether
  # W is categorical or continuous (matches rdrobust's $N convention).
  N.lev <- c(sum(T == 0), sum(T == 1))

  if (is.factor(W.covs)) {
    W.lev <- levels(W.covs)
    n.lev <- nlevels(W.covs)
    covs.cont <- FALSE
  } else {
    if (length(idx) == 1) {
      cond <- mean((W.covs == 0) | (W.covs == 1))
        if (cond == 1 ) {
          W.covs <- factor(W.covs)
          W.lev <- levels(W.covs)
          n.lev <- nlevels(W.covs)
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
      covs.cont <- FALSE
    }
  }
}




  # Initialize bandwidth storage
  h.vec  <- rep(0, N)
  h.lev  <- matrix(0, n.lev, 2)
  Nh.lev <- matrix(NA, n.lev, 2)

  ## Bandwidth selection. Each branch populates h.vec (per-observation
  ## kernel bandwidth) and h.lev (per-group L/R bandwidths). Nh.lev is
  ## then filled uniformly by .rdhte_fill_Nh(). The continuous /
  ## no-covs.hte case corresponds to n.lev == 1 || is.null(W.lev).
  if (!is.null(h)) {
    bwselect <- "Manual"
    if (n.lev > 1) {
      if (length(h) == 1) {
        h.vec <- matrix(h, N, 1)
        h.lev <- matrix(h, n.lev, 2)
      } else {
        if (abs(n.lev - length(h)) > 0)
          stop("check the number of bandwidths provided")
        for (l in seq_len(n.lev)) {
          h.lev[l, 1] <- h[l]
          h.lev[l, 2] <- h[l]
          h.vec[W.covs == W.lev[l]] <- h[l]
        }
      }
    } else {
      h.vec <- matrix(h, N, 1)
      h.lev <- matrix(h, n.lev, 2)
    }

  } else if (!is.null(h.l) & !is.null(h.r)) {
    bwselect <- "Manual"
    if (n.lev > 1) {
      if (length(h.l) == 1 & length(h.r) == 1) {
        h.vec[Xc <  0] <- h.l
        h.vec[Xc >= 0] <- h.r
        h.lev <- cbind(matrix(h.l, n.lev, 1), matrix(h.r, n.lev, 1))
      } else {
        if ((abs(n.lev - length(h.l)) > 0 | abs(n.lev - length(h.r)) > 0) &
            (length(h.l) > 1 | length(h.r) > 1))
          stop("check the number of bandwidths provided")
        for (l in seq_len(n.lev)) {
          h.lev[l, 1] <- h.l[l]
          h.lev[l, 2] <- h.r[l]
          h.vec[W.covs == W.lev[l] & T == 0] <- h.l[l]
          h.vec[W.covs == W.lev[l] & T == 1] <- h.r[l]
        }
      }
    } else {
      h.vec[Xc <  0] <- h.l
      h.vec[Xc >= 0] <- h.r
      h.lev[1, 1] <- h.l
      h.lev[1, 2] <- h.r
    }
  }

  # Bandwidth not provided
  if (is.null(h) & is.null(h.l) & is.null(h.r)) {
    if (isTRUE(bw.joint) | is.null(W.lev)) {
      # Joint bandwidth estimation
      rd.bw <- suppressWarnings(rdbwselect(y = Y, x = Xc, covs = covs, p = p, q = q,
                                            bwselect = bwselect, vce = vce,
                                            cluster = cluster, kernel = kernel,
                                            weights = weights))
      h.lev[, 1]    <- rep(rd.bw$bws[1], n.lev)
      h.lev[, 2]    <- rep(rd.bw$bws[2], n.lev)
      h.vec[Xc <  0] <- rd.bw$bws[1]
      h.vec[Xc >= 0] <- rd.bw$bws[2]
    } else {
      ## Perf note: per-cell rdbwselect is O(n.lev) calls into rdrobust;
      ## with cluster-robust V-fit each call is meaningfully more
      ## expensive than the joint case. When the cost looks high (many
      ## cells AND cluster vce), drop a one-time note so users know
      ## `bw.joint = TRUE` is available as a faster alternative.
      if (n.lev >= 4 && !is.null(cluster)) {
        message(sprintf(
          "rdhte: per-cell bandwidth selection across %d cells with cluster-robust vce; pass `bw.joint = TRUE` if a shared bandwidth is acceptable (typically ~%dx faster).",
          n.lev, n.lev))
      }
      ## Precompute the cell-subset masks ONCE; reusing `W.covs == W.lev[l]`
      ## inside the loop would recompute it (and again inside rdbwselect).
      cell_masks <- lapply(W.lev, function(lv) W.covs == lv)
      for (l in seq_len(n.lev)) {
        rd.bw <- suppressWarnings(rdbwselect(y = Y, x = Xc, covs = covs, p = p, q = q,
                                              bwselect = bwselect, vce = vce,
                                              cluster = cluster, kernel = kernel,
                                              weights = weights,
                                              subset = cell_masks[[l]]))
        h.lev[l, 1] <- rd.bw$bws[1]
        h.lev[l, 2] <- rd.bw$bws[2]
        m <- cell_masks[[l]]
        h.vec[m & T == 0] <- rd.bw$bws[1]
        h.vec[m & T == 1] <- rd.bw$bws[2]
      }
    }
  }

  ## Effective-N counting: one helper, six former branches.
  Nh.lev <- .rdhte_fill_Nh(Xc, T, h.lev, W.covs = W.covs, W.lev = W.lev)
  if (is.null(W.lev)) W.covs <- as.matrix(W.covs)


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

    ## ---- target.contrast: refit at h*_iota -------------------------------
    ## When the user supplies a target contrast, treat the fit above as a
    ## "pilot": extract (V_hat, B_hat) from it, compute the contrast-targeted
    ## MSE-optimal bandwidth, and refit both regressions at that h. The
    ## downstream per-level extraction then runs on the refit, so the
    ## reported Estimate / Estimate.bc / vcov reflect h*_iota.
    target.iota <- NULL
    h.pilot     <- NULL
    h.star      <- NULL
    V.hat       <- NULL
    B.hat       <- NULL
    if (!is.null(target.contrast)) {
      cn_p <- names(coef(rd_est))
      cn_q <- names(coef(rd_inf))
      varsigma_names <- c("T", cn_p[grepl("^T:W\\.covs", cn_p)])
      Xq_tag <- paste0("Xq", q)
      bias_names <- c(paste0("T:", Xq_tag),
                      cn_q[grepl(paste0("^T:", Xq_tag, ":W\\.covs"), cn_q)])
      n_rows <- length(varsigma_names)
      if (length(target.contrast) != n_rows)
        stop(sprintf("'target.contrast' must have length equal to the number of Estimate rows (%d).",
                     n_rows), call. = FALSE)
      if (length(bias_names) != n_rows)
        stop(sprintf(paste0("rdhte target.contrast: coefficient-name lookup failed ",
                            "(varsigma=%d, bias=%d, expected=%d). ",
                            "Refusing to silently mis-select the contrast."),
                     length(varsigma_names), length(bias_names), n_rows), call. = FALSE)
      ## Q maps Estimate rows -> varsigma. Factor: standard contrast pattern.
      ## Continuous: Estimate rows ARE varsigma entries, so Q is identity.
      Q <- if (isTRUE(covs.cont)) diag(n_rows) else .rdhte_Q_factor(n_rows)
      target.iota <- as.numeric(crossprod(Q, target.contrast))

      k0      <- .rdhte_kernel_k0(kernel)
      h.pilot <- mean(c(h.lev[1, 1], h.lev[1, 2]))
      n_full  <- sum(N.lev)
      VB      <- .rdhte_extract_VB_from_fits(rd_est, rd_inf, vce.hte,
                                             cluster_sub    = cluster[r.bw],
                                             varsigma_names = varsigma_names,
                                             bias_names     = bias_names,
                                             n              = n_full,
                                             h_pilot        = h.pilot,
                                             k0             = k0)
      V.hat <- VB$V
      B.hat <- VB$B
      hs    <- .rdhte_hstar_iota(V.hat, B.hat, target.iota, n_full)
      h.star <- hs$h_star

      ## Rebuild the same fit at h.star (joint, symmetric bandwidth)
      h.lev <- matrix(h.star, n.lev, 2)
      h.vec <- rep(h.star, length(Xc))
      r.bw  <- abs(Xc) <= h.star
      k.weights <- switch(kernel,
                          "tri" = pmax(1 - abs(Xc) / h.star, 0),
                          "epa" = pmax(1 - (Xc / h.star)^2, 0),
                          "uni" = as.numeric(abs(Xc) <= h.star),
                          as.numeric(r.bw))
      if (!is.null(weights)) k.weights <- k.weights * weights

      rd_est       <- lm(reg_formula_est, weights = k.weights, subset = r.bw)
      rd_inf       <- lm(reg_formula_inf, weights = k.weights, subset = r.bw)
      rd.vcov.full <- vcovCL(rd_inf, cluster = cluster[r.bw], type = vce.hte)

      if (n.lev > 1) {
        for (l in 1:n.lev) {
          ind <- W.covs == W.lev[l]
          Nh.lev[l, 1] <- sum((abs(Xc[ind]) <= h.star) & (T[ind] == 0))
          Nh.lev[l, 2] <- sum((abs(Xc[ind]) <= h.star) & (T[ind] == 1))
        }
      } else {
        Nh.lev <- matrix(c(sum(abs(Xc) <= h.star & T == 0),
                           sum(abs(Xc) <= h.star & T == 1)), 1, 2)
      }
      bwselect <- "msecontrast"
    }
    ## ----------------------------------------------------------------------

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
      if (isTRUE(covs.cont)) {
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
    if (!is.null(covs.hte) & isTRUE(covs.cont)) rdmodel <- "Sharp RD Heterogeneous Treatment Effects: Continuous."
    if (!is.null(covs.eff)) rdmodel <-paste(rdmodel, ", Covariate-Adjusted", sep="")
    if (!is.null(cluster))  rdmodel <-paste(rdmodel, ", Cluster-Adjusted", sep="")

  
    #if (cls == "formula") {
      #W.names = c("T", paste("T:",W.names, sep=""))
    #} else {
    #  W.names <- paste(covs.hte_chr,W.lev, sep="")  
    #}
    


  ## Slots populated only when target.contrast was used (NULL otherwise).
  ## target.contrast: the user-facing input (Estimate-row weights).
  ## target.iota:     the varsigma-space selector (Q' target.contrast).
  ## h.pilot, h.star: the pilot bandwidth and the contrast-targeted h.
  ## V, B:            the vector-first asymptotic objects from the pilot.
  ## Symbols only exist in the heterogeneous branch; default NULL so the
  ## ATE branch returns the same shape.
  if (!exists("target.iota", inherits = FALSE)) target.iota <- NULL
  if (!exists("h.pilot",     inherits = FALSE)) h.pilot     <- NULL
  if (!exists("h.star",      inherits = FALSE)) h.star      <- NULL
  if (!exists("V.hat",       inherits = FALSE)) V.hat       <- NULL
  if (!exists("B.hat",       inherits = FALSE)) B.hat       <- NULL

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
               W.names     = W.names,
               covs.hte_chr = covs.hte_chr,
               kernel      = Kernel,
               bwselect    = bwselect,
               vce         = vce_label,
               vce_select  = vce,
               c = c,
               h = h.lev,
               p = p,
               q = q,
               N = N.lev,
               Nh = Nh.lev,
               covs.cont = covs.cont,
               level = level,
               rdmodel = rdmodel,
               target.contrast = target.contrast,
               target.iota     = target.iota,
               h.pilot         = h.pilot,
               h.star          = h.star,
               V               = V.hat,
               B               = B.hat)
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

  invisible(x)
}



################################################################################
#' Summary method for \code{rdhte} objects.
#'
#' @param object Class \code{rdhte} object returned by \code{\link{rdhte}}.
#' @param ... ignored.
#'
#' @keywords internal
#' @return An S3 object inheriting from \code{rdhte} with class
#'   \code{summary.rdhte}; printed by \code{print.summary.rdhte}.
#' @export
summary.rdhte <- function(object, ...) {
  out <- object
  out$.summary <- TRUE
  class(out) <- c("summary.rdhte", class(object))
  out
}

#' Internal function.
#'
#' @param x Class \code{summary.rdhte} object returned by \code{summary.rdhte}.
#' @param ... ignored.
#'
#' @keywords internal
#' @return The input \code{x} returned invisibly; called for side effects.
#' @export
print.summary.rdhte <- function(x, ...) {
  cat(paste(x$rdmodel,"\n", sep=""))
  cat(paste("","\n", sep=""))

  name.long <- max(12, nchar(x$covs.hte_chr))
  llenght   <- name.long + 12 + 12 + 12 + 24 + 10 + 10 + 10 + 10

  cat(paste("Number of Obs.           ",  format(sum(x$N),   width=10, justify="right"),"\n", sep=""))
  cat(paste("BW type                  ",  format(x$bwselect, width=10, justify="right"),"\n", sep=""))
  cat(paste("Kernel                   ",  format(x$kernel,   width=10, justify="right"),"\n", sep=""))
  cat(paste("VCE method               ",  format(x$vce,      width=10, justify="right"),"\n", sep=""))
  cat("\n")

  cat(paste("Number of Obs.           ",  format(x$N[1], width=10, justify="right"),  "   ", format(x$N[2], width=10, justify="right"), "\n", sep=""))
  if (isTRUE(x$covs.cont)) {
    cat(paste("Eff. Number of Obs.      ",  format(x$Nh[1], width=10, justify="right"),  "   ", format(x$Nh[2], width=10, justify="right"), "\n", sep=""))
    llenght <- name.long + 12 + 12 + 12 + 24
  }
  cat(paste("Order est. (p)           ",  format(x$p, width=10, justify="right"),  "   ", format(x$p, width=10, justify="right"), "\n", sep=""))
  cat(paste("Order bias  (q)          ",  format(x$q, width=10, justify="right"),  "   ", format(x$q, width=10, justify="right"), "\n", sep=""))
  if (isTRUE(x$covs.cont)) {
    cat(paste("BW est. (h)              ",  format(sprintf("%10.3f", x$h[1])),  "   ", format(sprintf("%10.3f", x$h[2])), "\n", sep=""))
  }

  cat("\n")

  cat(paste(rep("=", llenght), collapse="")); cat("\n")
  cat(format(""                 , width=name.long, justify="right"))
  cat(format("Point"            , width=12, justify="right"))
  cat(format("Robust Inference" , width=24, justify="right"))
  cat("\n")

  if (isFALSE(x$covs.cont)) {
    cat(format(x$covs.hte_chr   , width=name.long, justify="right"))
    cat(format("Estimate"       , width=12, justify="right"))
  } else {
    cat(format(""               , width=name.long, justify="right"))
    cat(format("Estimate"       , width=12, justify="right"))
  }

  cat(format("z"        , width=12, justify="right"))
  cat(format("Pr(>|z|)" , width=12, justify="right"))
  cat(format(paste("[ ", x$level, "% ", "C.I. ]", sep=""), width=24, justify="centre"))

  if (isFALSE(x$covs.cont)) {
    cat(format("Nh-", width=10, justify="right"))
    cat(format("Nh+", width=10, justify="right"))
    cat(format("h-" , width=10, justify="right"))
    cat(format("h+" , width=10, justify="right"))
  }
  cat("\n")

  cat(paste(rep("-", llenght), collapse="")); cat("\n")

  for (i in 1:length(x$W.lev)) {
    if (isTRUE(x$covs.cont)) cat(format(x$W.names[i], width=name.long, justify="right"))
    else                     cat(format(x$W.lev[i],   width=name.long, justify="right"))
    cat(format(sprintf("%3.3f", x$Estimate[i]),                           width=12, justify="right"))
    cat(format(sprintf("%3.3f", x$t.rb[i]),                               width=12, justify="right"))
    cat(format(sprintf("%3.3f", x$pv.rb[i]),                              width=12, justify="right"))
    cat(format(paste("[", sprintf("%3.3f", x$ci.rb[i,1]), " , ", sep=""), width=12, justify="right"))
    cat(format(paste(     sprintf("%3.3f", x$ci.rb[i,2]), "]",   sep=""), width=12, justify="left"))

    if (isFALSE(x$covs.cont)) {
      cat(format(x$Nh[i,1],                 width=10, justify="right"))
      cat(format(x$Nh[i,2],                 width=10, justify="right"))
      cat(format(sprintf("%3.3f",x$h[i,1]), width=10, justify="right"))
      cat(format(sprintf("%3.3f",x$h[i,2]), width=10, justify="right"))
    }
    cat("\n")
  }

  cat(paste(rep("=", llenght), collapse="")); cat("\n")

  invisible(x)
}

#' Tidy a \code{rdhte} object
#'
#' broom-compatible \code{tidy()} method. Returns one row per heterogeneity
#' subgroup (factor case) or per series term (continuous case).
#'
#' @param x A \code{rdhte} object.
#' @param ... ignored.
#' @return A data frame with columns \code{term}, \code{estimate} (conventional),
#'   \code{std.error} (robust SE), \code{statistic}, \code{p.value},
#'   \code{conf.low}, \code{conf.high} (robust BC CI), \code{estimate.bc},
#'   \code{h.left}, \code{h.right}, \code{n.eff.left}, \code{n.eff.right}.
#' @keywords internal
#' @exportS3Method broom::tidy
tidy.rdhte <- function(x, ...) {
  is_cont <- isTRUE(x$covs.cont)
  terms   <- if (is_cont) as.character(x$W.names) else as.character(x$W.lev)
  data.frame(
    term        = terms,
    estimate    = as.numeric(x$Estimate),
    std.error   = as.numeric(x$se.rb),
    statistic   = as.numeric(x$t.rb),
    p.value     = as.numeric(x$pv.rb),
    conf.low    = as.numeric(x$ci.rb[, 1]),
    conf.high   = as.numeric(x$ci.rb[, 2]),
    estimate.bc = as.numeric(x$Estimate.bc),
    h.left      = as.numeric(x$h[, 1]),
    h.right     = as.numeric(x$h[, 2]),
    n.eff.left  = as.numeric(x$Nh[, 1]),
    n.eff.right = as.numeric(x$Nh[, 2]),
    row.names   = NULL,
    stringsAsFactors = FALSE
  )
}

#' Glance at a \code{rdhte} object
#'
#' broom-compatible \code{glance()} method. Returns a one-row summary of the
#' fit (sample sizes, polynomial orders, kernel, VCE, BW selector).
#'
#' @param x A \code{rdhte} object.
#' @param ... ignored.
#' @return A one-row data frame.
#' @keywords internal
#' @exportS3Method broom::glance
glance.rdhte <- function(x, ...) {
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
    level    = x$level,
    n.terms  = length(x$W.lev),
    covs.continuous = isTRUE(x$covs.cont),
    stringsAsFactors = FALSE
  )
}
