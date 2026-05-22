################################################################################
#' @title rdhte: RD Heterogeneous Treatment Effects Estimation and Inference
#'
#' @description Building on the recent developments in Calonico, Cattaneo, Farrell, Palomba, and Titiunik (2025), this package implements
#'  estimation and inference of heterogeneous treatment effects in RD designs.
#'  The package includes two main commands: \code{\link{rdhte}} conduct estimation and robust bias-corrected inference for conditional RD treatment effects,
#'  for a given choice of bandwidth parameter; and \code{\link{rdbwhte}} implements automatic bandwidth selection methods.
#'  We illustrate the methods implemented in the package \code{\link{rdhte}} using a canonical empirical application.
#'  We also demonstrate how the package \code{\link{rdhte}} complements, and in very specific cases recovers,
#'  the methods available in the packages rdrobust (Calonico, Cattaneo, Farrell, Titiunik (2017) and rdmulti, Cattaneo, Titiunik, VazquezBare (2020).
#'
#' Commands: \code{\link{rdhte}} for estimation and inference.
#'   \code{\link{rdbwhte}} for data-driven bandwidth selection.
#'
#' Related software packages useful for inference in regression discontinuity
#'   (RD) designs are described in the website:
#'   \url{https://rdpackages.github.io/}.
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
#'
#' @importFrom stats lm
#' @importFrom stats na.omit
#' @importFrom stats poly
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#' @importFrom stats model.matrix
#' @importFrom stats pchisq
#' @importFrom stats pf
#' @importFrom stats pt
#' @importFrom stats qt
#' @importFrom stats setNames
#' @importFrom stats as.formula 
#' @importFrom stats coef
#' @importFrom stats confint
#' @import rdrobust
#' @import sandwich
#' @import multcomp
#'
#' @aliases rdhte-package
"_PACKAGE"

## Silence R CMD check NOTE about the ggplot2 .data pronoun used in plot.rdhte().
utils::globalVariables(".data")
