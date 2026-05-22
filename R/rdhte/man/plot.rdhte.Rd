% Generated manually for plot.rdhte
\name{plot.rdhte}
\alias{plot.rdhte}
\title{Plot rdhte heterogeneous treatment effects}
\usage{
\method{plot}{rdhte}(x, sort = FALSE, point.size = 2.5,
                     errorbar.width = 0.2, zero.line = TRUE,
                     title = NULL, xlab = NULL, ylab = NULL, ...)
}
\arguments{
\item{x}{An object of class \code{rdhte} returned by \code{\link{rdhte}}.}

\item{sort}{Logical or character; if \code{TRUE} (or \code{"effect"}),
reorder groups along the x-axis by point estimate. \code{FALSE}
(default) keeps the original \code{W.lev} order.}

\item{point.size}{Numeric; size of the point markers. Default 2.5.}

\item{errorbar.width}{Numeric; width of the error-bar caps relative
to the x-axis discrete unit. Default 0.2.}

\item{zero.line}{Logical; if \code{TRUE} (default), draw a dashed
horizontal line at \code{y = 0}.}

\item{title, xlab, ylab}{Optional plot annotations. Defaults derived
from the rdhte object (\code{x$rdmodel}, \code{x$covs.hte_chr},
"Treatment effect").}

\item{...}{Currently unused.}
}
\value{
Invisibly, a \code{ggplot} object.
}
\description{
For an \code{rdhte} object fit with a categorical \code{covs.hte},
draw one point per group at the conventional point estimate
(\code{Estimate}) with the robust bias-corrected confidence interval
(\code{ci.rb}). A dashed horizontal line at zero gives a visual
reference for the null effect.
}
\details{
Continuous \code{covs.hte} (or no \code{covs.hte}) is not yet
supported -- the function errors with a clear message in those cases.

Requires the \pkg{ggplot2} package. The intervals shown are the same
robust bias-corrected CIs reported by \code{print(x)} and
\code{summary(x)}: they are centered on \code{Estimate.bc} (not
\code{Estimate}) and use the robust standard error \code{se.rb}.
Because the point and the CI center can differ slightly, the point
may sit just inside or just outside the bar; this is the rdrobust
convention and is not a plotting bug.
}
\examples{
\dontrun{
  set.seed(1)
  n <- 600
  x <- runif(n, -1, 1)
  W <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
  y <- 0.5 + (x >= 0) * (W == "A") * 1.0 +
             (x >= 0) * (W == "B") * 2.0 +
             (x >= 0) * (W == "C") * 0.3 + rnorm(n)
  m <- rdhte(y = y, x = x, covs.hte = W)
  plot(m)
  plot(m, sort = TRUE)        # reorder by effect size
}
}
