\name{rdhte_contrast}
\alias{rdhte_contrast}
\title{RD Heterogeneous Treatment Effects: contrast-targeted bandwidth}
\usage{
rdhte_contrast(model, contrast, refit = TRUE, level = 95, tol = 1e-10)
}
\arguments{
\item{model}{A fitted model returned by \code{\link{rdhte}}.}

\item{contrast}{Numeric vector with length equal to \code{nrow(model$Estimate)}.
The user-facing linear combination of Estimate rows whose MSE is optimized.}

\item{refit}{Logical. \code{TRUE} refits the model at the contrast-targeted
bandwidth. \code{FALSE} projects the existing fit and returns inference under
the model's current bandwidth.}

\item{level}{Confidence level in percentage form, in \code{[1, 100)}.}

\item{tol}{Tolerance below which the estimated contrast bias is treated as
zero and the contrast-targeted bandwidth is declared undefined.}
}
\value{
A list with the conventional estimate, robust bias-corrected estimate,
robust standard error, z-statistic, p-value, confidence interval, selected
bandwidth, pilot bandwidth, contrast diagnostics, and (when
\code{refit = TRUE}) the refitted \code{rdhte} object.
}
\description{
\code{rdhte_contrast} refits an \code{\link{rdhte}} model at the
MSE-optimal bandwidth for a single user-chosen linear combination of the
Estimate vector and returns the resulting point estimate, robust
bias-corrected estimate, standard error, p-value, and confidence interval.
}
\details{
Use \code{\link{rdhte_lincom}} instead when you want joint Wald inference on
several linear restrictions at one bandwidth.
}
\examples{
\dontrun{
set.seed(123)
n <- 1000
X <- runif(n, -1, 1)
W <- rbinom(n, 1, 0.5)
Y <- 3 + 2 * X + 1.5 * X^2 + 3 * W * (X >= 0) + rnorm(n)
m <- rdhte(y = Y, x = X, covs.hte = factor(W))
rdhte_contrast(m, contrast = c(-1, 1))
}
}
\seealso{\code{\link{rdhte}}, \code{\link{rdhte_lincom}}}
