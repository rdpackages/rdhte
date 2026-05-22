\name{rdhte_dataset}
\alias{rdhte_dataset}
\docType{data}
\title{Granzier-Pons-Tricaud Election Dataset for rdhte Illustrations}
\description{
A working extract of the data used in Granzier, Pons, and Tricaud
(2023, AEJ: Applied), distributed with the \pkg{rdhte} package as a
fixed-format example dataset for the package vignette and worked
examples. The running variable \code{x} is the candidate's
first-round margin against the threshold for advancing to the
runoff; \code{y} is an indicator for advancing.
}
\usage{data("rdhte_dataset")}
\format{
A data frame with 39,534 observations on the following variables.
\describe{
  \item{\code{y}}{Outcome: indicator (0/1) for advancing to the runoff.}
  \item{\code{x}}{Running variable: first-round margin against the
                  qualifying threshold (the cutoff is at zero).}
  \item{\code{cluster_var}}{Cluster identifier (district).}
  \item{\code{w_left}}{Binary: 1 if the candidate's party is left of
                       center, 0 otherwise.}
  \item{\code{w_ideology}}{Unordered categorical (4 levels): party
                            ideology bucket.}
  \item{\code{w_strength}}{Continuous: candidate's average prior
                            national-level vote share, used as a
                            proxy for ex-ante candidate strength.}
  \item{\code{w_strong}}{Binary: above-median strength
                         (\code{w_strength}).}
  \item{\code{w_strength_qrt}}{Ordered categorical (4 levels):
                                quartile bucket of \code{w_strength}.}
}
}
\source{
  Granzier, R., V. Pons, and C. Tricaud (2023). Coordination and
  Bandwagon Effects: How Past Rankings Shape the Behavior of Voters
  and Candidates. \emph{American Economic Journal: Applied Economics}
  15(4): 177-217.
}
\examples{
data(rdhte_dataset)
str(rdhte_dataset)
}
\keyword{datasets}
