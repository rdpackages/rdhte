#' Plot rdhte heterogeneous treatment effects
#'
#' For an \code{rdhte} object fit with a categorical \code{covs.hte},
#' draw one point per group at the conventional point estimate
#' (\code{Estimate}) with the robust bias-corrected confidence interval
#' (\code{ci.rb}). A dashed horizontal line at zero gives a visual
#' reference for the null effect.
#'
#' Continuous \code{covs.hte} (or no \code{covs.hte}) is not yet
#' supported -- the function errors with a clear message in those cases.
#'
#' @param x An object of class \code{rdhte} returned by \code{\link{rdhte}}.
#' @param sort Logical or character; if \code{TRUE} (or \code{"effect"}),
#'   reorder groups along the x-axis by point estimate. \code{FALSE}
#'   (default) keeps the original \code{W.lev} order.
#' @param point.size Numeric; size of the point markers. Default 2.5.
#' @param errorbar.width Numeric; width of the error-bar caps relative
#'   to the x-axis discrete unit. Default 0.2.
#' @param zero.line Logical; if \code{TRUE} (default), draw a dashed
#'   horizontal line at \code{y = 0}.
#' @param title,xlab,ylab Optional plot annotations. Defaults derived
#'   from the rdhte object (\code{x$rdmodel}, \code{x$covs.hte_chr},
#'   "Treatment effect").
#' @param ... Currently unused.
#'
#' @return Invisibly, a \code{ggplot} object.
#'
#' @details Requires the \pkg{ggplot2} package. The intervals shown are
#' the same robust bias-corrected CIs reported by \code{print(x)} and
#' \code{summary(x)}: they are centered on \code{Estimate.bc} (not
#' \code{Estimate}) and use the robust standard error \code{se.rb}.
#' Because the point and the CI center can differ slightly, the point
#' may sit just inside or just outside the bar; this is the rdrobust
#' convention and is not a plotting bug.
#'
#' @examples
#' \dontrun{
#'   set.seed(1)
#'   n <- 600
#'   x <- runif(n, -1, 1)
#'   W <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
#'   y <- 0.5 + (x >= 0) * (W == "A") * 1.0 +
#'              (x >= 0) * (W == "B") * 2.0 +
#'              (x >= 0) * (W == "C") * 0.3 + rnorm(n)
#'   m <- rdhte(y = y, x = x, covs.hte = W)
#'   plot(m)
#'   plot(m, sort = TRUE)        # reorder by effect size
#' }
#'
#' @export
plot.rdhte <- function(x, sort = FALSE,
                       point.size = 2.5,
                       errorbar.width = 0.2,
                       zero.line = TRUE,
                       title = NULL, xlab = NULL, ylab = NULL,
                       ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot.rdhte() requires the 'ggplot2' package. ",
         "Install it with install.packages(\"ggplot2\").",
         call. = FALSE)
  }

  if (isTRUE(x$covs.cont)) {
    stop("plot.rdhte() currently supports only categorical (factor) ",
         "covs.hte. The fitted model has continuous (or no) covs.hte.",
         call. = FALSE)
  }

  group_lab <- if (!is.null(x$W.lev)) as.character(x$W.lev) else as.character(x$W.names)
  if (is.null(group_lab) || length(group_lab) == 0L) {
    stop("plot.rdhte(): could not extract group labels from the rdhte object.",
         call. = FALSE)
  }

  est   <- as.numeric(x$Estimate)
  cilo  <- as.numeric(x$ci.rb[, 1])
  cihi  <- as.numeric(x$ci.rb[, 2])
  if (length(est) != length(group_lab) ||
      length(cilo) != length(group_lab) ||
      length(cihi) != length(group_lab)) {
    stop("plot.rdhte(): length of Estimate / ci.rb does not match number of groups.",
         call. = FALSE)
  }

  df <- data.frame(group = group_lab,
                   estimate = est,
                   ci_low = cilo,
                   ci_high = cihi,
                   stringsAsFactors = FALSE)

  if (isTRUE(sort) || identical(sort, "effect")) {
    df <- df[order(df$estimate), , drop = FALSE]
  }
  ## Preserve order on the x-axis using a factor with explicit levels.
  df$group <- factor(df$group, levels = df$group)

  if (is.null(title)) title <- x$rdmodel %||% "rdhte: heterogeneous treatment effects"
  if (is.null(xlab))  xlab  <- if (!is.null(x$covs.hte_chr) && nzchar(x$covs.hte_chr)) x$covs.hte_chr else "Group"
  if (is.null(ylab))  ylab  <- "Treatment effect"

  rotate_x <- max(nchar(as.character(df$group))) > 6 || nrow(df) > 6

  ## Build the ggplot object piecewise so each layer is testable.
  ##
  ## We deliberately reference ggplot2 functions via `ggplot2::` to keep
  ## ggplot2 a Suggests dependency.
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$group, y = .data$estimate))
  if (isTRUE(zero.line)) {
    p <- p + ggplot2::geom_hline(yintercept = 0,
                                 linetype = "dashed",
                                 colour = "grey40")
  }
  p <- p +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$ci_low,
                                        ymax = .data$ci_high),
                           width = errorbar.width) +
    ggplot2::geom_point(size = point.size) +
    ggplot2::labs(title = title, x = xlab, y = ylab) +
    ggplot2::theme_bw()
  if (isTRUE(rotate_x)) {
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                                hjust = 1))
  }

  print(p)
  invisible(p)
}

## Local null-default helper to keep the file self-contained without
## adding an rlang dependency.
`%||%` <- function(a, b) if (is.null(a)) b else a
