#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

script_path <- tryCatch(normalizePath(sys.frame(1)$ofile, mustWork = TRUE), error = function(e) NA_character_)
if (!is.na(script_path)) {
  repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
} else {
  repo_root <- normalizePath(".", mustWork = TRUE)
}

output <- if (length(args) >= 1) args[[1]] else file.path(repo_root, "docs", "audit", "baselines", "r-current.json")
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(rdrobust)
  library(sandwich)
  library(multcomp)
  library(jsonlite)
})

source(file.path(repo_root, "R", "rdhte", "R", "rdbwhte.R"))
source(file.path(repo_root, "R", "rdhte", "R", "rdhte.R"))
source(file.path(repo_root, "R", "rdhte", "R", "rdhte_lincom.R"))

clean_names <- function(x) {
  x <- ifelse(is.na(x) | x == "", "unnamed", x)
  gsub("[.]+", "_", make.names(x, unique = TRUE))
}

num <- function(x) {
  if (length(x) == 0 || is.null(x) || is.na(x) || !is.finite(x)) return(NULL)
  unname(as.numeric(x))
}

named_numbers <- function(x, names_hint = NULL) {
  vals <- as.numeric(x)
  nms <- names(x)
  if (is.null(nms) && !is.null(names_hint) && length(names_hint) == length(vals)) nms <- names_hint
  if (is.null(nms)) nms <- paste0("v", seq_along(vals))
  stats::setNames(lapply(vals, num), clean_names(nms))
}

matrix_numbers <- function(x, row_names = NULL, col_names = NULL) {
  m <- as.matrix(x)
  rn <- rownames(m)
  cn <- colnames(m)
  if (is.null(rn) && !is.null(row_names) && length(row_names) == nrow(m)) rn <- row_names
  if (is.null(cn) && !is.null(col_names) && length(col_names) == ncol(m)) cn <- col_names
  if (is.null(rn)) rn <- paste0("r", seq_len(nrow(m)))
  if (is.null(cn)) cn <- if (ncol(m) == 2) c("left", "right") else paste0("c", seq_len(ncol(m)))

  out <- list()
  for (i in seq_len(nrow(m))) {
    for (j in seq_len(ncol(m))) {
      out[[paste(clean_names(rn[i]), clean_names(cn[j]), sep = "_")]] <- num(m[i, j])
    }
  }
  out
}

numeric_matrix_or_null <- function(x, row_names = NULL, col_names = NULL) {
  if (!(is.numeric(x) || is.matrix(x) || is.data.frame(x))) return(NULL)
  matrix_numbers(x, row_names, col_names)
}

summarize_rdhte <- function(obj) {
  coef_names <- names(obj$Estimate)
  h_rows <- if (length(obj$W.lev) == nrow(as.matrix(obj$h))) obj$W.lev else NULL
  list(
    estimate = named_numbers(obj$Estimate),
    estimate_bc = named_numbers(obj$Estimate.bc, coef_names),
    se_rb = named_numbers(obj$se.rb, coef_names),
    t_rb = named_numbers(obj$t.rb, coef_names),
    pv_rb = named_numbers(obj$pv.rb, coef_names),
    ci_rb = matrix_numbers(obj$ci.rb, coef_names, c("lower", "upper")),
    h = matrix_numbers(obj$h, h_rows, c("left", "right")),
    n = named_numbers(obj$N),
    nh = matrix_numbers(obj$Nh, h_rows, c("left", "right")),
    vcov = numeric_matrix_or_null(obj$vcov, coef_names, coef_names),
    kernel = obj$kernel,
    vce = obj$vce,
    bwselect = obj$bwselect,
    rdmodel = obj$rdmodel
  )
}

summarize_rdbwhte <- function(obj) {
  h_rows <- if (length(obj$W.lev) == nrow(as.matrix(obj$h))) obj$W.lev else NULL
  list(
    h = matrix_numbers(obj$h, h_rows, c("left", "right")),
    n = named_numbers(obj$N),
    nh = matrix_numbers(obj$Nh, h_rows, c("left", "right")),
    kernel = obj$kernel,
    vce = obj$vce,
    bwselect = obj$bwselect,
    rdmodel = obj$rdmodel
  )
}

summarize_lincom <- function(obj) {
  individual <- obj$individual
  joint <- obj$joint
  list(
    individual = lapply(seq_len(nrow(individual)), function(i) {
      list(
        hypothesis = individual$hypothesis[i],
        estimate = num(individual$estimate[i]),
        t_stat = num(individual$t_stat[i]),
        p_value = num(individual$p_value[i]),
        conf_low = num(individual$conf.low[i]),
        conf_high = num(individual$conf.high[i])
      )
    }),
    joint = list(
      statistic_chi2 = num(joint$statistic.chi2[1]),
      df = num(joint$df[1]),
      p_value = num(joint$p_value[1])
    )
  )
}

data <- read.csv(file.path(repo_root, "R", "rdhte_dataset.csv"), stringsAsFactors = TRUE)
attach(data)
on.exit(detach(data), add = TRUE)

rd_left <- rdhte(y = y, x = x, covs.hte = factor(w_left), cluster = cluster_var)
rd_left_joint <- rdhte(y = y, x = x, covs.hte = w_left, cluster = cluster_var, bw.joint = TRUE)
bw_left <- rdbwhte(y = y, x = x, covs.hte = factor(w_left), cluster = cluster_var)
rd_ideology <- rdhte(y = y, x = x, covs.hte = factor(w_ideology), cluster = cluster_var)
rd_strength <- rdhte(y = y, x = x, covs.hte = w_strength, kernel = "uni", cluster = cluster_var)
rd_interaction <- rdhte(y = y, x = x, covs.hte = "w_left*w_strength", h = 0.1, cluster = cluster_var)
rd_average <- rdhte(y = y, x = x, h = 0.1, vce = "hc3")

cases <- list(
  binary_left = list(
    rdhte = summarize_rdhte(rd_left),
    lincom = summarize_lincom(
      rdhte_lincom(
        model = rd_left,
        linfct = c("`factor(w_left)1` - `factor(w_left)0` = 0"),
        digits = 12
      )
    )
  ),
  binary_left_joint = list(rdhte = summarize_rdhte(rd_left_joint)),
  binary_left_bw = list(rdbwhte = summarize_rdbwhte(bw_left)),
  categorical_ideology = list(rdhte = summarize_rdhte(rd_ideology)),
  continuous_strength = list(rdhte = summarize_rdhte(rd_strength)),
  interaction_strength = list(rdhte = summarize_rdhte(rd_interaction)),
  average_manual = list(rdhte = summarize_rdhte(rd_average))
)

result <- list(
  schema_version = 1,
  package = "rdhte",
  language = "r",
  source = "working-tree",
  timestamp_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  environment = list(
    r = R.version.string,
    platform = R.version$platform
  ),
  cases = cases
)

jsonlite::write_json(result, output, pretty = TRUE, auto_unbox = TRUE, digits = 16, null = "null")
cat("Wrote", output, "\n")
