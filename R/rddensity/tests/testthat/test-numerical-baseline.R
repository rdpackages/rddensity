baseline_data <- function() {
  index <- seq(0, 199)
  seq(-1.5, 1.5, length.out = 200) + 0.05 * sin(index)
}

masspoint_data <- function() {
  index <- seq(0, 239)
  round(seq(-1.4, 1.6, length.out = 240) + 0.08 * sin(index / 2), 1)
}

append_series <- function(rows, case_name, prefix, values) {
  values <- unlist(values)
  rbind(
    rows,
    data.frame(
      case = case_name,
      metric = paste(prefix, names(values), sep = "_"),
      value = as.numeric(values),
      stringsAsFactors = FALSE
    )
  )
}

summarize_rddensity <- function(case_name, fit) {
  rows <- data.frame(case = character(), metric = character(), value = numeric())
  for (prefix in c("hat", "sd_asy", "sd_jk", "test", "h", "N", "hat_p", "sd_asy_p", "sd_jk_p", "test_p")) {
    metric_prefix <- if (prefix == "N") "n" else prefix
    rows <- append_series(rows, case_name, metric_prefix, fit[[prefix]])
  }
  rows
}

summarize_rdbwdensity <- function(case_name, fit) {
  h <- as.data.frame(fit$h)
  rows <- data.frame(case = character(), metric = character(), value = numeric())
  for (row_name in row.names(h)) {
    for (col_name in colnames(h)) {
      rows <- rbind(
        rows,
        data.frame(
          case = case_name,
          metric = paste("h", row_name, col_name, sep = "_"),
          value = as.numeric(h[row_name, col_name]),
          stringsAsFactors = FALSE
        )
      )
    }
  }
  rows <- rbind(
    rows,
    data.frame(
      case = case_name,
      metric = paste("n", names(fit$N), sep = "_"),
      value = as.numeric(unlist(fit$N)),
      stringsAsFactors = FALSE
    )
  )
  rows
}

compute_baseline <- function() {
  x <- baseline_data()
  x_masspoints <- masspoint_data()
  rows <- list(
    summarize_rddensity("rddensity_fixed", rddensity(x, h = c(0.6, 0.6), bino = FALSE)),
    summarize_rddensity("rddensity_estimated", rddensity(x, bino = FALSE)),
    summarize_rddensity("rddensity_uniform", rddensity(x, h = c(0.7, 0.7), kernel = "uniform", bino = FALSE)),
    summarize_rddensity("rddensity_restricted", rddensity(x, h = c(0.8, 0.8), fitselect = "restricted", bino = FALSE)),
    summarize_rddensity("rddensity_all", rddensity(x, h = c(0.6, 0.6), all = TRUE, bino = FALSE)),
    summarize_rddensity("rddensity_epanechnikov_plugin", rddensity(x, h = c(0.75, 0.65), kernel = "epanechnikov", vce = "plugin", bino = FALSE)),
    summarize_rddensity("rddensity_nonzero_cutoff", rddensity(x, c = 0.15, h = c(0.5, 0.7), bino = FALSE)),
    summarize_rddensity("rddensity_masspoints_adjusted", rddensity(x_masspoints, h = c(0.8, 0.8), bino = FALSE)),
    summarize_rddensity("rddensity_masspoints_unadjusted", rddensity(x_masspoints, h = c(0.8, 0.8), massPoints = FALSE, bino = FALSE)),
    summarize_rdbwdensity("rdbwdensity_default", rdbwdensity(x)),
    summarize_rdbwdensity("rdbwdensity_uniform_plugin", rdbwdensity(x, kernel = "uniform", vce = "plugin")),
    summarize_rdbwdensity("rdbwdensity_nonzero_cutoff", rdbwdensity(x, c = 0.15, p = 3, kernel = "epanechnikov")),
    summarize_rdbwdensity("rdbwdensity_masspoints_adjusted", rdbwdensity(x_masspoints)),
    summarize_rdbwdensity("rdbwdensity_masspoints_unadjusted", rdbwdensity(x_masspoints, massPoints = FALSE))
  )
  do.call(rbind, rows)
}

if (identical(Sys.getenv("RDDENSITY_WRITE_BASELINE"), "true")) {
  write.csv(
    compute_baseline(),
    file.path("tests", "testthat", "fixtures", "numerical-baseline.csv"),
    row.names = FALSE,
    quote = FALSE,
    na = "NaN"
  )
} else {
  test_that("numerical baseline matches current contract", {
    expected <- read.csv(test_path("fixtures", "numerical-baseline.csv"), check.names = FALSE)
    observed <- compute_baseline()

    expect_equal(observed[c("case", "metric")], expected[c("case", "metric")], ignore_attr = TRUE)
    expect_equal(observed$value, expected$value, tolerance = 1e-10, ignore_attr = TRUE)
  })
}
