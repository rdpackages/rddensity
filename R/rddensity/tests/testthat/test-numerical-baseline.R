baseline_data <- function() {
  index <- seq(0, 199)
  seq(-1.5, 1.5, length.out = 200) + 0.05 * sin(index)
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
  rows <- list(
    summarize_rddensity("rddensity_fixed", rddensity(x, h = c(0.6, 0.6), bino = FALSE)),
    summarize_rddensity("rddensity_uniform", rddensity(x, h = c(0.7, 0.7), kernel = "uniform", bino = FALSE)),
    summarize_rddensity("rddensity_restricted", rddensity(x, h = c(0.8, 0.8), fitselect = "restricted", bino = FALSE)),
    summarize_rddensity("rddensity_all", rddensity(x, h = c(0.6, 0.6), all = TRUE, bino = FALSE)),
    summarize_rdbwdensity("rdbwdensity_default", rdbwdensity(x)),
    summarize_rdbwdensity("rdbwdensity_uniform_plugin", rdbwdensity(x, kernel = "uniform", vce = "plugin"))
  )
  do.call(rbind, rows)
}

test_that("numerical baseline matches current contract", {
  expected <- read.csv(test_path("fixtures", "numerical-baseline.csv"), check.names = FALSE)
  observed <- compute_baseline()

  expect_equal(observed[c("case", "metric")], expected[c("case", "metric")], ignore_attr = TRUE)
  expect_equal(observed$value, expected$value, tolerance = 1e-10, ignore_attr = TRUE)
})
