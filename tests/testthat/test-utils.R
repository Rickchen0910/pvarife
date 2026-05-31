test_that("lag_lead_matrix: lag-1 produces correct NA structure", {
  x <- matrix(1:8, nrow = 4, ncol = 2)
  res <- lag_lead_matrix(x, 1L, 1L)

  expect_equal(nrow(res), 4L)
  expect_equal(ncol(res), 2L)
  # First row must be NA (no lag-1 value at t=1)
  expect_true(all(is.na(res[1L, ])))
  # Rows 2..4 should equal x rows 1..3
  expect_equal(res[2:4, ], x[1:3, ], ignore_attr = TRUE)
})

test_that("lag_lead_matrix: lag-2 offsets correctly", {
  x <- matrix(seq_len(10), nrow = 5, ncol = 2)
  res <- lag_lead_matrix(x, 2L, 2L)
  expect_equal(ncol(res), 2L)
  expect_true(all(is.na(res[1:2, ])))
  expect_equal(res[3:5, ], x[1:3, ], ignore_attr = TRUE)
})

test_that("lag_lead_matrix: lead-1 produces correct NA structure", {
  x <- matrix(1:8, nrow = 4, ncol = 2)
  res <- lag_lead_matrix(x, -1L, -1L)
  # Last row must be NA (no lead value at t=T)
  expect_true(all(is.na(res[4L, ])))
  # Rows 1..3 should equal x rows 2..4
  expect_equal(res[1:3, ], x[2:4, ], ignore_attr = TRUE)
})

test_that("lag_lead_matrix: multiple lags give correct column count", {
  x <- matrix(1:12, nrow = 6, ncol = 2)
  res <- lag_lead_matrix(x, 1L, 3L)
  expect_equal(ncol(res), 6L)  # 2 cols * 3 lags
  expect_equal(nrow(res), 6L)
})

test_that("lag_lead_matrix: zero lag is identity", {
  x <- matrix(1:8, nrow = 4)
  res <- lag_lead_matrix(x, 0L, 0L)
  expect_equal(res, x, ignore_attr = TRUE)
})

test_that("ols_nan: clean data gives same result as lm()", {
  set.seed(1L)
  n <- 50L
  x <- cbind(1, matrix(stats::rnorm(n * 2L), nrow = n))
  b_true <- c(1.0, 2.0, -1.0)
  y <- matrix(x %*% b_true + stats::rnorm(n) * 0.5, ncol = 1L)

  res <- ols_nan(y, x)
  lm_b <- unname(stats::coef(stats::lm(y ~ x[, -1L])))

  expect_equal(as.numeric(res$beta), lm_b, tolerance = 1e-10)
})

test_that("ols_nan: rows with NA are correctly dropped", {
  set.seed(2L)
  n <- 20L
  x <- cbind(1, stats::rnorm(n))
  y <- matrix(x %*% c(1, 2) + stats::rnorm(n) * 0.3, ncol = 1L)

  # Introduce NAs
  y[c(3L, 7L, 15L)] <- NA
  res <- ols_nan(y, x)

  # Reference: lm on complete cases
  keep  <- !is.na(y)
  lm_b  <- unname(stats::coef(stats::lm(y[keep] ~ x[keep, -1L])))
  expect_equal(as.numeric(res$beta), lm_b, tolerance = 1e-10)
  expect_equal(nrow(res$y_clean), sum(keep))
})
