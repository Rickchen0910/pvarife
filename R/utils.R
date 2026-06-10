#' Build a matrix of lags and/or leads
#'
#' Creates a matrix whose columns are lagged or lead copies of the columns of
#' \code{x}, for integer offsets \code{from_lag} through \code{to_lag}. A
#' positive offset \eqn{k} produces a lag (rows \eqn{k+1} to \eqn{T} of the
#' output receive rows \eqn{1} to \eqn{T-k} of \code{x}); a negative offset
#' \eqn{k} produces a lead. Unfilled positions are set to \code{NA}.
#'
#' Faithful translation of \code{lagleadmatrix.m} from the replication package
#' of Tugan (2021).
#'
#' @param x A numeric matrix with \eqn{T} rows and \eqn{K} columns.
#' @param from_lag Integer. The smallest offset to include (may be negative for
#'   leads).
#' @param to_lag Integer. The largest offset to include (must be \eqn{\ge}
#'   \code{from_lag}).
#'
#' @return A matrix with \eqn{T} rows and \eqn{K \times (to\_lag - from\_lag +
#'   1)} columns. Column block \eqn{j} (1-indexed) corresponds to offset
#'   \eqn{from\_lag + j - 1}.
#'
#' @examples
#' x <- matrix(1:12, nrow = 4)
#' lag_lead_matrix(x, 1, 2)   # two lags
#' lag_lead_matrix(x, -1, 1)  # one lead, contemporaneous, one lag
#'
#' @references
#' Tugan, M. (2021). Panel VAR models with interactive fixed effects.
#' \emph{Econometrics Journal}, 24, 225--246.
#' \doi{10.1093/ectj/utaa021}
#'
#' @export
lag_lead_matrix <- function(x, from_lag, to_lag) {
  x <- as.matrix(x)
  nr <- nrow(x)
  nc <- ncol(x)
  lag_length <- to_lag - from_lag + 1L
  out <- matrix(NA_real_, nrow = nr, ncol = nc * lag_length)

  for (kk in seq(from_lag, to_lag)) {
    col_start <- (kk - from_lag) * nc + 1L
    col_end   <- (kk - from_lag + 1L) * nc

    if (kk >= 0L) {
      # lag: rows (kk+1):nr in output get rows 1:(nr-kk) of x
      if (kk < nr) {
        out[(kk + 1L):nr, col_start:col_end] <- x[1L:(nr - kk), , drop = FALSE]
      }
    } else {
      # lead: rows 1:(nr+kk) in output get rows (-kk+1):nr of x
      if (-kk < nr) {
        out[1L:(nr + kk), col_start:col_end] <- x[(-kk + 1L):nr, , drop = FALSE]
      }
    }
  }
  out
}


#' OLS with missing observations
#'
#' Computes ordinary least squares after dropping any row of
#' \code{cbind(y, x)} that contains at least one \code{NA}. Equivalent to
#' \code{ols_NaNincluded.m} in the Tugan (2021) replication package.
#'
#' @param y Numeric matrix or vector of dependent variables (\eqn{T \times M}).
#' @param x Numeric matrix of regressors (\eqn{T \times K}).
#'
#' @return A list with components:
#'   \describe{
#'     \item{beta}{OLS coefficient matrix (\eqn{K \times M}).}
#'     \item{residuals}{Residuals on the complete-case subsample.}
#'     \item{residuals_full}{Residuals on the full sample (NA where dropped).}
#'     \item{std_err}{Standard errors (only when \eqn{M = 1}).}
#'     \item{y_clean}{Dependent variable after dropping NA rows.}
#'     \item{x_clean}{Regressors after dropping NA rows.}
#'   }
#'
#' @noRd
ols_nan <- function(y, x) {
  y <- as.matrix(y)
  x <- as.matrix(x)
  data <- cbind(y, x)
  keep <- stats::complete.cases(data)
  y_clean <- y[keep, , drop = FALSE]
  x_clean <- x[keep, , drop = FALSE]
  beta <- solve(crossprod(x_clean), crossprod(x_clean, y_clean))
  resid <- y_clean - x_clean %*% beta
  resid_full <- y - x %*% beta

  result <- list(
    beta         = beta,
    residuals    = resid,
    residuals_full = resid_full,
    y_clean      = y_clean,
    x_clean      = x_clean
  )

  if (ncol(y) == 1L) {
    nn <- nrow(y_clean)
    pp <- nrow(beta)
    result$std_err <- sqrt(diag(solve(crossprod(x_clean)) *
                                  as.numeric(crossprod(resid) / (nn - pp))))
  }
  result
}


# ---------------------------------------------------------------------------
# Internal helpers for stacking / reshaping panel arrays
# ---------------------------------------------------------------------------

# Build Y_c and Z_c from the y[i,t,k] array.
# Returns a list: y_stack (NT*n_units x 1), z_stack (NT*n_units x K+K^2*n_lags),
#   y_c array (NT x 1 x n_units), z_c array (NT x K+K^2*n_lags x n_units),
#   and n_obs_c (n_units) with valid observation counts.
#' @noRd
.build_yz <- function(y_arr, n_lags) {
  n_units <- dim(y_arr)[1L]
  n_time  <- dim(y_arr)[2L]
  n_vars  <- dim(y_arr)[3L]

  n_rows <- n_vars * n_time          # NT rows per unit (before trimming kmax)
  n_cols_z <- n_vars + n_vars^2 * n_lags

  y_c <- array(NA_real_, dim = c(n_rows, 1L, n_units))
  z_c <- array(NA_real_, dim = c(n_rows, n_cols_z, n_units))

  for (ii in seq_len(n_units)) {
    # y_arr[ii, , ] is (n_time x n_vars)
    ymat <- y_arr[ii, , , drop = FALSE]
    ymat <- matrix(ymat, nrow = n_time, ncol = n_vars)

    # Period-level NA propagation, faithful to EmpiricalApplication.m
    # (NonStackedDVC(any(isnan(...),2),:,i) = NaN): if ANY variable is missing
    # at period t, the whole period is treated as missing. This also rules out
    # partially-observed periods downstream.
    incomplete_t <- !stats::complete.cases(ymat)
    if (any(incomplete_t)) ymat[incomplete_t, ] <- NA_real_

    lag_mat <- lag_lead_matrix(ymat, 1L, n_lags)  # T x K*n_lags

    yy_c <- matrix(NA_real_, nrow = n_vars * n_time, ncol = 1L)
    zz_c <- matrix(NA_real_, nrow = n_vars * n_time, ncol = n_cols_z)

    for (tt in seq_len(n_time)) {
      row_start <- n_vars * (tt - 1L) + 1L
      row_end   <- n_vars * tt

      yy_c[row_start:row_end, 1L] <- ymat[tt, ]

      zz_c[row_start:row_end, ] <- cbind(
        diag(n_vars),
        kronecker(diag(n_vars), matrix(lag_mat[tt, ], nrow = 1L))
      )
    }

    # Stack and remove the first n_lags*n_vars rows (no complete lags yet)
    dat <- cbind(yy_c, zz_c)
    dat[seq_len(n_lags * n_vars), ] <- NA_real_
    dat[apply(is.na(dat), 1L, any), ] <- NA_real_

    y_c[, 1L, ii] <- dat[, 1L]
    z_c[, , ii]   <- dat[, -1L]
  }

  # Build pooled stacks (dropping NA rows)
  y_stack <- do.call(rbind, lapply(seq_len(n_units), function(ii) {
    yy <- y_c[, 1L, ii]           # plain vector, length n_rows
    matrix(yy[!is.na(yy)], ncol = 1L)
  }))
  z_stack <- do.call(rbind, lapply(seq_len(n_units), function(ii) {
    zz   <- matrix(z_c[, , ii], nrow = n_rows, ncol = n_cols_z)
    keep <- !apply(is.na(zz), 1L, any)
    zz[keep, , drop = FALSE]
  }))

  list(
    y_c     = y_c,
    z_c     = z_c,
    y_stack = y_stack,
    z_stack = z_stack,
    n_units = n_units,
    n_time  = n_time,
    n_vars  = n_vars,
    n_lags  = n_lags,
    n_rows  = n_rows,
    n_cols_z = n_cols_z
  )
}
