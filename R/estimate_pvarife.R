#' Estimate a Panel VAR with Interactive Fixed Effects
#'
#' Jointly estimates VAR coefficients \eqn{\beta}, latent common factors
#' \eqn{F}, and factor loadings \eqn{\Lambda} for a panel vector autoregression
#' with interactive fixed effects, following the iterative algorithm of
#' Tugan (2021) based on Bai (2009).
#'
#' @details
#' The model is
#' \deqn{y_{i,t} = \sum_{l=1}^{\ell} \Theta_l y_{i,t-l} + F_t \lambda_i + e_{i,t},}
#' where \eqn{y_{i,t}} is a \eqn{K \times 1} vector of endogenous variables for
#' unit \eqn{i} at time \eqn{t}, \eqn{F_t} is an \eqn{r \times 1} vector of
#' unobservable common factors, \eqn{\lambda_i} is a unit-specific loading
#' vector, and \eqn{e_{i,t}} is an idiosyncratic error.
#'
#' The algorithm alternates between:
#' \enumerate{
#'   \item An \strong{inner loop} that extracts factors and loadings via PCA
#'     (principal components on the residual cross-product matrix) and imputes
#'     missing observations (EM step of Bai 2009).
#'   \item An \strong{outer loop} that updates \eqn{\beta} via least squares
#'     after projecting out the estimated factors (using \eqn{M_F = I - F(F'F)^{-1}F'}).
#' }
#'
#' @param y A numeric array of dimension \eqn{I \times T \times K}
#'   (units \eqn{\times} time \eqn{\times} variables). \code{NA} values are
#'   allowed for unbalanced panels.
#' @param n_lags Positive integer. Lag order \eqn{\ell}.
#' @param n_factors Positive integer. Number of interactive fixed effects \eqn{r}.
#' @param n_out Positive integer. Number of outer iterations (default 50).
#'   Corresponds to \code{out_number} in the MATLAB replication code.
#' @param n_in Positive integer. Number of inner PCA/EM iterations per outer
#'   step (default 10). Corresponds to \code{in_number} in the MATLAB code.
#' @param balanced_init Logical. If \code{TRUE} (default), the initial beta
#'   estimate is obtained from units that have at least 10 fully observed
#'   periods in the last window of the sample — matching the approach of the
#'   MATLAB \code{Initial_Step_in_Iteration.m} for unbalanced real data. Set
#'   to \code{FALSE} for balanced panels (e.g., Monte Carlo simulations) to
#'   skip this selection step and use all units directly.
#'
#' @return An object of class \code{"pvarife_result"}, which is a list with:
#'   \describe{
#'     \item{beta}{Coefficient vector of length \eqn{K + K^2 \ell}.
#'       The first \eqn{K} elements are equation-specific intercepts; the
#'       remaining \eqn{K^2 \ell} elements are VAR lag coefficients stacked
#'       as \eqn{[\Theta_1, \Theta_2, \ldots, \Theta_\ell]'} (column-major).}
#'     \item{ff}{Factor matrix of dimension \eqn{T \times r} (one row per
#'       time period; analogous to MATLAB \code{f}).}
#'     \item{factors_mat}{Block-diagonal factor matrix of dimension
#'       \eqn{TK \times Kr}, built as \eqn{I_K \otimes f_t'}.}
#'     \item{loadings}{Matrix of dimension \eqn{Kr \times I} (factor loadings
#'       per unit, stacked by variable).}
#'     \item{sigma}{Reduced-form error covariance matrix (\eqn{K \times K}).}
#'     \item{u_c}{Array of residuals \eqn{TK \times 1 \times I} (NA at
#'       unobserved positions).}
#'     \item{y_arr}{The original input array \eqn{I \times T \times K} (used,
#'       e.g., for initial conditions in \code{\link{bootstrap_irf_bands}}).}
#'     \item{y_c, z_c}{Stacked outcome/regressor arrays (\eqn{TK \times 1 \times I}
#'       and \eqn{TK \times (K + K^2\ell) \times I}).}
#'     \item{y_stack, z_stack}{Pooled outcome/regressor matrices (complete-case
#'       rows only).}
#'     \item{i_obs}{Integer matrix \eqn{TK \times I}: 1 = observed, 0 = missing.}
#'     \item{n_time_i}{Integer vector of length \eqn{I}: number of observed
#'       time periods per unit (analogous to MATLAB \code{TC}).}
#'     \item{tnc_i}{Integer vector of length \eqn{I}: \eqn{n\_time\_i \times K}
#'       (analogous to MATLAB \code{TNC}).}
#'     \item{n_lags, n_factors, n_vars, n_units, n_time}{Dimensions.}
#'   }
#'
#' @references
#' Tugan, M. (2021). Panel VAR models with interactive fixed effects.
#' \emph{Econometrics Journal}, 24, 225--246.
#' \doi{10.1093/ectj/utaa021}
#'
#' Bai, J. (2009). Panel data models with interactive fixed effects.
#' \emph{Econometrica}, 77(4), 1229--1279.
#'
#' @examples
#' sim <- sim_pvarife(n_units = 30, n_time = 20, n_vars = 2,
#'                    n_lags = 1, n_factors = 1, seed = 1)
#' fit <- pvarife(sim$y, n_lags = 1, n_factors = 1, n_out = 5, n_in = 3)
#' print(fit)
#'
#' @export
pvarife <- function(y, n_lags, n_factors, n_out = 50L, n_in = 10L,
                    balanced_init = TRUE) {
  y <- as.array(y)
  stopifnot(length(dim(y)) == 3L)
  n_out <- as.integer(n_out)
  n_in  <- as.integer(n_in)

  # Build stacked arrays
  yz <- .build_yz(y, n_lags)

  # Initial estimate
  init <- .initial_step(yz$y_stack, yz$z_stack, yz$y_c, yz$z_c, y, n_factors,
                        balanced_init = balanced_init)

  # Run inner/outer iteration
  res <- .inner_outer_iter(
    y_c      = yz$y_c,
    z_c      = yz$z_c,
    y_stack  = yz$y_stack,
    z_stack  = yz$z_stack,
    y_arr    = y,
    i_obs    = init$i_obs,
    beta_in  = init$beta_initial,
    w_c_in   = init$w_c,
    n_out    = n_out,
    n_in     = n_in,
    n_factors = n_factors
  )

  structure(
    c(res,
      list(
        y_arr   = y,
        y_c     = yz$y_c,
        z_c     = yz$z_c,
        y_stack = yz$y_stack,
        z_stack = yz$z_stack,
        i_obs   = init$i_obs,
        n_lags  = n_lags,
        n_factors = n_factors,
        n_vars  = yz$n_vars,
        n_units = yz$n_units,
        n_time  = yz$n_time
      )
    ),
    class = "pvarife_result"
  )
}


# Inner/outer EM iteration — faithful translation of Inner_Outer_Iteration.m
# from Tugan (2021).
#
# @noRd
.inner_outer_iter <- function(y_c, z_c, y_stack, z_stack, y_arr,
                               i_obs, beta_in, w_c_in,
                               n_out, n_in, n_factors) {

  n_units  <- dim(y_arr)[1L]
  n_time   <- dim(y_arr)[2L]
  n_vars   <- dim(y_arr)[3L]
  n_rows   <- dim(y_c)[1L]   # = n_vars * n_time

  beta_iter <- beta_in
  w_c_iter  <- w_c_in

  for (i_outer in seq_len(n_out)) {

    # -----------------------------------------------------------------------
    # Inner loop: update factors F and loadings lambda given current beta
    # -----------------------------------------------------------------------
    for (i_inner in seq_len(n_in)) {

      # Accumulate cross-product matrix: sum_i v_i v_i'
      vvprime <- matrix(0.0, nrow = n_time, ncol = n_time)

      for (ii in seq_len(n_units)) {
        w_vec <- w_c_iter[, 1L, ii]
        # Reshape (n_vars*n_time) vector -> T x K matrix
        # MATLAB: reshape(W_c', K, T)' => each column of W_c' is one variable
        v_mat <- matrix(w_vec, nrow = n_vars, ncol = n_time)
        v_mat <- t(v_mat)   # T x K
        vvprime <- vvprime + tcrossprod(v_mat)
      }

      # PCA: top n_factors eigenvectors, scaled by sqrt(T)
      scale_val <- as.numeric(n_time * n_vars * n_units)
      ev <- eigen(vvprime / scale_val, symmetric = TRUE)
      # eigen() returns descending order for symmetric matrices
      ff <- sqrt(n_time) * ev$vectors[, seq_len(n_factors), drop = FALSE]
      # ff is T x r

      # Build block-diagonal factors_mat: TK x Kr
      # Each T-block: I_K o f_t' (MATLAB: kron(eye(K), f(t,:)))
      factors_mat <- .build_factors_mat(ff, n_vars)

      # Update loadings: lambda = (I_I o (F'F)^{-1} F') W_stack
      # where W_stack = vec(W_c_1, ..., W_c_I)
      w_stack_full <- do.call(c, lapply(seq_len(n_units),
                                        function(ii) w_c_iter[, 1L, ii]))
      w_stack_full <- matrix(w_stack_full, ncol = 1L)

      # MATLAB: LAMDBA = kron(eye(C), inv(F'F)*F') * W_iteration
      #         lambda = reshape(LAMDBA', size(F,2), C)
      ff_inv <- solve(crossprod(factors_mat), t(factors_mat))  # (Kr x TK)
      # Split into per-unit blocks
      loadings <- matrix(NA_real_, nrow = ncol(factors_mat), ncol = n_units)
      for (ii in seq_len(n_units)) {
        row_start <- (ii - 1L) * n_rows + 1L
        row_end   <- ii * n_rows
        loadings[, ii] <- ff_inv %*% w_stack_full[row_start:row_end, , drop = FALSE]
      }
      # Note: MATLAB reshapes LAMDBA' so rows = r*K, cols = n_units
      # loadings here is (Kr x n_units): column ii = lambda_i

      # EM imputation: fill missing W_c entries with F * lambda_i
      # Observed entries: actual residual with current beta
      w_c_new <- array(0.0, dim = dim(w_c_iter))
      for (ii in seq_len(n_units)) {
        obs  <- i_obs[, ii] == 1L
        miss <- i_obs[, ii] == 0L

        if (any(miss)) {
          w_c_new[miss, 1L, ii] <-
            (factors_mat %*% loadings[, ii])[miss]
        }
        if (any(obs)) {
          yy <- y_c[obs, 1L, ii]
          zz <- z_c[obs, , ii, drop = FALSE]
          zz <- matrix(zz, nrow = sum(obs))
          w_c_new[obs, 1L, ii] <- yy - zz %*% beta_iter
        }
      }
      w_c_iter <- w_c_new
    }  # end inner loop

    # -----------------------------------------------------------------------
    # Outer update: least squares after M_F projection (quasi-differencing)
    # -----------------------------------------------------------------------
    m_f <- diag(n_rows) - factors_mat %*% solve(crossprod(factors_mat), t(factors_mat))

    sum_zprimemfz <- matrix(0.0, nrow = ncol(z_c), ncol = ncol(z_c))
    sum_zprimemfy <- matrix(0.0, nrow = ncol(z_c), ncol = 1L)

    for (ii in seq_len(n_units)) {
      obs_rows <- which(i_obs[, ii] == 1L)
      if (length(obs_rows) == 0L) next

      zz <- matrix(z_c[obs_rows, , ii], nrow = length(obs_rows))
      yy <- matrix(y_c[obs_rows, 1L, ii], ncol = 1L)
      mf_sub <- m_f[obs_rows, obs_rows, drop = FALSE]

      sum_zprimemfz <- sum_zprimemfz + crossprod(zz, mf_sub) %*% zz
      sum_zprimemfy <- sum_zprimemfy + crossprod(zz, mf_sub) %*% yy
    }

    beta_iter <- solve(sum_zprimemfz, sum_zprimemfy)

    # Update W_c with the new beta
    for (ii in seq_len(n_units)) {
      obs  <- i_obs[, ii] == 1L
      miss <- i_obs[, ii] == 0L
      if (any(miss)) {
        w_c_iter[miss, 1L, ii] <-
          (factors_mat %*% loadings[, ii])[miss]
      }
      if (any(obs)) {
        yy <- y_c[obs, 1L, ii]
        zz <- matrix(z_c[obs, , ii], nrow = sum(obs))
        w_c_iter[obs, 1L, ii] <- yy - zz %*% beta_iter
      }
    }

  }  # end outer loop

  # -----------------------------------------------------------------------
  # Compute final residuals u_c and Sigma
  # -----------------------------------------------------------------------
  u_c <- array(NA_real_, dim = dim(y_c))
  n_time_i <- integer(n_units)
  tnc_i    <- integer(n_units)
  sigma_acc <- matrix(0.0, nrow = n_vars, ncol = n_vars)

  for (ii in seq_len(n_units)) {
    obs_rows <- which(i_obs[, ii] == 1L)
    if (length(obs_rows) == 0L) next

    yy <- matrix(y_c[obs_rows, 1L, ii], ncol = 1L)
    zz <- matrix(z_c[obs_rows, , ii], nrow = length(obs_rows))
    ci <- (factors_mat %*% loadings[, ii])[obs_rows, , drop = FALSE]

    u_c[obs_rows, 1L, ii] <- as.numeric(yy - zz %*% beta_iter - ci)

    # Sigma from COMPLETE time periods only, MATLAB-style:
    # reshape the full TK residual vector (NA at unobserved rows) to T x K and
    # drop any period with a missing variable. Reshaping only the observed
    # rows would silently misalign whenever a period is partially observed
    # (observed row count not divisible by K).
    u_mat_full <- t(matrix(u_c[, 1L, ii], nrow = n_vars))   # T x K, NAs kept
    keep_t <- stats::complete.cases(u_mat_full)
    n_time_i[ii] <- sum(keep_t)
    tnc_i[ii]    <- sum(keep_t) * n_vars
    sigma_acc <- sigma_acc + crossprod(u_mat_full[keep_t, , drop = FALSE])
  }

  sigma <- sigma_acc / sum(n_time_i)

  list(
    beta        = beta_iter,
    ff          = ff,
    factors_mat = factors_mat,
    loadings    = loadings,
    sigma       = sigma,
    u_c         = u_c,
    n_time_i    = n_time_i,
    tnc_i       = tnc_i
  )
}


# Build block-diagonal TK x Kr factor matrix from T x r factor matrix ff.
# Block t (rows (t-1)*K+1 : t*K) = kron(I_K, f_t') = diag(K) o f[t,]
#' @noRd
.build_factors_mat <- function(ff, n_vars) {
  n_time_loc <- nrow(ff)
  n_fac <- ncol(ff)
  out <- matrix(0.0, nrow = n_time_loc * n_vars, ncol = n_vars * n_fac)
  for (tt in seq_len(n_time_loc)) {
    row_start <- (tt - 1L) * n_vars + 1L
    row_end   <- tt * n_vars
    out[row_start:row_end, ] <- kronecker(diag(n_vars), matrix(ff[tt, ], nrow = 1L))
  }
  out
}
