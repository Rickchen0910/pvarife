#' Compute impulse response functions under recursive (Cholesky) identification
#'
#' Computes structural impulse responses from an estimated \code{pvarife_result}
#' using a recursive (lower-triangular Cholesky) identification scheme. The
#' identification follows the short-run restriction approach of Tugan (2021).
#'
#' @details
#' The MA representation is computed recursively:
#' \deqn{B_0 = I_K, \quad B_h = \sum_{l=1}^{\ell} \Theta_l B_{h-l},}
#' with the convention \eqn{B_j = 0} for \eqn{j < 0}.
#'
#' \strong{Short-run identification} (default): The impact matrix is the
#' lower-triangular Cholesky factor of \eqn{\hat\Sigma}:
#' \eqn{A_0 = \mathrm{chol}(\hat\Sigma)^\top}.
#'
#' \strong{Long-run identification} (Blanchard-Quah type): The long-run
#' multiplier \eqn{C(1) = (I - \sum_\ell \Theta_\ell)^{-1} A_0} is constrained
#' to be lower-triangular. The impact matrix is
#' \deqn{A_0 = (I - \Theta)\,\mathrm{chol}(D)^\top,}
#' where \eqn{\Theta = \sum_{\ell=1}^{L} \Theta_\ell} and
#' \eqn{D = (I - \Theta)^{-1} \hat\Sigma (I - \Theta)^{-\top}}.
#' This identification restricts shock 1 to have no long-run effect on variable 2
#' (in a 2-variable system). Faithful to \code{IRs_to_Shocks_LR_Identification.m}
#' in the Monte Carlo replication code.
#'
#' \strong{Bias correction:} When \code{bias_correct = TRUE}, the impact matrix
#' is evaluated at the bias-corrected coefficient vector
#' \eqn{\hat\beta - \hat b} from \code{\link{asymptotic_var}}. The uncorrected
#' estimator is used by default (\code{bias_correct = FALSE}) for speed; users
#' who need confidence bands can rely on \code{\link{irf_bands}}, whose
#' median is already implicitly bias-corrected.
#'
#' The impulse response to shock \eqn{s} at horizon \eqn{h} is
#' \eqn{B_h A_0 e_s} where \eqn{e_s} is the \eqn{s}-th standard basis vector.
#'
#' @param fit An object of class \code{"pvarife_result"}.
#' @param n_periods Positive integer. Number of IRF horizons to compute.
#' @param shock Positive integer. Index of the structural shock (1 = first
#'   variable in the ordering). Default is 1.
#' @param diff_vars Integer vector. Indices of variables for which cumulative
#'   IRFs are reported (e.g., for variables entered in first differences).
#'   Default is \code{integer(0)} (no cumulation). For long-run identification
#'   \code{diff_vars = 1L} is common when variable 1 is treated as I(1).
#' @param identification Character. Either \code{"short_run"} (default,
#'   Cholesky of \eqn{\hat\Sigma}) or \code{"long_run"} (Blanchard-Quah:
#'   long-run multiplier lower-triangular). See Details.
#' @param bias_correct Logical. If \code{TRUE}, the IRF is evaluated at the
#'   bias-corrected coefficient vector \eqn{\hat\beta - \hat b} from
#'   \code{\link{asymptotic_var}}. Default \code{FALSE} (faster; use
#'   \code{\link{irf_bands}} for full bias-corrected inference).
#'
#' @return A matrix of dimension \eqn{K \times n\_periods}. Row \eqn{k} gives
#'   the response of variable \eqn{k} to the structural shock at horizons
#'   \eqn{1, \ldots, n\_periods}. The object has class
#'   \code{c("pvarife_irf", "matrix")}.
#'
#' @references
#' Tugan, M. (2021). Panel VAR models with interactive fixed effects.
#' \emph{Econometrics Journal}, 24, 225--246.
#' \doi{10.1093/ectj/utaa021}
#'
#' @examples
#' sim <- sim_pvarife(n_units = 30, n_time = 20, n_vars = 2,
#'                    n_lags = 1, n_factors = 1, seed = 1)
#' fit <- pvarife(sim$y, n_lags = 1, n_factors = 1, n_out = 5, n_in = 3)
#' ir  <- compute_irf(fit, n_periods = 8)
#' dim(ir)   # 2 x 8
#'
#' # Long-run identification
#' ir_lr <- compute_irf(fit, n_periods = 8, identification = "long_run",
#'                      diff_vars = 1L)
#'
#' @seealso \code{\link{irf_bands}}, \code{\link{bootstrap_irf_bands}}
#'
#' @export
compute_irf <- function(fit, n_periods, shock = 1L, diff_vars = integer(0),
                        identification = c("short_run", "long_run"),
                        bias_correct = FALSE) {
  stopifnot(inherits(fit, "pvarife_result"))
  identification <- match.arg(identification)
  n_periods <- as.integer(n_periods)
  shock     <- as.integer(shock)

  n_vars <- fit$n_vars
  n_lags <- fit$n_lags
  sigma  <- fit$sigma

  # Optionally apply bias correction to the coefficient vector
  beta <- as.numeric(fit$beta)
  if (bias_correct) {
    avar <- asymptotic_var(fit)
    beta <- beta - avar$bias
  }

  # Extract VAR coefficient matrices alpha_1, ..., alpha_n_lags
  # First n_vars elements of beta are intercepts (one per equation)
  alpha_vec <- beta[(n_vars + 1L):length(beta)]
  # Reshape: MATLAB alpha=(reshape(beta(N+1:end), N*kmax, N))'
  # so alpha is (K) x (K*n_lags), with alpha[,1:K] = Theta_1, etc.
  alpha <- matrix(alpha_vec, nrow = n_vars * n_lags, ncol = n_vars)
  alpha <- t(alpha)  # K x (K * n_lags)

  alpha_k <- array(0.0, dim = c(n_vars, n_vars, n_lags))
  for (ll in seq_len(n_lags)) {
    col_s <- (ll - 1L) * n_vars + 1L
    col_e <- ll * n_vars
    alpha_k[, , ll] <- alpha[, col_s:col_e, drop = FALSE]
  }

  # MA representation B_0, B_1, ..., B_{n_periods}
  # B_0 = I_K; B_h = sum_{l=1}^{min(h,n_lags)} alpha_l B_{h-l}
  b_arr <- array(0.0, dim = c(n_vars, n_vars, n_periods + 1L))
  b_arr[, , 1L] <- diag(n_vars)   # B_0

  for (hh in seq_len(n_periods)) {
    b_h <- matrix(0.0, n_vars, n_vars)
    for (ll in seq_len(min(hh, n_lags))) {
      b_h <- b_h + alpha_k[, , ll] %*% b_arr[, , hh - ll + 1L]
    }
    b_arr[, , hh + 1L] <- b_h
  }

  # Impact matrix A0 — depends on identification scheme
  if (identification == "short_run") {
    # Short-run: Cholesky of Sigma (lower-triangular)
    a0 <- t(chol(sigma))
  } else {
    # Long-run: Blanchard-Quah type (faithful to IRs_to_Shocks_LR_Identification.m)
    # Theta = sum of all lag matrices
    theta_lr <- Reduce("+", lapply(seq_len(n_lags), function(ll) alpha_k[, , ll]))
    # (I - Theta)^{-1}
    lr_inv <- solve(diag(n_vars) - theta_lr)
    # Long-run variance D = (I-Theta)^{-1} Sigma (I-Theta)^{-T}
    d_mat  <- lr_inv %*% sigma %*% t(lr_inv)
    # A0 = (I - Theta) * chol(D)'  =>  C(1) = lr_inv %*% A0 = chol(D)' (lower-triangular)
    a0 <- (diag(n_vars) - theta_lr) %*% t(chol(d_mat))
  }

  # Structural shock vector
  e_shock <- matrix(0.0, nrow = n_vars, ncol = 1L)
  e_shock[shock, 1L] <- 1.0

  # IRF: K x n_periods  (horizons 1,...,n_periods, analogous to MATLAB h=1:nir)
  ir <- matrix(NA_real_, nrow = n_vars, ncol = n_periods)
  for (hh in seq_len(n_periods)) {
    ir[, hh] <- as.numeric(b_arr[, , hh] %*% a0 %*% e_shock)
  }

  # Cumulate differenced variables
  if (length(diff_vars) > 0L) {
    ir[diff_vars, ] <- t(apply(ir[diff_vars, , drop = FALSE], 1L, cumsum))
  }

  structure(ir, class = c("pvarife_irf", "matrix"))
}


# Internal variant that accepts a pre-computed common component matrix
# (used in irf_bands() for the parametric simulation).
# @param beta_b  coefficient vector (numeric)
# @param y_c     NT x 1 x n_units array
# @param z_c     NT x n_cols_z x n_units array
# @param cc      NT x n_units matrix (common component: F*lambda)
# @param sigma_in K x K covariance matrix to use (from the random draw)
# @param n_vars, n_lags, n_time_i, i_obs — from the fit object
# @noRd
.irf_worker <- function(beta_b, y_c, z_c, cc, n_vars, n_lags,
                         n_time_i, i_obs, n_units, n_periods,
                         shock, diff_vars,
                         identification = "short_run") {
  # Compute residuals and Sigma for this draw (complete periods only;
  # full-vector reshape avoids misalignment under partial missingness)
  u_acc     <- matrix(0.0, nrow = n_vars, ncol = n_vars)
  sum_tc_b  <- 0L
  n_rows    <- dim(y_c)[1L]

  for (ii in seq_len(n_units)) {
    obs <- which(i_obs[, ii] == 1L)
    if (length(obs) == 0L) next
    yy   <- matrix(y_c[obs, 1L, ii], ncol = 1L)
    zz   <- matrix(z_c[obs, , ii], nrow = length(obs))
    cc_i <- matrix(cc[obs, ii], ncol = 1L)

    u_full <- rep(NA_real_, n_rows)
    u_full[obs] <- as.numeric(yy - zz %*% beta_b - cc_i)
    u_mat  <- t(matrix(u_full, nrow = n_vars))      # T x K, NAs kept
    keep_t <- stats::complete.cases(u_mat)
    sum_tc_b <- sum_tc_b + sum(keep_t)
    u_acc <- u_acc + crossprod(u_mat[keep_t, , drop = FALSE])
  }
  sigma_b <- u_acc / sum_tc_b

  # Build a minimal fake fit-like list to reuse compute_irf
  fake_fit <- list(
    beta      = beta_b,
    sigma     = sigma_b,
    n_vars    = n_vars,
    n_lags    = n_lags
  )
  class(fake_fit) <- "pvarife_result"

  tryCatch(
    compute_irf(fake_fit, n_periods = n_periods, shock = shock,
                diff_vars = diff_vars, identification = identification),
    error = function(e) NULL
  )
}
