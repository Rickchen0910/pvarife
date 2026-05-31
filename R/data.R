#' Simulated panel VAR dataset with interactive fixed effects
#'
#' A synthetic panel dataset generated from the DGP of Tugan (2021), Section
#' S10, via \code{\link{sim_pvarife}}. Provided as a compact example dataset
#' for vignettes and function examples.
#'
#' @format A list with the following components:
#' \describe{
#'   \item{y}{Numeric array of dimension \eqn{50 \times 30 \times 2}
#'     (50 units, 30 time periods, 2 variables).}
#'   \item{beta_true}{True coefficient vector of length 6 (2 intercepts + 4
#'     VAR lag coefficients).}
#'   \item{theta_true}{True VAR coefficient array of dimension
#'     \eqn{2 \times 2 \times 1}.}
#'   \item{sigma_true}{True reduced-form covariance matrix
#'     \eqn{2 \times 2}.}
#'   \item{factors_true}{True factor matrix of dimension
#'     \eqn{30 \times 1}.}
#'   \item{loadings_true}{True factor loadings matrix of dimension
#'     \eqn{2 \times 50}.}
#' }
#'
#' @details
#' Generated with \code{sim_pvarife(n_units = 50, n_time = 30, n_vars = 2,
#' n_lags = 1, n_factors = 1, seed = 42)}.
#'
#' @references
#' Tugan, M. (2021). Panel VAR models with interactive fixed effects.
#' \emph{Econometrics Journal}, 24, 225--246.
#' \doi{10.1093/ectj/utaa021}
#'
#' @source Simulated; see \code{\link{sim_pvarife}}.
"pvarife_sim"
