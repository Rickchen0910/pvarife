## Resubmission

This is a resubmission addressing the CRAN reviewer's comments:

1. **Explain all acronyms in the Description text.**
   The Description has been rewritten so that the only acronym, "VAR", is
   spelled out on first use as "vector autoregression (VAR)". The previous
   abbreviations "GLS" and "EM" have been removed; the estimation method is
   now described in full words ("an iterative algorithm that alternates
   between principal component estimation of the factors and least squares
   estimation of the VAR coefficients").

2. **Replace \dontrun{} with \donttest{} or unwrap fast examples.**
   - The example for `irf_bands()` has been unwrapped; it now runs in
     well under 5 seconds.
   - The example for `bootstrap_irf_bands()` re-estimates the model on many
     resamples and takes longer than 5 seconds, so its example is now wrapped
     in \donttest{} (it is fully executable, not \dontrun{}).

There are no remaining uses of \dontrun{} in the package.

## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

- macOS Tahoe 26.5, R 4.5.2 (local, aarch64-apple-darwin20): 0 errors, 0 warnings, 0 notes
- Windows Server 2022, R-devel and R-release (win-builder)

## Downstream dependencies

This is a new submission. There are no existing downstream dependencies.

## Notes on implementation

- The package implements the estimator of Tugan (2021)
  <doi:10.1093/ectj/utaa021>. The MATLAB replication code is used as a
  reference but is not distributed with the package.
- One intentional deviation from the MATLAB source is documented in
  asymptotic_var(): the bias term sums over all HAC lags g = 1, ..., G_bar,
  correcting a bug in the original MATLAB code.
