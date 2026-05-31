## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

- macOS Tahoe 26.5, R 4.5.2 (local, aarch64-apple-darwin20): 0 errors, 0 warnings, 0 notes
- Windows Server 2022, R 4.6.0, win-builder R-release (first submission): 0 errors, 0 warnings, 2 notes (see below)
- Windows Server 2022, R-devel, win-builder R-devel (this submission): result pending at time of submission; win-builder R-release FTP was unavailable

## Notes from Windows check (both resolved)

1. "New submission" — expected for a first CRAN submission.

2. "Possibly misspelled words: GLS, Tugan"
   - GLS is a standard statistical abbreviation (Generalised Least Squares),
     used throughout the econometrics literature.
   - Tugan is the surname of the paper author (Tugan 2021,
     <doi:10.1093/ectj/utaa021>) whose estimator this package implements.
   Both are added to inst/WORDLIST.

3. "Non-standard file: cran-comments.md" — added to .Rbuildignore so it is
   excluded from the built tarball; present only in the source directory.

## Downstream dependencies

This is a new submission. There are no existing downstream dependencies.

## Notes on implementation

- The package implements the estimator of Tugan (2021)
  <doi:10.1093/ectj/utaa021>. The MATLAB replication code is used
  as a reference but is not distributed with the package.
- One intentional deviation from the MATLAB source is documented in
  asymptotic_var(): the B_gamma bias term sums over all HAC lags
  g = 1, ..., G_bar, correcting a bug in the original MATLAB code
  (line 189 of Asymptotic_Distribution_of_beta.m).
