# pvarife

**Panel VAR Models with Interactive Fixed Effects**

An R implementation of the estimator proposed by Tuğan (2021, *Econometrics Journal*) for panel vector autoregression (VAR) models with interactive fixed effects. The package jointly estimates VAR coefficients, latent common factors, and unit-specific factor loadings via iterated GLS, following an inner/outer EM algorithm. Asymptotic inference is based on Theorem 2.3 of the paper.

---

## Installation

```r
# From CRAN (once accepted):
install.packages("pvarife")

# Development version from GitHub:
# install.packages("remotes")
remotes::install_github("Rickchen0910/pvarife")
```

---

## Quick Start

```r
library(pvarife)

# 1. Simulate a balanced panel (I=50 units, T=30 periods, K=2 variables)
sim <- sim_pvarife(n_units = 50, n_time = 30, n_vars = 2,
                   n_lags = 1, n_factors = 1, seed = 42)

# 2. Estimate
fit <- pvarife(sim$y, n_lags = 1, n_factors = 1)
print(fit)

# 3. Impulse response functions (short-run / Cholesky identification)
irf <- compute_irf(fit, n_periods = 21, shock = 1)
plot(irf)

# 4. 95% confidence bands (parametric simulation from asymptotic distribution)
bands <- irf_bands(fit, n_periods = 21, shock = 1, n_draw = 500, seed = 1)
plot(bands, normalise_by_h1 = TRUE)

# 5. Asymptotic bias and variance (Theorem 2.3)
avar <- asymptotic_var(fit)
summary(fit)
```

---

## Model

The model is

$$y_{i,t} = \sum_{\ell=1}^{L} \Theta_\ell\, y_{i,t-\ell} + F_t \lambda_i + e_{i,t}$$

where $y_{i,t}$ is a $K \times 1$ vector for unit $i$ at time $t$, $F_t$ is an $r \times 1$ vector of unobservable common factors, $\lambda_i$ is a unit-specific loading vector, and $e_{i,t}$ is idiosyncratic noise. The data array follows the convention `y[i, t, k]` (unit × time × variable). Unbalanced panels (`NA` entries) are supported.

---

## Key Functions

| Function | Description |
|----------|-------------|
| `pvarife()` | Main estimator — returns a `pvarife_result` object |
| `asymptotic_var()` | Theorem 2.3 bias and variance-covariance matrix |
| `compute_irf()` | Impulse response functions (short-run or long-run identification) |
| `irf_bands()` | Parametric 95% confidence bands from joint asymptotic distribution |
| `bootstrap_irf_bands()` | Classical residual bootstrap confidence bands |
| `sim_pvarife()` | Simulate data from the DGP of Tugan (2021, S10) |
| `extract_factors()` | Extract factors and loadings at an arbitrary coefficient vector |

---

## Identification Schemes

```r
# Short-run (recursive / Cholesky) — default
irf_sr <- compute_irf(fit, n_periods = 21, identification = "short_run")

# Long-run (Blanchard-Quah)
irf_lr <- compute_irf(fit, n_periods = 21, identification = "long_run",
                      diff_vars = 1L)
```

---

## Monte Carlo Replication

Scripts to replicate Table S.2 (bias/RMSE/coverage) and Figures S.1-S.2 (average IRFs) from the online appendix are in `scripts/`. They are designed for an SGE high-performance computing cluster and use the [MonteCarlo](https://CRAN.R-project.org/package=MonteCarlo) package.

```bash
# On HPC (SGE scheduler):
bash scripts/meta_mc.sh     # Table S.2
bash scripts/meta_irf.sh    # Figure S.1 (short-run)
```

---

## Citation

If you use this package, please cite the original paper:

> Tugan, M. (2021). Panel VAR models with interactive fixed effects. *Econometrics Journal*, 24, 225-246. <https://doi.org/10.1093/ectj/utaa021>

```bibtex
@article{tugan2021,
  author  = {Tu\u{g}an, Mustafa},
  title   = {Panel {VAR} models with interactive fixed effects},
  journal = {Econometrics Journal},
  year    = {2021},
  volume  = {24},
  pages   = {225--246},
  doi     = {10.1093/ectj/utaa021}
}
```

---

## Author

Binzhi Chen — University of Essex (<Binzhi.Chen9@gmail.com>)
