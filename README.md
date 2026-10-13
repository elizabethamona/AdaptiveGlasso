# AdaptiveGlasso

## Adaptive Graphical Lasso for Sparse Precision-Matrix and Network Estimation

**AdaptiveGlasso** is an R package implementing an **adaptive graphical lasso** approach for sparse precision-matrix and network estimation.

The method combines ridge stabilization, data-driven selection of the regularization parameter, an initial graphical lasso estimate, and adaptive edge-specific penalties. It was developed for settings in which stable network estimation may be difficult because the number of variables is large relative to the available sample size, while remaining applicable when the sample size is equal to or larger than the number of variables.

The package provides a single user-facing fitting function, `adaptive_glasso()`, together with methods for inspecting fitted models and extracting the estimated precision, covariance, adjacency, adaptive-weight, and penalty matrices.

---

## Why AdaptiveGlasso?

Estimating sparse conditional-dependence networks can be difficult in high-dimensional and small-sample settings, particularly when the number of variables approaches or exceeds the number of observations.

The standard graphical lasso estimates a sparse precision matrix using a common regularization structure. Adaptive penalization allows the amount of shrinkage to vary across potential edges according to an initial estimate of the precision matrix.

AdaptiveGlasso combines this idea with ridge stabilization of the sample covariance matrix. The approach therefore separates two roles:

1. **ridge stabilization** improves numerical stability of the covariance estimate; and
2. **adaptive penalization** applies edge-specific regularization in the final graphical lasso fit.

This structure is especially useful for sparse network estimation under challenging \(n<p\) conditions, but the method is not restricted to that setting.

---

## Method at a glance

For an \(n \times p\) data matrix \(X\), AdaptiveGlasso proceeds through the following steps.

### 1. Standardization

By default, variables are standardized before model fitting.

### 2. Ridge stabilization

The sample covariance matrix \(S\) is stabilized as

$$
S_{\delta} = S + \delta I,
$$

where \(\delta\) is controlled by the `ridge` argument.

The default is

```r
ridge = 0.10
```

### 3. Selection of the regularization parameter

The regularization parameter \(\lambda\) can be handled in two ways:

- `lambda = NULL` — select \(\lambda\) automatically by cross-validation using `CVglasso`;
- `lambda =` a positive numeric value — use the supplied value directly and skip cross-validation.

### 4. Initial graphical lasso estimate

An initial graphical lasso model is fitted to the ridge-stabilized covariance matrix to obtain an initial precision-matrix estimate,

$$
\widehat{\Omega}^{\,\mathrm{init}}.
$$

### 5. Adaptive weights

The edge-specific adaptive weights are

$$
w_{ij}
=
\frac{1}
{
\left|
\widehat{\omega}^{\,\mathrm{init}}_{ij}
\right| + d
},
$$

where \(d>0\) is a stabilization constant.

The default is

```r
d = 0.32
```

### 6. Adaptive penalties

The final edge-specific penalties are

$$
\rho_{ij}
=
\lambda w_{ij}.
$$

### 7. Final adaptive graphical lasso

A second graphical lasso model is fitted using the adaptive penalty matrix to obtain the final precision-matrix estimate.

The resulting fitted object contains the final precision and covariance matrices together with the quantities used throughout the adaptive fitting procedure.

---

## Installation

AdaptiveGlasso is currently available from GitHub and is being prepared for CRAN submission.

The package can be installed from GitHub using `remotes`:

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
}

remotes::install_github("elizabethamona/AdaptiveGlasso")
```

After installation:

```r
library(AdaptiveGlasso)
```

---

## Basic usage

The primary user-facing function is:

```r
adaptive_glasso()
```

The function accepts a numeric matrix or data frame with:

- observations in rows; and
- variables in columns.

For illustration:

```r
library(AdaptiveGlasso)

set.seed(123)

x <- matrix(
    rnorm(60 * 6),
    nrow = 60,
    ncol = 6
)
```

---

## Selecting lambda by cross-validation

The default behavior is to select the regularization parameter automatically.

```r
fit_cv <- adaptive_glasso(
    x = x,
    ridge = 0.10,
    d = 0.32,
    lambda = NULL
)
```

Because `lambda = NULL`, cross-validation is performed using `CVglasso`.

The selected value can be inspected using:

```r
fit_cv$lambda
```

and the fitted model can be printed directly:

```r
fit_cv
```

---

## Using a fixed lambda

Users may instead supply a positive numeric value for `lambda`.

For example:

```r
fit_fixed <- adaptive_glasso(
    x = x,
    ridge = 0.10,
    d = 0.32,
    lambda = 0.195237
)
```

When a numeric value is supplied, that value is used directly and cross-validation is not performed.

Thus:

```r
lambda = NULL
```

means **select lambda by cross-validation**, whereas

```r
lambda = 0.195237
```

means **use the supplied fixed value**.

---

## Default settings

The principal defaults are:

| Argument | Default | Purpose |
|---|---:|---|
| `ridge` | `0.10` | Ridge stabilization added to the covariance diagonal |
| `d` | `0.32` | Stabilization constant used in the adaptive weights |
| `lambda` | `NULL` | Select lambda by cross-validation |
| `lam.min.ratio` | `0.01` | Minimum-to-maximum lambda ratio used by `CVglasso` |
| `nfolds` | `5` | Number of cross-validation folds |
| `standardize` | `TRUE` | Standardize variables before fitting |
| `seed` | `123` | Seed used for reproducible cross-validation |
| `penalize_diagonal` | `TRUE` | Penalize diagonal elements in the graphical lasso fitting steps |
| `tol` | `1e-4` | Graphical lasso convergence tolerance |
| `maxit` | `10000` | Maximum number of graphical lasso iterations |

During automatic lambda selection, `CVglasso` uses its unpenalized-diagonal cross-validation formulation. The initial and final graphical lasso fits then use the value supplied through `penalize_diagonal`, which defaults to `TRUE`.

---

## Input requirements

The input supplied to `adaptive_glasso()` must be a numeric matrix or data frame.

The package requires:

- at least two observations;
- at least two variables;
- finite, non-missing values; and
- non-zero finite variable standard deviations when `standardize = TRUE`.

If column names are not supplied, AdaptiveGlasso automatically creates variable names of the form:

```r
V1, V2, V3, ...
```

By default:

```r
standardize = TRUE
```

so each variable is centered and scaled before covariance estimation.

---

## High-dimensional and conventional settings

AdaptiveGlasso can be fitted when:

$$
n < p,
$$

$$
n = p,
$$

or

$$
n > p.
$$

The method was motivated by high-dimensional network estimation where \(n<p\), but it is not restricted to that regime.

Automated package tests cover all three sample-size-to-dimension settings.

---

## Model output

An `adaptive_glasso()` fit returns an object of class:

```r
"adaptive_glasso_fit"
```

Printing the fitted object provides a concise summary:

```r
fit_cv
```

including:

- number of observations;
- number of variables;
- ridge stabilization constant;
- adaptive stabilization constant;
- selected or supplied lambda;
- estimated number of edges;
- diagonal-penalization setting; and
- final graphical lasso error flag.

Additional fitted quantities remain available directly from the object.

For example:

```r
fit_cv$lambda
fit_cv$edge_count
fit_cv$adaptive_weights
fit_cv$adaptive_penalty
fit_cv$precision
fit_cv$covariance
```

---

## Extracting fitted matrices

The `coef()` method provides a convenient interface for extracting the principal fitted matrices.

### Precision matrix

```r
coef(fit_cv, type = "precision")
```

### Adjacency matrix

```r
coef(fit_cv, type = "adjacency")
```

### Adaptive-weight matrix

```r
coef(fit_cv, type = "weights")
```

### Adaptive-penalty matrix

```r
coef(fit_cv, type = "penalty")
```

### Estimated covariance matrix

```r
coef(fit_cv, type = "covariance")
```

The adjacency matrix is constructed from the absolute values of the estimated off-diagonal precision-matrix entries, with its diagonal set to zero.

The signed conditional-association information remains available in the estimated precision matrix itself.

---

## Cross-validation output

When:

```r
lambda = NULL
```

the fitted object also stores the cross-validation results.

These include:

```r
fit_cv$cv
fit_cv$cv_best_index
fit_cv$cv_mean_error
fit_cv$lambda
```

When a numeric lambda is supplied directly, cross-validation is skipped.

---

## Downstream network analysis

AdaptiveGlasso estimates the sparse precision matrix and corresponding adjacency structure but does not prescribe a particular downstream network-analysis workflow.

The fitted matrices can be used with appropriate R tools for analyses such as:

- network visualization;
- community detection;
- node-strength analysis;
- closeness centrality;
- weighted betweenness centrality; and
- sensitivity or stability analyses.

Interpretation of individual edges and network statistics should account for the study design, sample size, estimation uncertainty, and scientific context.

---

## Data and reproducibility

The motivating application for AdaptiveGlasso uses a 41-variable mental-health symptom network with a small-sample, high-dimensional structure.

The original participant-level data are **not distributed with this package** because they are not public data.

For this reason, package documentation and examples use simulated data rather than distributing the motivating dataset.

The implementation has been checked against the motivating analysis to ensure that the package reproduces the intended fitting procedure.

---

## Testing

AdaptiveGlasso includes automated tests covering:

- fitted-object construction;
- ridge stabilization;
- adaptive-weight construction;
- adaptive-penalty construction;
- precision-matrix symmetry;
- adjacency-matrix construction;
- coefficient extraction;
- automatic cross-validation;
- invalid-input handling; and
- fitting under \(n<p\), \(n=p\), and \(n>p\) settings.

The package is also checked using standard R package checking tools before release.

---

## Documentation

Function-level documentation is available through:

```r
?adaptive_glasso
```

Package-level documentation is available through:

```r
?AdaptiveGlasso
```

The fitted-object methods can also be inspected using:

```r
?coef.adaptive_glasso_fit
?print.adaptive_glasso_fit
```

---

## Associated manuscript

AdaptiveGlasso accompanies the methodological manuscript:

**An Adaptive Graphical Lasso Approach to Modeling Symptom Networks of Common Mental Disorders in Eritrean Refugee Population**

The manuscript develops and evaluates the adaptive graphical lasso approach through methodological formulation, simulation studies, sensitivity analyses, and a real-world symptom-network application.

The manuscript is currently being finalized for submission. The complete publication citation will be added here when available.

---

## Software author and maintainer

**Elizabeth B. Amona, PhD**  
Biostatistics Core, Brown Cancer Center  
University of Louisville  
ORCID: https://orcid.org/0000-0002-7575-788X

---

## License

AdaptiveGlasso is released under the **MIT License**.

Copyright © 2026 Elizabeth B. Amona.