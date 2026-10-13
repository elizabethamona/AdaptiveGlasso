# AdaptiveGlasso

`AdaptiveGlasso` implements an adaptive graphical lasso method for sparse
precision-matrix and network estimation.

The method combines:

- ridge stabilization of the sample covariance matrix;
- cross-validated selection of the regularization parameter;
- an initial graphical lasso estimate;
- adaptive edge-specific penalty weights; and
- a final weighted graphical lasso estimate.

The adaptive weights are defined as

\[
w_{ij} =
\frac{1}{
|\widehat{\omega}^{\,\mathrm{init}}_{ij}| + d
},
\]

and the corresponding adaptive penalties are

\[
\rho_{ij} = \lambda w_{ij}.
\]

The package supports settings in which the number of observations may be
smaller than, equal to, or larger than the number of variables.

## Installation

The package is currently under development.

## Basic usage

```r
library(AdaptiveGlasso)

set.seed(123)

x <- matrix(
  rnorm(60 * 6),
  nrow = 60,
  ncol = 6
)

fit <- adaptive_glasso(
  x,
  ridge = 0.10,
  d = 0.32,
  lambda = 0.20
)

fit
