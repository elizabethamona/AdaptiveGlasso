# AdaptiveGlasso 1.0.0

## Initial release

- Implements adaptive graphical lasso estimation for sparse precision-matrix and network estimation.
- Supports automatic selection of the regularization parameter by cross-validation using `CVglasso`.
- Supports user-supplied fixed values of `lambda`.
- Uses ridge stabilization of the sample covariance matrix.
- Constructs adaptive edge-specific weights and penalty matrices from an initial graphical lasso estimate.
- Provides fitted precision, covariance, adjacency, adaptive-weight, and adaptive-penalty matrices.
- Includes methods for printing fitted objects and extracting fitted matrices with `coef()`.
- Supports settings with \(n < p\), \(n = p\), and \(n > p\).
- Includes automated unit tests and documentation for CRAN release.
