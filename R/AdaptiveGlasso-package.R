#' AdaptiveGlasso: Adaptive Graphical Lasso for Sparse Network Estimation
#'
#' @description
#' `AdaptiveGlasso` implements an adaptive graphical lasso method for sparse
#' precision-matrix and network estimation. The method combines ridge
#' stabilization, cross-validated regularization, an initial graphical lasso
#' estimate, and adaptive edge-specific penalties.
#'
#' @details
#' The main fitting function is [adaptive_glasso()].
#'
#' For a data matrix with observations in rows and variables in columns, the
#' method first optionally standardizes the variables and ridge-stabilizes the
#' sample covariance matrix. An initial graphical lasso estimate is then used
#' to construct adaptive weights
#'
#' \deqn{w_{ij} = \frac{1}{|\hat{\omega}^{init}_{ij}| + d},}
#'
#' where \eqn{d > 0} is a stabilization constant. The final edge-specific
#' penalties are
#'
#' \deqn{\rho_{ij} = \lambda w_{ij}.}
#'
#' If `lambda` is not supplied, it is selected by cross-validation using
#' `CVglasso`.
#'
#' The method can be applied when the sample size is smaller than, equal to,
#' or larger than the number of variables.
#'
#' @section Main functions:
#' - [adaptive_glasso()] fits the adaptive graphical lasso model.
#' - [coef.adaptive_glasso_fit()] extracts fitted matrices.
#'
#' @examples
#' set.seed(123)
#'
#' x <- matrix(
#'   rnorm(60 * 6),
#'   nrow = 60,
#'   ncol = 6
#' )
#'
#' fit <- adaptive_glasso(
#'   x,
#'   ridge = 0.10,
#'   d = 0.32,
#'   lambda = 0.20
#' )
#'
#' fit
#'
#' coef(fit, type = "precision")
#'
#' @keywords internal
"_PACKAGE"
