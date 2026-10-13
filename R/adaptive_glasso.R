#' Fit an adaptive graphical lasso model
#'
#' Fits an adaptive graphical lasso model for sparse precision-matrix
#' and network estimation. The procedure combines ridge stabilization,
#' an initial graphical lasso estimate, adaptive edge-specific penalties,
#' and a final weighted graphical lasso estimate.
#'
#' @details
#' The sample covariance matrix is first ridge-stabilized as
#' \deqn{S_{\delta} = S + \delta I,}
#' where `ridge` corresponds to \eqn{\delta}.
#'
#' An initial graphical lasso estimate is then used to construct the
#' adaptive weights
#' \deqn{w_{ij} = \frac{1}{|\hat{\omega}^{init}_{ij}| + d},}
#' followed by the edge-specific penalties
#' \deqn{\rho_{ij} = \lambda w_{ij}.}
#'
#' The value of `lambda` can be handled in two ways:
#'
#' - `lambda = NULL`: `lambda` is selected automatically by
#'   cross-validation using `CVglasso`.
#' - `lambda =` a positive numeric value: the supplied value is used
#'   directly and cross-validation is not performed.
#'
#' The method can be applied when the number of observations is smaller
#' than, equal to, or larger than the number of variables.
#'
#' @param x Numeric matrix or data frame with observations in rows and
#'   variables in columns.
#' @param ridge Non-negative ridge stabilization constant added to the
#'   diagonal of the sample covariance matrix. Default is `0.10`.
#' @param d Positive stabilization constant used in the adaptive weights
#'   \eqn{w_{ij} = 1/(|\hat{\omega}^{init}_{ij}| + d)}.
#'   Default is `0.32`.
#' @param lambda Optional positive regularization parameter. If `NULL`,
#'   `lambda` is selected automatically by cross-validation using
#'   `CVglasso`. If a positive numeric value is supplied, that fixed value
#'   is used and cross-validation is skipped.
#' @param lam.min.ratio Minimum-to-maximum lambda ratio used by
#'   `CVglasso` when `lambda = NULL`. Default is `0.01`.
#' @param nfolds Number of cross-validation folds used when
#'   `lambda = NULL`. Default is `5`.
#' @param standardize Logical; if `TRUE`, variables are standardized
#'   before fitting. Default is `TRUE`.
#' @param seed Optional integer seed used for reproducibility of
#'   cross-validation when `lambda = NULL`. Default is `123`.
#' @param penalize_diagonal Logical; whether diagonal precision-matrix
#'   elements are penalized. Default is `TRUE`, matching the fitting
#'   procedure used in the motivating analysis.
#' @param tol Positive convergence tolerance passed to the graphical
#'   lasso fitting routines. Default is `1e-4`.
#' @param maxit Maximum number of graphical lasso iterations.
#'   Default is `10000`.
#'
#' @return An object of class `"adaptive_glasso_fit"` containing the
#'   fitted precision and covariance matrices, adjacency matrix, initial
#'   precision estimate, adaptive weights, adaptive penalty matrix,
#'   selected or supplied `lambda`, ridge and adaptive stabilization
#'   constants, cross-validation results when applicable, and fitting
#'   diagnostics.
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
#' # Use a fixed lambda: cross-validation is not performed
#' fit <- adaptive_glasso(
#'   x,
#'   lambda = 0.195237
#' )
#'
#' fit
#' coef(fit, type = "precision")
#'
#' # Select lambda automatically by cross-validation
#' fit_cv <- adaptive_glasso(
#'   x,
#'   lambda = NULL,
#'   nfolds = 3
#' )
#'
#' fit_cv$lambda
#'
#' @export
adaptive_glasso <- function(
    x,
    ridge = 0.10,
    d = 0.32,
    lambda = NULL,
    lam.min.ratio = 0.01,
    nfolds = 5L,
    standardize = TRUE,
    seed = 123L,
    penalize_diagonal = TRUE,
    tol = 1e-4,
    maxit = 10000L
) {

  call <- match.call()

  ## ------------------------------------------------------------
  ## 1. Prepare and validate user data
  ## ------------------------------------------------------------

  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }

  if (!is.matrix(x) || !is.numeric(x)) {
    stop(
      "`x` must be a numeric matrix or data frame.",
      call. = FALSE
    )
  }

  if (nrow(x) < 2L || ncol(x) < 2L) {
    stop(
      "`x` must contain at least 2 observations and 2 variables.",
      call. = FALSE
    )
  }

  if (anyNA(x) || any(!is.finite(x))) {
    stop(
      "`x` must contain only finite, non-missing values.",
      call. = FALSE
    )
  }

  if (!is.logical(standardize) ||
      length(standardize) != 1L ||
      is.na(standardize)) {
    stop(
      "`standardize` must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  if (!is.null(lambda) &&
      (!is.numeric(lambda) ||
       length(lambda) != 1L ||
       !is.finite(lambda) ||
       lambda <= 0)) {
    stop(
      "`lambda` must be NULL or one positive finite number.",
      call. = FALSE
    )
  }

  if (!is.numeric(lam.min.ratio) ||
      length(lam.min.ratio) != 1L ||
      !is.finite(lam.min.ratio) ||
      lam.min.ratio <= 0 ||
      lam.min.ratio > 1) {
    stop(
      "`lam.min.ratio` must be in (0, 1].",
      call. = FALSE
    )
  }

  if (!is.numeric(nfolds) ||
      length(nfolds) != 1L ||
      !is.finite(nfolds) ||
      nfolds < 2) {
    stop(
      "`nfolds` must be an integer of at least 2.",
      call. = FALSE
    )
  }

  nfolds <- as.integer(nfolds)

  if (nfolds > nrow(x)) {
    stop(
      "`nfolds` cannot exceed the number of observations.",
      call. = FALSE
    )
  }

  if (is.null(colnames(x))) {
    colnames(x) <- paste0("V", seq_len(ncol(x)))
  }

  variable_names <- colnames(x)

  ## ------------------------------------------------------------
  ## 2. Standardize variables
  ## ------------------------------------------------------------

  if (standardize) {

    sds <- apply(x, 2L, stats::sd)

    if (any(!is.finite(sds)) || any(sds == 0)) {
      stop(
        paste0(
          "All variables must have non-zero finite ",
          "standard deviations before standardization."
        ),
        call. = FALSE
      )
    }

    x_fit <- scale(x)

    center <- attr(x_fit, "scaled:center")
    scale_values <- attr(x_fit, "scaled:scale")

    x_fit <- as.matrix(x_fit)

  } else {

    x_fit <- x
    center <- NULL
    scale_values <- NULL
  }

  ## ------------------------------------------------------------
  ## 3. Select lambda by cross-validation when needed
  ## ------------------------------------------------------------

  cv_fit <- NULL
  cv_best_index <- NA_integer_
  cv_mean_error <- NULL

  if (is.null(lambda)) {

    if (!is.null(seed)) {
      set.seed(seed)
    }

    cv_fit <- CVglasso::CVglasso(
      X = x_fit,
      lam.min.ratio = lam.min.ratio,
      nobs = nrow(x_fit),
      K = nfolds,
      diagonal = FALSE,
      tol = tol,
      maxit = maxit,
      trace = "none"
    )

    cv_mean_error <- rowMeans(cv_fit$CV.error)

    cv_best_index <- which.min(cv_mean_error)

    lambda <- cv_fit$Lambdas[cv_best_index]
  }

  ## ------------------------------------------------------------
  ## 4. Fit adaptive graphical lasso core
  ## ------------------------------------------------------------

  fit <- .adaptive_glasso_fit_matrix(
    x = x_fit,
    ridge = ridge,
    d = d,
    lambda = lambda,
    penalize_diagonal = penalize_diagonal,
    tol = tol,
    maxit = maxit
  )

  ## ------------------------------------------------------------
  ## 5. Add user-facing information
  ## ------------------------------------------------------------

  out <- c(
    list(
      call = call,
      standardized_data = x_fit,
      center = center,
      scale = scale_values,
      standardize = standardize,
      lam.min.ratio = lam.min.ratio,
      nfolds = nfolds,
      cv = cv_fit,
      cv_best_index = cv_best_index,
      cv_mean_error = cv_mean_error
    ),
    fit
  )

  out$variable_names <- variable_names

  class(out) <- "adaptive_glasso_fit"

  out
}

#' Print an adaptive graphical lasso fit
#'
#' @param x An object returned by `adaptive_glasso()`.
#' @param ... Additional arguments; currently unused.
#'
#' @return The fitted object, invisibly.
#'
#' @export
print.adaptive_glasso_fit <- function(x, ...) {

  cat("Adaptive graphical lasso fit\n")
  cat("-----------------------------\n")
  cat("Observations:", x$n, "\n")
  cat("Variables:", x$p, "\n")
  cat("Ridge stabilization:", x$ridge, "\n")
  cat("Adaptive stabilization d:", x$d, "\n")
  cat("Selected lambda:", signif(x$lambda, 6), "\n")
  cat("Estimated edges:", x$edge_count, "\n")
  cat("Penalize diagonal:", x$penalize_diagonal, "\n")
  cat("Final glasso error flag:", x$final_error_flag, "\n")

  invisible(x)
}

#' Extract matrices from an adaptive graphical lasso fit
#'
#' @param object An object returned by `adaptive_glasso()`.
#' @param type Matrix to return.
#' @param ... Additional arguments; currently unused.
#'
#' @return A matrix.
#'
#' @export
coef.adaptive_glasso_fit <- function(
    object,
    type = c(
      "precision",
      "adjacency",
      "weights",
      "penalty",
      "covariance"
    ),
    ...
) {

  type <- match.arg(type)

  switch(
    type,
    precision = object$precision,
    adjacency = object$adjacency,
    weights = object$adaptive_weights,
    penalty = object$adaptive_penalty,
    covariance = object$covariance
  )
}
