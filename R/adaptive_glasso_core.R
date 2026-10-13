.adaptive_glasso_fit_matrix <- function(
    x,
    ridge = 0.10,
    d = 0.32,
    lambda,
    penalize_diagonal = TRUE,
    tol = 1e-4,
    maxit = 10000L
) {

  ## ------------------------------------------------------------
  ## 1. Validate inputs
  ## ------------------------------------------------------------

  if (!is.matrix(x) || !is.numeric(x)) {
    stop("`x` must be a numeric matrix.", call. = FALSE)
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

  if (!is.numeric(ridge) ||
      length(ridge) != 1L ||
      !is.finite(ridge) ||
      ridge < 0) {
    stop(
      "`ridge` must be one non-negative finite number.",
      call. = FALSE
    )
  }

  if (!is.numeric(d) ||
      length(d) != 1L ||
      !is.finite(d) ||
      d <= 0) {
    stop(
      "`d` must be one positive finite number.",
      call. = FALSE
    )
  }

  if (!is.numeric(lambda) ||
      length(lambda) != 1L ||
      !is.finite(lambda) ||
      lambda <= 0) {
    stop(
      "`lambda` must be one positive finite number.",
      call. = FALSE
    )
  }

  if (!is.logical(penalize_diagonal) ||
      length(penalize_diagonal) != 1L ||
      is.na(penalize_diagonal)) {
    stop(
      "`penalize_diagonal` must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  if (!is.numeric(tol) ||
      length(tol) != 1L ||
      !is.finite(tol) ||
      tol <= 0) {
    stop(
      "`tol` must be one positive finite number.",
      call. = FALSE
    )
  }

  if (!is.numeric(maxit) ||
      length(maxit) != 1L ||
      !is.finite(maxit) ||
      maxit < 1) {
    stop(
      "`maxit` must be a positive integer.",
      call. = FALSE
    )
  }

  maxit <- as.integer(maxit)

  n <- nrow(x)
  p <- ncol(x)

  variable_names <- colnames(x)

  if (is.null(variable_names)) {
    variable_names <- paste0("V", seq_len(p))
    colnames(x) <- variable_names
  }

  ## ------------------------------------------------------------
  ## 2. Sample covariance and ridge stabilization
  ##
  ##    S_delta = S + ridge * I
  ## ------------------------------------------------------------

  sample_covariance <- stats::cov(x)

  ridge_covariance <-
    sample_covariance + ridge * diag(p)

  dimnames(sample_covariance) <-
    list(variable_names, variable_names)

  dimnames(ridge_covariance) <-
    list(variable_names, variable_names)

  ## ------------------------------------------------------------
  ## 3. Initial graphical lasso
  ## ------------------------------------------------------------

  initial_fit <- glasso::glasso(
    s = ridge_covariance,
    rho = lambda,
    nobs = n,
    thr = tol,
    maxit = maxit,
    penalize.diagonal = penalize_diagonal
  )

  initial_precision <- initial_fit$wi

  dimnames(initial_precision) <-
    list(variable_names, variable_names)

  ## ------------------------------------------------------------
  ## 4. Adaptive weights
  ##
  ##    w_ij = 1 / (|omega_init,ij| + d)
  ## ------------------------------------------------------------

  adaptive_weights <-
    1 / (abs(initial_precision) + d)

  if (!penalize_diagonal) {
    diag(adaptive_weights) <- 0
  }

  dimnames(adaptive_weights) <-
    list(variable_names, variable_names)

  ## ------------------------------------------------------------
  ## 5. Adaptive penalty matrix
  ##
  ##    rho_ij = lambda * w_ij
  ## ------------------------------------------------------------

  adaptive_penalty <-
    lambda * adaptive_weights

  dimnames(adaptive_penalty) <-
    list(variable_names, variable_names)

  ## ------------------------------------------------------------
  ## 6. Final adaptive graphical lasso
  ## ------------------------------------------------------------

  adaptive_fit <- glasso::glasso(
    s = ridge_covariance,
    rho = adaptive_penalty,
    nobs = n,
    thr = tol,
    maxit = maxit,
    penalize.diagonal = penalize_diagonal
  )

  precision <- adaptive_fit$wi
  covariance <- adaptive_fit$w

  dimnames(precision) <-
    list(variable_names, variable_names)

  dimnames(covariance) <-
    list(variable_names, variable_names)

  ## ------------------------------------------------------------
  ## 7. Network adjacency matrix
  ## ------------------------------------------------------------

  adjacency <- abs(precision)
  diag(adjacency) <- 0

  edge_count <-
    sum(adjacency[upper.tri(adjacency)] != 0)

  ## ------------------------------------------------------------
  ## 8. Return internal fit
  ## ------------------------------------------------------------

  list(
    n = n,
    p = p,
    variable_names = variable_names,
    sample_covariance = sample_covariance,
    ridge_covariance = ridge_covariance,
    initial_precision = initial_precision,
    adaptive_weights = adaptive_weights,
    adaptive_penalty = adaptive_penalty,
    precision = precision,
    covariance = covariance,
    adjacency = adjacency,
    edge_count = edge_count,
    ridge = ridge,
    d = d,
    lambda = lambda,
    penalize_diagonal = penalize_diagonal,
    initial_loglik = initial_fit$loglik,
    final_loglik = adaptive_fit$loglik,
    initial_iterations = initial_fit$niter,
    final_iterations = adaptive_fit$niter,
    initial_error_flag = initial_fit$errflag,
    final_error_flag = adaptive_fit$errflag
  )
}
