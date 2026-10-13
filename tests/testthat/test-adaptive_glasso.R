test_that("adaptive_glasso returns a valid fitted object", {

  set.seed(1)
  x <- matrix(rnorm(60 * 6), nrow = 60, ncol = 6)

  colnames(x) <- paste0("X", 1:6)

  fit <- adaptive_glasso(
    x,
    ridge = 0.10,
    d = 0.32,
    lambda = 0.20,
    standardize = TRUE,
    penalize_diagonal = TRUE
  )

  expect_s3_class(fit, "adaptive_glasso_fit")

  expect_equal(fit$n, 60)
  expect_equal(fit$p, 6)

  expect_equal(dim(fit$precision), c(6, 6))
  expect_equal(dim(fit$covariance), c(6, 6))
  expect_equal(dim(fit$adjacency), c(6, 6))

  expect_equal(
    fit$precision,
    t(fit$precision),
    tolerance = 1e-8
  )

  expect_equal(
    fit$adjacency,
    t(fit$adjacency),
    tolerance = 1e-8
  )

  expect_equal(
    unname(diag(fit$adjacency)),
    rep(0, 6)
  )
})


test_that("ridge stabilization is applied correctly", {

  set.seed(2)
  x <- matrix(rnorm(50 * 5), nrow = 50, ncol = 5)

  fit <- adaptive_glasso(
    x,
    ridge = 0.10,
    d = 0.32,
    lambda = 0.15
  )

  expected <-
    fit$sample_covariance + 0.10 * diag(5)

  expect_equal(
    fit$ridge_covariance,
    expected,
    tolerance = 1e-12
  )
})


test_that("adaptive weights are constructed correctly", {

  set.seed(3)
  x <- matrix(rnorm(70 * 6), nrow = 70, ncol = 6)

  fit <- adaptive_glasso(
    x,
    ridge = 0.10,
    d = 0.32,
    lambda = 0.20,
    penalize_diagonal = FALSE
  )

  expected_weights <-
    1 / (abs(fit$initial_precision) + 0.32)

  diag(expected_weights) <- 0

  expect_equal(
    fit$adaptive_weights,
    expected_weights,
    tolerance = 1e-12
  )

  expect_equal(
    unname(diag(fit$adaptive_weights)),
    rep(0, 6)
  )
})


test_that("adaptive penalty equals lambda times adaptive weights", {

  set.seed(4)
  x <- matrix(rnorm(60 * 6), nrow = 60, ncol = 6)

  fit <- adaptive_glasso(
    x,
    ridge = 0.10,
    d = 0.32,
    lambda = 0.20
  )

  expect_equal(
    fit$adaptive_penalty,
    fit$lambda * fit$adaptive_weights,
    tolerance = 1e-12
  )
})


test_that("coefficient extractor returns requested matrices", {

  set.seed(5)
  x <- matrix(rnorm(60 * 6), nrow = 60, ncol = 6)

  fit <- adaptive_glasso(
    x,
    lambda = 0.20
  )

  expect_equal(
    coef(fit, type = "precision"),
    fit$precision
  )

  expect_equal(
    coef(fit, type = "adjacency"),
    fit$adjacency
  )

  expect_equal(
    coef(fit, type = "weights"),
    fit$adaptive_weights
  )

  expect_equal(
    coef(fit, type = "penalty"),
    fit$adaptive_penalty
  )

  expect_equal(
    coef(fit, type = "covariance"),
    fit$covariance
  )
})


test_that("method works for n < p, n = p, and n > p", {

  set.seed(6)

  x_np <- matrix(rnorm(20 * 30), nrow = 20, ncol = 30)
  x_eq <- matrix(rnorm(20 * 20), nrow = 20, ncol = 20)
  x_ng <- matrix(rnorm(40 * 10), nrow = 40, ncol = 10)

  fit_np <- adaptive_glasso(x_np, lambda = 0.20)
  fit_eq <- adaptive_glasso(x_eq, lambda = 0.20)
  fit_ng <- adaptive_glasso(x_ng, lambda = 0.20)

  expect_s3_class(fit_np, "adaptive_glasso_fit")
  expect_s3_class(fit_eq, "adaptive_glasso_fit")
  expect_s3_class(fit_ng, "adaptive_glasso_fit")

  expect_equal(fit_np$final_error_flag, 0)
  expect_equal(fit_eq$final_error_flag, 0)
  expect_equal(fit_ng$final_error_flag, 0)

  expect_equal(dim(fit_np$precision), c(30, 30))
  expect_equal(dim(fit_eq$precision), c(20, 20))
  expect_equal(dim(fit_ng$precision), c(10, 10))
})


test_that("cross-validation selects a valid lambda", {

  set.seed(7)
  x <- matrix(rnorm(60 * 6), nrow = 60, ncol = 6)

  fit <- adaptive_glasso(
    x,
    lambda = NULL,
    nfolds = 3,
    seed = 123
  )

  expect_true(is.finite(fit$lambda))
  expect_gt(fit$lambda, 0)

  expect_false(is.null(fit$cv))
  expect_true(is.finite(fit$cv_best_index))

  expect_equal(
    fit$final_error_flag,
    0
  )
})


test_that("invalid inputs fail clearly", {

  x_constant <- matrix(
    1,
    nrow = 10,
    ncol = 3
  )

  expect_error(
    adaptive_glasso(
      x_constant,
      lambda = 0.20
    ),
    "non-zero finite standard deviations"
  )

  set.seed(8)

  x_missing <- matrix(
    rnorm(30),
    nrow = 10,
    ncol = 3
  )

  x_missing[1, 1] <- NA_real_

  expect_error(
    adaptive_glasso(
      x_missing,
      lambda = 0.20
    ),
    "finite, non-missing"
  )

  x_valid <- matrix(
    rnorm(40),
    nrow = 10,
    ncol = 4
  )

  expect_error(
    adaptive_glasso(
      x_valid,
      lambda = 0
    ),
    "positive finite"
  )

  expect_error(
    adaptive_glasso(
      x_valid,
      d = 0,
      lambda = 0.20
    ),
    "`d` must be one positive finite number"
  )

  expect_error(
    adaptive_glasso(
      x_valid,
      ridge = -0.10,
      lambda = 0.20
    ),
    "`ridge` must be one non-negative finite number"
  )
})
