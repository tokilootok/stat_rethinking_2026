test_that("A02 contracts reject invalid counts", {
  expect_error(validate_a02_input(-1, 2), "water")
  expect_error(validate_a02_input(0, 0), "At least one")
  expect_error(validate_a02_input(3, 11, 0), "future_tosses")
})

test_that("A02 grid posterior is normalized and matches Beta mean", {
  posterior <- globe_grid_posterior(3, 11, make_probability_grid(2001))
  exact <- globe_exact_posterior(3, 11)
  expect_equal(sum(posterior$posterior), 1)
  expect_equal(c(exact$alpha, exact$beta, exact$mean), c(4, 12, 0.25))
  expect_equal(sum(posterior$probability * posterior$posterior),
    exact$mean, tolerance = 1e-6)
})

test_that("A02 exact posterior predictive is a valid beta-binomial", {
  predictive <- globe_exact_predictive(3, 11, 5)
  expect_equal(predictive$future_water, 0:5)
  expect_equal(sum(predictive$probability), 1)
  expect_equal(sum(predictive$future_water * predictive$probability), 5 * 4 / 16)
})

test_that("A02 grid and exact posterior predictions agree", {
  posterior <- globe_grid_posterior(3, 11, make_probability_grid(2001))
  grid_prediction <- globe_grid_predictive(posterior, 5)
  exact_prediction <- globe_exact_predictive(3, 11, 5)
  expect_equal(grid_prediction$probability, exact_prediction$probability,
    tolerance = 1e-6)
})

test_that("A02 simulation is reproducible and approximates exact prediction", {
  first <- simulate_globe_predictive(3, 11, 5, draws = 50000, seed = 42)
  second <- simulate_globe_predictive(3, 11, 5, draws = 50000, seed = 42)
  exact <- globe_exact_predictive(3, 11, 5)
  empirical <- tabulate(first$future_water + 1L, nbins = 6) / nrow(first)
  expect_identical(first, second)
  expect_lt(max(abs(empirical - exact$probability)), 0.01)
})
