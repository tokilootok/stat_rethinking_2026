# A02: grid approximation and exact beta-binomial posterior prediction.

make_probability_grid <- function(size = 1001L) {
  size <- assert_scalar_count(size, "size", lower = 2L)
  seq(0, 1, length.out = size)
}

uniform_prior <- function(grid) {
  if (!is.numeric(grid) || any(!is.finite(grid)) || any(grid < 0 | grid > 1)) {
    stop("`grid` must contain finite probabilities between zero and one.", call. = FALSE)
  }
  rep.int(1, length(grid))
}

globe_grid_posterior <- function(water, land, grid, prior = uniform_prior(grid)) {
  input <- validate_a02_input(water, land)
  if (!is.numeric(grid) || length(grid) < 2L ||
      any(!is.finite(grid)) || any(grid < 0 | grid > 1)) {
    stop("`grid` must contain at least two finite probabilities in [0, 1].", call. = FALSE)
  }
  if (length(prior) != length(grid)) {
    stop("`prior` and `grid` must have the same length.", call. = FALSE)
  }
  log_likelihood <- dbinom(input$water, size = input$observations,
    prob = grid, log = TRUE)
  log_weight <- log_likelihood + log(prior)
  finite_weight <- is.finite(log_weight)
  if (!any(finite_weight)) {
    stop("The prior and likelihood have no overlapping support.", call. = FALSE)
  }
  relative_weight <- numeric(length(grid))
  relative_weight[finite_weight] <- exp(
    log_weight[finite_weight] - max(log_weight[finite_weight]))
  data.frame(
    probability = grid,
    prior = normalize_probability(prior, "prior"),
    likelihood = exp(log_likelihood),
    posterior = normalize_probability(relative_weight, "posterior weights")
  )
}

globe_exact_posterior <- function(water, land, alpha = 1, beta = 1) {
  input <- validate_a02_input(water, land)
  alpha <- assert_positive_scalar(alpha, "alpha")
  beta <- assert_positive_scalar(beta, "beta")
  alpha_post <- alpha + input$water
  beta_post <- beta + input$land
  list(
    alpha = alpha_post,
    beta = beta_post,
    mean = alpha_post / (alpha_post + beta_post),
    mode = if (alpha_post > 1 && beta_post > 1) {
      (alpha_post - 1) / (alpha_post + beta_post - 2)
    } else NA_real_
  )
}

beta_binomial_pmf <- function(successes, size, alpha, beta) {
  size <- assert_scalar_count(size, "size", lower = 0L)
  alpha <- assert_positive_scalar(alpha, "alpha")
  beta <- assert_positive_scalar(beta, "beta")
  if (!is.numeric(successes) || any(!is.finite(successes)) ||
      any(successes != floor(successes)) || any(successes < 0 | successes > size)) {
    stop("`successes` must contain integers between zero and `size`.", call. = FALSE)
  }
  exp(lchoose(size, successes) +
    lbeta(successes + alpha, size - successes + beta) -
    lbeta(alpha, beta))
}

globe_exact_predictive <- function(water, land, future_tosses = 5L,
                                   alpha = 1, beta = 1) {
  input <- validate_a02_input(water, land, future_tosses)
  posterior <- globe_exact_posterior(input$water, input$land, alpha, beta)
  future_water <- 0:input$future_tosses
  probability <- beta_binomial_pmf(future_water,
    size = input$future_tosses, alpha = posterior$alpha,
    beta = posterior$beta)
  data.frame(future_water = future_water,
    probability = normalize_probability(probability))
}

globe_grid_predictive <- function(grid_posterior, future_tosses = 5L) {
  future_tosses <- assert_scalar_count(future_tosses, "future_tosses", lower = 1L)
  required <- c("probability", "posterior")
  if (!is.data.frame(grid_posterior) || !all(required %in% names(grid_posterior))) {
    stop("`grid_posterior` must contain probability and posterior columns.", call. = FALSE)
  }
  future_water <- 0:future_tosses
  predictive <- vapply(future_water, function(k) {
    sum(dbinom(k, size = future_tosses,
      prob = grid_posterior$probability) * grid_posterior$posterior)
  }, numeric(1))
  data.frame(future_water = future_water,
    probability = normalize_probability(predictive))
}

simulate_globe_predictive <- function(water, land, future_tosses = 5L,
                                      draws = 10000L, alpha = 1, beta = 1,
                                      seed = NULL) {
  input <- validate_a02_input(water, land, future_tosses)
  draws <- assert_scalar_count(draws, "draws", lower = 1L)
  posterior <- globe_exact_posterior(input$water, input$land, alpha, beta)
  if (!is.null(seed)) set.seed(seed)
  probability_draws <- rbeta(draws, posterior$alpha, posterior$beta)
  future_water <- rbinom(draws, input$future_tosses, probability_draws)
  data.frame(draw = seq_len(draws), probability = probability_draws,
    future_water = future_water)
}

summarize_a02 <- function(grid_posterior, exact_posterior,
                          exact_predictive, grid_predictive) {
  list(
    posterior_family = sprintf("Beta(%g, %g)",
      exact_posterior$alpha, exact_posterior$beta),
    exact_mean = exact_posterior$mean,
    grid_mean = sum(grid_posterior$probability * grid_posterior$posterior),
    predictive_mean = sum(
      exact_predictive$future_water * exact_predictive$probability),
    maximum_grid_exact_difference = max(abs(
      exact_predictive$probability - grid_predictive$probability))
  )
}
