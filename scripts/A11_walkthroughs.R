# A11 — Monsters and Mixtures: original tidyverse teaching examples.
# Run from the repository root: Rscript scripts/A11_walkthroughs.R
# No book datasets or Stan compiler required. All data below are simulated.
# Grids approximate Bayesian posteriors; their range/resolution checks are printed.
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(tibble)
  library(readr); library(ggplot2)
})
options(width = 105, digits = 4)
out <- "courses/figures/A11"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
show_table <- function(x) print(as.data.frame(x), row.names = FALSE, digits = 4)
normalize_grid <- function(g) {
  g |> mutate(weight = exp(log_weight - max(log_weight)), weight = weight / sum(weight))
}
weighted_summary <- function(g, variables) {
  map_dfr(variables, function(v) {
    z <- g |> arrange(.data[[v]])
    q <- function(p) z[[v]][which(cumsum(z$weight) >= p)[1]]
    tibble(parameter = v, mean = sum(z[[v]] * z$weight), lo90 = q(.05), hi90 = q(.95))
  })
}
# A boundary check concerns finite support; refinement concerns grid resolution.
check_grid <- function(coarse, fine, variables, boundary) {
  a <- weighted_summary(coarse, variables)
  b <- weighted_summary(fine, variables)
  error <- max(abs(a$mean - b$mean))
  edge_mass <- sum(coarse$weight[boundary])
  stopifnot(error < .02, edge_mass < .001)
  tibble(max_mean_change = error, boundary_mass = edge_mass)
}
log_beta_binomial <- function(y, n, mu, kappa) {
  lchoose(n, y) + lbeta(y + mu*kappa, n - y + (1-mu)*kappa) - lbeta(mu*kappa, (1-mu)*kappa)
}
# BEGIN_6
cat("BEGIN_6\n")
# Problem 1A: each campaign cohort shares an unobserved conversion probability.
mu <- .20; n <- 30
bb_moments <- tibble(kappa = c(2, 10, 50, Inf)) |>
  mutate(mean = n*mu, variance_p = mu*(1-mu)/(kappa+1),
         variance_y = n*mu*(1-mu)*(1+(n-1)/(kappa+1)),
         within_cohort_correlation = 1/(kappa+1))
show_table(bb_moments)
# Same cohort-level P is shared by all 30 opportunities. New cohort, new P.
set.seed(1106)
bb_data <- tibble(cohort = 1:100, n = n, p_true = rbeta(100, mu*10, (1-mu)*10)) |>
  mutate(y = rbinom(n(), n, p_true))
show_table(bb_data |> summarise(cohorts = n(), mean_count = mean(y), variance_count = var(y)))
# END_6
cat("END_6\n")

# BEGIN_8
cat("BEGIN_8\n")
# Learn both the population mean and heterogeneity from repeated cohorts.
# Priors on transformed coordinates: logit(mu) ~ N(logit(.2),1), log(kappa) ~ N(log(10),1).
fit_bb <- function(step) {
  g <- expand_grid(a = seq(-4, .8, by = step), h = seq(-2, 7, by = step)) |>
    mutate(mu = plogis(a), kappa = exp(h),
           log_weight = dnorm(a, qlogis(.2), 1, log = TRUE) + dnorm(h, log(10), 1, log = TRUE))
  for (i in seq_len(nrow(bb_data))) {
    g$log_weight <- g$log_weight + log_beta_binomial(bb_data$y[i], bb_data$n[i], g$mu, g$kappa)
  }
  normalize_grid(g)
}
bb_fit <- fit_bb(.08); bb_fine <- fit_bb(.04)
bb_check <- check_grid(bb_fit, bb_fine, c("mu", "kappa"),
                      with(bb_fit, a < -3.8 | a > .6 | h < -1.8 | h > 6.8))
show_table(weighted_summary(bb_fit, c("mu", "kappa")))
show_table(bb_check)
write_csv(weighted_summary(bb_fit, c("mu", "kappa")), file.path(out, "beta_binomial_posterior.csv"))
bb_curve <- tibble(y = 0:30) |>
  mutate(Binomial = dbinom(y, 30, .2),
         `Beta-binomial` = exp(log_beta_binomial(y, 30, .2, 10))) |>
  pivot_longer(-y, names_to = "model", values_to = "probability")
ggsave(file.path(out, "beta_binomial.png"),
       ggplot(bb_curve, aes(y, probability, colour = model)) + geom_line() + geom_point() +
         theme_minimal() + labs(x = "Conversions out of 30", title = "Same mean: 6; different dispersion"),
       width = 8, height = 4, dpi = 150)
# END_8
cat("END_8\n")

# BEGIN_10
cat("BEGIN_10\n")
# Problem 1B: different days have different latent service-request rates.
# Gamma uses SHAPE and RATE, never an unspecified second parameter.
mu_nb <- 4
nb_moments <- tibble(phi = c(.5, 2, 10, Inf)) |>
  mutate(mean = mu_nb, variance_rate = mu_nb^2/phi, variance_count = mu_nb + mu_nb^2/phi,
         probability_zero = if_else(is.infinite(phi), exp(-mu_nb),
                                    dnbinom(0, size = phi, mu = mu_nb)))
show_table(nb_moments)
set.seed(1110)
nb_data <- tibble(day = 1:200, lambda_true = rgamma(200, shape = 2, rate = 2/4)) |>
  mutate(y = rpois(n(), lambda_true))
show_table(nb_data |> summarise(mean_count = mean(y), variance_count = var(y), zeros = sum(y == 0)))
# END_10
cat("END_10\n")

# BEGIN_12
cat("BEGIN_12\n")
# Learn the average rate and gamma mixing shape, then predict new days.
fit_nb <- function(step) {
  g <- expand_grid(a = seq(-.5, 3.5, by = step), h = seq(-3, 5, by = step)) |>
    mutate(mu = exp(a), phi = exp(h),
           log_weight = dnorm(a, log(4), 1, log = TRUE) + dnorm(h, log(2), 1, log = TRUE))
  for (y in nb_data$y) g$log_weight <- g$log_weight + dnbinom(y, mu = g$mu, size = g$phi, log = TRUE)
  normalize_grid(g)
}
nb_fit <- fit_nb(.08); nb_fine <- fit_nb(.04)
nb_check <- check_grid(nb_fit, nb_fine, c("mu", "phi"),
                      with(nb_fit, a < -.3 | a > 3.3 | h < -2.8 | h > 4.8))
show_table(weighted_summary(nb_fit, c("mu", "phi")))
show_table(nb_check)
set.seed(1112)
nb_draws <- nb_fit |> slice_sample(n = 3000, replace = TRUE, weight_by = weight)
nb_reps <- map2_dfr(nb_draws$mu, nb_draws$phi, function(mu, phi) {
  y <- rnbinom(nrow(nb_data), mu = mu, size = phi)
  tibble(mean = mean(y), variance = var(y), zero_fraction = mean(y == 0))
})
nb_ppc <- nb_reps |> pivot_longer(everything(), names_to = "statistic") |>
  group_by(statistic) |> summarise(lo90 = quantile(value, .05), hi90 = quantile(value, .95), .groups = "drop") |>
  left_join(tibble(statistic = c("mean", "variance", "zero_fraction"),
                   observed = c(mean(nb_data$y), var(nb_data$y), mean(nb_data$y == 0))), by = "statistic")
show_table(nb_ppc)
write_csv(nb_ppc, file.path(out, "negative_binomial_ppc.csv"))
nb_curve <- tibble(y = 0:30) |>
  mutate(Poisson = dpois(y, 4), `Gamma-Poisson` = dnbinom(y, size = 2, mu = 4)) |>
  pivot_longer(-y, names_to = "model", values_to = "probability")
ggsave(file.path(out, "gamma_poisson.png"),
       ggplot(nb_curve, aes(y, probability, colour = model)) + geom_line() + geom_point() +
         theme_minimal() + labs(x = "Requests per day", title = "Mean 4: variance 4 versus 12"),
       width = 8, height = 4, dpi = 150)
# END_12
cat("END_12\n")

# BEGIN_15
cat("BEGIN_15\n")
# Problem 2: some support desks are offline; active desks may also receive zero calls.
psi_true <- .30; lambda_true <- 2
zip_moments <- tibble(psi = psi_true, lambda = lambda_true) |>
  mutate(mean = (1-psi)*lambda, variance = (1-psi)*lambda*(1+psi*lambda),
         probability_zero = psi+(1-psi)*exp(-lambda),
         offline_given_zero = psi/probability_zero)
show_table(zip_moments)
set.seed(1115)
zip_data <- tibble(day = 1:300, offline_true = rbinom(300, 1, psi_true)) |>
  mutate(y = if_else(offline_true == 1, 0L, rpois(n(), lambda_true)))
show_table(zip_data |> summarise(days = n(), mean_count = mean(y), zero_fraction = mean(y == 0)))
# END_15
cat("END_15\n")

# BEGIN_17
cat("BEGIN_17\n")
# Marginal likelihood sums the two routes to zero; latent labels remain unobserved.
# Priors: logit(psi) ~ N(logit(.3),1), log(lambda) ~ N(log(2),1).
fit_zip <- function(step) {
  g <- expand_grid(a = seq(-5, 2, by = step), b = seq(-2, 3, by = step)) |>
    mutate(psi = plogis(a), lambda = exp(b),
           log_weight = dnorm(a, qlogis(.3), 1, log = TRUE) + dnorm(b, log(2), 1, log = TRUE))
  # Factoring out psi is avoided: log1p is accurate when the mixture is near 1.
  nz <- sum(zip_data$y == 0)
  pos <- zip_data$y[zip_data$y > 0]
  g$log_weight <- g$log_weight + nz * log1p((1-g$psi)*expm1(-g$lambda)) +
    length(pos)*log1p(-g$psi) + sum(pos)*log(g$lambda) - length(pos)*g$lambda - sum(lgamma(pos+1))
  normalize_grid(g) |> mutate(overall_mean = (1-psi)*lambda,
    probability_zero = psi+(1-psi)*exp(-lambda), offline_given_zero = psi/probability_zero)
}
zip_fit <- fit_zip(.05); zip_fine <- fit_zip(.025)
zip_check <- check_grid(zip_fit, zip_fine, c("psi", "lambda", "overall_mean"),
                       with(zip_fit, a < -4.8 | a > 1.8 | b < -1.8 | b > 2.8))
show_table(weighted_summary(zip_fit, c("psi", "lambda", "overall_mean", "offline_given_zero")))
show_table(zip_check)
# A Poisson calibrated to the observed mean has no separate zero-generating mechanism.
show_table(tibble(observed_zero_fraction = mean(zip_data$y == 0),
                 fitted_Poisson_zero_probability = exp(-mean(zip_data$y)),
                 posterior_ZIP_zero_probability = sum(zip_fit$weight * zip_fit$probability_zero)))
write_csv(weighted_summary(zip_fit, c("psi", "lambda", "overall_mean")), file.path(out, "zip_posterior.csv"))
zip_curve <- tibble(y = 0:12) |> mutate(
  Poisson = dpois(y, (1-psi_true)*lambda_true),
  ZIP = if_else(y == 0, psi_true, 0) + (1-psi_true)*dpois(y, lambda_true)) |>
  pivot_longer(-y, names_to = "model", values_to = "probability")
ggsave(file.path(out, "zero_inflation.png"),
       ggplot(zip_curve, aes(y, probability, fill = model)) + geom_col(position = "dodge") +
         theme_minimal() + labs(x = "Calls per day", title = "Both models have mean 1.4"),
       width = 8, height = 4, dpi = 150)
# END_17
cat("END_17\n")

# BEGIN_22
cat("BEGIN_22\n")
# Problem 3: low, medium, high satisfaction are ordered, not equally spaced measurements.
ordinal_prob <- function(eta, cuts) {
  stopifnot(length(cuts) >= 1, all(diff(cuts) > 0))
  diff(c(0, plogis(cuts - eta), 1))
}
ordinal_example <- map_dfr(c(0, 1), function(eta) {
  probabilities <- ordinal_prob(eta, c(-1, 1))
  tibble(eta = eta, category = c("Low", "Medium", "High"),
         probability = probabilities)
})
show_table(ordinal_example)
stopifnot(all(ordinal_example$probability >= 0))
# END_22
cat("END_22\n")

# BEGIN_25
cat("BEGIN_25\n")
# Fit TWO unknown ordered cutpoints and a slope; no free intercept (location identification).
# Ordered cutpoint prior: product N(0,1.5) densities restricted to c1<c2; beta ~ N(0,1).
set.seed(1125)
ordinal_data <- tibble(x = rep(0:1, each = 200)) |>
  mutate(y = map_int(x, function(x) sample(1:3, 1, prob = ordinal_prob(.8*x, c(-1, 1)))))
ordinal_counts <- ordinal_data |> count(x, y) |> complete(x = 0:1, y = 1:3, fill = list(n = 0))
fit_ordinal <- function(step) {
  g <- expand_grid(c1 = seq(-3, .5, by = step), c2 = seq(-.5, 3, by = step), beta = seq(-2, 3, by = step)) |>
    filter(c1 < c2) |>
    mutate(log_weight = dnorm(c1, 0, 1.5, log = TRUE) + dnorm(c2, 0, 1.5, log = TRUE) + dnorm(beta, 0, 1, log = TRUE))
  for (i in seq_len(nrow(ordinal_counts))) {
    cell <- ordinal_counts[i, ]
    f1 <- plogis(g$c1 - g$beta*cell$x); f2 <- plogis(g$c2 - g$beta*cell$x)
    p <- switch(as.character(cell$y), `1` = f1, `2` = f2-f1, `3` = 1-f2)
    g$log_weight <- g$log_weight + cell$n*log(p)
  }
  normalize_grid(g) |> mutate(high_difference = plogis(beta-c2)-plogis(-c2))
}
ordinal_fit <- fit_ordinal(.10); ordinal_fine <- fit_ordinal(.05)
ordinal_check <- check_grid(ordinal_fit, ordinal_fine, c("c1", "c2", "beta", "high_difference"),
  with(ordinal_fit, c1 < -2.8 | c1 > .3 | c2 < -.3 | c2 > 2.8 | beta < -1.8 | beta > 2.8))
show_table(ordinal_counts)
show_table(weighted_summary(ordinal_fit, c("c1", "c2", "beta", "high_difference")))
show_table(ordinal_check)
write_csv(weighted_summary(ordinal_fit, c("c1", "c2", "beta", "high_difference")), file.path(out, "ordinal_posterior.csv"))
ordinal_curve <- expand_grid(eta = seq(-2, 3, length.out = 101), category = 1:3) |>
  mutate(probability = map2_dbl(eta, category, function(e, k) ordinal_prob(e, c(-1,1))[k]))
ggsave(file.path(out, "ordinal.png"),
       ggplot(ordinal_curve, aes(eta, probability, colour = factor(category))) + geom_line() +
         theme_minimal() + labs(colour = "Category", title = "Positive eta shifts mass toward higher categories"),
       width = 8, height = 4, dpi = 150)
# END_25
cat("END_25\n")

# BEGIN_29
cat("BEGIN_29\n")
# Problem 4: ordered training levels with unequal increments.
delta <- c(.10, .60, .30); beta <- 1.5
monotonic_example <- tibble(level = 1:4, score = c(0, cumsum(delta))) |>
  mutate(eta = beta*score, probability_high = plogis(eta-1))
show_table(monotonic_example)
# Symmetric Dirichlet: expected increments are equal, individual vectors need not be.
set.seed(1129)
dirichlet_moments <- map_dfr(c(.5, 1, 2, 10), function(a) {
  z <- matrix(rgamma(30000, shape = a), ncol = 3)
  z <- z / rowSums(z)
  tibble(a = a, theoretical_mean = 1/3, theoretical_variance = 2/(9*(3*a+1)),
         simulated_mean = mean(z[,1]), simulated_variance = var(z[,1]))
})
show_table(dirichlet_moments)
# END_29
cat("END_29\n")

# BEGIN_31
cat("BEGIN_31\n")
# Conditional learning exercise: cutpoints are FIXED at (-1,1), as in an externally
# calibrated instrument. Learn beta and the three simplex increments, not cutpoints.
set.seed(1131)
mono_data <- tibble(level = rep(1:4, each = 150)) |>
  mutate(y = map_int(level, function(e) sample(1:3, 1, prob = ordinal_prob(monotonic_example$eta[e], c(-1,1)))))
mono_counts <- mono_data |> count(level, y) |> complete(level = 1:4, y = 1:3, fill = list(n = 0))
fit_monotonic <- function(M, step) {
  simplex <- expand_grid(j = seq_len(M-2), k = seq_len(M-2)) |> filter(j+k < M) |>
    transmute(d1 = j/M, d2 = k/M, d3 = 1-d1-d2)
  g <- crossing(simplex, beta = seq(-2, 4, by = step)) |>
    mutate(log_weight = dnorm(beta, 0, 1, log = TRUE) + log(d1)+log(d2)+log(d3))
  # Dirichlet(2,2,2) density is proportional to d1*d2*d3 on the simplex.
  for (i in seq_len(nrow(mono_counts))) {
    cell <- mono_counts[i, ]
    score <- switch(as.character(cell$level), `1` = 0, `2` = g$d1, `3` = g$d1+g$d2, `4` = 1)
    f1 <- plogis(-1-g$beta*score); f2 <- plogis(1-g$beta*score)
    p <- switch(as.character(cell$y), `1` = f1, `2` = f2-f1, `3` = 1-f2)
    g$log_weight <- g$log_weight + cell$n*log(p)
  }
  normalize_grid(g)
}
mono_fit <- fit_monotonic(30, .1); mono_fine <- fit_monotonic(60, .05)
# Simplex boundaries are genuine parameter boundaries, not artificial truncation.
mono_check <- check_grid(mono_fit, mono_fine, c("beta", "d1", "d2", "d3"),
                        with(mono_fit, beta < -1.8 | beta > 3.8))
show_table(weighted_summary(mono_fit, c("beta", "d1", "d2", "d3")))
show_table(mono_check)
mono_curve <- bind_rows(monotonic_example |> mutate(model = "Unequal increments"),
                       tibble(level = 1:4, score = (0:3)/3, eta = 1.5*(0:3)/3,
                              probability_high = plogis(eta-1), model = "Equal increments"))
ggsave(file.path(out, "monotonic.png"),
       ggplot(mono_curve, aes(level, probability_high, colour = model)) + geom_line() + geom_point() +
         theme_minimal() + labs(y = "Probability of high satisfaction", title = "Same endpoint effect, different intermediate levels"),
       width = 8, height = 4, dpi = 150)
write_csv(weighted_summary(mono_fit, c("beta", "d1", "d2", "d3")), file.path(out, "monotonic_posterior.csv"))
# END_31
cat("END_31\n")

# BEGIN_34
cat("BEGIN_34\n")
# Independent probability and moment checks, beyond posterior grid refinement.
stopifnot(abs(sum(exp(log_beta_binomial(0:30, 30, .2, 10))) - 1) < 1e-10)
y <- 0:30; p <- exp(log_beta_binomial(y, 30, .2, 10))
stopifnot(abs(sum(y*p)-6) < 1e-10, abs(sum((y-6)^2*p)-4.8*40/11) < 1e-9)
stopifnot(abs(sum(dnbinom(0:1000, size=2, mu=4))-1) < 1e-10)
for (eta in c(-10, 0, 10)) stopifnot(abs(sum(ordinal_prob(eta, c(-1,1)))-1) < 1e-12)
# A shift of BOTH cutpoints and eta must leave category probabilities unchanged.
stopifnot(max(abs(ordinal_prob(1,c(-1,1))-ordinal_prob(4,c(2,4)))) < 1e-12)
stopifnot(all(diff(monotonic_example$eta) >= 0), abs(sum(delta)-1) < 1e-12)
checks <- bind_rows(Beta_binomial = bb_check, Gamma_Poisson = nb_check, ZIP = zip_check,
                    Ordinal = ordinal_check, Monotonic_conditional = mono_check, .id = "model")
show_table(checks)
write_csv(checks, file.path(out, "numerical_checks.csv"))
cat("PASS: normalization, analytic moments, ordinal identification, monotonicity and grid checks.\n")
cat("No HMC, full SBC or empirical chapter-data fits were run.\n")
# END_34
cat("END_34\n")
