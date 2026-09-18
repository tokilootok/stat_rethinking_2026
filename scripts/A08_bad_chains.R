#!/usr/bin/env Rscript

# Deliberately broken examples for Lecture A08. They are tests of the
# diagnostic workflow, not templates for substantive analysis.

suppressPackageStartupMessages({
  library(rethinking)
  library(cmdstanr)
})

set.seed(202608)
dir.create("artifacts/A08", recursive = TRUE, showWarnings = FALSE)

# Experiment 1: weak priors create a badly conditioned posterior.
y <- c(-1, 1)
bad_weak_prior <- ulam(
  alist(
    y ~ dnorm(mu, sigma),
    mu <- alpha,
    alpha ~ dnorm(0, 1000),
    sigma ~ dexp(0.0001)
  ),
  data = list(y = y),
  chains = 4,
  cores = 4,
  iter = 1000,
  log_lik = FALSE
)

good_regularized <- ulam(
  alist(
    y ~ dnorm(mu, sigma),
    mu <- alpha,
    alpha ~ dnorm(0, 1),
    sigma ~ dexp(1)
  ),
  data = list(y = y),
  chains = 4,
  cores = 4,
  iter = 1000,
  log_lik = FALSE
)

print(precis(bad_weak_prior))
print(precis(good_regularized))

pdf("artifacts/A08/bad-versus-good-chains.pdf", width = 10, height = 8)
par(mfrow = c(2, 2))
traceplot(bad_weak_prior)
trankplot(bad_weak_prior)
traceplot(good_regularized)
trankplot(good_regularized)
dev.off()

# Experiment 2 is specified in the chapter: allowing signed discrimination
# produces two observationally equivalent IRT orientations. Its failure is
# structural label switching and must be repaired by identification, not by
# drawing more samples or increasing adapt_delta.

