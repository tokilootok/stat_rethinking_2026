#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(rethinking)
  library(cmdstanr)
})

args <- commandArgs(trailingOnly = TRUE)
csv_path <- if (length(args) >= 1L) args[[1L]] else "artifacts/A09/nwo_grants.csv"
quick <- identical(Sys.getenv("A09_QUICK"), "true")
d <- read.csv(csv_path, stringsAsFactors = FALSE)

dat <- list(
  A = d$awards,
  N = d$applications,
  G = d$gender_id,
  D = d$discipline_id
)

m_grants <- ulam(
  alist(
    A ~ dbinom(N, p),
    logit(p) <- a[G, D],
    matrix[2, 9]:a ~ dnorm(0, 1)
  ),
  data = dat,
  chains = if (quick) 2 else 4,
  cores = if (quick) 2 else 4,
  iter = if (quick) 500 else 2000,
  warmup = if (quick) 250 else 1000,
  control = list(adapt_delta = 0.95),
  log_lik = TRUE,
  cmdstan = TRUE
)

fit <- attr(m_grants, "cstanfit")
diagnostics <- fit$diagnostic_summary()
summary_all <- fit$summary()
print(diagnostics)
stopifnot(
  all(is.finite(summary_all$rhat)),
  max(summary_all$rhat, na.rm = TRUE) < if (quick) 1.10 else 1.01,
  sum(diagnostics$num_divergent) == 0L,
  sum(diagnostics$num_max_treedepth) == 0L
)

post <- extract.samples(m_grants)
p <- inv_logit(post$a)

apps_by_discipline <- tapply(d$applications, d$discipline_id, sum)
w_pool <- apps_by_discipline / sum(apps_by_discipline)
w_f <- as.numeric(tapply(
  d$applications[d$gender_id == 1], d$discipline_id[d$gender_id == 1], sum
))
w_f <- w_f / sum(w_f)
w_m <- as.numeric(tapply(
  d$applications[d$gender_id == 2], d$discipline_id[d$gender_id == 2], sum
))
w_m <- w_m / sum(w_m)

total_f <- as.vector(p[, 1, ] %*% w_f)
total_m <- as.vector(p[, 2, ] %*% w_m)
direct_f <- as.vector(p[, 1, ] %*% w_pool)
direct_m <- as.vector(p[, 2, ] %*% w_pool)

estimands <- data.frame(
  total_risk_difference = total_f - total_m,
  standardized_direct_risk_difference = direct_f - direct_m
)

print(t(vapply(estimands, function(x) {
  c(mean = mean(x), PI(x, prob = 0.89), Pr_gt_0 = mean(x > 0))
}, numeric(4))))

