#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(rethinking)
  library(cmdstanr)
})

args <- commandArgs(trailingOnly = TRUE)
csv_path <- if (length(args) >= 1L) args[[1L]] else "artifacts/A08/wines2012.csv"
quick <- identical(Sys.getenv("A08_QUICK"), "true")

d <- read.csv(csv_path, stringsAsFactors = FALSE)
judge_country_c <- vapply(
  seq_len(max(d$judge_id)),
  function(j) unique(d$judge_country_c[d$judge_id == j]),
  numeric(1)
)
wine_origin_c <- vapply(
  seq_len(max(d$wine_id)),
  function(w) unique(d$wine_origin_c[d$wine_id == w]),
  numeric(1)
)

dat <- list(
  S = d$score_z,
  W = d$wine_id,
  J = d$judge_id,
  A = judge_country_c,
  XW = wine_origin_c,
  NW = max(d$wine_id),
  NJ = max(d$judge_id),
  N = nrow(d)
)

f_irt <- alist(
  S ~ dnorm(mu, sigma_score),
  vector[N]:mu <- exp(a_D + b_D * A[J] + sigma_D * z_D[J]) *
        (q[W] - (a_H + b_H * A[J] + sigma_H * z_H[J])),

  vector[NW]:q ~ dnorm(b_O * XW, 1),
  vector[NJ]:z_H ~ dnorm(0, 1),
  vector[NJ]:z_D ~ dnorm(0, 1),

  a_H ~ dnorm(0, 0.5),
  b_H ~ dnorm(0, 0.5),
  sigma_H ~ dexp(2),

  a_D ~ dnorm(0, 0.35),
  b_D ~ dnorm(0, 0.35),
  sigma_D ~ dexp(2),

  b_O ~ dnorm(0, 0.5),
  sigma_score ~ dexp(1)
)

m_irt <- ulam(
  f_irt,
  data = dat,
  chains = if (quick) 2 else 4,
  cores = if (quick) 2 else 4,
  iter = if (quick) 500 else 2000,
  warmup = if (quick) 250 else 1000,
  control = list(adapt_delta = 0.99),
  log_lik = FALSE,
  cmdstan = TRUE
)

precis_irt <- precis(
  m_irt,
  depth = 2,
  pars = c("b_H", "b_D", "b_O", "sigma_H", "sigma_D", "sigma_score")
)
print(precis_irt)

diagnostics <- attr(m_irt, "cstanfit")$diagnostic_summary()
summary_all <- attr(m_irt, "cstanfit")$summary()
print(diagnostics)
stopifnot(
  all(is.finite(summary_all$rhat)),
  max(summary_all$rhat, na.rm = TRUE) < if (quick) 1.10 else 1.01,
  sum(diagnostics$num_divergent) == 0L,
  sum(diagnostics$num_max_treedepth) == 0L
)

post <- extract.samples(m_irt)
contrast <- data.frame(
  delta_harshness_us_fr = post$b_H,
  discrimination_ratio_us_fr = exp(post$b_D),
  delta_origin_us_fr = post$b_O
)

print(rbind(
  delta_harshness_us_fr = c(mean = mean(post$b_H), PI(post$b_H, prob = 0.89)),
  discrimination_ratio_us_fr = c(mean = mean(exp(post$b_D)), PI(exp(post$b_D), prob = 0.89)),
  delta_origin_us_fr = c(mean = mean(post$b_O), PI(post$b_O, prob = 0.89))
))
cat(sprintf("Pr(US judges harsher) = %.3f\n", mean(post$b_H > 0)))
cat(sprintf("Pr(US judges more discriminating) = %.3f\n", mean(post$b_D > 0)))
cat(sprintf("Pr(US wines higher latent quality) = %.3f\n", mean(post$b_O > 0)))

