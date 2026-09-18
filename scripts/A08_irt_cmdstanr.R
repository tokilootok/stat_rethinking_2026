#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
})

args <- commandArgs(trailingOnly = TRUE)
csv_path <- if (length(args) >= 1L) args[[1L]] else "artifacts/A08/wines2012.csv"
stan_path <- if (length(args) >= 2L) args[[2L]] else "scripts/A08_irt.stan"
output_dir <- if (length(args) >= 3L) args[[3L]] else "artifacts/A08"
quick <- identical(Sys.getenv("A08_QUICK"), "true")

d <- read.csv(csv_path, stringsAsFactors = FALSE)

stopifnot(
  nrow(d) == 180L,
  identical(sort(unique(d$judge_id)), seq_len(max(d$judge_id))),
  identical(sort(unique(d$wine_id)), seq_len(max(d$wine_id)))
)

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

dir.create(file.path(output_dir, "compiled"), recursive = TRUE, showWarnings = FALSE)

stan_data <- list(
  N = nrow(d),
  NW = max(d$wine_id),
  NJ = max(d$judge_id),
  score_z = d$score_z,
  wine_id = d$wine_id,
  judge_id = d$judge_id,
  wine_origin_c = wine_origin_c,
  judge_country_c = judge_country_c
)

model <- cmdstan_model(
  stan_path,
  exe_file = file.path(output_dir, "compiled", "A08_irt")
)
fit <- model$sample(
  data = stan_data,
  seed = 202608,
  chains = if (quick) 2 else 4,
  parallel_chains = if (quick) 2 else 4,
  iter_warmup = if (quick) 250 else 1000,
  iter_sampling = if (quick) 250 else 1000,
  adapt_delta = 0.99,
  refresh = if (quick) 0 else 250
)

key <- c(
  "delta_harshness",
  "delta_log_discrimination",
  "discrimination_ratio_us_fr",
  "delta_origin",
  "sigma_score"
)
print(fit$summary(variables = key))

diagnostics <- fit$diagnostic_summary()
print(diagnostics)

summary_all <- fit$summary()
stopifnot(
  all(is.finite(summary_all$rhat)),
  max(summary_all$rhat, na.rm = TRUE) < if (quick) 1.10 else 1.01,
  sum(diagnostics$num_divergent) == 0L,
  sum(diagnostics$num_max_treedepth) == 0L
)

dir.create(file.path(output_dir, "cmdstan"), recursive = TRUE, showWarnings = FALSE)
fit$save_output_files(dir = file.path(output_dir, "cmdstan"))
write.csv(fit$summary(), file.path(output_dir, "cmdstan_summary.csv"), row.names = FALSE)

