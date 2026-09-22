#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
})

args <- commandArgs(trailingOnly = TRUE)
csv_path <- if (length(args) >= 1L) args[[1L]] else "artifacts/A09/nwo_grants.csv"
stan_path <- if (length(args) >= 2L) args[[2L]] else "scripts/A09_grants.stan"
output_dir <- if (length(args) >= 3L) args[[3L]] else "artifacts/A09"
quick <- identical(Sys.getenv("A09_QUICK"), "true")

d <- read.csv(csv_path, stringsAsFactors = FALSE)
ND <- max(d$discipline_id)

weight_pooled <- as.numeric(tapply(d$applications, d$discipline_id, sum))
weight_pooled <- weight_pooled / sum(weight_pooled)
weight_female <- as.numeric(tapply(
  d$applications[d$gender_id == 1], d$discipline_id[d$gender_id == 1], sum
))
weight_female <- weight_female / sum(weight_female)
weight_male <- as.numeric(tapply(
  d$applications[d$gender_id == 2], d$discipline_id[d$gender_id == 2], sum
))
weight_male <- weight_male / sum(weight_male)

stan_data <- list(
  K = nrow(d),
  ND = ND,
  awards = d$awards,
  applications = d$applications,
  gender_id = d$gender_id,
  discipline_id = d$discipline_id,
  weight_pooled = weight_pooled,
  weight_female = weight_female,
  weight_male = weight_male
)

dir.create(file.path(output_dir, "compiled"), recursive = TRUE, showWarnings = FALSE)
model <- cmdstan_model(
  stan_path,
  exe_file = file.path(output_dir, "compiled", "A09_grants")
)
fit <- model$sample(
  data = stan_data,
  seed = 202609,
  chains = if (quick) 2 else 4,
  parallel_chains = if (quick) 2 else 4,
  iter_warmup = if (quick) 250 else 1000,
  iter_sampling = if (quick) 250 else 1000,
  adapt_delta = 0.95,
  refresh = if (quick) 0 else 250
)

key <- c(
  "total_female",
  "total_male",
  "total_risk_difference",
  "standardized_female",
  "standardized_male",
  "standardized_direct_risk_difference"
)
print(fit$summary(variables = key))

diagnostics <- fit$diagnostic_summary()
summary_all <- fit$summary()
print(diagnostics)
stopifnot(
  all(is.finite(summary_all$rhat)),
  max(summary_all$rhat, na.rm = TRUE) < if (quick) 1.10 else 1.01,
  sum(diagnostics$num_divergent) == 0L,
  sum(diagnostics$num_max_treedepth) == 0L
)

dir.create(file.path(output_dir, "cmdstan"), recursive = TRUE, showWarnings = FALSE)
fit$save_output_files(dir = file.path(output_dir, "cmdstan"))
write.csv(fit$summary(), file.path(output_dir, "cmdstan_summary.csv"), row.names = FALSE)

