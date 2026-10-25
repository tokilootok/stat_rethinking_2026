#!/usr/bin/env Rscript
# A10 — tidyverse + ulam companion. Run from the repository root:
#   Rscript scripts/A10_walkthroughs_tidyverse_ulam.R
# Syntax/short-chain smoke test (NOT final inference): append --quick.
# Dependencies: rethinking, cmdstanr, posterior, digest, dplyr, tidyr, purrr, tibble,
#               readr, ggplot2, plus an installed CmdStan C++ toolchain.
# Validated compact checkpoints are reused automatically; remove the checkpoints
# directory to request fresh fits. Insufficient mixing triggers longer-chain retries.
# Original grid route: scripts/A10_walkthroughs.R (retained for numerical checks).
# Outputs: courses/figures/A10_tidyverse_ulam/{full,quick}/.
# No real-admissions latent fit, full SBC, or causal-identification claim is made.

required <- c("rethinking", "cmdstanr", "posterior", "digest", "dplyr", "tidyr",
              "purrr", "tibble", "readr", "ggplot2")
missing_packages <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Install missing packages: ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(tibble)
  library(readr); library(ggplot2)
})
cmdstanr::check_cmdstan_toolchain(fix = FALSE, quiet = TRUE)
if (!file.exists("scripts/A10_walkthroughs.R")) stop("Run from the repository root.")
arguments <- commandArgs(trailingOnly = TRUE)
if (any(!arguments %in% "--quick")) stop("Only optional argument: --quick")
quick <- "--quick" %in% arguments
config <- list(chains = if (quick) 2L else 4L,
               warmup = if (quick) 250L else 1000L,
               sampling = if (quick) 250L else 1000L,
               adapt_delta = .95, max_treedepth = 12L)
out <- file.path("courses/figures/A10_tidyverse_ulam", if (quick) "quick" else "full")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
options(width = 110, digits = 4)
show_table <- function(x) print(as.data.frame(x), row.names = FALSE, digits = 4)
write_table <- function(x, filename) write_csv(x, file.path(out, filename))
summary90 <- function(tbl) {
  tbl |> summarise(mean = mean(value), lo90 = quantile(value, .05),
                   hi90 = quantile(value, .95), prob_positive = mean(value > 0), .groups = "drop")
}
# Preserve chain/iteration IDs: diagnostics must not treat shuffled draws as chains.
core_draws <- function(model, variables) {
  parameter_names <- names(model$core)[vapply(names(model$core), function(nm)
    any(nm == variables | startsWith(nm, paste0(variables, "["))), logical(1))]
  model$core |> select(.chain, .iteration, .draw, all_of(parameter_names))
}
contrast_diagnostics <- function(draws) {
  draws |> select(.chain, .iteration, .draw, department, value) |>
    pivot_wider(names_from = department, names_prefix = "department_", values_from = value) |>
    as.data.frame() |> posterior::as_draws_df() |>
    posterior::summarise_draws("mean", "sd", "rhat", "ess_bulk", "ess_tail", "mcse_mean") |>
    as_tibble() |> mutate(department = as.integer(sub("department_", "", variable)))
}
all_diagnostics <- list()
# Compact checkpoints retain core draws and diagnostics, not huge individual-u arrays.
# The fingerprint prevents reuse when model/data/settings/toolchain differ.
cache_dir <- file.path(out, "checkpoints")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
ignore_path <- file.path(dirname(out), ".gitignore")
ignore_rules <- if (file.exists(ignore_path)) readLines(ignore_path) else character()
if (!"checkpoints/" %in% ignore_rules) writeLines(c(ignore_rules, "checkpoints/"), ignore_path)
fit_ulam <- function(name, formula, data, seed, constraints = NULL) {
  fingerprint <- digest::digest(list(formula, data, config, seed, constraints,
    rethinking = as.character(packageVersion("rethinking")),
    posterior = as.character(packageVersion("posterior")),
    cmdstan = cmdstanr::cmdstan_version(), validation_version = 1L))
  cache_file <- file.path(cache_dir, paste0(name, ".rds"))
  if (file.exists(cache_file)) {
    cached <- readRDS(cache_file)
    if (identical(cached$fingerprint, fingerprint) && (quick || isTRUE(cached$report$pass))) {
      cat("Reusing matching", if (quick) "smoke-test" else "validated", "checkpoint:", name, "\n")
      all_diagnostics[[name]] <<- cached$report
      return(cached)
    }
  }
  attempts <- if (quick) 1L else 3L
  for (attempt in seq_len(attempts)) {
    warmup <- config$warmup * if (attempt == 1L) 1L else 2L
    sampling <- config$sampling * c(1L,4L,8L)[attempt]
    cat("\nFitting", name, "attempt", attempt, "with", config$chains, "chains,",
        warmup, "warmup and", sampling, "retained per chain; mode:",
        if (quick) "SMOKE TEST" else "FULL", "\n")
    fit_arguments <- list(
      flist = formula, data = data, chains = config$chains, cores = config$chains,
      iter = warmup + sampling, warmup = warmup,
      control = list(adapt_delta = config$adapt_delta), max_treedepth = config$max_treedepth,
      cmdstan = TRUE, log_lik = FALSE, seed = seed, refresh = 0, messages = FALSE
    )
    if (!is.null(constraints)) fit_arguments$constraints <- constraints
    model <- do.call(rethinking::ulam, fit_arguments)
    fit <- attr(model, "cstanfit")
    # Include all sampled parameters, including person-specific latent variables.
    diagnostics <- fit$summary() |> as_tibble()
    chains <- fit$diagnostic_summary() |> as_tibble()
    write_table(diagnostics, paste0(name, "_parameter_diagnostics.csv"))
    write_table(chains, paste0(name, "_chain_diagnostics.csv"))
    write_table(diagnostics, paste0(name, "_attempt_", attempt, "_parameter_diagnostics.csv"))
    write_table(chains, paste0(name, "_attempt_", attempt, "_chain_diagnostics.csv"))
    invisible(capture.output(stan_code <- rethinking::stancode(model)))
    writeLines(stan_code, file.path(out, paste0(name, ".stan")))
    report <- tibble(model = name, attempt = attempt, warmup = warmup, sampling = sampling,
      max_rhat = max(diagnostics$rhat, na.rm = TRUE),
      min_ess_bulk = min(diagnostics$ess_bulk, na.rm = TRUE),
      min_ess_tail = min(diagnostics$ess_tail, na.rm = TRUE),
      divergences = sum(chains$num_divergent), treedepth_hits = sum(chains$num_max_treedepth),
      min_ebfmi = min(chains$ebfmi),
      finite_diagnostics = all(is.finite(diagnostics$rhat)) &&
        all(is.finite(diagnostics$ess_bulk)) && all(is.finite(diagnostics$ess_tail))) |>
      mutate(pass = finite_diagnostics & max_rhat < 1.01 & min_ess_bulk >= 400 &
               min_ess_tail >= 400 & divergences == 0 & treedepth_hits == 0 & min_ebfmi >= .3)
    show_table(report)
    if (quick || report$pass || attempt == attempts) break
    # More draws address mixing/precision only. Geometry warnings need model investigation.
    if (report$divergences > 0 || report$treedepth_hits > 0 || report$min_ebfmi < .3) {
      warning(name, ": geometry warning; inspect diagnostics instead of retrying blindly.")
      break
    }
    cat("Agreement/precision insufficient; retrying with longer chains. Failed attempt is retained.\n")
  }
  all_diagnostics[[name]] <<- report
  core_names <- diagnostics$variable[grepl("^(a|delta|b|g)(\\[|$)", diagnostics$variable)]
  result <- list(core = fit$draws(variables = core_names, format = "draws_df") |> as_tibble(),
                 parameter_diagnostics = diagnostics, report = report, fingerprint = fingerprint)
  if (quick || report$pass) saveRDS(result, cache_file)
  result
}

# Section 6 — One known simulation; G labels are synthetic categories, not real genders.
cat("BEGIN_6\n")
set.seed(12)
N <- 2000L
sim <- tibble(id = seq_len(N), G = sample(1:2, N, replace = TRUE), U = rbinom(N, 1, .1)) |>
  mutate(D = 1L + rbinom(N, 1, if_else(G == 1L, as.double(U), .75)))
# The lookup keys are explicit; the RNG order matches the original walkthrough.
probability_lookup <- tribble(
  ~D, ~G, ~p_U0, ~p_U1,
  1L, 1L, .1, .3, 2L, 1L, .1, .3,
  1L, 2L, .1, .5, 2L, 2L, .3, .5)
sim <- sim |> left_join(probability_lookup, by = c("D", "G")) |>
  mutate(p_true = if_else(U == 1, p_U1, p_U0), A = rbinom(N, 1, p_true)) |>
  select(id, G, U, D, A, p_true)
counts <- sim |> group_by(G, D) |>
  summarise(admitted = sum(A), applied = n(), .groups = "drop") |>
  mutate(rate = admitted/applied)
stopifnot(nrow(counts) == 4, sum(counts$applied) == N, !anyNA(sim))
show_table(counts)
truth <- tibble(department = 1:2, observed_G1 = c(.1,.3), observed_G2 = c(.14,.32),
                intervention_G1 = .12, intervention_G2 = c(.14,.32)) |>
  mutate(observed_difference = observed_G1-observed_G2,
         controlled_difference = intervention_G1-intervention_G2)
show_table(truth)
cat("Known total contrast G1-G2:", .12-(.25*.14+.75*.32), "\n")
write_table(sim, "simulated_applicants.csv")
write_table(counts, "simulated_counts.csv")
# U and p_true are deliberately absent from the inference data.
data_observed <- list(A = as.integer(sim$A), G = as.integer(sim$G), D = as.integer(sim$D))
cat("END_6\n")

# Section 8 — Aggregate observed disparity: an ordinary Bernoulli ulam fit.
cat("BEGIN_8\n")
aggregate_fit <- fit_ulam("aggregate", alist(
  A ~ bernoulli(p), logit(p) <- a[G], a[G] ~ normal(0,1)
), data_observed[c("A", "G")], 1008)
aggregate_draws <- core_draws(aggregate_fit, "a") |>
  transmute(.chain, .iteration, .draw, department = 0L,
            value = plogis(`a[1]`)-plogis(`a[2]`))
show_table(aggregate_draws |> mutate(value = 100*value) |> summary90())
cat("Contrast direction throughout: G1 minus G2; printed intervals are central 90%.\nEND_8\n")

# Section 9 — Within-department associations, without correcting latent selection.
cat("BEGIN_9\n")
cell_fit <- fit_ulam("observed_cells", alist(
  A ~ bernoulli(p), logit(p) <- a[G,D], matrix[2,2]:a ~ normal(0,1)
), data_observed, 1009)
# ulam stores a[G,D]; convert named columns to explicit group/department keys.
cell_probabilities <- function(draws) {
  draws |> select(.chain, .iteration, .draw, starts_with("a[")) |>
    pivot_longer(starts_with("a["), names_to = "parameter", values_to = "a") |>
    extract(parameter, into = c("G", "department"), regex = "a\\[(\\d+),(\\d+)\\]", convert = TRUE)
}
cell_draws <- cell_probabilities(core_draws(cell_fit, "a")) |>
  mutate(p = plogis(a)) |> select(-a) |>
  pivot_wider(names_from = G, names_prefix = "G", values_from = p) |>
  mutate(value = G1-G2)
show_table(cell_draws |> mutate(value = 100*value) |> group_by(department) |> summary90())
cat("END_9\n")

# Section 10 — Seeing U would not create the missing combinations in this simulation.
cat("BEGIN_10\n")
show_table(sim |> count(G, D, U) |> complete(G = 1:2, D = 1:2, U = 0:1, fill = list(n = 0)))
cat("G1 has no U=1 in D1 and no U=0 in D2. Observing U alone does not restore overlap.\n")
cat("A common additive observed-U slope is also not the exact binary simulator.\nEND_10\n")

# Sections 12 and 14 — Interpret coefficients before assigning sensitivity values.
cat("BEGIN_12\n")
show_table(tibble(u = c(-1,0,1)) |>
  mutate(admission = plogis(qlogis(.1)+u), dept2_G1 = plogis(qlogis(.25)+u), dept2_G2 = .25))
cat("One unit of a logit coefficient multiplies odds by", exp(1), "\nEND_12\nBEGIN_14\n")
show_table(tibble(prior_quantile = c(.05,.5,.95)) |>
  mutate(coefficient = prior_quantile, odds_ratio = exp(coefficient),
         probability_from_10pct = plogis(qlogis(.1)+coefficient)))
cat("Uniform(0,1) coefficients restrict odds multipliers to [1,e], not all real effects.\nEND_14\n")

# Section 17 — Both observed processes share ONE latent u for EACH applicant.
cat("BEGIN_17\n")
data_joint <- c(data_observed, list(D2 = as.integer(sim$D == 2L), N = N))
fixed_formula <- alist(
  A ~ bernoulli(p),
  logit(p) <- a[G,D] + b[G]*u[i],
  matrix[2,2]:a ~ normal(0,1),
  D2 ~ bernoulli(q),
  logit(q) <- delta[G] + g[G]*u[i],
  delta[G] ~ normal(0,1),
  vector[N]:u ~ normal(0,1)
)
# Different from the binary generating U: these are explicitly Gaussian sensitivity models.
# Population targets integrate over a FRESH u ~ N(0,1), not fitted individual u[i].
# Selection loadings g affect the fitted a,b through the joint likelihood, but do
# not enter prediction directly once we hold the department fixed.
population_contrasts <- function(draws, b_fixed = NULL, Q = 601L) {
  coefficients <- if (is.null(b_fixed)) {
    draws |> select(.chain, .iteration, .draw, starts_with("b[")) |>
      pivot_longer(starts_with("b["), names_to = "G", values_to = "b") |>
      mutate(G = as.integer(sub("b\\[(\\d+)\\]", "\\1", G)))
  } else {
    crossing(draws |> select(.chain, .iteration, .draw), G = 1:2) |>
      mutate(b = b_fixed[G])
  }
  nodes <- qnorm((seq_len(Q)-.5)/Q)
  cells <- cell_probabilities(draws) |>
    left_join(coefficients, by = c(".chain", ".iteration", ".draw", "G")) |>
    mutate(p_population = map2_dbl(a, b, function(a,b) mean(plogis(a+b*nodes))))
  stopifnot(!anyNA(cells), nrow(cells) == 4*nrow(draws))
  cells |> select(.chain, .iteration, .draw, department, G, p_population) |>
    pivot_wider(names_from = G, names_prefix = "G", values_from = p_population) |>
    mutate(value = G1-G2)
}
cat("Normal(0,1) priors on a, delta and u; fixed b and g in each scenario.\n")
cat("The reported contrast holds department fixed and standardizes both groups to N(0,1).\nEND_17\n")

# Common plotting helper: summaries stay on the probability scale in CSV files.
plot_scenarios <- function(summary, levels) {
  summary |> mutate(across(c(mean,lo90,hi90), ~100*.x),
                    scenario = factor(scenario, levels = levels)) |>
    ggplot(aes(scenario, mean, ymin = lo90, ymax = hi90)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_pointrange(colour = "steelblue") + coord_flip() +
    facet_wrap(~department, labeller = label_both) + theme_minimal() +
    labs(x = NULL, y = "G1 minus G2 (percentage points)",
         title = "A10: sensitivity to latent selection",
         subtitle = if (quick) "SMOKE TEST: not final inference" else
           "ulam posterior means and central 90% intervals")
}

# Section 19 — Refit all five scenarios; changing coefficients without refitting is incorrect.
cat("BEGIN_19\n")
settings <- tribble(
  ~scenario, ~b1, ~b2, ~g1, ~g2,
  "No_latent", 0, 0, 0, 0,
  "Outcome_only", 1, 1, 0, 0,
  "Selection_only", 0, 0, 1, 0,
  "Official_fixed", 1, 1, 1, 0,
  "Stronger", 2, 2, 2, 0)
fixed_models <- set_names(map(seq_len(nrow(settings)), function(k) {
  s <- settings[k, ]
  fit_ulam(s$scenario, fixed_formula,
           c(data_joint, list(b = c(s$b1,s$b2), g = c(s$g1,s$g2))), 10190+k)
}), settings$scenario)
fixed_core <- map(fixed_models, core_draws, variables = c("a", "delta"))
fixed_contrasts <- imap_dfr(fixed_core, function(draws, name) {
  s <- filter(settings, scenario == name)
  population_contrasts(draws, c(s$b1,s$b2)) |> mutate(scenario = name)
})
scenario_summary <- fixed_contrasts |> group_by(scenario, department) |> summary90()
show_table(scenario_summary |> mutate(across(c(mean,lo90,hi90), function(x) 100*x)))
write_table(scenario_summary, "fixed_scenario_summary.csv")
# Explicit ggplot equivalent of the base-R two-panel loop in lesson 19.
fixed_plot_data <- scenario_summary |>
  mutate(across(c(mean, lo90, hi90), ~100*.x), # probability -> percentage points
         scenario = factor(scenario, levels = settings$scenario))
fixed_scenario_plot <- ggplot(fixed_plot_data, aes(x = mean, y = scenario)) +
  geom_segment(aes(x = lo90, xend = hi90, yend = scenario),
               colour = "steelblue", linewidth = .7) + # segments(): 90% intervals
  geom_point(size = 2) +                               # plot(..., pch = 19)
  geom_vline(xintercept = 0, linetype = "dashed") +      # abline(v = 0)
  facet_wrap(~department, nrow = 1, scales = "free_x",
             labeller = label_both) +                 # par(mfrow = c(1,2)) + loop
  scale_y_discrete(drop = FALSE) +                     # axis(): scenario labels
  theme_minimal() +                                    # automatic label margins
  labs(x = "G1 minus G2 (percentage points)", y = NULL,
       title = "A10: fixed sensitivity scenarios",
       subtitle = if (quick) "SMOKE TEST: not final inference" else
         "ulam posterior means and central 90% intervals")
if (interactive()) print(fixed_scenario_plot)
ggsave(file.path(out, "fixed_scenarios.png"), fixed_scenario_plot,
       width = 11, height = 5, dpi = 150)
# Check prediction integration at the SAME posterior parameter draws.
official_fine <- population_contrasts(fixed_core$Official_fixed, c(1,1), Q = 1201L)
integration_error <- max(abs(official_fine$value - filter(fixed_contrasts, scenario == "Official_fixed")$value))
stopifnot(integration_error < .001)
cat("Largest contrast change after doubling latent quadrature nodes:", integration_error, "\nEND_19\n")

# Section 20 — Joint posterior over sensitivity coefficients, not equal weighting of fixed scenarios.
cat("BEGIN_20\n")
prior_formula <- alist(
  A ~ bernoulli(p),
  logit(p) <- a[G,D] + b[G]*u[i],
  matrix[2,2]:a ~ normal(0,1),
  vector[2]:b ~ uniform(0,1),
  D2 ~ bernoulli(q),
  logit(q) <- delta[G] + g[G]*u[i],
  delta[G] ~ normal(0,1),
  vector[2]:g ~ uniform(0,1),
  vector[N]:u ~ normal(0,1)
)
prior_model <- fit_ulam("Prior_weighted", prior_formula, data_joint, 1020,
                        constraints = list(b = "lower=0,upper=1", g = "lower=0,upper=1"))
prior_core <- core_draws(prior_model, c("a", "delta", "b", "g"))
prior_contrasts <- population_contrasts(prior_core) |> mutate(scenario = "Prior_weighted")
show_table(prior_contrasts |> mutate(value = 100*value) |> group_by(department) |> summary90())
sensitivity_summary <- prior_model$parameter_diagnostics |>
  filter(grepl("^(b|g)\\[", variable))
show_table(sensitivity_summary)
write_table(sensitivity_summary, "sensitivity_coefficient_summary.csv")
cat("Uniform priors constrain all loadings to [0,1]. Results need not identify those loadings strongly.\n")
cat("The Stronger fixed scenario is outside this prior support; it is not part of this average.\nEND_20\n")

# Section 21 — Separate probability-scale tipping-point illustration (not the fitted logistic model).
cat("BEGIN_21\n")
# Same defaults and valid-probability restriction as the section 21 HTML simulator.
observed_difference <- .06
latent_risk_difference <- .20
threshold <- .02
p2_observed <- .10
q2 <- .10
tipping_curve <- tibble(imbalance = seq(0, .60, by = .001)) |>
  mutate(q1 = q2 + imbalance,
         low1 = p2_observed + observed_difference - latent_risk_difference*q1,
         low2 = p2_observed - latent_risk_difference*q2,
         high1 = low1 + latent_risk_difference,
         high2 = low2 + latent_risk_difference,
         corrected_difference = observed_difference - latent_risk_difference*imbalance,
         valid = q1 >= 0 & q1 <= 1 & low1 >= 0 & high1 <= 1 & low2 >= 0 & high2 <= 1)
tipping_table <- tibble(imbalance = c(0,.1,.2,.3,.4)) |>
  mutate(corrected_difference = observed_difference-latent_risk_difference*imbalance)
tipping_thresholds <- tibble(target = c("Zero contrast", "Decision threshold"),
                            contrast = c(0, threshold)) |>
  mutate(imbalance = (observed_difference-contrast)/latent_risk_difference)
show_table(tipping_table)
show_table(tipping_thresholds)
write_table(tipping_table, "tipping_table.csv")
write_table(tipping_curve, "tipping_curve.csv")
write_table(tipping_thresholds, "tipping_thresholds.csv")
tipping_plot <- tipping_curve |> filter(valid) |>
  ggplot(aes(100*imbalance, 100*corrected_difference)) +
  geom_line(colour = "steelblue", linewidth = 1) +
  geom_hline(yintercept = c(0, 100*threshold), linetype = "dashed") +
  geom_point(data = tipping_thresholds, aes(100*imbalance, 100*contrast), size = 3) +
  theme_minimal() +
  labs(x = "Latent-prevalence imbalance (percentage points)",
       y = "Standardized difference (percentage points)",
       title = "Section 21: when does the conclusion change?",
       subtitle = sprintf("Additive teaching model: observed gap %g points; hidden risk effect %g points",
                          100*observed_difference, 100*latent_risk_difference),
       caption = "Dots mark the decision threshold and zero crossing. Not a posterior fit.")
if (interactive()) print(tipping_plot)
ggsave(file.path(out, "tipping.png"), tipping_plot, width = 9, height = 5, dpi = 150)
stopifnot(any(tipping_curve$valid),
          max(abs(observed_difference-latent_risk_difference*tipping_thresholds$imbalance-
                  tipping_thresholds$contrast)) < 1e-12)
cat("END_21\n")

# Section 23 — Bring together both posterior analyses at the end of Problem 1.
cat("BEGIN_23\n")
all_contrasts <- bind_rows(fixed_contrasts, prior_contrasts)
summary_all <- all_contrasts |> group_by(scenario,department) |> summary90()
write_table(summary_all, "all_scenario_summary.csv")
scenario_plot <- plot_scenarios(summary_all, c(settings$scenario, "Prior_weighted"))
if (interactive()) print(scenario_plot)
ggsave(file.path(out, "scenarios.png"), scenario_plot, width = 11, height = 5, dpi = 150)
show_table(summary_all |> mutate(across(c(mean,lo90,hi90), ~100*.x)))
cat("Section 23 figures: fixed_scenarios.png (19) and scenarios.png (23).\n")
cat("Prior_weighted is a joint fit with independent Uniform(0,1) loadings, not an average of five rows.\nEND_23\n")

# Sections 25–26 — Real grouped counts; expanding rows cannot reconstruct missing ability.
cat("BEGIN_25\n")
ucb <- as_tibble(as.data.frame(UCBAdmissions))
show_table(ucb |> group_by(Admit) |> summarise(total = sum(Freq), .groups = "drop"))
expanded <- ucb |> uncount(Freq)
stopifnot(nrow(expanded) == 4526, all(xtabs(~Admit+Gender+Dept, expanded) == UCBAdmissions))
cat("Expanded rows:", nrow(expanded), "; retained columns:", paste(names(expanded), collapse=", "), "\n")
cat("No latent model is fitted to these real counts by this script.\nEND_25\nBEGIN_26\n")
softmax <- function(x) { w <- exp(x-max(x)); w/sum(w) }
show_table(tibble(department = LETTERS[1:6], baseline = softmax(rep(0,6)), shifted = softmax(c(1,rep(0,5)))))
stopifnot(max(abs(softmax(c(1,rep(0,5))) - softmax(c(1,rep(0,5))+7))) < 1e-12)
cat("END_26\n")

# Section 29 — Conjugate proxy teaching example, not measured UCB applicant proxies.
cat("BEGIN_29\n")
measurements <- c(1.2,.8,1)
show_table(tibble(number_of_proxies = c(1L,3L)) |>
  mutate(variance = 1/(1+number_of_proxies),
         mean = map2_dbl(number_of_proxies, variance, function(n,v) v*sum(measurements[seq_len(n)])),
         sd = sqrt(variance)))
cat("Loadings/noise SD are fixed to 1; conditional independence is assumed. Copies are not new measurements.\nEND_29\n")

# Section 32 — Non-collapsibility can occur without confounding.
cat("BEGIN_32\n")
noncollapsible <- expand_grid(X = 0:1, U = 0:1) |>
  mutate(p = plogis(-2+log(2)*X+2*U)) |> group_by(X) |> summarise(p = mean(p), .groups = "drop")
show_table(noncollapsible)
cat("Conditional OR: 2; marginal OR:", exp(diff(qlogis(noncollapsible$p))), "\nEND_32\n")

# Section 34 — Population PPC: new u, THEN selected department, THEN admission.
cat("BEGIN_34\n")
set.seed(1034)
posterior_rows <- fixed_core$Official_fixed |>
  slice_sample(n = min(500L, nrow(fixed_core$Official_fixed)))
ppc <- map_dfr(seq_len(nrow(posterior_rows)), function(s) {
  pars <- posterior_rows[s, ]
  map_dfr(1:2, function(group) {
    n <- sum(sim$G == group)
    u <- rnorm(n)
    delta <- pars[[paste0("delta[",group,"]")]]
    d <- 1L+rbinom(n, 1, plogis(delta+c(1,0)[group]*u))
    a1 <- pars[[paste0("a[",group,",1]")]]
    a2 <- pars[[paste0("a[",group,",2]")]]
    a <- rbinom(n, 1, plogis(if_else(d == 1L,a1,a2)+u))
    tibble(.replicate = s, G = group, department_2 = mean(d == 2L), admission = mean(a))
  })
})
observed <- sim |> group_by(G) |>
  summarise(department_2 = mean(D == 2), admission = mean(A), .groups = "drop") |>
  pivot_longer(-G, names_to = "statistic", values_to = "observed")
ppc_summary <- ppc |> pivot_longer(c(department_2,admission), names_to = "statistic") |>
  group_by(G, statistic) |> summary90() |> select(-prob_positive) |>
  left_join(observed, by = c("G", "statistic"))
show_table(ppc_summary)
write_table(ppc_summary, "population_ppc.csv")
cat("A good PPC for observables does not establish the latent causal story.\nEND_34\n")

# Section 35 — Prior prediction generates both processes from shared latent applicants.
cat("BEGIN_35\n")
set.seed(1035)
prior_prediction <- map_dfr(1:200, function(s) {
  # One joint parameter draw per synthetic dataset, not one new global prior per person.
  a <- matrix(rnorm(4), 2, 2); delta <- rnorm(2); b <- runif(2); g <- runif(2)
  tibble(G = rep(1:2, each = 250), u = rnorm(500)) |>
    mutate(D = 1L+rbinom(n(),1,plogis(delta[G]+g[G]*u)),
           A = rbinom(n(),1,plogis(a[cbind(G,D)]+b[G]*u))) |>
    group_by(G) |> summarise(department_2 = mean(D == 2L), admission = mean(A), .groups = "drop") |>
    mutate(.replicate = s)
})
show_table(prior_prediction |> pivot_longer(c(department_2,admission), names_to = "statistic") |>
             group_by(G, statistic) |> summary90() |> select(-prob_positive))
cat("These are predictions under the sensitivity-prior model, not its posterior.\nEND_35\n")

# Section 36 — Actual HMC diagnostics and comparison against the original grid route.
cat("BEGIN_36\n")
all_contrasts <- bind_rows(fixed_contrasts, prior_contrasts)
estimand_diagnostics <- all_contrasts |> group_split(scenario) |>
  map_dfr(function(x) contrast_diagnostics(x) |> mutate(scenario = unique(x$scenario)))
write_table(estimand_diagnostics, "contrast_diagnostics.csv")
# Evaluate ONLY the existing grid helper, not the whole script or its plotting side effects.
reference_env <- new.env(parent = globalenv())
reference_env$counts <- as.data.frame(counts)
reference_env$normalize <- function(z) { w <- exp(z-max(z)); w/sum(w) }
expressions <- parse("scripts/A10_walkthroughs.R")
helper <- keep(as.list(expressions), function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
                 identical(x[[2]], as.name("fit_group")))
stopifnot(length(helper) == 1L)
eval(helper[[1]], envir = reference_env)
set.seed(1036)
reference_means <- map_dfr(seq_len(nrow(settings)), function(k) {
  s <- settings[k, ]
  one <- reference_env$fit_group(1,s$b1,s$g1,K=401,Q=601,S=10)
  two <- reference_env$fit_group(2,s$b2,s$g2,K=401,Q=601,S=10)
  tibble(scenario = s$scenario, department = 1:2, grid_mean = one$mean-two$mean)
})
# Allow Monte Carlo error and residual quadrature error; do not demand identical random samples.
parity <- estimand_diagnostics |> inner_join(reference_means, by = c("scenario","department")) |>
  mutate(absolute_difference = abs(mean-grid_mean), tolerance = 5*mcse_mean+.001,
         pass = absolute_difference < tolerance)
show_table(parity |> select(scenario,department,mean,grid_mean,mcse_mean,absolute_difference,tolerance,pass))
write_table(parity, "grid_parity.csv")
# Ordinary cell posterior means supply an independent null-scenario check as well.
cell_reference <- counts |> mutate(p = map2_dbl(admitted, applied, function(y,n) {
  a <- seq(-8,8,length.out=1601)
  w <- reference_env$normalize(dbinom(y,n,plogis(a),log=TRUE)+dnorm(a,log=TRUE))
  sum(w*plogis(a))
})) |> select(G,D,p) |> pivot_wider(names_from=G,names_prefix="G",values_from=p) |>
  transmute(department=D,grid_mean=G1-G2)
cell_parity <- contrast_diagnostics(cell_draws) |> left_join(cell_reference,by="department") |>
  mutate(pass=abs(mean-grid_mean)<5*mcse_mean+.001)
stopifnot(all(parity$pass), all(cell_parity$pass))
chain_report <- bind_rows(all_diagnostics)
write_table(chain_report, "model_diagnostics.csv")
show_table(chain_report)
show_table(estimand_diagnostics |> select(scenario,department,rhat,ess_bulk,ess_tail,mcse_mean))
cat("END_36\n")

# Section 43 — Separate CRM decision calculation, with illustrative economic inputs.
cat("BEGIN_43\n")
show_table(tibble(scenario = c("Optimistic","Intermediate","Pessimistic"), uplift = c(.03,.015,.005)) |>
             mutate(incremental_profit_per_customer = 40*uplift-.50))
cat("Break-even conversion uplift:", .50/40, "\nEND_43\n")

# Save readable summaries, raw contrast draws and an explicit execution record.
summary_all <- all_contrasts |> group_by(scenario,department) |> summary90()
write_table(summary_all, "all_scenario_summary.csv")
write_table(all_contrasts, "contrast_draws.csv")
write_table(settings, "fixed_scenario_settings.csv")
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
writeLines(c(paste("mode:",if(quick) "quick (smoke only)" else "full"),
             paste("CmdStan:",cmdstanr::cmdstan_version()),
             paste("chains:",config$chains),paste("initial warmup per chain:",config$warmup),
             paste("initial retained per chain:",config$sampling),
             paste("max quadrature change:",integration_error),
             "Full SBC and real-UCB latent fits: not run.",
             "Generated Stan and diagnostic CSVs record actual fits and any longer-chain retries.",
             "Matching validated checkpoints preserve core draws and avoid unnecessary refits.",
             "Individual-u chain files use R's temporary directory; compact checkpoints are excluded from Git."),
           file.path(out,"execution.txt"))
if (!quick) {
  # Saving the results first leaves failures inspectable; never label failed diagnostics as passed.
  stopifnot(all(chain_report$pass), all(estimand_diagnostics$rhat < 1.01),
            all(estimand_diagnostics$ess_bulk >= 400), all(estimand_diagnostics$ess_tail >= 400))
  cat("PASS: full-run parameter/estimand diagnostics, quadrature refinement and fixed-scenario grid parity.\n")
} else {
  cat("SMOKE TEST completed. Read diagnostics; no final-inference certification.\n")
}
cat("Results:", out, "\n")
