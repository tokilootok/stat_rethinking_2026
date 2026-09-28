# A09 grant laboratory: tidyverse companion to A09_grants_walkthrough.R.
# Run from the repository root:
# Rscript scripts/A09_grants_walkthrough_tidyverse.R
# Packages: dplyr, tidyr, purrr, tibble, readr, ggplot2; rethinking for data only.
# Independent Normal-logit cell priors; numerical grid integration, not MCMC.
# Each section addresses the corresponding lesson in A09. Run in order.
# Outputs are separate from the original walkthrough's figures and tables.
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(readr)
  library(ggplot2)
})
output_dir <- "courses/figures/A09_grants_tidyverse"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Section 53 — What is one observation? One discipline × gender cell.
cat("BEGIN_53\n")
data("NWOGrants", package = "rethinking")
d <- as_tibble(NWOGrants) |>
  transmute(discipline_id = as.integer(discipline),
            discipline = as.character(discipline),
            gender_id = if_else(gender == "f", 1L, 2L),
            gender = if_else(gender == "f", "female", "male"),
            applications, awards) |>
  arrange(discipline_id, gender_id)
# Validate the keys before any joins: duplicates must never multiply rows.
validate_cells <- function(cells) {
  stopifnot(!anyNA(cells),
            nrow(distinct(cells, discipline_id, gender)) == nrow(cells),
            all(cells$applications > 0),
            all(cells$awards >= 0 & cells$awards <= cells$applications),
            all(cells$awards == floor(cells$awards)),
            all(cells$applications == floor(cells$applications)))
  pairs <- cells |> count(discipline_id)
  stopifnot(all(pairs$n == 2), setequal(cells$gender, c("female", "male")))
  invisible(TRUE)
}
validate_cells(d)
print(d |> summarise(cells = n(), applications = sum(applications), awards = sum(awards)))
cat("END_53\n")

# Section 55 — Awards need denominators; start with the observed gap.
cat("BEGIN_55\n")
totals <- d |> group_by(gender) |>
  summarise(across(c(awards, applications), sum), .groups = "drop") |>
  mutate(rate = awards / applications)
print(totals)
raw_gap <- totals |> select(gender, rate) |> pivot_wider(names_from = gender, values_from = rate) |>
  transmute(gap_pp = 100 * (female - male))
print(raw_gap)
print(tibble(awards = c(10, 10), applications = c(20, 100)) |> mutate(rate = awards / applications))
cat("END_55\n")

# Section 57 — Why grouped binomial counts retain the likelihood information.
cat("BEGIN_57\n")
first_cell <- slice(d, 1)
y <- c(rep(1, first_cell$awards), rep(0, first_cell$applications - first_cell$awards))
comparison <- tibble(p = c(.10, .25, .50)) |>
  mutate(log_individual = map_dbl(p, \(p) sum(dbinom(y, 1, p, log = TRUE))),
         log_grouped = dbinom(first_cell$awards, first_cell$applications, p, log = TRUE),
         gap = log_grouped - log_individual)
print(comparison)
stopifnot(max(abs(comparison$gap - lchoose(first_cell$applications, first_cell$awards))) < 1e-10)
cat("END_57\n")

# Section 59 — Odds ratios change under averaging even without confounding.
cat("BEGIN_59\n")
marginal <- expand_grid(x = 0:1, u = 0:1) |>
  mutate(p = plogis(-2 + log(2) * x + 2 * u)) |>
  group_by(x) |> summarise(p = mean(p), .groups = "drop")
print(tibble(conditional_OR = 2, marginal_OR = exp(diff(qlogis(marginal$p))),
             marginal_probability_difference = diff(marginal$p)))
cat("END_59\n")

# Sections 60–61 — Define the population before averaging probabilities.
cat("BEGIN_60\n")
weights <- d |> select(discipline_id, discipline, gender, applications) |>
  pivot_wider(names_from = gender, values_from = applications) |>
  mutate(w_f = female / sum(female), w_m = male / sum(male),
         w_pool = (female + male) / sum(female + male), w_equal = 1 / n()) |>
  select(discipline_id, discipline, starts_with("w_"))
# Join by ID, never by row position. Reject missing, extra or duplicate weights.
attach_weights <- function(probabilities, weight_table = weights) {
  stopifnot(!anyDuplicated(weight_table$discipline_id),
            setequal(probabilities$discipline_id, weight_table$discipline_id))
  left_join(probabilities, select(weight_table, -any_of("discipline")), by = "discipline_id")
}
# Input: .draw, discipline_id, gender, p. Output: one row per draw and target.
contrasts <- function(probabilities, weight_table = weights) {
  probabilities |> select(.draw, discipline_id, gender, p) |>
    pivot_wider(names_from = gender, values_from = p) |>
    attach_weights(weight_table) |>
    group_by(.draw) |>
    summarise(Total = sum(w_f * female) - sum(w_m * male),
              Standardized = sum(w_pool * (female - male)), .groups = "drop") |>
    pivot_longer(c(Total, Standardized), names_to = "estimand", values_to = "value")
}
raw_probabilities <- d |> transmute(.draw = 1L, discipline_id, gender, p = awards / applications)
print(contrasts(raw_probabilities) |> mutate(percentage_points = 100 * value))
cat("END_60\nBEGIN_61\n")
print(weights)
weight_sums <- weights |> summarise(across(starts_with("w_"), sum))
stopifnot(all(abs(unlist(weight_sums) - 1) < 1e-12))
cat("END_61\n")

# Section 62 — Audit the actual data before fitting.
cat("BEGIN_62\n")
validate_cells(d)
stopifnot(nrow(d) == 18, n_distinct(d$discipline_id) == 9,
          sum(d$applications) == 2823, sum(d$awards) == 467)
print(d, n = Inf)
cat("END_62\n")

# Section 64 — One independent posterior per cell: alpha ~ Normal(mu, sd),
# awards ~ Binomial(applications, plogis(alpha)). The grid is uniform in alpha,
# so normalized likelihood × prior weights approximate integration in alpha.
# No Jacobian is needed: we transform sampled alpha to p, not the integration measure.
cat("BEGIN_64\n")
fit_cells <- function(cells, prior_mean = 0, prior_sd = 1, step = .01, S = 60000) {
  validate_cells(cells)
  stopifnot(prior_sd > 0, step > 0, S > 0, S == floor(S))
  grid <- tibble(alpha = seq(-10, 10, by = step)) |> mutate(p = plogis(alpha))
  fitted <- cells |>
    mutate(grid = map2(awards, applications, \(a, n) {
      grid |> mutate(log_weight = dbinom(a, n, p, log = TRUE) +
                       dnorm(alpha, prior_mean, prior_sd, log = TRUE),
                     mass = exp(log_weight - max(log_weight)), mass = mass / sum(mass))
    }), boundary_mass = map_dbl(grid, \(g) sum(g$mass[abs(g$alpha) > 9])),
    p_mean = map_dbl(grid, \(g) sum(g$p * g$mass)))
  stopifnot(max(fitted$boundary_mass) < 1e-8)
  draws <- fitted |>
    mutate(samples = map(grid, \(g) tibble(.draw = seq_len(S),
                                         p = sample(g$p, S, replace = TRUE, prob = g$mass)))) |>
    select(discipline_id, gender, samples) |> unnest(samples)
  list(draws = draws, means = fitted |> select(discipline_id, gender, p_mean))
}
# Works on a value column and preserves any existing grouping (e.g. estimand).
summary90 <- function(tbl) {
  tbl |> summarise(mean = mean(value), lo90 = quantile(value, .05),
                   hi90 = quantile(value, .95), Pr_positive = mean(value > 0), .groups = "drop")
}
set.seed(9064)
fit_grid <- fit_cells(d)
p_draw <- fit_grid$draws
fine <- fit_cells(d, step = .005, S = 10)
refinement <- inner_join(fit_grid$means, fine$means, by = c("discipline_id", "gender"),
                         suffix = c("", "_fine")) |>
  summarise(max_change = max(abs(p_mean - p_mean_fine))) |> pull(max_change)
stopifnot(refinement < 1e-6)
print(tibble(cells = nrow(fit_grid$means), draws_per_cell = n_distinct(p_draw$.draw), refinement))
cat("END_64\n")

# Section 66 — Compare total and standardized female-minus-male differences.
cat("BEGIN_66\n")
posterior_contrasts <- contrasts(p_draw)
posterior_table <- posterior_contrasts |> group_by(estimand) |> summary90()
print(posterior_table)
probabilities <- p_draw |> attach_weights() |> group_by(.draw, gender) |>
  summarise(Observed_composition = sum(p * if_else(gender == "female", w_f, w_m)),
            Common_composition = sum(p * w_pool), .groups = "drop") |>
  group_by(gender) |> summarise(across(ends_with("composition"), mean), .groups = "drop")
print(probabilities)
write_csv(posterior_table, file.path(output_dir, "posterior_summary.csv"))
contrast_plot <- posterior_table |> mutate(across(c(mean, lo90, hi90), \(x) 100 * x)) |>
  ggplot(aes(x = estimand, y = mean, ymin = lo90, ymax = hi90)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_pointrange(colour = "steelblue") + coord_flip() + theme_minimal() +
  labs(title = "NWOGrants: changing the target composition", x = NULL,
       y = "Female minus male (percentage points)", subtitle = "Posterior means and central 90% intervals")
ggsave(file.path(output_dir, "contrasts.png"), contrast_plot, width = 9, height = 4.5, dpi = 150)
cat("END_66\n")

# Section 68 — Separate within-discipline differences from composition shifts.
cat("BEGIN_68\n")
cell_means <- fit_grid$means |> rename(p = p_mean) |> mutate(.draw = 1L)
decomposition <- cell_means |> select(-.draw) |>
  pivot_wider(names_from = gender, values_from = p) |> attach_weights() |>
  summarise(common = sum(w_pool * (female - male)),
            female_shift = sum((w_f - w_pool) * female),
            male_shift = -sum((w_m - w_pool) * male)) |>
  mutate(reconstructed = common + female_shift + male_shift)
correct <- contrasts(cell_means) |> filter(estimand == "Total") |> pull(value)
stopifnot(abs(decomposition$reconstructed - correct) < 1e-12)
print(decomposition |> mutate(across(everything(), \(x) 100 * x)))
cat("END_68\n")

# Section 69 — A row permutation must not silently change the answer.
cat("BEGIN_69\n")
# Deliberate BUG for teaching: keep IDs but give them other disciplines' weights.
wrong_weights <- weights |> mutate(across(c(w_f, w_m), rev))
wrong <- contrasts(cell_means, wrong_weights) |> filter(estimand == "Total") |> pull(value)
# Correct: reorder whole rows, keeping each discipline attached to its weights.
recovered <- contrasts(cell_means, arrange(weights, desc(discipline_id))) |>
  filter(estimand == "Total") |> pull(value)
stopifnot(abs(correct - recovered) < 1e-12)
print(tibble(correct_pp = 100 * correct, wrong_pp = 100 * wrong, recovered_pp = 100 * recovered))
cat("END_69\n")

# Section 70 — Could the prior generate plausible totals before observing awards?
cat("BEGIN_70\n")
set.seed(9070)
prior_counts <- function(mu, sd) {
  expand_grid(cell = seq_len(nrow(d)), .draw = seq_len(10000)) |>
    mutate(p = plogis(rnorm(n(), mu, sd)),
           awards = rbinom(n(), d$applications[cell], p)) |>
    group_by(.draw) |> summarise(value = sum(awards), .groups = "drop")
}
prior_default <- prior_counts(0, 1)
prior_centered <- prior_counts(qlogis(.17), .75)
prior_summary <- bind_rows(Normal_0_1 = prior_default, Rate_centered = prior_centered, .id = "prior") |>
  group_by(prior) |> summary90() |> select(-Pr_positive)
print(prior_summary)
cat("Observed awards:", sum(d$awards), "\n")
# The 0.17 center is data-informed sensitivity, not an independent prior validation.
cat("END_70\n")

# Section 71 — Expected counts versus realized counts in the smallest cell.
cat("BEGIN_71\n")
set.seed(9071)
smallest <- slice_min(d, applications, n = 1, with_ties = FALSE)
small_draws <- p_draw |> inner_join(smallest, by = c("discipline_id", "gender")) |>
  transmute(.draw, Expected_count = applications * p,
            Replicated_count = rbinom(n(), applications, p))
print(smallest)
print(small_draws |> pivot_longer(-.draw, names_to = "quantity", values_to = "value") |>
        group_by(quantity) |> summary90() |> select(-Pr_positive))
cat("END_71\n")

# Section 72 — Posterior predictive checks (PPC): compare replicated and observed data.
cat("BEGIN_72\n")
set.seed(9072)
replicated <- p_draw |> left_join(d, by = c("discipline_id", "gender")) |>
  mutate(y_rep = rbinom(n(), applications, p), expected = applications * p,
         variance = applications * p * (1 - p),
         T_obs = (awards - expected)^2 / variance, T_rep = (y_rep - expected)^2 / variance)
rep_totals <- bind_rows(
  replicated |> group_by(.draw, gender) |> summarise(value = sum(y_rep), .groups = "drop"),
  replicated |> group_by(.draw) |> summarise(value = sum(y_rep), .groups = "drop") |> mutate(gender = "All"))
observed <- bind_rows(totals |> transmute(gender, observed = awards),
                      tibble(gender = "All", observed = sum(d$awards)))
ppc_summary <- rep_totals |> group_by(gender) |> summary90() |> select(-Pr_positive) |>
  left_join(observed, by = "gender")
print(ppc_summary)
pearson <- replicated |> group_by(.draw) |>
  summarise(across(c(T_obs, T_rep), sum), .groups = "drop")
cat("Posterior predictive Pearson tail fraction:", mean(pearson$T_rep >= pearson$T_obs), "\n")
ppc_plot <- rep_totals |> filter(gender == "All") |>
  ggplot(aes(value)) + geom_histogram(bins = 35, fill = "lightblue", colour = "white") +
  geom_vline(xintercept = sum(d$awards), colour = "firebrick") + theme_minimal() +
  labs(title = "NWOGrants: replicated total awards", subtitle = "Red line: observed total",
       x = "Awards in a replicated cohort", y = "Posterior predictive draws")
ggsave(file.path(output_dir, "ppc.png"), ppc_plot, width = 8, height = 4.8, dpi = 150)
write_csv(ppc_summary, file.path(output_dir, "ppc_summary.csv"))
cat("END_72\n")

# Section 73 — Can we recover chosen truths with unequal composition?
# This single recovery exercise is NOT simulation-based calibration (SBC).
cat("BEGIN_73\n")
set.seed(9073)
truth <- d |> distinct(discipline_id) |>
  mutate(p = seq(.08, .30, length.out = n()))
simulated <- d |> left_join(truth, by = "discipline_id") |>
  mutate(awards = rbinom(n(), applications, p))
recovery <- fit_cells(simulated |> select(-p), S = 30000)
true_contrast <- simulated |> mutate(.draw = 1L) |> contrasts() |> rename(truth = value)
recovery_table <- contrasts(recovery$draws) |> group_by(estimand) |> summary90() |>
  left_join(select(true_contrast, estimand, truth), by = "estimand")
print(recovery_table)
stopifnot(abs(filter(true_contrast, estimand == "Standardized")$truth) < 1e-12)
cat("END_73\n")

# Section 75 — Distinguish changing the prior from changing the target population.
cat("BEGIN_75\n")
set.seed(9075)
alternative <- fit_cells(d, prior_mean = qlogis(.17), prior_sd = .75)
equal_weights <- weights |> mutate(w_pool = w_equal)
scenarios <- bind_rows(
  posterior_contrasts |> mutate(scenario = if_else(estimand == "Total", "Default_total", "Default_pooled")),
  contrasts(p_draw, equal_weights) |> filter(estimand == "Standardized") |> mutate(scenario = "Default_equal"),
  contrasts(alternative$draws) |> filter(estimand == "Standardized") |> mutate(scenario = "Centered_prior_pooled"))
sensitivity_table <- scenarios |> group_by(scenario) |> summary90()
print(sensitivity_table)
# Independent cell posteriors do not change when we omit a target discipline.
# Hence posterior means suffice for these linear, reweighted posterior means.
leave_out <- map_dfr(weights$discipline_id, \(omitted) {
  reduced_weights <- weights |> filter(discipline_id != omitted) |>
    mutate(across(c(w_f, w_m, w_pool), \(w) w / sum(w)))
  cell_means |> filter(discipline_id != omitted) |> contrasts(reduced_weights) |>
    filter(estimand == "Standardized") |> transmute(omitted, mean = value)
})
print(leave_out)
cat("This changes the target population; it does not assess prediction for a new discipline.\n")
write_csv(sensitivity_table, file.path(output_dir, "sensitivity_summary.csv"))
cat("END_75\n")

# Section 76 — A random shared cell probability creates overdispersion (toy only).
cat("BEGIN_76\n")
n <- 100; mu <- .17; kappa <- 20
# P ~ Beta(mu*kappa, (1-mu)*kappa); Y | P ~ Binomial(n, P).
print(tibble(model = c("Binomial", "Beta-binomial"), mean = n * mu,
             variance = n * mu * (1 - mu) * c(1, (n + kappa) / (1 + kappa))))
cat("END_76\n")

# Section 78 — What has this execution actually checked?
cat("BEGIN_78\n")
checks <- tibble(check = c("Complete valid cells", "Normalized weights", "Grid refinement",
                          "Decomposition", "Permutation defense", "Draw alignment"),
                passed = c(validate_cells(d), all(abs(unlist(weight_sums) - 1) < 1e-12),
                           refinement < 1e-6, abs(decomposition$reconstructed - correct) < 1e-12,
                           abs(correct - recovered) < 1e-12,
                           nrow(distinct(p_draw, .draw, discipline_id, gender)) == nrow(p_draw) &&
                             nrow(p_draw) == 60000 * nrow(d)))
print(checks)
stopifnot(all(checks$passed))
cat("PASS. Grid checks only: no HMC diagnostics, full SBC or cross-language validation.\n")
cat("Results saved in:", output_dir, "\nEND_78\n")
