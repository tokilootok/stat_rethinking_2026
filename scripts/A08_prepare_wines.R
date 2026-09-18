#!/usr/bin/env Rscript

# Build the one canonical analysis table used by the R, Stan, and PyMC
# implementations in Lecture A08.

suppressPackageStartupMessages(library(rethinking))

args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args) >= 1L) args[[1L]] else "artifacts/A08/wines2012.csv"

data(Wines2012)
d <- Wines2012

judge_levels <- levels(d$judge)
wine_levels <- levels(d$wine)
judge_id <- as.integer(d$judge)
wine_id <- as.integer(d$wine)

# These attributes are constant within judge and wine. Stop rather than silently
# accepting malformed data if that contract is ever violated.
judge_country <- vapply(
  seq_along(judge_levels),
  function(j) unique(d$judge.amer[judge_id == j]),
  numeric(1)
)
wine_origin <- vapply(
  seq_along(wine_levels),
  function(w) unique(d$wine.amer[wine_id == w]),
  numeric(1)
)

stopifnot(
  length(d$score) == 180L,
  length(judge_levels) == 9L,
  length(wine_levels) == 20L,
  all(lengths(lapply(seq_along(judge_levels), function(j) unique(d$judge.amer[judge_id == j]))) == 1L),
  all(lengths(lapply(seq_along(wine_levels), function(w) unique(d$wine.amer[wine_id == w]))) == 1L)
)

score_mean <- mean(d$score)
score_sd <- sd(d$score)

out <- data.frame(
  observation = seq_len(nrow(d)),
  score = d$score,
  score_z = as.numeric((d$score - score_mean) / score_sd),
  judge_id = judge_id,
  judge = as.character(d$judge),
  judge_american = as.integer(d$judge.amer),
  # Centered binary predictors make the coefficient the US-minus-France
  # contrast directly: (+0.5) - (-0.5) = 1.
  judge_country_c = d$judge.amer - 0.5,
  wine_id = wine_id,
  wine = as.character(d$wine),
  wine_american = as.integer(d$wine.amer),
  wine_origin_c = d$wine.amer - 0.5,
  flight = as.character(d$flight),
  stringsAsFactors = FALSE
)

dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
write.csv(out, output, row.names = FALSE)

cat(sprintf("Wrote %d rows, %d judges, and %d wines to %s\n",
            nrow(out), length(judge_levels), length(wine_levels), output))

