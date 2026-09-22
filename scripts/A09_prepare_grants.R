#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(rethinking))

args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args) >= 1L) args[[1L]] else "artifacts/A09/nwo_grants.csv"

data(NWOGrants)
d <- NWOGrants

out <- data.frame(
  cell_id = seq_len(nrow(d)),
  discipline_id = as.integer(d$discipline),
  discipline = as.character(d$discipline),
  gender_id = ifelse(d$gender == "f", 1L, 2L),
  gender = ifelse(d$gender == "f", "female", "male"),
  gender_c = ifelse(d$gender == "f", -0.5, 0.5),
  applications = as.integer(d$applications),
  awards = as.integer(d$awards),
  stringsAsFactors = FALSE
)

stopifnot(
  nrow(out) == 18L,
  length(unique(out$discipline_id)) == 9L,
  all(table(out$discipline_id) == 2L),
  all(out$awards >= 0L),
  all(out$awards <= out$applications),
  sum(out$applications) == 2823L,
  sum(out$awards) == 467L
)

dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
write.csv(out, output, row.names = FALSE)
cat(sprintf("Wrote %d cells and %d applications to %s\n",
            nrow(out), sum(out$applications), output))

