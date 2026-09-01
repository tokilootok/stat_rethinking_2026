# Shared validation helpers for the course pipelines.

assert_scalar_count <- function(x, name, lower = 0L, upper = Inf) {
  valid <- length(x) == 1L && is.numeric(x) && is.finite(x) &&
    x == floor(x) && x >= lower && x <= upper
  if (!valid) {
    stop(sprintf("`%s` must be one integer between %s and %s.",
      name, format(lower), format(upper)), call. = FALSE)
  }
  as.integer(x)
}

assert_positive_scalar <- function(x, name) {
  valid <- length(x) == 1L && is.numeric(x) && is.finite(x) && x > 0
  if (!valid) {
    stop(sprintf("`%s` must be one finite positive number.", name), call. = FALSE)
  }
  as.numeric(x)
}

normalize_probability <- function(weights, name = "weights") {
  if (!is.numeric(weights) || any(!is.finite(weights)) || any(weights < 0)) {
    stop(sprintf("`%s` must contain finite non-negative values.", name), call. = FALSE)
  }
  total <- sum(weights)
  if (total <= 0) {
    stop(sprintf("`%s` must have a positive sum.", name), call. = FALSE)
  }
  weights / total
}

validate_a01_input <- function(participants, claims) {
  participants <- assert_scalar_count(participants, "participants", lower = 1L)
  claims <- assert_scalar_count(claims, "claims", lower = 0L, upper = participants)
  list(participants = participants, claims = claims,
    nonclaims = participants - claims)
}

validate_a02_input <- function(water, land, future_tosses = 5L) {
  water <- assert_scalar_count(water, "water", lower = 0L)
  land <- assert_scalar_count(land, "land", lower = 0L)
  future_tosses <- assert_scalar_count(future_tosses, "future_tosses", lower = 1L)
  if (water + land == 0L) {
    stop("At least one globe-tossing observation is required.", call. = FALSE)
  }
  list(water = water, land = land, observations = water + land,
    future_tosses = future_tosses)
}
