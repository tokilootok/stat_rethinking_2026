# A01: Random Allocation Game.
# All participants flip a fair coin. Honest participants claim exactly on
# heads; dishonest participants always claim.

rag_ways_fixed <- function(honest, participants = 10L, claims = 8L) {
  input <- validate_a01_input(participants, claims)
  honest <- assert_scalar_count(honest, "honest", 0L, input$participants)
  if (honest < input$nonclaims) return(0)
  choose(honest, input$nonclaims) * 2^(input$participants - honest)
}

rag_ways_variable <- function(honest, participants = 10L, claims = 8L) {
  input <- validate_a01_input(participants, claims)
  honest <- assert_scalar_count(honest, "honest", 0L, input$participants)
  if (honest < input$nonclaims) return(0)
  choose(input$participants, input$nonclaims) *
    choose(input$claims, honest - input$nonclaims) *
    2^(input$participants - honest)
}

rag_claim_probability <- function(claims, honest, participants = 10L) {
  participants <- assert_scalar_count(participants, "participants", lower = 1L)
  claims <- assert_scalar_count(claims, "claims", 0L, participants)
  honest <- assert_scalar_count(honest, "honest", 0L, participants)
  honest_claims <- claims - (participants - honest)
  if (honest_claims < 0L || honest_claims > honest) return(0)
  dbinom(honest_claims, size = honest, prob = 0.5)
}

enumerate_rag_fixed <- function(honest, participants = 10L, claims = 8L) {
  input <- validate_a01_input(participants, claims)
  honest <- assert_scalar_count(honest, "honest", 0L, input$participants)
  if (input$participants > 20L) {
    stop("Brute-force enumeration is limited to 20 participants.", call. = FALSE)
  }
  paths <- 0:(2^input$participants - 1)
  claim_counts <- vapply(paths, function(path) {
    flips <- bitwAnd(bitwShiftR(path, 0:(input$participants - 1L)), 1L)
    honest_claims <- if (honest == 0L) 0L else sum(flips[seq_len(honest)])
    honest_claims + input$participants - honest
  }, integer(1))
  sum(claim_counts == input$claims)
}

simulate_rag <- function(honest, participants = 10L, repetitions = 1000L,
                         seed = NULL) {
  participants <- assert_scalar_count(participants, "participants", lower = 1L)
  honest <- assert_scalar_count(honest, "honest", 0L, participants)
  repetitions <- assert_scalar_count(repetitions, "repetitions", lower = 1L)
  if (!is.null(seed)) set.seed(seed)
  honest_claims <- if (honest == 0L) rep.int(0L, repetitions) else
    rbinom(repetitions, size = honest, prob = 0.5)
  honest_claims + participants - honest
}

rag_hypothesis_table <- function(participants = 10L, claims = 8L) {
  input <- validate_a01_input(participants, claims)
  hypotheses <- 0:input$participants
  fixed_ways <- vapply(hypotheses, rag_ways_fixed, numeric(1),
    participants = input$participants, claims = input$claims)
  variable_ways <- vapply(hypotheses, rag_ways_variable, numeric(1),
    participants = input$participants, claims = input$claims)
  likelihood <- vapply(hypotheses, rag_claim_probability, numeric(1),
    claims = input$claims, participants = input$participants)
  data.frame(
    honest = hypotheses,
    fixed_ways = fixed_ways,
    variable_ways = variable_ways,
    likelihood = likelihood,
    fixed_posterior = normalize_probability(likelihood),
    variable_microstate_posterior = normalize_probability(variable_ways),
    check.names = FALSE
  )
}

summarize_a01 <- function(table, participants = 10L, claims = 8L) {
  input <- validate_a01_input(participants, claims)
  required <- c("honest", "fixed_ways", "variable_ways", "likelihood")
  if (!all(required %in% names(table))) {
    stop("The A01 hypothesis table is missing required columns.", call. = FALSE)
  }
  at_maximum <- function(x) {
    abs(x - max(x)) <= sqrt(.Machine$double.eps) * max(1, abs(max(x)))
  }
  list(
    all_honest_ways = table$fixed_ways[table$honest == input$participants],
    five_honest_fixed_ways = table$fixed_ways[table$honest == 5L],
    five_honest_variable_ways = table$variable_ways[table$honest == 5L],
    fixed_count_modes = table$honest[at_maximum(table$fixed_ways)],
    variable_count_modes = table$honest[at_maximum(table$variable_ways)],
    likelihood_modes = table$honest[at_maximum(table$likelihood)],
    interpretation_note = paste(
      "Raw variable-identity path counts contain choose(N, H) possible identity",
      "assignments. They imply a different prior over H unless divided by that",
      "number of assignments. The normalized likelihood is reported separately."
    )
  )
}
