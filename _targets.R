library(targets)
library(tarchetypes)
library(quarto)

tar_option_set(packages = character(0), seed = 20260819, error = "stop")

source("R/contracts.R")
source("R/a01_rag.R")
source("R/a02_globe.R")

list(
  tar_target(a01_input,
    validate_a01_input(participants = 10L, claims = 8L)),
  tar_target(a01_hypotheses,
    rag_hypothesis_table(a01_input$participants, a01_input$claims)),
  tar_target(a01_summary,
    summarize_a01(a01_hypotheses, a01_input$participants, a01_input$claims)),

  tar_target(a02_input,
    validate_a02_input(water = 3L, land = 11L, future_tosses = 5L)),
  tar_target(a02_grid, make_probability_grid(size = 2001L)),
  tar_target(a02_grid_posterior,
    globe_grid_posterior(a02_input$water, a02_input$land, a02_grid)),
  tar_target(a02_exact_posterior,
    globe_exact_posterior(a02_input$water, a02_input$land)),
  tar_target(a02_exact_predictive,
    globe_exact_predictive(a02_input$water, a02_input$land,
      a02_input$future_tosses)),
  tar_target(a02_grid_predictive,
    globe_grid_predictive(a02_grid_posterior, a02_input$future_tosses)),
  tar_target(a02_predictive_simulation,
    simulate_globe_predictive(a02_input$water, a02_input$land,
      a02_input$future_tosses, draws = 50000L, seed = 20260202L)),
  tar_target(a02_summary,
    summarize_a02(a02_grid_posterior, a02_exact_posterior,
      a02_exact_predictive, a02_grid_predictive)),

  tar_quarto(a01_a02_report, path = "reports/A01_A02.qmd", quiet = TRUE)
)
