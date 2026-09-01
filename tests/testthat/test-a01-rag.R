test_that("A01 contracts reject impossible observations", {
  expect_error(validate_a01_input(0, 0), "participants")
  expect_error(validate_a01_input(10, 11), "claims")
  expect_error(rag_ways_fixed(11, 10, 8), "honest")
})

test_that("A01 analytical counts answer the homework", {
  expect_equal(rag_ways_fixed(10, 10, 8), 45)
  expect_equal(rag_ways_fixed(5, 10, 8), 320)
  expect_equal(rag_ways_variable(5, 10, 8), 80640)
  result <- summarize_a01(rag_hypothesis_table(10, 8), 10, 8)
  expect_equal(result$fixed_count_modes, c(3L, 4L))
  expect_equal(result$variable_count_modes, c(4L, 5L))
  expect_equal(result$likelihood_modes, c(3L, 4L))
})

test_that("A01 formula agrees with brute-force enumeration", {
  for (honest in 0:10) {
    expect_equal(rag_ways_fixed(honest, 10, 8),
      enumerate_rag_fixed(honest, 10, 8), info = paste("honest =", honest))
  }
})

test_that("A01 likelihood is normalized over possible claim counts", {
  for (honest in 0:10) {
    probabilities <- vapply(0:10, rag_claim_probability, numeric(1),
      honest = honest, participants = 10)
    expect_equal(sum(probabilities), 1)
  }
})

test_that("A01 simulation is reproducible", {
  first <- simulate_rag(5, repetitions = 20, seed = 42)
  second <- simulate_rag(5, repetitions = 20, seed = 42)
  expect_identical(first, second)
})
