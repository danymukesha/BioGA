library(BioGA)
library(testthat)

test_that("crossover_cpp returns matrix with correct dimensions and valid values", {
  genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
  population <- BioGA::initialize_population_cpp(genomic_data, population_size = 5)
  fitness <- BioGA::evaluate_fitness_cpp(genomic_data, population, c(1.0, 0.5))
  selected_parents <- BioGA::selection_cpp(population, fitness, num_parents = 2)
  offspring_size <- 2
  crossover_rate <- 0.9
  eta_c <- 20.0
  
  offspring <- BioGA::crossover_cpp(selected_parents, offspring_size, crossover_rate, eta_c)
  
  expect_equal(nrow(offspring), offspring_size)
  expect_equal(ncol(offspring), nrow(genomic_data))
  
  for (i in 1:nrow(offspring)) {
    for (j in 1:ncol(offspring)) {
      expect_true(offspring[i, j] >= min(selected_parents[, j]) &&
                    offspring[i, j] <= max(selected_parents[, j]))
    }
  }
  
  expect_no_warning(BioGA::crossover_cpp(selected_parents, offspring_size, crossover_rate, eta_c),
                    message = "crossover_cpp should not produce warnings")
  expect_no_error(BioGA::crossover_cpp(selected_parents, offspring_size, crossover_rate, eta_c),
                  message = "crossover_cpp should not produce errors")
})

# test_that("crossover_cpp handles edge cases", {
#   selected_parents <- matrix(rnorm(1), nrow = 1, ncol = 1)
#   offspring <- BioGA::crossover_cpp(selected_parents, 1, 0.0) # No crossover
#   expect_equal(nrow(offspring), 1)
#   expect_equal(ncol(offspring), 1)
#   expect_equal(offspring[1, 1], selected_parents[1, 1]) # No crossover should copy parent
# })
