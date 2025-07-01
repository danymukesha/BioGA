library(BioGA)
library(testthat)

test_that("mutation_cpp returns matrix with correct dimensions and valid mutations", {
  genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
  population <- BioGA::initialize_population_cpp(genomic_data, population_size = 5)
  fitness <- BioGA::evaluate_fitness_cpp(genomic_data, population, c(1.0, 0.5))
  selected_parents <- BioGA::selection_cpp(population, fitness, num_parents = 2)
  offspring <- BioGA::crossover_cpp(selected_parents, offspring_size = 2)
  mutation_rate <- 0.1
  iteration <- 1
  max_iterations <- 100
  
  mutated_offspring <- BioGA::mutation_cpp(offspring, mutation_rate, iteration, max_iterations)
  
  expect_equal(nrow(mutated_offspring), nrow(offspring))
  expect_equal(ncol(mutated_offspring), ncol(offspring))
  
  expect_true(any(mutated_offspring != offspring))
  
  expect_no_warning(BioGA::mutation_cpp(offspring, mutation_rate, iteration, max_iterations),
                    message = "mutation_cpp should not produce warnings")
  expect_no_error(BioGA::mutation_cpp(offspring, mutation_rate, iteration, max_iterations),
                  message = "mutation_cpp should not produce errors")
})

test_that("mutation_cpp handles zero mutation rate and network constraints", {
  offspring <- matrix(rnorm(20), nrow = 2, ncol = 10)
  network <- matrix(0, nrow = 10, ncol = 10) # Zero network (no constraints)
  
  mutated_offspring <- BioGA::mutation_cpp(offspring, 0.0, 1, 100, network)
  expect_equal(mutated_offspring, offspring)
  
  # test with network constraints
  network <- matrix(1.0, nrow = 10, ncol = 10) # Full constraints
  mutated_offspring <- BioGA::mutation_cpp(offspring, 0.1, 1, 100, network)
  expect_equal(nrow(mutated_offspring), nrow(offspring))
  expect_equal(ncol(mutated_offspring), ncol(offspring))
})