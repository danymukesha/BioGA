library(BioGA)
library(testthat)

test_that("evaluate_fitness_cpp returns matrix with correct dimensions and values", {
  
  genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
  population <- BioGA::initialize_population_cpp(genomic_data, population_size = 5)
  weights <- c(1.0, 0.5) 
  
  fitness <- BioGA::evaluate_fitness_cpp(genomic_data, population, weights)
  
  expect_equal(nrow(fitness), 5) 
  expect_equal(ncol(fitness), 2) 
  
  expect_true(all(fitness >= 0))
  
  expect_no_warning(BioGA::evaluate_fitness_cpp(genomic_data, population, weights),
                    message = "evaluate_fitness_cpp should not produce warnings")
  expect_no_error(BioGA::evaluate_fitness_cpp(genomic_data, population, weights),
                  message = "evaluate_fitness_cpp should not produce errors")
})

# test_that("evaluate_fitness_cpp handles edge cases", {
#   genomic_data <- matrix(0, nrow = 1, ncol = 1) 
#   population <- matrix(0, nrow = 1, ncol = 1) 
#   weights <- c(1.0)
#   
#   fitness <- BioGA::evaluate_fitness_cpp(genomic_data, population, weights)
#   expect_equal(nrow(fitness), 1)
#   expect_equal(ncol(fitness), 1)
#   expect_true(is.numeric(fitness[1, 1]))
# })
