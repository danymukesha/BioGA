library(BioGA)
library(testthat)

test_that("selection_cpp returns matrix with correct dimensions", {
  genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
  population <- BioGA::initialize_population_cpp(genomic_data, population_size = 5)
  fitness <- BioGA::evaluate_fitness_cpp(genomic_data, population, c(1.0, 0.5))
  num_parents <- 2
  
  selected_parents <- BioGA::selection_cpp(population, fitness, num_parents)
  
  expect_equal(nrow(selected_parents), num_parents)
  expect_equal(ncol(selected_parents), nrow(genomic_data))
  
  for (i in 1:nrow(selected_parents)) {
    expect_true(any(apply(population, 1, function(row) all(row == selected_parents[i, ]))))
  }
  
  expect_no_warning(BioGA::selection_cpp(population, fitness, num_parents),
                    message = "selection_cpp should not produce warnings")
  expect_no_error(BioGA::selection_cpp(population, fitness, num_parents),
                  message = "selection_cpp should not produce errors")
})

test_that("selection_cpp handles edge cases", {
  population <- matrix(rnorm(10), nrow = 1, ncol = 10)
  fitness <- matrix(1.0, nrow = 1, ncol = 1)
  selected_parents <- BioGA::selection_cpp(population, fitness, 1)
  expect_equal(nrow(selected_parents), 1)
  expect_equal(ncol(selected_parents), 10)
  expect_equal(selected_parents[1, ], population[1, ])
})
