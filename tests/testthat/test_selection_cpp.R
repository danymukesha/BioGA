test_that("selection_cpp returns a matrix with correct dimensions", {
  # Load The PACKAGE
  library(BioGA)
  
  # example of usage
  genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
  population <- BioGA::initialize_population_cpp(genomic_data,
                                                 population_size = 5)
  fitness <- BioGA::evaluate_fitness_cpp(genomic_data, population)
  selected_population <- BioGA::selection_cpp(population, fitness, 
                                              num_parents = 2)
  
  # Check dimensions
  expect_equal(2, nrow(selected_population))
  
})


test_that("selection_cpp returns a matrix with correct dimensions", {
  # Load The PACKAGE
  library(BioGA)
  
  # example of usage
  genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
  population <- BioGA::initialize_population_cpp(genomic_data,
                                                 population_size = 5)
  fitness <- BioGA::evaluate_fitness_cpp(genomic_data, population)
  selected_population <- BioGA::selection_cpp(population, fitness, 
                                              num_parents = 2)
  
  testthat::expect_no_warning(BioGA::evaluate_fitness_cpp(genomic_data, 
                                                          population),
                              message = "bananas")
  testthat::expect_no_error(BioGA::evaluate_fitness_cpp(genomic_data, 
                                                        population),
                            message = "bananas")
  
})
