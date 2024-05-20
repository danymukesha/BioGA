test_that("crossover_cpp not produce error or warning", {
  # Load The PACKAGE
  library(BioGA)
  
  # example of usage
  genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
  population <- BioGA::initialize_population_cpp(genomic_data,
                                                 population_size = 5)
  fitness <- BioGA::evaluate_fitness_cpp(genomic_data, population)
  selected_parents <- BioGA::selection_cpp(population, fitness,
                                           num_parents = 2)
  BioGA::crossover_cpp(selected_parents, offspring_size = 2)
  
  testthat::expect_no_warning(BioGA::crossover_cpp(selected_parents, 
                                                   offspring_size = 2),
                              message = "bananas")
  
  testthat::expect_no_error(BioGA::crossover_cpp(selected_parents, 
                                                   offspring_size = 2),
                              message = "bananas")
})