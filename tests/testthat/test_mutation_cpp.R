test_that("mutation_cpp not produce error or warning", {
  # example of usage
  genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
  population <- BioGA::initialize_population_cpp(genomic_data,
                                                 population_size = 5)
  fitness <- BioGA::evaluate_fitness_cpp(genomic_data, population)
  selected_parents <- BioGA::selection_cpp(population,
                                           fitness, num_parents = 2)
  offspring <- BioGA::crossover_cpp(selected_parents, offspring_size = 2)
  BioGA::mutation_cpp(offspring, mutation_rate = 0)
  
  testthat::expect_no_warning(BioGA::mutation_cpp(offspring, 
                                                  mutation_rate = 0),
                              message = "bananas")
  
  testthat::expect_no_error(BioGA::mutation_cpp(offspring, 
                                                mutation_rate = 0),
                            message = "bananas")
})