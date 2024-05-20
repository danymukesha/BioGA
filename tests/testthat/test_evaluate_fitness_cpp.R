test_that("evalutate_cpp not produce error or warning", {
  # example of usage
  genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
  population <- BioGA::initialize_population_cpp(genomic_data,
                                                 population_size = 5)
  
  output <- tryCatch(
    BioGA::evaluate_fitness_cpp(genomic_data, population), 
    error = function(e) return(NA)
  )
  
  testthat::expect_no_warning(BioGA::evaluate_fitness_cpp(genomic_data, 
                                                          population),
                              message = "bananas")
  testthat::expect_no_error(BioGA::evaluate_fitness_cpp(genomic_data, 
                                                          population),
                              message = "bananas")
})

