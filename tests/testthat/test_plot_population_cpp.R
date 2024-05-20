test_that("crossover_cpp not produce error or warning", {
  # Load The PACKAGE
  library(BioGA)
  
  # example of usage
  population <- matrix(runif(100), nrow = 10, ncol = 10)
  BioGA::plot_population(population)
  
  testthat::expect_no_warning(BioGA::plot_population(population),
                              message = "bananas")
})