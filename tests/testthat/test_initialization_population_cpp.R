test_that("initialize_population_cpp returns a matrix with correct dimensions", {
  # Load your package
  library(BioGA)

  # Generate some test data
  genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
  population_size <- 5

  # Call the function
  population <- initialize_population_cpp(genomic_data, population_size)

  # Check dimensions
  expect_equal(nrow(population), population_size)
  expect_equal(ncol(population), nrow(genomic_data))
})

test_that("initialize_population_cpp returns a matrix with values from genomic data", {
  # Load THE PACKAGE
  library(BioGA)

  # Generate some test data
  genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
  population_size <- 5

  # Call the function
  population <- initialize_population_cpp(genomic_data, population_size)

  # Check if all values in the population are present in genomic data
  for (i in 1:nrow(population)) {
    for (j in 1:ncol(population)) {
      expect_true(population[i, j] %in% genomic_data[j, ])
    }
  }
})
