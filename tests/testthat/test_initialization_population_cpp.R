library(BioGA)
library(testthat)

test_that("initialize_population_cpp returns matrix with correct dimensions", {
    genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
    population_size <- 5
    seed <- 123
    
    population <- BioGA::initialize_population_cpp(genomic_data, population_size, seed)
    
    expect_equal(nrow(population), population_size)
    expect_equal(ncol(population), nrow(genomic_data))
    
    population2 <- BioGA::initialize_population_cpp(genomic_data, population_size, seed)
    expect_equal(population, population2)
    
    expect_no_warning(BioGA::initialize_population_cpp(genomic_data, population_size, seed),
                      message = "initialize_population_cpp should not produce warnings")
    expect_no_error(BioGA::initialize_population_cpp(genomic_data, population_size, seed),
                    message = "initialize_population_cpp should not produce errors")
})

test_that("initialize_population_cpp handles clusters", {
    genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
    population_size <- 5
    clusters <- rep(1:2, each = 5) # Two clusters
    
    population <- BioGA::initialize_population_cpp(genomic_data, population_size, clusters = clusters)
    
    expect_equal(nrow(population), population_size)
    expect_equal(ncol(population), nrow(genomic_data))
    
    for (i in 1:nrow(population)) {
        for (j in 1:ncol(population)) {
            expect_true(population[i, j] %in% genomic_data[j, ])
        }
    }
})

test_that("initialize_population_cpp handles edge cases", {
    genomic_data <- matrix(1.0, nrow = 1, ncol = 1)
    population <- BioGA::initialize_population_cpp(genomic_data, 1)
    expect_equal(nrow(population), 1)
    expect_equal(ncol(population), 1)
    expect_equal(population[1, 1], 1.0)
})