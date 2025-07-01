library(BioGA)
library(testthat)

test_that("replacement_cpp returns matrix with correct dimensions and preserves elite", {
    genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
    population <- BioGA::initialize_population_cpp(genomic_data, population_size = 5)
    fitness <- BioGA::evaluate_fitness_cpp(genomic_data, population, c(1.0, 0.5))
    selected_parents <- BioGA::selection_cpp(population, fitness, num_parents = 2)
    offspring <- BioGA::crossover_cpp(selected_parents, offspring_size = 2)
    offspring_fitness <- BioGA::evaluate_fitness_cpp(genomic_data, offspring, c(1.0, 0.5))
    num_to_replace <- 1
    
    updated_population <- BioGA::replacement_cpp(population, offspring, fitness, offspring_fitness, num_to_replace)
    
    expect_equal(nrow(updated_population), nrow(population))
    expect_equal(ncol(updated_population), ncol(population))
    
    best_idx <- which.min(fitness[, 1])
    expect_true(any(apply(updated_population, 1, function(row) all(row == population[best_idx, ]))))
    
    expect_no_warning(BioGA::replacement_cpp(population, offspring, fitness, offspring_fitness, num_to_replace),
                      message = "replacement_cpp should not produce warnings")
    expect_no_error(BioGA::replacement_cpp(population, offspring, fitness, offspring_fitness, num_to_replace),
                    message = "replacement_cpp should not produce errors")
})

test_that("replacement_cpp handles edge cases", {
    population <- matrix(rnorm(10), nrow = 1, ncol = 10)
    offspring <- matrix(rnorm(10), nrow = 1, ncol = 10)
    fitness <- matrix(1.0, nrow = 1, ncol = 2)
    offspring_fitness <- matrix(2.0, nrow = 1, ncol = 2)
    updated_population <- BioGA::replacement_cpp(population, offspring, fitness, offspring_fitness, 0)
    expect_equal(updated_population, population)
})
