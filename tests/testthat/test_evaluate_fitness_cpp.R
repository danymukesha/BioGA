test_that("evaluate_fitness_cpp (real) returns two minimised objectives", {
    genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
    population <- BioGA::initialize_population_cpp(genomic_data, population_size = 5)
    weights <- c(1.0, 0.5)

    fitness <- BioGA::evaluate_fitness_cpp(genomic_data, population, weights)

    expect_equal(nrow(fitness), 5)
    expect_equal(ncol(fitness), 2)
    expect_true(all(fitness >= 0))
    expect_true(all(is.finite(fitness)))
})

test_that("evaluate_fitness_cpp accepts a single weight and recycles it", {
    genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
    population <- BioGA::initialize_population_cpp(genomic_data, population_size = 5)
    fitness <- BioGA::evaluate_fitness_cpp(genomic_data, population, weights = c(1.0))
    expect_equal(ncol(fitness), 2)
})

test_that("evaluate_fitness_cpp (binary) requires labels and scores sparsity", {
    genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
    labels <- rep(c(1L, 2L), each = 5)
    population <- BioGA::initialize_population_cpp(genomic_data, population_size = 5,
                                                   representation = "binary")

    expect_error(BioGA::evaluate_fitness_cpp(genomic_data, population, c(1, 1),
                                             representation = "binary"),
                 "labels")

    fitness <- BioGA::evaluate_fitness_cpp(genomic_data, population, c(1, 1),
                                           labels = labels, representation = "binary")
    expect_equal(ncol(fitness), 2)
    expect_true(all(fitness[, 1] <= 1e-12))
    expect_true(all(fitness[, 2] >= 0 & fitness[, 2] <= 1))
})

test_that("evaluate_fitness_cpp (binary) accepts character and factor labels", {
    genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
    population <- BioGA::initialize_population_cpp(genomic_data, population_size = 5,
                                                   representation = "binary")

    lab_char <- rep(c("Control", "Case"), each = 5)
    lab_fac <- factor(lab_char)

    fit_char <- BioGA::evaluate_fitness_cpp(genomic_data, population, c(1, 1),
                                            labels = lab_char, representation = "binary")
    fit_fac <- BioGA::evaluate_fitness_cpp(genomic_data, population, c(1, 1),
                                           labels = lab_fac, representation = "binary")
    expect_equal(dim(fit_char), c(5, 2))
    expect_equal(fit_char, fit_fac)
})

test_that("evaluate_fitness_cpp handles minimal input", {
    genomic_data <- matrix(1.0, nrow = 1, ncol = 1)
    population <- matrix(1.0, nrow = 1, ncol = 1)
    fitness <- BioGA::evaluate_fitness_cpp(genomic_data, population, weights = c(1.0))
    expect_equal(dim(fitness), c(1, 2))
    expect_true(is.numeric(fitness[1, 1]))
})
