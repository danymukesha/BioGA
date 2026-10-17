test_that("crossover_cpp returns matrix with correct dimensions", {
    genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
    population <- BioGA::initialize_population_cpp(genomic_data, population_size = 5)
    fitness <- BioGA::evaluate_fitness_cpp(genomic_data, population, c(1.0, 0.5))
    selected_parents <- BioGA::selection_cpp(population, fitness, num_parents = 2)

    offspring <- BioGA::crossover_cpp(selected_parents, offspring_size = 2)

    expect_equal(nrow(offspring), 2)
    expect_equal(ncol(offspring), nrow(genomic_data))
    expect_true(all(is.finite(offspring)))
})

test_that("SBX crossover preserves the parents' mean in expectation", {
    set.seed(1)
    parents <- rbind(c(1, 2, 3), c(3, 4, 5))
    offspring <- BioGA::crossover_cpp(parents, offspring_size = 5000,
                                      crossover_rate = 1.0, eta_c = 2.0)
    expect_equal(colMeans(offspring), colMeans(parents), tolerance = 0.05)
})

test_that("binary crossover keeps values in {0, 1}", {
    set.seed(1)
    parents <- matrix(c(1, 0, 1, 1, 1, 0), nrow = 2, byrow = TRUE)
    offspring <- BioGA::crossover_cpp(parents, offspring_size = 50,
                                      crossover_rate = 1.0, representation = "binary")
    expect_true(all(offspring %in% c(0, 1)))
})

test_that("crossover_cpp does not hang with a single parent", {
    parents <- matrix(rnorm(5), nrow = 1)
    offspring <- BioGA::crossover_cpp(parents, offspring_size = 3, crossover_rate = 0.9)
    expect_equal(offspring, matrix(parents[1, ], nrow = 3, ncol = 5, byrow = TRUE))
})

test_that("crossover_cpp copies parent when crossover_rate is zero", {
    parents <- rbind(c(1, 2, 3), c(4, 5, 6))
    offspring <- BioGA::crossover_cpp(parents, offspring_size = 10, crossover_rate = 0.0)
    for (i in seq_len(nrow(offspring))) {
        expect_true(any(apply(parents, 1, function(row) all(row == offspring[i, ]))))
    }
})
