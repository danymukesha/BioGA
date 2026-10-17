test_that("bioga_main_cpp returns a well-formed, valid result", {
    genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
    population_size <- 10
    num_generations <- 5

    result <- BioGA::bioga_main_cpp(
        genomic_data,
        population_size = population_size,
        num_generations = num_generations,
        crossover_rate = 0.9,
        eta_c = 20.0,
        mutation_rate = 0.1,
        num_parents = 5,
        num_offspring = 5,
        num_to_replace = 3,
        weights = c(1.0, 0.5),
        seed = 42
    )

    expect_true(is.list(result))
    expect_named(result, c("population", "fitness", "history", "best",
                           "best_fitness", "pareto_population", "pareto_fitness"))
    expect_equal(dim(result$population), c(population_size, nrow(genomic_data)))
    expect_equal(dim(result$fitness), c(population_size, 2))
    expect_length(result$history, num_generations)
    expect_equal(ncol(result$best), nrow(genomic_data))
    expect_true(all(diff(result$history) <= 1e-8))
})

test_that("bioga_main_cpp supports binary representation with labels", {
    genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
    labels <- rep(c(1L, 2L), each = 5)

    result <- BioGA::bioga_main_cpp(
        genomic_data,
        population_size = 10,
        num_generations = 5,
        crossover_rate = 0.9,
        eta_c = 20.0,
        mutation_rate = 0.1,
        num_parents = 5,
        num_offspring = 5,
        num_to_replace = 3,
        weights = c(1.0, 0.5),
        seed = 42,
        labels = labels,
        representation = "binary"
    )

    expect_true(all(result$population %in% c(0, 1)))
    expect_true(all(result$fitness[, 2] >= 0 & result$fitness[, 2] <= 1))
})

test_that("bioga_main_cpp handles clusters and network inputs", {
    genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
    clusters <- rep(1:2, each = 5)
    network <- matrix(0, nrow = 10, ncol = 10)

    result <- BioGA::bioga_main_cpp(
        genomic_data,
        population_size = 5,
        num_generations = 3,
        crossover_rate = 0.9,
        eta_c = 20.0,
        mutation_rate = 0.1,
        num_parents = 2,
        num_offspring = 2,
        num_to_replace = 1,
        weights = c(1.0, 0.5),
        clusters = clusters,
        network = network
    )

    expect_equal(dim(result$population), c(5, 10))
    expect_equal(dim(result$fitness), c(5, 2))
})

test_that("bioga_main_cpp handles minimal input without hanging", {
    genomic_data <- matrix(1.0, nrow = 1, ncol = 1)

    result <- BioGA::bioga_main_cpp(
        genomic_data,
        population_size = 1,
        num_generations = 1,
        crossover_rate = 0.0,
        eta_c = 20.0,
        mutation_rate = 0.0,
        num_parents = 1,
        num_offspring = 1,
        num_to_replace = 0,
        weights = c(1.0)
    )

    expect_equal(dim(result$population), c(1, 1))
    expect_equal(dim(result$fitness), c(1, 2))
})
