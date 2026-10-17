test_that("mutation_cpp (real) preserves dimensions and mutates", {
    set.seed(1)
    genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
    population <- BioGA::initialize_population_cpp(genomic_data, population_size = 5)
    fitness <- BioGA::evaluate_fitness_cpp(genomic_data, population, c(1.0, 0.5))
    selected_parents <- BioGA::selection_cpp(population, fitness, num_parents = 2)
    offspring <- BioGA::crossover_cpp(selected_parents, offspring_size = 2)

    mutated <- BioGA::mutation_cpp(offspring, mutation_rate = 0.5, iteration = 1,
                                   max_iterations = 100)

    expect_equal(dim(mutated), dim(offspring))
    expect_true(any(mutated != offspring))
})

test_that("mutation_cpp with zero rate returns the input unchanged", {
    offspring <- matrix(rnorm(20), nrow = 2, ncol = 10)
    expect_equal(BioGA::mutation_cpp(offspring, 0.0, 1, 100), offspring)
})

test_that("mutation_cpp (binary) flips bits and keeps values in {0, 1}", {
    set.seed(1)
    offspring <- matrix(rep(c(0, 1), 10), nrow = 2)
    mutated <- BioGA::mutation_cpp(offspring, mutation_rate = 1.0, iteration = 1,
                                   max_iterations = 10, representation = "binary")
    expect_true(all(mutated %in% c(0, 1)))
    expect_true(all(mutated == 1 - offspring))
})

test_that("mutation_cpp validates the network dimensions", {
    offspring <- matrix(rnorm(20), nrow = 2, ncol = 10)
    bad_network <- matrix(0, nrow = 5, ncol = 5)
    expect_error(BioGA::mutation_cpp(offspring, 0.1, 1, 100, network = bad_network),
                 "square matrix")
})

test_that("mutation_cpp damps hubs when a network is supplied", {
    set.seed(1)
    offspring <- matrix(0, nrow = 1, ncol = 3)
    network <- matrix(0, nrow = 3, ncol = 3)
    network[1, 2] <- network[2, 1] <- 1
    network[1, 3] <- network[3, 1] <- 1
    mutated <- BioGA::mutation_cpp(offspring, mutation_rate = 1.0, iteration = 1,
                                   max_iterations = 1, network = network)
    expect_equal(mutated[1, 1], 0, tolerance = 1e-12)
    expect_true(abs(mutated[1, 2]) > 0)
})
