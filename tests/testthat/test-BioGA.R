test_that("bioga_main_cpp runs and returns valid output", {
    set.seed(1)

    # Small synthetic genomic dataset
    genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)

    result <- BioGA::bioga_main_cpp(
        genomic_data        = genomic_data,
        population_size     = 20,
        num_generations     = 5,
        crossover_rate      = 0.9,
        eta_c               = 20.0,
        mutation_rate       = 0.1,
        num_parents         = 10,
        num_offspring       = 10,
        num_to_replace      = 5,
        weights             = c(1.0, 0.5),
        seed                = 123
    )

    # ---- Structure ----
    expect_type(result, "list")
    expect_named(result, c("population", "fitness"))

    # ---- Population ----
    expect_true(is.matrix(result$population))
    expect_equal(nrow(result$population), 20)

    # ---- Fitness ----
    expect_true(is.matrix(result$fitness))
    expect_equal(nrow(result$fitness), 20)
    expect_equal(ncol(result$fitness), length(c(1.0, 0.5)))

    # ---- Fitness sanity ----
    expect_true(all(is.finite(result$fitness)))
    expect_true(all(is.numeric(result$fitness)))
})


test_that("bioga_main_cpp is reproducible with a fixed seed", {
    genomic_data <- matrix(rnorm(100), nrow = 10)

    res1 <- BioGA::bioga_main_cpp(
        genomic_data, 20, 5, 0.9, 20.0, 0.1,
        10, 10, 5, c(1, 0.5),
        seed = 42
    )

    res2 <- BioGA::bioga_main_cpp(
        genomic_data, 20, 5, 0.9, 20.0, 0.1,
        10, 10, 5, c(1, 0.5),
        seed = 42
    )

    expect_equal(res1$population, res2$population)
    expect_equal(res1$fitness, res2$fitness)
})

test_that("bioga_main_cpp errors on invalid inputs", {
    genomic_data <- matrix(rnorm(20), nrow = 5)

    expect_error(
        BioGA::bioga_main_cpp(
            genomic_data,
            population_size = -1,
            num_generations = 5,
            crossover_rate  = 0.9,
            eta_c           = 20.0,
            mutation_rate   = 0.1,
            num_parents     = 5,
            num_offspring   = 5,
            num_to_replace  = 2,
            weights         = c(1, 0.5)
        )
    )
})
