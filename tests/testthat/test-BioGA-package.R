test_that("BioGA package loads and exports all functions", {
    expected_functions <- c(
        "bioga", "bioga_fitness", "initialize_population_cpp",
        "evaluate_fitness_cpp", "selection_cpp", "crossover_cpp",
        "mutation_cpp", "replacement_cpp", "bioga_main_cpp"
    )
    exported_functions <- getNamespaceExports("BioGA")
    expect_true(all(expected_functions %in% exported_functions))
})

test_that("BioGA handles a typical workflow without errors", {
    genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
    weights <- c(1.0, 0.5)

    expect_no_error({
        population <- BioGA::initialize_population_cpp(genomic_data, 5)
        fitness <- BioGA::evaluate_fitness_cpp(genomic_data, population, weights)
        selected_parents <- BioGA::selection_cpp(population, fitness, num_parents = 2)
        offspring <- BioGA::crossover_cpp(selected_parents, offspring_size = 2)
        mutated_offspring <- BioGA::mutation_cpp(offspring, mutation_rate = 0.1,
                                                 iteration = 1, max_iterations = 5)
        offspring_fitness <- BioGA::evaluate_fitness_cpp(genomic_data, mutated_offspring, weights)
        updated_population <- BioGA::replacement_cpp(
            population, mutated_offspring, fitness, offspring_fitness,
            num_to_replace = 1
        )
    })
})

test_that("bioga() provides a validated interface", {
    genomic_data <- matrix(rnorm(200), nrow = 20, ncol = 10)
    result <- BioGA::bioga(genomic_data, population_size = 10, num_generations = 3,
                           num_parents = 5, num_offspring = 5, num_to_replace = 2,
                           seed = 1)
    expect_named(result, c("population", "fitness", "history", "best",
                           "best_fitness", "pareto_population", "pareto_fitness"))

    expect_error(BioGA::bioga(genomic_data, representation = "binary"),
                 "labels")
    expect_error(BioGA::bioga("not a matrix"), "numeric matrix")
})

test_that("bioga() accepts a SummarizedExperiment", {
    skip_if_not_installed("SummarizedExperiment")
    genomic_data <- matrix(rnorm(200), nrow = 20, ncol = 10)
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(expr = genomic_data))
    result <- BioGA::bioga(se, population_size = 10, num_generations = 3,
                           num_parents = 5, num_offspring = 5, num_to_replace = 2,
                           seed = 1)
    expect_equal(ncol(result$population), nrow(genomic_data))
})
