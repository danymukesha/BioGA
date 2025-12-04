test_that("BioGA package loads and exports all functions", {
    expect_no_error(library(BioGA), message = "BioGA package should load without errors")

    expected_functions <- c(
        "initialize_population_cpp", "evaluate_fitness_cpp",
        "selection_cpp", "crossover_cpp", "mutation_cpp",
        "replacement_cpp", "bioga_main_cpp"
    )
    exported_functions <- ls("package:BioGA")
    expect_true(all(expected_functions %in% exported_functions), TRUE)
})

test_that("BioGA handles typical workflow without errors", {
    genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
    population_size <- 5
    num_generations <- 5
    crossover_rate <- 0.9
    mutation_rate <- 0.1
    num_parents <- 2
    num_offspring <- 2
    num_to_replace <- 1
    weights <- c(1.0, 0.5)

    expect_no_error(
        {
            population <- BioGA::initialize_population_cpp(genomic_data, population_size)
            fitness <- BioGA::evaluate_fitness_cpp(genomic_data, population, weights)
            selected_parents <- BioGA::selection_cpp(population, fitness, num_parents)
            offspring <- BioGA::crossover_cpp(selected_parents, num_offspring, crossover_rate)
            mutated_offspring <- BioGA::mutation_cpp(offspring, mutation_rate, 1, num_generations)
            offspring_fitness <- BioGA::evaluate_fitness_cpp(genomic_data, mutated_offspring, weights)
            updated_population <- BioGA::replacement_cpp(
                population, mutated_offspring,
                fitness, offspring_fitness, num_to_replace
            )
        },
        message = "Full BioGA workflow should run without errors"
    )
})
