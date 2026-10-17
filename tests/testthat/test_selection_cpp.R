test_that("selection_cpp returns matrix with correct dimensions", {
    genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
    population <- BioGA::initialize_population_cpp(genomic_data, population_size = 5)
    fitness <- BioGA::evaluate_fitness_cpp(genomic_data, population, c(1.0, 0.5))

    selected_parents <- BioGA::selection_cpp(population, fitness, num_parents = 2)

    expect_equal(nrow(selected_parents), 2)
    expect_equal(ncol(selected_parents), nrow(genomic_data))

    for (i in seq_len(nrow(selected_parents))) {
        expect_true(any(apply(population, 1, function(row) all(row == selected_parents[i, ]))))
    }
})

test_that("selection_cpp prefers non-dominated individuals", {
    set.seed(1)
    population <- matrix(rnorm(40), nrow = 4, ncol = 10)
    fitness <- rbind(c(0, 0), c(1, 1), c(2, 2), c(3, 3))

    sel <- replicate(200, as.numeric(BioGA::selection_cpp(population, fitness, 1)))
    best_count <- sum(apply(sel, 2, function(col) all(col == population[1, ])))
    worst_count <- sum(apply(sel, 2, function(col) all(col == population[4, ])))

    expect_gt(best_count, worst_count)
})

test_that("selection_cpp handles equal fitness without error", {
    population <- matrix(0, nrow = 4, ncol = 3)
    fitness <- matrix(1, nrow = 4, ncol = 2)
    selected_parents <- BioGA::selection_cpp(population, fitness, 2)
    expect_equal(nrow(selected_parents), 2)
})

test_that("selection_cpp handles single individual edge case", {
    population <- matrix(rnorm(10), nrow = 1, ncol = 10)
    fitness <- matrix(1.0, nrow = 1, ncol = 2)
    selected_parents <- BioGA::selection_cpp(population, fitness, 1)
    expect_equal(nrow(selected_parents), 1)
    expect_equal(selected_parents[1, ], population[1, ])
})
