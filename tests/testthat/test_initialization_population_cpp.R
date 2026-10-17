test_that("initialize_population_cpp returns matrix with correct dimensions", {
    genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
    population_size <- 5
    seed <- 123

    population <- BioGA::initialize_population_cpp(genomic_data, population_size, seed)

    expect_equal(nrow(population), population_size)
    expect_equal(ncol(population), nrow(genomic_data))

    population2 <- BioGA::initialize_population_cpp(genomic_data, population_size, seed)
    expect_equal(population, population2)
})

test_that("initialize_population_cpp (binary) returns indicators", {
    genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
    population <- BioGA::initialize_population_cpp(genomic_data, 5, seed = 1,
                                                   representation = "binary")
    expect_true(all(population %in% c(0, 1)))
})

test_that("initialize_population_cpp handles clusters", {
    genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
    clusters <- rep(1:2, each = 5)

    population <- BioGA::initialize_population_cpp(genomic_data, 5, clusters = clusters)

    expect_equal(dim(population), c(5, 10))
    for (i in seq_len(nrow(population))) {
        for (j in seq_len(ncol(population))) {
            expect_true(population[i, j] %in% genomic_data[j, ])
        }
    }
})

test_that("initialize_population_cpp validates cluster assignments", {
    genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
    expect_error(BioGA::initialize_population_cpp(genomic_data, 5, clusters = c(1, 2)),
                 "one assignment per gene")
    expect_error(BioGA::initialize_population_cpp(genomic_data, 5, clusters = c(-1, rep(1, 9))),
                 "non-negative")
})

test_that("initialize_population_cpp handles minimal input", {
    genomic_data <- matrix(1.0, nrow = 1, ncol = 1)
    population <- BioGA::initialize_population_cpp(genomic_data, 1)
    expect_equal(dim(population), c(1, 1))
    expect_equal(population[1, 1], 1.0)
})
