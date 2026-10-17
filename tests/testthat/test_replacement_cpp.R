test_that("replacement_cpp returns a population of the original size", {
    set.seed(1)
    population <- matrix(rnorm(40), nrow = 4, ncol = 10)
    offspring <- matrix(rnorm(30), nrow = 3, ncol = 10)
    fitness <- matrix(1.0, nrow = 4, ncol = 2)
    offspring_fitness <- matrix(2.0, nrow = 3, ncol = 2)

    updated <- BioGA::replacement_cpp(population, offspring, fitness, offspring_fitness, 2)

    expect_equal(dim(updated), dim(population))
})

test_that("replacement_cpp admits the best offspring", {
    set.seed(1)
    population <- matrix(rnorm(40), nrow = 4, ncol = 10)
    offspring <- matrix(rnorm(30), nrow = 3, ncol = 10)
    fitness <- matrix(c(10, 10, 9, 9, 8, 8, 7, 7), nrow = 4, byrow = TRUE)
    offspring_fitness <- matrix(c(0, 0, 1, 1, 2, 2), nrow = 3, byrow = TRUE)

    updated <- BioGA::replacement_cpp(population, offspring, fitness, offspring_fitness, 2)

    expect_true(any(apply(updated, 1, function(row) all(row == offspring[1, ]))))
})

test_that("replacement_cpp with zero returns the population unchanged", {
    population <- matrix(rnorm(40), nrow = 4, ncol = 10)
    offspring <- matrix(rnorm(30), nrow = 3, ncol = 10)
    fitness <- matrix(1.0, nrow = 4, ncol = 2)
    offspring_fitness <- matrix(2.0, nrow = 3, ncol = 2)

    expect_equal(BioGA::replacement_cpp(population, offspring, fitness, offspring_fitness, 0),
                 population)
})

test_that("replacement_cpp handles a single individual without out-of-bounds", {
    set.seed(1)
    population <- matrix(rnorm(5), nrow = 1)
    offspring <- matrix(rnorm(5), nrow = 1)
    fitness <- matrix(0.0, nrow = 1, ncol = 2)
    offspring_fitness <- matrix(0.0, nrow = 1, ncol = 2)

    updated <- BioGA::replacement_cpp(population, offspring, fitness, offspring_fitness, 1)
    expect_equal(dim(updated), c(1, 5))
})
