#' Plot Fitness Change Over Generations
#'
#' Plot the change in fitness values over generations.
#'
#' @param fitness_history A list containing fitness values for
#' each generation.
#'
#' @return Plot of fitness history
#'
#' @examples
#' # example of usage
#' fitness_history <- list(c(10, 8, 6, 4, 2), c(9, 7, 5, 3, 1))
#' plot_fitness_history(fitness_history)
#'
#' @export
plot_fitness_history <- function(fitness_history) {
    # Extract fitness values
    fitness_values <- unlist(fitness_history)

    # Create generation index
    generations <-
        rep(
            seq_along(fitness_history),
            vapply(fitness_history, length, integer(1))
        )

    # Create data frame
    df <- data.frame(
        Generation = generations,
        Fitness = fitness_values
    )

    # Plot
    ggplot2::ggplot(df, ggplot2::aes_string(
        x = "Generation",
        y = "Fitness"
    )) +
        geom_line() +
        labs(x = "Generation", y = "Fitness") +
        ggplot2::ggtitle("Fitness Change Over Generations")
}

#' Plot Fitness Values
#'
#' Plot the fitness values of the population over generations.
#'
#' @param fitness_values A numeric vector containing fitness values.
#'
#' @return Plot of fitness
#'
#' @examples
#' # example of usage
#' fitness_values <- c(10, 8, 6, 4, 2)
#' plot_fitness(fitness_values)
#'
#' @export
plot_fitness <- function(fitness_values) {
    generations <- seq_along(fitness_values)
    plot(generations, fitness_values,
        type = "l",
        xlab = "Generations", ylab = "Fitness",
        main = "Fitness Values Over Generations"
    )
}

#' Plot Population Distribution
#'
#' Plot the distribution of individuals in the population.
#'
#' @param population A numeric matrix containing the population data.
#'
#' @return Plot of population
#'
#'
#' @examples
#' # example of usage
#' population <- matrix(runif(100), nrow = 10, ncol = 10)
#' plot_population(population)
#'
#' @export
plot_population <- function(population) {
    par(mfrow = c(1, 2))
    boxplot(population, main = "Boxplot of Population")
    hist(population, main = "Histogram of Population")
}
