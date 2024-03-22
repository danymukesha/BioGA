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
