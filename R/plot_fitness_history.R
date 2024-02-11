#' Plot Fitness Change Over Generations
#'
#' Plot the change in fitness values over generations.
#'
#' @param fitness_history A list containing fitness values for each generation.
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
  generations <- rep(seq_along(fitness_history), sapply(fitness_history, length))

  # Create data frame
  df <- data.frame(Generation = generations, Fitness = fitness_values)

  # Plot
  library(ggplot2)
  ggplot(df, aes(x = Generation, y = Fitness)) +
    geom_line() +
    labs(x = "Generation", y = "Fitness") +
    ggtitle("Fitness Change Over Generations")
}

