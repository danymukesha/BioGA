#' Plot Fitness Change Over Generations
#'
#' Plot the change in fitness values over generations.
#'
#' @param fitness_history A list containing fitness values for each generation.
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
  generations <- rep(seq_along(fitness_history), vapply(fitness_history, length, integer(1)))

  # Create data frame
  df <- data.frame(Generation = generations, Fitness = fitness_values)

  # Plot
  ggplot2::ggplot(df, ggplot2::aes_string(x = "Generation", y = "Fitness")) +
    geom_line() +
    labs(x = "Generation", y = "Fitness") +
    ggplot2::ggtitle("Fitness Change Over Generations")
}

