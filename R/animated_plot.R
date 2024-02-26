#' Animate Fitness Change Over Generations
#'
#' Animate the change in fitness values over generations.
#'
#' @param fitness_history A list containing fitness values for each generation.
#'
#' @return Animation of fitness history
#'
#' @examples
#' # example of usage
#' fitness_history <- list(c(10, 8, 6, 4, 2), c(9, 7, 5, 3, 1))
#' animate_fitness_history(fitness_history)
#'
#' @export
animate_fitness_history <- function(fitness_history) {
  library(animation)

  # Function to generate each frame of the animation
  animate_plot <- function(iter) {
    for (i in seq_len(iter)) {
      temp <- data.frame(Generation = c(seq_len(i), seq_len(i)),
                         Variable = c(rep("mean", i), rep("best", i)),
                         Fitness = c(-fitness_history[[1]][seq_len(i)], -fitness_history[[2]][seq_len(i)]))
      pl <- ggplot2::ggplot(temp, ggplot2::aes(x = Generation, y = Fitness, group = Variable, colour = Variable)) +
        ggplot2::geom_line() +
        ggplot2::labs(x = "Generation", y = "Fitness") +
        ggplot2::ggtitle("Fitness Change Over Generations") +
        ggplot2::scale_x_continuous(limits = c(0, iter)) +
        ggplot2::scale_y_continuous(limits = c(0, 110)) +
        ggplot2::geom_hline(yintercept = max(temp$Fitness), linetype = "dashed") +
        ggplot2::annotate("text", x = 1, y = max(temp$Fitness) + 2, hjust = 0, size = 3, color = "black",
                          label = paste("Best solution:", max(temp$Fitness))) +
        ggplot2::scale_colour_brewer(palette = "Set1")

      print(pl)
    }
  }

  # Save the animation
  ani.options(interval = 0.1)
  saveGIF(animate_plot(length(fitness_history[[1]])), interval = 0.1, outdir = getwd())
}
