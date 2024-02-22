#' Plot Population Distribution
#'
#' Plot the distribution of individuals in the population.
#'
#' @param population A numeric matrix containing the population data.
#'
#' @return Plot of population
#'
#' @import graphics
#'
#' @examples
#' # example of usage
#' population <- matrix(runif(100), nrow = 10, ncol = 10)
#' plot_population(population)
#'
#' @export
plot_population <- function(population) {
  par(mfrow=c(1,2))
  boxplot(population, main = "Boxplot of Population")
  hist(population, main = "Histogram of Population")
}
