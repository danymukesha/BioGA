# Internal: coerce a genomic input to a numeric matrix. Accepts a plain
# numeric matrix or a SummarizedExperiment, from which the first assay is used.
.as_genomic_matrix <- function(genomic_data) {
  if (is(genomic_data, "SummarizedExperiment")) {
    genomic_data <- assay(genomic_data, 1)
  }
  if (!is.matrix(genomic_data) || !is.numeric(genomic_data)) {
    stop("genomic_data must be a numeric matrix or a SummarizedExperiment",
         call. = FALSE)
  }
  genomic_data
}

#' Run the BioGA genetic algorithm
#'
#' User-facing wrapper around [bioga_main_cpp()] with sensible defaults and
#' input validation. The algorithm minimises two objectives; in `"binary"`
#' mode pass class `labels` so that the first objective measures the predictive
#' signal of the selected genes.
#'
#' @param genomic_data Numeric matrix of genomic data (rows: genes,
#'   columns: samples).
#' @param population_size Number of individuals in the population.
#' @param num_generations Number of generations to run.
#' @param crossover_rate Probability of crossover.
#' @param mutation_rate Base probability of mutation.
#' @param num_parents Number of parents selected per generation.
#' @param num_offspring Number of offspring generated per generation.
#' @param num_to_replace Maximum number of offspring admitted per generation.
#' @param weights Numeric vector of objective weights (length 1 or 2).
#' @param labels Optional class labels, one per sample (integer, numeric,
#'   factor, character or logical). Required when
#'   `representation = "binary"`.
#' @param clusters Optional vector of gene cluster assignments.
#' @param network Optional symmetric gene network adjacency matrix.
#' @param seed Optional integer random seed for reproducibility.
#' @param eta_c SBX distribution index (default: 20).
#' @param representation Either `"real"` or `"binary"`.
#'
#' @return A list with `population`, `fitness`, `history`, `best`,
#'   `best_fitness`, `pareto_population` and `pareto_fitness`.
#' @seealso [bioga_fitness()], [bioga_main_cpp()]
#' @examples
#' genomic_data <- matrix(rnorm(200), nrow = 20, ncol = 10)
#' result <- bioga(genomic_data, population_size = 20, num_generations = 5)
#' @export
bioga <- function(genomic_data,
                  population_size = 50,
                  num_generations = 20,
                  crossover_rate = 0.9,
                  mutation_rate = 0.1,
                  num_parents = 20,
                  num_offspring = 20,
                  num_to_replace = 10,
                  weights = c(1.0, 0.5),
                  labels = NULL,
                  clusters = NULL,
                  network = NULL,
                  seed = NULL,
                  eta_c = 20.0,
                  representation = c("real", "binary")) {
  representation <- match.arg(representation)
  genomic_data <- .as_genomic_matrix(genomic_data)
  if (representation == "binary" && is.null(labels)) {
    stop("labels must be provided when representation = 'binary'", call. = FALSE)
  }

  bioga_main_cpp(
    genomic_data = genomic_data,
    population_size = as.integer(population_size),
    num_generations = as.integer(num_generations),
    crossover_rate = as.numeric(crossover_rate),
    eta_c = as.numeric(eta_c),
    mutation_rate = as.numeric(mutation_rate),
    num_parents = as.integer(num_parents),
    num_offspring = as.integer(num_offspring),
    num_to_replace = as.integer(num_to_replace),
    weights = as.numeric(weights),
    seed = seed,
    clusters = clusters,
    network = network,
    labels = labels,
    representation = representation
  )
}

#' Evaluate BioGA fitness
#'
#' User-facing wrapper around [evaluate_fitness_cpp()].
#'
#' @param genomic_data Numeric matrix of genomic data (rows: genes,
#'   columns: samples).
#' @param population Numeric matrix of individuals (rows: individuals,
#'   columns: genes).
#' @param weights Numeric vector of objective weights (length 1 or 2).
#' @param labels Optional class labels, one per sample (integer, numeric,
#'   factor, character or logical). Required when
#'   `representation = "binary"`.
#' @param representation Either `"real"` or `"binary"`.
#'
#' @return Numeric matrix of fitness scores (rows: individuals,
#'   columns: objectives).
#' @seealso [bioga()], [evaluate_fitness_cpp()]
#' @examples
#' genomic_data <- matrix(rnorm(200), nrow = 20, ncol = 10)
#' population <- initialize_population_cpp(genomic_data, population_size = 5)
#' bioga_fitness(genomic_data, population)
#' @export
bioga_fitness <- function(genomic_data,
                          population,
                          weights = c(1.0, 0.5),
                          labels = NULL,
                          representation = c("real", "binary")) {
  representation <- match.arg(representation)
  genomic_data <- .as_genomic_matrix(genomic_data)
  if (representation == "binary" && is.null(labels)) {
    stop("labels must be provided when representation = 'binary'", call. = FALSE)
  }
  evaluate_fitness_cpp(
    genomic_data = genomic_data,
    population = population,
    weights = as.numeric(weights),
    labels = labels,
    representation = representation
  )
}
