#' BioGA: Bioinformatics Genetic Algorithm
#'
#' @description
#' BioGA provides a genetic algorithm framework for analyzing and optimizing
#' high throughput genomic data. All core functions are implemented in C++
#' for performance via Rcpp and RcppParallel, with an easy-to-use R interface.
#'
#' @details
#' Two representations are supported. In `"real"` mode individuals are
#' continuous expression profiles; in `"binary"` mode individuals are gene
#' selection indicators and class `labels` are required. The workflow is:
#' \enumerate{
#'   \item \strong{Initialization:} \code{\link{initialize_population_cpp}()} —
#'     Creates an initial population, optionally guided by gene clusters.
#'   \item \strong{Fitness evaluation:} \code{\link{evaluate_fitness_cpp}()} —
#'     Evaluates two minimised objectives (data fidelity and parsimony in
#'     `"real"` mode; predictive signal and sparsity in `"binary"` mode),
#'     computed in parallel.
#'   \item \strong{Selection:} \code{\link{selection_cpp}()} —
#'     Selects parents using NSGA-II rank and crowding-distance tournament
#'     selection.
#'   \item \strong{Crossover:} \code{\link{crossover_cpp}()} —
#'     Generates offspring via simulated binary crossover (SBX) in `"real"`
#'     mode or uniform crossover in `"binary"` mode.
#'   \item \strong{Mutation:} \code{\link{mutation_cpp}()} —
#'     Applies adaptive Gaussian (or bit-flip) mutation with optional
#'     degree-based network damping.
#'   \item \strong{Replacement:} \code{\link{replacement_cpp}()} —
#'     Admits the best offspring by NSGA-II environmental selection.
#'   \item \strong{Main loop:} \code{\link{bioga_main_cpp}()} —
#'     Orchestrates the full pipeline; \code{\link{bioga}()} is the
#'     user-facing entry point.
#' }
#'
#' @section Citation:
#' If you use BioGA in your research, please cite:
#' \preformatted{
#' Mukesha D (2024). "Scalable and Interpretable Genetic Algorithm Framework
#' for High-Throughput Genomic Data." Zenodo. \doi{10.5281/zenodo.15801072}
#' }
#'
#' @author
#' \strong{Maintainer:} Dany Mukesha \email{danymukesha@gmail.com}
#'
#' @seealso
#' Useful links:
#' \itemize{
#'   \item \url{https://danymukesha.github.io/BioGA/}
#'   \item \url{https://bioconductor.org/packages/BioGA/}
#'   \item Report bugs at \url{https://github.com/danymukesha/BioGA/issues}
#' }
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom SummarizedExperiment SummarizedExperiment assay
#' @importFrom methods is
#' @importFrom Rcpp evalCpp sourceCpp
#' @import RcppParallel
#' @useDynLib BioGA, .registration = TRUE
## usethis namespace: end
NULL
