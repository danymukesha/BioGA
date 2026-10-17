#include <Rcpp.h>
#include <vector>
#include <limits>
#include "BioGA.h"
#include "nsga2.h"
using namespace Rcpp;

//' Main genetic algorithm loop for genomic data optimisation
//'
//' Runs a generational multi-objective genetic algorithm combining
//' initialization, fitness evaluation, NSGA-II selection, crossover,
//' mutation and NSGA-II replacement.
//'
//' @param genomic_data Numeric matrix of genomic data (rows: genes,
//' columns: samples).
//' @param population_size Number of individuals in the population.
//' @param num_generations Number of generations to run.
//' @param crossover_rate Probability of crossover.
//' @param eta_c SBX distribution index (real representation).
//' @param mutation_rate Base probability of mutation.
//' @param num_parents Number of parents selected per generation.
//' @param num_offspring Number of offspring generated per generation.
//' @param num_to_replace Maximum number of offspring admitted per generation.
//' @param weights Numeric vector of objective weights.
//' @param seed Optional random seed for reproducibility.
//' @param clusters Optional vector of gene cluster assignments.
//' @param network Optional symmetric gene network adjacency matrix.
//' @param labels Optional class labels (required when
//' `representation = "binary"`). Accepted as an integer, numeric, factor,
//' character or logical vector.
//' @param representation Either `"real"` or `"binary"`.
//' @return A list with `population` (final population), `fitness` (final
//' fitness matrix), `history` (best first objective per generation), `best`
//' (best individual), `best_fitness`, `pareto_population` and
//' `pareto_fitness` (the first non-dominated front of the final population).
//' @family main genetic algorithm loop
//' @seealso [initialize_population_cpp()], [evaluate_fitness_cpp()],
//'   [selection_cpp()], [crossover_cpp()], [mutation_cpp()],
//'   [replacement_cpp()]
//' @references
//'   Deb K, Pratap A, Agarwal S, Meyarivan T (2002).
//'   "A fast and elitist multiobjective genetic algorithm: NSGA-II."
//'   _IEEE Transactions on Evolutionary Computation_, 6(2), 182-197.
//'   \doi{10.1109/4235.996017}
//' @examples
//' genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
//' result <- BioGA::bioga_main_cpp(genomic_data, population_size = 50,
//'   num_generations = 10, crossover_rate = 0.9, eta_c = 20.0,
//'   mutation_rate = 0.1, num_parents = 20, num_offspring = 20,
//'   num_to_replace = 10, weights = c(1.0, 0.5), seed = 123)
//' @export
// [[Rcpp::export]]
List bioga_main_cpp(const NumericMatrix& genomic_data,
                    int population_size,
                    int num_generations,
                    double crossover_rate,
                    double eta_c,
                    double mutation_rate,
                    int num_parents,
                    int num_offspring,
                    int num_to_replace,
                    const NumericVector& weights,
                    Nullable<int> seed = R_NilValue,
                    Nullable<IntegerVector> clusters = R_NilValue,
                    Nullable<NumericMatrix> network = R_NilValue,
                    SEXP labels = R_NilValue,
                    std::string representation = "real") {
  if (genomic_data.nrow() == 0 || genomic_data.ncol() == 0)
    stop("genomic_data must have at least 1 row and 1 column");
  if (population_size < 1)
    stop("population_size must be >= 1");
  if (num_generations < 1)
    stop("num_generations must be >= 1");
  if (crossover_rate < 0.0 || crossover_rate > 1.0)
    stop("crossover_rate must be in [0, 1]");
  if (eta_c <= 0.0)
    stop("eta_c must be > 0");
  if (mutation_rate < 0.0 || mutation_rate > 1.0)
    stop("mutation_rate must be in [0, 1]");
  if (num_parents < 1 || num_parents > population_size)
    stop("num_parents must be in [1, population_size]");
  if (num_offspring < 1)
    stop("num_offspring must be >= 1");
  if (num_to_replace < 0)
    stop("num_to_replace must be >= 0");
  if (weights.length() == 0)
    stop("weights must not be empty");
  if (representation != "real" && representation != "binary")
    stop("representation must be either 'real' or 'binary'");

  NumericMatrix population =
    initialize_population_cpp(genomic_data, population_size, seed, clusters, representation);

  NumericVector history(num_generations);

  for (int generation = 0; generation < num_generations; ++generation) {
    NumericMatrix fitness =
      evaluate_fitness_cpp(genomic_data, population, weights, labels, representation);

    history[generation] = fitness(0, 0);
    for (int i = 1; i < population_size; ++i) {
      if (fitness(i, 0) < history[generation]) history[generation] = fitness(i, 0);
    }

    NumericMatrix parents = selection_cpp(population, fitness, num_parents);
    NumericMatrix offspring =
      crossover_cpp(parents, num_offspring, crossover_rate, eta_c, representation);
    NumericMatrix mutated =
      mutation_cpp(offspring, mutation_rate, generation, num_generations, network, representation);
    NumericMatrix offspring_fitness =
      evaluate_fitness_cpp(genomic_data, mutated, weights, labels, representation);

    population = replacement_cpp(population, mutated, fitness, offspring_fitness, num_to_replace);
  }

  NumericMatrix final_fitness =
    evaluate_fitness_cpp(genomic_data, population, weights, labels, representation);

  int best_idx = 0;
  for (int i = 1; i < population_size; ++i) {
    if (final_fitness(i, 0) < final_fitness(best_idx, 0)) best_idx = i;
  }

  std::vector<int> rank;
  std::vector<std::vector<int> > fronts;
  bioga::fast_nondominated_sort(final_fitness, rank, fronts);
  std::vector<int> first_front;
  if (!fronts.empty()) first_front = fronts[0];

  NumericMatrix pareto_population(first_front.size(), population.ncol());
  NumericMatrix pareto_fitness(first_front.size(), final_fitness.ncol());
  for (std::size_t a = 0; a < first_front.size(); ++a) {
    for (int j = 0; j < population.ncol(); ++j) {
      pareto_population(a, j) = population(first_front[a], j);
    }
    for (int k = 0; k < final_fitness.ncol(); ++k) {
      pareto_fitness(a, k) = final_fitness(first_front[a], k);
    }
  }

  NumericMatrix best(1, population.ncol());
  for (int j = 0; j < population.ncol(); ++j) best(0, j) = population(best_idx, j);

  NumericVector best_fitness(final_fitness.ncol());
  for (int k = 0; k < final_fitness.ncol(); ++k) best_fitness[k] = final_fitness(best_idx, k);

  return List::create(
    Named("population") = population,
    Named("fitness") = final_fitness,
    Named("history") = history,
    Named("best") = best,
    Named("best_fitness") = best_fitness,
    Named("pareto_population") = pareto_population,
    Named("pareto_fitness") = pareto_fitness);
}
