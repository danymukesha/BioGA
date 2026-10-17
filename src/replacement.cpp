#include <Rcpp.h>
#include <vector>
#include "BioGA.h"
#include "nsga2.h"
using namespace Rcpp;

//' Build the next population with NSGA-II environmental selection
//'
//' The `num_to_replace` best offspring (by NSGA-II rank and crowding distance)
//' replace the lowest-quality parents, keeping the best `population_size -
//' num_to_replace` parents. This preserves elitism while letting offspring
//' enter only on merit. If `num_to_replace <= 0` the population is unchanged.
//'
//' @param population Numeric matrix of current individuals.
//' @param offspring Numeric matrix of offspring.
//' @param fitness Numeric matrix of current population fitness scores.
//' @param offspring_fitness Numeric matrix of offspring fitness scores.
//' @param num_to_replace Maximum number of offspring to admit.
//' @return Numeric matrix of the updated population.
//' @family replacement
//' @seealso [initialize_population_cpp()], [evaluate_fitness_cpp()],
//'   [selection_cpp()], [crossover_cpp()], [mutation_cpp()],
//'   [bioga_main_cpp()]
//' @examples
//' genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
//' population <- BioGA::initialize_population_cpp(genomic_data, population_size = 5)
//' fitness <- BioGA::evaluate_fitness_cpp(genomic_data, population, c(1.0, 0.5))
//' offspring <- population
//' offspring_fitness <- fitness
//' BioGA::replacement_cpp(population, offspring, fitness, offspring_fitness,
//'   num_to_replace = 1)
//' @export
// [[Rcpp::export]]
NumericMatrix replacement_cpp(const NumericMatrix& population,
                              const NumericMatrix& offspring,
                              const NumericMatrix& fitness,
                              const NumericMatrix& offspring_fitness,
                              int num_to_replace) {
  if (population.nrow() == 0)
    stop("population must have at least 1 individual");
  if (offspring.nrow() == 0)
    stop("offspring must have at least 1 individual");
  if (population.ncol() != offspring.ncol())
    stop("population and offspring must have the same number of genes");
  if (fitness.nrow() != population.nrow())
    stop("fitness rows must match population rows");
  if (offspring_fitness.nrow() != offspring.nrow())
    stop("offspring_fitness rows must match offspring rows");
  if (fitness.ncol() != offspring_fitness.ncol())
    stop("fitness and offspring_fitness must have the same number of objectives");
  if (num_to_replace < 0)
    stop("num_to_replace must be >= 0");

  const int population_size = population.nrow();
  const int num_genes = population.ncol();
  const int num_offspring = offspring.nrow();

  if (num_to_replace == 0) {
    return clone(population);
  }

  int cap = num_to_replace;
  if (cap > num_offspring) cap = num_offspring;
  if (cap > population_size) cap = population_size;

  const int keep_parents = population_size - cap;

  std::vector<int> parent_ids;
  if (keep_parents > 0) {
    parent_ids = bioga::environmental_selection(fitness, keep_parents);
  }
  const std::vector<int> offspring_ids = bioga::environmental_selection(offspring_fitness, cap);

  NumericMatrix updated(population_size, num_genes);
  int row = 0;
  for (std::size_t a = 0; a < parent_ids.size(); ++a) {
    for (int j = 0; j < num_genes; ++j) updated(row, j) = population(parent_ids[a], j);
    ++row;
  }
  for (std::size_t a = 0; a < offspring_ids.size(); ++a) {
    for (int j = 0; j < num_genes; ++j) updated(row, j) = offspring(offspring_ids[a], j);
    ++row;
  }

  return updated;
}
