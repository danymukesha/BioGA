#include <Rcpp.h>
using namespace Rcpp;

//' Function to replace non-selected individuals in the population
//'
//' Replace non-selected individuals in the population
//' @param population Numeric matrix representing the population of individuals.
//' @param offspring Numeric matrix representing the offspring.
//' @param num_to_replace Number of individuals to replace.
//' @return Numeric matrix representing the updated population.
//' @export
// [[Rcpp::export]]
NumericMatrix replacement_cpp(const NumericMatrix& population, const NumericMatrix& offspring, int num_to_replace) {
  int population_size = population.nrow();
  int num_genes = population.ncol();

  // Allocate memory for the updated population matrix
  NumericMatrix updated_population = clone(population);

  // Replace non-selected individuals in the population with offspring
  for (int i = 0; i < num_to_replace; ++i) {
    int index_to_replace = rand() % population_size;
    for (int j = 0; j < num_genes; ++j) {
      updated_population(index_to_replace, j) = offspring(i, j);
    }
  }

  return updated_population;
}
