#include <Rcpp.h>
using namespace Rcpp;

//' Function to replace non-selected individuals in the population
//'
//' Replace non-selected individuals in the population
//' @param population Numeric matrix representing the population of
//' individuals.
//' @param offspring Numeric matrix representing the offspring.
//' @param num_to_replace Number of individuals to replace.
//' @return Numeric matrix representing the updated population.
//' @examples
//' # example of usage
//' genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
//' population <- BioGA::initialize_population_cpp(genomic_data,
//'                 population_size = 5)
//' fitness <- BioGA::evaluate_fitness_cpp(genomic_data, population)
//' selected_parents <- BioGA::selection_cpp(population, fitness,
//'                       num_parents = 2)
//' offspring <- BioGA::crossover_cpp(selected_parents, offspring_size = 2)
//' mutated_offspring <- BioGA::mutation_cpp(offspring, mutation_rate = 0)
//' BioGA::replacement_cpp(population, mutated_offspring, num_to_replace = 1)
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
