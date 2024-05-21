#include <Rcpp.h>
using namespace Rcpp;

//' Function to select individuals based on fitness scores
//'
//' @param population Numeric matrix representing the population of
//' individuals.
//' @param fitness Numeric vector of fitness scores for each individual.
//' @param num_parents Number of individuals to select.
//' @return Numeric matrix representing the selected individuals.
//' @examples
//' # example of usage
//' genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
//' population <- BioGA::initialize_population_cpp(genomic_data,
//'                 population_size = 5)
//' fitness <- BioGA::evaluate_fitness_cpp(genomic_data, population)
//' BioGA::selection_cpp(population, fitness, num_parents = 2)
//' @export
// [[Rcpp::export]]
NumericMatrix selection_cpp(const NumericMatrix& population, 
                            const NumericVector& fitness, 
                            int num_parents) {
  int population_size = population.nrow();
  int num_genes = population.ncol();

  // Allocate memory for the selected parents matrix
  NumericMatrix selected_parents(num_parents, num_genes);

  // Allocate memory for the selected flag array
  std::vector<bool> selected(population_size, false);

  // Select parents based on fitness scores
  for (int i = 0; i < num_parents; ++i) {
    // Find the individual with the minimum fitness score 
    // that has not been selected
    int min_fitness_index = -1;
    double min_fitness = std::numeric_limits<double>::infinity();
    for (int j = 0; j < population_size; ++j) {
      if (!selected[j] && fitness[j] < min_fitness) {
        min_fitness = fitness[j];
        min_fitness_index = j;
      }
    }

    // Copy the selected individual to the selected parents matrix
    for (int j = 0; j < num_genes; ++j) {
      selected_parents(i, j) = population(min_fitness_index, j);
    }

    // Mark the selected individual as selected
    selected[min_fitness_index] = true;
  }

  return selected_parents;
}
