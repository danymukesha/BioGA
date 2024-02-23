#include <Rcpp.h>
using namespace Rcpp;

//' Function to perform crossover between selected individuals
//'
//' @param selected_parents Numeric matrix representing the selected
//' individuals.
//' @param offspring_size Number of offspring to generate.
//' @return Numeric matrix representing the offspring.
//' @examples
//' # example of usage
//' genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
//' population <- BioGA::initialize_population_cpp(genomic_data,
//'                 population_size = 5)
//' fitness <- BioGA::evaluate_fitness_cpp(genomic_data, population)
//' selected_parents <- BioGA::selection_cpp(population, fitness,
//'                 num_parents = 2)
//' BioGA::crossover_cpp(selected_parents, offspring_size = 2)
//' @export
// [[Rcpp::export]]
NumericMatrix crossover_cpp(const NumericMatrix& selected_parents, int offspring_size) {
  int num_parents = selected_parents.nrow();
  int num_genes = selected_parents.ncol();

  // Allocate memory for the offspring matrix
  NumericMatrix offspring(offspring_size, num_genes);

  // Perform crossover between selected parents
  for (int i = 0; i < offspring_size; ++i) {
    int parent1_index = rand() % num_parents;
    int parent2_index = rand() % num_parents;
    for (int j = 0; j < num_genes; ++j) {
      offspring(i, j) = (selected_parents(parent1_index, j) + selected_parents(parent2_index, j)) / 2.0;
    }
  }

  return offspring;
}
