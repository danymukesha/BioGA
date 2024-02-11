#include <Rcpp.h>
using namespace Rcpp;

//' Function to perform crossover between selected individuals
//'
//' @param selected_parents Numeric matrix representing the selected individuals.
//' @param offspring_size Number of offspring to generate.
//' @return Numeric matrix representing the offspring.
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
