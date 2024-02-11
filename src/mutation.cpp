#include <Rcpp.h>
using namespace Rcpp;

//' Function to mutate the offspring
//'
//' @param offspring Numeric matrix representing the offspring.
//' @param mutation_rate Probability of mutation for each individual.
//' @return Numeric matrix representing the mutated offspring.
//' @export
// [[Rcpp::export]]
NumericMatrix mutation_cpp(const NumericMatrix& offspring, double mutation_rate) {
  int num_offspring = offspring.nrow();
  int num_genes = offspring.ncol();

  // Allocate memory for the mutated offspring matrix
  NumericMatrix mutated_offspring = clone(offspring);

  // Mutate offspring with a certain probability
  for (int i = 0; i < num_offspring; ++i) {
    for (int j = 0; j < num_genes; ++j) {
      if (R::runif(0.0, 1.0) < mutation_rate) {
        mutated_offspring(i, j) += R::rnorm(0.0, 0.1); // Add small random value
      }
    }
  }

  return mutated_offspring;
}
