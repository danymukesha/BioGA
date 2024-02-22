#include <Rcpp.h>
using namespace Rcpp;

//' Function to evaluate fitness using genomic data
//'
//' @param genomic_data Numeric matrix of genomic data where rows represent genes/features and columns represent samples.
//' @param population Numeric matrix representing the population of individuals.
//' @return Numeric vector of fitness scores for each individual.
//' @examples
//' # example of usage
//' genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
//' population <- BioGA::initialize_population_cpp(genomic_data, population_size = 5)
//' BioGA::evaluate_fitness_cpp(genomic_data, population)
//' @export
// [[Rcpp::export]]
NumericVector evaluate_fitness_cpp(const NumericMatrix& genomic_data, const NumericMatrix& population) {
  int num_genes = genomic_data.nrow();
  int num_samples = genomic_data.ncol();
  int population_size = population.nrow();

  // Allocate memory for the fitness vector
  NumericVector fitness(population_size);

  // Calculate fitness for each individual in the population
  for (int i = 0; i < population_size; ++i) {
    double sum_fitness = 0.0;
    for (int j = 0; j < num_genes; ++j) {
      double sum_gene_expression = 0.0;
      for (int k = 0; k < num_samples; ++k) {
        sum_gene_expression += pow(genomic_data(j, k) - population(i, j), 2);
      }
      sum_fitness += sum_gene_expression;
    }
    fitness[i] = sum_fitness;
  }

  return fitness;
}
