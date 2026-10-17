#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include "BioGA.h"
using namespace Rcpp;

//' Mutate offspring with an adaptive rate
//'
//' The mutation rate increases linearly with the generation index. Real-coded
//' mode applies Gaussian perturbation; binary mode flips selection bits.
//' When a gene network is supplied, the perturbation of each gene is damped by
//' its relative connectivity (degree), so network hubs change less.
//'
//' @param offspring Numeric matrix of offspring.
//' @param mutation_rate Base probability of mutation.
//' @param iteration Current GA iteration.
//' @param max_iterations Maximum number of GA iterations.
//' @param network Optional symmetric adjacency matrix
//' (rows: genes, columns: genes).
//' @param representation Either `"real"` or `"binary"`.
//' @return Numeric matrix of mutated offspring.
//' @family mutation
//' @seealso [initialize_population_cpp()], [evaluate_fitness_cpp()],
//'   [selection_cpp()], [crossover_cpp()], [replacement_cpp()],
//'   [bioga_main_cpp()]
//' @examples
//' genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
//' population <- BioGA::initialize_population_cpp(genomic_data, population_size = 5)
//' fitness <- BioGA::evaluate_fitness_cpp(genomic_data, population, c(1.0, 0.5))
//' parents <- BioGA::selection_cpp(population, fitness, num_parents = 2)
//' offspring <- BioGA::crossover_cpp(parents, offspring_size = 2)
//' BioGA::mutation_cpp(offspring, mutation_rate = 0.1, iteration = 1,
//'   max_iterations = 100)
//' @export
// [[Rcpp::export]]
NumericMatrix mutation_cpp(const NumericMatrix& offspring,
                           double mutation_rate,
                           int iteration,
                           int max_iterations,
                           Nullable<NumericMatrix> network = R_NilValue,
                           std::string representation = "real") {
  if (offspring.nrow() == 0)
    stop("offspring must have at least 1 individual");
  if (mutation_rate < 0.0 || mutation_rate > 1.0)
    stop("mutation_rate must be in [0, 1]");
  if (max_iterations < 1)
    stop("max_iterations must be >= 1");
  if (iteration < 0)
    stop("iteration must be >= 0");
  if (representation != "real" && representation != "binary")
    stop("representation must be either 'real' or 'binary'");

  const bool binary = (representation == "binary");
  const int num_offspring = offspring.nrow();
  const int num_genes = offspring.ncol();

  const double adaptive_rate =
    std::min(1.0, mutation_rate * (1.0 + 0.5 * static_cast<double>(iteration) / max_iterations));

  NumericMatrix mutated = clone(offspring);

  // Degree-based damping weights (1 = never mutate, 0 = full mutation).
  std::vector<double> damping(num_genes, 0.0);
  if (network.isNotNull()) {
    NumericMatrix net(network);
    if (net.nrow() != num_genes || net.ncol() != num_genes)
      stop("network must be a square matrix with one row/column per gene");
    double max_degree = 0.0;
    std::vector<double> degree(num_genes, 0.0);
    for (int j = 0; j < num_genes; ++j) {
      for (int k = 0; k < num_genes; ++k) {
        if (j != k) degree[j] += std::fabs(net(j, k));
      }
      if (degree[j] > max_degree) max_degree = degree[j];
    }
    if (max_degree > 0.0) {
      for (int j = 0; j < num_genes; ++j) damping[j] = degree[j] / max_degree;
    }
  }

  for (int i = 0; i < num_offspring; ++i) {
    for (int j = 0; j < num_genes; ++j) {
      if (R::runif(0.0, 1.0) >= adaptive_rate) continue;
      if (binary) {
        mutated(i, j) = (mutated(i, j) > 0.5) ? 0.0 : 1.0;
      } else {
        const double perturbation = R::rnorm(0.0, 0.1);
        mutated(i, j) += perturbation * (1.0 - damping[j]);
      }
    }
  }

  return mutated;
}
