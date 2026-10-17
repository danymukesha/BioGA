#include <Rcpp.h>
#include <cmath>
#include "BioGA.h"
using namespace Rcpp;

//' Generate offspring by crossover
//'
//' Real-coded mode uses simulated binary crossover (SBX); binary mode uses
//' uniform crossover. When only one parent is available the offspring are
//' exact copies of that parent.
//'
//' @param selected_parents Numeric matrix of selected individuals.
//' @param offspring_size Number of offspring to generate.
//' @param crossover_rate Probability of crossover per offspring (default: 0.9).
//' @param eta_c SBX distribution index (default: 20.0).
//' @param representation Either `"real"` or `"binary"`.
//' @return Numeric matrix of offspring.
//' @family crossover
//' @seealso [initialize_population_cpp()], [evaluate_fitness_cpp()],
//'   [selection_cpp()], [mutation_cpp()], [replacement_cpp()],
//'   [bioga_main_cpp()]
//' @references
//'   Deb K, Agrawal RB (1995).
//'   "Simulated binary crossover for continuous search space."
//'   _Complex Systems_, 9(2), 115-148.
//' @examples
//' genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
//' population <- BioGA::initialize_population_cpp(genomic_data, population_size = 5)
//' fitness <- BioGA::evaluate_fitness_cpp(genomic_data, population, c(1.0, 0.5))
//' parents <- BioGA::selection_cpp(population, fitness, num_parents = 2)
//' BioGA::crossover_cpp(parents, offspring_size = 2)
//' @export
// [[Rcpp::export]]
NumericMatrix crossover_cpp(const NumericMatrix& selected_parents,
                            int offspring_size,
                            double crossover_rate = 0.9,
                            double eta_c = 20.0,
                            std::string representation = "real") {
  if (selected_parents.nrow() == 0)
    stop("selected_parents must have at least 1 parent");
  if (offspring_size < 1)
    stop("offspring_size must be >= 1");
  if (crossover_rate < 0.0 || crossover_rate > 1.0)
    stop("crossover_rate must be in [0, 1]");
  if (eta_c <= 0.0)
    stop("eta_c must be > 0");
  if (representation != "real" && representation != "binary")
    stop("representation must be either 'real' or 'binary'");

  const bool binary = (representation == "binary");
  const int num_parents = selected_parents.nrow();
  const int num_genes = selected_parents.ncol();

  NumericMatrix offspring(offspring_size, num_genes);

  // With a single parent there is nothing to recombine.
  if (num_parents < 2) {
    for (int i = 0; i < offspring_size; ++i) {
      for (int j = 0; j < num_genes; ++j) {
        offspring(i, j) = selected_parents(0, j);
      }
    }
    return offspring;
  }

  for (int i = 0; i < offspring_size; ++i) {
    int parent1 = static_cast<int>(R::unif_rand() * num_parents);
    int parent2 = static_cast<int>(R::unif_rand() * num_parents);
    while (parent2 == parent1) {
      parent2 = static_cast<int>(R::unif_rand() * num_parents);
    }

    if (R::runif(0.0, 1.0) >= crossover_rate) {
      for (int j = 0; j < num_genes; ++j) {
        offspring(i, j) = selected_parents(parent1, j);
      }
      continue;
    }

    if (binary) {
      for (int j = 0; j < num_genes; ++j) {
        offspring(i, j) = (R::unif_rand() < 0.5) ? selected_parents(parent1, j)
                                                 : selected_parents(parent2, j);
      }
    } else {
      for (int j = 0; j < num_genes; ++j) {
        const double x = selected_parents(parent1, j);
        const double y = selected_parents(parent2, j);
        const double u = R::runif(0.0, 1.0);
        const double beta = (u <= 0.5)
          ? std::pow(2.0 * u, 1.0 / (eta_c + 1.0))
          : std::pow(1.0 / (2.0 * (1.0 - u)), 1.0 / (eta_c + 1.0));
        const double child1 = 0.5 * ((1.0 + beta) * x + (1.0 - beta) * y);
        const double child2 = 0.5 * ((1.0 - beta) * x + (1.0 + beta) * y);
        offspring(i, j) = (R::unif_rand() < 0.5) ? child1 : child2;
      }
    }
  }

  return offspring;
}
