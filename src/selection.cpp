#include <Rcpp.h>
#include <vector>
#include "BioGA.h"
#include "nsga2.h"
using namespace Rcpp;

//' Select parents using NSGA-II rank and crowding distance
//'
//' Ranks the population into non-dominated fronts, computes crowding
//' distances, and performs binary tournament selection that prefers lower
//' rank and, on ties, higher crowding distance.
//'
//' @param population Numeric matrix of individuals.
//' @param fitness Numeric matrix of fitness scores (rows: individuals,
//' columns: objectives, all minimised).
//' @param num_parents Number of individuals to select.
//' @return Numeric matrix of selected individuals.
//' @family selection
//' @seealso [initialize_population_cpp()], [evaluate_fitness_cpp()],
//'   [crossover_cpp()], [mutation_cpp()], [replacement_cpp()],
//'   [bioga_main_cpp()]
//' @references
//'   Deb K, Pratap A, Agarwal S, Meyarivan T (2002).
//'   "A fast and elitist multiobjective genetic algorithm: NSGA-II."
//'   _IEEE Transactions on Evolutionary Computation_, 6(2), 182-197.
//'   \doi{10.1109/4235.996017}
//' @examples
//' genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
//' population <- BioGA::initialize_population_cpp(genomic_data, population_size = 5)
//' fitness <- BioGA::evaluate_fitness_cpp(genomic_data, population, c(1.0, 0.5))
//' BioGA::selection_cpp(population, fitness, num_parents = 2)
//' @export
// [[Rcpp::export]]
NumericMatrix selection_cpp(const NumericMatrix& population,
                            const NumericMatrix& fitness,
                            int num_parents) {
  if (population.nrow() == 0)
    stop("population must have at least 1 individual");
  if (fitness.nrow() != population.nrow())
    stop("fitness rows must match population rows");
  if (fitness.ncol() < 1)
    stop("fitness must have at least one objective");
  if (num_parents < 1)
    stop("num_parents must be >= 1");
  if (num_parents > population.nrow())
    stop("num_parents cannot exceed population size");

  const int population_size = population.nrow();
  const int num_genes = population.ncol();

  std::vector<int> rank;
  std::vector<double> crowding;
  bioga::compute_rank_and_crowding(fitness, rank, crowding);

  NumericMatrix selected_parents(num_parents, num_genes);
  for (int i = 0; i < num_parents; ++i) {
    const int idx1 = static_cast<int>(R::unif_rand() * population_size);
    const int idx2 = static_cast<int>(R::unif_rand() * population_size);

    int winner = idx1;
    if (rank[idx2] < rank[idx1]) {
      winner = idx2;
    } else if (rank[idx2] == rank[idx1] && crowding[idx2] > crowding[idx1]) {
      winner = idx2;
    }

    for (int j = 0; j < num_genes; ++j) {
      selected_parents(i, j) = population(winner, j);
    }
  }

  return selected_parents;
}
