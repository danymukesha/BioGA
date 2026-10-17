#include <Rcpp.h>
#include "BioGA.h"
using namespace Rcpp;

//' Initialize a population, optionally guided by gene clusters
//'
//' Creates the starting population. In `"real"` mode each gene value is drawn
//' from an observed sample of that gene, so individuals preserve the marginal
//' distribution of the data. In `"binary"` mode each gene is encoded as
//' selected (1) or not selected (0).
//'
//' @param genomic_data Numeric matrix of genomic data (rows: genes,
//' columns: samples).
//' @param population_size Number of individuals in the population.
//' @param seed Optional random seed for reproducibility.
//' @param clusters Optional vector of gene cluster assignments (length equal
//' to the number of genes, non-negative integers).
//' @param representation Either `"real"` (continuous expression values) or
//' `"binary"` (gene selection indicators).
//' @return Numeric matrix of initialized population (rows: individuals,
//' columns: genes).
//' @family population initialization
//' @seealso [evaluate_fitness_cpp()], [selection_cpp()], [crossover_cpp()],
//'   [mutation_cpp()], [replacement_cpp()], [bioga_main_cpp()]
//' @examples
//' genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
//' BioGA::initialize_population_cpp(genomic_data, population_size = 5, seed = 123)
//' @export
// [[Rcpp::export]]
NumericMatrix initialize_population_cpp(const NumericMatrix& genomic_data,
                                        int population_size,
                                        Nullable<int> seed = R_NilValue,
                                        Nullable<IntegerVector> clusters = R_NilValue,
                                        std::string representation = "real") {
  if (genomic_data.nrow() == 0 || genomic_data.ncol() == 0)
    stop("genomic_data must have at least 1 row and 1 column");
  if (population_size < 1)
    stop("population_size must be >= 1");
  if (representation != "real" && representation != "binary")
    stop("representation must be either 'real' or 'binary'");

  const bool binary = (representation == "binary");
  const int num_genes = genomic_data.nrow();
  const int num_samples = genomic_data.ncol();

  // Optional reproducible seed.
  if (seed.isNotNull()) {
    Environment base = Environment::base_env();
    Function set_seed = base["set.seed"];
    set_seed(as<int>(seed));
  }

  NumericMatrix population(population_size, num_genes);

  if (clusters.isNotNull()) {
    IntegerVector cluster_assignments(clusters);
    if (cluster_assignments.size() != num_genes)
      stop("clusters must have one assignment per gene (length == nrow(genomic_data))");

    int max_cluster = -1;
    for (int j = 0; j < num_genes; ++j) {
      if (cluster_assignments[j] < 0)
        stop("clusters must contain non-negative integers");
      if (cluster_assignments[j] > max_cluster) max_cluster = cluster_assignments[j];
    }

    std::vector<std::vector<int> > cluster_indices(max_cluster + 1);
    for (int j = 0; j < num_genes; ++j) {
      cluster_indices[cluster_assignments[j]].push_back(j);
    }

    for (int i = 0; i < population_size; ++i) {
      for (std::size_t c = 0; c < cluster_indices.size(); ++c) {
        if (cluster_indices[c].empty()) continue;
        if (binary) {
          const double value = (R::unif_rand() < 0.5) ? 1.0 : 0.0;
          for (std::size_t a = 0; a < cluster_indices[c].size(); ++a) {
            population(i, cluster_indices[c][a]) = value;
          }
        } else {
          const int sample_idx = static_cast<int>(R::unif_rand() * num_samples);
          for (std::size_t a = 0; a < cluster_indices[c].size(); ++a) {
            population(i, cluster_indices[c][a]) = genomic_data(cluster_indices[c][a], sample_idx);
          }
        }
      }
    }
  } else {
    for (int i = 0; i < population_size; ++i) {
      for (int j = 0; j < num_genes; ++j) {
        if (binary) {
          population(i, j) = (R::unif_rand() < 0.5) ? 1.0 : 0.0;
        } else {
          const int sample_index = static_cast<int>(R::unif_rand() * num_samples);
          population(i, j) = genomic_data(j, sample_index);
        }
      }
    }
  }

  return population;
}
