#include <Rcpp.h>
using namespace Rcpp;

//' Function to initialize population with optional gene clustering
//'
//' @param genomic_data Numeric matrix of genomic data (rows: genes, columns: samples).
//' @param population_size Number of individuals in the population.
//' @param seed Optional random seed for reproducibility.
//' @param clusters Optional vector of gene cluster assignments.
//' @return Numeric matrix of initialized population.
//' @examples
//' genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
//' BioGA::initialize_population_cpp(genomic_data, population_size = 5, seed = 123)
//' @export
// [[Rcpp::export]]
NumericMatrix initialize_population_cpp(const NumericMatrix& genomic_data, 
                                      int population_size, 
                                      Nullable<int> seed = R_NilValue,
                                      Nullable<IntegerVector> clusters = R_NilValue) {
 int num_genes = genomic_data.nrow();
 int num_samples = genomic_data.ncol();
 
 // Set random seed if provided
 if (seed.isNotNull()) {
    int s = as<int>(seed);
    Environment base_env("package:base");
    Function set_seed = base_env["set.seed"];
    set_seed(s);
 }
 
 // Allocate memory for population matrix
 NumericMatrix population(population_size, num_genes);
 
 // Initialize with clustering if provided
 if (clusters.isNotNull()) {
    IntegerVector cluster_assignments(clusters);
    std::vector<std::vector<int>> cluster_indices(*std::max_element(cluster_assignments.begin(), 
                                                                    cluster_assignments.end()) + 1);
    for (int j = 0; j < num_genes; ++j) {
       cluster_indices[cluster_assignments[j]].push_back(j);
    }
    for (int i = 0; i < population_size; ++i) {
       for (size_t c = 0; c < cluster_indices.size(); ++c) {
          if (!cluster_indices[c].empty()) {
             int sample_idx = int(R::unif_rand() * num_samples);
             for (int j : cluster_indices[c]) {
                population(i, j) = genomic_data(j, sample_idx);
             }
          }
       }
    }
 } else {
    // Default random initialization
    for (int i = 0; i < population_size; ++i) {
       for (int j = 0; j < num_genes; ++j) {
          int sample_index = int(R::unif_rand() * num_samples);
          population(i, j) = genomic_data(j, sample_index);
       }
    }
 }
 
 return population;
}