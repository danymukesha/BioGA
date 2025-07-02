#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;

//' Function to select individuals using NSGA-II non-dominated sorting
//'
//' @param population Numeric matrix of individuals.
//' @param fitness Numeric matrix of fitness scores (columns: objectives).
//' @param num_parents Number of individuals to select.
//' @return Numeric matrix of selected individuals.
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
 int population_size = population.nrow();
 int num_genes = population.ncol();
 int num_objectives = fitness.ncol();
 
 // Allocate memory for selected parents
 NumericMatrix selected_parents(num_parents, num_genes);
 
 // Non-dominated sorting (simplified NSGA-II)
 std::vector<int> front(population_size, 0);
 std::vector<int> dominated_count(population_size, 0);
 std::vector<std::vector<int>> dominates(population_size);
 
 // Compute domination counts
 for (int i = 0; i < population_size; ++i) {
   for (int j = 0; j < population_size; ++j) {
     if (i == j) continue;
     bool dominates_i = true, dominates_j = true;
     for (int k = 0; k < num_objectives; ++k) {
       if (fitness(i, k) > fitness(j, k)) dominates_i = false;
       if (fitness(j, k) > fitness(i, k)) dominates_j = false;
     }
     if (dominates_i) {
       dominates[i].push_back(j);
       dominated_count[j]++;
     }
   }
 }
 
 // Assign/Build fronts
 std::vector<int> current_front;
 for (int i = 0; i < population_size; ++i) {
   if (dominated_count[i] == 0) {
     front[i] = 1;
     current_front.push_back(i);
   }
 }
 Rcpp::Rcout << "Current front size: " << current_front.size() << std::endl;
 
 if (current_front.empty()) {
    Rcpp::Rcout << "Warning: No non-dominated individuals found. \
    Using full population for selection.\n";
    // fallback: select from full population randomly
    for (int i = 0; i < population_size; ++i) {
       current_front.push_back(i);
    }
 }
 
 // if (current_front.empty()) {
 //    stop("No non-dominated individuals found in the current front.");
 // }
 
 // Tournament selection from first front
 for (int i = 0; i < num_parents; ++i) {
   int idx1 = current_front[int(R::unif_rand() * current_front.size())];
   int idx2 = current_front[int(R::unif_rand() * current_front.size())];
   int selected_idx = (fitness(idx1, 0) < fitness(idx2, 0)) ? idx1 : idx2; // Compare first objective
   for (int j = 0; j < num_genes; ++j) {
     selected_parents(i, j) = population(selected_idx, j);
   }
 }
 return selected_parents;
}