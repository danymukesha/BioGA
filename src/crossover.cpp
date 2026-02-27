#include <RcppParallel.h>
#include <Rcpp.h>
#include <vector>
#include <algorithm>

using namespace Rcpp;
using namespace RcppParallel;

//' Function to perform simulated binary crossover (SBX)
//'
//' @param selected_parents Numeric matrix of selected individuals.
//' @param offspring_size Number of offspring to generate.
//' @param crossover_rate Probability of crossover (default: 0.9).
//' @param eta_c SBX distribution index (default: 20.0).
//' @return Numeric matrix of offspring.
//' @examples
//' genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
//' population <- BioGA::initialize_population_cpp(genomic_data, 
//' population_size = 5)
//' fitness <- BioGA::evaluate_fitness_cpp(genomic_data, population, 
//' c(1.0, 0.5))
//' selected_parents <- BioGA::selection_cpp(population, fitness, 
//' num_parents = 2)
//' BioGA::crossover_cpp(selected_parents, offspring_size = 2)
//' @export
// [[Rcpp::export]]
NumericMatrix crossover_cpp(const NumericMatrix& selected_parents, 
                            int offspring_size, 
                            double crossover_rate = 0.9, 
                            double eta_c = 20.0) {
   int num_parents = selected_parents.nrow();
   int num_genes = selected_parents.ncol();
   NumericMatrix offspring(offspring_size, num_genes);
   
   double diversity = 0.0;
   for (int j = 0; j < num_genes; ++j) {
      double mean = 0.0, variance = 0.0;
      for (int i = 0; i < num_parents; ++i) mean += selected_parents(i, j);
      mean /= num_parents;
      for (int i = 0; i < num_parents; ++i) variance += pow(selected_parents(i, j) - mean, 2);
      diversity += sqrt(variance / num_parents);
   }
   diversity /= num_genes;
   double adaptive_crossover_rate = crossover_rate * (1.0 - std::min(1.0, diversity));
   
   for (int i = 0; i < offspring_size; ++i) {
      int p1_idx = int(R::unif_rand() * num_parents);
      int p2_idx = int(R::unif_rand() * num_parents);
      
      for (int j = 0; j < num_genes; ++j) {
         if (R::runif(0.0, 1.0) < adaptive_crossover_rate) {
            double u = R::runif(0.0, 1.0);
            double beta = (u <= 0.5) ? pow(2.0 * u, 1.0 / (eta_c + 1.0)) : 
               pow(1.0 / (2.0 * (1.0 - u)), 1.0 / (eta_c + 1.0));
            double p1 = selected_parents(p1_idx, j);
            double p2 = selected_parents(p2_idx, j);
            offspring(i, j) = 0.5 * ((p1 + p2) - beta * fabs(p2 - p1));
         } else {
            offspring(i, j) = selected_parents(p1_idx, j);
         }
      }
   }
   return offspring;
}
