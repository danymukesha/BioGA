#include <Rcpp.h>
using namespace Rcpp;

//' Function to replace population with elitism and diversity preservation
//'
//' @param population Numeric matrix of individuals.
//' @param offspring Numeric matrix of offspring.
//' @param fitness Numeric matrix of population fitness scores.
//' @param offspring_fitness Numeric matrix of offspring fitness scores.
//' @param num_to_replace Number of individuals to replace.
//' @return Numeric matrix of updated population.
//' @examples
//' genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
//' population <- BioGA::initialize_population_cpp(genomic_data, 
//' population_size = 5)
//' fitness <- BioGA::evaluate_fitness_cpp(genomic_data, population, 
//' c(1.0, 0.5))
//' selected_parents <- BioGA::selection_cpp(population, fitness, 
//' num_parents = 2)
//' offspring <- BioGA::crossover_cpp(selected_parents, offspring_size = 2)
//' offspring_fitness <- BioGA::evaluate_fitness_cpp(genomic_data, offspring, 
//' c(1.0, 0.5))
//' BioGA::replacement_cpp(population, offspring, fitness, offspring_fitness, 
//' num_to_replace = 1)
//' @export
// [[Rcpp::export]]
NumericMatrix replacement_cpp(const NumericMatrix& population, 
                             const NumericMatrix& offspring, 
                             const NumericMatrix& fitness, 
                             const NumericMatrix& offspring_fitness, 
                             int num_to_replace) {
 int population_size = population.nrow();
 int num_genes = population.ncol();
 int num_offspring = offspring.nrow();
 
 // Allocate memory for updated population
 NumericMatrix updated_population = clone(population);
 
 // Elitism: Copy best individual (based on first objective)
 int best_idx = 0;
 double best_fitness = fitness(0, 0);
 for (int i = 1; i < population_size; ++i) {
   if (fitness(i, 0) < best_fitness) {
     best_fitness = fitness(i, 0);
     best_idx = i;
   }
 }
 for (int j = 0; j < num_genes; ++j) {
   updated_population(0, j) = population(best_idx, j);
 }
 
 // Replace remaining individuals with diversity consideration
 for (int i = 0; i < num_to_replace && i < num_offspring; ++i) {
   int index_to_replace = 1 + int(R::unif_rand() * (population_size - 1)); // Preserve elite
   double diversity = 0.0;
   for (int j = 0; j < num_genes; ++j) {
     double diff = offspring(i, j) - updated_population(index_to_replace, j);
     diversity += diff * diff;
   }
   if (diversity > 1e-6) { // Replace only if sufficiently diverse
     for (int j = 0; j < num_genes; ++j) {
       updated_population(index_to_replace, j) = offspring(i, j);
     }
   }
 }
 
 return updated_population;
}