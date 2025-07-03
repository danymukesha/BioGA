#include <Rcpp.h>
using namespace Rcpp;

//' Function to mutate offspring with adaptive mutation and network constraints
//'
//' @param offspring Numeric matrix of offspring.
//' @param mutation_rate Base probability of mutation.
//' @param iteration Current GA iteration for adaptive mutation.
//' @param max_iterations Maximum number of GA iterations.
//' @param network Optional matrix of gene network constraints (rows: genes, 
//' cols: genes).
//' @return Numeric matrix of mutated offspring.
//' @examples
//' genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
//' population <- BioGA::initialize_population_cpp(genomic_data, 
//' population_size = 5)
//' fitness <- BioGA::evaluate_fitness_cpp(genomic_data, population, 
//' c(1.0, 0.5))
//' selected_parents <- BioGA::selection_cpp(population, fitness, 
//' num_parents = 2)
//' offspring <- BioGA::crossover_cpp(selected_parents, 
//' offspring_size = 2)
//' BioGA::mutation_cpp(offspring, mutation_rate = 0.1, iteration = 1, 
//' max_iterations = 100)
//' @export
// [[Rcpp::export]]
NumericMatrix mutation_cpp(const NumericMatrix& offspring, 
                          double mutation_rate, 
                          int iteration, 
                          int max_iterations, 
                          Nullable<NumericMatrix> network = R_NilValue) {
 int num_offspring = offspring.nrow();
 int num_genes = offspring.ncol();
 
 // Adaptive mutation rate (increases as iterations progress)
 double adaptive_mutation_rate = mutation_rate * (1.0 + 0.5 * (double)iteration / max_iterations);
 
 // Allocate memory for mutated offspring
 NumericMatrix mutated_offspring = clone(offspring);
 
 // Apply mutation with network constraints
 for (int i = 0; i < num_offspring; ++i) {
   for (int j = 0; j < num_genes; ++j) {
     if (R::runif(0.0, 1.0) < adaptive_mutation_rate) {
       double mutation = R::rnorm(0.0, 0.1);
       // Apply network constraint (if provided)
       if (network.isNotNull()) {
         NumericMatrix net(network);
         double constraint_factor = 0.0;
         for (int k = 0; k < num_genes; ++k) {
           constraint_factor += net(j, k) * mutated_offspring(i, k);
         }
         mutated_offspring(i, j) += mutation * (1.0 - constraint_factor);
       } else {
         mutated_offspring(i, j) += mutation;
       }
     }
   }
 }
 
 return mutated_offspring;
}