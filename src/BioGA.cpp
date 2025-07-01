#include <Rcpp.h>
#include "BioGA.h"
using namespace Rcpp;

//' Main genetic algorithm loop for genomic data optimization
//'
//' @param genomic_data Numeric matrix of genomic data (rows: genes, columns: samples).
//' @param population_size Number of individuals in the population.
//' @param num_generations Number of generations to run.
//' @param crossover_rate Probability of crossover.
//' @param mutation_rate Base probability of mutation.
//' @param num_parents Number of parents to select per generation.
//' @param num_offspring Number of offspring to generate per generation.
//' @param num_to_replace Number of individuals to replace per generation.
//' @param weights Numeric vector of weights for multi-objective fitness.
//' @param seed Optional random seed for reproducibility.
//' @param clusters Optional vector of gene cluster assignments.
//' @param network Optional matrix of gene network constraints.
//' @return List containing final population and fitness scores.
//' @examples
//' genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
//' result <- BioGA::bioga_main_cpp(genomic_data, population_size = 50, num_generations = 100,
//'                                 crossover_rate = 0.9, mutation_rate = 0.1,
//'                                 num_parents = 20, num_offspring = 20, num_to_replace = 10,
//'                                 weights = c(1.0, 0.5), seed = 123)
//' @export
// [[Rcpp::export]]
List bioga_main_cpp(const NumericMatrix& genomic_data,
                 int population_size,
                 int num_generations,
                 double crossover_rate,
                 double eta_c,
                 double mutation_rate,
                 int num_parents,
                 int num_offspring,
                 int num_to_replace,
                 const NumericVector& weights,
                 Nullable<int> seed = R_NilValue,
                 Nullable<IntegerVector> clusters = R_NilValue,
                 Nullable<NumericMatrix> network = R_NilValue) {

 NumericMatrix population = initialize_population_cpp(genomic_data, population_size, seed, clusters);
 
 // Main GA loop
 for (int generation = 0; generation < num_generations; ++generation) {
     
     NumericMatrix fitness = evaluate_fitness_cpp(genomic_data, population, weights);
     
     NumericMatrix parents = selection_cpp(population, fitness, num_parents);
     
     NumericMatrix offspring = crossover_cpp(parents, num_offspring, crossover_rate, eta_c);
     
     NumericMatrix mutated_offspring = mutation_cpp(offspring, mutation_rate, generation, num_generations, network);
     
     NumericMatrix offspring_fitness = evaluate_fitness_cpp(genomic_data, mutated_offspring, weights);
     
     population = replacement_cpp(population, mutated_offspring, fitness, offspring_fitness, num_to_replace);
 }
 
 NumericMatrix final_fitness = evaluate_fitness_cpp(genomic_data, population, weights);
 
 return List::create(Named("population") = population, Named("fitness") = final_fitness);
}