#include <RcppParallel.h>
#include <Rcpp.h>
#include <vector>
#include <algorithm>

using namespace Rcpp;
using namespace RcppParallel;

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
   int pop_size = population.nrow();
   int off_size = offspring.nrow();
   int num_genes = population.ncol();
   int num_objs = fitness.ncol();
   
   // 1. Combine parents and offspring into a pool
   int combined_size = pop_size + off_size;
   NumericMatrix pool_pop(combined_size, num_genes);
   NumericMatrix pool_fit(combined_size, num_objs);
   
   for(int i = 0; i < pop_size; ++i) {
      pool_pop(i, _) = population(i, _);
      pool_fit(i, _) = fitness(i, _);
   }
   for(int i = 0; i < off_size; ++i) {
      pool_pop(i + pop_size, _) = offspring(i, _);
      pool_fit(i + pop_size, _) = offspring_fitness(i, _);
   }
   
   // 2. Find Non-Dominated Elites (Front 1) of the combined pool
   std::vector<int> front1;
   for (int i = 0; i < combined_size; ++i) {
      bool dominated = false;
      for (int j = 0; j < combined_size; ++j) {
         if (i == j) continue;
         bool j_better_or_equal_all = true;
         bool j_strictly_better_one = false;
         for (int k = 0; k < num_objs; ++k) {
            if (pool_fit(j, k) > pool_fit(i, k)) j_better_or_equal_all = false;
            if (pool_fit(j, k) < pool_fit(i, k)) j_strictly_better_one = true;
         }
         if (j_better_or_equal_all && j_strictly_better_one) {
            dominated = true;
            break;
         }
      }
      if (!dominated) front1.push_back(i);
   }
   
   // 3. Construct the next generation
   NumericMatrix next_gen(pop_size, num_genes);
   int count = 0;
   
   // First, keep the entire Pareto Front (The Elites)
   for (int idx : front1) {
      if (count >= pop_size) break;
      next_gen(count++, _) = pool_pop(idx, _);
   }
   
   // Fill the remaining slots with children (or parents) randomly to maintain size
   while (count < pop_size) {
      int rand_idx = int(R::unif_rand() * combined_size);
      next_gen(count++, _) = pool_pop(rand_idx, _);
   }
   
   return next_gen;
}
