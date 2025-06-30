#include <RcppParallel.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace RcppParallel;

//' Function to evaluate fitness using genomic data with multi-objective support
//'
//' @param genomic_data Numeric matrix of genomic data (rows: genes, columns: samples).
//' @param population Numeric matrix representing the population of individuals.
//' @param weights Numeric vector of weights for multi-objective fitness (e.g., expression difference, sparsity).
//' @return Numeric matrix of fitness scores (columns: objectives, rows: individuals).
//' @examples
//' genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
//' population <- BioGA::initialize_population_cpp(genomic_data, population_size = 5)
//' weights <- c(1.0, 0.5) # Weight for expression difference and sparsity
//' BioGA::evaluate_fitness_cpp(genomic_data, population, weights)
//' @export
// [[Rcpp::export]]
NumericMatrix evaluate_fitness_cpp(const NumericMatrix& genomic_data, 
                                  const NumericMatrix& population,
                                  const NumericVector& weights) {
 int num_genes = genomic_data.nrow();
 int num_samples = genomic_data.ncol();
 int population_size = population.nrow();
 int num_objectives = weights.length(); // Number of objectives (e.g., 2: expression diff, sparsity)
 
 // Allocate memory for fitness matrix (rows: individuals, columns: objectives)
 NumericMatrix fitness(population_size, num_objectives);
 
 // Parallelized fitness calculation
 struct FitnessWorker : public Worker {
   const RMatrix<double> genomic_data;
   const RMatrix<double> population;
   const RVector<double> weights;
   RMatrix<double> fitness;
   
   FitnessWorker(const NumericMatrix& genomic_data, 
                 const NumericMatrix& population, 
                 const NumericVector& weights,
                 NumericMatrix& fitness)
     : genomic_data(genomic_data), population(population), 
       weights(weights), fitness(fitness) {}
   
   void operator()(std::size_t begin, std::size_t end) {
     for (std::size_t i = begin; i < end; ++i) {
       double sum_fitness = 0.0;
       double sparsity = 0.0;
       // Objective 1: Sum of squared differences (expression difference)
       for (int j = 0; j < genomic_data.nrow(); ++j) {
         double sum_gene_expression = 0.0;
         for (int k = 0; k < genomic_data.ncol(); ++k) {
           sum_gene_expression += pow(genomic_data(j, k) - population(i, j), 2);
         }
         sum_fitness += sum_gene_expression;
         // Objective 2: Sparsity (count non-zero genes)
         if (fabs(population(i, j)) > 1e-6) sparsity += 1.0;
       }
       fitness(i, 0) = weights[0] * sum_fitness; // Weighted expression difference
       fitness(i, 1) = weights[1] * sparsity / genomic_data.nrow(); // Weighted sparsity
     }
   }
 };
 
 FitnessWorker worker(genomic_data, population, weights, fitness);
 parallelFor(0, population_size, worker);
 
 return fitness;
}