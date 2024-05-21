#include <cstdlib> // Include the header for rand() function
#include <Rcpp.h>
using namespace Rcpp;

//' Function to initialize the population from genomic data
//'
//' @param genomic_data Numeric matrix of genomic data where rows represent
//' genes/features and columns represent samples.
//' @param population_size Number of individuals in the population.
//' @return Numeric matrix representing the initialized population.
//' @examples
//' # example of usage
//' genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
//' BioGA::initialize_population_cpp(genomic_data, population_size = 5)
//' @export
// [[Rcpp::export]]
 NumericMatrix initialize_population_cpp(const NumericMatrix& genomic_data, int population_size) {
   int num_genes = genomic_data.nrow();
   int num_samples = genomic_data.ncol();

   // Allocate memory for the population matrix
   NumericMatrix population(population_size, num_genes);

   // Seed the random number generator
   // srand(time(NULL));

   // Generate random population using genomic data
   for (int i = 0; i < population_size; ++i) {
     for (int j = 0; j < num_genes; ++j) {
       int sample_index = int(R::unif_rand()) % num_samples;
       population(i, j) = genomic_data(j, sample_index);
     }
   }

   return population;
 }
