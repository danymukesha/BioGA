#ifndef BIOGA_H
#define BIOGA_H

#include <Rcpp.h>
using namespace Rcpp;


// I decided to put here 4declaration only — no body here

NumericMatrix initialize_population_cpp(const NumericMatrix& genomic_data,
                                        int population_size,
                                        Nullable<int> seed = R_NilValue,
                                        Nullable<IntegerVector> clusters = R_NilValue);

NumericMatrix evaluate_fitness_cpp(const NumericMatrix& genomic_data,
                                   const NumericMatrix& population,
                                   const NumericVector& weights);

NumericMatrix selection_cpp(const NumericMatrix& population,
                            const NumericMatrix& fitness,
                            int num_parents);

NumericMatrix crossover_cpp(const NumericMatrix& parents,
                            int num_offspring,
                            double crossover_rate,
                            double eta_c);

NumericMatrix mutation_cpp(const NumericMatrix& offspring,
                           double mutation_rate,
                           int generation,
                           int num_generations,
                           Nullable<NumericMatrix> network = R_NilValue);

NumericMatrix replacement_cpp(const NumericMatrix& population,
                              const NumericMatrix& offspring,
                              const NumericMatrix& fitness,
                              const NumericMatrix& offspring_fitness,
                              int num_to_replace);

#endif 
