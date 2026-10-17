#ifndef BIOGA_H
#define BIOGA_H

#include <Rcpp.h>
#include <string>
using namespace Rcpp;

// Declarations only. Definitions (with default arguments) live in the
// corresponding .cpp files.

// Encode class labels as consecutive integers. Accepts integer, numeric,
// factor, character or logical vectors (mapped by unique level as with
// as.integer(as.factor(x))).
inline IntegerVector bioga_encode_labels(SEXP labels) {
  if (Rf_isNull(labels))
    stop("labels must be provided when representation = 'binary'");
  RObject x(labels);
  Function as_factor("as.factor");
  Function as_integer("as.integer");
  return as_integer(as_factor(x));
}

NumericMatrix initialize_population_cpp(const NumericMatrix& genomic_data,
                                        int population_size,
                                        Nullable<int> seed,
                                        Nullable<IntegerVector> clusters,
                                        std::string representation);

NumericMatrix evaluate_fitness_cpp(const NumericMatrix& genomic_data,
                                   const NumericMatrix& population,
                                   const NumericVector& weights,
                                   SEXP labels,
                                   std::string representation);

NumericMatrix selection_cpp(const NumericMatrix& population,
                            const NumericMatrix& fitness,
                            int num_parents);

NumericMatrix crossover_cpp(const NumericMatrix& selected_parents,
                            int offspring_size,
                            double crossover_rate,
                            double eta_c,
                            std::string representation);

NumericMatrix mutation_cpp(const NumericMatrix& offspring,
                           double mutation_rate,
                           int iteration,
                           int max_iterations,
                           Nullable<NumericMatrix> network,
                           std::string representation);

NumericMatrix replacement_cpp(const NumericMatrix& population,
                              const NumericMatrix& offspring,
                              const NumericMatrix& fitness,
                              const NumericMatrix& offspring_fitness,
                              int num_to_replace);

#endif
