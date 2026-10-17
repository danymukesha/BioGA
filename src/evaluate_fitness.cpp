#include <RcppParallel.h>
#include <Rcpp.h>
#include <map>
#include <vector>
#include <cmath>
#include "BioGA.h"
using namespace Rcpp;
using namespace RcppParallel;

//' Evaluate fitness with two competing objectives
//'
//' Objectives are always minimised. Weighted columns are returned so the
//' caller can inspect each objective (`weights` are recycled when shorter than
//' the number of objectives).
//'
//' In `"real"` mode individuals are continuous expression profiles:
//' \itemize{
//'   \item `f1`: squared Euclidean distance between the individual and every
//'     observed sample (data fidelity).
//'   \item `f2`: mean absolute gene value (L1 parsimony).
//' }
//'
//' In `"binary"` mode individuals are gene-selection indicators and `labels`
//' (one class label per sample) are required:
//' \itemize{
//'   \item `f1`: negated mean ANOVA F-statistic of the selected genes
//'     (predictive signal).
//'   \item `f2`: fraction of selected genes (cardinality / sparsity).
//' }
//'
//' @param genomic_data Numeric matrix of genomic data (rows: genes,
//' columns: samples).
//' @param population Numeric matrix of individuals (rows: individuals,
//' columns: genes).
//' @param weights Numeric vector of objective weights (recycled if shorter
//' than two).
//' @param labels Optional class labels (length == number of samples); an
//' integer, numeric, factor, character or logical vector is accepted and
//' mapped to integer class codes. Required for `representation = "binary"`.
//' @param representation Either `"real"` or `"binary"`.
//' @return Numeric matrix of fitness scores (rows: individuals,
//' columns: objectives).
//' @family fitness evaluation
//' @seealso [initialize_population_cpp()], [selection_cpp()], [crossover_cpp()],
//'   [mutation_cpp()], [replacement_cpp()], [bioga_main_cpp()]
//' @examples
//' genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
//' population <- BioGA::initialize_population_cpp(genomic_data, population_size = 5)
//' BioGA::evaluate_fitness_cpp(genomic_data, population, weights = c(1.0, 0.5))
//' @export
// [[Rcpp::export]]
NumericMatrix evaluate_fitness_cpp(const NumericMatrix& genomic_data,
                                   const NumericMatrix& population,
                                   const NumericVector& weights,
                                   SEXP labels = R_NilValue,
                                   std::string representation = "real") {
  if (genomic_data.nrow() == 0 || genomic_data.ncol() == 0)
    stop("genomic_data must have at least 1 row and 1 column");
  if (population.nrow() == 0)
    stop("population must have at least 1 individual");
  if (population.ncol() != genomic_data.nrow())
    stop("population columns must equal genomic_data rows (number of genes)");
  if (weights.length() == 0)
    stop("weights must not be empty");
  if (representation != "real" && representation != "binary")
    stop("representation must be either 'real' or 'binary'");

  const bool binary = (representation == "binary");
  const int num_genes = genomic_data.nrow();
  const int num_samples = genomic_data.ncol();
  const int population_size = population.nrow();
  const int num_objectives = 2;

  const double w1 = weights[0];
  const double w2 = (weights.length() > 1) ? weights[1] : weights[0];

  NumericMatrix fitness(population_size, num_objectives);

  if (!binary) {
    struct FitnessWorker : public Worker {
      const RMatrix<double> genomic_data;
      const RMatrix<double> population;
      RMatrix<double> fitness;
      const double w1;
      const double w2;
      const int num_genes;

      FitnessWorker(const NumericMatrix& genomic_data,
                    const NumericMatrix& population,
                    NumericMatrix& fitness,
                    double w1, double w2, int num_genes)
        : genomic_data(genomic_data), population(population),
          fitness(fitness), w1(w1), w2(w2), num_genes(num_genes) {}

      void operator()(std::size_t begin, std::size_t end) {
        const int n_samples = static_cast<int>(genomic_data.ncol());
        for (std::size_t i = begin; i < end; ++i) {
          double sse = 0.0;
          double l1 = 0.0;
          for (int j = 0; j < num_genes; ++j) {
            const double value = population(i, j);
            for (int k = 0; k < n_samples; ++k) {
              const double diff = genomic_data(j, k) - value;
              sse += diff * diff;
            }
            l1 += std::fabs(value);
          }
          fitness(i, 0) = w1 * sse;
          fitness(i, 1) = w2 * l1 / num_genes;
        }
      }
    };

    FitnessWorker worker(genomic_data, population, fitness, w1, w2, num_genes);
    parallelFor(0, population_size, worker);
    return fitness;
  }

  // Binary representation: supervised gene scores.
  IntegerVector y = bioga_encode_labels(labels);
  if (y.size() != num_samples)
    stop("labels length must equal the number of samples (columns of genomic_data)");

  std::map<int, std::vector<int> > groups;
  for (int k = 0; k < num_samples; ++k) {
    groups[y[k]].push_back(k);
  }
  const int num_groups = static_cast<int>(groups.size());
  if (num_groups < 2)
    stop("labels must contain at least 2 classes for binary representation");
  if (num_samples <= num_groups)
    stop("not enough samples to estimate within-class variance");

  NumericVector gene_scores(num_genes);
  for (int j = 0; j < num_genes; ++j) {
    double overall_mean = 0.0;
    for (int k = 0; k < num_samples; ++k) overall_mean += genomic_data(j, k);
    overall_mean /= num_samples;

    double ss_between = 0.0;
    double ss_within = 0.0;
    for (std::map<int, std::vector<int> >::const_iterator it = groups.begin();
         it != groups.end(); ++it) {
      const std::vector<int>& idx = it->second;
      double group_mean = 0.0;
      for (std::size_t a = 0; a < idx.size(); ++a) group_mean += genomic_data(j, idx[a]);
      group_mean /= idx.size();
      ss_between += idx.size() * (group_mean - overall_mean) * (group_mean - overall_mean);
      for (std::size_t a = 0; a < idx.size(); ++a) {
        const double diff = genomic_data(j, idx[a]) - group_mean;
        ss_within += diff * diff;
      }
    }
    gene_scores[j] = (ss_between / (num_groups - 1)) / (ss_within / (num_samples - num_groups) + 1e-12);
  }

  for (int i = 0; i < population_size; ++i) {
    int selected = 0;
    double score_sum = 0.0;
    for (int j = 0; j < num_genes; ++j) {
      if (population(i, j) > 0.5) {
        ++selected;
        score_sum += gene_scores[j];
      }
    }
    const double f1 = (selected > 0) ? -(score_sum / selected) : 0.0;
    const double f2 = static_cast<double>(selected) / num_genes;
    fitness(i, 0) = w1 * f1;
    fitness(i, 1) = w2 * f2;
  }

  return fitness;
}
