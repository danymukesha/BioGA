Working Example Fix
Here’s a revised and safe version of selection_cpp() that ensures:
    
    No crash if current_front is empty.

Falls back to selecting random individuals if no front is found.

Maintains basic NSGA-II selection logic.

```cpp
#include <Rcpp.h>
using namespace Rcpp;

//' NSGA-II selection function with safe fallback
//'
//' @param population NumericMatrix of individuals.
//' @param fitness NumericMatrix of fitness scores.
//' @param num_parents Integer number of parents to select.
//' @return NumericMatrix of selected individuals.
//' @export
// [[Rcpp::export]]
NumericMatrix selection_cpp(const NumericMatrix& population, 
                            const NumericMatrix& fitness, 
                            int num_parents) {

  int population_size = population.nrow();
  int num_genes = population.ncol();
  int num_objectives = fitness.ncol();

  // Output matrix for selected parents
  NumericMatrix selected_parents(num_parents, num_genes);

  // Non-dominated sorting: very simple version
  std::vector<int> dominated_count(population_size, 0);
  std::vector<std::vector<int>> dominates(population_size);
  std::vector<int> front;

  for (int i = 0; i < population_size; ++i) {
    for (int j = 0; j < population_size; ++j) {
      if (i == j) continue;
      bool i_dominates_j = true;
      bool j_dominates_i = true;

      for (int k = 0; k < num_objectives; ++k) {
        if (fitness(i, k) > fitness(j, k)) i_dominates_j = false;
        if (fitness(j, k) > fitness(i, k)) j_dominates_i = false;
      }

      if (i_dominates_j && !j_dominates_i) {
        dominates[i].push_back(j);
      } else if (j_dominates_i && !i_dominates_j) {
        dominated_count[i]++;
      }
    }
  }

  // Build the first front
  for (int i = 0; i < population_size; ++i) {
    if (dominated_count[i] == 0) {
      front.push_back(i);
    }
  }

  Rcpp::Rcout << "Current front size: " << front.size() << std::endl;

  // Fallback: use whole population if front is empty
  if (front.empty()) {
    Rcpp::Rcout << "Fallback: using random individuals from full population.\n";
    for (int i = 0; i < population_size; ++i) front.push_back(i);
  }

  // Tournament selection from front
  for (int i = 0; i < num_parents; ++i) {
    int idx1 = front[int(R::unif_rand() * front.size())];
    int idx2 = front[int(R::unif_rand() * front.size())];
    int selected = (fitness(idx1, 0) < fitness(idx2, 0)) ? idx1 : idx2;

    for (int j = 0; j < num_genes; ++j) {
      selected_parents(i, j) = population(selected, j);
    }
  }

  return selected_parents;
}
```

 Why this works
It ensures selection proceeds even when no individuals dominate others (common in random data).

Maintains evolutionary pressure via tournament selection (using first fitness objective).

Avoids crashing RStudio with out-of-range vector accesses.