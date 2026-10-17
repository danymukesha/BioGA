#ifndef BIOGA_NSGA2_H
#define BIOGA_NSGA2_H

#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <limits>

// Shared NSGA-II utilities used by selection, replacement and the main loop.
// All objectives are treated as minimisation objectives.
namespace bioga {

// i dominates j iff i is no worse on every objective and strictly better on
// at least one objective.
inline bool dominates(const Rcpp::NumericMatrix& fitness, int i, int j) {
  const int n_obj = fitness.ncol();
  bool strictly_better = false;
  for (int k = 0; k < n_obj; ++k) {
    const double fi = fitness(i, k);
    const double fj = fitness(j, k);
    if (fi > fj) return false;
    if (fi < fj) strictly_better = true;
  }
  return strictly_better;
}

// Layered non-dominated sort. Fills rank (0 = best front) and fronts.
inline void fast_nondominated_sort(const Rcpp::NumericMatrix& fitness,
                                   std::vector<int>& rank,
                                   std::vector<std::vector<int> >& fronts) {
  const int n = fitness.nrow();
  rank.assign(n, 0);
  fronts.clear();

  std::vector<std::vector<int> > dominated_by(n);
  std::vector<int> domination_count(n, 0);

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (i == j) continue;
      if (dominates(fitness, i, j)) {
        dominated_by[i].push_back(j);
      } else if (dominates(fitness, j, i)) {
        domination_count[i]++;
      }
    }
  }

  std::vector<int> current;
  for (int i = 0; i < n; ++i) {
    if (domination_count[i] == 0) current.push_back(i);
  }

  int level = 0;
  while (!current.empty()) {
    fronts.push_back(current);
    std::vector<int> next;
    for (std::size_t a = 0; a < current.size(); ++a) {
      const int i = current[a];
      for (std::size_t b = 0; b < dominated_by[i].size(); ++b) {
        const int j = dominated_by[i][b];
        if (--domination_count[j] == 0) {
          rank[j] = level + 1;
          next.push_back(j);
        }
      }
    }
    current = next;
    ++level;
  }
}

// Crowding distance for a single front, aligned with the front vector.
inline std::vector<double> crowding_distance(const Rcpp::NumericMatrix& fitness,
                                             const std::vector<int>& front) {
  const int l = static_cast<int>(front.size());
  const int n_obj = fitness.ncol();
  std::vector<double> distance(l, 0.0);
  if (l <= 2) {
    for (int a = 0; a < l; ++a) distance[a] = std::numeric_limits<double>::infinity();
    return distance;
  }

  for (int k = 0; k < n_obj; ++k) {
    std::vector<int> order(l);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
              [&](int a, int b) { return fitness(front[a], k) < fitness(front[b], k); });

    const double fmin = fitness(front[order.front()], k);
    const double fmax = fitness(front[order.back()], k);
    distance[order.front()] = std::numeric_limits<double>::infinity();
    distance[order.back()] = std::numeric_limits<double>::infinity();

    if (fmax - fmin > 1e-12) {
      for (int a = 1; a < l - 1; ++a) {
        distance[order[a]] +=
          (fitness(front[order[a + 1]], k) - fitness(front[order[a - 1]], k)) / (fmax - fmin);
      }
    }
  }
  return distance;
}

// Rank (0 = best) and crowding distance for every individual.
inline void compute_rank_and_crowding(const Rcpp::NumericMatrix& fitness,
                                      std::vector<int>& rank,
                                      std::vector<double>& crowding) {
  std::vector<std::vector<int> > fronts;
  fast_nondominated_sort(fitness, rank, fronts);
  crowding.assign(fitness.nrow(), 0.0);
  for (std::size_t f = 0; f < fronts.size(); ++f) {
    std::vector<double> cd = crowding_distance(fitness, fronts[f]);
    for (std::size_t a = 0; a < fronts[f].size(); ++a) {
      crowding[fronts[f][a]] = cd[a];
    }
  }
}

// Return the k best individuals (indices) by (rank ascending, crowding descending).
inline std::vector<int> environmental_selection(const Rcpp::NumericMatrix& fitness, int k) {
  const int n = fitness.nrow();
  std::vector<int> rank;
  std::vector<double> crowding;
  compute_rank_and_crowding(fitness, rank, crowding);

  std::vector<int> ids(n);
  std::iota(ids.begin(), ids.end(), 0);
  std::stable_sort(ids.begin(), ids.end(), [&](int a, int b) {
    if (rank[a] != rank[b]) return rank[a] < rank[b];
    return crowding[a] > crowding[b];
  });

  if (k < n) ids.resize(k);
  return ids;
}

}  // namespace bioga

#endif
