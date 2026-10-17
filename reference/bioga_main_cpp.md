# Main genetic algorithm loop for genomic data optimisation

Runs a generational multi-objective genetic algorithm combining
initialization, fitness evaluation, NSGA-II selection, crossover,
mutation and NSGA-II replacement.

## Usage

``` r
bioga_main_cpp(
  genomic_data,
  population_size,
  num_generations,
  crossover_rate,
  eta_c,
  mutation_rate,
  num_parents,
  num_offspring,
  num_to_replace,
  weights,
  seed = NULL,
  clusters = NULL,
  network = NULL,
  labels = NULL,
  representation = "real"
)
```

## Arguments

- genomic_data:

  Numeric matrix of genomic data (rows: genes, columns: samples).

- population_size:

  Number of individuals in the population.

- num_generations:

  Number of generations to run.

- crossover_rate:

  Probability of crossover.

- eta_c:

  SBX distribution index (real representation).

- mutation_rate:

  Base probability of mutation.

- num_parents:

  Number of parents selected per generation.

- num_offspring:

  Number of offspring generated per generation.

- num_to_replace:

  Maximum number of offspring admitted per generation.

- weights:

  Numeric vector of objective weights.

- seed:

  Optional random seed for reproducibility.

- clusters:

  Optional vector of gene cluster assignments.

- network:

  Optional symmetric gene network adjacency matrix.

- labels:

  Optional class labels (required when `representation = "binary"`).
  Accepted as an integer, numeric, factor, character or logical vector.

- representation:

  Either `"real"` or `"binary"`.

## Value

A list with `population` (final population), `fitness` (final fitness
matrix), `history` (best first objective per generation), `best` (best
individual), `best_fitness`, `pareto_population` and `pareto_fitness`
(the first non-dominated front of the final population).

## References

Deb K, Pratap A, Agarwal S, Meyarivan T (2002). "A fast and elitist
multiobjective genetic algorithm: NSGA-II." *IEEE Transactions on
Evolutionary Computation*, 6(2), 182-197.
[doi:10.1109/4235.996017](https://doi.org/10.1109/4235.996017)

## See also

[`initialize_population_cpp()`](https://danymukesha.github.io/BioGA/reference/initialize_population_cpp.md),
[`evaluate_fitness_cpp()`](https://danymukesha.github.io/BioGA/reference/evaluate_fitness_cpp.md),
[`selection_cpp()`](https://danymukesha.github.io/BioGA/reference/selection_cpp.md),
[`crossover_cpp()`](https://danymukesha.github.io/BioGA/reference/crossover_cpp.md),
[`mutation_cpp()`](https://danymukesha.github.io/BioGA/reference/mutation_cpp.md),
[`replacement_cpp()`](https://danymukesha.github.io/BioGA/reference/replacement_cpp.md)

## Examples

``` r
genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
result <- BioGA::bioga_main_cpp(genomic_data, population_size = 50,
  num_generations = 10, crossover_rate = 0.9, eta_c = 20.0,
  mutation_rate = 0.1, num_parents = 20, num_offspring = 20,
  num_to_replace = 10, weights = c(1.0, 0.5), seed = 123)
```
