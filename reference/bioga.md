# Run the BioGA genetic algorithm

User-facing wrapper around
[`bioga_main_cpp()`](https://danymukesha.github.io/BioGA/reference/bioga_main_cpp.md)
with sensible defaults and input validation. The algorithm minimises two
objectives; in `"binary"` mode pass class `labels` so that the first
objective measures the predictive signal of the selected genes.

## Usage

``` r
bioga(
  genomic_data,
  population_size = 50,
  num_generations = 20,
  crossover_rate = 0.9,
  mutation_rate = 0.1,
  num_parents = 20,
  num_offspring = 20,
  num_to_replace = 10,
  weights = c(1, 0.5),
  labels = NULL,
  clusters = NULL,
  network = NULL,
  seed = NULL,
  eta_c = 20,
  representation = c("real", "binary")
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

- mutation_rate:

  Base probability of mutation.

- num_parents:

  Number of parents selected per generation.

- num_offspring:

  Number of offspring generated per generation.

- num_to_replace:

  Maximum number of offspring admitted per generation.

- weights:

  Numeric vector of objective weights (length 1 or 2).

- labels:

  Optional class labels, one per sample (integer, numeric, factor,
  character or logical). Required when `representation = "binary"`.

- clusters:

  Optional vector of gene cluster assignments.

- network:

  Optional symmetric gene network adjacency matrix.

- seed:

  Optional integer random seed for reproducibility.

- eta_c:

  SBX distribution index (default: 20).

- representation:

  Either `"real"` or `"binary"`.

## Value

A list with `population`, `fitness`, `history`, `best`, `best_fitness`,
`pareto_population` and `pareto_fitness`.

## See also

[`bioga_fitness()`](https://danymukesha.github.io/BioGA/reference/bioga_fitness.md),
[`bioga_main_cpp()`](https://danymukesha.github.io/BioGA/reference/bioga_main_cpp.md)

## Examples

``` r
genomic_data <- matrix(rnorm(200), nrow = 20, ncol = 10)
result <- bioga(genomic_data, population_size = 20, num_generations = 5)
```
