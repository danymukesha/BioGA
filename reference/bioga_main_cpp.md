# Main genetic algorithm loop for genomic data optimization

Main genetic algorithm loop for genomic data optimization

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
  network = NULL
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

  SBX distribution index (default: 20.0).

- mutation_rate:

  Base probability of mutation.

- num_parents:

  Number of parents to select per generation.

- num_offspring:

  Number of offspring to generate per generation.

- num_to_replace:

  Number of individuals to replace per generation.

- weights:

  Numeric vector of weights for multi-objective fitness.

- seed:

  Optional random seed for reproducibility.

- clusters:

  Optional vector of gene cluster assignments.

- network:

  Optional matrix of gene network constraints.

## Value

List containing final population and fitness scores.

## Examples

``` r
genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
result <- BioGA::bioga_main_cpp(genomic_data,
population_size = 50, num_generations = 10,
        crossover_rate = 0.9, eta_c = 20.0, mutation_rate = 0.1,
        num_parents = 20, num_offspring = 20, num_to_replace = 10,
        weights = c(1.0, 0.5), seed = 123)
#> Current front size: 1
#> Current front size: 1
#> Current front size: 1
#> Current front size: 1
#> Current front size: 1
#> Current front size: 1
#> Current front size: 1
#> Current front size: 0
#> Warning: No non-dominated individuals found.     Using full population for selection.
#> Current front size: 1
#> Current front size: 1
```
