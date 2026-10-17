# Mutate offspring with an adaptive rate

The mutation rate increases linearly with the generation index.
Real-coded mode applies Gaussian perturbation; binary mode flips
selection bits. When a gene network is supplied, the perturbation of
each gene is damped by its relative connectivity (degree), so network
hubs change less.

## Usage

``` r
mutation_cpp(
  offspring,
  mutation_rate,
  iteration,
  max_iterations,
  network = NULL,
  representation = "real"
)
```

## Arguments

- offspring:

  Numeric matrix of offspring.

- mutation_rate:

  Base probability of mutation.

- iteration:

  Current GA iteration.

- max_iterations:

  Maximum number of GA iterations.

- network:

  Optional symmetric adjacency matrix (rows: genes, columns: genes).

- representation:

  Either `"real"` or `"binary"`.

## Value

Numeric matrix of mutated offspring.

## See also

[`initialize_population_cpp()`](https://danymukesha.github.io/BioGA/reference/initialize_population_cpp.md),
[`evaluate_fitness_cpp()`](https://danymukesha.github.io/BioGA/reference/evaluate_fitness_cpp.md),
[`selection_cpp()`](https://danymukesha.github.io/BioGA/reference/selection_cpp.md),
[`crossover_cpp()`](https://danymukesha.github.io/BioGA/reference/crossover_cpp.md),
[`replacement_cpp()`](https://danymukesha.github.io/BioGA/reference/replacement_cpp.md),
[`bioga_main_cpp()`](https://danymukesha.github.io/BioGA/reference/bioga_main_cpp.md)

## Examples

``` r
genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
population <- BioGA::initialize_population_cpp(genomic_data, population_size = 5)
fitness <- BioGA::evaluate_fitness_cpp(genomic_data, population, c(1.0, 0.5))
parents <- BioGA::selection_cpp(population, fitness, num_parents = 2)
offspring <- BioGA::crossover_cpp(parents, offspring_size = 2)
BioGA::mutation_cpp(offspring, mutation_rate = 0.1, iteration = 1,
  max_iterations = 100)
#>            [,1]      [,2]      [,3]      [,4]     [,5]      [,6]       [,7]
#> [1,] -0.6002596 0.5539177 0.5846137 0.7799651 1.253815 0.3796395 -0.2950715
#> [2,] -0.6002596 0.5539177 0.5846137 0.7799651 1.253815 0.3796395 -0.3186994
#>             [,8]      [,9]    [,10]
#> [1,] -0.04287046 0.8781335 1.360652
#> [2,] -0.04287046 0.8781335 1.360652
```
