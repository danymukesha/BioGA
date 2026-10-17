# Select parents using NSGA-II rank and crowding distance

Ranks the population into non-dominated fronts, computes crowding
distances, and performs binary tournament selection that prefers lower
rank and, on ties, higher crowding distance.

## Usage

``` r
selection_cpp(population, fitness, num_parents)
```

## Arguments

- population:

  Numeric matrix of individuals.

- fitness:

  Numeric matrix of fitness scores (rows: individuals, columns:
  objectives, all minimised).

- num_parents:

  Number of individuals to select.

## Value

Numeric matrix of selected individuals.

## References

Deb K, Pratap A, Agarwal S, Meyarivan T (2002). "A fast and elitist
multiobjective genetic algorithm: NSGA-II." *IEEE Transactions on
Evolutionary Computation*, 6(2), 182-197.
[doi:10.1109/4235.996017](https://doi.org/10.1109/4235.996017)

## See also

[`initialize_population_cpp()`](https://danymukesha.github.io/BioGA/reference/initialize_population_cpp.md),
[`evaluate_fitness_cpp()`](https://danymukesha.github.io/BioGA/reference/evaluate_fitness_cpp.md),
[`crossover_cpp()`](https://danymukesha.github.io/BioGA/reference/crossover_cpp.md),
[`mutation_cpp()`](https://danymukesha.github.io/BioGA/reference/mutation_cpp.md),
[`replacement_cpp()`](https://danymukesha.github.io/BioGA/reference/replacement_cpp.md),
[`bioga_main_cpp()`](https://danymukesha.github.io/BioGA/reference/bioga_main_cpp.md)

## Examples

``` r
genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
population <- BioGA::initialize_population_cpp(genomic_data, population_size = 5)
fitness <- BioGA::evaluate_fitness_cpp(genomic_data, population, c(1.0, 0.5))
BioGA::selection_cpp(population, fitness, num_parents = 2)
#>             [,1]      [,2]       [,3]      [,4]     [,5]      [,6]       [,7]
#> [1,]  1.09116354 0.0700649 -1.0061661 0.2912297 1.385545 1.4358605 -0.4794401
#> [2,] -0.04448748 0.6941418 -0.3148323 0.4334064 1.385545 0.6116994 -0.2336441
#>           [,8]       [,9]      [,10]
#> [1,] -1.585889 -0.3838501  0.9133362
#> [2,] -2.679716 -0.2075835 -0.7181257
```
