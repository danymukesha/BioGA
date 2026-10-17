# Build the next population with NSGA-II environmental selection

The `num_to_replace` best offspring (by NSGA-II rank and crowding
distance) replace the lowest-quality parents, keeping the best
`population_size - num_to_replace` parents. This preserves elitism while
letting offspring enter only on merit. If `num_to_replace <= 0` the
population is unchanged.

## Usage

``` r
replacement_cpp(
  population,
  offspring,
  fitness,
  offspring_fitness,
  num_to_replace
)
```

## Arguments

- population:

  Numeric matrix of current individuals.

- offspring:

  Numeric matrix of offspring.

- fitness:

  Numeric matrix of current population fitness scores.

- offspring_fitness:

  Numeric matrix of offspring fitness scores.

- num_to_replace:

  Maximum number of offspring to admit.

## Value

Numeric matrix of the updated population.

## See also

[`initialize_population_cpp()`](https://danymukesha.github.io/BioGA/reference/initialize_population_cpp.md),
[`evaluate_fitness_cpp()`](https://danymukesha.github.io/BioGA/reference/evaluate_fitness_cpp.md),
[`selection_cpp()`](https://danymukesha.github.io/BioGA/reference/selection_cpp.md),
[`crossover_cpp()`](https://danymukesha.github.io/BioGA/reference/crossover_cpp.md),
[`mutation_cpp()`](https://danymukesha.github.io/BioGA/reference/mutation_cpp.md),
[`bioga_main_cpp()`](https://danymukesha.github.io/BioGA/reference/bioga_main_cpp.md)

## Examples

``` r
genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
population <- BioGA::initialize_population_cpp(genomic_data, population_size = 5)
fitness <- BioGA::evaluate_fitness_cpp(genomic_data, population, c(1.0, 0.5))
offspring <- population
offspring_fitness <- fitness
BioGA::replacement_cpp(population, offspring, fitness, offspring_fitness,
  num_to_replace = 1)
#>            [,1]       [,2]       [,3]       [,4]       [,5]       [,6]
#> [1,]  1.2162502 -0.1395003  0.7395213 -0.5991232 -0.3460009 -0.8377676
#> [2,]  0.7915479  1.7971404  0.9512089  0.5384658  0.2118099 -0.5922097
#> [3,] -1.5701355  0.2197284  0.7395213 -1.1516032  0.7089544  0.4894620
#> [4,]  2.7081683  0.6403630 -1.2755434 -1.9111513 -0.3412932 -0.8377676
#> [5,]  1.2162502 -0.1395003  0.7395213 -0.5991232 -0.3460009 -0.8377676
#>            [,7]       [,8]       [,9]      [,10]
#> [1,]  1.3657766  0.8864217 -0.2919208  0.4890803
#> [2,]  1.3657766 -0.2840014 -0.4088598  0.1602164
#> [3,] -0.6104617  1.0553386 -0.2919208 -0.8241094
#> [4,] -1.1963517  0.4627379 -0.8047323  1.3168057
#> [5,]  1.3657766  0.8864217 -0.2919208  0.4890803
```
