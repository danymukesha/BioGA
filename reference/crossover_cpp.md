# Generate offspring by crossover

Real-coded mode uses simulated binary crossover (SBX); binary mode uses
uniform crossover. When only one parent is available the offspring are
exact copies of that parent.

## Usage

``` r
crossover_cpp(
  selected_parents,
  offspring_size,
  crossover_rate = 0.9,
  eta_c = 20,
  representation = "real"
)
```

## Arguments

- selected_parents:

  Numeric matrix of selected individuals.

- offspring_size:

  Number of offspring to generate.

- crossover_rate:

  Probability of crossover per offspring (default: 0.9).

- eta_c:

  SBX distribution index (default: 20.0).

- representation:

  Either `"real"` or `"binary"`.

## Value

Numeric matrix of offspring.

## References

Deb K, Agrawal RB (1995). "Simulated binary crossover for continuous
search space." *Complex Systems*, 9(2), 115-148.

## See also

[`initialize_population_cpp()`](https://danymukesha.github.io/BioGA/reference/initialize_population_cpp.md),
[`evaluate_fitness_cpp()`](https://danymukesha.github.io/BioGA/reference/evaluate_fitness_cpp.md),
[`selection_cpp()`](https://danymukesha.github.io/BioGA/reference/selection_cpp.md),
[`mutation_cpp()`](https://danymukesha.github.io/BioGA/reference/mutation_cpp.md),
[`replacement_cpp()`](https://danymukesha.github.io/BioGA/reference/replacement_cpp.md),
[`bioga_main_cpp()`](https://danymukesha.github.io/BioGA/reference/bioga_main_cpp.md)

## Examples

``` r
genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
population <- BioGA::initialize_population_cpp(genomic_data, population_size = 5)
fitness <- BioGA::evaluate_fitness_cpp(genomic_data, population, c(1.0, 0.5))
parents <- BioGA::selection_cpp(population, fitness, num_parents = 2)
BioGA::crossover_cpp(parents, offspring_size = 2)
#>           [,1]       [,2]      [,3]      [,4]       [,5]       [,6]      [,7]
#> [1,] 0.8606541 -0.8968653 0.1844905 0.4581388 -0.4174473 -0.9783576 0.2208403
#> [2,] 0.8606541 -0.8968653 0.1844905 0.4581388 -0.4174473 -0.9783576 0.2208403
#>           [,8]       [,9]    [,10]
#> [1,] 0.3136316 -0.1592713 0.684092
#> [2,] 0.3136316 -0.1592713 0.684092
```
