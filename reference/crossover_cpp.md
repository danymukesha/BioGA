# Function to perform simulated binary crossover (SBX)

Function to perform simulated binary crossover (SBX)

## Usage

``` r
crossover_cpp(
  selected_parents,
  offspring_size,
  crossover_rate = 0.9,
  eta_c = 20
)
```

## Arguments

- selected_parents:

  Numeric matrix of selected individuals.

- offspring_size:

  Number of offspring to generate.

- crossover_rate:

  Probability of crossover (default: 0.9).

- eta_c:

  SBX distribution index (default: 20.0).

## Value

Numeric matrix of offspring.

## Examples

``` r
genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
population <- BioGA::initialize_population_cpp(genomic_data, 
population_size = 5)
fitness <- BioGA::evaluate_fitness_cpp(genomic_data, population, 
c(1.0, 0.5))
selected_parents <- BioGA::selection_cpp(population, fitness, 
num_parents = 2)
#> Current front size: 1
BioGA::crossover_cpp(selected_parents, offspring_size = 2)
#>            [,1]      [,2]       [,3]       [,4]       [,5]      [,6]     [,7]
#> [1,] -0.8366399 0.4341165 -0.5434423 -0.1996363 -0.9343821 0.7124227 1.201867
#> [2,] -0.8366399 0.4341165 -0.5434423 -0.1996363 -0.9343821 0.7124227 1.201867
#>           [,8]      [,9]      [,10]
#> [1,] 0.2879157 0.1691907 -0.9783576
#> [2,] 0.2879157 0.1691907 -0.9783576
```
