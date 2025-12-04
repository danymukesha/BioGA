# Function to evaluate fitness using genomic data with multi-objective support

Function to evaluate fitness using genomic data with multi-objective
support

## Usage

``` r
evaluate_fitness_cpp(genomic_data, population, weights)
```

## Arguments

- genomic_data:

  Numeric matrix of genomic data (rows: genes, columns: samples).

- population:

  Numeric matrix representing the population of individuals.

- weights:

  Numeric vector of weights for multi-objective fitness (e.g.,
  expression difference, sparsity).

## Value

Numeric matrix of fitness scores (columns: objectives, rows:
individuals).

## Examples

``` r
genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
population <- BioGA::initialize_population_cpp(genomic_data, 
population_size = 5)
weights <- c(1.0, 0.5) # Weight for expression difference and sparsity
BioGA::evaluate_fitness_cpp(genomic_data, population, weights)
#>          [,1] [,2]
#> [1,] 174.1882  0.5
#> [2,] 185.6447  0.5
#> [3,] 282.9144  0.5
#> [4,] 260.4805  0.5
#> [5,] 239.6672  0.5
```
