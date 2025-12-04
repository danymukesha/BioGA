# Function to select individuals using NSGA-II non-dominated sorting

Function to select individuals using NSGA-II non-dominated sorting

## Usage

``` r
selection_cpp(population, fitness, num_parents)
```

## Arguments

- population:

  Numeric matrix of individuals.

- fitness:

  Numeric matrix of fitness scores (columns: objectives).

- num_parents:

  Number of individuals to select.

## Value

Numeric matrix of selected individuals.

## Examples

``` r
genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
population <- BioGA::initialize_population_cpp(genomic_data, 
     population_size = 5)
fitness <- BioGA::evaluate_fitness_cpp(genomic_data, population, 
     c(1.0, 0.5))
BioGA::selection_cpp(population, fitness, num_parents = 2)
#> Current front size: 1
#>            [,1]      [,2]     [,3]      [,4]      [,5]       [,6]      [,7]
#> [1,] -0.4356455 0.3461036 2.416773 -1.651049 0.1263159 0.05601673 0.5101325
#> [2,] -0.4356455 0.3461036 2.416773 -1.651049 0.1263159 0.05601673 0.5101325
#>            [,8]       [,9]     [,10]
#> [1,] -0.4503386 -0.3497542 0.7564064
#> [2,] -0.4503386 -0.3497542 0.7564064
```
