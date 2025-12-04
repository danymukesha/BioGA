# Function to replace population with elitism and diversity preservation

Function to replace population with elitism and diversity preservation

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

  Numeric matrix of individuals.

- offspring:

  Numeric matrix of offspring.

- fitness:

  Numeric matrix of population fitness scores.

- offspring_fitness:

  Numeric matrix of offspring fitness scores.

- num_to_replace:

  Number of individuals to replace.

## Value

Numeric matrix of updated population.

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
offspring <- BioGA::crossover_cpp(selected_parents, offspring_size = 2)
offspring_fitness <- BioGA::evaluate_fitness_cpp(genomic_data, offspring, 
c(1.0, 0.5))
BioGA::replacement_cpp(population, offspring, fitness, offspring_fitness, 
num_to_replace = 1)
#>            [,1]        [,2]       [,3]        [,4]       [,5]       [,6]
#> [1,]  0.8846505 -0.49929202 -0.7886220 -0.09031959  0.6843094 -1.3952743
#> [2,]  1.6509075  0.48545998 -1.6674751 -0.50219872  1.4960607 -0.8953634
#> [3,]  0.8846505 -0.49929202 -0.7886220 -0.09031959  0.6843094 -1.3952743
#> [4,] -0.6111659  0.08920722 -0.5739735  0.61798582 -0.3439172 -0.8953634
#> [5,] -0.3887799 -1.53290200 -1.6674751 -0.32468591  0.3860266 -1.3952743
#>            [,7]        [,8]       [,9]       [,10]
#> [1,]  0.8496430  0.05974994  0.1748027 -0.59461727
#> [2,]  0.1181445 -0.47624689  0.6007088  0.07455118
#> [3,]  0.8496430  0.05974994  0.1748027 -0.59461727
#> [4,] -1.3108015  0.13403865 -0.1009749 -1.75652740
#> [5,]  2.2930790  0.05974994  0.6007088 -0.71721816
```
