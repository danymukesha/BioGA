# Function to mutate offspring with adaptive mutation and network constraints

Function to mutate offspring with adaptive mutation and network
constraints

## Usage

``` r
mutation_cpp(
  offspring,
  mutation_rate,
  iteration,
  max_iterations,
  network = NULL
)
```

## Arguments

- offspring:

  Numeric matrix of offspring.

- mutation_rate:

  Base probability of mutation.

- iteration:

  Current GA iteration for adaptive mutation.

- max_iterations:

  Maximum number of GA iterations.

- network:

  Optional matrix of gene network constraints (rows: genes, cols:
  genes).

## Value

Numeric matrix of mutated offspring.

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
offspring <- BioGA::crossover_cpp(selected_parents, 
offspring_size = 2)
BioGA::mutation_cpp(offspring, mutation_rate = 0.1, iteration = 1, 
max_iterations = 100)
#>          [,1]      [,2]      [,3]       [,4]      [,5]       [,6]      [,7]
#> [1,] 1.096839 -0.640706 0.1813035 -0.1388914 0.3796395 -0.5023235 0.8951257
#> [2,] 1.096839 -0.640706 0.1813035 -0.1388914 0.3796395 -0.5259514 0.8951257
#>           [,8]      [,9]     [,10]
#> [1,] 0.6443765 0.5194072 0.6886403
#> [2,] 0.6443765 0.5194072 0.6886403
```
