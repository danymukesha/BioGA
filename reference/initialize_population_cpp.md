# Initialize a population, optionally guided by gene clusters

Creates the starting population. In `"real"` mode each gene value is
drawn from an observed sample of that gene, so individuals preserve the
marginal distribution of the data. In `"binary"` mode each gene is
encoded as selected (1) or not selected (0).

## Usage

``` r
initialize_population_cpp(
  genomic_data,
  population_size,
  seed = NULL,
  clusters = NULL,
  representation = "real"
)
```

## Arguments

- genomic_data:

  Numeric matrix of genomic data (rows: genes, columns: samples).

- population_size:

  Number of individuals in the population.

- seed:

  Optional random seed for reproducibility.

- clusters:

  Optional vector of gene cluster assignments (length equal to the
  number of genes, non-negative integers).

- representation:

  Either `"real"` (continuous expression values) or `"binary"` (gene
  selection indicators).

## Value

Numeric matrix of initialized population (rows: individuals, columns:
genes).

## See also

[`evaluate_fitness_cpp()`](https://danymukesha.github.io/BioGA/reference/evaluate_fitness_cpp.md),
[`selection_cpp()`](https://danymukesha.github.io/BioGA/reference/selection_cpp.md),
[`crossover_cpp()`](https://danymukesha.github.io/BioGA/reference/crossover_cpp.md),
[`mutation_cpp()`](https://danymukesha.github.io/BioGA/reference/mutation_cpp.md),
[`replacement_cpp()`](https://danymukesha.github.io/BioGA/reference/replacement_cpp.md),
[`bioga_main_cpp()`](https://danymukesha.github.io/BioGA/reference/bioga_main_cpp.md)

## Examples

``` r
genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
BioGA::initialize_population_cpp(genomic_data, population_size = 5, seed = 123)
#>             [,1]       [,2]     [,3]       [,4]        [,5]       [,6]
#> [1,]  0.06270301  1.6832769 1.463565 -0.1351125 -1.05541631 -0.2394497
#> [2,] -1.49079011 -0.4770160 2.558912 -0.1595557 -0.04804319 -1.9235357
#> [3,]  0.61886280 -0.4287279 2.558912  0.5639503 -1.71865548  0.5857998
#> [4,] -1.49079011 -0.3890261 2.558912 -1.6462332 -0.44972811  0.4242909
#> [5,]  0.24754278 -0.4770160 1.463565 -2.0166358 -0.04804319  0.3088376
#>            [,7]        [,8]       [,9]      [,10]
#> [1,] -0.4723847  0.71806477 -0.4393085 -0.6501188
#> [2,]  1.0104241  0.53409836 -0.1590448 -0.9046839
#> [3,] -0.4723847  1.49388563 -1.4871922  0.7206758
#> [4,] -0.6316866 -0.40966372 -0.1590448  0.2205800
#> [5,]  1.0104241 -0.05790095 -1.4871922  2.7254483
```
