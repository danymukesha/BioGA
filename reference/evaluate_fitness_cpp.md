# Evaluate fitness with two competing objectives

Objectives are always minimised. Weighted columns are returned so the
caller can inspect each objective (`weights` are recycled when shorter
than the number of objectives).

## Usage

``` r
evaluate_fitness_cpp(
  genomic_data,
  population,
  weights,
  labels = NULL,
  representation = "real"
)
```

## Arguments

- genomic_data:

  Numeric matrix of genomic data (rows: genes, columns: samples).

- population:

  Numeric matrix of individuals (rows: individuals, columns: genes).

- weights:

  Numeric vector of objective weights (recycled if shorter than two).

- labels:

  Optional class labels (length == number of samples); an integer,
  numeric, factor, character or logical vector is accepted and mapped to
  integer class codes. Required for `representation = "binary"`.

- representation:

  Either `"real"` or `"binary"`.

## Value

Numeric matrix of fitness scores (rows: individuals, columns:
objectives).

## Details

In `"real"` mode individuals are continuous expression profiles:

- `f1`: squared Euclidean distance between the individual and every
  observed sample (data fidelity).

- `f2`: mean absolute gene value (L1 parsimony).

In `"binary"` mode individuals are gene-selection indicators and
`labels` (one class label per sample) are required:

- `f1`: negated mean ANOVA F-statistic of the selected genes (predictive
  signal).

- `f2`: fraction of selected genes (cardinality / sparsity).

## See also

[`initialize_population_cpp()`](https://danymukesha.github.io/BioGA/reference/initialize_population_cpp.md),
[`selection_cpp()`](https://danymukesha.github.io/BioGA/reference/selection_cpp.md),
[`crossover_cpp()`](https://danymukesha.github.io/BioGA/reference/crossover_cpp.md),
[`mutation_cpp()`](https://danymukesha.github.io/BioGA/reference/mutation_cpp.md),
[`replacement_cpp()`](https://danymukesha.github.io/BioGA/reference/replacement_cpp.md),
[`bioga_main_cpp()`](https://danymukesha.github.io/BioGA/reference/bioga_main_cpp.md)

## Examples

``` r
genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
population <- BioGA::initialize_population_cpp(genomic_data, population_size = 5)
BioGA::evaluate_fitness_cpp(genomic_data, population, weights = c(1.0, 0.5))
#>          [,1]      [,2]
#> [1,] 211.7487 0.5084159
#> [2,] 171.4019 0.3824670
#> [3,] 190.9419 0.3457173
#> [4,] 212.2748 0.4887263
#> [5,] 228.8821 0.4212201
```
