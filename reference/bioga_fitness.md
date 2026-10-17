# Evaluate BioGA fitness

User-facing wrapper around
[`evaluate_fitness_cpp()`](https://danymukesha.github.io/BioGA/reference/evaluate_fitness_cpp.md).

## Usage

``` r
bioga_fitness(
  genomic_data,
  population,
  weights = c(1, 0.5),
  labels = NULL,
  representation = c("real", "binary")
)
```

## Arguments

- genomic_data:

  Numeric matrix of genomic data (rows: genes, columns: samples).

- population:

  Numeric matrix of individuals (rows: individuals, columns: genes).

- weights:

  Numeric vector of objective weights (length 1 or 2).

- labels:

  Optional class labels, one per sample (integer, numeric, factor,
  character or logical). Required when `representation = "binary"`.

- representation:

  Either `"real"` or `"binary"`.

## Value

Numeric matrix of fitness scores (rows: individuals, columns:
objectives).

## See also

[`bioga()`](https://danymukesha.github.io/BioGA/reference/bioga.md),
[`evaluate_fitness_cpp()`](https://danymukesha.github.io/BioGA/reference/evaluate_fitness_cpp.md)

## Examples

``` r
genomic_data <- matrix(rnorm(200), nrow = 20, ncol = 10)
population <- initialize_population_cpp(genomic_data, population_size = 5)
bioga_fitness(genomic_data, population)
#>          [,1]      [,2]
#> [1,] 304.2346 0.3832998
#> [2,] 323.1340 0.3914986
#> [3,] 220.9086 0.2619454
#> [4,] 287.5368 0.3194004
#> [5,] 347.6736 0.4414341
```
