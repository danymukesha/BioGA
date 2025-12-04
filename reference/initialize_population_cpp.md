# Function to initialize population with optional gene clustering

Function to initialize population with optional gene clustering

## Usage

``` r
initialize_population_cpp(
  genomic_data,
  population_size,
  seed = NULL,
  clusters = NULL
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

  Optional vector of gene cluster assignments.

## Value

Numeric matrix of initialized population.

## Examples

``` r
genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
BioGA::initialize_population_cpp(genomic_data, population_size = 5, 
seed = 123)
#>            [,1]       [,2]      [,3]        [,4]       [,5]        [,6]
#> [1,]  0.6644805 -0.3501541 0.2189152  0.07439935  0.3548044 -0.04733266
#> [2,] -0.3383459  0.5950387 0.3774680  0.27211081  2.5017651  0.65099328
#> [3,]  1.4264917  0.9967010 0.3774680 -0.06042597 -0.9732140  0.79157269
#> [4,] -0.3383459 -0.7006975 0.3774680 -0.08981425 -0.7795961  0.07504484
#> [5,] -0.6558381  0.5950387 0.2189152 -0.80241961  2.5017651  0.79416747
#>            [,7]       [,8]       [,9]      [,10]
#> [1,] -0.9944350  1.2929929  1.0936361 -1.0047212
#> [2,] -0.3395090  0.9344820 -0.2306096  0.6982138
#> [3,] -0.9944350 -0.4153374  0.2909621 -0.3092672
#> [4,] -0.4927276 -0.5939900 -0.2306096 -0.0550220
#> [5,] -0.3395090 -1.3970639  0.2909621 -1.1124601
```
