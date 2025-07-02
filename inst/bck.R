## Example of scenario

Consider an example scenario of using genetic algorithm optimization to select the best combination of genes for predicting a certain trait, such as disease susceptibility.

```{r load_library, message=FALSE, warning=FALSE}
library(BioGA)
library(SummarizedExperiment)
```

```{r parameters_setting}
# define parameters
num_genes <- 1000
num_samples <- 10

# parameters for genetic algorithm
population_size <- 100
generations <- 10
mutation_rate <- 0.5
```

```{r create_example}
# generate an example of genomic data using "SummarizedExperiment"
counts <- matrix(rpois(num_genes * num_samples, lambda = 10),
                 nrow = num_genes
)
rownames(counts) <- paste0("Gene", 1:num_genes)
colnames(counts) <- paste0("Sample", 1:num_samples)

# create "SummarizedExperiment" object
se <-
  SummarizedExperiment::SummarizedExperiment(assays = list(counts = counts))

# convert "SummarizedExperiment" to matrix for compatibility with `BioGA` 
genomic_data <- assay(se)
```

In the above example, `counts` is a matrix that represents the counts of gene expression levels across different samples. Each row corresponds to a gene, and each column corresponds to a sample. We used the `SummarizedExperiment` class to store this data, which is common Bioconductor class for representing rectangular feature x sample data, such as RNAseq count matrices or microarray data.

```{r}
head(genomic_data)
```

## Initialization

```{r}
# (select the number of canditate you wish `population`)
population <- 
  BioGA::initialize_population_cpp(genomic_data,
                                   population_size = 5
  )
```

As we already mentioned, the population represents a set of candidate combinations of genes that could be predictive of the trait. The population undergoes then genetic algorithm operations such as selection, crossover, mutation, and replacement to evolve towards individuals with higher predictive power for the trait.

## GA Optimization

```{r GA_optimization, eval=FALSE}
fitness_history <- list()

start_time <- Sys.time()
weights <- c(1.0, 0.5) 

generation <- 0
while (TRUE) {
  generation <- generation + 1
  
  # evaluate fitness
  fitness <- BioGA::evaluate_fitness_cpp(genomic_data, population, weights)
  fitness_history[[generation]] <- fitness
  
  # check termination condition
  if (generation == generations) { # defined number of generations
    break
  }
  
  # selection
  selected_parents <- BioGA::selection_cpp(population,
                                           fitness,
                                           num_parents = 2
  )
  
  # crossover and Mutation
  offspring <- BioGA::crossover_cpp(selected_parents, offspring_size = 4)
  # (no mutation in this example)
  mutated_offspring <- BioGA::mutation_cpp(offspring, 
                                           mutation_rate = mutation_rate,
                                           iteration = 1, max_iterations = 100)
  
  
  offspring_fitness <- BioGA::evaluate_fitness_cpp(genomic_data, offspring, c(1.0, 0.5))
  # replacement
  population <- BioGA::replacement_cpp(population, mutated_offspring, 
                                       fitness, offspring_fitness, num_to_replace = 1
  )
  
  # calculate time progress
  elapsed_time <- difftime(Sys.time(), start_time, units = "secs")
  
  cat(
    "\rGeneration:", generation, "- Elapsed Time:",
    format(elapsed_time, units = "secs\n"), "     "
  ) 
}
```

## Fitness calculation

The fitness calculation described in the provided code calculates a measure of dissimilarity between the gene expression profiles of individuals in the population and the genomic data. This measure of dissimilarity, or "fitness", quantifies how well the gene expression profile of an individual matches the genomic data.

Mathematically, the fitness calculation can be represented as follows:
  
  Let:
  
  -   $g_{ijk}$ be the gene expression level of gene $j$ in individual $i$ and sample $k$ from the genomic data.

-   $p_{ij}$ be the gene expression level of gene $j$ in individual $i$ from the population.

-   $N$ be the number of individuals in the population.

-   $G$ be the number of genes.

-   $S$ be the number of samples.

Then, the fitness $F_i$ for individual $i$ in the population can be calculated as the sum of squared differences between the gene expression levels of individual $i$ and the corresponding gene expression levels in the genomic data, across all genes and samples: $$
  F_i = \sum_{j=1}^{G} \sum_{k=1}^{S} (g_{ijk} - p_{ij})^2 
$$
  
  This fitness calculation aims to minimize the overall dissimilarity between the gene expression profiles of individuals in the population and the genomic data. Individuals with lower fitness scores are considered to have gene expression profiles that are more similar to the genomic data and are therefore more likely to be selected for further optimization in the genetic algorithm.

```{r plot_fitness, eval=FALSE}
# Plot fitness change over generations
BioGA::plot_fitness_history(fitness_history)
```

This showcases the integration of genetic algorithms with genomic data analysis and highlights the potential of genetic algorithms for feature selection in genomics.

Here's how `r Biocpkg("BioGA")` could work in the context of high throughput genomic data analysis:

1.  **Problem Definition**: `r Biocpkg("BioGA")` starts with a clear definition of the problem to be solved. This could include tasks such as identifying genetic markers associated with a particular disease, optimizing gene expression patterns, or clustering genomic data to identify patterns or groupings.

2.  **Representation**: Genomic data would need to be appropriately represented for use within the genetic algorithm framework. This might involve encoding the data in a suitable format, such as binary strings representing genes or chromosomes.

3.  **Fitness Evaluation**: `r Biocpkg("BioGA")` would define a fitness function that evaluates how well a particular solution performs with respect to the problem being addressed. In the context of genomic data analysis, this could involve measures such as classification accuracy, correlation with clinical outcomes, or fitness to a particular model.

4.  **Initialization**: The algorithm would initialize a population of candidate solutions, typically randomly or using some heuristic method. Each solution in the population represents a potential solution to the problem at hand.

5.  **Genetic Operations**: `r Biocpkg("BioGA")` would apply genetic operators such as selection, crossover, and mutation to evolve the population over successive generations. Selection identifies individuals with higher fitness to serve as parents for the next generation. Crossover combines genetic material from two parent solutions to produce offspring. Mutation introduces random changes to the offspring to maintain genetic diversity.

6.  **Termination Criteria**: The algorithm would continue iterating through generations until a termination criterion is met. This could be a maximum number of generations, reaching a satisfactory solution, or convergence of the population.

7.  **Result Analysis**: Once the algorithm terminates, `r Biocpkg("BioGA")` would analyze the final population to identify the best solution(s) found. This could involve further validation or interpretation of the results in the context of the original problem.

Other applications of `r Biocpkg("BioGA")` in genomic data analysis could include genome-wide association studies (GWAS), gene expression analysis, pathway analysis, and predictive modeling for personalized medicine, among others. By leveraging genetic algorithms, `r Biocpkg("BioGA")` offers a powerful approach to exploring complex genomic datasets and identifying meaningful patterns and associations.