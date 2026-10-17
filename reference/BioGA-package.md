# BioGA: Bioinformatics Genetic Algorithm

BioGA provides a genetic algorithm framework for analyzing and
optimizing high throughput genomic data. All core functions are
implemented in C++ for performance via Rcpp and RcppParallel, with an
easy-to-use R interface.

## Details

Two representations are supported. In `"real"` mode individuals are
continuous expression profiles; in `"binary"` mode individuals are gene
selection indicators and class `labels` are required. The workflow is:

1.  **Initialization:**
    [`initialize_population_cpp()`](https://danymukesha.github.io/BioGA/reference/initialize_population_cpp.md)
    — Creates an initial population, optionally guided by gene clusters.

2.  **Fitness evaluation:**
    [`evaluate_fitness_cpp()`](https://danymukesha.github.io/BioGA/reference/evaluate_fitness_cpp.md)
    — Evaluates two minimised objectives (data fidelity and parsimony in
    `"real"` mode; predictive signal and sparsity in `"binary"` mode),
    computed in parallel.

3.  **Selection:**
    [`selection_cpp()`](https://danymukesha.github.io/BioGA/reference/selection_cpp.md)
    — Selects parents using NSGA-II rank and crowding-distance
    tournament selection.

4.  **Crossover:**
    [`crossover_cpp()`](https://danymukesha.github.io/BioGA/reference/crossover_cpp.md)
    — Generates offspring via simulated binary crossover (SBX) in
    `"real"` mode or uniform crossover in `"binary"` mode.

5.  **Mutation:**
    [`mutation_cpp()`](https://danymukesha.github.io/BioGA/reference/mutation_cpp.md)
    — Applies adaptive Gaussian (or bit-flip) mutation with optional
    degree-based network damping.

6.  **Replacement:**
    [`replacement_cpp()`](https://danymukesha.github.io/BioGA/reference/replacement_cpp.md)
    — Admits the best offspring by NSGA-II environmental selection.

7.  **Main loop:**
    [`bioga_main_cpp()`](https://danymukesha.github.io/BioGA/reference/bioga_main_cpp.md)
    — Orchestrates the full pipeline;
    [`bioga()`](https://danymukesha.github.io/BioGA/reference/bioga.md)
    is the user-facing entry point.

## Citation

If you use BioGA in your research, please cite:


    Mukesha D (2024). "Scalable and Interpretable Genetic Algorithm Framework
    for High-Throughput Genomic Data." Zenodo. \doi{10.5281/zenodo.15801072}

## See also

Useful links:

- <https://danymukesha.github.io/BioGA/>

- <https://bioconductor.org/packages/BioGA/>

- Report bugs at <https://github.com/danymukesha/BioGA/issues>

## Author

**Maintainer:** Dany Mukesha <danymukesha@gmail.com>
