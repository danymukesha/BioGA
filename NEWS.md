# BioGA 0.99.13

*Commit -m "Optimized initialize_population.cpp"*

Improvements:

* Add option for biologically informed initialization using gene clusters.
* Improve random seed handling for reproducibility.

Changes:

* Added optional random seed for reproducibility.
* Included support for gene clustering to initialize biologically relevant populations.
* Improved code readability and documentation.

# BioGA 0.99.12

*Commit -m "Optimized replacement.cpp"*

indentation adjustments

# BioGA 0.99.11

*Commit -m "Optimized replacement.cpp"*

Improvements:

* Implement elitism to preserve best individuals.
* Add diversity-based replacement to avoid premature convergence.

Changes:

* Added elitism to preserve the best individual.
* Included diversity-based replacement to maintain population diversity.
* Updated to use multi-objective fitness matrix.

# BioGA 0.99.10

*Commit -m "Optimized mutation.cpp"*

Improvements:

* Add adaptive mutation rate based on iteration or fitness stagnation.
* Incorporate gene network constraints (placeholder for user-provided network).

Changes:

* Added adaptive mutation rate based on iteration progress.
* Included optional gene network constraints to ensure biologically relevant mutations.
* Maintained compatibility with existing functionality.

# BioGA 0.99.9

*Commit -m "Optimized selection.cpp"*

Improvements:

* Implement NSGA-II non-dominated sorting for multi-objective optimization.
* Add tournament selection for better diversity.

Changes:

* Implemented NSGA-II non-dominated sorting for multi-objective selection.
* Added tournament selection to maintain diversity.
* Updated to handle multi-objective fitness matrix.

# BioGA 0.99.8

*Commit -m "Optimized crossover.cpp"*
"
Improvements:

* Implement simulated binary crossover (SBX) for better exploration.
* Add adaptive crossover rate based on population diversity.

Changes:

* Replaced simple averaging with SBX crossover for better exploration.
* Added adaptive crossover rate based on population diversity.
* Included parameters for crossover rate and distribution index.

# BioGA 0.99.7

*Commit -m "Optimized evaluate_fitness.cpp"*

Improvements:

* Add multi-objective fitness evaluation (e.g., minimize expression difference and maximize gene sparsity).
* Use vectorized operations for faster computation.
* Add parallelization with RcppParallel for large datasets.

Changes:

* Added multi-objective support (expression difference and sparsity).
* Used `RcppParallel` for parallel computation.
* Returned a matrix of fitness scores for each objective.
* Added weights parameter for flexible objective prioritization.

# BioGA 0.99.6

* add the authors and date in the vignettes?
* no R version less than 4.4

# BioGA 0.99.5

# BioGA 0.99.4

* Fixed and added updates requested by Bioconductor Peer review

# BioGA 0.99.3

# BioGA 0.99.2

* Fixed ERROR: System files 'BioGA.Rproj' found that should not be Git tracked.

# BioGA 0.99.1

# BioGA 0.99.0

* Added a `NEWS.md` file to track changes to the package.
