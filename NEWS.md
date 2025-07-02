# BioGA 0.99.17

Added:
- Demo vignette for BioGA.
- Mathematical background vignette for BioGA.

Updated:
- BioGA vignettes with new examples and explanations.
- `selection_cpp` function for improved performance and bug fixes.

# BioGA 0.99.16

*UPDATE: test_crossover_cpp*

The test verifies dimensions, value ranges, and new parameters.

Changes:
- Added tests for new parameters (`crossover_rate`, `eta_c`).
- Verified that offspring values lie within parent ranges (due to SBX).
- Included edge case test for single parent and zero crossover rate.
- Simplified error/warning checks.

*UPDATE: test_selection_cpp*
The test checks dimensions, parent selection, and compatibility with multi-objective fitness.

Changes:
- Updated to handle multi-objective fitness matrix.
- Verified that selected parents are exact copies of population rows.
- Added edge case test for single individual and single objective.
- Removed redundant fitness evaluation test (already covered in `evaluate_fitness_cpp`).

*UPDATE: test_mutation_cpp*
The test checks dimensions, mutation effects, and new parameters.

Changes:
- Added tests for new parameters (`iteration`, `max_iterations`, `network`).
- Verified that zero mutation rate preserves offspring.
- Tested network constraints with extreme cases (zero and full constraints).
- Simplified error/warning checks.

*UPDATE/ADDED: test_replacement_cpp*

The test checks dimensions, elite preservation, and diversity.

Changes:
- Added test for new function signature (requires fitness matrices).
- Verified elite preservation based on first objective.
- Included edge case test for zero replacements.
- Added checks for dimensions and error-free execution.

*UPDATE: test_initialize_population_cpp*

The test checks dimensions, value validity, and new parameters.

Changes:
- Added test for `seed` parameter to ensure reproducibility.
- Included test for `clusters` parameter to verify correct initialization.
- Maintained tests for dimensions and value validity.
- Added edge case test for minimal input.

*UPDATE: test_bioga_main_cpp* {PROVISIONAL}

New test file is created to verify the main GA loop.

Feature:
- Tests the output structure (list with population and fitness).
- Verifies dimensions of population and fitness matrices.
- Checks optional parameters (`clusters`, `network`).
- Includes edge case test for minimal input.

*Test for Package (test-BioGA-package.R)*

Created a new version to test the overall package integrity and ensure all functions are exported correctly.

Feature:
- Tests package loading and function exports.
- Verifies a complete GA workflow without errors.
- Ensures compatibility with the optimized functions.

# BioGA 0.99.15

*UPDATE: test_evaluate_fitness_cpp* 
The test checks dimensions, values, and error-free execution.

Changes:

* Updated to test the new matrix output (multiple objectives).
* Added checks for non-negative fitness values (appropriate for expression difference and sparsity).
* Included edge case test for minimal input.
* Removed redundant `tryCatch` and simplified error/warning checks.

# BioGA 0.99.14

*Commit -m "New Feature - Main GA Loop"*

adding a main GA loop function that integrates all components 
and supports multi-objective optimization and parallelization.

Features:

* Integrates all optimized components into a cohesive GA loop.
* Supports multi-objective optimization, gene networks, and clustering.
* Returns both final population and fitness scores for analysis.

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
