# BioGA 0.99.18

## Major changes
- Rewrote core algorithm as a correct NSGA-II implementation: strict
  domination, full non-dominated front construction, crowding distance, and
  NSGA-II environmental selection for replacement (new `src/nsga2.h`).
- Added a hybrid representation via `representation = c("real", "binary")`.
  Real mode minimises squared error (`f1`) and mean absolute value (`f2`);
  binary mode minimises the negated mean ANOVA F-statistic of selected genes
  (`f1`) and the fraction of selected genes (`f2`) using class `labels`.
- Added idiomatic R wrappers `bioga()` and `bioga_fitness()`.
- `crossover_cpp` now generates both complementary SBX children and handles a
  single parent without hanging.
- `mutation_cpp` uses degree-based network damping (bounded, non-negative) for
  the real representation and bit-flip mutation for the binary representation.
- `replacement_cpp` now performs NSGA-II environmental selection;
  `num_to_replace = 0` returns the population unchanged.
- `bioga_main_cpp` returns `population`, `fitness`, `history`, `best`,
  `best_fitness`, `pareto_population`, and `pareto_fitness`, and accepts a
  `labels` argument for supervised objectives.

## Documentation
- Corrected scientific statements in the vignettes (removed unverified claims
  about tools, benchmarks, databases, and convergence guarantees; aligned
  formulae with the implementation).
- Lowered `Depends` to `R (>= 4.3)`.

## Testing
- Rewrote the test suite for the new API and behaviours.

# BioGA 0.99.17

## Added
- Demo vignette for BioGA.
- Mathematical background vignette for BioGA.

## Updated
- Vignettes with new examples and explanations.
- `selection_cpp` for improved performance and bug fixes.

# BioGA 0.99.16

## Testing improvements
- Tests for `crossover_cpp`: new parameters (`crossover_rate`, `eta_c`), SBX offspring range validation, edge cases.
- Tests for `selection_cpp`: multi-objective fitness matrix support, parent identity verification, single-individual edge case.
- Tests for `mutation_cpp`: new parameters (`iteration`, `max_iterations`, `network`), zero-rate preservation, network constraint extremes.
- Tests for `replacement_cpp`: new signature with fitness matrices, elite preservation, zero-replacement edge case.
- Tests for `initialize_population_cpp`: `seed` reproducibility, `clusters` parameter, minimal input edge case.
- Tests for `bioga_main_cpp`: output structure, dimension verification, optional parameters, minimal input.
- Package integration test: export verification, full workflow execution.

# BioGA 0.99.15

## Testing improvements
- Updated `evaluate_fitness_cpp` tests for multi-objective matrix output.
- Added non-negative fitness checks.
- Added minimal input edge case.
- Removed redundant `tryCatch`.

# BioGA 0.99.14

## New feature
- Main GA loop (`bioga_main_cpp`) integrating all components with multi-objective optimization and parallelization.

# BioGA 0.99.13

## Improvements
- Biologically informed initialization using gene clusters.
- Reproducible random seed handling.

# BioGA 0.99.12

## Maintenance
- Code cleanup and indentation in `replacement.cpp`.

# BioGA 0.99.11

## Improvements
- Elitism to preserve best individuals.
- Diversity-based replacement to avoid premature convergence.
- Multi-objective fitness matrix support.

# BioGA 0.99.10

## Improvements
- Adaptive mutation rate based on iteration progress.
- Optional gene network constraints for biologically relevant mutations.

# BioGA 0.99.9

## Improvements
- NSGA-II non-dominated sorting for multi-objective selection.
- Tournament selection for diversity.

# BioGA 0.99.8

## Improvements
- Simulated binary crossover (SBX) for better exploration.
- Adaptive crossover rate based on population diversity.
- Configurable crossover rate and distribution index.

# BioGA 0.99.7

## Improvements
- Multi-objective fitness evaluation (expression difference and sparsity).
- Parallelized computation via RcppParallel.
- Flexible objective weighting via `weights` parameter.

# BioGA 0.99.6

## Maintenance
- Added authors and date to vignettes.
- Set minimum R version to 4.4.

# BioGA 0.99.4

## Maintenance
- Updates requested by Bioconductor peer review.

# BioGA 0.99.2

## Bug fixes
- Fixed Git tracking of `.Rproj` system file.

# BioGA 0.99.0

## Initial release
- Added `NEWS.md` to track changes.
