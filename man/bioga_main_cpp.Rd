% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{bioga_main_cpp}
\alias{bioga_main_cpp}
\title{Main genetic algorithm loop for genomic data optimization}
\usage{
bioga_main_cpp(
  genomic_data,
  population_size,
  num_generations,
  crossover_rate,
  eta_c,
  mutation_rate,
  num_parents,
  num_offspring,
  num_to_replace,
  weights,
  seed = NULL,
  clusters = NULL,
  network = NULL
)
}
\arguments{
\item{genomic_data}{Numeric matrix of genomic data (rows: genes,
columns: samples).}

\item{population_size}{Number of individuals in the population.}

\item{num_generations}{Number of generations to run.}

\item{crossover_rate}{Probability of crossover.}

\item{eta_c}{SBX distribution index (default: 20.0).}

\item{mutation_rate}{Base probability of mutation.}

\item{num_parents}{Number of parents to select per generation.}

\item{num_offspring}{Number of offspring to generate per generation.}

\item{num_to_replace}{Number of individuals to replace per generation.}

\item{weights}{Numeric vector of weights for multi-objective fitness.}

\item{seed}{Optional random seed for reproducibility.}

\item{clusters}{Optional vector of gene cluster assignments.}

\item{network}{Optional matrix of gene network constraints.}
}
\value{
List containing final population and fitness scores.
}
\description{
Main genetic algorithm loop for genomic data optimization
}
\examples{
genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
result <- BioGA::bioga_main_cpp(genomic_data,
population_size = 50, num_generations = 10,
        crossover_rate = 0.9, eta_c = 20.0, mutation_rate = 0.1,
        num_parents = 20, num_offspring = 20, num_to_replace = 10,
        weights = c(1.0, 0.5), seed = 123)
}
