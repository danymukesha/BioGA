% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{crossover_cpp}
\alias{crossover_cpp}
\title{Function to perform simulated binary crossover (SBX)}
\usage{
crossover_cpp(
  selected_parents,
  offspring_size,
  crossover_rate = 0.9,
  eta_c = 20
)
}
\arguments{
\item{selected_parents}{Numeric matrix of selected individuals.}

\item{offspring_size}{Number of offspring to generate.}

\item{crossover_rate}{Probability of crossover (default: 0.9).}

\item{eta_c}{SBX distribution index (default: 20.0).}
}
\value{
Numeric matrix of offspring.
}
\description{
Function to perform simulated binary crossover (SBX)
}
\examples{
genomic_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
population <- BioGA::initialize_population_cpp(genomic_data, 
population_size = 5)
fitness <- BioGA::evaluate_fitness_cpp(genomic_data, population, 
c(1.0, 0.5))
selected_parents <- BioGA::selection_cpp(population, fitness, 
num_parents = 2)
BioGA::crossover_cpp(selected_parents, offspring_size = 2)
}
