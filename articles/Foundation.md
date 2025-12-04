# The Foundation of BioGA

## Introduction

The BioGA utilizes multi-objective Genetic Algorithm (GA) specifically
tailored for the optimization and analysis of genomic data. A detailed
mathematical and theoretical analysis of the algorithms implemented
within this package is described below; covering all major components of
the evolutionary process with formal mathematical notation, proofs, and
a discussion of key properties.

## The Genetic Algorithm Framework

The BioGA package employs a standard generational Genetic Algorithm
structure, modified for multi-objective optimization. This framework
involves the iterative evolutionary cycle of initialization, evaluation,
selection, reproduction (crossover and mutation), and replacement.

### Mathematical Representation

Let us define the core elements of the search space and the algorithm’s
state:

- $G = (V,E)$ is the gene network, where
  $V = \{ g_{1},g_{2},\ldots,g_{n}\}$ are $n$ genes.
- $X \in {\mathbb{R}}^{n \times m}$ is the genomic data matrix,
  consisting of $n$ genes across $m$ samples.
- $P_{t} \in {\mathbb{R}}^{p \times n}$ is the population matrix at
  generation $t$, containing $p$ individuals (rows) and $n$ gene values
  (columns).

The evolution from one generation to the next,
$\left. P_{t}\rightarrow P_{t + 1} \right.$, is described by the
following functional relationship:

$$P_{t + 1} = R\left( M\left( C\left( S\left( P_{t},f\left( P_{t} \right) \right),X \right) \right),P_{t},f\left( P_{t} \right) \right)$$

Where the functions represent the core operators:

- $f$: **Fitness Evaluation** function.
- $S$: **Selection** operator (NSGA-II inspired).
- $C$: **Crossover** operator (SBX).
- $M$: **Mutation** operator (Adaptive).
- $R$: **Replacement** operator (Elitism + Diversity Preservation).

## Population Initialization

The initialization phase creates the starting population, $P_{0}$, by
drawing random samples directly from the measured genomic data $X$.

### Mathematical Formulation

Given the genomic data $X \in {\mathbb{R}}^{n \times m}$ and a
population size $p$, the value of gene $j$ in the $i$-th individual of
the initial population is set as:

$$P_{0}\lbrack i,j\rbrack = X\lbrack j,k\rbrack\quad\text{where}\quad k \sim \text{Uniform}\{ 1,\ldots,m\}$$
This means that for each element $(i,j)$ in $P_{0}$, a sample $k$ is
chosen uniformly at random from the $m$ available samples, and the
corresponding gene expression value $X\lbrack j,k\rbrack$ is used.

**Initialization with Clustering:**

If a clustering breakdown is provided (e.g., $C$ is the set of
clusters), the initialization respects the structure:
$$P_{0}\lbrack i,j\rbrack = X\lbrack j,k\rbrack\quad\forall j \in c,\quad k \sim \text{Uniform}\{ 1,\ldots,m\}$$
for all genes $j$ belonging to the same cluster $c$.

### Statistical Properties

1.  **Distribution Preservation:** This method ensures that the initial
    population individuals maintain the original statistical
    distribution of expression values for each gene $j$.
2.  **Structural Integrity:** If clustering is used, the method
    preserves the pre-defined cluster structure within the initial
    population.
3.  **Expected Value:** The expected value of any element
    $P_{0}\lbrack i,j\rbrack$ over the initialization process is the
    mean of the observed expression for gene $j$:
    $${\mathbb{E}}\left\lbrack P_{0}\lbrack i,j\rbrack \right\rbrack = \mu_{j} = \frac{1}{m}\sum\limits_{k = 1}^{m}X\lbrack j,k\rbrack$$

## Fitness Evaluation

The `BioGA` package utilizes a multi-objective approach, simultaneously
optimizing two conflicting goals: data fidelity (Expression Difference)
and model parsimony (Sparsity).

### Objective 1: Expression Difference ($f_{1}$)

This function measures the squared Euclidean distance between the
individual $i$ and all $m$ observed data samples $X_{jk}$ for all $n$
genes. It quantifies how accurately the individual reflects the overall
expression patterns.

$$f_{1}(i) = \sum\limits_{j = 1}^{n}\sum\limits_{k = 1}^{m}\left( X_{jk} - P_{ij} \right)^{2}$$

#### Properties of $f_{1}$

1.  **Convexity:** $f_{1}$ is a convex function, which is a desirable
    property for optimization as it guarantees that any local minimum is
    also a global minimum in the absence of other constraints. The
    minimum occurs when $P_{ij}$ is the mean expression $\mu_{j}$.
2.  **Gradient (for theoretical analysis):** The partial derivative with
    respect to the gene value $P_{ij}$ is:
    $$\nabla_{P_{ij}}f_{1} = - 2\sum\limits_{k = 1}^{m}\left( X_{jk} - P_{ij} \right)$$

### Objective 2: Sparsity ($f_{2}$)

This objective promotes parsimonious solutions by penalizing non-zero
gene values. It measures the fraction of non-zero genes within an
individual $i$.

$$f_{2}(i) = \frac{\sum\limits_{j = 1}^{n}I\left( \left| P_{ij} \right| > \epsilon \right)}{n}$$

Here, $I$ is the indicator function, which equals 1 if the condition is
met, and $\epsilon$ is a small constant (typically $10^{- 6}$) used to
define “non-zero.”

#### Properties of $f_{2}$

1.  **Mathematical Form:** $f_{2}$ is inherently **Non-convex** and
    **Non-differentiable** due to the indicator function $I$. This makes
    standard gradient-based optimization methods unsuitable, justifying
    the use of a heuristic approach like the GA.
2.  **Goal:** It actively encourages sparse solutions, where many gene
    values are effectively zero.
3.  **Range:** The range of the function is restricted to
    $\lbrack 0,1\rbrack$, where $0$ represents a fully sparse individual
    ($P_{ij} = 0$ for all $j$) and $1$ represents a non-sparse
    individual.

### Combined Fitness

The overall fitness value $F(i)$ for an individual $i$ is a weighted
combination of the two objectives:

$$F(i) = w_{1}f_{1}(i) + w_{2}f_{2}(i)$$

Where $w_{1}$ and $w_{2}$ are user-specified weights that dictate the
trade-off between minimizing expression error ($f_{1}$) and maximizing
sparsity ($f_{2}$).

## Selection (NSGA-II Inspired)

The selection mechanism is based on the principles of the Non-dominated
Sorting Genetic Algorithm II (NSGA-II), ensuring that evolutionary
pressure drives the population toward the Pareto front.

### Domination Criteria

Individual $i$ dominates individual $j$ (denoted $i \succ j$) if and
only if:

1.  Individual $i$ is no worse than $j$ on all objectives, and
2.  Individual $i$ is strictly better than $j$ on at least one
    objective.

Formally:
$$i \succ j\quad\text{iff}\quad\forall k \in \{ 1,2\}:f_{k}(i) \leq f_{k}(j)\quad\text{and}\quad\exists k \in \{ 1,2\}:f_{k}(i) < f_{k}(j)$$

#### Proof of Partial Order

The two-objective domination relation $( \succ )$ defines a **partial
order** on the solution set because it satisfies the following
properties:

1.  **Irreflexive:** No individual dominates itself. If $i = j$, the
    second condition ($\exists k:f_{k}(i) < f_{k}(j)$) is false.
2.  **Antisymmetric Property:** If $i \succ j$, then $j$ cannot dominate
    $i$. If $i$ is strictly better on at least one objective, $j$ must
    be strictly worse ($\forall k:f_{k}(j) \nleq f_{k}(i)$).
3.  **Transitive:** If $i \succ j$ and $j \succ k$, then $i \succ k$.
    This follows directly from the transitivity of the combined $\leq$
    and $<$ relations.

### Front Construction

The population is sorted into non-dominated fronts:

1.  **Domination Metrics:** For every individual, we compute its
    **domination count** (number of individuals that dominate it) and
    its **dominated set** (set of individuals it dominates).
2.  **First Front ($\mathcal{F}_{1}$):** This front consists of all
    individuals whose domination count is $0$. These are the current
    non-dominated solutions.
3.  **Subsequent Fronts:** Individuals in $\mathcal{F}_{t}$ are
    temporarily removed. The domination counts of the individuals they
    dominated are reduced by one. The next front, $\mathcal{F}_{t + 1}$,
    consists of all remaining individuals whose updated domination count
    is $0$. This process continues until the entire population is
    sorted.

#### Theorem on Computational Complexity

The non-dominated sorting and front construction algorithm terminates in
quadratic time complexity relative to the population size.

**Theorem:** The front construction algorithm terminates in
$O\left( p^{2}o \right)$ time, where $p$ is the population size and $o$
is the number of objectives.

**Proof Sketch:**

- The pairwise domination check between any two individuals is an $O(o)$
  operation.
- The computation of the domination count and dominated set for all
  individuals requires $O\left( p^{2}o \right)$ time, as it involves
  $\left( \frac{p}{2} \right)$ comparisons.
- Once these metrics are computed, the iterative process of assigning
  individuals to fronts requires $O(p)$ time per front, leading to an
  overall complexity dominated by the $O\left( p^{2}o \right)$ initial
  comparison phase.

## Crossover (Simulated Binary Crossover - SBX)

The Simulated Binary Crossover (SBX) is a prominent operator for
real-coded GAs, designed to mimic the convergence properties of binary
crossover while operating directly on continuous variables.

Given two parent individuals, $x,y \in {\mathbb{R}}^{n}$, the SBX
operator generates a new offspring vector $z$. For each gene $j$:

With crossover probability $p_{c}$, the offspring gene $z_{j}$ is
determined: 1. A random number $u$ is drawn:
$u \sim \text{Uniform}(0,1)$. 2. The spreading factor $\beta$ is
computed based on the distribution index $\eta$: $$\beta = \begin{cases}
(2u)^{1/{(\eta + 1)}} & {{\text{if}\mspace{6mu}}u \leq 0.5} \\
\left( \frac{1}{2(1 - u)} \right)^{1/{(\eta + 1)}} & \text{otherwise}
\end{cases}$$ 3. The gene value $z_{j}$ is calculated:
$$z_{j} = 0.5\left\lbrack \left( x_{j} + y_{j} \right) - \beta\left| y_{j} - x_{j} \right| \right\rbrack$$

Otherwise (with probability $1 - p_{c}$):
$$z_{j} = x_{j}\quad\text{or}\quad z_{j} = y_{j}$$

### Properties of SBX

1.  **Mean Preservation:** The expected value of the offspring gene is
    average of the parents’ genes:
    ${\mathbb{E}}\left\lbrack z_{j} \right\rbrack = \frac{x_{j} + y_{j}}{2}$.
2.  **Distribution Control:** The distribution index $\eta$ controls the
    proximity of the offspring to the parents.
3.  **Low $\eta$ (e.g., $\left. \eta\rightarrow 0 \right.$):** Produces
    a wide spread and approaches uniform crossover, generating many
    extreme values.
4.  **High $\eta$ (e.g., $\left. \eta\rightarrow\infty \right.$):**
    Clusters offspring tightly around the parents, approaching no
    crossover.

## Mutation

The mutation operation introduces random small changes to the offspring
for exploration. The `BioGA` uses an adaptive mutation rate and
incorporates network constraints.

For each gene $j$, with adaptive probability $p_{m}(t)$:
$$p_{m}(t) = p_{0}\left( 1 + 0.5\frac{t}{T} \right)$$ where $p_{0}$ is
the initial rate, $t$ is the current generation, and $T$ is the total
number of generations. This rate increases linearly over time to enhance
exploration in later stages.

The change in the gene value, $z_{j}$, is applied as follows:

1.  A Gaussian perturbation is generated:
    $\Delta_{j} \sim N\left( 0,\sigma^{2} \right)$.

2.  **If a Gene Network $N$ is provided ($N_{jk}$ being the adjacency
    matrix):** The mutation is penalized by the network connectivity:
    $$\left. z_{j}\leftarrow z_{j} + \Delta_{j}\left( 1 - \sum\limits_{k}N_{jk}z_{k} \right) \right.$$

3.  **If no Network is provided:** Standard Gaussian mutation is
    applied: $$\left. z_{j}\leftarrow z_{j} + \Delta_{j} \right.$$

### Properties of Adaptive Network Mutation

1.  **Adaptive Rate:** The mutation rate monotonically increases with
    generation $t$, preventing premature convergence and maintaining
    diversity late in the run.
2.  **Network Constraint Impact:** The term
    $\left( 1 - \sum_{k}N_{jk}z_{k} \right)$ acts as a coefficient.
    Highly connected genes (large $\sum_{k}N_{jk}$) experience a reduced
    magnitude of mutation, preserving known network structures.
3.  **Expected Change:** Since
    ${\mathbb{E}}\left\lbrack \Delta_{j} \right\rbrack = 0$, the
    expected change remains zero:
    ${\mathbb{E}}\left\lbrack \Delta z_{j} \right\rbrack = 0$.
4.  **Variance Control:** For the network-constrained case, the variance
    of the change is:
    $$\text{Var}\left( \Delta z_{j} \right) = \sigma^{2}\left( 1 - \sum\limits_{k}N_{jk}z_{k} \right)^{2}$$
    This mathematical structure formalizes how network relationships
    dynamically control evolutionary exploration.

## Replacement

The replacement operator $R$ constructs the next population $P_{t + 1}$
from the current population $P_{t}$ and the offspring $O_{t}$. It
employs an Elitism strategy combined with a metric for diversity
preservation.

1.  **Elitism:** The best individual, $x^{*}$, defined as the solution
    minimizing the primary objective $f_{1}$, is always retained in the
    next generation:
    $$x^{*} = \text{argmin}_{x \in P_{t} \cup O_{t}}f_{1}(x)$$
2.  **Diversity-Preserving Replacement:** For the remaining $p - 1$
    slots in $P_{t + 1}$, a non-elitist replacement strategy is
    performed iteratively:
    - Randomly select an individual $x$ from the current population
      $P_{t}$.
    - Select an offspring $y$ from the newly generated offspring
      $O_{t}$.
    - Replace $x$ with $y$ only if the increase in diversity is
      significant:
      $$\left. {\text{Replace}\mspace{6mu}}x\leftarrow y\quad\text{if}\quad\text{diversity}(x,y) > \epsilon \right.$$
      Where the diversity metric is the squared Euclidean distance:
      $\text{diversity}(x,y) = \parallel x - y \parallel_{2}^{2}$.

### Theorem on Population Quality

**Theorem:** This replacement strategy preserves elitism while ensuring
that the population’s diversity does not collapse into a single point.

**Proof Sketch:**

1.  **Elitism Preservation (Quality):** By explicitly retaining $x^{*}$,
    the search process ensures that the best solution found so far
    across the primary objective $f_{1}$ is never lost. This guarantees
    monotonic, non-decreasing convergence quality on $f_{1}$.
2.  **Diversity Maintenance:** Replacements are conditional on
    $\text{diversity}(x,y) > \epsilon$. Since only replacements that
    increase or maintain the distance (and hence diversity) above a
    threshold are accepted, the expected diversity of the population is
    non-decreasing over the generations, effectively counteracting
    sampling drift.

## Convergence Analysis

The convergence of the BioGA algorithm, derived from established
Multi-objective Evolutionary Algorithm (MOEA) theory, guarantees
convergence to the optimal set under general conditions.

### Assumptions

For convergence proofs in MOEAs, the following are often assumed:

1.  **Finite Search Space:** While the gene values are continuous, they
    can be discretized in theory (or bounded in practice).
2.  **Strictly Positive Mutation Probability:** $p_{m}(t) > 0$ for all
    $t$, which ensures ergodicity.
3.  **Elitism:** The best solutions are preserved across generations.

### Theorem on Convergence

**Theorem:** Under the assumptions above, the BioGA algorithm converges
in probability to the global Pareto front of the multi-objective
optimization problem.

**Proof Sketch:**

1.  **Preservation of Optimality:** The NSGA-II selection and Elitism
    replacement ensure that once a solution is found on the Pareto
    front, it or a strictly better non-dominated alternative is
    retained.
2.  **Ergodicity:** The strictly positive mutation probability
    guarantees that any point in the search space can theoretically be
    reached from any other point in a finite number of steps.
3.  **Convergence Result:** Based on general MOEA convergence theorems
    (e.g., Rudolph 1998, or similar proofs for NSGA-II), the combination
    of selection pressure toward non-dominated solutions and exploratory
    power (mutation/crossover) ensures that the population will tend
    toward the true Pareto optimal set as
    $\left. t\rightarrow\infty \right.$.

## Computational Complexity

Analyzing the complexity is crucial for understanding the algorithm’s
performance on large genomic datasets. Key variables are defined:

- $p$ = population size
- $n$ = number of genes
- $m$ = number of samples
- $o$ = number of objectives ($o = 2$ in BioGA)
- $T$ = number of generations

### Component Complexities

| GA Component        | Complexity (Per Generation)              | Dominating Variables |
|:--------------------|:-----------------------------------------|:---------------------|
| Initialization      | $O(pn)$                                  | $p,n$                |
| Fitness Evaluation  | $O(pn) \times O(m)$$\rightarrow$$O(pmn)$ | $p,m,n$              |
| Selection (Sorting) | $O\left( p^{2}o \right)$                 | $p$ (quadratic)      |
| Crossover (SBX)     | $O(pn)$                                  | $p,n$                |
| Mutation (Adaptive) | $O(pn)$                                  | $p,n$                |
| Replacement         | $O(pn)$                                  | $p,n$                |

### Total Complexity

The total complexity for $T$ generations is dominated by the Fitness
Evaluation and the Selection phases:

$$\text{Total Complexity} = O\left( T \cdot \left( pmn + p^{2}o \right) \right)$$
Factoring $p$: $$\text{Total Complexity} = O\left( Tp(mn + po) \right)$$
This result highlights that scalability is limited by the quadratic
dependency on population size $p$ (due to sorting) and the linear
dependency on the product of genes/samples $mn$ (due to fitness
evaluation).

## Mathematical Optimization Interpretation

From a purely mathematical perspective, the BioGA algorithm serves as a
robust stochastic optimization solver for the following multi-objective
problem:

$$\begin{aligned}
{\text{minimize}\mspace{6mu}} & \left( f_{1}(P),f_{2}(P) \right) \\
{\text{subject to}\mspace{6mu}} & {P \in {\mathbb{R}}^{p \times n}}
\end{aligned}$$

**Core Characteristics that Favor GA:**

1.  **Multi-Objective Nature:** Explicitly seeking a set of trade-off
    solutions (the Pareto front).
2.  **High Dimensionality:** Search space size is $p \times n$, making
    deterministic grid search intractable.
3.  **Non-Differentiable Objectives:** The sparsity objective $f_{2}$ is
    non-differentiable, making standard calculus-based methods
    impossible.
4.  **Non-Convexity:** The fitness landscape, when combined, may be
    non-convex, meaning the GA’s global search capability is essential
    to avoid local optima.

## Biological Interpretation and Relevance

The mathematical operations within the BioGA framework have direct
analogues in the biological context of genomic data analysis, ensuring
the algorithm is biologically meaningful.

| Mathematical Component        | Biological Concept              | Function                                                                                                                                |
|:------------------------------|:--------------------------------|:----------------------------------------------------------------------------------------------------------------------------------------|
| **Population Initialization** | Initial Biological Variability  | Starting the search process by sampling observed expression levels.                                                                     |
| **Fitness ($f_{1},f_{2}$)**   | Functional Efficacy & Parsimony | $f_{1}$ ensures the optimized solution is relevant to observed data; $f_{2}$ enforces biological simplicity (Minimal Regulatory Model). |
| **Network Constraints**       | Gene-Gene Interactions          | Incorporating known biological constraints (e.g., regulatory or metabolic pathways) to guide mutation.                                  |
| **Clustering**                | Co-expressed Gene Modules       | Respecting known functional groupings during initialization.                                                                            |

## Conclusion

The `BioGA` package implements a theoretically sound and computationally
structured multi-objective evolutionary algorithm. Its mathematical
foundation—encompassing NSGA-II principles for multi-objective handling,
sophisticated SBX and adaptive network-constrained mutation for
exploration, and diversity-preserving replacement—provides a rigorous
framework for simultaneously achieving data fidelity and model sparsity
in genomic data optimization.

**Session Info**

``` r
sessioninfo::session_info()
#> ─ Session info ───────────────────────────────────────────────────────────────
#>  setting  value
#>  version  R version 4.5.2 (2025-10-31)
#>  os       Ubuntu 24.04.3 LTS
#>  system   x86_64, linux-gnu
#>  ui       X11
#>  language en
#>  collate  C.UTF-8
#>  ctype    C.UTF-8
#>  tz       UTC
#>  date     2025-12-04
#>  pandoc   3.1.11 @ /opt/hostedtoolcache/pandoc/3.1.11/x64/ (via rmarkdown)
#>  quarto   NA
#> 
#> ─ Packages ───────────────────────────────────────────────────────────────────
#>  package     * version date (UTC) lib source
#>  bookdown      0.45    2025-10-03 [1] RSPM
#>  bslib         0.9.0   2025-01-30 [1] RSPM
#>  cachem        1.1.0   2024-05-16 [1] RSPM
#>  cli           3.6.5   2025-04-23 [1] RSPM
#>  desc          1.4.3   2023-12-10 [1] RSPM
#>  digest        0.6.39  2025-11-19 [1] RSPM
#>  evaluate      1.0.5   2025-08-27 [1] RSPM
#>  fastmap       1.2.0   2024-05-15 [1] RSPM
#>  fs            1.6.6   2025-04-12 [1] RSPM
#>  htmltools     0.5.8.1 2024-04-04 [1] RSPM
#>  htmlwidgets   1.6.4   2023-12-06 [1] RSPM
#>  jquerylib     0.1.4   2021-04-26 [1] RSPM
#>  jsonlite      2.0.0   2025-03-27 [1] RSPM
#>  knitr         1.50    2025-03-16 [1] RSPM
#>  lifecycle     1.0.4   2023-11-07 [1] RSPM
#>  pkgdown       2.2.0   2025-11-06 [1] any (@2.2.0)
#>  R6            2.6.1   2025-02-15 [1] RSPM
#>  ragg          1.5.0   2025-09-02 [1] RSPM
#>  rlang         1.1.6   2025-04-11 [1] RSPM
#>  rmarkdown     2.30    2025-09-28 [1] RSPM
#>  sass          0.4.10  2025-04-11 [1] RSPM
#>  sessioninfo   1.2.3   2025-02-05 [1] RSPM
#>  systemfonts   1.3.1   2025-10-01 [1] RSPM
#>  textshaping   1.0.4   2025-10-10 [1] RSPM
#>  xfun          0.54    2025-10-30 [1] RSPM
#>  yaml          2.3.11  2025-11-28 [1] RSPM
#> 
#>  [1] /home/runner/work/_temp/Library
#>  [2] /opt/R/4.5.2/lib/R/site-library
#>  [3] /opt/R/4.5.2/lib/R/library
#> 
#> ──────────────────────────────────────────────────────────────────────────────
```
