# Simulate protein-level data with two groups

Draws protein abundances from a log-normal distribution and multiplies
them by group-specific fold changes for down-, up- and not-regulated
proteins; with `peptide = TRUE`, each protein gets a geometric number of
peptides. The `_config` variants add an
[AnalysisConfiguration](https://fgcz.github.io/prolfqua/reference/AnalysisConfiguration.md)
and missing values.

## Usage

``` r
sim_lfq_data(
  n_proteins = 20,
  n_per_group = 4,
  fc = list(A = c(D = -2, U = 2, N = 0), B = c(D = 1, U = -4)),
  prop = list(A = c(D = 10, U = 10), B = c(D = 5, U = 20)),
  mean_prot = 20,
  sdlog = log(1.2),
  probability_of_success = 0.3,
  peptide = FALSE,
  Nprot,
  N,
  PEPTIDE
)
```

## Arguments

- n_proteins:

  one positive whole number of proteins

- n_per_group:

  one positive whole number of samples per group

- fc:

  named list with one named numeric vector per group: fold changes of
  the down (D), up (U) and not (N) regulated proteins

- prop:

  named list with one named numeric vector per group: percentages of
  down (D), up (U) and not (N) regulated proteins

- mean_prot:

  one positive number, the mean protein abundance on the natural scale

- sdlog:

  one non-negative number, the standard deviation of the log-normal
  abundance distribution

- probability_of_success:

  one number in (0, 1\], the probability of success of the geometric
  distribution of peptide counts

- peptide:

  one TRUE/FALSE: TRUE simulates peptide-level data, FALSE protein-level
  data

- Nprot, N, PEPTIDE:

  deprecated; use `n_proteins`, `n_per_group` and `peptide`

## Value

A data frame containing simulated LFQ data.

## Examples

``` r
res <- sim_lfq_data(n_proteins = 10)
res <- sim_lfq_data(n_proteins = 10, peptide = TRUE)
res <- sim_lfq_data(n_proteins = 10, n_per_group = 4)
```
