# simulate protein level data with two groups

simulate protein level data with two groups

## Usage

``` r
sim_lfq_data(
  Nprot = 20,
  N = 4,
  fc = list(A = c(D = -2, U = 2, N = 0), B = c(D = 1, U = -4)),
  prop = list(A = c(D = 10, U = 10), B = c(D = 5, U = 20)),
  mean_prot = 20,
  sdlog = log(1.2),
  probability_of_success = 0.3,
  PEPTIDE = FALSE
)
```

## Arguments

- Nprot:

  number of proteins

- N:

  group size

- fc:

  D - down regulation N - matrix, U - regulation

- prop:

  proportion of down (D), up (U) and not regulated (N)

- mean_prot:

  mean protein abundance on the natural scale

- sdlog:

  standard deviation of the log-normal abundance distribution

- probability_of_success:

  probability of success for the geometric distribution of peptide
  counts

- PEPTIDE:

  if TRUE, simulate peptide-level data; if FALSE, protein-level

## Examples

``` r
res <- sim_lfq_data(Nprot = 10)
res <- sim_lfq_data(Nprot = 10, PEPTIDE = TRUE)
res <- sim_lfq_data(Nprot = 10, N=4)
```
