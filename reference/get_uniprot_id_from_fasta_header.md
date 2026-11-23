# Extract UniProt accessions from FASTA headers

Parses UniProt FASTA headers such as `sp|P12345|PROT1_HUMAN` in
`idcolumn` and adds the accession as the column `UniprotID`; identifiers
that are not FASTA headers get `NA`.

## Usage

``` r
get_uniprot_id_from_fasta_header(df, idcolumn = "protein_Id")
```

## Arguments

- df:

  data.frame

- idcolumn:

  character, column name to extract uniprot IDs from

## Value

data.frame

## See also

Other utilities:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://fgcz.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`is_contaminant()`](https://fgcz.github.io/prolfqua/reference/is_contaminant.md),
[`is_decoy()`](https://fgcz.github.io/prolfqua/reference/is_decoy.md),
[`table_facade()`](https://fgcz.github.io/prolfqua/reference/table_facade.md)

## Examples

``` r
df <- data.frame(
  protein_Id = c("sp|P12345|PROT1_HUMAN", "sp|Q67890|PROT2_HUMAN",
                 "tr|A0A0A0|PROT3_HUMAN", "NOFASTA_ID"),
  value = 1:4
)
result <- get_uniprot_id_from_fasta_header(df, idcolumn = "protein_Id")
stopifnot("UniprotID" %in% colnames(result))
stopifnot(sum(!is.na(result$UniprotID)) == 3)
```
