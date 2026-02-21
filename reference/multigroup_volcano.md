# plot volcano given multiple contrasts

plot volcano given multiple contrasts

## Usage

``` r
multigroup_volcano(
  .data,
  effect = "fc",
  significance = "p.adjust",
  contrast = "condition",
  colour = "colour",
  xintercept = c(-2, 2),
  yintercept = 0.05,
  label = NULL,
  size = 1,
  segment.size = 0.3,
  segment.alpha = 0.3,
  scales = "fixed",
  maxNrOfSignificantText = 20
)
```

## Arguments

- .data:

  data in long format

- effect:

  column containing effect sizes

- significance:

  column containing p-values, q.values etc

- contrast:

  column with contrast

- colour:

  colouring of points

- xintercept:

  fc thresholds

- yintercept:

  significance threshold

- label:

  column containing labels

- size:

  controls size of text

- segment.size:

  controls size of lines

- segment.alpha:

  controls visibility of lines

- scales:

  parameter to ggplot2::facet_wrap

- maxNrOfSignificantText:

  maximum number of significant labels to display

## See also

Other utilities:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`get_UniprotID_from_fasta_header()`](https://wolski.github.io/prolfqua/reference/get_UniprotID_from_fasta_header.md),
[`matrix_to_tibble()`](https://wolski.github.io/prolfqua/reference/matrix_to_tibble.md),
[`names_to_matrix()`](https://wolski.github.io/prolfqua/reference/names_to_matrix.md),
[`pairs_smooth()`](https://wolski.github.io/prolfqua/reference/pairs_smooth.md),
[`pairs_w_abline()`](https://wolski.github.io/prolfqua/reference/pairs_w_abline.md),
[`panel_cor()`](https://wolski.github.io/prolfqua/reference/panel_cor.md),
[`panel_hist()`](https://wolski.github.io/prolfqua/reference/panel_hist.md),
[`remove_NA_rows()`](https://wolski.github.io/prolfqua/reference/remove_NA_rows.md),
[`table_facade()`](https://wolski.github.io/prolfqua/reference/table_facade.md)

## Examples

``` r

show <- data.frame(logFC = rnorm(100, 0, 1), adj.P.Val = runif(100, 0, 1),
  Condition = sample(c("a","b")), colour = "forward", Name = paste0("n", 1:100))
prolfqua::multigroup_volcano( show,
effect="logFC",
significance = "adj.P.Val",
contrast="Condition",
colour=NULL,
label="Name",
maxNrOfSignificantText = 300)
```
