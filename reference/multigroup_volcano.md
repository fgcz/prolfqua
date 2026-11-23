# Volcano plot for multiple contrasts

Volcano plot for multiple contrasts

## Usage

``` r
multigroup_volcano(
  data,
  effect = "fc",
  significance = "p.adjust",
  contrast = "condition",
  colour = "colour",
  xintercept = c(-2, 2),
  yintercept = 0.05,
  label = NULL,
  size = 3,
  segment.size = 0.3,
  segment.alpha = 0.3,
  scales = "fixed",
  max_nr_significant_text = 20,
  .data
)
```

## Arguments

- data:

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

  label text size in mm (default 3, about 8.5 pt)

- segment.size:

  controls size of lines

- segment.alpha:

  controls visibility of lines

- scales:

  parameter to ggplot2::facet_wrap

- max_nr_significant_text:

  maximum number of significant labels to display

- .data:

  deprecated; use `data`

## Value

a ggplot object with one volcano panel per contrast; if `label` is
given, significant points are labelled

## Note

Deprecated in prolfqua 1.9.5; use `ModelInterface$coef_volcano()`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
