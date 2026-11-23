# Plot peptide intensities of protein as a function of the sample and factor

Plot peptide intensities of protein as a function of the sample and
factor

## Usage

``` r
plot_hierarchies_line(res, protein_name, lfqdata, show.legend = FALSE)
```

## Arguments

- res:

  data.frame

- protein_name:

  title of plot

- lfqdata:

  [LFQData](https://fgcz.github.io/prolfqua/reference/LFQData.md) object
  with the data and its configuration

- show.legend:

  logical, show legend in plot

## Value

a ggplot object with one line per feature across samples, coloured by
peptide, faceted by factor levels and titled `protein_name`

## Note

Deprecated in prolfqua 1.9.5; use `LFQData$get_Aggregator()$plot()`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
