# Decorates LFQData with methods to impute missing values

Decorates LFQData with methods to impute missing values

Decorates LFQData with methods to impute missing values

## Details

Use \`lfqdata\$get_Imputer()\` to create an instance. The imputer works
on a deep copy of the data, so the original LFQData is not modified.

Typical workflow:


    imp <- lfqdata$get_Imputer()
    imp$impute(method = "multRepl", lod = "global")
    imp$lfq$get_Plotter()$pca()

## See also

Other LFQData:
[`LFQData`](https://wolski.github.io/prolfqua/reference/LFQData.md),
[`LFQDataAggregator`](https://wolski.github.io/prolfqua/reference/LFQDataAggregator.md),
[`LFQDataPlotter`](https://wolski.github.io/prolfqua/reference/LFQDataPlotter.md),
[`LFQDataStats`](https://wolski.github.io/prolfqua/reference/LFQDataStats.md),
[`LFQDataSummariser`](https://wolski.github.io/prolfqua/reference/LFQDataSummariser.md),
[`LFQDataToSummarizedExperiment()`](https://wolski.github.io/prolfqua/reference/LFQDataToSummarizedExperiment.md)

## Public fields

- `lfq`:

  LFQData

## Methods

### Public methods

- [`LFQDataImp$new()`](#method-LFQDataImp-new)

- [`LFQDataImp$impute()`](#method-LFQDataImp-impute)

- [`LFQDataImp$clone()`](#method-LFQDataImp-clone)

------------------------------------------------------------------------

### Method `new()`

initialize

#### Usage

    LFQDataImp$new(lfqdata)

#### Arguments

- `lfqdata`:

  LFQData object

------------------------------------------------------------------------

### Method `impute()`

Impute missing values using zCompositions

#### Usage

    LFQDataImp$impute(
      method = c("multRepl", "GBM", "SQ", "BL", "CZM"),
      lod = c("global", "quantile")
    )

#### Arguments

- `method`:

  imputation method (default "multRepl")

- `lod`:

  limit of detection strategy (default "global")

#### Returns

invisible(self) for chaining

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    LFQDataImp$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
dd <- prolfqua::sim_lfq_data_peptide_config()
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
lfqdata <- LFQData$new(dd$data, dd$config)
imp <- LFQDataImp$new(lfqdata)
if (requireNamespace("zCompositions", quietly = TRUE)) {
  imp$impute(method = "multRepl", lod = "global")
  wide <- imp$lfq$to_wide(as.matrix = TRUE)
  stopifnot(!any(is.na(wide$data)))
}
```
