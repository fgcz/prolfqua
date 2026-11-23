# Encode intensities as a binary detection response

Encode intensities as a binary detection response

## Usage

``` r
encode_bin_resp(lfqdata, name = "bin_resp")
```

## Arguments

- lfqdata:

  [LFQData](https://fgcz.github.io/prolfqua/reference/LFQData.md) object
  with the data and its configuration

- name:

  column name for binary response

## Value

a tibble of the long-format data with an added integer column `name`: 1
if the response is observed, 0 if it is `NA`

## Note

Deprecated in prolfqua 1.9.5; use
`build_contrast_analysis(method = "firth")`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
