# copy script files and other from a package to workdir

copy script files and other from a package to workdir

## Usage

``` r
script_copy_helper_vec(runscripts, workdir = getwd(), packagename = "prolfqua")
```

## Examples

``` r
copied <- script_copy_helper_vec("extdata/metadata.csv", workdir = tempdir())
#> copy /home/runner/work/_temp/Library/prolfqua/extdata/metadata.csv to /tmp/RtmpIFz9yR/metadata.csv
#> your working directory now should contain: 1 new files:
basename(copied)
#> [1] "metadata.csv"
```
