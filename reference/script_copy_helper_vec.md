# Copy script and other files from a package to a working directory

Copy script and other files from a package to a working directory

## Usage

``` r
script_copy_helper_vec(
  runscripts,
  workdir = getwd(),
  packagename = "prolfqua",
  overwrite = TRUE
)
```

## Arguments

- runscripts:

  paths of the files relative to the package root (or `inst/`)

- workdir:

  directory to copy the files to

- packagename:

  package to copy the files from

- overwrite:

  overwrite existing files in `workdir`, default TRUE

## Value

paths of the copied files

## Examples

``` r
copied <- script_copy_helper_vec("CITATION", workdir = tempdir())
#> copy /home/runner/work/_temp/Library/prolfqua/CITATION to /tmp/RtmpxEIdvD/CITATION
#> your working directory now should contain: 1 new files:
basename(copied)
#> [1] "CITATION"
```
