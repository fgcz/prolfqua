# Find a file stored in a package

Looks for `file` in the installed package and, for a source package, in
its `inst/` folder.

## Usage

``` r
find_package_file(packagename, file)
```

## Arguments

- packagename:

  name of the R package

- file:

  relative path to the file within the package

## Value

the path to `file` in the installed package (or in its `inst/` folder of
a source package), or `NULL` if the file does not exist.

## Examples

``` r
find_package_file("prolfqua", "CITATION")
#> [1] "/home/runner/work/_temp/Library/prolfqua/CITATION"
```
