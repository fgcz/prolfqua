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

## Note

Deprecated in prolfqua 1.9.5; use
`prolfquapp::script_copy_helper_vec()`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
