# Registry of available contrast facade classes

A named list mapping short names to facade class names and their data
requirements. Each entry has:

- class:

  Character string naming the R6 facade class

- needs:

  One of `"aggregated"`, `"nested"`, or `"either"`

## Usage

``` r
FACADE_REGISTRY
```

## Format

An object of class `list` of length 8.
