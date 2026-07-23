# Find prior degrees of freedom for DEqMS

Uses
[`limma::trigammaInverse`](https://rdrr.io/pkg/limma/man/trigammainverse.html)
to estimate d0 from the mean residual variance after removing the
count-dependent trend.

## Usage

``` r
find_d0_deqms(mean_myfct)
```

## Arguments

- mean_myfct:

  mean of squared residuals minus trigamma correction

## Value

scalar prior degrees of freedom (d0)
