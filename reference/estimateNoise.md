# Estimation of noise variance for time series.

Given a time series signal, potentially with missing/unobserved values,
this function estimates the level of noise under the assumption of local
continuity/smoothness.

## Usage

``` r
estimateNoise(probs, read_pos, k, min_diffs)
```

## Arguments

- probs:

  Numeric vector of (observed) time series measurements.

- read_pos:

  Integer vector of measurement positions in the time series (same
  length as probs).

- k:

  Integer, maximum gap size tolerated when computing \\\Delta x\\.

- min_diffs:

  Integer, minimum number of lag-1 differences required to produce an
  estimate (default -1 = auto).

## Value

A vector with the following items:

1.  **Mean signal**

2.  **Total variance**

3.  **Noise variance estimate**

4.  **Number of adjacent positions used for the estimate**

## Details

\\\mathrm{Noise variance} \approx 0.5\\\mathrm{Var}(\Delta x)\\ where
\\\Delta x\\ are lag-1 differences that may skip up to *k* missing
values. This follows from error propagation and the assumption of low
varying x in adjacent measurements.

## Examples

``` r
estimateNoise(
    c(0.1, 0.25, 0.3, 0.45, 0.5, 0.7, 0.7),
    c(1L, 2L, 3L ,4L, 6L, 8L, 11L),
    2L, 1)
#>       mean      total  noise_raw     ndiffs 
#> 0.42857143 0.05154762 0.00300000 6.00000000 
```
