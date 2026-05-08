# Estimate a time series Signal to Noise Ratio (SNR) given total and noise variance components and (optionally) a background noise model.

Given time series variance components (total variance, raw noise
variance), compute the final noise, signal, and SNR. For noise
estimation use one of:

1.  **A raw noise variance** e.g the result of `estimateNoise`

2.  **A background noise model** based on linear predictors

3.  **A flooring rule**: final noise = max(raw, background)

The background noise model is a *general* linear noise baseline:
`baseline = sum_i beta[i] * feat[i]`.

## Usage

``` r
estimateSNR(totalVar, noiseRaw, eps, betas, features, noise_mode = "floor")
```

## Arguments

- totalVar:

  Total time series variance. Typically calculated with `estimateNoise`

- noiseRaw:

  noise variance estimate. Typically calculated with `estimateNoise`

- eps:

  Numeric. Lower bound applied to both noise and signal.

- betas, features:

  Coefficients (\\betas\\) and predictors (\\features\\) for the
  background noise model. Ignored when `noise_mode = "raw"`. Both
  vectors must have the same length when used. The predictors form a
  *design vector*: include a constant `1` if the model has an intercept.
  For example, for \\baseline = b0 + b1 \* mean\\, use
  `betas = c(b0, b1)`, `features = c(1, mean)`.

- noise_mode:

  Character scalar: one of `"raw"`, `"model"`, `"floor"`.

  - `"raw"` -\> use `noise_raw`

  - `"model"` -\> use `baseline = sum(betas * features)`

  - `"floor"` -\> use `max(noise_raw, baseline)`

## Value

Named numeric vector with elements: - snr (log2 scale) - signal - noise
final estimate used for signal and snr calculation - noise baseline
(sum(betas \* features); NA if not applied) - noise raw

## Examples

``` r
estimateSNR(totalVar = 1.0, noiseRaw = 1.0, eps = 0.01,
            betas = c(0.5, 2.0), features = c(2.0, 4.5),
            noise_mode = "model")
#>       snr    signal     noise  baseline       raw 
#> -9.965784  0.010000 10.000000 10.000000  1.000000 
```
