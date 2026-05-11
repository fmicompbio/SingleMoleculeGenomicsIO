# Get list of read names by sample from SummarizedExperiment object

Get list of read names by sample from SummarizedExperiment object

## Usage

``` r
getReadNamesBySample(se)
```

## Arguments

- se:

  A
  [`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object.

## Value

A named list with one entry per sample, containing the read names for
the respective sample.

## Author

Charlotte Soneson

## Examples

``` r
extractfiles <- system.file("extdata",
                            c("modkit_extract_rc_6mA_1.tsv.gz",
                              "modkit_extract_rc_6mA_2.tsv.gz"),
                            package = "SingleMoleculeGenomicsIO")
se <- readModkitExtract(extractfiles, modbase = "a", filter = "modkit",
                        BPPARAM = BiocParallel::SerialParam())
getReadNamesBySample(se)
#> $s1
#>  [1] "s1-233e48a7-f379-4dcf-9270-958231125563"
#>  [2] "s1-d52a5f6a-a60a-4f85-913e-eada84bfbfb9"
#>  [3] "s1-fc4646ce-66f9-401f-b968-e9b0cda14d61"
#>  [4] "s1-92e906ae-cddb-4347-a114-bf9137761a8d"
#>  [5] "s1-6cf74134-e550-4c02-bd2b-91385422ee25"
#>  [6] "s1-5d45d8d2-d5f5-47ff-a9fa-f3fd6b7bd3c7"
#>  [7] "s1-b6fea9db-c92d-4152-9d29-4d021bbc45e8"
#>  [8] "s1-49c1e21e-8cb0-415a-aba9-92912219c4bb"
#>  [9] "s1-b0b20f04-931f-4f60-b3e4-0ee1f5666a61"
#> [10] "s1-41ca0e97-11b3-454b-9741-bc373e29ef37"
#> 
#> $s2
#>  [1] "s2-daac487b-5406-42b5-b882-9020f1b03752"
#>  [2] "s2-034b625e-6230-4f8d-a713-3a32cd96c298"
#>  [3] "s2-274d50aa-f060-4bcf-901e-4cab771295f6"
#>  [4] "s2-3bee7d0d-d4cd-4b38-9839-10e1a0fce5de"
#>  [5] "s2-c0f0cdaa-216c-4c8a-8a2a-7e5c1e8445f5"
#>  [6] "s2-d03efe3b-a45b-430b-9cb6-7e5882e4faf8"
#>  [7] "s2-04784c5b-e31b-421f-a14c-7072ec62b50d"
#>  [8] "s2-9c27e4b2-d929-494a-8878-60c8fca73bf1"
#>  [9] "s2-7fc51790-2c6e-49e9-938c-985850ff85e7"
#> [10] "s2-573752c9-f768-46ff-9887-be555db141dc"
#> 
#> attr(,"source")
#> [1] "assay mod_prob"
```
