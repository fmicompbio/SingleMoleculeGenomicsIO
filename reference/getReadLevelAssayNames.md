# Get names of assays containing read-level data

The names of assays designated as containing read-level data are
extracted from `metadata(se)$readLevelData$assayNames`.

## Usage

``` r
getReadLevelAssayNames(se)
```

## Arguments

- se:

  A `SummarizedExperiment` object.

## Value

A (possibly empty) character vector with the names of the assays of se
containing read-level data.

## Author

Charlotte Soneson

## Examples

``` r
modbamfiles <- system.file("extdata", c("6mA_1_10reads.bam",
                                        "6mA_2_10reads.bam"),
                           package = "SingleMoleculeGenomicsIO")
se <- readModBam(bamfiles = modbamfiles, regions = "chr1:6940000-6955000",
                 modbase = "a", verbose = FALSE,
                 BPPARAM = BiocParallel::SerialParam())
se <- addReadStats(se, BPPARAM = BiocParallel::SerialParam())
#> Warning: Too few points to estimate noise floor (1); raw noise variances are used.
#> Warning: Too few points to estimate noise floor (0); raw noise variances are used.
getReadLevelAssayNames(se)
#> [1] "mod_prob"
```
