# Get names of colData columns containing read-level data

The names of `colData` column designated as containing read-level
annotations are extracted from
`metadata(se)$readLevelData$colDataColumns`.

## Usage

``` r
getReadLevelColDataNames(se)
```

## Arguments

- se:

  A `SummarizedExperiment` object.

## Value

A (possibly empty) character vector with the names of the columns of
colData(se) containing read-level data.

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
getReadLevelColDataNames(se)
#> [1] "readInfo" "QC"      
```
