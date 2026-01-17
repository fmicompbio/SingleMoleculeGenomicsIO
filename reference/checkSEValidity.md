# Check internal consistency of SummarizedExperiment object

All assays with read-level data must have the same number and order of
the reads, which must also agree with the order in `se$QC` if that
exists. All assays must have the same column names, which must also
agree with the column names of the object, and the `sample` column in
the `colData`.

## Usage

``` r
checkSEValidity(se, verbose = FALSE)
```

## Arguments

- se:

  A `SummarizedExperiment object`.

- verbose:

  A logical scalar. If `TRUE`, report on progress.

## Value

Silently returns `NULL`. If the object is not valid, an error will be
raised.

## Author

Charlotte Soneson

## Examples

``` r
library(GenomicRanges)
modbamfiles <- system.file("extdata", c("6mA_1_10reads.bam",
                                        "6mA_2_10reads.bam"),
                           package = "SingleMoleculeGenomicsIO")
se <- readModBam(bamfiles = modbamfiles, regions = "chr1:6940000-6955000",
                 modbase = "a", verbose = FALSE,
                 variantPositions = GPos(seqnames = "chr1",
                                         pos = c(6940000, 6940500)),
                 BPPARAM = BiocParallel::SerialParam())
checkSEValidity(se)
```
