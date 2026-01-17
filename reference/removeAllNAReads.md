# Remove all-NA reads

Remove reads (columns) that contain only NA values from a read-level
assay (`DataFrame` of `NaArray`s).

## Usage

``` r
removeAllNAReads(x, prune = TRUE)
```

## Arguments

- x:

  A
  [`DataFrame`](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)
  with
  [`NaArray`](https://rdrr.io/pkg/SparseArray/man/NaArray-class.html)
  objects in its columns (typically a read-level assay as returned by
  [`readModBam`](https://fmicompbio.github.io/SingleMoleculeGenomicsIO/reference/readModBam.md)).

- prune:

  A logical scalar. If `TRUE` (the default), samples (columns of the
  `DataFrame`) for which the NA-read filtering retains none of the reads
  will be completely removed. If `FALSE`, such samples are retained as a
  zero-column `NAMatrix`).

## Value

A `DataFrame` of `NaArray`s with all columns containing at least one
non-`NA` value.

## Examples

``` r
library(SummarizedExperiment)
modbamfile <- system.file("extdata", "6mA_1_10reads.bam",
                          package = "SingleMoleculeGenomicsIO")
se <- readModBam(bamfiles = modbamfile, regions = "chr1:6940000-6955000",
                 modbase = "a", verbose = FALSE,
                 BPPARAM = BiocParallel::SerialParam())
removeAllNAReads(assay(se[1:10, ], "mod_prob"))
#> DataFrame with 10 rows and 1 column
#>             s1
#>     <NaMatrix>
#> 1            0
#> 2            0
#> 3            0
#> 4            0
#> 5  0.275390625
#> 6  0.076171875
#> 7            0
#> 8            0
#> 9            0
#> 10           0
```
