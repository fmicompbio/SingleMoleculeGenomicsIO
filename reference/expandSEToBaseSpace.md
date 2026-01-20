# Expand the rows of a `RangedSummarizedExperiment` to single base resolution.

Expand the rows of a `RangedSummarizedExperiment` to single base
resolution.

## Usage

``` r
expandSEToBaseSpace(
  se,
  region = NULL,
  seqinfo = NULL,
  keepAssays = getReadLevelAssayNames(se),
  ignore.strand = TRUE
)
```

## Arguments

- se:

  [`RangedSummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/RangedSummarizedExperiment-class.html)
  object to be expanded to single base resolution.

- region:

  A
  [`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  object with a single region defining the range for expanding `se`.
  Alternatively, the region can be specified as a character scalar (e.g.
  "chr1:1200-1300") that can be coerced into a `GRanges` object. If
  `NULL` (the default), `region` is set to the range of the data in
  `se`.

- seqinfo:

  `NULL` or a
  [`Seqinfo`](https://rdrr.io/pkg/Seqinfo/man/Seqinfo-class.html) object
  containing information about the set of genomic sequences
  (chromosomes). Alternatively, a named numeric vector with genomic
  sequence names and lengths. Used to convert a character `region` to a
  `GRanges` object.

- keepAssays:

  Character vector indicating which (read-level) assays to expand to
  base space. Only these assays will be present in the returned object.

- ignore.strand:

  A logical scalar defining whether to ignore the strand of the
  `rowRanges(se)`.

## Value

A single-base resolution
[`RangedSummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/RangedSummarizedExperiment-class.html)
corresponding to `se`.

## Author

Charlotte Soneson, Michael Stadler

## Examples

``` r
modbamfile <- system.file("extdata", "6mA_1_10reads.bam",
                          package = "SingleMoleculeGenomicsIO")
se <- readModBam(bamfiles = modbamfile, regions = "chr1:6940000-6955000",
                 modbase = "a", verbose = TRUE,
                 BPPARAM = BiocParallel::SerialParam())
#> ℹ extracting base modifications from modBAM files
#> ℹ finding unique genomic positions...
#> ✔ finding unique genomic positions... [23ms]
#> 
#> ℹ collapsed 11300 positions to 4772 unique ones
#> ✔ collapsed 11300 positions to 4772 unique ones [161ms]
#> 
se_exp <- expandSEToBaseSpace(se)
dim(se)
#> [1] 4772    1
dim(se_exp)
#> [1] 15802     1
```
