# Expand the rows of a `RangedSummarizedExperiment` to single base resolution.

Expand the rows of a `RangedSummarizedExperiment` to single base
resolution.

## Usage

``` r
expandSEToBaseSpace(
  se,
  regions = NULL,
  flank = 0L,
  seqinfo = NULL,
  keepAssays = getReadLevelAssayNames(se),
  ignore.strand = TRUE,
  assayNameAmbig = "mod_prob",
  verbose = FALSE
)
```

## Arguments

- se:

  [`RangedSummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/RangedSummarizedExperiment-class.html)
  object to be expanded to single base resolution.

- regions:

  A
  [`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  object with one or more regions defining the ranges for expanding
  `se`. Alternatively, the regions can be specified as character scalars
  (e.g. "chr1:1200-1300") that can be coerced into a `GRanges` object.
  If `NULL` (the default), `regions` is set to the range of the data in
  `se`.

- flank:

  Integer scalar. If `regions=NULL` and `flank` is non-zero, `flank`
  positions will be added on both sides of each range in the derived
  `regions`. This is useful to add a buffer of NA values between
  observed values in non-overlapping regions.

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
  `rowRanges(se)`. If `FALSE`, `se` is first filtered (using
  `filterPositions(..., filters = "repeated.positions")`) to resolve any
  positions represented on both strands.

- assayNameAmbig:

  A character scalar indicating the assay to use to decide which row to
  retain if multiple rows represent the same genomic position (on
  different strands). The row with the largest row sum in this assay is
  retained. Only used if `ignore.strand` is `FALSE`.

- verbose:

  Logical scalar indicating whether to print out progress messages.

## Value

A single-base resolution
[`RangedSummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/RangedSummarizedExperiment-class.html)
corresponding to `se`. If `ignore.strand` is `TRUE`, the return object
is unstranded. If `ignore.strand` is `FALSE`, all positions present in
the input object will retain the assigned strand, while all other
positions will be unstranded.

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
#> ⠙ 0.000 Mio. genomic positions processed (0.001 Mio./s) [2ms]

#> ⠙ 0.000 Mio. genomic positions processed (0.001 Mio./s) [2ms]
#> ℹ finding unique genomic positions...
#> ✔ finding unique genomic positions... [49ms]
#> 
#> ⠙ 0.000 Mio. genomic positions processed (0.001 Mio./s) [2ms]

#> ℹ collapsed 11300 positions to 4772 unique ones
#> ✔ collapsed 11300 positions to 4772 unique ones [407ms]
#> 
#> ⠙ 0.000 Mio. genomic positions processed (0.001 Mio./s) [2ms]

se_exp <- expandSEToBaseSpace(se)
dim(se)
#> [1] 4772    1
dim(se_exp)
#> [1] 15802     1
```
