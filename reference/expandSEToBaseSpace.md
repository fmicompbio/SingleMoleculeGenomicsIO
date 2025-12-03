# Expand the rows of a `RangedSummarizedExperiment` to single base resolution.

Expand the rows of a `RangedSummarizedExperiment` to single base
resolution.

## Usage

``` r
expandSEToBaseSpace(
  se,
  region = NULL,
  seqinfo = NULL,
  keepAssays = .getReadLevelAssayNames(se),
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

  `NULL` or a `Seqinfo` object containing information about the set of
  genomic sequences (chromosomes). Alternatively, a named numeric vector
  with genomic sequence names and lengths. Used to convert a character
  `region` to a `GRanges` object.

- ignore.strand:

  A logical scalar defining whether to ignore the strand of the
  `rowRanges(se)`.

## Value

A single-base resolution
[`RangedSummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/RangedSummarizedExperiment-class.html)
corresponding to `se`.

## Author

Charlotte Soneson, Michael Stadler
