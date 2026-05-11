# Count pairs of modified bases by distance and modification state for bam files where modifications are indicated by sequence mismatches

For all pairs of bases with modification calls in a read (pair),
tabulate the number of pairs at a given distance with a given
modification state.

## Usage

``` r
countMismatchStatePairs(
  bamfile,
  bamFormat = "QuasR",
  regions = ".",
  sequenceContext = "GCH",
  nAlnsToSample = 0,
  seqnamesToSampleFrom = character(0),
  readBaseUnmod = "T",
  readBaseMod = "C",
  windowSize = 200,
  minMapQ = 0,
  minAlignedLength = 0,
  seqinfo = NULL,
  sequenceReference = NULL,
  BPPARAM = MulticoreParam(4L),
  verbose = FALSE
)
```

## Arguments

- bamfile:

  Character scalar giving the path to a `BAM` file, containing
  alignments with specific (e.g., C-to-T) mismatches. The `bamfile` must
  have an index.

- bamFormat:

  A character scalar giving the format of `BAM` files. Currently
  supported ar `BAM` files that have been created with `"QuasR"` (using
  `qAlign` in bisulfite mode) or `"Bismark"`.

- regions:

  A
  [`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  object specifying which genomic regions to extract the reads from.
  Alternatively, regions can be specified as a character vector (e.g.
  "chr1:1200-1300", "chr2:-6000", "chr1:10-", "chrM" or ".").

- sequenceContext:

  A character scalar with an odd number of characters, specifying the
  genomic context on the plus strand, for which to report (mis-)matches.
  The string may contain IUPAC ambiguity codes, e.g.
  `sequenceContext="GCH"` would correspond to all GpC dinucleotides,
  excluding C's that are in CpG dinucleotides. The call will be made by
  comparing the genomic base in the middle of `sequenceContext` to the
  aligned base in the read and interpreted according to `readBaseUnmod`
  and `readBaseMod`.

- nAlnsToSample:

  A numeric scalar. If non-zero, `regions` is ignored and approximately
  `nAlnsToSample` randomly selected alignments on `seqnamesToSampleFrom`
  are read from the `bamfile`. If `bamfile` contains paired alignments,
  they are included or excluded as pairs. Alignments are counted
  individually towards `nAlnsToSample`, thus for paired-end data this
  parameter has to be set to twice the number of pairs to be sampled. In
  order to make the results reproducible, make sure to set the random
  number seed using `set.seed`. Please note that secondary and
  supplementary alignments in `bamfiles` contribute to the total number
  of alignments but will not be sampled, thus the number of used
  alignments may be lower than `nAlnsToSample`.

- seqnamesToSampleFrom:

  A character vector with one or several sequence names (chromosomes)
  from which to sample alignments from (only used if `nAlnsToSample` is
  greater than zero).

- readBaseUnmod, readBaseMod:

  Character scalars defining the read bases that are interpreted as
  unmodified or modified, respectively, when aligned to the middle base
  of `sequenceContext`.

- windowSize:

  Numeric scalar giving the maximum window size covering pairs of
  modified bases to consider in pair-counting mode. A window size of 1
  corresponds to a single base, a size of 2 to directly adjacent bases,
  etc.

- minMapQ:

  Numeric scalar giving the minimal mapping quality to include
  alignments in pair-counting mode.

- minAlignedLength:

  Numeric scalar giving the minimal alignment length to include
  alignments in pair-counting mode.

- seqinfo:

  `NULL` or a
  [`Seqinfo`](https://rdrr.io/pkg/Seqinfo/man/Seqinfo-class.html) object
  containing information about the set of genomic sequences
  (chromosomes). Alternatively, a named numeric vector with genomic
  sequence names and lengths. Useful to set the sorting order of
  sequence names.

- sequenceReference:

  A [`BSgenome`](https://rdrr.io/pkg/BSgenome/man/BSgenome-class.html)
  object, or a character scalar giving the path to a fasta formatted
  file with reference sequences, or a
  [`DNAStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)
  object. The sequence context (see `sequenceContextWidth` argument)
  will be extracted from these sequences.

- BPPARAM:

  A
  [`BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  object that controls the number of parallel CPU threads to use for
  some of the steps in
  [`readModBam()`](https://fmicompbio.github.io/SingleMoleculeGenomicsIO/reference/readModBam.md).
  The default value is
  ([`MulticoreParam`](https://rdrr.io/pkg/BiocParallel/man/MulticoreParam-class.html)`(4L, RNGseed = 42L)`).
  If randomly sampling reads (`nAlnsToSample > 0`), make sure to set the
  `RNGseed` argument when constructing the `BPPARAM` object for
  reproducible results (see also
  [`vignette("Random_Numbers", package = "BiocParallel")`](https://bioconductor.org/packages/release/bioc/vignettes/BiocParallel/inst/doc/Random_Numbers.html)).

- verbose:

  A logical scalar. If `TRUE`, report on progress.

## Value

A `DataFrame` with `windowSize` rows and five columns. The `metadata`
slot contains header information from `bamfile`.

## Author

Charlotte Soneson, Michael Stadler

## Examples

``` r
bamfile <- system.file("extdata", "BisSeq_quasr_single.bam",
                       package = "SingleMoleculeGenomicsIO")
reffile <- system.file("extdata", "reference.fa.gz",
                       package = "SingleMoleculeGenomicsIO")
tbl <- countMismatchStatePairs(bamfile = bamfile, bamFormat = "QuasR",
                               regions = "chr1:6940000-6955000",
                               sequenceContext = "C",
                               readBaseUnmod = "T", readBaseMod = "C",
                               windowSize = 200,
                               sequenceReference = reffile,
                               BPPARAM = BiocParallel::SerialParam(),
                               verbose = TRUE)
#> ℹ finding positions with C
#> ℹ opening input file /Users/runner/work/_temp/Library/SingleMoleculeGenomicsIO/extdata/BisSeq_quasr_single.bam using 1 thread
#> ℹ finding positions with C

#> ℹ counting state-pairs for alignments
#> ℹ finding positions with C

#> ℹ read 184 alignments
#> ℹ finding positions with C

#> ✔ finding positions with C [64ms]
#> 
#> ⠙ 0.000 Mio. genomic positions processed (0.001 Mio./s) [2ms]

tbl
#> DataFrame with 200 rows and 5 columns
#>             S unmod_unmod unmod_mod mod_unmod   mod_mod
#>     <integer>   <numeric> <numeric> <numeric> <numeric>
#> 1           1        3242         0         0        68
#> 2           2         691        10         5         0
#> 3           3         561         9        10         1
#> 4           4         563        12        15         0
#> 5           5         533        12        13         0
#> ...       ...         ...       ...       ...       ...
#> 196       196           0         0         0         0
#> 197       197           0         0         0         0
#> 198       198           0         0         0         0
#> 199       199           0         0         0         0
#> 200       200           0         0         0         0
```
