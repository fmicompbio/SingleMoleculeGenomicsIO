# Read base modifications from "C-to-T" bam file(s)

Parse read-to-genome mismatches and return a
[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
object with information on base states. Base mismatches are typically
resulting from single molecule genomics experiments, representing the
accessibility of individual bases: For bisulfite-sequencing experiments,
C-to-T mismatches correspond to unmethylated (inaccessible) bases
(`readBaseUnmod="T", readBaseMod="C"`, see below), while in
deaminase-treatment based experiments, C-to-T mismatches represent
modified (accessible) bases (`readBaseUnmod="C", readBaseMod="T"`).
Secondary and supplementary alignments are ignored.

## Usage

``` r
readMismatchBam(
  bamfiles,
  bamFormat = "QuasR",
  regions = NULL,
  sequenceContext = "GCH",
  readBaseUnmod = "T",
  readBaseMod = "C",
  level = "read",
  overlapAggregation = "maxQscore",
  sampleAnnot = NULL,
  nAlnsToSample = 0,
  seqnamesToSampleFrom = "chr19",
  seqinfo = NULL,
  sequenceReference = NULL,
  variantPositions = NULL,
  trim = FALSE,
  maxCoverage = NULL,
  BPPARAM = MulticoreParam(4L, RNGseed = 42L),
  verbose = FALSE
)
```

## Arguments

- bamfiles:

  Character vector with one or several paths of `BAM` files, containing
  alignments with specific (e.g., C-to-T) mismatches. If `bamfiles` is a
  named vector, the names are used as sample names and prefixes for read
  names. Otherwise, the prefixes will be `s1`, ..., `sN`, where `N` is
  the length of `bamfiles`. All `bamfiles` must have an index.

- bamFormat:

  A character scalar giving the format of `BAM` files. Currently
  supported ar `BAM` files that have been created with `"QuasR"` (using
  `qAlign` in bisulfite mode) or `"Bismark"`.

- regions:

  A
  [`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  object specifying which genomic regions to extract the reads from.
  Alternatively, regions can be specified as a character vector (e.g.
  "chr1:1200-1300", "chr2:-6000" or "chrM") that will be coerced into a
  `GRanges` object. If the end coordinate is not provided (for example
  in "chr1:10-", which means "to the end of chr1"), it will be obtained
  from `seqinfo` or default to a large value if `seqinfo` is not
  provided. Note that unless `trim=TRUE`, the reads are not trimmed to
  the boundaries of the specified ranges. As a result, returned
  positions will typically extend out of the specified regions. If
  `nAlnsToSample` is set to a non-zero value, `regions` is ignored.

- sequenceContext:

  A character scalar with an odd number of characters, specifying the
  genomic context on the plus strand, for which to report (mis-)matches.
  The string may contain IUPAC ambiguity codes, e.g.
  `sequenceContext="GCH"` would correspond to all GpC dinucleotides,
  excluding C's that are in CpG dinucleotides. The call will be made by
  comparing the genomic base in the middle of `sequenceContext` to the
  aligned base in the read and interpreted according to `readBaseUnmod`
  and `readBaseMod`. The genomic sequence context will be returned in
  `rowData(x)$sequenceContext`.

- readBaseUnmod, readBaseMod:

  Character scalars defining the read bases that are interpreted as
  unmodified or modified, respectively, when aligned to the middle base
  of `sequenceContext`.

- level:

  Character scalar specifying the level of returned modification data.
  Supported values are:

  "read"

  :   : Extracts modification probabilities for individual reads into an
      assay called `"mod_prob"`, with each column (sample) consisting of
      a position-by-read
      [`NaMatrix`](https://rdrr.io/pkg/SparseArray/man/NaArray-class.html).
      This is the default.

  "quickread"

  :   : Like "read", but using pileup-based data rather than parsing
      individual alignments. This mode does not support variant labels
      or read sampling.

  "summary"

  :   : Counts the total and modified bases for each position and strand
      and returns them in assays named `"Nvalid"` and `"Nmod"`,
      respectively, as well as the `Nmod/Nvalid` ratio in the assay
      `"FracMod"`.

- overlapAggregation:

  Character scalar indicating how to aggregate modification
  probabilities for positions that are covered by both mates in a
  paired-end bam file. Only used if `level = "read"`. In other modes,
  aggregation corresponding to `overlapAggregation = "maxQscore"` is
  performed automatically. Currently supported values are `"maxQscore"`.

- sampleAnnot:

  A `data.frame` (or `NULL`) providing annotations for the samples. It
  must contain at least one column, named `"sample"`, which must contain
  all the values of `names(bamfiles)`. The provided annotations will be
  propagated to the returned `SummarizedExperiment` object.

- nAlnsToSample:

  A numeric scalar. If non-zero, `regions` is ignored and approximately
  `nAlnsToSample` randomly selected alignments on `seqnamesToSampleFrom`
  are read from each of the `bamfiles`. In order to make the results
  reproducible, make sure to set the `RNGseed` argument in the provided
  `BPPARAM` object (see below). Please note that secondary and
  supplementary alignments in `bamfiles` contribute to the total number
  of alignments but will not be sampled, thus the number of returned
  alignments may be lower than `nAlnsToSample`. Also, sampled reads that
  do not overlap a site with the indicated sequence context will not be
  returned, which may further reduce the number of returned alignments.
  Note that for paired-end bam files, individual reads are sampled and
  pairs will not be complete.

- seqnamesToSampleFrom:

  A character vector with one or several sequence names (chromosomes)
  from which to sample alignments from (only used if `nAlnsToSample` is
  greater than zero).

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

- variantPositions:

  An optional `GPos` object with seqnames and coordinates of single
  nucleotide variant positions, to be used to construct read labels for
  allele-specific analysis. Ignored if `NULL` or `nAlnsToSample > 0`
  (sampling-mode).

- trim:

  A logical scalar. If `TRUE`, the returned `SummarizedExperiment`
  object will only contain the positions overlapping the specified
  `regions`. If `FALSE` (default), the object will be extended to all
  positions covered by the reads overlapping `regions`. In both cases,
  only reads overlapping the specified `regions` are included.

- maxCoverage:

  Integer scalar used to increase the maximal coverage for which
  samtools pileup will allocate cache memory. Ignored if `level="read"`.
  A value of NULL will use the samtools default (at the time of writing
  8000). Large values will increase memory consumption. If the actual
  coverage is larger than this value, alignments may be silently
  ignored.

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

  Logical scalar. If `TRUE`, report on progress.

## Value

A
[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
object with genomic positions in rows and samples in columns. The assays
depend on the value of the `level` argument (see above).

## Author

Michael Stadler, Charlotte Soneson

## Examples

``` r
bamfile <- system.file("extdata", "BisSeq_quasr_single.bam",
                       package = "SingleMoleculeGenomicsIO")
reffile <- system.file("extdata", "reference.fa.gz",
                       package = "SingleMoleculeGenomicsIO")
se <- readMismatchBam(bamfiles = bamfile, regions = "chr1:6940000-6955000",
                      sequenceReference = reffile, sequenceContext = "NCG",
                      readBaseUnmod = "T", readBaseMod = "C",
                      verbose = TRUE, BPPARAM = BiocParallel::SerialParam())
#> ℹ finding positions with NCG
#> ℹ extracting base mismatches from BAM files
#> ℹ finding positions with NCG
#> ℹ opening input file /Users/runner/work/_temp/Library/SingleMoleculeGenomicsIO/extdata/BisSeq_quasr_single.bam using 1 thread
#> ℹ finding positions with NCG
#> ℹ reading alignments overlapping 1 region
#> ℹ finding positions with NCG
#> ℹ removed 0 unaligned (e.g. soft-masked) of 0 called bases
#> ℹ finding positions with NCG
#> ℹ read 184 alignments
#> ℹ finding positions with NCG
#> ✔ finding positions with NCG [102ms]
#> 
#> ⠙ 0.000 Mio. genomic positions processed (0.001 Mio./s) [3ms]
#> ℹ finding unique genomic positions...
#> ✔ finding unique genomic positions... [43ms]
#> 
#> ⠙ 0.000 Mio. genomic positions processed (0.001 Mio./s) [3ms]
#> ℹ collapsed 50 positions to 38 unique ones
#> ✔ collapsed 50 positions to 38 unique ones [37ms]
#> 
#> ⠙ 0.000 Mio. genomic positions processed (0.001 Mio./s) [3ms]
#> ℹ extracting sequence contexts
#> ✔ extracting sequence contexts [1.1s]
#> 
#> ⠙ 0.000 Mio. genomic positions processed (0.001 Mio./s) [3ms]
```
