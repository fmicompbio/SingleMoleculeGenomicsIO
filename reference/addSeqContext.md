# Add sequence context around positions of interest to a SummarizedExperiment

Convenience function to extract sequence context around positions of
interest (the `rowRanges` of a
[`RangedSummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/RangedSummarizedExperiment-class.html))
and add them the the
[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)'s
row data (`rowData(se)$sequenceContext`). The extracted sequences will
correspond to the regions defined as
`resize(rowRanges(x), width = sequenceContextWidth, fix = "center"`.
Sequence contexts are extracted using
[`extractSeqContext`](https://fmicompbio.github.io/SingleMoleculeGenomicsIO/reference/extractSeqContext.md).

## Usage

``` r
addSeqContext(x, sequenceContextWidth, sequenceReference)
```

## Arguments

- x:

  A
  [`RangedSummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/RangedSummarizedExperiment-class.html).

- sequenceContextWidth:

  A numeric scalar giving the width of the sequence context to be
  extracted from the reference (`sequenceReference` argument). This must
  be an odd number so that the sequence can be centered on the modified
  base. If `sequenceContextWidth = 0` (the default), no sequence context
  will be extracted.

- sequenceReference:

  A [`BSgenome`](https://rdrr.io/pkg/BSgenome/man/BSgenome-class.html)
  object, or a character scalar giving the path to a fasta formatted
  file with reference sequences, or a
  [`DNAStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)
  object. The sequence context (see `sequenceContextWidth` argument)
  will be extracted from these sequences.

## Value

A
[`RangedSummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/RangedSummarizedExperiment-class.html)
object with sequence contexts added as a
[`DNAStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)
object to `rowData(x)$sequenceContext`.

## See also

[`extractSeqContext`](https://fmicompbio.github.io/SingleMoleculeGenomicsIO/reference/extractSeqContext.md)

## Author

Michael Stadler

## Examples

``` r
# load package
library(SummarizedExperiment)

# file with sequence in fasta format of length 6957060
reffile <- system.file("extdata", "reference.fa.gz",
                       package = "SingleMoleculeGenomicsIO")

# define some regions at the end of the reference sequence
se <- SummarizedExperiment(
          assays = matrix(1:3, ncol=1),
          rowRanges = GRanges(
              "chr1", IRanges(start = 6957060 - c(4, 2, 0),
              width = 1, names = c("a", "b", "c")),
              strand = "-"))

# add sequence context (note the padding with N's)
rowRanges(se)
#> GRanges object with 3 ranges and 0 metadata columns:
#>     seqnames    ranges strand
#>        <Rle> <IRanges>  <Rle>
#>   a     chr1   6957056      -
#>   b     chr1   6957058      -
#>   c     chr1   6957060      -
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
se <- addSeqContext(se, 7, reffile)
rowRanges(se)
#> GRanges object with 3 ranges and 1 metadata column:
#>     seqnames    ranges strand | sequenceContext
#>        <Rle> <IRanges>  <Rle> |  <DNAStringSet>
#>   a     chr1   6957056      - |         CCCCTTT
#>   b     chr1   6957058      - |         NTCCCCT
#>   c     chr1   6957060      - |         NNNTCCC
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
