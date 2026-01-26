# Convert character region(s) to a `GRanges` object

This function takes a character vector with one or several strings
corresponding to genomic regions and converts them to a
[`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
object. The supported forms of strings are the same as the ones
supported by the `htslib` C API that is for example used in `samtools`:

- "REF" or "REF:":

  : All of seqname REF

- "REF:START":

  : Seqname REF from START to end of REF

- "REF:-END":

  : Seqname REF from 1 to END

- "REF:START-END":

  : Seqname REF from START to END

- ".":

  : All seqnames from 1 to their ends

Please note that the returned `GRanges` is generally parallel to the
elements of `regions`, with the exception of `"."`: This special region
is only allowed as a length-one `regions` argument and will possibly
result several returned regions corresponding to the sequences in
`seqinfo`.

## Usage

``` r
regionStringToGRanges(regions, seqinfo = NULL, maxend = .Machine$integer.max)
```

## Arguments

- regions:

  Character vector with region strings. See details for supported forms.

- seqinfo:

  One of `NULL`, an object for which a
  [`seqlengths`](https://rdrr.io/pkg/Seqinfo/man/seqinfo.html) method is
  available, such as a `BSgenome`, `Seqinfo` or `SummarizedExperiment`
  object or a named numeric vector with genomic sequence names and
  lengths. If not `NULL`, it will be used to obtain "END" for `regions`
  that do not specify it. Otherwise, `maxend` will be used for "END".

- maxend:

  Integer scalar giving the maximal value for END in cases where it is
  not given in `regions` (for example "REF:START") and was also not
  provided through other parameters.

## Value

A [`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
object.

## Author

Michael Stadler

## Examples

``` r
regionStringToGRanges("chr1:6940000-6955000")
#> GRanges object with 1 range and 0 metadata columns:
#>       seqnames          ranges strand
#>          <Rle>       <IRanges>  <Rle>
#>   [1]     chr1 6940000-6955000      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
regionStringToGRanges("chr1", seqinfo = c(chr1 = 999))
#> GRanges object with 1 range and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chr1     1-999      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
regionStringToGRanges("chr1")
#> GRanges object with 1 range and 0 metadata columns:
#>       seqnames       ranges strand
#>          <Rle>    <IRanges>  <Rle>
#>   [1]     chr1 1-2147483647      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
```
