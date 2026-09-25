# Get base coverage histogram for a BAM file

Get base coverage histogram for a BAM file

## Usage

``` r
getBaseCoverageForBam(
  bamfile,
  regions = NULL,
  maxDepth = 200L,
  method = "full",
  nThreads = 3L
)
```

## Arguments

- bamfile:

  A character scalar with the bam file name (and path).

- regions:

  Character vector specifying the region(s) for which to calculate
  coverage. Each region uses the grammar supported by the `htslib`
  function `sam_parse_region`, for example `"."` (all contigs), `"chr"`
  (whole contig), `"chr:START"`, `"chr:-END"` or `"chr:START-END"`. If
  `NULL`, the whole genome (`"."`) is used by default.

- maxDepth:

  An integer scalar defining the maximal depth to consider.

- method:

  Character scalar with the method used for coverage calculation.
  `"full"` (the default) considers CIGAR operations, thus not counting
  soft-clip bases and read-inserted positions as covered. `"simple"`
  ignores the CIGAR strings and covers the whole reference span between
  the first and the last aligned position of each alignment (using
  `htslib`'s `bam_endpos`); this is slightly faster but overestimates
  coverage in the soft-clipped ends and around indels.

- nThreads:

  A numeric scalar with the number of threads used for decompressing BAM
  records.

## Value

A named numeric vector of length `maxDepth + 1`, with values at index
`i` giving the number of positions that were overlapped by exactly `i-1`
alignments. Positions overlapped by more than `maxDepth` alignments are
also added to the value for `maxDepth` at index `maxDepth + 1`.

## Details

Secondary and supplementary alignments and unmapped reads are not
included.

## References

The algorithm was described in Pedersen BS and Quinlan AR. "Mosdepth:
quick coverage calculation for genomes and exomes". Bioinformatics.
2018; 34(5):867-868. <https://doi.org/10.1093/bioinformatics/btx699>

## Examples

``` r
modbamfile <- system.file("extdata", "6mA_1_10reads.bam", package = "SingleMoleculeGenomicsIO")
getBaseCoverageForBam(modbamfile, "chr1", 12L, "full")
#>         0         1         2         3         4         5         6         7 
#> 195138505      3841       811       563       973       596       675       463 
#>         8         9        10        11        12 
#>       420      1066      6366         0         0 
getBaseCoverageForBam(modbamfile, "chr1:6000000-7000000", 12L, "full")
#>      0      1      2      3      4      5      6      7      8      9     10 
#> 984227   3841    811    563    973    596    675    463    420   1066   6366 
#>     11     12 
#>      0      0 
getBaseCoverageForBam(modbamfile, "chr1:6000000-7000000", 12L, "simple")
#>      0      1      2      3      4      5      6      7      8      9     10 
#> 984199   3839    823    538    991    545    705    459    342    575   6985 
#>     11     12 
#>      0      0 
```
