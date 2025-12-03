# Add read segment annotation

Add segemnts for individual reads to the `colData` of a
`SummarizedExperiment` object with read-level data.

## Usage

``` r
annotateReadSegments(se, irlList, name)
```

## Arguments

- se:

  A
  [`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object.

- irlList:

  A list of
  [`IRangesList`](https://rdrr.io/pkg/IRanges/man/IRangesList-class.html)
  objects. Each list element corresponds to a sample (a columns in
  `se`), with the ranges in an `IRangesList` element giving the segments
  of an individual read (coordinates refer to the genome).

- name:

  A character scalar giving the name of the `colData(se)` column, in
  which the annotated segments should be stored. If `name` already
  exists in `colData(se)`, it will be overwritten.

## Value

A `SummarizedExperiment` object corresponding to `se`, with annotated
segments added to `colData(se)[[name]]`

## Examples

``` r
library(IRanges)
modbamfile <- system.file("extdata", "6mA_1_10reads.bam",
                          package = "SingleMoleculeGenomicsIO")
se <- readModBam(bamfiles = modbamfile, regions = "chr1:6935400-6936300",
                 modbase = "a", verbose = FALSE,
                 BPPARAM = BiocParallel::SerialParam())

segs <- list(s1 = IRangesList(
    "s1-d52a5f6a-a60a-4f85-913e-eada84bfbfb9" = IRanges(
        start = 6925834, width = 140),
    "s1-6cf74134-e550-4c02-bd2b-91385422ee25" = IRanges(
        start = c(6926000, 6926200), width = 140)))

se <- annotateReadSegments(se, segs, "nucl")
se$nucl
#> $s1
#> IRangesList object of length 10:
#> $`s1-233e48a7-f379-4dcf-9270-958231125563`
#> IRanges object with 0 ranges and 0 metadata columns:
#>        start       end     width
#>    <integer> <integer> <integer>
#> 
#> $`s1-d52a5f6a-a60a-4f85-913e-eada84bfbfb9`
#> IRanges object with 1 range and 0 metadata columns:
#>           start       end     width
#>       <integer> <integer> <integer>
#>   [1]   6925834   6925973       140
#> 
#> $`s1-fc4646ce-66f9-401f-b968-e9b0cda14d61`
#> IRanges object with 0 ranges and 0 metadata columns:
#>        start       end     width
#>    <integer> <integer> <integer>
#> 
#> ...
#> <7 more elements>
#> 
lengths(se$nucl$s1)
#> s1-233e48a7-f379-4dcf-9270-958231125563 s1-d52a5f6a-a60a-4f85-913e-eada84bfbfb9 
#>                                       0                                       1 
#> s1-fc4646ce-66f9-401f-b968-e9b0cda14d61 s1-92e906ae-cddb-4347-a114-bf9137761a8d 
#>                                       0                                       0 
#> s1-6cf74134-e550-4c02-bd2b-91385422ee25 s1-5d45d8d2-d5f5-47ff-a9fa-f3fd6b7bd3c7 
#>                                       2                                       0 
#> s1-b6fea9db-c92d-4152-9d29-4d021bbc45e8 s1-49c1e21e-8cb0-415a-aba9-92912219c4bb 
#>                                       0                                       0 
#> s1-b0b20f04-931f-4f60-b3e4-0ee1f5666a61 s1-41ca0e97-11b3-454b-9741-bc373e29ef37 
#>                                       0                                       0 
```
