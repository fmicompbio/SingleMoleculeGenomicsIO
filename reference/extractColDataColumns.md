# Extract one or more (possibly nested) columns from the colData of a SummarizedExperiment object

Extract one or more (possibly nested) columns from the colData of a
SummarizedExperiment object

## Usage

``` r
extractColDataColumns(
  se,
  colNames,
  alignWith = ifelse(any(grepl(":", colNames)), "read", "sample")
)
```

## Arguments

- se:

  A
  [`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object.

- colNames:

  A character vector corresponding to the names of annotation
  (`colData`) columns to extract and consolidate into a `data.frame`.
  The names can be either columns of `colData(se)` itself, or columns in
  a nested (read-level) annotation column of `colData(se)`. In the
  latter case, the `colNames` should be of the form
  `outerColName:innerColName`, where `outerColName` is a column name in
  `colData(se)`, and `innerColName` is a column name in each element of
  `colData(se)[[outerColName]]`.

- alignWith:

  A character scalar indicating whether the rows in the returned
  `data.frame` correspond to reads (`"read"`) or samples (`"sample"`).
  If any element of `colNames` refers to a nested column, `alignWith`
  must be `"read"`.

## Value

A `data.frame` with columns `sample`, `read` (if `alignWith = "read"`)
and `make.names(colNames)`.

## Author

Michael Stadler, Charlotte Soneson

## Examples

``` r
extractfiles <- system.file("extdata",
                            c("modkit_extract_rc_6mA_1.tsv.gz",
                              "modkit_extract_rc_6mA_2.tsv.gz"),
                            package = "SingleMoleculeGenomicsIO")
se <- readModkitExtract(extractfiles, modbase = "a", filter = "modkit",
                        BPPARAM = BiocParallel::SerialParam())
df <- extractColDataColumns(se = se,
                            colNames = c("modbase", "readInfo:ref_strand"),
                            alignWith = "read")
head(df)
#>     sample                                    read readInfo.ref_strand modbase
#> s11     s1 s1-233e48a7-f379-4dcf-9270-958231125563                   -       a
#> s12     s1 s1-d52a5f6a-a60a-4f85-913e-eada84bfbfb9                   -       a
#> s13     s1 s1-fc4646ce-66f9-401f-b968-e9b0cda14d61                   +       a
#> s14     s1 s1-92e906ae-cddb-4347-a114-bf9137761a8d                   -       a
#> s15     s1 s1-6cf74134-e550-4c02-bd2b-91385422ee25                   +       a
#> s16     s1 s1-5d45d8d2-d5f5-47ff-a9fa-f3fd6b7bd3c7                   +       a
dim(df)
#> [1] 20  4
df <- extractColDataColumns(se = se,
                            colNames = "modbase",
                            alignWith = "sample")
head(df)
#>   sample modbase
#> 1     s1       a
#> 2     s2       a
dim(df)
#> [1] 2 2
```
