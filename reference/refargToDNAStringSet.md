# Convert a supported reference sequence argument to a DNAStringSet

This function takes a supported argument value specifying the reference
sequence (e.g. genome) and returns it as a
[`DNAStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)
object. Currently supported argument values are:

- "BSgenome:":

  : A BSgenome object

- "DNAStringSet":

  : A DNAStringSet object

- "character":

  : A path to a fasta file

## Usage

``` r
refargToDNAStringSet(sequenceReference)
```

## Arguments

- sequenceReference:

  A supported value defining the reference sequence.

## Value

A
[`DNAStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)
object.

## Author

Michael Stadler, Charlotte Soneson

## Examples

``` r
ref <- refargToDNAStringSet(system.file("extdata", "reference.fa.gz",
                                        package = "SingleMoleculeGenomicsIO"))
ref
#> DNAStringSet object of length 1:
#>       width seq                                             names               
#> [1] 6957060 NNNNNNNNNNNNNNNNNNNNNN...GACAGTAATCCAGGAAAGGGGA chr1
```
