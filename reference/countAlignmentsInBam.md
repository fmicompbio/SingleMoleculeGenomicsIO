# Count alignments in an indexed bam file

Report the number of mapped and unmapped alignments in a sorted and
indexed bam file.

## Usage

``` r
countAlignmentsInBam(bamfile)
```

## Arguments

- bamfile:

  Character scalar giving the name of the bam file to count alignments
  in. The file must have an associated index.

## Value

A named list with elements `"mapped"` and `"unmapped"`, where the
`"mapped"` element is a named vector with the number of alignments per
chromosome, and the `"unmapped"` element is a numeric scalar giving the
number of unmapped reads in the bam file.

## Examples

``` r
countAlignmentsInBam(system.file("extdata/6mA_1_10reads.bam",
                                 package = "SingleMoleculeGenomicsIO"))
#> $mapped
#>       chr1       chr2       chr3       chr4       chr5       chr6       chr7 
#>         10          0          0          0          0          0          0 
#>       chr8       chr9      chr10      chr11      chr12      chr13      chr14 
#>          0          0          0          0          0          0          0 
#>      chr15      chr16      chr17      chr18      chr19       chrX       chrY 
#>          0          0          0          0          0          0          0 
#>       chrM GL456210.1 GL456211.1 GL456212.1 GL456219.1 GL456221.1 GL456233.2 
#>          0          0          0          0          0          0          0 
#> GL456239.1 GL456354.1 GL456359.1 GL456360.1 GL456366.1 GL456367.1 GL456368.1 
#>          0          0          0          0          0          0          0 
#> GL456370.1 GL456372.1 GL456378.1 GL456379.1 GL456381.1 GL456382.1 GL456383.1 
#>          0          0          0          0          0          0          0 
#> GL456385.1 GL456387.1 GL456389.1 GL456390.1 GL456392.1 GL456394.1 GL456396.1 
#>          0          0          0          0          0          0          0 
#> JH584295.1 JH584296.1 JH584297.1 JH584298.1 JH584299.1 JH584300.1 JH584301.1 
#>          0          0          0          0          0          0          0 
#> JH584302.1 JH584303.1 JH584304.1 MU069434.1 MU069435.1 
#>          0          0          0          0          0 
#> 
#> $unmapped
#> [1] 0
#> 
```
