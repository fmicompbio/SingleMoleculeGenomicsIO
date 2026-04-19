# SingleMoleculeGenomicsIO 0.1.4

* Bug fix in padding of extracted sequence context for regions on the minus strand

# SingleMoleculeGenomicsIO 0.1.3

* Return bam header information in the metadata of results from `readModBam`, `readMismatchBam`, `countStatePairs` and `countMismatchStatePairs`

# SingleMoleculeGenomicsIO 0.1.2

* Remove `removeAllNAReads` function (move to `footprintR`)

# SingleMoleculeGenomicsIO 0.1.1

* Add sampling mode to `countStatePairs` and `countMismatchStatePairs`

# SingleMoleculeGenomicsIO 0.1.0

* Add `readMismatchBam` and `countMismatchStatePairs` supporting BAM files
  from Bisulfite-seq or deaminase-based experiments.
* Initialize reading functions based on footprintR v0.4.0.
