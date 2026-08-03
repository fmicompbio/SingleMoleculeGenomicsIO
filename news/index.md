# Changelog

## SingleMoleculeGenomicsIO 0.1.6

- Fix reporting of total number of processed modifications
- Avoid processing of prior modification data upon removing of
  soft-masked bases (may improve performance of readModBam several-fold)

## SingleMoleculeGenomicsIO 0.1.5

- Breaking change: Require specification of column names for read
  regrouping by to include the name of the nested `colData` column where
  they should be found
- Do not allow column names in `colData(se)` containing ‘:’
- Add `addReadLevelAssay` function to simplify addition of read-level
  assays to a `SummarizedExperiment` object
- Add `extractColDataColumns` and `getReadNamesBySample` utility
  functions

## SingleMoleculeGenomicsIO 0.1.4

- Bug fix in padding of extracted sequence context for regions on the
  minus strand

## SingleMoleculeGenomicsIO 0.1.3

- Return bam header information in the metadata of results from
  `readModBam`, `readMismatchBam`, `countStatePairs` and
  `countMismatchStatePairs`

## SingleMoleculeGenomicsIO 0.1.2

- Remove `removeAllNAReads` function (move to `footprintR`)

## SingleMoleculeGenomicsIO 0.1.1

- Add sampling mode to `countStatePairs` and `countMismatchStatePairs`

## SingleMoleculeGenomicsIO 0.1.0

- Add `readMismatchBam` and `countMismatchStatePairs` supporting BAM
  files from Bisulfite-seq or deaminase-based experiments.
- Initialize reading functions based on footprintR v0.4.0.
