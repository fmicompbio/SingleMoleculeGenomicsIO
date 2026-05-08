test_that("addReadLevelAssay works", {
    modbamfile <- system.file("extdata", c("6mA_1_10reads.bam", "6mA_1_10reads.bam"),
                              package = "SingleMoleculeGenomicsIO")
    se <- readModBam(bamfiles = modbamfile, regions = "chr1:6940000-6955000",
                     modbase = "a", verbose = FALSE,
                     BPPARAM = BiocParallel::SerialParam())
    adf <- SummarizedExperiment::assay(se, "mod_prob")

    expect_error(addReadLevelAssay(se = 1, assayDF = adf, assayName = "new_assay"),
                 ".se. must be of class .SummarizedExperiment.")
    expect_error(addReadLevelAssay(se = se, assayDF = adf[[1]], assayName = "new_assay"),
                 ".assayDF. must be of class .DataFrame.")
    expect_error(addReadLevelAssay(se = se, assayDF = cbind(adf, adf), assayName = "new_assay"),
                 "The dimensions of .se.")
    expect_error(addReadLevelAssay(se = se, assayDF = adf, assayName = 1),
                 ".assayName. must be of class .character.")
    expect_error(addReadLevelAssay(se = se, assayDF = adf, assayName = "mod_prob",
                                   replaceExisting = FALSE),
                 ".se. already has an assay named mod_prob")
    expect_error(addReadLevelAssay(se = se, assayDF = adf, assayName = "new_assay",
                                   replaceExisting = 1),
                 ".replaceExisting. must be of class .logical.")
    expect_error(addReadLevelAssay(se = se, assayDF = adf, assayName = "new_assay",
                                   replaceExisting = FALSE, verbose = 1),
                 ".verbose. must be of class .logical.")

    # add read-level assay
    se2 <- addReadLevelAssay(se = se, assayDF = adf, assayName = "new_assay",
                             replaceExisting = FALSE, verbose = FALSE)
    expect_s4_class(se2, "RangedSummarizedExperiment")
    expect_identical(dim(se), dim(se2))
    expect_identical(SummarizedExperiment::assayNames(se2),
                     c("mod_prob", "new_assay"))
    expect_identical(S4Vectors::metadata(se2)$readLevelData$assayNames,
                     c("mod_prob", "new_assay"))

    # flatten based on each read-level assay
    se1a <- flattenReadLevelAssay(se, assayName = "mod_prob", keepReads = FALSE)
    se1b <- flattenReadLevelAssay(se2, assayName = "new_assay", keepReads = FALSE)
    expect_identical(se1a, se1b)
})
