test_that("countStatePairs works", {
    modbamfile <- system.file("extdata", "6mA_1_10reads.bam",
                              package = "SingleMoleculeGenomicsIO")

    # fail with wrong arguments
    expect_error(countStatePairs(),
                 ".bamfile. is missing")
    expect_error(countStatePairs(bamfile = "error"),
                 ".bamfile. does not exist")
    expect_error(countStatePairs(bamfile = modbamfile),
                 ".modbase. is missing")
    expect_error(countStatePairs(bamfile = modbamfile, modbase = "X"),
                 "invalid .modbase.")
    expect_error(countStatePairs(bamfile = modbamfile, regions = 1L,
                                 modbase = "a"),
                 "must be of class .character.")
    expect_error(countStatePairs(bamfile = modbamfile, regions = character(0),
                                 modbase = "a"),
                 ".regions. must contain at least one genomic range")
    expect_error(countStatePairs(bamfile = modbamfile, regions = NULL,
                                 modbase = "a", nAlnsToSample = -1),
                 ".nAlnsToSample. must be between 0 and Inf")
    expect_error(countStatePairs(bamfile = modbamfile, regions = NULL,
                                 modbase = "a", nAlnsToSample = 5,
                                 seqnamesToSampleFrom = NULL),
                 "must not be .NULL.")
    expect_error(countStatePairs(bamfile = modbamfile, regions = NULL,
                                 modbase = "a", nAlnsToSample = 5,
                                 seqnamesToSampleFrom = "chr2"),
                 "Cannot sample 5 alignments from a total of 0")
    expect_error(countStatePairs(bamfile = modbamfile, regions = ".",
                                 modbase = "a", threshUnmod = -1),
                 "must be between 0 and 1")
    expect_error(countStatePairs(bamfile = modbamfile, regions = ".",
                                 modbase = "a", threshMod = -1),
                 "must be between 0 and 1")
    expect_error(countStatePairs(bamfile = modbamfile, regions = ".",
                                 modbase = "a", threshUnmod = 0.8, threshMod = 0.3),
                 "must be less than or equal to")
    expect_error(countStatePairs(bamfile = modbamfile, regions = ".",
                                 modbase = "a", windowSize = "error"),
                 "must be of class .numeric.")
    expect_error(countStatePairs(bamfile = modbamfile, regions = ".",
                                 modbase = "a", minMapQ = "error"),
                 "must be of class .numeric.")
    expect_error(countStatePairs(bamfile = modbamfile, regions = ".",
                                 modbase = "a", minAlignedLength = -1),
                 "must be between 0 and Inf")
    expect_error(countStatePairs(bamfile = modbamfile, regions = ".",
                                 modbase = "a", BPPARAM = "error"),
                 "must be of class .BiocParallelParam.")
    expect_error(countStatePairs(bamfile = modbamfile, regions = ".",
                                 modbase = "a", verbose = "error",
                                 BPPARAM = BiocParallel::SerialParam()),
                 "must be of class .logical.")

    # expected results
    res1 <- countStatePairs(bamfile = modbamfile,
                            regions = GenomicRanges::GRanges("chr1:1-100000000"),
                            modbase = "a", windowSize = 300,
                            BPPARAM = BiocParallel::SerialParam())
    res2 <- countStatePairs(bamfile = modbamfile, regions = ".",
                            modbase = "a", windowSize = 300,
                            BPPARAM = BiocParallel::SerialParam())
    res3 <- countStatePairs(bamfile = modbamfile, regions = ".",
                            modbase = "h", windowSize = 300,
                            BPPARAM = BiocParallel::SerialParam())

    # expected results (sampling)
    expect_warning(
        expect_warning(
            res4 <- countStatePairs(bamfile = modbamfile, regions = ".",
                                    modbase = "a", windowSize = 300,
                                    nAlnsToSample = 10,
                                    seqnamesToSampleFrom = c("chr1", "error"),
                                    BPPARAM = BiocParallel::SerialParam()),
            "Ignoring .regions."),
        "Ignoring unknown target name: error")
    set.seed(42L)
    suppressMessages(expect_message(
        res5 <- countStatePairs(bamfile = modbamfile, regions = NULL,
                                modbase = "a", windowSize = 300,
                                nAlnsToSample = 5,
                                seqnamesToSampleFrom = "chr1",
                                BPPARAM = BiocParallel::SerialParam(),
                                verbose = TRUE)
    ))

    resL <- list(res1, res2, res3, res4, res5)

    for (i in seq_along(resL)) {
        expect_s4_class(resL[[i]], "DataFrame")
        expect_identical(colnames(resL[[i]]), c("S", "unmod_unmod", "unmod_mod",
                                                "mod_unmod", "mod_mod"))
        expect_identical(dim(resL[[i]]), c(300L, 5L))
    }

    expect_identical(res1, res2)
    expect_identical(sum(as.matrix(res3[, -1])), 0)

    expect_identical(res1, res4)
    expect_identical(colSums(as.matrix(res5)),
                     c(S = 45150, unmod_unmod = 414021, unmod_mod = 41737,
                       mod_unmod = 40658, mod_mod = 10917))
})
