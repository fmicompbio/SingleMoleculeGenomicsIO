test_that("countMismatchStatePairs works", {
    ref <- system.file("extdata", "reference.fa.gz", package = "SingleMoleculeGenomicsIO")
    bamfiles <- system.file("extdata",
                            c("BisSeq_quasr_paired.bam",
                              "BisSeq_bismark_paired.bam",
                              "BisSeq_quasr_single.bam",
                              "BisSeq_quasr_paired_discordant.bam",
                              "BisSeq_quasr_single_indels.bam"),
                            package = "SingleMoleculeGenomicsIO")
    names(bamfiles) <- c("quasr", "bismark", "quasrsingle", "quasrdiscordant",
                         "quasrindels")

    # fails with incorrect arguments
    expect_error(countMismatchStatePairs(bamfile = "error",
                                         BPPARAM = BiocParallel::SerialParam()),
                 ".bamfile. does not exist")
    expect_error(countMismatchStatePairs(bamfile = bamfiles,
                                         BPPARAM = BiocParallel::SerialParam()),
                 ".bamfile. must have length 1")
    expect_error(countMismatchStatePairs(bamfile = 1,
                                         BPPARAM = BiocParallel::SerialParam()),
                 ".bamfile. must be of class .character.")
    expect_error(countMismatchStatePairs(bamfile = bamfiles[1],
                                         bamFormat = "error",
                                         BPPARAM = BiocParallel::SerialParam()),
                 ".bamFormat. must be one of")
    expect_error(countMismatchStatePairs(bamfile = bamfiles[1],
                                         bamFormat = "Bismark",
                                         regions = "chr1:6940000-6955000",
                                         sequenceReference = ref,
                                         BPPARAM = BiocParallel::SerialParam()),
                 "Invalid Bismark bam format")
    expect_error(countMismatchStatePairs(bamfile = bamfiles[2],
                                         bamFormat = "QuasR",
                                         regions = "chr1:6940000-6955000",
                                         sequenceReference = ref,
                                         BPPARAM = BiocParallel::SerialParam()),
                 "Invalid QuasR bam format")
    expect_error(countMismatchStatePairs(bamfile = bamfiles[1],
                                         bamFormat = "QuasR",
                                         sequenceContext = c("ACT", "GTA"),
                                         BPPARAM = BiocParallel::SerialParam()),
                 ".sequenceContext. must have length 1")
    expect_error(countMismatchStatePairs(bamfile = bamfiles[1],
                                         bamFormat = "QuasR",
                                         sequenceContext = "NNN",
                                         BPPARAM = BiocParallel::SerialParam()),
                 "The central base of .sequenceContext.")
    expect_error(countMismatchStatePairs(bamfile = bamfiles[1],
                                         bamFormat = "QuasR",
                                         sequenceContext = "NN",
                                         BPPARAM = BiocParallel::SerialParam()),
                 "needs to have an odd number of characters")
    expect_error(countMismatchStatePairs(bamfile = bamfiles[1],
                                         bamFormat = "QuasR",
                                         nAlnsToSample = -1,
                                         seqinfo = c(chr1 = 1000000),
                                         sequenceReference = ref,
                                         BPPARAM = BiocParallel::SerialParam()),
                 ".nAlnsToSample. must be between 0 and Inf")
    expect_error(
        expect_warning(countMismatchStatePairs(bamfile = bamfiles[1],
                                               bamFormat = "QuasR",
                                               nAlnsToSample = 1000,
                                               seqnamesToSampleFrom = "chr1",
                                               seqinfo = c(chr1 = 1000000),
                                               sequenceReference = ref,
                                               BPPARAM = BiocParallel::SerialParam()),
                       "Ignoring .regions. because .nAlnsToSample. is greater than zero"),
                 "Cannot sample 1000 alignments from a total of ")
    expect_error(
        expect_warning(countMismatchStatePairs(bamfile = bamfiles[1],
                                               regions = NULL,
                                               bamFormat = "QuasR",
                                               nAlnsToSample = 1000,
                                               seqnamesToSampleFrom = "error",
                                               seqinfo = c(chr1 = 1000000),
                                               sequenceReference = ref,
                                               BPPARAM = BiocParallel::SerialParam()),
                       "Ignoring unknown target name: error"),
        "Cannot sample 1000 alignments from a total of 0")
    expect_error(countMismatchStatePairs(bamfile = bamfiles[1],
                                         bamFormat = "QuasR",
                                         readBaseUnmod = "N",
                                         BPPARAM = BiocParallel::SerialParam()),
                 "All values in .readBaseUnmod. must be one of")
    expect_error(countMismatchStatePairs(bamfile = bamfiles[1],
                                         bamFormat = "QuasR",
                                         readBaseMod = "N",
                                         BPPARAM = BiocParallel::SerialParam()),
                 "All values in .readBaseMod. must be one of")
    expect_error(countMismatchStatePairs(bamfile = bamfiles[1],
                                         bamFormat = "QuasR",
                                         readBaseUnmod = "C",
                                         readBaseMod = "C",
                                         BPPARAM = BiocParallel::SerialParam()),
                 "must not include common bases")
    expect_error(countMismatchStatePairs(bamfile = bamfiles[1],
                                         bamFormat = "QuasR",
                                         windowSize = "100",
                                         seqinfo = c(chr1 = 1000000),
                                         sequenceReference = ref,
                                         BPPARAM = BiocParallel::SerialParam()),
                 ".windowSize. must be of class .numeric.")
    expect_error(countMismatchStatePairs(bamfile = bamfiles[1],
                                         bamFormat = "QuasR",
                                         windowSize = -1,
                                         seqinfo = c(chr1 = 1000000),
                                         sequenceReference = ref,
                                         BPPARAM = BiocParallel::SerialParam()),
                 ".windowSize. must be between 1 and Inf")
    expect_error(countMismatchStatePairs(bamfile = bamfiles[1],
                                         bamFormat = "QuasR",
                                         minMapQ = "100",
                                         seqinfo = c(chr1 = 1000000),
                                         sequenceReference = ref,
                                         BPPARAM = BiocParallel::SerialParam()),
                 ".minMapQ. must be of class .numeric.")
    expect_error(countMismatchStatePairs(bamfile = bamfiles[1],
                                         bamFormat = "QuasR",
                                         minMapQ = c(1, 20),
                                         seqinfo = c(chr1 = 1000000),
                                         sequenceReference = ref,
                                         BPPARAM = BiocParallel::SerialParam()),
                 ".minMapQ. must have length 1")
    expect_error(countMismatchStatePairs(bamfile = bamfiles[1],
                                         bamFormat = "QuasR",
                                         minAlignedLength = "100",
                                         seqinfo = c(chr1 = 1000000),
                                         sequenceReference = ref,
                                         BPPARAM = BiocParallel::SerialParam()),
                 ".minAlignedLength. must be of class .numeric.")
    expect_error(countMismatchStatePairs(bamfile = bamfiles[1],
                                         bamFormat = "QuasR",
                                         minAlignedLength = c(1, 20),
                                         seqinfo = c(chr1 = 1000000),
                                         sequenceReference = ref,
                                         BPPARAM = BiocParallel::SerialParam()),
                 ".minAlignedLength. must have length 1")
    expect_error(countMismatchStatePairs(bamfile = bamfiles[1],
                                         bamFormat = "QuasR",
                                         regions = "error",
                                         BPPARAM = BiocParallel::SerialParam()),
                 ".sequenceReference. must be either a")
    expect_error(countMismatchStatePairs(bamfile = bamfiles[1],
                                         bamFormat = "QuasR",
                                         regions = "error",
                                         sequenceContext = "NCG",
                                         sequenceReference = ref,
                                         BPPARAM = BiocParallel::SerialParam()),
                 "Failed to get bam iterator")
    expect_error(countMismatchStatePairs(bamfile = bamfiles[1],
                                         bamFormat = "QuasR",
                                         regions = GenomicRanges::GRanges(),
                                         sequenceContext = "NCG",
                                         sequenceReference = ref,
                                         BPPARAM = BiocParallel::SerialParam()),
                 ".regions. must contain at least one genomic range")
    expect_error(countMismatchStatePairs(bamfile = bamfiles[1],
                                         bamFormat = "QuasR",
                                         regions = as("chr1:6940000-6955000", "GRanges"),
                                         seqinfo = "error",
                                         BPPARAM = BiocParallel::SerialParam()),
                 ".seqinfo. must be .NULL., a .Seqinfo. object")
    expect_error(countMismatchStatePairs(bamfile = bamfiles[1],
                                         bamFormat = "QuasR",
                                         regions = "chr1:6940000-6955000",
                                         seqinfo = "error",
                                         BPPARAM = BiocParallel::SerialParam()),
                 ".seqinfo. must be .NULL.")
    expect_error(countMismatchStatePairs(bamfile = bamfiles[1],
                                         bamFormat = "QuasR",
                                         regions = "chr1:6940000-6955000",
                                         BPPARAM = -1, sequenceReference = ref),
                 ".BPPARAM. must be of class .BiocParallelParam.")

    # works with correct arguments
    # ... single-end
    res <- countMismatchStatePairs(bamfile = bamfiles[3],
                                   bamFormat = "QuasR",
                                   regions = "chr1:6940000-6955000",
                                   sequenceContext = "GCH",
                                   windowSize = 200,
                                   sequenceReference = ref,
                                   BPPARAM = BiocParallel::SerialParam())
    expect_s4_class(res, "DFrame")
    expect_identical(dim(res), c(200L, 5L))
    expect_named(res, c("S", "unmod_unmod", "unmod_mod", "mod_unmod", "mod_mod"))
    ## compare to e.g. nomeR::get_ctable_from_SE(readMismatchBam(bamfiles = bamfiles[3], regions = "chr1:6940000-6955000", sequenceReference = ref, sequenceContext = "GCH"), "mod_prob", 0.5, 0.5, 0, 0, 200, FALSE, 1, FALSE)
    expect_identical(unname(unlist(res[1, ])), c(1, 513, 0, 0, 2))
    expect_identical(res[2:12, 2], c(0, 0, 17, 11, 10, 21, 19, 10, 11, 16, 9))
    expect_identical(unname(colSums(as.data.frame(res))), c(20100, 1292, 2, 2, 2))

    # ... single-end, smaller window size
    res <- countMismatchStatePairs(bamfile = bamfiles[3],
                                   bamFormat = "QuasR",
                                   regions = "chr1:6940000-6955000",
                                   sequenceContext = "GCH",
                                   windowSize = 20,
                                   sequenceReference = ref,
                                   BPPARAM = BiocParallel::SerialParam())
    expect_s4_class(res, "DFrame")
    expect_identical(dim(res), c(20L, 5L))
    expect_named(res, c("S", "unmod_unmod", "unmod_mod", "mod_unmod", "mod_mod"))
    ## compare to e.g. nomeR::get_ctable_from_SE(readMismatchBam(bamfiles = bamfiles[3], regions = "chr1:6940000-6955000", sequenceReference = ref, sequenceContext = "GCH"), "mod_prob", 0.5, 0.5, 0, 0, 20, FALSE, 1, FALSE)
    expect_identical(unname(unlist(res[1, ])), c(1, 513, 0, 0, 2))
    expect_identical(res[2:12, 2], c(0, 0, 17, 11, 10, 21, 19, 10, 11, 16, 9))
    expect_identical(unname(colSums(as.data.frame(res))), c(210, 733, 2, 0, 2))

    # ... paired-end, QuasR, one unpaired bam record
    res <- countMismatchStatePairs(bamfile = bamfiles[1],
                                   bamFormat = "QuasR",
                                   regions = "chr1:6925411-6925964",
                                   sequenceContext = "GCH",
                                   windowSize = 200,
                                   sequenceReference = ref,
                                   BPPARAM = BiocParallel::SerialParam())
    expect_s4_class(res, "DFrame")
    expect_identical(dim(res), c(200L, 5L))
    expect_named(res, c("S", "unmod_unmod", "unmod_mod", "mod_unmod", "mod_mod"))
    ## compare to e.g. nomeR::get_ctable_from_SE(readMismatchBam(bamfiles = bamfiles[1], regions = "chr1:6925411-6925964", sequenceReference = ref, sequenceContext = "GCH"), "mod_prob", 0.5, 0.5, 0, 0, 200, FALSE, 1, FALSE)
    expect_identical(unname(unlist(res[1, ])), c(1, 37, 0, 0, 1))
    expect_identical(res[2:12, 2], c(0, 0, 6, 2, 0, 5, 4, 2, 5, 0, 0))
    expect_identical(res[179:185, 2], c(1, 0, 0, 0, 0, 0, 1))
    expect_identical(unname(colSums(as.data.frame(res))), c(20100, 231, 1, 0, 1))

    # ... paired-end, Bismark
    res <- countMismatchStatePairs(bamfile = bamfiles[2],
                                   bamFormat = "Bismark",
                                   regions = "chr1:6925411-6925964",
                                   sequenceContext = "GCH",
                                   windowSize = 200,
                                   sequenceReference = ref,
                                   BPPARAM = BiocParallel::SerialParam())
    expect_s4_class(res, "DFrame")
    expect_identical(dim(res), c(200L, 5L))
    expect_named(res, c("S", "unmod_unmod", "unmod_mod", "mod_unmod", "mod_mod"))
    ## compare to e.g. nomeR::get_ctable_from_SE(readMismatchBam(bamfiles = bamfiles[2], bamFormat = "Bismark", regions = "chr1:6925411-6925964", sequenceReference = ref, sequenceContext = "GCH"), "mod_prob", 0.5, 0.5, 0, 0, 200, FALSE, 1, FALSE)
    expect_identical(unname(unlist(res[1, ])), c(1, 37, 0, 0, 1))
    expect_identical(res[2:12, 2], c(0, 0, 6, 2, 0, 5, 4, 2, 5, 0, 0))
    expect_identical(res[179:185, 2], c(1, 0, 0, 0, 0, 0, 1))
    expect_identical(unname(colSums(as.data.frame(res))), c(20100, 231, 1, 0, 1))

    # ... read with mismatch between mates (chr1:6925435)
    res <- countMismatchStatePairs(bamfile = bamfiles[4],
                                   bamFormat = "QuasR",
                                   regions = "chr1:6925411-6925964",
                                   sequenceContext = "C",
                                   windowSize = 200,
                                   sequenceReference = ref,
                                   BPPARAM = BiocParallel::SerialParam())
    expect_s4_class(res, "DFrame")
    expect_identical(dim(res), c(200L, 5L))
    expect_named(res, c("S", "unmod_unmod", "unmod_mod", "mod_unmod", "mod_mod"))
    ## compare to e.g. nomeR::get_ctable_from_SE(readMismatchBam(bamfiles = bamfiles[4], bamFormat = "QuasR", regions = "chr1:6925411-6925964", sequenceReference = ref, sequenceContext = "C"), "mod_prob", 0.5, 0.5, 0, 0, 200, FALSE, 1, FALSE)
    expect_identical(unname(unlist(res[1, ])), c(1, 12, 0, 0, 4))
    expect_identical(res[2:12, 2], c(3, 3, 1, 3, 1, 4, 2, 1, 1, 2, 3))
    expect_identical(res[33, 3], 2)
    expect_identical(unname(colSums(as.data.frame(res))), c(20100, 78, 32, 16, 10))

    # ... reads with indels
    res1 <- countMismatchStatePairs(bamfile = bamfiles[5],
                                   bamFormat = "QuasR",
                                   regions = "chr1:1-7000000",
                                   sequenceContext = "C",
                                   windowSize = 200,
                                   sequenceReference = ref,
                                   BPPARAM = BiocParallel::SerialParam())
    expect_s4_class(res1, "DFrame")
    expect_identical(dim(res1), c(200L, 5L))
    expect_named(res1, c("S", "unmod_unmod", "unmod_mod", "mod_unmod", "mod_mod"))
    ## compare to e.g. nomeR::get_ctable_from_SE(readMismatchBam(bamfiles = bamfiles[5], bamFormat = "QuasR", regions = "chr1:6925411-6925964", sequenceReference = ref, sequenceContext = "C"), "mod_prob", 0.5, 0.5, 0, 0, 200, FALSE, 1, FALSE)
    expect_identical(unname(unlist(res1[1, ])), c(1, 68, 0, 0, 0))
    expect_identical(res1[2:12, 2], c(21, 17, 9, 15, 9, 18, 14, 8, 14, 15, 20))
    expect_identical(unname(colSums(as.data.frame(res1))), c(20100, 805, 0, 0, 0))

    # ... with sampling
    # ... ... single-end QuasR (complete)
    res2 <- countMismatchStatePairs(bamfile = bamfiles[5],
                                    bamFormat = "QuasR",
                                    regions = NULL,
                                    sequenceContext = "C",
                                    nAlnsToSample = 3,
                                    seqnamesToSampleFrom = "chr1",
                                    windowSize = 200,
                                    sequenceReference = ref,
                                    BPPARAM = BiocParallel::SerialParam())
    expect_identical(res1, res2)
    # ... ... paired-end Bismark (complete)
    res3 <-  countMismatchStatePairs(bamfile = bamfiles[2],
                                     bamFormat = "Bismark",
                                     regions = "chr1",
                                     sequenceContext = "GCH",
                                     windowSize = 200,
                                     sequenceReference = ref,
                                     BPPARAM = BiocParallel::SerialParam())
    suppressMessages(expect_message(
        res4 <-  countMismatchStatePairs(bamfile = bamfiles[2],
                                         bamFormat = "Bismark",
                                         regions = NULL,
                                         sequenceContext = "GCH",
                                         nAlnsToSample = 788,
                                         seqnamesToSampleFrom = "chr1",
                                         windowSize = 200,
                                         sequenceReference = ref,
                                         BPPARAM = BiocParallel::SerialParam(),
                                         verbose = TRUE)
    ))
    expect_identical(res3, res4)
    # ... ... paired-end Bismark (subsample)
    set.seed(43L)
    res5 <-  countMismatchStatePairs(bamfile = bamfiles[2],
                                     bamFormat = "Bismark",
                                     regions = NULL,
                                     sequenceContext = "GCH",
                                     nAlnsToSample = 216,
                                     seqnamesToSampleFrom = "chr1",
                                     windowSize = 200,
                                     sequenceReference = ref,
                                     BPPARAM = BiocParallel::SerialParam())
    expect_true(all(as.vector(as.matrix(res4[, -1]) - as.matrix(res5[, -1])) >= 0))
    expect_true(any(as.vector(as.matrix(res4[, -1]) - as.matrix(res5[, -1])) > 0))
})
