suppressPackageStartupMessages({
    library(GenomicRanges)
    library(Rsamtools)
    library(Biostrings)
})

## -------------------------------------------------------------------------- ##
## Checks, readMismatchBam
## -------------------------------------------------------------------------- ##
test_that("readMismatchBam works", {
    # example data
    ref <- system.file("extdata", "reference.fa.gz", package = "SingleMoleculeGenomicsIO")
    bamfiles <- system.file("extdata",
                            c("BisSeq_quasr_paired.bam",
                              "BisSeq_bismark_paired.bam",
                              "BisSeq_quasr_single.bam"),
                            package = "SingleMoleculeGenomicsIO")
    names(bamfiles) <- c("quasr", "bismark", "quasrsingle")
    sample_annot <- data.frame(sample = c("quasr", "bismark", "quasrsingle"),
                               group = c("group1", "group1", "group2"),
                               condition = c("cond2", "cond1", "cond3"))

    # invalid arguments
    expect_error(readMismatchBam(bamfiles = character(0), bamFormat = "QuasR",
                                 regions = "chr1:6940000-6955000",
                                 BPPARAM = BiocParallel::SerialParam()),
                 "length of .bamfiles. must be between 1 and Inf")
    expect_error(readMismatchBam(bamfiles = "error", bamFormat = "QuasR",
                                 regions = "chr1:6940000-6955000",
                                 nAlnsToSample = 0,
                                 BPPARAM = BiocParallel::SerialParam()),
                 "not all .bamfiles. exist")
    expect_error(readMismatchBam(bamfiles = bamfiles[1], bamFormat = "Bismark",
                                 regions = "chr1:6940000-6955000",
                                 sequenceReference = ref,
                                 BPPARAM = BiocParallel::SerialParam()),
                 "Invalid Bismark bam format")
    expect_error(readMismatchBam(bamfiles = bamfiles[2], bamFormat = "QuasR",
                                 regions = "chr1:6940000-6955000",
                                 sequenceReference = ref,
                                 BPPARAM = BiocParallel::SerialParam()),
                 "Invalid QuasR bam format")
    expect_error(readMismatchBam(bamfiles = bamfiles, bamFormat = "QuasR",
                                 regions = NULL, nAlnsToSample = 0,
                                 sequenceContext = "NNN",
                                 BPPARAM = BiocParallel::SerialParam()),
                 "The central base of .sequenceContext.")
    expect_error(readMismatchBam(bamfiles = bamfiles, bamFormat = "QuasR",
                                 regions = NULL, nAlnsToSample = 0,
                                 sequenceContext = "ACGT",
                                 BPPARAM = BiocParallel::SerialParam()),
                 "needs to have an odd number of characters")
    expect_error(readMismatchBam(bamfiles = bamfiles, bamFormat = "QuasR",
                                 regions = NULL, nAlnsToSample = 0,
                                 sequenceContext = "NCG",
                                 readBaseUnmod = "C", readBaseMod = "C",
                                 BPPARAM = BiocParallel::SerialParam()),
                 "must not include common bases")
    expect_error(readMismatchBam(bamfiles = bamfiles, bamFormat = "QuasR",
                                 regions = NULL, nAlnsToSample = 0,
                                 sequenceContext = "NCG",
                                 BPPARAM = BiocParallel::SerialParam()),
                 ".regions. must contain at least one genomic range if not in sampling mode")
    expect_error(readMismatchBam(bamfiles = stats::setNames(unname(bamfiles),
                                                            c("s1", "s1")),
                                 bamFormat = "QuasR", regions = "chr1:6940000-6955000",
                                 nAlnsToSample = 0,
                                 BPPARAM = BiocParallel::SerialParam()),
                 "are not unique")
    expect_error(readMismatchBam(bamfiles = bamfiles, bamFormat = "QuasR",
                                 regions = "error", nAlnsToSample = 0,
                                 sequenceContext = "NCG",
                                 sequenceReference = NULL,
                                 BPPARAM = BiocParallel::SerialParam()),
                 "must be either a")
    expect_error(readMismatchBam(bamfiles = bamfiles, bamFormat = "QuasR",
                                 regions = "error", nAlnsToSample = 0,
                                 sequenceContext = "NCG",
                                 sequenceReference = ref,
                                 BPPARAM = BiocParallel::SerialParam()),
                 "Failed to get bam iterator")
    expect_error(readMismatchBam(bamfiles = bamfiles, bamFormat = "QuasR",
                                 regions = "chr1:6940000-6955000",
                                 nAlnsToSample = -1, sequenceContext = "NCG",
                                 sequenceReference = ref,
                                 BPPARAM = BiocParallel::SerialParam()),
                 "must be between 0 and Inf")
    expect_error(readMismatchBam(bamfiles = bamfiles, bamFormat = "QuasR",
                                 regions = "chr1:6940000-6955000",
                                 nAlnsToSample = "error",
                                 BPPARAM = BiocParallel::SerialParam()),
                 "must be of class .numeric.")
    expect_error(readMismatchBam(bamfiles = bamfiles, bamFormat = "QuasR",
                                 regions = "chr1:6940000-6955000",
                                 nAlnsToSample = 0,
                                 level = 1,
                                 BPPARAM = BiocParallel::SerialParam()),
                 "must be of class .character.")
    expect_error(readMismatchBam(bamfiles = bamfiles, bamFormat = "QuasR",
                                 regions = as("chr1:6940000-6955000", "GRanges"),
                                 seqinfo = "error",
                                 BPPARAM = BiocParallel::SerialParam()),
                 ".seqinfo. must be .NULL., a .Seqinfo. object")
    expect_error(readMismatchBam(bamfiles = bamfiles, bamFormat = "QuasR",
                                 regions = "chr1:6940000-6955000",
                                 nAlnsToSample = 0,
                                 level = "error",
                                 BPPARAM = BiocParallel::SerialParam()),
                 "must be one of")
    expect_error(
        expect_warning(
            expect_warning(readMismatchBam(bamfiles = bamfiles, bamFormat = "QuasR",
                                           regions = "chr1:6940000-6955000",
                                           nAlnsToSample = 10,
                                           seqnamesToSampleFrom = "error",
                                           sequenceReference = ref,
                                           BPPARAM = BiocParallel::SerialParam()),
                           "Ignoring .regions."),
            "Ignoring unknown target name"),
        "Cannot sample 10 alignments from a total of 0")
    expect_error(readMismatchBam(bamfiles = bamfiles, bamFormat = "QuasR",
                                 regions = "chr1:6940000-6955000",
                                 nAlnsToSample = 0,
                                 seqnamesToSampleFrom = "chr1", seqinfo = "error",
                                 BPPARAM = BiocParallel::SerialParam()),
                 ".seqinfo. must be .NULL.")
    expect_error(readMismatchBam(bamfiles = bamfiles, bamFormat = "QuasR",
                                 regions = "chr1:6940000-6955000",
                                 BPPARAM = -1, sequenceReference = ref),
                 ".BPPARAM. must be of class .BiocParallelParam.")
    expect_error(readMismatchBam(bamfiles = bamfiles, bamFormat = "QuasR",
                                 regions = "chr1:6940000-6955000",
                                 sequenceReference = ref,
                                 BPPARAM = BiocParallel::SerialParam, trim = 1),
                 ".trim. must be of class .logical.")
    expect_error(readMismatchBam(bamfiles = bamfiles, bamFormat = "QuasR",
                                 regions = "chr1:6940000-6955000",
                                 sequenceReference = ref,
                                 sampleAnnot = sample_annot[1, ],
                                 BPPARAM = BiocParallel::SerialParam()),
                 "Annotation information missing")
    expect_error(readMismatchBam(
        bamfiles = bamfiles, bamFormat = "QuasR",
        regions = "chr1:6940000-6955000",
        sampleAnnot = sample_annot[, c("group", "condition")],
        BPPARAM = BiocParallel::SerialParam()),
        ".sampleAnnot. must have at least a column")
    expect_error(readMismatchBam(bamfiles = bamfiles, bamFormat = "QuasR",
                                 regions = "chr1:6940000-6955000",
                                 nAlnsToSample = 5,
                                 level = "summary",
                                 BPPARAM = BiocParallel::SerialParam()),
                 "Read sampling is not supported")

    # query without reads
    seL <- list(readMismatchBam(bamfiles = bamfiles[1], bamFormat = "QuasR",
                                regions = "chr1:1-10",
                                sequenceReference = ref,
                                level = "read", verbose = FALSE,
                                BPPARAM = BiocParallel::SerialParam()),
                readMismatchBam(bamfiles = bamfiles[1], bamFormat = "QuasR",
                                regions = "chr1:1-10",
                                sequenceReference = ref,
                                level = "quickread", verbose = FALSE,
                                maxCoverage = 10000,
                                BPPARAM = BiocParallel::SerialParam()),
                readMismatchBam(bamfiles = bamfiles[1], bamFormat = "QuasR",
                                regions = "chr1:1-10",
                                sequenceReference = ref,
                                level = "summary", verbose = FALSE,
                                BPPARAM = BiocParallel::SerialParam()))
    expect_true(all(vapply(seL, is, logical(1), class2 = "SummarizedExperiment")))
    expect_true(all(vapply(seL, function(x) identical(dim(x), c(0L, 1L)), logical(1))))

    # expected results
    reg1 <- c("chr1:6940000-6955000", "chr1:6929000-6929500")
    reg2 <- GRanges("chr1", IRanges(start = 6940000, end = 6955000))
    reg3 <- rep("chr1:6940000-6955000", 3)
    reg3b <- "chr1:6925000-6926000"
    reg4 <- "chr1:6940000-6955000"
    reg5 <- c("chr1:6937000-6937001", "chr1:6928000-6928001")
    suppressMessages({
        expect_message(
            se1 <- readMismatchBam(bamfiles = bamfiles[1], bamFormat = "QuasR",
                                   regions = reg1, level = "read", nAlnsToSample = 0,
                                   sequenceContext = "GCH", sequenceReference = ref,
                                   seqnamesToSampleFrom = "chr1", verbose = TRUE,
                                   BPPARAM = BiocParallel::SerialParam())
        )
    })
    se1sum <- readMismatchBam(bamfiles = bamfiles[1], bamFormat = "QuasR",
                              regions = reg1, level = "summary", nAlnsToSample = 0,
                              sequenceContext = "GCH", sequenceReference = ref,
                              seqnamesToSampleFrom = "chr1", verbose = FALSE,
                              maxCoverage = 10000,
                              BPPARAM = BiocParallel::SerialParam())
    se1quick <- readMismatchBam(bamfiles = bamfiles[1], bamFormat = "QuasR",
                                regions = reg1, level = "quickread", nAlnsToSample = 0,
                                sequenceContext = "GCH", sequenceReference = ref,
                                seqnamesToSampleFrom = "chr1", verbose = FALSE,
                                BPPARAM = BiocParallel::SerialParam())
    # ... same as se1, but repeat the sample twice
    se1dupl <- readMismatchBam(bamfiles = c(s1 = bamfiles[1], s2 = bamfiles[1]),
                               bamFormat = "QuasR",
                               regions = reg1, level = "read", nAlnsToSample = 0,
                               sequenceContext = "GCH", sequenceReference = ref,
                               seqnamesToSampleFrom = "chr1", verbose = FALSE,
                               BPPARAM = BiocParallel::SerialParam())
    # ... automatic names
    se2 <- readMismatchBam(bamfiles = unname(bamfiles[1]), bamFormat = "QuasR",
                           regions = reg2, sequenceReference = ref,
                           nAlnsToSample = 0, seqnamesToSampleFrom = "chr1",
                           BPPARAM = BiocParallel::MulticoreParam(workers = 2L),
                           verbose = FALSE)
    se2sum <- readMismatchBam(bamfiles = unname(bamfiles[1]), bamFormat = "QuasR",
                              regions = reg2, sequenceReference = ref,
                              level = "summary",
                              nAlnsToSample = 0, seqnamesToSampleFrom = "chr1",
                              BPPARAM = BiocParallel::MulticoreParam(workers = 2L),
                              verbose = FALSE)
    se2quick <- readMismatchBam(bamfiles = unname(bamfiles[1]), bamFormat = "QuasR",
                                regions = reg2, sequenceReference = ref,
                                level = "quickread",
                                nAlnsToSample = 0, seqnamesToSampleFrom = "chr1",
                                BPPARAM = BiocParallel::MulticoreParam(workers = 2L),
                                verbose = FALSE)
    # ... Bismark
    se3 <- readMismatchBam(bamfiles = bamfiles[2], bamFormat = "Bismark",
                           regions = reg3, sequenceReference = ref,
                           nAlnsToSample = 0, seqnamesToSampleFrom = "chr1",
                           BPPARAM = BiocParallel::SerialParam(),
                           verbose = FALSE)
    se3b <- readMismatchBam(bamfiles = bamfiles[2], bamFormat = "Bismark",
                            regions = reg3b, sequenceReference = ref,
                            sequenceContext = "C",
                            nAlnsToSample = 0, seqnamesToSampleFrom = "chr1",
                            BPPARAM = BiocParallel::SerialParam(),
                            verbose = FALSE)
    se3sum <- readMismatchBam(bamfiles = bamfiles[2], bamFormat = "Bismark",
                              regions = reg3, sequenceReference = ref,
                              level = "summary",
                              nAlnsToSample = 0, seqnamesToSampleFrom = "chr1",
                              BPPARAM = BiocParallel::SerialParam(),
                              verbose = FALSE)
    se3quick <- readMismatchBam(bamfiles = bamfiles[2], bamFormat = "Bismark",
                                regions = reg3, sequenceReference = ref,
                                level = "quickread",
                                nAlnsToSample = 0, seqnamesToSampleFrom = "chr1",
                                BPPARAM = BiocParallel::SerialParam(),
                                verbose = FALSE)
    # ... compare to manually read alignments
    se4a <- readMismatchBam(bamfiles = bamfiles[1], bamFormat = "QuasR",
                            regions = reg5[1], sequenceReference = ref,
                            nAlnsToSample = 0, seqnamesToSampleFrom = "chr1",
                            BPPARAM = BiocParallel::SerialParam(),
                            verbose = FALSE)
    se4b <- readMismatchBam(bamfiles = bamfiles[1], bamFormat = "QuasR",
                            regions = reg5[1:2], sequenceReference = ref,
                            nAlnsToSample = 0, seqnamesToSampleFrom = "chr1",
                            BPPARAM = BiocParallel::SerialParam(),
                            verbose = FALSE)
    aln4a <- Rsamtools::scanBam(file = bamfiles[1],
                                param = Rsamtools::ScanBamParam(
                                    what = "qname",
                                    which = GRanges(reg5[1])
                                ))
    aln4b <- Rsamtools::scanBam(file = bamfiles[1],
                                param = Rsamtools::ScanBamParam(
                                    what = "qname",
                                    which = GRanges(reg5[1:2])
                                ))
    # ... variant names, sample reads
    expect_warning(
        se5a  <- readMismatchBam(bamfiles = bamfiles[1], bamFormat = "QuasR",
                                 regions = NULL, sequenceReference = ref,
                                 nAlnsToSample = 5, seqnamesToSampleFrom = "chr1",
                                 variantPositions = GPos("chr1", pos = 63000000),
                                 BPPARAM = BiocParallel::MulticoreParam(2L, RNGseed = 55L),
                                 verbose = FALSE)
    )
    suppressMessages({
        expect_message(
            se5b  <- readMismatchBam(bamfiles = bamfiles[1], bamFormat = "QuasR",
                                     regions = NULL, sequenceReference = ref,
                                     nAlnsToSample = 5, seqnamesToSampleFrom = "chr1",
                                     BPPARAM = BiocParallel::SerialParam(RNGseed = 55L),
                                     verbose = TRUE)
        )
    })
    # ... with sample annotation
    se6  <- readMismatchBam(bamfiles = bamfiles[1], bamFormat = "QuasR",
                            regions = reg5[1:2], sequenceReference = ref,
                            sampleAnnot = sample_annot,
                            BPPARAM = BiocParallel::SerialParam(RNGseed = 55L),
                            verbose = FALSE)
    se6sum  <- readMismatchBam(bamfiles = bamfiles[1], bamFormat = "QuasR",
                               regions = reg5[1:2], sequenceReference = ref,
                               level = "summary",
                               sampleAnnot = sample_annot,
                               BPPARAM = BiocParallel::SerialParam(RNGseed = 55L),
                               verbose = FALSE)
    se6quick  <- readMismatchBam(bamfiles = bamfiles[1], bamFormat = "QuasR",
                                 regions = reg5[1:2], sequenceReference = ref,
                                 level = "quickread",
                                 sampleAnnot = sample_annot,
                                 BPPARAM = BiocParallel::SerialParam(RNGseed = 55L),
                                 verbose = FALSE)
    # ... trim = TRUE
    suppressMessages({
        expect_message(
            se7 <- readMismatchBam(bamfiles = bamfiles[1], bamFormat = "QuasR",
                                   regions = reg1, level = "read", nAlnsToSample = 0,
                                   sequenceContext = "GCH", sequenceReference = ref,
                                   seqnamesToSampleFrom = "chr1", verbose = TRUE,
                                   trim = TRUE,
                                   BPPARAM = BiocParallel::SerialParam())
        )
    })
    se7sum <- readMismatchBam(bamfiles = bamfiles[1], bamFormat = "QuasR",
                              regions = reg1, level = "summary", nAlnsToSample = 0,
                              sequenceContext = "GCH", sequenceReference = ref,
                              seqnamesToSampleFrom = "chr1", verbose = FALSE,
                              trim = TRUE,
                              BPPARAM = BiocParallel::SerialParam())
    se7quick <- readMismatchBam(bamfiles = bamfiles[1], bamFormat = "QuasR",
                                regions = reg1, level = "quickread", nAlnsToSample = 0,
                                sequenceContext = "GCH", sequenceReference = ref,
                                seqnamesToSampleFrom = "chr1", verbose = FALSE,
                                trim = TRUE,
                                BPPARAM = BiocParallel::SerialParam())
    # ... trim = TRUE, single-end
    suppressMessages({
        expect_message(
            se8 <- readMismatchBam(bamfiles = bamfiles[3], bamFormat = "QuasR",
                                   regions = reg1, level = "read", nAlnsToSample = 0,
                                   sequenceContext = "GCH", sequenceReference = ref,
                                   seqnamesToSampleFrom = "chr1", verbose = TRUE,
                                   trim = TRUE,
                                   BPPARAM = BiocParallel::SerialParam())
        )
    })
    se8sum <- readMismatchBam(bamfiles = bamfiles[3], bamFormat = "QuasR",
                              regions = reg1, level = "summary", nAlnsToSample = 0,
                              sequenceContext = "GCH", sequenceReference = ref,
                              seqnamesToSampleFrom = "chr1", verbose = FALSE,
                              trim = TRUE,
                              BPPARAM = BiocParallel::SerialParam())
    se8quick <- readMismatchBam(bamfiles = bamfiles[3], bamFormat = "QuasR",
                                regions = reg1, level = "quickread", nAlnsToSample = 0,
                                sequenceContext = "GCH", sequenceReference = ref,
                                seqnamesToSampleFrom = "chr1", verbose = FALSE,
                                trim = TRUE,
                                BPPARAM = BiocParallel::SerialParam())

    seL <- list(se1, se1dupl, se2, se3, se3b, se4a, se4b, se5a, se5b)
    seLsum <- list(se1sum, se2sum, se3sum)
    seLquick <- list(se1quick, se2quick, se3quick)

    # ... structure
    expected_coldata_names <- c("sample", "n_reads", "readInfo")
    expected_coldata_names_summary <- c("sample")
    expected_read_info_names <- c("qscore", "read_length", "aligned_length",
                                  "variant_label", "ref_strand",
                                  "aligned_fraction")
    for (se in c(seL, list(se6, se7, se8), seLquick,
                 list(se6quick, se7quick, se8quick))) {
        expect_s4_class(se, "RangedSummarizedExperiment")
        expect_s4_class(rowRanges(se), "GPos")
        expect_s4_class(colData(se)$readInfo, "SimpleList")
        res_se <- lapply(colData(se)$readInfo, function(df) {
            expect_s4_class(df, "DataFrame")
            expect_named(df, expected_read_info_names)
        })
        expect_identical(assayNames(se), "mod_prob")
        expect_s4_class(assay(se, "mod_prob"), "DFrame")
        expect_s4_class(assay(se, "mod_prob")[[1]], "NaMatrix")
    }
    for (se in c(seL, list(se6), seLquick, list(se6quick))) {
        # se7/8 (trim=TRUE) may have some reads filtered out
        expect_equal(vapply(assay(se, "mod_prob"), ncol, 0), se$n_reads,
                     ignore_attr = TRUE)
    }
    for (se in c(seL, list(se7, se8), seLquick, list(se7quick, se8quick))) {
        expect_identical(colnames(colData(se)), expected_coldata_names)
    }
    expect_identical(colnames(colData(se6)), c(expected_coldata_names,
                                               "group", "condition"))
    expect_identical(colnames(colData(se6quick)), c(expected_coldata_names,
                                                    "group", "condition"))
    ## ... ... summary-level
    for (se in c(seLsum, list(se6sum, se7sum, se8sum))) {
        expect_s4_class(se, "RangedSummarizedExperiment")
        expect_s4_class(rowRanges(se), "GPos")
        expect_identical(assayNames(se), c("Nmod", "Nvalid", "FracMod"))
        expect_type(assay(se, "Nmod"), "double")
    }
    for (se in c(seLsum, list(se7sum, se8sum))) {
        expect_identical(colnames(colData(se)), expected_coldata_names_summary)
    }
    expect_identical(colnames(colData(se6sum)), c(expected_coldata_names_summary,
                                                  "group", "condition"))

    expect_identical(colnames(se1), names(bamfiles)[1])
    expect_identical(colnames(se1sum), names(bamfiles)[1])
    expect_identical(colnames(se1quick), names(bamfiles)[1])
    expect_identical(colnames(se1dupl), c("s1.quasr", "s2.quasr"))
    expect_identical(colnames(se2), c("s1"))
    expect_identical(colnames(se2sum), c("s1"))
    expect_identical(colnames(se2quick), c("s1"))
    expect_identical(colnames(se3), names(bamfiles)[2])
    expect_identical(colnames(se3b), names(bamfiles)[2])
    expect_identical(colnames(se3sum), names(bamfiles)[2])
    expect_identical(colnames(se3quick), names(bamfiles)[2])
    expect_identical(colnames(se4a), names(bamfiles)[1])
    expect_identical(colnames(se4b), names(bamfiles)[1])
    expect_identical(colnames(se5a), names(bamfiles)[1])
    expect_identical(colnames(se5b), names(bamfiles)[1])
    expect_identical(colnames(se6), names(bamfiles)[1])
    expect_identical(colnames(se6sum), names(bamfiles)[1])
    expect_identical(colnames(se6quick), names(bamfiles)[1])
    expect_identical(colnames(se7), names(bamfiles)[1])
    expect_identical(colnames(se7sum), names(bamfiles)[1])
    expect_identical(colnames(se7quick), names(bamfiles)[1])
    expect_identical(colnames(se8), names(bamfiles)[3])
    expect_identical(colnames(se8sum), names(bamfiles)[3])
    expect_identical(colnames(se8quick), names(bamfiles)[3])

    # ... se1 identical to each of the samples in se1dupl
    expect_identical(nnavals(assay(se1, "mod_prob")[[1]]),
                     nnavals(assay(se1dupl[, 1], "mod_prob")[[1]]))
    expect_identical(nnawhich(assay(se1, "mod_prob")[[1]]),
                     nnawhich(assay(se1dupl[, 1], "mod_prob")[[1]]))
    expect_identical(nnavals(assay(se1, "mod_prob")[[1]]),
                     nnavals(assay(se1dupl[, 2], "mod_prob")[[1]]))
    expect_identical(nnawhich(assay(se1, "mod_prob")[[1]]),
                     nnawhich(assay(se1dupl[, 2], "mod_prob")[[1]]))

    # ... content se1
    expect_identical(unname(se1$n_reads), 182L)
    expect_identical(dim(assay(se1, "mod_prob")[[1]]), c(454L, 182L))
    expect_identical(dim(se1), c(454L, 1L))
    expect_identical(colnames(se1), names(bamfiles)[1])
    expect_identical(lapply(se1$readInfo, rownames),
                     lapply(assay(se1, "mod_prob"), colnames))
    expect_equal(lapply(se1$readInfo, "[[", "qscore")[[1]][
        match(c("quasr-NB501735:53:HNNGGBGX3:2:23108:11466:4470",
                "quasr-NB501735:53:HNNGGBGX3:1:12106:19631:19994",
                "quasr-NB501735:42:HKFLKBGX3:2:13304:21790:14430",
                "quasr-NB501735:42:HKFLKBGX3:4:12405:23760:5729"),
              rownames(se1$readInfo[[1]]))
    ], c(35.68, 34.491935484, 32.5078125, 34.26))
    expect_identical(lapply(se1$readInfo, "[[", "read_length")[[1]][
        match(c("quasr-NB501735:53:HNNGGBGX3:2:23108:11466:4470",
                "quasr-NB501735:53:HNNGGBGX3:1:12106:19631:19994",
                "quasr-NB501735:42:HKFLKBGX3:2:13304:21790:14430",
                "quasr-NB501735:42:HKFLKBGX3:4:12405:23760:5729"),
              rownames(se1$readInfo[[1]]))
    ], c(150L, 124L, 128L, 300L))
    expect_identical(lapply(se1$readInfo, "[[", "aligned_length")[[1]][
        match(c("quasr-NB501735:53:HNNGGBGX3:2:23108:11466:4470",
                "quasr-NB501735:53:HNNGGBGX3:1:12106:19631:19994",
                "quasr-NB501735:42:HKFLKBGX3:2:13304:21790:14430",
                "quasr-NB501735:42:HKFLKBGX3:4:12405:23760:5729"),
              rownames(se1$readInfo[[1]]))
    ], c(150L, 124L, 128L, 300L))
    expect_identical(lapply(se1$readInfo, "[[", "variant_label"),
                     lapply(setNames(se1$n_reads, colnames(se1)), function(n) rep(NA_character_, n)))
    # ... compare to se1sum
    expect_identical(rownames(se1), rownames(se1sum))
    se1tmp <- flattenReadLevelAssay(se1)
    expect_identical(assay(se1tmp, "Nmod"), assay(se1sum, "Nmod"))
    expect_identical(assay(se1tmp, "Nvalid"), assay(se1sum, "Nvalid"))
    expect_identical(assay(se1tmp, "FracMod"), assay(se1sum, "FracMod"))
    expect_identical(colData(se1)[, c("sample")],
                     colData(se1sum)[, c("sample")])
    expect_identical(rowRanges(se1), rowRanges(se1sum))
    # ... compare to se1quick
    # ... ... there is no guarantee that the reads have to be in the
    #         same order
    se1quickreordered <- subsetReads(se1quick, colnames(SummarizedExperiment::assay(se1, "mod_prob")[[1]]))
    expect_identical(se1, se1quickreordered)
    expect_identical(rownames(se1), rownames(se1quickreordered))
    expect_identical(rowRanges(se1), rowRanges(se1quickreordered))
    expect_identical(assay(se1, "mod_prob"),
                     assay(se1quickreordered, "mod_prob"))
    expect_identical(metadata(se1), metadata(se1quickreordered))
    expect_identical(colData(se1)[, c("sample", "n_reads")],
                     colData(se1quickreordered)[, c("sample", "n_reads")])
    expect_identical(rownames(colData(se1)$readInfo$quasr),
                     rownames(colData(se1quickreordered)$readInfo$quasr))

    # ... content se2
    expect_identical(unname(se2$n_reads), 175L)
    expect_identical(dim(assay(se2, "mod_prob")[[1]]), c(432L, 175L))
    expect_identical(dim(se2), c(432L, 1L))
    expect_identical(colnames(se2), "s1")
    expect_identical(lapply(se2$readInfo, rownames),
                     lapply(assay(se2, "mod_prob"), colnames))
    expect_equal(lapply(se2$readInfo, "[[", "qscore")[[1]][
        match(c("s1-NB501735:53:HNNGGBGX3:4:12601:24136:1750",
                "s1-NB501735:42:HKFLKBGX3:2:12305:18518:19215",
                "s1-NB501735:53:HNNGGBGX3:2:21205:13860:2188",
                "s1-NB501735:53:HNNGGBGX3:3:23509:16843:7045"),
              rownames(se2$readInfo[[1]]))
    ], c(30.2, 35.69, 35.0430631, 34.990991))
    expect_identical(lapply(se2$readInfo, "[[", "read_length")[[1]][
        match(c("s1-NB501735:53:HNNGGBGX3:4:12601:24136:1750",
                "s1-NB501735:42:HKFLKBGX3:2:12305:18518:19215",
                "s1-NB501735:53:HNNGGBGX3:2:21205:13860:2188",
                "s1-NB501735:53:HNNGGBGX3:3:23509:16843:7045"),
              rownames(se2$readInfo[[1]]))
    ], c(90L, 200L, 149L, 149L))
    expect_identical(lapply(se2$readInfo, "[[", "aligned_length")[[1]][
        match(c("s1-NB501735:53:HNNGGBGX3:4:12601:24136:1750",
                "s1-NB501735:42:HKFLKBGX3:2:12305:18518:19215",
                "s1-NB501735:53:HNNGGBGX3:2:21205:13860:2188",
                "s1-NB501735:53:HNNGGBGX3:3:23509:16843:7045"),
              rownames(se2$readInfo[[1]]))
    ], c(90L, 200L, 149L, 149L))
    expect_identical(lapply(se2$readInfo, "[[", "variant_label"),
                     lapply(setNames(se2$n_reads, colnames(se2)), function(n) rep(NA_character_, n)))
    # ... compare to se2sum
    expect_identical(rownames(se2), rownames(se2sum))
    se2tmp <- flattenReadLevelAssay(se2)
    expect_identical(assay(se2tmp, "Nmod"), assay(se2sum, "Nmod"))
    expect_identical(assay(se2tmp, "Nvalid"), assay(se2sum, "Nvalid"))
    expect_identical(assay(se2tmp, "FracMod"), assay(se2sum, "FracMod"))
    expect_identical(colData(se2)[, c("sample")],
                     colData(se2sum)[, c("sample")])
    expect_identical(rowRanges(se2), rowRanges(se2sum))
    # ... compare to se2quick
    se2quickreordered <- subsetReads(se2quick, colnames(SummarizedExperiment::assay(se2, "mod_prob")[[1]]))
    expect_identical(se2, se2quickreordered)

    # ... content se3
    expect_identical(unname(se3$n_reads), 171L)
    expect_identical(dim(assay(se3, "mod_prob")[[1]]), c(426L, 171L))
    expect_identical(dim(se3), c(426L, 1L))
    expect_identical(lapply(se3$readInfo, rownames),
                     lapply(assay(se3, "mod_prob"), colnames))
    expect_equal(lapply(se3$readInfo, "[[", "qscore")[[1]][
        match(c("bismark-NB501735:53:HNNGGBGX3:4:12601:24136:1750_1:N:0:GCCAAT",
                "bismark-NB501735:42:HKFLKBGX3:2:12305:18518:19215_1:N:0:GCCAAT",
                "bismark-NB501735:53:HNNGGBGX3:2:21205:13860:2188_1:N:0:GCCAAT",
                "bismark-NB501735:53:HNNGGBGX3:3:23509:16843:7045_1:N:0:GCCAAT"),
              rownames(se3$readInfo[[1]]))
    ], c(30.2, 35.69, 35.0430631, 34.990991))
    expect_identical(lapply(se3$readInfo, "[[", "read_length")[[1]][
        match(c("bismark-NB501735:53:HNNGGBGX3:4:12601:24136:1750_1:N:0:GCCAAT",
                "bismark-NB501735:42:HKFLKBGX3:2:12305:18518:19215_1:N:0:GCCAAT",
                "bismark-NB501735:53:HNNGGBGX3:2:21205:13860:2188_1:N:0:GCCAAT",
                "bismark-NB501735:53:HNNGGBGX3:3:23509:16843:7045_1:N:0:GCCAAT"),
              rownames(se3$readInfo[[1]]))
    ], c(90L, 200L, 149L, 149L))
    expect_identical(lapply(se3$readInfo, "[[", "aligned_length")[[1]][
        match(c("bismark-NB501735:53:HNNGGBGX3:4:12601:24136:1750_1:N:0:GCCAAT",
                "bismark-NB501735:42:HKFLKBGX3:2:12305:18518:19215_1:N:0:GCCAAT",
                "bismark-NB501735:53:HNNGGBGX3:2:21205:13860:2188_1:N:0:GCCAAT",
                "bismark-NB501735:53:HNNGGBGX3:3:23509:16843:7045_1:N:0:GCCAAT"),
              rownames(se3$readInfo[[1]]))
    ], c(90L, 200L, 149L, 149L))
    # ... compare to se3sum
    expect_identical(rownames(se3), rownames(se3sum))
    se3tmp <- flattenReadLevelAssay(se3)
    expect_identical(assay(se3tmp, "Nmod"), assay(se3sum, "Nmod"))
    expect_identical(assay(se3tmp, "Nvalid"), assay(se3sum, "Nvalid"))
    expect_identical(assay(se3tmp, "FracMod"), assay(se3sum, "FracMod"))
    expect_identical(colData(se3)[, c("sample")],
                     colData(se3sum)[, c("sample")])
    expect_identical(rowRanges(se3), rowRanges(se3sum))
    # ... compare to se3quick
    se3quickreordered <- subsetReads(se3quick, colnames(SummarizedExperiment::assay(se3, "mod_prob")[[1]]))
    expect_identical(se3, se3quickreordered)

    # ... content of se3b
    true_readlevel_bismark <- data.frame(
        pos = c(6925411, 6925411, 6925417, 6925417, 6925426, 6925426, 6925435,
                6925435, 6925860, 6925866, 6925872, 6925872, 6925875, 6925875,
                6925963, 6925964),
        strand = c("+", "+", "+", "+", "+", "+", "+", "+", "-", "-", "-", "-",
                   "-", "-", "+", "-"),
        read_id = c("NB501735:42:HKFLKBGX3:3:11605:15919:14063_1:N:0:GCCAAT",
                    "NB501735:53:HNNGGBGX3:1:13103:18758:18374_1:N:0:GCCAAT",
                    "NB501735:42:HKFLKBGX3:3:11605:15919:14063_1:N:0:GCCAAT",
                    "NB501735:53:HNNGGBGX3:1:13103:18758:18374_1:N:0:GCCAAT",
                    "NB501735:42:HKFLKBGX3:3:11605:15919:14063_1:N:0:GCCAAT",
                    "NB501735:53:HNNGGBGX3:1:13103:18758:18374_1:N:0:GCCAAT",
                    "NB501735:42:HKFLKBGX3:3:11605:15919:14063_1:N:0:GCCAAT",
                    "NB501735:53:HNNGGBGX3:1:13103:18758:18374_1:N:0:GCCAAT",
                    "NB501735:42:HKFLKBGX3:3:22505:15192:19556_1:N:0:GCCAAT",
                    "NB501735:42:HKFLKBGX3:3:22505:15192:19556_1:N:0:GCCAAT",
                    "NB501735:42:HKFLKBGX3:3:22505:15192:19556_1:N:0:GCCAAT",
                    "NB501735:53:HNNGGBGX3:2:12308:22855:7543_1:N:0:GCCAAT",
                    "NB501735:42:HKFLKBGX3:3:22505:15192:19556_1:N:0:GCCAAT",
                    "NB501735:53:HNNGGBGX3:2:12308:22855:7543_1:N:0:GCCAAT",
                    "NB501735:42:HKFLKBGX3:1:23103:22334:7953_1:N:0:GCCAAT",
                    "NB501735:42:HKFLKBGX3:3:22505:15192:19556_1:N:0:GCCAAT"),
        mod_prob = c(0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1)
    )
    for (i in seq_len(nrow(true_readlevel_bismark))) {
        ridx <- BiocGenerics::which(
            start(se3b) == true_readlevel_bismark[i, "pos"] &
                strand(se3b) == true_readlevel_bismark[i, "strand"])
        cidx <- BiocGenerics::which(
            colnames(assay(se3b, "mod_prob")[[1]]) ==
                paste0("bismark-", true_readlevel_bismark[i, "read_id"]))
        expect_length(ridx, 1L)
        expect_length(cidx, 1L)
        expect_identical(unname(assay(se3b, "mod_prob")[[1]][ridx, cidx]),
                         true_readlevel_bismark[i, "mod_prob"])
    }

    # ... content of se4a and se4b (se4a should be a subset of se4b)
    # ... ... check ground truth
    expect_identical(names(aln4a), names(aln4b)[1])
    expect_identical(aln4a[[1]], aln4b[[1]])
    expect_length(aln4a[[1]]$qname, 3L)
    expect_length(intersect(aln4a[[1]]$qname, aln4b[[2]]$qname), 0L)
    # ... ... check return values
    mp4a <- assay(se4a, "mod_prob")
    mp4b <- assay(se4b, "mod_prob")
    expect_true(all(paste0("quasr-", aln4a[[1]]$qname) %in% colnames(mp4a$quasr)))
    expect_true(all(paste0("quasr-", aln4b[[1]]$qname) %in% colnames(mp4b$quasr)))
    idx <- intersect(rownames(se4a), rownames(se4b))
    expect_identical(idx, rownames(se4a))
    expect_identical(mp4a[idx, "quasr"][, paste0("quasr-", aln4a[[1]]$qname)],
                     mp4b[idx, "quasr"][, paste0("quasr-", aln4b[[1]]$qname)])

    # ... content of se5a and se5b
    expect_identical(se5a, se5b)

    # ... content of se6
    mp6 <- assay(se6, "mod_prob")
    expect_true(all(paste0("quasr-", aln4b[[1]]$qname) %in% colnames(mp6$quasr)))
    idx <- rownames(se6)
    expect_identical(mp6[idx, "quasr"][, paste0("quasr-", aln4b[[1]]$qname)],
                     mp4b[idx, "quasr"][, paste0("quasr-", aln4b[[1]]$qname)])
    # ... compare to se6sum
    expect_identical(rownames(se6), rownames(se6sum))
    se6tmp <- flattenReadLevelAssay(se6)
    expect_identical(assay(se6tmp, "Nmod"), assay(se6sum, "Nmod"))
    expect_identical(assay(se6tmp, "Nvalid"), assay(se6sum, "Nvalid"))
    expect_identical(assay(se6tmp, "FracMod"), assay(se6sum, "FracMod"))
    expect_identical(colData(se6)[, c("sample")],
                     colData(se6sum)[, c("sample")])
    expect_identical(rowRanges(se6), rowRanges(se6sum))
    # ... compare to se6quick
    expect_identical(se6, se6quick)

    # ... content se7 (like se1, but with trim=TRUE)
    expect_identical(unname(se7$n_reads), 182L)
    expect_identical(dim(assay(se7, "mod_prob")[[1]]), c(432L, 181L))
    expect_identical(dim(se7), c(432L, 1L))
    modprob1 <- as.matrix(assay(se1, "mod_prob"))
    modprob7 <- as.matrix(assay(se7, "mod_prob"))
    shared_rows <- intersect(rownames(modprob7), rownames(modprob1))
    shared_cols <- intersect(colnames(modprob7), colnames(modprob1))
    expect_identical(modprob1[shared_rows, shared_cols],
                     modprob7[shared_rows, shared_cols])

    # ... content se8 (trim, single-end)
    expect_identical(unname(se8$n_reads), 180L)
    expect_identical(dim(assay(se8, "mod_prob")[[1]]), c(394L, 179L))
    expect_identical(dim(se8), c(394L, 1L))
    expect_identical(colnames(se8), names(bamfiles)[3])
    expect_identical(lapply(se8$readInfo, rownames),
                     lapply(assay(se8, "mod_prob"), colnames))
    expect_equal(lapply(se8$readInfo, "[[", "qscore")[[1]][
        match(c("quasrsingle-NB501735:53:HNNGGBGX3:2:23108:11466:4470",
                "quasrsingle-NB501735:53:HNNGGBGX3:1:12106:19631:19994",
                "quasrsingle-NB501735:42:HKFLKBGX3:2:13304:21790:14430",
                "quasrsingle-NB501735:42:HKFLKBGX3:4:12405:23760:5729"),
              rownames(se8$readInfo[[1]]))
    ], c(35.7333333, 34.5322581, 30.234375, 35.34000))
    expect_identical(lapply(se8$readInfo, "[[", "read_length")[[1]][
        match(c("quasrsingle-NB501735:53:HNNGGBGX3:2:23108:11466:4470",
                "quasrsingle-NB501735:53:HNNGGBGX3:1:12106:19631:19994",
                "quasrsingle-NB501735:42:HKFLKBGX3:2:13304:21790:14430",
                "quasrsingle-NB501735:42:HKFLKBGX3:4:12405:23760:5729"),
              rownames(se8$readInfo[[1]]))
    ], c(75L, 62L, 64L, 150L))
    expect_identical(lapply(se8$readInfo, "[[", "aligned_length")[[1]][
        match(c("quasrsingle-NB501735:53:HNNGGBGX3:2:23108:11466:4470",
                "quasrsingle-NB501735:53:HNNGGBGX3:1:12106:19631:19994",
                "quasrsingle-NB501735:42:HKFLKBGX3:2:13304:21790:14430",
                "quasrsingle-NB501735:42:HKFLKBGX3:4:12405:23760:5729"),
              rownames(se8$readInfo[[1]]))
    ], c(75L, 62L, 64L, 150L))
    expect_identical(lapply(se8$readInfo, "[[", "variant_label"),
                     lapply(setNames(se8$n_reads-1, colnames(se8)), function(n) rep(NA_character_, n)))
    # ... compare to se8sum
    expect_identical(rownames(se8), rownames(se8sum))
    se8tmp <- flattenReadLevelAssay(se8)
    expect_identical(assay(se8tmp, "Nmod"), assay(se8sum, "Nmod"))
    expect_identical(assay(se8tmp, "Nvalid"), assay(se8sum, "Nvalid"))
    expect_identical(assay(se8tmp, "FracMod"), assay(se8sum, "FracMod"))
    expect_identical(colData(se8)[, c("sample")],
                     colData(se8sum)[, c("sample")])
    expect_identical(rowRanges(se8), rowRanges(se8sum))
    # ... compare to se8quick
    # ... ... there is no guarantee that the reads have to be in the
    #         same order
    se8quickreordered <- subsetReads(se8quick, colnames(SummarizedExperiment::assay(se8, "mod_prob")[[1]]))
    expect_identical(se8, se8quickreordered)
})

test_that("readMismatchBam correctly labels reads", {
    # example data
    bamfile <- system.file("extdata", "BisSeq_quasr_single.bam",
                           package = "SingleMoleculeGenomicsIO")
    ref <- system.file("extdata", "reference.fa.gz", package = "SingleMoleculeGenomicsIO")

    # extract alignments
    alns <- Rsamtools::scanBam(
        file = bamfile,
        param = ScanBamParam(which = GRanges("chr1", IRanges(6937000, 6939000)),
                             what = c("qname", "rname", "strand", "pos",
                                      "qwidth", "cigar", "seq")))[[1]]

    # extract first 2 and last 3 positions of each alignment
    rnames <- rep(as.character(alns$rname), each = 5L)
    # ... this was calculated using: GenomicAlignments::cigarWidthAlongReferenceSpace(alns$cigar)
    alnwidth <- c(145L, 149L, 150L, 48L, 148L, 89L, 150L, 59L, 102L, 63L, 40L,
                  84L, 87L, 90L, 150L, 74L, 79L, 116L, 71L, 74L, 148L, 150L, 65L,
                  121L)
    rpos <- unlist(lapply(seq_along(alns$qname), function(i) {
        alns$pos[i] + c(0:1, alnwidth[i] - 3:1)
    }))
    # ... add three positions that don't overlap any read (one before, two after)
    rnames <- c("chr1", rnames, "chr2", "chr1")
    rpos <- c(6925829L, rpos, 6941630L, 6941639L)
    # ... convert to GPos
    varpos <- GPos(seqnames = rnames, pos = rpos,
                   names = c("miss1", rep(alns$qname, each = 5), "miss2", "miss3"))

    # calculate expected labels
    softmaskStart <- suppressWarnings(
        ifelse(grepl("^[0-9]+S", alns$cigar),
               as.integer(sub("^([0-9]+)S.+$", "\\1", alns$cigar)),
               0L))
    softmaskEnd <- suppressWarnings(
        ifelse(grepl("[0-9]+S$", alns$cigar),
               as.integer(sub("^.+?([0-9]+)S$", "\\1", alns$cigar)),
               0L))
    readLabelParts <- paste0(subseq(x = alns$seq, start = softmaskStart + 1L, width = 2L),
                             subseq(x = alns$seq, end = width(alns$seq) - softmaskEnd, width = 3L))

    # run readMismatchBam
    se <- readMismatchBam(bamfiles = bamfile, bamFormat = "QuasR",
                          regions = varpos, sequenceContext = "NCN",
                          sequenceReference = ref,
                          variantPositions = varpos,
                          BPPARAM = BiocParallel::SerialParam())
    varposToSortedIdx <- match(varpos, metadata(se)$variantPositions)

    # compare to expected labels
    extractLabelParts <- unlist(lapply(seq_along(alns$qname), function(i) {
        rid <- alns$qname[i]
        idx <- varposToSortedIdx[mcols(varpos)$names == rid]
        paste(
            strsplit(se$readInfo$s1$variant_label[match(paste0("s1-", rid),
                                                        rownames(se$readInfo$s1))],
                     "")[[1]][idx],
            collapse = "")
    }))
    expect_identical(sort(sort(varpos), ignore.strand = TRUE), metadata(se)$variantPositions)
    expect_identical(extractLabelParts, readLabelParts)
    expect_true(all(grepl("^-.*--$", se$readInfo$s1$variant_label))) # missed positions

    # positions that are known to be variable across reads
    varpos2 <- GPos(seqnames = "chr1", pos = c(6937731, 6937788, 6937843, 6937857,
                                               6937873, 6937931, 6937932, 6938070,
                                               6938109))

    se2 <- readMismatchBam(bamfiles = bamfile, bamFormat = "QuasR",
                           sequenceContext = "NCN", sequenceReference = ref,
                           regions = varpos2,
                           variantPositions = varpos2,
                           BPPARAM = BiocParallel::SerialParam())
    bases <- c("A", "C", "G", "T", "-")
    expCnt <- matrix(
        as.integer(c(0, 2, 0, 1, 4,
                     0, 1, 0, 0, 6,
                     2, 0, 1, 0, 4,
                     3, 0, 0, 0, 4,
                     1, 0, 2, 0, 4,
                     0, 0, 1, 0, 6,
                     0, 0, 1, 0, 6,
                     0, 0, 0, 0, 7,
                     1, 0, 0, 0, 6)),
        ncol = length(bases), byrow = TRUE, dimnames = list(NULL, bases))
    obsCnt <- do.call(rbind, lapply(seq.int(9), function(i) {
        f <- factor(
            unlist(lapply(se2$readInfo$s1$variant_label, substr, i, i)),
            levels = bases
        )
        unclass(table(f))
    }))
    expect_identical(obsCnt, expCnt)
})
