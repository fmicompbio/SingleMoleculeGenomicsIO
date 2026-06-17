test_that("expandSEToBaseSpace works", {
    modbamfiles <- system.file("extdata",
                               c("6mA_1_10reads.bam", "6mA_2_10reads.bam"),
                               package = "SingleMoleculeGenomicsIO")
    names(modbamfiles) <- c("sample1", "sample2")
    se <- readModBam(bamfiles = modbamfiles,
                     regions = "chr1:6940000-6941000",
                     modbase = "a",
                     trim = TRUE,
                     BPPARAM = BiocParallel::SerialParam())
    se <- flattenReadLevelAssay(se)

    # expected fails
    expect_error(expandSEToBaseSpace(se = se, region = "chr1:6940100-6941000"),
                 "Not all positions in .se. are within .region.")
    expect_error(expandSEToBaseSpace(se = se, keepAssays = "doesnt-exist"),
                 "All values in .keepAssays. must be one of: mod_prob")
    expect_error(expandSEToBaseSpace(se = se, ignore.strand = FALSE),
                 "All values in .ignore.strand. must be one of: TRUE")

    # empty input
    seExp0 <- expandSEToBaseSpace(se = se[numeric(0), ])
    expect_identical(dim(seExp0), c(0L, ncol(se)))
    expect_identical(SummarizedExperiment::assayNames(seExp0),
                     getReadLevelAssayNames(se))
    expect_identical(lapply(SummarizedExperiment::assay(se, "mod_prob"), colnames),
                     lapply(SummarizedExperiment::assay(seExp0, "mod_prob"), colnames))

    # expected return value
    seExp1 <- expandSEToBaseSpace(se = se)
    expect_s4_class(seExp1, "RangedSummarizedExperiment")
    expect_identical(dim(seExp1),
                     c(width(range(SummarizedExperiment::rowRanges(se),
                                   ignore.strand = TRUE)),
                       ncol(se)))
    expect_identical(SummarizedExperiment::assayNames(seExp1),
                     getReadLevelAssayNames(se))
    expect_identical(lapply(SummarizedExperiment::assay(se, "mod_prob"), colnames),
                     lapply(SummarizedExperiment::assay(seExp1, "mod_prob"), colnames))
    for (s in colnames(se)) {
        assayOrig <- SummarizedExperiment::assay(se, "mod_prob")[[s]]
        assayExp <- SummarizedExperiment::assay(seExp1, "mod_prob")[[s]]

        valsOrig <- SparseArray::nnavals(assayOrig)
        valsExp <- SparseArray::nnavals(assayExp)
        expect_identical(valsOrig, valsExp)

        posOrig <- GenomicRanges::start(se)[SparseArray::nnawhich(assayOrig)]
        posExp <- GenomicRanges::start(seExp1)[SparseArray::nnawhich(assayExp)]
        expect_identical(posOrig, posExp)
    }
    expect_identical(as.character(unique(GenomicRanges::strand(seExp1))), "*")

    reg2 <- GenomicRanges::GRanges("chr1",
                                   IRanges::IRanges(6939000, 6942000))
    seExp2 <- expandSEToBaseSpace(se = se, region = reg2)
    expect_s4_class(seExp2, "RangedSummarizedExperiment")
    expect_identical(dim(seExp2), c(width(reg2), ncol(se)))
    expect_identical(SummarizedExperiment::assayNames(seExp2),
                     getReadLevelAssayNames(se))
    expect_identical(lapply(assay(se, "mod_prob"), colnames),
                     lapply(assay(seExp2, "mod_prob"), colnames))
    for (s in colnames(se)) {
        assayOrig <- assay(se, "mod_prob")[[s]]
        assayExp <- assay(seExp2, "mod_prob")[[s]]

        valsOrig <- SparseArray::nnavals(assayOrig)
        valsExp <- SparseArray::nnavals(assayExp)
        expect_identical(valsOrig, valsExp)

        posOrig <- GenomicRanges::start(se)[SparseArray::nnawhich(assayOrig)]
        posExp <- GenomicRanges::start(seExp2)[SparseArray::nnawhich(assayExp)]
        expect_identical(posOrig, posExp)
    }
    expect_identical(as.character(unique(GenomicRanges::strand(seExp2))), "*")
})
