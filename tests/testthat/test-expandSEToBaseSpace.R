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

    # change the chromosome for some of the positions
    set.seed(637128L)
    se1 <- subsetReads(se, randomSubset = 1)
    se2 <- subsetReads(se, reads = unlist(lapply(assay(se1, "mod_prob"), colnames)),
                       invert = TRUE)
    Seqinfo::seqlevels(SummarizedExperiment::rowRanges(se2)) <- "chr2"
    suppressWarnings(gpos <- c(rowRanges(se1), rowRanges(se2)))
    modmat <- make_zero_col_DFrame(nrow = length(gpos))
    for (nm in colnames(se)) {
        x1 <- assay(se1, "mod_prob")[[nm]]
        x2 <- assay(se2, "mod_prob")[[nm]]
        namat <- NaArray(dim = c(length(gpos), ncol(x1) + ncol(x2)),
                         dimnames = list(NULL, c(colnames(x1), colnames(x2))),
                         type = "double")
        i1 <- match(rowRanges(se1), gpos)
        j1 <- seq.int(ncol(x1))
        namat[cbind(rep(i1, length(j1)), rep(j1, each = length(i1)))] <- c(as.matrix(x1))
        i2 <- match(rowRanges(se2), gpos)
        j2 <- (ncol(x1) + seq.int(ncol(x2)))
        namat[cbind(rep(i2, length(j2)), rep(j2, each = length(i2)))] <- c(as.matrix(x2))
        modmat[[nm]] <- namat
    }
    seTwoChrom <- SummarizedExperiment(
        assays = list(mod_prob = modmat),
        colData = colData(se),
        metadata = metadata(se),
        rowRanges = gpos
    )

    # expected fails
    # expect_error(expandSEToBaseSpace(se = se, regions = "chr1:6940100-6941000"),
    #              "Not all positions in .se. are within .region.")
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

    # expanded region
    reg2 <- GenomicRanges::GRanges("chr1",
                                   IRanges::IRanges(6939000, 6942000))
    seExp2 <- expandSEToBaseSpace(se = se, regions = reg2)
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

    # multiple subregions
    seTwoChromExp <- expandSEToBaseSpace(se = seTwoChrom, regions = "chr1:6940000-6941000")
    expect_s4_class(seTwoChromExp, "RangedSummarizedExperiment")
    expect_identical(dim(seTwoChromExp), c(1001L, ncol(seTwoChrom)))
    expect_identical(SummarizedExperiment::assayNames(seTwoChromExp),
                     getReadLevelAssayNames(seTwoChrom))
    expect_identical(lapply(assay(seTwoChrom, "mod_prob"), colnames),
                     lapply(assay(seTwoChromExp, "mod_prob"), colnames))
    for (s in colnames(seTwoChrom)) {
        assayOrig <- assay(seTwoChrom[seqnames(seTwoChrom) == "chr1", ], "mod_prob")[[s]]
        assayExp <- assay(seTwoChromExp, "mod_prob")[[s]]

        valsOrig <- SparseArray::nnavals(assayOrig)
        valsExp <- SparseArray::nnavals(assayExp)
        expect_identical(valsOrig, valsExp)

        posOrig <- GenomicRanges::start(seTwoChrom[seqnames(seTwoChrom) == "chr1", ])[SparseArray::nnawhich(assayOrig)]
        posExp <- GenomicRanges::start(seTwoChromExp)[SparseArray::nnawhich(assayExp)]
        expect_identical(posOrig, posExp)
    }
    expect_identical(as.character(unique(GenomicRanges::strand(seTwoChromExp))), "*")

    # multiple chromosomes
    reg3 <- GRanges(seqnames = c("chr1", "chr2"),
                    ranges = IRanges(start = c(6940500, 6940200), width = 300))
    seTwoChromExp <- expandSEToBaseSpace(se = seTwoChrom, regions = reg3)
    expect_s4_class(seTwoChromExp, "RangedSummarizedExperiment")
    expect_identical(dim(seTwoChromExp), c(sum(width(reg3)), ncol(seTwoChrom)))
    expect_identical(SummarizedExperiment::assayNames(seTwoChromExp),
                     getReadLevelAssayNames(seTwoChrom))
    expect_identical(lapply(assay(seTwoChrom, "mod_prob"), colnames),
                     lapply(assay(seTwoChromExp, "mod_prob"), colnames))
    for (s in colnames(seTwoChrom)) {
        assayOrig <- assay(subsetByOverlaps(seTwoChrom, reg3), "mod_prob")[[s]]
        assayExp <- assay(seTwoChromExp, "mod_prob")[[s]]

        valsOrig <- SparseArray::nnavals(assayOrig)
        valsExp <- SparseArray::nnavals(assayExp)
        expect_identical(valsOrig, valsExp)

        posOrig <- GenomicRanges::start(subsetByOverlaps(seTwoChrom, reg3))[SparseArray::nnawhich(assayOrig)]
        posExp <- GenomicRanges::start(seTwoChromExp)[SparseArray::nnawhich(assayExp)]
        expect_identical(posOrig, posExp)
    }
    expect_identical(as.character(unique(GenomicRanges::strand(seTwoChromExp))), "*")
})
