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
    suppressWarnings(gpos <- c(SummarizedExperiment::rowRanges(se1),
                               SummarizedExperiment::rowRanges(se2)))
    modmat <- make_zero_col_DFrame(nrow = length(gpos))
    for (nm in colnames(se)) {
        x1 <- assay(se1, "mod_prob")[[nm]]
        x2 <- assay(se2, "mod_prob")[[nm]]
        colnames(x2) <- paste0("n-", colnames(x2))
        namat <- NaArray(dim = c(length(gpos), ncol(x1) + ncol(x2)),
                         dimnames = list(NULL, c(colnames(x1), colnames(x2))),
                         type = "double")
        i1 <- match(SummarizedExperiment::rowRanges(se1), gpos)
        j1 <- seq.int(ncol(x1))
        namat[cbind(rep(i1, length(j1)), rep(j1, each = length(i1)))] <- c(as.matrix(x1))
        i2 <- match(SummarizedExperiment::rowRanges(se2), gpos)
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

    # repeat some positions
    se1 <- se
    se2 <- se1[2:5, ]
    rr <- SummarizedExperiment::rowRanges(se2)
    strand(rr) <- ifelse(strand(rr) == "+", "-", "+")
    SummarizedExperiment::rowRanges(se2) <- rr
    suppressWarnings(gpos <- c(SummarizedExperiment::rowRanges(se1),
                               SummarizedExperiment::rowRanges(se2)))
    gpos <- gpos[order(start(gpos))]
    modmat <- make_zero_col_DFrame(nrow = length(gpos))
    for (nm in colnames(se)) {
        x1 <- assay(se1, "mod_prob")[[nm]]
        x2 <- assay(se2, "mod_prob")[[nm]] / 2
        colnames(x2) <- paste0("n-", colnames(x2))
        namat <- NaArray(dim = c(length(gpos), ncol(x1) + ncol(x2)),
                         dimnames = list(NULL, c(colnames(x1), colnames(x2))),
                         type = "double")
        i1 <- match(SummarizedExperiment::rowRanges(se1), gpos)
        j1 <- seq.int(ncol(x1))
        namat[cbind(rep(i1, length(j1)), rep(j1, each = length(i1)))] <- c(as.matrix(x1))
        i2 <- match(SummarizedExperiment::rowRanges(se2), gpos)
        j2 <- (ncol(x1) + seq.int(ncol(x2)))
        namat[cbind(rep(i2, length(j2)), rep(j2, each = length(i2)))] <- c(as.matrix(x2))
        modmat[[nm]] <- namat
    }
    seDupPos <- SummarizedExperiment(
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
    # expect_error(expandSEToBaseSpace(se = se, ignore.strand = FALSE),
    #              "All values in .ignore.strand. must be one of: TRUE")

    # empty input
    seExp0 <- expandSEToBaseSpace(se = se[numeric(0), ])
    expect_identical(dim(seExp0), c(0L, ncol(se)))
    expect_identical(SummarizedExperiment::assayNames(seExp0),
                     getReadLevelAssayNames(se))
    expect_identical(lapply(SummarizedExperiment::assay(se, "mod_prob"), colnames),
                     lapply(SummarizedExperiment::assay(seExp0, "mod_prob"), colnames))

    # expected return value - ignore.strand=TRUE
    seExp1 <- expandSEToBaseSpace(se = se, ignore.strand = TRUE)
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
    seExp2 <- expandSEToBaseSpace(se = se, regions = reg2, ignore.strand = TRUE)
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
    seTwoChromExp <- expandSEToBaseSpace(se = seTwoChrom,
                                         regions = "chr1:6940000-6941000",
                                         ignore.strand = TRUE)
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
    seTwoChromExp <- expandSEToBaseSpace(se = seTwoChrom, regions = reg3,
                                         ignore.strand = TRUE)
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

    # multiple chromosomes - ignore.strand=FALSE
    reg3 <- GRanges(seqnames = c("chr1", "chr2"),
                    ranges = IRanges(start = c(6940500, 6940200), width = 300))
    seTwoChromExp <- expandSEToBaseSpace(se = seTwoChrom, regions = reg3,
                                         ignore.strand = FALSE)
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
    expect_identical(GenomicRanges::strand(subsetByOverlaps(seTwoChrom, seTwoChromExp)),
                     GenomicRanges::strand(seTwoChromExp)[GenomicRanges::strand(seTwoChromExp) != "*"])

    # duplicated positions - ignore.strand=FALSE
    expect_message(
        seDupPosExp <- expandSEToBaseSpace(se = seDupPos,
                                           ignore.strand = FALSE, verbose = TRUE),
        "4 positions filtered out")
    expect_s4_class(seDupPosExp, "RangedSummarizedExperiment")
    expect_identical(dim(seDupPosExp),
                     c(width(range(SummarizedExperiment::rowRanges(seDupPos),
                                   ignore.strand = TRUE)),
                       ncol(seDupPos)))
    expect_identical(SummarizedExperiment::assayNames(seDupPosExp),
                     getReadLevelAssayNames(seDupPos))
    expect_identical(lapply(SummarizedExperiment::assay(seDupPos, "mod_prob"), colnames),
                     lapply(SummarizedExperiment::assay(seDupPosExp, "mod_prob"), colnames))
    for (s in colnames(seDupPos)) {
        assayOrig <- SummarizedExperiment::assay(seDupPos, "mod_prob")[[s]]
        assayOrig <- assayOrig[, seq_len(ncol(assayOrig) / 2)]
        assayExp <- SummarizedExperiment::assay(seDupPosExp, "mod_prob")[[s]]

        valsOrig <- SparseArray::nnavals(assayOrig)
        valsExp <- SparseArray::nnavals(assayExp)
        expect_identical(valsOrig, valsExp)

        posOrig <- GenomicRanges::start(seDupPos)[SparseArray::nnawhich(assayOrig)]
        posExp <- GenomicRanges::start(seDupPosExp)[SparseArray::nnawhich(assayExp)]
        expect_identical(posOrig, posExp)
    }
    expect_identical(GenomicRanges::strand(subsetByOverlaps(seDupPos, seDupPosExp)),
                     GenomicRanges::strand(seDupPosExp)[GenomicRanges::strand(seDupPosExp) != "*"])
})
