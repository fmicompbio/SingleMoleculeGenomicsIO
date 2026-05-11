suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(GenomicRanges)
    library(Biostrings)
    library(BSgenome)
})

## -------------------------------------------------------------------------- ##
## Checks, extractSeqContext
## -------------------------------------------------------------------------- ##
test_that("extractSeqContext works", {
    # example data
    ref <- system.file("extdata", "reference.fa.gz", package = "SingleMoleculeGenomicsIO")
    gnm <- Biostrings::readDNAStringSet(ref)
    gnm2 <- DNAStringSet(c(chr1 = "CATGTCCCCT"))
    regions <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(start = 6957060 - c(4, 2, 0),
                                  width = 1, names = c("x", "y", "z")))
    regions2 <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(start = 1 + c(0, 2, 4, 8),
                                  width = 3, names = c("a", "b", "c", "d")))
    regions3 <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(start = 1 + c(0, 2, 4, 8),
                                  width = 3, names = c("a", "b", "c", "d")),
        strand = "-")
    regions4 <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(start = 1 + c(0, 2, 4, 8),
                                  width = 3, names = c("a", "b", "c", "d")),
        strand = "+")
    # regions not overlapping with data range
    regions_wrong <- GenomicRanges::GRanges(
        seqnames = "chr2",
        ranges = IRanges::IRanges(start = 6957060 - c(4, 2, 0),
                                  width = 1, names = c("x", "y", "z")))
    se <- SummarizedExperiment(assays = matrix(1:3, ncol = 1), rowRanges = regions)

    # invalid arguments
    expect_error(extractSeqContext(x = "error"))
    expect_error(extractSeqContext(x = regions, sequenceContextWidth = -1))
    expect_error(extractSeqContext(x = regions, sequenceContextWidth = "error"))
    expect_error(extractSeqContext(x = regions, sequenceContextWidth = 7))
    expect_error(extractSeqContext(x = regions, sequenceContextWidth = 7, sequenceReference = "error"))
    expect_error(extractSeqContext(x = regions_wrong, sequenceContextWidth = 3,
                                   sequenceReference = ref),
                 "Not all chromosome names from .x. are present")


    # expected results
    expect_warning(s1 <- extractSeqContext(x = regions, sequenceContextWidth = 6, sequenceReference = ref))
    s2 <- extractSeqContext(x = regions, sequenceContextWidth = 7, sequenceReference = gnm)
    s4 <- extractSeqContext(x = unname(regions), sequenceContextWidth = 7, sequenceReference = gnm)
    s5 <- extractSeqContext(x = regions2, sequenceContextWidth = 7, sequenceReference = gnm)
    s6 <- extractSeqContext(x = resize(regions2, width = 1L, fix = "center"), sequenceContextWidth = 7, sequenceReference = gnm)
    s7 <- extractSeqContext(x = se, sequenceContextWidth = 7, sequenceReference = gnm)
    s8 <- extractSeqContext(x = regions2, sequenceContextWidth = 7, sequenceReference = gnm2)
    s9 <- extractSeqContext(x = regions3, sequenceContextWidth = 7, sequenceReference = gnm2)
    s10 <- extractSeqContext(x = c(regions4[seq(1, 2)], regions3[seq(1, 2)],
                                   regions4[seq(3, 4)], regions3[seq(3, 4)]),
                             sequenceContextWidth = 7, sequenceReference = gnm2)
    expect_s4_class(s1, "DNAStringSet")
    expect_s4_class(s2, "DNAStringSet")
    expect_s4_class(s4, "DNAStringSet")
    expect_s4_class(s5, "DNAStringSet")
    expect_s4_class(s6, "DNAStringSet")
    expect_s4_class(s7, "DNAStringSet")
    expect_s4_class(s8, "DNAStringSet")
    expect_s4_class(s9, "DNAStringSet")
    expect_identical(as.character(s1), c(x="AAAGGGG", y="AGGGGAN", z="GGGANNN"))
    expect_identical(s1, s2)
    expect_identical(unname(s1), s4)
    expect_identical(as.character(s5), c(a="NNNNNNN", b="NNNNNNN", c="NNNNNNN", d="NNNNNNN"))
    expect_identical(s5, s6)
    expect_identical(s7, s1)
    expect_identical(as.character(s8), c(a="NNCATGT", b="CATGTCC", c="TGTCCCC", d="CCCTNNN"))
    expect_identical(as.character(s9), c(a="ACATGNN", b="GGACATG", c="GGGGACA", d="NNNAGGG"))
    expect_identical(as.character(s10), c(a="NNCATGT", b="CATGTCC", a="ACATGNN", b="GGACATG",
                                          c="TGTCCCC", d="CCCTNNN", c="GGGGACA", d="NNNAGGG"))
})

## -------------------------------------------------------------------------- ##
## Checks, addSeqContext
## -------------------------------------------------------------------------- ##
test_that("addSeqContext works", {
    # example data
    ref <- system.file("extdata", "reference.fa.gz", package = "SingleMoleculeGenomicsIO")
    regions <- GenomicRanges::GRanges(seqnames = "chr1",
                                      ranges = IRanges::IRanges(start = 6957060 - c(4, 2, 0),
                                                                width = 1, names = c("x", "y", "z")))
    se <- SummarizedExperiment(assays = matrix(1:3, ncol = 1), rowRanges = regions)

    # invalid arguments
    expect_error(addSeqContext(x = "error"))
    expect_error(addSeqContext(x = regions, sequenceContextWidth = -1))
    expect_error(addSeqContext(x = regions, sequenceContextWidth = 7, sequenceReference = "error"))

    # expected results
    expect_warning(se1 <- addSeqContext(x = se, sequenceContextWidth = 6, sequenceReference = ref))
    expect_s4_class(se1, "RangedSummarizedExperiment")
    expect_identical(dim(se), dim(se1))
    expect_true("sequenceContext" %in% colnames(rowData(se1)))
    expect_identical(as.character(rowData(se1)$sequenceContext),
                     c(x = "AAAGGGG", y = "AGGGGAN", z = "GGGANNN"))
})
