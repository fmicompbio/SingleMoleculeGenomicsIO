suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(GenomicRanges)
})

## -------------------------------------------------------------------------- ##
## Checks, readModkitExtract
## -------------------------------------------------------------------------- ##
test_that("readModkitExtract works", {
    # example data
    ref <- system.file("extdata", "reference.fa.gz", package = "SingleMoleculeGenomicsIO")
    fnames <- c(
        s1_5mC = system.file("extdata", "modkit_extract_rc_5mC_1.tsv.gz",
                             package = "SingleMoleculeGenomicsIO"),
        s2_5mC = system.file("extdata", "modkit_extract_rc_5mC_2.tsv.gz",
                             package = "SingleMoleculeGenomicsIO"),
        s1_6mA = system.file("extdata", "modkit_extract_rc_6mA_1.tsv.gz",
                             package = "SingleMoleculeGenomicsIO"),
        s2_6mA = system.file("extdata", "modkit_extract_rc_6mA_2.tsv.gz",
                             package = "SingleMoleculeGenomicsIO")
    )
    sample_annot <- data.frame(
        sample = c("s1_5mC", "s2_5mC", "s1_6mA", "s2_6mA"),
        group = c("g1", "g1", "g2", "g2"))

    # invalid arguments
    expect_error(readModkitExtract(
        fnames = "error",
        BPPARAM = BiocParallel::SerialParam()),
        "not all .fnames. exist")
    expect_error(readModkitExtract(
        fnames = stats::setNames(c(fnames[[1]], fnames[[2]]), c("s1", "s1")),
        BPPARAM = BiocParallel::SerialParam()),
        ".names\\(fnames\\). are not unique")
    expect_error(readModkitExtract(
        fnames = fnames, modbase = NULL,
        BPPARAM = BiocParallel::SerialParam()),
        ".modbase. must not be .NULL.")
    expect_error(readModkitExtract(
        fnames = fnames, modbase = 1,
        BPPARAM = BiocParallel::SerialParam()),
        ".modbase. must be of class .character.")
    expect_error(readModkitExtract(
        fnames = fnames, modbase = c("m", "a"),
        BPPARAM = BiocParallel::SerialParam()),
        ".modbase. must have length 4")
    expect_error(readModkitExtract(
        fnames = fnames, modbase = "x",
        BPPARAM = BiocParallel::SerialParam()),
        "invalid .modbase. values")
    expect_error(readModkitExtract(
        fnames = fnames, modbase = c(s1 = "m", s2 = "m", s3 = "a", s4 = "a"),
        BPPARAM = BiocParallel::SerialParam()),
        "names of .modbase. and .fnames. don't agree")
    expect_error(readModkitExtract(
        fnames = fnames, modbase = c("m", "m", "a", "a"),
        BPPARAM = BiocParallel::SerialParam(), filter = "error"),
        "All values in .filter. must be one of")
    expect_error(readModkitExtract(
        fnames = fnames, modbase = c("m", "m", "a", "a"),
        BPPARAM = BiocParallel::SerialParam(), filter = c(0.1, 0.2)),
        ".filter. must be a named vector")
    expect_error(readModkitExtract(
        fnames = fnames, modbase = c("m", "m", "a", "a"),
        BPPARAM = BiocParallel::SerialParam(), filter = c(c = 0.1)),
        "a filter threshold needs to be supplied")
    expect_error(readModkitExtract(
        fnames = fnames, modbase = c("m", "m", "a", "a"),
        BPPARAM = BiocParallel::SerialParam(), nrows = -1),
        ".nrows. must be between")
    expect_error(readModkitExtract(
        fnames = fnames, modbase = c("m", "m", "a", "a"),
        BPPARAM = BiocParallel::SerialParam(), nrows = "error"),
        ".nrows. must be of class .numeric.")
    expect_error(readModkitExtract(
        fnames = fnames, modbase = c("m", "m", "a", "a"),
        BPPARAM = BiocParallel::SerialParam(), seqinfo = c(100)),
        ".seqinfo. must be ")
    expect_error(readModkitExtract(
        fnames = fnames, modbase = c("m", "m", "a", "a"),
        BPPARAM = BiocParallel::SerialParam(), seqinfo = c(chr2 = 1000)),
        ".seqnames. contains sequence names with no entries")
    expect_error(readModkitExtract(
        fnames = fnames, modbase = c("m", "m", "a", "a"), BPPARAM = "error"),
        ".BPPARAM. must be of class .BiocParallelParam.")
    expect_error(readModkitExtract(
        fnames = fnames, modbase = c("m", "m", "a", "a"),
        BPPARAM = BiocParallel::SerialParam(), verbose = "error"),
        ".verbose. must be of class .logical.")
    expect_error(readModkitExtract(
        fnames = fnames, modbase = c("m", "m", "a", "a"),
        BPPARAM = BiocParallel::SerialParam(), verbose = c(TRUE, FALSE)),
        ".verbose. must have length 1")
    expect_error(readModkitExtract(
        fnames = fnames, modbase = c("m", "m", "a", "a"),
        BPPARAM = BiocParallel::SerialParam(),
        sampleAnnot = sample_annot[1:2, ]),
        "Annotation information missing")
    expect_error(readModkitExtract(
        fnames = fnames, modbase = c("m", "m", "a", "a"),
        BPPARAM = BiocParallel::SerialParam(),
        sampleAnnot = sample_annot[, -1, drop = FALSE]),
        ".sampleAnnot. must have at least a column")

    # expected results
    # ... single file, no filtering
    suppressMessages(
        expect_message(
            rme <- readModkitExtract(fnames = fnames["s1_5mC"], modbase = "m",
                                     filter = NULL, nrows = Inf, seqinfo = NULL,
                                     sequenceContextWidth = 1, sequenceReference = ref,
                                     BPPARAM = BiocParallel::SerialParam(),
                                     verbose = TRUE)
    ))
    expect_s4_class(rme, "RangedSummarizedExperiment")
    expect_identical(dim(rme), c(6432L, 1L)) ## number of unique positions
    expect_identical(colnames(rme), "s1_5mC")
    expect_length(SummarizedExperiment::assays(rme), 1L)
    expect_named(SummarizedExperiment::assays(rme), "mod_prob")
    expect_s4_class(SummarizedExperiment::assay(rme, "mod_prob"), "DataFrame")
    expect_identical(dim(SummarizedExperiment::assay(rme, "mod_prob")[[1]]), c(6432L, 10L))
    expect_s4_class(SummarizedExperiment::assay(rme)[[1]], "NaMatrix")
    expect_false(is.null(colnames(rme)))
    expect_false(is.null(colnames(SummarizedExperiment::assay(rme)[[1]])))
    expect_length(S4Vectors::metadata(rme), 3)
    expect_named(S4Vectors::metadata(rme), c("modkit_threshold",
                                             "filter_threshold",
                                             "readLevelData"))
    expect_equal(S4Vectors::metadata(rme)$modkit_threshold,
                 list(s1_5mC = c(`m` = 0.7988281, `-` = 0.9082031)),
                 ignore_attr = TRUE)
    expect_equal(S4Vectors::metadata(rme)$filter_threshold,
                 list(s1_5mC = NULL),
                 ignore_attr = TRUE)
    expect_equal(sum(SummarizedExperiment::assay(rme)[[1]], na.rm = TRUE), 2297.13868)
    expect_identical(SparseArray::nnacount(SummarizedExperiment::assay(rme)[[1]]), 18531L) ## number of rows in the original file
    expect_equal(unclass(table(as.character(SummarizedExperiment::rowData(rme)$sequenceContext))),
                 c(A = 82L, C = 6159L, G = 107L, T = 84L), ignore_attr = TRUE)

    # ... single file, manual filtering
    rme <- readModkitExtract(fnames = fnames[["s1_5mC"]], modbase = "m",
                             filter = c(`m` = 0.6, `-` = 0.5),
                             nrows = Inf, seqinfo = NULL,
                             BPPARAM = BiocParallel::MulticoreParam(2L), verbose = FALSE)
    expect_s4_class(rme, "RangedSummarizedExperiment")
    expect_identical(dim(rme), c(6415L, 1L)) ## number of unique positions
    expect_identical(colnames(rme), "s1")
    expect_length(SummarizedExperiment::assays(rme), 1L)
    expect_named(SummarizedExperiment::assays(rme), "mod_prob")
    expect_s4_class(SummarizedExperiment::assay(rme, "mod_prob"), "DataFrame")
    expect_identical(dim(SummarizedExperiment::assay(rme, "mod_prob")[[1]]), c(6415L, 10L))
    expect_s4_class(SummarizedExperiment::assay(rme)[[1]], "NaMatrix")
    expect_false(is.null(colnames(rme)))
    expect_false(is.null(colnames(SummarizedExperiment::assay(rme)[[1]])))
    expect_length(S4Vectors::metadata(rme), 3)
    expect_named(S4Vectors::metadata(rme), c("modkit_threshold",
                                             "filter_threshold",
                                             "readLevelData"))
    expect_equal(S4Vectors::metadata(rme)$modkit_threshold,
                 list(s1_5mC = c(`m` = 0.7988281, `-` = 0.9082031)),
                 ignore_attr = TRUE)
    expect_equal(S4Vectors::metadata(rme)$filter_threshold,
                 list(s1_5mC = c(`m` = 0.6, `-` = 0.5)),
                 ignore_attr = TRUE)
    expect_equal(sum(SummarizedExperiment::assay(rme)[[1]], na.rm = TRUE), 2243.82032)
    expect_identical(SparseArray::nnacount(SummarizedExperiment::assay(rme)[[1]]), 18434L) ## number of rows in the original file

    # ... single file, automatic filtering
    rme <- readModkitExtract(fnames = fnames["s1_5mC"], modbase = "m",
                             filter = "modkit",
                             nrows = Inf, seqinfo = NULL,
                             BPPARAM = BiocParallel::SerialParam(),
                             verbose = FALSE)
    expect_s4_class(rme, "RangedSummarizedExperiment")
    expect_identical(dim(rme), c(5893L, 1L)) ## number of unique positions
    expect_identical(colnames(rme), "s1_5mC")
    expect_length(SummarizedExperiment::assays(rme), 1)
    expect_named(SummarizedExperiment::assays(rme), "mod_prob")
    expect_s4_class(SummarizedExperiment::assay(rme, "mod_prob"), "DataFrame")
    expect_identical(dim(SummarizedExperiment::assay(rme, "mod_prob")[[1]]), c(5893L, 10L))
    expect_s4_class(SummarizedExperiment::assay(rme)[[1]], "NaMatrix")
    expect_false(is.null(colnames(rme)))
    expect_false(is.null(colnames(SummarizedExperiment::assay(rme)[[1]])))
    expect_length(S4Vectors::metadata(rme), 3)
    expect_named(S4Vectors::metadata(rme), c("modkit_threshold",
                                             "filter_threshold",
                                             "readLevelData"))
    expect_equal(S4Vectors::metadata(rme)$modkit_threshold,
                 list(s1_5mC = c(`m` = 0.7988281, `-` = 0.9082031)),
                 ignore_attr = TRUE)
    expect_equal(S4Vectors::metadata(rme)$filter_threshold,
                 list(s1_5mC = c(`m` = 0.7988281, `-` = 0.9082031)),
                 ignore_attr = TRUE)
    expect_equal(sum(SummarizedExperiment::assay(rme)[[1]], na.rm = TRUE), 1649.92774)
    expect_identical(SparseArray::nnacount(SummarizedExperiment::assay(rme)[[1]]), 15325L) ## number of rows in the original file

    # ... multiple files, no filtering
    rme <- readModkitExtract(fnames = fnames[c("s1_5mC", "s2_5mC",
                                               "s1_6mA")],
                             modbase = c("m", "m", "a"),
                             filter = NULL, nrows = Inf, seqinfo = NULL,
                             BPPARAM = BiocParallel::SerialParam(),
                             verbose = FALSE)
    expect_s4_class(rme, "RangedSummarizedExperiment")
    expect_identical(dim(rme), c(18655L, 3L)) ## number of unique positions
    expect_identical(colnames(rme), c("s1_5mC", "s2_5mC", "s1_6mA"))
    expect_identical(colnames(SummarizedExperiment::colData(rme)),
                     c("sample", "modbase", "readInfo"))
    expect_equal(rme$modbase, c("m", "m", "a"), ignore_attr = TRUE)
    expect_length(SummarizedExperiment::assays(rme), 1L)
    expect_named(SummarizedExperiment::assays(rme), "mod_prob")
    expect_s4_class(SummarizedExperiment::assay(rme, "mod_prob"), "DataFrame")
    expect_identical(dim(SummarizedExperiment::assay(rme, "mod_prob")[[1]]), c(18655L, 10L))
    expect_identical(dim(SummarizedExperiment::assay(rme, "mod_prob")[[2]]), c(18655L, 10L))
    expect_identical(dim(SummarizedExperiment::assay(rme, "mod_prob")[[3]]), c(18655L, 10L))
    expect_identical(dim(as.matrix(SummarizedExperiment::assay(rme, "mod_prob"))), c(18655L, 30L))
    expect_s4_class(SummarizedExperiment::assay(rme)[[1]], "NaMatrix")
    expect_s4_class(SummarizedExperiment::assay(rme)[[2]], "NaMatrix")
    expect_s4_class(SummarizedExperiment::assay(rme)[[3]], "NaMatrix")
    expect_false(is.null(colnames(rme)))
    expect_length(S4Vectors::metadata(rme), 3)
    expect_named(S4Vectors::metadata(rme), c("modkit_threshold",
                                             "filter_threshold",
                                             "readLevelData"))
    expect_equal(S4Vectors::metadata(rme)$modkit_threshold,
                 list(s1_5mC = c(`m` = 0.7988281, `-` = 0.9082031),
                      s2_5mC = c(`m` = 0.7988281, `-` = 0.9003906),
                      s1_6mA = c(`a` = -Inf, `-` = 0.8964844)),
                 ignore_attr = TRUE)
    expect_equal(S4Vectors::metadata(rme)$filter_threshold,
                 list(s1_5mC = NULL, s2_5mC = NULL, s1_6mA = NULL),
                 ignore_attr = TRUE)
    expect_equal(sum(as.matrix(SummarizedExperiment::assay(rme)), na.rm = TRUE), 8236.457)
    expect_identical(SparseArray::nnacount(as.matrix(SummarizedExperiment::assay(rme))), 71750L) ## total number of rows in the original files

    # ... multiple files, no filtering, with sample annotation
    rme <- readModkitExtract(fnames = fnames[c("s1_5mC", "s2_5mC",
                                               "s1_6mA")],
                             modbase = c("m", "m", "a"),
                             sampleAnnot = sample_annot,
                             filter = NULL, nrows = Inf, seqinfo = NULL,
                             BPPARAM = BiocParallel::SerialParam(),
                             verbose = FALSE)
    expect_s4_class(rme, "RangedSummarizedExperiment")
    expect_identical(dim(rme), c(18655L, 3L)) ## number of unique positions
    expect_identical(colnames(rme), c("s1_5mC", "s2_5mC", "s1_6mA"))
    expect_identical(colnames(SummarizedExperiment::colData(rme)),
                     c("sample", "modbase", "readInfo", "group"))
    expect_equal(rme$modbase, c("m", "m", "a"), ignore_attr = TRUE)
    expect_equal(rme$group, c("g1", "g1", "g2"), ignore_attr = TRUE)
    expect_length(SummarizedExperiment::assays(rme), 1)
    expect_named(SummarizedExperiment::assays(rme), "mod_prob")
    expect_s4_class(SummarizedExperiment::assay(rme, "mod_prob"), "DataFrame")
    expect_identical(dim(SummarizedExperiment::assay(rme, "mod_prob")[[1]]), c(18655L, 10L))
    expect_identical(dim(SummarizedExperiment::assay(rme, "mod_prob")[[2]]), c(18655L, 10L))
    expect_identical(dim(SummarizedExperiment::assay(rme, "mod_prob")[[3]]), c(18655L, 10L))
    expect_identical(dim(as.matrix(SummarizedExperiment::assay(rme, "mod_prob"))), c(18655L, 30L))
    expect_s4_class(SummarizedExperiment::assay(rme)[[1]], "NaMatrix")
    expect_s4_class(SummarizedExperiment::assay(rme)[[2]], "NaMatrix")
    expect_s4_class(SummarizedExperiment::assay(rme)[[3]], "NaMatrix")
    expect_false(is.null(colnames(rme)))
    expect_length(S4Vectors::metadata(rme), 3)
    expect_named(S4Vectors::metadata(rme), c("modkit_threshold",
                                             "filter_threshold",
                                             "readLevelData"))
    expect_equal(S4Vectors::metadata(rme)$modkit_threshold,
                 list(s1_5mC = c(`m` = 0.7988281, `-` = 0.9082031),
                      s2_5mC = c(`m` = 0.7988281, `-` = 0.9003906),
                      s1_6mA = c(`a` = -Inf, `-` = 0.8964844)),
                 ignore_attr = TRUE)
    expect_equal(S4Vectors::metadata(rme)$filter_threshold,
                 list(s1_5mC = NULL, s2_5mC = NULL, s1_6mA = NULL),
                 ignore_attr = TRUE)
    expect_equal(sum(as.matrix(SummarizedExperiment::assay(rme)), na.rm = TRUE), 8236.457)
    expect_identical(SparseArray::nnacount(as.matrix(SummarizedExperiment::assay(rme))), 71750L) ## total number of rows in the original files

    # ... multiple files, manual filtering
    rme <- readModkitExtract(fnames = fnames[c("s1_5mC", "s2_5mC",
                                               "s1_6mA")],
                             modbase = c(s1_6mA = "a", s1_5mC = "m", s2_5mC = "m"),
                             filter = c(`m` = 0.6, `a` = 0.4, `-` = 0.3),
                             nrows = Inf, seqinfo = NULL,
                             BPPARAM = BiocParallel::SerialParam(),
                             verbose = FALSE)
    expect_s4_class(rme, "RangedSummarizedExperiment")
    expect_identical(dim(rme), c(18615L, 3L)) ## number of unique positions
    expect_identical(colnames(rme), c("s1_5mC", "s2_5mC", "s1_6mA"))
    expect_identical(rme$modbase, c("m", "m", "a"), ignore_attr = TRUE)
    expect_length(SummarizedExperiment::assays(rme), 1L)
    expect_named(SummarizedExperiment::assays(rme), "mod_prob")
    expect_s4_class(SummarizedExperiment::assay(rme, "mod_prob"), "DataFrame")
    expect_identical(dim(SummarizedExperiment::assay(rme, "mod_prob")[[1]]), c(18615L, 10L))
    expect_identical(dim(SummarizedExperiment::assay(rme, "mod_prob")[[2]]), c(18615L, 10L))
    expect_identical(dim(SummarizedExperiment::assay(rme, "mod_prob")[[3]]), c(18615L, 10L))
    expect_identical(dim(as.matrix(SummarizedExperiment::assay(rme, "mod_prob"))), c(18615L, 30L))
    expect_s4_class(SummarizedExperiment::assay(rme)[[1]], "NaMatrix")
    expect_s4_class(SummarizedExperiment::assay(rme)[[2]], "NaMatrix")
    expect_s4_class(SummarizedExperiment::assay(rme)[[3]], "NaMatrix")
    expect_false(is.null(colnames(rme)))
    expect_length(S4Vectors::metadata(rme), 3)
    expect_named(S4Vectors::metadata(rme), c("modkit_threshold",
                                             "filter_threshold",
                                             "readLevelData"))
    expect_equal(S4Vectors::metadata(rme)$modkit_threshold,
                 list(s1_5mC = c(`m` = 0.7988281, `-` = 0.9082031),
                      s2_5mC = c(`m` = 0.7988281, `-` = 0.9003906),
                      s1_6mA = c(`a` = -Inf, `-` = 0.8964844)),
                 ignore_attr = TRUE)
    expect_equal(S4Vectors::metadata(rme)$filter_threshold,
                 list(s1_5mC = c(`m` = 0.6, `a` = 0.4, `-` = 0.3),
                      s2_5mC = c(`m` = 0.6, `a` = 0.4, `-` = 0.3),
                      s1_6mA = c(`m` = 0.6, `a` = 0.4, `-` = 0.3)),
                 ignore_attr = TRUE)
    expect_equal(sum(as.matrix(SummarizedExperiment::assay(rme)), na.rm = TRUE), 8081.2012)
    expect_identical(SparseArray::nnacount(as.matrix(SummarizedExperiment::assay(rme))), 71467L) ## total number of rows in the original files

    # ... multiple files, automatic filtering
    rme <- readModkitExtract(fnames = fnames[c("s1_5mC", "s1_6mA", "s2_5mC")],
                             modbase = c("m", "a", "m"),
                             filter = "modkit",
                             nrows = Inf, seqinfo = NULL,
                             BPPARAM = BiocParallel::SerialParam(),
                             verbose = FALSE)
    expect_s4_class(rme, "RangedSummarizedExperiment")
    expect_identical(dim(rme), c(17459L, 3L)) ## number of unique positions
    expect_identical(colnames(rme), c("s1_5mC", "s1_6mA", "s2_5mC"))
    expect_equal(rme$modbase, c("m", "a", "m"), ignore_attr = TRUE)
    expect_length(SummarizedExperiment::assays(rme), 1L)
    expect_named(SummarizedExperiment::assays(rme), "mod_prob")
    expect_s4_class(SummarizedExperiment::assay(rme, "mod_prob"), "DataFrame")
    expect_identical(dim(SummarizedExperiment::assay(rme, "mod_prob")[[1]]), c(17459L, 10L))
    expect_identical(dim(SummarizedExperiment::assay(rme, "mod_prob")[[2]]), c(17459L, 10L))
    expect_identical(dim(SummarizedExperiment::assay(rme, "mod_prob")[[3]]), c(17459L, 10L))
    expect_identical(dim(as.matrix(SummarizedExperiment::assay(rme, "mod_prob"))), c(17459L, 30L))
    expect_s4_class(SummarizedExperiment::assay(rme)[[1]], "NaMatrix")
    expect_s4_class(SummarizedExperiment::assay(rme)[[2]], "NaMatrix")
    expect_s4_class(SummarizedExperiment::assay(rme)[[3]], "NaMatrix")
    expect_false(is.null(colnames(rme)))
    expect_length(S4Vectors::metadata(rme), 3)
    expect_named(S4Vectors::metadata(rme), c("modkit_threshold",
                                             "filter_threshold",
                                             "readLevelData"))
    expect_equal(S4Vectors::metadata(rme)$modkit_threshold,
                 list(s1_5mC = c(`m` = 0.7988281, `-` = 0.9082031),
                      s1_6mA = c(`a` = -Inf, `-` = 0.8964844),
                      s2_5mC = c(`m` = 0.7988281, `-` = 0.9003906)),
                 ignore_attr = TRUE)
    expect_equal(S4Vectors::metadata(rme)$filter_threshold,
                 list(s1_5mC = c(`m` = 0.7988281, `-` = 0.9082031),
                      s1_6mA = c(`a` = -Inf, `-` = 0.8964844),
                      s2_5mC = c(`m` = 0.7988281, `-` = 0.9003906)),
                 ignore_attr = TRUE)
    expect_equal(sum(as.matrix(SummarizedExperiment::assay(rme)), na.rm = TRUE), 5854.98049)
    expect_identical(SparseArray::nnacount(as.matrix(SummarizedExperiment::assay(rme))), 61228L) ## total number of rows in the original files
})
