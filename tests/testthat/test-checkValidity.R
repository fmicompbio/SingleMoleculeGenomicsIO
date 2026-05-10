test_that("validity checks work", {
    # example data
    fnames <- c(s1_5mC = system.file("extdata", "modkit_extract_rc_5mC_1.tsv.gz",
                                     package = "SingleMoleculeGenomicsIO"),
                s2_5mC = system.file("extdata", "modkit_extract_rc_5mC_2.tsv.gz",
                                     package = "SingleMoleculeGenomicsIO"),
                s1_6mA = system.file("extdata", "modkit_extract_rc_6mA_1.tsv.gz",
                                     package = "SingleMoleculeGenomicsIO"),
                s2_6mA = system.file("extdata", "modkit_extract_rc_6mA_2.tsv.gz",
                                     package = "SingleMoleculeGenomicsIO")
    )
    rme <- readModkitExtract(fnames = fnames[c("s1_5mC", "s2_5mC", "s1_6mA")],
                             modbase = c(s1_6mA = "a", s1_5mC = "m", s2_5mC = "m"),
                             filter = c(`m` = 0.6, `a` = 0.4, `-` = 0.3),
                             nrows = Inf, seqinfo = NULL,
                             BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
    suppressWarnings(expect_warning(
        rme <- addReadStats(rme, BPPARAM = BiocParallel::SerialParam()),
        "Too few points"))
    rme_withreads <- flattenReadLevelAssay(rme)
    rme_withoutreads <- flattenReadLevelAssay(rme, keepReads = FALSE)

    ## Test getReadLevelAssayNames
    expect_identical(getReadLevelAssayNames(rme_withreads), "mod_prob")
    expect_identical(getReadLevelAssayNames(rme_withoutreads), character(0))

    ## Test getReadNamesBySample
    expect_identical(getReadNamesBySample(rme_withreads),
                     structure(
                         list(s1_5mC = colnames(assay(rme_withreads, "mod_prob")[[1]]),
                              s2_5mC = colnames(assay(rme_withreads, "mod_prob")[[2]]),
                              s1_6mA = colnames(assay(rme_withreads, "mod_prob")[[3]])),
                         source = "assay mod_prob"))
    expect_equal(getReadNamesBySample(rme_withreads),
                 getReadNamesBySample(rme_withoutreads),
                 ignore_attr = TRUE)
    expect_identical(lengths(getReadNamesBySample(rme_withreads)),
                     c(s1_5mC = 10L, s2_5mC = 10L, s1_6mA = 10L))
    tmp <- rme_withoutreads
    tmp$QC <- tmp$readInfo <- NULL
    expect_error(getReadNamesBySample(tmp),
                 "does not contain any read-level assays or colData columns")

    ## Test checkSEValidity
    expect_no_error(checkSEValidity(rme_withreads))
    expect_no_error(checkSEValidity(rme_withoutreads))

    suppressMessages(expect_message(checkSEValidity(rme_withreads, verbose = TRUE)))
    suppressMessages(expect_message(checkSEValidity(rme_withoutreads, verbose = TRUE)))

    rme1 <- rme_withreads
    SummarizedExperiment::assay(rme1, "test") <- SummarizedExperiment::assay(rme1, "mod_prob")
    metadata(rme1)$readLevelData$assayNames <- c(metadata(rme1)$readLevelData$assayNames,
                                                 "test")
    suppressMessages(expect_message(checkSEValidity(rme1, verbose = TRUE)))

    rme1 <- rme_withreads
    rownames(rme1) <- as.character(rowRanges(rme1))
    rownames(rme1)[2] <- rownames(rme1)[1]
    expect_error(checkSEValidity(rme1),
                 "anyDuplicated(rownames(se)) == 0L is not TRUE", fixed = TRUE)

    rme1 <- rme_withreads
    SummarizedExperiment::assayNames(rme1) <- c("", "", "", "")
    expect_error(checkSEValidity(rme1),
                 '!is.null(assayNames(se)) && all(assayNames(se) != "") && anyDuplicated(assayNames(se)) ==  .... is not TRUE', fixed = TRUE)

    rme1 <- rme_withreads
    expect_length(assays(rme1), 4L)
    assays(rme1) <- list(assays(rme1)[[1]], assays(rme1)[[2]], assays(rme1)[[3]],
                         assays(rme1)[[4]])
    expect_null(assayNames(rme1))
    expect_error(checkSEValidity(rme1),
                 '!is.null(assayNames(se)) && all(assayNames(se) != "") && anyDuplicated(assayNames(se)) ==  .... is not TRUE', fixed = TRUE)

    rme1 <- rme_withreads
    rme1$QC <- rme1$QC[c(3, 1, 2)]
    expect_error(checkSEValidity(rme1),
                 "colnames(se) are not all TRUE", fixed = TRUE)

    rme1 <- rme_withreads
    SummarizedExperiment::colData(rme1) <- SummarizedExperiment::colData(rme1)[c(3, 1, 2), ]
    expect_error(checkSEValidity(rme1),
                 "colnames(assay(se, an, withDimnames = FALSE)) == colnames(se) are not all TRUE", fixed = TRUE)

    rme1 <- rme_withreads
    colnames(rme1) <- colnames(rme1)[c(3, 1, 2)]
    expect_error(checkSEValidity(rme1),
                 "colnames(assay(se, an, withDimnames = FALSE)) == colnames(se) are not all TRUE", fixed = TRUE)

    rme1 <- rme_withreads
    SummarizedExperiment::assay(rme1, "mod_prob", withDimnames = FALSE) <-
        SummarizedExperiment::assay(rme1, "mod_prob")[, c(3, 1, 2)]
    expect_error(checkSEValidity(rme1),
                 "colnames(assay(se, an, withDimnames = FALSE)) == colnames(se) are not all TRUE", fixed = TRUE)

    rme1 <- rme_withreads
    colnames(SummarizedExperiment::colData(rme1))[2] <- "one:two"
    expect_error(checkSEValidity(rme1),
                 "Column names in .colData.se.. can not contain")

    rme1 <- rme_withreads
    rme1$QC[[2]] <- rme1$QC[[2]][1:5, ]
    expect_error(checkSEValidity(rme1),
                 "Mismatching reads for assay mod_prob and colData column QC, sample s2_5mC")

    rme1 <- rme_withreads
    SummarizedExperiment::assay(rme1, "mod_prob")[[1]] <-
        SummarizedExperiment::assay(rme1, "mod_prob")[[1]][, 1:5]
    expect_error(checkSEValidity(rme1),
                 "Mismatching reads for assay mod_prob and colData column readInfo, sample s1_5mC")

    rme1 <- rme_withreads
    SummarizedExperiment::assay(rme1, "test") <- SummarizedExperiment::assay(rme1, "mod_prob")
    metadata(rme1)$readLevelData$assayNames <- c(metadata(rme1)$readLevelData$assayNames,
                                                 "test")
    SummarizedExperiment::assay(rme1, "mod_prob")[[1]] <-
        SummarizedExperiment::assay(rme1, "mod_prob")[[1]][, 1:5]
    expect_error(checkSEValidity(rme1),
                 "Mismatching reads for assay mod_prob and test, sample s1_5mC")

    rme1 <- rme_withreads
    SummarizedExperiment::assay(rme1, "test") <- SummarizedExperiment::assay(rme1, "mod_prob")
    metadata(rme1)$readLevelData$assayNames <- c(metadata(rme1)$readLevelData$assayNames,
                                                 "test")
    N <- ncol(SummarizedExperiment::assay(rme1, "mod_prob")[[1]])
    set.seed(123L)
    SummarizedExperiment::assay(rme1, "mod_prob")[[1]] <-
        SummarizedExperiment::assay(rme1, "mod_prob")[[1]][, sample.int(N, N)]
    expect_error(checkSEValidity(rme1),
                 "Mismatching reads for assay mod_prob and test, sample s1_5mC")
})
