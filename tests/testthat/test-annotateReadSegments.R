test_that("annotateReadSegments works", {
    # example data
    modbamfile <- system.file("extdata", "6mA_1_10reads.bam",
                               package = "SingleMoleculeGenomicsIO")
    se <- readModBam(bamfiles = modbamfile, regions = "chr1:6935400-6936300",
                     modbase = "a", verbose = FALSE,
                     BPPARAM = BiocParallel::SerialParam())
    segs <- list(s1 = IRanges::IRangesList(
         "s1-d52a5f6a-a60a-4f85-913e-eada84bfbfb9" = IRanges::IRanges(
             start = 6925834, width = 140),
         "s1-6cf74134-e550-4c02-bd2b-91385422ee25" = IRanges::IRanges(
             start = c(6926000, 6926200), width = 140)))

    # valid inputs
    expect_error(annotateReadSegments(se = "error"),
                 "must be of class .RangedSummarizedExperiment.")
    expect_error(annotateReadSegments(se = se, irlList = "error"),
                 "must be of class .list.")
    expect_error(annotateReadSegments(se = se, irlList = stats::setNames(segs, "wrong")),
                 "names of .irlList. must be identical to the column names of .se.")
    expect_error(annotateReadSegments(se = se, irlList = list(s1 = "error")),
                 ".irlList..i... must be of class .IRangesList.")
    se2 <- flattenReadLevelAssay(se, keepReads = FALSE)
    expect_error(annotateReadSegments(se = se2, irlList = segs),
                 ".se. must contain at least one read-level assay")
    expect_error(annotateReadSegments(se = se, irlList = list(
        s1 = IRanges::IRangesList("does-not-exist" = IRanges::IRanges(
            start = 6926000, width = 10)))),
                 "Not all read names")
    expect_error(annotateReadSegments(se = se, irlList = segs, name = 1L),
                 ".name. must be of class .character.")

    # expected results
    se3 <- annotateReadSegments(se, segs, "nucl")
    expect_s4_class(se3, "RangedSummarizedExperiment")
    expect_identical(dim(se3), dim(se))
    expect_true("nucl" %in% colnames(colData(se3)))
    expect_type(se3$nucl, "list")
    expect_named(se3$nucl, colnames(se))
    expect_s4_class(se3$nucl$s1, "IRangesList")
    expect_length(se3$nucl$s1, se3$n_reads)
    expect_identical(unname(lengths(se3$nucl$s1)),
                     c(0L, 1L, 0L, 0L, 2L, 0L, 0L, 0L, 0L, 0L))
    expect_identical(se3$nucl$s1[[2]], segs$s1[[1]])
    expect_identical(se3$nucl$s1[[5]], segs$s1[[2]])
})
