test_that("read regrouping works", {
    # get example data
    modbamfiles <- system.file("extdata", c("6mA_1_10reads.bam",
                                            "6mA_2_10reads.bam"),
                               package = "SingleMoleculeGenomicsIO")
    se <- readModBam(bamfiles = modbamfiles, regions = "chr1:6940000-6955000",
                     modbase = "a", verbose = FALSE,
                     sampleAnnot = data.frame(sample = c("s1", "s2"),
                                              group = c("g1", "g1")),
                     variantPositions = GPos(seqnames = "chr1",
                                             pos = c(6940000, 6940500)),
                     BPPARAM = BiocParallel::SerialParam())
    expect_warning(expect_warning(
        se <- addReadStats(se, name = "QC", BPPARAM = BiocParallel::SerialParam()),
        "Too few points"), "Too few points")
    segs <- list(s1 = IRanges::IRangesList(
        "s1-233e48a7-f379-4dcf-9270-958231125563" = IRanges::IRanges(
            start = 6925834, width = 140),
        "s1-92e906ae-cddb-4347-a114-bf9137761a8d" = IRanges::IRanges(
            start = c(6926000, 6926200), width = 140)),
        s2 = IRanges::IRangesList())
    se <- annotateReadSegments(se, segs, "nucl")
    # define read groups
    groups <- list(g1 = c("s1-233e48a7-f379-4dcf-9270-958231125563",
                          "s2-d03efe3b-a45b-430b-9cb6-7e5882e4faf8"),
                   g2 = "s1-92e906ae-cddb-4347-a114-bf9137761a8d",
                   g3 = c("s2-034b625e-6230-4f8d-a713-3a32cd96c298",
                          "s1-d52a5f6a-a60a-4f85-913e-eada84bfbfb9"))
    expectedOrder <- match(
        unlist(groups, use.names = FALSE),
        unlist(lapply(assay(se, "mod_prob"), colnames), use.names = FALSE))

    # test that functions fail with wrong input
    expect_error(regroupReads(se = "1", readGroups = groups),
                 "must be of class .RangedSummarizedExperiment.")
    expect_error(regroupReads(se = se, readGroups = 1),
                 "must be of class .list.")
    expect_error(regroupReads(se = se, readGroups = unname(groups)),
                 "must not be")
    expect_error(regroupReadsByColData(se = se, colNames = 1, withinSample = TRUE),
                 "must be of class .character.")
    expect_error(regroupReadsByColData(se = se, colNames = "missing", withinSample = TRUE),
                 "must be one of")
    expect_error(regroupReadsByColData(se = se, colNames = "variant_label",
                                       withinSample = 1),
                 "must be of class .logical.")
    expect_error(regroupReadsByColData(se = se, colNames = "variant_label",
                                       withinSample = c(TRUE, FALSE)),
                 "must have length 1")
    expect_error(regroupReadsByColData(se = se, colNames = "QC"),
                 "is not atomic and can not be used for read regrouping")

    # regroup reads based on predefined grouping
    sere <- regroupReads(se, readGroups = groups)
    expect_identical(lapply(assay(sere, "mod_prob"), ncol),
                     list(g1 = 2L, g2 = 1L, g3 = 2L))
    expect_identical(colnames(as.matrix(assay(sere, "mod_prob"))),
                     paste0(rep(names(groups), lengths(groups)), "-",
                            unlist(groups)))
    expect_identical(assayNames(sere), "mod_prob")
    expect_identical(unname(as.matrix(assay(sere, "mod_prob"))),
                     unname(as.matrix(assay(se, "mod_prob"))[, expectedOrder]))
    tmp <- do.call(rbind, sere$readInfo)
    rownames(tmp) <- sub("^g[0-9]-", "", rownames(tmp))
    expect_identical(tmp,
                     do.call(rbind, se$readInfo)[expectedOrder, ])
    tmp <- do.call(rbind, sere$QC)
    rownames(tmp) <- sub("^g[0-9]-", "", rownames(tmp))
    expect_identical(tmp,
                     do.call(rbind, se$QC)[expectedOrder, ])
    expect_identical(colnames(sere), names(groups))
    expect_identical(rowRanges(se), rowRanges(sere))
    tmp <- do.call(c, unname(se$nucl))[expectedOrder]
    names(tmp) <- paste0("g", c(1, 1, 2, 3, 3), "-", names(tmp))
    expect_identical(tmp, do.call(c, unname(sere$nucl)))
    expect_identical(lengths(sere$nucl), c(g1 = 2L, g2 = 1L, g3 = 2L))

    # ... identical results (with warning) if nonexistent reads are provided
    groups2 <- groups
    groups2$g1 <- c(groups2$g1, "missing1")
    groups2$g3 <- c(groups2$g3, "missing2")
    expect_warning({
        sere2 <- regroupReads(se, readGroups = groups2)
    }, "The following reads were not found")
    expect_identical(sere, sere2)

    # ... identical results also if non-read-level assays are present
    se2 <- flattenReadLevelAssay(se)
    sere2 <- regroupReads(se2, readGroups = groups)
    expect_identical(sere, sere2)

    # ... fail if there are no read-level assays
    se2 <- flattenReadLevelAssay(se, keepReads = FALSE)
    expect_error(regroupReads(se2, readGroups = groups),
                 "does not contain any read-level assays")

    # empty sample
    expect_warning({
        sere2 <- regroupReads(se, readGroups = c(groups, list(g4 = "missing")))
    }, "The following reads were not found")
    expect_identical(sere, sere2)

    # inconsistent modbase
    se2 <- se
    se2$modbase[2] <- "m"
    expect_error(regroupReads(se2, readGroups = groups),
                 "Some read groups correspond to reads with different modbases")

    # regroup by read annotation (variant label), across samples
    sere <- regroupReadsByColData(se, colNames = "variant_label",
                                  withinSample = FALSE)
    groups2 <- split(x = rownames(do.call(rbind, se$readInfo)),
                     f = do.call(rbind, se$readInfo)$variant_label)
    expectedOrder2 <- match(
        unlist(groups2, use.names = FALSE),
        unlist(lapply(assay(se, "mod_prob"), colnames), use.names = FALSE))
    expect_identical(lapply(assay(sere, "mod_prob"), ncol),
                     list(`G-` = 3L, GT = 2L))
    expect_identical(colnames(as.matrix(assay(sere, "mod_prob"))),
                     paste0(rep(names(groups2), lengths(groups2)), "-",
                            unlist(groups2)))
    expect_identical(assayNames(sere), "mod_prob")
    expect_identical(unname(as.matrix(assay(sere, "mod_prob"))),
                     unname(as.matrix(assay(se, "mod_prob"))[, expectedOrder2]))
    expect_identical(colnames(sere), names(groups2))
    expect_identical(rowRanges(se), rowRanges(sere))

    # ... within sample
    sere <- regroupReadsByColData(se, colNames = "variant_label",
                                  withinSample = TRUE)
    groups2 <- split(x = rownames(do.call(rbind, se$readInfo)),
                     f = paste0(rep(colnames(se), se$n_reads), "-",
                                do.call(rbind, se$readInfo)$variant_label))
    expectedOrder2 <- match(
        unlist(groups2, use.names = FALSE),
        unlist(lapply(assay(se, "mod_prob"), colnames), use.names = FALSE))
    expect_identical(lapply(assay(sere, "mod_prob"), ncol),
                     list(`s1-G-` = 1L, `s1-GT` = 2L, `s2-G-` = 2L))
    expect_identical(colnames(as.matrix(assay(sere, "mod_prob"))),
                     paste0(rep(names(groups2), lengths(groups2)), "-",
                            unlist(groups2)))
    expect_identical(assayNames(sere), "mod_prob")
    expect_identical(unname(as.matrix(assay(sere, "mod_prob"))),
                     unname(as.matrix(assay(se, "mod_prob"))[, expectedOrder2]))
    expect_identical(colnames(sere), names(groups2))
    expect_identical(rowRanges(se), rowRanges(sere))

    # multiple annotation columns
    se2 <- se
    se2$readInfo <- lapply(se2$readInfo, function(ri) {
        ri$label2 <- ri$variant_label
        ri
    })
    # ... across samples
    sere <- regroupReadsByColData(se2, colNames = c("variant_label", "label2"),
                                  withinSample = FALSE)
    groups2 <- split(x = rownames(do.call(rbind, se2$readInfo)),
                     f = paste0(do.call(rbind, se2$readInfo)$variant_label, "-",
                                do.call(rbind, se2$readInfo)$label2))
    expectedOrder2 <- match(
        unlist(groups2, use.names = FALSE),
        unlist(lapply(assay(se2, "mod_prob"), colnames), use.names = FALSE))
    expect_identical(lapply(assay(sere, "mod_prob"), ncol),
                     list(`G--G-` = 3L, `GT-GT` = 2L))
    expect_identical(colnames(as.matrix(assay(sere, "mod_prob"))),
                     paste0(rep(names(groups2), lengths(groups2)), "-",
                            unlist(groups2)))
    expect_identical(assayNames(sere), "mod_prob")
    expect_identical(unname(as.matrix(assay(sere, "mod_prob"))),
                     unname(as.matrix(assay(se2, "mod_prob"))[, expectedOrder2]))
    expect_identical(colnames(sere), names(groups2))
    expect_identical(rowRanges(se2), rowRanges(sere))

    # ... within sample
    sere <- regroupReadsByColData(se2, colNames = c("variant_label", "label2"),
                                  withinSample = TRUE)
    groups2 <- split(x = rownames(do.call(rbind, se$readInfo)),
                     f = paste0(rep(colnames(se), se2$n_reads), "-",
                                do.call(rbind, se2$readInfo)$variant_label, "-",
                                do.call(rbind, se2$readInfo)$label2))
    expectedOrder2 <- match(
        unlist(groups2, use.names = FALSE),
        unlist(lapply(assay(se2, "mod_prob"), colnames), use.names = FALSE))
    expect_identical(lapply(assay(sere, "mod_prob"), ncol),
                     list(`s1-G--G-` = 1L, `s1-GT-GT` = 2L, `s2-G--G-` = 2L))
    expect_identical(colnames(as.matrix(assay(sere, "mod_prob"))),
                     paste0(rep(names(groups2), lengths(groups2)), "-",
                            unlist(groups2)))
    expect_identical(assayNames(sere), "mod_prob")
    expect_identical(unname(as.matrix(assay(sere, "mod_prob"))),
                     unname(as.matrix(assay(se2, "mod_prob"))[, expectedOrder2]))
    expect_identical(colnames(sere), names(groups2))
    expect_identical(rowRanges(se2), rowRanges(sere))

    # colData column
    sere <- regroupReadsByColData(se, colNames = "group",
                                  withinSample = FALSE)
    expect_identical(dim(sere), c(nrow(se), 1L))
    expect_identical(dim(assay(sere, "mod_prob")[[1]]), c(nrow(se), 5L))
    expect_identical(colnames(sere), "g1")
    expect_identical(colnames(assay(sere, "mod_prob")[[1]]),
                     paste0("g1-", colnames(as.matrix(assay(se, "mod_prob")))))

    # ... works also if there is no modbase column
    setmp <- se
    setmp$modbase <- NULL
    sere2 <- regroupReadsByColData(setmp, colNames = "group",
                                   withinSample = FALSE)
    expect_identical(dim(sere2), c(nrow(setmp), 1L))
    expect_identical(dim(assay(sere2, "mod_prob")[[1]]), c(nrow(setmp), 5L))
    expect_identical(colnames(sere2), "g1")
    expect_identical(colnames(assay(sere2, "mod_prob")[[1]]),
                     paste0("g1-", colnames(as.matrix(assay(setmp, "mod_prob")))))
    expect_null(sere2$modbase)
    expect_identical(assays(sere), assays(sere2))

    # group + readInfo column
    sere1 <- regroupReadsByColData(se, colNames = "variant_label",
                                   withinSample = FALSE)
    sere2 <- regroupReadsByColData(se, colNames = c("variant_label", "group"),
                                   withinSample = FALSE)
    expect_identical(dim(sere1), dim(sere2))
    expect_identical(colnames(sere1), c("G-", "GT"))
    expect_identical(colnames(sere2), c("G--g1", "GT-g1"))
    expect_identical(nnavals(assay(sere1, "mod_prob")[[1]]),
                     nnavals(assay(sere2, "mod_prob")[[1]]))
})
