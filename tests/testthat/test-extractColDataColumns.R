test_that("extractColDataColumns works", {
    extractfiles <- system.file("extdata",
                                c("modkit_extract_rc_6mA_1.tsv.gz",
                                  "modkit_extract_rc_6mA_2.tsv.gz"),
                                package = "SingleMoleculeGenomicsIO")
    se <- readModkitExtract(extractfiles, modbase = "a", filter = "modkit",
                            BPPARAM = BiocParallel::SerialParam())
    suppressWarnings({
        se <- addReadStats(se, name = "QC")
    })

    expect_error(extractColDataColumns(se = 1, colNames = "modbase"),
                 ".se. must be of class .SummarizedExperiment.")
    expect_error(extractColDataColumns(se = se, colNames = 1),
                 ".colNames. must be of class .character.")
    expect_error(extractColDataColumns(se = se, colNames = "modbase",
                                       alignWith = 1),
                 ".alignWith. must be of class .character.")
    expect_error(extractColDataColumns(se = se, colNames = "modbase",
                                       alignWith = c("read", "sample")),
                 ".alignWith. must have length 1")
    expect_error(extractColDataColumns(se = se, colNames = "modbase",
                                       alignWith = "missing"),
                 "must be one of")
    expect_error(extractColDataColumns(se = se, colNames = "missing:ref_strand"),
                 "must be one of")
    expect_error(extractColDataColumns(se = se, colNames = "missing"),
                 "must be one of")
    expect_error(extractColDataColumns(se = se, colNames = "readInfo:ref_strand",
                                       alignWith = "sample"),
                 "is not supported if any value in .colNames.")
    expect_error(extractColDataColumns(se = se, colNames = "QC:PACModProb"),
                 "is not atomic")
    expect_error(extractColDataColumns(se = se, colNames = "QC"),
                 "is not atomic")

    df <- extractColDataColumns(se = se,
                                colNames = c("modbase", "readInfo:ref_strand"),
                                alignWith = "read")
    expect_s3_class(df, "data.frame")
    expect_named(df, c("sample", "read", "readInfo.ref_strand", "modbase"))
    expect_identical(dim(df), c(20L, 4L))
    expect_identical(df$readInfo.ref_strand,
                     do.call(rbind, se$readInfo)[["ref_strand"]])
    expect_identical(df$sample,
                     rep(colnames(se),
                         vapply(SummarizedExperiment::assay(se, "mod_prob"), ncol, 0L)))
    expect_identical(df$modbase,
                     unname(
                         rep(se$modbase,
                             vapply(SummarizedExperiment::assay(se, "mod_prob"), ncol, 0L))))
    expect_identical(df$read,
                     unname(unlist(lapply(SummarizedExperiment::assay(se, "mod_prob"),
                                          colnames))))

    df <- extractColDataColumns(se = se,
                                colNames = "modbase",
                                alignWith = "sample")
    expect_s3_class(df, "data.frame")
    expect_named(df, c("sample", "modbase"))
    expect_identical(dim(df), c(2L, 2L))
    expect_identical(df$sample, colnames(se))
    expect_identical(df$modbase, unname(se$modbase))
})
