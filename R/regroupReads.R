#' Regroup reads
#'
#' Regroup reads into new "samples", defined by \code{readGroups}.
#'
#' @param se A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object,
#'     e.g. obtained by \code{\link{readModBam}}. Rows represent positions, and
#'     columns samples. The object should have at least one read-level assay.
#' @param readGroups A named list, where each element corresponds to a new
#'     group of reads (a "sample").
#'
#' @returns A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'     with columns corresponding to the names of \code{readGroups}. Only the
#'     read-level assays will be returned (as summary-level assays are likely
#'     no longer consistent with the regrouped read-level assays).
#'
#' @author Charlotte Soneson
#' @export
#' @name regroupReads
#'
#' @examples
#' library(SummarizedExperiment)
#' library(GenomicRanges)
#' modbamfiles <- system.file("extdata", c("6mA_1_10reads.bam",
#'                                         "6mA_2_10reads.bam"),
#'                            package = "SingleMoleculeGenomicsIO")
#' se <- readModBam(bamfiles = modbamfiles, regions = "chr1:6940000-6955000",
#'                  modbase = "a", verbose = TRUE,
#'                  variantPositions = GPos(seqnames = "chr1",
#'                                          pos = c(6940000, 6940500)),
#'                  BPPARAM = BiocParallel::SerialParam())
#' # number of reads per sample
#' lapply(assay(se, "mod_prob"), ncol)
#' # define groups
#' groups <- list(g1 = c("s1-233e48a7-f379-4dcf-9270-958231125563",
#'                       "s2-d03efe3b-a45b-430b-9cb6-7e5882e4faf8"),
#'                g2 = "s1-92e906ae-cddb-4347-a114-bf9137761a8d",
#'                g3 = c("s2-034b625e-6230-4f8d-a713-3a32cd96c298",
#'                       "s1-d52a5f6a-a60a-4f85-913e-eada84bfbfb9"))
#' # regroup reads
#' sere <- regroupReads(se, readGroups = groups)
#' lapply(assay(sere, "mod_prob"), ncol)
#'
#' # regroup by read annotation (variant label)
#' # ... across samples
#' sere <- regroupReadsByColData(se, colNames = "readInfo:variant_label",
#'                               withinSample = FALSE)
#' sere
#' lapply(assay(sere, "mod_prob"), ncol)
#' # ... within sample
#' sere <- regroupReadsByColData(se, colNames = "readInfo:variant_label",
#'                               withinSample = TRUE)
#' sere
#' lapply(assay(sere, "mod_prob"), ncol)
#'
#' @importFrom S4Vectors metadata make_zero_col_DFrame
#' @importFrom SummarizedExperiment assayNames assay colData
#' @importFrom cli cli_abort cli_warn
#' @importFrom stats setNames
#'
regroupReads <- function(se, readGroups) {
    .assertVector(x = se, type = "RangedSummarizedExperiment")
    .assertVector(x = readGroups, type = "list")
    .assertVector(x = names(readGroups), type = "character")
    checkSEValidity(se = se, verbose = FALSE)

    # list read-level assays and colData columns
    rlAssays <- intersect(metadata(se)$readLevelData$assayNames, assayNames(se))
    if (length(rlAssays) == 0) {
        cli_abort("{.arg se} does not contain any read-level assays")
    }
    rlCols <- intersect(metadata(se)$readLevelData$colDataColumns,
                        colnames(colData(se)))

    # exclude any reads that are not found in the SE
    seReads <- colnames(as.matrix(assay(se, rlAssays[1])))
    msng <- setdiff(unlist(readGroups), seReads)
    if (length(msng) > 0) {
        cli_warn(paste0("The following reads were not found in {.arg se} ",
                        "and will be ignored: {msng}"))
        readGroups <- lapply(readGroups, function(rg) intersect(rg, seReads))
    }
    # exclude groups without reads
    readGroups <- readGroups[lengths(readGroups) > 0]

    # check that the modbase is consistent for each read group if it exists
    if (!is.null(se$modbase)) {
        mbmap <- setNames(
            rep(se$modbase, vapply(assay(se, rlAssays[1]), ncol, 0L)),
            unlist(lapply(assay(se, rlAssays[1]), colnames))
        )
        modbase <- lapply(readGroups, function(rg) unique(mbmap[rg]))
        if (any(lengths(modbase) > 1)) {
            cli_abort("Some read groups correspond to reads with different modbases")
        }
        cdata <- DataFrame(sample = names(readGroups),
                           modbase = unlist(modbase),
                           n_reads = lengths(readGroups))
    } else {
        cdata <- DataFrame(sample = names(readGroups),
                           n_reads = lengths(readGroups))
    }

    # generate regrouped assays
    aList <- lapply(setNames(rlAssays, rlAssays), function(rla) {
        mat <- as.matrix(assay(se, rla))
        mat <- lapply(readGroups, function(rg) mat[, rg, drop = FALSE])
        df <- make_zero_col_DFrame(nrow = nrow(se))
        for (nm in names(mat)) {
            tmp <- mat[[nm]]
            colnames(tmp) <- paste0(nm, "-", colnames(tmp))
            df[[nm]] <- tmp
        }
        df
    })

    # generate regrouped colData columns
    for (rlc in rlCols) {
        if (is(colData(se)[[rlc]][[1]], "data.frame") ||
            is(colData(se)[[rlc]][[1]], "DataFrame")) {
            tmp <- do.call(rbind, colData(se)[[rlc]])
            cdata[[rlc]] <- S4Vectors::SimpleList(
                lapply(setNames(names(readGroups), names(readGroups)),
                       function(nm) {
                           tmp2 <- tmp[readGroups[[nm]], , drop = FALSE]
                           rownames(tmp2) <- paste0(nm, "-", rownames(tmp2))
                           tmp2
                       }))
        } else if (is(colData(se)[[rlc]][[1]], "IRangesList")) {
            irl <- do.call(c, unname(colData(se)[[rlc]]))
            cdata[[rlc]] <- lapply(setNames(names(readGroups), names(readGroups)),
                                   function(nm) {
                                       tmp2 <- irl[readGroups[[nm]]]
                                       names(tmp2) <- paste0(nm, "-", names(tmp2))
                                       tmp2
                                   })
        }
    }

    # generate SummarizedExperiment object
    sere <- SummarizedExperiment(
        rowRanges = rowRanges(se),
        assays = aList,
        colData = cdata,
        metadata = metadata(se)
    )

    checkSEValidity(se = sere, verbose = FALSE)
    sere
}

#' Regroup reads by annotation column
#'
#' @rdname regroupReads
#' @param colNames A character vector corresponding to the names of annotation
#'     (\code{colData}) columns, the combination of which represent the
#'     desired grouping of the reads. The names can be either columns of
#'     \code{colData(se)} itself, or columns in a nested (read-level)
#'     annotation column of \code{colData(se)}. In the latter case, the
#'     \code{colNames} should be of the form \code{outerColName:innerColName},
#'     where \code{outerColName} is a column name in \code{colData(se)},
#'     and \code{innerColName} is a column name in each element of
#'     \code{colData(se)[[outerColName]]}.
#' @param withinSample A logical scalar, indicating whether the regrouping
#'     should be done within each current sample (column of \code{se}) or not.
#'     If \code{FALSE} (default), reads are pooled across samples before
#'     being regrouped.
#'
#' @author Charlotte Soneson
#' @export
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom cli cli_abort
regroupReadsByColData <- function(se, colNames, withinSample = FALSE) {
    checkSEValidity(se, verbose = FALSE)
    rlAssays <- getReadLevelAssayNames(se)
    if (length(rlAssays) == 0) {
        cli_abort("{.arg se} does not contain any read-level assays")
    }
    .assertVector(x = colNames, type = "character")
    colNamesTopLevel <- colNames[!grepl(":", colNames, fixed = TRUE)]
    colNamesNested <- colNames[grepl(":", colNames, fixed = TRUE)]
    matches <- regexec("^([^:]+):(.+)$", colNamesNested)
    tmp <- regmatches(colNamesNested, matches)
    colNamesNested <- setNames(vapply(tmp, "[", 3, FUN.VALUE = ""),
                               vapply(tmp, "[", 2, FUN.VALUE = ""))
    for (cn in colNamesTopLevel) {
        .assertScalar(x = cn, type = "character", validValues = colnames(colData(se)))
    }
    for (i in seq_along(colNamesNested)) {
        .assertScalar(x = names(colNamesNested)[i],
                      type = "character", validValues = getReadLevelColDataNames(se))
        .assertScalar(x = colNamesNested[i], type = "character",
                      validValues = colnames(se[[names(colNamesNested)[i]]][[1]]))
    }
    .assertScalar(x = withinSample, type = "logical")

    nReadsPerSample <- vapply(assay(se, rlAssays[1]), ncol, 0L)
    readInfo <- data.frame(
        sample = rep(se$sample, nReadsPerSample),
        read_id = unlist(lapply(assay(se, rlAssays[1]), colnames))
    )
    for (i in seq_along(colNamesNested)) {
        if (!is.atomic(se[[names(colNamesNested)[i]]][[1]][[colNamesNested[i]]])) {
            cli_abort("The {.var {names(colNamesNested)[i]}:{colNamesNested[i]}} column is not atomic and can not be used for read regrouping")
        }
        tmp <- unlist(lapply(se[[names(colNamesNested)[i]]], "[[", colNamesNested[i]))
        readInfo[[paste0(names(colNamesNested)[i], "_", colNamesNested[i])]] <- tmp
    }
    # additional columns from colData(se)
    for (m in colNamesTopLevel) {
        if (!is.atomic(se[[m]])) {
            cli_abort("The {.var {m}} column is not atomic and can not be used for read regrouping")
        }
        readInfo[[m]] <- rep(se[[m]], nReadsPerSample)
    }
    # generate read groups
    if (withinSample) {
        readGroups <- split(x = readInfo$read_id,
                            f = apply(readInfo[, c("sample",
                                                   paste0(names(colNamesNested),
                                                          rep("_", length(colNamesNested)),
                                                          colNamesNested),
                                                   colNamesTopLevel), drop = FALSE],
                                      1, paste, collapse = "-"))
    } else {
        readGroups <- split(x = readInfo$read_id,
                            f = apply(readInfo[, c(paste0(names(colNamesNested),
                                                          rep("_", length(colNamesNested)),
                                                          colNamesNested),
                                                   colNamesTopLevel), drop = FALSE],
                                      1, paste, collapse = "-"))
    }
    regroupReads(se = se, readGroups = readGroups)
}
