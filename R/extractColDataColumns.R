#' Extract one or more (possibly nested) columns from the colData of a SummarizedExperiment object
#'
#' @param se A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object.
#' @param colNames A character vector corresponding to the names of annotation
#'     (\code{colData}) columns to extract and consolidate into a
#'     \code{data.frame}. The names can be either columns of
#'     \code{colData(se)} itself, or columns in a nested (read-level)
#'     annotation column of \code{colData(se)}. In the latter case, the
#'     \code{colNames} should be of the form \code{outerColName:innerColName},
#'     where \code{outerColName} is a column name in \code{colData(se)},
#'     and \code{innerColName} is a column name in each element of
#'     \code{colData(se)[[outerColName]]}.
#' @param alignWith A character scalar indicating whether the rows in the
#'     returned \code{data.frame} correspond to reads (\code{"read"}) or
#'     samples (\code{"sample"}). If any element of \code{colNames} refers to
#'     a nested column, \code{alignWith} must be \code{"read"}.
#'
#' @returns
#' A \code{data.frame} with columns \code{sample}, \code{read} (if
#' \code{alignWith = "read"}) and \code{make.names(colNames)}.
#'
#' @examples
#' extractfiles <- system.file("extdata",
#'                             c("modkit_extract_rc_6mA_1.tsv.gz",
#'                               "modkit_extract_rc_6mA_2.tsv.gz"),
#'                             package = "SingleMoleculeGenomicsIO")
#' se <- readModkitExtract(extractfiles, modbase = "a", filter = "modkit",
#'                         BPPARAM = BiocParallel::SerialParam())
#' df <- extractColDataColumns(se = se,
#'                             colNames = c("modbase", "readInfo:ref_strand"),
#'                             alignWith = "read")
#' head(df)
#' dim(df)
#' df <- extractColDataColumns(se = se,
#'                             colNames = "modbase",
#'                             alignWith = "sample")
#' head(df)
#' dim(df)
#'
#' @export
#' @author Michael Stadler, Charlotte Soneson
#'
#' @importFrom cli cli_abort
#' @importFrom SummarizedExperiment colData
extractColDataColumns <- function(
        se, colNames,
        alignWith = ifelse(any(grepl(":", colNames)), "read", "sample")) {
    .assertVector(x = se, type = "SummarizedExperiment")
    checkSEValidity(se, verbose = FALSE)
    .assertVector(x = colNames, type = "character")
    .assertScalar(x = alignWith, type = "character", validValues = c("sample", "read"))
    if (identical(alignWith, "sample") && any(grepl(":", colNames))) {
        cli_abort("Setting {.arg alignWith} to {alignWith} is not supported if any value in {.arg colNames} contains ':'")
    }

    if (alignWith == "read") {
        readsBySample <- getReadNamesBySample(se)
        nReadsPerSample <- lengths(readsBySample)
        readInfo <- data.frame(
            sample = rep(names(readsBySample), nReadsPerSample),
            read = unlist(readsBySample)
        )
    } else {
        readInfo <- data.frame(
            sample = se$sample
        )
        nReadsPerSample <- 1
    }

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

    for (i in seq_along(colNamesNested)) {
        if (!is.atomic(se[[names(colNamesNested)[i]]][[1]][[colNamesNested[i]]])) {
            cli_abort("The {.var {names(colNamesNested)[i]}:{colNamesNested[i]}} column is not atomic")
        }
        tmp <- unlist(lapply(se[[names(colNamesNested)[i]]], "[[", colNamesNested[i]))
        readInfo[[make.names(paste0(names(colNamesNested)[i], ":", colNamesNested[i]))]] <- tmp
    }
    for (m in colNamesTopLevel) {
        if (!is.atomic(se[[m]])) {
            cli_abort("The {.var {m}} column is not atomic")
        }
        readInfo[[make.names(m)]] <- rep(se[[m]], nReadsPerSample)
    }

    readInfo
}
