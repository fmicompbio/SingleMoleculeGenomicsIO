#' Get names of assays containing read-level data
#'
#' The names of assays designated as containing read-level data are
#' extracted from \code{metadata(se)$readLevelData$assayNames}.
#'
#' @export
#'
#' @param se A \code{SummarizedExperiment} object.
#'
#' @author Charlotte Soneson
#'
#' @return A (possibly empty) character vector with the names of the assays of
#' se containing read-level data.
#'
#' @examples
#' modbamfiles <- system.file("extdata", c("6mA_1_10reads.bam",
#'                                         "6mA_2_10reads.bam"),
#'                            package = "SingleMoleculeGenomicsIO")
#' se <- readModBam(bamfiles = modbamfiles, regions = "chr1:6940000-6955000",
#'                  modbase = "a", verbose = FALSE,
#'                  BPPARAM = BiocParallel::SerialParam())
#' se <- addReadStats(se, BPPARAM = BiocParallel::SerialParam())
#' getReadLevelAssayNames(se)
#'
#' @importFrom SummarizedExperiment assayNames
getReadLevelAssayNames <- function(se) {
    .assertVector(x = se, type = "SummarizedExperiment")
    intersect(metadata(se)$readLevelData$assayNames,
              assayNames(se))
}

#' Get names of colData columns containing read-level data
#'
#' The names of \code{colData} column designated as containing read-level
#' annotations are extracted from \code{metadata(se)$readLevelData$colDataColumns}.
#'
#' @export
#'
#' @param se A \code{SummarizedExperiment} object.
#'
#' @author Charlotte Soneson
#'
#' @return A (possibly empty) character vector with the names of the columns of
#' colData(se) containing read-level data.
#'
#' @examples
#' modbamfiles <- system.file("extdata", c("6mA_1_10reads.bam",
#'                                         "6mA_2_10reads.bam"),
#'                            package = "SingleMoleculeGenomicsIO")
#' se <- readModBam(bamfiles = modbamfiles, regions = "chr1:6940000-6955000",
#'                  modbase = "a", verbose = FALSE,
#'                  BPPARAM = BiocParallel::SerialParam())
#' se <- addReadStats(se, BPPARAM = BiocParallel::SerialParam())
#' getReadLevelColDataNames(se)
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom BiocGenerics colnames
#' @importFrom S4Vectors metadata
getReadLevelColDataNames <- function(se) {
    .assertVector(x = se, type = "SummarizedExperiment")
    intersect(metadata(se)$readLevelData$colDataColumns,
              colnames(colData(se)))
}

#' Get list of read names by sample from SummarizedExperiment object
#'
#' @param se A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object.
#'
#' @export
#' @author Charlotte Soneson
#'
#' @returns A named list with one entry per sample, containing the read names
#' for the respective sample.
#'
#' @examples
#' extractfiles <- system.file("extdata",
#'                             c("modkit_extract_rc_6mA_1.tsv.gz",
#'                               "modkit_extract_rc_6mA_2.tsv.gz"),
#'                             package = "SingleMoleculeGenomicsIO")
#' se <- readModkitExtract(extractfiles, modbase = "a", filter = "modkit",
#'                         BPPARAM = BiocParallel::SerialParam())
#' getReadNamesBySample(se)
#'
#' @importFrom SummarizedExperiment assay colData
#' @importFrom cli cli_abort
getReadNamesBySample <- function(se) {
    .assertVector(x = se, type = "SummarizedExperiment")
    rlAssays <- getReadLevelAssayNames(se)
    rlColNames <- getReadLevelColDataNames(se)
    if (length(rlAssays) > 0) {
        readsBySample <- lapply(assay(se, rlAssays[1]), colnames)
        attr(readsBySample, "source") <- paste0("assay ", rlAssays[1])
    } else if (length(rlColNames) > 0) {
        readsBySample <- lapply(colData(se)[[rlColNames[1]]], rownames)
        attr(readsBySample, "source") <- paste0("colData column ", rlColNames[1])
    } else {
        cli_abort("{.arg se} does not contain any read-level assays or colData columns")
    }
    return(readsBySample)
}

#' Check internal consistency of SummarizedExperiment object
#'
#' All assays with read-level data must have the same number and order of
#' the reads, which must also agree with the order in \code{se$QC} if that
#' exists. All assays must have the same column names, which must also
#' agree with the column names of the object, and the \code{sample} column
#' in the \code{colData}.
#'
#' @export
#'
#' @param se A \code{SummarizedExperiment object}.
#' @param verbose A logical scalar. If \code{TRUE}, report on progress.
#'
#' @author Charlotte Soneson
#'
#' @return Silently returns \code{NULL}. If the object is not valid, an error
#' will be raised.
#'
#' @importFrom SummarizedExperiment colData assayNames assay
#' @importFrom BiocGenerics nrow
#' @importFrom cli cli_abort
#' @importFrom S4Vectors metadata
#'
#' @examples
#' library(GenomicRanges)
#' modbamfiles <- system.file("extdata", c("6mA_1_10reads.bam",
#'                                         "6mA_2_10reads.bam"),
#'                            package = "SingleMoleculeGenomicsIO")
#' se <- readModBam(bamfiles = modbamfiles, regions = "chr1:6940000-6955000",
#'                  modbase = "a", verbose = FALSE,
#'                  variantPositions = GPos(seqnames = "chr1",
#'                                          pos = c(6940000, 6940500)),
#'                  BPPARAM = BiocParallel::SerialParam())
#' checkSEValidity(se)
#'
checkSEValidity <- function(se, verbose = FALSE) {
    stopifnot(is(se, "SummarizedExperiment"))

    .message("Checking assay names")
    stopifnot(!is.null(assayNames(se)) &&
                  all(assayNames(se) != "") &&
                  anyDuplicated(assayNames(se)) == 0L)

    if (nrow(se) > 0) {
        .message("Checking row names")
        if (!is.null(rownames(se))) {
            stopifnot(anyDuplicated(rownames(se)) == 0L)
        }
    }

    stopifnot(!is.null(metadata(se)$readLevelData) &&
                  is.list(metadata(se)$readLevelData) &&
                  all(c("assayNames", "colDataColumns") %in%
                          names(metadata(se)$readLevelData)))

    .message("Checking consistency of sample names")
    stopifnot("sample" %in% colnames(colData(se)))

    for (an in assayNames(se)) {
        stopifnot(colnames(assay(
            se, an, withDimnames = FALSE)) == colnames(se))
    }
    for (cn in getReadLevelColDataNames(se)) {
        stopifnot(names(se[[cn]]) == colnames(se))
    }
    if (any(grepl(":", colnames(colData(se)), fixed = TRUE))) {
        cli_abort("Column names in {.var colData(se)} can not contain ':'")
    }

    rlAssays <- getReadLevelAssayNames(se)
    if (length(rlAssays) > 0) {
        .message("Read-level assay found")
        ## Choose one assay as the reference to compare to
        refAssay <- rlAssays[1]
        refReads <- lapply(assay(se, refAssay), colnames)
        for (an in setdiff(rlAssays, refAssay)) {
            .message("Comparing {refAssay} and {an}")
            for (sn in colnames(se)) {
                if (!all(colnames(assay(se, an)[[sn]]) == refReads[[sn]])) {
                    cli_abort(paste0(
                        "Mismatching reads for assays {refAssay} and {an}, ",
                        "sample {sn}"))
                }
            }
        }
        for (cn in getReadLevelColDataNames(se)) {
            .message("Read-level column data found, checking consistency")
            for (sn in colnames(se)) {
                if (!all(rownames(se[[cn]][[sn]]) == refReads[[sn]])) {
                    cli_abort(paste0(
                        "Mismatching reads for assay {refAssay} and ",
                         "colData column {cn}, sample {sn}"))
                }
            }
        }
    }
}
