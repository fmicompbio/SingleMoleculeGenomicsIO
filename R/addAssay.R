#' Add a read-level assay to a SummarizedExperiment object
#'
#' @param se A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object.
#' @param assayDF A \code{\link[S4Vectors]{DataFrame}} with the same
#'     dimensionality as \code{se}. Each column should be a read-level
#'     \code{\link[SparseArray]{NaMatrix}} for a single sample.
#' @param assayName A character scalar giving the name of the new assay.
#' @param replaceExisting A logical scalar indicating whether to replace an
#'     existing assay with the same name or not.
#' @param verbose Logical scalar. If \code{TRUE}, report on progress.
#'
#' @export
#' @author Charlotte Soneson
#'
#' @returns
#' A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object with the
#' new read-level assay added.
#'
#' @examples
#' library(SummarizedExperiment)
#' library(S4Vectors)
#' modbamfile <- system.file("extdata", "6mA_1_10reads.bam",
#'                           package = "SingleMoleculeGenomicsIO")
#' se <- readModBam(bamfiles = modbamfile, regions = "chr1:6940000-6955000",
#'                  modbase = "a", verbose = TRUE,
#'                  BPPARAM = BiocParallel::SerialParam())
#' # duplicate the 'mod_prob' assay into a new assay named 'new_assay'
#' se <- addReadLevelAssay(se, assayDF = assay(se, "mod_prob"),
#'                         assayName = "new_assay")
#' assayNames(se)
#' metadata(se)$readLevelData
#'
#' @importFrom SummarizedExperiment assayNames assay
#' @importFrom cli cli_abort
#' @importFrom S4Vectors metadata
addReadLevelAssay <- function(se, assayDF, assayName,
                              replaceExisting = FALSE, verbose = FALSE) {
    .assertVector(x = se, type = "SummarizedExperiment")
    .assertVector(x = assayDF, type = "DataFrame")
    .assertScalar(x = assayName, type = "character")
    .assertScalar(x = replaceExisting, type = "logical")
    .assertScalar(x = verbose, type = "logical")

    if (!replaceExisting && assayName %in% assayNames(se)) {
        cli_abort("{.var se} already has an assay named {assayName}. Set {.arg replaceExisting} to TRUE if you want to replace it.")
    }
    if (!all(dim(se) == dim(assayDF))) {
        cli_abort("The dimensions of {.var se} ({dim(se)}) are not the same as those of {.var assayDF} ({dim(assayDF)})")
    }

    suppressWarnings({
        # currently, assigning to assays triggers a deprecation warning
        # (introduced in https://github.com/Bioconductor/IRanges/commit/b4e9e7e8530a822980259c37cef186c652ba8be5)
        # see issue at https://github.com/Bioconductor/SummarizedExperiment/issues/74
        assay(se, assayName) <- assayDF
    })
    metadata(se)$readLevelData$assayNames <-
        union(metadata(se)$readLevelData$assayNames, assayName)
    checkSEValidity(se, verbose = verbose)
    se
}
